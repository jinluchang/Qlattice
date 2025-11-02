#include <qlat/compressed-eigen-io.h>

namespace qlat
{  //

template <class T>
static void caxpy_single(ComplexT<T>* res, const ComplexT<T>& ca,
                         const ComplexT<T>* x, const ComplexT<T>* y,
                         const Int c_size)
{
  for (Int i = 0; i < c_size; i++) {
    res[i] = ca * x[i] + y[i];
  }
}

void read_floats(Vector<RealF> out, const Vector<uint8_t> fp_data)
{
  Qassert(out.data_size() == fp_data.size());
  memcpy(out.data(), fp_data.data(), fp_data.size());
  to_from_little_endian(out);
}

void read_floats_fp16(RealF* out, const uint8_t* ptr, const int64_t n,
                      const Int nsc)
{
  const int64_t nsites = n / nsc;
  Qassert(n % nsc == 0);
  const unsigned short* in = (const unsigned short*)ptr;
  // do for each site
  for (int64_t site = 0; site < nsites; site++) {
    RealF* ev = &out[site * nsc];
    const unsigned short* bptr = &in[site * (nsc + 1)];
    const unsigned short exp = bptr[0];
    const RealF max = unmap_fp16_exp(exp);
    const RealF min = -max;
    for (Int i = 0; i < nsc; i++) {
      ev[i] = fp_unmap(bptr[1 + i], min, max, USHRT_MAX);
    }
  }
}

void read_floats_fp16(Vector<RealF> out, const Vector<uint8_t> fp_data,
                      const Int nsc)
{
  const Long size = fp_16_size(out.size(), nsc);
  Qassert(fp_data.size() == size);
  std::vector<uint16_t> buffer(size / sizeof(uint16_t));
  memcpy(buffer.data(), fp_data.data(), fp_data.size());
  to_from_little_endian(get_data(buffer));
  read_floats_fp16(out.data(), (const uint8_t*)buffer.data(), out.size(), nsc);
}

CompressedEigenSystemInfo populate_eigen_system_info(
    const CompressedEigenSystemDenseInfo& cesdi,
    const std::vector<crc32_t>& crcs)
{
  CompressedEigenSystemInfo cesi;
  cesi.initialized = true;
  cesi.crcs = crcs;
  //
  cesi.total_site = cesdi.total_site;
  cesi.node_site = cesdi.node_site;
  cesi.block_site = cesdi.block_site;
  cesi.ls = cesdi.ls;
  cesi.neig = cesdi.neig;
  cesi.nkeep = cesdi.nkeep;
  cesi.nkeep_single = cesdi.nkeep_single;
  cesi.FP16_COEF_EXP_SHARE_FLOATS = cesdi.FP16_COEF_EXP_SHARE_FLOATS;
  //
  Qassert(cesi.nkeep >= cesi.nkeep_single);
  cesi.nkeep_fp16 = cesi.nkeep - cesi.nkeep_single;
  cesi.total_node = cesi.total_site / cesi.node_site;
  cesi.total_block = cesi.total_site / cesi.block_site;
  cesi.node_block = cesi.node_site / cesi.block_site;
  Qassert(cesi.total_site == cesi.total_node * cesi.node_site);
  Qassert(cesi.total_site == cesi.total_block * cesi.block_site);
  Qassert(cesi.node_site == cesi.node_block * cesi.block_site);
  Qassert((int)cesi.crcs.size() == product(cesi.total_node));
  return cesi;
}

CompressedEigenSystemInfo read_compressed_eigen_system_info(
    const std::string& root)
{
  TIMER_VERBOSE("read_compressed_eigen_system_info");
  CompressedEigenSystemDenseInfo cesdi;
  std::vector<crc32_t> crcs;
  if (get_id_node() == 0) {
    std::vector<std::string> lines = qgetlines(root + "/metadata.txt");
    for (Int i = 0; i < 4; ++i) {
      reads(cesdi.node_site[i], info_get_prop(lines, ssprintf("s[%d] = ", i)));
      reads(cesdi.block_site[i], info_get_prop(lines, ssprintf("b[%d] = ", i)));
      Int nb;
      reads(nb, info_get_prop(lines, ssprintf("nb[%d] = ", i)));
      Qassert(cesdi.node_site[i] == nb * cesdi.block_site[i]);
    }
    reads(cesdi.ls, info_get_prop(lines, ssprintf("s[4] = ")));
    Int b5, nb5;
    reads(b5, info_get_prop(lines, ssprintf("b[4] = ")));
    reads(nb5, info_get_prop(lines, ssprintf("nb[4] = ")));
    Qassert(b5 == cesdi.ls);
    Qassert(nb5 == 1);
    reads(cesdi.neig, info_get_prop(lines, "neig = "));
    reads(cesdi.nkeep, info_get_prop(lines, "nkeep = "));
    reads(cesdi.nkeep_single, info_get_prop(lines, "nkeep_single = "));
    Int blocks;
    reads(blocks, info_get_prop(lines, "blocks = "));
    Qassert(blocks == product(cesdi.node_site / cesdi.block_site));
    reads(cesdi.FP16_COEF_EXP_SHARE_FLOATS,
          info_get_prop(lines, "FP16_COEF_EXP_SHARE_FLOATS = "));
    for (Long idx = 0; true; ++idx) {
      const std::string value =
          info_get_prop(lines, ssprintf("crc32[%d] = ", idx));
      if (value == "") {
        break;
      }
      crcs.push_back(read_crc32(value));
    }
    cesdi.total_site = Coordinate();
    for (Int i = 0; i < 4; ++i) {
      const std::string str_gs = info_get_prop(lines, ssprintf("gs[%d] = ", i));
      if (str_gs != "") {
        reads(cesdi.total_site[i], str_gs);
      }
    }
    const Long global_volume = crcs.size() * product(cesdi.node_site);
    if (cesdi.total_site == Coordinate()) {
      if (64 * 64 * 64 * 128 == global_volume) {
        cesdi.total_site = Coordinate(64, 64, 64, 128);
      } else if (48 * 48 * 48 * 96 == global_volume) {
        cesdi.total_site = Coordinate(48, 48, 48, 96);
      } else if (48 * 48 * 48 * 64 == global_volume) {
        cesdi.total_site = Coordinate(48, 48, 48, 64);
      } else if (24 * 24 * 24 * 64 == global_volume) {
        cesdi.total_site = Coordinate(24, 24, 24, 64);
      } else if (32 * 32 * 32 * 64 == global_volume) {
        cesdi.total_site = Coordinate(32, 32, 32, 64);
      } else if (16 * 16 * 16 * 32 == global_volume) {
        cesdi.total_site = Coordinate(16, 16, 16, 32);
      } else {
        cesdi.total_site = cesdi.node_site;
        cesdi.total_site[3] *= crcs.size();
      }
      displayln_info(0, fname + ssprintf(": using guessed total_site=%s",
                                         show(cesdi.total_site).c_str()));
    }
    Qassert(product(cesdi.total_site) == global_volume);
    const std::string str_gs5 = info_get_prop(lines, ssprintf("gs[4] = "));
    if (str_gs5 != "") {
      Int gs5;
      reads(gs5, str_gs5);
      Qassert(cesdi.ls == gs5);
    }
  }
  bcast(Vector<Char>((Char*)&cesdi, sizeof(CompressedEigenSystemDenseInfo)));
  crcs.resize(product(cesdi.total_site / cesdi.node_site));
  bcast(get_data(crcs));
  return populate_eigen_system_info(cesdi, crcs);
}

void write_compressed_eigen_system_info(const CompressedEigenSystemInfo& cesi,
                                        const std::string& root)
{
  TIMER_VERBOSE("write_compressed_eigen_system_info");
  if (get_id_node() == 0) {
    QFile fp = qfopen(root + "/metadata.txt.partial", "w");
    Qassert(not fp.null());
    for (Int i = 0; i < 4; ++i) {
      qwrite_data(ssprintf("s[%d] = %d\n", i, cesi.node_site[i]), fp);
    }
    qwrite_data(ssprintf("s[4] = %d\n", cesi.ls), fp);
    for (Int i = 0; i < 4; ++i) {
      qwrite_data(ssprintf("b[%d] = %d\n", i, cesi.block_site[i]), fp);
    }
    qwrite_data(ssprintf("b[4] = %d\n", cesi.ls), fp);
    for (Int i = 0; i < 4; ++i) {
      qwrite_data(ssprintf("nb[%d] = %d\n", i, cesi.node_block[i]), fp);
    }
    qwrite_data(ssprintf("nb[4] = 1\n"), fp);
    qwrite_data(ssprintf("neig = %d\n", cesi.neig), fp);
    qwrite_data(ssprintf("nkeep = %d\n", cesi.nkeep), fp);
    qwrite_data(ssprintf("nkeep_single = %d\n", cesi.nkeep_single), fp);
    qwrite_data(ssprintf("blocks = %d\n", product(cesi.node_block)), fp);
    qwrite_data(ssprintf("FP16_COEF_EXP_SHARE_FLOATS = %d\n",
                         cesi.FP16_COEF_EXP_SHARE_FLOATS),
                fp);
    for (Int i = 0; i < (int)cesi.crcs.size(); ++i) {
      qwrite_data(ssprintf("crc32[%d] = %08X\n", i, cesi.crcs[i]), fp);
    }
    for (Int i = 0; i < 4; ++i) {
      qwrite_data(ssprintf("gs[%d] = %d\n", i, cesi.total_site[i]), fp);
    }
    qwrite_data(ssprintf("gs[4] = %d\n", cesi.ls), fp);
    qfclose(fp);
    qrename(root + "/metadata.txt.partial", root + "/metadata.txt");
  }
}

CompressedEigenSystemInfo resize_compressed_eigen_system_info(
    const CompressedEigenSystemInfo& cesi, const Coordinate& new_size_node)
{
  CompressedEigenSystemDenseInfo cesdi;
  cesdi.total_site = cesi.total_site;
  cesdi.block_site = cesi.block_site;
  cesdi.ls = cesi.ls;
  cesdi.neig = cesi.neig;
  cesdi.nkeep = cesi.nkeep;
  cesdi.nkeep_single = cesi.nkeep_single;
  cesdi.FP16_COEF_EXP_SHARE_FLOATS = cesi.FP16_COEF_EXP_SHARE_FLOATS;
  Qassert(cesdi.total_site % new_size_node == Coordinate());
  cesdi.node_site = cesdi.total_site / new_size_node;
  return populate_eigen_system_info(
      cesdi, std::vector<crc32_t>(product(new_size_node), 0));
}

void init_half_vector(HalfVector& hv, const Geometry& geo, const Int ls)
// eo = odd
{
  TIMER("init_half_vector");
  hv.ls = ls;
  const Geometry geo_odd = geo_eo(geo_resize(geo), 1);
  const Int multiplicity = ls * HalfVector::c_size;
  hv.init();
  hv.init(geo_odd, multiplicity);
}

Geometry get_geo_from_cesi(const CompressedEigenSystemInfo& cesi,
                           const Int id_node,
                           const Coordinate& new_size_node)
{
  GeometryNode geon;
  Coordinate node_site;
  if (new_size_node == Coordinate()) {
    geon.size_node = cesi.total_site / cesi.node_site;
    node_site = cesi.node_site;
  } else {
    geon.size_node = new_size_node;
    node_site = cesi.total_site / new_size_node;
    Qassert(cesi.total_site == new_size_node * node_site);
  }
  geon.id_node = id_node;
  geon.num_node = product(geon.size_node);
  geon.coor_node = coordinate_from_index(geon.id_node, geon.size_node);
  Geometry geo_full;
  geo_full.init(geon, node_site);
  return geo_full;
}

Geometry block_geometry(const Geometry& geo_full, const Coordinate& block_site)
// new multiplicity = 1
{
  Qassert(geo_full.node_site % block_site == Coordinate());
  const Coordinate node_block = geo_full.node_site / block_site;
  Geometry geo;
  geo.init(geo_full.geon, node_block);
  return geo;
}

void init_compressed_eigen_system_data(
    CompressedEigenSystemData& cesd, const CompressedEigenSystemInfo& cesi,
    const Int id_node, const Coordinate& new_size_node)
{
  TIMER("init_compressed_eigen_system_data");
  cesd.cesi = cesi;
  cesd.block_vol_eo = product(cesi.block_site) / 2;
  cesd.ls = cesi.ls;
  cesd.basis_c_size = cesd.block_vol_eo * cesd.ls * HalfVector::c_size;
  cesd.coef_c_size = cesi.nkeep;
  cesd.basis_size_single = cesd.basis_c_size * sizeof(ComplexF);
  cesd.basis_size_fp16 = fp_16_size(cesd.basis_c_size * 2, 24);
  cesd.coef_size_single = cesi.nkeep_single * sizeof(ComplexF);
  cesd.coef_size_fp16 =
      fp_16_size(cesi.nkeep_fp16 * 2, cesi.FP16_COEF_EXP_SHARE_FLOATS);
  cesd.coef_size = cesd.coef_size_single + cesd.coef_size_fp16;
  cesd.bases_size_single = cesi.nkeep_single * cesd.basis_size_single;
  cesd.bases_size_fp16 = cesi.nkeep_fp16 * cesd.basis_size_fp16;
  cesd.bases_offset_single = 0;
  cesd.bases_offset_fp16 = cesd.bases_offset_single + cesd.bases_size_single;
  cesd.coefs_offset = cesd.bases_offset_fp16 + cesd.bases_size_fp16;
  cesd.end_offset = cesd.coefs_offset + cesi.neig * cesd.coef_size;
  const Geometry geo_full = get_geo_from_cesi(cesi, id_node, new_size_node);
  const Geometry geo = block_geometry(geo_full, cesi.block_site);
  cesd.init(geo, cesd.end_offset);
}

void init_compressed_eigen_system_bases(CompressedEigenSystemBases& cesb,
                                        const Int n_basis,
                                        const Geometry& geo_full,
                                        const Coordinate& block_site,
                                        const Int ls)
{
  if (cesb.initialized) {
    return;
  }
  TIMER("init_compressed_eigen_system_bases");
  cesb.n_basis = n_basis;
  cesb.block_vol_eo = product(block_site) / 2;
  cesb.ls = ls;
  cesb.c_size_vec = cesb.block_vol_eo * cesb.ls * HalfVector::c_size;
  const Geometry geo = block_geometry(geo_full, block_site);
  cesb.geo_full = geo_full;
  cesb.block_site = block_site;
  cesb.init();
  cesb.init(geo, n_basis * cesb.c_size_vec);
}

void init_compressed_eigen_system_bases(
    CompressedEigenSystemBases& cesb, const CompressedEigenSystemInfo& cesi,
    const Int id_node, const Coordinate& new_size_node)
{
  const Geometry geo_full = get_geo_from_cesi(cesi, id_node, new_size_node);
  init_compressed_eigen_system_bases(cesb, cesi.nkeep, geo_full,
                                     cesi.block_site, cesi.ls);
}

void init_compressed_eigen_system_coefs(CompressedEigenSystemCoefs& cesc,
                                        const Int n_vec, const Int n_basis,
                                        const Geometry& geo_full,
                                        const Coordinate& block_site,
                                        const Int ls)
{
  if (cesc.initialized) {
    return;
  }
  TIMER("init_compressed_eigen_system_coefs");
  cesc.n_vec = n_vec;
  cesc.n_basis = n_basis;
  cesc.c_size_vec = n_basis;
  const Geometry geo = block_geometry(geo_full, block_site);
  cesc.geo_full = geo_full;
  cesc.block_site = block_site;
  cesc.ls = ls;
  cesc.init();
  cesc.init(geo, n_vec * cesc.c_size_vec);
}

void init_compressed_eigen_system_coefs(
    CompressedEigenSystemCoefs& cesc, const CompressedEigenSystemInfo& cesi,
    const Int id_node, const Coordinate& new_size_node)
{
  const Geometry geo_full = get_geo_from_cesi(cesi, id_node, new_size_node);
  init_compressed_eigen_system_coefs(cesc, cesi.neig, cesi.nkeep, geo_full,
                                     cesi.block_site, cesi.ls);
}

VBFile vbopen(const std::string& fn, const std::string& mode)
{
  VBFile fp;
  fp.fn = fn;
  fp.mode = mode;
  fp.fp = qfopen(fp.fn, fp.mode);
  Qassert(not fp.fp.null());
  return fp;
}

void vbflush(VBFile& fp)
{
  if (fp.entry_total_size == 0) {
    return;
  }
  fp.buffer.resize(fp.entry_total_size);
  if (fp.mode == "r") {
    qread_data(get_data(fp.buffer), fp.fp);
    Long pos = 0;
    for (Long i = 0; i < (Long)fp.entries.size(); ++i) {
      Vector<uint8_t> d = fp.entries[i];
      memcpy(d.data(), &fp.buffer[pos], d.size());
      pos += d.size();
    }
  } else if (fp.mode == "w") {
    Long pos = 0;
    for (Long i = 0; i < (Long)fp.entries.size(); ++i) {
      const Vector<uint8_t> d = fp.entries[i];
      memcpy(&fp.buffer[pos], d.data(), d.size());
      pos += d.size();
    }
    qwrite_data(get_data(fp.buffer), fp.fp);
  } else {
    Qassert(false);
  }
  fp.buffer.clear();
  fp.entry_total_size = 0;
  fp.entries.clear();
}

void vbclose(VBFile& fp)
{
  vbflush(fp);
  qfclose(fp.fp);
}

void vbread_data(const Vector<uint8_t>& v, VBFile& fp)
{
  Qassert(fp.mode == "r");
  if (fp.entry_total_size + v.size() >= fp.buffer_limit) {
    vbflush(fp);
  }
  if (v.size() >= fp.buffer_limit) {
    qread_data(v, fp.fp);
  } else {
    fp.entry_total_size += v.size();
    fp.entries.push_back(v);
  }
}

void vbwrite_data(const Vector<uint8_t>& v, VBFile& fp)
{
  Qassert(fp.mode == "w");
  if (fp.entry_total_size + v.size() >= fp.buffer_limit) {
    vbflush(fp);
  }
  if (v.size() >= fp.buffer_limit) {
    qwrite_data(v, fp.fp);
  } else {
    fp.entry_total_size += v.size();
    fp.entries.push_back(v);
  }
}

Int vbseek(VBFile& fp, const Long offset, const Int whence)
{
  TIMER("vbseek");
  vbflush(fp);
  return qfseek(fp.fp, offset, whence);
}

void set_vfile_size(VFile& fp)
{
  TIMER("set_vfile_size");
  if (fp.size < 0) {
    QFile fpr = qfopen(fp.fn, fp.mode);
    Qassert(not fpr.null());
    qfseek(fpr, 0, SEEK_END);
    fp.size = qftell(fpr);
    qfclose(fpr);
  }
}

VFile vopen(const std::string& fn, const std::string& mode)
{
  VFile fp;
  fp.fn = fn;
  fp.mode = mode;
  fp.pos = 0;
  fp.size = -1;
  fp.read_entries.clear();
  if (fp.mode == "r") {
    fp.size = -1;
  } else if (fp.mode == "w") {
    fp.size = 0;
  } else {
    Qassert(false);
  }
  return fp;
}

void vclose(VFile& fp)
{
  Long n_seek = 0;
  if (not fp.read_entries.empty()) {
    TIMER_VERBOSE_FLOPS("vclose-read");
    Qassert(fp.mode == "r");
    set_vfile_size(fp);
    VBFile fpr = vbopen(fp.fn, "r");
    Long pos = 0;
    for (std::multimap<Long, Vector<uint8_t>>::const_iterator it =
             fp.read_entries.begin();
         it != fp.read_entries.end(); ++it) {
      const Long new_pos = it->first;
      Vector<uint8_t> data = it->second;
      if (new_pos != pos) {
        TIMER_FLOPS("vclose-fseek");
        timer.flops = std::abs(new_pos - pos);
        pos = new_pos;
        vbseek(fpr, pos, SEEK_SET);
        n_seek += 1;
      }
      vbread_data(data, fpr);
      pos += data.size();
      qassert(0 <= pos and pos <= fp.size);
      timer.flops += data.size();
    }
    vbclose(fpr);
    displayln(fname + ssprintf(": n_seek = %d.", n_seek));
  } else if (not fp.write_entries.empty()) {
    TIMER_VERBOSE_FLOPS("vclose-write");
    Qassert(fp.mode == "w");
    VBFile fpr = vbopen(fp.fn + ".partial", "w");
    Long pos = 0;
    for (std::multimap<Long, Vector<uint8_t>>::const_iterator it =
             fp.write_entries.begin();
         it != fp.write_entries.end(); ++it) {
      const Long new_pos = it->first;
      Vector<uint8_t> data = it->second;
      if (new_pos != pos) {
        TIMER_FLOPS("vclose-fseek");
        timer.flops = std::abs(new_pos - pos);
        pos = new_pos;
        vbseek(fpr, pos, SEEK_SET);
        n_seek += 1;
      }
      vbwrite_data(data, fpr);
      pos += data.size();
      timer.flops += data.size();
    }
    vbclose(fpr);
    qrename(fp.fn + ".partial", fp.fn);
    displayln(fname + ssprintf(": n_seek = %d.", n_seek));
  }
  fp.pos = 0;
  fp.read_entries.clear();
}

Int vseek(VFile& fp, const Long offset, const Int whence)
{
  if (whence == SEEK_END) {
    set_vfile_size(fp);
  }
  if (whence == SEEK_SET) {
    fp.pos = offset;
  } else if (whence == SEEK_CUR) {
    fp.pos += offset;
  } else if (whence == SEEK_END) {
    fp.pos = fp.size + offset;
  } else {
    Qassert(false);
  }
  Qassert(0 <= fp.pos);
  return 0;
}

Long vtell(const VFile& fp) { return fp.pos; }

void load_block_data(CompressedEigenSystemData& cesd, const Coordinate& xl,
                     const CompressedEigenSystemInfo& cesi, VFile& fp,
                     const Coordinate& xl_file)
{
  TIMER("load_block_data");
  Vector<uint8_t> data = cesd.get_elems(xl);
  const Int block_idx = index_from_coordinate(xl_file, cesi.node_block);
  const Int block_size = product(cesi.node_block);
  vseek(fp, 0, SEEK_END);
  const Long file_size = vtell(fp);
  Qassert(file_size == block_size * cesd.end_offset);
  const Long bases_offset_single = block_size * cesd.bases_offset_single;
  const Long bases_offset_fp16 = block_size * cesd.bases_offset_fp16;
  const Long coefs_offset = block_size * cesd.coefs_offset;
  {
    vseek(fp, bases_offset_single + block_idx * cesd.bases_size_single,
          SEEK_SET);
    Vector<uint8_t> chunk(&data[cesd.bases_offset_single],
                          cesd.bases_size_single);
    vread_data(chunk, fp);
  }
  {
    vseek(fp, bases_offset_fp16 + block_idx * cesd.bases_size_fp16, SEEK_SET);
    Vector<uint8_t> chunk(&data[cesd.bases_offset_fp16], cesd.bases_size_fp16);
    vread_data(chunk, fp);
  }
  for (Int i = 0; i < cesi.neig; ++i) {
    vseek(fp,
          coefs_offset + i * block_size * cesd.coef_size +
              block_idx * cesd.coef_size,
          SEEK_SET);
    Vector<uint8_t> chunk(&data[cesd.coefs_offset + i * cesd.coef_size],
                          cesd.coef_size);
    vread_data(chunk, fp);
  }
}

void save_block_data(const CompressedEigenSystemData& cesd,
                     const Coordinate& xl,
                     const CompressedEigenSystemInfo& cesi, VFile& fp)
{
  TIMER("load_block_data");
  const Vector<uint8_t> data = cesd.get_elems_const(xl);
  const Int block_idx = index_from_coordinate(xl, cesi.node_block);
  const Int block_size = product(cesi.node_block);
  const Long bases_offset_single = block_size * cesd.bases_offset_single;
  const Long bases_offset_fp16 = block_size * cesd.bases_offset_fp16;
  const Long coefs_offset = block_size * cesd.coefs_offset;
  {
    vseek(fp, bases_offset_single + block_idx * cesd.bases_size_single,
          SEEK_SET);
    Vector<uint8_t> chunk(&data[cesd.bases_offset_single],
                          cesd.bases_size_single);
    vwrite_data(chunk, fp);
  }
  {
    vseek(fp, bases_offset_fp16 + block_idx * cesd.bases_size_fp16, SEEK_SET);
    Vector<uint8_t> chunk(&data[cesd.bases_offset_fp16], cesd.bases_size_fp16);
    vwrite_data(chunk, fp);
  }
  for (Int i = 0; i < cesi.neig; ++i) {
    vseek(fp,
          coefs_offset + i * block_size * cesd.coef_size +
              block_idx * cesd.coef_size,
          SEEK_SET);
    Vector<uint8_t> chunk(&data[cesd.coefs_offset + i * cesd.coef_size],
                          cesd.coef_size);
    vwrite_data(chunk, fp);
  }
}

crc32_t block_data_crc(const CompressedEigenSystemData& cesd,
                       const Coordinate& xl,
                       const CompressedEigenSystemInfo& cesi,
                       const Coordinate& xl_file)
{
  TIMER("block_data_crc");
  crc32_t crc = 0;
  const Vector<uint8_t> data = cesd.get_elems_const(xl);
  const Int block_idx = index_from_coordinate(xl_file, cesi.node_block);
  const Int block_size = product(cesi.node_block);
  const Long file_size = block_size * cesd.end_offset;
  const Long bases_offset_single = block_size * cesd.bases_offset_single;
  const Long bases_offset_fp16 = block_size * cesd.bases_offset_fp16;
  const Long coefs_offset = block_size * cesd.coefs_offset;
  {
    const Vector<uint8_t> chunk(&data[cesd.bases_offset_single],
                                cesd.bases_size_single);
    crc ^= crc32_combine(
        crc32(chunk), 0,
        file_size -
            (bases_offset_single + (block_idx + 1) * cesd.bases_size_single));
  }
  {
    const Vector<uint8_t> chunk(&data[cesd.bases_offset_fp16],
                                cesd.bases_size_fp16);
    crc ^= crc32_combine(crc32(chunk), 0,
                         file_size - (bases_offset_fp16 +
                                      (block_idx + 1) * cesd.bases_size_fp16));
  }
  for (Int i = 0; i < cesi.neig; ++i) {
    const Vector<uint8_t> chunk(&data[cesd.coefs_offset + i * cesd.coef_size],
                                cesd.coef_size);
    crc ^= crc32_combine(
        crc32(chunk), 0,
        file_size - (coefs_offset + i * block_size * cesd.coef_size +
                     (block_idx + 1) * cesd.coef_size));
  }
  return crc;
}

std::vector<crc32_t> load_node_data(CompressedEigenSystemData& cesd,
                                    const CompressedEigenSystemInfo& cesi,
                                    const std::string& path)
// interface
// cesd need to be initialized beforehand (or the machine layout will be used)
{
  TIMER_VERBOSE_FLOPS("load_node_data");
  if (not cesd.initialized) {
    displayln_info(
        0,
        "initialize compressed eigen system data with current machine layout");
    init_compressed_eigen_system_data(cesd, cesi, get_id_node(),
                                      get_size_node());
  }
  std::vector<crc32_t> crcs(product(cesi.total_node), 0);
  const Geometry& geo = cesd.geo();
  const Int idx_size = product(cesi.total_node);
  std::vector<VFile> fps(idx_size);
  for (Int idx = 0; idx < idx_size; ++idx) {
    const Int dir_idx = compute_dist_file_dir_id(idx, idx_size);
    fps[idx] =
        vopen(path + ssprintf("/%02d/%010d.compressed", dir_idx, idx), "r");
  }
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const Coordinate xl_file = xg % cesi.node_block;
    const Int idx =
        index_from_coordinate(xg / cesi.node_block, cesi.total_node);
    qassert(0 <= idx && idx < idx_size);
    displayln_info(2, "load: fn='" + fps[idx].fn + "' ; index=" + show(index) +
                          "/" + show(geo.local_volume()) + " ; xl=" + show(xl) +
                          "/" + show(geo.node_site));
    load_block_data(cesd, xl, cesi, fps[idx], xl_file);
  }
  for (Int idx = 0; idx < idx_size; ++idx) {
    vclose(fps[idx]);
  }
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const Coordinate xl_file = xg % cesi.node_block;
    const Int idx =
        index_from_coordinate(xg / cesi.node_block, cesi.total_node);
    qassert(0 <= idx && idx < idx_size);
    crc32_t crc = block_data_crc(cesd, xl, cesi, xl_file);
#pragma omp critical
    crcs[idx] ^= crc;
  }
  if (cesd.geo().geon.size_node == cesi.total_node) {
    const Int id_node = cesd.geo().geon.id_node;
    if (crcs[id_node] == cesi.crcs[id_node]) {
      displayln(fname + ssprintf(": crc check successfull."));
    } else {
      displayln(
          fname +
          ssprintf(
              ": ERROR: crc check failed id_node=%d read=%08X computed=%08X.",
              id_node, cesi.crcs[id_node], crcs[id_node]));
      ssleep(1.0);
      Qassert(false);
    }
  }
  return crcs;
}

crc32_t save_node_data(const CompressedEigenSystemData& cesd,
                       const CompressedEigenSystemInfo& cesi,
                       const std::string& path)
// interface
// cesd need to be initialized beforehand (or the machine layout will be used)
{
  TIMER_VERBOSE_FLOPS("save_node_data");
  const Geometry& geo = cesd.geo();
  crc32_t crc_node = 0;
  const Int idx_size = product(geo.geon.size_node);
  const Int idx = geo.geon.id_node;
  const Int dir_idx = compute_dist_file_dir_id(idx, idx_size);
  qmkdir(path);
  qmkdir(path + ssprintf("/%02d", dir_idx));
  VFile fp =
      vopen(path + ssprintf("/%02d/%010d.compressed", dir_idx, idx), "w");
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    displayln_info(0, "save: fn='" + fp.fn + "' ; index=" + show(index) + "/" +
                          show(geo.local_volume()) + " ; xl=" + show(xl) + "/" +
                          show(geo.node_site));
    save_block_data(cesd, xl, cesi, fp);
  }
  vclose(fp);
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    crc32_t crc = block_data_crc(cesd, xl, cesi, xl);
#pragma omp critical
    crc_node ^= crc;
  }
  return crc_node;
}

void load_block(CompressedEigenSystemBases& cesb,
                CompressedEigenSystemCoefs& cesc,
                const CompressedEigenSystemData& cesd, const Coordinate& xl,
                const CompressedEigenSystemInfo& cesi)
{
  TIMER_VERBOSE("load_block");
  Vector<ComplexF> bases = cesb.get_elems(xl);
  Vector<ComplexF> coefs = cesc.get_elems(xl);
  const Vector<uint8_t> data = cesd.get_elems_const(xl);
  {
    const Long c_size = cesi.nkeep_single * cesb.c_size_vec;
    const Long d_size = c_size * sizeof(ComplexF);
    read_floats(Vector<RealF>((RealF*)&bases[0], c_size * 2),
                Vector<uint8_t>(&data[cesd.bases_offset_single], d_size));
  }
  {
    const Long c_size = cesi.nkeep_fp16 * cesb.c_size_vec;
    const Long d_size = fp_16_size(c_size * 2, 24);
    read_floats_fp16(
        Vector<RealF>((RealF*)&bases[cesi.nkeep_single * cesb.c_size_vec],
                      c_size * 2),
        Vector<uint8_t>(&data[cesd.bases_offset_fp16], d_size), 24);
  }
#pragma omp parallel for
  for (Int i = 0; i < cesb.n_basis; ++i) {
    Qassert(cesb.c_size_vec ==
            cesb.ls * cesb.block_vol_eo * HalfVector::c_size);
    Vector<ComplexF> basis(&bases[i * cesb.c_size_vec], cesb.c_size_vec);
    std::vector<ComplexF> buffer(cesb.c_size_vec);
    for (Int s = 0; s < cesb.ls; ++s) {
      for (Long index = 0; index < cesb.block_vol_eo; ++index) {
        memcpy(&buffer[(index * cesb.ls + s) * HalfVector::c_size],
               &basis[(s * cesb.block_vol_eo + index) * HalfVector::c_size],
               HalfVector::c_size * sizeof(ComplexF));
      }
    }
    assign(basis, get_data(buffer));
  }
  for (Int i = 0; i < cesc.n_vec; ++i) {
    Vector<RealF> coef((RealF*)&coefs[i * cesc.c_size_vec],
                       cesc.c_size_vec * 2);
    Vector<uint8_t> dc(&data[cesd.coefs_offset + i * cesd.coef_size],
                       cesd.coef_size);
    read_floats(Vector<RealF>(&coef[0], cesi.nkeep_single * 2),
                Vector<uint8_t>(&dc[0], cesi.nkeep_single * sizeof(ComplexF)));
    read_floats_fp16(
        Vector<RealF>(&coef[cesi.nkeep_single * 2], cesi.nkeep_fp16 * 2),
        Vector<uint8_t>(
            &dc[cesi.nkeep_single * sizeof(ComplexF)],
            fp_16_size(cesi.nkeep_fp16 * 2, cesi.FP16_COEF_EXP_SHARE_FLOATS)),
        cesi.FP16_COEF_EXP_SHARE_FLOATS);
  }
}

std::vector<crc32_t> load_node(CompressedEigenSystemBases& cesb,
                               CompressedEigenSystemCoefs& cesc,
                               const CompressedEigenSystemInfo& cesi,
                               const std::string& path)
// interface
// cesb and cesc need to be initialized beforehand (or the machine layout will
// be used)
{
  TIMER_VERBOSE("load_node");
  CompressedEigenSystemData cesd;
  if (cesb.initialized) {
    init_compressed_eigen_system_data(cesd, cesi, cesb.geo().geon.id_node,
                                      cesb.geo().geon.size_node);
  }
  std::vector<crc32_t> crcs = load_node_data(cesd, cesi, path);
  if (not cesb.initialized) {
    displayln_info(
        0,
        "initialize compressed eigen system bases and coefs with current "
        "machine layout");
    init_compressed_eigen_system_bases(cesb, cesi, get_id_node(),
                                       get_size_node());
    init_compressed_eigen_system_coefs(cesc, cesi, get_id_node(),
                                       get_size_node());
  }
  Qassert(cesb.geo() == cesc.geo());
  const Geometry& geo = cesb.geo();
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    load_block(cesb, cesc, cesd, xl, cesi);
  }
  return crcs;
}

std::vector<RealD> read_eigen_values(const std::string& path)
{
  TIMER_VERBOSE("read_eigen_values");
  Long n_eigen_values = 0;
  std::vector<RealD> vals;
  if (0 == get_id_node()) {
    const std::string filename = path + "/eigen-values.txt";
    QFile file = qfopen(filename, "r");
    Qassert(not file.null());
    n_eigen_values = read_long(qgetline(file));
    glb_sum(n_eigen_values);
    vals.resize(n_eigen_values, 0.0);
    displayln(fname + ssprintf("Reading %d eigen-values.\n", n_eigen_values));
    for (Int k = 0; k < n_eigen_values; k++) {
      vals[k] = read_double(qgetline(file));
      displayln(ssprintf("%d %24.17E", k, sqrt(vals[k])));
    }
    qfclose(file);
  } else {
    glb_sum(n_eigen_values);
    vals.resize(n_eigen_values, 0.0);
  }
  bcast(get_data(vals));
  return vals;
}

void init_blocked_half_vector(BlockedHalfVector& bhv, const Geometry& geo_full,
                              const Coordinate& block_site, const Int ls)
{
  TIMER("init_blocked_half_vector");
  bhv.block_vol_eo = product(block_site) / 2;
  bhv.ls = ls;
  const Geometry geo = block_geometry(geo_full, block_site);
  bhv.geo_full = geo_full;
  bhv.block_site = block_site;
  bhv.init();
  bhv.init(geo, bhv.block_vol_eo * ls * HalfVector::c_size);
}

void decompress_eigen_system(std::vector<BlockedHalfVector>& bhvs,
                             const CompressedEigenSystemBases& cesb,
                             const CompressedEigenSystemCoefs& cesc)
// interface
{
  TIMER_VERBOSE("decompress_eigen_system");
  const Geometry geo_full = geo_resize(cesb.geo_full);
  const Coordinate& block_site = cesb.block_site;
  const Int ls = cesb.ls;
  Qassert(cesb.geo() == cesc.geo());
  Qassert(geo_full == cesb.geo_full);
  Qassert(geo_full == cesc.geo_full);
  Qassert(block_site == cesc.block_site);
  Qassert(ls == cesc.ls);
  const Int n_vec = cesc.n_vec;
  if (n_vec == 0) {
    return;
  }
  {
    TIMER_VERBOSE("decompress_eigen_system-init");
    bhvs.resize(n_vec);
    for (Int i = 0; i < n_vec; ++i) {
      init_blocked_half_vector(bhvs[i], geo_full, block_site, ls);
    }
  }
  const Long block_size = cesb.block_vol_eo * ls * HalfVector::c_size;
  const Long n_basis = cesb.n_basis;
  Qassert(n_basis == cesc.n_basis);
  const Geometry& geo = bhvs[0].geo();
  Qassert(block_size == bhvs[0].multiplicity);
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Vector<ComplexF> bases = cesb.get_elems_const(index);
    const Vector<ComplexF> coefs = cesc.get_elems_const(index);
    std::vector<Vector<ComplexF>> blocks(n_vec);
    for (Int i = 0; i < n_vec; ++i) {
      blocks[i] = bhvs[i].get_elems(index);
    }
    TIMER_VERBOSE_FLOPS("decompress_eigen_system-block");
    timer.flops += n_vec * n_basis * block_size * 8;
#pragma omp parallel for
    for (Int i = 0; i < n_vec; ++i) {
      for (Int j = 0; j < n_basis; ++j) {
        caxpy_single(blocks[i].data(), coefs[i * n_basis + j],
                     &bases[j * block_size], blocks[i].data(), block_size);
      }
    }
  }
}

void convert_half_vector(BlockedHalfVector& bhv, const HalfVector& hv,
                         const Coordinate& block_site)
// interface
{
  TIMER("convert_half_vector");
  const Int ls = hv.ls;
  init_blocked_half_vector(bhv, geo_resize(hv.geo()), block_site, ls);
  const Geometry& geo = hv.geo();
  const Coordinate node_block = geo.node_site / block_site;
  Qassert(geo.is_only_local);
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    qassert((xl[0] + xl[1] + xl[2] + xl[3]) % 2 == 2 - geo.eo);
    const Coordinate bxl = xl / block_site;
    const Long bindex = index_from_coordinate(bxl, node_block);
    const Vector<ComplexF> site = hv.get_elems_const(index);
    qassert(site.size() == ls * HalfVector::c_size);
    const Long bidx = index_from_coordinate(xl % block_site, block_site) / 2;
    memcpy(&bhv.get_elems(bindex)[bidx * site.size()], site.data(),
           site.data_size());
  }
}

void convert_half_vector(HalfVector& hv, const BlockedHalfVector& bhv)
// interface
{
  TIMER("convert_half_vector");
  const Coordinate& block_site = bhv.block_site;
  const Int ls = bhv.ls;
  init_half_vector(hv, bhv.geo_full, ls);
  const Geometry& geo = hv.geo();
  const Coordinate node_block = geo.node_site / block_site;
  Qassert(geo.is_only_local);
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    qassert((xl[0] + xl[1] + xl[2] + xl[3]) % 2 == 2 - geo.eo);
    const Coordinate bxl = xl / block_site;
    const Long bindex = index_from_coordinate(bxl, node_block);
    Vector<ComplexF> site = hv.get_elems(index);
    qassert(site.size() == ls * HalfVector::c_size);
    const Long bidx = index_from_coordinate(xl % block_site, block_site) / 2;
    memcpy(site.data(), &bhv.get_elems_const(bindex)[bidx * site.size()],
           site.data_size());
  }
}

void convert_half_vectors(std::vector<HalfVector>& hvs,
                          std::vector<BlockedHalfVector>& bhvs)
// interface
// will clear bhvs to save space
{
  TIMER_VERBOSE("convert_half_vectors");
  clear(hvs);
  hvs.resize(bhvs.size());
  for (Int i = 0; i < (int)hvs.size(); ++i) {
    convert_half_vector(hvs[i], bhvs[i]);
    Qassert(hvs[i].geo().is_only_local);
    bhvs[i].init();
  }
  clear(bhvs);
}

void convert_half_vector_bfm_format(Vector<ComplexF> bfm_data,
                                    const HalfVector& hv)
// interface
// bfm will have t_dir simd layout
{
  TIMER("convert_half_vector_bfm_format");
  const Long size = bfm_data.size();
  Qassert(hv.geo().is_only_local);
  Qassert((Long)hv.field.size() == size);
#pragma omp parallel for
  for (Long m = 0; m < size / 2; ++m) {
    bfm_data[m * 2] = hv.field[m];
    bfm_data[m * 2 + 1] = hv.field[size / 2 + m];
  }
}

Long load_compressed_eigen_vectors(vector<RealD>& eigen_values,
                                   CompressedEigenSystemInfo& cesi,
                                   CompressedEigenSystemBases& cesb,
                                   CompressedEigenSystemCoefs& cesc,
                                   const std::string& path)
// interface
// cesb and cesc will be reinitialized
// geometry will be same as machine geometry
{
  if (!does_file_exist_sync_node(path + "/metadata.txt")) {
    displayln_info(0,
                   ssprintf("load_compressed_eigen_vectors: '%s' do not exist.",
                            path.c_str()));
    return 0;
  }
  TIMER_VERBOSE_FLOPS("load_compressed_eigen_vectors");
  Long total_bytes = 0;
  eigen_values = read_eigen_values(path);
  cesi = read_compressed_eigen_system_info(path);
  cesb.init();
  cesc.init();
  std::vector<crc32_t> crcs;
  const Int n_cycle = std::max(1, get_num_node() / dist_read_par_limit());
  {
    Long bytes = 0;
    for (Int i = 0; i < n_cycle; ++i) {
      TIMER_VERBOSE_FLOPS("load_compressed_eigen_vectors-load-cycle");
      if (get_id_node() % n_cycle == i) {
        crcs = load_node(cesb, cesc, cesi, path);
        bytes = get_data(cesb).data_size() + get_data(cesc).data_size();
      } else {
        bytes = 0;
      }
      glb_sum(bytes);
      displayln_info(
          0, fname + ssprintf(": cycle / n_cycle = %4d / %4d", i + 1, n_cycle));
      timer.flops += bytes;
      total_bytes += bytes;
    }
    timer.flops += total_bytes;
  }
  glb_sum(get_data_char(crcs));
  for (Int j = 0; j < (int)crcs.size(); ++j) {
    if (crcs[j] != cesi.crcs[j]) {
      qwarn(ssprintf("file-idx=%d loaded=%08X metadata=%08X", j, crcs[j],
                     cesi.crcs[j]));
      Qassert(false);
    }
  }
  return total_bytes;
}

crc32_t save_half_vectors(const std::vector<HalfVector>& hvs,
                          const std::string& fn,
                          const bool is_saving_crc,
                          const bool is_bfm_format)
// if is_bfm_format then will apply t_dir simd
// always big endianness
{
  TIMER_VERBOSE_FLOPS("save_half_vectors");
  Qassert(hvs.size() > 0);
  crc32_t crc = 0;
  QFile fp = qfopen(fn + ".partial", "w");
  Qassert(not fp.null());
  std::vector<ComplexF> buffer;
  for (Int i = 0; i < (int)hvs.size(); ++i) {
    TIMER_FLOPS("save_half_vectors-iter");
    const Long size = hvs[i].field.size();
    buffer.resize(size);
    timer.flops += get_data(buffer).data_size();
    if (is_bfm_format) {
      convert_half_vector_bfm_format(get_data(buffer), hvs[i]);
    }
    to_from_big_endian(get_data(buffer));  // always save data in big endianness
    crc = crc32_par(crc, get_data(buffer));
    qwrite_data(get_data(buffer), fp);
  }
  qfclose(fp);
  if (is_saving_crc) {
    qtouch(fn + ".crc32", ssprintf("%08X\n", crc));
  }
  qrename(fn + ".partial", fn);
  return crc;
}

Long decompress_eigen_vectors_node(
    const std::string& old_path, const CompressedEigenSystemInfo& cesi,
    const std::string& new_path, const Int idx,
    const Coordinate& new_size_node)
// interface
// new_size_node can be Coordinate()
// single node code
{
  TIMER_VERBOSE("decompress_eigen_vectors_node");
  Coordinate size_node = new_size_node;
  if (size_node == Coordinate()) {
    size_node = cesi.total_node;
  }
  const Int idx_size = product(size_node);
  const Int dir_idx = compute_dist_file_dir_id(idx, idx_size);
  qmkdir(new_path);
  qmkdir(new_path + ssprintf("/%02d", dir_idx));
  const std::string new_fn = new_path + ssprintf("/%02d/%010d", dir_idx, idx);
  CompressedEigenSystemBases cesb;
  CompressedEigenSystemCoefs cesc;
  init_compressed_eigen_system_bases(cesb, cesi, idx, size_node);
  init_compressed_eigen_system_coefs(cesc, cesi, idx, size_node);
  std::vector<crc32_t> crcs = load_node(cesb, cesc, cesi, old_path);
  to_from_big_endian(get_data(crcs));
  QFile fp = qfopen(new_fn + ".orig-crc32", "w");
  Qassert(not fp.null());
  qwrite_data(get_data(crcs), fp);
  qfclose(fp);
  std::vector<BlockedHalfVector> bhvs;
  decompress_eigen_system(bhvs, cesb, cesc);
  std::vector<HalfVector> hvs;
  convert_half_vectors(hvs, bhvs);
  save_half_vectors(hvs, new_fn, true, true);
  return 0;
}

crc32_t resize_compressed_eigen_vectors_node(
    std::vector<crc32_t>& crcs_acc, const std::string& old_path,
    const CompressedEigenSystemInfo& cesi, const std::string& new_path,
    const Int idx, const Coordinate& size_node)
// interface
// single node code
{
  TIMER_VERBOSE("resize_compressed_eigen_vectors_node");
  CompressedEigenSystemData cesd;
  init_compressed_eigen_system_data(cesd, cesi, idx, size_node);
  std::vector<crc32_t> crcs = load_node_data(cesd, cesi, old_path);
  if (crcs_acc.size() == 0) {
    crcs_acc.resize(crcs.size(), 0);
  } else if (crcs.size() != crcs_acc.size()) {
    Qassert(false);
  }
  for (Int i = 0; i < (int)crcs.size(); ++i) {
    crcs_acc[i] ^= crcs[i];
  }
  return save_node_data(
      cesd, resize_compressed_eigen_system_info(cesi, size_node), new_path);
}

void combine_crc32(const std::string& path, const Int idx_size,
                   const CompressedEigenSystemInfo& cesi)
{
  TIMER_VERBOSE("combine_crc32");
  if (0 == get_id_node()) {
    {
      // check orig-crc32
      const std::vector<crc32_t>& crcs_orig = cesi.crcs;
      const Int size = crcs_orig.size();
      std::vector<crc32_t> crcs_acc(size, 0);
      for (Int idx = 0; idx < idx_size; ++idx) {
        std::vector<crc32_t> crcs(size, 0);
        const Int dir_idx = compute_dist_file_dir_id(idx, idx_size);
        QFile fp = qfopen(
            path + ssprintf("/%02d/%010d.orig-crc32", dir_idx, idx), "r");
        if (not fp.null()) {
          qread_data(get_data(crcs), fp);
          to_from_big_endian(get_data(crcs));
          qfclose(fp);
        } else {
          displayln(fname +
                    ssprintf(": ERROR: %02d/%010d orig-crc32 do not exist",
                             dir_idx, idx));
          continue;
        }
        for (Int i = 0; i < size; ++i) {
          crcs_acc[i] ^= crcs[i];
        }
      }
      for (Int i = 0; i < size; ++i) {
        if (crcs_acc[i] != crcs_orig[i]) {
          displayln(fname +
                    ssprintf(": ERROR: mismatch %d/%d orig=%08X acc=%08X", i,
                             size, crcs_orig[i], crcs_acc[i]));
        } else {
          displayln(fname + ssprintf(": match %d/%d orig=%08X acc=%08X", i,
                                     size, crcs_orig[i], crcs_acc[i]));
        }
      }
    }
    const std::string fn = path + "/checksums.txt";
    if (not does_file_exist(fn)) {
      std::vector<crc32_t> crcs(idx_size, 0);
      for (Int idx = 0; idx < idx_size; ++idx) {
        const Int dir_idx = compute_dist_file_dir_id(idx, idx_size);
        const std::string path_data =
            path + ssprintf("/%02d/%010d", dir_idx, idx);
        if (does_file_exist(path_data + ".crc32")) {
          crcs[idx] = read_crc32(qcat(path_data + ".crc32"));
          displayln(
              ssprintf("reading crc %7d/%d %08X", idx, idx_size, crcs[idx]));
        } else {
          crcs[idx] = compute_crc32(path_data);
          qtouch(path_data + ".crc32", ssprintf("%08X\n", crcs[idx]));
          displayln(
              ssprintf("computing crc %7d/%d %08X", idx, idx_size, crcs[idx]));
        }
      }
      crc32_t crc = dist_crc32(crcs);
      QFile fp = qfopen(fn + ".partial", "w");
      Qassert(not fp.null());
      qwrite_data(ssprintf("%08X\n", crc), fp);
      qwrite_data("\n", fp);
      for (size_t i = 0; i < crcs.size(); ++i) {
        qwrite_data(ssprintf("%08X\n", crcs[i]), fp);
      }
      qfclose(fp);
      qrename(fn + ".partial", fn);
    }
    qtouch(path + "/checkpoint");
    // for (Int idx = 0; idx < idx_size; ++idx) {
    //   const Int dir_idx = compute_dist_file_dir_id(idx, idx_size);
    //   qremove(path + ssprintf("/%02d/%010d.crc32", dir_idx, idx));
    //   qremove(path + ssprintf("/%02d/%010d.orig-crc32", dir_idx, idx));
    // }
  }
}

void decompress_eigen_vectors(const std::string& old_path,
                              const std::string& new_path,
                              const Coordinate& new_size_node)
// interface
{
  if (does_file_exist_sync_node(new_path + "/checkpoint")) {
    return;
  }
  TIMER_VERBOSE("decompress_eigen_vectors");
  displayln_info(0, fname + ssprintf(": old_path: '") + old_path + "'");
  displayln_info(0, fname + ssprintf(": new_path: '") + new_path + "'");
  qmkdir_info(new_path);
  CompressedEigenSystemInfo cesi;
  cesi = read_compressed_eigen_system_info(old_path);
  Coordinate size_node = new_size_node;
  if (size_node == Coordinate()) {
    size_node = cesi.total_node;
  }
  Long idx_size = product(size_node);
  if (0 == get_id_node()) {
    const std::string eigen_values = qcat(old_path + "/eigen-values.txt");
    qtouch(new_path + "/eigen-values.txt", eigen_values);
  }
  displayln_info(0, fname + ssprintf(": idx_size=%d", idx_size));
  std::vector<int> avails(idx_size, 0);
  if (0 == get_id_node()) {
    for (Int idx = 0; idx < idx_size; ++idx) {
      const Int dir_idx = compute_dist_file_dir_id(idx, idx_size);
      if (!does_file_exist(new_path + ssprintf("/%02d/%010d", dir_idx, idx))) {
        avails[idx] = 1;
      }
    }
  }
  bcast(get_data(avails));
  Int n_avails = 0;
  for (Int idx = 0; idx < idx_size; ++idx) {
    if (avails[idx] != 0) {
      n_avails += 1;
    }
  }
  displayln_info(0, fname + ssprintf(": n_avails=%d", n_avails));
  if (n_avails == 0) {
    if (obtain_lock(new_path + "/lock")) {
      combine_crc32(new_path, idx_size, cesi);
      release_lock();
    }
    return;
  }
  std::vector<int> order(n_avails, -1);
  {
    TIMER_VERBOSE("decompress_eigen_vectors-shuffle");
    RngState rs(get_global_rng_state(), fname);
    split_rng_state(rs, rs, show(get_time()));
    split_rng_state(rs, rs, get_id_node());
    for (Int i = 0; i < n_avails; ++i) {
      Int idx = -1;
      do {
        idx = rand_gen(rs) % idx_size;
      } while (avails[idx] == 0);
      avails[idx] = 0;
      order[i] = idx;
    }
  }
  ssleep(2.0 * get_id_node());
  Int num_done = 0;
  for (Int i = 0; i < (int)order.size(); ++i) {
    displayln(fname + ssprintf(": %5d/%d order[%03d/%05d/%d]=%010d",
                               get_id_node(), get_num_node(), num_done, i,
                               order.size(), order[i]));
    const Int idx = order[i];
    const Int dir_idx = compute_dist_file_dir_id(idx, idx_size);
    const std::string path_data =
        new_path + ssprintf("/%02d/%010d", dir_idx, idx);
    const std::string path_lock =
        new_path + ssprintf("/%02d-%010d-lock", dir_idx, idx);
    if (!does_file_exist(path_data) and obtain_lock_all_node(path_lock)) {
      decompress_eigen_vectors_node(old_path, cesi, new_path, idx,
                                    new_size_node);
      release_lock_all_node();
      num_done += 1;
      displayln(fname + ssprintf(": %5d/%d order[%03d/%05d/%d]=%010d finished",
                                 get_id_node(), get_num_node(), num_done, i,
                                 order.size(), order[i]));
      Timer::display();
    }
  }
  displayln(
      fname +
      ssprintf(": process %5d/%d finish '", get_id_node(), get_num_node()) +
      old_path + "'");
}

bool check_compressed_eigen_vectors(const std::string& path)
// interface
{
  TIMER_VERBOSE("check_compressed_eigen_vectors");
  const CompressedEigenSystemInfo cesi =
      read_compressed_eigen_system_info(path);
  const Coordinate size_node = cesi.total_site / cesi.node_site;
  Long is_failed = 0;
  for (Int idx = 0; idx < product(size_node); ++idx) {
    if (idx % get_num_node() == get_id_node()) {
      const crc32_t crc = compute_crc32(
          path + ssprintf("/%02d/%010d.compressed",
                          compute_dist_file_dir_id(idx, product(size_node)),
                          idx));
      displayln(fname + ssprintf(": idx=%d computed=%08X previous=%08X", idx,
                                 crc, cesi.crcs[idx]));
      if (cesi.crcs[idx] != crc) {
        is_failed += 1;
        displayln(
            fname +
            ssprintf(": WARNING: mismatch idx=%d computed=%08X previous=%08X",
                     idx, crc, cesi.crcs[idx]));
        qrename(path + "/metadata.txt", path + "/metadata.txt.wrong_crc");
      }
    }
  }
  glb_sum(is_failed);
  return is_failed > 0;
}

bool resize_compressed_eigen_vectors(const std::string& old_path,
                                     const std::string& new_path,
                                     const Coordinate& size_node)
// interface
{
  TIMER_VERBOSE("resize_compressed_eigen_vectors");
  displayln_info(0, fname + ssprintf(": old_path: '") + old_path + "'");
  displayln_info(0, fname + ssprintf(": new_path: '") + new_path + "'");
  if (does_file_exist_sync_node(new_path)) {
    displayln_info(
        0, fname + ssprintf(": new_path: '%s' exists.", new_path.c_str()));
    return false;
  }
  qmkdir_info(new_path);
  CompressedEigenSystemInfo cesi;
  cesi = read_compressed_eigen_system_info(old_path);
  Long idx_size = product(size_node);
  if (0 == get_id_node()) {
    const std::string eigen_values = qcat(old_path + "/eigen-values.txt");
    qtouch(new_path + "/eigen-values.txt", eigen_values);
  }
  displayln_info(0, fname + ssprintf(": idx_size=%d", idx_size));
  const CompressedEigenSystemInfo cesi_old =
      read_compressed_eigen_system_info(old_path);
  CompressedEigenSystemInfo cesi_new =
      resize_compressed_eigen_system_info(cesi_old, size_node);
  std::vector<crc32_t> crcs_acc;
  for (Int idx = 0; idx < product(size_node); ++idx) {
    if (idx % get_num_node() == get_id_node()) {
      cesi_new.crcs[idx] = resize_compressed_eigen_vectors_node(
          crcs_acc, old_path, cesi_old, new_path, idx, size_node);
      displayln(fname + ssprintf(": resized %d/%d", idx, product(size_node)));
      Timer::display();
    }
  }
  glb_sum(get_data_char(cesi_new.crcs));
  glb_sum(get_data_char(crcs_acc));
  Qassert(crcs_acc.size() == cesi_old.crcs.size());
  for (Int j = 0; j < (int)cesi_old.crcs.size(); ++j) {
    if (crcs_acc[j] != cesi_old.crcs[j]) {
      qwarn(ssprintf("file-idx=%d loaded=%08X metadata=%08X", j, crcs_acc[j],
                     cesi.crcs[j]));
      Qassert(false);
    }
  }
  displayln_info(0, fname + ssprintf(": loaded data checksum matched."));
  write_compressed_eigen_system_info(cesi_new, new_path);
  displayln_info(0, fname + ssprintf(": checking saved data checksum"));
  return check_compressed_eigen_vectors(new_path);
}

void decompressed_eigen_vectors_check_crc32(const std::string& path)
// interface
{
  if (does_file_exist_sync_node(path + "/checksums-check.txt") or
      not does_file_exist_sync_node(path + "/checkpoint")) {
    return;
  }
  TIMER_VERBOSE("decompressed_eigen_vectors_check_crc32");
  if (not obtain_lock(path + "/lock")) {
    return;
  }
  displayln_info(0, fname + ": " + path);
  Qassert(does_file_exist_sync_node(path + "/checksums.txt"));
  Long idx_size = 0;
  std::vector<crc32_t> crcs_load;
  if (get_id_node() == 0) {
    const std::string fn = path + "/checksums.txt";
    const std::vector<std::string> lines = qgetlines(fn);
    for (size_t i = 0; i < lines.size() - 2; ++i) {
      const std::string& l = lines[i + 2];
      if (l.size() > 0 and l[0] != '\n') {
        crcs_load.push_back(read_crc32(l));
      }
    }
    idx_size = crcs_load.size();
    displayln(fname + ssprintf(": idx_size=%d", idx_size));
    glb_sum(idx_size);
    const crc32_t crc_l0 = read_crc32(lines[0]);
    const crc32_t crc_c0 = dist_crc32(crcs_load);
    displayln(
        ssprintf("summary crc stored=%08X computed=%08X", crc_l0, crc_c0));
    Qassert(crc_l0 == crc_c0);
  } else {
    glb_sum(idx_size);
    crcs_load.resize(idx_size);
  }
  bcast(get_data(crcs_load));
  std::vector<crc32_t> crcs(idx_size, 0);
  // #pragma omp parallel for
  Long num_failure = 0;
  for (Int idx = get_id_node(); idx < idx_size; idx += get_num_node()) {
    const Int dir_idx = compute_dist_file_dir_id(idx, idx_size);
    const std::string path_data = path + ssprintf("/%02d/%010d", dir_idx, idx);
    crcs[idx] = compute_crc32(path_data);
    displayln(ssprintf("reading crc %7d/%d stored=%08X computed=%08X", idx,
                       idx_size, crcs_load[idx], crcs[idx]));
    if (crcs_load[idx] != crcs[idx]) {
      displayln(ssprintf(
          "reading crc %7d/%d stored=%08X computed=%08X mismatch ('%s')", idx,
          idx_size, crcs_load[idx], crcs[idx], path.c_str()));
      qremove(path_data);
      qremove(path_data + ".crc32");
      qremove(path + "/checkpoint");
      num_failure += 1;
    }
  }
  glb_sum(get_data_char(crcs));
  glb_sum(num_failure);
  if (0 == get_id_node() and num_failure == 0) {
    crc32_t crc = dist_crc32(crcs);
    const std::string fn = path + "/checksums-check.txt";
    QFile fp = qfopen(fn + ".partial", "w");
    Qassert(not fp.null());
    qwrite_data(ssprintf("%08X\n", crc), fp);
    qwrite_data("\n", fp);
    for (size_t i = 0; i < crcs.size(); ++i) {
      qwrite_data(ssprintf("%08X\n", crcs[i]), fp);
    }
    qfclose(fp);
    qrename(fn + ".partial", fn);
  }
  release_lock();
}

bool eigen_system_repartition(const Coordinate& new_size_node,
                              const std::string& path,
                              const std::string& new_path)
// interface_function
{
  TIMER_VERBOSE("eigen_system_repartition");
  bool is_failed = false;
  const std::string npath = remove_trailing_slashes(path);
  if (std::string(npath, npath.length() - 4, 4) == ".tmp") {
    return true;
  }
  const std::string new_npath = remove_trailing_slashes(new_path);
  if (not does_file_exist_sync_node(npath + "/metadata.txt")) {
    displayln_info(
        0, ssprintf("repartition: WARNING: not a folder to partition: '%s'.",
                    npath.c_str()));
    return true;
  }
  CompressedEigenSystemInfo cesi;
  cesi = read_compressed_eigen_system_info(npath);
  if (cesi.total_node == new_size_node and
      (new_npath == "" or new_npath == npath)) {
    displayln_info(
        0, fname + ssprintf(": size_node=%s ; no need to repartition '%s'.",
                            show(cesi.total_node).c_str(), npath.c_str()));
    return true;
  } else if (new_npath == npath or new_npath == "") {
    Qassert(not does_file_exist_sync_node(npath + "-repartition-new.tmp"));
    Qassert(not does_file_exist_sync_node(npath + "-repartition-old.tmp"));
    is_failed = resize_compressed_eigen_vectors(
        npath, npath + "-repartition-new.tmp", new_size_node);
    SYNC_NODE();
    if (does_file_exist_sync_node(npath +
                                  "-repartition-new.tmp/metadata.txt")) {
      qrename_info(npath, npath + "-repartition-old.tmp");
      qrename_info(npath + "-repartition-new.tmp", npath);
      qremove_all_info(npath + "-repartition-old.tmp");
    }
  } else {
    is_failed =
        resize_compressed_eigen_vectors(npath, new_npath, new_size_node);
    SYNC_NODE();
  }
  return is_failed;
}

}  // namespace qlat
