#pragma once

#include <qlat/field-serial-io.h>
#include <qlat/selected-field.h>

namespace qlat
{  //

inline void mk_grid_field_selection(FieldM<int64_t, 1>& f_rank,
                                    const Coordinate& total_site,
                                    const long n_per_tslice_,
                                    const RngState& rs)
// each time slice has "n_per_tslice = spatial_vol / ratio" points been
// selected. not selected points has value = -1; selected points has value from
// 0 to "n_per_tslice - 1" in random order (different per time slice)
{
  TIMER_VERBOSE("mk_grid_field_selection");
  const long spatial_vol = total_site[0] * total_site[1] * total_site[2];
  const long n_per_tslice = n_per_tslice_ == -1 ? spatial_vol : n_per_tslice_;
  qassert(n_per_tslice > 0);
  qassert(spatial_vol % n_per_tslice == 0);
  const long ratio = spatial_vol / n_per_tslice;
  qassert(ratio >= 0);
  RngState rs_shift(rs, "random_shift");
  const Coordinate random_shift = mod(
      Coordinate(rand_gen(rs_shift), rand_gen(rs_shift), rand_gen(rs_shift), 0),
      total_site);
  Geometry geo;
  geo.init(total_site, 1);
  f_rank.init();
  f_rank.init(geo);
  qassert(f_rank.geo().is_only_local);
  qthread_for(index, geo.local_volume(), {
    f_rank.get_elem(index) = -1;
  });
  std::vector<Field<int64_t> > fs;
  const Coordinate new_size_node = get_default_serial_new_size_node(geo);
  shuffle_field(fs, f_rank, new_size_node);
  qassert(fs.size() <= 1);
  if (fs.size() == 1) {
    // require each tslice is on one node.
    Field<int64_t>& nf = fs[0];
    const Geometry& ngeo = nf.geo();
    qassert(ngeo.multiplicity == 1);
    qassert(new_size_node == ngeo.geon.size_node);
    const int t_start = ngeo.node_site[3] * ngeo.geon.coor_node[3];
    const int t_end = t_start + ngeo.node_site[3];
#pragma omp parallel for
    for (int t = t_start; t < t_end; ++t) {
      RngState rst = rs.split(t);
      std::vector<long> ranks(n_per_tslice);
      for (long i = 0; i < n_per_tslice; ++i) {
        ranks[i] = i;
      }
      random_permute(ranks, rst);
      long idx = 0;
      for (long index = 0; index < ngeo.local_volume(); ++index) {
        const Coordinate xl = ngeo.coordinate_from_index(index);
        const Coordinate xg = ngeo.coordinate_g_from_l(xl);
        if (xg[3] != t) {
          continue;
        }
        const Coordinate x = mod(xg + random_shift, total_site);
        bool check = true;
        switch (ratio) {
          case 1:
            check = true;
            break;
          case 2:
            check = check and (x[0] + x[1] + x[2]) % 2 == 0;
            break;
          case 4:
            check = check and (x[0] + x[1]) % 2 == 0;
            check = check and (x[1] + x[2]) % 2 == 0;
            check = check and (x[0] + x[2]) % 2 == 0;
            break;
          case 8:
            check = check and x[0] % 2 == 0;
            check = check and x[1] % 2 == 0;
            check = check and x[2] % 2 == 0;
            break;
          case 16:
            check = check and x[0] % 2 == 0;
            check = check and x[1] % 2 == 0;
            check = check and x[2] % 2 == 0;
            check = check and (x[0] + x[1] + x[2]) % 4 == 0;
            break;
          case 32:
            check = check and x[0] % 2 == 0;
            check = check and x[1] % 2 == 0;
            check = check and x[2] % 2 == 0;
            check = check and (x[0] + x[1]) % 4 == 0;
            check = check and (x[1] + x[2]) % 4 == 0;
            check = check and (x[0] + x[2]) % 4 == 0;
            break;
          case 64:
            check = check and x[0] % 4 == 0;
            check = check and x[1] % 4 == 0;
            check = check and x[2] % 4 == 0;
            break;
          default:
            displayln_info(fname + ssprintf(": ERROR: ratio=%d.", ratio));
            qassert(false);
            break;
        }
        if (check) {
          qassert(idx < n_per_tslice);
          nf.get_elem(index) = ranks[idx];
          idx += 1;
        } else {
          nf.get_elem(index) = -1;
        }
      }
      qassert(idx == n_per_tslice);
    }
  }
  shuffle_field_back(f_rank, fs, new_size_node);
}

inline void mk_field_selection(FieldM<int64_t, 1>& f_rank,
                               const Coordinate& total_site,
                               const long n_per_tslice, const RngState& rs)
// interface function
// not selected points has value = -1;
// random select n_per_tslice points based on ranks from mk_grid_field_selection
{
  TIMER_VERBOSE("mk_field_selection(n_per_tslice,rs)");
  const long spatial_vol = total_site[0] * total_site[1] * total_site[2];
  mk_grid_field_selection(f_rank, total_site, spatial_vol, rs);
  if (n_per_tslice == -1 or n_per_tslice == spatial_vol) {
    return;
  }
  set_n_per_tslice(f_rank, n_per_tslice);
}

inline void set_grid_field_selection(FieldSelection& fsel,
                                     const Coordinate& total_site,
                                     const long n_per_tslice,
                                     const RngState& rs)
{
  TIMER_VERBOSE("set_grid_field_selection(fsel,total_site,n_per_tslice,rs)");
  fsel.init();
  mk_grid_field_selection(fsel.f_rank, total_site, n_per_tslice, rs);
  update_field_selection(fsel);
  update_field_selection(fsel, n_per_tslice);
}

inline void set_field_selection(FieldSelection& fsel,
                                const Coordinate& total_site,
                                const long n_per_tslice, const RngState& rs)
{
  TIMER_VERBOSE("set_field_selection(fsel,total_site,n_per_tslice,rs)");
  fsel.init();
  mk_field_selection(fsel.f_rank, total_site, n_per_tslice, rs);
  update_field_selection(fsel);
  update_field_selection(fsel, n_per_tslice);
}

inline void set_field_selection(FieldSelection& fsel,
                                const Coordinate& total_site,
                                const long n_per_tslice, const RngState& rs,
                                const PointSelection& psel)
{
  TIMER_VERBOSE("set_field_selection(fsel,total_site,n_per_tslice,rs,psel)");
  fsel.init();
  mk_field_selection(fsel.f_rank, total_site, n_per_tslice, rs);
  add_field_selection(fsel.f_rank, psel);
  update_field_selection(fsel);
  update_field_selection(fsel, n_per_tslice);
}

inline long write_field_selection(const FieldSelection& fsel,
                                  const std::string& path)
{
  TIMER_VERBOSE("write_field_selection");
  return write_field_64(fsel.f_rank, path);
}

inline long read_field_selection(FieldSelection& fsel, const std::string& path,
                                 const long n_per_tslice)
{
  TIMER_VERBOSE("read_field_selection");
  fsel.init();
  FieldM<int64_t, 1> f_rank;
  const long total_bytes = read_field_64(f_rank, path);
  if (total_bytes > 0) {
    set_field_selection(fsel, f_rank, n_per_tslice);
  }
  return total_bytes;
}

template <class M>
crc32_t field_crc32(const SelectedField<M>& sf, const FieldSelection& fsel,
                    const Coordinate& new_size_node_ = Coordinate())
{
  TIMER_VERBOSE_FLOPS("field_crc32(sf)");
  const Geometry& geo = sf.geo();
  qassert(geo.is_only_local);
  qassert(fsel.f_rank.geo() == geo_remult(geo));
  // const Coordinate total_site = geo.total_site();
  const Coordinate new_size_node = new_size_node_ != Coordinate()
                                       ? new_size_node_
                                       : get_default_serial_new_size_node(geo);
  qassert(new_size_node[0] == 1);
  qassert(new_size_node[1] == 1);
  qassert(new_size_node[2] == 1);
  std::vector<FieldSelection> fsels;
  const ShufflePlan sp = make_shuffle_plan(fsels, fsel, new_size_node);
  std::vector<SelectedField<M> > sfs;
  shuffle_field(sfs, sf, sp);
  qassert(fsels.size() == sfs.size());
  const int new_num_node = product(new_size_node);
  crc32_t crc = 0;
  for (int i = 0; i < (int)sfs.size(); ++i) {
    const int new_id_node = sfs[i].geo().geon.id_node;
    qassert(sfs[i].geo().geon.num_node == new_num_node);
    const Vector<M> v = get_data(sfs[i].field);
    crc ^= crc32_shift(crc32_par(v),
                       (new_num_node - new_id_node - 1) * v.data_size());
  }
  glb_sum_byte(crc);
  return crc;
}

inline std::string make_selected_field_header(const Geometry& geo,
                                              const long n_per_tslice,
                                              const int sizeof_M,
                                              const crc32_t crc32)
{
  const Coordinate total_site = geo.total_site();
  std::ostringstream out;
  // const std::string todo = "NOT yet implemented";
  out << "BEGIN_SELECTED_FIELD_HEADER" << std::endl;
  out << "selected_field_version = 1.0" << std::endl;
  out << "total_site[0] = " << total_site[0] << std::endl;
  out << "total_site[1] = " << total_site[1] << std::endl;
  out << "total_site[2] = " << total_site[2] << std::endl;
  out << "total_site[3] = " << total_site[3] << std::endl;
  out << "n_per_tslice = " << n_per_tslice << std::endl;
  out << "multiplicity = " << geo.multiplicity << std::endl;
  out << "sizeof(M) = " << sizeof_M << std::endl;
  out << ssprintf("selected_field_crc32 = %08X", crc32) << std::endl;
  out << "END_HEADER" << std::endl;
  return out.str();
}

template <class M>
long write_selected_field(const SelectedField<M>& sf, const std::string& path,
                          const FieldSelection& fsel,
                          const Coordinate& new_size_node_ = Coordinate())
{
  TIMER_VERBOSE_FLOPS("write_selected_field");
  displayln_info(fname + ssprintf(": fn='%s'.", path.c_str()));
  const Geometry& geo = sf.geo();
  qassert(geo.is_only_local);
  qassert(fsel.f_rank.geo() == geo_remult(geo));
  const Coordinate total_site = geo.total_site();
  const Coordinate new_size_node = new_size_node_ != Coordinate()
                                       ? new_size_node_
                                       : get_default_serial_new_size_node(geo);
  qassert(new_size_node[0] == 1);
  qassert(new_size_node[1] == 1);
  qassert(new_size_node[2] == 1);
  std::vector<FieldSelection> fsels;
  const ShufflePlan sp = make_shuffle_plan(fsels, fsel, new_size_node);
  std::vector<SelectedField<M> > sfs;
  shuffle_field(sfs, sf, sp);
  qassert(fsels.size() == sfs.size());
  long check_n_per_tslice = 0;
  for (int i = 0; i < (int)sfs.size(); ++i) {
    const Coordinate& node_site = sfs[i].geo().node_site;
    qassert(node_site[0] == total_site[0]);
    qassert(node_site[1] == total_site[1]);
    qassert(node_site[2] == total_site[2]);
    if (sfs[i].n_elems != node_site[3] * fsel.n_per_tslice) {
      check_n_per_tslice += 1;
    }
  }
  glb_sum(check_n_per_tslice);
  if (check_n_per_tslice > 0) {
    displayln_info(
        fname +
        ssprintf(
            ": WARNING fsels.n_per_tslice=%d do not match with data. n_fail=%d",
            fsel.n_per_tslice, check_n_per_tslice));
  }
  const int new_num_node = product(new_size_node);
  crc32_t crc = 0;
  for (int i = 0; i < (int)sfs.size(); ++i) {
    const int new_id_node = sfs[i].geo().geon.id_node;
    qassert(sfs[i].geo().geon.num_node == new_num_node);
    const Vector<M> v = get_data(sfs[i].field);
    crc ^= crc32_shift(crc32_par(v),
                       (new_num_node - new_id_node - 1) * v.data_size());
  }
  glb_sum_byte(crc);
  if (get_force_field_write_sizeof_M() == 0) {
    qtouch_info(path + ".partial", make_selected_field_header(
                                       geo, fsel.n_per_tslice, sizeof(M), crc));
  } else {
    const int sizeof_M = get_force_field_write_sizeof_M();
    qassert((geo.multiplicity * sizeof(M)) % sizeof_M == 0);
    const int multiplicity = (geo.multiplicity * sizeof(M)) / sizeof_M;
    qtouch_info(path + ".partial",
                make_selected_field_header(geo_remult(geo, multiplicity),
                                           fsel.n_per_tslice, sizeof_M, crc));
    get_force_field_write_sizeof_M() = 0;
  }
  const int mpi_tag = 8;
  if (get_id_node() == 0) {
    qassert(sfs.size() > 0);
    Vector<M> v = get_data(sfs[0].field);
    QFile fp = qfopen(path + ".partial", "a");
    qassert(not fp.null());
    const int num_node = get_num_node();
    for (int new_id_node = 0; new_id_node < new_num_node; ++new_id_node) {
      const int id_node =
          get_id_node_from_new_id_node(new_id_node, new_num_node, num_node);
      if (0 == id_node) {
        assign(v, get_data(sfs[new_id_node].field));
      } else {
        mpi_recv(v.data(), v.data_size(), MPI_BYTE, id_node, mpi_tag,
                 get_comm(), MPI_STATUS_IGNORE);
      }
      qwrite_data(v, fp);
    }
    qfclose(fp);
  } else {
    for (size_t i = 0; i < sfs.size(); ++i) {
      const Vector<M> v = get_data(sfs[i].field);
      mpi_send((void*)v.data(), v.data_size(), MPI_BYTE, 0, mpi_tag,
               get_comm());
    }
  }
  qrename_info(path + ".partial", path);
  const long total_bytes =
      fsel.n_per_tslice * total_site[3] * geo.multiplicity * sizeof(M);
  timer.flops += total_bytes;
  return total_bytes;
}

inline long read_selected_geo_info(Coordinate& total_site, int& multiplicity,
                                   long& n_per_tslice, int& sizeof_M,
                                   crc32_t& crc, const std::string& path)
{
  TIMER("read_selected_geo_info");
  long pos = 0;
  if (get_id_node() == 0) {
    QFile fp = qfopen(path, "r");
    if (not fp.null()) {
      const std::string header = "BEGIN_SELECTED_FIELD_HEADER\n";
      std::vector<char> check_line(header.size(), 0);
      if (1 == qfread(check_line.data(), header.size(), 1, fp)) {
        if (std::string(check_line.data(), check_line.size()) == header) {
          std::vector<std::string> infos;
          infos.push_back(header);
          while (infos.back() != "END_HEADER\n" && infos.back() != "") {
            infos.push_back(qgetline(fp));
          }
          for (int m = 0; m < 4; ++m) {
            reads(total_site[m],
                  info_get_prop(infos, ssprintf("total_site[%d] = ", m)));
          }
          reads(multiplicity, info_get_prop(infos, "multiplicity = "));
          reads(n_per_tslice, info_get_prop(infos, "n_per_tslice = "));
          reads(sizeof_M, info_get_prop(infos, "sizeof(M) = "));
          crc = read_crc32(info_get_prop(infos, "selected_field_crc32 = "));
        }
      }
      pos = qftell(fp);
    }
    qfclose(fp);
  }
  bcast(get_data_one_elem(pos));
  bcast(get_data_one_elem(total_site));
  bcast(get_data_one_elem(multiplicity));
  bcast(get_data_one_elem(n_per_tslice));
  bcast(get_data_one_elem(sizeof_M));
  bcast(get_data_one_elem(crc));
  return pos;
}

template <class M>
long read_selected_field(SelectedField<M>& sf, const std::string& path,
                         const FieldSelection& fsel,
                         const Coordinate& new_size_node_ = Coordinate())
{
  TIMER_VERBOSE_FLOPS("read_selected_field");
  displayln_info(fname + ssprintf(": fn='%s'.", path.c_str()));
  sf.init();
  Coordinate total_site;
  int multiplicity = 0;
  long n_per_tslice = 0;
  int sizeof_M = 0;
  crc32_t crc_info = 0;
  const long pos = read_selected_geo_info(
      total_site, multiplicity, n_per_tslice, sizeof_M, crc_info, path);
  if (total_site == Coordinate() or multiplicity == 0) {
    displayln_info(fname +
                   ssprintf(": fn='%s' can not be parsed.", path.c_str()));
    return 0;
  }
  qassert(fsel.n_per_tslice == n_per_tslice);
  get_incorrect_field_read_sizeof_M() = 0;
  if (sizeof_M != sizeof(M)) {
    get_incorrect_field_read_sizeof_M() = sizeof_M;
    displayln(fname + ssprintf(": WARNING: sizeof(M) do not match. "
                               "Expected %d, Actual file %d",
                               sizeof(M), sizeof_M));
    qassert((multiplicity * sizeof_M) % sizeof(M) == 0);
    multiplicity = (multiplicity * sizeof_M) / sizeof(M);
  }
  Geometry geo;
  geo.init(total_site, multiplicity);
  qassert(fsel.f_rank.geo() == geo_remult(geo));
  const Coordinate new_size_node = new_size_node_ != Coordinate()
                                       ? new_size_node_
                                       : get_default_serial_new_size_node(geo);
  qassert(new_size_node[0] == 1);
  qassert(new_size_node[1] == 1);
  qassert(new_size_node[2] == 1);
  std::vector<FieldSelection> fsels;
  const ShufflePlan sp = make_shuffle_plan(fsels, fsel, new_size_node);
  if (not is_no_shuffle(sp)) {
    qassert(fsels.size() == sp.geos_recv.size());
  }
  std::vector<SelectedField<M> > sfs;
  sfs.resize(fsels.size());
  for (size_t i = 0; i < sfs.size(); ++i) {
    sfs[i].init(fsels[i], geo.multiplicity);
  }
  long check_n_per_tslice = 0;
  for (int i = 0; i < (int)sfs.size(); ++i) {
    const Coordinate& node_site = sfs[i].geo().node_site;
    qassert(node_site[0] == total_site[0]);
    qassert(node_site[1] == total_site[1]);
    qassert(node_site[2] == total_site[2]);
    if (sfs[i].n_elems != node_site[3] * fsel.n_per_tslice) {
      check_n_per_tslice += 1;
    }
  }
  glb_sum(check_n_per_tslice);
  if (check_n_per_tslice > 0) {
    displayln_info(
        fname +
        ssprintf(
            ": WARNING fsels.n_per_tslice=%d do not match with data. n_fail=%d",
            fsel.n_per_tslice, check_n_per_tslice));
  }
  const int new_num_node = product(new_size_node);
  crc32_t crc = 0;
  if (sfs.size() > 0) {
    QFile fp = qfopen(path, "r");
    qassert(not fp.null());
    qfseek(fp,
           pos + sfs[0].geo().geon.id_node * get_data(sfs[0].field).data_size(),
           SEEK_SET);
    for (int i = 0; i < (int)sfs.size(); ++i) {
      const int new_id_node = sfs[i].geo().geon.id_node;
      qassert(sfs[i].geo().geon.num_node == new_num_node);
      const Vector<M> v = get_data(sfs[i].field);
      qread_data(v, fp);
      crc ^= crc32_shift(crc32_par(v),
                         (new_num_node - new_id_node - 1) * v.data_size());
    }
    qfclose(fp);
  }
  glb_sum_byte(crc);
  if (crc != crc_info) {
    displayln_info(
        fname +
        ssprintf(": crc of data = %08X ; crc of header = %08X", crc, crc_info));
    qassert(false);
  }
  sf.init(fsel, geo.multiplicity);
  shuffle_field_back(sf, sfs, sp);
  const long total_bytes =
      fsel.n_per_tslice * total_site[3] * geo.multiplicity * sizeof(M);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
long write_selected_field_64(const SelectedField<M>& f, const std::string& path,
                             const FieldSelection& fsel,
                             const Coordinate& new_size_node_ = Coordinate())
{
  TIMER_VERBOSE_FLOPS("write_selected_field_64");
  SelectedField<M> ff;
  ff = f;
  to_from_big_endian_64(get_data(ff));
  const long total_bytes = write_selected_field(ff, path, fsel, new_size_node_);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
long read_selected_field_64(SelectedField<M>& sf, const std::string& path,
                            const FieldSelection& fsel,
                            const Coordinate& new_size_node_ = Coordinate())
{
  TIMER_VERBOSE_FLOPS("read_selected_field_64");
  const long total_bytes = read_selected_field(sf, path, fsel, new_size_node_);
  if (total_bytes == 0) {
    return 0;
  } else {
    to_from_big_endian_64(get_data(sf));
    timer.flops += total_bytes;
    return total_bytes;
  }
}

template <class M>
long write_selected_field_double(
    const SelectedField<M>& f, const std::string& path,
    const FieldSelection& fsel, const Coordinate& new_size_node_ = Coordinate())
{
  TIMER_VERBOSE_FLOPS("write_selected_field_double");
  SelectedField<M> ff;
  ff = f;
  to_from_big_endian_64(get_data(ff));
  const long total_bytes = write_selected_field(ff, path, fsel, new_size_node_);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
long read_selected_field_double(SelectedField<M>& sf, const std::string& path,
                                const FieldSelection& fsel,
                                const Coordinate& new_size_node_ = Coordinate())
{
  TIMER_VERBOSE_FLOPS("read_selected_field_double");
  const long total_bytes = read_selected_field(sf, path, fsel, new_size_node_);
  if (total_bytes == 0) {
    return 0;
  } else {
    to_from_big_endian_64(get_data(sf));
    timer.flops += total_bytes;
    return total_bytes;
  }
}

template <class M, class N>
void convert_field_float_from_double(SelectedField<N>& ff,
                                     const SelectedField<M>& f)
// interface_function
{
  TIMER("convert_field_float_from_double(sf)");
  qassert(f.geo().is_only_local);
  qassert(sizeof(M) % sizeof(double) == 0);
  qassert(sizeof(N) % sizeof(float) == 0);
  qassert(f.geo().multiplicity * sizeof(M) / 2 % sizeof(N) == 0);
  const int multiplicity = f.geo().multiplicity * sizeof(M) / 2 / sizeof(N);
  const Geometry geo = geo_remult(f.geo(), multiplicity);
  const long n_elems = f.n_elems;
  ff.init(geo, n_elems, multiplicity);
  const Vector<M> fdata = get_data(f);
  const Vector<double> fd((double*)fdata.data(),
                          fdata.data_size() / sizeof(double));
  Vector<N> ffdata = get_data(ff);
  Vector<float> ffd((float*)ffdata.data(), ffdata.data_size() / sizeof(float));
  qassert(ffd.size() == fd.size());
  qacc_for(i, ffd.size(), {
    ffd[i] = fd[i];
  });
}

template <class M, class N>
void convert_field_double_from_float(SelectedField<N>& ff,
                                     const SelectedField<M>& f)
// interface_function
{
  TIMER("convert_field_double_from_float(sf)");
  qassert(f.geo().is_only_local);
  qassert(sizeof(M) % sizeof(float) == 0);
  qassert(sizeof(N) % sizeof(double) == 0);
  qassert(f.geo().multiplicity * sizeof(M) * 2 % sizeof(N) == 0);
  const int multiplicity = f.geo().multiplicity * sizeof(M) * 2 / sizeof(N);
  const Geometry geo = geo_remult(f.geo(), multiplicity);
  const long n_elems = f.n_elems;
  ff.init(geo, n_elems, multiplicity);
  const Vector<M> fdata = get_data(f);
  const Vector<float> fd((float*)fdata.data(),
                         fdata.data_size() / sizeof(float));
  Vector<N> ffdata = get_data(ff);
  Vector<double> ffd((double*)ffdata.data(),
                     ffdata.data_size() / sizeof(double));
  qassert(ffd.size() == fd.size());
  qacc_for(i, ffd.size(), {
    ffd[i] = fd[i];
  });
}

template <class M>
long write_selected_field_float_from_double(
    const SelectedField<M>& f, const std::string& path,
    const FieldSelection& fsel, const Coordinate& new_size_node_ = Coordinate())
{
  TIMER_VERBOSE_FLOPS("write_selected_field_float_from_double");
  SelectedField<float> ff;
  convert_field_float_from_double(ff, f);
  to_from_big_endian_32(get_data(ff));
  const long total_bytes = write_selected_field(ff, path, fsel, new_size_node_);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
long read_selected_field_double_from_float(
    SelectedField<M>& f, const std::string& path, const FieldSelection& fsel,
    const Coordinate& new_size_node_ = Coordinate())
{
  TIMER_VERBOSE_FLOPS("read_selected_field_double_from_float");
  SelectedField<float> ff;
  const long total_bytes = read_selected_field(ff, path, fsel, new_size_node_);
  if (total_bytes == 0) {
    return 0;
  } else {
    to_from_big_endian_32(get_data(ff));
    convert_field_double_from_float(f, ff);
    timer.flops += total_bytes;
    return total_bytes;
  }
}

inline bool is_selected_field(const std::string& path)
{
  TIMER("is_selected_field");
  long nfile = 0;
  if (get_id_node() == 0) {
    QFile fp = qfopen(path, "r");
    if (not fp.null()) {
      const std::string header = "BEGIN_SELECTED_FIELD_HEADER\n";
      std::vector<char> check_line(header.size(), 0);
      if (1 == qfread(check_line.data(), header.size(), 1, fp)) {
        if (std::string(check_line.data(), check_line.size()) == header) {
          nfile = 1;
        }
      }
    }
    qfclose(fp);
  }
  bcast(get_data(nfile));
  return nfile > 0;
}

template <class M>
long write_selected_field(const Field<M>& f, const std::string& path,
                          const FieldSelection& fsel,
                          const Coordinate& new_size_node_ = Coordinate())
{
  TIMER_VERBOSE_FLOPS("write_selected_field(f)");
  SelectedField<M> sf;
  set_selected_field(sf, f, fsel);
  return write_selected_field(sf, path, fsel, new_size_node_);
}

template <class M>
long write_selected_field_64(const Field<M>& f, const std::string& path,
                             const FieldSelection& fsel,
                             const Coordinate& new_size_node_ = Coordinate())
{
  TIMER_VERBOSE_FLOPS("write_selected_field_64(f)");
  SelectedField<M> sf;
  set_selected_field(sf, f, fsel);
  return write_selected_field_64(sf, path, fsel, new_size_node_);
}

template <class M>
long write_selected_field_double(
    const Field<M>& f, const std::string& path, const FieldSelection& fsel,
    const Coordinate& new_size_node_ = Coordinate())
{
  TIMER_VERBOSE_FLOPS("write_selected_field_64(f)");
  SelectedField<M> sf;
  set_selected_field(sf, f, fsel);
  return write_selected_field_double(sf, path, fsel, new_size_node_);
}

template <class M>
long write_selected_field_float_from_double(
    const Field<M>& f, const std::string& path, const FieldSelection& fsel,
    const Coordinate& new_size_node_ = Coordinate())
{
  TIMER_VERBOSE_FLOPS("write_selected_field_64(f)");
  SelectedField<M> sf;
  set_selected_field(sf, f, fsel);
  return write_selected_field_float_from_double(sf, path, fsel, new_size_node_);
}

template <class M>
long read_selected_field(Field<M>& f, const std::string& path,
                         const FieldSelection& fsel,
                         const Coordinate& new_size_node_ = Coordinate())
{
  TIMER_VERBOSE_FLOPS("read_selected_field(f)");
  SelectedField<M> sf;
  const long total_bytes = read_selected_field(sf, path, fsel, new_size_node_);
  if (total_bytes > 0) {
    set_field_selected(f, sf, fsel);
  }
  return total_bytes;
}

template <class M>
long read_selected_field_64(Field<M>& f, const std::string& path,
                            const FieldSelection& fsel,
                            const Coordinate& new_size_node_ = Coordinate())
{
  TIMER_VERBOSE_FLOPS("read_selected_field_64(f)");
  SelectedField<M> sf;
  const long total_bytes =
      read_selected_field_64(sf, path, fsel, new_size_node_);
  if (total_bytes > 0) {
    set_field_selected(f, sf, fsel);
  }
  return total_bytes;
}

template <class M>
long read_selected_field_double(Field<M>& f, const std::string& path,
                                const FieldSelection& fsel,
                                const Coordinate& new_size_node_ = Coordinate())
{
  TIMER_VERBOSE_FLOPS("read_selected_field_double(f)");
  SelectedField<M> sf;
  const long total_bytes =
      read_selected_field_double(sf, path, fsel, new_size_node_);
  if (total_bytes > 0) {
    set_field_selected(f, sf, fsel);
  }
  return total_bytes;
}

template <class M>
long read_selected_field_double_from_float(
    Field<M>& f, const std::string& path, const FieldSelection& fsel,
    const Coordinate& new_size_node_ = Coordinate())
{
  TIMER_VERBOSE_FLOPS("read_selected_field_double_from_float(f)");
  SelectedField<M> sf;
  const long total_bytes =
      read_selected_field_double_from_float(sf, path, fsel, new_size_node_);
  if (total_bytes > 0) {
    set_field_selected(f, sf, fsel);
  }
  return total_bytes;
}

}  // namespace qlat
