#pragma once

#include <qlat/field-base-io.h>
#include <qlat/field-dist-io.h>
#include <qlat/field-serial-io.h>

namespace qlat
{  //

template <class M>
crc32_t field_simple_checksum(const Field<M>& f)
// call with data endianness native to the machine
{
  TIMER("field_simple_checksum");
  qassert(f.geo().is_only_local);
  crc32_t ret = 0;
  const Vector<M> v = get_data(f);
  Vector<crc32_t> vc((crc32_t*)v.data(), v.data_size() / sizeof(crc32_t));
  for (Long i = 0; i < vc.size(); ++i) {
    ret += vc[i];
  }
  Long sum = ret;
  glb_sum(sum);
  ret = (crc32_t)sum;
  return ret;
}

template <class M>
crc32_t field_crc32_shuffle(const Field<M>& f)
{
  TIMER_VERBOSE_FLOPS("field_crc32_shuffle");
  const Geometry& geo = f.geo();
  const Coordinate new_size_node = get_default_serial_new_size_node(geo);
  std::vector<Field<M>> fs;
  shuffle_field(fs, f, new_size_node);
  crc32_t ret = 0;
  for (Int i = 0; i < (int)fs.size(); ++i) {
    const Int new_id_node = fs[i].geo().geon.id_node;
    const Int new_num_node = fs[i].geo().geon.num_node;
    const Vector<M> v = get_data(fs[i]);
    ret ^= crc32_shift(crc32_par(v),
                       (new_num_node - new_id_node - 1) * v.data_size());
  }
  glb_sum(get_data_char(ret));
  timer.flops += get_data(f).data_size() * geo.geon.num_node;
  return ret;
}

template <class M>
crc32_t field_crc32_sites(const Field<M>& f)
{
  TIMER_VERBOSE_FLOPS("field_crc32_sites");
  const Geometry& geo = f.geo();
  const Int multiplicity = f.multiplicity;
  const Long total_volume = geo.total_volume();
  const Long data_size_site = multiplicity * sizeof(M);
  const Int v_limit = omp_get_max_threads();
  std::vector<crc32_t> crcs(v_limit, 0);
#pragma omp parallel
  {
    crc32_t crc = 0;
#pragma omp for
    for (Long index = 0; index < geo.local_volume(); ++index) {
      const Coordinate xl = geo.coordinate_from_index(index);
      const Coordinate xg = geo.coordinate_g_from_l(xl);
      const Long gindex = geo.g_index_from_g_coordinate(xg);
      const Long offset = data_size_site * (total_volume - gindex - 1);
      const Vector<M> v = f.get_elems_const(xl);
      crc ^= crc32_shift(crc32(v), offset);
    }
    const Int id = omp_get_thread_num();
    crcs[id] = crc;
  }
  crc32_t ret = 0;
  for (Int i = 0; i < v_limit; ++i) {
    ret ^= crcs[i];
  }
  glb_sum(get_data_char(ret));
  timer.flops += get_data(f).data_size() * geo.geon.num_node;
  return ret;
}

template <class M>
crc32_t field_crc32(const Field<M>& f)
{
  return field_crc32_shuffle(f);
}

// ----------------------

template <class M>
Long dist_write_field(const Field<M>& f, const Coordinate& new_size_node,
                      const std::string& path)
// interface_function
{
  TIMER_VERBOSE_FLOPS("dist_write_field");
  displayln_info(fname + ssprintf(": fn='%s'.", path.c_str()));
  std::vector<Field<M>> fs;
  shuffle_field(fs, f, new_size_node);
  Long total_bytes = dist_write_fields(fs, product(new_size_node), path);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
Long dist_read_field(Field<M>& f, const std::string& path)
// interface_function
{
  TIMER_VERBOSE_FLOPS("dist_read_field");
  displayln_info(fname + ssprintf(": fn='%s'.", path.c_str()));
  Geometry geo;
  Int multiplicity;
  std::vector<Field<M>> fs;
  Coordinate new_size_node;
  const Long total_bytes = dist_read_fields(fs, geo, multiplicity, new_size_node, path);
  if (total_bytes == 0) {
    return 0;
  } else {
    f.init(geo, multiplicity);
    qassert(f.geo() == geo);
    qassert(f.multiplicity == multiplicity);
    shuffle_field_back(f, fs, new_size_node);
    timer.flops += total_bytes;
    return total_bytes;
  }
}

template <class M, class N>
void convert_field_float_from_double(Field<N>& ff, const Field<M>& f)
// interface_function
{
  TIMER("convert_field_float_from_double");
  const Geometry& geo = f.geo();
  qassert(geo.is_only_local);
  qassert(sizeof(M) % sizeof(RealD) == 0);
  qassert(sizeof(N) % sizeof(RealF) == 0);
  qassert(f.multiplicity * sizeof(M) / 2 % sizeof(N) == 0);
  const Int multiplicity = f.multiplicity * sizeof(M) / 2 / sizeof(N);
  ff.init(geo, multiplicity);
  const Vector<M> fdata = get_data(f);
  const Vector<double> fd((double*)fdata.data(),
                          fdata.data_size() / sizeof(RealD));
  Vector<N> ffdata = get_data(ff);
  Vector<RealF> ffd((RealF*)ffdata.data(), ffdata.data_size() / sizeof(RealF));
  qassert(ffd.size() == fd.size());
  qacc_for(i, ffd.size(), { ffd[i] = fd[i]; });
}

template <class M, class N>
void convert_field_double_from_float(Field<N>& ff, const Field<M>& f)
// interface_function
{
  TIMER("convert_field_double_from_float");
  const Geometry& geo = f.geo();
  const Int multiplicity = f.multiplicity * sizeof(M) * 2 / sizeof(N);
  qassert(geo.is_only_local);
  qassert(sizeof(M) % sizeof(RealF) == 0);
  qassert(sizeof(N) % sizeof(RealD) == 0);
  qassert(f.multiplicity * sizeof(M) * 2 % sizeof(N) == 0);
  ff.init(geo, multiplicity);
  const Vector<M> fdata = get_data(f);
  const Vector<RealF> fd((RealF*)fdata.data(),
                         fdata.data_size() / sizeof(RealF));
  Vector<N> ffdata = get_data(ff);
  Vector<double> ffd((double*)ffdata.data(),
                     ffdata.data_size() / sizeof(RealD));
  qassert(ffd.size() == fd.size());
  qacc_for(i, ffd.size(), { ffd[i] = fd[i]; });
}

template <class M>
Long dist_write_field_float_from_double(const Field<M>& f,
                                        const std::string& path)
// interface_function
{
  TIMER_VERBOSE_FLOPS("dist_write_field_float_from_double");
  Field<RealF> ff;
  convert_field_float_from_double(ff, f);
  to_from_big_endian(get_data(ff));
  const Long total_bytes = dist_write_field(ff, path);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
Long dist_write_field_float_from_double(const Field<M>& f,
                                        const Coordinate& new_size_node,
                                        const std::string& path)
// interface_function
{
  TIMER_VERBOSE_FLOPS("dist_write_field_float_from_double");
  Field<RealF> ff;
  convert_field_float_from_double(ff, f);
  to_from_big_endian(get_data(ff));
  const Long total_bytes = dist_write_field(ff, new_size_node, path);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
Long dist_read_field_double_from_float(Field<M>& f, const std::string& path)
// interface_function
{
  TIMER_VERBOSE_FLOPS("dist_read_field_double_from_float");
  Field<RealF> ff;
  const Long total_bytes = dist_read_field(ff, path);
  if (total_bytes == 0) {
    return 0;
  } else {
    to_from_big_endian(get_data(ff));
    convert_field_double_from_float(f, ff);
    timer.flops += total_bytes;
    return total_bytes;
  }
}

template <class M>
Long dist_write_field_double(const Field<M>& f, const std::string& path)
// interface_function
{
  TIMER_VERBOSE_FLOPS("dist_write_field_double");
  Field<M> ff;
  ff.init(f);
  to_from_big_endian(get_data(ff));
  const Long total_bytes = dist_write_field(ff, path);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
Long dist_write_field_double(const Field<M>& f, const Coordinate& new_size_node,
                             const std::string& path)
// interface_function
{
  TIMER_VERBOSE_FLOPS("dist_write_field_double");
  Field<M> ff;
  ff.init(f);
  to_from_big_endian(get_data(ff));
  const Long total_bytes = dist_write_field(ff, new_size_node, path);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
Long dist_read_field_double(Field<M>& f, const std::string& path)
// interface_function
{
  TIMER_VERBOSE_FLOPS("dist_read_field_double");
  const Long total_bytes = dist_read_field(f, path);
  if (total_bytes == 0) {
    return 0;
  } else {
    to_from_big_endian(get_data(f));
    timer.flops += total_bytes;
    return total_bytes;
  }
}

// ----------------------

inline std::string make_field_header(const Geometry& geo, const Int multiplicity, const Int sizeof_M,
                                     const crc32_t crc32)
{
  const Coordinate total_site = geo.total_site();
  std::ostringstream out;
  // const std::string todo = "NOT yet implemented";
  out << "BEGIN_FIELD_HEADER" << std::endl;
  out << "field_version = 1.0" << std::endl;
  out << "total_site[0] = " << total_site[0] << std::endl;
  out << "total_site[1] = " << total_site[1] << std::endl;
  out << "total_site[2] = " << total_site[2] << std::endl;
  out << "total_site[3] = " << total_site[3] << std::endl;
  out << "multiplicity = " << multiplicity << std::endl;
  out << "sizeof(M) = " << sizeof_M << std::endl;
  out << ssprintf("field_crc32 = %08X", crc32) << std::endl;
  out << "END_HEADER" << std::endl;
  return out.str();
}

inline void read_geo_info(Coordinate& total_site, Int& multiplicity,
                          Int& sizeof_M, crc32_t& crc, const std::string& path)
{
  TIMER("read_geo_info");
  if (get_id_node() == 0) {
    QFile qfile = qfopen(path, "r");
    if (not qfile.null()) {
      const std::string header = "BEGIN_FIELD_HEADER\n";
      std::vector<char> check_line(header.size(), 0);
      if (1 == qfread(check_line.data(), header.size(), 1, qfile)) {
        if (std::string(check_line.data(), check_line.size()) == header) {
          std::vector<std::string> infos;
          infos.push_back(header);
          while (infos.back() != "END_HEADER\n" && infos.back() != "") {
            infos.push_back(qgetline(qfile));
          }
          for (Int m = 0; m < 4; ++m) {
            reads(total_site[m],
                  info_get_prop(infos, ssprintf("total_site[%d] = ", m)));
          }
          reads(multiplicity, info_get_prop(infos, "multiplicity = "));
          reads(sizeof_M, info_get_prop(infos, "sizeof(M) = "));
          crc = read_crc32(info_get_prop(infos, "field_crc32 = "));
        }
      }
    }
    qfclose(qfile);
  }
  bcast(Vector<Coordinate>(&total_site, 1));
  bcast(Vector<int>(&multiplicity, 1));
  bcast(Vector<int>(&sizeof_M, 1));
  bcast(Vector<crc32_t>(&crc, 1));
}

template <class M>
Long write_field(const Field<M>& f, const std::string& path,
                 const Coordinate& new_size_node = Coordinate())
// if new_size_node != Coordinate() then use dist_write_field
{
  if (new_size_node != Coordinate()) {
    return dist_write_field(f, new_size_node, path);
  }
  TIMER_VERBOSE_FLOPS("write_field");
  displayln_info(fname + ssprintf(": fn='%s'.", path.c_str()));
  qassert(is_initialized(f));
  const Geometry& geo = f.geo();
  const crc32_t crc32 = field_crc32(f);
  if (get_force_field_write_sizeof_M() == 0) {
    qtouch_info(path + ".partial", make_field_header(geo, f.multiplicity, sizeof(M), crc32));
  } else {
    const Int sizeof_M = get_force_field_write_sizeof_M();
    qassert((f.multiplicity * sizeof(M)) % sizeof_M == 0);
    const Int multiplicity = (f.multiplicity * sizeof(M)) / sizeof_M;
    qtouch_info(
        path + ".partial",
        make_field_header(geo, multiplicity, sizeof_M, crc32));
    get_force_field_write_sizeof_M() = 0;
  }
  const Long file_size = serial_write_field(
      f, path + ".partial",
      get_default_serial_new_size_node(geo, dist_write_par_limit()));
  qrename_info(path + ".partial", path);
  timer.flops += file_size;
  return file_size;
}

template <class M>
Long read_field(Field<M>& f, const std::string& path,
                const Coordinate& new_size_node_ = Coordinate())
// assume new_size_node is properly choosen so that concatenate the new fields
// would be correct. eg. new_size_node = Coordinate(1,1,1,2)
{
  if (does_file_exist_qar_sync_node(path + "/geo-info.txt")) {
    return dist_read_field(f, path);
  }
  TIMER_VERBOSE_FLOPS("read_field");
  displayln_info(fname + ssprintf(": fn='%s'.", path.c_str()));
  Coordinate total_site;
  Int multiplicity = 0;
  Int sizeof_M = 0;
  crc32_t crc = 0;
  read_geo_info(total_site, multiplicity, sizeof_M, crc, path);
  if (total_site == Coordinate() or multiplicity == 0) {
    return 0;
  }
  get_incorrect_field_read_sizeof_M() = 0;
  if (sizeof_M != sizeof(M)) {
    get_incorrect_field_read_sizeof_M() = sizeof_M;
    displayln_info(fname + ssprintf(": WARNING: sizeof(M) do not match. "
                                    "Expected %d, Actual file %d",
                                    sizeof(M), sizeof_M));
    qassert((multiplicity * sizeof_M) % sizeof(M) == 0);
    multiplicity = (multiplicity * sizeof_M) / sizeof(M);
  }
  Geometry geo;
  geo.init(total_site);
  f.init(geo, multiplicity);
  const Long data_size =
      geo.geon.num_node * geo.local_volume() * multiplicity * sizeof(M);
  const Coordinate new_size_node =
      new_size_node_ == Coordinate()
          ? get_default_serial_new_size_node(geo, dist_read_par_limit())
          : new_size_node_;
  const Long file_size =
      serial_read_field_par(f, path, new_size_node, -data_size, SEEK_END);
  if (file_size != data_size) {
    displayln_info(
        fname +
        ssprintf(": file size do not match read_size=%ld vs data_size=%ld",
                 file_size, data_size));
    qassert(false);
  }
  const crc32_t f_crc = field_crc32(f);
  const bool is_checking = is_checksum_mismatch();
  is_checksum_mismatch() = false;
  if (crc != f_crc) {
    displayln_info(fname + ssprintf(": WARNING: file crc32 do not match "
                                    "read_info_crc=%06X vs read_data_crc=%06X",
                                    crc, f_crc));
    is_checksum_mismatch() = true;
  }
  if (is_checking) {
    qassert(is_checksum_mismatch() == false);
  }
  timer.flops += file_size;
  return file_size;
}

template <class M>
Long write_field_float_from_double(
    const Field<M>& f, const std::string& path,
    const Coordinate& new_size_node = Coordinate())
// interface_function
{
  TIMER_VERBOSE_FLOPS("write_field_float_from_double");
  Field<RealF> ff;
  convert_field_float_from_double(ff, f);
  to_from_big_endian(get_data(ff));
  const Long total_bytes = write_field(ff, path, new_size_node);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
Long read_field_double_from_float(
    Field<M>& f, const std::string& path,
    const Coordinate& new_size_node = Coordinate())
// interface_function
{
  TIMER_VERBOSE_FLOPS("read_field_double_from_float");
  Field<RealF> ff;
  const Long total_bytes = read_field(ff, path, new_size_node);
  if (total_bytes == 0) {
    return 0;
  } else {
    to_from_big_endian(get_data(ff));
    convert_field_double_from_float(f, ff);
    timer.flops += total_bytes;
    return total_bytes;
  }
}

template <class M>
Long write_field_64(const Field<M>& f, const std::string& path,
                    const Coordinate& new_size_node = Coordinate())
// interface_function
{
  TIMER_VERBOSE_FLOPS("write_field_64");
  Field<M> ff;
  ff.init(f);
  to_from_big_endian(get_data(ff));
  const Long total_bytes = write_field(ff, path, new_size_node);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
Long read_field_64(Field<M>& f, const std::string& path,
                   const Coordinate& new_size_node = Coordinate())
// interface_function
{
  TIMER_VERBOSE_FLOPS("read_field_64");
  const Long total_bytes = read_field(f, path, new_size_node);
  if (total_bytes == 0) {
    return 0;
  } else {
    to_from_big_endian(get_data(f));
    timer.flops += total_bytes;
    return total_bytes;
  }
}

template <class M>
Long write_field_double(const Field<M>& f, const std::string& path,
                        const Coordinate& new_size_node = Coordinate())
// interface_function
{
  TIMER_VERBOSE_FLOPS("write_field_double");
  Field<M> ff;
  ff.init(f);
  to_from_big_endian(get_data(ff));
  const Long total_bytes = write_field(ff, path, new_size_node);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
Long read_field_double(Field<M>& f, const std::string& path,
                       const Coordinate& new_size_node = Coordinate())
// interface_function
{
  TIMER_VERBOSE_FLOPS("read_field_double");
  const Long total_bytes = read_field(f, path, new_size_node);
  if (total_bytes == 0) {
    return 0;
  } else {
    to_from_big_endian(get_data(f));
    timer.flops += total_bytes;
    return total_bytes;
  }
}

// --------------------

bool is_dist_field(const std::string& path);

bool is_field(const std::string& path);

bool is_d_field(const std::string& path);

bool dist_repartition(const Coordinate& new_size_node, const std::string& path,
                      const std::string& new_path = "");

// --------------------

#ifdef QLAT_INSTANTIATE_FIELD_IO
#define QLAT_EXTERN
#else
#define QLAT_EXTERN extern
#endif

#define QLAT_EXTERN_TEMPLATE(TYPENAME)                           \
                                                                 \
  QLAT_EXTERN template crc32_t field_crc32<TYPENAME>(            \
      const Field<TYPENAME>& f);                                 \
                                                                 \
  QLAT_EXTERN template Long dist_write_field<TYPENAME>(          \
      const Field<TYPENAME>& f, const Coordinate& new_size_node, \
      const std::string& path);                                  \
                                                                 \
  QLAT_EXTERN template Long dist_read_field<TYPENAME>(           \
      Field<TYPENAME> & f, const std::string& path);             \
                                                                 \
  QLAT_EXTERN template Long write_field<TYPENAME>(               \
      const Field<TYPENAME>& f, const std::string& path,         \
      const Coordinate& new_size_node);                          \
                                                                 \
  QLAT_EXTERN template Long read_field<TYPENAME>(                \
      Field<TYPENAME> & f, const std::string& path,              \
      const Coordinate& new_size_node)

QLAT_CALL_WITH_TYPES(QLAT_EXTERN_TEMPLATE);
#undef QLAT_EXTERN_TEMPLATE

#undef QLAT_EXTERN

}  // namespace qlat
