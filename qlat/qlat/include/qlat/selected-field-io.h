#pragma once

#include <qlat/field-io.h>
#include <qlat/selected-field.h>

namespace qlat
{  //

std::string make_selected_field_header(const Geometry& geo, Int multiplicity,
                                       const Int sizeof_M, const crc32_t crc32);

Long read_selected_geo_info(Coordinate& total_site, Int& multiplicity,
                            Int& sizeof_M, crc32_t& crc,
                            const std::string& path);

bool is_selected_field(const std::string& path);

// ---------------------------------------

template <class M>
crc32_t field_crc32(const SelectedField<M>& sf, const FieldSelection& fsel,
                    const Coordinate& new_size_node_ = Coordinate())
{
  TIMER_VERBOSE_FLOPS("field_crc32(sf)");
  const Geometry& geo = sf.geo();
  qassert(geo.is_only_local);
  qassert(fsel.f_rank.geo() == geo);
  // const Coordinate total_site = geo.total_site();
  const Coordinate new_size_node = new_size_node_ != Coordinate()
                                       ? new_size_node_
                                       : get_default_serial_new_size_node(geo);
  qassert(new_size_node[0] == 1);
  qassert(new_size_node[1] == 1);
  qassert(new_size_node[2] == 1);
  std::vector<FieldSelection> fsels;
  const ShufflePlan sp = make_shuffle_plan(fsels, fsel, new_size_node);
  std::vector<SelectedField<M>> sfs;
  shuffle_field(sfs, sf, sp);
  qassert(fsels.size() == sfs.size());
  const Int new_num_node = product(new_size_node);
  crc32_t crc = 0;
  for (Int i = 0; i < (int)sfs.size(); ++i) {
    const Int new_id_node = sfs[i].geo().geon.id_node;
    qassert(sfs[i].geo().geon.num_node == new_num_node);
    const Vector<M> v = get_data(sfs[i].field);
    crc ^= crc32_shift(crc32_par(v),
                       (new_num_node - new_id_node - 1) * v.data_size());
  }
  glb_sum(get_data_char(crc));
  return crc;
}

template <class M>
Long write_selected_field(const SelectedField<M>& sf, const std::string& path,
                          const FieldSelection& fsel,
                          const Coordinate& new_size_node_ = Coordinate())
// interface function
{
  TIMER_VERBOSE_FLOPS("write_selected_field");
  displayln_info(fname + ssprintf(": fn='%s'.", path.c_str()));
  const Geometry& geo = sf.geo();
  const Int multiplicity = sf.multiplicity;
  qassert(geo.is_only_local);
  qassert(fsel.f_rank.geo() == geo);
  const Coordinate new_size_node =
      new_size_node_ != Coordinate()
          ? new_size_node_
          : get_default_serial_new_size_node(geo, dist_write_par_limit());
  qassert(new_size_node[0] == 1);
  qassert(new_size_node[1] == 1);
  qassert(new_size_node[2] == 1);
  std::vector<FieldSelection> fsels;
  const ShufflePlan sp = make_shuffle_plan(fsels, fsel, new_size_node);
  std::vector<SelectedField<M>> sfs;
  shuffle_field(sfs, sf, sp);
  qassert(fsels.size() == sfs.size());
  const Int new_num_node = product(new_size_node);
  std::vector<Long> n_elems_vec(new_num_node, 0);
  for (size_t i = 0; i < sfs.size(); ++i) {
    const Int id_node = fsels[i].f_rank.geo().geon.id_node;
    n_elems_vec[id_node] = fsels[i].n_elems;
    sfs[i].init(fsels[i], multiplicity);
  }
  glb_sum(n_elems_vec);
  std::vector<Long> data_offset_vec(new_num_node + 1, 0);
  Long total_bytes = 0;
  for (Int i = 0; i < new_num_node; ++i) {
    data_offset_vec[i] = total_bytes;
    total_bytes += n_elems_vec[i] * multiplicity * sizeof(M);
  }
  data_offset_vec[new_num_node] = total_bytes;
  crc32_t crc = 0;
  for (Int i = 0; i < (int)sfs.size(); ++i) {
    const Int new_id_node = sfs[i].geo().geon.id_node;
    qassert(sfs[i].geo().geon.num_node == new_num_node);
    const Vector<M> v = get_data(sfs[i].field);
    crc ^= crc32_shift(crc32_par(v),
                       (total_bytes - data_offset_vec[new_id_node + 1]));
  }
  glb_sum(get_data_char(crc));
  if (get_force_field_write_sizeof_M() == 0) {
    qtouch_info(path + ".partial",
                make_selected_field_header(geo, multiplicity, sizeof(M), crc));
  } else {
    const Int sizeof_M = get_force_field_write_sizeof_M();
    qassert((multiplicity * sizeof(M)) % sizeof_M == 0);
    const Int multiplicity_new = (multiplicity * sizeof(M)) / sizeof_M;
    qtouch_info(path + ".partial",
                make_selected_field_header(geo, multiplicity_new, sizeof_M, crc));
    get_force_field_write_sizeof_M() = 0;
  }
  const Int mpi_tag = 8;
  if (get_id_node() == 0) {
    qassert(sfs.size() > 0);
    QFile fp = qfopen(path + ".partial", "a");
    qassert(not fp.null());
    const Int num_node = get_num_node();
    for (Int new_id_node = 0; new_id_node < new_num_node; ++new_id_node) {
      vector<M> vec(n_elems_vec[new_id_node] * multiplicity);
      Vector<M> v = get_data(vec);
      const Int id_node =
          get_id_node_from_new_id_node(new_id_node, new_num_node, num_node);
      if (0 == id_node) {
        assign(v, get_data(sfs[new_id_node].field));
      } else {
        mpi_recv(v.data(), v.data_size(), MPI_BYTE, id_node, mpi_tag,
                 get_comm());
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
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
Long read_selected_field(SelectedField<M>& sf, const std::string& path,
                         const FieldSelection& fsel,
                         const Coordinate& new_size_node_ = Coordinate())
// interface function
{
  TIMER_VERBOSE_FLOPS("read_selected_field");
  displayln_info(fname + ssprintf(": fn='%s'.", path.c_str()));
  sf.init();
  Coordinate total_site;
  Int multiplicity = 0;
  Int sizeof_M = 0;
  crc32_t crc_info = 0;
  const Long pos = read_selected_geo_info(total_site, multiplicity, sizeof_M,
                                          crc_info, path);
  if (total_site == Coordinate() or multiplicity == 0) {
    displayln_info(fname +
                   ssprintf(": fn='%s' can not be parsed.", path.c_str()));
    return 0;
  }
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
  geo.init(total_site);
  qassert(fsel.f_rank.geo() == geo);
  const Coordinate new_size_node =
      new_size_node_ != Coordinate()
          ? new_size_node_
          : get_default_serial_new_size_node(geo, dist_read_par_limit());
  qassert(new_size_node[0] == 1);
  qassert(new_size_node[1] == 1);
  qassert(new_size_node[2] == 1);
  std::vector<FieldSelection> fsels;
  const ShufflePlan sp = make_shuffle_plan(fsels, fsel, new_size_node);
  if (not is_no_shuffle(sp)) {
    qassert(fsels.size() == sp.geos_recv.size());
  }
  std::vector<SelectedField<M>> sfs;
  sfs.resize(fsels.size());
  const Int new_num_node = product(new_size_node);
  std::vector<Long> n_elems_vec(new_num_node, 0);
  for (size_t i = 0; i < sfs.size(); ++i) {
    const Int id_node = fsels[i].f_rank.geo().geon.id_node;
    n_elems_vec[id_node] = fsels[i].n_elems;
    sfs[i].init(fsels[i], multiplicity);
  }
  glb_sum(n_elems_vec);
  std::vector<Long> data_offset_vec(new_num_node + 1, 0);
  Long total_bytes = 0;
  for (Int i = 0; i < new_num_node; ++i) {
    data_offset_vec[i] = total_bytes;
    total_bytes += n_elems_vec[i] * multiplicity * sizeof(M);
  }
  data_offset_vec[new_num_node] = total_bytes;
  crc32_t crc = 0;
  if (sfs.size() > 0) {
    QFile fp = qfopen(path, "r");
    qassert(not fp.null());
    qfseek(fp, pos + data_offset_vec[sfs[0].geo().geon.id_node], SEEK_SET);
    for (Int i = 0; i < (int)sfs.size(); ++i) {
      const Int new_id_node = sfs[i].geo().geon.id_node;
      qassert(sfs[i].geo().geon.num_node == new_num_node);
      const Vector<M> v = get_data(sfs[i].field);
      qread_data(v, fp);
      crc ^= crc32_shift(crc32_par(v),
                         (total_bytes - data_offset_vec[new_id_node + 1]));
    }
    qfclose(fp);
  }
  glb_sum(get_data_char(crc));
  if (crc != crc_info) {
    displayln_info(
        fname +
        ssprintf(": crc of data = %08X ; crc of header = %08X", crc, crc_info));
    qassert(false);
  }
  sf.init(fsel, multiplicity);
  shuffle_field_back(sf, sfs, sp);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
Long write_selected_field_64(const SelectedField<M>& f, const std::string& path,
                             const FieldSelection& fsel,
                             const Coordinate& new_size_node_ = Coordinate())
{
  TIMER_VERBOSE_FLOPS("write_selected_field_64");
  SelectedField<M> ff;
  ff = f;
  to_from_big_endian(get_data(ff));
  const Long total_bytes = write_selected_field(ff, path, fsel, new_size_node_);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
Long read_selected_field_64(SelectedField<M>& sf, const std::string& path,
                            const FieldSelection& fsel,
                            const Coordinate& new_size_node_ = Coordinate())
{
  TIMER_VERBOSE_FLOPS("read_selected_field_64");
  const Long total_bytes = read_selected_field(sf, path, fsel, new_size_node_);
  if (total_bytes == 0) {
    return 0;
  } else {
    to_from_big_endian(get_data(sf));
    timer.flops += total_bytes;
    return total_bytes;
  }
}

template <class M>
Long write_selected_field_double(
    const SelectedField<M>& f, const std::string& path,
    const FieldSelection& fsel, const Coordinate& new_size_node_ = Coordinate())
{
  TIMER_VERBOSE_FLOPS("write_selected_field_double");
  SelectedField<M> ff;
  ff = f;
  to_from_big_endian(get_data(ff));
  const Long total_bytes = write_selected_field(ff, path, fsel, new_size_node_);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
Long read_selected_field_double(SelectedField<M>& sf, const std::string& path,
                                const FieldSelection& fsel,
                                const Coordinate& new_size_node_ = Coordinate())
{
  TIMER_VERBOSE_FLOPS("read_selected_field_double");
  const Long total_bytes = read_selected_field(sf, path, fsel, new_size_node_);
  if (total_bytes == 0) {
    return 0;
  } else {
    to_from_big_endian(get_data(sf));
    timer.flops += total_bytes;
    return total_bytes;
  }
}

template <class M>
Long write_selected_field_float_from_double(
    const SelectedField<M>& f, const std::string& path,
    const FieldSelection& fsel, const Coordinate& new_size_node_ = Coordinate())
{
  TIMER_VERBOSE_FLOPS("write_selected_field_float_from_double");
  SelectedField<float> ff;
  convert_field_float_from_double(ff, f);
  to_from_big_endian(get_data(ff));
  const Long total_bytes = write_selected_field(ff, path, fsel, new_size_node_);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
Long read_selected_field_double_from_float(
    SelectedField<M>& f, const std::string& path, const FieldSelection& fsel,
    const Coordinate& new_size_node_ = Coordinate())
{
  TIMER_VERBOSE_FLOPS("read_selected_field_double_from_float");
  SelectedField<float> ff;
  const Long total_bytes = read_selected_field(ff, path, fsel, new_size_node_);
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
Long write_selected_field(const Field<M>& f, const std::string& path,
                          const FieldSelection& fsel,
                          const Coordinate& new_size_node_ = Coordinate())
{
  TIMER_VERBOSE_FLOPS("write_selected_field(f)");
  SelectedField<M> sf;
  set_selected_field(sf, f, fsel);
  return write_selected_field(sf, path, fsel, new_size_node_);
}

template <class M>
Long write_selected_field_64(const Field<M>& f, const std::string& path,
                             const FieldSelection& fsel,
                             const Coordinate& new_size_node_ = Coordinate())
{
  TIMER_VERBOSE_FLOPS("write_selected_field_64(f)");
  SelectedField<M> sf;
  set_selected_field(sf, f, fsel);
  return write_selected_field_64(sf, path, fsel, new_size_node_);
}

template <class M>
Long write_selected_field_double(
    const Field<M>& f, const std::string& path, const FieldSelection& fsel,
    const Coordinate& new_size_node_ = Coordinate())
{
  TIMER_VERBOSE_FLOPS("write_selected_field_64(f)");
  SelectedField<M> sf;
  set_selected_field(sf, f, fsel);
  return write_selected_field_double(sf, path, fsel, new_size_node_);
}

template <class M>
Long write_selected_field_float_from_double(
    const Field<M>& f, const std::string& path, const FieldSelection& fsel,
    const Coordinate& new_size_node_ = Coordinate())
{
  TIMER_VERBOSE_FLOPS("write_selected_field_64(f)");
  SelectedField<M> sf;
  set_selected_field(sf, f, fsel);
  return write_selected_field_float_from_double(sf, path, fsel, new_size_node_);
}

template <class M>
Long read_selected_field(Field<M>& f, const std::string& path,
                         const FieldSelection& fsel,
                         const Coordinate& new_size_node_ = Coordinate())
{
  TIMER_VERBOSE_FLOPS("read_selected_field(f)");
  SelectedField<M> sf;
  const Long total_bytes = read_selected_field(sf, path, fsel, new_size_node_);
  if (total_bytes > 0) {
    set_field_selected(f, sf, fsel);
  }
  return total_bytes;
}

template <class M>
Long read_selected_field_64(Field<M>& f, const std::string& path,
                            const FieldSelection& fsel,
                            const Coordinate& new_size_node_ = Coordinate())
{
  TIMER_VERBOSE_FLOPS("read_selected_field_64(f)");
  SelectedField<M> sf;
  const Long total_bytes =
      read_selected_field_64(sf, path, fsel, new_size_node_);
  if (total_bytes > 0) {
    set_field_selected(f, sf, fsel);
  }
  return total_bytes;
}

template <class M>
Long read_selected_field_double(Field<M>& f, const std::string& path,
                                const FieldSelection& fsel,
                                const Coordinate& new_size_node_ = Coordinate())
{
  TIMER_VERBOSE_FLOPS("read_selected_field_double(f)");
  SelectedField<M> sf;
  const Long total_bytes =
      read_selected_field_double(sf, path, fsel, new_size_node_);
  if (total_bytes > 0) {
    set_field_selected(f, sf, fsel);
  }
  return total_bytes;
}

template <class M>
Long read_selected_field_double_from_float(
    Field<M>& f, const std::string& path, const FieldSelection& fsel,
    const Coordinate& new_size_node_ = Coordinate())
{
  TIMER_VERBOSE_FLOPS("read_selected_field_double_from_float(f)");
  SelectedField<M> sf;
  const Long total_bytes =
      read_selected_field_double_from_float(sf, path, fsel, new_size_node_);
  if (total_bytes > 0) {
    set_field_selected(f, sf, fsel);
  }
  return total_bytes;
}

// --------------------

#ifdef QLAT_INSTANTIATE_SELECTED_FIELD
#define QLAT_EXTERN
#else
#define QLAT_EXTERN extern
#endif

#define QLAT_EXTERN_TEMPLATE(TYPENAME)                                  \
                                                                        \
  QLAT_EXTERN template SelectedField<TYPENAME>& operator+= <TYPENAME>(  \
      SelectedField<TYPENAME> & f, const SelectedField<TYPENAME>& f1);  \
                                                                        \
  QLAT_EXTERN template SelectedField<TYPENAME>& operator-= <TYPENAME>(  \
      SelectedField<TYPENAME> & f, const SelectedField<TYPENAME>& f1);  \
                                                                        \
  QLAT_EXTERN template SelectedField<TYPENAME>& operator*=              \
      <TYPENAME>(SelectedField<TYPENAME> & f, const double factor);     \
                                                                        \
  QLAT_EXTERN template SelectedField<TYPENAME>& operator*=              \
      <TYPENAME>(SelectedField<TYPENAME> & f, const ComplexD factor);   \
                                                                        \
  QLAT_EXTERN template void only_keep_selected_points<TYPENAME>(        \
      Field<TYPENAME> & f, const FieldSelection& fsel);                 \
                                                                        \
  QLAT_EXTERN template double qnorm<TYPENAME>(                          \
      const SelectedField<TYPENAME>& sp);                               \
                                                                        \
  QLAT_EXTERN template void qnorm_field<TYPENAME>(                      \
      SelectedField<double> & sp, const SelectedField<TYPENAME>& sp1);  \
                                                                        \
  QLAT_EXTERN template void set_selected_field<TYPENAME>(               \
      SelectedField<TYPENAME> & sf, const Field<TYPENAME>& f,           \
      const FieldSelection& fsel);                                      \
                                                                        \
  QLAT_EXTERN template void set_selected_field<TYPENAME>(               \
      SelectedField<TYPENAME> & sf, const SelectedField<TYPENAME>& sf0, \
      const FieldSelection& fsel, const FieldSelection& fsel0,          \
      const bool is_keeping_data);                                      \
                                                                        \
  QLAT_EXTERN template void set_selected_field<TYPENAME>(               \
      SelectedField<TYPENAME> & sf, const SelectedPoints<TYPENAME>& sp, \
      const FieldSelection& fsel, const PointsSelection& psel,          \
      const bool is_keeping_data);                                      \
                                                                        \
  QLAT_EXTERN template void set_selected_points<TYPENAME>(              \
      SelectedPoints<TYPENAME> & sp, const SelectedField<TYPENAME>& sf, \
      const PointsSelection& psel, const FieldSelection& fsel,          \
      const bool is_keeping_data);                                      \
                                                                        \
  QLAT_EXTERN template void set_field_selected<TYPENAME>(               \
      Field<TYPENAME> & f, const SelectedField<TYPENAME>& sf,           \
      const FieldSelection& fsel, const bool is_keeping_data);          \
                                                                        \
  QLAT_EXTERN template void acc_field<TYPENAME>(                        \
      Field<TYPENAME> & f, const SelectedField<TYPENAME>& sf,           \
      const FieldSelection& fsel);                                      \
                                                                        \
  QLAT_EXTERN template void field_glb_sum_tslice<TYPENAME>(             \
      SelectedPoints<TYPENAME> & sp, const SelectedField<TYPENAME>& sf, \
      const FieldSelection& fsel, const Int t_dir);                     \
                                                                        \
  QLAT_EXTERN template crc32_t field_crc32<TYPENAME>(                   \
      const SelectedField<TYPENAME>& sf, const FieldSelection& fsel,    \
      const Coordinate& new_size_node_);                                \
                                                                        \
  QLAT_EXTERN template Long write_selected_field<TYPENAME>(             \
      const SelectedField<TYPENAME>& sf, const std::string& path,       \
      const FieldSelection& fsel, const Coordinate& new_size_node_);    \
                                                                        \
  QLAT_EXTERN template Long read_selected_field<TYPENAME>(              \
      SelectedField<TYPENAME> & sf, const std::string& path,            \
      const FieldSelection& fsel, const Coordinate& new_size_node_);    \
                                                                        \
  QLAT_EXTERN template void set_u_rand<TYPENAME>(                       \
      SelectedField<TYPENAME> & sf, const FieldSelection& fsel,         \
      const RngState& rs, const RealD upper, const RealD lower);        \
                                                                        \
  QLAT_EXTERN template void set_g_rand<TYPENAME>(                       \
      SelectedField<TYPENAME> & sf, const FieldSelection& fsel,         \
      const RngState& rs, const RealD center, const RealD sigma)

QLAT_CALL_WITH_TYPES(QLAT_EXTERN_TEMPLATE);
#undef QLAT_EXTERN_TEMPLATE

#undef QLAT_EXTERN

}  // namespace qlat
