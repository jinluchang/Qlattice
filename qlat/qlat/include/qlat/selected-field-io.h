#pragma once

#include <qlat/field-serial-io.h>
#include <qlat/selected-field.h>

namespace qlat
{  //

std::string make_selected_field_header(const Geometry& geo,
                                       const long n_per_tslice,
                                       const int sizeof_M, const crc32_t crc32);

long read_selected_geo_info(Coordinate& total_site, int& multiplicity,
                            long& n_per_tslice, int& sizeof_M, crc32_t& crc,
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
