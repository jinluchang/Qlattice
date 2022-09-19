#pragma once

#include <qlat/field-dist-io.h>
#include <qlat/field-io.h>

#include <stdio.h>
#include <ctime>

#include <fstream>
#include <iostream>

namespace qlat
{  //

inline Coordinate get_default_serial_new_size_node(const Geometry& geo, const int max_num_ = 0)
{
  const int num_node = geo.geon.num_node;
  const int max_num = max_num_ <= 0 or max_num_ > num_node ? num_node : max_num_;
  const Coordinate total_site = geo.total_site();
  Coordinate new_size_node = Coordinate(1, 1, 1, total_site[3]);
  while (max_num < new_size_node[3]) {
    if (new_size_node[3] % 2 == 0) {
      new_size_node[3] /= 2;
    } else if (new_size_node[3] % 3 == 0) {
      new_size_node[3] /= 3;
    } else if (new_size_node[3] % 5 == 0) {
      new_size_node[3] /= 5;
    } else if (new_size_node[3] % 7 == 0) {
      new_size_node[3] /= 7;
    } else if (new_size_node[3] % 11 == 0) {
      new_size_node[3] /= 11;
    } else if (new_size_node[3] % 13 == 0) {
      new_size_node[3] /= 13;
    } else {
      new_size_node[3] = 1;
    }
  }
  qassert(total_site % new_size_node == Coordinate());
  return new_size_node;
}

template <class M>
long serial_write_field(const Field<M>& f, const std::string& path,
                        const Coordinate& new_size_node)
// will append to the file
// assume new_size_node is properly chosen so that concatenate the new fields
// would be correct. eg. new_size_node = Coordinate(1,1,1,2)
{
  TIMER_VERBOSE_FLOPS("serial_write_field");
  std::vector<Field<M> > fs;
  shuffle_field(fs, f, new_size_node);
  const int mpi_tag = 6;
  if (get_id_node() == 0) {
    qassert(fs.size() > 0);
    Field<M> f;
    f = fs[0];
    Vector<M> v = get_data(f);
    const int num_node = get_num_node();
    const int new_num_node = product(new_size_node);
    QFile qfile = qfopen(path, "a");
    qassert(not qfile.null());
    for (int new_id_node = 0; new_id_node < new_num_node; ++new_id_node) {
      const int id_node =
          get_id_node_from_new_id_node(new_id_node, new_num_node, num_node);
      if (0 == id_node) {
        assign(v, get_data(fs[new_id_node]));
      } else {
        mpi_recv(v.data(), v.data_size(), MPI_BYTE, id_node, mpi_tag,
                 get_comm(), MPI_STATUS_IGNORE);
      }
      qwrite_data(v, qfile);
    }
    qfile.close();
  } else {
    for (size_t i = 0; i < fs.size(); ++i) {
      const Vector<M> v = get_data(fs[i]);
      mpi_send((void*)v.data(), v.data_size(), MPI_BYTE, 0, mpi_tag,
               get_comm());
    }
  }
  const long file_size = get_data(f).data_size() * f.geo().geon.num_node;
  timer.flops += file_size;
  return file_size;
}

template <class M>
long serial_read_field(Field<M>& f, const std::string& path,
                       const Coordinate& new_size_node, const long offset = 0,
                       const int whence = SEEK_SET)
// will read from offset relative to whence
// assume new_size_node is properly choosen so that concatenate the new fields
// would be correct. eg. new_size_node = Coordinate(1,1,1,2)
{
  TIMER_VERBOSE_FLOPS("serial_read_field");
  if (not does_file_exist_qar_sync_node(path)) {
    displayln_info(fname +
                   ssprintf(": file does not exist: '%s'", path.c_str()));
    return 0;
  }
  const Geometry& geo = f.geo();
  std::vector<Field<M> > fs;
  const std::vector<Geometry> new_geos =
      make_dist_io_geos(geo.total_site(), geo.multiplicity, new_size_node);
  fs.resize(new_geos.size());
  for (size_t i = 0; i < fs.size(); ++i) {
    fs[i].init(new_geos[i]);
  }
  const int mpi_tag = 7;
  if (get_id_node() == 0) {
    qassert(fs.size() > 0);
    Field<M> f;
    f.init(fs[0].geo());
    Vector<M> v = get_data(f);
    const int num_node = get_num_node();
    const int new_num_node = product(new_size_node);
    QFile qfile = qfopen(path, "r");
    qassert(not qfile.null());
    qfseek(qfile, offset, whence);
    for (int new_id_node = 0; new_id_node < new_num_node; ++new_id_node) {
      const int id_node =
          get_id_node_from_new_id_node(new_id_node, new_num_node, num_node);
      qread_data(v, qfile);
      if (0 == id_node) {
        assign(get_data(fs[new_id_node]), v);
      } else {
        mpi_send((void*)v.data(), v.data_size(), MPI_BYTE, id_node, mpi_tag,
                 get_comm());
      }
    }
    qfile.close();
  } else {
    for (size_t i = 0; i < fs.size(); ++i) {
      Vector<M> v = get_data(fs[i]);
      mpi_recv(v.data(), v.data_size(), MPI_BYTE, 0, mpi_tag, get_comm(),
               MPI_STATUS_IGNORE);
    }
  }
  shuffle_field_back(f, fs, new_size_node);
  const long file_size = get_data(f).data_size() * f.geo().geon.num_node;
  timer.flops += file_size;
  return file_size;
}

template <class M>
long serial_read_field_par(Field<M>& f, const std::string& path,
                           const Coordinate& new_size_node,
                           const long offset = 0, const int whence = SEEK_SET)
// will read from offset relative to whence
// assume new_size_node is properly choosen so that concatenate the new fields
// would be correct. eg. new_size_node = Coordinate(1,1,1,2)
{
  TIMER_VERBOSE_FLOPS("serial_read_field_par");
  if (not does_file_exist_qar_sync_node(path)) {
    displayln_info(fname +
                   ssprintf(": file does not exist: '%s'", path.c_str()));
    return 0;
  }
  const Geometry& geo = f.geo();
  std::vector<Field<M> > fs;
  const std::vector<Geometry> new_geos =
      make_dist_io_geos(geo.total_site(), geo.multiplicity, new_size_node);
  fs.resize(new_geos.size());
  for (size_t i = 0; i < fs.size(); ++i) {
    fs[i].init(new_geos[i]);
  }
  if (fs.size() > 0) {
    QFile qfile = qfopen(path, "r");
    qassert(not qfile.null());
    qfseek(qfile,
           offset + fs[0].geo().geon.id_node * get_data(fs[0]).data_size(),
           whence);
    for (size_t i = 0; i < fs.size(); ++i) {
      Vector<M> v = get_data(fs[i]);
      qread_data(v, qfile);
    }
    qfile.close();
  }
  shuffle_field_back(f, fs, new_size_node);
  sync_node();
  const long file_size = get_data(f).data_size() * f.geo().geon.num_node;
  timer.flops += file_size;
  return file_size;
}

template <class M>
long serial_write_field(const Field<M>& f, const std::string& path)
// interface_function
{
  return serial_write_field(
      f, path, get_default_serial_new_size_node(f.geo(), dist_write_par_limit()));
}

template <class M>
long serial_read_field(Field<M>& f, const std::string& path,
                       const long offset = 0, const int whence = SEEK_SET)
// interface_function
{
  return serial_read_field(
      f, path, get_default_serial_new_size_node(f.geo(), dist_read_par_limit()),
      offset, whence);
}

template <class M>
long serial_read_field_par(Field<M>& f, const std::string& path,
                           const long offset = 0, const int whence = SEEK_SET)
// interface_function
{
  return serial_read_field_par(
      f, path, get_default_serial_new_size_node(f.geo(), dist_read_par_limit()),
      offset, whence);
}

template <class M>
crc32_t field_simple_checksum(const Field<M>& f)
{
  TIMER("field_simple_checksum");
  qassert(f.geo().is_only_local());
  crc32_t ret = 0;
  const Vector<M> v = get_data(f);
  Vector<crc32_t> vc((crc32_t*)v.data(), v.data_size() / sizeof(crc32_t));
  for (long i = 0; i < vc.size(); ++i) {
    ret += vc[i];
  }
  long sum = ret;
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
  std::vector<Field<M> > fs;
  shuffle_field(fs, f, new_size_node);
  crc32_t ret = 0;
  for (int i = 0; i < (int)fs.size(); ++i) {
    const int new_id_node = fs[i].geo().geon.id_node;
    const int new_num_node = fs[i].geo().geon.num_node;
    const Vector<M> v = get_data(fs[i]);
    ret ^= crc32_shift(crc32_par(v),
                       (new_num_node - new_id_node - 1) * v.data_size());
  }
  glb_sum_byte(ret);
  timer.flops += get_data(f).data_size() * geo.geon.num_node;
  return ret;
}

template <class M>
crc32_t field_crc32_sites(const Field<M>& f)
{
  TIMER_VERBOSE_FLOPS("field_crc32_sites");
  const Geometry& geo = f.geo();
  const long total_volume = geo.total_volume();
  const long data_size_site = geo.multiplicity * sizeof(M);
  const int v_limit = omp_get_max_threads();
  std::vector<crc32_t> crcs(v_limit, 0);
#pragma omp parallel
  {
    crc32_t crc = 0;
#pragma omp for
    for (long index = 0; index < geo.local_volume(); ++index) {
      const Coordinate xl = geo.coordinate_from_index(index);
      const Coordinate xg = geo.coordinate_g_from_l(xl);
      const long gindex = geo.g_index_from_g_coordinate(xg);
      const long offset = data_size_site * (total_volume - gindex - 1);
      const Vector<M> v = f.get_elems_const(xl);
      crc ^= crc32_shift(crc32(v), offset);
    }
    const int id = omp_get_thread_num();
    crcs[id] = crc;
  }
  crc32_t ret = 0;
  for (int i = 0; i < v_limit; ++i) {
    ret ^= crcs[i];
  }
  glb_sum_byte(ret);
  timer.flops += get_data(f).data_size() * geo.geon.num_node;
  return ret;
}

template <class M>
crc32_t field_crc32(const Field<M>& f)
{
  return field_crc32_shuffle(f);
}

inline std::string make_field_header(const Geometry& geo, const int sizeof_M,
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
  out << "multiplicity = " << geo.multiplicity << std::endl;
  out << "sizeof(M) = " << sizeof_M << std::endl;
  out << ssprintf("field_crc32 = %08X", crc32) << std::endl;
  out << "END_HEADER" << std::endl;
  return out.str();
}

template <class M>
long write_field(const Field<M>& f, const std::string& path,
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
    qtouch_info(path + ".partial", make_field_header(geo, sizeof(M), crc32));
  } else {
    const int sizeof_M = get_force_field_write_sizeof_M();
    qassert((geo.multiplicity * sizeof(M)) % sizeof_M == 0);
    const int multiplicity = (geo.multiplicity * sizeof(M)) / sizeof_M;
    qtouch_info(
        path + ".partial",
        make_field_header(geo_remult(geo, multiplicity), sizeof_M, crc32));
    get_force_field_write_sizeof_M() = 0;
  }
  const long file_size = serial_write_field(
      f, path + ".partial",
      get_default_serial_new_size_node(geo, dist_write_par_limit()));
  qrename_info(path + ".partial", path);
  timer.flops += file_size;
  return file_size;
}

inline void read_geo_info(Coordinate& total_site, int& multiplicity, int& sizeof_M,
                          crc32_t& crc, const std::string& path)
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
          for (int m = 0; m < 4; ++m) {
            reads(total_site[m],
                  info_get_prop(infos, ssprintf("total_site[%d] = ", m)));
          }
          reads(multiplicity, info_get_prop(infos, "multiplicity = "));
          reads(sizeof_M, info_get_prop(infos, "sizeof(M) = "));
          crc = read_crc32(info_get_prop(infos, "field_crc32 = "));
        }
      }
    }
    qfile.close();
  }
  bcast(Vector<Coordinate>(&total_site, 1));
  bcast(Vector<int>(&multiplicity, 1));
  bcast(Vector<int>(&sizeof_M, 1));
  bcast(Vector<crc32_t>(&crc, 1));
}

template <class M>
long read_field(Field<M>& f, const std::string& path,
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
  int multiplicity = 0;
  int sizeof_M = 0;
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
  geo.init(total_site, multiplicity);
  f.init(geo, geo.multiplicity);
  const long data_size =
      geo.geon.num_node * geo.local_volume() * geo.multiplicity * sizeof(M);
  const Coordinate new_size_node =
      new_size_node_ == Coordinate()
          ? get_default_serial_new_size_node(geo, dist_read_par_limit())
          : new_size_node_;
  const long file_size =
      serial_read_field_par(f, path, new_size_node, -data_size, SEEK_END);
  if (file_size != data_size) {
    displayln_info(
        fname +
        ssprintf(": file size do not match read_size=%ld vs data_size=%ld",
                 file_size, data_size));
    qassert(false);
  }
  const crc32_t f_crc = field_crc32(f);
  is_checksum_missmatch() = false;
  if (crc != f_crc) {
    displayln_info(fname + ssprintf(": WARNING: file crc32 do not match "
                                    "read_info_crc=%06X vs read_data_crc=%06X",
                                    crc, f_crc));
    is_checksum_missmatch() = true;
  }
  timer.flops += file_size;
  return file_size;
}

template <class M>
long write_field_float_from_double(
    const Field<M>& f, const std::string& path,
    const Coordinate& new_size_node = Coordinate())
// interface_function
{
  TIMER_VERBOSE_FLOPS("write_field_float_from_double");
  Field<float> ff;
  convert_field_float_from_double(ff, f);
  to_from_big_endian_32(get_data(ff));
  const long total_bytes = write_field(ff, path, new_size_node);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
long read_field_double_from_float(
    Field<M>& f, const std::string& path,
    const Coordinate& new_size_node = Coordinate())
// interface_function
{
  TIMER_VERBOSE_FLOPS("read_field_double_from_float");
  Field<float> ff;
  const long total_bytes = read_field(ff, path, new_size_node);
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
long write_field_64(const Field<M>& f, const std::string& path,
                    const Coordinate& new_size_node = Coordinate())
// interface_function
{
  TIMER_VERBOSE_FLOPS("write_field_64");
  Field<M> ff;
  ff.init(f);
  to_from_big_endian_64(get_data(ff));
  const long total_bytes = write_field(ff, path, new_size_node);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
long read_field_64(Field<M>& f, const std::string& path,
                   const Coordinate& new_size_node = Coordinate())
// interface_function
{
  TIMER_VERBOSE_FLOPS("read_field_64");
  const long total_bytes = read_field(f, path, new_size_node);
  if (total_bytes == 0) {
    return 0;
  } else {
    to_from_big_endian_64(get_data(f));
    timer.flops += total_bytes;
    return total_bytes;
  }
}

template <class M>
long write_field_double(const Field<M>& f, const std::string& path,
                        const Coordinate& new_size_node = Coordinate())
// interface_function
{
  TIMER_VERBOSE_FLOPS("write_field_double");
  Field<M> ff;
  ff.init(f);
  to_from_big_endian_64(get_data(ff));
  const long total_bytes = write_field(ff, path, new_size_node);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
long read_field_double(Field<M>& f, const std::string& path,
                       const Coordinate& new_size_node = Coordinate())
// interface_function
{
  TIMER_VERBOSE_FLOPS("read_field_double");
  const long total_bytes = read_field(f, path, new_size_node);
  if (total_bytes == 0) {
    return 0;
  } else {
    to_from_big_endian_64(get_data(f));
    timer.flops += total_bytes;
    return total_bytes;
  }
}

inline bool is_field(const std::string& path)
{
  TIMER("is_field");
  long nfile = 0;
  if (get_id_node() == 0) {
    QFile qfile = qfopen(path, "r");
    if (not qfile.null()) {
      const std::string header = "BEGIN_FIELD_HEADER\n";
      std::vector<char> check_line(header.size(), 0);
      if (1 == qfread(check_line.data(), header.size(), 1, qfile)) {
        if (std::string(check_line.data(), check_line.size()) == header) {
          nfile = 1;
        }
      }
    }
    qfile.close();
  }
  bcast(get_data(nfile));
  return nfile > 0;
}

inline bool is_d_field(const std::string& path)
{
  return is_dist_field(path) or is_field(path);
}

inline bool dist_repartition(const Coordinate& new_size_node,
                             const std::string& path,
                             const std::string& new_path = "")
// interface_function
{
  bool is_failed = false;
  const std::string npath = remove_trailing_slashes(path);
  if (std::string(npath, npath.length() - 4, 4) == ".tmp") {
    return true;
  }
  if (does_file_exist_sync_node(npath + "-lock")) {
    return true;
  }
  const std::string new_npath = remove_trailing_slashes(new_path);
  const bool is_dir = is_directory_sync_node(path);
  if (is_dir and new_size_node != Coordinate(1, 1, 1, 1) and
      (new_npath == npath or new_npath == "")) {
    if (is_dist_field(npath)) {
      Geometry geo;
      int sizeof_M;
      Coordinate size_node;
      dist_read_geo_info(geo, sizeof_M, size_node, npath);
      if (size_node == new_size_node and
          (new_path == "" or new_npath == npath)) {
        displayln_info(
            ssprintf("repartition: size_node=%s ; no need to repartition '%s'.",
                     show(size_node).c_str(), npath.c_str()));
        return true;
      }
    } else {
      displayln_info(
          ssprintf("repartition: WARNING: not a folder to partition: '%s'.",
                   npath.c_str()));
    }
  }
  if (not is_dir and new_size_node == Coordinate(1, 1, 1, 1) and
      (new_npath == npath or new_npath == "")) {
    displayln_info(
        ssprintf("repartition: size_node=%s ; no need to repartition '%s'.",
                 show(new_size_node).c_str(), npath.c_str()));
    return true;
  }
  if (not obtain_lock(npath + "-lock")) {
    return true;
  }
  TIMER_VERBOSE("dist_repartition");
  Field<float> f;
  read_field(f, npath);
  if (get_incorrect_field_read_sizeof_M() != 0) {
    get_force_field_write_sizeof_M() = get_incorrect_field_read_sizeof_M();
    get_incorrect_field_read_sizeof_M() = 0;
  }
  if (new_npath == npath or new_npath == "") {
    qassert(not does_file_exist_sync_node(npath + "-repartition-new.tmp"));
    qassert(not does_file_exist_sync_node(npath + "-repartition-old.tmp"));
    if (new_size_node == Coordinate(1, 1, 1, 1)) {
      write_field(f, npath + "-repartition-new.tmp");
    } else {
      dist_write_field(f, new_size_node, npath + "-repartition-new.tmp");
    }
    qrename_info(npath, npath + "-repartition-old.tmp");
    qrename_info(npath + "-repartition-new.tmp", npath);
    qremove_all_info(npath + "-repartition-old.tmp");
  } else {
    if (new_size_node == Coordinate(1, 1, 1, 1)) {
      write_field(f, new_npath);
    } else {
      dist_write_field(f, new_size_node, new_npath);
    }
  }
  release_lock();
  return is_failed;
}

}  // namespace qlat
