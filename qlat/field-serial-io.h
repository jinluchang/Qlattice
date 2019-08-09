#pragma once

#include <qlat/field-dist-io.h>
#include <qlat/field-io.h>

#include <stdio.h>
#include <ctime>

#include <fstream>
#include <iostream>

namespace qlat
{  //

inline Coordinate get_default_serial_new_size_node(const Geometry& geo)
{
  const int num_node = geo.geon.num_node;
  const Coordinate total_site = geo.total_site();
  Coordinate new_size_node = Coordinate(1, 1, 1, total_site[3]);
  while (num_node < new_size_node[3]) {
    if (new_size_node[3] % 5 == 0) {
      new_size_node[3] /= 5;
    } else if (new_size_node[3] % 3 == 0) {
      new_size_node[3] /= 3;
    } else if (new_size_node[3] % 2 == 0) {
      new_size_node[3] /= 2;
    } else {
      break;
    }
  }
  qassert(total_site % new_size_node == Coordinate());
  return new_size_node;
}

template <class M>
long serial_write_field(const Field<M>& f, const std::string& path,
                        const Coordinate& new_size_node)
// will append to the file
// assume new_size_node is properly choosen so that concatenate the new fields
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
    FILE* fp = qopen(path, "a");
    qassert(fp != NULL);
    for (int new_id_node = 0; new_id_node < new_num_node; ++new_id_node) {
      const int id_node =
          get_id_node_from_new_id_node(new_id_node, new_num_node, num_node);
      if (0 == id_node) {
        assign(v, get_data(fs[new_id_node]));
      } else {
        MPI_Recv(v.data(), v.data_size(), MPI_BYTE, id_node, mpi_tag,
                 get_comm(), MPI_STATUS_IGNORE);
      }
      qwrite_data(v, fp);
    }
    qclose(fp);
  } else {
    for (size_t i = 0; i < fs.size(); ++i) {
      const Vector<M> v = get_data(fs[i]);
      MPI_Send((void*)v.data(), v.data_size(), MPI_BYTE, 0, mpi_tag,
               get_comm());
    }
  }
  const long file_size = get_data(f).data_size() * f.geo.geon.num_node;
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
  if (not does_file_exist_sync_node(path)) {
    displayln_info(fname +
                   ssprintf(": file does not exist: '%s'", path.c_str()));
    return 0;
  }
  const Geometry& geo = f.geo;
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
    f.init(fs[0].geo);
    Vector<M> v = get_data(f);
    const int num_node = get_num_node();
    const int new_num_node = product(new_size_node);
    FILE* fp = qopen(path, "r");
    qassert(fp != NULL);
    fseek(fp, offset, whence);
    for (int new_id_node = 0; new_id_node < new_num_node; ++new_id_node) {
      const int id_node =
          get_id_node_from_new_id_node(new_id_node, new_num_node, num_node);
      qread_data(v, fp);
      if (0 == id_node) {
        assign(get_data(fs[new_id_node]), v);
      } else {
        MPI_Send((void*)v.data(), v.data_size(), MPI_BYTE, id_node, mpi_tag,
                 get_comm());
      }
    }
    qclose(fp);
  } else {
    for (size_t i = 0; i < fs.size(); ++i) {
      Vector<M> v = get_data(fs[i]);
      MPI_Recv(v.data(), v.data_size(), MPI_BYTE, 0, mpi_tag, get_comm(),
               MPI_STATUS_IGNORE);
    }
  }
  shuffle_field_back(f, fs, new_size_node);
  const long file_size = get_data(f).data_size() * f.geo.geon.num_node;
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
  if (not does_file_exist_sync_node(path)) {
    displayln_info(fname +
                   ssprintf(": file does not exist: '%s'", path.c_str()));
    return 0;
  }
  const Geometry& geo = f.geo;
  std::vector<Field<M> > fs;
  const std::vector<Geometry> new_geos =
      make_dist_io_geos(geo.total_site(), geo.multiplicity, new_size_node);
  fs.resize(new_geos.size());
  for (size_t i = 0; i < fs.size(); ++i) {
    fs[i].init(new_geos[i]);
  }
  if (fs.size() > 0) {
    FILE* fp = qopen(path, "r");
    qassert(fp != NULL);
    fseek(fp, offset + fs[0].geo.geon.id_node * get_data(fs[0]).data_size(),
          whence);
    for (size_t i = 0; i < fs.size(); ++i) {
      Vector<M> v = get_data(fs[i]);
      qread_data(v, fp);
    }
    qclose(fp);
  }
  shuffle_field_back(f, fs, new_size_node);
  sync_node();
  const long file_size = get_data(f).data_size() * f.geo.geon.num_node;
  timer.flops += file_size;
  return file_size;
}

template <class M>
long serial_write_field(const Field<M>& f, const std::string& path)
// interface_function
{
  return serial_write_field(f, path, get_default_serial_new_size_node(f.geo));
}

template <class M>
long serial_read_field(Field<M>& f, const std::string& path,
                       const long offset = 0, const int whence = SEEK_SET)
// interface_function
{
  return serial_read_field(f, path, get_default_serial_new_size_node(f.geo),
                           offset, whence);
}

template <class M>
long serial_read_field_par(Field<M>& f, const std::string& path,
                           const long offset = 0, const int whence = SEEK_SET)
// interface_function
{
  return serial_read_field_par(f, path, get_default_serial_new_size_node(f.geo),
                               offset, whence);
}

template <class M>
crc32_t field_simple_checksum(const Field<M>& f)
{
  TIMER("field_simple_checksum");
  qassert(f.geo.is_only_local());
  crc32_t ret = 0;
  const Vector<M> v = get_data(f);
  Vector<crc32_t> vc((crc32_t*)v.data(), v.data_size() / sizeof(crc32_t));
  for (long i = 0; i < vc.size(); ++i) {
    ret += vc[i];
  }
  glb_sum_byte(ret);
  return ret;
}

template <class M>
crc32_t field_crc32_slow(const Field<M>& f)
{
  TIMER_VERBOSE_FLOPS("field_crc32_slow");
  const Geometry& geo = f.geo;
  const Coordinate new_size_node = get_default_serial_new_size_node(geo);
  const int new_num_node = product(geo.total_site() / new_size_node);
  const int num_node = get_num_node();
  const int id_node = get_id_node();
  std::vector<Field<M> > fs;
  shuffle_field(fs, f, new_size_node);
  crc32_t ret = 0;
  for (int i = 0; i < std::min(num_node, new_num_node); ++i) {
    if (i == id_node) {
      for (int k = 0; k < (int)fs.size(); ++k) {
        ret = crc32_par(ret, get_data(fs[k]));
      }
    } else {
      ret = 0;
    }
    glb_sum_byte(ret);
  }
  timer.flops += get_data(f).data_size() * geo.geon.num_node;
  return ret;
}

template <class M>
crc32_t field_crc32(const Field<M>& f)
{
  TIMER_VERBOSE_FLOPS("field_crc32");
  const Geometry& geo = f.geo;
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

// template<class M>
// void field_import_serial(qlat::Field<M>& destination,
//     const std::string& read_addr,
//     const long offset = 0,
//     const int whence = SEEK_SET,
//     const int num_of_reading_threads = 0)
// {
//   TIMER_VERBOSE("field_import_serial");
//
//   Geometry geo_only_local = geo_resize(destination.geo, 0);;
//
//   Field<M> field_recv;
//   field_recv.init(geo_only_local);
//
//   Field<M> field_send;
//   field_send.init(geo_only_local);
//
//   Field<M> field_rslt;
//   field_rslt.init(geo_only_local);
//
//   Coordinate total_site = geo_only_local.total_site();
//
//   long range_low = geo_only_local.local_volume() * get_id_node();
//   long range_high = range_low + geo_only_local.local_volume();
//
//   // for every node:
//   //
//   FILE *inputFile = NULL;
//
//   // Well as you can see this is not really serial reading anymore. The
//   sertial reading speed is unbearablly slow.
//   // Anyway it is tested. And it seems to be right.
//   inputFile = fopen(read_addr.c_str(), "r");
//   qassert(inputFile != NULL);
//   qassert(!ferror(inputFile));
//   fseek(inputFile, offset, whence);
//   qassert(!feof(inputFile));
//
//   sync_node();
//   int cycle_limit = 0;
//   if(num_of_reading_threads > 0)
//     cycle_limit = (int)ceil((double)get_num_node() / num_of_reading_threads);
//   else
//     cycle_limit = 1;
//   for(int cycle = 0; cycle < cycle_limit; cycle++){
//     if(get_id_node() % cycle_limit == cycle){
//       std::cout << "Reading STARTED:  Node Number =\t"
//         << get_id_node() << std::endl;
//       M *ptr = get_data(field_send).data();
//       long size = sizeof(M) * geo_only_local.local_volume()
//         * geo_only_local.multiplicity;
//       assert(!fseek(inputFile, size * get_id_node(), SEEK_CUR));
//       timer_fread((char*)ptr, size, inputFile);
//       std::cout << "Reading FINISHED: Node Number =\t"
//         << get_id_node() << std::endl;
//       fclose(inputFile);
//     }
//     sync_node();
//   }
//
//   // 	if(get_id_node() == 0){
//   //         	// input.open(read_addr.c_str());
//   // 		inputFile = fopen(read_addr.c_str(), "r");
//   // 		char line[1000];
//   // 		char indicator[] = "END_HEADER";
//   //
//   // 		int pos_ = -1; fpos_t pos;
//   // 		rewind(inputFile);
//   // 		while(fgets(line, 1000, inputFile) != NULL)
//   // 		{if(strstr(line, indicator) != NULL){
//   // 			fgetpos(inputFile, &pos); pos_ = 1;  break;
//   // 		}}
//   // 		qassert(pos_ > -1); qassert(!feof(inputFile));
//   // 	}
//   //
//   // 	for(int i = 0; i < get_num_node(); i++){
//   //
//   // 		if(get_id_node() == 0){
//   // 			std::cout << "Reading Cycle: " << i <<
//   std::endl;
//   // 			M *ptr = get_data(field_send).data();
//   //                         long size = sizeof(M) *
//   geo_only_local.local_volume() * geo_only_local.multiplicity;
//   //                         timer_fread((char*)ptr, size, inputFile);
//   // 		// 	fflush(inputFile);
//   // 		}
//   // 		sync_node();
//   // 		get_data_dir(get_data(field_recv), get_data(field_send),
//   0);
//   //                 qswap(field_recv, field_send);
//   // 	}
//   //
//   // 	if(get_id_node() == 0) fclose(inputFile);
//
//   field_rslt = field_send;
//
//   for(int i = 0; i < get_num_node(); i++){
//
//     int id_send_node = (get_id_node() + i) % get_num_node();
//
//     Coordinate coor_send_node = qlat::coordinate_from_index(id_send_node,
//     geo_only_local.geon.size_node);
// #pragma omp parallel for
//     for(int index = 0; index < geo_only_local.local_volume(); index++){
//       Coordinate local_coor = geo_only_local.coordinate_from_index(index);
//       Coordinate global_coor;
//       for (int mu = 0; mu < 4; mu++) {
//         global_coor[mu] = local_coor[mu] + coor_send_node[mu] *
//         geo_only_local.node_site[mu];
//       }
//       long global_index = index_from_coordinate(global_coor, total_site);
//       if(global_index >= range_low && global_index < range_high)
//       {
//         Coordinate local_coor_read =
//         geo_only_local.coordinate_from_index(global_index - range_low);
//         assign(field_send.get_elems(local_coor),
//         field_rslt.get_elems_const(local_coor_read));
//       }
//     }
//
//     get_data_dir(get_data(field_recv), get_data(field_send), 0);
//     qswap(field_recv, field_send);
//     if(get_id_node() == 0) std::cout << "Shuffling CYCLE:\t" << i <<
//     std::endl;
//   }
//
//   destination = field_send;
//
//   sync_node();
// }

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
                 const Coordinate& new_size_node_ = Coordinate())
// assume new_size_node is properly choosen so that concatenate the new fields
// would be correct. eg. new_size_node = Coordinate(1,1,1,2)
{
  TIMER_VERBOSE_FLOPS("write_field");
  qassert(is_initialized(f));
  const Geometry& geo = f.geo;
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
  const Coordinate new_size_node = new_size_node_ == Coordinate()
                                       ? get_default_serial_new_size_node(geo)
                                       : new_size_node_;
  const long file_size =
      serial_write_field(f, path + ".partial", new_size_node);
  qrename_info(path + ".partial", path);
  timer.flops += file_size;
  return file_size;
}

inline void read_geo_info(Coordinate& total_site, int& multiplicity, int& sizeof_M,
                          crc32_t& crc, const std::string& path)
{
  TIMER("read_geo_info");
  if (get_id_node() == 0) {
    FILE* fp = qopen(path, "r");
    if (fp != NULL) {
      const std::string header = "BEGIN_FIELD_HEADER\n";
      std::vector<char> check_line(header.size(), 0);
      if (1 == fread(check_line.data(), header.size(), 1, fp)) {
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
          reads(sizeof_M, info_get_prop(infos, "sizeof(M) = "));
          crc = read_crc32(info_get_prop(infos, "field_crc32 = "));
        }
      }
    }
    qclose(fp);
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
  if (does_file_exist_sync_node(path + "/geo-info.txt")) {
    return dist_read_field(f, path);
  }
  TIMER_VERBOSE_FLOPS("read_field");
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
    displayln(fname + ssprintf(": WARNING: sizeof(M) do not match. "
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
  const Coordinate new_size_node = new_size_node_ == Coordinate()
                                       ? get_default_serial_new_size_node(geo)
                                       : new_size_node_;
  const long file_size =
      serial_read_field_par(f, path, new_size_node, -data_size, SEEK_END);
  qassert(crc == field_crc32(f));
  qassert(file_size == data_size);
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
    FILE* fp = qopen(path, "r");
    if (fp != NULL) {
      const std::string header = "BEGIN_FIELD_HEADER\n";
      std::vector<char> check_line(header.size(), 0);
      if (1 == fread(check_line.data(), header.size(), 1, fp)) {
        if (std::string(check_line.data(), check_line.size()) == header) {
          nfile = 1;
        }
      }
    }
    qclose(fp);
  }
  bcast(get_data(nfile));
  return nfile > 0;
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
  return is_failed;
}

}  // namespace qlat
