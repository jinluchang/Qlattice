#pragma once

#include <qlat/field-dist-io.h>
#include <qlat/field-io.h>

#include <stdio.h>
#include <ctime>

#include <fstream>
#include <iostream>

QLAT_START_NAMESPACE

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
crc32_t field_crc32(const Field<M>& f)
{
  TIMER_VERBOSE_FLOPS("field_crc32");
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
      for (int k = 0; k < fs.size(); ++k) {
        ret = crc32_par(ret, get_data(fs[k]));
      }
    }
    glb_sum_byte(ret);
  }
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

QLAT_END_NAMESPACE
