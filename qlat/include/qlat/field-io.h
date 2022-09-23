#pragma once

#include <qlat/config.h>
#include <qlat/geometry.h>
#include <qlat/mpi.h>


#include <stdio.h>
#include <ctime>

#include <fstream>
#include <iostream>

namespace qlat
{  //

typedef array<Complex, 6> MatrixTruncatedSU3;
typedef array<Complex, 9> MatrixSU3;

class rePort
{
 public:
  std::ostream *os;
  rePort() { os = &std::cout; }
};

template <class T>
const rePort &operator<<(const rePort &p, const T &data)
{
  if (get_id_node() == 0) *(p.os) << data;
  return p;
}

inline const rePort &operator<<(const rePort &p,
                                std::ostream &(*func)(std::ostream &))
{
  if (get_id_node() == 0) *(p.os) << func;
  return p;
}

// static const rePort report;

inline std::string str_printf(const char *format, ...)
{
  char cstr[512];
  va_list args;
  va_start(args, format);
  vsnprintf(cstr, sizeof(cstr), format, args);
  return std::string(cstr);
}

inline int Printf(const char *format, ...)
{
  if (!get_id_node()) {
    va_list args;
    va_start(args, format);
    return vprintf(format, args);
  } else {
    return 0;
  }
}

inline FILE *Fopen(const char *filename, const char *mode)
{
  if (!get_id_node()) {
    return fopen(filename, mode);
  } else {
    return NULL;
  }
}

inline int Fprintf(FILE *pFile, const char *format, ...)
{
  if (!get_id_node()) {
    va_list args;
    va_start(args, format);
    return vfprintf(pFile, format, args);
  } else {
    return 0;
  }
}

inline int Fflush(FILE *pFile)
{
  if (!get_id_node()) {
    return fflush(pFile);
  } else {
    return 0;
  }
}

inline bool Is_not_null(FILE *pFile)
{
  if (!get_id_node()) {
    return pFile != NULL;
  } else {
    return false;
  }
}

template <class M, class N>
void castTruncated(M &x, const N &y)
{
  qassert(sizeof(M) <= sizeof(N));
  memcpy(&x, &y, sizeof(M));
}

template <class M, class N>
void fieldCastTruncated(Field<M> &dest, const Field<N> &src)
{
  TIMER("fieldCastTruncated");
  const Geometry &geo = src.geo();
  dest.init(geo);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    Coordinate xl = geo.coordinate_from_index(index);
    const Vector<N> s = src.get_elems_const(xl);
    Vector<M> d = dest.get_elems(xl);
    for (int m = 0; m < geo.multiplicity; ++m) {
      castTruncated(d[m], s[m]);
    }
  }
}

template <class M>
uint32_t fieldChecksumSum32(const Field<M> &f)
{
  TIMER("fieldChecksumSum32");

  Geometry geo = geo_resize(f.geo(), 0);
  qassert(geo.is_only_local());

  Field<M> f_local;
  f_local.init(geo);
  f_local = f;

  qassert(sizeof(M) % sizeof(uint32_t) == 0);
  long sum = 0;
  const uint32_t *data = (const uint32_t *)f_local.field.data();
  const long size = f_local.field.size() * sizeof(M) / sizeof(uint32_t);
  for (long i = 0; i < size; ++i) {
    sum += data[i];
  }
#ifdef USE_MULTI_NODE
  glb_sum(sum);
#endif
  uint32_t cs = sum;
  displayln_info("qlat::" + fname + "(): " + ssprintf("check sum = %x", cs));
  return cs;
}

template <class M>
std::string field_hash_crc32(const qlat::Field<M> &origin)
{
  // somehow this checksum function does not agree with CPS's one.
  // Do not know why. But I am not sure what algorithm CPS uses.

  TIMER("field_hash_crc32");

  Geometry geo_only_local = geo_resize(origin.geo(), 0);
  crc32_t hash;
  for (int id_node = 0; id_node < get_num_node(); id_node++) {
    if (get_id_node() == id_node) {
      crc32(hash, (void *)get_data(origin).data(),
            get_data(origin).size() * sizeof(M));
    }
    sync_node();
    MPI_Bcast((void *)&hash, 4, MPI_BYTE, id_node, get_comm());
  }
  return ssprintf("%08X", hash);
}

inline void timer_fwrite(char *ptr, long size, FILE *outputFile)
{
  TIMER("timer_fwrite");
  fwrite(ptr, size, 1, outputFile);
}

template <class M>
void sophisticated_make_to_order(Field<M> &result, const Field<M> &origin)
{
  TIMER("sophisticated_make_to_order");

  Geometry geo_only_local = geo_resize(origin.geo(), 0);
  ;

  Field<M> field_recv;
  field_recv.init(geo_only_local);

  Field<M> field_send;
  field_send.init(geo_only_local);
  field_send = origin;

  Field<M> field_rslt;
  field_rslt.init(geo_only_local);

  Coordinate total_site = geo_only_local.total_site();

  long range_low = geo_only_local.local_volume() * get_id_node();
  long range_high = range_low + geo_only_local.local_volume();

  for (int i = 0; i < get_num_node(); i++) {
    // std::cout << "Shuffle loop: " << i << std::endl;

    int id_send_node = (get_id_node() + i) % get_num_node();

    Coordinate coor_send_node = qlat::coordinate_from_index(
        id_send_node, geo_only_local.geon.size_node);
#pragma omp parallel for
    for (int index = 0; index < geo_only_local.local_volume(); index++) {
      Coordinate local_coor = geo_only_local.coordinate_from_index(index);
      Coordinate global_coor;
      for (int mu = 0; mu < 4; mu++) {
        global_coor[mu] =
            local_coor[mu] + coor_send_node[mu] * geo_only_local.node_site[mu];
      }
      long global_index = index_from_coordinate(global_coor, total_site);
      if (global_index >= range_low && global_index < range_high) {
        Coordinate local_coor_write =
            geo_only_local.coordinate_from_index(global_index - range_low);
        assign(field_rslt.get_elems(local_coor_write),
               field_send.get_elems_const(local_coor));
      }
    }

    get_data_dir(get_data(field_recv), get_data(field_send), 0);
    qswap(field_recv, field_send);
  }
  result.init(geo_only_local);
  result = field_rslt;
}

template <class M>
void sophisticated_serial_write(const qlat::Field<M> &origin,
                                const std::string &write_addr,
                                const bool is_append = false)
{
  TIMER("sophisticated_serial_write");

  Geometry geo_only_local = geo_resize(origin.geo(), 0);

  Field<M> field_recv;
  field_recv.init(geo_only_local);

  Field<M> field_send;
  field_send.init(geo_only_local);
  field_send = origin;

  FILE *outputFile = NULL;
  if (get_id_node() == 0) {
    if (is_append)
      outputFile = fopen(write_addr.c_str(), "a");
    else
      outputFile = fopen(write_addr.c_str(), "w");
  }

  for (int i = 0; i < get_num_node(); i++) {
    if (get_id_node() == 0) {
      M *ptr = get_data(field_send).data();
      qassert(ptr != NULL);
      long size = sizeof(M) * geo_only_local.local_volume() *
                  geo_only_local.multiplicity;
      std::cout << "Writing CYCLE: " << i << "\tSIZE = " << size << std::endl;
      timer_fwrite((char *)ptr, size, outputFile);
      fflush(outputFile);
    }

    get_data_dir(get_data(field_recv), get_data(field_send), 0);
    qswap(field_recv, field_send);
  }

  if (get_id_node() == 0) fclose(outputFile);

  displayln("Export file CLOSED");

  sync_node();
}

// std::string cps_Matrix_header_generator(const qlat::Field<cps::Matrix>
// &origin,
//							const bool
//does_skip_third = false){ 	NOT yet implemented :( 	qassert(false);
//return "NOT
// IMPLEMENTED.";
// }

inline void timer_fread(char *ptr, long size, FILE *inputFile)
{
  TIMER_VERBOSE("timer_fread");
  fread(ptr, size, 1, inputFile);
}

template <class M>
void sophisticated_serial_read(qlat::Field<M> &destination,
                               const std::string &import, int pos,
                               const int num_of_reading_threads = 0)
{
  // Blindly read binary data into field. All checking should be handled by
  // wrapper functions.

  TIMER_VERBOSE("sophisticated_serial_read");

  Geometry geo_only_local = geo_resize(destination.geo(), 0);
  ;

  Field<M> field_recv;
  field_recv.init(geo_only_local);

  Field<M> field_send;
  field_send.init(geo_only_local);

  Field<M> field_rslt;
  field_rslt.init(geo_only_local);

  Coordinate total_site = geo_only_local.total_site();

  // for every node:
  //
  FILE *input = NULL;

  // Well as you can see this is not really serial reading anymore. The sertial
  // reading speed is unbearablly slow. Anyway it is tested. And it seems to be
  // right.
  input = fopen(import.c_str(), "rb");
  qassert(input != NULL);
  qassert(!ferror(input));

  fseek(input, pos, SEEK_SET);

#ifndef USE_SINGLE_NODE

  long range_low = geo_only_local.local_volume() * get_id_node();
  long range_high = range_low + geo_only_local.local_volume();

  sync_node();
  int cycle_limit = 0;
  if (num_of_reading_threads > 0)
    cycle_limit = (int)ceil((double)get_num_node() / num_of_reading_threads);
  else
    cycle_limit = 1;
  for (int cycle = 0; cycle < cycle_limit; cycle++) {
    if (get_id_node() % cycle_limit == cycle) {
      std::cout << "Reading STARTED:  Node Number =\t" << get_id_node()
                << std::endl;
      M *ptr = field_send.field.data();
      long size = sizeof(M) * geo_only_local.local_volume() *
                  geo_only_local.multiplicity;
      qassert(!fseek(input, size * get_id_node(), SEEK_CUR));
      timer_fread((char *)ptr, size, input);
      std::cout << "Reading FINISHED: Node Number =\t" << get_id_node()
                << std::endl;
      fclose(input);
    }
    sync_node();
  }

  //	if(get_id_node() == 0){
  //			// input.open(read_addr.c_str());
  //		inputFile = fopen(read_addr.c_str(), "r");
  //		char line[1000];
  //		char indicator[] = "END_HEADER";
  //
  //		int pos_ = -1; fpos_t pos;
  //		rewind(inputFile);
  //		while(fgets(line, 1000, inputFile) != NULL)
  //		{if(strstr(line, indicator) != NULL){
  //			fgetpos(inputFile, &pos); pos_ = 1;  break;
  //		}}
  //		qassert(pos_ > -1); qassert(!feof(inputFile));
  //	}
  //
  //	for(int i = 0; i < get_num_node(); i++){
  //
  //		if(get_id_node() == 0){
  //			std::cout << "Reading Cycle: " << i << std::endl;
  //			M *ptr = get_data(field_send).data();
  //                         long size = sizeof(M) *
  //                         geo_only_local.local_volume() *
  //                         geo_only_local.multiplicity;
  //                         timer_fread((char*)ptr, size, inputFile);
  //		//	fflush(inputFile);
  //		}
  //		sync_node();
  //		get_data_dir(get_data(field_recv), get_data(field_send), 0);
  //                 qswap(field_recv, field_send);
  //	}
  //
  //	if(get_id_node() == 0) fclose(inputFile);

  field_rslt = field_send;

  for (int i = 0; i < get_num_node(); i++) {
    int id_send_node = (get_id_node() + i) % get_num_node();

    Coordinate coor_send_node = qlat::coordinate_from_index(
        id_send_node, geo_only_local.geon.size_node);
#pragma omp parallel for
    for (int index = 0; index < geo_only_local.local_volume(); index++) {
      Coordinate local_coor = geo_only_local.coordinate_from_index(index);
      Coordinate global_coor;
      for (int mu = 0; mu < 4; mu++) {
        global_coor[mu] =
            local_coor[mu] + coor_send_node[mu] * geo_only_local.node_site[mu];
      }
      long global_index = index_from_coordinate(global_coor, total_site);
      if (global_index >= range_low && global_index < range_high) {
        Coordinate local_coor_read =
            geo_only_local.coordinate_from_index(global_index - range_low);
        assign(field_send.get_elems(local_coor),
               field_rslt.get_elems_const(local_coor_read));
      }
    }

    get_data_dir(get_data(field_recv), get_data(field_send), 0);
    qswap(field_recv, field_send);
    if (get_id_node() == 0) std::cout << "Shuffling CYCLE:\t" << i << std::endl;
  }

  destination = field_send;

  sync_node();

#else

  M *ptr = field_rslt.field.data();
  long size = geo_only_local.local_volume() * geo_only_local.multiplicity;
  //	printf("read = %d\n", fread((char*)ptr, sizeof(M), size, input));
  fread((char *)ptr, sizeof(M), size, input);

  destination = field_rslt;

#endif
}

}  // namespace qlat
