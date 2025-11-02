// utils_geo_cache.h
// Gen Wang
// Sep. 2025

#ifndef UTILS_GEO_CACHE_H
#define UTILS_GEO_CACHE_H
#pragma once

#include <string.h>
#include <sys/resource.h>
#include <mpi.h>
#include <time.h>
#include <typeinfo>
#include <iterator>
#include <cstdarg>
#include <qlat-utils/mat-vec.h>
#include <qlat-utils/eigen.h>
#include <qlat/qcd.h>

#include <qlat/vector_utils/utils_macro_defines.h>

#ifndef QLAT_NO_SYSINFO
#include <sys/sysinfo.h>
#endif

namespace qlat{

#define QLAT_MAX_GEO_CACHE 512

bool Compare_geo(const Geometry& g0, const Geometry& g1, const bool compare_bytes = true);
bool Compare_geo_less(const Geometry& g0, const Geometry& g1, const bool compare_bytes = true);

// compare elements
inline bool compare_geo(const Geometry& g0, const Geometry& g1){
  return Compare_geo(g0, g1, false);
}

void Get_geo_local(const qlat::Geometry& geo, Geometry& geo_l);
Geometry& get_geo_cache(const Coordinate& tot, box<Geometry>& geo_BOX);
Geometry& get_geo_resize(const Geometry& geo, const Coordinate& gl, const Coordinate& gr, box<Geometry>& geo_BOX);
Geometry& get_geo_local(const Geometry& geo, box<Geometry>& geo_BOX);

// read geo into vectors
void geo_to_nv(const Geometry& geo, std::vector<Int >& nv, std::vector<Int >& Nv, std::vector<Int>& mv);
void geo_to_nv(const Geometry& geo, vector<Int >& nv, vector<Int >& Nv, vector<Int >& mv);

// buffers for global geo box, will be replaced latter
struct Gbox_key {
  Geometry geo;
  Gbox_key(const Geometry& geo_)
  {
    const MemType gmem  = check_mem_type(&geo_);
    Qassert(gmem == MemType::Uvm or gmem == MemType::Cpu);
    copy_mem(&geo, MemType::Cpu, &geo_, gmem, sizeof(Geometry));
  }
};

bool operator<(const Gbox_key& x, const Gbox_key& y);

Geometry& get_geo_cache(const Coordinate& Lat);
Geometry& get_geo_cache(const Geometry& geo);
Geometry& get_geo_cache(const Geometry& geo, const Coordinate& gl, const Coordinate& gr);
Geometry& get_geo_local(const Geometry& geo);
void clear_geo_cache();

}

#endif
