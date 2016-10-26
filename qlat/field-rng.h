#pragma once

#include <qlat/config.h>
#include <qlat/field.h>

#include <rng-state.h>

#include <cassert>

QLAT_START_NAMESPACE

template <class T>
void split_rng_state(RngState& rs, const RngState& rs0, const T& s)
{
  splitRngState(rs, rs0, s);
}

inline uint64_t rand_gen(RngState& rs)
{
  return randGen(rs);
}

inline double u_rand_gen(RngState& rs, const double upper = 1.0, const double lower = 0.0)
{
  return uRandGen(rs, upper, lower);
}

inline double g_rand_gen(RngState& rs, const double sigma = 1.0, const double center = 0.0)
{
  return gRandGen(rs, sigma, center);
}

struct RngField : FieldM<RngState,1>
{
  virtual const char* cname()
  {
    return "RngField";
  }
  //
  virtual void init()
  {
    FieldM<RngState,1>::init();
  }
  virtual void init(const Geometry& geo_, const RngState& rs)
  {
    FieldM<RngState,1>::init(geo_);
    Coordinate total_site;
    for (int mu = 0; mu < DIM; ++mu) {
      total_site[mu] = geo.total_site(mu);
    }
#pragma omp parallel for
    for (long index = 0; index < geo.local_volume(); ++index) {
      Coordinate x; geo.coordinate_from_index(x, index);
      Coordinate xg; geo.coordinate_g_from_l(xg, x);
      long gindex = index_from_coordinate(xg, total_site);
      split_rng_state(get_elem(x), rs, gindex);
    }
  }
  //
  RngField()
  {
    FieldM<RngState,1>::init();
  }
  RngField(const RngField& rf)
  {
    assert(false);
  }
};

QLAT_END_NAMESPACE
