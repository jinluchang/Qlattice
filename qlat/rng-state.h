#pragma once

#include <qlat/config.h>

#include <rng-state.h>

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

inline double g_rand_gen(RngState& rs, const double center = 0.0, const double sigma = 1.0)
{
  return gRandGen(rs, center, sigma);
}

inline RngState& get_global_rng_state()
{
  return getGlobalRngState();
}

QLAT_END_NAMESPACE
