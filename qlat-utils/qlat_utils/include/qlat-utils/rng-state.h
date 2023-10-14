#pragma once

#include <string>
#include <cstdint>

namespace qlat
{  //

struct API RngState;

void reset(RngState& rs);

void reset(RngState& rs, const std::string& seed);

void reset(RngState& rs, const long seed);

void split_rng_state(RngState& rs, const RngState& rs0,
                     const std::string& sindex);

void split_rng_state(RngState& rs, const RngState& rs0, const long sindex = 0);

void set_type(RngState& rs, const unsigned long type);

uint64_t rand_gen(RngState& rs);

double u_rand_gen(RngState& rs, const double upper = 1.0,
                  const double lower = 0.0);

double g_rand_gen(RngState& rs, const double center = 0.0,
                  const double sigma = 1.0);

void compute_hash_with_input(uint32_t hash[8], const RngState& rs,
                             const std::string& input);

struct API RngState {
  uint64_t numBytes;
  uint32_t hash[8];
  unsigned long type;
  unsigned long index;
  //
  uint64_t cache[3];
  double gaussian;
  int cacheAvail;
  bool gaussianAvail;
  //
  inline void init() { reset(*this); }
  //
  RngState() { init(); }
  RngState(const std::string& seed) { reset(*this, seed); }
  RngState(const long seed) { reset(*this, seed); }
  RngState(const RngState& rs0, const std::string& sindex)
  {
    std::memset((void*)this, 0, sizeof(RngState));
    split_rng_state(*this, rs0, sindex);
  }
  RngState(const RngState& rs0, const long sindex)
  {
    std::memset((void*)this, 0, sizeof(RngState));
    split_rng_state(*this, rs0, sindex);
  }
  //
  RngState split(const std::string& sindex) const
  {
    return RngState(*this, sindex);
  }
  RngState split(const long sindex) const { return RngState(*this, sindex); }
  //
  RngState newtype(const unsigned long type) const
  {
    RngState rs(*this);
    set_type(rs, type);
    return rs;
  }
};

API inline RngState& get_global_rng_state()
{
  static RngState rs;
  return rs;
}

}  // namespace qlat
