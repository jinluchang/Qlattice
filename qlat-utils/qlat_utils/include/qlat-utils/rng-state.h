#pragma once

#include <string>
#include <vector>
#include <cstdint>
#include <qlat-utils/show.h>

namespace qlat
{  //

struct API RngState;

void reset(RngState& rs);

void reset(RngState& rs, const std::string& seed);

void reset(RngState& rs, const Long seed);

void split_rng_state(RngState& rs, const RngState& rs0,
                     const std::string& sindex);

void split_rng_state(RngState& rs, const RngState& rs0, const Long sindex = 0);

void set_type(RngState& rs, const uint64_t type);

uint64_t rand_gen(RngState& rs);

RealD u_rand_gen(RngState& rs, const RealD upper = 1.0,
                  const RealD lower = 0.0);

RealD g_rand_gen(RngState& rs, const RealD center = 0.0,
                  const RealD sigma = 1.0);

void compute_hash_with_input(uint32_t hash[8], const RngState& rs,
                             const std::string& input);

struct API RngState {
  uint64_t numBytes;
  uint32_t hash[8];
  uint64_t type;
  uint64_t index;
  //
  uint64_t cache[3];
  RealD gaussian;
  Int cacheAvail;
  bool gaussianAvail;
  //
  inline void init() { reset(*this); }
  //
  RngState() { init(); }
  RngState(const std::string& seed) { reset(*this, seed); }
  RngState(const Long seed) { reset(*this, seed); }
  RngState(const RngState& rs0, const std::string& sindex)
  {
    memset((void*)this, 0, sizeof(RngState));
    split_rng_state(*this, rs0, sindex);
  }
  RngState(const RngState& rs0, const Long sindex)
  {
    memset((void*)this, 0, sizeof(RngState));
    split_rng_state(*this, rs0, sindex);
  }
  //
  RngState split(const std::string& sindex) const
  {
    return RngState(*this, sindex);
  }
  RngState split(const Long sindex) const { return RngState(*this, sindex); }
  //
  RngState newtype(const uint64_t type) const
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

template <class M>
void random_permute(std::vector<M>& vec, const RngState& rs_)
{
  RngState rs = rs_;
  const Long size = (Long)vec.size();
  M tmp;
  for (Long k = 0; k < size; ++k) {
    const Long kk = rand_gen(rs) % (size - k);
    tmp = vec[k];
    vec[k] = vec[k + kk];
    vec[k + kk] = tmp;
  }
}

}  // namespace qlat
