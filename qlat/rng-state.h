#pragma once

#include <qlat/config.h>

#include <hash-cpp/sha256.h>
#include <timer.h>

#include <cstring>
#include <cmath>

QLAT_START_NAMESPACE

const int RngStateSize = 8 * 8;

struct RngState
{
  uint8_t state[RngStateSize];
  //
  uint64_t cache[4];
  double gaussion;
  int cacheAvail;
  bool gaussionAvail;
  //
  RngState()
  {
    init();
  }
  //
  RngState(const uint64_t seed, const uint64_t shift = 0)
  {
    init(seed, shift);
  }
  //
  void init()
  {
    std::memset(state, 0, sizeof(state));
    cacheAvail = 0;
    gaussionAvail = false;
  }
  //
  void init(const uint64_t seed, const uint64_t shift = 0)
  {
    init(seed, 0, 0, 0, shift);
  }
  //
  void init(const uint64_t seed, const uint64_t type, const uint64_t traj, const uint64_t index = 0, const uint64_t shift = 0)
  {
    init();
    uint64_t ss[RngStateSize/8];
    std::memset(ss, 0, sizeof(ss));
    ss[RngStateSize/8-1] = seed;
    ss[RngStateSize/8-2] = type;
    ss[RngStateSize/8-3] = traj;
    ss[RngStateSize/8-4] = index;
    ss[0] = shift;
    // std::memcpy(state, ss, sizeof(ss));
    for (int i = 0; i < RngStateSize/8; ++i) {
      for (int j = i * 8; j < i * 8 + 8; ++j) {
        state[j] = ss[i] % 256;
        ss[i] /= 256;
      }
    }
  }
};

inline void shiftRngState(RngState& rs, const uint64_t shift = 1)
{
  uint64_t pre = shift;
  int i = 0;
  do {
    pre += rs.state[i];
    rs.state[i] = pre % 256;
    pre /= 256;
    i += 1;
  } while (pre > 0 && i < RngStateSize);
}

inline uint64_t randGen(RngState& rs)
{
  if (rs.cacheAvail > 0) {
    rs.cacheAvail -= 1;
    return rs.cache[rs.cacheAvail];
  } else {
    SHA256 sha;
    sha.add(rs.state, RngStateSize * sizeof(uint8_t));
    uint8_t hash[4 * 8];
    sha.getHash(hash);
    for (int i = 0; i < 4; ++i) {
      rs.cache[i] = 0;
      uint64_t unit = 1;
      for (int j = i * 8; j < i * 8 + 8; ++j) {
        rs.cache[i] += hash[j] * unit;
        unit *= 256;
      }
    }
    // sha.getHash((unsigned char*)rs.cache);
    shiftRngState(rs);
    rs.cacheAvail = 3;
    return rs.cache[3];
  }
}

inline double uRandGen(RngState& rs, const double upper = 1.0, const double lower = 0.0)
{
  uint64_t u = randGen(rs);
  const double fac = 1.0 / (256.0 * 256.0 * 256.0 * 256.0) / (256.0 * 256.0 * 256.0 * 256.0);
  return u * fac * (upper - lower) + lower;
}

inline double gRandGen(RngState& rs, const double sigma = 1.0, const double center = 0.0)
{
  if (rs.gaussionAvail) {
    rs.gaussionAvail = false;
    return rs.gaussion * sigma + center;
  } else {
    const char* cname = "";
    const char* fname = "gRandGen()";
    int num_try = 1;
    double v1, v2, rsq;
    do {
      v1 = uRandGen(rs, 1.0, -1.0);
      v2 = uRandGen(rs, 1.0, -1.0);
      if ((num_try % 1000)==0) {
        Display(cname,fname,"num_try=%d v1=%e v2=%e\n",num_try,v1,v2);
      }
      rsq = v1*v1 + v2*v2;
      num_try++;
    } while ((num_try < 10000) && (rsq >= 1.0 || rsq == 0));
    if (num_try > 9999) {
      Display(cname, fname, "failed after 10000 tries (corrupted RNG?), returning ridiculous numbers (1e+10)\n");
      return 1e+10;
    }
    // pick 2 uniform numbers in the square extending from
    // -1 to 1 in each direction, see if they are in the
    // unit circle, and try again if they are not.
    //
    double fac = std::sqrt(-2.0 * std::log(rsq)/rsq);
    rs.gaussion = v1 * fac;
    rs.gaussionAvail = true;
    return v2 * fac * sigma + center;
  }
}

QLAT_END_NAMESPACE
