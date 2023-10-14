// vim: set ts=2 sw=2 expandtab:

// Copyright (c) 2022 Luchang Jin
// All rights reserved.

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <qlat-utils/assert.h>
#include <qlat-utils/config.h>
#include <qlat-utils/rng-state.h>
#include <qlat-utils/sha256.h>
#include <qlat-utils/show.h>
#include <stdint.h>

#include <climits>
#include <cmath>
#include <cstring>
#include <istream>
#include <ostream>
#include <string>
#include <vector>

namespace qlat
{  //

void set_type(RngState& rs, const unsigned long type)
{
  qassert(ULONG_MAX == rs.type);
  qassert(ULONG_MAX != type);
  rs.type = type;
  rs.cacheAvail = 0;
  rs.gaussianAvail = false;
}

const size_t RNG_STATE_NUM_OF_INT32 = 2 + 8 + 2 + 2 + 3 * 2 + 2 + 1 + 1;

uint64_t patchTwoUint32(const uint32_t a, const uint32_t b)
{
  return (uint64_t)a << 32 | (uint64_t)b;
}

void splitTwoUint32(uint32_t& a, uint32_t& b, const uint64_t x)
{
  b = (uint32_t)x;
  a = (uint32_t)(x >> 32);
  qassert(x == patchTwoUint32(a, b));
}

void exportRngState(uint32_t* v, const RngState& rs)
{
  qassert(24 == RNG_STATE_NUM_OF_INT32);
  splitTwoUint32(v[0], v[1], rs.numBytes);
  for (int i = 0; i < 8; ++i) {
    v[2 + i] = rs.hash[i];
  }
  splitTwoUint32(v[10], v[11], rs.type);
  splitTwoUint32(v[12], v[13], rs.index);
  for (int i = 0; i < 3; ++i) {
    splitTwoUint32(v[14 + i * 2], v[14 + i * 2 + 1], rs.cache[i]);
  }
  union {
    double d;
    uint64_t l;
  } g;
  g.d = rs.gaussian;
  splitTwoUint32(v[20], v[21], g.l);
  v[22] = rs.cacheAvail;
  v[23] = rs.gaussianAvail;
}

void importRngState(RngState& rs, const uint32_t* v)
{
  qassert(24 == RNG_STATE_NUM_OF_INT32);
  rs.numBytes = patchTwoUint32(v[0], v[1]);
  for (int i = 0; i < 8; ++i) {
    rs.hash[i] = v[2 + i];
  }
  rs.type = patchTwoUint32(v[10], v[11]);
  rs.index = patchTwoUint32(v[12], v[13]);
  for (int i = 0; i < 3; ++i) {
    rs.cache[i] = patchTwoUint32(v[14 + i * 2], v[14 + i * 2 + 1]);
  }
  union {
    double d;
    uint64_t l;
  } g;
  g.l = patchTwoUint32(v[20], v[21]);
  rs.gaussian = g.d;
  rs.cacheAvail = v[22];
  rs.gaussianAvail = v[23];
}

void exportRngState(std::vector<uint32_t>& v, const RngState& rs)
{
  v.resize(RNG_STATE_NUM_OF_INT32);
  exportRngState(v.data(), rs);
}

void importRngState(RngState& rs, const std::vector<uint32_t>& v)
{
  qassert(RNG_STATE_NUM_OF_INT32 == v.size());
  importRngState(rs, v.data());
}

std::ostream& operator<<(std::ostream& os, const RngState& rs)
{
  std::vector<uint32_t> v(RNG_STATE_NUM_OF_INT32);
  exportRngState(v, rs);
  for (size_t i = 0; i < v.size() - 1; ++i) {
    os << v[i] << " ";
  }
  os << v.back();
  return os;
}

std::istream& operator>>(std::istream& is, RngState& rs)
{
  std::vector<uint32_t> v(RNG_STATE_NUM_OF_INT32);
  for (size_t i = 0; i < v.size(); ++i) {
    is >> v[i];
  }
  importRngState(rs, v);
  return is;
}

std::string show(const RngState& rs) { return shows(rs); }

bool operator==(const RngState& rs1, const RngState& rs2)
{
  return 0 == std::memcmp(&rs1, &rs2, sizeof(RngState));
}

bool operator!=(const RngState& rs1, const RngState& rs2)
{
  return !(rs1 == rs2);
}

void reset(RngState& rs)
{
  std::memset((void*)&rs, 0, sizeof(RngState));
  rs.numBytes = 0;
  sha256::setInitialHash(rs.hash);
  rs.type = ULONG_MAX;
  rs.index = 0;
  rs.cache[0] = 0;
  rs.cache[1] = 0;
  rs.cache[2] = 0;
  rs.gaussian = 0.0;
  rs.cacheAvail = 0;
  rs.gaussianAvail = false;
}

void reset(RngState& rs, const std::string& seed)
{
  reset(rs);
  split_rng_state(rs, rs, seed);
}

void reset(RngState& rs, const long seed) { reset(rs, show(seed)); }

void split_rng_state(RngState& rs, const RngState& rs0,
                            const std::string& sindex)
// produce a new rng ``rs'' uniquely identified by ``rs0'' and ``sindex''
// will not affect old rng ``rs0''
// the function should behave correctly even if ``rs'' is actually ``rs0''
{
  std::string data;
  if (ULONG_MAX == rs0.type) {
    data = ssprintf("[%lu] {%s}", rs0.index, sindex.c_str());
  } else {
    data = ssprintf("[%lu,%lu] {%s}", rs0.type, rs0.index, sindex.c_str());
  }
  const int nBlocks = (data.length() - 1) / 64 + 1;
  data.resize(nBlocks * 64, ' ');
  sha256::processBlock(rs.hash, rs0.hash, (const uint8_t*)data.c_str());
  for (int i = 1; i < nBlocks; ++i) {
    sha256::processBlock(rs.hash, rs.hash,
                         (const uint8_t*)data.c_str() + i * 64);
  }
  rs.numBytes = rs0.numBytes + nBlocks * 64;
  rs.type = ULONG_MAX;
  rs.index = 0;
  rs.cache[0] = 0;
  rs.cache[1] = 0;
  rs.cache[2] = 0;
  rs.gaussian = 0.0;
  rs.cacheAvail = 0;
  rs.gaussianAvail = false;
}

void split_rng_state(RngState& rs, const RngState& rs0,
                            const long sindex)
{
  split_rng_state(rs, rs0, show(sindex));
}

void compute_hash_with_input(uint32_t hash[8], const RngState& rs,
                                    const std::string& input)
{
  sha256::processInput(hash, rs.hash, rs.numBytes,
                       (const uint8_t*)input.c_str(), input.length());
}

uint64_t rand_gen(RngState& rs)
{
  qassert(0 <= rs.cacheAvail && rs.cacheAvail <= 3);
  rs.index += 1;
  if (rs.cacheAvail > 0) {
    rs.cacheAvail -= 1;
    uint64_t r = rs.cache[rs.cacheAvail];
    rs.cache[rs.cacheAvail] = 0;
    return r;
  } else {
    uint32_t hash[8];
    if (ULONG_MAX == rs.type) {
      compute_hash_with_input(hash, rs, ssprintf("[%lu]", rs.index));
    } else {
      compute_hash_with_input(hash, rs,
                              ssprintf("[%lu,%lu]", rs.type, rs.index));
    }
    rs.cache[0] = patchTwoUint32(hash[0], hash[1]);
    rs.cache[1] = patchTwoUint32(hash[2], hash[3]);
    rs.cache[2] = patchTwoUint32(hash[4], hash[5]);
    rs.cacheAvail = 3;
    return patchTwoUint32(hash[6], hash[7]);
  }
}

double u_rand_gen(RngState& rs, const double upper, const double lower)
{
  uint64_t u = rand_gen(rs);
  const double fac =
      1.0 / (256.0 * 256.0 * 256.0 * 256.0) / (256.0 * 256.0 * 256.0 * 256.0);
  return u * fac * (upper - lower) + lower;
}

double g_rand_gen(RngState& rs, const double center, const double sigma)
{
  rs.index += 1;
  if (rs.gaussianAvail) {
    rs.gaussianAvail = false;
    return rs.gaussian * sigma + center;
  } else {
    // pick 2 uniform numbers in the square extending from
    // -1 to 1 in each direction, see if they are in the
    // unit circle, and try again if they are not.
    int num_try = 1;
    double v1, v2, rsq;
    do {
      v1 = u_rand_gen(rs, 1.0, -1.0);
      v2 = u_rand_gen(rs, 1.0, -1.0);
      if ((num_try % 1000) == 0) {
        printf("g_rand_gen : WARNING num_try=%d v1=%e v2=%e\n", num_try, v1,
               v2);
      }
      rsq = v1 * v1 + v2 * v2;
      num_try++;
    } while ((num_try < 10000) && (rsq >= 1.0 || rsq == 0));
    if (num_try > 9999) {
      printf(
          "g_rand_gen : WARNING failed after 10000 tries (corrupted RNG?), "
          "returning ridiculous numbers (1e+10)\n");
      return 1e+10;
    }
    double fac = std::sqrt(-2.0 * std::log(rsq) / rsq);
    rs.gaussian = v1 * fac;
    rs.gaussianAvail = true;
    return v2 * fac * sigma + center;
  }
}

}  // namespace qlat
