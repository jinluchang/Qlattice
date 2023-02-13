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

#pragma once

#include <qlat-utils/config.h>
#include <qlat-utils/qcd-setting.h>
#include <qlat-utils/qacc.h>
#include <qlat-utils/assert.h>
#include <qlat-utils/array.h>
#include <qlat-utils/rng-state.h>
#include <qlat-utils/core.h>

#include <vector>
#include <iostream>
#include <sstream>

namespace qlat
{

// --------------------

qacc bool qisnan(const float& arg) { return std::isnan(arg); }

qacc bool qisnan(const double& arg) { return std::isnan(arg); }

template <class T>
bool qisnan(const std::complex<T>& arg)
{
  return qisnan(arg.real()) or qisnan(arg.imag());
}

template <class M, unsigned long N>
qacc bool qisnan(const array<M, N>& arr)
{
  for (int i = 0; i < (int)N; ++i) {
    if (qisnan(arr[i])) {
      return true;
    }
  }
  return false;
}

template <class M>
bool qisnan(const std::vector<M>& arr)
{
  for (size_t i = 0; i < arr.size(); ++i) {
    if (qisnan(arr[i])) {
      return true;
    }
  }
  return false;
}

// --------------------

qacc long modl(const long x, const long len)
{
  qassert(0 < len);
  const int m = x % len;
  if (0 <= m) {
    return m;
  } else {
    return m + len;
  }
}

qacc int mod(const int x, const int len)
{
  qassert(0 < len);
  const int m = x % len;
  if (0 <= m) {
    return m;
  } else {
    return m + len;
  }
}

qacc double mod(const double x, const double len)
{
  qassert(0 < len);
  const double m = x - trunc(x / len) * len;
  if (0 <= m) {
    return m;
  } else {
    return m + len;
  }
}

qacc int smod(const int x, const int len)
{
  qassert(0 < len);
  const int m = mod(x, len);
  if (m * 2 < len) {
    return m;
  } else {
    return m - len;
  }
}

qacc double smod(const double x, const double len)
{
  qassert(0 < len);
  const double m = mod(x, len);
  if (m * 2 < len) {
    return m;
  } else {
    return m - len;
  }
}

qacc double smod_sym(const double x, const double len,
                     const double eps = 1.0e-8)
{
  const double m = smod(x, len);
  if (std::abs(std::abs(m * 2) - len) < eps) {
    return 0;
  } else {
    return m;
  }
}

qacc int middle_mod(const int x, const int y, const int len)
{
  qassert(0 < len);
  const int xm = mod(x, len);
  const int ym = mod(y, len);
  if (xm <= ym) {
    const int r = smod(ym - xm, len);
    return mod(xm + r / 2, len);
  } else {
    const int r = smod(xm - ym, len);
    return mod(ym + r / 2, len);
  }
}

qacc double middle_mod(const double x, const double y, const double len)
{
  qassert(0 < len);
  const double xm = mod(x, len);
  const double ym = mod(y, len);
  if (xm <= ym) {
    const double r = smod(ym - xm, len);
    return mod(xm + r / 2, len);
  } else {
    const double r = smod(xm - ym, len);
    return mod(ym + r / 2, len);
  }
}

// --------------------

struct API Coordinate : public array<int, DIMN> {
  qacc Coordinate() { array<int, DIMN>::fill(0); }
  qacc Coordinate(int first, int second, int third, int fourth)
  {
    int* p = data();
    p[0] = first;
    p[1] = second;
    p[2] = third;
    p[3] = fourth;
  }
};

inline std::string show(const Coordinate& x)
{
  return ssprintf("%dx%dx%dx%d", x[0], x[1], x[2], x[3]);
}

qacc long product(const Coordinate& coor)
{
  long ret = 1;
  for (int i = 0; i < (int)coor.size(); i++) {
    ret *= coor[i];
  }
  return ret;
}

qacc Coordinate coordinate_from_index(long index, const Coordinate& size)
{
  Coordinate x;
  x[0] = index % size[0];
  index /= size[0];
  x[1] = index % size[1];
  index /= size[1];
  x[2] = index % size[2];
  index /= size[2];
  x[3] = index % size[3];
  return x;
}

qacc long index_from_coordinate(const Coordinate& x, const Coordinate& size)
{
  return ((((long)x[3] * (long)size[2]) + (long)x[2]) * (long)size[1] +
          (long)x[1]) *
             (long)size[0] +
         (long)x[0];
}

qacc int eo_from_coordinate(const Coordinate& xl)
{
  return 2 - (xl[0] + xl[1] + xl[2] + xl[3] + 16 * 1024 * 1024) % 2;
}

qacc bool operator==(const Coordinate& c1, const Coordinate& c2)
{
  return c1[0] == c2[0] and c1[1] == c2[1] and c1[2] == c2[2] and
         c1[3] == c2[3];
}

qacc bool operator!=(const Coordinate& c1, const Coordinate& c2)
{
  return !(c1 == c2);
}

qacc Coordinate operator+(const Coordinate& coor1, const Coordinate& coor2)
{
  return Coordinate(coor1[0] + coor2[0], coor1[1] + coor2[1],
                    coor1[2] + coor2[2], coor1[3] + coor2[3]);
}

qacc Coordinate operator-(const Coordinate& coor1, const Coordinate& coor2)
{
  return Coordinate(coor1[0] - coor2[0], coor1[1] - coor2[1],
                    coor1[2] - coor2[2], coor1[3] - coor2[3]);
}

qacc Coordinate operator-(const Coordinate& coor)
{
  return Coordinate(-coor[0], -coor[1], -coor[2], -coor[3]);
}

qacc Coordinate operator*(const int integer, const Coordinate& coor)
{
  return Coordinate(integer * coor[0], integer * coor[1], integer * coor[2],
                    integer * coor[3]);
}

qacc Coordinate operator*(const Coordinate& coor, const int integer)
{
  return integer * coor;
}

qacc Coordinate operator*(const Coordinate& coor1, const Coordinate& coor2)
{
  return Coordinate(coor1[0] * coor2[0], coor1[1] * coor2[1],
                    coor1[2] * coor2[2], coor1[3] * coor2[3]);
}

qacc Coordinate operator%(const Coordinate& coor1, const Coordinate& coor2)
{
  return Coordinate(coor1[0] % coor2[0], coor1[1] % coor2[1],
                    coor1[2] % coor2[2], coor1[3] % coor2[3]);
}

qacc Coordinate operator%(const Coordinate& coor, const int integer)
{
  return Coordinate(coor[0] % integer, coor[1] % integer, coor[2] % integer,
                    coor[3] % integer);
}

qacc Coordinate mod(const Coordinate& x, const Coordinate& size)
{
  Coordinate ret;
  ret[0] = mod(x[0], size[0]);
  ret[1] = mod(x[1], size[1]);
  ret[2] = mod(x[2], size[2]);
  ret[3] = mod(x[3], size[3]);
  return ret;
}

qacc Coordinate smod(const Coordinate& x, const Coordinate& size)
{
  Coordinate ret;
  ret[0] = smod(x[0], size[0]);
  ret[1] = smod(x[1], size[1]);
  ret[2] = smod(x[2], size[2]);
  ret[3] = smod(x[3], size[3]);
  return ret;
}

qacc Coordinate middle_mod(const Coordinate& x, const Coordinate& y,
                           const Coordinate& size)
{
  Coordinate ret;
  ret[0] = middle_mod(x[0], y[0], size[0]);
  ret[1] = middle_mod(x[1], y[1], size[1]);
  ret[2] = middle_mod(x[2], y[2], size[2]);
  ret[3] = middle_mod(x[3], y[3], size[3]);
  return ret;
}

qacc Coordinate c_rand_gen(RngState& rs, const Coordinate& size)
{
  const long total_vol = product(size);
  const long ri = rand_gen(rs) % total_vol;
  return coordinate_from_index(ri, size);
}

// --------------------

struct API CoordinateD : public array<double, DIMN> {
  qacc CoordinateD() { memset(this, 0, sizeof(CoordinateD)); }
  qacc CoordinateD(const array<double, DIMN>& arr)
  {
    CoordinateD& c = *this;
    c = arr;
    qassert(false == qisnan(c));
  }
  qacc CoordinateD(const double x0, const double x1, const double x2,
                   const double x3)
  {
    qassert(DIMN == 4);
    CoordinateD& c = *this;
    c[0] = x0;
    c[1] = x1;
    c[2] = x2;
    c[3] = x3;
    qassert(false == qisnan(c));
  }
  qacc CoordinateD(const Coordinate& x)
  {
    CoordinateD& c = *this;
    for (int i = 0; i < DIMN; ++i) {
      c[i] = x[i];
    }
  }
};

inline std::string show(const CoordinateD& c)
{
  return ssprintf("(%23.16e,%23.16e,%23.16e,%23.16e)", c[0], c[1], c[2], c[3]);
}

// --------------------

qacc long sqr(const Coordinate& xg)
{
  return sqr((long)xg[0]) + sqr((long)xg[1]) + sqr((long)xg[2]) +
         sqr((long)xg[3]);
}

qacc Coordinate operator/(const Coordinate& coor, const int integer)
{
  return Coordinate(coor[0] / integer, coor[1] / integer, coor[2] / integer,
                    coor[3] / integer);
}

qacc Coordinate operator/(const Coordinate& coor1, const Coordinate& coor2)
{
  return Coordinate(coor1[0] / coor2[0], coor1[1] / coor2[1],
                    coor1[2] / coor2[2], coor1[3] / coor2[3]);
}

qacc int sum(const Coordinate& coor)
{
  int ret = 0;
  for (int i = 0; i < (int)coor.size(); i++) {
    ret += coor[i];
  }
  return ret;
}

qacc int parity(const Coordinate& coor)
{
  return 2 - sum(coor) % 2;  // 2 for even, 1 for odd
}

inline void regularize(Coordinate& coor, const Coordinate& regularizer)
{
#ifndef SUPPRESS_REG_COOR
  // #ifndef USE_SINGLE_NODE // SINGLE_NODE mode (or !USE_MULTI_NODE) does NOT
  // work (yet?) so this is meaningless. Thus it is commented out.
  warn("use regular_coordinate");
#endif
  for (int mu = 0; mu < DIMN; mu++) {
    coor[mu] = (coor[mu] % regularizer[mu] + regularizer[mu]) % regularizer[mu];
  }
}

qacc Coordinate coordinate_shifts(const Coordinate& x) { return x; }

qacc Coordinate coordinate_shifts(const Coordinate& x, const int dir)
{
  Coordinate xsh = x;
  qassert(-DIMN <= dir && dir < DIMN);
  if (0 <= dir) {
    xsh[dir] += 1;
  } else {
    xsh[-dir - 1] -= 1;
  }
  return xsh;
}

qacc Coordinate coordinate_shifts(const Coordinate& x, const int dir1,
                                  const int dir2)
{
  Coordinate xsh = x;
  qassert(-DIMN <= dir1 && dir1 < DIMN);
  qassert(-DIMN <= dir2 && dir2 < DIMN);
  if (0 <= dir1) {
    xsh[dir1] += 1;
  } else {
    xsh[-dir1 - 1] -= 1;
  }
  if (0 <= dir2) {
    xsh[dir2] += 1;
  } else {
    xsh[-dir2 - 1] -= 1;
  }
  return xsh;
}

qacc Coordinate coordinate_shifts(const Coordinate& x, const int dir1,
                                  const int dir2, const int dir3)
{
  Coordinate xsh = x;
  qassert(-DIMN <= dir1 && dir1 < DIMN);
  qassert(-DIMN <= dir2 && dir2 < DIMN);
  qassert(-DIMN <= dir3 && dir3 < DIMN);
  if (0 <= dir1) {
    xsh[dir1] += 1;
  } else {
    xsh[-dir1 - 1] -= 1;
  }
  if (0 <= dir2) {
    xsh[dir2] += 1;
  } else {
    xsh[-dir2 - 1] -= 1;
  }
  if (0 <= dir3) {
    xsh[dir3] += 1;
  } else {
    xsh[-dir3 - 1] -= 1;
  }
  return xsh;
}

qacc Coordinate coordinate_shifts(const Coordinate& x, const int dir1,
                                  const int dir2, const int dir3,
                                  const int dir4)
{
  Coordinate xsh = x;
  qassert(-DIMN <= dir1 && dir1 < DIMN);
  qassert(-DIMN <= dir2 && dir2 < DIMN);
  qassert(-DIMN <= dir3 && dir3 < DIMN);
  qassert(-DIMN <= dir4 && dir4 < DIMN);
  if (0 <= dir1) {
    xsh[dir1] += 1;
  } else {
    xsh[-dir1 - 1] -= 1;
  }
  if (0 <= dir2) {
    xsh[dir2] += 1;
  } else {
    xsh[-dir2 - 1] -= 1;
  }
  if (0 <= dir3) {
    xsh[dir3] += 1;
  } else {
    xsh[-dir3 - 1] -= 1;
  }
  if (0 <= dir4) {
    xsh[dir4] += 1;
  } else {
    xsh[-dir4 - 1] -= 1;
  }
  return xsh;
}

inline Coordinate coordinate_shifts(const Coordinate& x,
                                    const std::vector<int>& path)
{
  Coordinate ret = x;
  for (int i = 0; i < (int)path.size(); ++i) {
    const int dir = path[i];
    qassert(-DIMN <= dir && dir < DIMN);
    ret = coordinate_shifts(ret, dir);
  }
  return ret;
}

template <class CharT, class Traits>
std::basic_ostream<CharT, Traits>& operator<<(
    std::basic_ostream<CharT, Traits>& os, const Coordinate& coor)
{
  std::basic_ostringstream<CharT, Traits> os_;
  os_.flags(os.flags());
  os_.imbue(os.getloc());
  os_.precision(os.precision());
  os_.setf(std::ios::showpos);
  os_ << "(" << coor[0] << ", " << coor[1] << ", " << coor[2] << ", " << coor[3]
      << ")";
  return os << os_.str();
}

}  // namespace qlat
