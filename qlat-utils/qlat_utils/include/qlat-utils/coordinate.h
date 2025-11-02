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

#include <qlat-utils/array.h>
#include <qlat-utils/qassert.h>
#include <qlat-utils/config.h>
#include <qlat-utils/core.h>
#include <qlat-utils/qacc.h>
#include <qlat-utils/qcd-setting.h>
#include <qlat-utils/rng-state.h>
#include <qlat-utils/timer.h>

#include <iostream>
#include <sstream>
#include <vector>

namespace qlat
{

qacc Long modl(const Long x, const Long len)
{
  qassert(0 < len);
  const Int m = x % len;
  if (0 <= m) {
    return m;
  } else {
    return m + len;
  }
}

qacc Int mod(const Int x, const Int len)
{
  qassert(0 < len);
  const Int m = x % len;
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

qacc Int smod(const Int x, const Int len)
{
  qassert(0 < len);
  const Int m = mod(x, len);
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

qacc Int smod_sym(const Int x, const Int len)
{
  qassert(0 < len);
  const Int m = smod(x, len);
  if (abs(m * 2) == len) {
    return 0;
  } else {
    return m;
  }
}

qacc double smod_sym(const double x, const double len,
                     const double eps = 1.0e-8)
{
  const double m = smod(x, len);
  if (abs(abs(m * 2) - len) < eps) {
    return 0;
  } else {
    return m;
  }
}

qacc Int middle_mod(const Int x, const Int y, const Int len)
{
  qassert(0 < len);
  const Int xm = mod(x, len);
  const Int ym = mod(y, len);
  if (xm <= ym) {
    const Int r = smod(ym - xm, len);
    return mod(xm + r / 2, len);
  } else {
    const Int r = smod(xm - ym, len);
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

inline std::string show(const Coordinate& x)
{
  return ssprintf("%dx%dx%dx%d", x[0], x[1], x[2], x[3]);
}

qacc Long product(const Coordinate& coor)
{
  Long ret = 1;
  for (Int i = 0; i < (Int)coor.size(); i++) {
    ret *= coor[i];
  }
  return ret;
}

qacc Coordinate coordinate_from_index(Long index, const Coordinate& size)
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

qacc Long index_from_coordinate(const Coordinate& x, const Coordinate& size)
{
  return ((((Long)x[3] * (Long)size[2]) + (Long)x[2]) * (Long)size[1] +
          (Long)x[1]) *
             (Long)size[0] +
         (Long)x[0];
}

qacc Int eo_from_coordinate(const Coordinate& xl)
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

qacc Coordinate operator*(const Int integer, const Coordinate& coor)
{
  return Coordinate(integer * coor[0], integer * coor[1], integer * coor[2],
                    integer * coor[3]);
}

qacc Coordinate operator*(const Coordinate& coor, const Int integer)
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

qacc Coordinate operator%(const Coordinate& coor, const Int integer)
{
  return Coordinate(coor[0] % integer, coor[1] % integer, coor[2] % integer,
                    coor[3] % integer);
}

qacc Coordinate operator/(const Coordinate& coor, const Int integer)
{
  return Coordinate(coor[0] / integer, coor[1] / integer, coor[2] / integer,
                    coor[3] / integer);
}

qacc Coordinate operator/(const Coordinate& coor1, const Coordinate& coor2)
{
  return Coordinate(coor1[0] / coor2[0], coor1[1] / coor2[1],
                    coor1[2] / coor2[2], coor1[3] / coor2[3]);
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

qacc Coordinate smod_sym(const Coordinate& x, const Coordinate& size)
{
  Coordinate ret;
  ret[0] = smod_sym(x[0], size[0]);
  ret[1] = smod_sym(x[1], size[1]);
  ret[2] = smod_sym(x[2], size[2]);
  ret[3] = smod_sym(x[3], size[3]);
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
  const Long total_vol = product(size);
  const Long ri = rand_gen(rs) % total_vol;
  return coordinate_from_index(ri, size);
}

// --------------------

inline std::string show(const CoordinateD& c)
{
  return ssprintf("(%23.16e,%23.16e,%23.16e,%23.16e)", c[0], c[1], c[2], c[3]);
}

// --------------------

qacc Long sqr(const Coordinate& xg)
{
  return sqr((Long)xg[0]) + sqr((Long)xg[1]) + sqr((Long)xg[2]) +
         sqr((Long)xg[3]);
}

qacc Int sum(const Coordinate& coor)
{
  Int ret = 0;
  for (Int i = 0; i < (Int)coor.size(); i++) {
    ret += coor[i];
  }
  return ret;
}

qacc Int parity(const Coordinate& coor)
{
  return 2 - sum(coor) % 2;  // 2 for even, 1 for odd
}

inline void regularize(Coordinate& coor, const Coordinate& regularizer)
{
#ifndef SUPPRESS_REG_COOR
  warn("use regular_coordinate");
#endif
  for (Int mu = 0; mu < DIMN; mu++) {
    coor[mu] = (coor[mu] % regularizer[mu] + regularizer[mu]) % regularizer[mu];
  }
}

qacc Coordinate coordinate_shifts(const Coordinate& x) { return x; }

qacc Coordinate coordinate_shifts(const Coordinate& x, const Int dir)
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

qacc Coordinate coordinate_shifts(const Coordinate& x, const Int dir1,
                                  const Int dir2)
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

qacc Coordinate coordinate_shifts(const Coordinate& x, const Int dir1,
                                  const Int dir2, const Int dir3)
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

qacc Coordinate coordinate_shifts(const Coordinate& x, const Int dir1,
                                  const Int dir2, const Int dir3,
                                  const Int dir4)
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
                                    const std::vector<Int>& path)
{
  Coordinate ret = x;
  for (Int i = 0; i < (Int)path.size(); ++i) {
    const Int dir = path[i];
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
