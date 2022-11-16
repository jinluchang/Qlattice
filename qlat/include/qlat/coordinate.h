#pragma once

#include <qlat/setup.h>

#include <iostream>
#include <sstream>

namespace qlat
{  //

qacc long sqr(const qlat::Coordinate& xg)
{
  return sqr((long)xg[0]) + sqr((long)xg[1]) + sqr((long)xg[2]) +
         sqr((long)xg[3]);
}

qacc Coordinate operator/(const Coordinate &coor, const int integer)
{
  return Coordinate(coor[0] / integer, coor[1] / integer, coor[2] / integer,
                    coor[3] / integer);
}

qacc Coordinate operator/(const Coordinate &coor1, const Coordinate &coor2)
{
  return Coordinate(coor1[0] / coor2[0], coor1[1] / coor2[1],
                    coor1[2] / coor2[2], coor1[3] / coor2[3]);
}

qacc int sum(const Coordinate &coor)
{
  int ret = 0;
  for (int i = 0; i < (int)coor.size(); i++) {
    ret += coor[i];
  }
  return ret;
}

qacc int parity(const Coordinate &coor)
{
  return 2 - sum(coor) % 2;  // 2 for even, 1 for odd
}

inline void regularize(Coordinate &coor, const Coordinate &regularizer)
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

qacc Coordinate coordinate_shifts(const Coordinate &x) { return x; }

qacc Coordinate coordinate_shifts(const Coordinate &x, const int dir)
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

qacc Coordinate coordinate_shifts(const Coordinate &x, const int dir1,
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

qacc Coordinate coordinate_shifts(const Coordinate &x, const int dir1,
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

qacc Coordinate coordinate_shifts(const Coordinate &x, const int dir1,
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

inline Coordinate coordinate_shifts(const Coordinate &x,
                                    const std::vector<int> path)
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
std::basic_ostream<CharT, Traits> &operator<<(
    std::basic_ostream<CharT, Traits> &os, const Coordinate coor)
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
