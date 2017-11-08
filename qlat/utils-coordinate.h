// vim: set ts=2 sw=2 expandtab:

#pragma once

#include <qlat/utils.h>
#include <qlat/coordinate.h>
#include <qlat/coordinate-d.h>

QLAT_START_NAMESPACE

inline int mod(const int x, const int len)
{
  qassert(0 < len);
  const int m = x % len;
  if (0 <= m) {
    return m;
  } else {
    return m + len;
  }
}

inline double mod(const double x, const double len)
{
  qassert(0 < len);
  const double m = x - trunc(x / len) * len;
  if (0 <= m) {
    return m;
  } else {
    return m + len;
  }
}

inline int smod(const int x, const int len)
{
  qassert(0 < len);
  const int m = mod(x, len);
  if (m * 2 < len) {
    return m;
  } else {
    return m - len;
  }
}

inline double smod(const double x, const double len)
{
  qassert(0 < len);
  const double m = mod(x, len);
  if (m * 2 < len) {
    return m;
  } else {
    return m - len;
  }
}

inline int middle_mod(const int x, const int y, const int len)
{
  qassert(0 < len);
  const int xm = mod(x, len);
  const int ym = mod(y, len);
  if (xm <= ym) {
    const int r = smod(ym - xm, len);
    return mod(xm + r/2, len);
  } else {
    const int r = smod(xm - ym, len);
    return mod(ym + r/2, len);
  }
}

inline double middle_mod(const double x, const double y, const double len)
{
  qassert(0 < len);
  const double xm = mod(x, len);
  const double ym = mod(y, len);
  if (xm <= ym) {
    const double r = smod(ym - xm, len);
    return mod(xm + r/2, len);
  } else {
    const double r = smod(xm - ym, len);
    return mod(ym + r/2, len);
  }
}

inline Coordinate mod(const Coordinate& x, const Coordinate& size)
{
  Coordinate ret;
  ret[0] = mod(x[0], size[0]);
  ret[1] = mod(x[1], size[1]);
  ret[2] = mod(x[2], size[2]);
  ret[3] = mod(x[3], size[3]);
  return ret;
}

inline Coordinate smod(const Coordinate& x, const Coordinate& size)
{
  Coordinate ret;
  ret[0] = smod(x[0], size[0]);
  ret[1] = smod(x[1], size[1]);
  ret[2] = smod(x[2], size[2]);
  ret[3] = smod(x[3], size[3]);
  return ret;
}

inline CoordinateD smod(const CoordinateD& x, const CoordinateD& size)
{
  CoordinateD ret;
  ret[0] = smod(x[0], size[0]);
  ret[1] = smod(x[1], size[1]);
  ret[2] = smod(x[2], size[2]);
  ret[3] = smod(x[3], size[3]);
  return ret;
}

inline Coordinate middle_mod(const Coordinate& x, const Coordinate& y, const Coordinate& size)
{
  Coordinate ret;
  ret[0] = middle_mod(x[0], y[0], size[0]);
  ret[1] = middle_mod(x[1], y[1], size[1]);
  ret[2] = middle_mod(x[2], y[2], size[2]);
  ret[3] = middle_mod(x[3], y[3], size[3]);
  return ret;
}

inline CoordinateD middle_mod(const CoordinateD& x, const CoordinateD& y, const CoordinateD& size)
{
  CoordinateD ret;
  ret[0] = middle_mod(x[0], y[0], size[0]);
  ret[1] = middle_mod(x[1], y[1], size[1]);
  ret[2] = middle_mod(x[2], y[2], size[2]);
  ret[3] = middle_mod(x[3], y[3], size[3]);
  return ret;
}

inline bool is_reaching_edge_coordinate(const Coordinate& x, const Coordinate& size)
{
  return std::abs(x[0]) * 2 == size[0]
    || std::abs(x[1]) * 2 == size[1]
    || std::abs(x[2]) * 2 == size[2]
    || std::abs(x[3]) * 2 == size[3];
}

inline bool is_outside_coordinate(const Coordinate& x, const Coordinate& size)
{
  return std::abs(x[0]) * 2 > size[0]
    || std::abs(x[1]) * 2 > size[1]
    || std::abs(x[2]) * 2 > size[2]
    || std::abs(x[3]) * 2 > size[3];
}

inline long distance_sq_relative_coordinate_g(const Coordinate& xg)
{
  return sqr((long)xg[0]) + sqr((long)xg[1]) + sqr((long)xg[2]) + sqr((long)xg[3]);
}

inline double distance_relative_coordinate_g(const Coordinate& xg)
{
  return sqrt(distance_sq_relative_coordinate_g(xg));
}

inline Coordinate coordinate_from_index(long index, const Coordinate& size)
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

inline long index_from_coordinate(const Coordinate& x, const Coordinate& size)
{
  return (((x[3] * size[2]) + x[2]) * size[1] + x[1]) * size[0] + x[0];
}

inline Coordinate regular_coordinate(const Coordinate& x, const Coordinate& size)
{
  return mod(x, size);
}

inline Coordinate relative_coordinate(const Coordinate& x, const Coordinate& size)
{
  return smod(x, size);
}

inline void regularize_coordinate(Coordinate& x, const Coordinate& size)
{
  x = regular_coordinate(x, size);
}

inline CoordinateD relative_coordinate(const CoordinateD& x, const CoordinateD& size)
{
  return smod(x, size);
}

inline Coordinate middle_coordinate(const Coordinate& x, const Coordinate& y, const Coordinate& size)
{
  return middle_mod(x, y, size);
}

inline CoordinateD middle_coordinate(const CoordinateD& x, const CoordinateD& y, const CoordinateD& size)
{
  return middle_mod(x, y, size);
}

QLAT_END_NAMESPACE

namespace qshow {

inline std::string show(const qlat::Coordinate& x) {
  return ssprintf("%dx%dx%dx%d", x[0], x[1], x[2], x[3]);
}

}

#ifndef USE_NAMESPACE
using namespace qshow;
#endif
