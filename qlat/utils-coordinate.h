// vim: set ts=2 sw=2 expandtab:

#pragma once

#include <qlat/coordinate-d.h>
#include <qlat/coordinate.h>
#include <qlat/utils.h>

QLAT_START_NAMESPACE

inline std::string show(const qlat::Coordinate& x)
{
  return ssprintf("%dx%dx%dx%d", x[0], x[1], x[2], x[3]);
}

inline long sqr(const qlat::Coordinate& xg)
{
  return sqr((long)xg[0]) + sqr((long)xg[1]) + sqr((long)xg[2]) +
         sqr((long)xg[3]);
}

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

inline double smod_sym(const double x, const double len,
                       const double eps = 1.0e-8)
{
  const double m = smod(x, len);
  if (std::abs(std::abs(m * 2) - len) < eps) {
    return 0;
  } else {
    return m;
  }
}

inline int middle_mod(const int x, const int y, const int len)
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

inline double middle_mod(const double x, const double y, const double len)
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

inline Coordinate mod(const Coordinate& x, const Coordinate& size)
{
  Coordinate ret;
  ret[0] = mod(x[0], size[0]);
  ret[1] = mod(x[1], size[1]);
  ret[2] = mod(x[2], size[2]);
  ret[3] = mod(x[3], size[3]);
  return ret;
}

inline CoordinateD mod(const CoordinateD& x, const CoordinateD& size)
{
  CoordinateD ret;
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

inline Coordinate middle_mod(const Coordinate& x, const Coordinate& y,
                             const Coordinate& size)
{
  Coordinate ret;
  ret[0] = middle_mod(x[0], y[0], size[0]);
  ret[1] = middle_mod(x[1], y[1], size[1]);
  ret[2] = middle_mod(x[2], y[2], size[2]);
  ret[3] = middle_mod(x[3], y[3], size[3]);
  return ret;
}

inline CoordinateD middle_mod(const CoordinateD& x, const CoordinateD& y,
                              const CoordinateD& size)
{
  CoordinateD ret;
  ret[0] = middle_mod(x[0], y[0], size[0]);
  ret[1] = middle_mod(x[1], y[1], size[1]);
  ret[2] = middle_mod(x[2], y[2], size[2]);
  ret[3] = middle_mod(x[3], y[3], size[3]);
  return ret;
}

inline bool is_reaching_edge_coordinate(const Coordinate& x,
                                        const Coordinate& size)
{
  return std::abs(x[0]) * 2 == size[0] || std::abs(x[1]) * 2 == size[1] ||
         std::abs(x[2]) * 2 == size[2] || std::abs(x[3]) * 2 == size[3];
}

inline bool is_outside_coordinate(const Coordinate& x, const Coordinate& size)
{
  return std::abs(x[0]) * 2 > size[0] || std::abs(x[1]) * 2 > size[1] ||
         std::abs(x[2]) * 2 > size[2] || std::abs(x[3]) * 2 > size[3];
}

inline long distance_sq_relative_coordinate_g(const Coordinate& xg)
{
  return sqr(xg);
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
  return ((((long)x[3] * (long)size[2]) + (long)x[2]) * (long)size[1] +
          (long)x[1]) *
             (long)size[0] +
         (long)x[0];
}

inline Coordinate regular_coordinate(const Coordinate& x,
                                     const Coordinate& size)
{
  return mod(x, size);
}

inline Coordinate relative_coordinate(const Coordinate& x,
                                      const Coordinate& size)
{
  return smod(x, size);
}

inline void regularize_coordinate(Coordinate& x, const Coordinate& size)
{
  x = regular_coordinate(x, size);
}

inline CoordinateD relative_coordinate(const CoordinateD& x,
                                       const CoordinateD& size)
{
  return smod(x, size);
}

inline Coordinate middle_coordinate(const Coordinate& x, const Coordinate& y,
                                    const Coordinate& size)
{
  return middle_mod(x, y, size);
}

inline CoordinateD middle_coordinate(const CoordinateD& x, const CoordinateD& y,
                                     const CoordinateD& size)
{
  return middle_mod(x, y, size);
}

struct EpsilonTensorTable {
  int tensor[4][4][4][4];
  //
  EpsilonTensorTable() { init(); }
  //
  void init()
  {
    std::memset(this, 0, sizeof(tensor));
    setv(0, 1, 2, 3);
    setv(0, 2, 3, 1);
    setv(0, 3, 1, 2);
  }
  //
  void setv(const int a, const int b, const int c, const int d)
  {
    set(a, b, c, d, 1);
    set(a, b, d, c, -1);
  }
  void set(const int a, const int b, const int c, const int d, const int val)
  {
    tensor[a][b][c][d] = val;
    tensor[b][c][d][a] = -val;
    tensor[c][d][a][b] = val;
    tensor[d][a][b][c] = -val;
  }
};

inline int epsilon_tensor(const int a, const int b, const int c, const int d)
{
  static EpsilonTensorTable table;
  return table.tensor[a][b][c][d];
}

inline int epsilon_tensor(const int i, const int j, const int k)
{
  return epsilon_tensor(i, j, k, 3);
}

inline Coordinate read_coordinate(const std::string& str)
{
  long x = 0;
  long y = 0;
  long z = 0;
  long t = 0;
  long cur = 0;
  char c;
  assert(parse_long(x, cur, str));
  assert(parse_char(c, cur, str) and (c == 'x' or c == ','));
  assert(parse_long(y, cur, str));
  assert(parse_char(c, cur, str) and (c == 'x' or c == ','));
  assert(parse_long(z, cur, str));
  assert(parse_char(c, cur, str) and (c == 'x' or c == ','));
  assert(parse_long(t, cur, str));
  return Coordinate(x, y, z, t);
}

QLAT_END_NAMESPACE
