// vim: set ts=2 sw=2 expandtab:

#pragma once

#include <qlat-utils/coordinate.h>
#include <qlat/coordinate-d.h>

namespace qlat
{  //

qacc bool is_reaching_edge_coordinate(const Coordinate& x,
                                      const Coordinate& size)
{
  return std::abs(x[0]) * 2 == size[0] || std::abs(x[1]) * 2 == size[1] ||
         std::abs(x[2]) * 2 == size[2] || std::abs(x[3]) * 2 == size[3];
}

qacc bool is_outside_coordinate(const Coordinate& x, const Coordinate& size)
{
  return std::abs(x[0]) * 2 > size[0] || std::abs(x[1]) * 2 > size[1] ||
         std::abs(x[2]) * 2 > size[2] || std::abs(x[3]) * 2 > size[3];
}

qacc long distance_sq_relative_coordinate_g(const Coordinate& xg)
{
  return sqr(xg);
}

qacc double distance_relative_coordinate_g(const Coordinate& xg)
{
  return sqrt(distance_sq_relative_coordinate_g(xg));
}

qacc Coordinate regular_coordinate(const Coordinate& x, const Coordinate& size)
{
  return mod(x, size);
}

qacc Coordinate relative_coordinate(const Coordinate& x, const Coordinate& size)
{
  return smod(x, size);
}

qacc void regularize_coordinate(Coordinate& x, const Coordinate& size)
{
  x = regular_coordinate(x, size);
}

qacc CoordinateD relative_coordinate(const CoordinateD& x,
                                     const CoordinateD& size)
{
  return smod(x, size);
}

qacc Coordinate middle_coordinate(const Coordinate& x, const Coordinate& y,
                                  const Coordinate& size)
{
  return middle_mod(x, y, size);
}

qacc CoordinateD middle_coordinate(const CoordinateD& x, const CoordinateD& y,
                                   const CoordinateD& size)
{
  return middle_mod(x, y, size);
}

struct API EpsilonTensorTable {
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

API inline int epsilon_tensor(const int a, const int b, const int c,
                              const int d)
{
  static EpsilonTensorTable table;
  return table.tensor[a][b][c][d];
}

inline int epsilon_tensor(const int i, const int j, const int k)
{
  return epsilon_tensor(i, j, k, 3);
}

qacc int epsilon_tensor_acc(const int i, const int j, const int k)
{
  if (i == 0 and j == 1 and k == 2) {
    return 1;
  } else if (i == 1 and j == 2 and k == 0) {
    return 1;
  } else if (i == 2 and j == 0 and k == 1) {
    return 1;
  } else if (i == 2 and j == 1 and k == 0) {
    return -1;
  } else if (i == 1 and j == 0 and k == 2) {
    return -1;
  } else if (i == 0 and j == 2 and k == 1) {
    return -1;
  }
  return 0;
}

qacc int epsilon_tensor_acc(const int a, const int b, const int c, const int d)
{
  if (d == 3) {
    return epsilon_tensor_acc(a, b, c);
  } else if (c == 3) {
    return -epsilon_tensor_acc(a, b, d);
  } else if (b == 3) {
    return epsilon_tensor_acc(a, c, d);
  } else if (a == 3) {
    return -epsilon_tensor_acc(b, c, d);
  }
  return 0;
}

inline Coordinate read_coordinate(const std::string& str)
{
  long x = 0;
  long y = 0;
  long z = 0;
  long t = 0;
  long cur = 0;
  char c;
  qassert(parse_long(x, cur, str));
  qassert(parse_char(c, cur, str) and (c == 'x' or c == ','));
  qassert(parse_long(y, cur, str));
  qassert(parse_char(c, cur, str) and (c == 'x' or c == ','));
  qassert(parse_long(z, cur, str));
  qassert(parse_char(c, cur, str) and (c == 'x' or c == ','));
  qassert(parse_long(t, cur, str));
  return Coordinate(x, y, z, t);
}

}  // namespace qlat
