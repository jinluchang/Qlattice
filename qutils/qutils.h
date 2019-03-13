#pragma once

#ifndef OLD_CPP
#include <array>
#endif

#include <vector>
#include <complex>

namespace qutils
{  //

typedef std::complex<double> Complex;

typedef std::complex<float> ComplexF;

const double PI = 3.141592653589793;

const Complex ii(0, 1);

template <class T>
T sqr(const T& x)
{
  return x * x;
}

template <class M>
void clear(std::vector<M>& vec)
{
  std::vector<M> empty;
  swap(empty, vec);
}

template <class M>
std::array<M, 0> make_array()
{
  std::array<M, 0> arr;
  return arr;
}

template <class M>
std::array<M, 1> make_array(const M& x)
{
  std::array<M, 1> arr;
  arr[0] = x;
  return arr;
}

template <class M>
std::array<M, 2> make_array(const M& x, const M& x1)
{
  std::array<M, 2> arr;
  arr[0] = x;
  arr[1] = x1;
  return arr;
}

template <class M>
std::array<M, 3> make_array(const M& x, const M& x1, const M& x2)
{
  std::array<M, 3> arr;
  arr[0] = x;
  arr[1] = x1;
  arr[2] = x2;
  return arr;
}

template <class M>
std::array<M, 4> make_array(const M& x, const M& x1, const M& x2, const M& x3)
{
  std::array<M, 4> arr;
  arr[0] = x;
  arr[1] = x1;
  arr[2] = x2;
  arr[3] = x3;
  return arr;
}

template <class M>
std::array<M, 5> make_array(const M& x, const M& x1, const M& x2, const M& x3,
                            const M& x4)
{
  std::array<M, 5> arr;
  arr[0] = x;
  arr[1] = x1;
  arr[2] = x2;
  arr[3] = x3;
  arr[4] = x4;
  return arr;
}

template <class M>
std::array<M, 6> make_array(const M& x, const M& x1, const M& x2, const M& x3,
                            const M& x4, const M& x5)
{
  std::array<M, 6> arr;
  arr[0] = x;
  arr[1] = x1;
  arr[2] = x2;
  arr[3] = x3;
  arr[4] = x4;
  arr[5] = x5;
  return arr;
}

template <class M>
std::array<M, 7> make_array(const M& x, const M& x1, const M& x2, const M& x3,
                            const M& x4, const M& x5, const M& x6)
{
  std::array<M, 7> arr;
  arr[0] = x;
  arr[1] = x1;
  arr[2] = x2;
  arr[3] = x3;
  arr[4] = x4;
  arr[5] = x5;
  arr[6] = x6;
  return arr;
}

template <class M>
std::array<M, 8> make_array(const M& x, const M& x1, const M& x2, const M& x3,
                            const M& x4, const M& x5, const M& x6, const M& x7)
{
  std::array<M, 8> arr;
  arr[0] = x;
  arr[1] = x1;
  arr[2] = x2;
  arr[3] = x3;
  arr[4] = x4;
  arr[5] = x5;
  arr[6] = x6;
  arr[7] = x7;
  return arr;
}

template <class M>
std::array<M, 9> make_array(const M& x, const M& x1, const M& x2, const M& x3,
                            const M& x4, const M& x5, const M& x6, const M& x7,
                            const M& x8)
{
  std::array<M, 9> arr;
  arr[0] = x;
  arr[1] = x1;
  arr[2] = x2;
  arr[3] = x3;
  arr[4] = x4;
  arr[5] = x5;
  arr[6] = x6;
  arr[7] = x7;
  arr[8] = x8;
  return arr;
}

template <class M>
std::array<M, 10> make_array(const M& x, const M& x1, const M& x2, const M& x3,
                             const M& x4, const M& x5, const M& x6, const M& x7,
                             const M& x8, const M& x9)
{
  std::array<M, 10> arr;
  arr[0] = x;
  arr[1] = x1;
  arr[2] = x2;
  arr[3] = x3;
  arr[4] = x4;
  arr[5] = x5;
  arr[6] = x6;
  arr[7] = x7;
  arr[8] = x8;
  arr[9] = x9;
  return arr;
}

}  // namespace qutils
