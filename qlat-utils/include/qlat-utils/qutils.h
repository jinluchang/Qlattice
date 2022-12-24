#pragma once

#include <qlat-utils/core.h>

#ifdef QLAT_USE_MACHINE_ENDIAN_H
#include <machine/endian.h>
#else
#include <endian.h>
#endif

#include <vector>
#include <map>
#include <set>

namespace qlat
{  //

// -------------------

qacc void set_zero(char& x) { x = 0; }

qacc void set_zero(int8_t& x) { x = 0; }

qacc void set_zero(long& x) { x = 0; }

qacc void set_zero(long long& x) { x = 0; }

qacc void set_zero(double& x) { x = 0; }

qacc void set_zero(float& x) { x = 0; }

qacc void set_zero(Complex& x) { x = 0; }

qacc void set_zero(ComplexF& x) { x = 0; }

template <class M, unsigned long N>
qacc void set_zero(array<M, N>& arr)
{
  long size = N * sizeof(M);
  std::memset((void*)arr.data(), 0, size);
}

template <class M>
void set_zero(std::vector<M>& vec)
{
  long size = vec.size() * sizeof(M);
  std::memset((void*)vec.data(), 0, size);
}

// -------------------

qacc void set_unit(char& x, const long& coef = 1) { x = coef; }

qacc void set_unit(char& x, const Complex& coef = 1) { x = coef.real(); }

qacc void set_unit(int8_t& x, const long& coef = 1) { x = coef; }

qacc void set_unit(int8_t& x, const Complex& coef = 1) { x = coef.real(); }

qacc void set_unit(long& x, const long& coef = 1) { x = coef; }

qacc void set_unit(long& x, const Complex& coef) { x = coef.real(); }

qacc void set_unit(long long& x, const long& coef = 1) { x = coef; }

qacc void set_unit(long long& x, const Complex& coef) { x = coef.real(); }

qacc void set_unit(float& x, const double& coef = 1.0) { x = coef; }

qacc void set_unit(float& x, const Complex& coef) { x = coef.real(); }

qacc void set_unit(double& x, const double& coef = 1.0) { x = coef; }

qacc void set_unit(double& x, const Complex& coef) { x = coef.real(); }

qacc void set_unit(Complex& x, const Complex& coef = 1.0) { x = coef; }

qacc void set_unit(ComplexF& x, const Complex& coef = 1.0) { x = coef; }

// -------------------

qacc double qnorm(const double& x) { return x * x; }

qacc double qnorm(const double& x, const double& y) { return x * y; }

template <class T, size_t N>
qacc double qnorm(const array<T, N>& mm)
{
  double sum = 0.0;
  for (size_t i = 0; i < N; ++i) {
    sum += qnorm(mm[i]);
  }
  return sum;
}

template <class T>
double qnorm(const std::vector<T>& mm)
{
  double sum = 0.0;
  for (size_t i = 0; i < mm.size(); ++i) {
    sum += qnorm(mm[i]);
  }
  return sum;
}

// -------------------

template <class M>
void clear(std::vector<M>& vec)
{
  std::vector<M> empty;
  std::swap(empty, vec);
}

// -------------------

qacc bool is_big_endian()
{
#if defined(__BYTE_ORDER) && (__BYTE_ORDER != 0) && \
    (__BYTE_ORDER == __BIG_ENDIAN)
  return true;
#else
  return false;
#endif
}

qacc bool is_little_endian() { return not is_big_endian(); }

// -------------------

}  // namespace qlat
