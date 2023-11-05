#pragma once

#include <qlat-utils/core.h>

#include <vector>
#include <map>
#include <set>

namespace qlat
{  //

// -------------------

qacc void set_zero(char& x) { x = 0; }

qacc void set_zero(int8_t& x) { x = 0; }

qacc void set_zero(int& x) { x = 0; }

qacc void set_zero(Long& x) { x = 0; }

qacc void set_zero(double& x) { x = 0; }

qacc void set_zero(float& x) { x = 0; }

qacc void set_zero(ComplexD& x) { x = 0; }

qacc void set_zero(ComplexF& x) { x = 0; }

template <class M, unsigned long N>
qacc void set_zero(array<M, N>& arr)
{
  Long size = N * sizeof(M);
  std::memset((void*)arr.data(), 0, size);
}

template <class M>
void set_zero(std::vector<M>& vec)
{
  Long size = vec.size() * sizeof(M);
  std::memset((void*)vec.data(), 0, size);
}

// -------------------

qacc void set_unit(char& x, const Long& coef = 1) { x = coef; }

qacc void set_unit(char& x, const ComplexD& coef = 1) { x = coef.real(); }

qacc void set_unit(int8_t& x, const Long& coef = 1) { x = coef; }

qacc void set_unit(int8_t& x, const ComplexD& coef = 1) { x = coef.real(); }

qacc void set_unit(int& x, const Long& coef = 1) { x = coef; }

qacc void set_unit(int& x, const ComplexD& coef) { x = coef.real(); }

qacc void set_unit(Long& x, const Long& coef = 1) { x = coef; }

qacc void set_unit(Long& x, const ComplexD& coef) { x = coef.real(); }

qacc void set_unit(float& x, const double& coef = 1.0) { x = coef; }

qacc void set_unit(float& x, const ComplexD& coef) { x = coef.real(); }

qacc void set_unit(double& x, const double& coef = 1.0) { x = coef; }

qacc void set_unit(double& x, const ComplexD& coef) { x = coef.real(); }

qacc void set_unit(ComplexD& x, const ComplexD& coef = 1.0) { x = coef; }

qacc void set_unit(ComplexF& x, const ComplexD& coef = 1.0) { x = coef; }

// -------------------

qacc double qnorm(const char x) { return x * x; }

qacc double qnorm(const int8_t x) { return x * x; }

qacc double qnorm(const int x) { return x * x; }

qacc double qnorm(const Long x) { return x * x; }

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

}  // namespace qlat
