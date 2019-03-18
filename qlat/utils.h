// vim: set ts=2 sw=2 expandtab:

#pragma once

#include <qlat/config.h>

#include <endian.h>

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

namespace qlat
{  //

inline void set_zero(double& x) { x = 0; }

inline void set_zero(Complex& x) { x = 0; }

template <class M, unsigned long N>
void set_zero(std::array<M, N>& arr)
{
  long size = N * sizeof(M);
  std::memset(arr.data(), 0, size);
}

inline void set_unit(double& x, const double& coef = 1.0) { x = coef; }

inline void set_unit(Complex& x, const Complex& coef = 1.0) { x = coef; }

template <class M>
void set_zero(std::vector<M>& vec)
{
  long size = vec.size() * sizeof(M);
  std::memset(vec.data(), 0, size);
}

inline double qnorm(const double& x) { return x * x; }

inline double qnorm(const Complex& x) { return std::norm(x); }

template <class T, size_t N>
inline double qnorm(const std::array<T, N>& mm)
{
  double sum = 0.0;
  for (size_t i = 0; i < N; ++i) {
    sum += qnorm(mm[i]);
  }
  return sum;
}

template <class Vec>
bool is_equal_vec(const Vec& v1, const Vec& v2)
{
  const bool b = v1.size() == v2.size();
  if (not b) {
    return false;
  } else {
    const long s = v1.size();
    for (long i = 0; i < s; ++i) {
      if (not(v1[i] == v2[i])) {
        return false;
      }
    }
    return true;
  }
}

template <class M>
bool operator==(const std::vector<M>& v1, const std::vector<M>& v2)
{
  return is_equal_vec(v1, v2);
}

template <class M, int N>
bool operator==(const std::array<M, N>& v1, const std::array<M, N>& v2)
{
  return is_equal_vec(v1, v2);
}

template <class M>
void set_zero(Vector<M> vec)
{
  std::memset(vec.data(), 0, vec.data_size());
}

template <class M, int N>
struct Array {
  M* p;
  //
  Array<M, N>() { p = NULL; }
  Array<M, N>(const Array<M, N>& arr) { p = arr.p; }
  Array<M, N>(const Vector<M>& vec)
  {
    qassert(N == vec.size());
    p = vec.p;
  }
  Array<M, N>(const std::array<M, N>& arr) { p = (M*)arr.data(); }
  Array<M, N>(const M* p_) { p = (M*)p_; }
  Array<M, N>(const M& x)
  {
    qassert(N == 1);
    p = (M*)&x;
  }
  //
  const M& operator[](int i) const
  {
    qassert(0 <= i && i < N);
    return p[i];
  }
  M& operator[](int i)
  {
    qassert(0 <= i && i < N);
    return p[i];
  }
  //
  M* data() { return p; }
  const M* data() const { return p; }
  //
  int size() const { return N; }
  //
  long data_size() const { return N * sizeof(M); }
  //
  const Array<M, N>& operator=(const Array<M, N>& v)
  {
    p = v.p;
    return *this;
  }
  const Array<M, N>& operator=(const Vector<M>& v)
  {
    qassert(N == v.size());
    p = v.p;
    return *this;
  }
};

template <class M, int N>
void set_zero(Array<M, N> arr)
{
  long size = N * sizeof(M);
  std::memset(arr.data(), 0, size);
}

template <class M, int N>
Vector<M> get_data(Array<M, N> arr)
{
  return Vector<M>(arr.data(), arr.size());
}

template <class M>
Vector<M> get_data(Vector<M> vec)
{
  return vec;
}

template <class M, unsigned long N>
Vector<M> get_data(const std::array<M, N>& vec)
{
  return Vector<M>((M*)vec.data(), vec.size());
}

template <class M>
Vector<M> get_data(const std::vector<M>& vec)
{
  return Vector<M>((M*)vec.data(), vec.size());
}

inline Vector<char> get_data(const std::string& str)
{
  return Vector<char>(&str[0], str.length());
}

template <class M>
Vector<M> get_data(const Handle<M>& h)
{
  return Vector<M>(h.p, 1);
}

template <class M>
Vector<M> get_data(const ConstHandle<M>& h)
{
  return Vector<M>(h.p, 1);
}

inline Vector<long> get_data(const long& x) { return Vector<long>(&x, 1); }

inline Vector<double> get_data(const double& x)
{
  return Vector<double>(&x, 1);
}

inline Vector<int> get_data(const int& x) { return Vector<int>(&x, 1); }

inline Vector<float> get_data(const float& x) { return Vector<float>(&x, 1); }

template <class M>
Vector<double> get_data_double(const M& v)
{
  return Vector<double>(&v, sizeof(M) / sizeof(double));
}

template <class M>
Vector<long> get_data_long(const M& v)
{
  return Vector<long>(&v, sizeof(M) / sizeof(long));
}

template <class M>
long get_data_size(const M& x)
{
  return get_data(x).data_size();
}

template <class M, int N>
void assign(std::array<M, N>& vec, const Array<M, N>& src)
{
  memcpy(vec.data(), src.data(), src.data_size());
}

template <class M, int N>
void assign(std::array<M, N>& vec, const Vector<M>& src)
{
  qassert(N == src.size());
  memcpy(vec.data(), src.data(), src.data_size());
}

template <class M, int N>
void assign(std::vector<M>& vec, const Array<M, N>& src)
{
  vec.resize(src.size());
  memcpy(vec.data(), src.data(), src.data_size());
}

template <class M>
void assign(std::vector<M>& vec, const Vector<M>& src)
{
  vec.resize(src.size());
  memcpy(vec.data(), src.data(), src.data_size());
}

template <class M>
void assign(Vector<M> vec, const Vector<M>& src)
{
  qassert(vec.size() == src.size());
  memcpy(vec.data(), src.data(), src.data_size());
}

template <class M, int N>
void assign(Vector<M> vec, const Array<M, N>& src)
{
  qassert(vec.size() == N);
  memcpy(vec.data(), src.data(), src.data_size());
}

template <class M, int N>
void assign(Array<M, N> vec, const Array<M, N>& src)
{
  memcpy(vec.data(), src.data(), src.data_size());
}

template <class M, int N>
void assign(Array<M, N> vec, const Vector<M>& src)
{
  qassert(src.size() == N);
  memcpy(vec.data(), src.data(), src.data_size());
}

template <class M, class N>
void assign_truncate(M& x, const N& y)
{
  if (sizeof(M) <= sizeof(N)) {
    memcpy(&x, &y, sizeof(M));
  } else {
    // if M has a larger size, than leave the extra space untouched
    memcpy(&x, &y, sizeof(N));
  }
}

inline bool is_integer(const double& x)
{
  const double diff = x - (long)x;
  return 1e-6 > diff || diff > 1 - 1e-6;
}

template <class M>
inline bool is_integer(const std::vector<M>& v)
{
  for (int i = 0; i < (int)v.size(); ++i) {
    if (!is_integer(v[i])) {
      return false;
    }
  }
  return true;
}

template <class M, unsigned long N>
inline bool is_integer(const std::array<M, N>& v)
{
  for (int i = 0; i < N; ++i) {
    if (!is_integer(v[i])) {
      return false;
    }
  }
  return true;
}

template <class M>
inline void random_permute(std::vector<M>& vec, RngState& rs)
{
  const long size = (long)vec.size();
  M tmp;
  for (long k = 0; k < size; ++k) {
    const long kk = rand_gen(rs) % (size - k);
    tmp = vec[k];
    vec[k] = vec[k + kk];
    vec[k + kk] = tmp;
  }
}

inline std::string show(const Complex& x)
{
  return ssprintf("(%24.17E + %24.17E j)", x.real(), x.imag());
}

}  // namespace qlat
