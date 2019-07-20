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

inline void from_big_endian_32(char* str, const size_t len)
// obsolete
{
  to_from_big_endian_32(str, len);
}

inline void from_big_endian_64(char* str, const size_t len)
// obsolete
{
  to_from_big_endian_64(str, len);
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
    std::memcpy((void*)&x, (void*)&y, sizeof(M));
  } else {
    // if M has a larger size, than leave the extra space untouched
    std::memcpy((void*)&x, (void*)&y, sizeof(N));
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
