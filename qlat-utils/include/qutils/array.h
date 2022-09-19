#pragma once

#include <cassert>

#include <qutils/qacc.h>

namespace qlat
{  //

template <class M, unsigned long N>
struct array {
  M v[N];
  //
  qacc unsigned long size() const { return N; }
  //
  qacc void fill(const M& x)
  {
    for (unsigned long i = 0; i < N; ++i) {
      v[i] = x;
    }
  }
  //
  qacc M* data() { return v; }
  qacc const M* data() const { return v; }
  //
  qacc M& operator[](unsigned long k) { return v[k]; };
  qacc const M& operator[](unsigned long k) const { return v[k]; };
};

template <class M, unsigned long N>
qacc bool operator<(const array<M, N>& a1, const array<M, N>& a2)
{
  for (unsigned long i = 0; i < N; ++i) {
    if (a1[i] < a2[i]) {
      return true;
    } else if (a1[i] > a2[i]) {
      return false;
    }
  }
  return false;
}

template <class M, unsigned long N>
qacc bool operator<=(const array<M, N>& a1, const array<M, N>& a2)
{
  for (unsigned long i = 0; i < N; ++i) {
    if (a1[i] < a2[i]) {
      return true;
    } else if (a1[i] > a2[i]) {
      return false;
    }
  }
  return true;
}

}  // namespace qlat
