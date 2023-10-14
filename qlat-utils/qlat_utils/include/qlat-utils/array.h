#pragma once

#include <cassert>

#include <qlat-utils/show.h>
#include <qlat-utils/qacc.h>

namespace qlat
{  //

template <class M, unsigned long N>
struct API array {
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

template <class M>
struct API array<M, 0> {
  qacc unsigned long size() const { return 0; }
  //
  qacc void fill(const M& x) { (void)x; }
  qacc M* data() { return NULL; }
  //
  qacc const M* data() const { return NULL; }
  //
  qacc M& operator[](unsigned long k)
  {
    (void)k;
    assert(false);
    static M x;
    return x;
  };
  qacc const M& operator[](unsigned long k) const
  {
    (void)k;
    assert(false);
    static M x;
    return x;
  };
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
