#pragma once

#include <qlat-utils/qacc.h>
#include <qlat-utils/show.h>

#include <cassert>

namespace qlat
{  //

template <class M, std::size_t N>
struct API array {
  M v[N];
  //
  qacc std::size_t size() const { return N; }
  //
  qacc void fill(const M& x)
  {
    for (std::size_t i = 0; i < N; ++i) {
      v[i] = x;
    }
  }
  //
  qacc M* data() { return v; }
  qacc const M* data() const { return v; }
  //
  qacc M& operator[](std::size_t k)
  {
    assert(k < N);
    return v[k];
  };
  qacc const M& operator[](std::size_t k) const
  {
    assert(k < N);
    return v[k];
  };
};

template <class M>
struct API array<M, 0> {
  qacc std::size_t size() const { return 0; }
  //
  qacc void fill(const M& x) { (void)x; }
  qacc M* data() { return NULL; }
  //
  qacc const M* data() const { return NULL; }
  //
  qacc M& operator[](std::size_t k)
  {
    (void)k;
    assert(false);
    static M x;
    return x;
  };
  qacc const M& operator[](std::size_t k) const
  {
    (void)k;
    assert(false);
    static M x;
    return x;
  };
};

template <class M, std::size_t N>
qacc bool operator<(const array<M, N>& a1, const array<M, N>& a2)
{
  for (std::size_t i = 0; i < N; ++i) {
    if (a1[i] < a2[i]) {
      return true;
    } else if (a1[i] > a2[i]) {
      return false;
    }
  }
  return false;
}

template <class M, std::size_t N>
qacc bool operator<=(const array<M, N>& a1, const array<M, N>& a2)
{
  for (std::size_t i = 0; i < N; ++i) {
    if (a1[i] < a2[i]) {
      return true;
    } else if (a1[i] > a2[i]) {
      return false;
    }
  }
  return true;
}

}  // namespace qlat
