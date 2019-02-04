#pragma once

#include <cassert>

namespace std
{
template <class M, unsigned long N>
struct array {
  M v[N];
  //
  unsigned long size() const { return N; }
  //
  void fill(const M& x)
  {
    for (int i = 0; i < N; ++i) {
      v[i] = x;
    }
  }
  //
  M* data() { return v; }
  const M* data() const { return v; }
  //
  M& operator[](int k) { return v[k]; };
  const M& operator[](int k) const { return v[k]; };
};

template <class M, unsigned long N>
inline bool operator<(const array<M, N>& a1, const array<M, N>& a2)
{
  for (int i = 0; i < N; ++i) {
    if (a1[i] < a2[i]) {
      return true;
    } else if (a1[i] > a2[i]) {
      return false;
    }
  }
  return false;
}

template <class M, unsigned long N>
inline bool operator<=(const array<M, N>& a1, const array<M, N>& a2)
{
  for (int i = 0; i < N; ++i) {
    if (a1[i] < a2[i]) {
      return true;
    } else if (a1[i] > a2[i]) {
      return false;
    }
  }
  return true;
}

}  // namespace std
