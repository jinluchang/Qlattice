#pragma once

#include <qlat-utils/qacc.h>
#include <qlat-utils/array.h>
#include <qlat-utils/assert.h>

namespace qlat
{  //

template <class M>
struct API Handle {
  M* p;
  //
  qacc Handle<M>() { init(); }
  qacc Handle<M>(M& obj) { init(obj); }
  //
  qacc void init() { p = NULL; }
  qacc void init(M& obj) { p = (M*)&obj; }
  //
  qacc bool null() const { return p == NULL; }
  //
  qacc M& operator()() const
  {
    qassert(NULL != p);
    return *p;
  }
};

template <class M>
struct API ConstHandle {
  const M* p;
  //
  qacc ConstHandle<M>() { init(); }
  qacc ConstHandle<M>(const M& obj) { init(obj); }
  qacc ConstHandle<M>(const Handle<M>& h) { init(h()); }
  //
  qacc void init() { p = NULL; }
  qacc void init(const M& obj) { p = (M*)&obj; }
  //
  qacc bool null() const { return p == NULL; }
  //
  qacc const M& operator()() const
  {
    qassert(NULL != p);
    return *p;
  }
};

template <class M>
struct API Vector {
  M* p;
  long n;
  //
  qacc Vector<M>()
  {
    p = NULL;
    n = 0;
  }
  qacc Vector<M>(const Vector<M>& vec)
  {
    p = vec.p;
    n = vec.n;
  }
  template <int N>
  qacc Vector<M>(const array<M, N>& arr)
  {
    p = (M*)arr.data();
    n = arr.size();
  }
  Vector<M>(const std::vector<M>& vec)
  {
    p = (M*)vec.data();
    n = vec.size();
  }
  qacc Vector<M>(const M* p_, const long n_)
  {
    p = (M*)p_;
    n = n_;
  }
  qacc Vector<M>(const M& x)
  {
    p = (M*)&x;
    n = 1;
  }
  //
  qacc const M& operator[](const long i) const
  {
    if (not(0 <= i && i < n)) {
#ifndef QLAT_IN_ACC
      displayln(
          ssprintf("ERROR: expect: 0 <= i && i < n but: i=%d n=%d sizeof(M)=%d",
                   i, n, sizeof(M)));
#endif
      qassert(false);
    }
    return p[i];
  }
  qacc M& operator[](const long i)
  {
    if (not(0 <= i && i < n)) {
#ifndef QLAT_IN_ACC
      displayln(
          ssprintf("ERROR: expect: 0 <= i && i < n but: i=%d n=%d sizeof(M)=%d",
                   i, n, sizeof(M)));
#endif
      qassert(false);
    }
    return p[i];
  }
  //
  qacc M* data() { return p; }
  qacc const M* data() const { return p; }
  //
  qacc long size() const { return n; }
  //
  qacc long data_size() const { return n * sizeof(M); }
  //
  qacc const Vector<M>& operator=(const Vector<M>& v)
  {
    n = v.n;
    p = v.p;
    return *this;
  }
};

template <class M, int N>
struct API Array {
  M* p;
  //
  qacc Array<M, N>() { p = NULL; }
  qacc Array<M, N>(const Array<M, N>& arr) { p = arr.p; }
  qacc Array<M, N>(const Vector<M>& vec)
  {
    qassert(N == vec.size());
    p = vec.p;
  }
  qacc Array<M, N>(const array<M, N>& arr) { p = (M*)arr.data(); }
  qacc Array<M, N>(const M* p_) { p = (M*)p_; }
  qacc Array<M, N>(const M& x)
  {
    qassert(N == 1);
    p = (M*)&x;
  }
  //
  qacc const M& operator[](int i) const
  {
    qassert(0 <= i && i < N);
    return p[i];
  }
  qacc M& operator[](int i)
  {
    qassert(0 <= i && i < N);
    return p[i];
  }
  //
  qacc M* data() { return p; }
  qacc const M* data() const { return p; }
  //
  qacc int size() const { return N; }
  //
  qacc long data_size() const { return N * sizeof(M); }
  //
  qacc const Array<M, N>& operator=(const Array<M, N>& v)
  {
    p = v.p;
    return *this;
  }
  qacc const Array<M, N>& operator=(const Vector<M>& v)
  {
    qassert(N == v.size());
    p = v.p;
    return *this;
  }
};

}  // namespace qlat
