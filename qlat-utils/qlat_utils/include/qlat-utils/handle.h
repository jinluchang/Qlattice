#pragma once

#include <qlat-utils/array.h>
#include <qlat-utils/assert.h>
#include <qlat-utils/qacc.h>

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
  //
  qacc M& val() const
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
  //
  qacc const M& val() const
  {
    qassert(NULL != p);
    return *p;
  }
};

template <class M>
struct API Vector {
  M* p;
  Long n;
  //
  qacc Vector<M>()
  {
    p = NULL;
    n = 0;
  }
  qacc Vector<M>(const Vector<M>& vec)
  {
    if (vec.p == NULL) {
      qassert(vec.n == 0);
    }
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
  qacc Vector<M>(const M* p_, const Long n_)
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
  template <class N>
  qacc void set_cast(const Vector<N>& vec)
  {
    if (vec.p == NULL) {
      qassert(vec.n == 0);
      p = NULL;
      n = 0;
    } else {
      Long total_size = vec.n * sizeof(N);
      p = (M*)vec.p;
      n = total_size / (Long)sizeof(M);
      if (n * (Long)sizeof(M) != total_size) {
#ifndef QLAT_IN_ACC
        qerr(
            ssprintf("Vector::set_cast: n=%ld ; sizeof(M)=%ld ; "
                     "total_size=%ld=%ld*%ld .",
                     (long)n, (long)sizeof(M), (long)total_size, (long)vec.n,
                     (long)sizeof(N)));
#else
        qassert(false);
#endif
      }
    }
  }
  //
  qacc const M& operator[](const Long i) const
  {
    if (not(0 <= i && i < n)) {
#ifndef QLAT_IN_ACC
      qerr(
          ssprintf("ERROR: expect: 0 <= i && i < n but: i=%d n=%d sizeof(M)=%d",
                   i, n, sizeof(M)));
#else
      printf("%ld %ld %ld", (long)i, (long)n, (long)sizeof(M));
      qassert(false);
#endif
    }
    return p[i];
  }
  qacc M& operator[](const Long i)
  {
    if (not(0 <= i && i < n)) {
#ifndef QLAT_IN_ACC
      qerr(
          ssprintf("ERROR: expect: 0 <= i && i < n but: i=%d n=%d sizeof(M)=%d",
                   i, n, sizeof(M)));
#else
      printf("%ld %ld %ld", (long)i, (long)n, (long)sizeof(M));
      qassert(false);
#endif
    }
    return p[i];
  }
  //
  qacc M* data() { return p; }
  qacc const M* data() const { return p; }
  //
  qacc Long size() const { return n; }
  //
  qacc Long data_size() const { return n * sizeof(M); }
  //
  qacc Vector<M>& operator=(const Vector<M>& vec)
  {
    if (vec.p == NULL) {
      qassert(vec.n == 0);
    }
    n = vec.n;
    p = vec.p;
    return *this;
  }
};

template <class M, Long N>
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
  qacc const M& operator[](Long i) const
  {
    qassert(0 <= i && i < N);
    return p[i];
  }
  qacc M& operator[](Long i)
  {
    qassert(0 <= i && i < N);
    return p[i];
  }
  //
  qacc M* data() { return p; }
  qacc const M* data() const { return p; }
  //
  qacc Long size() const { return N; }
  //
  qacc Long data_size() const { return N * sizeof(M); }
  //
  qacc Array<M, N>& operator=(const Array<M, N>& v)
  {
    p = v.p;
    return *this;
  }
  qacc Array<M, N>& operator=(const Vector<M>& v)
  {
    qassert(N == v.size());
    p = v.p;
    return *this;
  }
};

}  // namespace qlat
