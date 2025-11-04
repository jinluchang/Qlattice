#pragma once

#include <qlat-utils/array.h>
#include <qlat-utils/qassert.h>
#include <qlat-utils/qacc.h>
#include <qlat-utils/timer.h>

namespace qlat
{  //

template <class M>
struct API Handle {
  mutable M* p;
  //
  qacc Handle() { init(); }
  qacc Handle(M& obj) { init(obj); }
  //
  qacc void init() { p = NULL; }
  qacc void init(M& obj) { p = &obj; }
  qacc void init(Handle<M>& h) { p = h.p; }
  //
  qacc bool null() const { return p == NULL; }
  //
  qacc const M& operator()() const
  {
    qassert(NULL != p);
    return *p;
  }
  qacc M& operator()()
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
  qacc M& val()
  {
    qassert(NULL != p);
    return *p;
  }
  //
  qacc M& cast_const() const
  {
    qassert(NULL != p);
    return *p;
  }
};

template <class M>
struct API ConstHandle {
  mutable const M* p;
  //
  qacc ConstHandle() { init(); }
  qacc ConstHandle(const M& obj) { init(obj); }
  qacc ConstHandle(const Handle<M>& h) { init(h()); }
  //
  qacc void init() { p = NULL; }
  qacc void init(const M& obj) { p = &obj; }
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
  mutable M* p;
  Long n;
  //
  qacc Vector()
  {
    p = NULL;
    n = 0;
  }
  qacc Vector(const Vector<M>& vec)
  {
    if (vec.p == NULL) {
      qassert(vec.n == 0);
    }
    p = vec.p;
    n = vec.n;
  }
  template <std::size_t N>
  qacc Vector(const array<M, N>& arr)
  {
    p = (M*)arr.data();
    n = arr.size();
  }
  Vector(const std::vector<M>& vec)
  {
    p = (M*)vec.data();
    n = vec.size();
  }
  qacc Vector(const M* p_, const Long n_)
  {
    p = (M*)p_;
    n = n_;
  }
  qacc Vector(const M& x)
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
      qassert(n * (Long)sizeof(M) == total_size);
    }
  }
  //
  qacc const M& operator[](const Long i) const
  {
    qassert(0 <= i && i < n);
    return p[i];
  }
  qacc M& operator[](const Long i)
  {
    qassert(0 <= i && i < n);
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

template <class M, std::size_t N>
struct API Array {
  mutable M* p;
  //
  qacc Array() { p = NULL; }
  qacc Array(const Array<M, N>& arr) { p = arr.p; }
  qacc Array(const Vector<M>& vec)
  {
    qassert(N == vec.size());
    p = vec.p;
  }
  qacc Array(const array<M, N>& arr) { p = (M*)arr.data(); }
  qacc Array(const M* p_) { p = (M*)p_; }
  qacc Array(const M& x)
  {
    qassert(N == 1);
    p = (M*)&x;
  }
  //
  qacc const M& operator[](Long i) const
  {
    qassert(0 <= i && i < (Long)N);
    return p[i];
  }
  qacc M& operator[](Long i)
  {
    qassert(0 <= i && i < (Long)N);
    return p[i];
  }
  //
  qacc M* data() { return p; }
  qacc const M* data() const { return p; }
  //
  qacc Long size() const { return (Long)N; }
  //
  qacc Long data_size() const { return (Long)N * sizeof(M); }
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
