#pragma once

#include "qutils.h"

namespace qlat
{  //

template <class M>
std::array<M, 0> make_array()
{
  std::array<M, 0> arr;
  return arr;
}

template <class M>
std::array<M, 1> make_array(const M& x)
{
  std::array<M, 1> arr;
  arr[0] = x;
  return arr;
}

template <class M>
std::array<M, 2> make_array(const M& x, const M& x1)
{
  std::array<M, 2> arr;
  arr[0] = x;
  arr[1] = x1;
  return arr;
}

template <class M>
std::array<M, 3> make_array(const M& x, const M& x1, const M& x2)
{
  std::array<M, 3> arr;
  arr[0] = x;
  arr[1] = x1;
  arr[2] = x2;
  return arr;
}

template <class M>
std::array<M, 4> make_array(const M& x, const M& x1, const M& x2, const M& x3)
{
  std::array<M, 4> arr;
  arr[0] = x;
  arr[1] = x1;
  arr[2] = x2;
  arr[3] = x3;
  return arr;
}

template <class M>
std::array<M, 5> make_array(const M& x, const M& x1, const M& x2, const M& x3,
                            const M& x4)
{
  std::array<M, 5> arr;
  arr[0] = x;
  arr[1] = x1;
  arr[2] = x2;
  arr[3] = x3;
  arr[4] = x4;
  return arr;
}

template <class M>
std::array<M, 6> make_array(const M& x, const M& x1, const M& x2, const M& x3,
                            const M& x4, const M& x5)
{
  std::array<M, 6> arr;
  arr[0] = x;
  arr[1] = x1;
  arr[2] = x2;
  arr[3] = x3;
  arr[4] = x4;
  arr[5] = x5;
  return arr;
}

template <class M>
std::array<M, 7> make_array(const M& x, const M& x1, const M& x2, const M& x3,
                            const M& x4, const M& x5, const M& x6)
{
  std::array<M, 7> arr;
  arr[0] = x;
  arr[1] = x1;
  arr[2] = x2;
  arr[3] = x3;
  arr[4] = x4;
  arr[5] = x5;
  arr[6] = x6;
  return arr;
}

template <class M>
std::array<M, 8> make_array(const M& x, const M& x1, const M& x2, const M& x3,
                            const M& x4, const M& x5, const M& x6, const M& x7)
{
  std::array<M, 8> arr;
  arr[0] = x;
  arr[1] = x1;
  arr[2] = x2;
  arr[3] = x3;
  arr[4] = x4;
  arr[5] = x5;
  arr[6] = x6;
  arr[7] = x7;
  return arr;
}

template <class M>
std::array<M, 9> make_array(const M& x, const M& x1, const M& x2, const M& x3,
                            const M& x4, const M& x5, const M& x6, const M& x7,
                            const M& x8)
{
  std::array<M, 9> arr;
  arr[0] = x;
  arr[1] = x1;
  arr[2] = x2;
  arr[3] = x3;
  arr[4] = x4;
  arr[5] = x5;
  arr[6] = x6;
  arr[7] = x7;
  arr[8] = x8;
  return arr;
}

template <class M>
std::array<M, 10> make_array(const M& x, const M& x1, const M& x2, const M& x3,
                             const M& x4, const M& x5, const M& x6, const M& x7,
                             const M& x8, const M& x9)
{
  std::array<M, 10> arr;
  arr[0] = x;
  arr[1] = x1;
  arr[2] = x2;
  arr[3] = x3;
  arr[4] = x4;
  arr[5] = x5;
  arr[6] = x6;
  arr[7] = x7;
  arr[8] = x8;
  arr[9] = x9;
  return arr;
}

template <class M>
struct Handle {
  M* p;
  //
  Handle<M>() { init(); }
  Handle<M>(M& obj) { init(obj); }
  //
  void init() { p = NULL; }
  void init(M& obj) { p = (M*)&obj; }
  //
  bool null() const { return p == NULL; }
  //
  M& operator()() const
  {
    qassert(NULL != p);
    return *p;
  }
};

template <class M>
struct ConstHandle {
  const M* p;
  //
  ConstHandle<M>() { init(); }
  ConstHandle<M>(const M& obj) { init(obj); }
  ConstHandle<M>(const Handle<M>& h) { init(h()); }
  //
  void init() { p = NULL; }
  void init(const M& obj) { p = (M*)&obj; }
  //
  bool null() const { return p == NULL; }
  //
  const M& operator()() const
  {
    qassert(NULL != p);
    return *p;
  }
};

template <class M>
struct Vector {
  M* p;
  long n;
  //
  Vector<M>()
  {
    p = NULL;
    n = 0;
  }
  Vector<M>(const Vector<M>& vec)
  {
    p = vec.p;
    n = vec.n;
  }
  template <int N>
  Vector<M>(const std::array<M, N>& arr)
  {
    p = (M*)arr.data();
    n = arr.size();
  }
  Vector<M>(const std::vector<M>& vec)
  {
    p = (M*)vec.data();
    n = vec.size();
  }
  Vector<M>(const M* p_, const long n_)
  {
    p = (M*)p_;
    n = n_;
  }
  Vector<M>(const M& x)
  {
    p = (M*)&x;
    n = 1;
  }
  //
  const M& operator[](long i) const
  {
    qassert(0 <= i && i < n);
    return p[i];
  }
  M& operator[](long i)
  {
    qassert(0 <= i && i < n);
    return p[i];
  }
  //
  M* data() { return p; }
  const M* data() const { return p; }
  //
  long size() const { return n; }
  //
  long data_size() const { return n * sizeof(M); }
  //
  const Vector<M>& operator=(const Vector<M>& v)
  {
    n = v.n;
    p = v.p;
    return *this;
  }
};

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
void to_from_little_endian_16(Vector<M> v)
{
  to_from_little_endian_16((void*)v.data(), v.data_size());
}

template <class M>
void to_from_little_endian_32(Vector<M> v)
{
  to_from_little_endian_32((void*)v.data(), v.data_size());
}

template <class M>
void to_from_little_endian_64(Vector<M> v)
{
  to_from_little_endian_64((void*)v.data(), v.data_size());
}

template <class M>
void to_from_big_endian_16(Vector<M> v)
{
  to_from_big_endian_16((void*)v.data(), v.data_size());
}

template <class M>
void to_from_big_endian_32(Vector<M> v)
{
  to_from_big_endian_32((void*)v.data(), v.data_size());
}

template <class M>
void to_from_big_endian_64(Vector<M> v)
{
  to_from_big_endian_64((void*)v.data(), v.data_size());
}

}  // namespace qlat
