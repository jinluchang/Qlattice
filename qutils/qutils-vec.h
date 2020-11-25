#pragma once

#include <qutils/qutils.h>
#include <qutils/rng-state.h>

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
    if (not(0 <= i && i < n)) {
      displayln(
          ssprintf("ERROR: expect: 0 <= i && i < n but: i=%d n=%d sizeof(M)=%d",
                   i, n, sizeof(M)));
      qassert(false);
    }
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
  std::memset((void*)vec.data(), 0, vec.data_size());
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
  std::memset((void*)arr.data(), 0, size);
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

template <class M>
Vector<M> get_data_one_elem(const M& x)
{
  return Vector<M>(&x, 1);
}

inline Vector<long> get_data(const long& x) { return get_data_one_elem(x); }

inline Vector<double> get_data(const double& x) { return get_data_one_elem(x); }

inline Vector<int> get_data(const int& x) { return get_data_one_elem(x); }

inline Vector<float> get_data(const float& x) { return get_data_one_elem(x); }

template <class T>
double qnorm(const Vector<T>& mm)
{
  double sum = 0.0;
  const long size = mm.size();
  for (long i = 0; i < size; ++i) {
    sum += qnorm(mm[i]);
  }
  return sum;
}

template <class T>
double qnorm(const Vector<T>& m1, const Vector<T>& m2)
{
  double sum = 0.0;
  const long size = m1.size();
  qassert(size == (long)m2.size());
  for (long i = 0; i < size; ++i) {
    sum += qnorm(m1[i], m2[i]);
  }
  return sum;
}

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

template <class T>
double qnorm_double(const Vector<T>& m1, const Vector<T>& m2)
{
  const Vector<double> dm1((double*)m1.data(), m1.data_size() / sizeof(double));
  const Vector<double> dm2((double*)m2.data(), m2.data_size() / sizeof(double));
  return qnorm(dm1, dm2);
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
  std::memcpy((void*)vec.data(), (void*)src.data(), src.data_size());
}

template <class M, int N>
void assign(std::array<M, N>& vec, const Vector<M>& src)
{
  qassert(N == src.size());
  std::memcpy((void*)vec.data(), (void*)src.data(), src.data_size());
}

template <class M, int N>
void assign(std::vector<M>& vec, const Array<M, N>& src)
{
  vec.resize(src.size());
  std::memcpy((void*)vec.data(), (void*)src.data(), src.data_size());
}

template <class M>
void assign(std::vector<M>& vec, const Vector<M>& src)
{
  vec.resize(src.size());
  std::memcpy((void*)vec.data(), (void*)src.data(), src.data_size());
}

template <class M>
void assign(Vector<M> vec, const Vector<M>& src)
{
  qassert(vec.size() == src.size());
  std::memcpy((void*)vec.data(), (void*)src.data(), src.data_size());
}

template <class M, class N>
void assign(Vector<M> vec, const Vector<N>& src)
{
  qassert(vec.data_size() == src.data_size());
  std::memcpy((void*)vec.data(), (void*)src.data(), src.data_size());
}

template <class M, int N>
void assign(Vector<M> vec, const Array<M, N>& src)
{
  qassert(vec.size() == N);
  std::memcpy((void*)vec.data(), (void*)src.data(), src.data_size());
}

template <class M, int N>
void assign(Array<M, N> vec, const Array<M, N>& src)
{
  std::memcpy((void*)vec.data(), (void*)src.data(), src.data_size());
}

template <class M, int N>
void assign(Array<M, N> vec, const Vector<M>& src)
{
  qassert(src.size() == N);
  std::memcpy((void*)vec.data(), (void*)src.data(), src.data_size());
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

template <class M, int N>
Array<M, N> operator+=(Array<M, N> v, const Array<M, N> v1)
{
  for (int i = 0; i < N; ++i) {
    v.p[i] += v1.p[i];
  }
  return v;
}

template <class M, int N>
Array<M, N> operator-=(Array<M, N> v, const Array<M, N> v1)
{
  for (int i = 0; i < N; ++i) {
    v.p[i] -= v1.p[i];
  }
  return v;
}

template <class M, int N>
Array<M, N> operator*=(Array<M, N> v, const double factor)
{
  for (int i = 0; i < N; ++i) {
    v.p[i] *= factor;
  }
  return v;
}

template <class M, int N>
Array<M, N> operator*=(Array<M, N> v, const Complex factor)
{
  for (int i = 0; i < N; ++i) {
    v.p[i] *= factor;
  }
  return v;
}

template <class M>
Vector<M> operator+=(Vector<M> v, const Vector<M> v1)
{
  qassert(v.size() == v1.size());
  for (long i = 0; i < v.size(); ++i) {
    v.p[i] += v1.p[i];
  }
  return v;
}

template <class M>
Vector<M> operator-=(Vector<M> v, const Vector<M> v1)
{
  qassert(v.size() == v1.size());
  for (long i = 0; i < v.size(); ++i) {
    v.p[i] -= v1.p[i];
  }
  return v;
}

template <class M>
Vector<M> operator*=(Vector<M> v, const double factor)
{
  for (long i = 0; i < v.size(); ++i) {
    v.p[i] *= factor;
  }
  return v;
}

template <class M>
Vector<M> operator*=(Vector<M> v, const Complex factor)
{
  for (long i = 0; i < v.size(); ++i) {
    v.p[i] *= factor;
  }
  return v;
}

template <class M, unsigned long N>
std::array<M, N> operator+(const std::array<M, N>& v1,
                           const std::array<M, N>& v2)
{
  std::array<M, N> ret;
  for (unsigned long i = 0; i < N; ++i) {
    ret[i] = v1[i] + v2[i];
  }
  return ret;
}

template <class M, unsigned long N>
std::array<M, N> operator-(const std::array<M, N>& v1,
                           const std::array<M, N>& v2)
{
  std::array<M, N> ret;
  for (unsigned long i = 0; i < N; ++i) {
    ret[i] = v1[i] - v2[i];
  }
  return ret;
}

template <class M>
inline void random_permute(std::vector<M>& vec, const RngState& rs_)
{
  RngState rs = rs_;
  const long size = (long)vec.size();
  M tmp;
  for (long k = 0; k < size; ++k) {
    const long kk = rand_gen(rs) % (size - k);
    tmp = vec[k];
    vec[k] = vec[k + kk];
    vec[k + kk] = tmp;
  }
}

template <class M>
void set_u_rand_double(Vector<M> v, const RngState& rs, const double upper = 1.0,
                      const double lower = -1.0)
{
  RngState rsi = rs;
  Vector<double> dv((double*)v.data(), v.data_size() / sizeof(double));
  for (int m = 0; m < dv.size(); ++m) {
    dv[m] = u_rand_gen(rsi, 1.0, -1.0);
  }
}

template <class M>
void set_u_rand_float(Vector<M> v, const RngState& rs, const double upper = 1.0,
                      const double lower = -1.0)
{
  RngState rsi = rs;
  Vector<float> dv((float*)v.data(), v.data_size() / sizeof(float));
  for (int m = 0; m < dv.size(); ++m) {
    dv[m] = u_rand_gen(rsi, 1.0, -1.0);
  }
}

template <class M>
bool does_vec_have(const std::vector<M>& vec, const M& val)
{
  for (long i = 0; i < (long)vec.size(); ++i) {
    if (vec[i] == val) {
      return true;
    }
  }
  return false;
}

bool does_string_have_tag(const std::string& str, const std::string& tag)
{
  const std::vector<std::string> tags = split_line_with_spaces(str);
  return does_vec_have(tags, tag);
}

template <class M>
void vector_append(std::vector<M>& v, const std::vector<M>& v1)
{
  v.insert(v.end(), v1.begin(), v1.end());
}

}  // namespace qlat
