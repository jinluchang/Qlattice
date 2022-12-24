#pragma once

#include <qlat-utils/qutils-extra.h>
#include <qlat-utils/rng-state.h>

namespace qlat
{  //

template <class M>
qacc array<M, 0> make_array()
{
  array<M, 0> arr;
  return arr;
}

template <class M>
qacc array<M, 1> make_array(const M& x)
{
  array<M, 1> arr;
  arr[0] = x;
  return arr;
}

template <class M>
qacc array<M, 2> make_array(const M& x, const M& x1)
{
  array<M, 2> arr;
  arr[0] = x;
  arr[1] = x1;
  return arr;
}

template <class M>
qacc array<M, 3> make_array(const M& x, const M& x1, const M& x2)
{
  array<M, 3> arr;
  arr[0] = x;
  arr[1] = x1;
  arr[2] = x2;
  return arr;
}

template <class M>
qacc array<M, 4> make_array(const M& x, const M& x1, const M& x2, const M& x3)
{
  array<M, 4> arr;
  arr[0] = x;
  arr[1] = x1;
  arr[2] = x2;
  arr[3] = x3;
  return arr;
}

template <class M>
qacc array<M, 5> make_array(const M& x, const M& x1, const M& x2, const M& x3,
                            const M& x4)
{
  array<M, 5> arr;
  arr[0] = x;
  arr[1] = x1;
  arr[2] = x2;
  arr[3] = x3;
  arr[4] = x4;
  return arr;
}

template <class M>
qacc array<M, 6> make_array(const M& x, const M& x1, const M& x2, const M& x3,
                            const M& x4, const M& x5)
{
  array<M, 6> arr;
  arr[0] = x;
  arr[1] = x1;
  arr[2] = x2;
  arr[3] = x3;
  arr[4] = x4;
  arr[5] = x5;
  return arr;
}

template <class M>
qacc array<M, 7> make_array(const M& x, const M& x1, const M& x2, const M& x3,
                            const M& x4, const M& x5, const M& x6)
{
  array<M, 7> arr;
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
qacc array<M, 8> make_array(const M& x, const M& x1, const M& x2, const M& x3,
                            const M& x4, const M& x5, const M& x6, const M& x7)
{
  array<M, 8> arr;
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
qacc array<M, 9> make_array(const M& x, const M& x1, const M& x2, const M& x3,
                            const M& x4, const M& x5, const M& x6, const M& x7,
                            const M& x8)
{
  array<M, 9> arr;
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
qacc array<M, 10> make_array(const M& x, const M& x1, const M& x2, const M& x3,
                             const M& x4, const M& x5, const M& x6, const M& x7,
                             const M& x8, const M& x9)
{
  array<M, 10> arr;
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
qacc void set_zero(Vector<M> vec)
{
  std::memset((void*)vec.data(), 0, vec.data_size());
}

template <class M, int N>
qacc void set_zero(Array<M, N> arr)
{
  long size = N * sizeof(M);
  std::memset((void*)arr.data(), 0, size);
}

template <class M, int N>
qacc Vector<M> get_data(Array<M, N> arr)
{
  return Vector<M>(arr.data(), arr.size());
}

template <class M>
qacc Vector<M> get_data(Vector<M> vec)
{
  return vec;
}

template <class M>
qacc Vector<M> get_data(Vector<M> vec, const long size)
// only get a portion of the vec
// vec should be at least size long
{
  qassert(vec.size() >= size);
  return Vector<M>(vec.data(), size);
}

template <class M, unsigned long N>
qacc Vector<M> get_data(const array<M, N>& vec)
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
qacc Vector<M> get_data(const Handle<M>& h)
{
  return Vector<M>(h.p, 1);
}

template <class M>
qacc Vector<M> get_data(const ConstHandle<M>& h)
{
  return Vector<M>(h.p, 1);
}

template <class M>
qacc Vector<M> get_data_one_elem(const M& x)
{
  return Vector<M>(&x, 1);
}

qacc Vector<long> get_data(const long& x) { return get_data_one_elem(x); }

qacc Vector<double> get_data(const double& x) { return get_data_one_elem(x); }

qacc Vector<int> get_data(const int& x) { return get_data_one_elem(x); }

qacc Vector<float> get_data(const float& x) { return get_data_one_elem(x); }

template <class T>
qacc double qnorm(const Vector<T>& mm)
{
  double sum = 0.0;
  const long size = mm.size();
  for (long i = 0; i < size; ++i) {
    sum += qnorm(mm[i]);
  }
  return sum;
}

template <class T>
qacc double qnorm(const Vector<T>& m1, const Vector<T>& m2)
{
  double sum = 0.0;
  const long size = m1.size();
  qassert(size == (long)m2.size());
  for (long i = 0; i < size; ++i) {
    sum += qnorm(m1[i], m2[i]);
  }
  return sum;
}

template <class N, class M>
qacc Vector<N> get_data_as(const Vector<M>& v)
{
  const long n = v.data_size() / sizeof(N);
  Vector<N> v1((N*)v.data(), n);
  qassert(v1.data_size() == v.data_size());
  return v1;
}

template <class M>
qacc Vector<double> get_data_double(const Vector<M>& v)
{
  return get_data_as<double>(v);
}

template <class M>
qacc Vector<Complex> get_data_complex(const Vector<M>& v)
{
  return get_data_as<Complex>(v);
}

template <class M>
qacc Vector<long> get_data_long(const Vector<M>& v)
{
  return get_data_as<long>(v);
}

template <class M>
qacc long get_data_size(const M& x)
{
  return get_data(x).data_size();
}

template <class T>
qacc double qnorm_double(const Vector<T>& m1, const Vector<T>& m2)
{
  const Vector<double> dm1((double*)m1.data(), m1.data_size() / sizeof(double));
  const Vector<double> dm2((double*)m2.data(), m2.data_size() / sizeof(double));
  return qnorm(dm1, dm2);
}

template <class M>
qacc void to_from_little_endian_16(Vector<M> v)
{
  to_from_little_endian_16((void*)v.data(), v.data_size());
}

template <class M>
qacc void to_from_little_endian_32(Vector<M> v)
{
  to_from_little_endian_32((void*)v.data(), v.data_size());
}

template <class M>
qacc void to_from_little_endian_64(Vector<M> v)
{
  to_from_little_endian_64((void*)v.data(), v.data_size());
}

template <class M>
qacc void to_from_big_endian_16(Vector<M> v)
{
  to_from_big_endian_16((void*)v.data(), v.data_size());
}

template <class M>
qacc void to_from_big_endian_32(Vector<M> v)
{
  to_from_big_endian_32((void*)v.data(), v.data_size());
}

template <class M>
qacc void to_from_big_endian_64(Vector<M> v)
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
qacc void assign(array<M, N>& vec, const Array<M, N>& src)
{
  std::memcpy((void*)vec.data(), (void*)src.data(), src.data_size());
}

template <class M, int N>
qacc void assign(array<M, N>& vec, const Vector<M>& src)
{
  qassert(N == src.size());
  std::memcpy((void*)vec.data(), (void*)src.data(), src.data_size());
}

template <class M, int N>
qacc void assign(std::vector<M>& vec, const Array<M, N>& src)
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
qacc void assign(Vector<M> vec, const Vector<M>& src)
{
  qassert(vec.size() == src.size());
  std::memcpy((void*)vec.data(), (void*)src.data(), src.data_size());
}

template <class M, class N>
qacc void assign(Vector<M> vec, const Vector<N>& src)
{
  qassert(vec.data_size() == src.data_size());
  std::memcpy((void*)vec.data(), (void*)src.data(), src.data_size());
}

template <class M, int N>
qacc void assign(Vector<M> vec, const Array<M, N>& src)
{
  qassert(vec.size() == N);
  std::memcpy((void*)vec.data(), (void*)src.data(), src.data_size());
}

template <class M, int N>
qacc void assign(Array<M, N> vec, const Array<M, N>& src)
{
  std::memcpy((void*)vec.data(), (void*)src.data(), src.data_size());
}

template <class M, int N>
qacc void assign(Array<M, N> vec, const Vector<M>& src)
{
  qassert(src.size() == N);
  std::memcpy((void*)vec.data(), (void*)src.data(), src.data_size());
}

template <class M, class N>
qacc void assign_truncate(M& x, const N& y)
{
  if (sizeof(M) <= sizeof(N)) {
    std::memcpy((void*)&x, (void*)&y, sizeof(M));
  } else {
    // if M has a larger size, than leave the extra space untouched
    std::memcpy((void*)&x, (void*)&y, sizeof(N));
  }
}

qacc bool is_integer(const double& x)
{
  const double diff = x - (long)x;
  return 1e-6 > diff || diff > 1 - 1e-6;
}

template <class M>
qacc bool is_integer(const std::vector<M>& v)
{
  for (int i = 0; i < (int)v.size(); ++i) {
    if (!is_integer(v[i])) {
      return false;
    }
  }
  return true;
}

template <class M, unsigned long N>
qacc bool is_integer(const array<M, N>& v)
{
  for (int i = 0; i < N; ++i) {
    if (!is_integer(v[i])) {
      return false;
    }
  }
  return true;
}

template <class M, int N>
qacc Array<M, N> operator+=(Array<M, N> v, const Array<M, N> v1)
{
  for (int i = 0; i < N; ++i) {
    v.p[i] += v1.p[i];
  }
  return v;
}

template <class M, int N>
qacc Array<M, N> operator-=(Array<M, N> v, const Array<M, N> v1)
{
  for (int i = 0; i < N; ++i) {
    v.p[i] -= v1.p[i];
  }
  return v;
}

template <class M, int N>
qacc Array<M, N> operator*=(Array<M, N> v, const double factor)
{
  for (int i = 0; i < N; ++i) {
    v.p[i] *= factor;
  }
  return v;
}

template <class M, int N>
qacc Array<M, N> operator*=(Array<M, N> v, const Complex factor)
{
  for (int i = 0; i < N; ++i) {
    v.p[i] *= factor;
  }
  return v;
}

template <class M>
qacc Vector<M> operator+=(Vector<M> v, const Vector<M> v1)
{
  qassert(v.size() == v1.size());
  for (long i = 0; i < v.size(); ++i) {
    v.p[i] += v1.p[i];
  }
  return v;
}

template <class M>
qacc Vector<M> operator-=(Vector<M> v, const Vector<M> v1)
{
  qassert(v.size() == v1.size());
  for (long i = 0; i < v.size(); ++i) {
    v.p[i] -= v1.p[i];
  }
  return v;
}

template <class M>
qacc Vector<M> operator*=(Vector<M> v, const double factor)
{
  for (long i = 0; i < v.size(); ++i) {
    v.p[i] *= factor;
  }
  return v;
}

template <class M>
qacc Vector<M> operator*=(Vector<M> v, const Complex factor)
{
  for (long i = 0; i < v.size(); ++i) {
    v.p[i] *= factor;
  }
  return v;
}

template <class M, unsigned long N>
qacc array<M, N> operator+(const array<M, N>& v1, const array<M, N>& v2)
{
  array<M, N> ret;
  for (unsigned long i = 0; i < N; ++i) {
    ret[i] = v1[i] + v2[i];
  }
  return ret;
}

template <class M, unsigned long N>
qacc array<M, N> operator-(const array<M, N>& v1, const array<M, N>& v2)
{
  array<M, N> ret;
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
void set_u_rand_double(Vector<M> v, const RngState& rs,
                       const double upper = 1.0, const double lower = -1.0)
{
  RngState rsi = rs;
  Vector<double> dv((double*)v.data(), v.data_size() / sizeof(double));
  for (int m = 0; m < dv.size(); ++m) {
    dv[m] = u_rand_gen(rsi, upper, lower);
  }
}

template <class M>
void set_u_rand_float(Vector<M> v, const RngState& rs, const double upper = 1.0,
                      const double lower = -1.0)
{
  RngState rsi = rs;
  Vector<float> dv((float*)v.data(), v.data_size() / sizeof(float));
  for (int m = 0; m < dv.size(); ++m) {
    dv[m] = u_rand_gen(rsi, upper, lower);
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

inline bool does_string_have_tag(const std::string& str, const std::string& tag)
{
  const std::vector<std::string> tags = split_line_with_spaces(str);
  return does_vec_have(tags, tag);
}

template <class M>
void vector_append(std::vector<M>& v, const std::vector<M>& v1)
{
  v.insert(v.end(), v1.begin(), v1.end());
}

template <class M>
std::vector<long> vector_map_size(const std::vector<std::vector<M> >& datatable)
{
  std::vector<long> row_sizes(datatable.size());
  for (size_t i = 0; i < datatable.size(); ++i) {
    const std::vector<M>& row = datatable[i];
    row_sizes[i] = row.size();
  }
  return row_sizes;
}

template <class M>
std::vector<M> vector_concat(const std::vector<std::vector<M> >& datatable)
{
  size_t total_size = 0;
  for (size_t i = 0; i < datatable.size(); ++i) {
    const std::vector<M>& row = datatable[i];
    total_size += row.size();
  }
  std::vector<M> data(total_size);
  size_t count = 0;
  for (size_t i = 0; i < datatable.size(); ++i) {
    const std::vector<M>& row = datatable[i];
    for (size_t j = 0; j < row.size(); ++j) {
      data[count] = row[j];
      count += 1;
    }
  }
  return data;
}

template <class M>
std::vector<std::vector<M> > vector_split(const std::vector<M>& data,
                                          const std::vector<long>& row_sizes)
{
  std::vector<std::vector<M> > datatable;
  datatable.resize(row_sizes.size());
  size_t count = 0;
  for (size_t i = 0; i < datatable.size(); ++i) {
    std::vector<M>& row = datatable[i];
    row.resize(row_sizes[i]);
    for (size_t j = 0; j < row.size(); ++j) {
      row[j] = data[count];
      count += 1;
    }
  }
  return datatable;
}

}  // namespace qlat
