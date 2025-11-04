#pragma once

#include <qlat-utils/utils-extra.h>
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

// template <class T>
// qacc RealD qnorm(const Vector<T>& mm)
// {
//   RealD sum = 0.0;
//   const Long size = mm.size();
//   for (Long i = 0; i < size; ++i) {
//     sum += qnorm(mm[i]);
//   }
//   return sum;
// }

template <class T>
qacc RealD qnorm(const Vector<T>& m1, const Vector<T>& m2)
{
  RealD sum = 0.0;
  const Long size = m1.size();
  qassert(size == (Long)m2.size());
  for (Long i = 0; i < size; ++i) {
    sum += qnorm(m1[i], m2[i]);
  }
  return sum;
}

template <class T>
qacc RealD qnorm_double(const Vector<T>& m1, const Vector<T>& m2)
{
  const Vector<RealD> dm1((RealD*)m1.data(), m1.data_size() / sizeof(RealD));
  const Vector<RealD> dm2((RealD*)m2.data(), m2.data_size() / sizeof(RealD));
  return qnorm(dm1, dm2);
}

template <class M, class N>
qacc void assign_direct(M& x, const N& y)
{
  x = y;
}

template <class M, class N>
qacc void iadd_direct(M& x, const N& y)
{
  x += y;
}

template <class M, class N>
qacc void isub_direct(M& x, const N& y)
{
  x -= y;
}

template <class M, class N>
qacc void imul_direct(M& x, const N& y)
{
  x *= y;
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

qacc bool is_integer(const RealD& x)
{
  const RealD diff = x - (Long)x;
  return 1e-6 > diff || diff > 1 - 1e-6;
}

template <class M>
qacc bool is_integer(const std::vector<M>& v)
{
  for (Int i = 0; i < (Int)v.size(); ++i) {
    if (!is_integer(v[i])) {
      return false;
    }
  }
  return true;
}

template <class M, std::size_t N>
qacc bool is_integer(const array<M, N>& v)
{
  for (Int i = 0; i < (Int)N; ++i) {
    if (!is_integer(v[i])) {
      return false;
    }
  }
  return true;
}

template <class M, std::size_t N>
qacc Array<M, N> operator+=(Array<M, N> v, const Array<M, N> v1)
{
  for (Int i = 0; i < (Int)N; ++i) {
    v.p[i] += v1.p[i];
  }
  return v;
}

template <class M, std::size_t N>
qacc Array<M, N> operator-=(Array<M, N> v, const Array<M, N> v1)
{
  for (Int i = 0; i < (Int)N; ++i) {
    v.p[i] -= v1.p[i];
  }
  return v;
}

template <class M, std::size_t N>
qacc Array<M, N> operator*=(Array<M, N> v, const RealD factor)
{
  for (Int i = 0; i < (Int)N; ++i) {
    v.p[i] *= factor;
  }
  return v;
}

template <class M, std::size_t N>
qacc Array<M, N> operator*=(Array<M, N> v, const ComplexD factor)
{
  for (Int i = 0; i < (Int)N; ++i) {
    v.p[i] *= factor;
  }
  return v;
}

template <class M>
qacc Vector<M> operator+=(Vector<M> v, const Vector<M> v1)
{
  qassert(v.size() == v1.size());
  for (Long i = 0; i < v.size(); ++i) {
    v.p[i] += v1.p[i];
  }
  return v;
}

template <class M>
qacc Vector<M> operator-=(Vector<M> v, const Vector<M> v1)
{
  qassert(v.size() == v1.size());
  for (Long i = 0; i < v.size(); ++i) {
    v.p[i] -= v1.p[i];
  }
  return v;
}

template <class M>
qacc Vector<M> operator*=(Vector<M> v, const RealD factor)
{
  for (Long i = 0; i < v.size(); ++i) {
    v.p[i] *= factor;
  }
  return v;
}

template <class M>
qacc Vector<M> operator*=(Vector<M> v, const ComplexD factor)
{
  for (Long i = 0; i < v.size(); ++i) {
    v.p[i] *= factor;
  }
  return v;
}

template <class M, std::size_t N>
qacc array<M, N> operator+(const array<M, N>& v1, const array<M, N>& v2)
{
  array<M, N> ret;
  for (std::size_t i = 0; i < N; ++i) {
    ret[i] = v1[i] + v2[i];
  }
  return ret;
}

template <class M, std::size_t N>
qacc array<M, N> operator-(const array<M, N>& v1, const array<M, N>& v2)
{
  array<M, N> ret;
  for (std::size_t i = 0; i < N; ++i) {
    ret[i] = v1[i] - v2[i];
  }
  return ret;
}

template <class M>
void set_u_rand(Vector<M> v, const RngState& rs, const RealD upper = 1.0,
                const RealD lower = -1.0)
{
  if (not is_composed_of_real<M>()) {
    qassert(is_composed_of_real<M>());
    return;
  }
  using Real = typename IsDataValueType<M>::ElementaryType;
  RngState rsi = rs;
  Vector<Real> dv((Real*)v.data(), v.data_size() / sizeof(Real));
  for (Int m = 0; m < dv.size(); ++m) {
    dv[m] = u_rand_gen(rsi, upper, lower);
  }
}

template <class M>
bool does_vec_have(const std::vector<M>& vec, const M& val)
{
  for (Long i = 0; i < (Long)vec.size(); ++i) {
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
std::vector<Long> vector_map_size(const std::vector<std::vector<M> >& datatable)
{
  std::vector<Long> row_sizes(datatable.size());
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
std::vector<std::vector<M>> vector_split(const vector<M>& data,
                                         const vector<Long>& row_sizes)
{
  std::vector<std::vector<M>> datatable;
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

void convert_double_from_float(Vector<RealD> vd, const Vector<RealF> vf);

void convert_float_from_double(Vector<RealF> vf, const Vector<RealD> vd);

}  // namespace qlat
