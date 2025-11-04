#pragma once

#include <qlat-utils/types.h>

#include <map>
#include <set>
#include <vector>

namespace qlat
{  //

// -------------------

template <class M>
struct IsDataVectorType {
  using DataType = typename IsDataValueType<M>::DataType;
  using BasicDataType = typename IsDataValueType<DataType>::BasicDataType;
  using ElementaryType = typename IsDataValueType<DataType>::ElementaryType;
  static constexpr bool value = is_data_value_type<DataType>();
};

template <class M>
struct IsDataVectorType<Vector<M>> {
  using DataType = M;
  using BasicDataType = typename IsDataValueType<DataType>::BasicDataType;
  using ElementaryType = typename IsDataValueType<DataType>::ElementaryType;
  static constexpr bool value = is_data_value_type<DataType>();
};

template <class M>
struct IsDataVectorType<std::vector<M>> {
  using DataType = M;
  using BasicDataType = typename IsDataValueType<DataType>::BasicDataType;
  using ElementaryType = typename IsDataValueType<DataType>::ElementaryType;
  static constexpr bool value = is_data_value_type<DataType>();
};

template <class M, std::size_t N>
struct IsDataVectorType<Array<M, N>> {
  using DataType = M;
  using BasicDataType = typename IsDataValueType<DataType>::BasicDataType;
  using ElementaryType = typename IsDataValueType<DataType>::ElementaryType;
  static constexpr bool value = is_data_value_type<DataType>();
};

template <class M>
struct IsDataVectorType<Handle<M>> {
  using DataType = M;
  using BasicDataType = typename IsDataValueType<DataType>::BasicDataType;
  using ElementaryType = typename IsDataValueType<DataType>::ElementaryType;
  static constexpr bool value = is_data_value_type<DataType>();
};

template <class M>
struct IsDataVectorType<ConstHandle<M>> {
  using DataType = M;
  using BasicDataType = typename IsDataValueType<DataType>::BasicDataType;
  using ElementaryType = typename IsDataValueType<DataType>::ElementaryType;
  static constexpr bool value = is_data_value_type<DataType>();
};

template <>
struct IsDataVectorType<std::string> {
  using DataType = char;
  using BasicDataType = typename IsDataValueType<DataType>::BasicDataType;
  using ElementaryType = typename IsDataValueType<DataType>::ElementaryType;
  static constexpr bool value = is_data_value_type<DataType>();
};

// -------------------

template <class M>
qacc constexpr bool is_data_vector_type()
// data vector types
// data value types + all kinds of vector of data value types
{
  return IsDataVectorType<M>::value;
}

// -------------------

template <class M>
struct IsGetDataType {
  using DataType = typename IsDataVectorType<M>::DataType;
  using BasicDataType = typename IsDataValueType<DataType>::BasicDataType;
  using ElementaryType = typename IsDataValueType<DataType>::ElementaryType;
  static constexpr bool value = is_data_vector_type<M>();
};

// -------------------

template <class M>
qacc constexpr bool is_get_data_type()
// get data types
// basically data vector types
// support get_data function
{
  return IsGetDataType<M>::value;
}

// -------------------

template <class M, QLAT_ENABLE_IF(is_basic_data_type<M>())>
qacc Vector<M> get_data(const M& x)
{
  return Vector<M>(&x, 1);
}

template <class M, std::size_t N, QLAT_ENABLE_IF(is_data_value_type<M>())>
qacc Vector<M> get_data(const array<M, N>& arr)
{
  return Vector<M>(arr.data(), arr.size());
}

template <class M, std::size_t N, QLAT_ENABLE_IF(is_data_value_type<M>())>
qacc Vector<M> get_data(const std::array<M, N>& arr)
{
  return Vector<M>(arr.data(), arr.size());
}

template <class M, QLAT_ENABLE_IF(is_data_value_type<M>())>
qacc Vector<M> get_data(Vector<M> vec)
{
  return vec;
}

template <class M, std::size_t N, QLAT_ENABLE_IF(is_data_value_type<M>())>
qacc Vector<M> get_data(Array<M, N> arr)
{
  return Vector<M>(arr.data(), arr.size());
}

template <class M, QLAT_ENABLE_IF(is_data_value_type<M>())>
qacc Vector<M> get_data(const std::vector<M>& vec)
{
  return Vector<M>(vec.data(), vec.size());
}

template <class M, QLAT_ENABLE_IF(is_data_value_type<M>())>
qacc Vector<M> get_data(Handle<M> h)
{
  return Vector<M>(h.p, 1);
}

template <class M, QLAT_ENABLE_IF(is_data_value_type<M>())>
qacc Vector<M> get_data(ConstHandle<M> h)
{
  return Vector<M>(h.p, 1);
}

template <class M, QLAT_ENABLE_IF(is_data_value_type<M>())>
qacc Vector<M> get_data(Vector<M> vec, const Long size)
// only get a portion of the vec
// vec should be at least size Long
{
  qassert(vec.size() >= size);
  return Vector<M>(vec.data(), size);
}

Vector<char> get_data(const std::string& str);

// -------------------

template <class M>
qacc Long get_data_size(const M& x)
{
  return get_data(x).data_size();
}

// -------------------

template <class T, QLAT_ENABLE_IF(is_get_data_type<T>())>
qacc Vector<Char> get_data_char(const T& xx)
{
  using M = typename IsGetDataType<T>::DataType;
  const Vector<M> vec = get_data(xx);
  const Long n = vec.data_size();
  Vector<Char> v1((Char*)vec.data(), n);
  return v1;
}

template <class N, class T,
          class E1 = typename IsGetDataType<N>::ElementaryType,
          class E2 = typename IsGetDataType<T>::ElementaryType,
          QLAT_ENABLE_IF(is_get_data_type<T>())>
qacc Vector<N> get_data_as(const T& xx)
{
  static_assert(is_same<typename IsGetDataType<T>::ElementaryType,
                        typename IsGetDataType<N>::ElementaryType>(),
                "get_data_as type error");
  using M = typename IsGetDataType<T>::DataType;
  const Vector<M> vec = get_data(xx);
  const Long n = vec.data_size() / sizeof(N);
  Vector<N> v1((N*)vec.data(), n);
  qassert(v1.data_size() == vec.data_size());
  return v1;
}

template <class T, QLAT_ENABLE_IF(is_get_data_type<T>() and
                                  is_composed_of_real_d<
                                      typename IsGetDataType<T>::DataType>())>
qacc Vector<RealD> get_data_real_d(const T& xx)
{
  return get_data_as<RealD>(xx);
}

template <class T, QLAT_ENABLE_IF(is_get_data_type<T>() and
                                  is_composed_of_complex_d<
                                      typename IsGetDataType<T>::DataType>())>
qacc Vector<ComplexD> get_data_complex_d(const T& xx)
{
  return get_data_as<ComplexD>(xx);
}

template <class T, QLAT_ENABLE_IF(is_get_data_type<T>() and
                                  is_composed_of_long<
                                      typename IsGetDataType<T>::DataType>())>
qacc Vector<Long> get_data_long(const T& xx)
{
  return get_data_as<Long>(xx);
}

// -------------------

template <class M>
qacc Vector<M> get_data_one_elem(const M& x)
{
  return Vector<M>(&x, 1);
}

template <class M>
qacc Vector<Char> get_data_char(const Vector<M> v)
{
  Vector<Char> v1((Char*)v.data(), v.data_size());
  return v1;
}

template <class N, class M>
qacc Vector<N> get_data_as(const Vector<M> v)
{
  const Long n = v.data_size() / sizeof(N);
  Vector<N> v1((N*)v.data(), n);
  qassert(v1.data_size() == v.data_size());
  return v1;
}

template <class M>
qacc Vector<RealD> get_data_real_d(const Vector<M> v)
{
  return get_data_as<RealD>(v);
}

template <class M>
qacc Vector<ComplexD> get_data_complex_d(const Vector<M> v)
{
  return get_data_as<ComplexD>(v);
}

template <class M>
qacc Vector<Long> get_data_long(const Vector<M> v)
{
  return get_data_as<Long>(v);
}

// -------------------

template <class T, class E = typename IsGetDataType<T>::ElementaryType,
          QLAT_ENABLE_IF(is_get_data_type<T>())>
qacc Vector<E> get_data_in_elementary_type(const T& xx)
{
  using M = typename IsGetDataType<T>::DataType;
  static_assert(sizeof(M) % sizeof(E) == 0, "get_data_in_elementary_type");
  constexpr Int m = sizeof(M) / sizeof(E);
  const Vector<M> vec = get_data(xx);
  return Vector<E>((E*)vec.p, vec.n * m);
}

// -------------------

template <class T, QLAT_ENABLE_IF(is_get_data_type<T>())>
qacc void set_zero(T& xx)
{
  Vector<Char> vec = get_data_char(xx);
  memset(vec.data(), 0, vec.size());
}

template <class M>
qacc void set_zero(Vector<M> xx)
{
  Vector<Char> vec = get_data_char(xx);
  memset(vec.data(), 0, vec.size());
}

template <class M, std::size_t N>
qacc void set_zero(Array<M, N> xx)
{
  Vector<Char> vec = get_data_char(xx);
  memset(vec.data(), 0, vec.size());
}

// -------------------

qacc void assign(Vector<Char> xx, const Vector<Char>& yy)
{
  qassert(xx.size() == yy.size());
#ifndef QLAT_IN_ACC
  std::memcpy((void*)xx.data(), (void*)yy.data(), xx.size());
#else
  const Long num = xx.size() / sizeof(Long);
  const Long offset = num * sizeof(Long);
  const Long rem = xx.size() - offset;
  Long* px = (Long*)xx.data();
  const Long* py = (Long*)yy.data();
  for (Long i = 0; i < num; ++i) {
    px[i] = py[i];
  }
  if (rem > 0) {
    Char* px = xx.data();
    const Char* py = yy.data();
    for (Long i = offset; i < xx.size(); ++i) {
      px[i] = py[i];
    }
  }
#endif
}

qacc void assign(Vector<Long> xx, const Vector<Long>& yy)
{
  qassert(xx.size() == yy.size());
#ifndef QLAT_IN_ACC
  std::memcpy((void*)xx.data(), (void*)yy.data(), sizeof(Long) * xx.size());
#else
  const Long num = xx.size();
  Long* px = xx.data();
  const Long* py = yy.data();
  for (Long i = 0; i < num; ++i) {
    px[i] = py[i];
  }
#endif
}

template <class T1, class T2,
          class E1 = typename IsGetDataType<T1>::ElementaryType,
          class E2 = typename IsGetDataType<T2>::ElementaryType,
          QLAT_ENABLE_IF(is_get_data_type<T1>() and
                         is_get_data_type<T2>())>
qacc void assign(T1& xx, const T2& yy)
{
  // static_assert(is_same<E1, E2>(), "assign type error");
  // Not used as function may be used in python binding.
  if (not is_same<E1, E2>()) {
    qerr(ssprintf("assign type mismatch: %s %s",
                  IsBasicDataType<E1>::get_type_name().c_str(),
                  IsBasicDataType<E2>::get_type_name().c_str()));
  }
  Vector<Char> vx = get_data_char(xx);
  const Vector<Char> vy = get_data_char(yy);
  assign(vx, vy);
}

template <
    class M1, class T2, class E1 = typename IsGetDataType<M1>::ElementaryType,
    class E2 = typename IsGetDataType<T2>::ElementaryType,
    QLAT_ENABLE_IF((is_data_value_type<M1>() and is_get_data_type<T2>()))>
qacc void assign(Vector<M1> xx, const T2& yy)
{
  return assign<Vector<M1>, T2>(xx, yy);
}

// -------------------

template <class M, QLAT_ENABLE_IF(is_number<M>())>
qacc void set_unit(M& x, const Long& coef = 1)
{
  x = coef;
}

template <class M, QLAT_ENABLE_IF(is_number<M>())>
qacc void set_unit(M& x, const RealD& coef = 1.0)
{
  x = coef;
}

template <class M, QLAT_ENABLE_IF(is_complex<M>())>
qacc void set_unit(M& x, const ComplexD& coef = 1.0)
{
  x = coef;
}

template <class M, QLAT_ENABLE_IF(is_integer<M>() or is_real<M>())>
qacc void set_unit(M& x, const ComplexD& coef = 1.0)
{
  x = coef.real();
}

// -------------------

template <class M, QLAT_ENABLE_IF(is_integer<M>())>
qacc RealD qnorm(const M& x)
{
  return x * x;
}

template <class T, QLAT_ENABLE_IF(is_get_data_type<T>() and
                                  not is_basic_data_type<T>())>
qacc RealD qnorm(const T& xx)
{
  using M = typename IsGetDataType<T>::DataType;
  const Vector<M> vec = get_data(xx);
  RealD sum = 0.0;
  for (Long i = 0; i < vec.size(); ++i) {
    sum += qnorm(vec[i]);
  }
  return sum;
}

// -------------------

template <class M>
void clear(std::vector<M>& vec)
{
  std::vector<M> empty;
  std::swap(empty, vec);
}

// -------------------

template <class T,
          QLAT_ENABLE_IF(is_get_data_type<T>() and (not is_number<T>()))>
qacc bool operator==(const T& x1, const T& x2)
{
  const Vector<Char> v1 = get_data_char(x1);
  const Vector<Char> v2 = get_data_char(x2);
  if (v1.size() != v2.size()) {
    return false;
  }
  const Int cmp = std::memcmp(v1.data(), v2.data(), v1.size());
  return cmp == 0;
}

template <class T,
          QLAT_ENABLE_IF(is_get_data_type<T>() and (not is_number<T>()))>
qacc bool operator!=(const T& x1, const T& x2)
{
  return not(x1 == x2);
}

// -------------------

}  // namespace qlat
