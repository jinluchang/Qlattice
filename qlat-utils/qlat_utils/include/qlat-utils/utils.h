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

template <class M, Long N>
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

template <class M, size_t N, QLAT_ENABLE_IF(is_data_value_type<M>())>
qacc Vector<M> get_data(const array<M, N>& arr)
{
  return Vector<M>(arr.data(), arr.size());
}

template <class M, QLAT_ENABLE_IF(is_data_value_type<M>())>
qacc Vector<M> get_data(Vector<M> vec)
{
  return vec;
}

template <class M, Long N, QLAT_ENABLE_IF(is_data_value_type<M>())>
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
// vec should be at least size long
{
  qassert(vec.size() >= size);
  return Vector<M>(vec.data(), size);
}

inline Vector<char> get_data(const std::string& str)
{
  return Vector<char>(&str[0], str.length());
}

// -------------------

template <class T, QLAT_ENABLE_IF(is_get_data_type<T>())>
qacc Vector<char> get_data_char(const T& xx)
{
  using M = typename IsGetDataType<T>::DataType;
  const Vector<M> vec = get_data(xx);
  const Long n = vec.data_size();
  Vector<char> v1((char*)vec.data(), n);
  return v1;
}

template <class N, class T, QLAT_ENABLE_IF(is_get_data_type<T>())>
qacc Vector<N> get_data_as(const T& xx)
{
  static_assert(is_same<typename IsDataVectorType<T>::ElementaryType,
                        typename IsDataValueType<N>::ElementaryType>(),
                "get_data_as type error");
  using M = typename IsGetDataType<T>::DataType;
  const Vector<M> vec = get_data(xx);
  const Long n = vec.data_size() / sizeof(N);
  Vector<N> v1((N*)vec.data(), n);
  qassert(v1.data_size() == vec.data_size());
  return v1;
}

template <class T,
          QLAT_ENABLE_IF(is_get_data_type<T>() and
                         is_composed_of_real_d<typename IsGetDataType<T>::DataType>())>
qacc Vector<RealD> get_data_real_d(const T& xx)
{
  return get_data_as<RealD>(xx);
}

template <class T, QLAT_ENABLE_IF(
                       is_get_data_type<T>() and
                       is_composed_of_complex_d<typename IsGetDataType<T>::DataType>())>
qacc Vector<ComplexD> get_data_complex_d(const T& xx)
{
  return get_data_as<ComplexD>(xx);
}

template <class T,
          QLAT_ENABLE_IF(is_get_data_type<T>() and
                         is_composed_of_long<typename IsGetDataType<T>::DataType>())>
qacc Vector<Long> get_data_long(const T& xx)
{
  return get_data_as<Long>(xx);
}

// -------------------

template <class T, class E = typename IsDataVectorType<T>::ElementaryType,
          QLAT_ENABLE_IF(is_get_data_type<T>())>
qacc Vector<E> get_data_in_elementary_type(const T& xx)
{
  using M = typename IsGetDataType<T>::DataType;
  static_assert(sizeof(M) % sizeof(E) == 0, "get_data_in_elementary_type");
  constexpr int m = sizeof(M) / sizeof(E);
  const Vector<M> vec = get_data(xx);
  return Vector<E>((E*)vec.p, vec.n * m);
}

// -------------------

template <class T, QLAT_ENABLE_IF(is_get_data_type<T>())>
qacc void set_zero(T& xx)
{
  using M = typename IsGetDataType<T>::DataType;
  Vector<M> vec = get_data(xx);
  Long size = vec.size() * sizeof(M);
  std::memset((void*)vec.data(), 0, size);
}

template <class T1, class T2,
          class E1 = typename IsDataVectorType<T1>::ElementaryType,
          class E2 = typename IsDataVectorType<T2>::ElementaryType,
          QLAT_ENABLE_IF(is_data_vector_type<T1>() and
                         is_data_vector_type<T2>())>
qacc void assign(T1& xx, const T2& yy)
{
  // static_assert(is_same<E1, E2>(), "assign type error");
  // Not used as function may be used in python binding.
  if (not is_same<E1, E2>()) {
    qerr(ssprintf("assign type mismatch: %s %s",
                  IsBasicDataType<E1>::get_type_name().c_str(),
                  IsBasicDataType<E2>::get_type_name().c_str()));
  }
  Vector<E1> vx = get_data_in_elementary_type(xx);
  const Vector<E2> vy = get_data_in_elementary_type(yy);
  qassert(vx.size() == vy.size());
  std::memcpy((void*)vx.data(), (void*)vy.data(), vx.size() * sizeof(E1));
}

template <
    class M1, class T2, class E1 = typename IsDataValueType<M1>::ElementaryType,
    class E2 = typename IsDataVectorType<T2>::ElementaryType,
    QLAT_ENABLE_IF((is_data_value_type<M1>() and is_data_vector_type<T2>()))>
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

template <class T, QLAT_ENABLE_IF(is_data_vector_type<T>() and
                                  not is_basic_data_type<T>())>
RealD qnorm(const T& xx)
{
  using M = typename IsDataVectorType<T>::DataType;
  const Vector<M> vec = get_data(xx);
  RealD sum = 0.0;
  for (Long i = 0; i < vec.size(); ++i) {
    sum += qnorm(vec[i]);
  }
  return sum;
}

// -------------------

template <class M, class N, QLAT_ENABLE_IF(is_real<M>() and is_real<N>())>
qacc RealD qnorm(const M& x, const N& y)
{
  return x * y;
}

// -------------------

template <class M>
void clear(std::vector<M>& vec)
{
  std::vector<M> empty;
  std::swap(empty, vec);
}

// -------------------

}  // namespace qlat
