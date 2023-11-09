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
  static constexpr bool value = is_data_value_type<M>();
  using DataType = typename IsDataValueType<M>::DataType;
};

template <class M>
struct IsDataVectorType<Vector<M>> {
  static constexpr bool value = is_data_value_type<M>();
  using DataType = M;
};

template <class M>
struct IsDataVectorType<std::vector<M>> {
  static constexpr bool value = is_data_value_type<M>();
  using DataType = M;
};

template <class M, Long N>
struct IsDataVectorType<Array<M, N>> {
  static constexpr bool value = is_data_value_type<M>();
  using DataType = M;
};

template <class M>
struct IsDataVectorType<Handle<M>> {
  static constexpr bool value = is_data_value_type<M>();
  using DataType = M;
};

template <class M>
struct IsDataVectorType<ConstHandle<M>> {
  static constexpr bool value = is_data_value_type<M>();
  using DataType = M;
};

// -------------------

template <class M>
qacc constexpr bool is_data_vector_type()
// data vector types
// data value types + all kinds of vector of data value types
{
  return IsDataVectorType<M>::value;
};

// -------------------

template <class M>
struct IsGetDataType {
  static constexpr bool value = is_data_vector_type<M>();
  using DataType = typename IsDataVectorType<M>::DataType;
};

// -------------------

template <class M>
qacc constexpr bool is_get_data_type()
// get data types
// data vector types + many other containers of data value types
// support get_data function
{
  return IsGetDataType<M>::value;
};

// -------------------

template <class M, QLAT_ENABLE_IF(is_data_value_type<M>())>
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
qacc Vector<M> get_data(Handle<M>& h)
{
  return Vector<M>(h.p, 1);
}

template <class M, QLAT_ENABLE_IF(is_data_value_type<M>())>
qacc Vector<M> get_data(ConstHandle<M>& h)
{
  return Vector<M>(h.p, 1);
}

// -------------------

template <class T, QLAT_ENABLE_IF(is_get_data_type<T>())>
qacc void set_zero(T& xx)
{
  using M = typename IsGetDataType<T>::DataType;
  const Vector<M> vec = get_data(xx);
  Long size = vec.size() * sizeof(M);
  std::memset((void*)vec.data(), 0, size);
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
