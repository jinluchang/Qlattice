#pragma once

#include <qlat-utils/types.h>

#include <map>
#include <set>
#include <vector>

namespace qlat
{  //

// -------------------

template <class M>
struct HasGetData {
  static constexpr bool value = is_data_value_type<M>();
  using DataType = M;
};

template <class M>
struct HasGetData<Vector<M>> {
  static constexpr bool value = is_data_value_type<M>();
  using DataType = M;
};

template <class M>
struct HasGetData<std::vector<M>> {
  static constexpr bool value = is_data_value_type<M>();
  using DataType = M;
};

template <class M, Long N>
struct HasGetData<Array<M, N>> {
  static constexpr bool value = is_data_value_type<M>();
  using DataType = M;
};

template <class M, size_t N>
struct HasGetData<array<M, N>> {
  static constexpr bool value = is_data_value_type<M>();
  using DataType = M;
};

template <class M>
struct HasGetData<Handle<M>> {
  static constexpr bool value = is_data_value_type<M>();
  using DataType = M;
};

template <class M>
struct HasGetData<ConstHandle<M>> {
  static constexpr bool value = is_data_value_type<M>();
  using DataType = M;
};

template <class M>
qacc constexpr bool has_get_data()
{
  return HasGetData<M>::value;
};

// -------------------

template <class M, QLAT_ENABLE_IF(is_data_value_type<M>())>
qacc Vector<M> get_data(const M& x)
{
  return Vector<M>(&x, 1);
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

template <class M, size_t N, QLAT_ENABLE_IF(is_data_value_type<M>())>
qacc Vector<M> get_data(const array<M, N>& arr)
{
  return Vector<M>(arr.data(), arr.size());
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

template <class T, QLAT_ENABLE_IF(has_get_data<T>())>
void set_zero(T& xx)
{
  using M = typename HasGetData<T>::DataType;
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

template <class T,
          QLAT_ENABLE_IF(has_get_data<T>() and not is_data_value_type<T>())>
RealD qnorm(const T& xx)
{
  using M = typename HasGetData<T>::DataType;
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
