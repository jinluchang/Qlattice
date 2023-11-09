#pragma once

#include <qlat-utils/mat-vec.h>
#include <qlat-utils/rng-state.h>

#include <cstdint>
#include <vector>

namespace qlat
{  //

using DataTable = std::vector<std::vector<RealD>>;

using crc32_t = uint32_t;

// -------------------------------------------------------------------------

template <class M, size_t N>
qacc bool qisnan(const array<M, N>& arr)
{
  for (int i = 0; i < (int)N; ++i) {
    if (qisnan(arr[i])) {
      return true;
    }
  }
  return false;
}

template <class M>
bool qisnan(const std::vector<M>& arr)
{
  for (size_t i = 0; i < arr.size(); ++i) {
    if (qisnan(arr[i])) {
      return true;
    }
  }
  return false;
}

// --------------------

struct API Coordinate : public array<Int, DIMN> {
  qacc Coordinate() { array<Int, DIMN>::fill(0); }
  qacc Coordinate(Int first, Int second, Int third, Int fourth)
  {
    Int* p = data();
    p[0] = first;
    p[1] = second;
    p[2] = third;
    p[3] = fourth;
  }
};

struct API CoordinateD : public array<RealD, DIMN> {
  qacc CoordinateD() { memset(this, 0, sizeof(CoordinateD)); }
  qacc CoordinateD(const array<RealD, DIMN>& arr)
  {
    CoordinateD& c = *this;
    c = arr;
    qassert(false == qisnan(c));
  }
  qacc CoordinateD(const RealD x0, const RealD x1, const RealD x2,
                   const RealD x3)
  {
    qassert(DIMN == 4);
    CoordinateD& c = *this;
    c[0] = x0;
    c[1] = x1;
    c[2] = x2;
    c[3] = x3;
    qassert(false == qisnan(c));
  }
  qacc CoordinateD(const Coordinate& x)
  {
    CoordinateD& c = *this;
    for (int i = 0; i < DIMN; ++i) {
      c[i] = x[i];
    }
  }
};

// -------------------------------------------------------------------------

template <class M>
std::string get_type_name()
{
  std::string ret = "unknown";
  if (is_same<M, char>()) {
    ret = "char";
  } else if (is_same<M, signed char>()) {
    ret = "signed_char";
  } else if (is_same<M, unsigned char>()) {
    ret = "unsigned_char";
  } else if (is_same<M, int8_t>()) {
    ret = "Int8t";
  } else if (is_same<M, int16_t>()) {
    ret = "Int16t";
  } else if (is_same<M, int32_t>()) {
    ret = "Int32t";
  } else if (is_same<M, int64_t>()) {
    ret = "Int64t";
  } else if (is_same<M, uint8_t>()) {
    ret = "UInt8t";
  } else if (is_same<M, uint16_t>()) {
    ret = "UInt16t";
  } else if (is_same<M, uint32_t>()) {
    ret = "UInt32t";
  } else if (is_same<M, uint64_t>()) {
    ret = "UInt64t";
  } else if (is_same<M, RealF>()) {
    ret = "RealF";
  } else if (is_same<M, RealD>()) {
    ret = "RealD";
  } else if (is_same<M, ComplexD>()) {
    ret = "ComplexD";
  } else if (is_same<M, ComplexF>()) {
    ret = "ComplexF";
  } else if (is_same<M, ColorMatrixD>()) {
    ret = "ColorMatrixD";
  } else if (is_same<M, ColorMatrixF>()) {
    ret = "ColorMatrixF";
  } else if (is_same<M, WilsonMatrixD>()) {
    ret = "WilsonMatrixD";
  } else if (is_same<M, WilsonMatrixF>()) {
    ret = "WilsonMatrixF";
  } else if (is_same<M, SpinMatrixD>()) {
    ret = "SpinMatrixD";
  } else if (is_same<M, SpinMatrixF>()) {
    ret = "SpinMatrixF";
  } else if (is_same<M, NonRelWilsonMatrixD>()) {
    ret = "NonRelWilsonMatrixD";
  } else if (is_same<M, NonRelWilsonMatrixF>()) {
    ret = "NonRelWilsonMatrixF";
  } else if (is_same<M, IsospinMatrixD>()) {
    ret = "IsospinMatrixD";
  } else if (is_same<M, IsospinMatrixF>()) {
    ret = "IsospinMatrixF";
  } else if (is_same<M, WilsonVectorD>()) {
    ret = "WilsonVectorD";
  } else if (is_same<M, WilsonVectorF>()) {
    ret = "WilsonVectorF";
  } else if (is_same<M, SpinVectorD>()) {
    ret = "SpinVectorD";
  } else if (is_same<M, SpinVectorF>()) {
    ret = "SpinVectorF";
  }
  return ret;
}

// -------------------------------------------------------------------------

template <class M>
qacc constexpr bool is_signed_integer()
{
  bool ret = false;
  if (is_same<M, int8_t>()) {
    ret = true;
  } else if (is_same<M, int16_t>()) {
    ret = true;
  } else if (is_same<M, int32_t>()) {
    ret = true;
  } else if (is_same<M, int64_t>()) {
    ret = true;
  }
  return ret;
}

template <class M>
qacc constexpr bool is_unsigned_integer()
{
  bool ret = false;
  if (is_same<M, uint8_t>()) {
    ret = true;
  } else if (is_same<M, uint16_t>()) {
    ret = true;
  } else if (is_same<M, uint32_t>()) {
    ret = true;
  } else if (is_same<M, uint64_t>()) {
    ret = true;
  }
  return ret;
}

template <class M>
qacc constexpr bool is_integer()
{
  return is_signed_integer<M>() or is_unsigned_integer<M>();
}

template <class M>
qacc constexpr bool is_number()
{
  return is_integer<M>() or is_real<M>() or is_complex<M>();
}

template <class M>
qacc constexpr bool is_char()
{
  bool ret = false;
  if (is_same<M, char>()) {
    ret = true;
  } else if (is_same<M, signed char>()) {
    ret = true;
  } else if (is_same<M, unsigned char>()) {
    ret = true;
  }
  return ret;
}

template <class M>
qacc constexpr bool is_composed_of_complex_d()
{
  bool ret = false;
  if (is_same<M, ComplexD>()) {
    ret = true;
  } else if (is_same<M, ColorMatrixD>()) {
    ret = true;
  } else if (is_same<M, WilsonMatrixD>()) {
    ret = true;
  } else if (is_same<M, SpinMatrixD>()) {
    ret = true;
  } else if (is_same<M, NonRelWilsonMatrixD>()) {
    ret = true;
  } else if (is_same<M, IsospinMatrixD>()) {
    ret = true;
  } else if (is_same<M, WilsonVectorD>()) {
    ret = true;
  } else if (is_same<M, SpinVectorD>()) {
    ret = true;
  }
  return ret;
}

template <class M>
qacc constexpr bool is_composed_of_real_d()
{
  bool ret = is_composed_of_complex_d<M>();
  if (is_same<M, RealD>()) {
    ret = true;
  } else if (is_same<M, CoordinateD>()) {
    ret = true;
  } else if (is_same<M, AdjointColorMatrixD>()) {
    ret = true;
  }
  return ret;
}

template <class M>
qacc constexpr bool is_composed_of_complex_f()
{
  bool ret = false;
  if (is_same<M, ComplexF>()) {
    ret = true;
  } else if (is_same<M, ColorMatrixF>()) {
    ret = true;
  } else if (is_same<M, WilsonMatrixF>()) {
    ret = true;
  } else if (is_same<M, SpinMatrixF>()) {
    ret = true;
  } else if (is_same<M, NonRelWilsonMatrixF>()) {
    ret = true;
  } else if (is_same<M, IsospinMatrixF>()) {
    ret = true;
  } else if (is_same<M, WilsonVectorF>()) {
    ret = true;
  } else if (is_same<M, SpinVectorF>()) {
    ret = true;
  }
  return ret;
}

template <class M>
qacc constexpr bool is_composed_of_real_f()
{
  bool ret = is_composed_of_complex_f<M>();
  if (is_same<M, RealF>()) {
    ret = true;
  } else if (is_same<M, AdjointColorMatrixF>()) {
    ret = true;
  }
  return ret;
}

template <class M>
qacc constexpr bool is_composed_of_long()
{
  bool ret = false;
  if (is_same<M, Long>()) {
    ret = true;
  }
  return ret;
}

template <class M>
qacc constexpr bool is_composed_of_int()
{
  bool ret = false;
  if (is_same<M, Int>()) {
    ret = true;
  } else if (is_same<M, Coordinate>()) {
    ret = true;
  }
  return ret;
}

template <class M>
qacc constexpr bool is_composed_of_char()
// Char is int8_t
// NOTE: not char, not signed char, not unsigned char
{
  bool ret = false;
  if (is_same<M, Char>()) {
    ret = true;
  }
  return ret;
}

template <class M>
qacc constexpr int element_size_of()
// for example: size for convert endianness
{
  int ret = 0;
  if (is_integer<M>() or is_real<M>() or is_char<M>()) {
    ret = sizeof(M);
  } else if (is_composed_of_long<M>()) {
    ret = sizeof(Long);
  } else if (is_composed_of_int<M>()) {
    ret = sizeof(Int);
  } else if (is_composed_of_char<M>()) {
    ret = sizeof(Char);
  } else if (is_composed_of_real_d<M>()) {
    ret = sizeof(RealD);
  } else if (is_composed_of_real_f<M>()) {
    ret = sizeof(RealF);
  } else {
    ret = 0;
  }
  return ret;
}

// -------------------------------------------------------------------------

template <class M>
qacc constexpr bool is_basic_data_type()
// basic data types
// Long, RealD, ComplexD, ColorMatrixD, etc
{
  return element_size_of<M>() > 0;
}

// -------------------------------------------------------------------------

template <class M>
struct IsDataValueType {
  static constexpr bool value = is_basic_data_type<M>();
  using DataType = M;
};

template <class M, size_t N>
struct IsDataValueType<array<M, N>> {
  static constexpr bool value = IsDataValueType<M>::value;
  using DataType = typename IsDataValueType<M>::DataType;
};

template <>
struct IsDataValueType<RngState> {
  static constexpr bool value = true;
  using DataType = RngState;
};

// -------------------------------------------------------------------------

template <class M>
qacc constexpr bool is_data_value_type()
// data value types
// basic data types + qlat::array of basic data types + RngState + etc
// data value types should contain its data within itself (no pointers)
{
  return IsDataValueType<M>::value;
}

// -------------------------------------------------------------------------

#define MAXTYPE 128
#define FLOATIND 30

enum DATA_TYPE {
  CHAR_TYPE = 0 + MAXTYPE * sizeof(char),
  UCHAR_TYPE = 1 + MAXTYPE * sizeof(unsigned char),
  SHORT_TYPE = 2 + MAXTYPE * sizeof(short),
  USHORT_TYPE = 3 + MAXTYPE * sizeof(unsigned short),
  //
  INT_TYPE = 4 + MAXTYPE * sizeof(int),
  UINT_TYPE = 5 + MAXTYPE * sizeof(unsigned int),
  LONG_TYPE = 6 + MAXTYPE * sizeof(Long),
  ULONG_TYPE = 7 + MAXTYPE * sizeof(unsigned long),
  LONGL_TYPE = 8 + MAXTYPE * sizeof(long long),
  ULONGL_TYPE = 9 + MAXTYPE * sizeof(unsigned long long),
  INT8_TYPE = 10 + MAXTYPE * sizeof(std::int8_t),
  UINT8_TYPE = 11 + MAXTYPE * sizeof(std::uint8_t),
  INT16_TYPE = 12 + MAXTYPE * sizeof(std::int16_t),
  UINT16_TYPE = 13 + MAXTYPE * sizeof(std::uint16_t),
  INT32_TYPE = 14 + MAXTYPE * sizeof(std::int32_t),
  UINT32_TYPE = 15 + MAXTYPE * sizeof(std::uint32_t),
  INT64_TYPE = 16 + MAXTYPE * sizeof(std::int64_t),
  UINT64_TYPE = 17 + MAXTYPE * sizeof(std::uint64_t),
  //
  DOUBLE_TYPE = FLOATIND + 0 + MAXTYPE * sizeof(double),
  FLOAT_TYPE = FLOATIND + 1 + MAXTYPE * sizeof(float),
  ComplexD_TYPE = FLOATIND + 2 + MAXTYPE * sizeof(ComplexD),
  ComplexF_TYPE = FLOATIND + 3 + MAXTYPE * sizeof(ComplexF),
  //
  ColorMatrix_TYPE = FLOATIND + 4 + MAXTYPE * sizeof(ColorMatrixT<double>),
  ColorMatrixF_TYPE = FLOATIND + 5 + MAXTYPE * sizeof(ColorMatrixT<float>),
  WilsonMatrix_TYPE = FLOATIND + 6 + MAXTYPE * sizeof(WilsonMatrixT<double>),
  WilsonMatrixF_TYPE = FLOATIND + 7 + MAXTYPE * sizeof(WilsonMatrixT<float>),
  SpinMatrix_TYPE = FLOATIND + 8 + MAXTYPE * sizeof(SpinMatrixT<double>),
  SpinMatrixF_TYPE = FLOATIND + 9 + MAXTYPE * sizeof(SpinMatrixT<float>),
  WilsonVector_TYPE = FLOATIND + 10 + MAXTYPE * sizeof(WilsonVectorT<double>),
  WilsonVectorF_TYPE = FLOATIND + 11 + MAXTYPE * sizeof(WilsonVectorT<float>),
  //
  NonRelWilsonMatrix_TYPE =
      FLOATIND + 12 + MAXTYPE * sizeof(NonRelWilsonMatrixT<double>),
  NonRelWilsonMatrixF_TYPE =
      FLOATIND + 13 + MAXTYPE * sizeof(NonRelWilsonMatrixT<float>),
  INVALID_TYPE = 9999999
};

template <class M>
qacc DATA_TYPE get_data_type()
{
  return INVALID_TYPE;
}
template <>
qacc DATA_TYPE get_data_type<char>()
{
  return CHAR_TYPE;
}
template <>
qacc DATA_TYPE get_data_type<unsigned char>()
{
  return UCHAR_TYPE;
}
template <>
qacc DATA_TYPE get_data_type<short>()
{
  return SHORT_TYPE;
}
template <>
qacc DATA_TYPE get_data_type<unsigned short>()
{
  return USHORT_TYPE;
}
template <>
qacc DATA_TYPE get_data_type<int>()
{
  return INT_TYPE;
}
template <>
qacc DATA_TYPE get_data_type<unsigned int>()
{
  return UINT_TYPE;
}
template <>
qacc DATA_TYPE get_data_type<Long>()
{
  return LONG_TYPE;
}
template <>
qacc DATA_TYPE get_data_type<unsigned long>()
{
  return ULONG_TYPE;
}
template <>
qacc DATA_TYPE get_data_type<std::int8_t>()
{
  return INT8_TYPE;
}
// template<> qacc DATA_TYPE get_data_type<std::int8_t         >(){return
// INT8_TYPE          ; } template<> qacc DATA_TYPE get_data_type<std::uint8_t
// >(){return  UINT8_TYPE          ; } template<> qacc DATA_TYPE
// get_data_type<std::int16_t        >(){return   INT16_TYPE         ; }
// template<> qacc DATA_TYPE get_data_type<std::uint16_t       >(){return
// UINT16_TYPE         ; } template<> qacc DATA_TYPE get_data_type<std::int32_t
// >(){return   INT32_TYPE         ; } template<> qacc DATA_TYPE
// get_data_type<std::uint32_t       >(){return  UINT32_TYPE         ; }
// template<> qacc DATA_TYPE get_data_type<std::int64_t        >(){return
// INT64_TYPE         ; } template<> qacc DATA_TYPE get_data_type<std::uint64_t
// >(){return  UINT64_TYPE         ; }
template <>
qacc DATA_TYPE get_data_type<double>()
{
  return DOUBLE_TYPE;
}
template <>
qacc DATA_TYPE get_data_type<float>()
{
  return FLOAT_TYPE;
}
template <>
qacc DATA_TYPE get_data_type<ComplexD>()
{
  return ComplexD_TYPE;
}
template <>
qacc DATA_TYPE get_data_type<ComplexF>()
{
  return ComplexF_TYPE;
}
template <>
qacc DATA_TYPE get_data_type<ColorMatrixT<double>>()
{
  return ColorMatrix_TYPE;
}
template <>
qacc DATA_TYPE get_data_type<ColorMatrixT<float>>()
{
  return ColorMatrixF_TYPE;
}
template <>
qacc DATA_TYPE get_data_type<WilsonMatrixT<double>>()
{
  return WilsonMatrix_TYPE;
}
template <>
qacc DATA_TYPE get_data_type<WilsonMatrixT<float>>()
{
  return WilsonMatrixF_TYPE;
}
template <>
qacc DATA_TYPE get_data_type<SpinMatrixT<double>>()
{
  return SpinMatrix_TYPE;
}
template <>
qacc DATA_TYPE get_data_type<SpinMatrixT<float>>()
{
  return SpinMatrixF_TYPE;
}
template <>
qacc DATA_TYPE get_data_type<WilsonVectorT<double>>()
{
  return WilsonVector_TYPE;
}
template <>
qacc DATA_TYPE get_data_type<WilsonVectorT<float>>()
{
  return WilsonVectorF_TYPE;
}
template <>
qacc DATA_TYPE get_data_type<NonRelWilsonMatrixT<double>>()
{
  return NonRelWilsonMatrix_TYPE;
}
template <>
qacc DATA_TYPE get_data_type<NonRelWilsonMatrixT<float>>()
{
  return NonRelWilsonMatrixF_TYPE;
}

template <class M>
qacc bool get_data_type_is_double()
{
  DATA_TYPE cur = get_data_type<M>();
  if (cur < FLOATIND or cur == INVALID_TYPE) {
    if (get_id_node() == 0) {
      printf("Given type not float/double %d \n", cur);
    }
    qassert(false);
  }
  if (cur % 2 == 0) {
    return true;
  }
  if (cur % 2 == 1) {
    return false;
  }
  return true;
}

// -------------------------------------------------------------------------

}  // namespace qlat
