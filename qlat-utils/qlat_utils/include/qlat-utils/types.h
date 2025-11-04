#pragma once

#include <qlat-utils/mat-vec.h>
#include <qlat-utils/rng-state.h>

#include <cstdint>
#include <vector>
#include <array>

namespace qlat
{  //

using DataTable = std::vector<std::vector<RealD>>;

using crc32_t = uint32_t;

// -------------------------------------------------------------------------

template <class M, std::size_t N>
qacc bool qisnan(const array<M, N>& arr)
{
  for (Int i = 0; i < (Int)N; ++i) {
    if (qisnan(arr[i])) {
      return true;
    }
  }
  return false;
}

template <class M>
bool qisnan(const std::vector<M>& arr)
{
  for (std::size_t i = 0; i < arr.size(); ++i) {
    if (qisnan(arr[i])) {
      return true;
    }
  }
  return false;
}

// --------------------

struct API Coordinate : public array<Int, DIMN> {
  qacc Coordinate() { init(); }
  qacc Coordinate(Int first, Int second, Int third, Int fourth)
  {
    Int* p = data();
    p[0] = first;
    p[1] = second;
    p[2] = third;
    p[3] = fourth;
  }
  qacc void init() { array<Int, DIMN>::fill(0); }
};

struct API CoordinateD : public array<RealD, DIMN> {
  qacc CoordinateD() { init(); }
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
    for (Int i = 0; i < DIMN; ++i) {
      c[i] = x[i];
    }
  }
  qacc void init() { array<RealD, DIMN>::fill(0); }
};

// -------------------------------------------------------------------------

template <class M>
struct IsBasicDataType {
  static constexpr bool value = false;
  static constexpr bool is_complex = false;
  static const std::string get_type_name() { return "unknown_type"; }
  using ElementaryType = M;
};

template <>
struct IsBasicDataType<char> {
  static constexpr bool value = true;
  static constexpr bool is_complex = false;
  static const std::string get_type_name() { return "Char"; }
  using ElementaryType = Char;
};

template <>
struct IsBasicDataType<int8_t> {
  static constexpr bool value = true;
  static constexpr bool is_complex = false;
  static const std::string get_type_name() { return "Int8t"; }
  using ElementaryType = int8_t;
};

template <>
struct IsBasicDataType<int16_t> {
  static constexpr bool value = true;
  static constexpr bool is_complex = false;
  static const std::string get_type_name() { return "Int16t"; }
  using ElementaryType = int16_t;
};

template <>
struct IsBasicDataType<int32_t> {
  static constexpr bool value = true;
  static constexpr bool is_complex = false;
  static const std::string get_type_name() { return "Int32t"; }
  using ElementaryType = int32_t;
};

template <>
struct IsBasicDataType<int64_t> {
  static constexpr bool value = true;
  static constexpr bool is_complex = false;
  static const std::string get_type_name() { return "Int64t"; }
  using ElementaryType = int64_t;
};

template <>
struct IsBasicDataType<uint8_t> {
  static constexpr bool value = true;
  static constexpr bool is_complex = false;
  static const std::string get_type_name() { return "UInt8t"; }
  using ElementaryType = uint8_t;
};

template <>
struct IsBasicDataType<uint16_t> {
  static constexpr bool value = true;
  static constexpr bool is_complex = false;
  static const std::string get_type_name() { return "UInt16t"; }
  using ElementaryType = uint16_t;
};

template <>
struct IsBasicDataType<uint32_t> {
  static constexpr bool value = true;
  static constexpr bool is_complex = false;
  static const std::string get_type_name() { return "UInt32t"; }
  using ElementaryType = uint32_t;
};

template <>
struct IsBasicDataType<uint64_t> {
  static constexpr bool value = true;
  static constexpr bool is_complex = false;
  static const std::string get_type_name() { return "UInt64t"; }
  using ElementaryType = uint64_t;
};

template <>
struct IsBasicDataType<RealF> {
  static constexpr bool value = true;
  static constexpr bool is_complex = false;
  static const std::string get_type_name() { return "RealF"; }
  using ElementaryType = RealF;
};

template <>
struct IsBasicDataType<RealD> {
  static constexpr bool value = true;
  static constexpr bool is_complex = false;
  static const std::string get_type_name() { return "RealD"; }
  using ElementaryType = RealD;
};

template <>
struct IsBasicDataType<RealDD> {
  static constexpr bool value = true;
  static constexpr bool is_complex = false;
  static const std::string get_type_name() { return "RealDD"; }
  using ElementaryType = RealDD;
};

template <>
struct IsBasicDataType<ComplexF> {
  static constexpr bool value = true;
  static constexpr bool is_complex = true;
  static const std::string get_type_name() { return "ComplexF"; }
  using ElementaryType = RealF;
};

template <>
struct IsBasicDataType<ComplexD> {
  static constexpr bool value = true;
  static constexpr bool is_complex = true;
  static const std::string get_type_name() { return "ComplexD"; }
  using ElementaryType = RealD;
};

template <>
struct IsBasicDataType<ComplexT<RealDD>> {
  static constexpr bool value = true;
  static constexpr bool is_complex = true;
  static const std::string get_type_name() { return "ComplexDD"; }
  using ElementaryType = RealDD;
};

template <Int DIMN, class T>
struct IsBasicDataType<MatrixT<DIMN, T>> {
  static constexpr bool value = IsBasicDataType<T>::value;
  static constexpr bool is_complex = IsBasicDataType<T>::is_complex;
  static const std::string get_type_name() { return "MatrixT_" + IsBasicDataType<T>::get_type_name(); }
  using ElementaryType = typename IsBasicDataType<T>::ElementaryType;
};

template <Int DIMN, class T>
struct IsBasicDataType<MvectorT<DIMN, T>> {
  static constexpr bool value = IsBasicDataType<T>::value;
  static constexpr bool is_complex = IsBasicDataType<T>::is_complex;
  static const std::string get_type_name() { return "MvectorT_" + IsBasicDataType<T>::get_type_name(); }
  using ElementaryType = typename IsBasicDataType<T>::ElementaryType;
};

template <>
struct IsBasicDataType<ColorMatrixT<RealDD>> {
  static constexpr bool value = true;
  static constexpr bool is_complex = true;
  static const std::string get_type_name() { return "ColorMatrixDD"; }
  using ElementaryType = RealDD;
};

template <>
struct IsBasicDataType<ColorMatrixD> {
  static constexpr bool value = true;
  static constexpr bool is_complex = true;
  static const std::string get_type_name() { return "ColorMatrixD"; }
  using ElementaryType = RealD;
};

template <>
struct IsBasicDataType<ColorMatrixF> {
  static constexpr bool value = true;
  static constexpr bool is_complex = true;
  static const std::string get_type_name() { return "ColorMatrixF"; }
  using ElementaryType = RealF;
};

template <>
struct IsBasicDataType<WilsonMatrixD> {
  static constexpr bool value = true;
  static constexpr bool is_complex = true;
  static const std::string get_type_name() { return "WilsonMatrixD"; }
  using ElementaryType = RealD;
};

template <>
struct IsBasicDataType<WilsonMatrixF> {
  static constexpr bool value = true;
  static constexpr bool is_complex = true;
  static const std::string get_type_name() { return "WilsonMatrixF"; }
  using ElementaryType = RealF;
};

template <>
struct IsBasicDataType<SpinMatrixD> {
  static constexpr bool value = true;
  static constexpr bool is_complex = true;
  static const std::string get_type_name() { return "SpinMatrixD"; }
  using ElementaryType = RealD;
};

template <>
struct IsBasicDataType<SpinMatrixF> {
  static constexpr bool value = true;
  static constexpr bool is_complex = true;
  static const std::string get_type_name() { return "SpinMatrixF"; }
  using ElementaryType = RealF;
};

template <>
struct IsBasicDataType<NonRelWilsonMatrixD> {
  static constexpr bool value = true;
  static constexpr bool is_complex = true;
  static const std::string get_type_name() { return "NonRelWilsonMatrixD"; }
  using ElementaryType = RealD;
};

template <>
struct IsBasicDataType<NonRelWilsonMatrixF> {
  static constexpr bool value = true;
  static constexpr bool is_complex = true;
  static const std::string get_type_name() { return "NonRelWilsonMatrixF"; }
  using ElementaryType = RealF;
};

template <>
struct IsBasicDataType<IsospinMatrixD> {
  static constexpr bool value = true;
  static constexpr bool is_complex = true;
  static const std::string get_type_name() { return "IsospinMatrixD"; }
  using ElementaryType = RealD;
};

template <>
struct IsBasicDataType<IsospinMatrixF> {
  static constexpr bool value = true;
  static constexpr bool is_complex = true;
  static const std::string get_type_name() { return "IsospinMatrixF"; }
  using ElementaryType = RealF;
};

template <>
struct IsBasicDataType<AdjointColorMatrixD> {
  static constexpr bool value = true;
  static constexpr bool is_complex = false;
  static const std::string get_type_name() { return "AdjointColorMatrixD"; }
  using ElementaryType = RealD;
};

template <>
struct IsBasicDataType<AdjointColorMatrixF> {
  static constexpr bool value = true;
  static constexpr bool is_complex = false;
  static const std::string get_type_name() { return "AdjointColorMatrixF"; }
  using ElementaryType = RealF;
};

template <>
struct IsBasicDataType<WilsonVectorD> {
  static constexpr bool value = true;
  static constexpr bool is_complex = true;
  static const std::string get_type_name() { return "WilsonVectorD"; }
  using ElementaryType = RealD;
};

template <>
struct IsBasicDataType<WilsonVectorF> {
  static constexpr bool value = true;
  static constexpr bool is_complex = true;
  static const std::string get_type_name() { return "WilsonVectorF"; }
  using ElementaryType = RealF;
};

template <>
struct IsBasicDataType<SpinVectorD> {
  static constexpr bool value = true;
  static constexpr bool is_complex = true;
  static const std::string get_type_name() { return "SpinVectorD"; }
  using ElementaryType = RealD;
};

template <>
struct IsBasicDataType<SpinVectorF> {
  static constexpr bool value = true;
  static constexpr bool is_complex = true;
  static const std::string get_type_name() { return "SpinVectorF"; }
  using ElementaryType = RealF;
};

template <>
struct IsBasicDataType<Coordinate> {
  static constexpr bool value = true;
  static constexpr bool is_complex = false;
  static const std::string get_type_name() { return "Coordinate"; }
  using ElementaryType = Int;
};

template <>
struct IsBasicDataType<CoordinateD> {
  static constexpr bool value = true;
  static constexpr bool is_complex = false;
  static const std::string get_type_name() { return "CoordinateD"; }
  using ElementaryType = RealD;
};

template <>
struct IsBasicDataType<RngState> {
  static constexpr bool value = true;
  static constexpr bool is_complex = false;
  static const std::string get_type_name() { return "RngState"; }
  using ElementaryType = Char;
};

// -------------------------------------------------------------------------

template <class M>
qacc constexpr bool is_basic_data_type()
// basic data types
// Long, RealD, ComplexD, ColorMatrixD, etc
{
  return IsBasicDataType<M>::value;
}

template <class M>
std::string get_type_name()
{
  return IsBasicDataType<M>::get_type_name();
}

template <class M>
qacc constexpr bool is_elementary_type()
// elementary types
// Long, RealD, etc
{
  return is_same<typename IsBasicDataType<M>::ElementaryType, M>();
}

// -------------------------------------------------------------------------

template <class M>
struct IsDataValueType {
  using DataType = M;
  using BasicDataType = DataType;
  using ElementaryType =
      typename IsBasicDataType<BasicDataType>::ElementaryType;
  static constexpr bool value = is_basic_data_type<BasicDataType>();
  static constexpr bool is_complex = IsBasicDataType<BasicDataType>::is_complex;
};

template <class M, std::size_t N>
struct IsDataValueType<array<M, N>> {
  using DataType = M;
  using BasicDataType = typename IsDataValueType<DataType>::BasicDataType;
  using ElementaryType = typename IsBasicDataType<DataType>::ElementaryType;
  static constexpr bool value = IsDataValueType<DataType>::value;
  static constexpr bool is_complex = IsDataValueType<DataType>::is_complex;
};

template <class M, std::size_t N>
struct IsDataValueType<std::array<M, N>> {
  using DataType = M;
  using BasicDataType = typename IsDataValueType<DataType>::BasicDataType;
  using ElementaryType = typename IsBasicDataType<DataType>::ElementaryType;
  static constexpr bool value = IsDataValueType<DataType>::value;
  static constexpr bool is_complex = IsDataValueType<DataType>::is_complex;
};

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
  if (is_same<M, Char>()) {
    ret = true;
  }
  return ret;
}

// -------------------------------------------------------------------------

template <class M>
qacc constexpr bool is_composed_of_real()
{
  return is_real<typename IsDataValueType<M>::ElementaryType>();
}

template <class M>
qacc constexpr bool is_composed_of_complex()
{
  return IsDataValueType<M>::is_complex and is_composed_of_real<M>();
}

template <class M>
qacc constexpr bool is_composed_of_integer()
{
  return is_integer<typename IsDataValueType<M>::ElementaryType>();
}

// -------------------------------------------------------------------------

template <class M>
qacc constexpr bool is_composed_of_real_d()
{
  return is_same<typename IsDataValueType<M>::ElementaryType, RealD>();
}

template <class M>
qacc constexpr bool is_composed_of_real_f()
{
  return is_same<typename IsDataValueType<M>::ElementaryType, RealF>();
}

template <class M>
qacc constexpr bool is_composed_of_complex_d()
{
  return IsDataValueType<M>::is_complex and is_composed_of_real_d<M>();
}

template <class M>
qacc constexpr bool is_composed_of_complex_f()
{
  return IsDataValueType<M>::is_complex and is_composed_of_real_f<M>();
}

template <class M>
qacc constexpr bool is_composed_of_long()
{
  return is_same<typename IsDataValueType<M>::ElementaryType, Long>();
}

template <class M>
qacc constexpr bool is_composed_of_int()
{
  return is_same<typename IsDataValueType<M>::ElementaryType, Int>();
}

template <class M>
qacc constexpr bool is_composed_of_char()
{
  return is_same<typename IsDataValueType<M>::ElementaryType, Char>();
}

// -------------------------------------------------------------------------

template <class M>
qacc constexpr Int element_size_of()
// for example: size for convert endianness
{
  return sizeof(typename IsDataValueType<M>::ElementaryType);
}

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
  LONG_TYPE = 6 + MAXTYPE * sizeof(long),
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
  DOUBLE_TYPE = FLOATIND + 0 + MAXTYPE * sizeof(RealD),
  FLOAT_TYPE = FLOATIND + 1 + MAXTYPE * sizeof(float),
  ComplexD_TYPE = FLOATIND + 2 + MAXTYPE * sizeof(ComplexD),
  ComplexF_TYPE = FLOATIND + 3 + MAXTYPE * sizeof(ComplexF),
  //
  ColorMatrix_TYPE = FLOATIND + 4 + MAXTYPE * sizeof(ColorMatrixT<RealD>),
  ColorMatrixF_TYPE = FLOATIND + 5 + MAXTYPE * sizeof(ColorMatrixT<float>),
  WilsonMatrix_TYPE = FLOATIND + 6 + MAXTYPE * sizeof(WilsonMatrixT<RealD>),
  WilsonMatrixF_TYPE = FLOATIND + 7 + MAXTYPE * sizeof(WilsonMatrixT<float>),
  SpinMatrix_TYPE = FLOATIND + 8 + MAXTYPE * sizeof(SpinMatrixT<RealD>),
  SpinMatrixF_TYPE = FLOATIND + 9 + MAXTYPE * sizeof(SpinMatrixT<float>),
  WilsonVector_TYPE = FLOATIND + 10 + MAXTYPE * sizeof(WilsonVectorT<RealD>),
  WilsonVectorF_TYPE = FLOATIND + 11 + MAXTYPE * sizeof(WilsonVectorT<float>),
  //
  NonRelWilsonMatrix_TYPE =
      FLOATIND + 12 + MAXTYPE * sizeof(NonRelWilsonMatrixT<RealD>),
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
qacc DATA_TYPE get_data_type<Int>()
{
  return INT_TYPE;
}
template <>
qacc DATA_TYPE get_data_type<uint32_t>()
{
  return UINT_TYPE;
}
template <>
qacc DATA_TYPE get_data_type<Long>()
{
  return LONG_TYPE;
}
template <>
qacc DATA_TYPE get_data_type<uint64_t>()
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
qacc DATA_TYPE get_data_type<RealD>()
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
qacc DATA_TYPE get_data_type<ColorMatrixT<RealD>>()
{
  return ColorMatrix_TYPE;
}
template <>
qacc DATA_TYPE get_data_type<ColorMatrixT<float>>()
{
  return ColorMatrixF_TYPE;
}
template <>
qacc DATA_TYPE get_data_type<WilsonMatrixT<RealD>>()
{
  return WilsonMatrix_TYPE;
}
template <>
qacc DATA_TYPE get_data_type<WilsonMatrixT<float>>()
{
  return WilsonMatrixF_TYPE;
}
template <>
qacc DATA_TYPE get_data_type<SpinMatrixT<RealD>>()
{
  return SpinMatrix_TYPE;
}
template <>
qacc DATA_TYPE get_data_type<SpinMatrixT<float>>()
{
  return SpinMatrixF_TYPE;
}
template <>
qacc DATA_TYPE get_data_type<WilsonVectorT<RealD>>()
{
  return WilsonVector_TYPE;
}
template <>
qacc DATA_TYPE get_data_type<WilsonVectorT<float>>()
{
  return WilsonVectorF_TYPE;
}
template <>
qacc DATA_TYPE get_data_type<NonRelWilsonMatrixT<RealD>>()
{
  return NonRelWilsonMatrix_TYPE;
}
template <>
qacc DATA_TYPE get_data_type<NonRelWilsonMatrixT<float>>()
{
  return NonRelWilsonMatrixF_TYPE;
}

template <class M>
qacc Int get_data_type_is_Double()
{
  return -1;
}

template <>
qacc Int get_data_type_is_Double<RealF>()
{
  return  0;
}

template <>
qacc Int get_data_type_is_Double<RealD>()
{
  return  1;
}

template <>
qacc Int get_data_type_is_Double<RealDD>()
{
  return  2;
}

template <class M>
qacc Int Is_data_double()
{
  using D = typename IsBasicDataType<M>::ElementaryType;
  Int type = get_data_type_is_Double<D>();
  qassert(type != -1);
  return type;
}

template <class M>
qacc bool get_data_type_is_double()
{
  DATA_TYPE cur = get_data_type<M>();
  if ((cur < FLOATIND) or (cur == INVALID_TYPE)) {
    if (get_id_node() == 0) {
      printf("Given type not float/double %d \n", (Int)cur);
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
