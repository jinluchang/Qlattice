#pragma once

#include <qlat-utils/matrix.h>
#include <qlat-utils/endian.h>

#include <cstdint>
#include <vector>

namespace qlat
{  //

typedef std::vector<std::vector<RealD>> DataTable;

typedef uint32_t crc32_t;

// -------------------------------------------------------------------------

template <class M>
qacc constexpr bool is_signed_integer()
{
  return false;
}

template <>
qacc constexpr bool is_signed_integer<int8_t>()
{
  return true;
}

template <>
qacc constexpr bool is_signed_integer<int16_t>()
{
  return true;
}

template <>
qacc constexpr bool is_signed_integer<int32_t>()
{
  return true;
}

template <>
qacc constexpr bool is_signed_integer<int64_t>()
{
  return true;
}

// -------------------------------------------------------------------------

template <class M>
qacc constexpr bool is_unsigned_integer()
{
  return false;
}

template <>
qacc constexpr bool is_unsigned_integer<uint8_t>()
{
  return true;
}

template <>
qacc constexpr bool is_unsigned_integer<uint16_t>()
{
  return true;
}

template <>
qacc constexpr bool is_unsigned_integer<uint32_t>()
{
  return true;
}

template <>
qacc constexpr bool is_unsigned_integer<int64_t>()
{
  return true;
}

// -------------------------------------------------------------------------

template <class M>
qacc constexpr bool is_integer()
{
  return is_signed_integer<M>() or is_unsigned_integer<M>();
}

// -------------------------------------------------------------------------

template <class M>
std::string get_type_name()
{
  return "unknown";
}

template <>
inline std::string get_type_name<int8_t>()
{
  return "Int8t";
}

template <>
inline std::string get_type_name<int16_t>()
{
  return "Int16t";
}

template <>
inline std::string get_type_name<int32_t>()
{
  return "Int32t";
}

template <>
inline std::string get_type_name<int64_t>()
{
  return "Int64t";
}

template <>
inline std::string get_type_name<uint8_t>()
{
  return "UInt8t";
}

template <>
inline std::string get_type_name<uint16_t>()
{
  return "UInt16t";
}

template <>
inline std::string get_type_name<uint32_t>()
{
  return "UInt32t";
}

template <>
inline std::string get_type_name<uint64_t>()
{
  return "UInt64t";
}

template <>
inline std::string get_type_name<RealF>()
{
  return "RealF";
}

template <>
inline std::string get_type_name<RealD>()
{
  return "RealD";
}

template <>
inline std::string get_type_name<ComplexF>()
{
  return "ComplexF";
}

template <>
inline std::string get_type_name<ComplexD>()
{
  return "ComplexD";
}

template <>
inline std::string get_type_name<WilsonVectorD>()
{
  return "WilsonVectorD";
}

template <>
inline std::string get_type_name<SpinMatrixD>()
{
  return "SpinMatrixD";
}

template <>
inline std::string get_type_name<NonRelWilsonMatrixD>()
{
  return "NonRelWilsonMatrixD";
}

template <>
inline std::string get_type_name<WilsonMatrixD>()
{
  return "WilsonMatrixD";
}

template <>
inline std::string get_type_name<ColorMatrixD>()
{
  return "ColorMatrixD";
}

template <>
inline std::string get_type_name<IsospinMatrixD>()
{
  return "IsospinMatrixD";
}

template <>
inline std::string get_type_name<WilsonVectorF>()
{
  return "WilsonVectorF";
}

template <>
inline std::string get_type_name<SpinMatrixF>()
{
  return "SpinMatrixF";
}

template <>
inline std::string get_type_name<NonRelWilsonMatrixF>()
{
  return "NonRelWilsonMatrixF";
}

template <>
inline std::string get_type_name<WilsonMatrixF>()
{
  return "WilsonMatrixF";
}

template <>
inline std::string get_type_name<ColorMatrixF>()
{
  return "ColorMatrixF";
}

template <>
inline std::string get_type_name<IsospinMatrixF>()
{
  return "IsospinMatrixF";
}

// -------------------------------------------------------------------------

template <class M>
qacc constexpr bool is_composed_of_complex_d()
{
  return false;
}

template <>
qacc constexpr bool is_composed_of_complex_d<ComplexD>()
{
  return true;
}

template <>
qacc constexpr bool is_composed_of_complex_d<ColorMatrixD>()
{
  return true;
}

template <>
qacc constexpr bool is_composed_of_complex_d<WilsonMatrixD>()
{
  return true;
}

template <>
qacc constexpr bool is_composed_of_complex_d<SpinMatrixD>()
{
  return true;
}

template <>
qacc constexpr bool is_composed_of_complex_d<NonRelWilsonMatrixD>()
{
  return true;
}

template <>
qacc constexpr bool is_composed_of_complex_d<IsospinMatrixD>()
{
  return true;
}

template <>
qacc constexpr bool is_composed_of_complex_d<AdjointColorMatrixD>()
{
  return true;
}

template <>
qacc constexpr bool is_composed_of_complex_d<WilsonVectorD>()
{
  return true;
}

template <>
qacc constexpr bool is_composed_of_complex_d<SpinVectorD>()
{
  return true;
}

// -------------------------------------------------------------------------

template <class M>
qacc constexpr bool is_composed_of_real_d()
{
  return is_composed_of_complex_d<M>();
}

template <>
qacc constexpr bool is_composed_of_real_d<RealD>()
{
  return true;
}

// -------------------------------------------------------------------------

template <class M>
qacc constexpr bool is_composed_of_complex_f()
{
  return false;
}

template <>
qacc constexpr bool is_composed_of_complex_f<ComplexF>()
{
  return true;
}

template <>
qacc constexpr bool is_composed_of_complex_f<ColorMatrixF>()
{
  return true;
}

template <>
qacc constexpr bool is_composed_of_complex_f<WilsonMatrixF>()
{
  return true;
}

template <>
qacc constexpr bool is_composed_of_complex_f<SpinMatrixF>()
{
  return true;
}

template <>
qacc constexpr bool is_composed_of_complex_f<NonRelWilsonMatrixF>()
{
  return true;
}

template <>
qacc constexpr bool is_composed_of_complex_f<IsospinMatrixF>()
{
  return true;
}

template <>
qacc constexpr bool is_composed_of_complex_f<AdjointColorMatrixF>()
{
  return true;
}

template <>
qacc constexpr bool is_composed_of_complex_f<WilsonVectorF>()
{
  return true;
}

template <>
qacc constexpr bool is_composed_of_complex_f<SpinVectorF>()
{
  return true;
}

// -------------------------------------------------------------------------

template <class M>
qacc constexpr bool is_composed_of_real_f()
{
  return is_composed_of_complex_f<M>();
}

template <>
qacc constexpr bool is_composed_of_real_f<float>()
{
  return true;
}

// -------------------------------------------------------------------------

template <class M>
qacc constexpr bool is_composed_of_int64()
{
  return false;
}

template <>
qacc constexpr bool is_composed_of_int64<int64_t>()
{
  return true;
}

// -------------------------------------------------------------------------

template <class M>
qacc constexpr bool is_composed_of_int32()
{
  return false;
}

template <>
qacc constexpr bool is_composed_of_int32<int32_t>()
{
  return true;
}

// -------------------------------------------------------------------------

template <class M>
qacc constexpr bool is_composed_of_int8()
{
  return false;
}

template <>
qacc constexpr bool is_composed_of_int8<int8_t>()
{
  return true;
}

// -------------------------------------------------------------------------

template <class M>
qacc constexpr int element_size_of()
// for example: size for convert endianness
{
  int ret = 0;
  if (is_integer<M>() or is_real<M>()) {
    ret = sizeof(M);
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
qacc constexpr bool is_data_value_type()  // for example size for convert endian
{
  return element_size_of<M>() > 0;
}

// -------------------------------------------------------------------------

template <class M, QLAT_ENABLE_IF(is_data_value_type<M>())>
qacc void to_from_little_endian(Vector<M> v)
{
  constexpr int size = element_size_of<M>();
  to_from_little_endian<size>((void*)v.data(), v.data_size());
}

template <class M, QLAT_ENABLE_IF(is_data_value_type<M>())>
qacc void to_from_big_endian(Vector<M> v)
{
  constexpr int size = element_size_of<M>();
  to_from_big_endian<size>((void*)v.data(), v.data_size());
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
