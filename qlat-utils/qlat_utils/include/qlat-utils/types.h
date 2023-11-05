#pragma once

#include <qlat-utils/matrix.h>
#include <qlat-utils/endian.h>

#include <cstdint>
#include <vector>

namespace qlat
{  //

typedef std::vector<std::vector<RealD>> DataTable;

typedef uint32_t crc32_t;

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
inline std::string get_type_name<WilsonVector>()
{
  return "WilsonVector";
}

template <>
inline std::string get_type_name<SpinMatrix>()
{
  return "SpinMatrix";
}

template <>
inline std::string get_type_name<NonRelWilsonMatrix>()
{
  return "NonRelWilsonMatrix";
}

template <>
inline std::string get_type_name<WilsonMatrix>()
{
  return "WilsonMatrix";
}

template <>
inline std::string get_type_name<ColorMatrix>()
{
  return "ColorMatrix";
}

template <>
inline std::string get_type_name<IsospinMatrix>()
{
  return "IsospinMatrix";
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
qacc constexpr bool is_composed_of_complex_d<ColorMatrix>()
{
  return true;
}

template <>
qacc constexpr bool is_composed_of_complex_d<WilsonMatrix>()
{
  return true;
}

template <>
qacc constexpr bool is_composed_of_complex_d<SpinMatrix>()
{
  return true;
}

template <>
qacc constexpr bool is_composed_of_complex_d<WilsonVector>()
{
  return true;
}

template <>
qacc constexpr bool is_composed_of_complex_d<IsospinMatrix>()
{
  return true;
}

template <>
qacc constexpr bool is_composed_of_complex_d<NonRelWilsonMatrix>()
{
  return true;
}

// -------------------------------------------------------------------------

template <class M>
qacc constexpr bool is_composed_of_real_d()
{
  return false;
}

template <>
qacc constexpr bool is_composed_of_real_d<RealD>()
{
  return true;
}

template <>
qacc constexpr bool is_composed_of_real_d<ComplexD>()
{
  return true;
}

template <>
qacc constexpr bool is_composed_of_real_d<ColorMatrix>()
{
  return true;
}

template <>
qacc constexpr bool is_composed_of_real_d<WilsonMatrix>()
{
  return true;
}

template <>
qacc constexpr bool is_composed_of_real_d<SpinMatrix>()
{
  return true;
}

template <>
qacc constexpr bool is_composed_of_real_d<WilsonVector>()
{
  return true;
}

template <>
qacc constexpr bool is_composed_of_real_d<IsospinMatrix>()
{
  return true;
}

template <>
qacc constexpr bool is_composed_of_real_d<NonRelWilsonMatrix>()
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

// -------------------------------------------------------------------------

template <class M>
qacc constexpr bool is_composed_of_real_f()
{
  return false;
}

template <>
qacc constexpr bool is_composed_of_real_f<float>()
{
  return true;
}

template <>
qacc constexpr bool is_composed_of_real_f<ComplexF>()
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
qacc constexpr int element_size_of();  // for example size for convert endian

template <>
qacc constexpr int element_size_of<int8_t>()
{
  return 1;
}

template <>
qacc constexpr int element_size_of<int16_t>()
{
  return 2;
}

template <>
qacc constexpr int element_size_of<int32_t>()
{
  return 4;
}

template <>
qacc constexpr int element_size_of<int64_t>()
{
  return 8;
}

template <>
qacc constexpr int element_size_of<uint8_t>()
{
  return 1;
}

template <>
qacc constexpr int element_size_of<uint16_t>()
{
  return 2;
}

template <>
qacc constexpr int element_size_of<uint32_t>()
{
  return 4;
}

template <>
qacc constexpr int element_size_of<uint64_t>()
{
  return 8;
}

template <>
qacc constexpr int element_size_of<RealF>()
{
  return 4;
}

template <>
qacc constexpr int element_size_of<RealD>()
{
  return 8;
}

template <>
qacc constexpr int element_size_of<ComplexF>()
{
  return 4;
}

template <>
qacc constexpr int element_size_of<ComplexD>()
{
  return 8;
}

template <>
qacc constexpr int element_size_of<WilsonVector>()
{
  return 8;
}

template <>
qacc constexpr int element_size_of<SpinMatrix>()
{
  return 8;
}

template <>
qacc constexpr int element_size_of<NonRelWilsonMatrix>()
{
  return 8;
}

template <>
qacc constexpr int element_size_of<WilsonMatrix>()
{
  return 8;
}

template <>
qacc constexpr int element_size_of<ColorMatrix>()
{
  return 8;
}

template <>
qacc constexpr int element_size_of<IsospinMatrix>()
{
  return 8;
}

// -------------------------------------------------------------------------

template <class M>
qacc void to_from_little_endian(Vector<M> v)
{
  constexpr int size = element_size_of<M>();
  to_from_little_endian<size>((void*)v.data(), v.data_size());
}

template <class M>
qacc void to_from_big_endian(Vector<M> v)
{
  constexpr int size = element_size_of<M>();
  to_from_big_endian<size>((void*)v.data(), v.data_size());
}

// -------------------------------------------------------------------------

}  // namespace qlat
