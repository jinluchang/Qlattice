// vim: set ts=2 sw=2 expandtab:

#pragma once

#include <cstring>

#include <endian.h>

#include <zlib.h>
#include <show.h>

QLAT_START_NAMESPACE

typedef uint32_t crc32_t;

// inline crc32_t crc32(const crc32_t initial, const void* data, const long len)
// {
//   // crc32 of zlib was incorrect for very large sizes, so do it block-wise
//   const size_t step = 1024*1024*1024;
//   const unsigned char* blk = (const unsigned char*)data;
//   long remain_len = len;
//   crc32_t crc = initial;
//   while (remain_len > step) {
//     crc = crc32_z(crc, blk, step);
//     blk += step;
//     remain_len -= step;
//   }
//   crc = crc32_z(crc, blk, remain_len);
//   return crc;
// }

inline crc32_t crc32(const crc32_t initial, const void* data, const long len)
{
  return crc32_z(initial, (const unsigned char*)data, len);
}

inline crc32_t crc32(const void* smessage, long nBytes)
{
  return crc32(0, smessage, nBytes);
}

template <class M>
inline crc32_t crc32(const crc32_t initial, const Vector<M> v)
{
  return crc32(initial, v.data(), v.data_size());
}

template <class M>
inline crc32_t crc32(const Vector<M> v)
{
  return crc32(0, v.data(), v.data_size());
}

inline void crc32_check()
{
  const char* test = "123456789";
  const crc32_t CHECK_VALUE = 0xCBF43926;
  displayln_info(ssprintf("The check value for the %s standard is 0x%X", "CRC32", 0xCBF43926));
  displayln_info(ssprintf("The crc32() of \"123456789\" is 0x%X", crc32(test, std::strlen(test))));
  qassert(CHECK_VALUE == crc32(test, std::strlen(test)));
}

QLAT_END_NAMESPACE

#ifndef USE_NAMESPACE
using namespace qshow;
#endif
