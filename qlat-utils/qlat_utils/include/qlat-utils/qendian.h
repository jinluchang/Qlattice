#pragma once

#ifdef QLAT_USE_MACHINE_ENDIAN_H
#include <machine/endian.h>
#else
#include <endian.h>
#endif

#include <qlat-utils/handle.h>
#include <qlat-utils/qacc-func.h>
#include <qlat-utils/types.h>

namespace qlat
{  //

// -------------------

qacc bool is_big_endian()
{
#if defined(__BYTE_ORDER) && (__BYTE_ORDER != 0) && \
    (__BYTE_ORDER == __BIG_ENDIAN)
  return true;
#else
  return false;
#endif
}

qacc bool is_little_endian() { return not is_big_endian(); }

// -------------------

qacc uint16_t flip_endian_16(uint16_t x) { return ((x >> 8)) | ((x << 8)); }

qacc uint32_t flip_endian_32(uint32_t x)
{
  return ((x >> 24)) | ((x >> 8) & 0x0000FF00) | ((x << 8) & 0x00FF0000) |
         ((x << 24));
}

qacc uint64_t flip_endian_64(uint64_t x)
{
  return ((x >> 56)) | ((x >> 40) & 0xFF00) | ((x >> 24) & 0xFF0000) |
         ((x >> 8) & 0xFF000000) | ((x << 8) & 0xFF00000000) |
         ((x << 24) & 0xFF0000000000) | ((x << 40) & 0xFF000000000000) |
         ((x << 56));
}

qacc void flip_endian_16(void* str, const size_t len, const bool is_par)
{
  qassert(0 == len % 2);
  uint16_t* p = (uint16_t*)str;
  if (is_par) {
    qthread_for(i, (Long)(len / 2), { p[i] = flip_endian_16(p[i]); });
  } else {
    qfor(i, (Long)(len / 2), { p[i] = flip_endian_16(p[i]); });
  }
}

qacc void flip_endian_32(void* str, const size_t len, const bool is_par)
{
  qassert(0 == len % 4);
  uint32_t* p = (uint32_t*)str;
  if (is_par) {
    qthread_for(i, (Long)(len / 4), { p[i] = flip_endian_32(p[i]); });
  } else {
    qfor(i, (Long)(len / 4), { p[i] = flip_endian_32(p[i]); });
  }
}

qacc void flip_endian_64(void* str, const size_t len, const bool is_par)
{
  qassert(0 == len % 8);
  uint64_t* p = (uint64_t*)str;
  if (is_par) {
    qthread_for(i, (Long)(len / 8), { p[i] = flip_endian_64(p[i]); });
  } else {
    qfor(i, (Long)(len / 8), { p[i] = flip_endian_64(p[i]); });
  }
}

template <Int SIZE>
qacc void flip_endian(void* str, const size_t len, const bool is_par);

template <>
qacc void flip_endian<1>(void* str, const size_t len, const bool is_par)
{
  (void)str;
  (void)len;
  (void)is_par;
}

template <>
qacc void flip_endian<2>(void* str, const size_t len, const bool is_par)
{
  flip_endian_16(str, len, is_par);
}

template <>
qacc void flip_endian<4>(void* str, const size_t len, const bool is_par)
{
  flip_endian_32(str, len, is_par);
}

template <>
qacc void flip_endian<8>(void* str, const size_t len, const bool is_par)
{
  flip_endian_64(str, len, is_par);
}

// -------------------

template <Int SIZE>
qacc void to_from_little_endian(void* str, const size_t len, const bool is_par)
{
  qassert(0 == len % SIZE);
  if (is_big_endian()) {
    flip_endian<SIZE>(str, len, is_par);
  }
}

template <Int SIZE>
qacc void to_from_big_endian(void* str, const size_t len, const bool is_par)
{
  qassert(0 == len % SIZE);
  if (is_little_endian()) {
    flip_endian<SIZE>(str, len, is_par);
  }
}

// -------------------

template <class M, QLAT_ENABLE_IF(is_data_value_type<M>())>
qacc void to_from_little_endian(Vector<M> v, const bool is_par = false)
{
  constexpr Int size = element_size_of<M>();
  to_from_little_endian<size>((void*)v.data(), v.data_size(), is_par);
}

template <class M, QLAT_ENABLE_IF(is_data_value_type<M>())>
qacc void to_from_big_endian(Vector<M> v, const bool is_par = false)
{
  constexpr Int size = element_size_of<M>();
  to_from_big_endian<size>((void*)v.data(), v.data_size(), is_par);
}

// -------------------

}  // namespace qlat
