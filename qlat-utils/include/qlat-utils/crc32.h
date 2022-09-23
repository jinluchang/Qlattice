// vim: set ts=2 sw=2 expandtab:

#pragma once

#include <cstring>

#include <zlib.h>

#include <qlat-utils/show.h>
#include <qlat-utils/qutils.h>
#include <qlat-utils/qutils-vec.h>

namespace qlat
{  //

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

inline crc32_t crc32(const void* smessage, const long nBytes)
{
  return crc32(0, smessage, nBytes);
}

inline crc32_t crc32_shift(const crc32_t initial, const long offset)
// shift initial left by offset length
// if offset == 0 then return initial
// offset should be the length of the part after the initial part (which
// evaluate to crc initial) xor all the results gives the final crc32
{
  return crc32_combine(initial, 0, offset);
}

template <class M>
crc32_t crc32(const crc32_t initial, const Vector<M> v)
{
  return crc32(initial, v.data(), v.data_size());
}

template <class M>
crc32_t crc32(const Vector<M> v)
{
  return crc32(0, v.data(), v.data_size());
}

inline crc32_t crc32_par(const void* smessage, const long nBytes)
{
  const uint8_t* ptr = (const uint8_t*)smessage;
  const int v_limit = omp_get_max_threads();
  std::vector<crc32_t> crcs(v_limit, 0);
  crc32_t ret = 0;
#pragma omp parallel
  {
    const int nthreads = omp_get_num_threads();
    const int id = omp_get_thread_num();
    long start = 0, len = 0;
    split_work(start, len, nBytes, nthreads, id);
    const crc32_t crc = crc32(ptr + start, len);
    crcs[id] = crc32_shift(crc, nBytes - start - len);
#pragma omp barrier
    if (0 == id) {
      for (int i = 0; i < nthreads; ++i) {
        ret ^= crcs[i];
      }
    }
  }
  return ret;
}

template <class M>
crc32_t crc32_par(const Vector<M> v)
{
  return crc32_par(v.data(), v.data_size());
}

template <class M>
crc32_t crc32_par(const crc32_t initial, const Vector<M> v)
{
  return crc32_shift(initial, v.data_size()) ^ crc32_par(v);
}

inline void crc32_check()
{
  const char* test = "123456789";
  const crc32_t check_value = 0xCBF43926;
  displayln(ssprintf("The check value for the %s standard is 0x%X", "CRC32",
                     check_value));
  const crc32_t v1 = crc32(test, std::strlen(test));
  displayln(ssprintf("The crc32() of \"123456789\" is 0x%X", v1));
  const crc32_t v2 = crc32(crc32(test, 3), test + 3, 6);
  displayln(ssprintf("The crc32() of \"123456789\" is 0x%X (concat)", v2));
  const crc32_t v3 = crc32_shift(crc32(test, 3), 6) ^ crc32(test + 3, 6);
  displayln(ssprintf("The crc32() of \"123456789\" is 0x%X (crc32_shift)", v3));
  qassert(check_value == v1);
  qassert(check_value == v2);
  qassert(check_value == v3);
  const int nthreads = omp_get_max_threads();
  displayln(ssprintf("nthreads=%d", nthreads));
  const long limit = 1024 * 1024;
  std::vector<uint8_t> test_data(limit, 0);
  for (long i = 0; i < limit; ++i) {
    test_data[i] = i * 7 + 23;
  }
  for (int i = 0; i < 128; ++i) {
    const Vector<uint8_t> v(&test_data[0], i);
    qassert(v.data_size() <= limit);
    if (crc32_par(check_value, v) != crc32(check_value, v)) {
      displayln(ssprintf("i=%d, v.data_size()=%ld", i, v.data_size()));
      qassert(false);
    }
  }
  for (int i = 0; i < 128; ++i) {
    const Vector<double> v((const double*)&test_data[0], i);
    qassert(v.data_size() <= limit);
    if (crc32_par(check_value, v) != crc32(check_value, v)) {
      displayln(ssprintf("i=%d, v.data_size()=%ld", i, v.data_size()));
      qassert(false);
    }
  }
  for (int i = 0; i < 4; ++i) {
    const Vector<uint8_t> v((const uint8_t*)&test_data[0], i * 16 * 1024 + 37);
    qassert(v.data_size() <= limit);
    if (crc32_par(check_value, v) != crc32(check_value, v)) {
      displayln(ssprintf("i=%d, v.data_size()=%ld", i, v.data_size()));
      qassert(false);
    }
  }
}

inline crc32_t read_crc32(const std::string& s)
{
  crc32_t crc32;
  std::sscanf(s.c_str(), "%X", &crc32);
  return crc32;
}

}  // namespace qlat
