// vim: set ts=2 sw=2 expandtab:

#pragma once

#include <cstring>

#include <endian.h>

#include <hash-cpp/crc32.h>
#include <show.h>

QLAT_START_NAMESPACE

typedef uint32_t crc32_t;

inline crc32_t crc32(const crc32_t initial, const void* smessage, long nBytes) {
  qassert(sizeof(CRC32) == sizeof(crc32_t));
  CRC32 c;
  std::memcpy((void*)&c, (void*)&initial, sizeof(crc32_t));
  c.add(smessage, nBytes);
  crc32_t ret;
  std::memcpy((void*)&ret, (void*)&c, sizeof(crc32_t));
  return ret;
}

inline crc32_t crc32(const void* smessage, int nBytes) {
  return crc32(0, smessage, nBytes);
}

inline void crc32_check() {
  const char* test = "123456789";
  const crc32_t CHECK_VALUE = 0xCBF43926;
  displayln_info(ssprintf("The check value for the %s standard is 0x%X", "CRC32", 0xCBF43926));
  displayln_info(ssprintf("The crc32() of \"123456789\" is 0x%X", crc32(test, std::strlen(test))));
  qassert(CHECK_VALUE == crc32(test, std::strlen(test)));
}

QLAT_END_NAMESPACE

namespace qshow {

inline std::string show(const qlat::crc32_t x) {
  return ssprintf("%08X", x);
}

}

#ifndef USE_NAMESPACE
using namespace qshow;
#endif
