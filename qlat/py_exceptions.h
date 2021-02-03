#pragma once

#include <Python.h>
#include <qlat/qlat.h>

#define STRX(x) #x

#define STR(x) STRX(x)

#define pqassert(x)                                              \
  {                                                              \
    if (!(x))                                                    \
      throw std::string("Assert " #x " failed in file " __FILE__ \
                        ":" STR(__LINE__));                      \
  };

#define pqerr(...)                                                    \
  {                                                                   \
    const std::string msg =                                           \
        qlat::ssprintf(__VA_ARGS__) +                                 \
        qlat::ssprintf(" in from '%s' line %d ", __FILE__, __LINE__); \
    throw std::string(msg);                                           \
  };

