#pragma once

#define field_dispatch(p_ret, fname, ctype, ...)                            \
  {                                                                         \
    if ("ColorMatrix" == (ctype)) {                                         \
      (p_ret) = fname##_##ctype<ColorMatrix>(__VA_ARGS__);                  \
    } else if ("WilsonMatrix" == (ctype)) {                                 \
      (p_ret) = fname##_##ctype<WilsonMatrix>(__VA_ARGS__);                 \
    } else if ("WilsonVector" == (ctype)) {                                 \
      (p_ret) = fname##_##ctype<WilsonVector>(__VA_ARGS__);                 \
    } else if ("double" == (ctype)) {                                       \
      (p_ret) = fname##_##ctype<double>(__VA_ARGS__);                       \
    } else if ("Complex" == (ctype)) {                                      \
      (p_ret) = fname##_##ctype<Complex>(__VA_ARGS__);                      \
    } else {                                                                \
      pqerr("%s %s='%s' does not exist.", #fname, #ctype, (ctype).c_str()); \
      (p_ret) = NULL;                                                       \
    }                                                                       \
  }
