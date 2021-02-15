#pragma once

#define FIELD_DISPATCH(p_ret, fname, ctype, ...)                            \
  {                                                                         \
    if ("ColorMatrix" == (ctype)) {                                         \
      (p_ret) = fname<ColorMatrix>(__VA_ARGS__);                            \
    } else if ("WilsonMatrix" == (ctype)) {                                 \
      (p_ret) = fname<WilsonMatrix>(__VA_ARGS__);                           \
    } else if ("WilsonVector" == (ctype)) {                                 \
      (p_ret) = fname<WilsonVector>(__VA_ARGS__);                           \
    } else if ("double" == (ctype)) {                                       \
      (p_ret) = fname<double>(__VA_ARGS__);                                 \
    } else if ("Complex" == (ctype)) {                                      \
      (p_ret) = fname<Complex>(__VA_ARGS__);                                \
    } else if ("long" == (ctype)) {                                         \
      (p_ret) = fname<long>(__VA_ARGS__);                                   \
    } else if ("int64_t" == (ctype)) {                                      \
      (p_ret) = fname<int64_t>(__VA_ARGS__);                                \
    } else {                                                                \
      pqerr("%s %s='%s' does not exist.", #fname, #ctype, (ctype).c_str()); \
      (p_ret) = NULL;                                                       \
    }                                                                       \
  }
