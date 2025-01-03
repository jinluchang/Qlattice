#pragma once

#if __GNUC__ >= 13
#define PUSH_DIAGNOSTIC_DISABLE_DANGLING_REF \
  _Pragma("GCC diagnostic push")             \
      _Pragma("GCC diagnostic ignored \"-Wdangling-reference\"")
#define DIAGNOSTIC_POP _Pragma("GCC diagnostic pop")
#else
#define PUSH_DIAGNOSTIC_DISABLE_DANGLING_REF
#define DIAGNOSTIC_POP
#endif
