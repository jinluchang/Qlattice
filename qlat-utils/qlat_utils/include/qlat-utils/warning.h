#pragma once

#if __GNUC__ >= 13
#define QLAT_PUSH_DIAGNOSTIC_DISABLE_DANGLING_REF \
  _Pragma("GCC diagnostic push")                  \
      _Pragma("GCC diagnostic ignored \"-Wdangling-reference\"")
#define QLAT_DIAGNOSTIC_POP _Pragma("GCC diagnostic pop")
#else
#define QLAT_PUSH_DIAGNOSTIC_DISABLE_DANGLING_REF
#define QLAT_DIAGNOSTIC_POP
#endif

/*
 *
Example:
 *
QLAT_PUSH_DIAGNOSTIC_DISABLE_DANGLING_REF;
const ShufflePlan& sp = get_shuffle_plan(f.geo().total_site(), new_size_node);
QLAT_DIAGNOSTIC_POP;
 *
 */
