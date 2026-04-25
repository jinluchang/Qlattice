// vim: set ts=2 sw=2 expandtab:

// Copyright (c) 2022 Luchang Jin
// All rights reserved.

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <qlat-utils/qacc.h>

#include <cstdio>

#if defined QLAT_IN_ACC

#define PRINT_ERR_MSG(tag, str)                                             \
  std::printf("%s: %s from '%s' line %d. (IN_ACC)", (tag), (str), __FILE__, \
              __LINE__)

#define qqwarn(str)               \
  {                               \
    PRINT_ERR_MSG("qwarn", #str); \
  }

#define qqerr(str)               \
  {                              \
    PRINT_ERR_MSG("qerr", #str); \
    assert(false);               \
  }

#define qqassert(x)                 \
  {                                 \
    if (not(x)) {                   \
      PRINT_ERR_MSG("qassert", #x); \
      assert(false);                \
    }                               \
  }

#define qqassert_info(x, ...)            \
  {                                      \
    if (not(x)) {                        \
      {__VA_ARGS__};                     \
      PRINT_ERR_MSG("qassert_info", #x); \
      assert(false);                     \
    }                                    \
  }

#else

#define PRINT_ERR_MSG(tag, str)                                               \
  {                                                                           \
    std::string msg =                                                         \
        qlat::ssprintf("%s: %s from '%s' line %d. (id_node=%d id_thread=%d)", \
                       qlat::get_c_str(tag), qlat::get_c_str(str), __FILE__,  \
                       __LINE__, qlat::get_id_node(), qlat::get_id_thread()); \
    qlat::displayln_c_stdout(msg);                                            \
    qlat::Timer::display_stack_always();                                      \
  }

#define qqwarn(str)              \
  {                              \
    PRINT_ERR_MSG("qwarn", str); \
  }

#define qqerr(str)              \
  {                             \
    PRINT_ERR_MSG("qerr", str); \
    throw std::string(str);     \
  }

#define qqassert(x)                 \
  {                                 \
    if (not(x)) {                   \
      PRINT_ERR_MSG("qassert", #x); \
      throw std::string(#x);        \
    }                               \
  }

#define qqassert_info(x, ...)            \
  {                                      \
    if (not(x)) {                        \
      {__VA_ARGS__};                     \
      PRINT_ERR_MSG("qassert_info", #x); \
      throw std::string(#x);             \
    }                                    \
  }

#endif

// Use this version if you want to keep the assert even if SKIP_ASSERT is
// defined
#define Qwarn(str) qqwarn(str)
#define Qerr(str) qqerr(str)
#define Qassert(x) qqassert(x)
#define Qassert_info(x, ...) qqassert_info(x, {__VA_ARGS__})

// compile without assert message and vector in qacc lengh checks
#ifdef SKIP_ASSERT

#define qwarn(str) assert(true)
#define qerr(str) assert(false)
#define qassert(x) assert(true)
#define qassert_info(x, ...) assert(true)

#else

#define qwarn(str) qqwarn(str)
#define qerr(str) qqerr(str)
#define qassert(x) qqassert(x)
#define qassert_info(x, ...) qqassert_info(x, {__VA_ARGS__})

#endif
