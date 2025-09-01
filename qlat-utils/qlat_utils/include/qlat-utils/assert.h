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

#else

#define MK_ERR_MSG(tag, str)                                            \
  qlat::ssprintf("%s: %s from '%s' line %d. (id_node=%d id_thread=%d)", \
                 qlat::get_c_str(tag), qlat::get_c_str(str), __FILE__,  \
                 __LINE__, qlat::get_id_node(), qlat::get_id_thread())

#define qqwarn(str)                             \
  {                                             \
    std::string msg = MK_ERR_MSG("qwarn", str); \
    qlat::displayln_c_stdout(msg);              \
    qlat::Timer::display_stack_always();        \
  }

#define qqerr(str)                             \
  {                                            \
    std::string msg = MK_ERR_MSG("qerr", str); \
    qlat::displayln_c_stdout(msg);             \
    qlat::Timer::display_stack_always();       \
    throw std::string(msg);                    \
  };

#define qqassert(x)                                \
  {                                                \
    if (not(x)) {                                  \
      std::string msg = MK_ERR_MSG("qassert", #x); \
      qlat::displayln_c_stdout(msg);               \
      qlat::Timer::display_stack_always();         \
      throw std::string(msg);                      \
    }                                              \
  }

#endif

#define Qassert(x) qqassert(x)

// compile without assert message and vector in qacc lengh checks
#ifdef SKIP_ASSERT

#define qassert(x) assert(true)
#define qwarn(str) assert(true)
#define qerr(str) assert(false)

#else

#define qassert(x) qqassert(x)
#define qwarn(str) qqwarn(str)
#define qerr(str) qqerr(str)

#endif
