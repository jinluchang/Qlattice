#pragma once

#include <qlat-utils/show.h>
#include <qlat-utils/timer.h>
#include <qlat-utils/qacc.h>
#include <qlat-utils/qacc-func.h>
#include <qlat-utils/array.h>
#include <qlat-utils/complex.h>

#include <unistd.h>
#include <cassert>

#define qwarn(str)                                                     \
  {                                                                    \
    std::string msg = qlat::ssprintf(                                  \
        "qwarn: %s from '%s' line %d. (id_node=%d thread_num=%d "      \
        "id_node_in_shuffle=%d)",                                      \
        qlat::get_c_str(str), __FILE__, __LINE__, qlat::get_id_node(), \
        qlat::get_thread_num(), qlat::get_id_node_in_shuffle());       \
    qlat::displayln(msg);                                              \
    qlat::Timer::display_stack_always();                               \
  }

#define qqassert(x)                      \
  {                                      \
    if (not(x)) {                        \
      qwarn("qassert: " #x);             \
      usleep((useconds_t)(0.1 * 1.0e6)); \
      throw std::string("qassert");      \
    }                                    \
  }

// #define SKIP_ASSERT

#ifdef SKIP_ASSERT
#define qassert(x) assert(true)
#elif defined QLAT_IN_ACC
#define qassert(x) assert(x)
#else
#define qassert(x) qqassert(x)
#endif

namespace qlat
{  //

const double PI = 3.141592653589793;

template <class T>
qacc T sqr(const T& x)
{
  return x * x;
}

// ----------------

qacc long modl(const long x, const long len)
{
  qassert(0 < len);
  const int m = x % len;
  if (0 <= m) {
    return m;
  } else {
    return m + len;
  }
}

qacc int mod(const int x, const int len)
{
  qassert(0 < len);
  const int m = x % len;
  if (0 <= m) {
    return m;
  } else {
    return m + len;
  }
}

qacc double mod(const double x, const double len)
{
  qassert(0 < len);
  const double m = x - trunc(x / len) * len;
  if (0 <= m) {
    return m;
  } else {
    return m + len;
  }
}

qacc int smod(const int x, const int len)
{
  qassert(0 < len);
  const int m = mod(x, len);
  if (m * 2 < len) {
    return m;
  } else {
    return m - len;
  }
}

qacc double smod(const double x, const double len)
{
  qassert(0 < len);
  const double m = mod(x, len);
  if (m * 2 < len) {
    return m;
  } else {
    return m - len;
  }
}

qacc double smod_sym(const double x, const double len,
                     const double eps = 1.0e-8)
{
  const double m = smod(x, len);
  if (std::abs(std::abs(m * 2) - len) < eps) {
    return 0;
  } else {
    return m;
  }
}

qacc int middle_mod(const int x, const int y, const int len)
{
  qassert(0 < len);
  const int xm = mod(x, len);
  const int ym = mod(y, len);
  if (xm <= ym) {
    const int r = smod(ym - xm, len);
    return mod(xm + r / 2, len);
  } else {
    const int r = smod(xm - ym, len);
    return mod(ym + r / 2, len);
  }
}

qacc double middle_mod(const double x, const double y, const double len)
{
  qassert(0 < len);
  const double xm = mod(x, len);
  const double ym = mod(y, len);
  if (xm <= ym) {
    const double r = smod(ym - xm, len);
    return mod(xm + r / 2, len);
  } else {
    const double r = smod(xm - ym, len);
    return mod(ym + r / 2, len);
  }
}

// ----------------

}  // namespace qlat
