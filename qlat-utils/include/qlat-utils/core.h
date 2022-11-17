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

#pragma once

#include <qlat-utils/config.h>

#include <qlat-utils/show.h>
#include <qlat-utils/timer.h>
#include <qlat-utils/qacc.h>
#include <qlat-utils/qacc-func.h>
#include <qlat-utils/complex.h>
#include <qlat-utils/array.h>
#include <qlat-utils/vector.h>

#include <unistd.h>
#include <cassert>

#ifdef SKIP_ASSERT
#define qassert(x) assert(true)
#elif defined QLAT_IN_ACC
#define qassert(x) assert(x)
#else
#define qassert(x) qqassert(x)
#endif

// #define SKIP_ASSERT

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
