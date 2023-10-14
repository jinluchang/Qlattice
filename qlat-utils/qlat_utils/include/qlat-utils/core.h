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

#include <qlat-utils/array.h>
#include <qlat-utils/assert.h>
#include <qlat-utils/complex.h>
#include <qlat-utils/config.h>
#include <qlat-utils/qacc-func.h>
#include <qlat-utils/qacc.h>
#include <qlat-utils/show.h>
#include <qlat-utils/timer.h>
#include <qlat-utils/vector.h>
#include <qlat-utils/version.h>
#include <unistd.h>

namespace qlat
{  //

const double PI = 3.141592653589793;

template <class T>
qacc T sqr(const T& x)
{
  return x * x;
}

inline void warn(const std::string& str = "")
{
  if (str != "") {
    displayln_info(ssprintf("WARNING: %s", str.c_str()));
  }
}

}  // namespace qlat
