// vim: set ts=2 sw=2 expandtab:

#pragma once

#include <qlat/config.h>

#include <endian.h>

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

namespace qlat
{  //

inline std::string show(const Complex& x)
{
  return ssprintf("(%24.17E + %24.17E j)", x.real(), x.imag());
}

}  // namespace qlat
