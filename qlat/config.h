// vim: set ts=2 sw=2 expandtab:

#pragma once

#include <complex>
#include <array>

#define QLAT_START_NAMESPACE namespace qlat {
#define QLAT_END_NAMESPACE }

#define USE_MULTI_NODE

#include <timer.h>


QLAT_START_NAMESPACE

const int DIM = 4;

const int NUM_COLOR = 3;

typedef std::complex<double> Complex;

const char* const cname = "Qlat";

inline void warn(const std::string& str = "")
{
  if (str != "") {
    displayln_info(ssprintf("WARNING: %s", str.c_str()), stderr);
  }
}

QLAT_END_NAMESPACE

#include <qlat/coordinate.h>
