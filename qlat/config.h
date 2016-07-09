#pragma once

#include <complex>
#include <array>

#define QLAT_START_NAMESPACE namespace qlat {
#define QLAT_END_NAMESPACE }

#define USE_MULTI_NODE

QLAT_START_NAMESPACE

const int DIM = 4;

typedef std::complex<double> Complex;

using Coordinate = std::array<int,DIM>;

const char* const cname = "QLAT";

QLAT_END_NAMESPACE
