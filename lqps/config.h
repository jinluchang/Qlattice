#pragma once

#include <complex>
#include <array>

#define LQPS_START_NAMESPACE namespace lqps {
#define LQPS_END_NAMESPACE }

#define USE_MULTI_NODE

LQPS_START_NAMESPACE

using Complex = std::complex<double>;

const int DIM = 4;

using Coordinate = std::array<int,DIM>;

LQPS_END_NAMESPACE
