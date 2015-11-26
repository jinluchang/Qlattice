#pragma once

#include <complex>
#include <array>

#define LQPS_START_NAMESPACE namespace lqps {
#define LQPS_END_NAMESPACE }

#define USE_MULTI_NODE

LQPS_START_NAMESPACE

const int DIM = 4;

using Coordinate = std::array<int,DIM>;

using Complex = std::complex<double>;

LQPS_END_NAMESPACE
