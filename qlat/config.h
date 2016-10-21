// vim: set ts=2 sw=2 expandtab:

#pragma once

#include <complex>
#include <array>

#define QLAT_START_NAMESPACE namespace qlat {
#define QLAT_END_NAMESPACE }

#define USE_MULTI_NODE

QLAT_START_NAMESPACE

const int DIM = 4;

const int NUM_COLOR = 3;

typedef std::complex<double> Complex;

const char* const cname = "Qlat";

QLAT_END_NAMESPACE

#include <qlat/coordinate.h>
