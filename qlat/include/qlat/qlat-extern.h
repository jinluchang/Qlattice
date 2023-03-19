#pragma once

#include <qlat-utils/mat.h>
#include <qlat-utils/matrix-hmc.h>

#define QLAT_CALL_WITH_TYPES(FUNC) \
  FUNC(ColorMatrix)                \
  FUNC(WilsonMatrix)               \
  FUNC(NonRelWilsonMatrix)         \
  FUNC(IsospinMatrix)              \
  FUNC(SpinMatrix)                 \
  FUNC(WilsonVector)               \
  FUNC(ComplexD)                   \
  FUNC(ComplexF)                   \
  FUNC(double)                     \
  FUNC(float)                      \
  FUNC(int64_t)                    \
  FUNC(char)                       \
  FUNC(int8_t)
