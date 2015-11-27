#pragma once

#include <lqps/config.h>
#include <lqps/utils.h>
#include <lqps/mpi.h>
#include <lqps/geometry.h>

#include <vector>

LQPS_START_NAMESPACE

template <class M>
struct Field
{
  bool initialized;
  Geometry geo;
  std::vector<M> field;
};

LQPS_END_NAMESPACE
