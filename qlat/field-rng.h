#pragma once

#include <qlat/config.h>
#include <qlat/field.h>

#include <rng-state.h>

#include <cassert>

QLAT_START_NAMESPACE

struct RngField : FieldM<RngState,1>
{
  virtual const char* cname()
  {
    return "RngField";
  }
  //
  virtual void init()
  {
    FieldM<RngState,1>::init();
  }
  virtual void init(const Geometry& geo_, const RngState& rs)
  {
    FieldM<RngState,1>::init(geo_);
    Coordinate totalSite;
    for (int mu = 0; mu < DIM; ++mu) {
      totalSite[mu] = geo.totalSite(mu);
    }
#pragma omp parallel for
    for (long index = 0; index < geo.localVolume(); ++index) {
      Coordinate x; geo.coordinateFromIndex(x, index);
      Coordinate xg; geo.coordinateGfL(xg, x);
      long gindex = indexFromCoordinate(xg, totalSite);
      splitRngState(getElem(x), rs, gindex);
    }
  }
  //
  RngField()
  {
    FieldM<RngState,1>::init();
  }
  RngField(const RngField& rf)
  {
    assert(false);
  }
};

QLAT_END_NAMESPACE
