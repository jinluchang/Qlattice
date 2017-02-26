#pragma once

#include <qlat/config.h>
#include <qlat/field.h>

QLAT_START_NAMESPACE

struct RngField : FieldM<RngState,1>
{
  virtual const std::string& cname()
  {
    static const std::string s = "RngField";
    return s;
  }
  //
  virtual void init()
  {
    FieldM<RngState,1>::init();
  }
  virtual void init(const Geometry& geo_, const RngState& rs)
  {
    FieldM<RngState,1>::init(geo_);
    Coordinate total_site = geo.total_site();
#pragma omp parallel for
    for (long index = 0; index < geo.local_volume(); ++index) {
      Coordinate x = geo.coordinate_from_index(index);
      Coordinate xg = geo.coordinate_g_from_l(x);
      long gindex = index_from_coordinate(xg, total_site);
      split_rng_state(get_elem(x), rs, gindex);
    }
  }
  //
  RngField()
  {
    FieldM<RngState,1>::init();
  }
  RngField(const RngField& rf)
  {
    qassert(false);
  }
};

QLAT_END_NAMESPACE
