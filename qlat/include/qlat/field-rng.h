#pragma once

#include <qlat/setup.h>
#include <qlat/field.h>

namespace qlat
{  //

struct RngField : FieldM<RngState, 1> {
  using FieldM<RngState, 1>::init;
  void init(const Geometry& geo_, const RngState& rs)
  {
    FieldM<RngState, 1>::init(geo_);
    Coordinate total_site = geo().total_site();
#pragma omp parallel for
    for (long index = 0; index < geo().local_volume(); ++index) {
      Coordinate x = geo().coordinate_from_index(index);
      Coordinate xg = geo().coordinate_g_from_l(x);
      long gindex = index_from_coordinate(xg, total_site);
      split_rng_state(get_elem(x), rs, gindex);
    }
  }
  //
  RngField() { FieldM<RngState, 1>::init(); }
  RngField(const RngField& rf)
  {
    (void)rf;
    qassert(false);
  }
};

}  // namespace qlat
