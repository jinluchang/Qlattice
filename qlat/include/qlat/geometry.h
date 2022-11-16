#pragma once

#include <qlat/setup.h>
#include <qlat/mpi.h>

namespace qlat
{  //

qacc bool is_initialized(const Geometry& geo) { return geo.initialized; }

inline bool check_matching_geo(const Geometry& geo1, const Geometry& geo2)
{
  if (is_matching_geo(geo1, geo2)) {
    return true;
  } else {
    displayln("geo1 =\n" + show(geo1));
    displayln("geo2 =\n" + show(geo2));
    return false;
  }
}

inline bool check_matching_geo_mult(const Geometry& geo1, const Geometry& geo2)
{
  if (is_matching_geo_mult(geo1, geo2)) {
    return true;
  } else {
    displayln("geo1 =\n" + show(geo1));
    displayln("geo2 =\n" + show(geo2));
    return false;
  }
}

}  // namespace qlat
