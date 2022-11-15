#pragma once

#include <qlat/setup.h>
#include <qlat/mpi.h>

namespace qlat
{  //

qacc bool operator==(const Geometry& geo1, const Geometry& geo2)
{
  return geo1.initialized == geo2.initialized && geo1.eo == geo2.eo &&
         geo1.geon == geo2.geon && geo1.multiplicity == geo2.multiplicity &&
         geo1.node_site == geo2.node_site &&
         geo1.expansion_left == geo2.expansion_left &&
         geo1.expansion_right == geo2.expansion_right &&
         geo1.node_site_expanded == geo2.node_site_expanded &&
         geo1.is_only_local == geo2.is_only_local;
}

qacc bool operator!=(const Geometry& geo1, const Geometry& geo2)
{
  return !(geo1 == geo2);
}

qacc Geometry geo_resize(const Geometry& geo_, const int thick = 0)
{
  Geometry geo = geo_;
  geo.resize(thick);
  return geo;
}

qacc Geometry geo_resize(const Geometry& geo_,
                         const Coordinate& expansion_left_,
                         const Coordinate& expansion_right_)
{
  Geometry geo = geo_;
  geo.resize(expansion_left_, expansion_right_);
  return geo;
}

qacc Geometry geo_remult(const Geometry& geo_, const int multiplicity_ = 1)
{
  Geometry geo = geo_;
  geo.remult(multiplicity_);
  return geo;
}

qacc Geometry geo_reform(const Geometry& geo_, const int multiplicity_ = 1,
                         const int thick = 0)
// do not change eo
{
  Geometry geo = geo_;
  geo.remult(multiplicity_);
  geo.resize(thick);
  return geo;
}

qacc Geometry geo_reform(const Geometry& geo_, const int multiplicity_,
                         const Coordinate& expansion_left_,
                         const Coordinate& expansion_right_)
// do not change eo
{
  Geometry geo = geo_;
  geo.remult(multiplicity_);
  geo.resize(expansion_left_, expansion_right_);
  return geo;
}

qacc Geometry geo_eo(const Geometry& geo_, const int eo_ = 0)
// 0:regular; 1:odd; 2:even
{
  Geometry geo = geo_;
  geo.eo = eo_;
  return geo;
}

qacc bool is_initialized(const Geometry& geo) { return geo.initialized; }

inline std::string show(const qlat::Geometry& geo)
{
  std::string s;
  s += ssprintf("{ initialized   = %s\n", show(geo.initialized).c_str());
  s += ssprintf(", eo            = %s\n", show(geo.eo).c_str());
  s += ssprintf(", geon          =\n%s\n", show(geo.geon).c_str());
  s += ssprintf(", multiplicity  = %s\n", show(geo.multiplicity).c_str());
  s += ssprintf(", node_site     = %s\n", show(geo.node_site).c_str());
  s += ssprintf(", expan_left    = %s\n", show(geo.expansion_left).c_str());
  s += ssprintf(", expan_right   = %s\n", show(geo.expansion_right).c_str());
  s += ssprintf(", node_site_exp = %s\n", show(geo.node_site_expanded).c_str());
  s += ssprintf(", is_only_local = %s }", show(geo.is_only_local).c_str());
  return s;
}

qacc bool is_matching_geo(const Geometry& geo1, const Geometry& geo2)
{
  return geo1.initialized == geo2.initialized && geo1.geon == geo2.geon &&
         geo1.node_site == geo2.node_site;
}

qacc bool is_matching_geo_mult(const Geometry& geo1, const Geometry& geo2)
{
  return is_matching_geo(geo1, geo2) && geo1.eo == geo2.eo &&
         geo1.multiplicity == geo2.multiplicity;
}

qacc bool is_matching_geo_included(const Geometry& geo1, const Geometry& geo2)
{
  bool include = is_matching_geo_mult(geo1, geo2);
  for(int i=0;i<4;i++) {
    if(geo2.expansion_left[i] < geo1.expansion_left[i]) {
      include = false;
    }
  }
  for(int i=0;i<4;i++) {
    if(geo2.expansion_right[i] < geo1.expansion_right[i]) {
      include = false;
    }
  }
  return include;
}

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
