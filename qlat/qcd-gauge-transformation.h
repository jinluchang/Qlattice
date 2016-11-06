#pragma once

#include <qlat/qcd.h>

QLAT_START_NAMESPACE

struct GaugeTransform : FieldM<ColorMatrix,1>
{
  virtual const char* cname()
  {
    return "GaugeTransform";
  }
};

inline void gt_apply_gauge_transform(GaugeTransform& gt0, const GaugeTransform& gt1)
  // gt can be the same as gt0
{
  TIMER("gt_apply_gauge_transform");
  assert(is_matching_geo_mult(gt0.geo, gt1.geo));
  const Geometry& geo = gt0.geo;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const ColorMatrix& t1 = gt1.get_elem(xl);
    ColorMatrix& t0 = gt0.get_elem(xl);
    t0 = t1 * t0;
  }
}

inline void gf_apply_gauge_transform_no_comm(GaugeField& gf, const GaugeField& gf0, const GaugeTransform& gt)
  // gf can be the same as gf0
  // assuming comm for gt is done
{
  TIMER("gf_apply_gauge_transform_no_comm");
  assert(is_matching_geo_mult(gf.geo, gf0.geo));
  assert(is_matching_geo(gf.geo, gt.geo));
  const Geometry& geo = gf.geo;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    Coordinate xl = geo.coordinate_from_index(index);
    Vector<ColorMatrix> v = gf.get_elems(xl);
    const Vector<ColorMatrix> v0 = gf0.get_elems_const(xl);
    const ColorMatrix& t0 = gt.get_elem(xl);
    for (int m = 0; m < DIM; ++m) {
      xl[m] += 1;
      const ColorMatrix& t1 = gt.get_elem(xl);
      v[m] = t0 * v0[m] * t1.adjoint();
      xl[m] -= 1;
    }
  }
}

inline void gf_apply_gauge_transform(GaugeField& gf, const GaugeField& gf0, const GaugeTransform& gt)
{
  TIMER("gf_apply_gauge_transform");
  GaugeTransform gt1;
  gt1.init(geo_resize(gt.geo, 1));
  gt1 = gt;
  refresh_expanded(gt1);
  gf_apply_gauge_transform_no_comm(gf, gf0, gt1);
}


QLAT_END_NAMESPACE
