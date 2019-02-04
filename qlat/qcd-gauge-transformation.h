#pragma once

#include <qlat/fermion-action.h>
#include <qlat/qcd-prop.h>
#include <qlat/qcd-utils.h>
#include <qlat/qcd.h>

#include <fftw3.h>

QLAT_START_NAMESPACE

struct GaugeTransform : FieldM<ColorMatrix, 1> {
  virtual const std::string& cname()
  {
    static const std::string s = "GaugeTransform";
    return s;
  }
};

struct U1GaugeTransform : FieldM<ComplexF, 1> {
  virtual const std::string& cname()
  {
    static const std::string s = "U1GaugeTransform";
    return s;
  }
};

inline void gt_apply_gauge_transformation(GaugeTransform& gt0,
                                          const GaugeTransform& gt1)
// gt0 can be the same as gt1
// gt0 <- gt1 * gt0
{
  TIMER("gt_apply_gauge_transformation");
  qassert(is_matching_geo_mult(gt0.geo, gt1.geo));
  const Geometry& geo = gt0.geo;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const ColorMatrix& t1 = gt1.get_elem(xl);
    ColorMatrix& t0 = gt0.get_elem(xl);
    t0 = t1 * t0;
  }
}

inline void gf_apply_gauge_transformation_no_comm(GaugeField& gf,
                                                  const GaugeField& gf0,
                                                  const GaugeTransform& gt)
// gf can be the same as gf0
// assuming comm for gt is done
// gf <- gt * gf0
{
  TIMER("gf_apply_gauge_transformation_no_comm");
  qassert(is_matching_geo(gf0.geo, gt.geo));
  const Geometry& geo = gf0.geo;
  gf.init(geo_resize(geo, 0));
  qassert(is_matching_geo(gf.geo, gf0.geo));
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    Coordinate xl = geo.coordinate_from_index(index);
    Vector<ColorMatrix> v = gf.get_elems(xl);
    const Vector<ColorMatrix> v0 = gf0.get_elems_const(xl);
    const ColorMatrix& t0 = gt.get_elem(xl);
    for (int m = 0; m < DIMN; ++m) {
      xl[m] += 1;
      const ColorMatrix& t1 = gt.get_elem(xl);
      v[m] = t0 * v0[m] * matrix_adjoint(t1);
      xl[m] -= 1;
    }
  }
}

inline void gf_apply_gauge_transformation(GaugeField& gf, const GaugeField& gf0,
                                          const GaugeTransform& gt)
{
  TIMER("gf_apply_gauge_transformation");
  qassert(is_matching_geo(gf0.geo, gt.geo));
  GaugeTransform gt1;
  gt1.init(geo_resize(gt.geo, 1));
  gt1 = gt;
  refresh_expanded(gt1);
  gf_apply_gauge_transformation_no_comm(gf, gf0, gt1);
}

inline void gt_inverse(GaugeTransform& gt, const GaugeTransform& gt0)
{
  TIMER("gt_inverse");
  gt.init(geo_resize(gt0.geo));
  const Geometry& geo = gt.geo;
  qassert(is_matching_geo_mult(gt.geo, gt0.geo));
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const ColorMatrix& t = gt.get_elem(xl);
    gt.get_elem(xl) = matrix_adjoint(gt0.get_elem(xl));
  }
}

inline void ff_apply_gauge_transformation(FermionField4d& ff,
                                          const FermionField4d& ff0,
                                          const GaugeTransform& gt)
{
  TIMER("ff_apply_gauge_transformation");
  qassert(is_matching_geo(ff0.geo, gt.geo));
  const Geometry& geo = ff0.geo;
  ff.init(geo_resize(geo));
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<WilsonVector> v = ff.get_elems(xl);
    const Vector<WilsonVector> v0 = ff0.get_elems_const(xl);
    const ColorMatrix& t = gt.get_elem(xl);
    for (int m = 0; m < v0.size(); ++m) {
      v[m] = t * v0[m];
    }
  }
}

inline void prop_apply_gauge_transformation(Propagator4d& prop,
                                            const Propagator4d& prop0,
                                            const GaugeTransform& gt)
{
  TIMER("prop_apply_gauge_transformation");
  qassert(is_matching_geo(prop0.geo, gt.geo));
  const Geometry& geo = prop0.geo;
  prop.init(geo_resize(geo));
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<WilsonMatrix> v = prop.get_elems(xl);
    const Vector<WilsonMatrix> v0 = prop0.get_elems_const(xl);
    const ColorMatrix& t = gt.get_elem(xl);
    for (int m = 0; m < v0.size(); ++m) {
      v[m] = t * v0[m];
    }
  }
}

inline void gf_apply_rand_gauge_transformation(GaugeField& gf,
                                               const GaugeField& gf0,
                                               const RngState& rs)
{
  const Geometry geo = geo_reform(gf0.geo);
  GaugeTransform gt;
  gt.init(geo);
  set_g_rand_color_matrix_field(gt, rs, 1.0);
  gf_apply_gauge_transformation(gf, gf0, gt);
}

inline void make_temporal_gauge_transformation(GaugeTransform& gt,
                                               const GaugeField& gf,
                                               const int tgref = 0,
                                               const int dir = 3)
// after tranform: ``gf.get_elem(xl, dir) = unit'' is true from ``xg[dir] =
// tgref'' until as far as possible
// ``gt.get_elem(xl) = unit'' if ``xg[dir] = tgref''
{
  TIMER("make_temporal_gauge_transformation");
  const Geometry geo = geo_reform(gf.geo, 0);
  gt.init(geo);
  assert(is_matching_geo(gt.geo, gf.geo));
  Coordinate expension_left, expension_right;
  set_zero(expension_left);
  set_zero(expension_right);
  expension_left[dir] = 1;
  const Geometry geo1 = geo_resize(geo, expension_left, expension_right);
  GaugeField gf1;
  gf1.init(geo1);
  gf1 = gf;
  refresh_expanded(gf1);
  GaugeTransform gt1;
  gt1.init(geo1);
  set_unit(gt1);
  const Coordinate total_site = geo.total_site();
  for (int tgrel = 1; tgrel < total_site[dir]; ++tgrel) {
    refresh_expanded(gt1);
    const int tg = mod(tgref + tgrel, total_site[dir]);
#pragma omp parallel for
    for (long index = 0; index < geo.local_volume(); ++index) {
      Coordinate xl = geo.coordinate_from_index(index);
      Coordinate xg = geo.coordinate_g_from_l(xl);
      if (tg == xg[dir]) {
        ColorMatrix& t1 = gt1.get_elem(xl);
        xl[dir] -= 1;
        const ColorMatrix& v = gf1.get_elem(xl, dir);
        const ColorMatrix& t0 = gt1.get_elem(xl);
        t1 = t0 * v;
      }
    }
  }
  gt = gt1;
}

inline void make_tree_gauge_transformation(
    GaugeTransform& gt, const GaugeField& gf,
    const Coordinate& xgref = Coordinate(0, 0, 0, 0),
    const Coordinate& dirs = Coordinate(0, 1, 2, 3))
{
  TIMER("make_tree_gauge_transformation");
  const Geometry& geo = geo_reform(gf.geo);
  if (false == is_initialized(gt)) {
    gt.init(geo);
  }
  assert(is_matching_geo(gt.geo, gf.geo));
  set_unit(gt);
  GaugeTransform gt_dir;
  gt_dir.init(geo);
  GaugeField gft;
  gft.init(geo);
  gft = gf;
  for (int m = 0; m < DIMN; ++m) {
    make_temporal_gauge_transformation(gt_dir, gft, xgref[dirs[m]], dirs[m]);
    gf_apply_gauge_transformation(gft, gft, gt_dir);
    gt_apply_gauge_transformation(gt, gt_dir);
  }
}

template <class Inverter>
struct GaugeTransformInverter
// gt_inv should be: gt_inverse(gt_inv, gt);
// the result should be the same as inverse with gf_fix where
// gf_fix is: gf_apply_gauge_transformation(gf_fix, gf, gt);
{
  Geometry geo;
  FermionAction fa;
  GaugeField gf;
  //
  ConstHandle<Inverter> inv;
  GaugeTransform gt, gt_inv;
  //
  GaugeTransformInverter() { init(); }
  GaugeTransformInverter(const Inverter& inv_, const GaugeTransform& gt_)
  {
    init(inv_, gt_);
  }
  //
  void init()
  {
    geo.init();
    fa.init();
    gf.init();
    inv.init();
    gt.init();
    gt_inv.init();
  }
  //
  void init(const Inverter& inv_, const GaugeTransform& gt_)
  {
    inv.init(inv_);
    gt = gt_;
    gt_inverse(gt_inv, gt);
    geo = inv().geo;
    fa = inv().fa;
    gf = inv().gf;
  }
};

template <class Inverter>
inline void inverse(FermionField4d& out, const FermionField4d& in,
                    const GaugeTransformInverter<Inverter>& gtinv)
{
  TIMER_VERBOSE("inverse(out,in,gt_inv)");
  const Inverter& inv = gtinv.inv();
  const Geometry geo = geo_reform(inv.geo);
  FermionField4d src;
  ff_apply_gauge_transformation(src, in, gtinv.gt_inv);
  inverse(out, src, inv);
  ff_apply_gauge_transformation(out, out, gtinv.gt);
}

// -------------------------------------------------------------------------

template <class Inverter>
void set_wall_src_propagator(Propagator4d& prop, const int tslice,
                             const CoordinateD& lmom, const Inverter& inv,
                             const GaugeTransform& gt,
                             const GaugeTransform& gt_inv)
// gt_inv should be: gt_inverse(gt_inv, gt);
// the result should be the same as inverse with gf_fix where
// gf_fix is: gf_apply_gauge_transformation(gf_fix, gf, gt);
{
  TIMER_VERBOSE("set_wall_src_propagator");
  warn("obsolete");
  const Geometry geo = geo_reform(inv.geo);
  prop.init(geo);
  qassert(prop.geo == geo);
  FermionField4d src, sol;
  src.init(geo);
  sol.init(geo);
  for (int cs = 0; cs < 4 * NUM_COLOR; ++cs) {
    set_tslice_mom_src_fermion_field(src, tslice, lmom, cs);
    ff_apply_gauge_transformation(src, src, gt_inv);
    set_zero(sol);
    inverse(sol, src, inv);
    ff_apply_gauge_transformation(sol, sol, gt);
    set_propagator_col_from_fermion_field(prop, cs, sol);
  }
}

QLAT_END_NAMESPACE
