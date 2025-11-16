#pragma once

#include <fftw3.h>
#include <qlat/fermion-action.h>
#include <qlat/qcd-prop.h>
#include <qlat/qcd-utils.h>
#include <qlat/qcd.h>
#include <qlat/selected-points.h>

namespace qlat
{  //

void gt_apply_gauge_transformation(GaugeTransform& gt0,
                                   const GaugeTransform& gt1);

void gt_apply_gauge_transformation(GaugeTransform& gt,
                                   const GaugeTransform& gt0,
                                   const GaugeTransform& gt1);

void gf_apply_gauge_transformation(GaugeField& gf, const GaugeField& gf0,
                                   const GaugeTransform& gt,
                                   const bool is_dagger = false);

void gt_invert(GaugeTransform& gt, const GaugeTransform& gt0);

void ff_apply_gauge_transformation(FermionField4d& ff,
                                   const FermionField4d& ff0,
                                   const GaugeTransform& gt);

void prop_apply_gauge_transformation(Propagator4d& prop,
                                     const Propagator4d& prop0,
                                     const GaugeTransform& gt);

void prop_apply_gauge_transformation(SelectedField<WilsonMatrix>& prop,
                                     const SelectedField<WilsonMatrix>& prop0,
                                     const GaugeTransform& gt,
                                     const FieldSelection& fsel);

void prop_apply_gauge_transformation(vector<WilsonMatrix>& prop,
                                     const vector<WilsonMatrix>& prop0,
                                     const GaugeTransform& gt,
                                     const std::vector<Coordinate>& pcs);

void prop_apply_gauge_transformation(SelectedPoints<WilsonMatrix>& prop,
                                     const SelectedPoints<WilsonMatrix>& prop0,
                                     const GaugeTransform& gt,
                                     const PointsSelection& psel);

void gf_apply_rand_gauge_transformation(GaugeField& gf, const GaugeField& gf0,
                                        const RngState& rs);

void make_temporal_gauge_transformation(GaugeTransform& gt,
                                        const GaugeField& gf,
                                        const Int tgref = 0, const Int dir = 3);

void make_tree_gauge_transformation(GaugeTransform& gt, const GaugeField& gf,
                                    const Coordinate& xgref = Coordinate(0, 0,
                                                                         0, 0),
                                    const Coordinate& dirs = Coordinate(0, 1, 2,
                                                                        3));

template <class Inverter>
struct GaugeTransformInverter
// gt_inv should be: gt_invert(gt_inv, gt);
// the result should be the same as invert with gf_fix where
// gf_fix is: gf_apply_gauge_transformation(gf_fix, gf, gt);
{
  box<Geometry> geo;
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
    gt_invert(gt_inv, gt);
    geo.set(inv().geo());
    fa = inv().fa;
    gf = inv().gf;
  }
};

template <class Inverter>
inline void invert(FermionField4d& out, const FermionField4d& in,
                   const GaugeTransformInverter<Inverter>& gtinv)
{
  TIMER_VERBOSE("invert(out,in,gt_inv)");
  const Inverter& inv = gtinv.inv();
  FermionField4d src;
  ff_apply_gauge_transformation(src, in, gtinv.gt_inv);
  invert(out, src, inv);
  ff_apply_gauge_transformation(out, out, gtinv.gt);
}

// -------------------------------------------------------------------------

template <class Inverter>
void set_wall_src_propagator(Propagator4d& prop, const Int tslice,
                             const CoordinateD& lmom, const Inverter& inv,
                             const GaugeTransform& gt,
                             const GaugeTransform& gt_inv)
// gt_inv should be: gt_invert(gt_inv, gt);
// the result should be the same as invert with gf_fix where
// gf_fix is: gf_apply_gauge_transformation(gf_fix, gf, gt);
{
  TIMER_VERBOSE("set_wall_src_propagator");
  warn("obsolete");
  const Geometry geo = geo_resize(inv.geo());
  prop.init(geo);
  qassert(prop.geo() == geo);
  FermionField4d src, sol;
  src.init(geo);
  sol.init(geo);
  for (Int cs = 0; cs < 4 * NUM_COLOR; ++cs) {
    set_tslice_mom_src_fermion_field(src, tslice, lmom, cs);
    ff_apply_gauge_transformation(src, src, gt_inv);
    set_zero(sol);
    invert(sol, src, inv);
    ff_apply_gauge_transformation(sol, sol, gt);
    set_propagator_col_from_fermion_field(prop, cs, sol);
  }
}

}  // namespace qlat
