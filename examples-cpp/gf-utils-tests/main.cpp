#include <qlat/qlat.h>

#include <cmath>
#include <complex>
#include <iostream>

using namespace qlat;

void simple_tests()
{
  TIMER_VERBOSE("simple_tests");
  RngState rs(get_global_rng_state(), fname);
  // const Coordinate total_site(16, 16, 16, 32);
  // const Coordinate total_site(8, 8, 8, 8);
  const Coordinate total_site(4, 4, 4, 8);
  Geometry geo;
  geo.init(total_site);
  GaugeField gf;
  gf.init(geo);
  set_unit(gf);
  gf_show_info(gf);
  set_g_rand_color_matrix_field(gf, RngState(rs, "gf-0.1"), 0.1);
  gf_show_info(gf, 1);
  set_g_rand_color_matrix_field(gf, RngState(rs, "gf-0.3"), 0.3);
  gf_show_info(gf, 1);
  GaugeTransform gt;
  gt.init(geo);
  set_g_rand_color_matrix_field(gt, RngState(rs, "gt-1.0"), 1.0);
  gf_apply_gauge_transformation(gf, gf, gt);
  gf_show_info(gf, 1);
  gf_apply_gauge_transformation(gf, gf, gt);
  gf_show_info(gf);
  unitarize(gf);
  gf_show_info(gf, 1);
  make_tree_gauge_transformation(gt, gf);
  gf_apply_gauge_transformation(gf, gf, gt);
  gf_show_info(gf);
  GaugeField gfs;
  gf_ape_smear(gfs, gf, 0.1);
  gf_show_info(gfs, 1);
  gf_hyp_smear(gfs, gf, 0.1, 0.0, 0.3);
  gf_show_info(gfs);
  gf_hyp_smear(gfs, gf, 0.75, 0.6, 0.3);
  gf_show_info(gfs, 1);
  //
  gf_ape_smear_spatial(gfs, gf, 0.5);
  GaugeField gfs1;
  set_left_expanded_gauge_field(gfs1, gfs);
  Propagator4d prop;
  prop.init(geo);
  prop_smear(prop, gfs1, 0.5, 10);
  //
  {
    TIMER_VERBOSE("test-gf-stout-smear");
    GaugeField gf_orig;
    gf_orig = gf;
    GaugeField gfs_stout, gfs_block, gfs_stout_alias;
    gf_stout_smear(gfs_stout, gf, 0.1, 1);
    gf_stout_smear_block(gfs_block, gf, Coordinate(), 0.1, 1);
    gfs_block -= gfs_stout;
    displayln_info(
        ssprintf("CHECK: gf_stout_smear vs gf_stout_smear_block "
                 "(no blocking): qnorm(diff) = %.2E",
                 qnorm(gfs_block)));
    gfs_stout_alias = gf;
    gf_stout_smear(gfs_stout_alias, gfs_stout_alias, 0.1, 1);
    gfs_stout_alias -= gfs_stout;
    displayln_info(
        ssprintf("CHECK: gf_stout_smear with gf == gf0: "
                 "qnorm(diff) = %.2E",
                 qnorm(gfs_stout_alias)));
    GaugeField gfs_stout_1, gfs_stout_2, gfs_stout_12, gfs_stout_zero;
    gf_stout_smear(gfs_stout_1, gf, 0.1, 1);
    gf_stout_smear(gfs_stout_2, gfs_stout_1, 0.1, 1);
    gf_stout_smear(gfs_stout_12, gf, 0.1, 2);
    gfs_stout_12 -= gfs_stout_2;
    displayln_info(
        ssprintf("CHECK: gf_stout_smear num_step = 2 vs two num_step = 1: "
                 "qnorm(diff) = %.2E",
                 qnorm(gfs_stout_12)));
    gf_stout_smear(gfs_stout_zero, gf, 0.1, 0);
    gfs_stout_zero -= gf;
    displayln_info(
        ssprintf("CHECK: gf_stout_smear with num_step = 0: "
                 "qnorm(diff) = %.2E",
                 qnorm(gfs_stout_zero)));
    GaugeField gfs_no_proj, gfs_no_proj_2, gfs_no_proj_zero, gfs_no_proj_alias;
    gf_ape_smear_no_proj(gfs_no_proj, gf, 0.1, 1);
    gf_ape_smear_no_proj(gfs_no_proj_2, gfs_no_proj, 0.1, 1);
    gf_ape_smear_no_proj(gfs_no_proj_alias, gf, 0.1, 2);
    gfs_no_proj_2 -= gfs_no_proj_alias;
    displayln_info(
        ssprintf("CHECK: gf_ape_smear_no_proj num_step = 2 vs two "
                 "num_step = 1: qnorm(diff) = %.2E",
                 qnorm(gfs_no_proj_2)));
    // The alpha = 0 update is the identity, but the unconditional rescale to
    // <N^2> = 3 still multiplies the whole field by one real scalar, so the
    // output is that scalar times the input.  The input is scaled link by
    // link, with a spread in the link norms, so that the scalar differs
    // noticeably from one and neither a missing rescale nor a per-link one
    // can pass; the scalar is derived from the input data as the one bringing
    // <N^2> to 3 instead of being assumed.
    GaugeField gf_scaled;
    gf_scaled = gf;
    qacc_for(index, gf_scaled.geo().local_volume(), {
      const Geometry& geo_scaled = gf_scaled.geo();
      const Coordinate xl = geo_scaled.coordinate_from_index(index);
      Vector<ColorMatrix> v = gf_scaled.get_elems(xl);
      for (Int mu = 0; mu < DIMN; ++mu) {
        const Long k = xl[0] + 2 * xl[1] + 3 * xl[2] + 4 * xl[3] + mu;
        v[mu] *= 0.5 * (1.0 + 0.02 * (RealD)(k % 11));
      }
    });
    const RealD s_expected = std::sqrt(
        3.0 / (qnorm(gf_scaled) / ((RealD)geo.total_volume() * (RealD)DIMN)));
    gfs_no_proj_zero = gf_scaled;
    gf_ape_smear_no_proj(gfs_no_proj_zero, gfs_no_proj_zero, 0.0, 1);
    GaugeField gf_zero_expected, gf_zero_resid;
    gf_zero_expected = gf_scaled;
    gf_zero_expected *= s_expected;
    gf_zero_resid = gfs_no_proj_zero;
    gf_zero_resid -= gf_zero_expected;
    const RealD zero_rel_dev = qnorm(gf_zero_resid) / qnorm(gfs_no_proj_zero);
    displayln_info(
        ssprintf("CHECK: gf_ape_smear_no_proj with alpha = 0 is a uniform real "
                 "rescale = %s",
                 zero_rel_dev < 1e-12 ? "True" : "False"));
    GaugeField gfs_no_proj_alias2;
    gfs_no_proj_alias2 = gf;
    gf_ape_smear_no_proj(gfs_no_proj_alias2, gfs_no_proj_alias2, 0.1, 1);
    gfs_no_proj_alias2 -= gfs_no_proj;
    displayln_info(
        ssprintf("CHECK: gf_ape_smear_no_proj with gf == gf0: "
                 "qnorm(diff) = %.2E",
                 qnorm(gfs_no_proj_alias2)));
    GaugeField gfs_hyp_1, gfs_hyp_2, gfs_hyp_12, gfs_hyp_zero;
    gf_hyp_smear(gfs_hyp_1, gf, 0.75, 0.6, 0.3, 1);
    gf_hyp_smear(gfs_hyp_2, gfs_hyp_1, 0.75, 0.6, 0.3, 1);
    gf_hyp_smear(gfs_hyp_12, gf, 0.75, 0.6, 0.3, 2);
    gfs_hyp_12 -= gfs_hyp_2;
    displayln_info(
        ssprintf("CHECK: gf_hyp_smear num_step = 2 vs two num_step = 1: "
                 "qnorm(diff) = %.2E",
                 qnorm(gfs_hyp_12)));
    gf_hyp_smear(gfs_hyp_zero, gf, 0.75, 0.6, 0.3, 0);
    gfs_hyp_zero -= gf;
    displayln_info(
        ssprintf("CHECK: gf_hyp_smear with num_step = 0: "
                 "qnorm(diff) = %.2E",
                 qnorm(gfs_hyp_zero)));
    GaugeField gfs_block_zero, gfs_local_zero;
    gf_stout_smear_block(gfs_block_zero, gf, Coordinate(), 0.1, 0);
    gfs_block_zero -= gf;
    displayln_info(
        ssprintf("CHECK: gf_stout_smear_block with num_step = 0: "
                 "qnorm(diff) = %.2E",
                 qnorm(gfs_block_zero)));
    gf_stout_smear_local(gfs_local_zero, gf, Coordinate(), 0.1, 0);
    gfs_local_zero -= gf;
    displayln_info(
        ssprintf("CHECK: gf_stout_smear_local with num_step = 0: "
                 "qnorm(diff) = %.2E",
                 qnorm(gfs_local_zero)));
    GaugeField gf_after;
    gf_after = gf;
    gf_after -= gf_orig;
    displayln_info(
        ssprintf("CHECK: smear functions leave the input gauge "
                 "field unchanged: qnorm(diff) = %.2E",
                 qnorm(gf_after)));
  }
  {
    TIMER_VERBOSE("test-field_shift");
    GaugeField gf1, gf2;
    const Coordinate c1(1, 2, 3, 4);
    const Coordinate c2(7, 3, 1, 8);
    field_shift(gf1, gf, c1);
    field_shift(gf2, gf1, c2);
    field_shift(gf1, gf, c1 + c2);
    gf1 -= gf2;
    displayln_info(
        ssprintf("CHECK: Consistency: orig qnorm: %.10E ; shift qnorm %.10E ; "
                 "diff qnorm: %.2E",
                 qnorm(gf), qnorm(gf2), qnorm(gf1)));
    field_shift_direct(gf1, gf, c1 + c2);
    gf1 -= gf2;
    displayln_info(
        ssprintf("CHECK: Reference (field_shift_direct): orig qnorm: %.10E ; "
                 "shift qnorm %.10E ; diff qnorm: %.2E",
                 qnorm(gf), qnorm(gf2), qnorm(gf1)));
    field_shift_steps(gf1, gf, c1 + c2);
    gf1 -= gf2;
    displayln_info(
        ssprintf("CHECK: Reference (field_shift_steps): orig qnorm: %.10E ; "
                 "shift qnorm %.10E ; diff qnorm: %.2E",
                 qnorm(gf), qnorm(gf2), qnorm(gf1)));
  }
}

void show_matrix()
{
  const qlat::array<SpinMatrix, 16>& cps_gms =
      SpinMatrixConstants::get_cps_gms();
  const qlat::array<SpinMatrix, 4>& gammas =
      SpinMatrixConstants::get_cps_gammas();
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  for (int i = 0; i < 16; ++i) {
    displayln_info(ssprintf("cps_gms[%d] =\n", i) + show(cps_gms[i]));
  }
  displayln_info(ssprintf("gamma_x * gamma_y * gamma_z * gamma_t =\n") +
                 show(gammas[0] * gammas[1] * gammas[2] * gammas[3]));
  displayln_info(ssprintf("gamma_5 =\n") + show(gamma5));
}

int main(int argc, char* argv[])
{
  begin(&argc, &argv);
  get_global_rng_state() = RngState(get_global_rng_state(), "qcd-utils-tests");
  show_matrix();
  simple_tests();
  displayln_info("CHECK: finished successfully.");
  Timer::display();
  end();
  return 0;
}
