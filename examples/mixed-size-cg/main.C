#include <qlat/qlat.h>

namespace qlat
{  //

inline void set_rand_gauge_field(GaugeField& gf, const Geometry& geo,
                          const RngState& rs)
{
  TIMER_VERBOSE("set_rand_gauge_field");
  gf.init();
  gf.init(geo);
  set_g_rand_color_matrix_field(gf, RngState(rs, fname), 1.0, 10);
  for (int i = 0; i < 40; ++i) {
    gf_ape_smear(gf, gf, 0.1);
  }
}

inline void set_rand_fermion_field(FermionField4d& ff, const Geometry& geo,
                            const RngState& rs)
{
  TIMER_VERBOSE("set_rand_gauge_field");
  ff.init();
  ff.init(geo);
  set_u_rand_double(ff, RngState(rs, fname));
  displayln_info(fname + ssprintf(": qnorm=%24.17E", qnorm(ff)));
}

inline void simple_cg(FermionField4d& ff_sol, const GaugeField& gf,
                      const FermionField4d& ff_src, const FermionAction& fa)
{
  TIMER_VERBOSE("simple_cg");
  InverterDomainWall inv;
  setup_inverter(inv, gf, fa);
  inv.max_mixed_precision_cycle() = 100;
  inv.max_num_iter() = 300;
  inv.stop_rsd() = 1e-10;
  ff_sol.init();
  ff_sol.init(ff_src.geo());
  invert(ff_sol, ff_src, inv);
}

inline void set_smooth_gauge_field(GaugeField& sgf, const GaugeField& gf)
{
  TIMER_VERBOSE("set_smooth_gauge_field");
  GaugeField sgf0;
  sgf0 = gf;
  gf_ape_smear(sgf0, sgf0, 0.6);
  gf_ape_smear(sgf0, sgf0, 0.6);
  const Coordinate expansion_left(2, 2, 2, 2);
  const Coordinate expansion_right(1, 1, 1, 1);
  const Geometry geo1 = geo_resize(sgf0.geo(), expansion_left, expansion_right);
  sgf.init(geo1);
  sgf = sgf0;
  refresh_expanded(sgf); // TODO
}

inline void reduce_half_gauge_field(GaugeField& hgf, const GaugeField& gf)
// xl = coordinate_shifts(hxl * 2, 0)
{
  TIMER_VERBOSE("reduce_half_gauge_field");
  const Geometry& geo = gf.geo();
  Geometry hgeo;
  hgeo.init(geo.geon, geo.node_site / 2, geo.multiplicity);
  const Coordinate expansion_left(1, 1, 1, 1);
  const Coordinate expansion_right(0, 0, 0, 0);
  hgeo = geo_resize(hgeo, expansion_left, expansion_right);
  hgf.init(hgeo);
#pragma omp parallel for
  for (long hindex = 0; hindex < hgeo.local_volume(); ++hindex) {
    const Coordinate hxl = hgeo.coordinate_from_index(hindex);
    const Coordinate xl = coordinate_shifts(hxl * 2, 0);
    for (int m = 0; m < hgeo.multiplicity; ++m) {
      hgf.get_elem(hxl, m) =
          gf.get_elem(xl, m) * gf.get_elem(coordinate_shifts(xl, m), m);
    }
  }
  refresh_expanded_1(hgf);
}

inline FermionAction get_half_fermion_action(const FermionAction& fa)
{
  FermionAction hfa = fa;
  hfa.mass /= 2.0;
  return hfa;
}

struct InverterDomainWallMixedSize : InverterDomainWall {
  InverterDomainWall hinv;
  GaugeField sgf;
  //
  InverterDomainWallMixedSize() { init(); }
  //
  void init() { InverterDomainWall::init(); }
  //
  void setup() { InverterDomainWall::setup(); }
  //
  void setup(const GaugeField& gf_, const FermionAction& fa_)
  {
    InverterDomainWall::setup(gf_, fa_);
    update_hf();
  }
  //
  void setup(const GaugeField& gf_, const FermionAction& fa_, LowModes& lm_)
  {
    InverterDomainWall::setup(gf_, fa_, lm_);
    update_hf();
  }
  //
  bool check_local_volume()
  {
    return geo.node_site % 4 == Coordinate();
  }
  //
  void update_hf()
  {
    qassert(check_local_volume());
    GaugeField hgf;
    FermionAction hfa;
    set_smooth_gauge_field(sgf, gf);
    reduce_half_gauge_field(hgf, sgf);
    hfa = get_half_fermion_action(fa);
    hinv.setup(hgf, hfa);
  }
};

inline void multiply_hermop_sym2(FermionField5d& out, const FermionField5d& in,
                                 const InverterDomainWallMixedSize& inv)
// odd <- odd (works for even <- even as well)
{
  multiply_hermop_sym2(out, in, inv.gf, inv.fa);
}

inline void reduce_half_fermion_field(FermionField5d& hff,
                                   const FermionField5d& ff)
// xl = coordinate_shifts(hxl * 2, 0)
{
  TIMER_VERBOSE("reduce_half_fermion_field");
  const Geometry& geo = ff.geo();
  Geometry hgeo;
  hgeo.init(geo.geon, geo.node_site / 2, geo.multiplicity);
  hgeo.eo = geo.eo;
  hff.init(hgeo);
#pragma omp parallel for
  for (long hindex = 0; hindex < hgeo.local_volume(); ++hindex) {
    const Coordinate hxl = hgeo.coordinate_from_index(hindex);
    const Coordinate xl = coordinate_shifts(hxl * 2, 0);
    Vector<WilsonVector> hv = hff.get_elems(hxl);
    const Vector<WilsonVector> v = ff.get_elems_const(xl);
    for (int m = 0; m < hgeo.multiplicity; ++m) {
      hv[m] = v[m];
    }
  }
  // TODO
}

inline void extend_half_fermion_field(FermionField5d& ff, const FermionField5d& hff)
// xl = coordinate_shifts(hxl * 2, 0)
{
  TIMER_VERBOSE("extend_half_fermion_field");
  const Geometry& hgeo = hff.geo();
  Geometry geo;
  geo.init(hgeo.geon, hgeo.node_site * 2, hgeo.multiplicity);
  geo.eo = hgeo.eo;
  ff.init(geo);
  set_zero(ff);
#pragma omp parallel for
  for (long hindex = 0; hindex < hgeo.local_volume(); ++hindex) {
    const Coordinate hxl = hgeo.coordinate_from_index(hindex);
    const Coordinate xl = coordinate_shifts(hxl * 2, 0);
    const Vector<WilsonVector> hv = hff.get_elems_const(hxl);
    Vector<WilsonVector> v = ff.get_elems(xl);
    for (int m = 0; m < hgeo.multiplicity; ++m) {
      v[m] += (WilsonVector)((Complex)16.0 * hv[m]);
    }
  }
  // TODO
}

inline long cg_with_herm_sym_2_half(FermionField5d& sol,
                                    const FermionField5d& src,
                                    const InverterDomainWallMixedSize& inv,
                                    const double stop_rsd = 1e-8,
                                    const long max_num_iter = 50000)
{
  TIMER_VERBOSE_FLOPS("cg_with_herm_sym_2_half(5d,5d,inv)");
  FermionField5d hsol, hsrc;
  reduce_half_fermion_field(hsrc, src);
  reduce_half_fermion_field(hsol, sol);
  const long half_iter =
      cg_with_f(hsol, hsrc, inv.hinv, multiply_hermop_sym2, stop_rsd, max_num_iter * 16);
  timer.flops += 5500 * half_iter * inv.fa.ls * inv.geo().local_volume() / 16;
  extend_half_fermion_field(sol, hsol);
  // set_zero(sol); // ADJUST ME
  const long iter =
      cg_with_f(sol, src, inv, multiply_hermop_sym2, stop_rsd, max_num_iter);
  timer.flops += 5500 * iter * inv.fa.ls * inv.geo().local_volume();
  return iter;
}

inline long cg_with_herm_sym_2(FermionField5d& sol, const FermionField5d& src,
                               const InverterDomainWallMixedSize& inv,
                               const double stop_rsd = 1e-8,
                               const long max_num_iter = 50000)
{
  TIMER_VERBOSE_FLOPS("cg_with_herm_sym_2(5d,5d,inv)");
  const long iter =
      cg_with_f(sol, src, inv, multiply_hermop_sym2, stop_rsd, max_num_iter);
  timer.flops += 5500 * iter * inv.fa.ls * inv.geo().local_volume();
  return iter;
}

long invert_with_cg_with_guess_half(
    FermionField5d& out, const FermionField5d& in,
    const InverterDomainWallMixedSize& inv,
    long cg(FermionField5d&, const FermionField5d&, const InverterDomainWall&,
            const double, const long),
    double stop_rsd = -1, long max_num_iter = -1,
    long max_mixed_precision_cycle = -1, bool dminus_multiplied_already = false)
{
  TIMER_VERBOSE("invert_with_cg_with_guess_half");
  if (stop_rsd < 0) {
    stop_rsd = inv.stop_rsd();
  }
  if (max_num_iter < 0) {
    max_num_iter = inv.max_num_iter();
  }
  if (max_mixed_precision_cycle < 0) {
    max_mixed_precision_cycle = inv.max_mixed_precision_cycle();
  }
  FermionField5d dm_in;
  if (not dminus_multiplied_already and inv.fa.is_multiplying_dminus) {
    multiply_d_minus(dm_in, in, inv);
  } else {
    dm_in.init(geo_resize(in.geo()));
    dm_in = in;
  }
  const double qnorm_dm_in = qnorm(dm_in);
  FermionField5d tmp;
  multiply_m(tmp, out, inv);
  dm_in -= tmp;
  const double qnorm_dm_in_sub = qnorm(dm_in);
  const double ratio = sqrt(qnorm_dm_in / qnorm_dm_in_sub);
  stop_rsd *= ratio;
  displayln_info(fname + ssprintf(": guess ratio = %.3E", ratio));
  FermionField5d hout, hin;
  reduce_half_fermion_field(hin, dm_in);
  const long total_iter =
      invert_with_cg(hout, hin, inv.hinv, cg, stop_rsd, max_num_iter,
                     max_mixed_precision_cycle, true);
  extend_half_fermion_field(tmp, hout);
  // set_zero(tmp); // ADJUST ME
  out += tmp;
  return total_iter;
}

inline long invert_mix(FermionField5d& out, const FermionField5d& in,
                       const InverterDomainWallMixedSize& inv)
{
  TIMER_VERBOSE("invert_mix_prec(5d,5d,inv-ms)");
  const long total_iter = invert_with_cg(out, in, inv, cg_with_herm_sym_2_half);
  return total_iter;
}

inline long invert_mix_prec(FermionField5d& out, const FermionField5d& in,
                            const InverterDomainWallMixedSize& inv)
{
  TIMER_VERBOSE("invert_mix_prec(5d,5d,inv-ms)");
  long total_iter = 0;
  FermionField5d dm_in;
  if (inv.fa.is_multiplying_dminus) {
    multiply_d_minus(dm_in, in, inv);
  } else {
    dm_in.init(geo_resize(in.geo()));
    dm_in = in;
  }
  total_iter += invert_with_cg(out, dm_in, inv, cg_with_herm_sym_2,
                               inv.stop_rsd(), inv.max_num_iter(), 1, true);
  int cycle;
  for (cycle = 1; cycle <= inv.max_mixed_precision_cycle(); ++cycle) {
    invert_with_cg_with_guess_half(out, dm_in, inv, cg_with_herm_sym_2, inv.stop_rsd(),
                                   inv.max_num_iter() * 8, 1, true);
    const long iter =
        invert_with_cg_with_guess(out, dm_in, inv, cg_with_herm_sym_2,
                                  inv.stop_rsd(), inv.max_num_iter(), 1, true);
    total_iter += iter;
    if (iter < inv.max_num_iter()) {
      break;
    }
  }
  displayln_info(fname + ssprintf(": total_iter=%ld cycle=%d stop_rsd=%.3E",
                                  total_iter, cycle, inv.stop_rsd()));
  return total_iter;
}

inline long invert(FermionField5d& out, const FermionField5d& in,
                   const InverterDomainWallMixedSize& inv)
{
  TIMER_VERBOSE("invert(5d,5d,inv-ms)");
  return invert_mix(out, in, inv);
  // return invert_mix_prec(out, in, inv);
}

inline long invert(FermionField4d& out, const FermionField4d& in,
                    const InverterDomainWallMixedSize& inv)
{
  return invert_dwf(out, in, inv);
}

inline void mixed_size_cg(FermionField4d& ff_sol, const GaugeField& gf,
                          const FermionField4d& ff_src, const FermionAction& fa)
{
  TIMER_VERBOSE("mixed_size_cg");
  InverterDomainWallMixedSize inv;
  setup_inverter(inv, gf, fa);
  inv.max_mixed_precision_cycle() = 100;
  inv.max_num_iter() = 300;
  inv.stop_rsd() = 1e-10;
  ff_sol.init();
  ff_sol.init(ff_src.geo());
  invert(ff_sol, ff_src, inv);
}

}  // namespace qlat

int main(int argc, char* argv[])
{
  using namespace qlat;
  begin(&argc, &argv);
  const Coordinate total_site(4,4,4,8);
  // const Coordinate total_site(8,8,8,16);
  // const Coordinate total_site(16,16,16,32);
  const RngState rs("seed-mixed-size-cg");
  const FermionAction fa(0.01, 16, 1.8, 1.0, true, true);
  Geometry geo;
  geo.init(total_site, 1);
  GaugeField gf;
  set_rand_gauge_field(gf, geo, RngState(rs, "set-field"));
  // set_unit(gf);
  // gf.init(geo);
  // const std::string path = get_env("HOME") + "/qcdarchive/DWF_iwa_nf2p1/16c32/2plus1_16nt32_IWASAKI_b2p13_ls16_M1p8_ms0p04_mu0p01_rhmc_multi_timescale_ukqcd/ckpoint_lat.IEEE64BIG.1000";
  // load_gauge_field(gf, path);
  gf_show_info(gf);
  FermionField4d ff_src, ff_sol, ff_sol_mixed;
  set_rand_fermion_field(ff_src, geo, rs);
  mixed_size_cg(ff_sol_mixed, gf, ff_src, fa);
  simple_cg(ff_sol, gf, ff_src, fa);
  ff_sol_mixed -= ff_sol;
  displayln_info(ssprintf(": sol qnorm=%24.17E", qnorm(ff_sol)));
  displayln_info(ssprintf(": diff qnorm=%24.17E", qnorm(ff_sol_mixed)));
  displayln_info("CHECK: finished successfully.");
  Timer::display();
  end();
  return 0;
}
