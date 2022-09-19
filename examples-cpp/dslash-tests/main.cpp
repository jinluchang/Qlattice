#include <qlat/qlat.h>

#include <complex>
#include <iostream>

using namespace qlat;
using namespace std;

void simple_dwf_tests()
{
  TIMER_VERBOSE("simple_dwf_tests");
  RngState rs(get_global_rng_state(), fname);
  // const Coordinate total_site(16, 16, 16, 32);
  // const Coordinate total_site(8, 8, 8, 16);
  const Coordinate total_site(4, 4, 4, 8);
  FermionAction fa(0.1, 12, 1.8);
  Geometry geo;
  geo.init(total_site, 1);
  GaugeField gf;
  gf.init(geo);
  set_unit(gf);
  gf_show_info(gf);
  set_g_rand_color_matrix_field(gf, RngState(rs, "gf-0.1"), 0.1);
  gf_show_info(gf, 1);
  InverterDomainWall inv;
  inv.stop_rsd() = 1e-12;
  setup_inverter(inv, gf, fa);
  FermionField4d ff4din;
  ff4din.init(geo);
  set_zero(ff4din);
  FermionField5d ff5din, ff5dout;
  FermionField5d ff5din1, ff5dout1;
  FermionField5d ff5d;
  ff5din.init(geo_remult(geo, fa.ls));
  ff5dout.init(geo_remult(geo, fa.ls));
  ff5din1.init(geo_remult(geo, fa.ls));
  ff5dout1.init(geo_remult(geo, fa.ls));
  set_zero(ff5din);
  set_zero(ff5din1);
  set_zero(ff5dout);
  set_zero(ff5dout1);
  const Coordinate xg =
      mod(Coordinate(rand_gen(rs), rand_gen(rs), rand_gen(rs), rand_gen(rs)),
          geo.total_site());
  const int cs = mod(rand_gen(rs), 12);
  set_point_src_fermion_field(ff4din, xg, cs);
  displayln_info("testing multiply_m_full and multiply_m");
  fermion_field_5d_from_4d(ff5din, ff4din, 0, fa.ls - 1);
  ff5din1 = ff5din;
  for (int i = 0; i < 20; ++i) {
    // project_eo(ff5din, 1);
    // project_eo(ff5din1, 1);
    multiply_m_full(ff5dout, ff5din, inv);
    // multiply_m(ff5dout1, ff5din1, inv);
    // multiply_m_with_prec_sym2(ff5dout1, ff5din1, inv);
    multiply_m_dwf(ff5dout1, ff5din1, inv);
    // project_eo(ff5dout, 1);
    // project_eo(ff5dout1, 1);
    ff5d = ff5dout;
    ff5d -= ff5dout1;
    displayln_info(ssprintf("%E, %E qnorm(diff) = %E", qnorm(ff5dout),
                            qnorm(ff5dout1), qnorm(ff5d)));
    ff5dout *= 1.0 / sqrt(qnorm(ff5dout));
    ff5dout1 *= 1.0 / sqrt(qnorm(ff5dout1));
    ff5d = ff5dout;
    ff5d -= ff5dout1;
    displayln_info(ssprintf("qnorm(diff after qnormalization) = %E", qnorm(ff5d)));
    ff5din = ff5dout;
    ff5din1 = ff5dout1;
  }
}

void simple_tests()
{
  TIMER_VERBOSE("simple_tests");
  RngState rs(get_global_rng_state(), fname);
  // const Coordinate total_site(16, 16, 16, 32);
  // const Coordinate total_site(8, 8, 8, 16);
  const Coordinate total_site(4, 4, 4, 8);
  FermionAction fa(0.1, 12, 1.8, 2.5 + 1.5);
  fa.is_using_zmobius = true;
  if (true) {
    std::vector<Complex> omega(12, 0);
    omega[0] = 1.0903256131299373;
    omega[1] = 0.9570283702230611;
    omega[2] = 0.7048886040934104;
    omega[3] = 0.48979921782791747;
    omega[4] = 0.328608311201356;
    omega[5] = 0.21664245377015995;
    omega[6] = 0.14121112711957107;
    omega[7] = 0.0907785101745156;
    omega[8] = Complex(0.05608303440064219, -0.007537158177840385);
    omega[9] = Complex(0.05608303440064219, 0.007537158177840385);
    omega[10] = Complex(0.0365221637144842, -0.03343945161367745);
    omega[11] = Complex(0.0365221637144842, 0.03343945161367745);
    for (int i = 0; i < fa.ls; ++i) {
      fa.bs[i] = 0.5 * (1.0 / omega[i] + 1.0);
      fa.cs[i] = fa.bs[i] - 1.0;
    }
  }
  Geometry geo;
  geo.init(total_site, 1);
  GaugeField gf;
  gf.init(geo);
  set_unit(gf);
  gf_show_info(gf);
  set_g_rand_color_matrix_field(gf, RngState(rs, "gf-0.1"), 0.1);
  gf_show_info(gf, 1);
  InverterDomainWall inv;
  inv.stop_rsd() = 1e-12;
  setup_inverter(inv, gf, fa);
  find_max_eigen_value_hermop_sym2(
      inv, RngState(rs, "find_max_eigen_value_hermop_sym2"));
  FermionField4d ff4din;
  ff4din.init(geo);
  set_zero(ff4din);
  FermionField5d ff5din, ff5dout;
  FermionField5d ff5din1, ff5dout1;
  FermionField5d ff5d;
  ff5din.init(geo_remult(geo, fa.ls));
  ff5dout.init(geo_remult(geo, fa.ls));
  ff5din1.init(geo_remult(geo, fa.ls));
  ff5dout1.init(geo_remult(geo, fa.ls));
  set_zero(ff5din);
  set_zero(ff5din1);
  set_zero(ff5dout);
  set_zero(ff5dout1);
  const Coordinate xg =
      mod(Coordinate(rand_gen(rs), rand_gen(rs), rand_gen(rs), rand_gen(rs)),
          geo.total_site());
  const int cs = mod(rand_gen(rs), 12);
  set_point_src_fermion_field(ff4din, xg, cs);
  displayln_info("testing multiply_m_full and multiply_m");
  fermion_field_5d_from_4d(ff5din, ff4din, 0, fa.ls - 1);
  ff5din1 = ff5din;
  for (int i = 0; i < 20; ++i) {
    // project_eo(ff5din, 1);
    // project_eo(ff5din1, 1);
    multiply_m_full(ff5dout, ff5din, inv);
    multiply_m(ff5dout1, ff5din1, inv);
    // multiply_m_with_prec_sym2(ff5dout1, ff5din1, inv);
    // multiply_m_dwf(ff5dout1, ff5din1, inv);
    // project_eo(ff5dout, 1);
    // project_eo(ff5dout1, 1);
    ff5d = ff5dout;
    ff5d -= ff5dout1;
    displayln_info(ssprintf("%E, %E qnorm(diff) = %E", qnorm(ff5dout),
                            qnorm(ff5dout1), qnorm(ff5d)));
    ff5dout *= 1.0 / sqrt(qnorm(ff5dout));
    ff5dout1 *= 1.0 / sqrt(qnorm(ff5dout1));
    ff5d = ff5dout;
    ff5d -= ff5dout1;
    displayln_info(ssprintf("qnorm(diff after qnormalization) = %E", qnorm(ff5d)));
    ff5din = ff5dout;
    ff5din1 = ff5dout1;
  }
  displayln_info("testing multiply_m_full and multiply_m_with_prec_sym2");
  fermion_field_5d_from_4d(ff5din, ff4din, 0, fa.ls - 1);
  ff5din1 = ff5din;
  for (int i = 0; i < 20; ++i) {
    // project_eo(ff5din, 1);
    // project_eo(ff5din1, 1);
    multiply_m_full(ff5dout, ff5din, inv);
    // multiply_m(ff5dout1, ff5din1, inv);
    multiply_m_with_prec_sym2(ff5dout1, ff5din1, inv);
    // multiply_m_dwf(ff5dout1, ff5din1, inv);
    // project_eo(ff5dout, 1);
    // project_eo(ff5dout1, 1);
    ff5d = ff5dout;
    ff5d -= ff5dout1;
    displayln_info(ssprintf("%E, %E qnorm(diff) = %E", qnorm(ff5dout),
                            qnorm(ff5dout1), qnorm(ff5d)));
    ff5dout *= 1.0 / sqrt(qnorm(ff5dout));
    ff5dout1 *= 1.0 / sqrt(qnorm(ff5dout1));
    ff5d = ff5dout;
    ff5d -= ff5dout1;
    displayln_info(ssprintf("qnorm(diff after qnormalization) = %E", qnorm(ff5d)));
    ff5din = ff5dout;
    ff5din1 = ff5dout1;
  }
  FermionField5d ffeven, ffodd;
  get_half_fermion(ffeven, ff5dout, 2);
  get_half_fermion(ffodd, ff5dout, 1);
  multiply_m_e_e(ffeven, ffeven, fa);
  multiply_m_e_e_inv(ffeven, ffeven, fa);
  multiply_m_e_e(ffodd, ffodd, fa);
  multiply_m_e_e_inv(ffodd, ffodd, fa);
  set_half_fermion(ff5d, ffeven, 2);
  set_half_fermion(ff5d, ffodd, 1);
  ff5d -= ff5dout;
  displayln_info(ssprintf("m_e_e qnorm(diff) = %E", qnorm(ff5d)));
  get_half_fermion(ffeven, ff5dout, 2);
  get_half_fermion(ffodd, ff5dout, 1);
  multiply_mdag_e_e(ffeven, ffeven, fa);
  multiply_mdag_e_e_inv(ffeven, ffeven, fa);
  multiply_mdag_e_e(ffodd, ffodd, fa);
  multiply_mdag_e_e_inv(ffodd, ffodd, fa);
  set_half_fermion(ff5d, ffeven, 2);
  set_half_fermion(ff5d, ffodd, 1);
  ff5d -= ff5dout;
  displayln_info(ssprintf("mdag_e_e qnorm(diff) = %E", qnorm(ff5d)));
  displayln_info(ssprintf("qnorm = %E, dot = %E", qnorm(ffeven),
                          dot_product(ffeven, ffeven).real()));
  FermionField5d tmp;
  multiply_m_e_e(tmp, ffeven, fa);
  displayln_info(ssprintf("dot = %s", show(dot_product(ffeven, tmp)).c_str()));
  multiply_mdag_e_e(tmp, ffeven, fa);
  displayln_info(ssprintf("dot = %s", show(dot_product(tmp, ffeven)).c_str()));
  multiply_m_e_e_inv(tmp, ffeven, fa);
  displayln_info(ssprintf("dot = %s", show(dot_product(ffeven, tmp)).c_str()));
  multiply_mdag_e_e_inv(tmp, ffeven, fa);
  displayln_info(ssprintf("dot = %s", show(dot_product(tmp, ffeven)).c_str()));
  tmp.init();
  multiply_m(tmp, ff5dout, inv);
  displayln_info(ssprintf("dot = %s", show(dot_product(tmp, ff5dout)).c_str()));
  multiply_mdag(tmp, ff5dout, inv);
  displayln_info(ssprintf("dot = %s", show(dot_product(ff5dout, tmp)).c_str()));
  //
  displayln_info(ssprintf("qnorm = %E, dot = %E", qnorm(ffodd),
                          dot_product(ffodd, ffodd).real()));
  tmp.init();
  multiply_hermop_sym2(tmp, ffodd, inv);
  FermionField5d ffodd2;
  displayln_info(ssprintf("multiply qnorm = %E", qnorm(tmp)));
  cg_with_f(ffodd2, tmp, inv, multiply_hermop_sym2);
  displayln_info(ssprintf("qnorm = %E", qnorm(ffodd2)));
  ffodd2 -= ffodd;
  displayln_info(ssprintf("diff qnorm = %E", qnorm(ffodd2)));
  //
  tmp.init();
  displayln_info(ssprintf("orig qnorm = %E", qnorm(ff5dout)));
  multiply_m(tmp, ff5dout, inv);
  displayln_info(ssprintf("tmp qnorm = %E", qnorm(tmp)));
  FermionField5d sol;
  invert(sol, tmp, inv);
  multiply_d_minus(tmp, tmp, inv);
  multiply_m(sol, sol, inv);
  sol -= tmp;
  displayln_info(ssprintf("invert diff qnorm = %E", qnorm(sol)));
}

int main(int argc, char* argv[])
{
  begin(&argc, &argv);
  get_global_rng_state() = RngState(get_global_rng_state(), "qcd-utils-tests");
  simple_dwf_tests();
  simple_tests();
  displayln_info("CHECK: finished successfully.");
  Timer::display();
  end();
  return 0;
}
