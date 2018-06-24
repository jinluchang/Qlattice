#include <qlat/qlat.h>

#include <iostream>
#include <complex>

using namespace qlat;
using namespace std;

void simple_tests()
{
  TIMER_VERBOSE("simple_tests");
  RngState rs(get_global_rng_state(), fname);
  // const Coordinate total_site(16, 16, 16, 32);
  // const Coordinate total_site(8, 8, 8, 8);
  const Coordinate total_site(4, 4, 4, 8);
  const FermionAction fa(0.1, 8, 1.8, 2.0);
  Geometry geo;
  geo.init(total_site, 1);
  GaugeField gf;
  gf.init(geo);
  set_unit(gf);
  gf_show_info(gf);
  set_g_rand_color_matrix_field(gf, RngState(rs, "gf-0.1"), 0.1);
  gf_show_info(gf, 1);
  InverterDomainWall inv;
  setup_inverter(inv, gf, fa);
  FermionField4d ff4din;
  ff4din.init(geo);
  FermionField5d ff5din, ff5dout;
  FermionField5d ff5din1, ff5dout1;
  FermionField5d ff5d;
  ff5din.init(geo_remult(geo, fa.ls));
  ff5dout.init(geo_remult(geo, fa.ls));
  ff5din1.init(geo_remult(geo, fa.ls));
  ff5dout1.init(geo_remult(geo, fa.ls));
  const Coordinate xg = mod(Coordinate(rand_gen(rs), rand_gen(rs), rand_gen(rs), rand_gen(rs)), geo.total_site());
  const int cs = mod(rand_gen(rs), 12);
  set_fermion_field_point_src(ff4din, xg, cs);
  fermion_field_5d_from_4d(ff5din, ff4din, 0, fa.ls-1);
  ff5din1 = ff5din;
  for (int i = 0; i < 20; ++i) {
    // project_eo(ff5din, 1);
    // project_eo(ff5din1, 1);
    multiply_m(ff5dout, ff5din, inv);
    // multiply_m_dwf(ff5dout1, ff5din1, inv);
    multiply_m_from_eo(ff5dout1, ff5din1, inv);
    // project_eo(ff5dout, 1);
    // project_eo(ff5dout1, 1);
    ff5d = ff5dout;
    ff5d -= ff5dout1;
    displayln_info(ssprintf("%E, %E norm(diff) = %E", norm(ff5dout), norm(ff5dout1), norm(ff5d)));
    ff5dout *= 1.0 / sqrt(norm(ff5dout));
    ff5dout1 *= 1.0 / sqrt(norm(ff5dout1));
    ff5d = ff5dout;
    ff5d -= ff5dout1;
    displayln_info(ssprintf("norm(diff after normalization) = %E", norm(ff5d)));
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
  displayln_info(ssprintf("norm(diff) = %E", norm(ff5d)));
  displayln_info(ssprintf("norm = %E, dot = %E", norm(ffeven), dot_product(ffeven, ffeven).real()));
  FermionField5d tmp;
  multiply_m_e_e(tmp, ffeven, fa);
  displayln_info(ssprintf("dot = %s", show(dot_product(ffeven, tmp)).c_str()));
  multiply_mdag_e_e(tmp, ffeven, fa);
  displayln_info(ssprintf("dot = %s", show(dot_product(tmp, ffeven)).c_str()));
}

int main(int argc, char* argv[])
{
  Timer::max_function_name_length_shown() = 50;
  begin(&argc, &argv);
  get_global_rng_state() = RngState(get_global_rng_state(), "qcd-utils-tests");
  simple_tests();
  end();
  Timer::display();
  return 0;
}
