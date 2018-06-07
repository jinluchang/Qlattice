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
  const FermionAction fa(0.1, 8, 1.8);
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
  ff5din.init(geo_remult(geo, fa.ls));
  ff5dout.init(geo_remult(geo, fa.ls));
  ff5din1.init(geo_remult(geo, fa.ls));
  ff5dout1.init(geo_remult(geo, fa.ls));
  const Coordinate xg = mod(Coordinate(rand_gen(rs), rand_gen(rs), rand_gen(rs), rand_gen(rs)), geo.total_site());
  const int cs = mod(rand_gen(rs), 12);
  set_fermion_field_point_src(ff4din, xg, cs);
  fermion_field_5d_from_4d(ff5din, ff4din, 0, fa.ls-1);
  ff5din1 = ff5din;
  ff5dout = ff5din;
  ff5dout1 = ff5din1;
  for (int i = 0; i < 20; ++i) {
    multiply_m_dwf(ff5dout1, ff5dout1, inv);
    ff5din = ff5dout;
    multiply_m(ff5dout, ff5din, inv);
    ff5dout *= 1.0 / norm(ff5dout);
    ff5dout1 *= 1.0 / norm(ff5dout1);
    FermionField5d ff5d;
    ff5d = ff5dout;
    ff5d -= ff5dout1;
    displayln_info(ssprintf("norm(diff) = %E", norm(ff5d)));
  }
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
