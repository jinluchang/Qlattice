#include <qlat/qlat.h>

#include <complex>
#include <iostream>

using namespace qlat;
using namespace std;

void coordinateHalf(Coordinate& xh, const Coordinate& x)
{
  for (int mu = 0; mu < DIMN; ++mu) {
    xh[mu] = x[mu] / 2;
  }
}

void test1()
{
  TIMER("test1");
  // Coordinate total_site(16, 16, 16, 32);
  Coordinate total_site(4, 4, 4, 8);
  Geometry geo;
  geo.init(total_site, 1);
  Coordinate total_siteHalf;
  coordinateHalf(total_siteHalf, total_site);
  Geometry geoHalf;
  geoHalf.init(total_siteHalf, 1);
  const int seed = 1231;
  const int type = 1;
  const int traj = 1;
  RngState rs(seed);
  splitRngState(rs, rs, type);
  splitRngState(rs, rs, traj);
  RngField rf;
  rf.init(geo, rs);
  FieldM<Complex, 1> af;
  FieldM<double, 1> sumf;
  FieldM<double, 1> sigma2f;
  af.init(geoHalf);
  sumf.init(geoHalf);
  sigma2f.init(geoHalf);
  double gsum = 0.0;
  double gsigma2 = 0.0;
  const int Ni = 2 * 2 * 2 * 2;
  const int Ntake = 1;
  const int Nb = geo.total_volume() / Ni;
  const int Ntraj = 16;
  for (long traj = 0; traj < Ntraj; ++traj) {
    set_zero(af);
    set_zero(sumf);
    set_zero(sigma2f);
    for (long index = 0; index < geo.local_volume(); ++index) {
      Coordinate x = geo.coordinate_from_index(index);
      Coordinate xh;
      coordinateHalf(xh, x);
      RngState& rs = rf.get_elem(x);
      af.get_elem(xh) += polar(1.0, u_rand_gen(rs, PI, -PI));
    }
    for (long index = 0; index < geoHalf.local_volume(); ++index) {
      Coordinate x = geoHalf.coordinate_from_index(index);
      Complex& a = af.get_elem(x);
      sumf.get_elem(x) += qnorm(a);
      sigma2f.get_elem(x) += sqr(qnorm(a));
    }
    const double sum = field_glb_sum_double(sumf)[0];
    const double sigma2 = field_glb_sum_double(sigma2f)[0];
    gsum += sum / Nb;
    gsigma2 += sqr(sum / Nb);
    if (0 == get_id_node()) {
      cout << "traj     : " << traj << endl;
      cout << "Expected : " << Ni * Ntake << endl;
      cout << "Mean     : " << sum / Nb << endl;
      cout << "Var      : " << sqrt(sigma2 / Nb - sqr(sum / Nb)) / sqrt(Nb - 1)
           << endl;
    }
  }
  if (0 == get_id_node()) {
    cout << "# Final" << endl;
    cout << "Expected : " << Ni * Ntake << endl;
    cout << "Mean     : " << gsum / Ntraj << endl;
    cout << "Var      : "
         << sqrt(gsigma2 / Ntraj - sqr(gsum / Ntraj)) / sqrt(Ntraj - 1) << endl;
  }
}

int main(int argc, char* argv[])
{
  begin(&argc, &argv);
  test1();
  displayln_info("CHECK: finished successfully.");
  Timer::display();
  end();
  return 0;
}
