#include <qlat/qlat.h>

const char* cname = "Main";

namespace qlat
{  //

void setField(Field<ComplexD>& f)
{
  TIMER("setField");
  const Geometry& geo = f.geo();
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    Coordinate x = geo.coordinate_from_index(index);
    Vector<ComplexD> fx = f.get_elems(x);
    for (int m = 0; m < f.multiplicity; ++m) {
      fx[m] =
          geo.geon.id_node * sqrt(2) + index * sqrt(3) * ii + m * sqrt(5);
    }
  }
}

SpinMatrix projPositiveState(const SpinMatrix& x)
{
  const SpinMatrix psm = (SpinMatrixConstants::get_unit() +
                                SpinMatrixConstants::get_cps_gamma(3)) /
                               (ComplexD)2.0;
  return psm * x * psm;
}

SpinMatrix lblMuonLine(const int tsnk, const int tsrc,
                             const QedGaugeField& egf1,
                             const QedGaugeField& egf2,
                             const QedGaugeField& egf3, const RealD mass,
                             const CoordinateD& momtwist)
{
  TIMER("lblMuonLine");
  const Geometry& geo = egf1.geo();
  SpinPropagator4d sol;
  sol.init(geo);
  SpinPropagator4d src;
  src.init(geo);
  set_zero(src);
  set_wall_source_plusm(src, 1.0, tsrc);
  free_invert(sol, src, mass, 1.0, momtwist);
  //
  set_zero(src);
  sequential_photon_spin_propagator_plusm(src, ii, egf1, sol);
  free_invert(sol, src, mass, 1.0, momtwist);
  //
  set_zero(src);
  sequential_photon_spin_propagator_plusm(src, ii, egf2, sol);
  free_invert(sol, src, mass, 1.0, momtwist);
  //
  set_zero(src);
  sequential_photon_spin_propagator_plusm(src, ii, egf3, sol);
  free_invert(src, src, mass, 1.0, momtwist);
  //
  SpinPropagator4d snk;
  snk.init(geo);
  set_zero(snk);
  set_wall_source_plusm(snk, 1.0, tsnk);
  return contract_spin_propagator4d(snk, src);
}

SpinMatrix lblMuonLineC(const int tsnk, const int tsrc,
                              const QedGaugeField& egf1,
                              const QedGaugeField& egf2,
                              const QedGaugeField& egf3, const RealD mass,
                              const CoordinateD& momtwist)
{
  TIMER("lblMuonLineC");
  SpinMatrix sm;
  set_zero(sm);
  sm += lblMuonLine(tsnk, tsrc, egf1, egf2, egf3, mass, momtwist);
  sm += lblMuonLine(tsnk, tsrc, egf2, egf3, egf1, mass, momtwist);
  sm += lblMuonLine(tsnk, tsrc, egf3, egf1, egf2, mass, momtwist);
  sm += lblMuonLine(tsnk, tsrc, egf3, egf2, egf1, mass, momtwist);
  sm += lblMuonLine(tsnk, tsrc, egf2, egf1, egf3, mass, momtwist);
  sm += lblMuonLine(tsnk, tsrc, egf1, egf3, egf2, mass, momtwist);
  sm *= 1.0 / 6.0;
  return sm;
}

SpinMatrix lblMuonPartPointSrc(const Geometry& geo, const int tsnk,
                                     const int tsrc,
                                     const Coordinate& xg1, const int mu1,
                                     const Coordinate& xg2, const int mu2,
                                     const Coordinate& xg3, const int mu3,
                                     const RealD mass,
                                     const CoordinateD& momtwist)
{
  TIMER("lblMuonPartPointSrc");
  QedGaugeField egf1;
  egf1.init(geo);
  QedGaugeField egf2;
  egf2.init(geo);
  QedGaugeField egf3;
  egf3.init(geo);
  set_zero(egf1);
  set_zero(egf2);
  set_zero(egf3);
  set_point_source_plusm(egf1, 1.0, xg1, mu1);
  set_point_source_plusm(egf2, 1.0, xg2, mu2);
  set_point_source_plusm(egf3, 1.0, xg3, mu3);
  // ADJUST ME
  prop_photon_invert(egf1, momtwist);
  prop_photon_invert(egf2, momtwist);
  prop_photon_invert(egf3, momtwist);
  // ADJUST ME
  // return lblMuonLine(tsnk, tsrc, egf1, egf2, egf3, mass, momtwist);
  return lblMuonLineC(tsnk, tsrc, egf1, egf2, egf3, mass, momtwist);
}

void lblMagneticMomentSpinMatrix(Array<SpinMatrix, 3> bs,
                                 const Geometry& geo, const int tsnk,
                                 const int tsrc, const RealD mass,
                                 const CoordinateD& momtwist)
// pretend to the operator to be
// \Sigma_i * mass / 2
{
  TIMER("lblMagneticMomentSpinMatrix");
  SpinPropagator4d snk;
  snk.init(geo);
  SpinPropagator4d src;
  src.init(geo);
  set_zero(snk);
  set_zero(src);
  set_wall_source_plusm(snk, 1.0, tsnk);
  set_wall_source_plusm(src, 1.0, tsrc);
  free_invert(snk, snk, mass, 1.0, momtwist);
  free_invert(src, src, mass, 1.0, momtwist);
  const int top =
      mod(tsrc + mod(tsnk - tsrc, geo.total_site()[3]) / 2,
                geo.total_site()[3]);
  Coordinate xgop(0, 0, 0, top);
  Coordinate xlop = geo.coordinate_l_from_g(xgop);
  set_zero(bs);
  if (geo.is_local(xlop)) {
    Display(cname, fname.c_str(), "src =\n%s\n",
                  show(src.get_elem(xlop)).c_str());
    Display(cname, fname.c_str(), "snk =\n%s\n",
                  show(snk.get_elem(xlop)).c_str());
    for (int i = 0; i < 3; ++i) {
      bs[i] = SpinMatrixConstants::get_gamma5() *
              matrix_adjoint(snk.get_elem(xlop)) *
              SpinMatrixConstants::get_gamma5();
      bs[i] *= SpinMatrixConstants::get_cps_cap_sigmas()[i] *
               (ComplexD)(mass / 2.0);
      bs[i] *= src.get_elem(xlop);
      bs[i] = projPositiveState(bs[i]);
    }
  }
  glb_sum(Vector<RealD>((RealD*)bs.data(),
                                    get_data_size(bs) / sizeof(RealD)));
}

ComplexD linearFit(const SpinMatrix& x,
                         const SpinMatrix& base)
{
  using namespace qlat;
  const int size = 4 * 4;
  ComplexD s1 = 0.0;
  ComplexD s2 = 0.0;
  for (int i = 0; i < size; ++i) {
    s1 += qconj(base.data()[i]) * x.data()[i];
    s2 += qconj(base.data()[i]) * base.data()[i];
  }
  return s1 / s2;
}

void lblShowMuonPartPointSrc(const Geometry& geo, const int tsnk,
                             const int tsrc,
                             const Array<SpinMatrix, 3>& bs,
                             const Coordinate& xg1, const int mu1,
                             const Coordinate& xg2, const int mu2,
                             const Coordinate& xg3, const int mu3,
                             const RealD mass,
                             const CoordinateD& momtwist)
{
  TIMER("lblShowMuonPartPointSrc");
  DisplayInfo(cname, fname.c_str(), "mass = %.2f\n", mass);
  DisplayInfo(cname, fname.c_str(), "xg1 = %s ; xg2 = %s ; xg3 = %s .\n",
                    show(xg1).c_str(), show(xg2).c_str(),
                    show(xg3).c_str());
  DisplayInfo(cname, fname.c_str(), "mu1 = %d ; mu2 = %d ; mu3 = %d .\n",
                    mu1, mu2, mu3);
  SpinMatrix muonline = lblMuonPartPointSrc(
      geo, tsnk, tsrc, xg1, mu1, xg2, mu2, xg3, mu3, mass, momtwist);
  muonline = projPositiveState(muonline);
  DisplayInfo(cname, fname.c_str(), "qnorm(muonline) = %.16e\n",
                    qnorm(muonline));
  if (qnorm(muonline) < 1.0e-30) {
    return;
  }
  DisplayInfo(cname, fname.c_str(), "muonline =\n%s\n",
                    show(muonline).c_str());
  DisplayInfo(cname, fname.c_str(), "linearFit[0] = %s\n",
                    show(linearFit(muonline, bs[0])).c_str());
  DisplayInfo(cname, fname.c_str(), "linearFit[1] = %s\n",
                    show(linearFit(muonline, bs[1])).c_str());
  DisplayInfo(cname, fname.c_str(), "linearFit[2] = %s\n",
                    show(linearFit(muonline, bs[2])).c_str());
  DisplayInfo(
      cname, fname.c_str(), "linearFitUni = %s\n",
      show(
          linearFit(muonline,
                    projPositiveState(SpinMatrixConstants::get_unit())))
          .c_str());
  displayln_info("CHECK: " + fname +
                       ssprintf(": linearFit[0] * 1e9 = %10.6f",
                                      1e9 * linearFit(muonline, bs[0]).real()));
  displayln_info("CHECK: " + fname +
                       ssprintf(": linearFit[1] * 1e9 = %10.6f",
                                      1e9 * linearFit(muonline, bs[1]).real()));
  displayln_info("CHECK: " + fname +
                       ssprintf(": linearFit[2] * 1e9 = %10.6f",
                                      1e9 * linearFit(muonline, bs[2]).real()));
  displayln_info(
      "CHECK: " + fname +
      ssprintf(
          ": linearFitUni * 1e9 = %10.6f",
          1e9 * linearFit(muonline, projPositiveState(
                                        SpinMatrixConstants::get_unit()))
                    .imag()));
}

void lblMuonPart()
{
  TIMER("lblMuonPart");
  // ADJUST ME
  // Coordinate total_site(16, 16, 16, 32);
  Coordinate total_site(4, 4, 4, 8);
  // Coordinate total_site(32, 32, 32, 128);
  Geometry geo;
  geo.init(total_site);
  DisplayInfo(cname, fname.c_str(), "geo =\n%s\n",
                    show(geo).c_str());
  CoordinateD momtwist;
  momtwist[0] = 0.0;
  momtwist[1] = 0.0;
  momtwist[2] = 0.0;
  momtwist[3] = 0.0;
  const RealD mass = 0.10;
  const int tsnk = total_site[3] / 4 * 3;
  const int tsrc = total_site[3] / 4;
  array<SpinMatrix, 3> bs;
  lblMagneticMomentSpinMatrix(bs, geo, tsnk, tsrc, mass, momtwist);
  DisplayInfo(cname, fname.c_str(), "bs[0] =\n%s\n",
                    show(bs[0]).c_str());
  DisplayInfo(cname, fname.c_str(), "bs[1] =\n%s\n",
                    show(bs[1]).c_str());
  DisplayInfo(cname, fname.c_str(), "bs[2] =\n%s\n",
                    show(bs[2]).c_str());
  // ADJUST ME
  // Coordinate xg1(0, 0, 0, total_site[3]/2 + 0);
  // Coordinate xg2(0, 0, 0, total_site[3]/2 + 0);
  // Coordinate xg3(0, 0, 0, total_site[3]/2 + 0);
  //
  Coordinate xg1(1, 1, 1, total_site[3] / 2 + 2);
  Coordinate xg2(2, 3, 1, total_site[3] / 2 - 1);
  Coordinate xg3(2, 2, 2, total_site[3] / 2 + 0);
  //
  // Coordinate xg1(2, 1, 5, total_site[3]/2 - 4);
  // Coordinate xg2(3, 4, 0, total_site[3]/2 + 0);
  // Coordinate xg3(1, 5, 2, total_site[3]/2 + 4);
  //
  for (int mu1 = 0; mu1 < 4; ++mu1) {
    for (int mu2 = 0; mu2 < 4; ++mu2) {
      for (int mu3 = 0; mu3 < 4; ++mu3) {
        lblShowMuonPartPointSrc(geo, tsnk, tsrc, bs, xg1, mu1, xg2, mu2, xg3,
                                mu3, mass, momtwist);
        Timer::display();
        return;
      }
    }
  }
}

void displayGammas()
{
  TIMER("displayGammas");
  DisplayInfo(
      cname, fname.c_str(), "gamma5 =\n%s\n",
      show(SpinMatrixConstants::get_gamma5()).c_str());
  DisplayInfo(
      cname, fname.c_str(), "gamma0 * gamma1 * gamma2 * gamma3 =\n%s\n",
      show(
          (SpinMatrix)(SpinMatrixConstants::get_cps_gamma(0) *
                             SpinMatrixConstants::get_cps_gamma(1) *
                             SpinMatrixConstants::get_cps_gamma(2) *
                             SpinMatrixConstants::get_cps_gamma(3)))
          .c_str());
  DisplayInfo(
      cname, fname.c_str(), "gamma3 =\n%s\n",
      show(SpinMatrixConstants::get_cps_gamma(3)).c_str());
  DisplayInfo(cname, fname.c_str(), "ii * gamma0 =\n%s\n",
                    show((SpinMatrix)(
                                   (ComplexD)ii *
                                   SpinMatrixConstants::get_cps_gamma(0)))
                        .c_str());
  DisplayInfo(cname, fname.c_str(), "ii * gamma1 =\n%s\n",
                    show((SpinMatrix)(
                                   (ComplexD)ii *
                                   SpinMatrixConstants::get_cps_gamma(1)))
                        .c_str());
  DisplayInfo(cname, fname.c_str(), "ii * gamma2 =\n%s\n",
                    show((SpinMatrix)(
                                   (ComplexD)ii *
                                   SpinMatrixConstants::get_cps_gamma(2)))
                        .c_str());
}

void displaySpinPropagator4d()
{
  TIMER("displaySpinPropagator4d");
  // Coordinate total_site(16, 16, 16, 32);
  Coordinate total_site(4, 4, 4, 8);
  Geometry geo;
  geo.init(total_site);
  DisplayInfo(cname, fname.c_str(), "geo =\n%s\n",
                    show(geo).c_str());
  CoordinateD momtwist;
  momtwist[0] = 0.0;
  momtwist[1] = 0.0;
  momtwist[2] = 0.0;
  momtwist[3] = 0.0;
  const RealD mass = 0.1;
  SpinPropagator4d prop;
  prop.init(geo);
  set_zero(prop);
  Coordinate xgsrc(0, 0, 0, 0);
  Coordinate xlsrc = geo.coordinate_l_from_g(xgsrc);
  if (geo.is_local(xlsrc)) {
    set_unit(prop.get_elem(xlsrc));
  }
  free_invert(prop, prop, mass, 1.0, momtwist);
  Coordinate xgsnk(0, 0, 0, 0);
  Coordinate xlsnk = geo.coordinate_l_from_g(xgsnk);
  DisplayInfo(cname, fname.c_str(), "xgsnk = %s .\n",
                    show(xgsnk).c_str());
  if (geo.is_local(xlsnk)) {
    Display(cname, fname.c_str(), "prop[xgsnk] =\n%s\n",
                  show(prop.get_elem(xlsnk)).c_str());
  }
}

}  // namespace qlat

int main(int argc, char* argv[])
{
  using namespace qlat;
  begin(&argc, &argv);
  displayGammas();
  displaySpinPropagator4d();
  Timer::display();
  lblMuonPart();
  displayln_info("CHECK: finished successfully.");
  Timer::display();
  end();
  return 0;
}
