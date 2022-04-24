#include <qlat/qlat.h>

#include <iostream>

const char* cname = "Main";

void setField(qlat::Field<qlat::Complex>& f)
{
  TIMER("setField");
  const qlat::Geometry& geo = f.geo();
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    qlat::Coordinate x = geo.coordinate_from_index(index);
    qlat::Vector<qlat::Complex> fx = f.get_elems(x);
    for (int m = 0; m < geo.multiplicity; ++m) {
      fx[m] =
          geo.geon.id_node * sqrt(2) + index * sqrt(3) * qlat::ii + m * sqrt(5);
    }
  }
}

qlat::SpinMatrix projPositiveState(const qlat::SpinMatrix& x)
{
  const qlat::SpinMatrix psm = (qlat::SpinMatrixConstants::get_unit() +
                                qlat::SpinMatrixConstants::get_gamma(3)) /
                               (qlat::ComplexT)2.0;
  return psm * x * psm;
}

qlat::SpinMatrix lblMuonLine(const int tsnk, const int tsrc,
                             const qlat::QedGaugeField& egf1,
                             const qlat::QedGaugeField& egf2,
                             const qlat::QedGaugeField& egf3, const double mass,
                             const qlat::array<double, qlat::DIMN>& momtwist)
{
  TIMER("lblMuonLine");
  const qlat::Geometry& geo = egf1.geo();
  qlat::SpinPropagator4d sol;
  sol.init(geo);
  qlat::SpinPropagator4d src;
  src.init(geo);
  qlat::set_zero(src);
  qlat::set_wall_source_plusm(src, 1.0, tsrc);
  qlat::prop_spin_propagator4d(src, mass, momtwist);
  //
  sol = src;
  qlat::set_zero(src);
  qlat::sequential_photon_spin_propagator_plusm(src, qlat::ii, egf1, sol);
  qlat::prop_spin_propagator4d(src, mass, momtwist);
  //
  sol = src;
  qlat::set_zero(src);
  qlat::sequential_photon_spin_propagator_plusm(src, qlat::ii, egf2, sol);
  qlat::prop_spin_propagator4d(src, mass, momtwist);
  //
  sol = src;
  qlat::set_zero(src);
  qlat::sequential_photon_spin_propagator_plusm(src, qlat::ii, egf3, sol);
  qlat::prop_spin_propagator4d(src, mass, momtwist);
  //
  qlat::SpinPropagator4d snk;
  snk.init(geo);
  qlat::set_zero(snk);
  qlat::set_wall_source_plusm(snk, 1.0, tsnk);
  return qlat::contract_spin_propagator4d(snk, src);
}

qlat::SpinMatrix lblMuonLineC(const int tsnk, const int tsrc,
                              const qlat::QedGaugeField& egf1,
                              const qlat::QedGaugeField& egf2,
                              const qlat::QedGaugeField& egf3,
                              const double mass,
                              const qlat::array<double, qlat::DIMN>& momtwist)
{
  TIMER("lblMuonLineC");
  qlat::SpinMatrix sm;
  qlat::set_zero(sm);
  sm += lblMuonLine(tsnk, tsrc, egf1, egf2, egf3, mass, momtwist);
  sm += lblMuonLine(tsnk, tsrc, egf2, egf3, egf1, mass, momtwist);
  sm += lblMuonLine(tsnk, tsrc, egf3, egf1, egf2, mass, momtwist);
  sm += lblMuonLine(tsnk, tsrc, egf3, egf2, egf1, mass, momtwist);
  sm += lblMuonLine(tsnk, tsrc, egf2, egf1, egf3, mass, momtwist);
  sm += lblMuonLine(tsnk, tsrc, egf1, egf3, egf2, mass, momtwist);
  sm *= 1.0 / 6.0;
  return sm;
}

qlat::SpinMatrix lblMuonPartPointSrc(
    const qlat::Geometry& geo, const int tsnk, const int tsrc,
    const qlat::Coordinate& xg1, const int mu1, const qlat::Coordinate& xg2,
    const int mu2, const qlat::Coordinate& xg3, const int mu3,
    const double mass, const qlat::array<double, qlat::DIMN>& momtwist)
{
  TIMER("lblMuonPartPointSrc");
  qlat::QedGaugeField egf1;
  egf1.init(geo);
  qlat::QedGaugeField egf2;
  egf2.init(geo);
  qlat::QedGaugeField egf3;
  egf3.init(geo);
  qlat::set_zero(egf1);
  qlat::set_zero(egf2);
  qlat::set_zero(egf3);
  qlat::set_point_source_plusm(egf1, 1.0, xg1, mu1);
  qlat::set_point_source_plusm(egf2, 1.0, xg2, mu2);
  qlat::set_point_source_plusm(egf3, 1.0, xg3, mu3);
  // ADJUST ME
  qlat::prop_photon_invert(egf1, momtwist);
  qlat::prop_photon_invert(egf2, momtwist);
  qlat::prop_photon_invert(egf3, momtwist);
  // ADJUST ME
  // return lblMuonLine(tsnk, tsrc, egf1, egf2, egf3, mass, momtwist);
  return lblMuonLineC(tsnk, tsrc, egf1, egf2, egf3, mass, momtwist);
}

void lblMagneticMomentSpinMatrix(qlat::Array<qlat::SpinMatrix, 3> bs,
                                 const qlat::Geometry& geo, const int tsnk,
                                 const int tsrc, const double mass,
                                 const qlat::array<double, qlat::DIMN>& momtwist)
// pretend to the operator to be
// \Sigma_i * mass / 2
{
  TIMER("lblMagneticMomentSpinMatrix");
  qlat::SpinPropagator4d snk;
  snk.init(geo);
  qlat::SpinPropagator4d src;
  src.init(geo);
  qlat::set_zero(snk);
  qlat::set_zero(src);
  qlat::set_wall_source_plusm(snk, 1.0, tsnk);
  qlat::set_wall_source_plusm(src, 1.0, tsrc);
  qlat::prop_spin_propagator4d(snk, mass, momtwist);
  qlat::prop_spin_propagator4d(src, mass, momtwist);
  const int top =
      qlat::mod(tsrc + qlat::mod(tsnk - tsrc, geo.total_site()[3]) / 2,
                geo.total_site()[3]);
  qlat::Coordinate xgop(0, 0, 0, top);
  qlat::Coordinate xlop = geo.coordinate_l_from_g(xgop);
  qlat::set_zero(bs);
  if (geo.is_local(xlop)) {
    qlat::Display(cname, fname.c_str(), "src =\n%s\n",
            qlat::show(src.get_elem(xlop)).c_str());
    qlat::Display(cname, fname.c_str(), "snk =\n%s\n",
            qlat::show(snk.get_elem(xlop)).c_str());
    for (int i = 0; i < 3; ++i) {
      bs[i] = qlat::SpinMatrixConstants::get_gamma5() *
              matrix_adjoint(snk.get_elem(xlop)) *
              qlat::SpinMatrixConstants::get_gamma5();
      bs[i] *= qlat::SpinMatrixConstants::get_cap_sigma(i) * (qlat::ComplexT)(mass / 2.0);
      bs[i] *= src.get_elem(xlop);
      bs[i] = projPositiveState(bs[i]);
    }
  }
  qlat::glb_sum(qlat::Vector<double>((double*)bs.data(),
                                     get_data_size(bs) / sizeof(double)));
}

qlat::Complex linearFit(const qlat::SpinMatrix& x, const qlat::SpinMatrix& base)
{
  return (Eigen::conj(base.em().array()) * x.em().array()).sum() /
         (Eigen::conj(base.em().array()) * base.em().array()).sum();
}

void lblShowMuonPartPointSrc(const qlat::Geometry& geo, const int tsnk,
                             const int tsrc,
                             const qlat::Array<qlat::SpinMatrix, 3>& bs,
                             const qlat::Coordinate& xg1, const int mu1,
                             const qlat::Coordinate& xg2, const int mu2,
                             const qlat::Coordinate& xg3, const int mu3,
                             const double mass,
                             const qlat::array<double, qlat::DIMN>& momtwist)
{
  TIMER("lblShowMuonPartPointSrc");
  qlat::DisplayInfo(cname, fname.c_str(), "mass = %.2f\n", mass);
  qlat::DisplayInfo(cname, fname.c_str(), "xg1 = %s ; xg2 = %s ; xg3 = %s .\n",
                    qlat::show(xg1).c_str(), qlat::show(xg2).c_str(),
                    qlat::show(xg3).c_str());
  qlat::DisplayInfo(cname, fname.c_str(), "mu1 = %d ; mu2 = %d ; mu3 = %d .\n", mu1,
                    mu2, mu3);
  qlat::SpinMatrix muonline = lblMuonPartPointSrc(
      geo, tsnk, tsrc, xg1, mu1, xg2, mu2, xg3, mu3, mass, momtwist);
  muonline = projPositiveState(muonline);
  qlat::DisplayInfo(cname, fname.c_str(), "qnorm(muonline) = %.16e\n",
                    qlat::qnorm(muonline));
  if (qlat::qnorm(muonline) < 1.0e-30) {
    return;
  }
  qlat::DisplayInfo(cname, fname.c_str(), "muonline =\n%s\n",
                    qlat::show(muonline).c_str());
  qlat::DisplayInfo(cname, fname.c_str(), "linearFit[0] = %s\n",
                    qlat::show(linearFit(muonline, bs[0])).c_str());
  qlat::DisplayInfo(cname, fname.c_str(), "linearFit[1] = %s\n",
                    qlat::show(linearFit(muonline, bs[1])).c_str());
  qlat::DisplayInfo(cname, fname.c_str(), "linearFit[2] = %s\n",
                    qlat::show(linearFit(muonline, bs[2])).c_str());
  qlat::DisplayInfo(
      cname, fname.c_str(), "linearFitUni = %s\n",
      qlat::show(
          linearFit(muonline,
                    projPositiveState(qlat::SpinMatrixConstants::get_unit())))
          .c_str());
  qlat::DisplayInfo(cname, fname.c_str(), "linearFit[0] * 1e9 = %10.2f\n",
                    1e9 * linearFit(muonline, bs[0]).real());
  qlat::DisplayInfo(cname, fname.c_str(), "linearFit[1] * 1e9 = %10.2f\n",
                    1e9 * linearFit(muonline, bs[1]).real());
  qlat::DisplayInfo(cname, fname.c_str(), "linearFit[2] * 1e9 = %10.2f\n",
                    1e9 * linearFit(muonline, bs[2]).real());
}

void lblMuonPart()
{
  TIMER("lblMuonPart");
  // ADJUST ME
  // qlat::Coordinate total_site(16, 16, 16, 32);
  qlat::Coordinate total_site(4, 4, 4, 8);
  // qlat::Coordinate total_site(32, 32, 32, 128);
  qlat::Geometry geo;
  geo.init(total_site, 1);
  qlat::DisplayInfo(cname, fname.c_str(), "geo =\n%s\n", qlat::show(geo).c_str());
  qlat::array<double, qlat::DIMN> momtwist;
  momtwist[0] = 0.0;
  momtwist[1] = 0.0;
  momtwist[2] = 0.0;
  momtwist[3] = 0.0;
  const double mass = 0.10;
  const int tsnk = total_site[3] / 4 * 3;
  const int tsrc = total_site[3] / 4;
  qlat::array<qlat::SpinMatrix, 3> bs;
  lblMagneticMomentSpinMatrix(bs, geo, tsnk, tsrc, mass, momtwist);
  qlat::DisplayInfo(cname, fname.c_str(), "bs[0] =\n%s\n", qlat::show(bs[0]).c_str());
  qlat::DisplayInfo(cname, fname.c_str(), "bs[1] =\n%s\n", qlat::show(bs[1]).c_str());
  qlat::DisplayInfo(cname, fname.c_str(), "bs[2] =\n%s\n", qlat::show(bs[2]).c_str());
  // ADJUST ME
  // qlat::Coordinate xg1(0, 0, 0, total_site[3]/2 + 0);
  // qlat::Coordinate xg2(0, 0, 0, total_site[3]/2 + 0);
  // qlat::Coordinate xg3(0, 0, 0, total_site[3]/2 + 0);
  //
  qlat::Coordinate xg1(1, 1, 1, total_site[3] / 2 + 2);
  qlat::Coordinate xg2(2, 3, 1, total_site[3] / 2 - 1);
  qlat::Coordinate xg3(2, 2, 2, total_site[3] / 2 + 0);
  //
  // qlat::Coordinate xg1(2, 1, 5, total_site[3]/2 - 4);
  // qlat::Coordinate xg2(3, 4, 0, total_site[3]/2 + 0);
  // qlat::Coordinate xg3(1, 5, 2, total_site[3]/2 + 4);
  //
  for (int mu1 = 0; mu1 < 4; ++mu1) {
    for (int mu2 = 0; mu2 < 4; ++mu2) {
      for (int mu3 = 0; mu3 < 4; ++mu3) {
        lblShowMuonPartPointSrc(geo, tsnk, tsrc, bs, xg1, mu1, xg2, mu2, xg3,
                                mu3, mass, momtwist);
        qlat::Timer::display();
        return;
      }
    }
  }
}

void displayGammas()
{
  TIMER("displayGammas");
  qlat::DisplayInfo(cname, fname.c_str(), "gamma5 =\n%s\n",
              qlat::show(qlat::SpinMatrixConstants::get_gamma5()).c_str());
  qlat::DisplayInfo(
      cname, fname.c_str(), "gamma0 * gamma1 * gamma2 * gamma3 =\n%s\n",
      qlat::show((qlat::SpinMatrix)(qlat::SpinMatrixConstants::get_gamma(0) *
                                    qlat::SpinMatrixConstants::get_gamma(1) *
                                    qlat::SpinMatrixConstants::get_gamma(2) *
                                    qlat::SpinMatrixConstants::get_gamma(3)))
          .c_str());
  qlat::DisplayInfo(cname, fname.c_str(), "gamma3 =\n%s\n",
              qlat::show(qlat::SpinMatrixConstants::get_gamma(3)).c_str());
  qlat::DisplayInfo(
      cname, fname.c_str(), "ii * gamma0 =\n%s\n",
      qlat::show((qlat::SpinMatrix)((qlat::ComplexT)qlat::ii *
                                    qlat::SpinMatrixConstants::get_gamma(0)))
          .c_str());
  qlat::DisplayInfo(
      cname, fname.c_str(), "ii * gamma1 =\n%s\n",
      qlat::show((qlat::SpinMatrix)((qlat::ComplexT)qlat::ii *
                                    qlat::SpinMatrixConstants::get_gamma(1)))
          .c_str());
  qlat::DisplayInfo(
      cname, fname.c_str(), "ii * gamma2 =\n%s\n",
      qlat::show((qlat::SpinMatrix)((qlat::ComplexT)qlat::ii *
                                    qlat::SpinMatrixConstants::get_gamma(2)))
          .c_str());
}

void displaySpinPropagator4d()
{
  TIMER("displaySpinPropagator4d");
  // qlat::Coordinate total_site(16, 16, 16, 32);
  qlat::Coordinate total_site(4, 4, 4, 8);
  qlat::Geometry geo;
  geo.init(total_site, 1);
  qlat::DisplayInfo(cname, fname.c_str(), "geo =\n%s\n", qlat::show(geo).c_str());
  qlat::array<double, qlat::DIMN> momtwist;
  momtwist[0] = 0.0;
  momtwist[1] = 0.0;
  momtwist[2] = 0.0;
  momtwist[3] = 0.0;
  const double mass = 0.1;
  qlat::SpinPropagator4d prop;
  prop.init(geo);
  set_zero(prop);
  qlat::Coordinate xgsrc(0, 0, 0, 0);
  qlat::Coordinate xlsrc = geo.coordinate_l_from_g(xgsrc);
  if (geo.is_local(xlsrc)) {
    qlat::set_unit(prop.get_elem(xlsrc));
  }
  qlat::prop_spin_propagator4d(prop, mass, momtwist);
  qlat::Coordinate xgsnk(0, 0, 0, 0);
  qlat::Coordinate xlsnk = geo.coordinate_l_from_g(xgsnk);
  qlat::DisplayInfo(cname, fname.c_str(), "xgsnk = %s .\n", qlat::show(xgsnk).c_str());
  if (geo.is_local(xlsnk)) {
    qlat::Display(cname, fname.c_str(), "prop[xgsnk] =\n%s\n",
                  qlat::show(prop.get_elem(xlsnk)).c_str());
  }
}

int main(int argc, char* argv[])
{
  qlat::begin(&argc, &argv);
  displayGammas();
  displaySpinPropagator4d();
  qlat::Timer::display();
  lblMuonPart();
  qlat::displayln_info("CHECK: finished successfully.");
  qlat::Timer::display();
  qlat::end();
  return 0;
}
