#pragma once

#include <qlat/matrix.h>
#include <qlat/field.h>
#include <qlat/field-fft.h>

#include <Eigen/Eigen>

#include <cmath>

QLAT_START_NAMESPACE

struct QedGaugeField : FieldM<Complex,DIM>
{
  virtual const std::string& cname()
  {
    static const std::string s = "QedGaugeField";
    return s;
  }
};

struct ComplexScalerField : FieldM<Complex,1>
{
  virtual const std::string& cname()
  {
    static const std::string s = "ComplexScalerField";
    return s;
  }
};

struct SpinPropagator4d : FieldM<SpinMatrix,1>
{
  virtual const std::string& cname()
  {
    static const std::string s = "SpinPropagator4d";
    return s;
  }
};

inline void prop_mom_photon_invert(QedGaugeField& egf, const std::array<double,DIM>& momtwist)
  // Feynman Gauge
  // All spatial zero mode removed.
  // egf in momentum space.
{
  TIMER("prop_mom_photon_invert");
  const Geometry& geo = egf.geo;
  for (long index = 0; index < geo.local_volume(); ++index) {
    Coordinate kl = geo.coordinate_from_index(index);
    Coordinate kg = geo.coordinate_g_from_l(kl);
    std::array<double,DIM> kk;
    // FIXME unused 'ks'
	// std::array<double,DIM> ks;
    double s2 = 0.0;
    for (int i = 0; i < DIM; i++) {
      const Coordinate total_site = geo.total_site();
      kg[i] = smod(kg[i], total_site[i]);
      kk[i] = 2.0 * PI * (kg[i] + momtwist[i]) / (double)total_site[i];
      s2 += 4.0 * sqr(std::sin(kk[i] / 2.0));
    }
    if (0.0 == kk[0] && 0.0 == kk[1] && 0.0 == kk[2]) {
      for (int mu = 0; mu < geo.multiplicity; ++mu) {
        egf.get_elem(kl, mu) *= 0.0;
      }
    } else {
      double s2inv = 1.0 / s2;
      for (int mu = 0; mu < geo.multiplicity; ++mu) {
        egf.get_elem(kl, mu) *= s2inv;
      }
    }
  }
}

inline void prop_photon_invert(QedGaugeField& egf, const std::array<double,DIM>& momtwist)
  // Feynman Gauge
  // All spatial zero mode removed.
  // egf in coordinate space.
{
  TIMER_VERBOSE("prop_photon_invert");
  const Geometry& geo = egf.geo;
  fft_complex_field(egf, true);
  prop_mom_photon_invert(egf, momtwist);
  fft_complex_field(egf, false);
  egf *= 1.0 / geo.total_volume();
}

inline void prop_mom_complex_scaler_invert(ComplexScalerField& csf, const double mass, const std::array<double,DIM>& momtwist)
{
  TIMER("prop_mom_complex_scaler_invert");
  const Geometry& geo = csf.geo;
  for (long index = 0; index < geo.local_volume(); ++index) {
    Coordinate kl = geo.coordinate_from_index(index);
    Coordinate kg = geo.coordinate_g_from_l(kl);
    std::array<double,DIM> kk;
    // FIXME my compiler says unused variable 'ks'
	// std::array<double,DIM> ks;
    double s2 = 0.0;
    for (int i = 0; i < DIM; i++) {
      Coordinate total_site = geo.total_site();
      kg[i] = smod(kg[i], total_site[i]);
      kk[i] = 2.0 * PI * (kg[i] + momtwist[i]) / (double)total_site[i];
      s2 += 4.0 * sqr(std::sin(kk[i] / 2.0));
    }
    double s2inv = 1.0 / s2;
    for (int mu = 0; mu < geo.multiplicity; ++mu) {
      csf.get_elem(kl, mu) *= s2inv;
    }
  }
}

inline void prop_complex_scaler_invert(ComplexScalerField& csf, const double mass, const std::array<double,DIM>& momtwist)
{
  TIMER_VERBOSE("prop_complex_scaler_invert");
  const Geometry& geo = csf.geo;
  fft_complex_field(csf, true);
  prop_mom_complex_scaler_invert(csf, mass, momtwist);
  fft_complex_field(csf, false);
  csf *= 1.0 / geo.total_volume();
}

inline void prop_mom_spin_propagator4d(SpinPropagator4d& sp4d, const double mass, const std::array<double,DIM>& momtwist)
  // DWF infinite L_s
  // M_5 = 1.0
{
  TIMER("prop_mom_spin_propagator4d");
  const Geometry& geo = sp4d.geo;
  const double m5 = 1.0;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    Coordinate kl = geo.coordinate_from_index(index);
    Coordinate kg = geo.coordinate_g_from_l(kl);
    std::array<double,DIM> kk, ks;
    double p2 = 0.0;
    double wp = 1.0 - m5;
    SpinMatrix pg; set_zero(pg);
    for (int i = 0; i < DIM; ++i) {
      Coordinate total_site = geo.total_site();
      kg[i] = smod(kg[i], total_site[i]);
      kk[i] = 2.0 * PI * (kg[i] + momtwist[i]) / (double)total_site[i];
      ks[i] = sin(kk[i]);
      pg += SpinMatrixConstants::get_gamma(i) * (Complex)ks[i];
      p2 += sqr(ks[i]);
      wp += 2.0 * sqr(sin(kk[i]/2.0));
    }
    const double calpha = (1.0 + sqr(wp) + p2) / 2.0 / wp;
    const double alpha = std::acosh(calpha);
    const double lwa = 1.0 - wp * exp(-alpha);
    SpinMatrix m; set_unit(m, mass * lwa);
    SpinMatrix ipgm = pg;
    ipgm *= -ii;
    ipgm += m;
    ipgm *= lwa / (p2 + sqr(mass * lwa));
    if (1.0e-10 > p2 && 1.0e-10 > lwa) {
      // if (0.0 != norm(ipgm)) {
      //   Display(cname, fname, "kg = %s\n", show(kg).c_str());
      //   Display(cname, fname, "p2         = %13.5E\n", p2);
      //   Display(cname, fname, "wp         = %13.5E\n", wp);
      //   Display(cname, fname, "alpha      = %13.5E\n", alpha);
      //   Display(cname, fname, "lwa        = %13.5E\n", lwa);
      //   Display(cname, fname, "norm(ipgm) = %13.5E\n", norm(ipgm));
      // }
      set_zero(sp4d.get_elem(kl));
      continue;
    }
    ipgm *= sp4d.get_elem(kl);
    sp4d.get_elem(kl) = ipgm;
  }
}

inline void prop_spin_propagator4d(SpinPropagator4d& sp4d, const double mass, const std::array<double,DIM>& momtwist)
{
  TIMER_VERBOSE("prop_spin_propagator4d");
  const Geometry& geo = sp4d.geo;
  fft_complex_field(sp4d, true);
  prop_mom_spin_propagator4d(sp4d, mass, momtwist);
  fft_complex_field(sp4d, false);
  sp4d *= 1.0 / geo.total_volume();
}

inline void set_point_source_plusm(QedGaugeField& f, const Complex& coef, const Coordinate& xg, const int mu)
{
  TIMER("set_point_source_plusm");
  const Geometry& geo = f.geo;
  Coordinate xl = geo.coordinate_l_from_g(xg);
  if (geo.is_local(xl)) {
    f.get_elem(xl,mu) += coef;
  }
}

inline void set_box_source_plusm(SpinPropagator4d& f, const Complex& coef, const Coordinate& xg1, const Coordinate& xg2)
  // 0 <= xg1 <= xg2 <= total_site
  // FIXME: Do not handle the cross boundary case very well.
{
  TIMER("set_box_source_plusm");
  const Geometry& geo = f.geo;
  SpinMatrix sm; set_unit(sm, coef);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    Coordinate xl = geo.coordinate_from_index(index);
    Coordinate xg = geo.coordinate_g_from_l(xl);
    if (xg1 <= xg && xg < xg2) {
      f.get_elem(xl) += sm;
    }
  }
}

inline void set_wall_source_plusm(SpinPropagator4d& f, const Complex& coef, const int t)
{
  TIMER("set_wall_source_plusm");
  const Geometry& geo = f.geo;
  Coordinate total_site = geo.total_site();
  qassert(0 <= t && t < total_site[3]);
  const Coordinate xg1(0, 0, 0, t);
  const Coordinate xg2(total_site[0], total_site[1], total_site[2], t+1);
  set_box_source_plusm(f, coef, xg1, xg2);
}

inline void sequential_photon_spin_propagator_plusm(
    SpinPropagator4d& src, const Complex coef,
    const QedGaugeField& egf, const SpinPropagator4d& sol)
{
  TIMER("sequential_photon_spin_propagator_plusm");
  const Geometry& geo = sol.geo;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    Coordinate xl = geo.coordinate_from_index(index);
    for (int mu = 0; mu < DIM; mu++) {
      // a = A_\mu(x)
      Complex a = egf.get_elem(xl, mu);
      // tmp = \gamma_\mu \psi(x)
      SpinMatrix tmp = SpinMatrixConstants::get_gamma(mu) * sol.get_elem(xl);
      // tmp = coef * \gamma_\mu A_\mu(x) \psi(x)
      tmp *= a * coef;
      src.get_elem(xl) += tmp;
    }
  }
}

inline SpinMatrix contract_spin_propagator4d(const SpinPropagator4d& snk, const SpinPropagator4d& src)
{
  TIMER("contractSpinPropagator");
  const Geometry& geo = src.geo;
  SpinMatrix sum; set_zero(sum);
#pragma omp parallel
  {
    SpinMatrix psum; set_zero(psum);
#pragma omp for nowait
    for (long index = 0; index < geo.local_volume(); ++index) {
      Coordinate xl = geo.coordinate_from_index(index);
      psum += matrix_adjoint(snk.get_elem(xl)) * src.get_elem(xl);
    }
    for (int i = 0; i < omp_get_num_threads(); ++i) {
#pragma omp barrier
      if (omp_get_thread_num() == i) {
          sum += psum;
      }
    }
  }
  glb_sum_double(sum);
  return sum;
}

QLAT_END_NAMESPACE
