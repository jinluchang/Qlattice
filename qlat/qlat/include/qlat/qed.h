#pragma once

#include <qlat-utils/matrix.h>
#include <qlat/core.h>
#include <qlat/field-fft.h>
#include <qlat/field.h>

#include <cmath>

namespace qlat
{  //

template <class T = Real>
struct QedGaugeFieldT : FieldM<ComplexT<T>, DIMN> {
};

template <class T = Real>
struct ComplexScalerFieldT : FieldM<ComplexT<T>, 1> {
};

template <class T = Real>
struct SpinPropagator4dT : FieldM<SpinMatrixT<T>, 1> {
};

#ifndef QLAT_NO_DEFAULT_TYPE

using QedGaugeField = QedGaugeFieldT<>;

using ComplexScalerField = ComplexScalerFieldT<>;

using SpinPropagator4d = SpinPropagator4dT<>;

using SpinProp = SpinPropagator4d;

#endif

template <class T>
inline void take_real_part_and_multiply_sqrt2(Field<T>& f)
{
  TIMER("take_real_part_and_multiply_sqrt2");
  const Geometry& geo = f.geo();
  // const Int multiplicity = f.multiplicity;
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<T> fv = f.get_elems(xl);
    for (Int m = 0; m < f.multiplicity; ++m) {
      fv[m] = sqrt(2.0) * fv[m].real();
    }
  }
}

template <class T>
inline void set_mom_stochastic_qed_field_feynman(Field<T>& f,
                                                 const Geometry& geo,
                                                 const Int multiplicity,
                                                 const RngState& rs)
// use QED_L scheme: all spatial zero mode removed.
{
  TIMER("set_mom_stochastic_qed_field_feynman");
  f.init(geo, multiplicity);
  const RealD total_volume = geo.total_volume();
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate kl = geo.coordinate_from_index(index);
    const Coordinate kg = geo.coordinate_g_from_l(kl);
    const Long gindex = geo.g_index_from_g_coordinate(kg);
    RngState rst = rs.newtype(gindex);
    RealD s2 = 0.0;
    array<RealD, DIMN> kk;
    for (Int i = 0; i < DIMN; i++) {
      const Coordinate total_site = geo.total_site();
      kk[i] = 2.0 * PI * smod(kg[i], total_site[i]) / (RealD)total_site[i];
      s2 += 4.0 * sqr(std::sin(kk[i] / 2.0));
    }
    Vector<ComplexD> fv = f.get_elems(kl);
    if (0.0 == kk[0] && 0.0 == kk[1] && 0.0 == kk[2]) {
      for (Int m = 0; m < multiplicity; ++m) {
        fv[m] = 0.0;
      }
    } else {
      const RealD sigma = std::sqrt(1.0 / (2.0 * total_volume * s2));
      for (Int m = 0; m < multiplicity; ++m) {
        const RealD re = g_rand_gen(rst, 0.0, sigma);
        const RealD im = g_rand_gen(rst, 0.0, sigma);
        fv[m] = T(re, im);
      }
    }
  }
}

template <class T>
inline void set_stochastic_qed_field_feynman(Field<T>& f, const Geometry& geo,
                                             const Int multiplicity,
                                             const RngState& rs)
{
  TIMER("set_stochastic_qed_field_feynman");
  set_mom_stochastic_qed_field_feynman(f, geo, multiplicity, rs);
  fft_complex_field(f, false);
  take_real_part_and_multiply_sqrt2(f);
}

template <class T>
inline void set_mom_stochastic_qed_field_mass(Field<T>& f, const Geometry& geo,
                                              const Int multiplicity,
                                              const RealD mass,
                                              const RngState& rs)
// use mass scheme, all zero modes are kept
{
  TIMER("set_mom_stochastic_qed_field_mass");
  f.init(geo, multiplicity);
  const RealD total_volume = geo.total_volume();
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate kl = geo.coordinate_from_index(index);
    const Coordinate kg = geo.coordinate_g_from_l(kl);
    const Long gindex = geo.g_index_from_g_coordinate(kg);
    RngState rst = rs.newtype(gindex);
    RealD s2 = sqr(mass);
    array<RealD, DIMN> kk;
    for (Int i = 0; i < DIMN; i++) {
      const Coordinate total_site = geo.total_site();
      kk[i] = 2.0 * PI * smod(kg[i], total_site[i]) / (RealD)total_site[i];
      s2 += 4.0 * sqr(std::sin(kk[i] / 2.0));
    }
    Vector<ComplexD> fv = f.get_elems(kl);
    const RealD sigma = std::sqrt(1.0 / (2.0 * total_volume * s2));
    for (Int m = 0; m < multiplicity; ++m) {
      const RealD re = g_rand_gen(rst, 0.0, sigma);
      const RealD im = g_rand_gen(rst, 0.0, sigma);
      fv[m] = T(re, im);
    }
  }
}

template <class T>
inline void set_stochastic_qed_field_mass(Field<T>& f, const Geometry& geo,
                                          const RealD mass, const RngState& rs)
{
  TIMER("set_stochastic_qed_field_mass");
  set_mom_stochastic_qed_field_mass(f, geo, mass, rs);
  fft_complex_field(f, false);
  take_real_part_and_multiply_sqrt2(f);
}

inline void prop_mom_photon_invert(QedGaugeField& egf,
                                   const CoordinateD& momtwist)
// Feynman Gauge
// All spatial zero mode removed.
// egf in momentum space.
{
  TIMER("prop_mom_photon_invert");
  const Geometry& geo = egf.geo();
  const Int multiplicity = egf.multiplicity;
  for (Long index = 0; index < geo.local_volume(); ++index) {
    Coordinate kl = geo.coordinate_from_index(index);
    Coordinate kg = geo.coordinate_g_from_l(kl);
    array<RealD, DIMN> kk;
    RealD s2 = 0.0;
    for (Int i = 0; i < DIMN; i++) {
      const Coordinate total_site = geo.total_site();
      kg[i] = smod(kg[i], total_site[i]);
      kk[i] = 2.0 * PI * (kg[i] + momtwist[i]) / (RealD)total_site[i];
      s2 += 4.0 * sqr(std::sin(kk[i] / 2.0));
    }
    if (0.0 == kk[0] && 0.0 == kk[1] && 0.0 == kk[2]) {
      for (Int mu = 0; mu < multiplicity; ++mu) {
        egf.get_elem(kl, mu) *= 0.0;
      }
    } else {
      RealD s2inv = 1.0 / s2;
      for (Int mu = 0; mu < multiplicity; ++mu) {
        egf.get_elem(kl, mu) *= s2inv;
      }
    }
  }
}

inline void prop_photon_invert(QedGaugeField& egf, const CoordinateD& momtwist)
// Feynman Gauge
// All spatial zero mode removed.
// egf in coordinate space.
{
  TIMER_VERBOSE("prop_photon_invert");
  const Geometry geo = egf.geo();
  fft_complex_field(egf, true);
  prop_mom_photon_invert(egf, momtwist);
  fft_complex_field(egf, false);
  egf *= 1.0 / geo.total_volume();
}

inline void prop_mom_complex_scalar_invert(ComplexScalerField& csf,
                                           const RealD mass,
                                           const CoordinateD& momtwist)
{
  TIMER("prop_mom_complex_scalar_invert");
  const Geometry& geo = csf.geo();
  const Int multiplicity = csf.multiplicity;
  for (Long index = 0; index < geo.local_volume(); ++index) {
    Coordinate kl = geo.coordinate_from_index(index);
    Coordinate kg = geo.coordinate_g_from_l(kl);
    array<RealD, DIMN> kk;
    RealD s2 = sqr(mass);
    for (Int i = 0; i < DIMN; i++) {
      Coordinate total_site = geo.total_site();
      kg[i] = smod(kg[i], total_site[i]);
      kk[i] = 2.0 * PI * (kg[i] + momtwist[i]) / (RealD)total_site[i];
      s2 += 4.0 * sqr(std::sin(kk[i] / 2.0));
    }
    RealD s2inv = 1.0 / s2;
    for (Int mu = 0; mu < multiplicity; ++mu) {
      csf.get_elem(kl, mu) *= s2inv;
    }
  }
}

inline void prop_complex_scalar_invert(ComplexScalerField& csf,
                                       const RealD mass,
                                       const CoordinateD& momtwist)
{
  TIMER_VERBOSE("prop_complex_scalar_invert");
  const Geometry& geo = csf.geo();
  fft_complex_field(csf, true);
  prop_mom_complex_scalar_invert(csf, mass, momtwist);
  fft_complex_field(csf, false);
  csf *= 1.0 / geo.total_volume();
}

inline RealD acosh(const RealD x)
{
  return std::log(x + std::sqrt(x + 1.0) * std::sqrt(x - 1.0));
}

inline void prop_free_scalar_invert(Field<ComplexD>& f, const RealD mass,
                                    const CoordinateD& momtwist)
{
  TIMER("prop_free_scalar_invert");
  const Geometry& geo = f.geo();
  const Coordinate total_site = geo.total_site();
  const RealD m_pi_sq = 4.0 * sqr(std::sinh(mass / 2.0));
  qacc_for(index, geo.local_volume(), {
    const Coordinate kl = geo.coordinate_from_index(index);
    Coordinate kg = geo.coordinate_g_from_l(kl);
    CoordinateD kk, ks;
    RealD s2 = 0.0;
    for (Int i = 0; i < DIMN; i++) {
      kg[i] = smod(kg[i], total_site[i]);
      kk[i] = 2.0 * PI * (kg[i] + momtwist[i]) / (RealD)total_site[i];
      ks[i] = 2.0 * std::sin(kk[i] / 2.0);
      s2 += sqr(ks[i]);
    }
    const RealD fac = 1.0 / (m_pi_sq + s2);
    Vector<ComplexD> v = f.get_elems(kl);
    for (Int i = 0; i < v.size(); ++i) {
      v[i] *= fac;
    }
  })
}

template <class T>
void prop_mom_spin_propagator4d(SpinPropagator4dT<T>& sp4d, const RealD mass,
                                const RealD m5, const CoordinateD& momtwist)
// DWF infinite L_s
// M_5 <= 1.0
{
  TIMER("prop_mom_spin_propagator4d");
  const Geometry& geo = sp4d.geo();
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    Coordinate kl = geo.coordinate_from_index(index);
    Coordinate kg = geo.coordinate_g_from_l(kl);
    array<RealD, DIMN> kk, ks;
    RealD p2 = 0.0;
    RealD wp = 1.0 - m5;
    SpinMatrixT<T> pg;
    set_zero(pg);
    for (Int i = 0; i < DIMN; ++i) {
      Coordinate total_site = geo.total_site();
      kg[i] = smod(kg[i], total_site[i]);
      kk[i] = 2.0 * PI * (kg[i] + momtwist[i]) / (RealD)total_site[i];
      ks[i] = sin(kk[i]);
      pg += SpinMatrixConstantsT<T>::get_cps_gamma(i) * (ComplexT<T>)ks[i];
      p2 += sqr(ks[i]);
      wp += 2.0 * sqr(sin(kk[i] / 2.0));
    }
    const RealD calpha = (1.0 + sqr(wp) + p2) / 2.0 / wp;
    const RealD alpha = acosh(calpha);
    const RealD lwa = 1.0 - wp * exp(-alpha);
    SpinMatrixT<T> m;
    set_unit(m, mass * lwa);
    SpinMatrixT<T> ipgm = pg;
    ipgm *= (ComplexT<T>)(-ii);
    ipgm += m;
    ipgm *= lwa / (p2 + sqr(mass * lwa));
    if (1.0e-10 > p2 && 1.0e-10 > lwa) {
      // if (0.0 != qnorm(ipgm)) {
      //   Display(cname, fname, "kg = %s\n", show(kg).c_str());
      //   Display(cname, fname, "p2         = %13.5E\n", p2);
      //   Display(cname, fname, "wp         = %13.5E\n", wp);
      //   Display(cname, fname, "alpha      = %13.5E\n", alpha);
      //   Display(cname, fname, "lwa        = %13.5E\n", lwa);
      //   Display(cname, fname, "qnorm(ipgm) = %13.5E\n", qnorm(ipgm));
      // }
      set_zero(sp4d.get_elem(kl));
      continue;
    }
    ipgm *= sp4d.get_elem(kl);
    sp4d.get_elem(kl) = ipgm;
  }
}

template <class T>
void prop_spin_propagator4d(SpinPropagator4dT<T>& sp4d, const RealD mass,
                            const RealD m5, const CoordinateD& momtwist)
{
  TIMER_VERBOSE("prop_spin_propagator4d");
  const Geometry geo = sp4d.geo();
  fft_complex_field(sp4d, true);
  prop_mom_spin_propagator4d(sp4d, mass, m5, momtwist);
  fft_complex_field(sp4d, false);
  sp4d *= 1.0 / geo.total_volume();
}

inline void set_point_source_plusm(QedGaugeField& f, const ComplexD& coef,
                                   const Coordinate& xg, const Int mu)
{
  TIMER("set_point_source_plusm");
  const Geometry& geo = f.geo();
  Coordinate xl = geo.coordinate_l_from_g(xg);
  if (geo.is_local(xl)) {
    f.get_elem(xl, mu) += coef;
  }
}

template <class T>
void set_box_source_plusm(SpinPropagator4dT<T>& f, const ComplexD& coef,
                          const Coordinate& xg1, const Coordinate& xg2)
// 0 <= xg1 <= xg2 <= total_site
// FIXME: Do not handle the cross boundary case very well.
{
  TIMER("set_box_source_plusm");
  const Geometry& geo = f.geo();
  SpinMatrixT<T> sm;
  set_unit(sm, coef);
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    Coordinate xl = geo.coordinate_from_index(index);
    Coordinate xg = geo.coordinate_g_from_l(xl);
    if (xg1 <= xg && xg < xg2) {
      f.get_elem(xl) += sm;
    }
  }
}

template <class T>
void set_wall_source_plusm(SpinPropagator4dT<T>& f, const ComplexD& coef,
                           const Int t)
{
  TIMER("set_wall_source_plusm");
  const Geometry& geo = f.geo();
  Coordinate total_site = geo.total_site();
  qassert(0 <= t && t < total_site[3]);
  const Coordinate xg1(0, 0, 0, t);
  const Coordinate xg2(total_site[0], total_site[1], total_site[2], t + 1);
  set_box_source_plusm(f, coef, xg1, xg2);
}

template <class T>
void sequential_photon_spin_propagator_plusm(SpinPropagator4dT<T>& src,
                                             const ComplexD coef,
                                             const QedGaugeField& egf,
                                             const SpinPropagator4dT<T>& sol)
{
  TIMER("sequential_photon_spin_propagator_plusm");
  const Geometry& geo = sol.geo();
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    Coordinate xl = geo.coordinate_from_index(index);
    for (Int mu = 0; mu < DIMN; mu++) {
      // a = A_\mu(x)
      ComplexD a = egf.get_elem(xl, mu);
      // tmp = \gamma_\mu \psi(x)
      SpinMatrixT<T> tmp =
          SpinMatrixConstantsT<T>::get_cps_gamma(mu) * sol.get_elem(xl);
      // tmp = coef * \gamma_\mu A_\mu(x) \psi(x)
      tmp *= (ComplexT<T>)(a * coef);
      src.get_elem(xl) += tmp;
    }
  }
}

template <class T>
SpinMatrixT<T> contract_spin_propagator4d(const SpinPropagator4dT<T>& snk,
                                          const SpinPropagator4dT<T>& src)
{
  TIMER("contractSpinPropagator");
  const Geometry& geo = src.geo();
  SpinMatrixT<T> sum;
  set_zero(sum);
#pragma omp parallel
  {
    SpinMatrixT<T> psum;
    set_zero(psum);
#pragma omp for nowait
    for (Long index = 0; index < geo.local_volume(); ++index) {
      Coordinate xl = geo.coordinate_from_index(index);
      psum += matrix_adjoint(snk.get_elem(xl)) * src.get_elem(xl);
    }
    for (Int i = 0; i < omp_get_num_threads(); ++i) {
#pragma omp barrier
      if (omp_get_thread_num() == i) {
        sum += psum;
      }
    }
  }
  glb_sum(sum);
  return sum;
}

void set_left_expanded_gauge_field(Field<ComplexD>& gf1,
                                   const Field<ComplexD>& gf);

void free_invert(SpinProp& sp_sol, SpinProp& sp_src, const RealD mass,
                 const RealD m5 = 1.0,
                 const CoordinateD& momtwist = CoordinateD());

void invert_qed(SpinProp& sp_sol, const SpinProp& sp_src,
                const Field<ComplexD>& gf1, const RealD mass, const RealD m5,
                const Int ls, const bool is_dagger = false,
                const RealD stop_rsd = 1e-8, const Long max_num_iter = 50000);

void fermion_field_4d_from_5d_qed(Field<ComplexD>& ff4d,
                                  const Field<ComplexD>& ff5d, const Int ls,
                                  const Int upper, const Int lower);

void fermion_field_5d_from_4d_qed(Field<ComplexD>& ff5d,
                                  const Field<ComplexD>& ff4d, const Int ls,
                                  const Int upper, const Int lower);

Long invert_dwf_qed(Field<ComplexD>& f_out4d, const Field<ComplexD>& f_in4d,
                    const Field<ComplexD>& gf1, const RealD mass,
                    const RealD m5, const Int ls, const bool is_dagger = false,
                    const RealD stop_rsd = 1e-8,
                    const Long max_num_iter = 50000);

ComplexD dot_product(const Field<ComplexD>& ff1, const Field<ComplexD>& ff2);

Long cg_with_m_dwf_qed(Field<ComplexD>& f_out5d, const Field<ComplexD>& f_in5d,
                       const Field<ComplexD>& gf1, const RealD mass,
                       const RealD m5, const Int ls, const bool is_dagger = false,
                       const RealD stop_rsd = 1e-8,
                       const Long max_num_iter = 50000);

void multiply_m_dwf_qed(Field<ComplexD>& f_out5d, const Field<ComplexD>& f_in5d,
                        const Field<ComplexD>& gf1, const RealD mass,
                        const RealD m5, const Int ls, const bool is_dagger = false);

}  // namespace qlat
