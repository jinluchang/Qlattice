#pragma once

#include <qlat/qcd.h>
#include <qlat/qcd-utils.h>
#include <qlat/fermion-action.h>

QLAT_START_NAMESPACE

struct InverterParams
{
  double stop_rsd;
  long max_num_iter;
  long max_mixed_precision_cycle;
  //
  void init()
  {
    stop_rsd = 1.0e-8;
    max_num_iter = 50000;
    max_mixed_precision_cycle = 100;
  }
  //
  InverterParams()
  {
    init();
  }
};

struct LowModes
{
  // TODO
};

struct InverterDomainWall
{
  Geometry geo;
  FermionAction fa;
  GaugeField gf;
  InverterParams ip;
  //
  InverterDomainWall()
  {
    init();
  }
  //
  void init()
  {
    ip.init();
  }
  //
  void setup()
  {
    TIMER_VERBOSE("Inv::setup");
  }
  void setup(const GaugeField& gf_, const FermionAction& fa_)
  {
    geo = geo_reform(gf_.geo);
    gf.init();
    set_left_expanded_gauge_field(gf, gf_);
    fa = fa_;
    setup();
  }
  void setup(const GaugeField& gf_, const FermionAction& fa_, const LowModes& lm_)
  {
    setup(gf_, fa_);
    // TODO
  }
  //
  double& stop_rsd()
  {
    return ip.stop_rsd;
  }
  //
  long& max_num_iter()
  {
    return ip.max_num_iter;
  }
  //
  long& max_mixed_precision_cycle()
  {
    return ip.max_mixed_precision_cycle;
  }
};

inline void load_or_compute_low_modes(LowModes& lm, const std::string& path,
    const GaugeField& gf, const FermionAction& fa, const LancArg& la)
{
  TIMER_VERBOSE("load_or_compute_low_modes");
  // TODO
}

inline void setup_inverter(InverterDomainWall& inv)
{
  inv.setup();
}

inline void setup_inverter(InverterDomainWall& inv, const GaugeField& gf, const FermionAction& fa)
{
  inv.setup(gf, fa);
}

inline void setup_inverter(InverterDomainWall& inv, const GaugeField& gf, const FermionAction& fa, const LowModes& lm)
{
  inv.setup(gf, fa, lm);
}

inline void multiply_m_dwf_no_comm(FermionField5d& out, const FermionField5d& in, const InverterDomainWall& inv)
{
  TIMER("multiply_m_dwf_no_comm(5d,5d,Inv)");
  const Geometry geo = geo_resize(in.geo);
  qassert(is_matching_geo(inv.geo, geo));
  out.init(geo);
  set_zero(out);
  const FermionAction& fa = inv.fa;
  const GaugeField& gf = inv.gf;
  qassert(in.geo.multiplicity == fa.ls);
  qassert(out.geo.multiplicity == fa.ls);
  qassert(fa.mobius_scale == 1.0);
  qassert(fa.bs.size() == fa.ls);
  qassert(fa.cs.size() == fa.ls);
  const std::array<SpinMatrix,4>& gammas = SpinMatrixConstants::get_cps_gammas();
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  const SpinMatrix& unit = SpinMatrixConstants::get_unit();
  const SpinMatrix p_p = 0.5 * (unit + gamma5);
  const SpinMatrix p_m = 0.5 * (unit - gamma5);
  std::array<SpinMatrix,4> p_mu_p;
  std::array<SpinMatrix,4> p_mu_m;
  for (int mu = 0; mu < 4; ++mu) {
    p_mu_p[mu] = 0.5 * (unit + gammas[mu]);
    p_mu_m[mu] = 0.5 * (unit - gammas[mu]);
  }
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<WilsonVector> v = out.get_elems(xl);
    {
      const Vector<WilsonVector> iv = in.get_elems_const(xl);
      for (int m = 0; m < fa.ls; ++m) {
        v[m] = (5.0 - fa.m5) * iv[m];
        v[m] -= p_m * (m < fa.ls-1 ? iv[m+1] : (WilsonVector)(-fa.mass * iv[0]));
        v[m] -= p_p * (m > 0 ? iv[m-1] : (WilsonVector)(-fa.mass * iv[fa.ls-1]));
      }
    }
    for (int mu = 0; mu < 4; ++mu) {
      const Coordinate xl_p = coordinate_shifts(xl, mu);
      const Coordinate xl_m = coordinate_shifts(xl, -mu-1);
      const ColorMatrix u_p = gf.get_elem(xl, mu);
      const ColorMatrix u_m = matrix_adjoint(gf.get_elem(xl_m, mu));
      const Vector<WilsonVector> iv_p = in.get_elems_const(xl_p);
      const Vector<WilsonVector> iv_m = in.get_elems_const(xl_m);
      for (int m = 0; m < fa.ls; ++m) {
        v[m] -= u_p * (p_mu_m[mu] * iv_p[m]);
        v[m] -= u_m * (p_mu_p[mu] * iv_m[m]);
      }
    }
  }
}

inline void multiply_m_dwf(FermionField5d& out, const FermionField5d& in, const InverterDomainWall& inv)
  // out can be the same object as in
{
  TIMER("multiply_m_dwf(5d,5d,Inv)");
  const Geometry geo1 = geo_resize(in.geo, 1);
  FermionField5d in1;
  in1.init(geo1);
  in1 = in;
  refresh_expanded_1(in1);
  multiply_m_dwf_no_comm(out, in1, inv);
}

inline void multiply_wilson_d_no_comm(FermionField5d& out, const FermionField5d& in, const GaugeField& gf, const double mass)
  // gf.init();
  // set_left_expanded_gauge_field(gf, gf_);
  // in.geo = geo_reform(geo, 1, ls);
  // refresh_expanded_1(in);
{
  TIMER("multiply_wilson_d_no_comm(5d,5d,Inv)");
  const Geometry geo = geo_resize(in.geo);
  qassert(is_matching_geo(gf.geo, geo));
  out.init(geo);
  set_zero(out);
  const int ls = in.geo.multiplicity;
  qassert(out.geo.multiplicity == ls);
  const std::array<SpinMatrix,4>& gammas = SpinMatrixConstants::get_cps_gammas();
  const SpinMatrix& unit = SpinMatrixConstants::get_unit();
  std::array<SpinMatrix,4> p_mu_p;
  std::array<SpinMatrix,4> p_mu_m;
  for (int mu = 0; mu < 4; ++mu) {
    p_mu_p[mu] = 0.5 * (unit + gammas[mu]);
    p_mu_m[mu] = 0.5 * (unit - gammas[mu]);
  }
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<WilsonVector> v = out.get_elems(xl);
    {
      const Vector<WilsonVector> iv = in.get_elems_const(xl);
      for (int m = 0; m < ls; ++m) {
        v[m] = (4.0 + mass) * iv[m];
      }
    }
    for (int mu = 0; mu < 4; ++mu) {
      const Coordinate xl_p = coordinate_shifts(xl, mu);
      const Coordinate xl_m = coordinate_shifts(xl, -mu-1);
      const ColorMatrix u_p = gf.get_elem(xl, mu);
      const ColorMatrix u_m = matrix_adjoint(gf.get_elem(xl_m, mu));
      const Vector<WilsonVector> iv_p = in.get_elems_const(xl_p);
      const Vector<WilsonVector> iv_m = in.get_elems_const(xl_m);
      for (int m = 0; m < ls; ++m) {
        v[m] -= u_p * (p_mu_m[mu] * iv_p[m]);
        v[m] -= u_m * (p_mu_p[mu] * iv_m[m]);
      }
    }
  }
}

inline void multiply_d_minus(FermionField5d& out, const FermionField5d& in, const InverterDomainWall& inv)
{
  TIMER("multiply_d_minus(5d,5d,Inv)");
  const Geometry geo = geo_resize(in.geo);
  out.init(geo);
  set_zero(out);
  const FermionAction& fa = inv.fa;
  qassert(is_matching_geo(inv.geo, in.geo));
  qassert(is_matching_geo(inv.geo, out.geo));
  qassert(in.geo.multiplicity == fa.ls);
  qassert(out.geo.multiplicity == fa.ls);
  qassert(fa.bs.size() == fa.ls);
  qassert(fa.cs.size() == fa.ls);
  const GaugeField& gf = inv.gf;
  const Geometry geo1 = geo_resize(in.geo, 1);
  FermionField5d in1;
  in1.init(geo1);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<WilsonVector> v = out.get_elems(xl);
    Vector<WilsonVector> v1 = in1.get_elems(xl);
    const Vector<WilsonVector> iv = in.get_elems_const(xl);
    for (int m = 0; m < fa.ls; ++m) {
      const Complex& c = fa.cs[m];
      v1[m] = c * iv[m];
      v[m] -= iv[m];
    }
  }
  refresh_expanded_1(in1);
  FermionField5d out1;
  multiply_wilson_d_no_comm(out1, in1, gf, -fa.m5);
  out += out1;
}

inline void multiply_m(FermionField5d& out, const FermionField5d& in, const InverterDomainWall& inv)
  // out can be the same object as in
{
  TIMER("multiply_m(5d,5d,Inv)");
  const Geometry geo = geo_resize(in.geo);
  out.init(geo);
  set_zero(out);
  const FermionAction& fa = inv.fa;
  qassert(is_matching_geo(inv.geo, in.geo));
  qassert(is_matching_geo(inv.geo, out.geo));
  qassert(geo.multiplicity == fa.ls);
  qassert(in.geo.multiplicity == fa.ls);
  qassert(out.geo.multiplicity == fa.ls);
  qassert(fa.bs.size() == fa.ls);
  qassert(fa.cs.size() == fa.ls);
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  const SpinMatrix& unit = SpinMatrixConstants::get_unit();
  const SpinMatrix p_p = 0.5 * (unit + gamma5);
  const SpinMatrix p_m = 0.5 * (unit - gamma5);
  const GaugeField& gf = inv.gf;
  const Geometry geo1 = geo_resize(in.geo, 1);
  FermionField5d in1, fftmp;
  in1.init(geo1);
  fftmp.init(geo);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<WilsonVector> iv = in.get_elems_const(xl);
    Vector<WilsonVector> v = fftmp.get_elems(xl);
    Vector<WilsonVector> v1 = in1.get_elems(xl);
    for (int m = 0; m < fa.ls; ++m) {
      const Complex& b = fa.bs[m];
      const Complex& c = fa.cs[m];
      v1[m] = iv[m];
      v[m] = b * iv[m];
      const WilsonVector tmp =
        (p_m * (m < fa.ls-1 ? iv[m+1] : (WilsonVector)(-fa.mass * iv[0]))) +
        (p_p * (m > 0 ? iv[m-1] : (WilsonVector)(-fa.mass * iv[fa.ls-1])));
      v1[m] += c * tmp;
      v[m] -= tmp;
    }
  }
  refresh_expanded_1(in1);
  multiply_wilson_d_no_comm(out, in1, gf, -fa.m5);
  out += fftmp;
}

inline void get_half_fermion(FermionField5d& half, const FermionField5d& ff, const int eo)
{
  TIMER("get_half_fermion");
  Geometry geoh = geo_resize(ff.geo);
  geoh.eo = eo;
  half.init(geoh);
  qassert(half.geo.eo == eo);
  qassert(ff.geo.eo == 0);
  qassert(is_matching_geo(ff.geo, half.geo));
#pragma omp parallel for
  for (long index = 0; index < geoh.local_volume(); ++index) {
    const Coordinate xl = geoh.coordinate_from_index(index);
    assign(half.get_elems(xl), ff.get_elems_const(xl));
  }
}

inline void set_half_fermion(FermionField5d& ff, const FermionField5d& half, const int eo)
{
  TIMER("set_half_fermion");
  const Geometry geoh = half.geo;
  Geometry geo = geo_resize(geoh);
  geo.eo = 0;
  ff.init(geo);
  qassert(half.geo.eo == eo);
  qassert(ff.geo.eo == 0);
  qassert(is_matching_geo(ff.geo, half.geo));
#pragma omp parallel for
  for (long index = 0; index < geoh.local_volume(); ++index) {
    const Coordinate xl = geoh.coordinate_from_index(index);
    assign(ff.get_elems(xl), half.get_elems_const(xl));
  }
}

inline void multiply_m_eo_eo_no_comm(FermionField5d& out, const FermionField5d& in, const InverterDomainWall& inv,
    const int eo_out, const int eo_in)
  // out need to be initialized with correct geo and eo
{
  TIMER_VERBOSE("half_multiply_m_no_comm");
  qassert(is_matching_geo(out.geo, in.geo));
  out.geo.eo == eo_out;
  in.geo.eo == eo_in;
  // TODO
}

inline bool& is_checking_inverse()
  // qlat parameter
{
  static bool b = false;
  return b;
}

QLAT_END_NAMESPACE
