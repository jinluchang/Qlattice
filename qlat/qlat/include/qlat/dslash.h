#pragma once

#include <qlat/fermion-action.h>
#include <qlat/qcd-utils.h>
#include <qlat/qcd.h>
#include <qlat/compressed-eigen-io.h>

namespace qlat
{  //

API inline bool& is_checking_invert()
// qlat parameter
{
  static bool b = true;
  return b;
}

API inline bool& is_cg_verbose()
// qlat parameter
{
  static bool b = false;
  return b;
}

struct LowModesInfo {
  bool initialized;
  std::string path;
  GaugeField gf;
  FermionAction fa;
  LancArg la;
  //
  LowModesInfo() { init(); }
  //
  void init()
  {
    initialized = false;
    path = "";
    gf.init();
    fa.init();
    la.init();
  }
};

struct LowModes {
  bool initialized;
  LowModesInfo lmi;
  vector<double> eigen_values;
  CompressedEigenSystemInfo cesi;
  CompressedEigenSystemBases cesb;
  CompressedEigenSystemCoefs cesc;
  //
  LowModes() { init(); }
  //
  void init()
  {
    initialized = false;
    lmi.init();
    eigen_values.init();
    cesi.init();
    cesb.init();
    cesc.init();
  }
  void init(const Geometry& geo, const Int multiplicity,
            const Coordinate& block_site, const Long neig, const Long nkeep)
  // multiplicity = ls
  {
    initialized = true;
    lmi.init();
    eigen_values.resize(neig);
    CompressedEigenSystemDenseInfo cesdi;
    cesdi.total_site = geo.total_site();
    cesdi.block_site = block_site;
    cesdi.ls = multiplicity;
    cesdi.neig = neig;
    cesdi.nkeep = nkeep;
    cesdi.nkeep_single = nkeep;
    cesdi.FP16_COEF_EXP_SHARE_FLOATS = 10;
    qassert(cesdi.total_site % geo.geon.size_node == Coordinate());
    cesdi.node_site = cesdi.total_site / geo.geon.size_node;
    cesi = populate_eigen_system_info(
        cesdi, std::vector<crc32_t>(product(geo.geon.size_node), 0));
    init_compressed_eigen_system_bases(cesb, cesi, geo.geon.id_node);
    init_compressed_eigen_system_coefs(cesc, cesi, geo.geon.id_node);
  }
};

Long load_low_modes(LowModes& lm, const std::string& path);

Long load_or_compute_low_modes(LowModes& lm, const std::string& path,
                               const GaugeField& gf, const FermionAction& fa,
                               const LancArg& la);

void load_low_modes_delay(LowModes& lm, const std::string& path);

void load_or_compute_low_modes_delay(LowModes& lm, const std::string& path,
                                     const GaugeField& gf,
                                     const FermionAction& fa,
                                     const LancArg& la);

Long force_low_modes(LowModes& lm);

void set_u_rand(LowModes& lm, const RngState& rs);

Long save_low_modes_decompress(LowModes& lm, const std::string& path);

void deflate(HalfVector& hv_out, const HalfVector& hv_in, LowModes& lm);

void deflate(FermionField5d& out, const FermionField5d& in, LowModes& lm);

void benchmark_deflate(const Geometry& geo, const Int ls,
                       const Coordinate& block_site, const Long neig,
                       const Long nkeep, const RngState& rs);

struct InverterParams {
  double stop_rsd;
  Long max_num_iter;
  Long max_mixed_precision_cycle;
  Int solver_type;  // 0 -> CG, 1-> EIGCG, 2->MSPCG
  Int higher_precision;
  Int lower_precision;
  //
  void init()
  {
    stop_rsd = 1.0e-8;
    max_num_iter = 200;
    max_mixed_precision_cycle = 300;
    solver_type = 0;
    higher_precision = 8;
    lower_precision = 8;
  }
  //
  InverterParams() { init(); }
};

struct InverterDomainWall {
  box<Geometry> geo;
  FermionAction fa;
  GaugeField gf;
  InverterParams ip;
  Handle<LowModes> lm;
  //
  InverterDomainWall() { init(); }
  //
  void init()
  {
    geo.init();
    fa.init();
    gf.init();
    ip.init();
    lm.init();
  }
  //
  void setup() {}
  void setup(const GaugeField& gf_, const FermionAction& fa_)
  {
    TIMER_VERBOSE("Inv::setup(gf,fa)");
    geo.set(geo_resize(gf_.geo()));
    gf.init();
    set_left_expanded_gauge_field(gf, gf_);
    fa = fa_;
    lm.init();
  }
  //
  void setup(const GaugeField& gf_, const FermionAction& fa_,
             const InverterParams& ip_)
  {
    TIMER_VERBOSE("Inv::setup(gf,fa)");
    geo.set(geo_resize(gf_.geo()));
    gf.init();
    set_left_expanded_gauge_field(gf, gf_);
    fa = fa_;
    ip = ip_;
    lm.init();
  }
  //
  void setup(const GaugeField& gf_, const FermionAction& fa_, LowModes& lm_)
  {
    TIMER_VERBOSE("Inv::setup(gf,fa,lm)");
    geo.set(geo_resize(gf_.geo()));
    gf.init();
    set_left_expanded_gauge_field(gf, gf_);
    fa = fa_;
    lm.init(lm_);
  }
  //
  double& stop_rsd() { return ip.stop_rsd; }
  const double& stop_rsd() const { return ip.stop_rsd; }
  //
  Long& max_num_iter() { return ip.max_num_iter; }
  const Long& max_num_iter() const { return ip.max_num_iter; }
  //
  Long& max_mixed_precision_cycle() { return ip.max_mixed_precision_cycle; }
  const Long& max_mixed_precision_cycle() const
  {
    return ip.max_mixed_precision_cycle;
  }
};

template <class Inv>
void setup_inverter(Inv& inv)
{
  inv.setup();
}

template <class Inv>
void setup_inverter(Inv& inv, const GaugeField& gf, const FermionAction& fa)
{
  inv.setup(gf, fa);
}

template <class Inv>
void setup_inverter(Inv& inv, const GaugeField& gf, const FermionAction& fa,
                    LowModes& lm)
{
  inv.setup(gf, fa, lm);
}

void multiply_m_dwf_no_comm(FermionField5d& out, const FermionField5d& in,
                            const InverterDomainWall& inv);

void multiply_m_dwf(FermionField5d& out, const FermionField5d& in,
                    const InverterDomainWall& inv);

void multiply_wilson_d_no_comm(FermionField5d& out, const FermionField5d& in,
                               const GaugeField& gf, const double mass);

void multiply_d_minus(FermionField5d& out, const FermionField5d& in,
                      const GaugeField& gf, const FermionAction& fa);

void multiply_d_minus(FermionField5d& out, const FermionField5d& in,
                      const InverterDomainWall& inv);

void multiply_m_full(FermionField5d& out, const FermionField5d& in,
                     const InverterDomainWall& inv);

void get_half_fermion(FermionField5d& half, const FermionField5d& ff,
                      const Int eo);

void set_half_fermion(FermionField5d& ff, const FermionField5d& half,
                      const Int eo);

void project_eo(FermionField5d& ff, const Int eo);

void multiply_m_e_e(FermionField5d& out, const FermionField5d& in,
                    const FermionAction& fa);

void multiply_mdag_e_e(FermionField5d& out, const FermionField5d& in,
                       const FermionAction& fa);

void multiply_m_e_e_inv(FermionField5d& out, const FermionField5d& in,
                        const FermionAction& fa);

void multiply_mdag_e_e_inv(FermionField5d& out, const FermionField5d& in,
                           const FermionAction& fa);

void multiply_wilson_d_e_o_no_comm(FermionField5d& out,
                                   const FermionField5d& in,
                                   const GaugeField& gf);

void multiply_wilson_ddag_e_o_no_comm(FermionField5d& out,
                                      const FermionField5d& in,
                                      const GaugeField& gf);

void multiply_m_e_o(FermionField5d& out, const FermionField5d& in,
                    const GaugeField& gf, const FermionAction& fa);

void multiply_mdag_e_o(FermionField5d& out, const FermionField5d& in,
                       const GaugeField& gf, const FermionAction& fa);

void multiply_m(FermionField5d& out, const FermionField5d& in,
                const GaugeField& gf, const FermionAction& fa);

void multiply_mdag(FermionField5d& out, const FermionField5d& in,
                   const GaugeField& gf, const FermionAction& fa);

void multiply_mpc_sym2(FermionField5d& out, const FermionField5d& in,
                       const GaugeField& gf, const FermionAction& fa);

void multiply_mpcdag_sym2(FermionField5d& out, const FermionField5d& in,
                          const GaugeField& gf, const FermionAction& fa);

void multiply_hermop_sym2(FermionField5d& out, const FermionField5d& in,
                          const GaugeField& gf, const FermionAction& fa);

void multiply_m_e_e(FermionField5d& out, const FermionField5d& in,
                    const InverterDomainWall& inv);

void multiply_mdag_e_e(FermionField5d& out, const FermionField5d& in,
                       const InverterDomainWall& inv);

void multiply_m_e_e_inv(FermionField5d& out, const FermionField5d& in,
                        const InverterDomainWall& inv);

void multiply_mdag_e_e_inv(FermionField5d& out, const FermionField5d& in,
                           const InverterDomainWall& inv);

void multiply_m_e_o(FermionField5d& out, const FermionField5d& in,
                    const InverterDomainWall& inv);

void multiply_mdag_e_o(FermionField5d& out, const FermionField5d& in,
                       const InverterDomainWall& inv);

void multiply_m_eo_eo(FermionField5d& out, const FermionField5d& in,
                      const InverterDomainWall& inv, const Int eo_out,
                      const Int eo_in);

void multiply_mdag_eo_eo(FermionField5d& out, const FermionField5d& in,
                         const InverterDomainWall& inv, const Int eo_out,
                         const Int eo_in);

void multiply_m(FermionField5d& out, const FermionField5d& in,
                const InverterDomainWall& inv);

void multiply_mdag(FermionField5d& out, const FermionField5d& in,
                   const InverterDomainWall& inv);

void multiply_mpc_sym2(FermionField5d& out, const FermionField5d& in,
                       const InverterDomainWall& inv);

void multiply_mpcdag_sym2(FermionField5d& out, const FermionField5d& in,
                          const InverterDomainWall& inv);

void multiply_hermop_sym2(FermionField5d& out, const FermionField5d& in,
                          const InverterDomainWall& inv);

void multiply_m_with_prec_sym2(FermionField5d& out, const FermionField5d& in,
                               const InverterDomainWall& inv);

ComplexD dot_product(const FermionField5d& ff1, const FermionField5d& ff2);

template <class Inv>
inline Long cg_with_f(
    FermionField5d& out, const FermionField5d& in, const Inv& inv,
    void f(FermionField5d&, const FermionField5d&, const Inv&),
    const double stop_rsd = 1e-8, const Long max_num_iter = 50000)
// f(out, in, inv);
{
  TIMER("cg_with_f");
  qassert(&out != &in);
  const Geometry geo = geo_resize(in.geo());
  if (not is_initialized(out)) {
    out.init(geo, in.multiplicity);
    set_zero(out);
  } else {
    out.init(geo, in.multiplicity);
  }
  if (max_num_iter == 0) {
    return 0;
  }
  FermionField5d r, p, tmp, ap;
  r.init(geo, in.multiplicity);
  p.init(geo, in.multiplicity);
  tmp.init(geo, in.multiplicity);
  r = in;
  f(tmp, out, inv);
  r -= tmp;
  p = r;
  const double qnorm_in = qnorm(in);
  displayln_info(
      fname +
      ssprintf(
          ": start max_num_iter=%4ld        sqrt(qnorm_in)=%.3E stop_rsd=%.3E",
          max_num_iter, sqrt(qnorm_in), stop_rsd));
  double qnorm_r = qnorm(r);
  for (Long iter = 1; iter <= max_num_iter; ++iter) {
    f(ap, p, inv);
    const double alpha = qnorm_r / dot_product(p, ap).real();
    tmp = p;
    tmp *= alpha;
    out += tmp;
    tmp = ap;
    tmp *= alpha;
    r -= tmp;
    const double new_qnorm_r = qnorm(r);
    if (is_cg_verbose()) {
      displayln_info(
          fname +
          ssprintf(": iter=%4ld sqrt(qnorm_r/qnorm_in)=%.3E stop_rsd=%.3E",
                   iter, sqrt(new_qnorm_r / qnorm_in), stop_rsd));
    }
    if (new_qnorm_r <= qnorm_in * sqr(stop_rsd)) {
      displayln_info(
          fname +
          ssprintf(
              ": final iter=%4ld sqrt(qnorm_r/qnorm_in)=%.3E stop_rsd=%.3E",
              iter, sqrt(new_qnorm_r / qnorm_in), stop_rsd));
      return iter;
    }
    const double beta = new_qnorm_r / qnorm_r;
    p *= beta;
    p += r;
    qnorm_r = new_qnorm_r;
  }
  displayln_info(fname + ssprintf(": final iter=%4ld (exceeded max_num_iter) "
                                  "sqrt(qnorm_r/qnorm_in)=%.3E stop_rsd=%.3E",
                                  max_num_iter + 1, sqrt(qnorm_r / qnorm_in),
                                  stop_rsd));
  return max_num_iter + 1;
}

template <class Inv>
void set_odd_prec_field_sym2(FermionField5d& in_o_p, FermionField5d& out_e_p,
                             const FermionField5d& in, const Inv& inv)
{
  TIMER_VERBOSE("set_odd_prec_field_sym2");
  qassert(&in_o_p != &out_e_p);
  qassert(&in_o_p != &in);
  qassert(&out_e_p != &in);
  FermionField5d in_e;
  get_half_fermion(in_e, in, 2);
  get_half_fermion(in_o_p, in, 1);
  FermionField5d tmp;
  multiply_m_e_e_inv(out_e_p, in_e, inv);
  multiply_m_e_o(tmp, out_e_p, inv);
  in_o_p -= tmp;
  multiply_mpcdag_sym2(in_o_p, in_o_p, inv);
}

template <class Inv>
void restore_field_from_odd_prec_sym2(FermionField5d& out,
                                      const FermionField5d& out_o_p,
                                      const FermionField5d& out_e_p,
                                      const Inv& inv)
{
  TIMER_VERBOSE("restore_field_from_odd_prec_sym2");
  qassert(&out != &out_o_p);
  qassert(&out != &out_e_p);
  FermionField5d out_e, out_o;
  out_e = out_e_p;
  out_o = out_o_p;
  multiply_m_e_e_inv(out_o, out_o, inv);
  FermionField5d tmp;
  multiply_m_e_o(tmp, out_o, inv);
  multiply_m_e_e_inv(tmp, tmp, inv);
  out_e -= tmp;
  set_half_fermion(out, out_e, 2);
  set_half_fermion(out, out_o, 1);
}

template <class Inv>
Long invert_with_cg(FermionField5d& out, const FermionField5d& in,
                    const Inv& inv,
                    Long cg(FermionField5d&, const FermionField5d&, const Inv&,
                            const double, const Long),
                    double stop_rsd = -1, Long max_num_iter = -1,
                    Long max_mixed_precision_cycle = -1,
                    bool dminus_multiplied_already = false)
// cg(out, in, inv, stop_rsd, max_iter)
{
  TIMER_VERBOSE_FLOPS("invert_with_cg(5d,5d,inv,cg)");
  qassert(&out != &in);
  if (stop_rsd < 0) {
    stop_rsd = inv.stop_rsd();
  }
  if (max_num_iter < 0) {
    max_num_iter = inv.max_num_iter();
  }
  if (max_mixed_precision_cycle < 0) {
    max_mixed_precision_cycle = inv.max_mixed_precision_cycle();
  }
  if (not is_initialized(out)) {
    out.init(geo_resize(in.geo()), in.multiplicity);
    set_zero(out);
  } else {
    out.init(geo_resize(in.geo()), in.multiplicity);
  }
  FermionField5d dm_in;
  if (not dminus_multiplied_already and inv.fa.is_multiplying_dminus) {
    multiply_d_minus(dm_in, in, inv);
  } else {
    dm_in.init(geo_resize(in.geo()), in.multiplicity);
    dm_in = in;
  }
  const double dm_in_qnorm = qnorm(dm_in);
  displayln_info(fname +
                 ssprintf(": dm_in sqrt(qnorm) = %E", sqrt(dm_in_qnorm)));
  if (dm_in_qnorm == 0.0) {
    displayln_info(fname + ssprintf(": WARNING: dm_in qnorm is zero."));
    out.init(geo_resize(in.geo()), in.multiplicity);
    set_zero(out);
    return 0;
  }
  Long total_iter = 0;
  if (inv.fa.is_using_zmobius == true and inv.fa.cg_diagonal_mee == 2) {
    FermionField5d in_o_p, out_e_p, out_o_p;
    set_odd_prec_field_sym2(in_o_p, out_e_p, dm_in, inv);
    out_o_p.init(in_o_p.geo(), in_o_p.multiplicity);
    set_zero(out_o_p);
    const double qnorm_in_o_p = qnorm(in_o_p);
    FermionField5d tmp, itmp;
    tmp.init(in_o_p.geo(), in_o_p.multiplicity);
    itmp = in_o_p;
    double qnorm_itmp = qnorm_in_o_p;
    Int cycle;
    for (cycle = 1; cycle <= max_mixed_precision_cycle; ++cycle) {
      if (not inv.lm.null() and inv.lm().initialized) {
        deflate(tmp, itmp, inv.lm.cast_const());
      } else {
        set_zero(tmp);
      }
      const Long iter =
          cg(tmp, itmp, inv, stop_rsd * sqrt(qnorm_in_o_p / qnorm_itmp),
             max_num_iter);
      total_iter += iter;
      out_o_p += tmp;
      if (iter <= max_num_iter) {
        itmp.init();
        break;
      }
      multiply_hermop_sym2(itmp, out_o_p, inv);
      itmp *= -1.0;
      itmp += in_o_p;
      qnorm_itmp = qnorm(itmp);
    }
    timer.flops += 5500 * total_iter * inv.fa.ls * inv.geo().local_volume();
    displayln_info(fname + ssprintf(": total_iter=%ld cycle=%d stop_rsd=%.3E",
                                    total_iter, cycle, stop_rsd));
    restore_field_from_odd_prec_sym2(out, out_o_p, out_e_p, inv);
  } else {
    qassert(false);
  }
  if (is_checking_invert()) {
    FermionField5d tmp;
    multiply_m(tmp, out, inv);
    tmp -= dm_in;
    displayln_info(fname + ssprintf(": checking %E from %E", sqrt(qnorm(tmp)),
                                    sqrt(qnorm(dm_in))));
  }
  return total_iter;
}

template <class Inv>
Long invert_with_cg_with_guess(FermionField5d& out, const FermionField5d& in,
                               const Inv& inv,
                               Long cg(FermionField5d&, const FermionField5d&,
                                       const Inv&, const double, const Long),
                               double stop_rsd = -1, Long max_num_iter = -1,
                               Long max_mixed_precision_cycle = -1,
                               bool dminus_multiplied_already = false)
// cg(out, in, inv, stop_rsd, max_iter)
{
  TIMER_VERBOSE("invert_with_cg_with_guess");
  qassert(&out != &in);
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
    dm_in.init(geo_resize(in.geo()), in.multiplicity);
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
  const Long total_iter =
      invert_with_cg(tmp, dm_in, inv, cg, stop_rsd, max_num_iter,
                     max_mixed_precision_cycle, true);
  out += tmp;
  return total_iter;
}

Long cg_with_herm_sym_2(FermionField5d& sol, const FermionField5d& src,
                        const InverterDomainWall& inv,
                        const double stop_rsd = 1e-8,
                        const Long max_num_iter = 50000);

Long invert(FermionField5d& out, const FermionField5d& in,
            const InverterDomainWall& inv);

Long invert(FermionField4d& out, const FermionField4d& in,
            const InverterDomainWall& inv);

double find_max_eigen_value_hermop_sym2(const InverterDomainWall& inv,
                                        const RngState& rs,
                                        const Long max_iter = 100);

}  // namespace qlat
