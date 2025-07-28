#pragma once

#include <qlat/fermion-action.h>
#include <qlat/qcd-utils.h>
#include <qlat/qcd.h>

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
  void init(const Geometry& geo, const Int multiplicity, const Coordinate& block_site, const Long neig,
            const Long nkeep)
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

inline Long load_low_modes(LowModes& lm, const std::string& path)
{
  TIMER_VERBOSE("load_low_modes");
  lm.initialized = false;
  if (path == "/dev/null" or path == "") {
    return 0;
  } else {
    const Long total_bytes = load_compressed_eigen_vectors(
        lm.eigen_values, lm.cesi, lm.cesb, lm.cesc, path);
    if (0 != total_bytes) {
      lm.initialized = true;
    }
    return total_bytes;
  }
}

inline Long load_or_compute_low_modes(LowModes& lm, const std::string& path,
                                      const GaugeField& gf,
                                      const FermionAction& fa,
                                      const LancArg& la)
// TODO: currently only load low modes
{
  TIMER_VERBOSE("load_or_compute_low_modes");
  (void)gf;
  (void)fa;
  (void)la;
  Long total_bytes = load_low_modes(lm, path);
  return total_bytes;
}

inline void load_low_modes_delay(LowModes& lm, const std::string& path)
{
  TIMER_VERBOSE("load_low_modes_delay");
  lm.initialized = true;
  lm.lmi.init();
  lm.lmi.initialized = true;
  lm.lmi.path = path;
}

inline void load_or_compute_low_modes_delay(LowModes& lm,
                                            const std::string& path,
                                            const GaugeField& gf,
                                            const FermionAction& fa,
                                            const LancArg& la)
{
  TIMER_VERBOSE("load_or_compute_low_modes_delay");
  lm.initialized = true;
  lm.lmi.init();
  lm.lmi.initialized = true;
  lm.lmi.path = path;
  lm.lmi.gf = gf;
  lm.lmi.fa = fa;
  lm.lmi.la = la;
}

inline Long force_low_modes(LowModes& lm)
{
  TIMER("force_low_modes");
  if (lm.lmi.initialized) {
    Long total_bytes = load_or_compute_low_modes(lm, lm.lmi.path, lm.lmi.gf,
                                                 lm.lmi.fa, lm.lmi.la);
    lm.lmi.init();
    return total_bytes;
  }
  return 0;
}

inline void set_u_rand(LowModes& lm, const RngState& rs)
{
  TIMER_VERBOSE("set_u_rand(lm,rs)");
  qassert(lm.initialized);
  force_low_modes(lm);
  set_u_rand(get_data(lm.eigen_values), rs.split("eigen_values"));
  set_u_rand(lm.cesb, rs.split("cesb"));
  set_u_rand(lm.cesc, rs.split("cesc"));
}

inline Long save_low_modes_decompress(LowModes& lm, const std::string& path)
{
  TIMER_VERBOSE("save_low_modes_decompress");
  qassert(lm.initialized);
  force_low_modes(lm);
  std::vector<BlockedHalfVector> bhvs;
  decompress_eigen_system(bhvs, lm.cesb, lm.cesc);
  std::vector<HalfVector> hvs;
  convert_half_vectors(hvs, bhvs);
  const int id_node = get_id_node();
  const int num_node = get_num_node();
  const int idx = id_node;
  const int idx_size = num_node;
  const int dir_idx = compute_dist_file_dir_id(idx, idx_size);
  qmkdir(path);
  qmkdir(path + ssprintf("/%02d", dir_idx));
  const std::string fn = path + ssprintf("/%02d/%010d", dir_idx, idx);
  Long total_bytes = 0;
  const int n_cycle = std::max(1, num_node / dist_write_par_limit());
  std::vector<crc32_t> crcs(num_node, 0);
  for (int i = 0; i < n_cycle; i++) {
    Long bytes = 0;
    if (id_node % n_cycle == i) {
      qassert(hvs.size() >= 1);
      bytes = hvs.size() * get_data(hvs[0]).data_size();
      crcs[id_node] = save_half_vectors(hvs, fn, false, true);
    }
    glb_sum(bytes);
    total_bytes += bytes;
    displayln_info(
        ssprintf("qlat::%s: cycle / n_cycle = %4d / %4d ; total_bytes = %15ld",
                 fname.c_str(), i + 1, n_cycle, total_bytes));
  }
  glb_sum_byte_vec(get_data(crcs));
  const crc32_t crc = dist_crc32(crcs);
  if (get_id_node() == 0) {
    const std::string fn = path + "/checksums.txt";
    QFile fp = qfopen(fn, "w");
    qassert(not fp.null());
    qwrite_data(ssprintf("%08X\n", crc), fp);
    qwrite_data("\n", fp);
    for (size_t i = 0; i < crcs.size(); ++i) {
      qwrite_data(ssprintf("%08X\n", crcs[i]), fp);
    }
    qfclose(fp);
  }
  if (get_id_node() == 0) {
    const std::string fn = path + "/eigen-values.txt";
    QFile fp = qfopen(fn, "w");
    qassert(not fp.null());
    qwrite_data(ssprintf("%ld\n", lm.eigen_values.size()), fp);
    for (Long i = 0; i < lm.eigen_values.size(); ++i) {
      qwrite_data(ssprintf("%.20lE\n", lm.eigen_values[i]), fp);
    }
    qfclose(fp);
  }
  qtouch_info(path + "/checkpoint");
  timer.flops += total_bytes;
  return total_bytes;
}

inline void deflate(HalfVector& hv_out, const HalfVector& hv_in, LowModes& lm)
// hv_out can be the same as hv_in
{
  force_low_modes(lm);
  if (not lm.initialized) {
    hv_out.init(geo_resize(hv_in.geo()), hv_in.multiplicity);
    set_zero(hv_out);
    return;
  }
  TIMER("deflate(hv,hv,lm)");
  const int ls = lm.cesi.ls;
  const Long block_size = lm.cesb.block_vol_eo * ls * HalfVector::c_size;
  const Long n_basis = lm.cesb.n_basis;
  const Long n_vec = lm.cesc.n_vec;
  qassert(n_vec == (Long)lm.eigen_values.size());
  BlockedHalfVector bhv;
  convert_half_vector(bhv, hv_in, lm.cesi.block_site);
  const Geometry geo = geo_resize(bhv.geo());
  Field<ComplexD> chv, phv;
  chv.init(geo, n_basis);
  phv.init(geo, n_vec);
  set_zero(chv);
  set_zero(phv);
  {
    TIMER("deflate-project");
#pragma omp parallel for
    for (Long index = 0; index < geo.local_volume(); ++index) {
      // const Coordinate xl = geo.coordinate_from_index(index);
      const Vector<ComplexF> vb = bhv.get_elems_const(index);
      const Vector<ComplexF> vbs = lm.cesb.get_elems_const(index);
      const Vector<ComplexF> vcs = lm.cesc.get_elems_const(index);
      // project to coarse grid
      Vector<ComplexD> vc = chv.get_elems(index);
      for (int j = 0; j < n_basis; ++j) {
        ComplexD& vc_j = vc.p[j];
        const Vector<ComplexF> vbs_j(vbs.p + j * block_size, block_size);
        for (Long k = 0; k < block_size; ++k) {
          const ComplexF& vb_k = vb.p[k];
          vc_j += qconj(vbs_j.p[k]) * vb_k;
        }
      }
      // compute inner products
      Vector<ComplexD> vp = phv.get_elems(index);
      for (int i = 0; i < n_vec; ++i) {
        ComplexD& vp_i = vp.p[i];
        const Vector<ComplexF> vcs_i(vcs.p + i * n_basis, n_basis);
        for (int j = 0; j < n_basis; ++j) {
          const ComplexD& vc_j = vc.p[j];
          vp_i += (ComplexD)qconj(vcs_i.p[j]) * vc_j;
        }
      }
    }
  }
  std::vector<ComplexD> phv_sum(n_vec, 0.0);
  {
    TIMER("deflate-glbsum");
    // glb sum inner products
    for (Long index = 0; index < geo.local_volume(); ++index) {
      // const Coordinate xl = geo.coordinate_from_index(index);
      Vector<ComplexD> vp = phv.get_elems(index);
#pragma omp parallel for
      for (int i = 0; i < n_vec; ++i) {
        phv_sum[i] += vp[i];
      }
    }
    phv.init();
    glb_sum_double_vec(get_data(phv_sum));
    // scale by eigen values
#pragma omp parallel for
    for (int i = 0; i < n_vec; ++i) {
      phv_sum[i] /= lm.eigen_values[i];
    }
  }
  {
    TIMER("deflate-produce");
    // producing coarse space vector
    set_zero(chv);
#pragma omp parallel for
    for (Long index = 0; index < geo.local_volume(); ++index) {
      // const Coordinate xl = geo.coordinate_from_index(index);
      const Vector<ComplexF> vbs = lm.cesb.get_elems_const(index);
      const Vector<ComplexF> vcs = lm.cesc.get_elems_const(index);
      // compute inner products
      Vector<ComplexD> vc = chv.get_elems(index);
      for (int i = 0; i < n_vec; ++i) {
        const ComplexD& phv_sum_i = phv_sum[i];
        const Vector<ComplexF> vcs_i(vcs.p + i * n_basis, n_basis);
        for (int j = 0; j < n_basis; ++j) {
          ComplexD& vc_j = vc.p[j];
          vc_j += (ComplexD)(vcs_i.p[j]) * phv_sum_i;
        }
      }
      // project to fine grid
      Vector<ComplexF> vb = bhv.get_elems(index);
      for (int j = 0; j < n_basis; ++j) {
        const ComplexF& vc_j = (ComplexF)vc[j];
        const Vector<ComplexF> vbs_j(vbs.p + j * block_size, block_size);
        for (Long k = 0; k < block_size; ++k) {
          ComplexF& vb_k = vb.p[k];
          vb_k += vbs_j.p[k] * vc_j;
        }
      }
    }
  }
  convert_half_vector(hv_out, bhv);
}

inline void deflate(FermionField5d& out, const FermionField5d& in, LowModes& lm)
// out can be the same as in
{
  force_low_modes(lm);
  if (not lm.initialized) {
    out.init(geo_resize(in.geo()), in.multiplicity);
    set_zero(out);
    return;
  }
  TIMER_VERBOSE("deflate(5d,5d,lm)");
  const Geometry& geo = geo_resize(in.geo());
  qassert(geo.eo == 1 or geo.eo == 2);
  const int ls = in.multiplicity;
  qassert(ls == lm.cesi.ls);
  HalfVector hv;
  init_half_vector(hv, geo, ls);
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    qassert((xl[0] + xl[1] + xl[2] + xl[3]) % 2 == 2 - geo.eo);
    Vector<ComplexF> vhv = hv.get_elems(index);
    const Vector<WilsonVector> vin = in.get_elems_const(index);
    qassert(vhv.size() ==
            vin.size() * (Long)sizeof(WilsonVector) / (Long)sizeof(ComplexD));
    const Vector<ComplexD> vff((const ComplexD*)vin.data(), vhv.size());
    for (int m = 0; m < vhv.size(); ++m) {
      vhv[m] = vff[m];
    }
  }
  deflate(hv, hv, lm);
  if (is_initialized(out) and out.geo().eo == 3 - geo.eo) {
    out.geo().eo = geo.eo;
  }
  out.init(geo, in.multiplicity);
  qassert(out.geo().eo == geo.eo);
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    qassert((xl[0] + xl[1] + xl[2] + xl[3]) % 2 == 2 - geo.eo);
    const Vector<ComplexF> vhv = hv.get_elems(index);
    Vector<WilsonVector> vout = out.get_elems(index);
    qassert(vhv.size() ==
            vout.size() * (Long)sizeof(WilsonVector) / (Long)sizeof(ComplexD));
    Vector<ComplexD> vff((ComplexD*)vout.data(), vhv.size());
    for (int m = 0; m < vhv.size(); ++m) {
      vff[m] = vhv[m];
    }
  }
}

inline void benchmark_deflate(const Geometry& geo, const int ls,
                              const Coordinate& block_site, const Long neig,
                              const Long nkeep, const RngState& rs)
{
  TIMER_VERBOSE("benchmark_deflate");
  displayln_info(ssprintf("geo = %s", show(geo).c_str()));
  displayln_info(ssprintf("block_site = %s", show(block_site).c_str()));
  displayln_info(ssprintf("ls = %d, neig = %d, nkeep = %d", ls, neig, nkeep));
  LowModes lm;
  lm.init(geo, ls, block_site, neig, nkeep);
  set_u_rand(lm, rs.split("lm"));
  FermionField5d in, out;
  in.init(geo_eo(geo, 1), ls);
  out.init(geo_eo(geo, 1), ls);
  set_u_rand(in, rs.split("in"));
  sync_node();
  {
    TIMER_VERBOSE("benchmark_deflate-deflate");
    for (int i = 0; i < 4; ++i) {
      set_zero(out);
      deflate(out, in, lm);
    }
  }
  sync_node();
}

struct InverterParams {
  double stop_rsd;
  Long max_num_iter;
  Long max_mixed_precision_cycle;
  int solver_type;  // 0 -> CG, 1-> EIGCG, 2->MSPCG
  int higher_precision;
  int lower_precision;
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

inline void multiply_m_dwf_no_comm(FermionField5d& out,
                                   const FermionField5d& in,
                                   const InverterDomainWall& inv)
{
  TIMER("multiply_m_dwf_no_comm(5d,5d,Inv)");
  qassert(&out != &in);
  const Geometry geo = geo_resize(in.geo());
  qassert(is_matching_geo(inv.geo(), geo));
  out.init(geo, in.multiplicity);
  set_zero(out);
  const FermionAction& fa = inv.fa;
  const GaugeField& gf = inv.gf;
  qassert(in.multiplicity == fa.ls);
  qassert(out.multiplicity == fa.ls);
  qassert(fa.mobius_scale == 1.0);
  qassert((int)fa.bs.size() == fa.ls);
  qassert((int)fa.cs.size() == fa.ls);
  const array<SpinMatrix, 4>& gammas =
      SpinMatrixConstants::get_cps_gammas();
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  const SpinMatrix& unit = SpinMatrixConstants::get_unit();
  const SpinMatrix p_p = (ComplexD)0.5 * (unit + gamma5);
  const SpinMatrix p_m = (ComplexD)0.5 * (unit - gamma5);
  array<SpinMatrix, 4> p_mu_p;
  array<SpinMatrix, 4> p_mu_m;
  for (int mu = 0; mu < 4; ++mu) {
    p_mu_p[mu] = (ComplexD)0.5 * (unit + gammas[mu]);
    p_mu_m[mu] = (ComplexD)0.5 * (unit - gammas[mu]);
  }
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<WilsonVector> v = out.get_elems(xl);
    {
      const Vector<WilsonVector> iv = in.get_elems_const(xl);
      for (int m = 0; m < fa.ls; ++m) {
        v[m] = (ComplexD)(5.0 - fa.m5) * iv[m];
        v[m] -= p_m * (m < fa.ls - 1
                           ? iv[m + 1]
                           : (WilsonVector)((ComplexD)(-fa.mass) * iv[0]));
        v[m] -= p_p *
                (m > 0 ? iv[m - 1]
                       : (WilsonVector)((ComplexD)(-fa.mass) * iv[fa.ls - 1]));
      }
    }
    for (int mu = 0; mu < 4; ++mu) {
      const Coordinate xl_p = coordinate_shifts(xl, mu);
      const Coordinate xl_m = coordinate_shifts(xl, -mu - 1);
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

inline void multiply_m_dwf(FermionField5d& out, const FermionField5d& in,
                           const InverterDomainWall& inv)
// out can be the same object as in
{
  TIMER("multiply_m_dwf(5d,5d,Inv)");
  const Geometry geo1 = geo_resize(in.geo(), 1);
  FermionField5d in1;
  in1.init(geo1, in.multiplicity);
  in1 = in;
  refresh_expanded_1(in1);
  multiply_m_dwf_no_comm(out, in1, inv);
}

inline void multiply_wilson_d_no_comm(FermionField5d& out,
                                      const FermionField5d& in,
                                      const GaugeField& gf, const double mass)
// set_left_expanded_gauge_field(gf, gf_);
// in.geo() = geo_resize(geo, 1);
// in.multiplicity == ls;
// refresh_expanded_1(in);
{
  TIMER("multiply_wilson_d_no_comm(5d,5d,gf,mass)");
  qassert(&out != &in);
  const Geometry geo = geo_resize(in.geo());
  qassert(is_matching_geo(gf.geo(), geo));
  out.init(geo, in.multiplicity);
  set_zero(out);
  const int ls = in.multiplicity;
  qassert(out.multiplicity == ls);
  const array<SpinMatrix, 4>& gammas =
      SpinMatrixConstants::get_cps_gammas();
  const SpinMatrix& unit = SpinMatrixConstants::get_unit();
  array<SpinMatrix, 4> p_mu_p;
  array<SpinMatrix, 4> p_mu_m;
  for (int mu = 0; mu < 4; ++mu) {
    p_mu_p[mu] = (ComplexD)0.5 * (unit + gammas[mu]);
    p_mu_m[mu] = (ComplexD)0.5 * (unit - gammas[mu]);
  }
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<WilsonVector> v = out.get_elems(xl);
    {
      const Vector<WilsonVector> iv = in.get_elems_const(xl);
      for (int m = 0; m < ls; ++m) {
        v[m] = (ComplexD)(4.0 + mass) * iv[m];
      }
    }
    for (int mu = 0; mu < 4; ++mu) {
      const Coordinate xl_p = coordinate_shifts(xl, mu);
      const Coordinate xl_m = coordinate_shifts(xl, -mu - 1);
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

inline void multiply_d_minus(FermionField5d& out, const FermionField5d& in,
                             const GaugeField& gf, const FermionAction& fa)
// out can be the same object as in
{
  TIMER("multiply_d_minus(5d,5d,gf,fa)");
  const Geometry geo = geo_resize(in.geo());
  qassert(is_matching_geo(gf.geo(), in.geo()));
  qassert(in.multiplicity == fa.ls);
  qassert((int)fa.bs.size() == fa.ls);
  qassert((int)fa.cs.size() == fa.ls);
  const Geometry geo1 = geo_resize(in.geo(), 1);
  FermionField5d in1, out1;
  in1.init(geo1, in.multiplicity);
  out1.init(geo, in.multiplicity);
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<WilsonVector> v = out1.get_elems(xl);
    Vector<WilsonVector> v1 = in1.get_elems(xl);
    const Vector<WilsonVector> iv = in.get_elems_const(xl);
    for (int m = 0; m < fa.ls; ++m) {
      const ComplexD& c = fa.cs[m];
      v1[m] = (ComplexD)(-c) * iv[m];
      v[m] = iv[m];
    }
  }
  refresh_expanded_1(in1);
  out.init(geo, in.multiplicity);
  qassert(out.multiplicity == fa.ls);
  qassert(is_matching_geo(gf.geo(), out.geo()));
  multiply_wilson_d_no_comm(out, in1, gf, -fa.m5);
  out += out1;
}

inline void multiply_d_minus(FermionField5d& out, const FermionField5d& in,
                             const InverterDomainWall& inv)
// out can be the same object as in
{
  TIMER("multiply_d_minus(5d,5d,Inv)");
  multiply_d_minus(out, in, inv.gf, inv.fa);
}

inline void multiply_m_full(FermionField5d& out, const FermionField5d& in,
                            const InverterDomainWall& inv)
{
  TIMER("multiply_m_full(5d,5d,Inv)");
  qassert(&out != &in);
  const Geometry geo = geo_resize(in.geo());
  out.init(geo, in.multiplicity);
  set_zero(out);
  const FermionAction& fa = inv.fa;
  qassert(is_matching_geo(inv.geo(), in.geo()));
  qassert(is_matching_geo(inv.geo(), out.geo()));
  qassert(in.multiplicity == fa.ls);
  qassert(in.multiplicity == fa.ls);
  qassert(out.multiplicity == fa.ls);
  qassert((int)fa.bs.size() == fa.ls);
  qassert((int)fa.cs.size() == fa.ls);
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  const SpinMatrix& unit = SpinMatrixConstants::get_unit();
  const SpinMatrix p_p = (ComplexD)0.5 * (unit + gamma5);
  const SpinMatrix p_m = (ComplexD)0.5 * (unit - gamma5);
  const GaugeField& gf = inv.gf;
  const Geometry geo1 = geo_resize(in.geo(), 1);
  FermionField5d in1, fftmp;
  in1.init(geo1, in.multiplicity);
  fftmp.init(geo, in.multiplicity);
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<WilsonVector> iv = in.get_elems_const(xl);
    Vector<WilsonVector> v = fftmp.get_elems(xl);
    Vector<WilsonVector> v1 = in1.get_elems(xl);
    for (int m = 0; m < fa.ls; ++m) {
      const ComplexD& b = fa.bs[m];
      const ComplexD& c = fa.cs[m];
      v1[m] = (ComplexD)b * iv[m];
      v[m] = iv[m];
      const WilsonVector tmp =
          (p_m * (m < fa.ls - 1
                      ? iv[m + 1]
                      : (WilsonVector)((ComplexD)(-fa.mass) * iv[0]))) +
          (p_p * (m > 0
                      ? iv[m - 1]
                      : (WilsonVector)((ComplexD)(-fa.mass) * iv[fa.ls - 1])));
      v1[m] += (ComplexD)c * tmp;
      v[m] -= tmp;
    }
  }
  refresh_expanded_1(in1);
  multiply_wilson_d_no_comm(out, in1, gf, -fa.m5);
  out += fftmp;
}

inline void get_half_fermion(FermionField5d& half, const FermionField5d& ff,
                             const int eo)
// 2:even 1:odd
{
  TIMER("get_half_fermion");
  qassert(&half != &ff);
  Geometry geoh = geo_resize(ff.geo());
  geoh.eo = eo;
  half.init(geoh, ff.multiplicity);
  qassert(half.geo().eo == eo);
  qassert(ff.geo().eo == 0);
  qassert(is_matching_geo(ff.geo(), half.geo()));
#pragma omp parallel for
  for (Long index = 0; index < geoh.local_volume(); ++index) {
    const Coordinate xl = geoh.coordinate_from_index(index);
    assign(half.get_elems(xl), ff.get_elems_const(xl));
  }
}

inline void set_half_fermion(FermionField5d& ff, const FermionField5d& half,
                             const int eo)
// 2:even 1:odd
{
  TIMER("set_half_fermion");
  qassert(&half != &ff);
  const Geometry geoh = half.geo();
  Geometry geo = geo_resize(geoh);
  geo.eo = 0;
  ff.init(geo, half.multiplicity);
  qassert(half.geo().eo == eo);
  qassert(ff.geo().eo == 0);
  qassert(is_matching_geo(ff.geo(), half.geo()));
#pragma omp parallel for
  for (Long index = 0; index < geoh.local_volume(); ++index) {
    const Coordinate xl = geoh.coordinate_from_index(index);
    assign(ff.get_elems(xl), half.get_elems_const(xl));
  }
}

inline void project_eo(FermionField5d& ff, const int eo)
{
  TIMER("project_eo");
  qassert(eo == 1 or eo == 2);
  FermionField5d half;
  get_half_fermion(half, ff, eo);
  set_zero(ff);
  set_half_fermion(ff, half, eo);
}

inline void multiply_m_e_e(FermionField5d& out, const FermionField5d& in,
                           const FermionAction& fa)
// out can be the same object as in
// works for _o_o as well
{
  TIMER("multiply_m_e_e");
  if (is_initialized(out) and out.geo().eo == 3 - in.geo().eo) {
    out.geo().eo = in.geo().eo;
  }
  out.init(geo_resize(in.geo()), in.multiplicity);
  qassert(is_matching_geo(out.geo(), in.geo()));
  qassert(in.geo().eo == 1 or in.geo().eo == 2);
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  const SpinMatrix& unit = SpinMatrixConstants::get_unit();
  const SpinMatrix p_p = (ComplexD)0.5 * (unit + gamma5);
  const SpinMatrix p_m = (ComplexD)0.5 * (unit - gamma5);
  FermionField5d in_copy;
  ConstHandle<FermionField5d> hin;
  if (&out != &in) {
    hin.init(in);
  } else {
    in_copy.init(geo_resize(in.geo()), in.multiplicity);
    in_copy = in;
    hin.init(in_copy);
  }
  const Geometry& geo = out.geo();
  std::vector<ComplexD> bee(fa.ls), cee(fa.ls);
  for (int m = 0; m < fa.ls; ++m) {
    bee[m] = 1.0 + fa.bs[m] * (4.0 - fa.m5);
    cee[m] = 1.0 - fa.cs[m] * (4.0 - fa.m5);
  }
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<WilsonVector> iv = hin().get_elems_const(xl);
    Vector<WilsonVector> v = out.get_elems(xl);
    for (int m = 0; m < fa.ls; ++m) {
      v[m] = (ComplexD)bee[m] * iv[m];
      const WilsonVector tmp =
          (p_m * (m < fa.ls - 1
                      ? iv[m + 1]
                      : (WilsonVector)((ComplexD)(-fa.mass) * iv[0]))) +
          (p_p * (m > 0
                      ? iv[m - 1]
                      : (WilsonVector)((ComplexD)(-fa.mass) * iv[fa.ls - 1])));
      v[m] -= (ComplexD)cee[m] * tmp;
    }
  }
}

inline void multiply_mdag_e_e(FermionField5d& out, const FermionField5d& in,
                              const FermionAction& fa)
// out can be the same object as in
// works for _o_o as well
{
  TIMER("multiply_mdag_e_e");
  if (is_initialized(out) and out.geo().eo == 3 - in.geo().eo) {
    out.geo().eo = in.geo().eo;
  }
  out.init(geo_resize(in.geo()), in.multiplicity);
  qassert(is_matching_geo(out.geo(), in.geo()));
  qassert(in.geo().eo == 1 or in.geo().eo == 2);
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  const SpinMatrix& unit = SpinMatrixConstants::get_unit();
  const SpinMatrix p_p = (ComplexD)0.5 * (unit + gamma5);
  const SpinMatrix p_m = (ComplexD)0.5 * (unit - gamma5);
  FermionField5d in_copy;
  ConstHandle<FermionField5d> hin;
  if (&out != &in) {
    hin.init(in);
  } else {
    in_copy.init(geo_resize(in.geo()), in.multiplicity);
    in_copy = in;
    hin.init(in_copy);
  }
  const Geometry& geo = out.geo();
  std::vector<ComplexD> bee(fa.ls), cee(fa.ls);
  for (int m = 0; m < fa.ls; ++m) {
    bee[m] = qconj(1.0 + fa.bs[m] * (4.0 - fa.m5));
    cee[m] = qconj(1.0 - fa.cs[m] * (4.0 - fa.m5));
  }
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<WilsonVector> iv = hin().get_elems_const(xl);
    Vector<WilsonVector> v = out.get_elems(xl);
    for (int m = 0; m < fa.ls; ++m) {
      v[m] = (ComplexD)bee[m] * iv[m];
      const WilsonVector tmp =
          (p_p * (m < fa.ls - 1
                      ? (WilsonVector)((ComplexD)cee[m + 1] * iv[m + 1])
                      : (WilsonVector)(
                            (ComplexD)(-qconj((ComplexD)fa.mass) * cee[0]) *
                            iv[0]))) +
          (p_m * (m > 0
                      ? (WilsonVector)((ComplexD)cee[m - 1] * iv[m - 1])
                      : (WilsonVector)((ComplexD)(-qconj((ComplexD)fa.mass) *
                                                  cee[fa.ls - 1]) *
                                       iv[fa.ls - 1])));
      v[m] -= tmp;
    }
  }
}

inline void multiply_m_e_e_inv(FermionField5d& out, const FermionField5d& in,
                               const FermionAction& fa)
// out can be the same object as in
// works for _o_o as well
{
  TIMER("multiply_m_e_e_inv");
  if (is_initialized(out) and out.geo().eo == 3 - in.geo().eo) {
    out.geo().eo = in.geo().eo;
  }
  out.init(geo_resize(in.geo()), in.multiplicity);
  qassert(is_matching_geo(out.geo(), in.geo()));
  qassert(out.geo().eo == in.geo().eo);
  qassert(in.geo().eo == 1 or in.geo().eo == 2);
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  const SpinMatrix& unit = SpinMatrixConstants::get_unit();
  const SpinMatrix p_p = (ComplexD)0.5 * (unit + gamma5);
  const SpinMatrix p_m = (ComplexD)0.5 * (unit - gamma5);
  const Geometry& geo = out.geo();
  std::vector<ComplexD> bee(fa.ls), cee(fa.ls);
  for (int m = 0; m < fa.ls; ++m) {
    bee[m] = 1.0 + fa.bs[m] * (4.0 - fa.m5);
    cee[m] = 1.0 - fa.cs[m] * (4.0 - fa.m5);
  }
  std::vector<ComplexD> lee(fa.ls - 1), leem(fa.ls - 1);
  for (int m = 0; m < fa.ls - 1; ++m) {
    lee[m] = -cee[m + 1] / bee[m];
    leem[m] = m == 0 ? fa.mass * cee[fa.ls - 1] / bee[0]
                     : leem[m - 1] * cee[m - 1] / bee[m];
  }
  std::vector<ComplexD> dee(fa.ls, 0.0);
  dee[fa.ls - 1] = fa.mass * cee[fa.ls - 1];
  for (int m = 0; m < fa.ls - 1; ++m) {
    dee[fa.ls - 1] *= cee[m] / bee[m];
  }
  for (int m = 0; m < fa.ls; ++m) {
    dee[m] += bee[m];
  }
  std::vector<ComplexD> uee(fa.ls - 1), ueem(fa.ls - 1);
  for (int m = 0; m < fa.ls - 1; ++m) {
    uee[m] = -cee[m] / bee[m];
    ueem[m] =
        m == 0 ? fa.mass * cee[0] / bee[0] : ueem[m - 1] * cee[m] / bee[m];
  }
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<WilsonVector> iv = in.get_elems_const(xl);
    Vector<WilsonVector> v = out.get_elems(xl);
    if (v.data() != iv.data()) {
      std::memcpy(v.data(), iv.data(), iv.data_size());
    }
    WilsonVector tmp;
    // {L^m_{ee}}^{-1}
    set_zero(tmp);
    for (int m = 0; m < fa.ls - 1; ++m) {
      tmp += (ComplexD)(-leem[m]) * v[m];
    }
    v[fa.ls - 1] += p_m * tmp;
    // {L'_{ee}}^{-1}
    for (int m = 1; m < fa.ls; ++m) {
      v[m] += (ComplexD)(-lee[m - 1]) * (p_p * v[m - 1]);
    }
    // {D_{ee}}^{-1}
    for (int m = 0; m < fa.ls; ++m) {
      v[m] *= (ComplexD)(1.0 / dee[m]);
    }
    // {U^'_{ee}}^{-1}
    for (int m = fa.ls - 2; m >= 0; --m) {
      v[m] += (ComplexD)(-uee[m]) * (p_m * v[m + 1]);
    }
    // {U^m_{ee}}^{-1}
    for (int m = 0; m < fa.ls - 1; ++m) {
      v[m] += (ComplexD)(-ueem[m]) * (p_p * v[fa.ls - 1]);
    }
  }
}

inline void multiply_mdag_e_e_inv(FermionField5d& out, const FermionField5d& in,
                                  const FermionAction& fa)
// out can be the same object as in
// works for _o_o as well
{
  TIMER("multiply_mdag_e_e_inv");
  if (is_initialized(out) and out.geo().eo == 3 - in.geo().eo) {
    out.geo().eo = in.geo().eo;
  }
  out.init(geo_resize(in.geo()), in.multiplicity);
  qassert(is_matching_geo(out.geo(), in.geo()));
  qassert(out.geo().eo == in.geo().eo);
  qassert(in.geo().eo == 1 or in.geo().eo == 2);
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  const SpinMatrix& unit = SpinMatrixConstants::get_unit();
  const SpinMatrix p_p = (ComplexD)0.5 * (unit + gamma5);
  const SpinMatrix p_m = (ComplexD)0.5 * (unit - gamma5);
  const Geometry& geo = out.geo();
  std::vector<ComplexD> bee(fa.ls), cee(fa.ls);
  for (int m = 0; m < fa.ls; ++m) {
    bee[m] = 1.0 + fa.bs[m] * (4.0 - fa.m5);
    cee[m] = 1.0 - fa.cs[m] * (4.0 - fa.m5);
  }
  std::vector<ComplexD> lee(fa.ls - 1), leem(fa.ls - 1);
  for (int m = 0; m < fa.ls - 1; ++m) {
    lee[m] = -cee[m + 1] / bee[m];
    leem[m] = m == 0 ? fa.mass * cee[fa.ls - 1] / bee[0]
                     : leem[m - 1] * cee[m - 1] / bee[m];
  }
  std::vector<ComplexD> dee(fa.ls, 0.0);
  dee[fa.ls - 1] = fa.mass * cee[fa.ls - 1];
  for (int m = 0; m < fa.ls - 1; ++m) {
    dee[fa.ls - 1] *= cee[m] / bee[m];
  }
  for (int m = 0; m < fa.ls; ++m) {
    dee[m] += bee[m];
  }
  std::vector<ComplexD> uee(fa.ls - 1), ueem(fa.ls - 1);
  for (int m = 0; m < fa.ls - 1; ++m) {
    uee[m] = -cee[m] / bee[m];
    ueem[m] =
        m == 0 ? fa.mass * cee[0] / bee[0] : ueem[m - 1] * cee[m] / bee[m];
  }
  for (int m = 0; m < fa.ls; ++m) {
    bee[m] = qconj(bee[m]);
    cee[m] = qconj(cee[m]);
    dee[m] = qconj(dee[m]);
  }
  for (int m = 0; m < fa.ls - 1; ++m) {
    lee[m] = qconj(lee[m]);
    leem[m] = qconj(leem[m]);
    uee[m] = qconj(uee[m]);
    ueem[m] = qconj(ueem[m]);
  }
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<WilsonVector> iv = in.get_elems_const(xl);
    Vector<WilsonVector> v = out.get_elems(xl);
    if (v.data() != iv.data()) {
      std::memcpy(v.data(), iv.data(), iv.data_size());
    }
    WilsonVector tmp;
    // {U^m_{ee}}^\dagger^{-1}
    set_zero(tmp);
    for (int m = 0; m < fa.ls - 1; ++m) {
      tmp += (ComplexD)(-ueem[m]) * v[m];
    }
    v[fa.ls - 1] += p_p * tmp;
    // {U^'_{ee}}^\dagger^{-1}
    for (int m = 1; m < fa.ls; ++m) {
      v[m] += (ComplexD)(-uee[m - 1]) * (p_m * v[m - 1]);
    }
    // {D_{ee}}^\dagger^{-1}
    for (int m = 0; m < fa.ls; ++m) {
      v[m] *= (ComplexD)(1.0 / dee[m]);
    }
    // {L'_{ee}}^\dagger^{-1}
    for (int m = fa.ls - 2; m >= 0; --m) {
      v[m] += (ComplexD)(-lee[m]) * (p_p * v[m + 1]);
    }
    // {L^m_{ee}}^\dagger^{-1}
    for (int m = 0; m < fa.ls - 1; ++m) {
      v[m] += (ComplexD)(-leem[m]) * (p_m * v[fa.ls - 1]);
    }
  }
}

inline void multiply_wilson_d_e_o_no_comm(FermionField5d& out,
                                          const FermionField5d& in,
                                          const GaugeField& gf)
// set_left_expanded_gauge_field(gf, gf_);
// in.geo() = geo_resize(geo, 1);
// in.multiplicity = ls;
// refresh_expanded_1(in);
{
  TIMER("multiply_wilson_d_e_o_no_comm(5d,5d,gf)");
  qassert(&out != &in);
  qassert(is_matching_geo(gf.geo(), in.geo()));
  qassert(in.geo().eo == 1 or in.geo().eo == 2);
  Geometry geo = geo_resize(in.geo());
  geo.eo = 3 - in.geo().eo;
  if (is_initialized(out) and out.geo().eo == 3 - geo.eo) {
    out.geo().eo = geo.eo;
  }
  out.init(geo, in.multiplicity);
  set_zero(out);
  const int ls = in.multiplicity;
  qassert(out.multiplicity == ls);
  qassert(is_matching_geo(out.geo(), in.geo()));
  qassert(out.geo().eo != in.geo().eo);
  qassert(out.geo().eo == 1 or out.geo().eo == 2);
  const array<SpinMatrix, 4>& gammas =
      SpinMatrixConstants::get_cps_gammas();
  const SpinMatrix& unit = SpinMatrixConstants::get_unit();
  array<SpinMatrix, 4> p_mu_p;
  array<SpinMatrix, 4> p_mu_m;
  for (int mu = 0; mu < 4; ++mu) {
    p_mu_p[mu] = (ComplexD)0.5 * (unit + gammas[mu]);
    p_mu_m[mu] = (ComplexD)0.5 * (unit - gammas[mu]);
  }
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<WilsonVector> v = out.get_elems(xl);
    for (int mu = 0; mu < 4; ++mu) {
      const Coordinate xl_p = coordinate_shifts(xl, mu);
      const Coordinate xl_m = coordinate_shifts(xl, -mu - 1);
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

inline void multiply_wilson_ddag_e_o_no_comm(FermionField5d& out,
                                             const FermionField5d& in,
                                             const GaugeField& gf)
// set_left_expanded_gauge_field(gf, gf_);
// in.geo() = geo_resize(geo, 1);
// in.multiplicity = ls;
// refresh_expanded_1(in);
{
  TIMER("multiply_wilson_ddag_e_o_no_comm(5d,5d,gf)");
  qassert(&out != &in);
  qassert(is_matching_geo(gf.geo(), in.geo()));
  qassert(in.geo().eo == 1 or in.geo().eo == 2);
  Geometry geo = geo_resize(in.geo());
  geo.eo = 3 - in.geo().eo;
  if (is_initialized(out) and out.geo().eo == 3 - geo.eo) {
    out.geo().eo = geo.eo;
  }
  out.init(geo, in.multiplicity);
  set_zero(out);
  const int ls = in.multiplicity;
  qassert(out.multiplicity == ls);
  qassert(is_matching_geo(out.geo(), in.geo()));
  qassert(out.geo().eo != in.geo().eo);
  qassert(out.geo().eo == 1 or out.geo().eo == 2);
  const array<SpinMatrix, 4>& gammas =
      SpinMatrixConstants::get_cps_gammas();
  const SpinMatrix& unit = SpinMatrixConstants::get_unit();
  array<SpinMatrix, 4> p_mu_p;
  array<SpinMatrix, 4> p_mu_m;
  for (int mu = 0; mu < 4; ++mu) {
    p_mu_p[mu] = (ComplexD)0.5 * (unit + gammas[mu]);
    p_mu_m[mu] = (ComplexD)0.5 * (unit - gammas[mu]);
  }
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<WilsonVector> v = out.get_elems(xl);
    for (int mu = 0; mu < 4; ++mu) {
      const Coordinate xl_p = coordinate_shifts(xl, mu);
      const Coordinate xl_m = coordinate_shifts(xl, -mu - 1);
      const ColorMatrix u_p = gf.get_elem(xl, mu);
      const ColorMatrix u_m = matrix_adjoint(gf.get_elem(xl_m, mu));
      const Vector<WilsonVector> iv_p = in.get_elems_const(xl_p);
      const Vector<WilsonVector> iv_m = in.get_elems_const(xl_m);
      for (int m = 0; m < ls; ++m) {
        v[m] -= u_p * (p_mu_p[mu] * iv_p[m]);
        v[m] -= u_m * (p_mu_m[mu] * iv_m[m]);
      }
    }
  }
}

inline void multiply_m_e_o(FermionField5d& out, const FermionField5d& in,
                           const GaugeField& gf, const FermionAction& fa)
// out can be the same object as in
// works for _o_e as well
{
  TIMER("multiply_m_e_o(5d,5d,gf,fa)");
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  const SpinMatrix& unit = SpinMatrixConstants::get_unit();
  const SpinMatrix p_p = (ComplexD)0.5 * (unit + gamma5);
  const SpinMatrix p_m = (ComplexD)0.5 * (unit - gamma5);
  const int in_geo_eo = in.geo().eo;
  FermionField5d in1;
  in1.init(geo_resize(in.geo(), 1), in.multiplicity);
  const Geometry& geo = in.geo();
  std::vector<ComplexD> beo(fa.ls), ceo(fa.ls);
  for (int m = 0; m < fa.ls; ++m) {
    beo[m] = fa.bs[m];
    ceo[m] = -fa.cs[m];
  }
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<WilsonVector> iv = in.get_elems_const(xl);
    Vector<WilsonVector> v = in1.get_elems(xl);
    for (int m = 0; m < fa.ls; ++m) {
      v[m] = (ComplexD)beo[m] * iv[m];
      const WilsonVector tmp =
          (p_m * (m < fa.ls - 1
                      ? iv[m + 1]
                      : (WilsonVector)((ComplexD)(-fa.mass) * iv[0]))) +
          (p_p * (m > 0
                      ? iv[m - 1]
                      : (WilsonVector)((ComplexD)(-fa.mass) * iv[fa.ls - 1])));
      v[m] -= (ComplexD)ceo[m] * tmp;
    }
  }
  refresh_expanded_1(in1);
  multiply_wilson_d_e_o_no_comm(out, in1, gf);
  qassert(is_matching_geo(out.geo(), in.geo()));
  qassert(out.geo().eo != in_geo_eo);
  qassert(in.geo().eo == 1 or in.geo().eo == 2);
  qassert(out.geo().eo == 1 or out.geo().eo == 2);
}

inline void multiply_mdag_e_o(FermionField5d& out, const FermionField5d& in,
                              const GaugeField& gf, const FermionAction& fa)
// out can be the same object as in
// works for _o_e as well
{
  TIMER("multiply_mdag_e_o(5d,5d,gf,fa)");
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  const SpinMatrix& unit = SpinMatrixConstants::get_unit();
  const SpinMatrix p_p = (ComplexD)0.5 * (unit + gamma5);
  const SpinMatrix p_m = (ComplexD)0.5 * (unit - gamma5);
  const int in_geo_eo = in.geo().eo;
  qassert(is_matching_geo(gf.geo(), in.geo()));
  qassert(in.geo().eo == 1 or in.geo().eo == 2);
  Geometry geo = geo_resize(in.geo());
  geo.eo = 3 - in.geo().eo;
  std::vector<ComplexD> beo(fa.ls), ceo(fa.ls);
  for (int m = 0; m < fa.ls; ++m) {
    beo[m] = qconj(fa.bs[m]);
    ceo[m] = qconj(-fa.cs[m]);
  }
  FermionField5d in1;
  in1.init(geo_resize(in.geo(), 1), in.multiplicity);
  in1 = in;
  refresh_expanded_1(in1);
  FermionField5d out1;
  multiply_wilson_ddag_e_o_no_comm(out1, in1, gf);
  in1.init();
  if (is_initialized(out) and out.geo().eo == 3 - geo.eo) {
    out.geo().eo = geo.eo;
  }
  out.init(geo, in.multiplicity);
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<WilsonVector> iv = out1.get_elems_const(xl);
    Vector<WilsonVector> v = out.get_elems(xl);
    for (int m = 0; m < fa.ls; ++m) {
      v[m] = (ComplexD)beo[m] * iv[m];
      const WilsonVector tmp =
          (p_p * (m < fa.ls - 1
                      ? (WilsonVector)((ComplexD)ceo[m + 1] * iv[m + 1])
                      : (WilsonVector)(
                            (ComplexD)(-qconj((ComplexD)fa.mass) * ceo[0]) *
                            iv[0]))) +
          (p_m * (m > 0
                      ? (WilsonVector)((ComplexD)ceo[m - 1] * iv[m - 1])
                      : (WilsonVector)((ComplexD)(-qconj((ComplexD)fa.mass) *
                                                  ceo[fa.ls - 1]) *
                                       iv[fa.ls - 1])));
      v[m] -= tmp;
    }
  }
  qassert(is_matching_geo(out.geo(), in.geo()));
  qassert(out.geo().eo != in_geo_eo);
  qassert(in.geo().eo == 1 or in.geo().eo == 2);
  qassert(out.geo().eo == 1 or out.geo().eo == 2);
}

inline void multiply_m(FermionField5d& out, const FermionField5d& in,
                       const GaugeField& gf, const FermionAction& fa)
// out can be the same object as in
{
  TIMER("multiply_m");
  FermionField5d in_e, in_o;
  get_half_fermion(in_e, in, 2);
  get_half_fermion(in_o, in, 1);
  FermionField5d tmp_e, tmp_o;
  FermionField5d out_e, out_o;
  multiply_m_e_o(tmp_e, in_o, gf, fa);
  multiply_m_e_e(out_e, in_e, fa);
  multiply_m_e_e(tmp_o, in_o, fa);
  multiply_m_e_o(out_o, in_e, gf, fa);
  out_e += tmp_e;
  out_o += tmp_o;
  set_half_fermion(out, out_e, 2);
  set_half_fermion(out, out_o, 1);
}

inline void multiply_mdag(FermionField5d& out, const FermionField5d& in,
                          const GaugeField& gf, const FermionAction& fa)
// out can be the same object as in
{
  TIMER("multiply_mdag");
  FermionField5d in_e, in_o;
  get_half_fermion(in_e, in, 2);
  get_half_fermion(in_o, in, 1);
  FermionField5d tmp_e, tmp_o;
  FermionField5d out_e, out_o;
  multiply_mdag_e_o(tmp_e, in_o, gf, fa);
  multiply_mdag_e_e(out_e, in_e, fa);
  multiply_mdag_e_e(tmp_o, in_o, fa);
  multiply_mdag_e_o(out_o, in_e, gf, fa);
  out_e += tmp_e;
  out_o += tmp_o;
  set_half_fermion(out, out_e, 2);
  set_half_fermion(out, out_o, 1);
}

inline void multiply_mpc_sym2(FermionField5d& out, const FermionField5d& in,
                              const GaugeField& gf, const FermionAction& fa)
// out can be the same object as in
// odd <- odd (works for even <- even as well)
{
  TIMER("multiply_mpc_sym2");
  FermionField5d tmp;
  multiply_m_e_e_inv(tmp, in, fa);
  multiply_m_e_o(tmp, tmp, gf, fa);
  multiply_m_e_e_inv(tmp, tmp, fa);
  multiply_m_e_o(tmp, tmp, gf, fa);
  if (is_initialized(out) and out.geo().eo == 3 - in.geo().eo) {
    out.geo().eo = in.geo().eo;
  }
  out.init(geo_resize(in.geo()), in.multiplicity);
  out = in;
  out -= tmp;
}

inline void multiply_mpcdag_sym2(FermionField5d& out, const FermionField5d& in,
                                 const GaugeField& gf, const FermionAction& fa)
// out can be the same object as in
// odd <- odd (works for even <- even as well)
{
  TIMER("multiply_mpcdag_sym2");
  FermionField5d tmp;
  multiply_mdag_e_o(tmp, in, gf, fa);
  multiply_mdag_e_e_inv(tmp, tmp, fa);
  multiply_mdag_e_o(tmp, tmp, gf, fa);
  multiply_mdag_e_e_inv(tmp, tmp, fa);
  if (is_initialized(out) and out.geo().eo == 3 - in.geo().eo) {
    out.geo().eo = in.geo().eo;
  }
  out.init(geo_resize(in.geo()), in.multiplicity);
  out = in;
  out -= tmp;
}

inline void multiply_hermop_sym2(FermionField5d& out, const FermionField5d& in,
                                 const GaugeField& gf, const FermionAction& fa)
// out can be the same object as in
// odd <- odd (works for even <- even as well)
{
  TIMER_FLOPS("multiply_hermop_sym2");
  multiply_mpc_sym2(out, in, gf, fa);
  multiply_mpcdag_sym2(out, out, gf, fa);
  timer.flops += 5500 * fa.ls * gf.geo().local_volume();
}

inline void multiply_m_e_e(FermionField5d& out, const FermionField5d& in,
                           const InverterDomainWall& inv)
// out can be the same object as in
{
  multiply_m_e_e(out, in, inv.fa);
}

inline void multiply_mdag_e_e(FermionField5d& out, const FermionField5d& in,
                              const InverterDomainWall& inv)
// out can be the same object as in
{
  multiply_mdag_e_e(out, in, inv.fa);
}

inline void multiply_m_e_e_inv(FermionField5d& out, const FermionField5d& in,
                               const InverterDomainWall& inv)
// out can be the same object as in
{
  multiply_m_e_e_inv(out, in, inv.fa);
}

inline void multiply_mdag_e_e_inv(FermionField5d& out, const FermionField5d& in,
                                  const InverterDomainWall& inv)
// out can be the same object as in
{
  multiply_mdag_e_e_inv(out, in, inv.fa);
}

inline void multiply_m_e_o(FermionField5d& out, const FermionField5d& in,
                           const InverterDomainWall& inv)
// out can be the same object as in
{
  multiply_m_e_o(out, in, inv.gf, inv.fa);
}

inline void multiply_mdag_e_o(FermionField5d& out, const FermionField5d& in,
                              const InverterDomainWall& inv)
// out can be the same object as in
{
  multiply_mdag_e_o(out, in, inv.gf, inv.fa);
}

inline void multiply_m_eo_eo(FermionField5d& out, const FermionField5d& in,
                             const InverterDomainWall& inv, const int eo_out,
                             const int eo_in)
// out can be the same object as in
// out need to be initialized with correct geo and eo
{
  TIMER("multiply_m_eo_eo");
  Geometry geo = geo_resize(in.geo());
  geo.eo = eo_out;
  out.init(geo, in.multiplicity);
  qassert(is_matching_geo(out.geo(), in.geo()));
  qassert(out.geo().eo == eo_out);
  qassert(in.geo().eo == eo_in);
  if (eo_out == eo_in) {
    multiply_m_e_e(out, in, inv.fa);
  } else {
    multiply_m_e_o(out, in, inv.gf, inv.fa);
  }
}

inline void multiply_mdag_eo_eo(FermionField5d& out, const FermionField5d& in,
                                const InverterDomainWall& inv, const int eo_out,
                                const int eo_in)
// out can be the same object as in
// out need to be initialized with correct geo and eo
{
  TIMER("multiply_mdag_eo_eo");
  Geometry geo = geo_resize(in.geo());
  geo.eo = eo_out;
  out.init(geo, in.multiplicity);
  qassert(is_matching_geo(out.geo(), in.geo()));
  qassert(out.geo().eo == eo_out);
  qassert(in.geo().eo == eo_in);
  if (eo_out == eo_in) {
    multiply_mdag_e_e(out, in, inv.fa);
  } else {
    multiply_mdag_e_o(out, in, inv.gf, inv.fa);
  }
}

inline void multiply_m(FermionField5d& out, const FermionField5d& in,
                       const InverterDomainWall& inv)
// out can be the same object as in
{
  multiply_m(out, in, inv.gf, inv.fa);
}

inline void multiply_mdag(FermionField5d& out, const FermionField5d& in,
                          const InverterDomainWall& inv)
// out can be the same object as in
{
  multiply_mdag(out, in, inv.gf, inv.fa);
}

inline void multiply_mpc_sym2(FermionField5d& out, const FermionField5d& in,
                              const InverterDomainWall& inv)
// out can be the same object as in
// odd <- odd (works for even <- even as well)
{
  multiply_mpc_sym2(out, in, inv.gf, inv.fa);
}

inline void multiply_mpcdag_sym2(FermionField5d& out, const FermionField5d& in,
                                 const InverterDomainWall& inv)
// out can be the same object as in
// odd <- odd (works for even <- even as well)
{
  multiply_mpcdag_sym2(out, in, inv.gf, inv.fa);
}

inline void multiply_hermop_sym2(FermionField5d& out, const FermionField5d& in,
                                 const InverterDomainWall& inv)
// out can be the same object as in
// odd <- odd (works for even <- even as well)
{
  multiply_hermop_sym2(out, in, inv.gf, inv.fa);
}

inline void multiply_m_with_prec_sym2(FermionField5d& out,
                                      const FermionField5d& in,
                                      const InverterDomainWall& inv)
// out can be the same object as in
{
  TIMER("multiply_m_with_prec_sym2");
  FermionField5d in_e, in_o;
  FermionField5d out_e, out_o;
  get_half_fermion(in_e, in, 2);
  get_half_fermion(in_o, in, 1);
  //
  multiply_m_e_o(out_e, in_o, inv);
  multiply_m_e_e_inv(out_e, out_e, inv);
  out_e += in_e;
  multiply_m_e_e(out_o, in_o, inv);
  //
  multiply_m_e_e(out_e, out_e, inv);
  multiply_mpc_sym2(out_o, out_o, inv);
  //
  FermionField5d tmp;
  multiply_m_e_e_inv(tmp, out_e, inv);
  multiply_m_e_o(tmp, tmp, inv);
  out_o += tmp;
  //
  set_half_fermion(out, out_e, 2);
  set_half_fermion(out, out_o, 1);
}

inline ComplexD dot_product(const FermionField5d& ff1, const FermionField5d& ff2)
// return ff1^dag * ff2
{
  TIMER("dot_product");
  qassert(is_matching_geo(ff1.geo(), ff2.geo()));
  qassert(ff1.geo().eo == ff2.geo().eo);
  const Geometry& geo = ff1.geo();
  ComplexD sum = 0.0;
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<WilsonVector> v1 = ff1.get_elems_const(xl);
    const Vector<WilsonVector> v2 = ff2.get_elems_const(xl);
    const Vector<ComplexD> cv1((const ComplexD*)v1.data(),
                               v1.data_size() / sizeof(ComplexD));
    const Vector<ComplexD> cv2((const ComplexD*)v2.data(),
                               v2.data_size() / sizeof(ComplexD));
    qassert(cv1.size() == cv2.size());
    for (int k = 0; k < cv1.size(); ++k) {
      sum += qconj(cv1[k]) * cv2[k];
    }
  }
  glb_sum(sum);
  return sum;
}

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
  displayln_info(
      fname +
      ssprintf(
          ": final max_num_iter=%4ld sqrt(qnorm_r/qnorm_in)=%.3E stop_rsd=%.3E",
          max_num_iter, sqrt(qnorm_r / qnorm_in), stop_rsd));
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
    int cycle;
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

inline Long cg_with_herm_sym_2(FermionField5d& sol, const FermionField5d& src,
                               const InverterDomainWall& inv,
                               const double stop_rsd = 1e-8,
                               const Long max_num_iter = 50000)
{
  TIMER_VERBOSE_FLOPS("cg_with_herm_sym_2(5d,5d,inv)");
  qassert(&sol != &src);
  const Long iter =
      cg_with_f(sol, src, inv, multiply_hermop_sym2, stop_rsd, max_num_iter);
  timer.flops += 5500 * iter * inv.fa.ls * inv.geo().local_volume();
  return iter;
}

inline Long invert(FermionField5d& out, const FermionField5d& in,
                   const InverterDomainWall& inv)
{
  qassert(&out != &in);
  return invert_with_cg(out, in, inv, cg_with_herm_sym_2);
}

inline Long invert(FermionField4d& out, const FermionField4d& in,
                   const InverterDomainWall& inv)
{
  qassert(&out != &in);
  return invert_dwf(out, in, inv);
}

inline double find_max_eigen_value_hermop_sym2(const InverterDomainWall& inv,
                                               const RngState& rs,
                                               const Long max_iter = 100)
{
  TIMER_VERBOSE("find_max_eigen_value_hermop_sym2");
  Geometry geo = geo_resize(inv.geo());
  geo.eo = 1;
  FermionField5d ff;
  ff.init(geo, inv.fa.ls);
  set_u_rand(ff, rs);
  ff *= 1.0 / sqrt(qnorm(ff));
  double sqrt_qnorm_ratio = 1.0;
  for (Long i = 0; i < max_iter; ++i) {
    multiply_hermop_sym2(ff, ff, inv);
    sqrt_qnorm_ratio = sqrt(qnorm(ff));
    displayln_info(fname + ssprintf(": %5d: sqrt_qnorm_ratio =%24.17E."
                                    " (max eigen value of hermop_sym2)",
                                    i + 1, sqrt_qnorm_ratio));
    ff *= 1.0 / sqrt_qnorm_ratio;
  }
  return sqrt_qnorm_ratio;
}

}  // namespace qlat
