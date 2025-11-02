#define QLAT_INSTANTIATE_SMEAR

#include <qlat/compressed-eigen-io.h>
#include <qlat/dslash.h>

namespace qlat
{  //

Long load_low_modes(LowModes& lm, const std::string& path)
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

Long load_or_compute_low_modes(LowModes& lm, const std::string& path,
                               const GaugeField& gf, const FermionAction& fa,
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

void load_low_modes_delay(LowModes& lm, const std::string& path)
{
  TIMER_VERBOSE("load_low_modes_delay");
  lm.initialized = true;
  lm.lmi.init();
  lm.lmi.initialized = true;
  lm.lmi.path = path;
}

void load_or_compute_low_modes_delay(LowModes& lm, const std::string& path,
                                     const GaugeField& gf,
                                     const FermionAction& fa, const LancArg& la)
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

Long force_low_modes(LowModes& lm)
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

void set_u_rand(LowModes& lm, const RngState& rs)
{
  TIMER_VERBOSE("set_u_rand(lm,rs)");
  Qassert(lm.initialized);
  force_low_modes(lm);
  set_u_rand(get_data(lm.eigen_values), rs.split("eigen_values"));
  set_u_rand(lm.cesb, rs.split("cesb"));
  set_u_rand(lm.cesc, rs.split("cesc"));
}

Long save_low_modes_decompress(LowModes& lm, const std::string& path)
{
  TIMER_VERBOSE("save_low_modes_decompress");
  Qassert(lm.initialized);
  force_low_modes(lm);
  std::vector<BlockedHalfVector> bhvs;
  decompress_eigen_system(bhvs, lm.cesb, lm.cesc);
  std::vector<HalfVector> hvs;
  convert_half_vectors(hvs, bhvs);
  const Int id_node = get_id_node();
  const Int num_node = get_num_node();
  const Int idx = id_node;
  const Int idx_size = num_node;
  const Int dir_idx = compute_dist_file_dir_id(idx, idx_size);
  qmkdir(path);
  qmkdir(path + ssprintf("/%02d", dir_idx));
  const std::string fn = path + ssprintf("/%02d/%010d", dir_idx, idx);
  Long total_bytes = 0;
  const Int n_cycle = std::max(1, num_node / dist_write_par_limit());
  std::vector<crc32_t> crcs(num_node, 0);
  for (Int i = 0; i < n_cycle; i++) {
    Long bytes = 0;
    if (id_node % n_cycle == i) {
      Qassert(hvs.size() >= 1);
      bytes = hvs.size() * get_data(hvs[0]).data_size();
      crcs[id_node] = save_half_vectors(hvs, fn, false, true);
    }
    glb_sum(bytes);
    total_bytes += bytes;
    displayln_info(
        ssprintf("qlat::%s: cycle / n_cycle = %4d / %4d ; total_bytes = %15ld",
                 fname.c_str(), i + 1, n_cycle, total_bytes));
  }
  glb_sum(get_data_char(crcs));
  const crc32_t crc = dist_crc32(crcs);
  if (get_id_node() == 0) {
    const std::string fn = path + "/checksums.txt";
    QFile fp = qfopen(fn, "w");
    Qassert(not fp.null());
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
    Qassert(not fp.null());
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

void deflate(HalfVector& hv_out, const HalfVector& hv_in, LowModes& lm)
// hv_out can be the same as hv_in
{
  force_low_modes(lm);
  if (not lm.initialized) {
    hv_out.init(geo_resize(hv_in.geo()), hv_in.multiplicity);
    set_zero(hv_out);
    return;
  }
  TIMER("deflate(hv,hv,lm)");
  const Int ls = lm.cesi.ls;
  const Long block_size = lm.cesb.block_vol_eo * ls * HalfVector::c_size;
  const Long n_basis = lm.cesb.n_basis;
  const Long n_vec = lm.cesc.n_vec;
  Qassert(n_vec == (Long)lm.eigen_values.size());
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
      for (Int j = 0; j < n_basis; ++j) {
        ComplexD& vc_j = vc.p[j];
        const Vector<ComplexF> vbs_j(vbs.p + j * block_size, block_size);
        for (Long k = 0; k < block_size; ++k) {
          const ComplexF& vb_k = vb.p[k];
          vc_j += qconj(vbs_j.p[k]) * vb_k;
        }
      }
      // compute inner products
      Vector<ComplexD> vp = phv.get_elems(index);
      for (Int i = 0; i < n_vec; ++i) {
        ComplexD& vp_i = vp.p[i];
        const Vector<ComplexF> vcs_i(vcs.p + i * n_basis, n_basis);
        for (Int j = 0; j < n_basis; ++j) {
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
      for (Int i = 0; i < n_vec; ++i) {
        phv_sum[i] += vp[i];
      }
    }
    phv.init();
    glb_sum(phv_sum);
    // scale by eigen values
#pragma omp parallel for
    for (Int i = 0; i < n_vec; ++i) {
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
      for (Int i = 0; i < n_vec; ++i) {
        const ComplexD& phv_sum_i = phv_sum[i];
        const Vector<ComplexF> vcs_i(vcs.p + i * n_basis, n_basis);
        for (Int j = 0; j < n_basis; ++j) {
          ComplexD& vc_j = vc.p[j];
          vc_j += (ComplexD)(vcs_i.p[j]) * phv_sum_i;
        }
      }
      // project to fine grid
      Vector<ComplexF> vb = bhv.get_elems(index);
      for (Int j = 0; j < n_basis; ++j) {
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

void deflate(FermionField5d& out, const FermionField5d& in, LowModes& lm)
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
  Qassert(geo.eo == 1 or geo.eo == 2);
  const Int ls = in.multiplicity;
  Qassert(ls == lm.cesi.ls);
  HalfVector hv;
  init_half_vector(hv, geo, ls);
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    (void) xl;
    qassert((xl[0] + xl[1] + xl[2] + xl[3]) % 2 == 2 - geo.eo);
    Vector<ComplexF> vhv = hv.get_elems(index);
    const Vector<WilsonVector> vin = in.get_elems_const(index);
    qassert(vhv.size() ==
            vin.size() * (Long)sizeof(WilsonVector) / (Long)sizeof(ComplexD));
    const Vector<ComplexD> vff((const ComplexD*)vin.data(), vhv.size());
    for (Int m = 0; m < vhv.size(); ++m) {
      vhv[m] = vff[m];
    }
  }
  deflate(hv, hv, lm);
  if (is_initialized(out) and out.geo().eo == 3 - geo.eo) {
    out.geo().eo = geo.eo;
  }
  out.init(geo, in.multiplicity);
  Qassert(out.geo().eo == geo.eo);
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    (void) xl;
    qassert((xl[0] + xl[1] + xl[2] + xl[3]) % 2 == 2 - geo.eo);
    const Vector<ComplexF> vhv = hv.get_elems(index);
    Vector<WilsonVector> vout = out.get_elems(index);
    qassert(vhv.size() ==
            vout.size() * (Long)sizeof(WilsonVector) / (Long)sizeof(ComplexD));
    Vector<ComplexD> vff((ComplexD*)vout.data(), vhv.size());
    for (Int m = 0; m < vhv.size(); ++m) {
      vff[m] = vhv[m];
    }
  }
}

void benchmark_deflate(const Geometry& geo, const Int ls,
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
  SYNC_NODE();
  {
    TIMER_VERBOSE("benchmark_deflate-deflate");
    for (Int i = 0; i < 4; ++i) {
      set_zero(out);
      deflate(out, in, lm);
    }
  }
  SYNC_NODE();
}

void multiply_m_dwf_no_comm(FermionField5d& out, const FermionField5d& in,
                            const InverterDomainWall& inv)
{
  TIMER("multiply_m_dwf_no_comm(5d,5d,Inv)");
  Qassert(&out != &in);
  const Geometry geo = geo_resize(in.geo());
  Qassert(is_matching_geo(inv.geo(), geo));
  out.init(geo, in.multiplicity);
  set_zero(out);
  const FermionAction& fa = inv.fa;
  const GaugeField& gf = inv.gf;
  Qassert(in.multiplicity == fa.ls);
  Qassert(out.multiplicity == fa.ls);
  Qassert(fa.mobius_scale == 1.0);
  Qassert((int)fa.bs.size() == fa.ls);
  Qassert((int)fa.cs.size() == fa.ls);
  const array<SpinMatrix, 4>& gammas = SpinMatrixConstants::get_cps_gammas();
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  const SpinMatrix& unit = SpinMatrixConstants::get_unit();
  const SpinMatrix p_p = (ComplexD)0.5 * (unit + gamma5);
  const SpinMatrix p_m = (ComplexD)0.5 * (unit - gamma5);
  array<SpinMatrix, 4> p_mu_p;
  array<SpinMatrix, 4> p_mu_m;
  for (Int mu = 0; mu < 4; ++mu) {
    p_mu_p[mu] = (ComplexD)0.5 * (unit + gammas[mu]);
    p_mu_m[mu] = (ComplexD)0.5 * (unit - gammas[mu]);
  }
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<WilsonVector> v = out.get_elems(xl);
    {
      const Vector<WilsonVector> iv = in.get_elems_const(xl);
      for (Int m = 0; m < fa.ls; ++m) {
        v[m] = (ComplexD)(5.0 - fa.m5) * iv[m];
        v[m] -= p_m * (m < fa.ls - 1
                           ? iv[m + 1]
                           : (WilsonVector)((ComplexD)(-fa.mass) * iv[0]));
        v[m] -= p_p *
                (m > 0 ? iv[m - 1]
                       : (WilsonVector)((ComplexD)(-fa.mass) * iv[fa.ls - 1]));
      }
    }
    for (Int mu = 0; mu < 4; ++mu) {
      const Coordinate xl_p = coordinate_shifts(xl, mu);
      const Coordinate xl_m = coordinate_shifts(xl, -mu - 1);
      const ColorMatrix u_p = gf.get_elem(xl, mu);
      const ColorMatrix u_m = matrix_adjoint(gf.get_elem(xl_m, mu));
      const Vector<WilsonVector> iv_p = in.get_elems_const(xl_p);
      const Vector<WilsonVector> iv_m = in.get_elems_const(xl_m);
      for (Int m = 0; m < fa.ls; ++m) {
        v[m] -= u_p * (p_mu_m[mu] * iv_p[m]);
        v[m] -= u_m * (p_mu_p[mu] * iv_m[m]);
      }
    }
  }
}

void multiply_m_dwf(FermionField5d& out, const FermionField5d& in,
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

void multiply_wilson_d_no_comm(FermionField5d& out, const FermionField5d& in,
                               const GaugeField& gf, const double mass)
// set_left_expanded_gauge_field(gf, gf_);
// in.geo() = geo_resize(geo, 1);
// in.multiplicity == ls;
// refresh_expanded_1(in);
{
  TIMER("multiply_wilson_d_no_comm(5d,5d,gf,mass)");
  Qassert(&out != &in);
  const Geometry geo = geo_resize(in.geo());
  Qassert(is_matching_geo(gf.geo(), geo));
  out.init(geo, in.multiplicity);
  set_zero(out);
  const Int ls = in.multiplicity;
  Qassert(out.multiplicity == ls);
  const array<SpinMatrix, 4>& gammas = SpinMatrixConstants::get_cps_gammas();
  const SpinMatrix& unit = SpinMatrixConstants::get_unit();
  array<SpinMatrix, 4> p_mu_p;
  array<SpinMatrix, 4> p_mu_m;
  for (Int mu = 0; mu < 4; ++mu) {
    p_mu_p[mu] = (ComplexD)0.5 * (unit + gammas[mu]);
    p_mu_m[mu] = (ComplexD)0.5 * (unit - gammas[mu]);
  }
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<WilsonVector> v = out.get_elems(xl);
    {
      const Vector<WilsonVector> iv = in.get_elems_const(xl);
      for (Int m = 0; m < ls; ++m) {
        v[m] = (ComplexD)(4.0 + mass) * iv[m];
      }
    }
    for (Int mu = 0; mu < 4; ++mu) {
      const Coordinate xl_p = coordinate_shifts(xl, mu);
      const Coordinate xl_m = coordinate_shifts(xl, -mu - 1);
      const ColorMatrix u_p = gf.get_elem(xl, mu);
      const ColorMatrix u_m = matrix_adjoint(gf.get_elem(xl_m, mu));
      const Vector<WilsonVector> iv_p = in.get_elems_const(xl_p);
      const Vector<WilsonVector> iv_m = in.get_elems_const(xl_m);
      for (Int m = 0; m < ls; ++m) {
        v[m] -= u_p * (p_mu_m[mu] * iv_p[m]);
        v[m] -= u_m * (p_mu_p[mu] * iv_m[m]);
      }
    }
  }
}

void multiply_d_minus(FermionField5d& out, const FermionField5d& in,
                      const GaugeField& gf, const FermionAction& fa)
// out can be the same object as in
{
  TIMER("multiply_d_minus(5d,5d,gf,fa)");
  const Geometry geo = geo_resize(in.geo());
  Qassert(is_matching_geo(gf.geo(), in.geo()));
  Qassert(in.multiplicity == fa.ls);
  Qassert((int)fa.bs.size() == fa.ls);
  Qassert((int)fa.cs.size() == fa.ls);
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
    for (Int m = 0; m < fa.ls; ++m) {
      const ComplexD& c = fa.cs[m];
      v1[m] = (ComplexD)(-c) * iv[m];
      v[m] = iv[m];
    }
  }
  refresh_expanded_1(in1);
  out.init(geo, in.multiplicity);
  Qassert(out.multiplicity == fa.ls);
  Qassert(is_matching_geo(gf.geo(), out.geo()));
  multiply_wilson_d_no_comm(out, in1, gf, -fa.m5);
  out += out1;
}

void multiply_d_minus(FermionField5d& out, const FermionField5d& in,
                      const InverterDomainWall& inv)
// out can be the same object as in
{
  TIMER("multiply_d_minus(5d,5d,Inv)");
  multiply_d_minus(out, in, inv.gf, inv.fa);
}

void multiply_m_full(FermionField5d& out, const FermionField5d& in,
                     const InverterDomainWall& inv)
{
  TIMER("multiply_m_full(5d,5d,Inv)");
  Qassert(&out != &in);
  const Geometry geo = geo_resize(in.geo());
  out.init(geo, in.multiplicity);
  set_zero(out);
  const FermionAction& fa = inv.fa;
  Qassert(is_matching_geo(inv.geo(), in.geo()));
  Qassert(is_matching_geo(inv.geo(), out.geo()));
  Qassert(in.multiplicity == fa.ls);
  Qassert(in.multiplicity == fa.ls);
  Qassert(out.multiplicity == fa.ls);
  Qassert((int)fa.bs.size() == fa.ls);
  Qassert((int)fa.cs.size() == fa.ls);
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
    for (Int m = 0; m < fa.ls; ++m) {
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

void get_half_fermion(FermionField5d& half, const FermionField5d& ff,
                      const Int eo)
// 2:even 1:odd
{
  TIMER("get_half_fermion");
  Qassert(&half != &ff);
  Geometry geoh = geo_resize(ff.geo());
  geoh.eo = eo;
  half.init(geoh, ff.multiplicity);
  Qassert(half.geo().eo == eo);
  Qassert(ff.geo().eo == 0);
  Qassert(is_matching_geo(ff.geo(), half.geo()));
#pragma omp parallel for
  for (Long index = 0; index < geoh.local_volume(); ++index) {
    const Coordinate xl = geoh.coordinate_from_index(index);
    assign(half.get_elems(xl), ff.get_elems_const(xl));
  }
}

void set_half_fermion(FermionField5d& ff, const FermionField5d& half,
                      const Int eo)
// 2:even 1:odd
{
  TIMER("set_half_fermion");
  Qassert(&half != &ff);
  const Geometry geoh = half.geo();
  Geometry geo = geo_resize(geoh);
  geo.eo = 0;
  ff.init(geo, half.multiplicity);
  Qassert(half.geo().eo == eo);
  Qassert(ff.geo().eo == 0);
  Qassert(is_matching_geo(ff.geo(), half.geo()));
#pragma omp parallel for
  for (Long index = 0; index < geoh.local_volume(); ++index) {
    const Coordinate xl = geoh.coordinate_from_index(index);
    assign(ff.get_elems(xl), half.get_elems_const(xl));
  }
}

void project_eo(FermionField5d& ff, const Int eo)
{
  TIMER("project_eo");
  Qassert(eo == 1 or eo == 2);
  FermionField5d half;
  get_half_fermion(half, ff, eo);
  set_zero(ff);
  set_half_fermion(ff, half, eo);
}

void multiply_m_e_e(FermionField5d& out, const FermionField5d& in,
                    const FermionAction& fa)
// out can be the same object as in
// works for _o_o as well
{
  TIMER("multiply_m_e_e");
  if (is_initialized(out) and out.geo().eo == 3 - in.geo().eo) {
    out.geo().eo = in.geo().eo;
  }
  out.init(geo_resize(in.geo()), in.multiplicity);
  Qassert(is_matching_geo(out.geo(), in.geo()));
  Qassert(in.geo().eo == 1 or in.geo().eo == 2);
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
  for (Int m = 0; m < fa.ls; ++m) {
    bee[m] = 1.0 + fa.bs[m] * (4.0 - fa.m5);
    cee[m] = 1.0 - fa.cs[m] * (4.0 - fa.m5);
  }
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<WilsonVector> iv = hin().get_elems_const(xl);
    Vector<WilsonVector> v = out.get_elems(xl);
    for (Int m = 0; m < fa.ls; ++m) {
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

void multiply_mdag_e_e(FermionField5d& out, const FermionField5d& in,
                       const FermionAction& fa)
// out can be the same object as in
// works for _o_o as well
{
  TIMER("multiply_mdag_e_e");
  if (is_initialized(out) and out.geo().eo == 3 - in.geo().eo) {
    out.geo().eo = in.geo().eo;
  }
  out.init(geo_resize(in.geo()), in.multiplicity);
  Qassert(is_matching_geo(out.geo(), in.geo()));
  Qassert(in.geo().eo == 1 or in.geo().eo == 2);
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
  for (Int m = 0; m < fa.ls; ++m) {
    bee[m] = qconj(1.0 + fa.bs[m] * (4.0 - fa.m5));
    cee[m] = qconj(1.0 - fa.cs[m] * (4.0 - fa.m5));
  }
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<WilsonVector> iv = hin().get_elems_const(xl);
    Vector<WilsonVector> v = out.get_elems(xl);
    for (Int m = 0; m < fa.ls; ++m) {
      v[m] = (ComplexD)bee[m] * iv[m];
      const WilsonVector tmp =
          (p_p * (m < fa.ls - 1
                      ? (WilsonVector)((ComplexD)cee[m + 1] * iv[m + 1])
                      : (WilsonVector)((ComplexD)(-qconj((ComplexD)fa.mass) *
                                                  cee[0]) *
                                       iv[0]))) +
          (p_m * (m > 0 ? (WilsonVector)((ComplexD)cee[m - 1] * iv[m - 1])
                        : (WilsonVector)((ComplexD)(-qconj((ComplexD)fa.mass) *
                                                    cee[fa.ls - 1]) *
                                         iv[fa.ls - 1])));
      v[m] -= tmp;
    }
  }
}

void multiply_m_e_e_inv(FermionField5d& out, const FermionField5d& in,
                        const FermionAction& fa)
// out can be the same object as in
// works for _o_o as well
{
  TIMER("multiply_m_e_e_inv");
  if (is_initialized(out) and out.geo().eo == 3 - in.geo().eo) {
    out.geo().eo = in.geo().eo;
  }
  out.init(geo_resize(in.geo()), in.multiplicity);
  Qassert(is_matching_geo(out.geo(), in.geo()));
  Qassert(out.geo().eo == in.geo().eo);
  Qassert(in.geo().eo == 1 or in.geo().eo == 2);
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  const SpinMatrix& unit = SpinMatrixConstants::get_unit();
  const SpinMatrix p_p = (ComplexD)0.5 * (unit + gamma5);
  const SpinMatrix p_m = (ComplexD)0.5 * (unit - gamma5);
  const Geometry& geo = out.geo();
  std::vector<ComplexD> bee(fa.ls), cee(fa.ls);
  for (Int m = 0; m < fa.ls; ++m) {
    bee[m] = 1.0 + fa.bs[m] * (4.0 - fa.m5);
    cee[m] = 1.0 - fa.cs[m] * (4.0 - fa.m5);
  }
  std::vector<ComplexD> lee(fa.ls - 1), leem(fa.ls - 1);
  for (Int m = 0; m < fa.ls - 1; ++m) {
    lee[m] = -cee[m + 1] / bee[m];
    leem[m] = m == 0 ? fa.mass * cee[fa.ls - 1] / bee[0]
                     : leem[m - 1] * cee[m - 1] / bee[m];
  }
  std::vector<ComplexD> dee(fa.ls, 0.0);
  dee[fa.ls - 1] = fa.mass * cee[fa.ls - 1];
  for (Int m = 0; m < fa.ls - 1; ++m) {
    dee[fa.ls - 1] *= cee[m] / bee[m];
  }
  for (Int m = 0; m < fa.ls; ++m) {
    dee[m] += bee[m];
  }
  std::vector<ComplexD> uee(fa.ls - 1), ueem(fa.ls - 1);
  for (Int m = 0; m < fa.ls - 1; ++m) {
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
    for (Int m = 0; m < fa.ls - 1; ++m) {
      tmp += (ComplexD)(-leem[m]) * v[m];
    }
    v[fa.ls - 1] += p_m * tmp;
    // {L'_{ee}}^{-1}
    for (Int m = 1; m < fa.ls; ++m) {
      v[m] += (ComplexD)(-lee[m - 1]) * (p_p * v[m - 1]);
    }
    // {D_{ee}}^{-1}
    for (Int m = 0; m < fa.ls; ++m) {
      v[m] *= (ComplexD)(1.0 / dee[m]);
    }
    // {U^'_{ee}}^{-1}
    for (Int m = fa.ls - 2; m >= 0; --m) {
      v[m] += (ComplexD)(-uee[m]) * (p_m * v[m + 1]);
    }
    // {U^m_{ee}}^{-1}
    for (Int m = 0; m < fa.ls - 1; ++m) {
      v[m] += (ComplexD)(-ueem[m]) * (p_p * v[fa.ls - 1]);
    }
  }
}

void multiply_mdag_e_e_inv(FermionField5d& out, const FermionField5d& in,
                           const FermionAction& fa)
// out can be the same object as in
// works for _o_o as well
{
  TIMER("multiply_mdag_e_e_inv");
  if (is_initialized(out) and out.geo().eo == 3 - in.geo().eo) {
    out.geo().eo = in.geo().eo;
  }
  out.init(geo_resize(in.geo()), in.multiplicity);
  Qassert(is_matching_geo(out.geo(), in.geo()));
  Qassert(out.geo().eo == in.geo().eo);
  Qassert(in.geo().eo == 1 or in.geo().eo == 2);
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  const SpinMatrix& unit = SpinMatrixConstants::get_unit();
  const SpinMatrix p_p = (ComplexD)0.5 * (unit + gamma5);
  const SpinMatrix p_m = (ComplexD)0.5 * (unit - gamma5);
  const Geometry& geo = out.geo();
  std::vector<ComplexD> bee(fa.ls), cee(fa.ls);
  for (Int m = 0; m < fa.ls; ++m) {
    bee[m] = 1.0 + fa.bs[m] * (4.0 - fa.m5);
    cee[m] = 1.0 - fa.cs[m] * (4.0 - fa.m5);
  }
  std::vector<ComplexD> lee(fa.ls - 1), leem(fa.ls - 1);
  for (Int m = 0; m < fa.ls - 1; ++m) {
    lee[m] = -cee[m + 1] / bee[m];
    leem[m] = m == 0 ? fa.mass * cee[fa.ls - 1] / bee[0]
                     : leem[m - 1] * cee[m - 1] / bee[m];
  }
  std::vector<ComplexD> dee(fa.ls, 0.0);
  dee[fa.ls - 1] = fa.mass * cee[fa.ls - 1];
  for (Int m = 0; m < fa.ls - 1; ++m) {
    dee[fa.ls - 1] *= cee[m] / bee[m];
  }
  for (Int m = 0; m < fa.ls; ++m) {
    dee[m] += bee[m];
  }
  std::vector<ComplexD> uee(fa.ls - 1), ueem(fa.ls - 1);
  for (Int m = 0; m < fa.ls - 1; ++m) {
    uee[m] = -cee[m] / bee[m];
    ueem[m] =
        m == 0 ? fa.mass * cee[0] / bee[0] : ueem[m - 1] * cee[m] / bee[m];
  }
  for (Int m = 0; m < fa.ls; ++m) {
    bee[m] = qconj(bee[m]);
    cee[m] = qconj(cee[m]);
    dee[m] = qconj(dee[m]);
  }
  for (Int m = 0; m < fa.ls - 1; ++m) {
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
    for (Int m = 0; m < fa.ls - 1; ++m) {
      tmp += (ComplexD)(-ueem[m]) * v[m];
    }
    v[fa.ls - 1] += p_p * tmp;
    // {U^'_{ee}}^\dagger^{-1}
    for (Int m = 1; m < fa.ls; ++m) {
      v[m] += (ComplexD)(-uee[m - 1]) * (p_m * v[m - 1]);
    }
    // {D_{ee}}^\dagger^{-1}
    for (Int m = 0; m < fa.ls; ++m) {
      v[m] *= (ComplexD)(1.0 / dee[m]);
    }
    // {L'_{ee}}^\dagger^{-1}
    for (Int m = fa.ls - 2; m >= 0; --m) {
      v[m] += (ComplexD)(-lee[m]) * (p_p * v[m + 1]);
    }
    // {L^m_{ee}}^\dagger^{-1}
    for (Int m = 0; m < fa.ls - 1; ++m) {
      v[m] += (ComplexD)(-leem[m]) * (p_m * v[fa.ls - 1]);
    }
  }
}

void multiply_wilson_d_e_o_no_comm(FermionField5d& out,
                                   const FermionField5d& in,
                                   const GaugeField& gf)
// set_left_expanded_gauge_field(gf, gf_);
// in.geo() = geo_resize(geo, 1);
// in.multiplicity = ls;
// refresh_expanded_1(in);
{
  TIMER("multiply_wilson_d_e_o_no_comm(5d,5d,gf)");
  Qassert(&out != &in);
  Qassert(is_matching_geo(gf.geo(), in.geo()));
  Qassert(in.geo().eo == 1 or in.geo().eo == 2);
  Geometry geo = geo_resize(in.geo());
  geo.eo = 3 - in.geo().eo;
  if (is_initialized(out) and out.geo().eo == 3 - geo.eo) {
    out.geo().eo = geo.eo;
  }
  out.init(geo, in.multiplicity);
  set_zero(out);
  const Int ls = in.multiplicity;
  Qassert(out.multiplicity == ls);
  Qassert(is_matching_geo(out.geo(), in.geo()));
  Qassert(out.geo().eo != in.geo().eo);
  Qassert(out.geo().eo == 1 or out.geo().eo == 2);
  const array<SpinMatrix, 4>& gammas = SpinMatrixConstants::get_cps_gammas();
  const SpinMatrix& unit = SpinMatrixConstants::get_unit();
  array<SpinMatrix, 4> p_mu_p;
  array<SpinMatrix, 4> p_mu_m;
  for (Int mu = 0; mu < 4; ++mu) {
    p_mu_p[mu] = (ComplexD)0.5 * (unit + gammas[mu]);
    p_mu_m[mu] = (ComplexD)0.5 * (unit - gammas[mu]);
  }
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<WilsonVector> v = out.get_elems(xl);
    for (Int mu = 0; mu < 4; ++mu) {
      const Coordinate xl_p = coordinate_shifts(xl, mu);
      const Coordinate xl_m = coordinate_shifts(xl, -mu - 1);
      const ColorMatrix u_p = gf.get_elem(xl, mu);
      const ColorMatrix u_m = matrix_adjoint(gf.get_elem(xl_m, mu));
      const Vector<WilsonVector> iv_p = in.get_elems_const(xl_p);
      const Vector<WilsonVector> iv_m = in.get_elems_const(xl_m);
      for (Int m = 0; m < ls; ++m) {
        v[m] -= u_p * (p_mu_m[mu] * iv_p[m]);
        v[m] -= u_m * (p_mu_p[mu] * iv_m[m]);
      }
    }
  }
}

void multiply_wilson_ddag_e_o_no_comm(FermionField5d& out,
                                      const FermionField5d& in,
                                      const GaugeField& gf)
// set_left_expanded_gauge_field(gf, gf_);
// in.geo() = geo_resize(geo, 1);
// in.multiplicity = ls;
// refresh_expanded_1(in);
{
  TIMER("multiply_wilson_ddag_e_o_no_comm(5d,5d,gf)");
  Qassert(&out != &in);
  Qassert(is_matching_geo(gf.geo(), in.geo()));
  Qassert(in.geo().eo == 1 or in.geo().eo == 2);
  Geometry geo = geo_resize(in.geo());
  geo.eo = 3 - in.geo().eo;
  if (is_initialized(out) and out.geo().eo == 3 - geo.eo) {
    out.geo().eo = geo.eo;
  }
  out.init(geo, in.multiplicity);
  set_zero(out);
  const Int ls = in.multiplicity;
  Qassert(out.multiplicity == ls);
  Qassert(is_matching_geo(out.geo(), in.geo()));
  Qassert(out.geo().eo != in.geo().eo);
  Qassert(out.geo().eo == 1 or out.geo().eo == 2);
  const array<SpinMatrix, 4>& gammas = SpinMatrixConstants::get_cps_gammas();
  const SpinMatrix& unit = SpinMatrixConstants::get_unit();
  array<SpinMatrix, 4> p_mu_p;
  array<SpinMatrix, 4> p_mu_m;
  for (Int mu = 0; mu < 4; ++mu) {
    p_mu_p[mu] = (ComplexD)0.5 * (unit + gammas[mu]);
    p_mu_m[mu] = (ComplexD)0.5 * (unit - gammas[mu]);
  }
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<WilsonVector> v = out.get_elems(xl);
    for (Int mu = 0; mu < 4; ++mu) {
      const Coordinate xl_p = coordinate_shifts(xl, mu);
      const Coordinate xl_m = coordinate_shifts(xl, -mu - 1);
      const ColorMatrix u_p = gf.get_elem(xl, mu);
      const ColorMatrix u_m = matrix_adjoint(gf.get_elem(xl_m, mu));
      const Vector<WilsonVector> iv_p = in.get_elems_const(xl_p);
      const Vector<WilsonVector> iv_m = in.get_elems_const(xl_m);
      for (Int m = 0; m < ls; ++m) {
        v[m] -= u_p * (p_mu_p[mu] * iv_p[m]);
        v[m] -= u_m * (p_mu_m[mu] * iv_m[m]);
      }
    }
  }
}

void multiply_m_e_o(FermionField5d& out, const FermionField5d& in,
                    const GaugeField& gf, const FermionAction& fa)
// out can be the same object as in
// works for _o_e as well
{
  TIMER("multiply_m_e_o(5d,5d,gf,fa)");
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  const SpinMatrix& unit = SpinMatrixConstants::get_unit();
  const SpinMatrix p_p = (ComplexD)0.5 * (unit + gamma5);
  const SpinMatrix p_m = (ComplexD)0.5 * (unit - gamma5);
  const Int in_geo_eo = in.geo().eo;
  FermionField5d in1;
  in1.init(geo_resize(in.geo(), 1), in.multiplicity);
  const Geometry& geo = in.geo();
  std::vector<ComplexD> beo(fa.ls), ceo(fa.ls);
  for (Int m = 0; m < fa.ls; ++m) {
    beo[m] = fa.bs[m];
    ceo[m] = -fa.cs[m];
  }
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<WilsonVector> iv = in.get_elems_const(xl);
    Vector<WilsonVector> v = in1.get_elems(xl);
    for (Int m = 0; m < fa.ls; ++m) {
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
  Qassert(is_matching_geo(out.geo(), in.geo()));
  Qassert(out.geo().eo != in_geo_eo);
  Qassert(in.geo().eo == 1 or in.geo().eo == 2);
  Qassert(out.geo().eo == 1 or out.geo().eo == 2);
}

void multiply_mdag_e_o(FermionField5d& out, const FermionField5d& in,
                       const GaugeField& gf, const FermionAction& fa)
// out can be the same object as in
// works for _o_e as well
{
  TIMER("multiply_mdag_e_o(5d,5d,gf,fa)");
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  const SpinMatrix& unit = SpinMatrixConstants::get_unit();
  const SpinMatrix p_p = (ComplexD)0.5 * (unit + gamma5);
  const SpinMatrix p_m = (ComplexD)0.5 * (unit - gamma5);
  const Int in_geo_eo = in.geo().eo;
  Qassert(is_matching_geo(gf.geo(), in.geo()));
  Qassert(in.geo().eo == 1 or in.geo().eo == 2);
  Geometry geo = geo_resize(in.geo());
  geo.eo = 3 - in.geo().eo;
  std::vector<ComplexD> beo(fa.ls), ceo(fa.ls);
  for (Int m = 0; m < fa.ls; ++m) {
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
    for (Int m = 0; m < fa.ls; ++m) {
      v[m] = (ComplexD)beo[m] * iv[m];
      const WilsonVector tmp =
          (p_p * (m < fa.ls - 1
                      ? (WilsonVector)((ComplexD)ceo[m + 1] * iv[m + 1])
                      : (WilsonVector)((ComplexD)(-qconj((ComplexD)fa.mass) *
                                                  ceo[0]) *
                                       iv[0]))) +
          (p_m * (m > 0 ? (WilsonVector)((ComplexD)ceo[m - 1] * iv[m - 1])
                        : (WilsonVector)((ComplexD)(-qconj((ComplexD)fa.mass) *
                                                    ceo[fa.ls - 1]) *
                                         iv[fa.ls - 1])));
      v[m] -= tmp;
    }
  }
  Qassert(is_matching_geo(out.geo(), in.geo()));
  Qassert(out.geo().eo != in_geo_eo);
  Qassert(in.geo().eo == 1 or in.geo().eo == 2);
  Qassert(out.geo().eo == 1 or out.geo().eo == 2);
}

void multiply_m(FermionField5d& out, const FermionField5d& in,
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

void multiply_mdag(FermionField5d& out, const FermionField5d& in,
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

void multiply_mpc_sym2(FermionField5d& out, const FermionField5d& in,
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

void multiply_mpcdag_sym2(FermionField5d& out, const FermionField5d& in,
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

void multiply_hermop_sym2(FermionField5d& out, const FermionField5d& in,
                          const GaugeField& gf, const FermionAction& fa)
// out can be the same object as in
// odd <- odd (works for even <- even as well)
{
  TIMER_FLOPS("multiply_hermop_sym2");
  multiply_mpc_sym2(out, in, gf, fa);
  multiply_mpcdag_sym2(out, out, gf, fa);
  timer.flops += 5500 * fa.ls * gf.geo().local_volume();
}

void multiply_m_e_e(FermionField5d& out, const FermionField5d& in,
                    const InverterDomainWall& inv)
// out can be the same object as in
{
  multiply_m_e_e(out, in, inv.fa);
}

void multiply_mdag_e_e(FermionField5d& out, const FermionField5d& in,
                       const InverterDomainWall& inv)
// out can be the same object as in
{
  multiply_mdag_e_e(out, in, inv.fa);
}

void multiply_m_e_e_inv(FermionField5d& out, const FermionField5d& in,
                        const InverterDomainWall& inv)
// out can be the same object as in
{
  multiply_m_e_e_inv(out, in, inv.fa);
}

void multiply_mdag_e_e_inv(FermionField5d& out, const FermionField5d& in,
                           const InverterDomainWall& inv)
// out can be the same object as in
{
  multiply_mdag_e_e_inv(out, in, inv.fa);
}

void multiply_m_e_o(FermionField5d& out, const FermionField5d& in,
                    const InverterDomainWall& inv)
// out can be the same object as in
{
  multiply_m_e_o(out, in, inv.gf, inv.fa);
}

void multiply_mdag_e_o(FermionField5d& out, const FermionField5d& in,
                       const InverterDomainWall& inv)
// out can be the same object as in
{
  multiply_mdag_e_o(out, in, inv.gf, inv.fa);
}

void multiply_m_eo_eo(FermionField5d& out, const FermionField5d& in,
                      const InverterDomainWall& inv, const Int eo_out,
                      const Int eo_in)
// out can be the same object as in
// out need to be initialized with correct geo and eo
{
  TIMER("multiply_m_eo_eo");
  Geometry geo = geo_resize(in.geo());
  geo.eo = eo_out;
  out.init(geo, in.multiplicity);
  Qassert(is_matching_geo(out.geo(), in.geo()));
  Qassert(out.geo().eo == eo_out);
  Qassert(in.geo().eo == eo_in);
  if (eo_out == eo_in) {
    multiply_m_e_e(out, in, inv.fa);
  } else {
    multiply_m_e_o(out, in, inv.gf, inv.fa);
  }
}

void multiply_mdag_eo_eo(FermionField5d& out, const FermionField5d& in,
                         const InverterDomainWall& inv, const Int eo_out,
                         const Int eo_in)
// out can be the same object as in
// out need to be initialized with correct geo and eo
{
  TIMER("multiply_mdag_eo_eo");
  Geometry geo = geo_resize(in.geo());
  geo.eo = eo_out;
  out.init(geo, in.multiplicity);
  Qassert(is_matching_geo(out.geo(), in.geo()));
  Qassert(out.geo().eo == eo_out);
  Qassert(in.geo().eo == eo_in);
  if (eo_out == eo_in) {
    multiply_mdag_e_e(out, in, inv.fa);
  } else {
    multiply_mdag_e_o(out, in, inv.gf, inv.fa);
  }
}

void multiply_m(FermionField5d& out, const FermionField5d& in,
                const InverterDomainWall& inv)
// out can be the same object as in
{
  multiply_m(out, in, inv.gf, inv.fa);
}

void multiply_mdag(FermionField5d& out, const FermionField5d& in,
                   const InverterDomainWall& inv)
// out can be the same object as in
{
  multiply_mdag(out, in, inv.gf, inv.fa);
}

void multiply_mpc_sym2(FermionField5d& out, const FermionField5d& in,
                       const InverterDomainWall& inv)
// out can be the same object as in
// odd <- odd (works for even <- even as well)
{
  multiply_mpc_sym2(out, in, inv.gf, inv.fa);
}

void multiply_mpcdag_sym2(FermionField5d& out, const FermionField5d& in,
                          const InverterDomainWall& inv)
// out can be the same object as in
// odd <- odd (works for even <- even as well)
{
  multiply_mpcdag_sym2(out, in, inv.gf, inv.fa);
}

void multiply_hermop_sym2(FermionField5d& out, const FermionField5d& in,
                          const InverterDomainWall& inv)
// out can be the same object as in
// odd <- odd (works for even <- even as well)
{
  multiply_hermop_sym2(out, in, inv.gf, inv.fa);
}

void multiply_m_with_prec_sym2(FermionField5d& out, const FermionField5d& in,
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

ComplexD dot_product(const FermionField5d& ff1, const FermionField5d& ff2)
// return ff1^dag * ff2
{
  TIMER("dot_product");
  Qassert(is_matching_geo(ff1.geo(), ff2.geo()));
  Qassert(ff1.geo().eo == ff2.geo().eo);
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
    for (Int k = 0; k < cv1.size(); ++k) {
      sum += qconj(cv1[k]) * cv2[k];
    }
  }
  glb_sum(sum);
  return sum;
}

Long cg_with_herm_sym_2(FermionField5d& sol, const FermionField5d& src,
                        const InverterDomainWall& inv, const double stop_rsd,
                        const Long max_num_iter)
{
  TIMER_VERBOSE_FLOPS("cg_with_herm_sym_2(5d,5d,inv)");
  Qassert(&sol != &src);
  const Long iter =
      cg_with_f(sol, src, inv, multiply_hermop_sym2, stop_rsd, max_num_iter);
  timer.flops += 5500 * iter * inv.fa.ls * inv.geo().local_volume();
  return iter;
}

Long invert(FermionField5d& out, const FermionField5d& in,
            const InverterDomainWall& inv)
{
  Qassert(&out != &in);
  return invert_with_cg(out, in, inv, cg_with_herm_sym_2);
}

Long invert(FermionField4d& out, const FermionField4d& in,
            const InverterDomainWall& inv)
{
  Qassert(&out != &in);
  return invert_dwf(out, in, inv);
}

double find_max_eigen_value_hermop_sym2(const InverterDomainWall& inv,
                                        const RngState& rs, const Long max_iter)
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
