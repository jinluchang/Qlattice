#include <qlat-utils/matrix.h>
#include <qlat-utils/qacc-func.h>
#include <qlat-utils/types.h>
#include <qlat/core.h>
#include <qlat/qed.h>
#include <qlat/dslash.h>

namespace qlat
{  //

void set_left_expanded_gauge_field(Field<ComplexD>& gf1,
                                   const Field<ComplexD>& gf)
{
  TIMER("set_left_expanded_gauge_field(qed-gf)");
  const Coordinate expansion_left(1, 1, 1, 1);
  const Coordinate expansion_right(0, 0, 0, 0);
  const Geometry geo1 = geo_resize(gf.geo(), expansion_left, expansion_right);
  gf1.init(geo1, 4);
  Qassert(gf1.geo() == geo1);
  gf1 = gf;
  refresh_expanded_1(gf1);
}

void free_invert(SpinProp& sp_sol, SpinProp& sp_src, const RealD mass,
                 const RealD m5, const CoordinateD& momtwist)
{
  TIMER("free_invert(qed-ff)");
  sp_sol.init(sp_src);
  prop_spin_propagator4d(sp_sol, mass, m5, momtwist);
}

void invert_qed(SpinProp& sp_sol, const SpinProp& sp_src, const Field<ComplexD>& gf1,
                const RealD mass, const RealD m5, const Int ls,
                const bool is_dagger)
{
  TIMER("invert_qed(qed-sp)");
  Qassert(sp_src.multiplicity == 1);
  const Geometry& geo = sp_src.geo();
  Qassert(geo.is_only_local);
  Qassert(&sp_sol != &sp_src);
  sp_sol.init(geo, sp_src.multiplicity);
  Qassert(sp_sol.geo() == geo);
  Qassert(sp_sol.multiplicity == 1);
  Field<ComplexD> sp_sol_c, sp_src_c;
  qswap_cast(sp_sol_c, sp_sol);
  sp_src_c.set_view_cast(sp_src);
  Qassert(sp_sol_c.multiplicity == 4 * 4);
  Qassert(sp_src_c.multiplicity == 4 * 4);
  Field<ComplexD> sol4d, src4d;
  sol4d.init(geo, 4);
  src4d.init(geo, 4);
  for (Int i = 0; i < 4; ++i) {
    for (Int j = 0; j < 4; ++j) {
      set_field_m(src4d, sp_src_c, j, j * 4 + i);
    }
    invert_dwf_qed(sol4d, src4d, gf1, mass, m5, ls, is_dagger);
    for (Int j = 0; j < 4; ++j) {
      set_field_m(sp_sol_c, sol4d, j * 4 + i, j);
    }
  }
  qswap_cast(sp_sol, sp_sol_c);
}

void fermion_field_4d_from_5d_qed(Field<ComplexD>& ff4d,
                                  const Field<ComplexD>& ff5d, const Int ls,
                                  const Int upper, const Int lower)
// upper componets are right handed
// lower componets are left handed
{
  TIMER("fermion_field_4d_from_5d_qed(qed-ff)");
  const Geometry& geo = ff5d.geo();
  Qassert(geo.is_only_local);
  Qassert(4 * ls == ff5d.multiplicity);
  ff4d.init(geo, 4);
  set_zero(ff4d);
  qacc_for(index, geo.local_volume(), {
    const Vector<ComplexD> iv = ff5d.get_elems_const(index);
    Vector<ComplexD> v = ff4d.get_elems(index);
    v[0] = iv[upper * 4 + 0];
    v[1] = iv[upper * 4 + 1];
    v[2] = iv[lower * 4 + 2];
    v[3] = iv[lower * 4 + 3];
  });
}

void fermion_field_5d_from_4d_qed(Field<ComplexD>& ff5d,
                                  const Field<ComplexD>& ff4d, const Int ls,
                                  const Int upper, const Int lower)
// upper componets are right handed
// lower componets are left handed
{
  TIMER("fermion_field_5d_from_4d_qed");
  const Geometry& geo = ff4d.geo();
  Qassert(geo.is_only_local);
  ff5d.init(geo, 4 * ls);
  set_zero(ff5d);
  qacc_for(index, geo.local_volume(), {
    const Vector<ComplexD> iv = ff4d.get_elems_const(index);
    Vector<ComplexD> v = ff5d.get_elems(index);
    v[upper * 4 + 0] = iv[0];
    v[upper * 4 + 1] = iv[1];
    v[lower * 4 + 2] = iv[2];
    v[lower * 4 + 3] = iv[3];
  });
}

Long invert_dwf_qed(Field<ComplexD>& out, const Field<ComplexD>& in,
                    const Field<ComplexD>& gf1, const RealD mass,
                    const RealD m5, const Int ls, const bool is_dagger)
// properly project to 4d fermion field
// if is_dagger is false (default), then M out = in
// if is_dagger is true, then M^dag out = in
{
  TIMER("invert_dwf_qed(qed-ff)");
  const Geometry& geo = in.geo();
  Qassert(geo.is_only_local);
  Field<ComplexD> sol5d, src5d;
  sol5d.init(geo, 4 * ls);
  src5d.init(geo, 4 * ls);
  fermion_field_5d_from_4d_qed(src5d, in, ls, 0, ls - 1);
  set_zero(sol5d);
  multiply_m_dwf_qed(src5d, src5d, gf1, mass, m5, ls, is_dagger);
  const Long iter = cg_with_m_dwf_qed(sol5d, src5d, gf1, mass, m5, ls, is_dagger);
  fermion_field_4d_from_5d_qed(out, sol5d, ls, ls - 1, 0);
  return iter;
}

ComplexD dot_product(const Field<ComplexD>& ff1, const Field<ComplexD>& ff2)
// return ff1^dag * ff2
{
  TIMER("dot_product(qed-ff)");
  Qassert(ff1.geo().is_only_local);
  Qassert(ff2.geo().is_only_local);
  Qassert(ff1.geo() == ff2.geo());
  Qassert(ff1.multiplicity == ff2.multiplicity);
  const Geometry& geo = ff1.geo();
  Field<ComplexD> fsum;
  fsum.init(geo, 1);
  qacc_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<ComplexD> v1 = ff1.get_elems_const(xl);
    const Vector<ComplexD> v2 = ff2.get_elems_const(xl);
    Qassert(v1.size() == v2.size());
    ComplexD s = 0.0;
    for (Int k = 0; k < v1.size(); ++k) {
      s += qconj(v1[k]) * v2[k];
    }
    fsum.get_elem(index) = s;
  });
  const std::vector<ComplexD> sum = field_glb_sum(fsum);
  return sum[0];
}

Long cg_with_m_dwf_qed(Field<ComplexD>& out, const Field<ComplexD>& in,
                       const Field<ComplexD>& gf1, const RealD mass,
                       const RealD m5, const Int ls, const bool is_dagger,
                       const RealD stop_rsd, const Long max_num_iter)
// if is_dagger is false, then M^dag M out = in
// if is_dagger is true, then M M^dag out = in
{
  TIMER("cg_with_m_dwf_qed");
  Qassert(mass > 0.0);
  Qassert(ls > 0);
  Qassert(stop_rsd >= 0.0);
  Qassert(max_num_iter >= 0);
  Qassert(in.geo().is_only_local);
  // Treat out as initial guess (zero if not provided)
  out.init_zero(in.geo(), in.multiplicity);
  Qassert(out.geo().is_only_local);
  Qassert(out.multiplicity == in.multiplicity);
  Qassert(gf1.multiplicity == 4);
  Qassert(is_matching_geo(gf1.geo(), in.geo()));
  if (max_num_iter == 0) {
    return 0;
  }
  // implement conjugate gradient with normal equation solver
  Field<ComplexD> r, p, tmp, ap;
  r = in;
  multiply_m_dwf_qed(tmp, out, gf1, mass, m5, ls, is_dagger);
  multiply_m_dwf_qed(tmp, out, gf1, mass, m5, ls, not is_dagger);
  r -= tmp;
  p = r;
  const RealD qnorm_in = qnorm(in);
  displayln_info(
      fname +
      ssprintf(
          ": start max_num_iter=%4ld        sqrt(qnorm_in)=%.3E stop_rsd=%.3E",
          max_num_iter, sqrt(qnorm_in), stop_rsd));
  RealD qnorm_r = qnorm(r);
  for (Long iter = 1; iter <= max_num_iter; ++iter) {
    multiply_m_dwf_qed(ap, p, gf1, mass, m5, ls, is_dagger);
    multiply_m_dwf_qed(ap, p, gf1, mass, m5, ls, not is_dagger);
    const RealD alpha = qnorm_r / dot_product(p, ap).real();
    tmp = p;
    tmp *= alpha;
    out += tmp;
    tmp = ap;
    tmp *= alpha;
    r -= tmp;
    const RealD new_qnorm_r = qnorm(r);
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
    const RealD beta = new_qnorm_r / qnorm_r;
    p *= beta;
    p += r;
    qnorm_r = new_qnorm_r;
  }
  displayln_info(
      fname +
      ssprintf(
          ": final max_num_iter=%4ld sqrt(qnorm_r/qnorm_in)=%.3E stop_rsd=%.3E",
          max_num_iter + 1, sqrt(qnorm_r / qnorm_in), stop_rsd));
  return max_num_iter + 1;
}



void multiply_m_dwf_qed(Field<ComplexD>& out, const Field<ComplexD>& in,
                        const Field<ComplexD>& gf1, const RealD mass,
                        const RealD m5, const Int ls, const bool is_dagger)
// set_left_expanded_gauge_field(gf1, gf);
// in.geo() should not be expanded.
// out can be the same object as in.
// mass is the masss of the fermion, for example, its value can be `0.1`.
// m5 should typically to 1.0
// ls should typically to a large even integer, such as 64.
{
  TIMER("multiply_m_dwf_qed");
  Qassert(mass > 0.0);
  Qassert(ls > 0);
  const Geometry& geo = in.geo();
  const Geometry& geo_gf1 = gf1.geo();
  Qassert(is_matching_geo(geo, geo_gf1));
  const Coordinate size_node = geo_gf1.geon.size_node;
  for (Int mu = 0; mu < DIMN; ++mu) {
    if (size_node[mu] > 1) {
      Qassert(geo_gf1.expansion_left[mu] == 1);
      Qassert(geo_gf1.expansion_right[mu] == 0);
    }
  }
  Qassert(gf1.multiplicity == 4);
  Qassert(in.multiplicity == 4 * ls);
  Qassert(geo.is_only_local == true);
  const Coordinate expansion_left(1, 1, 1, 1);
  const Coordinate expansion_right(1, 1, 1, 1);
  const Geometry geo1 = geo_resize(geo, expansion_left, expansion_right);
  Qassert(geo.eo == 0);
  Qassert(geo1.eo == 0);
  Qassert(geo_gf1.eo == 0);
  const box<SpinMatrixConstants>& smc = get_spin_matrix_constants();
  Field<ComplexD> in1;
  in1.init(geo1, 4 * ls);
  in1 = in;
  refresh_expanded_1(in1);
  out.init(geo, 4 * ls);
  Qassert(out.geo().is_only_local);
  set_zero(out);
  out = in1;
  const RealD dagger_factor = is_dagger ? -1.0 : 1.0;
  qacc_for(index, geo.local_volume(), {
    const SpinMatrix& gamma_5 = smc().gamma5;
    const array<SpinMatrix, 4>& gammas = smc().cps_gammas;
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<ComplexD> v = out.get_elems(index);
    {
      const Vector<ComplexD> iv = in1.get_elems_const(xl);
      vec_plusm(v, (ComplexD)(5.0 - m5), iv);
      const Vector<ComplexD> iv_m(iv.p + 4, 4 * (ls - 1));
      const Vector<ComplexD> iv_p(iv.p, 4 * (ls - 1));
      const Vector<ComplexD> iv_m0(iv.p, 4);
      const Vector<ComplexD> iv_p0(iv.p + 4 * (ls - 1), 4);
      Vector<ComplexD> v_m(v.p, 4 * (ls - 1));
      Vector<ComplexD> v_p(v.p + 4, 4 * (ls - 1));
      Vector<ComplexD> v_m0(v.p + 4 * (ls - 1), 4);
      Vector<ComplexD> v_p0(v.p, 4);
      vec_plusm(v_m, (ComplexD)(-0.5), iv_m);
      vec_plusm(v_p, (ComplexD)(-0.5), iv_p);
      vec_plusm(v_m0, (ComplexD)(0.5 * mass), iv_m0);
      vec_plusm(v_p0, (ComplexD)(0.5 * mass), iv_p0);
      mat_mul_multi_vec_plusm(v_m, (ComplexD)(dagger_factor * 0.5), gamma_5,
                              iv_m);
      mat_mul_multi_vec_plusm(v_p, (ComplexD)(-dagger_factor * 0.5), gamma_5,
                              iv_p);
      mat_mul_multi_vec_plusm(v_m0, (ComplexD)(-dagger_factor * 0.5 * mass),
                              gamma_5, iv_m0);
      mat_mul_multi_vec_plusm(v_p0, (ComplexD)(dagger_factor * 0.5 * mass),
                              gamma_5, iv_p0);
    }
    for (Int mu = 0; mu < 4; ++mu) {
      const Coordinate xl_p = coordinate_shifts(xl, mu);
      const Coordinate xl_m = coordinate_shifts(xl, -mu - 1);
      ComplexD u_p = gf1.get_elem(xl, mu);
      ComplexD u_m = 1.0 / gf1.get_elem(xl_m, mu);
      if (is_dagger) {
        u_p = qconj(1.0 / u_p);
        u_m = qconj(1.0 / u_m);
      }
      const Vector<ComplexD> iv_p = in1.get_elems_const(xl_p);
      const Vector<ComplexD> iv_m = in1.get_elems_const(xl_m);
      vec_plusm(v, (ComplexD)(-0.5 * u_p), iv_p);
      vec_plusm(v, (ComplexD)(-0.5 * u_m), iv_m);
      mat_mul_multi_vec_plusm(v, (ComplexD)(dagger_factor * 0.5 * u_p),
                              gammas[mu], iv_p);
      mat_mul_multi_vec_plusm(v, (ComplexD)(-dagger_factor * 0.5 * u_m),
                              gammas[mu], iv_m);
    }
  });
}

}  // namespace qlat
