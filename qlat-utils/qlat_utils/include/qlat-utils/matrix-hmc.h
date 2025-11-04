#pragma once

#include <qlat-utils/matrix.h>
#include <qlat-utils/rng-state.h>
#include <qlat-utils/utils-vec.h>

namespace qlat
{  //

template <class T, std::size_t N>
qacc void normalize_array_real(Array<T, N> p)
//------------------------------------------------------------------
// Few routines that are needed by the Unitarize routine
//------------------------------------------------------------------
{
  RealD norm = 0;
  for (Long i = 0; i < (Long)N; i++) {
    norm += p[i] * p[i];
  }
  if (not(norm == 1.0)) {
    norm = 1.0 / sqrt(norm);
    for (Long i = 0; i < (Long)N; i++) {
      p[i] *= norm;
    }
  }
}

template <class T, std::size_t N>
qacc void normalize_array_complex(Array<ComplexT<T>, N> p)
//------------------------------------------------------------------
// Few routines that are needed by the Unitarize routine
//------------------------------------------------------------------
{
  Array<T, N * 2> p_r((T*)p.data());
  normalize_array_real(p_r);
}

template <class T, std::size_t N>
qacc void orthogonalize_array_complex(Array<T, N> p2, const Array<T, N> p1)
//	v2' = v2 - v1 * (v1^*, v2)
// 	then	(v1^*, v2') = 0
{
  ComplexD c = 0.0;
  for (Long i = 0; i < (Long)N; ++i) {
    c += qconj(p1[i]) * p2[i];
  }
  if (not(c == 0.0)) {
    for (Long i = 0; i < (Long)N; ++i) {
      p2[i] -= c * p1[i];
    }
  }
}

template <class T>
qacc void cross_product_conj(Array<T, 3> v3, const Array<T, 3> v1,
                             const Array<T, 3> v2)
// v3 = ( v1 x v2 )^*
{
  v3[0] = qconj(v1[1] * v2[2] - v1[2] * v2[1]);
  v3[1] = qconj(v1[2] * v2[0] - v1[0] * v2[2]);
  v3[2] = qconj(v1[0] * v2[1] - v1[1] * v2[0]);
}

template <class T>
qacc void unitarize(ColorMatrixT<T>& cm)
{
  Array<ComplexT<T>, 3> p1(cm.data());
  Array<ComplexT<T>, 3> p2(cm.data() + 3);
  Array<ComplexT<T>, 3> p3(cm.data() + 6);
  normalize_array_complex(p1);
  orthogonalize_array_complex(p2, p1);
  normalize_array_complex(p2);
  cross_product_conj(p3, p1, p2);
}

qacc ColorMatrix make_anti_hermitian_matrix(const array<RealD, 8>& basis)
{
  qassert(3 == NUM_COLOR);
  ColorMatrix m;
  Array<RealD, 18> p(m.d());
  const RealD s3 = (1.0 / std::sqrt(3.0)) * basis[7];
  p[0] = 0.0;
  p[8] = 0.0;
  p[16] = 0.0;
  p[1] = basis[2] + s3;
  p[9] = -basis[2] + s3;
  p[17] = -2.0 * s3;
  p[2] = basis[1];
  p[3] = basis[0];
  p[4] = basis[4];
  p[5] = basis[3];
  p[6] = -basis[1];
  p[7] = basis[0];
  p[10] = basis[6];
  p[11] = basis[5];
  p[12] = -basis[4];
  p[13] = basis[3];
  p[14] = -basis[6];
  p[15] = basis[5];
  return m;
}

qacc array<RealD, 8> basis_projection_anti_hermitian_matrix(
    const ColorMatrix& m)
// m is a tr_less_anti_hermitian_matrix
{
  const Array<RealD, 18> p(m.d());
  array<RealD, 8> basis;
  basis[0] = p[3];
  basis[1] = p[2];
  basis[2] = 0.5 * (p[1] - p[9]);
  basis[3] = p[5];
  basis[4] = p[4];
  basis[5] = p[11];
  basis[6] = p[10];
  basis[7] = std::sqrt(3.0) * 0.5 * (p[1] + p[9]);
  return basis;
}

inline ColorMatrix make_g_rand_anti_hermitian_matrix(RngState& rs,
                                                     const RealD sigma)
//  Creates an anti-hermitian 3x3 complex matrix with each complex
//  element drawn at random from a Gaussian distribution with zero mean.
//  Hence the matrices are distributed according to
//
//  exp[- Tr(mat^2)/(2 sigma**2)]
{
  const RealD s = sigma / std::sqrt(2);
  array<RealD, 8> a;
  for (Int i = 0; i < 8; ++i) {
    a[i] = g_rand_gen(rs, 0.0, s);
  }
  return make_anti_hermitian_matrix(a);
}

qacc RealD neg_half_tr_square(const ColorMatrix& m)
// ret == (basis**2).sum()
// basis = basis_projection_anti_hermitian_matrix(m)
{
  const Array<RealD, 18> p(m.d());
  return sqr(p[3]) + sqr(p[2]) + sqr(p[1] - p[9]) * 0.25 + sqr(p[5]) +
         sqr(p[4]) + sqr(p[11]) + sqr(p[10]) +
         sqr(p[17]) * 0.75;  // 0.75 = (sqrt(3)/2)^2
}

template <class M>
qacc M make_matrix_exp(const M& a, const Int max_order = 20)
{
  M unit;
  set_unit(unit);
  M t2 = a;
  M t3;
  for (Int j = max_order; j > 1; --j) {
    t3 = unit + (1.0 / j) * t2;
    t2 = a * t3;
  }
  t3 = unit + t2;
  return t3;
}

template <class T>
qacc ColorMatrixT<T> make_color_matrix_exp(const ColorMatrixT<T>& a,
                                           const Int max_order = 9)
{
  ColorMatrixT<T> ret = make_matrix_exp(a, max_order);
  unitarize(ret);
  return ret;
}

template <class T>
qacc ColorMatrixT<T> matrix_evolve(const ColorMatrixT<T>& gf_cm,
                                   const ColorMatrixT<T>& gm_cm,
                                   const RealD step_size)
// return exp(gm_cm * step_size) * gf_cm
{
  const ColorMatrixT<T> t = (ComplexT<T>)step_size * gm_cm;
  return make_matrix_exp(t) * gf_cm;
}

template <class T>
qacc ColorMatrixT<T> matrix_evolve_dual(const ColorMatrixT<T>& gf_cm,
                                        const ColorMatrixT<T>& gm_cm,
                                        const RealD step_size)
// return gf_cm * exp(-gm_cm * step_size)
{
  const ColorMatrixT<T> t = (ComplexT<T>)(-step_size) * gm_cm;
  return gf_cm * make_matrix_exp(t);
}

qacc ColorMatrix make_tr_less_anti_herm_matrix(const ColorMatrix& m)
// (m - m^\dagger) / 2 - Tr(m - m^\dagger) / 6
{
  ColorMatrix ret = m - matrix_adjoint(m);
  ret *= 0.5;
  Array<RealD, 18> p(ret.d());
  const RealD c = (p[1] + p[9] + p[17]) / 3.0;
  p[1] -= c;
  p[9] -= c;
  p[17] -= c;
  return ret;
}

struct API ColorMatrixConstants;

qacc AdjointColorMatrix make_adjoint_representation(
    const ColorMatrix& m, const ColorMatrixConstants& cmcs);

qacc AdjointColorMatrix make_diff_exp_map(const AdjointColorMatrix& am,
                                          const Int max_order = 20);

qacc AdjointColorMatrix make_diff_exp_map(const ColorMatrix& m,
                                          const ColorMatrixConstants& cmcs,
                                          const Int max_order = 20);

qacc AdjointColorMatrix make_diff_exp_map_diff_ref(
    const ColorMatrix& m, const Int a, const ColorMatrixConstants& cmcs);

qacc AdjointColorMatrix make_diff_exp_map_diff(const AdjointColorMatrix& am,
                                               const Int a,
                                               const ColorMatrixConstants& cmcs,
                                               const Int max_order = 20);

qacc AdjointColorMatrix make_diff_exp_map_diff(const ColorMatrix& m,
                                               const Int a,
                                               const ColorMatrixConstants& cmcs,
                                               const Int max_order = 20);

struct API ColorMatrixConstants {
  //
  // Luscher's convention: https://arxiv.org/abs/0907.5491v3
  //
  // except the norm of T^a matrices: tr(T^a T^b) = -2 \delta^{a,b}
  //
  ColorMatrix unit;
  array<ColorMatrix, 8> ts;
  AdjointColorMatrix aunit;
  array<AdjointColorMatrix, 8> f;
  //
  qacc ColorMatrixConstants() { init(); }
  //
  qacc void init()
  {
    set_unit(unit);
    for (Int a = 0; a < 8; ++a) {
      array<RealD, 8> basis;
      set_zero(basis);
      basis[a] = 1.0;
      ts[a] = make_anti_hermitian_matrix(basis);
    }
    set_unit(aunit);
    for (Int a = 0; a < 8; ++a) {
      for (Int b = 0; b < 8; ++b) {
        const ColorMatrix m = ts[a] * ts[b] - ts[b] * ts[a];
        array<RealD, 8> basis = basis_projection_anti_hermitian_matrix(m);
        for (Int c = 0; c < 8; ++c) {
          f[c](a, b) = basis[c];
        }
      }
    }
  }
  //
  void check() const;
  //
  API static const box<ColorMatrixConstants>& get_instance_box()
  {
    static box<ColorMatrixConstants> cmcs =
        box<ColorMatrixConstants>(ColorMatrixConstants(), MemType::Uvm);
    return cmcs;
  }
  //
  API static const ColorMatrixConstants& get_instance()
  {
    return get_instance_box()();
  }
  //
  API static const ColorMatrix& get_unit() { return get_instance().unit; }
  //
  API static const array<ColorMatrix, 8>& get_ts() { return get_instance().ts; }
  //
  API static const AdjointColorMatrix& get_aunit()
  {
    return get_instance().aunit;
  }
  //
  API static const array<AdjointColorMatrix, 8>& get_f()
  {
    return get_instance().f;
  }
};

qacc AdjointColorMatrix make_adjoint_representation(
    const ColorMatrix& m, const ColorMatrixConstants& cmcs)
{
  const array<AdjointColorMatrix, 8>& f = cmcs.f;
  const array<RealD, 8> basis = basis_projection_anti_hermitian_matrix(m);
  AdjointColorMatrix am;
  set_zero(am);
  for (Int a = 0; a < 8; ++a) {
    am += (-basis[a]) * f[a];
  }
  return am;
}

qacc AdjointColorMatrix make_diff_exp_map(const AdjointColorMatrix& am,
                                          const Int max_order)
{
  AdjointColorMatrix aunit;
  set_unit(aunit);
  AdjointColorMatrix t2 = -am;
  AdjointColorMatrix t3;
  for (Int j = max_order; j > 1; --j) {
    t3 = aunit + (1.0 / (j + 1)) * t2;
    t2 = -am * t3;
  }
  t3 = aunit + 0.5 * t2;
  return t3;
}

qacc AdjointColorMatrix make_diff_exp_map(const ColorMatrix& m,
                                          const ColorMatrixConstants& cmcs,
                                          const Int max_order)
{
  return make_diff_exp_map(make_adjoint_representation(m, cmcs), max_order);
}

qacc AdjointColorMatrix make_diff_exp_map_diff_ref(
    const ColorMatrix& m, const Int a, const ColorMatrixConstants& cmcs)
{
  const array<ColorMatrix, 8>& ts = cmcs.ts;
  const RealD eps = 1e-4;
  const ColorMatrix m_p = m + (ComplexD)eps * ts[a];
  const ColorMatrix m_n = m - (ComplexD)eps * ts[a];
  const AdjointColorMatrix j_p = make_diff_exp_map(m_p, cmcs);
  const AdjointColorMatrix j_n = make_diff_exp_map(m_n, cmcs);
  return (1.0 / (2.0 * eps)) * (j_p - j_n);
}

qacc AdjointColorMatrix make_diff_exp_map_diff(const AdjointColorMatrix& am,
                                               const Int a,
                                               const ColorMatrixConstants& cmcs,
                                               const Int max_order)
{
  const array<AdjointColorMatrix, 8>& f = cmcs.f;
  const AdjointColorMatrix& aunit = cmcs.aunit;
  AdjointColorMatrix t2 = -am;
  AdjointColorMatrix dt2 = f[a];
  AdjointColorMatrix t3;
  AdjointColorMatrix dt3;
  for (Int j = max_order; j > 1; --j) {
    t3 = aunit + (1.0 / (j + 1)) * t2;
    dt3 = (1.0 / (j + 1)) * dt2;
    t2 = -am * t3;
    dt2 = f[a] * t3 - am * dt3;
  }
  dt3 = 0.5 * dt2;
  return dt3;
}

qacc AdjointColorMatrix make_diff_exp_map_diff(const ColorMatrix& m,
                                               const Int a,
                                               const ColorMatrixConstants& cmcs,
                                               const Int max_order)
{
  return make_diff_exp_map_diff(make_adjoint_representation(m, cmcs), a, cmcs,
                                max_order);
}

}  // namespace qlat
