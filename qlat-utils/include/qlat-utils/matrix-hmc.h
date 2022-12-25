#pragma once

#include <qlat-utils/matrix.h>
#include <qlat-utils/rng-state.h>
#include <qlat-utils/qutils-vec.h>

namespace qlat
{  //

template <class T>
qacc void unitarize(ColorMatrixT<T>& cm)
{
  cm.em().row(0).normalize();
  cm.em().row(2) =
      cm.em().row(1) - cm.em().row(0).dot(cm.em().row(1)) * cm.em().row(0);
  cm.em().row(1) = cm.em().row(2).normalized();
  cm.em().row(2) = cm.em().row(0).cross(cm.em().row(1));
}

qacc ColorMatrixT<> make_anti_hermitian_matrix(const array<double, 8>& basis)
{
  qassert(3 == NUM_COLOR);
  ColorMatrixT<> m;
  Array<double, 18> p(m.d());
  const double s3 = (1.0 / std::sqrt(3.0)) * basis[7];
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

qacc array<double, 8> basis_projection_anti_hermitian_matrix(
    const ColorMatrix& m)
// m is a tr_less_anti_hermitian_matrix
{
  const Array<double, 18> p(m.d());
  array<double, 8> basis;
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

inline ColorMatrixT<> make_g_rand_anti_hermitian_matrix(RngState& rs,
                                                        const double sigma)
//  Creates an antihermitian 3x3 complex matrix with each complex
//  element drawn at random from a gaussian distribution with zero mean.
//  Hence the matrices are distributed according to
//
//  exp[- Tr(mat^2)/(2 sigma**2)]
{
  const double s = sigma / std::sqrt(2);
  array<double, 8> a;
  for (int i = 0; i < 8; ++i) {
    a[i] = g_rand_gen(rs, 0.0, s);
  }
  return make_anti_hermitian_matrix(a);
}

qacc double neg_half_tr_square(const ColorMatrixT<>& m)
{
  const Array<double, 18> p(m.d());
  return sqr(p[3]) + sqr(p[2]) + sqr(p[1] - p[9]) * 0.25 + sqr(p[5]) +
         sqr(p[4]) + sqr(p[11]) + sqr(p[10]) +
         sqr(p[17]) * 0.75;  // 0.75 = (sqrt(3)/2)^2
}

template <class M>
qacc M make_matrix_exp(const M& a, const int max_order = 20)
{
  M unit;
  set_unit(unit);
  M t2 = a;
  M t3;
  for (int j = max_order; j > 1; --j) {
    t3 = unit + (1.0 / j) * t2;
    t2 = a * t3;
  }
  t3 = unit + t2;
  return t3;
}

template <class T>
qacc ColorMatrixT<T> make_color_matrix_exp(const ColorMatrixT<T>& a,
                                           const int max_order = 9)
{
  ColorMatrixT<T> ret = make_matrix_exp(a, max_order);
  unitarize(ret);
  return ret;
}

template <class T>
qacc ColorMatrixT<T> matrix_evolve(const ColorMatrixT<T>& gf_cm,
                                   const ColorMatrixT<T>& gm_cm,
                                   const double step_size)
{
  const ColorMatrixT<T> t = (ComplexT<T>)step_size * gm_cm;
  return make_matrix_exp(t) * gf_cm;
}

qacc ColorMatrixT<> make_tr_less_anti_herm_matrix(const ColorMatrixT<>& m)
// (m - m^\dagger) / 2 - Tr(m - m^\dagger) / 6
{
  ColorMatrixT<> ret = m - matrix_adjoint(m);
  ret *= 0.5;
  Array<double, 18> p(ret.d());
  const double c = (p[1] + p[9] + p[17]) / 3.0;
  p[1] -= c;
  p[9] -= c;
  p[17] -= c;
  return ret;
}

struct API AdjointColorMatrix : MatrixT<8, double> {
  qacc AdjointColorMatrix() {}
  qacc AdjointColorMatrix(const MatrixT<8, double>& m) { *this = m; }
  //
  qacc const AdjointColorMatrix& operator=(const MatrixT<8, double>& m)
  {
    *this = (const AdjointColorMatrix&)m;
    return *this;
  }
};

struct API ColorMatrixConstants;

qacc AdjointColorMatrix make_adjoint_representation(
    const ColorMatrix& m, const ColorMatrixConstants& cmcs);

qacc AdjointColorMatrix make_diff_exp_map(const AdjointColorMatrix& am,
                                          const int max_order = 20);

qacc AdjointColorMatrix make_diff_exp_map(const ColorMatrix& m,
                                          const ColorMatrixConstants& cmcs,
                                          const int max_order = 20);

qacc AdjointColorMatrix make_diff_exp_map_diff_ref(
    const ColorMatrix& m, const int a, const ColorMatrixConstants& cmcs);

qacc AdjointColorMatrix make_diff_exp_map_diff(const AdjointColorMatrix& am,
                                               const int a,
                                               const ColorMatrixConstants& cmcs,
                                               const int max_order = 20);

qacc AdjointColorMatrix make_diff_exp_map_diff(const ColorMatrix& m,
                                               const int a,
                                               const ColorMatrixConstants& cmcs,
                                               const int max_order = 20);

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
    for (int a = 0; a < 8; ++a) {
      array<double, 8> basis;
      set_zero(basis);
      basis[a] = 1.0;
      ts[a] = make_anti_hermitian_matrix(basis);
    }
    set_unit(aunit);
    for (int a = 0; a < 8; ++a) {
      for (int b = 0; b < 8; ++b) {
        const ColorMatrix m = ts[a] * ts[b] - ts[b] * ts[a];
        array<double, 8> basis = basis_projection_anti_hermitian_matrix(m);
        for (int c = 0; c < 8; ++c) {
          f[c](a, b) = basis[c];
        }
      }
    }
  }
  //
  void check() const
  {
    for (int a = 0; a < 8; ++a) {
      ColorMatrix m = make_tr_less_anti_herm_matrix(ts[a]);
      qassert(qnorm(m - ts[a]) < 1e-20);
    }
    for (int a = 0; a < 8; ++a) {
      array<double, 8> basis;
      set_zero(basis);
      basis[a] = 1.0;
      qassert(qnorm(make_anti_hermitian_matrix(basis) - ts[a]) < 1e-20);
      qassert(qnorm(basis - basis_projection_anti_hermitian_matrix(ts[a])) <
              1e-20);
    }
    for (int a = 0; a < 8; ++a) {
      for (int b = 0; b < 8; ++b) {
        Complex tr = matrix_trace(ts[a], ts[b]);
        if (a == b) {
          tr -= -2;
        }
        if (qnorm(tr) >= 1e-20) {
          displayln_info(
              ssprintf("a=%d b=%d tr-dif=%s", a, b, show(tr).c_str()));
          qassert(false);
        }
      }
    }
    RngState rs("ColorMatrixConstants::check");
    ColorMatrix x = make_g_rand_anti_hermitian_matrix(rs, 1.0);
    qassert(qnorm(x - make_tr_less_anti_herm_matrix(x)) < 1e-20);
    AdjointColorMatrix adx = make_adjoint_representation(x, *this);
    array<double, 8> basis = basis_projection_anti_hermitian_matrix(x);
    for (int a = 0; a < 8; ++a) {
      for (int b = 0; b < 8; ++b) {
        Complex sum = 0.0;
        for (int c = 0; c < 8; ++c) {
          sum += -f[c](a, b) * basis[c];
        }
        qassert(qnorm(adx(a, b) - sum) < 1e-20);
      }
    }
    const double coef = 0.2;
    const AdjointColorMatrix exp_adx = make_matrix_exp(coef * adx);
    for (int b = 0; b < 8; ++b) {
      array<double, 8> basis_b;
      for (int a = 0; a < 8; ++a) {
        basis_b[a] = exp_adx(a, b);
      }
      const ColorMatrix mx_b = make_anti_hermitian_matrix(basis_b);
      const ColorMatrix exp_x = make_matrix_exp((Complex)coef * x);
      const ColorMatrix exp_n_x = make_matrix_exp(-(Complex)coef * x);
      if (qnorm(mx_b - exp_x * ts[b] * exp_n_x) >= 1e-20) {
        displayln_info("adx_b: " + show(qnorm(mx_b)));
        displayln_info("exp_x b exp_n_x: " +
                       show(qnorm(exp_x * ts[b] * exp_n_x)));
        displayln_info("diff: " + show(qnorm(mx_b - exp_x * ts[b] * exp_n_x)));
        qassert(false);
      }
    }
    const AdjointColorMatrix j_x = make_diff_exp_map((Complex)coef * x, *this);
    const AdjointColorMatrix j_n_x =
        make_diff_exp_map(-(Complex)coef * x, *this);
    qassert(qnorm(j_x - matrix_adjoint(j_n_x)) < 1e-20);
    qassert(qnorm(j_n_x - exp_adx * j_x) < 1e-20);
    for (int a = 0; a < 8; ++a) {
      const AdjointColorMatrix am0 = make_diff_exp_map_diff_ref(x, a, *this);
      const AdjointColorMatrix am1 = make_diff_exp_map_diff(x, a, *this);
      qassert(qnorm(am0 - am1) < 1e-8);
    }
  }
  //
  API static const box_acc<ColorMatrixConstants>& get_instance_box()
  {
    static box_acc<ColorMatrixConstants> cmcs =
        box_acc<ColorMatrixConstants>(ColorMatrixConstants());
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
  API static const AdjointColorMatrix& get_aunit() { return get_instance().aunit; }
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
  const array<double, 8> basis = basis_projection_anti_hermitian_matrix(m);
  AdjointColorMatrix am;
  set_zero(am);
  for (int a = 0; a < 8; ++a) {
    am += (-basis[a]) * f[a];
  }
  return am;
}

qacc AdjointColorMatrix make_diff_exp_map(const AdjointColorMatrix& am,
                                          const int max_order)
{
  AdjointColorMatrix aunit;
  set_unit(aunit);
  AdjointColorMatrix t2 = -am;
  AdjointColorMatrix t3;
  for (int j = max_order; j > 1; --j) {
    t3 = aunit + (1.0 / (j + 1)) * t2;
    t2 = -am * t3;
  }
  t3 = aunit + 0.5 * t2;
  return t3;
}

qacc AdjointColorMatrix make_diff_exp_map(const ColorMatrix& m,
                                          const ColorMatrixConstants& cmcs,
                                          const int max_order)
{
  return make_diff_exp_map(make_adjoint_representation(m, cmcs), max_order);
}

qacc AdjointColorMatrix make_diff_exp_map_diff_ref(
    const ColorMatrix& m, const int a, const ColorMatrixConstants& cmcs)
{
  const array<ColorMatrix, 8>& ts = cmcs.ts;
  const double eps = 1e-4;
  const ColorMatrix m_p = m + (Complex)eps * ts[a];
  const ColorMatrix m_n = m - (Complex)eps * ts[a];
  const AdjointColorMatrix j_p = make_diff_exp_map(m_p, cmcs);
  const AdjointColorMatrix j_n = make_diff_exp_map(m_n, cmcs);
  return (1.0 / (2.0 * eps)) * (j_p - j_n);
}

qacc AdjointColorMatrix make_diff_exp_map_diff(const AdjointColorMatrix& am,
                                               const int a,
                                               const ColorMatrixConstants& cmcs,
                                               const int max_order)
{
  const array<AdjointColorMatrix, 8>& f = cmcs.f;
  const AdjointColorMatrix& aunit = cmcs.aunit;
  AdjointColorMatrix t2 = -am;
  AdjointColorMatrix dt2 = f[a];
  AdjointColorMatrix t3;
  AdjointColorMatrix dt3;
  for (int j = max_order; j > 1; --j) {
    t3 = aunit + (1.0 / (j + 1)) * t2;
    dt3 = (1.0 / (j + 1)) * dt2;
    t2 = -am * t3;
    dt2 = f[a] * t3 - am * dt3;
  }
  dt3 = 0.5 * dt2;
  return dt3;
}

qacc AdjointColorMatrix make_diff_exp_map_diff(const ColorMatrix& m,
                                               const int a,
                                               const ColorMatrixConstants& cmcs,
                                               const int max_order)
{
  return make_diff_exp_map_diff(make_adjoint_representation(m, cmcs), a, cmcs,
                                max_order);
}

}  // namespace qlat
