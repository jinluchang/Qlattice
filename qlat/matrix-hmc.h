#pragma once

#include <qlat/matrix.h>

namespace qlat
{  //

template <class T>
void unitarize(ColorMatrixT<T>& cm)
{
  cm.em().row(0).normalize();
  cm.em().row(2) =
      cm.em().row(1) - cm.em().row(0).dot(cm.em().row(1)) * cm.em().row(0);
  cm.em().row(1) = cm.em().row(2).normalized();
  cm.em().row(2) = cm.em().row(0).cross(cm.em().row(1));
}

inline ColorMatrixT<> make_anti_hermitian_matrix(
    const std::array<double, 8>& basis)
{
  qassert(3 == NUM_COLOR);
  ColorMatrixT<> m;
  Array<double, 18> p(m.d());
  const double s3 = 1.0 / std::sqrt(3.0) * basis[7];
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

inline std::array<double, 8> basis_projection_anti_hermitian_matrix(
    const ColorMatrix& m)
// m is a tr_less_anti_hermitian_matrix
{
  const Array<double, 18> p(m.d());
  std::array<double, 8> basis;
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
  std::array<double, 8> a;
  for (int i = 0; i < 8; ++i) {
    a[i] = g_rand_gen(rs, 0.0, s);
  }
  return make_anti_hermitian_matrix(a);
}

inline double neg_half_tr_square(const ColorMatrixT<>& m)
{
  const Array<double, 18> p(m.d());
  return sqr(p[3]) + sqr(p[2]) + sqr(p[1] - p[9]) * 0.25 + sqr(p[5]) +
         sqr(p[4]) + sqr(p[11]) + sqr(p[10]) +
         sqr(p[17]) * 0.75;  // 0.75 = (sqrt(3)/2)^2
}

template <class M>
M make_matrix_exp(const M& a, const int max_order = 20)
{
  M unit;
  set_unit(unit);
  M t2 = a;
  M t3;
  for (int j = max_order; j > 1; --j) {
    t3 = unit + (Complex)(1.0 / j) * t2;
    t2 = a * t3;
  }
  t3 = unit + t2;
  return t3;
}

template <class T>
ColorMatrixT<T> make_color_matrix_exp_no_unitarize(const ColorMatrixT<T>& a)
{
  return make_matrix_exp(a);
}

template <class T>
ColorMatrixT<T> make_color_matrix_exp(const ColorMatrixT<T>& a)
{
  ColorMatrixT<T> ret = make_color_matrix_exp_no_unitarize(a);
  unitarize(ret);
  return ret;
}

template <class T>
ColorMatrixT<T> matrix_evolve(const ColorMatrixT<T>& gf_cm,
                              const ColorMatrixT<T>& gm_cm,
                              const double step_size)
{
  const ColorMatrixT<T> t = (T)step_size * gm_cm;
  return make_color_matrix_exp_no_unitarize(t) * gf_cm;
}

inline ColorMatrixT<> make_tr_less_anti_herm_matrix(const ColorMatrixT<>& m)
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

struct AdjointColorMatrix : MatrixT<8, ComplexT> {
  AdjointColorMatrix() {}
  AdjointColorMatrix(const MatrixT<8, ComplexT>& m) { *this = m; }
  //
  const AdjointColorMatrix& operator=(const MatrixT<8, ComplexT>& m)
  {
    *this = (const AdjointColorMatrix&)m;
    return *this;
  }
};

inline AdjointColorMatrix make_adjoint_representation(const ColorMatrix& cm);

inline AdjointColorMatrix make_diff_exp_map(const AdjointColorMatrix& am,
                                            const int max_order = 20);

inline AdjointColorMatrix make_diff_exp_map(const ColorMatrix& m,
                                            const int max_order = 20);

struct ColorMatrixConstants {
  //
  // Luscher's convention: https://arxiv.org/abs/0907.5491v3
  //
  // except the norm of T^a matrices: tr(T^a T^b) = -2 \delta^{a,b}
  //
  ColorMatrix unit;
  std::array<ColorMatrix, 8> ts;
  AdjointColorMatrix aunit;
  std::array<AdjointColorMatrix, 8> f;
  //
  ColorMatrixConstants() { init(); }
  //
  void init()
  {
    set_unit(unit);
    for (int a = 0; a < 8; ++a) {
      std::array<double, 8> basis;
      set_zero(basis);
      basis[a] = 1.0;
      ts[a] = make_anti_hermitian_matrix(basis);
    }
    set_unit(aunit);
    for (int a = 0; a < 8; ++a) {
      for (int b = 0; b < 8; ++b) {
        const ColorMatrix m = ts[a] * ts[b] - ts[b] * ts[a];
        std::array<double, 8> basis = basis_projection_anti_hermitian_matrix(m);
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
      std::array<double, 8> basis;
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
    AdjointColorMatrix adx = make_adjoint_representation(x);
    std::array<double, 8> basis = basis_projection_anti_hermitian_matrix(x);
    for (int a = 0; a < 8; ++a) {
      for (int b = 0; b < 8; ++b) {
        Complex sum = 0.0;
        for (int c = 0; c < 8; ++c) {
          sum += -f[c](a, b) * basis[c];
        }
        qassert(qnorm(adx(a, b) - sum) < 1e-20);
      }
    }
    const Complex coef = 0.2;
    const AdjointColorMatrix exp_adx = make_matrix_exp(coef * adx);
    for (int b = 0; b < 8; ++b) {
      std::array<double, 8> basis_b;
      for (int a = 0; a < 8; ++a) {
        basis_b[a] = exp_adx(a, b).real();
      }
      const ColorMatrix mx_b = make_anti_hermitian_matrix(basis_b);
      const ColorMatrix exp_x = make_matrix_exp(coef * x);
      const ColorMatrix exp_n_x = make_matrix_exp(-coef * x);
      if (qnorm(mx_b - exp_x * ts[b] * exp_n_x) >= 1e-20) {
        displayln_info("adx_b: " + show(qnorm(mx_b)));
        displayln_info("exp_x b exp_n_x: " +
                       show(qnorm(exp_x * ts[b] * exp_n_x)));
        displayln_info("diff: " + show(qnorm(mx_b - exp_x * ts[b] * exp_n_x)));
        qassert(false);
      }
    }
    const AdjointColorMatrix j_x = make_diff_exp_map(coef * x);
    const AdjointColorMatrix j_n_x = make_diff_exp_map(-coef * x);
    qassert(qnorm(j_x - matrix_adjoint(j_n_x)) < 1e-20);
    qassert(qnorm(j_n_x - exp_adx * j_x) < 1e-20);
  }
  //
  static const ColorMatrixConstants& get_instance()
  {
    static ColorMatrixConstants cmcs;
    return cmcs;
  }
  //
  static const ColorMatrix& get_unit() { return get_instance().unit; }
  //
  static const std::array<ColorMatrix, 8>& get_ts()
  {
    return get_instance().ts;
  }
  //
  static const AdjointColorMatrix& get_aunit() { return get_instance().aunit; }
  //
  static const std::array<AdjointColorMatrix, 8>& get_f()
  {
    return get_instance().f;
  }
};

inline AdjointColorMatrix make_adjoint_representation(const ColorMatrix& m)
{
  const std::array<AdjointColorMatrix, 8>& f = ColorMatrixConstants::get_f();
  const std::array<double, 8> basis = basis_projection_anti_hermitian_matrix(m);
  AdjointColorMatrix am;
  set_zero(am);
  for (int a = 0; a < 8; ++a) {
    am += (Complex)(-basis[a]) * f[a];
  }
  return am;
}

inline AdjointColorMatrix make_diff_exp_map(const AdjointColorMatrix& am,
                                            const int max_order)
{
  const AdjointColorMatrix& unit = ColorMatrixConstants::get_aunit();
  AdjointColorMatrix t2 = -am;
  AdjointColorMatrix t3;
  for (int j = max_order; j > 1; --j) {
    t3 = unit + (Complex)(1.0 / (j + 1)) * t2;
    t2 = -am * t3;
  }
  t3 = unit + (Complex)0.5 * t2;
  return t3;
}

inline AdjointColorMatrix make_diff_exp_map(const ColorMatrix& m,
                                            const int max_order)
{
  return make_diff_exp_map(make_adjoint_representation(m), max_order);
}

}  // namespace qlat
