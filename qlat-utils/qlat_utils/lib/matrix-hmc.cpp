#include <qlat-utils/matrix-hmc.h>

namespace qlat
{  //

void ColorMatrixConstants::check() const
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
      ComplexD tr(matrix_trace(ts[a], ts[b]));
      if (a == b) {
        tr -= -2;
      }
      if (qnorm(tr) >= 1e-20) {
        displayln_info(ssprintf("a=%d b=%d tr-dif=%s", a, b, show(tr).c_str()));
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
      ComplexD sum(0.0);
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
    const ColorMatrix exp_x = make_matrix_exp((ComplexD)coef * x);
    const ColorMatrix exp_n_x = make_matrix_exp(-(ComplexD)coef * x);
    if (qnorm(mx_b - exp_x * ts[b] * exp_n_x) >= 1e-20) {
      displayln_info("adx_b: " + show(qnorm(mx_b)));
      displayln_info("exp_x b exp_n_x: " +
                     show(qnorm(exp_x * ts[b] * exp_n_x)));
      displayln_info("diff: " + show(qnorm(mx_b - exp_x * ts[b] * exp_n_x)));
      qassert(false);
    }
  }
  const AdjointColorMatrix j_x = make_diff_exp_map((ComplexD)coef * x, *this);
  const AdjointColorMatrix j_n_x =
      make_diff_exp_map(-(ComplexD)coef * x, *this);
  qassert(qnorm(j_x - matrix_adjoint(j_n_x)) < 1e-20);
  qassert(qnorm(j_n_x - exp_adx * j_x) < 1e-20);
  for (int a = 0; a < 8; ++a) {
    const AdjointColorMatrix am0 = make_diff_exp_map_diff_ref(x, a, *this);
    const AdjointColorMatrix am1 = make_diff_exp_map_diff(x, a, *this);
    qassert(qnorm(am0 - am1) < 1e-8);
  }
}

}  // namespace qlat
