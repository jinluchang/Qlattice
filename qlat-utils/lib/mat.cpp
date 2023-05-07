#include <qlat-utils/mat.h>
#include <qlat-utils/matrix.h>
#include <qlat-utils/vector.h>
#include <qlat-utils/rng-state.h>

namespace qlat
{  //

const SpinMatrixT<>& get_gamma_matrix(const int mu)
// CPS's convention gamma matrices
{
  return SpinMatrixConstantsT<>::get_cps_gamma(mu);
}

template <class M>
void set_rand_complex_mat(M& m, RngState& rs)
{
  TIMER("set_rand_complex_mat")
  Vector<Complex> v = get_data(m);
  for (long k = 0; k < v.size(); ++k) {
    v[k].real(u_rand_gen(rs));
    v[k].imag(u_rand_gen(rs));
  }
}

void benchmark_matrix_functions(const long count)
{
  TIMER_VERBOSE("benchmark_matrix_functions");
  const long n_mat = 8;
  vector_acc<WilsonMatrix> wm(n_mat);
  vector_acc<SpinMatrix> sm(n_mat);
  for (long i = 0; i < n_mat; ++i) {
    RngState rs(fname + ssprintf(": %d/%d", i, n_mat));
    set_rand_complex_mat(wm[i], rs);
    set_rand_complex_mat(sm[i], rs);
  }
  {
    TIMER_VERBOSE_FLOPS("benchmark_matrix_functions-wm*wm");
    for (long i = 0; i < count; ++i) {
      const WilsonMatrix& wm1 = wm[(i * 2) % n_mat];
      const WilsonMatrix& wm2 = wm[(i * 2 + 2) % n_mat];
      WilsonMatrix& wm3 = wm[(i * 2 + 1) % n_mat];
      wm3 = wm1 * wm2;
    }
    timer.flops += count * 13536;
  }
  {
    TIMER_VERBOSE_FLOPS(fname + "benchmark_matrix_functions-sm*wm");
    for (long i = 0; i < count; ++i) {
      const SpinMatrix& sm1 = sm[(i * 2) % n_mat];
      const WilsonMatrix& wm2 = wm[(i * 2 + 2) % n_mat];
      WilsonMatrix& wm3 = wm[(i *2 + 1) % n_mat];
      wm3 = sm1 * wm2;
    }
    timer.flops += count * 4320;
  }
  {
    TIMER_VERBOSE_FLOPS(fname + "benchmark_matrix_functions-wm*sm");
    for (long i = 0; i < count; ++i) {
      const WilsonMatrix& wm1 = wm[(i * 2) % n_mat];
      const SpinMatrix& sm2 = sm[(i * 2 + 2) % n_mat];
      WilsonMatrix& wm3 = wm[(i * 2 + 1) % n_mat];
      wm3 = wm1 * sm2;
    }
    timer.flops += count * 4320;
  }
  {
    TIMER_VERBOSE_FLOPS(fname + "benchmark_matrix_functions-sm*sm");
    for (long i = 0; i < count; ++i) {
      const SpinMatrix& sm1 = sm[(i * 2) % n_mat];
      const SpinMatrix& sm2 = sm[(i * 2 + 2) % n_mat];
      SpinMatrix& sm3 = sm[(i * 2 + 1) % n_mat];
      sm3 = sm1 * sm2;
    }
    timer.flops += count * 480;
  }
}

}  // namespace qlat
