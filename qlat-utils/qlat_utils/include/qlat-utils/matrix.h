#pragma once

#include <qlat-utils/mat-vec.h>
#include <qlat-utils/mvector.h>
#include <qlat-utils/vector.h>

#include <cmath>

namespace qlat
{  //

template <int DIMN, class T>
qacc MatrixT<DIMN, T> operator+(const MatrixT<DIMN, T>& x,
                                const MatrixT<DIMN, T>& y)
{
  MatrixT<DIMN, T> ret;
  for (int i = 0; i < DIMN * DIMN; ++i) {
    ret.p[i] = x.p[i] + y.p[i];
  }
  return ret;
}

template <int DIMN, class T>
qacc MatrixT<DIMN, T> operator-(const MatrixT<DIMN, T>& x,
                                const MatrixT<DIMN, T>& y)
{
  MatrixT<DIMN, T> ret;
  for (int i = 0; i < DIMN * DIMN; ++i) {
    ret.p[i] = x.p[i] - y.p[i];
  }
  return ret;
}

template <int DIMN, class T>
qacc MatrixT<DIMN, T> operator-(const MatrixT<DIMN, T>& x)
{
  MatrixT<DIMN, T> ret;
  for (int i = 0; i < DIMN * DIMN; ++i) {
    ret.p[i] = -x.p[i];
  }
  return ret;
}

template <int DIMN, class T>
qacc MatrixT<DIMN, T> operator*(const MatrixT<DIMN, T>& x,
                                const MatrixT<DIMN, T>& y)
{
  MatrixT<DIMN, T> ret;
  MatrixT<DIMN, T> yt = matrix_transpose(y);
  for (int i = 0; i < DIMN; ++i) {
    const int id = i * DIMN;
    for (int j = 0; j < DIMN; ++j) {
      const int jd = j * DIMN;
      T s = 0.0;
      for (int k = 0; k < DIMN; ++k) {
        s += x.p[id + k] * yt.p[jd + k];
      }
      ret.p[id + j] = s;
    }
  }
  return ret;
}

template <int DIMN, class T>
qacc MatrixT<DIMN, T> operator*(const double x, const MatrixT<DIMN, T>& y)
{
  MatrixT<DIMN, T> ret;
  for (int i = 0; i < DIMN * DIMN; ++i) {
    ret.p[i] = x * y.p[i];
  }
  return ret;
}

template <int DIMN, class T>
qacc MatrixT<DIMN, T> operator*(const MatrixT<DIMN, T>& x, const double y)
{
  MatrixT<DIMN, T> ret;
  for (int i = 0; i < DIMN * DIMN; ++i) {
    ret.p[i] = x.p[i] * y;
  }
  return ret;
}

template <int DIMN, class T>
qacc MatrixT<DIMN, T> operator*(const ComplexD& x, const MatrixT<DIMN, T>& y)
{
  MatrixT<DIMN, T> ret;
  for (int i = 0; i < DIMN * DIMN; ++i) {
    ret.p[i] = x * y.p[i];
  }
  return ret;
}

template <int DIMN, class T>
qacc MatrixT<DIMN, T> operator*(const MatrixT<DIMN, T>& x, const ComplexD& y)
{
  MatrixT<DIMN, T> ret;
  for (int i = 0; i < DIMN * DIMN; ++i) {
    ret.p[i] = x.p[i] * y;
  }
  return ret;
}

template <int DIMN, class T>
qacc MatrixT<DIMN, T> operator/(const MatrixT<DIMN, T>& x, const T& y)
{
  MatrixT<DIMN, T> ret;
  for (int i = 0; i < DIMN * DIMN; ++i) {
    ret.p[i] = x.p[i] / y;
  }
  return ret;
}

template <int DIMN, class T>
qacc MvectorT<DIMN, T> operator*(const MatrixT<DIMN, T>& x,
                                 const MvectorT<DIMN, T>& y)
{
  MvectorT<DIMN, T> ret;
  for (int i = 0; i < DIMN; ++i) {
    const int id = i * DIMN;
    T s = 0.0;
    for (int j = 0; j < DIMN; ++j) {
      s += x.p[id + j] * y.p[j];
    }
    ret.p[i] = s;
  }
  return ret;
}

template <int DIMN, class T>
qacc void set_zero(MatrixT<DIMN, T>& m)
{
  memset((void*)&m, 0, sizeof(MatrixT<DIMN, T>));
}

template <int DIMN>
qacc void set_unit(MatrixT<DIMN, RealF>& m, const ComplexD& coef = 1.0)
{
  set_zero(m);
  for (int i = 0; i < DIMN; ++i) {
    m(i, i) = coef.real();
  }
}

template <int DIMN>
qacc void set_unit(MatrixT<DIMN, RealD>& m, const ComplexD& coef = 1.0)
{
  set_zero(m);
  for (int i = 0; i < DIMN; ++i) {
    m(i, i) = coef.real();
  }
}

template <int DIMN, class T>
qacc void set_unit(MatrixT<DIMN, T>& m, const ComplexD& coef = 1.0)
{
  set_zero(m);
  for (int i = 0; i < DIMN; ++i) {
    m(i, i) = coef;
  }
}

template <int DIMN, class T>
qacc RealD qnorm(const MatrixT<DIMN, T>& m)
{
  RealD s = 0.0;
  for (int i = 0; i < DIMN * DIMN; ++i) {
    s += qnorm(m.p[i]);
  }
  return s;
}

template <int DIMN, class T>
qacc RealD qnorm(const MatrixT<DIMN, T>& m1, const MatrixT<DIMN, T>& m2)
{
  RealD s = 0.0;
  for (int i = 0; i < DIMN * DIMN; ++i) {
    s += qnorm(m1.p[i], m2.p[i]);
  }
  return s;
}

template <int DIMN, class T>
qacc ComplexD matrix_trace(const MatrixT<DIMN, T>& x)
{
  ComplexD s = 0.0;
  for (int i = 0; i < DIMN; ++i) {
    s += x.p[i * (DIMN + 1)];
  }
  return s;
}

template <int DIMN, class T>
qacc ComplexD matrix_trace(const MatrixT<DIMN, T>& x, const MatrixT<DIMN, T>& y)
{
  MatrixT<DIMN, T> yt = matrix_transpose(y);
  ComplexD s = 0.0;
  for (int i = 0; i < DIMN; ++i) {
    const int id = i * DIMN;
    for (int j = 0; j < DIMN; ++j) {
      s += x.p[id + j] * yt.p[id + j];
    }
  }
  return s;
}

template <int DIMN, class T>
qacc MatrixT<DIMN, T> matrix_adjoint(const MatrixT<DIMN, T>& x)
{
  MatrixT<DIMN, T> ret;
  for (int i = 0; i < DIMN; ++i) {
    const int id = i * DIMN;
    for (int j = 0; j < DIMN; ++j) {
      const int jd = j * DIMN;
      ret.p[jd + i] = qconj(x.p[id + j]);
    }
  }
  return ret;
}

template <int DIMN, class T>
qacc MatrixT<DIMN, T> matrix_transpose(const MatrixT<DIMN, T>& x)
{
  MatrixT<DIMN, T> ret;
  for (int i = 0; i < DIMN; ++i) {
    const int id = i * DIMN;
    for (int j = 0; j < DIMN; ++j) {
      const int jd = j * DIMN;
      ret.p[jd + i] = x.p[id + j];
    }
  }
  return ret;
}

template <int DIMN, class T>
qacc MatrixT<DIMN, T> matrix_conjugate(const MatrixT<DIMN, T>& x)
{
  MatrixT<DIMN, T> ret;
  for (int i = 0; i < DIMN * DIMN; ++i) {
    ret.p[i] = qconj(x.p[i]);
  }
  return ret;
}

template <int DIMN, class T>
qacc void matrix_ludcmp(MatrixT<DIMN, T>& a,
                        array<int, (unsigned long)DIMN>& indx, T& d)
// https://www.astro.umd.edu/~ricotti/NEWWEB/teaching/ASTR415/InClassExamples/NR3/legacy/nr2/C_211/progs.htm
// https://www.astro.umd.edu/~ricotti/NEWWEB/teaching/ASTR415/InClassExamples/NR3/legacy/nr2/C_211/recipes/ludcmp.c
//
// a: input and output
// indx, d: output
{
  MvectorT<DIMN, T> vv;
  d = 1.0;
  for (int i = 0; i < DIMN; ++i) {
    double big = 0.0;
    for (int j = 0; j < DIMN; ++j) {
      const double temp = std::abs(a(i, j));
      if (temp > big) {
        big = temp;
      }
    }
    if (big == 0.0) {
      // matrix singular
      qassert(false);
    }
    vv(i) = 1.0 / big;
  }
  for (int j = 0; j < DIMN; ++j) {
    for (int i = 0; i < j; ++i) {
      T sum = a(i, j);
      for (int k = 0; k < i; ++k) {
        sum -= a(i, k) * a(k, j);
      }
      a(i, j) = sum;
    }
    double big = 0.0;
    int imax = -1;
    for (int i = j; i < DIMN; ++i) {
      T sum = a(i, j);
      for (int k = 0; k < j; ++k) {
        sum -= a(i, k) * a(k, j);
      }
      a(i, j) = sum;
      const T dum = vv(i) * std::abs(sum);
      if (dum >= big) {
        big = dum;
        imax = i;
      }
    }
    if (j != imax) {
      for (int k = 0; k < DIMN; ++k) {
        const T dum = a(imax, k);
        a(imax, k) = a(j, k);
        a(j, k) = dum;
      }
      d = -d;
      vv(imax) = vv(j);
    }
    indx[j] = imax;
    if (a(j, j) == 0.0) {
      // matrix singular
      qassert(false);
    }
    if (j != DIMN - 1) {
      const T dum = 1.0 / a(j, j);
      for (int i = j + 1; i < DIMN; ++i) {
        a(i, j) *= dum;
      }
    }
  }
}

template <int DIMN, class T>
qacc void matrix_lubksb(MvectorT<DIMN, T>& b, const MatrixT<DIMN, T>& a,
                        const array<int, (unsigned long)DIMN>& indx)
// https://www.astro.umd.edu/~ricotti/NEWWEB/teaching/ASTR415/InClassExamples/NR3/legacy/nr2/C_211/progs.htm
// https://www.astro.umd.edu/~ricotti/NEWWEB/teaching/ASTR415/InClassExamples/NR3/legacy/nr2/C_211/recipes/lubksb.c
//
// b: input and output
{
  int ii = -1;
  for (int i = 0; i < DIMN; ++i) {
    const int ip = indx[i];
    T sum = b(ip);
    b(ip) = b(i);
    if (ii >= 0) {
      for (int j = ii; j < i; j++) {
        sum -= a(i, j) * b(j);
      }
    } else if (sum != 0.0) {
      ii = i;
    }
    b(i) = sum;
  }
  for (int i = DIMN - 1; i >= 0; --i) {
    T sum = b(i);
    for (int j = i + 1; j < DIMN; ++j) {
      sum -= a(i, j) * b(j);
    }
    b(i) = sum / a(i, i);
  }
}

template <int DIMN, class T>
qacc MatrixT<DIMN, T> matrix_inverse(const MatrixT<DIMN, T>& x)
{
  MatrixT<DIMN, T> ret;
  MatrixT<DIMN, T> a = x;
  array<int, DIMN> indx;
  T d;
  matrix_ludcmp(a, indx, d);
  MvectorT<DIMN, T> col;
  for (int j = 0; j < DIMN; ++j) {
    set_zero(col);
    col(j) = 1.0;
    matrix_lubksb(col, a, indx);
    for (int i = 0; i < DIMN; ++i) {
      ret(i, j) = col(i);
    }
  }
  return ret;
}

template <int DIMN, class T>
qacc T matrix_determinant(const MatrixT<DIMN, T>& x)
{
  MatrixT<DIMN, T> a = x;
  array<int, DIMN> indx;
  T d;
  matrix_ludcmp(a, indx, d);
  for (int j = 0; j < DIMN; ++j) {
    d *= a(j, j);
  }
  return d;
}

template <class T = Real>
struct API SpinMatrixConstantsT {
  SpinMatrixT<T> unit;
  array<SpinMatrixT<T>, 4>
      gammas;  // Not using CPS's convention, but a more standard one.
  array<SpinMatrixT<T>, 4> cps_gammas;  // CPS's convention gamma matrices
  SpinMatrixT<T> gamma5;                // Same as CPS's gamma5
  array<SpinMatrixT<T>, 3> cap_sigmas;
  array<SpinMatrixT<T>, 3> cps_cap_sigmas;  // CPS's convention sigmas
  array<SpinMatrixT<T>, 16> gms;
  array<SpinMatrixT<T>, 16> cps_gms;
  //
  qacc SpinMatrixConstantsT() { init(); }
  //
  qacc void init()
  {
    ComplexD ii(0.0, 1.0);
    // TIMER_VERBOSE("SpinMatrixConstants::init()");
    unit << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;
    // gamma_x
    gammas[0] << 0, 0, 0, 1, 0, 0, 1, 0, 0, -1, 0, 0, -1, 0, 0, 0;
    // gamma_y
    gammas[1] << 0, 0, 0, (ComplexT<T>)-ii, 0, 0, (ComplexT<T>)ii, 0, 0,
        (ComplexT<T>)ii, 0, 0, (ComplexT<T>)-ii, 0, 0, 0;
    // gamma_z
    gammas[2] << 0, 0, 1, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 1, 0, 0;
    gammas[0] *= (ComplexT<T>)(-ii);
    gammas[1] *= (ComplexT<T>)(-ii);
    gammas[2] *= (ComplexT<T>)(-ii);
    // gamma_t
    gammas[3] << 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0;
    //
    cps_gammas[0] = -gammas[0];
    cps_gammas[1] = gammas[1];
    cps_gammas[2] = -gammas[2];
    cps_gammas[3] = gammas[3];
    // gamma_5 = gamma_x * gamma_y * gamma_z * gamma_t;
    gamma5 << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1;
    // Sigma_x
    cap_sigmas[0] << 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0;
    // Sigma_y
    cap_sigmas[1] << 0, (ComplexT<T>)-ii, 0, 0, (ComplexT<T>)ii, 0, 0, 0, 0, 0,
        0, (ComplexT<T>)-ii, 0, 0, (ComplexT<T>)ii, 0;
    // Sigma_z
    cap_sigmas[2] << 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1;
    //
    cps_cap_sigmas[0] = -cap_sigmas[0];
    cps_cap_sigmas[1] = cap_sigmas[1];
    cps_cap_sigmas[2] = -cap_sigmas[2];
    //
    for (int a = 0; a < 2; ++a) {
      const SpinMatrixT<T> ma = a == 0 ? unit : gammas[0];
      for (int b = 0; b < 2; ++b) {
        const SpinMatrixT<T> mb =
            b == 0 ? ma : (SpinMatrixT<T>)(ma * gammas[1]);
        for (int c = 0; c < 2; ++c) {
          const SpinMatrixT<T> mc =
              c == 0 ? mb : (SpinMatrixT<T>)(mb * gammas[2]);
          for (int d = 0; d < 2; ++d) {
            const SpinMatrixT<T> md =
                d == 0 ? mc : (SpinMatrixT<T>)(mc * gammas[3]);
            const int idx = a + 2 * b + 4 * c + 8 * d;
            gms[idx] = md;
          }
        }
      }
    }
    for (int a = 0; a < 2; ++a) {
      const SpinMatrixT<T> ma = a == 0 ? unit : cps_gammas[0];
      for (int b = 0; b < 2; ++b) {
        const SpinMatrixT<T> mb =
            b == 0 ? ma : (SpinMatrixT<T>)(ma * cps_gammas[1]);
        for (int c = 0; c < 2; ++c) {
          const SpinMatrixT<T> mc =
              c == 0 ? mb : (SpinMatrixT<T>)(mb * cps_gammas[2]);
          for (int d = 0; d < 2; ++d) {
            const SpinMatrixT<T> md =
                d == 0 ? mc : (SpinMatrixT<T>)(mc * cps_gammas[3]);
            const int idx = a + 2 * b + 4 * c + 8 * d;
            cps_gms[idx] = md;
          }
        }
      }
    }
  }
  //
  API static const box<SpinMatrixConstantsT<T>>& get_instance_box()
  {
    static box<SpinMatrixConstantsT<T>> smcs =
        box<SpinMatrixConstantsT<T>>(SpinMatrixConstantsT<T>(), MemType::Uvm);
    return smcs;
  }
  //
  API static const SpinMatrixConstantsT<T>& get_instance()
  {
    return get_instance_box()();
  }
  //
  API static const SpinMatrixT<T>& get_unit() { return get_instance().unit; }
  API static const SpinMatrixT<T>& get_gamma(int mu)
  {
    if (0 <= mu and mu < 4) {
      const SpinMatrix& gamma_mu = get_instance().gammas[mu];
      return gamma_mu;
    } else {
      qassert(mu == 5);
      const SpinMatrix& gamma5 = get_instance().gamma5;
      return gamma5;
    }
  }
  API static const SpinMatrixT<T>& get_cps_gamma(int mu)
  {
    if (0 <= mu and mu < 4) {
      const SpinMatrix& gamma_mu = get_instance().cps_gammas[mu];
      return gamma_mu;
    } else {
      qassert(mu == 5);
      const SpinMatrix& gamma5 = get_instance().gamma5;
      return gamma5;
    }
  }
  API static const array<SpinMatrixT<T>, 4>& get_gammas()
  {
    return get_instance().gammas;
  }
  API static const array<SpinMatrixT<T>, 4>& get_cps_gammas()
  {
    return get_instance().cps_gammas;
  }
  API static const array<SpinMatrixT<T>, 16>& get_gms()
  {
    return get_instance().gms;
  }
  API static const array<SpinMatrixT<T>, 16>& get_cps_gms()
  {
    return get_instance().cps_gms;
  }
  API static const SpinMatrixT<T>& get_gamma5()
  {
    return get_instance().gamma5;
  }
  API static const SpinMatrixT<T>& get_cap_sigma(int i)
  {
    return get_instance().cap_sigmas[i];
  }
  API static const array<SpinMatrixT<T>, 3>& get_cap_sigmas()
  {
    return get_instance().cap_sigmas;
  }
  API static const array<SpinMatrixT<T>, 3>& get_cps_cap_sigmas()
  {
    return get_instance().cps_cap_sigmas;
  }
  API static const SpinMatrixT<T> get_cps_sigmas(int i, int j)
  {
    const SpinMatrixConstantsT<T>& instance = get_instance();
    return (instance.cps_gammas[i] * instance.cps_gammas[j] -
            instance.cps_gammas[j] * instance.cps_gammas[i]) /
           2.;
  }
};

inline const box<SpinMatrixConstantsT<> >& get_spin_matrix_constants()
{
  return SpinMatrixConstantsT<>::get_instance_box();
}

template <class T>
qacc WilsonMatrixT<T> operator*(const ColorMatrixT<T>& cm,
                                const WilsonMatrixT<T>& m)
{
  WilsonMatrixT<T> ret;
  set_zero(ret);
  for (int s1 = 0; s1 < 4; ++s1) {
    for (int s2 = 0; s2 < 4; ++s2) {
      for (int c1 = 0; c1 < NUM_COLOR; ++c1) {
        for (int c2 = 0; c2 < NUM_COLOR; ++c2) {
          const int s1_s2_c2 =
              s1 * (NUM_COLOR * 4 * NUM_COLOR) + s2 * NUM_COLOR + c2;
          for (int c3 = 0; c3 < NUM_COLOR; ++c3) {
            ret.p[s1_s2_c2 + c1 * (4 * NUM_COLOR)] +=
                cm.p[c1 * NUM_COLOR + c3] *
                m.p[s1_s2_c2 + c3 * (4 * NUM_COLOR)];
          }
        }
      }
    }
  }
  return ret;
}

template <class T>
qacc WilsonMatrixT<T> operator*(const WilsonMatrixT<T>& m,
                                const ColorMatrixT<T>& cm)
{
  WilsonMatrixT<T> ret;
  set_zero(ret);
  for (int s1 = 0; s1 < 4; ++s1) {
    for (int s2 = 0; s2 < 4; ++s2) {
      for (int c1 = 0; c1 < NUM_COLOR; ++c1) {
        for (int c2 = 0; c2 < NUM_COLOR; ++c2) {
          const int s1_c1_s2 =
              (s1 * NUM_COLOR + c1) * (4 * NUM_COLOR) + s2 * NUM_COLOR;
          for (int c3 = 0; c3 < NUM_COLOR; ++c3) {
            ret.p[s1_c1_s2 + c2] +=
                m.p[s1_c1_s2 + c3] * cm.p[c3 * NUM_COLOR + c2];
          }
        }
      }
    }
  }
  return ret;
}

template <class T>
qacc WilsonMatrixT<T> operator*(const SpinMatrixT<T>& sm,
                                const WilsonMatrixT<T>& m)
{
  WilsonMatrixT<T> ret;
  set_zero(ret);
  for (int s1 = 0; s1 < 4; ++s1) {
    for (int s3 = 0; s3 < 4; ++s3) {
      const ComplexT<T>& sm_s1_s3 = sm.p[s1 * 4 + s3];
      if (sm_s1_s3 == 0.0) {
        continue;
      }
      for (int s2 = 0; s2 < 4; ++s2) {
        for (int c1 = 0; c1 < NUM_COLOR; ++c1) {
          const int s1_c1_s2 =
              (s1 * NUM_COLOR + c1) * (4 * NUM_COLOR) + s2 * NUM_COLOR;
          const int s3_c1_s2 =
              (s3 * NUM_COLOR + c1) * (4 * NUM_COLOR) + s2 * NUM_COLOR;
          for (int c2 = 0; c2 < NUM_COLOR; ++c2) {
            ret.p[s1_c1_s2 + c2] += sm_s1_s3 * m.p[s3_c1_s2 + c2];
          }
        }
      }
    }
  }
  return ret;
}

template <class T>
qacc WilsonMatrixT<T> operator*(const WilsonMatrixT<T>& m,
                                const SpinMatrixT<T>& sm)
{
  WilsonMatrixT<T> ret;
  set_zero(ret);
  for (int s1 = 0; s1 < 4; ++s1) {
    for (int s3 = 0; s3 < 4; ++s3) {
      const ComplexT<T>& sm_s3_s1 = sm.p[s3 * 4 + s1];
      if (sm_s3_s1 == 0.0) {
        continue;
      }
      for (int s2 = 0; s2 < 4; ++s2) {
        for (int c2 = 0; c2 < NUM_COLOR; ++c2) {
          const int s2_c2_s1 =
              (s2 * NUM_COLOR + c2) * (4 * NUM_COLOR) + s1 * NUM_COLOR;
          const int s2_c2_s3 =
              (s2 * NUM_COLOR + c2) * (4 * NUM_COLOR) + s3 * NUM_COLOR;
          for (int c1 = 0; c1 < NUM_COLOR; ++c1) {
            ret.p[s2_c2_s1 + c1] += m.p[s2_c2_s3 + c1] * sm_s3_s1;
          }
        }
      }
    }
  }
  return ret;
}

template <class T>
qacc ComplexD matrix_trace(const SpinMatrixT<T>& sm, const WilsonMatrixT<T>& m)
{
  ComplexD ret = 0;
  for (int s1 = 0; s1 < 4; ++s1) {
    for (int s3 = 0; s3 < 4; ++s3) {
      const ComplexT<T>& sm_s1_s3 = sm.p[s1 * 4 + s3];
      if (sm_s1_s3 == 0.0) {
        continue;
      }
      const int s3_s1 = s3 * (NUM_COLOR * 4 * NUM_COLOR) + s1 * NUM_COLOR;
      for (int c1 = 0; c1 < NUM_COLOR; ++c1) {
        ret += sm_s1_s3 * m.p[s3_s1 + c1 * (4 * NUM_COLOR + 1)];
      }
    }
  }
  return ret;
}

template <class T>
qacc ComplexD matrix_trace(const WilsonMatrixT<T>& m, const SpinMatrixT<T>& sm)
{
  return matrix_trace(sm, m);
}

template <class T>
qacc ComplexD matrix_trace(const ColorMatrixT<T>& cm, const WilsonMatrixT<T>& m)
{
  ComplexD ret = 0;
  for (int c1 = 0; c1 < NUM_COLOR; ++c1) {
    for (int c3 = 0; c3 < NUM_COLOR; ++c3) {
      const ComplexT<T>& cm_c1_c3 = cm.p[c1 * 3 + c3];
      if (cm_c1_c3 == 0.0) {
        continue;
      }
      const int c3_c1 = c3 * (4 * NUM_COLOR) + c1;
      for (int s1 = 0; s1 < 4; ++s1) {
        ret += cm_c1_c3 * m.p[c3_c1 + s1 * (NUM_COLOR * 4 * NUM_COLOR + NUM_COLOR)];
      }
    }
  }
  return ret;
}

template <class T>
qacc ComplexD matrix_trace(const WilsonMatrixT<T>& m, const ColorMatrixT<T>& cm)
{
  return matrix_trace(cm, m);
}

template <class T>
qacc WilsonVectorT<T> operator*(const ColorMatrixT<T>& cm,
                                const WilsonVectorT<T>& m)
{
  WilsonVectorT<T> ret;
  set_zero(ret);
  for (int s1 = 0; s1 < 4; ++s1) {
    for (int c1 = 0; c1 < NUM_COLOR; ++c1) {
      for (int c2 = 0; c2 < NUM_COLOR; ++c2) {
        ret.p[s1 * NUM_COLOR + c1] +=
            cm.p[c1 * NUM_COLOR + c2] * m.p[s1 * NUM_COLOR + c2];
      }
    }
  }
  return ret;
}

template <class T>
qacc WilsonVectorT<T> operator*(const SpinMatrixT<T>& sm,
                                const WilsonVectorT<T>& m)
{
  WilsonVectorT<T> ret;
  set_zero(ret);
  for (int s1 = 0; s1 < 4; ++s1) {
    for (int s2 = 0; s2 < 4; ++s2) {
      const ComplexT<T>& sm_s1_s2 = sm.p[s1 * 4 + s2];
      if (sm_s1_s2 == 0.0) {
        continue;
      }
      for (int c1 = 0; c1 < NUM_COLOR; ++c1) {
        ret.p[s1 * NUM_COLOR + c1] +=
            sm.p[s1 * 4 + s2] * m.p[s2 * NUM_COLOR + c1];
      }
    }
  }
  return ret;
}

template <class T>
qacc void convert_mspincolor_from_wm(WilsonMatrixT<T>& msc,
                                     const WilsonMatrixT<T>& wm)
{
  for (int s1 = 0; s1 < 4; ++s1) {
    for (int s2 = 0; s2 < 4; ++s2) {
      for (int c1 = 0; c1 < NUM_COLOR; ++c1) {
        for (int c2 = 0; c2 < NUM_COLOR; ++c2) {
          msc.p[(s1 * 4 + s2) * NUM_COLOR * NUM_COLOR + c1 * NUM_COLOR + c2] =
              wm(s1 * NUM_COLOR + c1, s2 * NUM_COLOR + c2);
        }
      }
    }
  }
}

template <class T>
qacc void convert_wm_from_mspincolor(WilsonMatrixT<T>& wm,
                                     const WilsonMatrixT<T>& msc)
{
  for (int s1 = 0; s1 < 4; ++s1) {
    for (int s2 = 0; s2 < 4; ++s2) {
      for (int c1 = 0; c1 < NUM_COLOR; ++c1) {
        for (int c2 = 0; c2 < NUM_COLOR; ++c2) {
          wm(s1 * NUM_COLOR + c1, s2 * NUM_COLOR + c2) =
              msc.p[(s1 * 4 + s2) * NUM_COLOR * NUM_COLOR + c1 * NUM_COLOR +
                    c2];
        }
      }
    }
  }
}

template <class T>
qacc ComplexT<T> epsilon_contraction(const int v_s1, const int b_s1, const int v_s2,
                           const int b_s2, const int v_s3, const int b_s3,
                           const WilsonMatrixT<T>& wm1,
                           const WilsonMatrixT<T>& wm2,
                           const WilsonMatrixT<T>& wm3)
{
  array<array<int, 4>, 6> eps_table;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      eps_table[i][j] = (i + j) % 3;
      eps_table[i + 3][j] = (i + 3 - j) % 3;
    }
    eps_table[i][3] = 1;
    eps_table[i + 3][3] = -1;
  }
  const int n_color = NUM_COLOR;
  const int n_spin_color = 4 * NUM_COLOR;
  const int n_color_spin_color = NUM_COLOR * 4 * NUM_COLOR;
  ComplexT<T> sum = 0;
  for (int v = 0; v < 6; ++v) {
    const int v_c1 = eps_table[v][0];
    const int v_c2 = eps_table[v][1];
    const int v_c3 = eps_table[v][2];
    const int v_sign = eps_table[v][3];
    for (int b = 0; b < 6; ++b) {
      const int b_c1 = eps_table[b][0];
      const int b_c2 = eps_table[b][1];
      const int b_c3 = eps_table[b][2];
      const int b_sign = eps_table[b][3];
      const int sign = v_sign * b_sign;
      const ComplexT<T>& val1 =
          wm1.p[v_s1 * n_color_spin_color + v_c1 * n_spin_color +
                b_s1 * n_color + b_c1];
      const ComplexT<T>& val2 =
          wm2.p[v_s2 * n_color_spin_color + v_c2 * n_spin_color +
                b_s2 * n_color + b_c2];
      const ComplexT<T>& val3 =
          wm3.p[v_s3 * n_color_spin_color + v_c3 * n_spin_color +
                b_s3 * n_color + b_c3];
      sum += sign * val1 * val2 * val3;
    }
  }
  return sum;
}

#ifndef QLAT_NO_DEFAULT_TYPE

typedef SpinMatrixConstantsT<> SpinMatrixConstants;

#endif

template <int DIMN, class T>
std::string show(const MatrixT<DIMN, T>& m)
{
  std::ostringstream out;
  out << "[ ";
  for (int i = 0; i < DIMN; ++i) {
    const int id = i * DIMN;
    out << "[ ";
    for (int j = 0; j < DIMN; ++j) {
      out << m.p[id + j] << ", ";
    }
    out << "], ";
  }
  out << "]";
  return out.str();
}

}  // namespace qlat
