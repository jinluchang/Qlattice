#pragma once

#include <qlat/config.h>

#include <cmath>

QLAT_START_NAMESPACE

template <int DIMN, class T = ComplexT>
struct MatrixT {
  T p[DIMN * DIMN];
  //
  // convert to double array
  double* d() { return (double*)p; }
  const double* d() const { return (const double*)p; }
  //
  // convert to Eigen MatrixT
  Eigen::Matrix<T, DIMN, DIMN, Eigen::RowMajor>& em()
  {
    return *((Eigen::Matrix<T, DIMN, DIMN, Eigen::RowMajor>*)this);
  }
  const Eigen::Matrix<T, DIMN, DIMN, Eigen::RowMajor>& em() const
  {
    return *((Eigen::Matrix<T, DIMN, DIMN, Eigen::RowMajor>*)this);
  }
  //
  T& operator()(int i, int j)
  {
    qassert(0 <= i && i < DIMN);
    qassert(0 <= j && j < DIMN);
    return p[i * DIMN + j];
  }
  const T& operator()(int i, int j) const
  {
    qassert(0 <= i && i < DIMN);
    qassert(0 <= j && j < DIMN);
    return p[i * DIMN + j];
  }
  //
  const MatrixT& operator+=(const MatrixT& x)
  {
    *this = *this + x;
    return *this;
  }
  //
  const MatrixT& operator-=(const MatrixT& x)
  {
    *this = *this - x;
    return *this;
  }
  //
  const MatrixT& operator*=(const MatrixT& x)
  {
    *this = *this * x;
    return *this;
  }
  //
  const MatrixT& operator*=(const T& x)
  {
    *this = *this * x;
    return *this;
  }
  //
  const MatrixT& operator/=(const T& x)
  {
    *this = *this / x;
    return *this;
  }
};

template <int DIMN, class T>
MatrixT<DIMN, T> operator+(const MatrixT<DIMN, T>& x, const MatrixT<DIMN, T>& y)
{
  MatrixT<DIMN, T> ret;
  ret.em() = x.em() + y.em();
  return ret;
}

template <int DIMN, class T>
MatrixT<DIMN, T> operator-(const MatrixT<DIMN, T>& x, const MatrixT<DIMN, T>& y)
{
  MatrixT<DIMN, T> ret;
  ret.em() = x.em() - y.em();
  return ret;
}

template <int DIMN, class T>
MatrixT<DIMN, T> operator-(const MatrixT<DIMN, T>& x)
{
  MatrixT<DIMN, T> ret;
  ret.em() = -x.em();
  return ret;
}

template <int DIMN, class T>
MatrixT<DIMN, T> operator*(const MatrixT<DIMN, T>& x, const MatrixT<DIMN, T>& y)
{
  MatrixT<DIMN, T> ret;
  ret.em() = x.em() * y.em();
  return ret;
}

template <int DIMN, class T>
MatrixT<DIMN, T> operator*(const T& x, const MatrixT<DIMN, T>& y)
{
  MatrixT<DIMN, T> ret;
  ret.em() = x * y.em();
  return ret;
}

template <int DIMN, class T>
MatrixT<DIMN, T> operator*(const MatrixT<DIMN, T>& x, const T& y)
{
  MatrixT<DIMN, T> ret;
  ret.em() = x.em() * y;
  return ret;
}

template <int DIMN, class T>
MatrixT<DIMN, T> operator/(const MatrixT<DIMN, T>& x, const T& y)
{
  MatrixT<DIMN, T> ret;
  ret.em() = x.em() / y;
  return ret;
}

template <int DIMN, class T>
MvectorT<DIMN, T> operator*(const MatrixT<DIMN, T>& x,
                            const MvectorT<DIMN, T>& y)
{
  MvectorT<DIMN, T> ret;
  ret.em() = x.em() * y.em();
  return ret;
}

template <int DIMN, class T>
void set_zero(MatrixT<DIMN, T>& m)
{
  memset(&m, 0, sizeof(MatrixT<DIMN, T>));
}

template <int DIMN, class T>
void set_unit(MatrixT<DIMN, T>& m, const Complex& coef = 1.0)
{
  set_zero(m);
  for (int i = 0; i < DIMN; ++i) {
    m(i, i) = coef;
  }
}

template <int DIMN, class T>
double qnorm(const MatrixT<DIMN, T>& m)
{
  return m.em().squaredNorm();
}

template <int DIMN, class T>
Complex matrix_trace(const MatrixT<DIMN, T>& x)
{
  return x.em().trace();
}

template <int DIMN, class T>
MatrixT<DIMN, T> matrix_adjoint(const MatrixT<DIMN, T>& x)
{
  MatrixT<DIMN, T> ret;
  ret.em() = x.em().adjoint();
  return ret;
}

template <int DIMN, class T>
MatrixT<DIMN, T> matrix_transpose(const MatrixT<DIMN, T>& x)
{
  MatrixT<DIMN, T> ret;
  ret.em() = x.em().transpose();
  return ret;
}

template <int DIMN, class T>
MatrixT<DIMN, T> matrix_conjugate(const MatrixT<DIMN, T>& x)
{
  MatrixT<DIMN, T> ret;
  ret.em() = x.em().conjugate();
  return ret;
}

template <class T = ComplexT>
struct ColorMatrixT : MatrixT<NUM_COLOR, T> {
  ColorMatrixT() {}
  ColorMatrixT(const MatrixT<NUM_COLOR, T>& m) { *this = m; }
  //
  const ColorMatrixT& operator=(const MatrixT<NUM_COLOR, T>& m)
  {
    *this = (const ColorMatrixT&)m;
    return *this;
  }
};

template <class T>
inline void unitarize(ColorMatrixT<T>& cm)
{
  cm.em().row(0).normalize();
  cm.em().row(2) =
      cm.em().row(1) - cm.em().row(0).dot(cm.em().row(1)) * cm.em().row(0);
  cm.em().row(1) = cm.em().row(2).normalized();
  cm.em().row(2) = cm.em().row(0).cross(cm.em().row(1));
}

inline ColorMatrixT<> make_anti_hermitian_matrix(const Array<double, 8> a)
{
  qassert(3 == NUM_COLOR);
  ColorMatrixT<> m;
  Array<double, 18> p(m.d());
  const double s3 =
      0.5773502691896258 * a[7];  // 1/sqrt(3) = 0.5773502691896258;
  p[0] = 0.0;
  p[8] = 0.0;
  p[16] = 0.0;
  p[1] = a[2] + s3;
  p[9] = -a[2] + s3;
  p[17] = -2.0 * s3;
  p[2] = a[1];
  p[3] = a[0];
  p[4] = a[4];
  p[5] = a[3];
  p[6] = -a[1];
  p[7] = a[0];
  p[10] = a[6];
  p[11] = a[5];
  p[12] = -a[4];
  p[13] = a[3];
  p[14] = -a[6];
  p[15] = a[5];
  return m;
}

inline ColorMatrixT<> make_g_rand_anti_hermitian_matrix(RngState& rs,
                                                        const double sigma)
{
  const double s = sigma / std::sqrt(2);
  std::array<double, 8> a;
  for (int i = 0; i < 8; ++i) {
    a[i] = g_rand_gen(rs, 0.0, s);
  }
  return make_anti_hermitian_matrix(a);
}

template <class T>
inline ColorMatrixT<T> make_color_matrix_exp(const ColorMatrixT<T>& a)
{
  ColorMatrixT<T> t2 = a;
  ColorMatrixT<T> t3 = a;
  ColorMatrixT<T> unit;
  set_unit(unit);
  for (int j = 9; j > 1; --j) {
    t3 = unit + (T)(1.0 / j) * t2;
    t2 = a * t3;
  }
  t3 = unit + t2;
  unitarize(t3);
  return t3;
}

template <class T = ComplexT>
struct WilsonMatrixT : MatrixT<4 * NUM_COLOR, T> {
  WilsonMatrixT() {}
  WilsonMatrixT(const MatrixT<4 * NUM_COLOR, T>& m) { *this = m; }
  //
  const WilsonMatrixT& operator=(const MatrixT<4 * NUM_COLOR, T>& m)
  {
    *this = (const WilsonMatrixT&)m;
    return *this;
  }
};

template <class T = ComplexT>
struct SpinMatrixT : MatrixT<4, T> {
  SpinMatrixT() {}
  SpinMatrixT(const MatrixT<4, T>& m) { *this = m; }
  //
  const SpinMatrixT& operator=(const MatrixT<4, T>& m)
  {
    *this = (const SpinMatrixT&)m;
    return *this;
  }
};

template <class T = ComplexT>
struct SpinMatrixConstantsT {
  SpinMatrixT<T> unit;
  std::array<SpinMatrixT<T>, 4>
      gammas;  // Not using CPS's convention, but a more standard one.
  std::array<SpinMatrixT<T>, 4> cps_gammas;  // CPS's convention gamma matrices
  SpinMatrixT<T> gamma5;                     // Same as CPS's gamma5
  std::array<SpinMatrixT<T>, 3> cap_sigmas;
  std::array<SpinMatrixT<T>, 3> cps_cap_sigmas;  // CPS's convention sigmas
  std::array<SpinMatrixT<T>, 16> gms;
  std::array<SpinMatrixT<T>, 16> cps_gms;
  //
  SpinMatrixConstantsT() { init(); }
  //
  void init()
  {
    TIMER_VERBOSE("SpinMatrixConstants::init()");
    unit.em() << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;
    // gamma_x
    gammas[0].em() << 0, 0, 0, 1, 0, 0, 1, 0, 0, -1, 0, 0, -1, 0, 0, 0;
    // gamma_y
    gammas[1].em() << 0, 0, 0, (T)-ii, 0, 0, (T)ii, 0, 0, (T)ii, 0, 0, (T)-ii, 0, 0, 0;
    // gamma_z
    gammas[2].em() << 0, 0, 1, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 1, 0, 0;
    gammas[0] *= (T)(-ii);
    gammas[1] *= (T)(-ii);
    gammas[2] *= (T)(-ii);
    // gamma_t
    gammas[3].em() << 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0;
    //
    cps_gammas[0] = -gammas[0];
    cps_gammas[1] = gammas[1];
    cps_gammas[2] = -gammas[2];
    cps_gammas[3] = gammas[3];
    // gamma_5 = gamma_x * gamma_y * gamma_z * gamma_t;
    gamma5.em() << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1;
    // Sigma_x
    cap_sigmas[0].em() << 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0;
    // Sigma_y
    cap_sigmas[1].em() << 0, (T)-ii, 0, 0, (T)ii, 0, 0, 0, 0, 0, 0, (T)-ii, 0, 0, (T)ii, 0;
    // Sigma_z
    cap_sigmas[2].em() << 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1;
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
  static const SpinMatrixConstantsT& get_instance()
  {
    static SpinMatrixConstantsT smcs;
    return smcs;
  }
  //
  static const SpinMatrixT<T>& get_unit() { return get_instance().unit; }
  static const SpinMatrixT<T>& get_gamma(int mu)
  {
    qassert(0 <= mu && mu < 4);
    return get_instance().gammas[mu];
  }
  static const std::array<SpinMatrixT<T>, 4>& get_gammas()
  {
    return get_instance().gammas;
  }
  static const std::array<SpinMatrixT<T>, 4>& get_cps_gammas()
  {
    return get_instance().cps_gammas;
  }
  static const std::array<SpinMatrixT<T>, 16>& get_gms()
  {
    return get_instance().gms;
  }
  static const std::array<SpinMatrixT<T>, 16>& get_cps_gms()
  {
    return get_instance().cps_gms;
  }
  static const SpinMatrixT<T>& get_gamma5() { return get_instance().gamma5; }
  static const SpinMatrixT<T>& get_cap_sigma(int i)
  {
    return get_instance().cap_sigmas[i];
  }
  static const std::array<SpinMatrixT<T>, 3>& get_cap_sigmas()
  {
    return get_instance().cap_sigmas;
  }
  static const std::array<SpinMatrixT<T>, 3>& get_cps_cap_sigmas()
  {
    return get_instance().cps_cap_sigmas;
  }
  static const SpinMatrixT<T> get_cps_sigmas(int i, int j)
  {
    return (get_instance().cps_gammas[i] * get_instance().cps_gammas[j] -
            get_instance().cps_gammas[j] * get_instance().cps_gammas[i]) /
           2.;
  }
};

template <class T>
WilsonMatrixT<T> operator*(const ColorMatrixT<T>& cm, const WilsonMatrixT<T>& m)
{
  WilsonMatrixT<T> ret;
  set_zero(ret);
  for (int s1 = 0; s1 < 4; ++s1) {
    for (int s2 = 0; s2 < 4; ++s2) {
      for (int c1 = 0; c1 < NUM_COLOR; ++c1) {
        for (int c2 = 0; c2 < NUM_COLOR; ++c2) {
          for (int c3 = 0; c3 < NUM_COLOR; ++c3) {
            ret(s1 * NUM_COLOR + c1, s2 * NUM_COLOR + c2) +=
                cm(c1, c3) * m(s1 * NUM_COLOR + c3, s2 * NUM_COLOR + c2);
          }
        }
      }
    }
  }
  return ret;
}

template <class T>
WilsonMatrixT<T> operator*(const WilsonMatrixT<T>& m, const ColorMatrixT<T>& cm)
{
  return matrix_adjoint((const ColorMatrixT<T>&)matrix_adjoint(cm) *
                        (const WilsonMatrixT<T>&)matrix_adjoint(m));
}

template <class T>
WilsonMatrixT<T> operator*(const SpinMatrixT<T>& sm, const WilsonMatrixT<T>& m)
{
  WilsonMatrixT<T> ret;
  set_zero(ret);
  for (int s1 = 0; s1 < 4; ++s1) {
    for (int s2 = 0; s2 < 4; ++s2) {
      for (int s3 = 0; s3 < 4; ++s3) {
        for (int c1 = 0; c1 < NUM_COLOR; ++c1) {
          for (int c2 = 0; c2 < NUM_COLOR; ++c2) {
            ret(s1 * NUM_COLOR + c1, s2 * NUM_COLOR + c2) +=
                sm(s1, s3) * m(s3 * NUM_COLOR + c1, s2 * NUM_COLOR + c2);
          }
        }
      }
    }
  }
  return ret;
}

template <class T>
WilsonMatrixT<T> operator*(const WilsonMatrixT<T>& m, const SpinMatrixT<T>& sm)
{
  return matrix_adjoint((const SpinMatrixT<T>&)matrix_adjoint(sm) *
                        (const WilsonMatrixT<T>&)matrix_adjoint(m));
}

template <class T>
WilsonVectorT<T> operator*(const ColorMatrixT<T>& cm, const WilsonVectorT<T>& m)
{
  WilsonVectorT<T> ret;
  set_zero(ret);
  for (int s1 = 0; s1 < 4; ++s1) {
    for (int c1 = 0; c1 < NUM_COLOR; ++c1) {
      for (int c2 = 0; c2 < NUM_COLOR; ++c2) {
        ret(s1 * NUM_COLOR + c1) += cm(c1, c2) * m(s1 * NUM_COLOR + c2);
      }
    }
  }
  return ret;
}

template <class T>
WilsonVectorT<T> operator*(const SpinMatrixT<T>& sm, const WilsonVectorT<T>& m)
{
  WilsonVectorT<T> ret;
  set_zero(ret);
  for (int s1 = 0; s1 < 4; ++s1) {
    for (int s2 = 0; s2 < 4; ++s2) {
      for (int c1 = 0; c1 < NUM_COLOR; ++c1) {
        ret(s1 * NUM_COLOR + c1) += sm(s1, s2) * m(s2 * NUM_COLOR + c1);
      }
    }
  }
  return ret;
}

#ifndef QLAT_NO_DEFAULT_TYPE

typedef ColorMatrixT<> ColorMatrix;

typedef WilsonMatrixT<> WilsonMatrix;

typedef SpinMatrixT<> SpinMatrix;

typedef SpinMatrixConstantsT<> SpinMatrixConstants;

#endif

QLAT_END_NAMESPACE

namespace qshow
{
template <int DIMN, class T>
std::string show(const qlat::MatrixT<DIMN, T>& m)
{
  std::ostringstream out;
  out << m.em();
  return out.str();
}

}  // namespace qshow

#ifndef USE_NAMESPACE
using namespace qshow;
#endif
