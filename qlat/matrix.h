#pragma once

#include <qlat/config.h>
#include <eigen3/Eigen/Eigen>

#include <cmath>

QLAT_START_NAMESPACE

template <int DIMN>
struct Matrix
{
  Complex p[DIMN * DIMN];
  //
  // convert to double array
  double* d()
  {
    return (double*)p;
  }
  const double* d() const
  {
    return (const double*)p;
  }
  //
  // convert to Eigen Matrix
  Eigen::Matrix<Complex,DIMN,DIMN,Eigen::RowMajor>& em()
  {
    return *((Eigen::Matrix<Complex,DIMN,DIMN,Eigen::RowMajor>*)this);
  }
  const Eigen::Matrix<Complex,DIMN,DIMN,Eigen::RowMajor>& em() const
  {
    return *((Eigen::Matrix<Complex,DIMN,DIMN,Eigen::RowMajor>*)this);
  }
  //
  Complex& operator()(int i, int j)
  {
    qassert(0 <= i && i < DIMN);
    qassert(0 <= j && j < DIMN);
    return p[i * DIMN + j];
  }
  const Complex& operator()(int i, int j) const
  {
    qassert(0 <= i && i < DIMN);
    qassert(0 <= j && j < DIMN);
    return p[i * DIMN + j];
  }
  //
  const Matrix& operator+=(const Matrix& x)
  {
    *this = *this + x;
    return *this;
  }
  //
  const Matrix& operator-=(const Matrix& x)
  {
    *this = *this - x;
    return *this;
  }
  //
  const Matrix& operator*=(const Matrix& x)
  {
    *this = *this * x;
    return *this;
  }
  //
  const Matrix& operator*=(const Complex& x)
  {
    *this = *this * x;
    return *this;
  }
  //
  const Matrix& operator/=(const Complex& x)
  {
    *this = *this / x;
    return *this;
  }
};

template <int DIMN>
Matrix<DIMN> operator+(const Matrix<DIMN>& x, const Matrix<DIMN>& y)
{
  Matrix<DIMN> ret;
  ret.em() = x.em() + y.em();
  return ret;
}

template <int DIMN>
Matrix<DIMN> operator-(const Matrix<DIMN>& x, const Matrix<DIMN>& y)
{
  Matrix<DIMN> ret;
  ret.em() = x.em() - y.em();
  return ret;
}

template <int DIMN>
Matrix<DIMN> operator-(const Matrix<DIMN>& x)
{
  Matrix<DIMN> ret;
  ret.em() = -x.em();
  return ret;
}

template <int DIMN>
Matrix<DIMN> operator*(const Matrix<DIMN>& x, const Matrix<DIMN>& y)
{
  Matrix<DIMN> ret;
  ret.em() = x.em() * y.em();
  return ret;
}

template <int DIMN>
Matrix<DIMN> operator*(const Complex& x, const Matrix<DIMN>& y)
{
  Matrix<DIMN> ret;
  ret.em() = x * y.em();
  return ret;
}

template <int DIMN>
Matrix<DIMN> operator*(const Matrix<DIMN>& x, const Complex& y)
{
  Matrix<DIMN> ret;
  ret.em() = x.em() * y;
  return ret;
}

template <int DIMN>
Matrix<DIMN> operator/(const Matrix<DIMN>& x, const Complex& y)
{
  Matrix<DIMN> ret;
  ret.em() = x.em() / y;
  return ret;
}

template <int DIMN>
Mvector<DIMN> operator*(const Matrix<DIMN>& x, const Mvector<DIMN>& y)
{
  Mvector<DIMN> ret;
  ret.em() = x.em() * y.em();
  return ret;
}

template <int DIMN>
void set_zero(Matrix<DIMN>& m)
{
  memset(&m, 0, sizeof(Matrix<DIMN>));
}

template <int DIMN>
void set_unit(Matrix<DIMN>& m, const Complex& coef = 1.0)
{
  set_zero(m);
  for (int i = 0; i < DIMN; ++i) {
    m(i,i) = coef;
  }
}

template <int DIMN>
double norm(const Matrix<DIMN>& m)
{
  return m.em().squaredNorm();
}

template <int DIMN>
Complex matrix_trace(const Matrix<DIMN>& x)
{
  return x.em().trace();
}

template <int DIMN>
Matrix<DIMN> matrix_adjoint(const Matrix<DIMN>& x)
{
  Matrix<DIMN> ret;
  ret.em() = x.em().adjoint();
  return ret;
}

struct ColorMatrix : Matrix<NUM_COLOR>
{
  ColorMatrix()
  {
  }
  ColorMatrix(const Matrix<NUM_COLOR>& m)
  {
    *this = m;
  }
  //
  const ColorMatrix& operator=(const Matrix<NUM_COLOR>& m)
  {
    *this = (const ColorMatrix&)m;
    return *this;
  }
};

inline void unitarize(ColorMatrix& cm)
{
  cm.em().row(0).normalize();
  cm.em().row(2) = cm.em().row(1) - cm.em().row(0).dot(cm.em().row(1)) * cm.em().row(0);
  cm.em().row(1) = cm.em().row(2).normalized();
  cm.em().row(2) = cm.em().row(0).cross(cm.em().row(1));
}

inline ColorMatrix make_anti_hermitian_matrix(const Array<double, 8> a)
{
  qassert(3 == NUM_COLOR);
  ColorMatrix m;
  Array<double,18> p(m.d());
  const double s3 = 0.5773502691896258 * a[7];       // 1/sqrt(3) = 0.5773502691896258;
  p[0 ] =  0.0;
  p[8 ] =  0.0;
  p[16] =  0.0;
  p[1 ] =  a[2 ] + s3;
  p[9 ] = -a[2 ] + s3;
  p[17] = -2.0   * s3;
  p[2 ] =  a[1 ];
  p[3 ] =  a[0 ];
  p[4 ] =  a[4 ];
  p[5 ] =  a[3 ];
  p[6 ] = -a[1 ];
  p[7 ] =  a[0 ];
  p[10] =  a[6 ];
  p[11] =  a[5 ];
  p[12] = -a[4 ];
  p[13] =  a[3 ];
  p[14] = -a[6 ];
  p[15] =  a[5 ];
  return m;
}

inline ColorMatrix make_g_rand_anti_hermitian_matrix(RngState& rs, const double sigma)
{
  const double s = sigma / std::sqrt(2);
  std::array<double, 8> a;
  for (int i = 0; i < 8; ++i) {
    a[i] = g_rand_gen(rs, 0.0, s);
  }
  return make_anti_hermitian_matrix(a);
}

inline ColorMatrix make_color_matrix_exp(const ColorMatrix& a)
{
  ColorMatrix t2 = a;
  ColorMatrix t3 = a;
  ColorMatrix unit;
  set_unit(unit);
  for(int j = 9; j > 1; --j) {
    t3 = unit + (1.0/j) * t2;
    t2 = a * t3;
  }
  t3 = unit + t2;
  unitarize(t3);
  return t3;
}

struct WilsonMatrix : Matrix<4*NUM_COLOR>
{
  WilsonMatrix()
  {
  }
  WilsonMatrix(const Matrix<4*NUM_COLOR>& m)
  {
    *this = m;
  }
  //
  const WilsonMatrix& operator=(const Matrix<4*NUM_COLOR>& m)
  {
    *this = (const WilsonMatrix&)m;
    return *this;
  }
};

struct SpinMatrix : Matrix<4>
{
  SpinMatrix()
  {
  }
  SpinMatrix(const Matrix<4>& m)
  {
    *this = m;
  }
  //
  const SpinMatrix& operator=(const Matrix<4>& m)
  {
    *this = (const SpinMatrix&)m;
    return *this;
  }
};

struct SpinMatrixConstants
{
  SpinMatrix unit;
  std::array<SpinMatrix,4> gammas; // Not using CPS's convention, but a more standard one.
  std::array<SpinMatrix,4> cps_gammas; // CPS's convention gamma matrices
  SpinMatrix gamma5; // Same as CPS's gamma5
  std::array<SpinMatrix,3> cap_sigmas;
  std::array<SpinMatrix,3> cps_cap_sigmas; // CPS's convention sigmas
  std::array<SpinMatrix,16> gms;
  std::array<SpinMatrix,16> cps_gms;
  //
  SpinMatrixConstants()
  {
    init();
  }
  //
  void init()
  {
    unit.em() <<
        1,   0,   0,   0,
        0,   1,   0,   0,
        0,   0,   1,   0,
        0,   0,   0,   1;
    // gamma_x
    gammas[0].em() <<
        0,   0,   0,   1,
        0,   0,   1,   0,
        0,  -1,   0,   0,
       -1,   0,   0,   0;
    // gamma_y
    gammas[1].em() <<
        0,   0,   0, -ii,
        0,   0,  ii,   0,
        0,  ii,   0,   0,
      -ii,   0,   0,   0;
    // gamma_z
    gammas[2].em() <<
        0,   0,   1,   0,
        0,   0,   0,  -1,
       -1,   0,   0,   0,
        0,   1,   0,   0;
    gammas[0] *= -ii;
    gammas[1] *= -ii;
    gammas[2] *= -ii;
    // gamma_t
    gammas[3].em() <<
        0,   0,   1,   0,
        0,   0,   0,   1,
        1,   0,   0,   0,
        0,   1,   0,   0;
    //
    cps_gammas[0] = -gammas[0];
    cps_gammas[1] = gammas[1];
    cps_gammas[2] = -gammas[2];
    cps_gammas[3] = gammas[3];
    // gamma_5
    gamma5.em() <<
        1,   0,   0,   0,
        0,   1,   0,   0,
        0,   0,  -1,   0,
        0,   0,   0,  -1;
    // Sigma_x
    cap_sigmas[0].em() <<
        0,   1,   0,   0,
        1,   0,   0,   0,
        0,   0,   0,   1,
        0,   0,   1,   0;
    // Sigma_y
    cap_sigmas[1].em() <<
        0, -ii,   0,   0,
       ii,   0,   0,   0,
        0,   0,   0, -ii,
        0,   0,  ii,   0;
    // Sigma_z
    cap_sigmas[2].em() <<
        1,   0,   0,   0,
        0,  -1,   0,   0,
        0,   0,   1,   0,
        0,   0,   0,  -1;
    //
    cps_cap_sigmas[0] = -cap_sigmas[0];
    cps_cap_sigmas[1] = cap_sigmas[1];
    cps_cap_sigmas[2] = -cap_sigmas[2];
    //
    for (int a = 0; a < 2; ++a) {
      const SpinMatrix ma = a == 0 ? unit : gammas[0];
      for (int b = 0; b < 2; ++b) {
        const SpinMatrix mb = b == 0 ? ma : (SpinMatrix)(ma * gammas[1]);
        for (int c = 0; c < 2; ++c) {
          const SpinMatrix mc = c == 0 ? mb : (SpinMatrix)(mb * gammas[2]);
          for (int d = 0; d < 2; ++d) {
            const SpinMatrix md = d == 0 ? mc : (SpinMatrix)(mc * gammas[3]);
            const int idx = a + 2 * b + 4 * c + 8 * d;
            gms[idx] = md;
          }
        }
      }
    }
    for (int a = 0; a < 2; ++a) {
      const SpinMatrix ma = a == 0 ? unit : cps_gammas[0];
      for (int b = 0; b < 2; ++b) {
        const SpinMatrix mb = b == 0 ? ma : (SpinMatrix)(ma * cps_gammas[1]);
        for (int c = 0; c < 2; ++c) {
          const SpinMatrix mc = c == 0 ? mb : (SpinMatrix)(mb * cps_gammas[2]);
          for (int d = 0; d < 2; ++d) {
            const SpinMatrix md = d == 0 ? mc : (SpinMatrix)(mc * cps_gammas[3]);
            const int idx = a + 2 * b + 4 * c + 8 * d;
            cps_gms[idx] = md;
          }
        }
      }
    }
  }
  //
  static const SpinMatrixConstants& get_instance()
  {
    static SpinMatrixConstants smcs;
    return smcs;
  }
  //
  static const SpinMatrix& get_unit()
  {
    return get_instance().unit;
  }
  static const SpinMatrix& get_gamma(int mu)
  {
    qassert(0 <= mu && mu < 4);
    return get_instance().gammas[mu];
  }
  static const std::array<SpinMatrix,4>& get_gammas()
  {
    return get_instance().gammas;
  }
  static const std::array<SpinMatrix,4>& get_cps_gammas()
  {
    return get_instance().cps_gammas;
  }
  static const std::array<SpinMatrix,16>& get_gms()
  {
    return get_instance().gms;
  }
  static const std::array<SpinMatrix,16>& get_cps_gms()
  {
    return get_instance().cps_gms;
  }
  static const SpinMatrix& get_gamma5()
  {
    return get_instance().gamma5;
  }
  static const SpinMatrix& get_cap_sigma(int i)
  {
    return get_instance().cap_sigmas[i];
  }
  static const std::array<SpinMatrix,3>& get_cap_sigmas()
  {
    return get_instance().cap_sigmas;
  }
  static const std::array<SpinMatrix,3>& get_cps_cap_sigmas()
  {
    return get_instance().cps_cap_sigmas;
  }
  static const SpinMatrix get_cps_sigmas(int i, int j)
  {
    return (get_instance().cps_gammas[i]*get_instance().cps_gammas[j]-get_instance().cps_gammas[j]*get_instance().cps_gammas[i])/2.;
  }
};

inline WilsonMatrix operator*(const ColorMatrix& cm, const WilsonMatrix& m)
{
  WilsonMatrix ret;
  set_zero(ret);
  for (int s1 = 0; s1 < 4; ++s1) {
    for (int s2 = 0; s2 < 4; ++s2) {
      for (int c1 = 0; c1 < NUM_COLOR; ++c1) {
        for (int c2 = 0; c2 < NUM_COLOR; ++c2) {
          for (int c3 = 0; c3 < NUM_COLOR; ++c3) {
            ret(s1*NUM_COLOR+c1, s2*NUM_COLOR+c2) += cm(c1,c3) * m(s1*NUM_COLOR+c3, s2*NUM_COLOR+c2);
          }
        }
      }
    }
  }
  return ret;
}

inline WilsonMatrix operator*(const WilsonMatrix& m, const ColorMatrix& cm)
{
  return matrix_adjoint(matrix_adjoint(cm) * matrix_adjoint(m));
}

inline WilsonMatrix operator*(const SpinMatrix& sm, const WilsonMatrix& m)
{
  WilsonMatrix ret;
  set_zero(ret);
  for (int s1 = 0; s1 < 4; ++s1) {
    for (int s2 = 0; s2 < 4; ++s2) {
      for (int s3 = 0; s3 < 4; ++s3) {
        for (int c1 = 0; c1 < NUM_COLOR; ++c1) {
          for (int c2 = 0; c2 < NUM_COLOR; ++c2) {
            ret(s1*NUM_COLOR+c1, s2*NUM_COLOR+c2) += sm(s1,s3) * m(s3*NUM_COLOR+c1, s2*NUM_COLOR+c2);
          }
        }
      }
    }
  }
  return ret;
}

inline WilsonMatrix operator*(const WilsonMatrix& m, const SpinMatrix& sm)
{
  return matrix_adjoint(matrix_adjoint(sm) * matrix_adjoint(m));
}

inline WilsonVector operator*(const ColorMatrix& cm, const WilsonVector& m)
{
  WilsonVector ret;
  set_zero(ret);
  for (int s1 = 0; s1 < 4; ++s1) {
    for (int c1 = 0; c1 < NUM_COLOR; ++c1) {
      for (int c2 = 0; c2 < NUM_COLOR; ++c2) {
        ret(s1*NUM_COLOR+c1) += cm(c1,c2) * m(s1*NUM_COLOR+c2);
      }
    }
  }
  return ret;
}

inline WilsonVector operator*(const SpinMatrix& sm, const WilsonVector& m)
{
  WilsonVector ret;
  set_zero(ret);
  for (int s1 = 0; s1 < 4; ++s1) {
    for (int s2 = 0; s2 < 4; ++s2) {
      for (int c1 = 0; c1 < NUM_COLOR; ++c1) {
        ret(s1*NUM_COLOR+c1) += sm(s1,s2) * m(s2*NUM_COLOR+c1);
      }
    }
  }
  return ret;
}

QLAT_END_NAMESPACE

namespace qshow {

template <int DIMN>
std::string show(const qlat::Matrix<DIMN>& m)
{
  std::ostringstream out;
  out << m.em();
  return out.str();
}

}

#ifndef USE_NAMESPACE
using namespace qshow;
#endif
