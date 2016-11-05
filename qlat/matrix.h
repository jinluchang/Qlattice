#pragma once

#include <qlat/field.h>
#include <qlat/field-fft.h>

#include <eigen3/Eigen/Eigen>

#include <cmath>

QLAT_START_NAMESPACE

typedef Eigen::Matrix<Complex,NUM_COLOR,NUM_COLOR,Eigen::RowMajor> ColorMatrix;

inline set_zero(ColorMatrix& m)
{
  m.setZero();
}

inline void set_unit(ColorMatrix& m, const Complex& coef = 1.0)
{
  set_zero(m);
  for (int i = 0; i < m.rows() && i < m.cols(); ++i) {
    m(i,i) = coef;
  }
}

inline double norm(const ColorMatrix& m)
{
  return m.squaredNorm();
}

inline std::string show(const ColorMatrix& m)
{
  std::ostringstream out;
  out << m;
  return out.str();
}

typedef Eigen::Matrix<Complex,4*NUM_COLOR,4*NUM_COLOR,Eigen::RowMajor> WilsonMatrix;

inline set_zero(WilsonMatrix& m)
{
  m.setZero();
}

inline void set_unit(WilsonMatrix& m, const Complex& coef = 1.0)
{
  set_zero(m);
  for (int i = 0; i < m.rows() && i < m.cols(); ++i) {
    m(i,i) = coef;
  }
}

inline double norm(const WilsonMatrix& m)
{
  return m.squaredNorm();
}

inline std::string show(const WilsonMatrix& m)
{
  std::ostringstream out;
  out << m;
  return out.str();
}

inline void unitarize(ColorMatrix& cm)
{
  cm.row(0).normalize();
  cm.row(2) = cm.row(1) - cm.row(0).dot(cm.row(1)) * cm.row(0);
  cm.row(1) = cm.row(2).normalized();
  cm.row(2) = cm.row(0).cross(cm.row(1));
}

typedef Eigen::Matrix<Complex,4,4,Eigen::RowMajor> SpinMatrix;

inline set_zero(SpinMatrix& m)
{
  m.setZero();
}

inline void set_unit(SpinMatrix& m, const Complex& coef = 1.0)
{
  set_zero(m);
  for (int i = 0; i < m.rows() && i < m.cols(); ++i) {
    m(i,i) = coef;
  }
}

inline double norm(const SpinMatrix& m)
{
  return m.squaredNorm();
}

inline std::string show(const SpinMatrix& m)
{
  std::ostringstream out;
  out << m;
  return out.str();
}

struct SpinMatrixConstants
{
  SpinMatrix unit;
  std::array<SpinMatrix,4> gammas;
  // Not using CPS's convention, but a more standard one.
  SpinMatrix gamma5;
  // Same as CPS's gamma5
  std::array<SpinMatrix,3> cap_sigmas;
  //
  SpinMatrixConstants()
  {
    init();
  }
  //
  void init()
  {
    unit <<
        1,   0,   0,   0,
        0,   1,   0,   0,
        0,   0,   1,   0,
        0,   0,   0,   1;
    // gamma_x
    gammas[0] <<
        0,   0,   0,   1,
        0,   0,   1,   0,
        0,  -1,   0,   0,
       -1,   0,   0,   0;
    // gamma_y
    gammas[1] <<
        0,   0,   0, -ii,
        0,   0,  ii,   0,
        0,  ii,   0,   0,
      -ii,   0,   0,   0;
    // gamma_z
    gammas[2] <<
        0,   0,   1,   0,
        0,   0,   0,  -1,
       -1,   0,   0,   0,
        0,   1,   0,   0;
    gammas[0] *= -ii;
    gammas[1] *= -ii;
    gammas[2] *= -ii;
    // gamma_t
    gammas[3] <<
        0,   0,   1,   0,
        0,   0,   0,   1,
        1,   0,   0,   0,
        0,   1,   0,   0;
    // gamma_5
    gamma5 <<
        1,   0,   0,   0,
        0,   1,   0,   0,
        0,   0,  -1,   0,
        0,   0,   0,  -1;
    // Sigma_x
    cap_sigmas[0] <<
        0,   1,   0,   0,
        1,   0,   0,   0,
        0,   0,   0,   1,
        0,   0,   1,   0;
    // Sigma_y
    cap_sigmas[1] <<
        0, -ii,   0,   0,
       ii,   0,   0,   0,
        0,   0,   0, -ii,
        0,   0,  ii,   0;
    // Sigma_z
    cap_sigmas[2] <<
        1,   0,   0,   0,
        0,  -1,   0,   0,
        0,   0,   1,   0,
        0,   0,   0,  -1;
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
  static const SpinMatrix& get_gamma5()
  {
    return get_instance().gamma5;
  }
  static const SpinMatrix& get_cap_sigma(int i)
  {
    return get_instance().cap_sigmas[i];
  }
};

QLAT_END_NAMESPACE
