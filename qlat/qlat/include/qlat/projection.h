#pragma once

#include "muon-line-config.h"

#include <qlat/qlat.h>

#include <iostream>
#include <cmath>
#include <cassert>

namespace qlat
{  //

template <class T, std::size_t N>
inline std::string show(const array<T, N>& mm)
{
  if (0 == N) {
    return "()";
  } else {
    std::ostringstream out;
    out << "(";
    out << show(mm[0]);
    for (std::size_t i = 1; i < N; ++i) {
      out << "," << show(mm[i]);
    }
    out << ")";
    return out.str();
  }
}

template <class T>
inline std::string show(const std::vector<T>& mm)
{
  if (0 == mm.size()) {
    return "[]";
  } else {
    std::ostringstream out;
    out << "[";
    out << show(mm[0]);
    for (std::size_t i = 1; i < mm.size(); ++i) {
      out << "," << show(mm[i]);
    }
    out << "]";
    return out.str();
  }
}

typedef array<RealD, 3> MagneticMoment;

typedef array<MagneticMoment, 4 * 4 * 4> ManyMagneticMoments;

qacc MagneticMoment operator*(const RealD x, const MagneticMoment& m)
{
  MagneticMoment ret;
  ret[0] = x * m[0];
  ret[1] = x * m[1];
  ret[2] = x * m[2];
  return ret;
}

qacc MagneticMoment operator*(const MagneticMoment& m, const RealD x)
{
  return x * m;
}

qacc MagneticMoment operator+(const MagneticMoment& m1,
                              const MagneticMoment& m2)
{
  MagneticMoment ret;
  ret[0] = m1[0] + m2[0];
  ret[1] = m1[1] + m2[1];
  ret[2] = m1[2] + m2[2];
  return ret;
}

qacc MagneticMoment operator-(const MagneticMoment& m1,
                              const MagneticMoment& m2)
{
  MagneticMoment ret;
  ret[0] = m1[0] - m2[0];
  ret[1] = m1[1] - m2[1];
  ret[2] = m1[2] - m2[2];
  return ret;
}

qacc ManyMagneticMoments operator*(const RealD a, const ManyMagneticMoments& m)
{
  ManyMagneticMoments ret;
  for (std::size_t i = 0; i < ret.size(); ++i) {
    ret[i] = a * m[i];
  }
  return ret;
}

qacc ManyMagneticMoments operator*(const ManyMagneticMoments& m, const RealD a)
{
  return a * m;
}

qacc ManyMagneticMoments& operator*=(ManyMagneticMoments& m, const RealD a)
{
  m = m * a;
  return m;
}

qacc ManyMagneticMoments operator+(const ManyMagneticMoments& m1,
                                   const ManyMagneticMoments& m2)
{
  ManyMagneticMoments ret;
  for (std::size_t i = 0; i < ret.size(); ++i) {
    ret[i] = m1[i] + m2[i];
  }
  return ret;
}

qacc ManyMagneticMoments operator-(const ManyMagneticMoments& m1,
                                   const ManyMagneticMoments& m2)
{
  ManyMagneticMoments ret;
  for (std::size_t i = 0; i < ret.size(); ++i) {
    ret[i] = m1[i] - m2[i];
  }
  return ret;
}

qacc ManyMagneticMoments& operator+=(ManyMagneticMoments& m1,
                                     const ManyMagneticMoments& m2)
{
  m1 = m1 + m2;
  return m1;
}

inline MagneticMoment computeProjection(const Int nu, const Int rho,
                                        const Int mu,
                                        const std::vector<RealD>& integral)
{
  assert(integral.size() == 5 * 5);
  const SpinMatrixConstants& smc = SpinMatrixConstants::get_instance();
  const SpinMatrix& unit = smc.unit;
  const array<SpinMatrix, 3>& cap_sigmas = smc.cap_sigmas;
  array<SpinMatrix, 5> gammas;
  for (Int i = 0; i < 4; ++i) {
    gammas[i] = smc.gammas[i];
  }
  gammas[4] = gammas[3] + unit;
  SpinMatrix proj =
      (Complex)0.5 * (unit + gammas[3]);  // proj = (1 + gamma_0) / 2
  SpinMatrix magMat;
  set_zero(magMat);
  for (Int i = 0; i < 5; ++i) {
    for (Int j = 0; j < 5; ++j) {
      magMat += -ii / (4.0 * sqr(PI)) * integral[i * 5 + j] * proj *
                gammas[nu] * gammas[i] * gammas[rho] * gammas[j] * gammas[mu] *
                proj;
    }
  }
  MagneticMoment ans;
  for (Int i = 0; i < 3; ++i) {
    ans[i] = 0.5 * matrix_trace(magMat * cap_sigmas[i])
                       .real();  // indeed should be 0.5, note the proj matrix
  }
  return ans;
}

inline ManyMagneticMoments computeProjections(
    const std::vector<RealD>& integral)
{
  ManyMagneticMoments ans;
  for (Int nu = 0; nu < 4; ++nu) {
    for (Int rho = 0; rho < 4; ++rho) {
      for (Int mu = 0; mu < 4; ++mu) {
        // ans[y,z,x]
        ans[nu * 16 + rho * 4 + mu] = computeProjection(nu, rho, mu, integral);
      }
    }
  }
  return ans;
}

inline ManyMagneticMoments permuteNuRhoMu(const ManyMagneticMoments& mmm,
                                          const Int i, const Int j, const Int k)
// i, j, k are different integer, range from 0, 1, 2
// represent different permutation
{
  assert(0 <= i && i < 3);
  assert(0 <= j && j < 3);
  assert(0 <= k && k < 3);
  assert(i != j && i != k && j != k);
  ManyMagneticMoments ans;
  array<int, 3> nms, oms;
  for (Int nu = 0; nu < 4; ++nu) {
    for (Int rho = 0; rho < 4; ++rho) {
      for (Int mu = 0; mu < 4; ++mu) {
        oms[0] = nu;
        oms[1] = rho;
        oms[2] = mu;
        nms[0] = oms[i];
        nms[1] = oms[j];
        nms[2] = oms[k];
        ans[nms[0] * 16 + nms[1] * 4 + nms[2]] =
            mmm[oms[0] * 16 + oms[1] * 4 + oms[2]];
      }
    }
  }
  return ans;
}

qacc const RealD& get_m_comp(const ManyMagneticMoments& mmm, const Int i,
                              const Int rho, const Int sigma, const Int nu)
{
  return mmm[16 * sigma + 4 * nu + rho][i];
}

qacc RealD& get_m_comp(ManyMagneticMoments& mmm, const Int i, const Int rho,
                        const Int sigma, const Int nu)
{
  return mmm[16 * sigma + 4 * nu + rho][i];
}

qacc ManyMagneticMoments permute_rho_sigma_nu(const ManyMagneticMoments& mmm,
                                              const Int i, const Int j,
                                              const Int k)
// i, j, k are different integer, range from 0, 1, 2
// represent different permutation
{
  assert(0 <= i && i < 3);
  assert(0 <= j && j < 3);
  assert(0 <= k && k < 3);
  assert(i != j && i != k && j != k);
  ManyMagneticMoments ans;
  array<int, 3> oms;
  for (Int nu = 0; nu < 4; ++nu) {
    for (Int rho = 0; rho < 4; ++rho) {
      for (Int mu = 0; mu < 4; ++mu) {
        oms[0] = nu;
        oms[1] = rho;
        oms[2] = mu;
        for (Int a = 0; a < 3; ++a) {
          get_m_comp(ans, a, oms[i], oms[j], oms[k]) =
              get_m_comp(mmm, a, oms[0], oms[1], oms[2]);
        }
      }
    }
  }
  return ans;
}

inline std::string showMagneticMoment(const MagneticMoment& mm)
{
  return ssprintf("(%9.6f,%9.6f,%9.6f)", mm[0], mm[1], mm[2]);
}

inline std::string showManyMagneticMoments(const ManyMagneticMoments& mmm)
{
  const RealD sqrt_qnorm = std::sqrt(qnorm(mmm));
  std::ostringstream out;
  out << "qnorm = " << show(sqrt_qnorm) << std::endl;
  for (Int nu = 0; nu < 4; ++nu) {
    for (Int rho = 0; rho < 4; ++rho) {
      for (Int mu = 0; mu < 4; ++mu) {
        const MagneticMoment& mm = mmm[nu * 16 + rho * 4 + mu];
        out << show(nu) + show(rho) + show(mu) << " "
            << ssprintf("%9.6f", std::sqrt(qnorm(mm)) / sqrt_qnorm) << " "
            << showMagneticMoment((1.0 / sqrt_qnorm) * mm) << std::endl;
      }
    }
  }
  return out.str();
}

inline ManyMagneticMoments averageManyMagneticMoments(
    const std::vector<ManyMagneticMoments>& mmms)
{
  const std::size_t size = mmms.size();
  ManyMagneticMoments mmm;
  set_zero(mmm);
  if (size > 0) {
    for (std::size_t i = 0; i < size; ++i) {
      for (size_t m = 0; m < 4 * 4 * 4; ++m) {
        for (Int k = 0; k < 3; ++k) {
          mmm[m][k] += mmms[i][m][k];
        }
      }
    }
    for (std::size_t m = 0; m < 4 * 4 * 4; ++m) {
      for (Int k = 0; k < 3; ++k) {
        mmm[m][k] /= (RealD)size;
      }
    }
  }
  return mmm;
}

inline void test_projection()
{
  TIMER_VERBOSE("test_projection");
  const SpinMatrixConstants& smc = SpinMatrixConstants::get_instance();
  const SpinMatrix& unit = smc.unit;
  const array<SpinMatrix, 4>& gammas = smc.gammas;
  const array<SpinMatrix, 3>& cap_sigmas = smc.cap_sigmas;
  SpinMatrix proj =
      (Complex)0.5 * (unit + gammas[3]);  // proj = (1 + gamma_0) / 2
  displayln_info(
      shows(0.5 * matrix_trace(proj * cap_sigmas[1] * proj * cap_sigmas[1])));
  std::vector<RealD> integral(25);
  set_zero(integral);
  RngState rs("test_projection");
  for (Int i = 0; i < 5; ++i) {
    for (Int j = 0; j < 5; ++j) {
      if (3 > i && 3 > j) {
        integral[i * 5 + j] = u_rand_gen(rs);
      }
    }
  }
  displayln_info(
      shows(0.5 * matrix_trace(proj * cap_sigmas[1] * proj * cap_sigmas[1])));
  MagneticMoment mm = computeProjection(3, 3, 3, integral);
  displayln_info(show(mm));
  ManyMagneticMoments mmm = computeProjections(integral);
  displayln_info(showManyMagneticMoments(mmm));
}

}  // namespace qlat
