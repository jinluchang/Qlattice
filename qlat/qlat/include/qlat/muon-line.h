#pragma once

#include "integration-multidimensional.h"
#include "compute-int-mult.h"
#include "projection.h"
#include "interpolation-bilinear.h"

#include <qlat/qlat.h>

#include <cmath>
#include <cstdlib>
#include <cassert>
#include <string>

// Eq. (9) of https://arxiv.org/pdf/2304.04423.pdf
//
// % \mathcal M_{i,\rho,\sigma,\lambda}(x,y,z)
// % = muonLineSym(m * (x - z), m * (y - z))[16 * sigma + 4 * lambda + rho][i]
// % = muon_line_sym_py(m * (x - z), m * (y - z))[i * 64 + rho * 16 + sigma * 4 + lambda]
//
// % interpolate saved data and extrapolate different interpolations
//
// % = get_muon_line_m_extra(m * x, m * y, m * z , tag)[16 * sigma  + 4 * lambda + rho][i]
//
// % = get_m_comp(get_muon_line_m_extra_lat(x, y, z, total_site, m, tag), i, rho, sigma, lambda)
//
// % tag = 0 : subtraction
// % tag = 1 : no subtraction
//
// % m is the muon mass in the lattice unit (or lattice spacing in muon mass unit)
// % total_site is the lattice size in lattice unit

namespace qlat
{  //

inline Int qrename_partial(const std::string& path)
{
  TIMER("qrename_partial");
  return qrename(path + ".partial", path);
}

inline Int qrename_partial_info(const std::string& path)
{
  TIMER("qrename_partial_info");
  return qrename_info(path + ".partial", path);
}

inline std::vector<std::string> get_lines(std::istream& is)
{
  std::vector<std::string> ret;
  while (is.good()) {
    std::string str;
    std::getline(is, str);
    ret.push_back(str);
  }
  return ret;
}

inline std::vector<std::string> lines(const std::string& s)
{
  std::istringstream is(s);
  return get_lines(is);
}

inline std::string pad_cut_string(const std::string& s, const Int width,
                                  const char c = ' ')
{
  std::string ret(s);
  ret.resize(width, c);
  return ret;
}

inline std::string compare_multiline_string(const std::string& s1,
                                            const std::string& s2,
                                            const Int width = 80)
{
  std::vector<std::string> ls1(lines(s1)), ls2(lines(s2));
  std::ostringstream out;
  const size_t size = std::max(ls1.size(), ls2.size());
  ls1.resize(size, "");
  ls2.resize(size, "");
  for (size_t i = 0; i < size; ++i) {
    out << pad_cut_string(ls1[i], width) << " | "
        << pad_cut_string(ls2[i], width) << std::endl;
  }
  return out.str();
}

typedef Eigen::Matrix<RealD, 3, 3, Eigen::RowMajor> SpatialO3Matrix;

inline SpatialO3Matrix makeRotationAroundX(const RealD phi)
{
  using namespace std;
  SpatialO3Matrix ret;
  ret << 1, 0, 0, 0, cos(phi), -sin(phi), 0, sin(phi), cos(phi);
  return ret;
}

inline SpatialO3Matrix makeRotationAroundZ(const RealD phi)
{
  using namespace std;
  SpatialO3Matrix ret;
  ret << cos(phi), -sin(phi), 0, sin(phi), cos(phi), 0, 0, 0, 1;
  return ret;
}

inline SpatialO3Matrix makeRandomSpatialO3Matrix(RngState& rs)
{
  const RealD reflex_x = u_rand_gen(rs) >= 0.5 ? 1.0 : -1.0;
  const RealD reflex_y = u_rand_gen(rs) >= 0.5 ? 1.0 : -1.0;
  const RealD reflex_z = u_rand_gen(rs) >= 0.5 ? 1.0 : -1.0;
  SpatialO3Matrix reflex;
  reflex << reflex_x, 0, 0, 0, reflex_y, 0, 0, 0, reflex_z;
  const RealD phi1 = u_rand_gen(rs, PI, -PI);
  const RealD phi2 = u_rand_gen(rs, PI, -PI);
  const RealD phi3 = u_rand_gen(rs, PI, -PI);
  return reflex * makeRotationAroundX(phi3) * makeRotationAroundZ(phi2) *
         makeRotationAroundX(phi1);
}

inline CoordinateD operator*(const SpatialO3Matrix& m, const CoordinateD& x)
{
  Eigen::Matrix<RealD, 3, 1> vec;
  vec << x[0], x[1], x[2];
  vec = m * vec;
  CoordinateD ret;
  ret[0] = vec[0];
  ret[1] = vec[1];
  ret[2] = vec[2];
  ret[3] = x[3];
  return ret;
}

typedef array<RealD, 92> ManyMagneticMomentsCompressed;

inline ManyMagneticMomentsCompressed operator*(
    const RealD a, const ManyMagneticMomentsCompressed& m)
{
  ManyMagneticMomentsCompressed ret;
  for (size_t i = 0; i < ret.size(); ++i) {
    ret[i] = a * m[i];
  }
  return ret;
}

inline ManyMagneticMomentsCompressed operator*(
    const ManyMagneticMomentsCompressed& m, const RealD a)
{
  return a * m;
}

inline ManyMagneticMomentsCompressed& operator*=(
    ManyMagneticMomentsCompressed& m, const RealD a)
{
  m = a * m;
  return m;
}

inline ManyMagneticMomentsCompressed operator+(
    const ManyMagneticMomentsCompressed& m1,
    const ManyMagneticMomentsCompressed& m2)
{
  ManyMagneticMomentsCompressed ret;
  for (size_t i = 0; i < ret.size(); ++i) {
    ret[i] = m1[i] + m2[i];
  }
  return ret;
}

inline ManyMagneticMomentsCompressed operator-(
    const ManyMagneticMomentsCompressed& m1,
    const ManyMagneticMomentsCompressed& m2)
{
  ManyMagneticMomentsCompressed ret;
  for (size_t i = 0; i < ret.size(); ++i) {
    ret[i] = m1[i] - m2[i];
  }
  return ret;
}

inline SpatialO3Matrix makeProperRotation(const CoordinateD& x,
                                          const CoordinateD& y)
{
  // TIMER_VERBOSE("makeProperRotation");
  // displayln(ssprintf("y=") + show(y));
  SpatialO3Matrix rot = makeRotationAroundX(0);
  SpatialO3Matrix rotx = makeRotationAroundX(0);
  SpatialO3Matrix rotz = makeRotationAroundX(0);
  SpatialO3Matrix xrotx = makeRotationAroundX(0);
  if (is_very_close(y[2], 0) && is_very_close(y[1], 0)) {
    rotx = makeRotationAroundX(0);
  } else {
    const RealD phi_x = std::atan2(y[2], y[1]);
    rotx = makeRotationAroundX(-phi_x);
  }
  const CoordinateD y1 = rotx * y;
  // displayln(ssprintf("y1=") + show(y1));
  Qassert(is_very_close(y1[2], 0));
  Qassert(y1[1] >= 0 || is_very_close(y1[1], 0));
  if (is_very_close(y1[1], 0) && is_very_close(y1[0], 0)) {
    rotz = makeRotationAroundZ(0);
  } else {
    const RealD phi_z = std::atan2(y1[1], y1[0]);
    rotz = makeRotationAroundZ(-phi_z);
  }
  const CoordinateD y2 = rotz * y1;
  // displayln(ssprintf("y2=") + show(y2));
  Qassert(is_very_close(y2[2], 0));
  Qassert(is_very_close(y2[1], 0));
  Qassert(y2[0] >= 0 || is_very_close(y2[0], 0));
  rot = rotz * rotx;
  const CoordinateD x2 = rot * x;
  if (is_very_close(x2[2], 0) && is_very_close(x2[1], 0)) {
    xrotx = makeRotationAroundX(0);
  } else {
    const RealD xphi_x = std::atan2(x2[2], x2[1]);
    xrotx = makeRotationAroundX(-xphi_x);
  }
  const CoordinateD x3 = xrotx * x2;
  // displayln(ssprintf("x3=") + show(x3));
  Qassert(is_very_close(x3[2], 0));
  Qassert(x3[1] >= 0 || is_very_close(x3[1], 0));
  rot = xrotx * rot;
  if (is_very_close(y2[0], 0) && is_very_close(y2[1], 0) &&
      is_very_close(y2[2], 0)) {
    // important
    // if y has no spatial component
    // then eta = 0
    // thus the spatial component of x need to be along y direction
    SpatialO3Matrix rotzz = makeRotationAroundZ(0);
    if (is_very_close(x3[1], 0) && is_very_close(x3[0], 0)) {
      rotzz = makeRotationAroundZ(0);
    } else {
      const RealD phi_z = std::atan2(x3[1], x3[0]);
      rotzz = makeRotationAroundZ(-phi_z + PI / 2.0);
    }
    const CoordinateD x4 = rotzz * x3;
    // displayln(ssprintf("y2=") + show(y2));
    Qassert(is_very_close(x4[2], 0));
    Qassert(is_very_close(x4[0], 0));
    Qassert(x4[1] >= 0 || is_very_close(x4[1], 0));
    rot = rotzz * rot;
  }
  return rot;
}

inline ManyMagneticMoments muonLine(const CoordinateD& x, const CoordinateD& y,
                                    const IntegrationEps& eps)
{
  TIMER("muonLine");
  std::vector<RealD> integral = integrateMuonLine(x, y, eps);
  return computeProjections(integral);
}

inline ManyMagneticMoments muonLineSymR(const CoordinateD& x,
                                        const CoordinateD& y,
                                        const IntegrationEps& eps)
{
  TIMER("muonLineSymR");
  // already performed in the integrand
  return muonLine(x, y, eps);
}

inline std::vector<RealD> std_vector_from_many_magnetic_moments(
    const ManyMagneticMoments& mmm)
// return ret;
//
// \mathcal M_{i,\rho,\sigma,\lambda}(x,y,z)
// = ret[i * 64 + rho * 16 + sigma * 4 + lambda]
// = mmm[16 * sigma + 4 * lambda + rho][i]
{
  std::vector<RealD> ret(3 * 4 * 4 * 4, 0.0);
  for (Int i = 0; i < 3; ++i) {
    for (Int rho = 0; rho < 4; ++rho) {
      for (Int sigma = 0; sigma < 4; ++sigma) {
        for (Int lambda = 0; lambda < 4; ++lambda) {
          ret[i * 64 + rho * 16 + sigma * 4 + lambda] =
              mmm[16 * sigma + 4 * lambda + rho][i];
        }
      }
    }
  }
  return ret;
}

// ------------------------------------------------------

inline ManyMagneticMoments muonLineSym(const CoordinateD& x,
                                       const CoordinateD& y,
                                       const IntegrationEps& eps)
// interface function
// important, return average instead of sum of all six permutations
// This function return the avg of six different muon-line diagram
// directly compute the result
//
{
  TIMER("muonLineSym");
  std::vector<ManyMagneticMoments> mmms(3);
  mmms[0] = muonLineSymR(x, y, eps);
  mmms[1] = permuteNuRhoMu(muonLineSymR(-x + y, -x, eps), 2, 0,
                           1);  // x',2 <- y,0 ; y',0 <- z,1 ; z',1 <- x,2
  mmms[2] = permuteNuRhoMu(muonLineSymR(-y, x - y, eps), 1, 2,
                           0);  // x',2 <- z,1 ; y',0 <- x,2 ; z',1 <- y,0
  return averageManyMagneticMoments(mmms);
}

inline std::vector<RealD> muon_line_sym_py(const CoordinateD& x,
                                           const CoordinateD& y,
                                           const IntegrationEps& eps)
// python interface function
// interface function
// important, return average instead of sum of all six permutations
// This function return the avg of six different muon-line diagram
// directly compute the result
//
// return ret;
//
// \mathcal M_{i,\rho,\sigma,\lambda}(x,y,z)
// = ret[i * 64 + rho * 16 + sigma * 4 + lambda]
{
  TIMER("muon_line_sym_py");
  const ManyMagneticMoments mmm = muonLineSym(x, y, eps);
  return std_vector_from_many_magnetic_moments(mmm);
}

// ------------------------------------------------------

inline void coordinatesFromParams(CoordinateD& x, CoordinateD& y,
                                  const std::vector<RealD>& params)
// if y[0] == 0 then theta = 0 and eta = PI/2
{
  assert(5 == params.size());
  const RealD d = DISTANCE_LIMIT * std::pow(params[0], 2.0);
  const RealD alpha = std::pow(params[1], 2.0);
  const RealD theta = params[2] * PI;
  const RealD phi = params[3] * PI;
  const RealD eta = params[4] * PI;
  x[0] = d * alpha * std::sin(phi) * std::cos(eta);
  x[1] = d * alpha * std::sin(phi) * std::sin(eta);
  x[2] = 0.0;
  x[3] = d * alpha * std::cos(phi);
  y[0] = d * std::sin(theta);
  y[1] = 0.0;
  y[2] = 0.0;
  y[3] = d * std::cos(theta);
  if (theta == PI || theta == 0.0) {
    y[0] = 0.0;
    x[0] = 0.0;
    x[1] = d * alpha * std::sin(phi);
  }
  if (phi == PI || phi == 0.0) {
    x[0] = 0.0;
    x[1] = 0.0;
  }
}

inline void paramsFromCoordinates(std::vector<RealD>& params,
                                  const CoordinateD& x, const CoordinateD& y)
{
  const RealD x_len = coordinate_len(x) + 1.0E-99;
  const RealD d = coordinate_len(y) + 1.0E-99;
  RealD alpha = x_len / d;
  if (alpha > 1.0) {
    Qassert(alpha < 1.0 + 1e-10);
    alpha = 1.0;
  }
  const RealD cos_theta = y[3] / d;
  const RealD cos_phi = x[3] / x_len;
  RealD cos_eta = (x[0] * y[0] + x[1] * y[1] + x[2] * y[2]) /
                   (std::sqrt((x[0] * x[0] + x[1] * x[1] + x[2] * x[2]) *
                              (y[0] * y[0] + y[1] * y[1] + y[2] * y[2])) +
                    1.0E-99);
  if (cos_eta > 1.0) {
    cos_eta = 1.0;
  } else if (cos_eta < -1.0) {
    cos_eta = -1.0;
  }
  params.resize(5);
  params[0] = std::pow(d / DISTANCE_LIMIT, 1.0 / 2.0);
  params[1] = std::pow(alpha, 1.0 / 2.0);
  params[2] = std::acos(cos_theta) / PI;
  params[3] = std::acos(cos_phi) / PI;
  params[4] = std::acos(cos_eta) / PI;
  if (qisnan(params)) {
    displayln_c_stdout(shows("paramsFromCoordinates ") + show(x));
    displayln_c_stdout(shows("paramsFromCoordinates ") + show(y));
    displayln_c_stdout(shows("paramsFromCoordinates ") + show(params));
    Qassert(false);
  }
}

inline ManyMagneticMomentsCompressed compressManyMagneticMoments(
    const ManyMagneticMoments& m)
{
  ManyMagneticMomentsCompressed ret;
  size_t index = 0;
  bool bi, bj, bk, bl;
  for (Int i = 0; i < 4; ++i) {
    bi = 2 == i;
    for (Int j = 0; j < 4; ++j) {
      bj = bi ^ (2 == j);
      for (Int k = 0; k < 4; ++k) {
        bk = bj ^ (2 == k);
        for (Int l = 0; l < 3; ++l) {
          bl = bk ^ (2 == l);
          if (bl) {
            ret[index] = m[i * 4 * 4 + j * 4 + k][l];
            index += 1;
          }
        }
      }
    }
  }
  assert(ret.size() == index);
  return ret;
}

inline ManyMagneticMoments uncompressManyMagneticMoments(
    const ManyMagneticMomentsCompressed& mc)
{
  ManyMagneticMoments ret;
  set_zero(ret);
  size_t index = 0;
  bool bi, bj, bk, bl;
  for (Int i = 0; i < 4; ++i) {
    bi = 2 == i;
    for (Int j = 0; j < 4; ++j) {
      bj = bi ^ (2 == j);
      for (Int k = 0; k < 4; ++k) {
        bk = bj ^ (2 == k);
        for (Int l = 0; l < 3; ++l) {
          bl = bk ^ (2 == l);
          if (bl) {
            ret[i * 4 * 4 + j * 4 + k][l] = mc[index];
            index += 1;
          }
        }
      }
    }
  }
  assert(mc.size() == index);
  return ret;
}

inline ManyMagneticMomentsCompressed muonLineSymParamsCompressed(
    const std::vector<RealD>& params, const IntegrationEps& eps)
// Target function for interpolation
{
  TIMER("muonLineSymParamsCompressed");
  CoordinateD x, y;
  coordinatesFromParams(x, y, params);
  const ManyMagneticMoments mmm = muonLineSym(x, y, eps);
  assert(false == qisnan(mmm));
  const ManyMagneticMomentsCompressed ret = compressManyMagneticMoments(mmm);
  // const RealD mmm_len = std::sqrt(qnorm(ret));
  return ret;
}

typedef InterpolationBilinearNd<ManyMagneticMomentsCompressed> MuonLineInterp;

API inline std::vector<MuonLineInterp>& get_muonline_interps()
{
  static std::vector<MuonLineInterp> interps(1);
  return interps;
}

API inline Long& get_default_muonline_interp_idx()
{
  static Long idx = IS_USING_MUON_LINE_INTERPOLATION ? 0 : -1;
  return idx;
}

inline MuonLineInterp& getMuonLineInterp(
    const Long idx = get_default_muonline_interp_idx())
{
  Qassert(idx >= 0);
  std::vector<MuonLineInterp>& interps = get_muonline_interps();
  if (idx >= (Long)interps.size()) {
    interps.resize(idx + 1);
  }
  return interps[idx];
}

inline ManyMagneticMomentsCompressed muonLineSymParamsCompressedInterpolate(
    const std::vector<RealD>& params,
    const Int b_interp = get_default_muonline_interp_idx())
{
  const MuonLineInterp& interpolation = getMuonLineInterp(b_interp);
  // displayln(ssprintf("muonLineSymParamsCompressedInterpolate: ") +
  // show(params) + " " + show(qnorm(interpolation(params))));
  return interpolation(params);
}

inline ManyMagneticMoments muonLineSymParams(
    const std::vector<RealD>& params, const IntegrationEps& eps,
    const Int b_interp = get_default_muonline_interp_idx())
{
  ManyMagneticMomentsCompressed mmm;
  if (b_interp >= 0) {
    mmm = muonLineSymParamsCompressedInterpolate(params, b_interp);
  } else {
    mmm = muonLineSymParamsCompressed(params, eps);
  }
  return uncompressManyMagneticMoments(mmm);
}

inline ManyMagneticMoments muonLineSymThroughParam(
    const CoordinateD& x, const CoordinateD& y, const IntegrationEps& eps,
    const Int b_interp = get_default_muonline_interp_idx())
{
  std::vector<RealD> params;
  paramsFromCoordinates(params, x, y);
  return muonLineSymParams(params, eps, b_interp);
}

inline MagneticMoment operator*(const SpatialO3Matrix& m,
                                const MagneticMoment& v)
{
  Eigen::Matrix<RealD, 3, 1> vec;
  vec << v[0], v[1], v[2];
  MagneticMoment ret;
  vec = m.determinant() * m * vec;
  ret[0] = vec[0];
  ret[1] = vec[1];
  ret[2] = vec[2];
  return ret;
}

inline ManyMagneticMoments operator*(const SpatialO3Matrix& m,
                                     const ManyMagneticMoments& v)
{
  // TIMER_VERBOSE("SpatialO3Matrix*ManyMagneticMoments");
  ManyMagneticMoments ret, tmp;
  for (Int i = 0; i < 4 * 4 * 4; ++i) {
    ret[i] = m * v[i];
  }
  // displayln(show(sqrt(qnorm(ret))));
  tmp = ret;
  for (Int i = 0; i < 4; ++i) {
    for (Int j = 0; j < 4; ++j) {
      const Int prefix = i * 4 * 4 + j * 4;
      for (Int k = 0; k < 3; ++k) {
        MagneticMoment& mm = ret[prefix + k];
        set_zero(mm);
        for (Int kk = 0; kk < 3; ++kk) {
          mm = mm + m(k, kk) * tmp[prefix + kk];
        }
      }
    }
  }
  // displayln(show(sqrt(qnorm(ret))));
  tmp = ret;
  for (Int i = 0; i < 4; ++i) {
    for (Int j = 0; j < 4; ++j) {
      const Int prefix = i * 4 * 4 + j;
      for (Int k = 0; k < 3; ++k) {
        MagneticMoment& mm = ret[prefix + k * 4];
        set_zero(mm);
        for (Int kk = 0; kk < 3; ++kk) {
          mm = mm + m(k, kk) * tmp[prefix + kk * 4];
        }
      }
    }
  }
  // displayln(show(sqrt(qnorm(ret))));
  tmp = ret;
  for (Int i = 0; i < 4; ++i) {
    for (Int j = 0; j < 4; ++j) {
      const Int prefix = i * 4 + j;
      for (Int k = 0; k < 3; ++k) {
        MagneticMoment& mm = ret[prefix + k * 4 * 4];
        set_zero(mm);
        for (Int kk = 0; kk < 3; ++kk) {
          mm = mm + m(k, kk) * tmp[prefix + kk * 4 * 4];
        }
      }
    }
  }
  // displayln(show(sqrt(qnorm(ret))));
  return ret;
}

inline ManyMagneticMoments muonLineSymRotate(
    const CoordinateD& x, const CoordinateD& y, const IntegrationEps& eps,
    const Int b_interp = get_default_muonline_interp_idx())
{
  const SpatialO3Matrix rot = makeProperRotation(x, y);
  CoordinateD xr = rot * x;
  CoordinateD yr = rot * y;
  // for (Int i = 0; i < 4; ++i) {
  //   if (is_very_close(xr[i], 0)) {
  //     xr[i] = 0.0;
  //   }
  //   if (is_very_close(yr[i], 0)) {
  //     yr[i] = 0.0;
  //   }
  // }
  const SpatialO3Matrix rott = rot.transpose();
  return rott * muonLineSymThroughParam(xr, yr, eps, b_interp);
}

inline ManyMagneticMoments muonLineSymPermute(
    const CoordinateD& x, const CoordinateD& y, const IntegrationEps& eps,
    const Int b_interp = get_default_muonline_interp_idx())
{
  const RealD xyl = coordinate_len(y - x);
  RealD xl = coordinate_len(x);
  RealD yl = coordinate_len(y);
  if (is_very_close(xl, xyl)) {
    xl = xyl;
  }
  if (is_very_close(yl, xyl)) {
    yl = xyl;
  }
  if (is_very_close(yl, xl)) {
    yl = xl;
  }
  std::vector<ManyMagneticMoments> mmms;
  if (xl <= xyl && xyl <= yl) {
    mmms.push_back(muonLineSymRotate(x, y, eps, b_interp));  // y z x
  }
  if (xyl <= yl && yl <= xl) {
    mmms.push_back(
        permuteNuRhoMu(muonLineSymRotate(-x + y, -x, eps, b_interp),
                       2, 0, 1));  // z x y
  }
  if (yl <= xl && xl <= xyl) {
    mmms.push_back(
        permuteNuRhoMu(muonLineSymRotate(-y, x - y, eps, b_interp),
                       1, 2, 0));  // x y z
  }
  if (yl <= xyl && xyl <= xl) {
    mmms.push_back(permuteNuRhoMu(
        muonLineSymRotate(y, x, eps, b_interp), 2, 1, 0));  // x z y
  }
  if (xl <= yl && yl <= xyl) {
    mmms.push_back(
        permuteNuRhoMu(muonLineSymRotate(-x, -x + y, eps, b_interp),
                       0, 2, 1));  // y x z
  }
  if (xyl <= xl && xl <= yl) {
    mmms.push_back(
        permuteNuRhoMu(muonLineSymRotate(x - y, -y, eps, b_interp),
                       1, 0, 2));  // z y x
  }
  Qassert(mmms.size() > 0);
  return averageManyMagneticMoments(mmms);
}

inline ManyMagneticMoments muonLineSymTransform(
    const CoordinateD& x, const CoordinateD& y, const IntegrationEps& eps,
    const Int b_interp = get_default_muonline_interp_idx())
// interface function
// This function return the avg of six different muon-line diagram
// proper rotation transformation is made in order to reduce to 5 parameters
// interpolation is performed if IS_USING_MUON_LINE_INTERPOLATION = true
//
// ManyMagneticMoments format [y-pol,z-pol,x-pol][mag-dir]
// argument x and y assume z = 0
{
  return muonLineSymPermute(x, y, eps, b_interp);
}

inline void paramsFromCoordinatesPermute(std::vector<RealD>& params,
                                         const CoordinateD& x,
                                         const CoordinateD& y)
{
  const RealD xl = coordinate_len(x);
  const RealD yl = coordinate_len(y);
  const RealD xyl = coordinate_len(y - x);
  if (xl <= xyl && xyl <= yl) {
    paramsFromCoordinates(params, x, y);
  } else if (xyl <= yl && yl <= xl) {
    paramsFromCoordinates(params, -x + y, -x);
  } else if (yl <= xl && xl <= xyl) {
    paramsFromCoordinates(params, -y, x - y);
  } else if (yl <= xyl && xyl <= xl) {
    paramsFromCoordinates(params, y, x);
  } else if (xl <= yl && yl <= xyl) {
    paramsFromCoordinates(params, -x, -x + y);
  } else if (xyl <= xl && xl <= yl) {
    paramsFromCoordinates(params, x - y, -y);
  } else {
    assert(false);
  }
}

inline void compare_many_magnetic_moments(const std::string& tag,
                                          const CoordinateD& x,
                                          const CoordinateD& y,
                                          const ManyMagneticMoments& mmm,
                                          const ManyMagneticMoments& mmmp)
{
  std::vector<RealD> params;
  paramsFromCoordinatesPermute(params, x, y);
  const RealD diff_percent = 100.0 * sqrt(qnorm(mmmp - mmm) / qnorm(mmm));
  const bool is_print = diff_percent > 0.0001 || true;
  if (is_print) {
#pragma omp critical
    {
      displayln_c_stdout(compare_multiline_string(
          showManyMagneticMoments(mmm), showManyMagneticMoments(mmmp), 48));
      displayln_c_stdout(tag + ": " +
                         ssprintf("CHECKING: %10.2e %10.2e %10.4f%%",
                                  sqrt(qnorm(mmm)), sqrt(qnorm(mmmp - mmm)),
                                  diff_percent));
      displayln_c_stdout(tag + ": " + shows("params= ") + show(params));
      displayln_c_stdout(
          tag + ": " +
          ssprintf(" x  = %8.4f %s", coordinate_len(x), show(x).c_str()));
      displayln_c_stdout(
          tag + ": " +
          ssprintf(" y  = %8.4f %s", coordinate_len(y), show(y).c_str()));
      displayln_c_stdout(tag + ": " +
                         ssprintf("y-x = %8.4f %s", coordinate_len(y - x),
                                  show(y - x).c_str()));
      displayln_c_stdout(tag + ": " + ssprintf("x  y  ") +
                         show(sqrt(qnorm(mmm))));
      displayln_c_stdout(
          tag + ": " +
          ssprintf("DATA: %24.17E %24.17E %24.17E %24.17E %24.17E   "
                   "%24.17E %24.17E %24.17E  %24.17E",
                   params[0], params[1], params[2], params[3], params[4],
                   sqrt(qnorm(mmm)), sqrt(qnorm(mmmp)), sqrt(qnorm(mmmp - mmm)),
                   sqrt(qnorm(mmmp - mmm) / qnorm(mmm))));
    }
  }
}

inline ManyMagneticMoments muonLineSymParamsCheck(const CoordinateD& x,
                                                  const CoordinateD& y,
                                                  const IntegrationEps& eps)
{
  TIMER_VERBOSE("muonLineSymParamsCheck");
  ManyMagneticMoments mmm = muonLineSym(x, y, eps);
  ManyMagneticMoments mmmp = muonLineSymTransform(x, y, eps);
  compare_many_magnetic_moments("params-check", x, y, mmm, mmmp);
  return mmm;
}

inline ManyMagneticMoments muonLineSymRotateCheck(const SpatialO3Matrix& rot,
                                                  const CoordinateD& x,
                                                  const CoordinateD& y,
                                                  const IntegrationEps& eps)
{
  TIMER_VERBOSE("muonLineSymRotateCheck");
  ManyMagneticMoments mmm = muonLineSymTransform(x, y, eps);
  const CoordinateD xr = rot * x;
  const CoordinateD yr = rot * y;
  const SpatialO3Matrix rott = rot.transpose();
  ManyMagneticMoments mmmp = rott * muonLineSymTransform(xr, yr, eps);
  compare_many_magnetic_moments("rotate-check", x, y, mmm, mmmp);
  return mmm;
}

inline void initializeMuonLineInterpolation(const std::vector<Int>& dims,
                                            const IntegrationEps& eps)
// computing the muon-line interpolation database
// take quite some time
{
  TIMER_VERBOSE("initializeMuonLineInterpolation");
  MuonLineInterp& interpolation = getMuonLineInterp();
  interpolation.init();
  assert(dims.size() == 5);
  interpolation.add_dimension(dims[0], 1.0, 0.0);  // d
  interpolation.add_dimension(dims[1], 1.0, 0.0);  // alpha
  interpolation.add_dimension(dims[2], 1.0, 0.0);  // theta
  interpolation.add_dimension(dims[3], 1.0, 0.0);  // phi
  interpolation.add_dimension(dims[4], 1.0, 0.0);  // eta
  const Long jobs_total = interpolation.size();
  const Long num_nodes = get_num_node();
  const Long jobs_per_nodes = jobs_total / num_nodes;
  const Long jobs_parallel = jobs_per_nodes * num_nodes;
  const Long jobs_left = jobs_total - jobs_parallel;
  const Long my_start = jobs_left + jobs_per_nodes * get_id_node();
  std::vector<ManyMagneticMomentsCompressed> workplace(jobs_per_nodes);
#pragma omp parallel for schedule(dynamic)
  for (Long i = 0; i < jobs_per_nodes; ++i) {
    const Long idx = my_start + i;
    TIMER_VERBOSE("interp-initial-iter");
#pragma omp critical
    {
      displayln_c_stdout(ssprintf("jobs-par: %5d %10ld %10ld/%ld",
                                  get_id_node(), my_start, i, jobs_per_nodes));
    }
    workplace[i] =
        muonLineSymParamsCompressed(interpolation.get_coor(idx), eps);
  }
  Vector<ManyMagneticMomentsCompressed> all_workspace(&interpolation[jobs_left],
                                                      jobs_parallel);
  all_gather(all_workspace, Vector<ManyMagneticMomentsCompressed>(workplace));
#pragma omp parallel for schedule(dynamic)
  for (Long idx = 0; idx < jobs_left; ++idx) {
    TIMER_VERBOSE("interp-initial-iter");
#pragma omp critical
    if (0 == get_id_node()) {
      displayln_c_stdout(ssprintf("jobs-left: %10ld/%ld", idx, jobs_left));
      displayln_c_stdout(show(idx));
    }
    interpolation[idx] =
        muonLineSymParamsCompressed(interpolation.get_coor(idx), eps);
  }
}

inline void saveMuonLineInterpolation(const std::string& path)
{
  TIMER_VERBOSE("saveMuonLineInterpolation");
  if (0 == get_id_node()) {
    MuonLineInterp& interp = getMuonLineInterp();
    qmkdir(path);
    FILE* fdims = qopen(path + "/dims.txt" + ".partial", "w");
    const std::vector<InterpolationDim>& dims = interp.dims;
    fprintf(fdims, "# i dims[i].n dims[i].xhigh dims[i].xlow\n");
    for (size_t i = 0; i < dims.size(); ++i) {
      fprintf(fdims, "%3ld %5d %24.17E %24.17E\n", i, dims[i].n, dims[i].xhigh,
              dims[i].xlow);
    }
    qfclose(fdims);
    qrename_partial(path + "/dims.txt");
    FILE* fdata = qopen(path + "/data.txt" + ".partial", "w");
    const std::vector<ManyMagneticMomentsCompressed>& data = interp.data;
    fprintf(fdata, "# idx params[0-4] ManyMagneticMomentsCompressed[0-91]\n");
    for (size_t idx = 0; idx < data.size(); ++idx) {
      fprintf(fdata, "%10ld", idx);
      const std::vector<RealD> params = interp.get_coor(idx);
      for (size_t i = 0; i < params.size(); ++i) {
        fprintf(fdata, " %24.17E", params[i]);
      }
      fprintf(fdata, "  ");
      const ManyMagneticMomentsCompressed& mmm = data[idx];
      for (size_t i = 0; i < mmm.size(); ++i) {
        fprintf(fdata, " %24.17E", mmm[i]);
      }
      fprintf(fdata, "\n");
    }
    qfclose(fdata);
    qrename_partial(path + "/data.txt");
    qtouch(path + "/checkpoint");
  }
  SYNC_NODE();
}

inline bool loadMuonLineInterpolation(const std::string& path,
                                      const size_t interp_idx = 0)
{
  TIMER_VERBOSE("loadMuonLineInterpolation");
  displayln_info(fname +
                 ssprintf(": load '%s' as idx=%ld", path.c_str(), interp_idx));
  if (!does_file_exist_qar_sync_node(path + "/checkpoint")) {
    return false;
  }
  MuonLineInterp& interp = getMuonLineInterp(interp_idx);
  interp.init();
  std::vector<std::vector<RealD>> dims;
  if (0 == get_id_node()) {
    Qassert(does_file_exist_qar(path + "/dims.txt"));
    dims = qload_datatable(path + "/dims.txt");
    for (size_t i = 0; i < dims.size(); ++i) {
      Qassert(i == (size_t)dims[i][0]);
      interp.add_dimension((int)dims[i][1], dims[i][2], dims[i][3]);
    }
  }
  // sync MuonLineInterp dims across all nodes
  bcast(dims);
  if (0 != get_id_node()) {
    for (size_t i = 0; i < dims.size(); ++i) {
      Qassert(i == (size_t)dims[i][0]);
      interp.add_dimension((int)dims[i][1], dims[i][2], dims[i][3]);
    }
  }
  Long limit = 0;
  if (0 == get_id_node()) {
    // Can be obtained from
    // split -da 10 -l 10000 data.txt data.txt.
    while (does_file_exist_qar(path + ssprintf("/data.txt.%010d", limit))) {
      limit += 1;
    }
  }
  glb_sum(limit);
  std::vector<DataTable> tables(limit);
  const Int num_node = get_num_node();
  const Int id_node = get_id_node();
  for (size_t i = id_node; i < (size_t)limit; i += num_node) {
    std::string fn = path + ssprintf("/data.txt.%010d", i);
    Qassert(does_file_exist_qar(fn));
    tables[i] = qload_datatable_par(fn);
  }
  if (limit == 0 && get_id_node() == 0) {
    Qassert(does_file_exist_qar(path + "/data.txt"));
    tables.push_back(qload_datatable_par(path + "/data.txt"));
  }
  DataTable data;
  for (size_t i = 0; i < tables.size(); ++i) {
    for (size_t j = 0; j < tables[i].size(); ++j) {
      data.push_back(tables[i][j]);
    }
  }
  set_zero(get_data(interp.data));
#pragma omp parallel for
  for (size_t i = 0; i < data.size(); ++i) {
    const std::vector<RealD>& data_vec = data[i];
    qassert(data_vec.size() == 1 + 5 + 92);
    const size_t idx = (size_t)data_vec[0];
    std::vector<RealD> params = interp.get_coor(idx);
    qassert(params.size() == 5);
    for (size_t k = 0; k < params.size(); ++k) {
      qassert(params[k] == data_vec[1 + k]);
    }
    ManyMagneticMomentsCompressed mmmc;
    qassert(mmmc.size() == 92);
    for (size_t k = 0; k < mmmc.size(); ++k) {
      mmmc[k] = data_vec[6 + k];
    }
    interp[idx] = mmmc;
  }
  // sync MuonLineInterp data across all nodes
  glb_sum(get_data(interp.data));
  return true;
}

inline void load_or_compute_muonline_interpolation(const std::string& path,
                                                   const std::vector<Int>& dims,
                                                   const IntegrationEps& eps)
{
  if (!loadMuonLineInterpolation(path)) {
    test_fCalc();
    initializeMuonLineInterpolation(dims, eps);
    saveMuonLineInterpolation(path);
  }
}

inline void save_part_muonline_interpolation_data(
    const std::string& path, const size_t fn_idx, const size_t start_idx,
    const Vector<ManyMagneticMomentsCompressed> data)
{
  TIMER_VERBOSE("save_part_muonline_interpolation_data");
  MuonLineInterp& interp = getMuonLineInterp();
  std::string fn = path + ssprintf("/data.txt.%010d", fn_idx);
  FILE* fdata = qopen(fn + ".partial", "w");
  fprintf(fdata, "# idx params[0-4] ManyMagneticMomentsCompressed[0-91]\n");
  for (size_t k = 0; k < (size_t)data.size(); ++k) {
    size_t idx = start_idx + k;
    fprintf(fdata, "%10ld", idx);
    const std::vector<RealD> params = interp.get_coor(idx);
    for (size_t i = 0; i < params.size(); ++i) {
      fprintf(fdata, " %24.17E", params[i]);
    }
    fprintf(fdata, "  ");
    const ManyMagneticMomentsCompressed& mmm = data[k];
    for (size_t i = 0; i < mmm.size(); ++i) {
      fprintf(fdata, " %24.17E", mmm[i]);
    }
    fprintf(fdata, "\n");
  }
  qfclose(fdata);
  qrename_partial(fn);
}

inline bool is_part_muonline_interpolation_data_done(const std::string& path,
                                                     const size_t fn_idx)
{
  std::string fn = path + ssprintf("/data.txt.%010d", fn_idx);
  return does_file_exist_qar(fn);
}

struct PointPairWeight {
  CoordinateD rxy, rxz;
  RealD weight;
};

inline Int boundary_multiplicity(const Coordinate& xg, const Coordinate& lb,
                                 const Coordinate& ub)
{
  Int ret = 1;
  for (Int i = 0; i < 4; ++i) {
    if (xg[i] == lb[i] || xg[i] == ub[i]) {
      ret *= 2;
    }
  }
  return ret;
}

inline std::vector<PointPairWeight> shift_lat_corr(const Coordinate& x,
                                                   const Coordinate& y,
                                                   const Coordinate& z,
                                                   const Coordinate& total_site,
                                                   const RealD a)
// x and y are closest
// rxy and rxz are returned
// z must be (1) within the box centered at x (2) within the box centered at y
// if at boundary of the box, need to be divided by the multiplicity
{
  const Coordinate rxy0 = relative_coordinate(y - x, total_site);
  std::vector<Coordinate> rxy_list;
  for (Int r0 = rxy0[0]; r0 <= total_site[0] / 2; r0 += total_site[0]) {
    for (Int r1 = rxy0[1]; r1 <= total_site[1] / 2; r1 += total_site[1]) {
      for (Int r2 = rxy0[2]; r2 <= total_site[2] / 2; r2 += total_site[2]) {
        for (Int r3 = rxy0[3]; r3 <= total_site[3] / 2; r3 += total_site[3]) {
          rxy_list.push_back(Coordinate(r0, r1, r2, r3));
        }
      }
    }
  }
  Qassert(rxy_list[0] == rxy0);
  std::vector<PointPairWeight> ret;
  for (Int i = 0; i < (int)rxy_list.size(); ++i) {
    const Coordinate rxy = rxy_list[i];
    const CoordinateD rdxy = a * CoordinateD(rxy);
    Coordinate lb, ub;
    for (Int m = 0; m < 4; ++m) {
      lb[m] = std::max(0, rxy[m]) - total_site[m] / 2;
      ub[m] = std::min(0, rxy[m]) + total_site[m] / 2;
    }
    Coordinate rxz = relative_coordinate(z - x, total_site);
    for (Int r0 = rxz[0]; r0 <= ub[0]; r0 += total_site[0]) {
      if (lb[0] > r0) {
        continue;
      }
      for (Int r1 = rxz[1]; r1 <= ub[1]; r1 += total_site[1]) {
        if (lb[1] > r1) {
          continue;
        }
        for (Int r2 = rxz[2]; r2 <= ub[2]; r2 += total_site[2]) {
          if (lb[2] > r2) {
            continue;
          }
          for (Int r3 = rxz[3]; r3 <= ub[3]; r3 += total_site[3]) {
            if (lb[3] > r3) {
              continue;
            }
            PointPairWeight ppw;
            ppw.rxy = rdxy;
            const Coordinate rxz(r0, r1, r2, r3);
            ppw.rxz = a * CoordinateD(rxz);
            const RealD l1 = coordinate_len(ppw.rxz);
            const RealD l2 = coordinate_len(ppw.rxz - ppw.rxy);
            const RealD l3 = coordinate_len(ppw.rxy);
            if (l1 >= DISTANCE_LIMIT - 0.01 or l2 >= DISTANCE_LIMIT - 0.01 or
                l3 >= DISTANCE_LIMIT - 0.01) {
              // one edge is too long, outside of the interpolation range
              continue;
            }
            ppw.weight = 1.0 / (RealD)rxy_list.size() /
                         (RealD)boundary_multiplicity(rxz, lb, ub);
            ret.push_back(ppw);
          }
        }
      }
    }
  }
  return ret;
}

// ------------------------------------------------------

inline void clear_muon_line_interpolations()
{
  get_muonline_interps().resize(1);
  getMuonLineInterp(0).init();
}

inline Long get_number_of_muon_line_interpolations()
{
  const Long size = get_muonline_interps().size();
  Qassert(size >= 1);
  Long count = 0;
  for (Long i = 0; i < size; ++i) {
    if (not getMuonLineInterp(i).empty()) {
      count += 1;
    }
  }
  return count;
}

inline bool compute_save_muonline_interpolation_cc(const std::string& path,
                                                   const std::vector<Int>& dims,
                                                   const IntegrationEps& eps)
// preferred way to generate interpolation field
{
  TIMER_VERBOSE("compute_save_muonline_interpolation_cc");
  std::vector<std::vector<RealD>> ldims;
  if (0 == get_id_node() && does_file_exist_qar(path + "/dims.txt")) {
    std::vector<std::vector<RealD>> ldims =
        qload_datatable(path + "/dims.txt");
    Qassert(ldims.size() == dims.size());
    for (Long i = 0; i < (Long)ldims.size(); ++i) {
      Qassert(i == (Long)ldims[i][0]);
      Qassert(dims[i] == (int)ldims[i][1]);
      Qassert(1.0 == ldims[i][2]);
      Qassert(0.0 == ldims[i][3]);
    }
  }
  assert(dims.size() == 5);
  MuonLineInterp& interp = getMuonLineInterp();
  interp.init();
  for (Long i = 0; i < (Long)dims.size(); ++i) {
    interp.add_dimension(dims[i], 1.0, 0.0);
  }
  set_zero(get_data(interp.data));
  if (!does_file_exist_qar_sync_node(path + "/checkpoint")) {
    qmkdir_info(path);
    if (!obtain_lock(path + "/lock")) {
      return false;
    }
    if (0 == get_id_node() && !does_file_exist_qar(path + "/dims.txt")) {
      FILE* fdims = qopen(path + "/dims.txt" + ".partial", "w");
      const std::vector<InterpolationDim>& dims = interp.dims;
      fprintf(fdims, "# i dims[i].n dims[i].xhigh dims[i].xlow\n");
      for (Long i = 0; i < (Long)dims.size(); ++i) {
        fprintf(fdims, "%3ld %5d %24.17E %24.17E\n", (long)i, dims[i].n,
                dims[i].xhigh, dims[i].xlow);
      }
      qfclose(fdims);
      qrename_partial(path + "/dims.txt");
    }
    test_fCalc();
    const Long jobs_total = interp.size();
    displayln_info(ssprintf("jobs_total=%ld", jobs_total));
    // ADJUST ME
    const Long job_chunk_size = 64;
    std::vector<Long> jobs;
    for (Long start = 0; start <= jobs_total - job_chunk_size;
         start += job_chunk_size) {
      jobs.push_back(start);
    }
    displayln_info(ssprintf("jobs.size()=%ld", jobs.size()));
    std::vector<array<ManyMagneticMomentsCompressed, job_chunk_size>> data;
    Qassert(get_num_node() > 1);
    if (0 == get_id_node()) {
      const Long size = jobs.size();
      data.resize(size);
      Long num_running_jobs = 0;
      Long idx = 0;
      for (Int dest = 1; dest < get_num_node(); ++dest) {
        while (is_part_muonline_interpolation_data_done(path, idx)) {
          idx += 1;
        }
        if (idx >= (Long)jobs.size()) {
          break;
        }
        displayln_info(
            ssprintf("send-job: %5d %10ld/%ld", dest, jobs[idx], jobs_total));
        send_job(idx, jobs[idx], dest);
        idx += 1;
        num_running_jobs += 1;
      }
      while (idx < (Long)jobs.size()) {
        check_time_limit();
        TIMER_VERBOSE("muonline_interpolation_cc-traj");
        Int source;
        int64_t flag;
        array<ManyMagneticMomentsCompressed, job_chunk_size> result;
        receive_result(source, flag, result);
        num_running_jobs -= 1;
        data[flag] = result;
        save_part_muonline_interpolation_data(path, flag, jobs[flag],
                                              get_data(result));
        while (is_part_muonline_interpolation_data_done(path, idx)) {
          idx += 1;
        }
        if (idx >= (Long)jobs.size()) {
          break;
        }
        displayln_info(
            ssprintf("send-job: %5d %10ld/%ld", source, jobs[idx], jobs_total));
        send_job(idx, jobs[idx], source);
        idx += 1;
        num_running_jobs += 1;
      }
      while (num_running_jobs > 0) {
        check_time_limit();
        TIMER_VERBOSE("muonline_interpolation_cc-traj");
        Int source;
        int64_t flag;
        array<ManyMagneticMomentsCompressed, job_chunk_size> result;
        receive_result(source, flag, result);
        data[flag] = result;
        save_part_muonline_interpolation_data(path, flag, jobs[flag],
                                              get_data(result));
        num_running_jobs -= 1;
      }
      for (Int dest = 1; dest < get_num_node(); ++dest) {
        Long job = 0;
        displayln_info(
            ssprintf("send-job: %5d %10ld/%ld", dest, (Long)-1, jobs_total));
        send_job(-1, job, dest);
      }
      idx = 0;
      for (Long i = 0; i < (Long)data.size(); ++i) {
        for (Long j = 0; j < job_chunk_size; ++j) {
          interp[idx] = data[i][j];
          idx += 1;
        }
      }
    } else {
      while (true) {
        int64_t flag;
        Long job;
        receive_job(flag, job);
        if (-1 == flag) {
          displayln(ssprintf("par: %5d done", get_id_node()));
          break;
        }
        array<ManyMagneticMomentsCompressed, job_chunk_size> result;
#pragma omp parallel for schedule(dynamic)
        for (Long i = 0; i < job_chunk_size; ++i) {
          const Long idx = job + i;
          result[i] =
              muonLineSymParamsCompressed(interp.get_coor(idx), eps);
#pragma omp critical
          {
            displayln_c_stdout(ssprintf("par: %5d %10ld/%ld %10ld/%ld",
                                        get_id_node(), idx, jobs_total, i,
                                        job_chunk_size));
          }
        }
        send_result(flag, result);
      }
    }
    if (0 == get_id_node()) {
      const Long last_flag = jobs_total / job_chunk_size;
      const Long last_start_idx = last_flag * job_chunk_size;
#pragma omp parallel for schedule(dynamic)
      for (Long idx = last_start_idx; idx < jobs_total; ++idx) {
        interp[idx] =
            muonLineSymParamsCompressed(interp.get_coor(idx), eps);
#pragma omp critical
        {
          displayln_c_stdout(
              ssprintf("compute_save_muonline_interpolation_cc: %10ld/%ld", idx,
                       jobs_total));
        }
      }
      save_part_muonline_interpolation_data(
          path, last_flag, last_start_idx,
          Vector<ManyMagneticMomentsCompressed>(&interp[last_start_idx],
                                                jobs_total - last_start_idx));
      qtouch(path + "/checkpoint");
    }
    release_lock();
  }
  // ADJUST ME
  return true;
  // return loadMuonLineInterpolation(path);
}

inline bool load_multiple_muonline_interpolations(
    const std::string& path, std::vector<Long> idx_list = std::vector<Long>())
// python interface function
// interface function
// used to load the data
{
  TIMER_VERBOSE("load_multiple_muonline_interpolations");
  Long limit = 0;
  if (0 == get_id_node()) {
    while (does_file_exist_qar(path + ssprintf("/%010d/checkpoint", limit))) {
      limit += 1;
    }
  }
  glb_sum(limit);
  displayln_info(fname + ssprintf(": limit=%ld", (long)limit));
  if (0 == limit) {
    return loadMuonLineInterpolation(path);
  }
  if (idx_list.size() == 0) {
    for (Long idx = 0; idx < limit; ++idx) {
      idx_list.push_back(idx);
    }
  }
  for (Long i = 0; i < (Long)idx_list.size(); ++i) {
    const Long idx = idx_list[i];
    loadMuonLineInterpolation(path + ssprintf("/%010d", idx), idx);
  }
  return true;
}

inline ManyMagneticMoments get_muon_line_m(
    const CoordinateD& x, const CoordinateD& y, const CoordinateD& z,
    const Int idx = get_default_muonline_interp_idx(),
    const IntegrationEps& eps = IntegrationEps())
{
  // ADJUST ME
  return muonLineSymTransform(x - z, y - z, eps, idx);
  // return muonLineSym(x - z, y - z, eps);
}

inline std::vector<RealD> get_muon_line_m_py(
    const CoordinateD& x, const CoordinateD& y, const CoordinateD& z,
    const Int idx = get_default_muonline_interp_idx(),
    const IntegrationEps& eps = IntegrationEps())
{
  const ManyMagneticMoments mmm = get_muon_line_m(x, y, z, idx, eps);
  return std_vector_from_many_magnetic_moments(mmm);
}

inline std::vector<std::vector<RealD>> get_muon_line_m_extra_weights_default()
{
  std::vector<std::vector<RealD>> ret;
  ret.resize(3);
  ret[0].resize(6, 0.0);
  ret[0][5] = 3.0476190476190476;
  ret[0][3] = -2.3142857142857143;
  ret[0][1] = 0.26666666666666667;
  ret[1].resize(12, 0.0);
  ret[1][11] = 3.0476190476190476;
  ret[1][9] = -2.3142857142857143;
  ret[1][7] = 0.26666666666666667;
  ret[2].resize(1, 0.0);
  ret[2][0] = 1.0;
  return ret;
}

API inline std::vector<std::vector<RealD>>& get_muon_line_m_extra_weights()
// interface function
// weight[tag][idx]
// ->
// weight for muon-line interpolation of idx for tag in `get_muon_line_m_extra`
{
  static std::vector<std::vector<RealD>> weight =
      get_muon_line_m_extra_weights_default();
  return weight;
}

inline void set_muon_line_m_extra_weights(
    const std::vector<std::vector<RealD>>& weights)
{
  get_muon_line_m_extra_weights() = weights;
}

inline ManyMagneticMoments get_muon_line_m_extra(const CoordinateD& x,
                                                 const CoordinateD& y,
                                                 const CoordinateD& z,
                                                 const Int tag)
// % \mathcal M_{i,\rho,\sigma,\lambda}(x,y,z)
// % interpolate saved data and extrapolate different interpolations
// % = get_muon_line_m_extra(m * x, m * y, m * z , tag)[16 * sigma  + 4 * lambda + rho][i]
// % = get_m_comp(get_muon_line_m_extra_lat(x, y, z, total_site, m, tag), i, rho, sigma, lambda)
// % tag = 0 : subtraction
// % tag = 1 : no subtraction
// % m is the muon mass in the lattice unit (or lattice spacing in muon mass unit)
// % total_site is the lattice size in lattice unit
// % Interpolations should be loaded for this function:
// % 0: 6^5 with-sub
// % 1: 8^5 with-sub (used for tag 0)
// % 2: 10^5 with-sub
// % 3: 12^5 with-sub (used for tag 0)
// % 4: 14^5 with-sub
// % 5: 16^5 with-sub (used for tag 0)
// % 6: 6^5 no-sub
// % 7: 8^5 no-sub (used for tag 1)
// % 8: 10^5 no-sub
// % 9: 12^5 no-sub (used for tag 1)
// % l0: 14^5 no-sub
// % 11: 16^5 no-sub (used for tag 1)
{
  // TIMER("get_muon_line_m_extra");
  const std::vector<std::vector<RealD>>& weights =
      get_muon_line_m_extra_weights();
  Qassert(0 <= tag and tag < (int)weights.size());
  const std::vector<RealD> ws = weights[tag];
  const Int size = ws.size();
  ManyMagneticMoments m;
  set_zero(m);
  for (Int i = 0; i < size; ++i) {
    const RealD w = ws[i];
    if (w == 0.0) {
      continue;
    }
    m += w * get_muon_line_m(x, y, z, i);
  }
  return m;
}

inline std::vector<RealD> get_muon_line_m_extra_py(const CoordinateD& x,
                                                    const CoordinateD& y,
                                                    const CoordinateD& z,
                                                    const Int tag)
{
  const ManyMagneticMoments mmm = get_muon_line_m_extra(x, y, z, tag);
  return std_vector_from_many_magnetic_moments(mmm);
}

inline ManyMagneticMoments get_muon_line_m_extra_lat(
    const Coordinate& x, const Coordinate& y, const Coordinate& z,
    const Coordinate& total_site, const RealD a, const Int tag)
// interface
// % \mathcal M_{i,\rho,\sigma,\lambda}(x,y,z)
// % interpolate saved data and extrapolate different interpolations
// % = get_muon_line_m_extra(m * x, m * y, m * z , tag)[16 * sigma  + 4 * lambda + rho][i]
// % = get_m_comp(get_muon_line_m_extra_lat(x, y, z, total_site, m, tag), i, rho, sigma, lambda)
// % tag = 0 : subtraction
// % tag = 1 : no subtraction
// % m is the muon mass in the lattice unit (or lattice spacing in muon mass unit)
// % total_site is the lattice size in lattice unit
{
  // TIMER("get_muon_line_m_extra_lat");
  const Coordinate rxy = relative_coordinate(y - x, total_site);
  const Coordinate ryz = relative_coordinate(z - y, total_site);
  const Coordinate rzx = relative_coordinate(x - z, total_site);
  const Long d2xy = distance_sq_relative_coordinate_g(rxy);
  const Long d2yz = distance_sq_relative_coordinate_g(ryz);
  const Long d2zx = distance_sq_relative_coordinate_g(rzx);
  ManyMagneticMoments mmm;
  set_zero(mmm);
  if (d2xy * sqr(a) > sqr(DISTANCE_LIMIT)) {
    return mmm;
  }
  if (d2yz * sqr(a) > sqr(DISTANCE_LIMIT)) {
    return mmm;
  }
  if (d2zx * sqr(a) > sqr(DISTANCE_LIMIT)) {
    return mmm;
  }
  Int ppws_counts = 0;
  std::vector<PointPairWeight> ppws_total;
  if (d2xy <= d2yz and d2xy <= d2zx) {
    std::vector<PointPairWeight> ppws = shift_lat_corr(x, y, z, total_site, a);
    vector_append(ppws_total, ppws);
    ppws_counts += 1;
  }
  if (d2yz <= d2xy and d2yz <= d2zx) {
    std::vector<PointPairWeight> ppws = shift_lat_corr(y, z, x, total_site, a);
    // y z : z-y -y
    // -z y-z : y z
    for (Int i = 0; i < (int)ppws.size(); ++i) {
      const PointPairWeight ppw = ppws[i];
      ppws[i].rxy = -ppw.rxz;
      ppws[i].rxz = ppw.rxy - ppw.rxz;
    }
    vector_append(ppws_total, ppws);
    ppws_counts += 1;
  }
  if (d2zx <= d2xy and d2zx <= d2yz) {
    std::vector<PointPairWeight> ppws = shift_lat_corr(z, x, y, total_site, a);
    // y z : -z y-z
    // z-y -y : y z
    for (Int i = 0; i < (int)ppws.size(); ++i) {
      const PointPairWeight ppw = ppws[i];
      ppws[i].rxy = ppw.rxz - ppw.rxy;
      ppws[i].rxz = -ppw.rxy;
    }
    vector_append(ppws_total, ppws);
    ppws_counts += 1;
  }
  for (Int i = 0; i < (int)ppws_total.size(); ++i) {
    const PointPairWeight& ppw = ppws_total[i];
    mmm += (ppw.weight / ppws_counts) *
           get_muon_line_m_extra(CoordinateD(), ppw.rxy, ppw.rxz, tag);
  }
  return mmm;
}

inline std::vector<RealD> get_muon_line_m_extra_lat_py(
    const Coordinate& x, const Coordinate& y, const Coordinate& z,
    const Coordinate& total_site, const RealD a, const Int tag)
// interface
// tag = 0 sub
// tag = 1 nosub
{
  const ManyMagneticMoments mmm =
      get_muon_line_m_extra_lat(x, y, z, total_site, a, tag);
  return std_vector_from_many_magnetic_moments(mmm);
}

// ------------------------------------------------------

inline void load_compute_save_muonline_interpolation(
    const std::string& path, const std::vector<Int>& dims,
    const IntegrationEps& eps)
{
  if (!load_multiple_muonline_interpolations(path)) {
    // ADJUST ME
    // load_or_compute_muonline_interpolation(path, dims, eps);
    compute_save_muonline_interpolation_cc(path, dims, eps);
    loadMuonLineInterpolation(path);
  }
}

inline void test_muonline_transformation()
{
  TIMER_VERBOSE_FLOPS("test_muonline_transformation");
  RngState rs(get_global_rng_state(), "test_muonline_transformation");
  const size_t size = 128 * 64 * 2;
  const RealD high = 3.0;
  const RealD low = -3.0;
  const size_t jobs_total = size;
  const size_t num_nodes = get_num_node();
  const size_t jobs_per_nodes = jobs_total / num_nodes;
  const size_t jobs_parallel = jobs_per_nodes * num_nodes;
  const size_t jobs_left = jobs_total - jobs_parallel;
  size_t my_start = jobs_left + jobs_per_nodes * get_id_node();
  const size_t my_end = my_start + jobs_per_nodes;
  if (0 == get_id_node()) {
    my_start = 0;
  }
#pragma omp parallel for schedule(dynamic)
  for (Int i = my_start; i < (int)my_end; ++i) {
    const Int a = i / 9;
    const Int b = i % 9;
    const Int bx = b / 3;
    const Int by = b % 3;
    RngState rsi(rs, a);
    CoordinateD x, y;
    for (Int m = 0; m < 4; ++m) {
      x[m] = u_rand_gen(rsi, high, low);
      y[m] = u_rand_gen(rsi, high, low);
    }
    for (Int m = 0; m < 4; ++m) {
      if (1 == bx) {
        x[m] /= 5.0;
      } else if (2 == bx) {
        x[m] /= 20.0;
      }
      if (1 == by) {
        y[m] /= 5.0;
      } else if (2 == by) {
        y[m] /= 20.0;
      }
    }
    muonLineSymParamsCheck(x, y, IntegrationEps());
  }
  timer.flops += size;
}

inline void test_muonline_transform_scaling()
{
  TIMER_VERBOSE("test_muonline_transform_scaling");
  RngState rs(get_global_rng_state(), "test_muonline_transform_scaling");
  const Int size = 1024;
  const RealD high = 1.0;
  const RealD low = -1.0;
  // const RealD ratio = 2.0;
  // const RealD ratio2 = 10.0;
  const size_t jobs_total = size;
  const size_t num_nodes = get_num_node();
  const size_t jobs_per_nodes = jobs_total / num_nodes;
  const size_t jobs_parallel = jobs_per_nodes * num_nodes;
  const size_t jobs_left = jobs_total - jobs_parallel;
  size_t my_start = jobs_left + jobs_per_nodes * get_id_node();
  const size_t my_end = my_start + jobs_per_nodes;
  if (0 == get_id_node()) {
    my_start = 0;
  }
#pragma omp parallel for schedule(dynamic)
  for (Int i = my_start; i < (int)my_end; ++i) {
    RngState rsi(rs, i);
    CoordinateD x, y;
    for (Int m = 0; m < 4; ++m) {
      x[m] = u_rand_gen(rsi, high, low);
      x[m] *= pow(u_rand_gen(rsi, 1.0, 0.0), 10);
      y[m] = u_rand_gen(rsi, high, low);
      y[m] *= pow(u_rand_gen(rsi, 1.0, 0.0), 10);
    }
    ManyMagneticMoments mmm = muonLineSym(x, y, IntegrationEps());
    ManyMagneticMoments mmmp = muonLineSymTransform(x, y, IntegrationEps(), false);
    // ManyMagneticMoments mmmp = ratio * muonLineSymTransform(x/ratio, y/ratio,
    // 1e-12, 1e-3, false); ManyMagneticMoments mmmpp = ratio2 *
    // muonLineSymTransform(x/ratio2, y/ratio2, 1e-12, 1e-3, false);
    // ManyMagneticMoments mmmppp = ratio*ratio2 *
    // muonLineSymTransform(x/ratio2/ratio, y/ratio2/ratio, 1e-12, 1e-3, false);
    compare_many_magnetic_moments("scaling", x, y, mmm, mmmp);
    // compare_many_magnetic_moments("scaling1", x, y, mmmp, mmmpp);
    // compare_many_magnetic_moments("scaling2", x, y, mmmpp, mmmppp);
  }
}

inline void test_muonline_interp()
{
  TIMER_VERBOSE("test_muonline_interp")
  RngState rs(get_global_rng_state(), "test_muonline_interp");
  const Int size = 1024;
  // const RealD high = 1.0;
  // const RealD low = -1.0;
  // const RealD ratio = 2.0;
  // const RealD ratio2 = 10.0;
  const size_t jobs_total = size;
  const size_t num_nodes = get_num_node();
  const size_t jobs_per_nodes = jobs_total / num_nodes;
  const size_t jobs_parallel = jobs_per_nodes * num_nodes;
  const size_t jobs_left = jobs_total - jobs_parallel;
  size_t my_start = jobs_left + jobs_per_nodes * get_id_node();
  const size_t my_end = my_start + jobs_per_nodes;
  if (0 == get_id_node()) {
    my_start = 0;
  }
  const MuonLineInterp& interp = getMuonLineInterp();
  const size_t interp_size = interp.size();
#pragma omp parallel for schedule(dynamic)
  for (Int i = my_start; i < (int)my_end; ++i) {
    RngState rsi(rs, i);
    const size_t idx = rand_gen(rsi) % interp_size;
    const std::vector<RealD> params = interp.get_coor(idx);
    CoordinateD x, y;
    coordinatesFromParams(x, y, params);
    for (Int i = 0; i < 4; ++i) {
      if (is_very_close(x[i], 0)) {
        x[i] = 0.0;
      }
      if (is_very_close(y[i], 0)) {
        y[i] = 0.0;
      }
    }
    ManyMagneticMoments mmm = muonLineSym(x, y, IntegrationEps());
    ManyMagneticMoments mmmp = muonLineSymTransform(x, y, IntegrationEps(), false);
    ManyMagneticMoments mmmpp = muonLineSymTransform(x, y, IntegrationEps(), true);
    {
      compare_many_magnetic_moments("checking", x, y, mmm, mmmp);
      compare_many_magnetic_moments("checking2", x, y, mmm, mmmpp);
      displayln_c_stdout(ssprintf("test_muonline_interp idx=%d params=%s", idx,
                                  show(params).c_str()));
    }
  }
}

inline void test_muonline_rotate()
{
  TIMER_VERBOSE_FLOPS("test_muonline_rotate");
  RngState rs(get_global_rng_state(), "test_muonline_rotate");
  const size_t size = 128 * 64 * 2;
  const RealD high = 3.0;
  const RealD low = -3.0;
  const size_t jobs_total = size;
  const size_t num_nodes = get_num_node();
  const size_t jobs_per_nodes = jobs_total / num_nodes;
  const size_t jobs_parallel = jobs_per_nodes * num_nodes;
  const size_t jobs_left = jobs_total - jobs_parallel;
  size_t my_start = jobs_left + jobs_per_nodes * get_id_node();
  const size_t my_end = my_start + jobs_per_nodes;
  if (0 == get_id_node()) {
    my_start = 0;
  }
#pragma omp parallel for schedule(dynamic)
  for (Int i = my_start; i < (int)my_end; ++i) {
    RngState rsi(rs, i);
    CoordinateD x, y;
    for (Int m = 0; m < 4; ++m) {
      x[m] = u_rand_gen(rsi, high, low);
      x[m] *= pow(u_rand_gen(rsi, 1.0, 0.0), 5);
      y[m] = u_rand_gen(rsi, high, low);
      y[m] *= pow(u_rand_gen(rsi, 1.0, 0.0), 5);
    }
    muonLineSymRotateCheck(makeRandomSpatialO3Matrix(rsi), x, y, IntegrationEps());
  }
  timer.flops += size;
}

inline void test_muonline_int()
{
  TIMER_VERBOSE_FLOPS("test_muonline_int");
  RngState rs(get_global_rng_state(), "test_muonline_int");
  const size_t size = 128 * 64 * 2;
  const RealD high = 4.0;
  const RealD low = -3.0;
  const size_t jobs_total = size;
  const size_t num_nodes = get_num_node();
  const size_t jobs_per_nodes = jobs_total / num_nodes;
  const size_t jobs_parallel = jobs_per_nodes * num_nodes;
  const size_t jobs_left = jobs_total - jobs_parallel;
  size_t my_start = jobs_left + jobs_per_nodes * get_id_node();
  const size_t my_end = my_start + jobs_per_nodes;
  if (0 == get_id_node()) {
    my_start = 0;
  }
#pragma omp parallel for schedule(dynamic)
  for (Int i = my_start; i < (int)my_end; ++i) {
    RngState rsi(rs, i);
    CoordinateD x, y;
    for (Int m = 0; m < 4; ++m) {
      x[m] = 0.1 * (int)u_rand_gen(rsi, high, low);
      y[m] = 0.1 * (int)u_rand_gen(rsi, high, low);
    }
    muonLineSymParamsCheck(x, y, IntegrationEps());
    // muonLineSymParamsCheck(x, y, 1e-10, 1e-4);
    // muonLineSymRotateCheck(makeRandomSpatialO3Matrix(rsi), x, y);
    // muonLineSymRotateCheck(makeRandomSpatialO3Matrix(rsi), x, y);
    // muonLineSymRotateCheck(makeRandomSpatialO3Matrix(rsi), x, y);
    // muonLineSymRotateCheck(makeRandomSpatialO3Matrix(rsi), x, y);
  }
  timer.flops += size;
}

inline void test_muonLine()
{
  SYNC_NODE();
  TIMER_VERBOSE("test_muonLine");
  if (IS_USING_MUON_LINE_INTERPOLATION) {
    load_compute_save_muonline_interpolation(
        "huge-data-muon-line-interpolation", std::vector<Int>(5, 12),
        IntegrationEps());
  }
  // const CoordinateD x(0.1, 0.2, 0.0, 0.5);
  // const CoordinateD y(0.3, 0.0, -0.2, 0.1);
  // ManyMagneticMoments mmm;
  // mmm = muonLineSym(x,y, 1e-8, 1e-3);
  // displayln(showManyMagneticMoments(mmm));
  // mmm = muonLineSym(x,y);
  // displayln(showManyMagneticMoments(mmm));
  // test_muonline_transform_scaling();
  test_muonline_interp();
  // test_muonline_int();
  // test_muonline_rotate();
  // test_muonline_transformation();
}

}  // namespace qlat
