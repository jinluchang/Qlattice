#pragma once

#include <qlat/field.h>
#include <qlat/field-fft.h>
#include <qlat/qed.h>
#include <qlat/field-expand.h>

#include <eigen3/Eigen/Eigen>

#include <cmath>

QLAT_START_NAMESPACE

typedef Eigen::Matrix<Complex,NUM_COLOR,NUM_COLOR,Eigen::RowMajor> ColorMatrix;

struct GaugeField : FieldM<ColorMatrix,4>
{
  virtual const char* cname()
  {
    return "GaugeField";
  }
};

typedef Eigen::Matrix<Complex,4*NUM_COLOR,4*NUM_COLOR,Eigen::RowMajor> WilsonMatrix;

struct Propagator4d : FieldM<WilsonMatrix,1>
{
  virtual const char* cname()
  {
    return "Propagator4d";
  }
};

inline void unitarize(ColorMatrix& cm)
{
  // fdisplayln(stdout, shows("unitarize before\n") + show(cm));
  cm.row(0).normalize();
  cm.row(1) = cm.row(1) - cm.row(1).dot(cm.row(0)) * cm.row(0);
  cm.row(1).normalize();
  cm.row(2) = cm.row(0).cross(cm.row(1));
  // fdisplayln(stdout, shows("unitarize after\n") + show(cm));
}

inline void unitarize(GaugeField& gf)
{
  TIMER_VERBOSE("unitarize(gf)");
  const Geometry& geo = gf.geo;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    Coordinate xl = geo.coordinate_from_index(index);
    Vector<ColorMatrix> v = gf.get_elems(xl);
    for (int m = 0; m < geo.multiplicity; ++m) {
      unitarize(v[m]);
    }
  }
}

inline double gf_avg_plaq(const GaugeField& gf)
  // assume proper communication is done
{
  TIMER("gf_avg_plaq");
  const Geometry& geo = gf.geo;
  std::vector<double> sums(omp_get_max_threads(), 0.0);
#pragma omp parallel
  {
    double sum_avg_plaq = 0.0;
#pragma omp for
    for (long index = 0; index < geo.local_volume(); ++index) {
      Coordinate xl = geo.coordinate_from_index(index);
      const Vector<ColorMatrix> v = gf.get_elems_const(xl);
      std::vector<Vector<ColorMatrix> > vms(DIM);
      for (int m = 0; m < DIM; ++m) {
        xl[m] += 1;
        vms[m] = gf.get_elems_const(xl);
        xl[m] -= 1;
      }
      double avg_plaq = 0.0;
      for (int m1 = 1; m1 < DIM; ++m1) {
        for (int m2 = 0; m2 < m1; ++m2) {
          ColorMatrix cm = v[m1] * vms[m1][m2] * (v[m2] * vms[m2][m1]).adjoint();
          avg_plaq += cm.trace().real() / NUM_COLOR;
          if (isnan(avg_plaq)) {
            fdisplayln(stdout, ssprintf("WARNING: isnan in gf_avg_plaq"));
            assert(false);
          }
        }
      }
      avg_plaq /= DIM * (DIM-1) / 2;
      sum_avg_plaq += avg_plaq;
    }
    sums[omp_get_thread_num()] = sum_avg_plaq;
  }
  double sum = 0.0;
  for (size_t i = 0; i < sums.size(); ++i) {
    sum += sums[i];
  }
  glb_sum(sum);
  sum /= geo.total_volume();
  return sum;
}

inline double gf_avg_plaq_with_comm(const GaugeField& gf)
{
  TIMER("gf_avg_plaq_with_comm");
  GaugeField gf1;
  gf1.init(geo_resize(gf.geo, 1));
  gf1 = gf;
  refresh_expanded(gf1);
  return gf_avg_plaq(gf1);
}

inline double gf_avg_link_trace(const GaugeField& gf)
{
  TIMER("gf_avg_link_trace");
  const Geometry& geo = gf.geo;
  std::vector<double> sums(omp_get_max_threads(), 0.0);
#pragma omp parallel
  {
    double sum_avg_link_trace = 0.0;
#pragma omp for
    for (long index = 0; index < geo.local_volume(); ++index) {
      Coordinate xl = geo.coordinate_from_index(index);
      const Vector<ColorMatrix> v = gf.get_elems_const(xl);
      double avg_link_trace = 0.0;
      for (int m = 0; m < DIM; ++m) {
        avg_link_trace += v[m].trace().real() / NUM_COLOR;
      }
      avg_link_trace /= DIM;
      sum_avg_link_trace += avg_link_trace;
    }
    sums[omp_get_thread_num()] = sum_avg_link_trace;
  }
  double sum = 0.0;
  for (size_t i = 0; i < sums.size(); ++i) {
    sum += sums[i];
  }
  glb_sum(sum);
  sum /= geo.total_volume();
  return sum;
}

inline void load_gauge_field(GaugeField& gf, const std::string& path)
  // assuming gf already initialized and have correct size;
{
  TIMER_VERBOSE("load_gauge_field");
  assert(is_initialized(gf));
  const Geometry& geo = gf.geo;
  FieldM<std::array<Complex, 6>, 4> gft;
  gft.init(geo);
  sophisticated_serial_read(gft, path, 1);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    Coordinate xl = geo.coordinate_from_index(index);
    Vector<std::array<Complex, 6> > vt = gft.get_elems(xl);
    from_big_endian_64((char*)vt.data(), vt.data_size());
    Vector<ColorMatrix> v = gf.get_elems(xl);
    for (int m = 0; m < geo.multiplicity; ++m) {
      assign_truncate(v[m], vt[m]);
      unitarize(v[m]);
    }
  }
}

QLAT_END_NAMESPACE
