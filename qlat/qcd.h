#pragma once

#include <qlat/matrix.h>
#include <qlat/field.h>
#include <qlat/field-fft.h>
#include <qlat/field-expand.h>

#include <eigen3/Eigen/Eigen>

#include <cmath>

QLAT_START_NAMESPACE

struct GaugeField : FieldM<ColorMatrix,4>
{
  virtual const std::string& cname()
  {
    static const std::string s = "GaugeField";
    return s;
  }
};

struct Propagator4d : FieldM<WilsonMatrix,1>
{
  virtual const std::string&  cname()
  {
    static const std::string s = "Propagator4d";
    return s;
  }
};

struct FermionField4d : FieldM<WilsonVector,1>
{
  virtual const std::string& cname()
  {
    static const std::string s = "FermionField4d";
    return s;
  }
};

inline void unitarize(Field<ColorMatrix>& gf)
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

inline double gf_avg_plaq_no_comm(const GaugeField& gf)
  // assume proper communication is done
{
  TIMER("gf_avg_plaq_no_comm");
  const Geometry& geo = gf.geo;
  std::vector<double> sums(omp_get_max_threads(), 0.0);
#pragma omp parallel
  {
    double sum_avg_plaq = 0.0;
#pragma omp for
    for (long index = 0; index < geo.local_volume(); ++index) {
      Coordinate xl = geo.coordinate_from_index(index);
      const Vector<ColorMatrix> v = gf.get_elems_const(xl);
      std::vector<Vector<ColorMatrix> > vms(DIMN);
      for (int m = 0; m < DIMN; ++m) {
        xl[m] += 1;
        vms[m] = gf.get_elems_const(xl);
        xl[m] -= 1;
      }
      double avg_plaq = 0.0;
      for (int m1 = 1; m1 < DIMN; ++m1) {
        for (int m2 = 0; m2 < m1; ++m2) {
          ColorMatrix cm = v[m1] * vms[m1][m2] * matrix_adjoint(v[m2] * vms[m2][m1]);
          avg_plaq += matrix_trace(cm).real() / NUM_COLOR;
          if (std::isnan(avg_plaq)) {
            fdisplayln(stdout, ssprintf("WARNING: isnan in gf_avg_plaq"));
            qassert(false);
          }
        }
      }
      avg_plaq /= DIMN * (DIMN-1) / 2;
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

inline double gf_avg_plaq(const GaugeField& gf)
{
  TIMER("gf_avg_plaq");
  GaugeField gf1;
  gf1.init(geo_resize(gf.geo, Coordinate(0,0,0,0), Coordinate(1,1,1,1)));
  gf1 = gf;
  refresh_expanded(gf1);
  return gf_avg_plaq_no_comm(gf1);
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
      for (int m = 0; m < DIMN; ++m) {
        avg_link_trace += matrix_trace(v[m]).real() / NUM_COLOR;
      }
      avg_link_trace /= DIMN;
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
  qassert(is_initialized(gf));
  const Geometry& geo = gf.geo;
  FieldM<std::array<Complex, 6>, 4> gft;
  gft.init(geo);
  field_import_serial(gft, path, -get_data_size(gft) * get_num_node(), SEEK_END);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    Coordinate xl = geo.coordinate_from_index(index);
    Vector<std::array<Complex, 6> > vt = gft.get_elems(xl);
    to_from_big_endian_64(get_data(vt));
    Vector<ColorMatrix> v = gf.get_elems(xl);
    for (int m = 0; m < geo.multiplicity; ++m) {
      assign_truncate(v[m], vt[m]);
      unitarize(v[m]);
    }
  }
}

inline WilsonMatrix make_wilson_matrix_from_vectors(const std::array<ConstHandle<WilsonVector>,4*NUM_COLOR>& cols)
{
  WilsonMatrix ret;
  for (int i = 0; i < 4*NUM_COLOR; ++i) {
    for (int j = 0; j < 4*NUM_COLOR; ++j) {
      ret(j, i) = cols[i]()(j);
    }
  }
  return ret;
}

inline void set_propagator_from_fermion_fields(Propagator4d& prop, const Array<FermionField4d,4*NUM_COLOR> ffs)
{
  TIMER_VERBOSE("set_propagator_from_fermion_fields");
  const Geometry geo = ffs[0].geo;
  for (int i = 0; i < 4*NUM_COLOR; ++i) {
    qassert(geo == ffs[i].geo);
  }
  prop.init(geo);
  qassert(prop.geo == geo);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    std::array<ConstHandle<WilsonVector>, 4*NUM_COLOR> cols;
    for (int k = 0; k < 4*NUM_COLOR; ++k) {
      cols[k].init(ffs[k].get_elem(xl));
    }
    prop.get_elem(xl) = make_wilson_matrix_from_vectors(cols);
  }
}

QLAT_END_NAMESPACE
