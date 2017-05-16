#pragma once

#include <qlat/matrix.h>
#include <qlat/field.h>
#include <qlat/field-fft.h>
#include <qlat/field-expand.h>

#include <eigen3/Eigen/Eigen>

#include <cmath>
#include <string>
#include <sstream>

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

struct GaugeFieldInfo
{
  std::string ensemble_id;
  std::string ensemble_label;
  std::string creator;
  std::string date;
  long sequence_num;
  double beta;
  double plaq, trace;
  crc32_t simple_checksum, crc32;
  Coordinate total_site;
  //
  GaugeFieldInfo()
  {
    ensemble_id = "42";
    ensemble_label = "default-ensemble";
    creator = "Qlat";
    time_t now = std::time(NULL);	
    date = shows(std::ctime(&now));
    sequence_num = 0;
    beta = 0.0;
    plaq = 1.0;
    trace = 0.0;
    simple_checksum = 0;
    crc32 = 0;
  }
};

std::string make_gauge_field_header(const GaugeFieldInfo& gfi = GaugeFieldInfo())
{
  std::ostringstream out;
  const std::string todo = "NOT yet implemented";
  out << "BEGIN_HEADER" << std::endl;
  out << "HDR_VERSION = 1.0" << std::endl;
  out << "DATATYPE = 4D_SU3_GAUGE" << std::endl;
  out << "DIMENSION_1 = " << gfi.total_site[0] << std::endl;
  out << "DIMENSION_2 = " << gfi.total_site[1] << std::endl;
  out << "DIMENSION_3 = " << gfi.total_site[2] << std::endl;
  out << "DIMENSION_4 = " << gfi.total_site[3] << std::endl;
  out << ssprintf("LINK_TRACE = %.12f", gfi.plaq) << std::endl;
  out << ssprintf("PLAQUETTE = %.12f", gfi.trace) << std::endl;
  out << ssprintf("CHECKSUM = %08x", gfi.simple_checksum) << std::endl;
  out << ssprintf("CRC32HASH = %08x", gfi.crc32) << std::endl;
  out << "CREATOR = " << gfi.creator << std::endl;
  out << "ARCHIVE_DATE = " << gfi.date << std::endl;
  out << "ENSEMBLE_ID = " << gfi.ensemble_id << std::endl;
  out << "ENSEMBLE_LABEL = " << gfi.ensemble_label << std::endl;
  out << ssprintf("BETA = %.12f", gfi.beta) << std::endl; 
  out << ssprintf("SEQUENCE_NUMBER = %ld", gfi.sequence_num) << std::endl;
  out << "FLOATING_POINT = IEEE64BIG" << std::endl;
  out << "END_HEADER" << std::endl;
  return out.str();
}

inline void save_gauge_field(const GaugeField& gf, const std::string& path, const GaugeFieldInfo& gfi_ = GaugeFieldInfo())
{
  TIMER_VERBOSE("save_gauge_field");
  qassert(is_initialized(gf));
  const Geometry& geo = gf.geo;
  FieldM<std::array<Complex, 6>, 4> gft;
  gft.init(geo);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<ColorMatrix> v = gf.get_elems_const(xl);
    Vector<std::array<Complex, 6> > vt = gft.get_elems(xl);
    for (int m = 0; m < geo.multiplicity; ++m) {
      assign_truncate(vt[m], v[m]);
    }
    to_from_big_endian_64(get_data(vt));
  }
  GaugeFieldInfo gfi = gfi_;
  gfi.plaq = gf_avg_plaq(gf);
  gfi.trace = gf_avg_link_trace(gf);
  gfi.simple_checksum = field_simple_checksum(gft);
  gfi.crc32 = field_crc32(gft);
  gfi.total_site = gf.geo.total_site();
  qtouch_info(path, make_gauge_field_header(gfi));
  serial_write_field(gft, path);
}

inline void load_gauge_field(GaugeField& gf, const std::string& path, const bool par_read = false)
  // assuming gf already initialized and have correct size;
{
  TIMER_VERBOSE("load_gauge_field");
  qassert(is_initialized(gf));
  const Geometry& geo = gf.geo;
  FieldM<std::array<Complex, 6>, 4> gft;
  gft.init(geo);
  if (par_read) {
    serial_read_field_par(gft, path, -get_data_size(gft) * get_num_node(), SEEK_END);
  } else {
    serial_read_field(gft, path, -get_data_size(gft) * get_num_node(), SEEK_END);
  }
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<std::array<Complex, 6> > vt = gft.get_elems(xl);
    to_from_big_endian_64(get_data(vt));
    Vector<ColorMatrix> v = gf.get_elems(xl);
    for (int m = 0; m < geo.multiplicity; ++m) {
      assign_truncate(v[m], vt[m]);
      unitarize(v[m]);
    }
  }
}

inline void load_gauge_field_milc(GaugeField& gf, const std::string& path, const bool par_read = false)
  // assuming gf already initialized and have correct size;
{
  TIMER_VERBOSE("load_gauge_field_milc");
  qassert(is_initialized(gf));
  const Geometry& geo = gf.geo;
  FieldM<std::array<std::complex<float>, 9>, 4> gft;
  gft.init(geo);
  // ADJUST ME
  if (par_read) {
    serial_read_field_par(gft, path, 0x730, SEEK_SET);
  } else {
    serial_read_field(gft, path, 0x730, SEEK_SET);
  }
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    Coordinate xl = geo.coordinate_from_index(index);
    Vector<std::array<std::complex<float>, 9> > vt = gft.get_elems(xl);
    to_from_big_endian_32((char*)vt.data(), vt.data_size());
    Vector<ColorMatrix> v = gf.get_elems(xl);
    for (int m = 0; m < geo.multiplicity; ++m) {
      // assign_truncate(v[m], vt[m]);
      v[m](0,0) = vt[m][0*3 + 0];
      v[m](0,1) = vt[m][0*3 + 1];
      v[m](0,2) = vt[m][0*3 + 2];
      v[m](1,0) = vt[m][1*3 + 0];
      v[m](1,1) = vt[m][1*3 + 1];
      v[m](1,2) = vt[m][1*3 + 2];
      v[m](2,0) = vt[m][2*3 + 0];
      v[m](2,1) = vt[m][2*3 + 1];
      v[m](2,2) = vt[m][2*3 + 2];
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
