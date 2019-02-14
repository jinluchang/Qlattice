#pragma once

#include <qlat/config.h>
#include <qlat/field.h>
#include <qlat/matrix.h>

#ifndef QLAT_FFTW_OFF
#include <qlat/field-fft.h>
#endif

#include <qlat/field-expand.h>

#include <cmath>
#include <sstream>
#include <string>

QLAT_START_NAMESPACE

template <class T = ComplexT>
struct GaugeFieldT : FieldM<ColorMatrixT<T>, 4> {
  virtual const std::string& cname()
  {
    static const std::string s = "GaugeField";
    return s;
  }
};

template <class T = ComplexT>
struct Propagator4dT : FieldM<WilsonMatrixT<T>, 1> {
  virtual const std::string& cname()
  {
    static const std::string s = "Propagator4d";
    return s;
  }
};

template <class T = ComplexT>
struct FermionField4dT : FieldM<WilsonVectorT<T>, 1> {
  virtual const std::string& cname()
  {
    static const std::string s = "FermionField4d";
    return s;
  }
};

template <class T = ComplexT>
struct FermionField5dT : Field<WilsonVectorT<T> > {
  virtual const std::string& cname()
  {
    static const std::string s = "FermionField5d";
    return s;
  }
};

#ifndef QLAT_NO_DEFAULT_TYPE

typedef GaugeFieldT<> GaugeField;

typedef Propagator4dT<> Propagator4d;

typedef FermionField4dT<> FermionField4d;

typedef FermionField5dT<> FermionField5d;

#endif

template <class T>
void unitarize(Field<ColorMatrixT<T> >& gf)
{
  TIMER_VERBOSE("unitarize(gf)");
  const Geometry& geo = gf.geo;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    Coordinate xl = geo.coordinate_from_index(index);
    Vector<ColorMatrixT<T> > v = gf.get_elems(xl);
    for (int m = 0; m < geo.multiplicity; ++m) {
      unitarize(v[m]);
    }
  }
}

template <class T>
double gf_avg_plaq_no_comm(const GaugeFieldT<T>& gf)
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
      const Vector<ColorMatrixT<T> > v = gf.get_elems_const(xl);
      std::vector<Vector<ColorMatrixT<T> > > vms(DIMN);
      for (int m = 0; m < DIMN; ++m) {
        xl[m] += 1;
        vms[m] = gf.get_elems_const(xl);
        xl[m] -= 1;
      }
      double avg_plaq = 0.0;
      for (int m1 = 1; m1 < DIMN; ++m1) {
        for (int m2 = 0; m2 < m1; ++m2) {
          ColorMatrixT<T> cm =
              v[m1] * vms[m1][m2] * matrix_adjoint(v[m2] * vms[m2][m1]);
          avg_plaq += matrix_trace(cm).real() / NUM_COLOR;
          if (std::isnan(avg_plaq)) {
            fdisplayln(stdout, ssprintf("WARNING: isnan in gf_avg_plaq"));
            qassert(false);
          }
        }
      }
      avg_plaq /= DIMN * (DIMN - 1) / 2;
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

template <class T>
double gf_avg_plaq(const GaugeFieldT<T>& gf)
{
  TIMER("gf_avg_plaq");
  GaugeFieldT<T> gf1;
  gf1.init(geo_resize(gf.geo, Coordinate(0, 0, 0, 0), Coordinate(1, 1, 1, 1)));
  gf1 = gf;
  refresh_expanded(gf1);
  return gf_avg_plaq_no_comm(gf1);
}

template <class T>
double gf_avg_spatial_plaq_no_comm(const GaugeFieldT<T>& gf)
// assume proper communication is done
{
  TIMER("gf_avg_spatial_plaq_no_comm");
  const Geometry& geo = gf.geo;
  std::vector<double> sums(omp_get_max_threads(), 0.0);
#pragma omp parallel
  {
    double sum_avg_plaq = 0.0;
#pragma omp for
    for (long index = 0; index < geo.local_volume(); ++index) {
      Coordinate xl = geo.coordinate_from_index(index);
      const Vector<ColorMatrixT<T> > v = gf.get_elems_const(xl);
      std::vector<Vector<ColorMatrixT<T> > > vms(DIMN - 1);
      for (int m = 0; m < DIMN - 1; ++m) {
        xl[m] += 1;
        vms[m] = gf.get_elems_const(xl);
        xl[m] -= 1;
      }
      double avg_plaq = 0.0;
      for (int m1 = 1; m1 < 3; ++m1) {
        for (int m2 = 0; m2 < m1; ++m2) {
          ColorMatrixT<T> cm =
              v[m1] * vms[m1][m2] * matrix_adjoint(v[m2] * vms[m2][m1]);
          avg_plaq += matrix_trace(cm).real() / NUM_COLOR;
          if (std::isnan(avg_plaq)) {
            fdisplayln(stdout, ssprintf("WARNING: isnan in gf_avg_plaq"));
            qassert(false);
          }
        }
      }
      avg_plaq /= (DIMN - 1) * (DIMN - 2) / 2;
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

template <class T>
double gf_avg_spatial_plaq(const GaugeFieldT<T>& gf)
{
  TIMER("gf_avg_spatial_plaq");
  GaugeFieldT<T> gf1;
  gf1.init(geo_resize(gf.geo, Coordinate(0, 0, 0, 0), Coordinate(1, 1, 1, 0)));
  gf1 = gf;
  refresh_expanded(gf1);
  return gf_avg_spatial_plaq_no_comm(gf1);
}

template <class T>
double gf_avg_link_trace(const GaugeFieldT<T>& gf)
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
      const Vector<ColorMatrixT<T> > v = gf.get_elems_const(xl);
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

struct GaugeFieldInfo {
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

inline std::string make_gauge_field_header(
    const GaugeFieldInfo& gfi = GaugeFieldInfo())
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
  out << "ARCHIVE_DATE = " << gfi.date;
  out << "ENSEMBLE_ID = " << gfi.ensemble_id << std::endl;
  out << "ENSEMBLE_LABEL = " << gfi.ensemble_label << std::endl;
  out << ssprintf("BETA = %.12f", gfi.beta) << std::endl;
  out << ssprintf("SEQUENCE_NUMBER = %ld", gfi.sequence_num) << std::endl;
  out << "FLOATING_POINT = IEEE64BIG" << std::endl;
  out << "END_HEADER" << std::endl;
  return out.str();
}

inline void save_gauge_field(const GaugeFieldT<Complex>& gf, const std::string& path,
                             const GaugeFieldInfo& gfi_ = GaugeFieldInfo())
{
  TIMER_VERBOSE_FLOPS("save_gauge_field");
  qassert(is_initialized(gf));
  const Geometry& geo = gf.geo;
  FieldM<std::array<Complex, 6>, 4> gft;
  gft.init(geo);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<ColorMatrixT<> > v = gf.get_elems_const(xl);
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
  timer.flops += get_data(gft).data_size() * gft.geo.geon.num_node;
}

inline long load_gauge_field(GaugeFieldT<Complex>& gf, const std::string& path)
// assuming gf already initialized and have correct size;
{
  TIMER_VERBOSE_FLOPS("load_gauge_field");
  displayln_info(fname + ssprintf(": '%s'.", path.c_str()));
  qassert(is_initialized(gf));
  const Geometry& geo = gf.geo;
  FieldM<std::array<Complex, 6>, 4> gft;
  gft.init(geo);
  const long file_size = serial_read_field(
      gft, path, -get_data_size(gft) * get_num_node(), SEEK_END);
  if (0 == file_size) {
    return 0;
  }
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<std::array<Complex, 6> > vt = gft.get_elems(xl);
    to_from_big_endian_64(get_data(vt));
    Vector<ColorMatrixT<> > v = gf.get_elems(xl);
    for (int m = 0; m < geo.multiplicity; ++m) {
      assign_truncate(v[m], vt[m]);
      unitarize(v[m]);
    }
  }
  timer.flops += file_size;
  return file_size;
}

inline long load_gauge_field_par(GaugeFieldT<Complex>& gf, const std::string& path)
// assuming gf already initialized and have correct size;
{
  TIMER_VERBOSE_FLOPS("load_gauge_field_par");
  displayln_info(fname + ssprintf(": '%s'.", path.c_str()));
  qassert(is_initialized(gf));
  const Geometry& geo = gf.geo;
  FieldM<std::array<Complex, 6>, 4> gft;
  gft.init(geo);
  const long file_size = serial_read_field_par(
      gft, path, -get_data_size(gft) * get_num_node(), SEEK_END);
  if (file_size == 0) {
    return 0;
  }
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<std::array<Complex, 6> > vt = gft.get_elems(xl);
    to_from_big_endian_64(get_data(vt));
    Vector<ColorMatrixT<> > v = gf.get_elems(xl);
    for (int m = 0; m < geo.multiplicity; ++m) {
      assign_truncate(v[m], vt[m]);
      unitarize(v[m]);
    }
  }
  timer.flops += file_size;
  return file_size;
}

inline long load_gauge_field_cps3x3(GaugeFieldT<Complex>& gf, const std::string& path)
// assuming gf already initialized and have correct size;
{
  TIMER_VERBOSE_FLOPS("load_gauge_field_cps3x3");
  displayln_info(fname + ssprintf(": '%s'.", path.c_str()));
  qassert(is_initialized(gf));
  const Geometry& geo = gf.geo;
  FieldM<std::array<Complex, 9>, 4> gft;
  gft.init(geo);
  const long file_size = serial_read_field(
      gft, path, -get_data_size(gft) * get_num_node(), SEEK_END);
  if (file_size == 0) {
    return 0;
  }
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<std::array<Complex, 9> > vt = gft.get_elems(xl);
    to_from_big_endian_64(get_data(vt));
    Vector<ColorMatrixT<> > v = gf.get_elems(xl);
    for (int m = 0; m < geo.multiplicity; ++m) {
      assign_truncate(v[m], vt[m]);
    }
  }
  timer.flops += file_size;
  return file_size;
}

inline long load_gauge_field_milc(GaugeFieldT<Complex>& gf, const std::string& path,
                                  const bool par_read = false)
// assuming gf already initialized and have correct size;
{
  TIMER_VERBOSE_FLOPS("load_gauge_field_milc");
  displayln_info(fname + ssprintf(": '%s'.", path.c_str()));
  qassert(is_initialized(gf));
  const Geometry& geo = gf.geo;
  FieldM<std::array<std::complex<float>, 9>, 4> gft;
  gft.init(geo);
  // ADJUST ME
  long file_size = 0;
  if (par_read) {
    file_size = serial_read_field_par(gft, path, 0x730, SEEK_SET);
  } else {
    file_size = serial_read_field(gft, path, 0x730, SEEK_SET);
  }
  if (0 == file_size) {
    return 0;
  }
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    Coordinate xl = geo.coordinate_from_index(index);
    Vector<std::array<std::complex<float>, 9> > vt = gft.get_elems(xl);
    to_from_big_endian_32((char*)vt.data(), vt.data_size());
    Vector<ColorMatrixT<> > v = gf.get_elems(xl);
    for (int m = 0; m < geo.multiplicity; ++m) {
      // assign_truncate(v[m], vt[m]);
      v[m](0, 0) = vt[m][0 * 3 + 0];
      v[m](0, 1) = vt[m][0 * 3 + 1];
      v[m](0, 2) = vt[m][0 * 3 + 2];
      v[m](1, 0) = vt[m][1 * 3 + 0];
      v[m](1, 1) = vt[m][1 * 3 + 1];
      v[m](1, 2) = vt[m][1 * 3 + 2];
      v[m](2, 0) = vt[m][2 * 3 + 0];
      v[m](2, 1) = vt[m][2 * 3 + 1];
      v[m](2, 2) = vt[m][2 * 3 + 2];
      unitarize(v[m]);
    }
  }
  timer.flops += file_size;
  return file_size;
}

template <class T>
void twist_boundary_at_boundary(GaugeFieldT<T>& gf, double mom, int mu)
{
  TIMER_VERBOSE_FLOPS("twist_boundary_at_boundary");
  const Geometry& geo = gf.geo;
  const double amp = 2.0 * PI * mom;
  const int len = geo.total_site()[mu];
  for (int index = 0; index < geo.local_volume(); index++) {
    Coordinate xl = geo.coordinate_from_index(index);
    Coordinate xg = geo.coordinate_g_from_l(xl);
    if (xg[mu] == len - 1) {
      ColorMatrixT<T>& mat = gf.get_elem(xl, mu);
      mat *= std::polar(1.0, amp);
    }
  }
}

QLAT_END_NAMESPACE
