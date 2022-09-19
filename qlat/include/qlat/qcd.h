#pragma once

#include <qlat/config.h>
#include <qlat/field-serial-io.h>
#include <qlat/field.h>
#include <qlat/matrix-hmc.h>
#include <qlat/matrix.h>

#ifndef QLAT_FFTW_OFF
#include <qlat/field-fft.h>
#endif

#include <qlat/field-expand.h>
#include <qlat/fields-io.h>

#include <cmath>
#include <sstream>
#include <string>

namespace qlat
{  //

template <class T = ComplexT>
struct GaugeFieldT : FieldM<ColorMatrixT<T>, 4> {
};

template <class T = ComplexT>
struct Propagator4dT : FieldM<WilsonMatrixT<T>, 1> {
};

template <class T = ComplexT>
struct FermionField4dT : FieldM<WilsonVectorT<T>, 1> {
};

template <class T = ComplexT>
struct FermionField5dT : Field<WilsonVectorT<T> > {
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
  const Geometry& geo = gf.geo();
  qacc_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<ColorMatrixT<T> > v = gf.get_elems(xl);
    for (int m = 0; m < geo.multiplicity; ++m) {
      unitarize(v[m]);
    }
  });
}

template <class T>
double gf_avg_plaq_no_comm(const GaugeFieldT<T>& gf)
// assume proper communication is done
{
  TIMER("gf_avg_plaq_no_comm");
  const Geometry& geo = gf.geo();
  std::vector<double> sums(omp_get_max_threads(), 0.0);
#pragma omp parallel
  {
    double sum_avg_plaq = 0.0;
#pragma omp for
    for (long index = 0; index < geo.local_volume(); ++index) {
      Coordinate xl = geo.coordinate_from_index(index);
      const Vector<ColorMatrixT<T> > v = gf.get_elems_const(xl);
      array<Vector<ColorMatrixT<T> >, DIMN> vms;
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
  gf1.init(geo_resize(gf.geo(), Coordinate(0, 0, 0, 0), Coordinate(1, 1, 1, 1)));
  gf1 = gf;
  refresh_expanded(gf1);
  return gf_avg_plaq_no_comm(gf1);
}

template <class T>
double gf_avg_spatial_plaq_no_comm(const GaugeFieldT<T>& gf)
// assume proper communication is done
{
  TIMER("gf_avg_spatial_plaq_no_comm");
  const Geometry& geo = gf.geo();
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
  gf1.init(geo_resize(gf.geo(), Coordinate(0, 0, 0, 0), Coordinate(1, 1, 1, 0)));
  gf1 = gf;
  refresh_expanded(gf1);
  return gf_avg_spatial_plaq_no_comm(gf1);
}

template <class T>
double gf_avg_link_trace(const GaugeFieldT<T>& gf)
{
  TIMER("gf_avg_link_trace");
  const Geometry& geo = gf.geo();
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
  std::string datatype;
  std::string floating_point;
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
    datatype = "4D_SU3_GAUGE";
    floating_point = "IEEE64BIG";
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
  // const std::string todo = "NOT yet implemented";
  out << "BEGIN_HEADER" << std::endl;
  out << "HDR_VERSION = 1.0" << std::endl;
  out << "DATATYPE = 4D_SU3_GAUGE" << std::endl;
  out << "DIMENSION_1 = " << gfi.total_site[0] << std::endl;
  out << "DIMENSION_2 = " << gfi.total_site[1] << std::endl;
  out << "DIMENSION_3 = " << gfi.total_site[2] << std::endl;
  out << "DIMENSION_4 = " << gfi.total_site[3] << std::endl;
  out << ssprintf("LINK_TRACE = %.12f", gfi.trace) << std::endl;
  out << ssprintf("PLAQUETTE = %.12f", gfi.plaq) << std::endl;
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

inline void read_gauge_field_header(GaugeFieldInfo& gfi,
                                    const std::string& path)
{
  TIMER("read_gauge_field_header");
  if (get_id_node() == 0) {
    FILE* fp = qopen(path, "r");
    if (fp != NULL) {
      const std::string header = "BEGIN_HEADER\n";
      std::vector<char> check_line(header.size(), 0);
      if (1 == fread(check_line.data(), header.size(), 1, fp)) {
        if (std::string(check_line.data(), check_line.size()) == header) {
          std::vector<std::string> infos;
          infos.push_back(header);
          while (infos.back() != "END_HEADER\n" && infos.back() != "") {
            infos.push_back(qgetline(fp));
          }
          for (int m = 0; m < 4; ++m) {
            reads(gfi.total_site[m],
                  info_get_prop(infos, ssprintf("DIMENSION_%d = ", m + 1)));
          }
          reads(gfi.trace, info_get_prop(infos, "LINK_TRACE = "));
          reads(gfi.plaq, info_get_prop(infos, "PLAQUETTE = "));
          std::string info;
          info = info_get_prop(infos, "CHECKSUM = ");
          if (info != "") {
            gfi.simple_checksum = read_crc32(info);
          }
          info = info_get_prop(infos, "CRC32HASH = ");
          if (info != "") {
            gfi.crc32 = read_crc32(info);
          }
          gfi.datatype = info_get_prop(infos, "DATATYPE = ");
          gfi.floating_point = info_get_prop(infos, "FLOATING_POINT = ");
          remove_trailing_newline(gfi.datatype);
          remove_trailing_newline(gfi.floating_point);
        }
      }
    }
    qclose(fp);
  }
  bcast(Vector<Coordinate>(&gfi.total_site, 1));
  bcast(Vector<double>(&gfi.trace, 1));
  bcast(Vector<double>(&gfi.plaq, 1));
  bcast(Vector<crc32_t>(&gfi.simple_checksum, 1));
  bcast(Vector<crc32_t>(&gfi.crc32, 1));
  bcast(gfi.floating_point);
  bcast(gfi.datatype);
}

template <class T>
long save_gauge_field(const GaugeFieldT<T>& gf, const std::string& path,
                      const GaugeFieldInfo& gfi_ = GaugeFieldInfo())
{
  TIMER_VERBOSE_FLOPS("save_gauge_field");
  qassert(is_initialized(gf));
  const Geometry& geo = gf.geo();
  FieldM<array<Complex, 6>, 4> gft;
  gft.init(geo);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<ColorMatrixT<T> > v = gf.get_elems_const(xl);
    Vector<array<Complex, 6> > vt = gft.get_elems(xl);
    for (int m = 0; m < geo.multiplicity; ++m) {
      vt[m][0] = v[m](0, 0);
      vt[m][1] = v[m](0, 1);
      vt[m][2] = v[m](0, 2);
      vt[m][3] = v[m](1, 0);
      vt[m][4] = v[m](1, 1);
      vt[m][5] = v[m](1, 2);
    }
  }
  GaugeFieldInfo gfi = gfi_;
  gfi.simple_checksum = field_simple_checksum(gft);
  to_from_big_endian_64(get_data(gft));
  gfi.plaq = gf_avg_plaq(gf);
  gfi.trace = gf_avg_link_trace(gf);
  gfi.crc32 = field_crc32(gft);
  gfi.total_site = gf.geo().total_site();
  qtouch_info(path + ".partial", make_gauge_field_header(gfi));
  const long file_size = serial_write_field(gft, path + ".partial");
  qrename_info(path + ".partial", path);
  timer.flops += file_size;
  return file_size;
}

template <class T = ComplexT>
long load_gauge_field(GaugeFieldT<T>& gf, const std::string& path)
// assuming gf already initialized and have correct size;
{
  TIMER_VERBOSE_FLOPS("load_gauge_field");
  displayln_info(fname + ssprintf(": '%s'.", path.c_str()));
  gf.init();
  GaugeFieldInfo gfi;
  read_gauge_field_header(gfi, path);
  const bool is_two_row = gfi.datatype == "4D_SU3_GAUGE";
  const bool is_three_row = gfi.datatype == "4D_SU3_GAUGE_3x3";
  const int n_complex_su3 = is_two_row ? 6 : ( is_three_row ? 9 : 0);
  if (n_complex_su3 == 0) {
    displayln(fname + ssprintf(": gfi.datatype '%s' id_node=%d.",
                               gfi.datatype.c_str(), get_id_node()));
    qassert(false);
  }
  Geometry geo;
  geo.init(gfi.total_site, 4);
  Field<Complex> gft;
  gft.init(geo_remult(geo, 4 * n_complex_su3));
  const long file_size = serial_read_field_par(
      gft, path, -get_data_size(gft) * get_num_node(), SEEK_END);
  if (0 == file_size) {
    return 0;
  }
  if (gfi.floating_point == "IEEE64BIG") {
    to_from_big_endian_64(get_data(gft));
  } else if (gfi.floating_point == "IEEE64LITTLE") {
    to_from_little_endian_64(get_data(gft));
  } else {
    qassert(false);
  }
  crc32_t simple_checksum = field_simple_checksum(gft);
  if (simple_checksum != gfi.simple_checksum) {
    if (get_id_node() == 0) {
      qwarn(fname +
            ssprintf(": WARNING: fn='%s' CHECKSUM= %08X (calc) %08X (read)",
                     path.c_str(), simple_checksum, gfi.simple_checksum));
    }
  }
  gf.init(geo);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<Complex> vt = gft.get_elems(xl);
    Vector<ColorMatrixT<T> > v = gf.get_elems(xl);
    for (int m = 0; m < geo.multiplicity; ++m) {
      v[m](0, 0) = vt[m * n_complex_su3 + 0];
      v[m](0, 1) = vt[m * n_complex_su3 + 1];
      v[m](0, 2) = vt[m * n_complex_su3 + 2];
      v[m](1, 0) = vt[m * n_complex_su3 + 3];
      v[m](1, 1) = vt[m * n_complex_su3 + 4];
      v[m](1, 2) = vt[m * n_complex_su3 + 5];
      if (is_three_row) {
        v[m](2, 0) = vt[m * n_complex_su3 + 6];
        v[m](2, 1) = vt[m * n_complex_su3 + 7];
        v[m](2, 2) = vt[m * n_complex_su3 + 8];
      } else {
        unitarize(v[m]);
      }
    }
  }
  timer.flops += file_size;
  return file_size;
}

inline long load_gauge_field_cps3x3(GaugeFieldT<Complex>& gf,
                                    const std::string& path)
// assuming gf already initialized and have correct size;
{
  TIMER_VERBOSE_FLOPS("load_gauge_field_cps3x3");
  displayln_info(fname + ssprintf(": '%s'.", path.c_str()));
  qassert(is_initialized(gf));
  const Geometry& geo = gf.geo();
  FieldM<array<Complex, 9>, 4> gft;
  gft.init(geo);
  const long file_size = serial_read_field_par(
      gft, path, -get_data_size(gft) * get_num_node(), SEEK_END);
  if (file_size == 0) {
    return 0;
  }
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<array<Complex, 9> > vt = gft.get_elems(xl);
    to_from_big_endian_64(get_data(vt));
    Vector<ColorMatrixT<Complex> > v = gf.get_elems(xl);
    for (int m = 0; m < geo.multiplicity; ++m) {
      assign_truncate(v[m], vt[m]);
    }
  }
  timer.flops += file_size;
  return file_size;
}

inline long load_gauge_field_milc(GaugeFieldT<Complex>& gf,
                                  const std::string& path,
                                  const bool par_read = false)
// assuming gf already initialized and have correct size;
{
  TIMER_VERBOSE_FLOPS("load_gauge_field_milc");
  displayln_info(fname + ssprintf(": '%s'.", path.c_str()));
  qassert(is_initialized(gf));
  const Geometry& geo = gf.geo();
  FieldM<array<ComplexF, 9>, 4> gft;
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
    Vector<array<ComplexF, 9> > vt = gft.get_elems(xl);
    to_from_big_endian_32((char*)vt.data(), vt.data_size());
    Vector<ColorMatrixT<Complex> > v = gf.get_elems(xl);
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
void twist_boundary_at_boundary(GaugeFieldT<T>& gf, double lmom, int mu)
{
  TIMER_VERBOSE_FLOPS("twist_boundary_at_boundary");
  const Geometry& geo = gf.geo();
  const double amp = 2.0 * PI * lmom;
  const int len = geo.total_site()[mu];
  for (int index = 0; index < geo.local_volume(); index++) {
    Coordinate xl = geo.coordinate_from_index(index);
    Coordinate xg = geo.coordinate_g_from_l(xl);
    if (xg[mu] == len - 1) {
      ColorMatrixT<T>& mat = gf.get_elem(xl, mu);
      mat *= ComplexT(std::polar(1.0, amp));
    }
  }
}

}  // namespace qlat
