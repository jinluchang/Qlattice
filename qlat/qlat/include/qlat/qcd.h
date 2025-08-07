#pragma once

#include <qlat-utils/matrix-hmc.h>
#include <qlat-utils/matrix.h>
#include <qlat/field-expand.h>
#include <qlat/field-io.h>
#include <qlat/fields-io.h>

#include <cmath>
#include <sstream>
#include <string>

namespace qlat
{  //

RealD gf_avg_spatial_plaq(const GaugeField& gf);

RealD gf_avg_plaq(const GaugeField& gf);

RealD gf_avg_link_trace(const GaugeField& gf);

void gf_plaq_field(Field<RealD>& f_plaq, const GaugeField& gf);

struct U1GaugeTransform : FieldM<ComplexF, 1> {
};

template <class T>
void unitarize(Field<ColorMatrixT<T> >& gf)
{
  TIMER_VERBOSE("unitarize(gf)");
  qacc_for(index, gf.geo().local_volume(), {
    const Geometry& geo = gf.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<ColorMatrixT<T>> v = gf.get_elems(xl);
    for (int m = 0; m < gf.multiplicity; ++m) {
      unitarize(v[m]);
    }
  });
}

// GaugeField IO

struct API GaugeFieldInfo {
  std::string ensemble_id;
  std::string ensemble_label;
  std::string creator;
  std::string date;
  std::string datatype;
  std::string floating_point;
  Long sequence_num;
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
    QFile qfile = qfopen(path, "r");
    if (not qfile.null()) {
      const std::string header = "BEGIN_HEADER\n";
      std::vector<char> check_line(header.size(), 0);
      if (1 == qfread(check_line.data(), header.size(), 1, qfile)) {
        if (std::string(check_line.data(), check_line.size()) == header) {
          std::vector<std::string> infos;
          infos.push_back(header);
          while (infos.back() != "END_HEADER\n" && infos.back() != "") {
            infos.push_back(qgetline(qfile));
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
          gfi.datatype = remove_trailing_newline(gfi.datatype);
          gfi.floating_point = info_get_prop(infos, "FLOATING_POINT = ");
          gfi.floating_point = remove_trailing_newline(gfi.floating_point);
        }
      }
    }
    qfclose(qfile);
  }
  bcast(gfi.total_site);
  bcast(gfi.trace);
  bcast(gfi.plaq);
  bcast(gfi.simple_checksum);
  bcast(gfi.crc32);
  bcast(gfi.floating_point);
  bcast(gfi.datatype);
}

template <class T>
Long save_gauge_field(const GaugeFieldT<T>& gf, const std::string& path,
                      const GaugeFieldInfo& gfi_ = GaugeFieldInfo())
{
  TIMER_VERBOSE_FLOPS("save_gauge_field");
  qassert(is_initialized(gf));
  const Geometry& geo = gf.geo();
  const Int multiplicity = gf.multiplicity;
  FieldM<ComplexD, 4 * 6> gft;
  gft.init(geo);
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<ColorMatrixT<T> > v = gf.get_elems_const(xl);
    Vector<ComplexD> vt = gft.get_elems(xl);
    for (int m = 0; m < multiplicity; ++m) {
      vt[6 * m + 0] = v[m](0, 0);
      vt[6 * m + 1] = v[m](0, 1);
      vt[6 * m + 2] = v[m](0, 2);
      vt[6 * m + 3] = v[m](1, 0);
      vt[6 * m + 4] = v[m](1, 1);
      vt[6 * m + 5] = v[m](1, 2);
    }
  }
  GaugeFieldInfo gfi = gfi_;
  gfi.simple_checksum = field_simple_checksum(gft); // before to_from_big_endian_64
  to_from_big_endian(get_data(gft));
  gfi.plaq = gf_avg_plaq(gf);
  gfi.trace = gf_avg_link_trace(gf);
  gfi.crc32 = field_crc32(gft);
  gfi.total_site = gf.geo().total_site();
  qtouch_info(path + ".partial", make_gauge_field_header(gfi));
  const Long file_size = serial_write_field(gft, path + ".partial");
  qrename_info(path + ".partial", path);
  timer.flops += file_size;
  return file_size;
}

template <class T = Real>
Long load_gauge_field(GaugeFieldT<T>& gf, const std::string& path)
{
  TIMER_VERBOSE_FLOPS("load_gauge_field");
  displayln_info(fname + ssprintf(": '%s'.", path.c_str()));
  gf.init();
  GaugeFieldInfo gfi;
  read_gauge_field_header(gfi, path);
  const bool is_two_row = gfi.datatype == "4D_SU3_GAUGE";
  const bool is_three_row = gfi.datatype == "4D_SU3_GAUGE_3x3";
  const int n_complex_su3 = is_two_row ? 6 : (is_three_row ? 9 : 0);
  if (n_complex_su3 == 0) {
    displayln(fname + ssprintf(": gfi.datatype '%s' id_node=%d.",
                               gfi.datatype.c_str(), get_id_node()));
    qassert(false);
  }
  Geometry geo;
  geo.init(gfi.total_site);
  const Int multiplicity = 4;
  Field<ComplexD> gft;
  gft.init(geo, multiplicity * n_complex_su3);
  const Long file_size = serial_read_field_par(
      gft, path, -get_data_size(gft) * get_num_node(), SEEK_END);
  if (0 == file_size) {
    return 0;
  }
  crc32_t crc32 = field_crc32(gft);
  if (crc32 != gfi.crc32) {
    if (get_id_node() == 0) {
      qwarn(fname + ssprintf(": WARNING: fn='%s' CHECKSUM= %08X (calc) %08X "
                             "(read) possibly missing CRC32HASH field",
                             path.c_str(), crc32, gfi.crc32));
    }
  }
  if (gfi.floating_point == "IEEE64BIG") {
    to_from_big_endian(get_data(gft));
  } else if (gfi.floating_point == "IEEE64LITTLE") {
    to_from_little_endian(get_data(gft));
  } else {
    qassert(false);
  }
  crc32_t simple_checksum =
      field_simple_checksum(gft);  // after endianness conversion
  if (simple_checksum != gfi.simple_checksum) {
    if (get_id_node() == 0) {
      qwarn(fname +
            ssprintf(": WARNING: fn='%s' CHECKSUM= %08X (calc) %08X (read)",
                     path.c_str(), simple_checksum, gfi.simple_checksum));
    }
  }
  gf.init(geo);
  qacc_for(index, geo.local_volume(), {
    const Geometry& geo = gf.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<ComplexD> vt = gft.get_elems(xl);
    Vector<ColorMatrixT<T>> v = gf.get_elems(xl);
    for (int m = 0; m < multiplicity; ++m) {
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
  });
  timer.flops += file_size;
  return file_size;
}

template <class T = Real>
inline Long load_gauge_field_cps3x3(GaugeFieldT<T>& gf,
                                    const std::string& path)
// assuming gf already initialized and have correct size;
{
  TIMER_VERBOSE_FLOPS("load_gauge_field_cps3x3");
  displayln_info(fname + ssprintf(": '%s'.", path.c_str()));
  qassert(is_initialized(gf));
  const Geometry& geo = gf.geo();
  FieldM<ComplexD, 4 * 9> gft;
  gft.init(geo);
  const Long file_size = serial_read_field_par(
      gft, path, -get_data_size(gft) * get_num_node(), SEEK_END);
  if (file_size == 0) {
    return 0;
  }
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<ComplexD> vt = gft.get_elems(xl);
    to_from_big_endian(vt);
    Vector<ColorMatrixT<T> > v = gf.get_elems(xl);
    for (int m = 0; m < gf.multiplicity; ++m) {
      assign_truncate(v[m], vt[m]);
    }
  }
  timer.flops += file_size;
  return file_size;
}

template <class T = Real>
inline Long load_gauge_field_milc(GaugeFieldT<T>& gf,
                                  const std::string& path,
                                  const bool par_read = false)
// assuming gf already initialized and have correct size;
{
  TIMER_VERBOSE_FLOPS("load_gauge_field_milc");
  displayln_info(fname + ssprintf(": '%s'.", path.c_str()));
  qassert(is_initialized(gf));
  const Geometry& geo = gf.geo();
  FieldM<ComplexF, 4 * 9> gft;
  gft.init(geo);
  // ADJUST ME
  Long file_size = 0;
  if (par_read) {
    file_size = serial_read_field_par(gft, path, 0x730, SEEK_SET);
  } else {
    file_size = serial_read_field(gft, path, 0x730, SEEK_SET);
  }
  if (0 == file_size) {
    return 0;
  }
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    Coordinate xl = geo.coordinate_from_index(index);
    Vector<ComplexF> vt = gft.get_elems(xl);
    to_from_big_endian(vt);
    Vector<ColorMatrixT<T> > v = gf.get_elems(xl);
    for (int m = 0; m < gf.multiplicity; ++m) {
      // assign_truncate(v[m], vt[m]);
      v[m](0, 0) = vt[9 * m + 0 * 3 + 0];
      v[m](0, 1) = vt[9 * m + 0 * 3 + 1];
      v[m](0, 2) = vt[9 * m + 0 * 3 + 2];
      v[m](1, 0) = vt[9 * m + 1 * 3 + 0];
      v[m](1, 1) = vt[9 * m + 1 * 3 + 1];
      v[m](1, 2) = vt[9 * m + 1 * 3 + 2];
      v[m](2, 0) = vt[9 * m + 2 * 3 + 0];
      v[m](2, 1) = vt[9 * m + 2 * 3 + 1];
      v[m](2, 2) = vt[9 * m + 2 * 3 + 2];
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
      mat *= ComplexT<T>(qpolar(1.0, amp));
    }
  }
}

// GaugeTransform IO

struct API GaugeTransformInfo {
  std::string hdr_version;
  std::string storage_format;
  Coordinate total_site;
  crc32_t simple_checksum;
  std::string floating_point;
  int data_per_site;
  std::string gf_type;
  double gf_accuracy;
  //
  GaugeTransformInfo() { init(); }
  //
  void init()
  {
    hdr_version = "1.0";
    storage_format = "1.0";
    total_site = Coordinate();
    simple_checksum = 0;
    floating_point = "IEEE64BIG";
    data_per_site = 18;
    gf_type = "COULOMB_T";
    gf_accuracy = 1e-14;
  }
};

inline std::string make_gauge_transform_header(
    const GaugeTransformInfo& info = GaugeTransformInfo())
{
  std::ostringstream out;
  // const std::string todo = "NOT yet implemented";
  out << "BEGIN_HEADER" << std::endl;
  out << "HDR_VERSION = " << info.hdr_version << std::endl;
  out << "STORAGE_FORMAT = " << info.storage_format << std::endl;
  out << "DIMENSION_1 = " << info.total_site[0] << std::endl;
  out << "DIMENSION_2 = " << info.total_site[1] << std::endl;
  out << "DIMENSION_3 = " << info.total_site[2] << std::endl;
  out << "DIMENSION_4 = " << info.total_site[3] << std::endl;
  out << "CHECKSUM = " << show_crc32(info.simple_checksum) << std::endl;
  out << "FLOATING_POINT = " << info.floating_point << std::endl;
  out << "DATA_PER_SITE = " << info.data_per_site << std::endl;
  out << "GF_TYPE = " << info.gf_type << std::endl;
  out << "GF_ACCURACY = " << info.gf_accuracy << std::endl;
  out << "END_HEADER" << std::endl;
  return out.str();
}

inline void read_gauge_transform_header(GaugeTransformInfo& info,
                                        const std::string& path)
{
  TIMER("read_gauge_transform_header");
  if (get_id_node() == 0) {
    QFile qfile = qfopen(path, "r");
    if (not qfile.null()) {
      const std::string header = "BEGIN_HEADER\n";
      std::vector<char> check_line(header.size(), 0);
      if (1 == qfread(check_line.data(), header.size(), 1, qfile)) {
        if (std::string(check_line.data(), check_line.size()) == header) {
          std::vector<std::string> infos;
          infos.push_back(header);
          while (infos.back() != "END_HEADER\n" && infos.back() != "") {
            infos.push_back(qgetline(qfile));
          }
          info.hdr_version =
              remove_trailing_newline(info_get_prop(infos, "HDR_VERSION = "));
          info.storage_format = remove_trailing_newline(
              info_get_prop(infos, "STORAGE_FORMAT = "));
          for (int m = 0; m < 4; ++m) {
            reads(info.total_site[m],
                  info_get_prop(infos, ssprintf("DIMENSION_%d = ", m + 1)));
          }
          info.simple_checksum =
              read_crc32(info_get_prop(infos, "CHECKSUM = "));
          info.floating_point = remove_trailing_newline(
              info_get_prop(infos, "FLOATING_POINT = "));
          info.data_per_site =
              read_long(info_get_prop(infos, "DATA_PER_SITE = "));
          info.gf_type =
              remove_trailing_newline(info_get_prop(infos, "GF_TYPE = "));
          info.gf_accuracy = read_double(
              info_get_prop(infos, "GF_ACCURACY = ", "GF_ACCRUACY = "));
        }
      }
    }
    qfclose(qfile);
  }
  bcast(info.hdr_version);
  bcast(info.storage_format);
  bcast(info.total_site);
  bcast(info.simple_checksum);
  bcast(info.floating_point);
  bcast(info.data_per_site);
  bcast(info.gf_type);
  bcast(info.gf_accuracy);
}

inline Long save_gauge_transform_cps(
    const GaugeTransform& gt, const std::string& path,
    const GaugeTransformInfo& info_ = GaugeTransformInfo())
{
  TIMER_VERBOSE_FLOPS("save_gauge_transform_cps");
  qassert(is_initialized(gt));
  GaugeTransform gt1;
  gt1 = gt;
  const Geometry& geo = gt1.geo();
  GaugeTransformInfo info = info_;
  info.total_site = geo.total_site();
  info.simple_checksum = field_simple_checksum(gt1); // before to_from_big_endian_64
  to_from_big_endian(get_data(gt1));
  qtouch_info(path + ".partial", make_gauge_transform_header(info));
  const Long file_size = serial_write_field(gt1, path + ".partial");
  qrename_info(path + ".partial", path);
  timer.flops += file_size;
  return file_size;
}

inline Long load_gauge_transform_cps(GaugeTransform& gt, const std::string& path)
// USE: read_field_double(gt, path) for qlat format GaugeTransform
{
  TIMER_VERBOSE_FLOPS("load_gauge_transform_cps");
  displayln_info(fname + ssprintf(": '%s'.", path.c_str()));
  gt.init();
  GaugeTransformInfo info;
  read_gauge_transform_header(info, path);
  qassert(info.data_per_site == 18);
  const Geometry geo(info.total_site);
  gt.init(geo, 1);
  const Long file_size = serial_read_field_par(
      gt, path, -get_data_size(gt) * get_num_node(), SEEK_END);
  if (0 == file_size) {
    displayln_info(fname + ssprintf(": failed to read any content."));
    gt.init();
    return 0;
  }
  if (info.floating_point == "IEEE64BIG") {
    to_from_big_endian(get_data(gt));
  } else if (info.floating_point == "IEEE64LITTLE") {
    to_from_little_endian(get_data(gt));
  } else {
    qassert(false);
  }
  crc32_t simple_checksum = field_simple_checksum(gt); // after endianness conversion
  if (simple_checksum != info.simple_checksum) {
    if (get_id_node() == 0) {
      qwarn(fname +
            ssprintf(": WARNING: fn='%s' CHECKSUM= %08X (calc) %08X (read)",
                     path.c_str(), simple_checksum, info.simple_checksum));
    }
  }
  timer.flops += file_size;
  return file_size;
}

}  // namespace qlat
