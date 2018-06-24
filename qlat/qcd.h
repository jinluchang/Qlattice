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
  virtual const std::string& cname()
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

struct FermionField5d : Field<WilsonVector>
{
  virtual const std::string& cname()
  {
    static const std::string s = "FermionField5d";
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

inline std::string make_gauge_field_header(const GaugeFieldInfo& gfi = GaugeFieldInfo())
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

inline void save_gauge_field(const GaugeField& gf, const std::string& path, const GaugeFieldInfo& gfi_ = GaugeFieldInfo())
{
  TIMER_VERBOSE_FLOPS("save_gauge_field");
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
  timer.flops += get_data(gft).data_size() * gft.geo.geon.num_node;
}

inline long load_gauge_field(GaugeField& gf, const std::string& path)
  // assuming gf already initialized and have correct size;
{
  TIMER_VERBOSE_FLOPS("load_gauge_field");
  displayln_info(fname + ssprintf(": '%s'.", path.c_str()));
  qassert(is_initialized(gf));
  const Geometry& geo = gf.geo;
  FieldM<std::array<Complex, 6>, 4> gft;
  gft.init(geo);
  const long file_size = serial_read_field(gft, path, -get_data_size(gft) * get_num_node(), SEEK_END);
  if (0 == file_size) {
    return 0;
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
  timer.flops += file_size;
  return file_size;
}

inline long load_gauge_field_par(GaugeField& gf, const std::string& path)
  // assuming gf already initialized and have correct size;
{
  TIMER_VERBOSE_FLOPS("load_gauge_field_par");
  displayln_info(fname + ssprintf(": '%s'.", path.c_str()));
  qassert(is_initialized(gf));
  const Geometry& geo = gf.geo;
  FieldM<std::array<Complex, 6>, 4> gft;
  gft.init(geo);
  const long file_size = serial_read_field_par(gft, path, -get_data_size(gft) * get_num_node(), SEEK_END);
  if (file_size == 0) {
    return 0;
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
  timer.flops += file_size;
  return file_size;
}

inline long load_gauge_field_cps3x3(GaugeField& gf, const std::string& path)
  // assuming gf already initialized and have correct size;
{
  TIMER_VERBOSE_FLOPS("load_gauge_field");
  displayln_info(fname + ssprintf(": '%s'.", path.c_str()));
  qassert(is_initialized(gf));
  const Geometry& geo = gf.geo;
  FieldM<std::array<Complex, 9>, 4> gft;
  gft.init(geo);
  const long file_size = serial_read_field(gft, path, -get_data_size(gft) * get_num_node(), SEEK_END);
  if (file_size == 0) {
    return 0;
  }
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<std::array<Complex, 9> > vt = gft.get_elems(xl);
    to_from_big_endian_64(get_data(vt));
    Vector<ColorMatrix> v = gf.get_elems(xl);
    for (int m = 0; m < geo.multiplicity; ++m) {
      assign_truncate(v[m], vt[m]);
    }
  }
  timer.flops += file_size;
  return file_size;
}

inline long load_gauge_field_milc(GaugeField& gf, const std::string& path, const bool par_read = false)
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
  timer.flops += file_size;
  return file_size;
}

inline void twist_boundary_at_boundary(GaugeField& gf, double mom, int mu)
{
  TIMER_VERBOSE_FLOPS("twist_boundary_at_boundary");
  const Geometry& geo = gf.geo;
  const double amp = 2.0 * PI * mom;
  const int len = geo.total_site()[mu];
  for (int index = 0; index < geo.local_volume(); index++) {
    Coordinate xl = geo.coordinate_from_index(index);
    Coordinate xg = geo.coordinate_g_from_l(xl);
    if (xg[mu] == len - 1) {
      ColorMatrix& mat = gf.get_elem(xl, mu);
      mat *= std::polar(1.0, amp);
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
  const Geometry& geo = ffs[0].geo;
  for (int i = 0; i < 4*NUM_COLOR; ++i) {
    qassert(geo == ffs[i].geo);
  }
  prop.init(geo_reform(geo));
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

inline void set_wilson_matrix_col_from_vector(WilsonMatrix& wm, const int idx, const WilsonVector& col)
{
  for (int j = 0; j < 4*NUM_COLOR; ++j) {
    wm(j, idx) = col(j);
  }
}

inline void set_propagator_col_from_fermion_field(Propagator4d& prop, const int idx, const FermionField4d& ff)
{
  TIMER("set_propagator_col_from_fermion_field");
  const Geometry& geo = ff.geo;
  qassert(geo == prop.geo);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    set_wilson_matrix_col_from_vector(prop.get_elem(xl), idx, ff.get_elem(xl));
  }
}

inline void set_wilson_vector_from_matrix_col(WilsonVector& col, const WilsonMatrix& wm, const int idx)
{
  for (int j = 0; j < 4*NUM_COLOR; ++j) {
    col(j) = wm(j, idx);
  }
}

inline void set_fermion_field_from_propagator_col(FermionField4d& ff, const Propagator4d& prop, const int idx)
{
  TIMER("set_fermion_field_from_propagator_col");
  const Geometry& geo = prop.geo;
  ff.init(geo_reform(geo));
  qassert(geo == ff.geo);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    set_wilson_vector_from_matrix_col(ff.get_elem(xl), prop.get_elem(xl), idx);
  }
}

inline void fermion_field_5d_from_4d(FermionField5d& ff5d, const FermionField4d& ff4d, const int upper, const int lower)
  // ff5d need to be initialized
  // upper componets are right handed
  // lower componets are left handed
{
  TIMER("fermion_field_5d_from_4d");
  const Geometry& geo = ff5d.geo;
  set_zero(ff5d);
  const int sizewvh = sizeof(WilsonVector) / 2;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    Coordinate x = geo.coordinate_from_index(index);
    memcpy((char*)&(ff5d.get_elem(x, upper)),
        (const char*)&(ff4d.get_elem(x)),
        sizewvh);
    memcpy((char*)&(ff5d.get_elem(x, lower)) + sizewvh,
        (const char*)&(ff4d.get_elem(x)) + sizewvh,
        sizewvh);
  }
}

inline void fermion_field_4d_from_5d(FermionField4d& ff4d, const FermionField5d& ff5d, const int upper, const int lower)
  // upper componets are right handed
  // lower componets are left handed
{
  TIMER("fermion_field_4d_from_5d");
  const Geometry& geo = ff5d.geo;
  ff4d.init(geo_reform(geo));
  set_zero(ff4d);
  const int sizewvh = sizeof(WilsonVector) / 2;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    Coordinate x = geo.coordinate_from_index(index);
    memcpy((char*)&(ff4d.get_elem(x)),
        (const char*)&(ff5d.get_elem(x, upper)),
        sizewvh);
    memcpy((char*)&(ff4d.get_elem(x)) + sizewvh,
        (const char*)&(ff5d.get_elem(x, lower)) + sizewvh,
        sizewvh);
  }
}

template <class Inverter>
inline void inverse_dwf(FermionField4d& sol, const FermionField4d& src, const Inverter& inv)
  // sol do not need to be initialized
  // inv.geo must be the geometry of the fermion field
  // inverse(sol5d, src5d, inv) perform the inversion
{
  TIMER_VERBOSE("inverse_dwf(4d,4d,inv)");
  const Geometry& geo = src.geo;
  sol.init(geo);
  const Geometry geo_ls = geo_reform(inv.geo, inv.fa.ls, 0);
  const int ls = geo_ls.multiplicity;
  FermionField5d sol5d, src5d;
  sol5d.init(geo_ls);
  src5d.init(geo_ls);
  fermion_field_5d_from_4d(src5d, src, 0, ls-1);
  fermion_field_5d_from_4d(sol5d, sol, ls-1, 0);
  inverse(sol5d, src5d, inv);
  fermion_field_4d_from_5d(sol, sol5d, ls-1, 0);
}

template <class Inverter>
inline void inverse(Propagator4d& sol, const Propagator4d& src, const Inverter& inv)
  // sol do not need to be initialized
  // inv.geo must be the geometry of the fermion field
  // inverse(4d, 4d, inv) perform the inversion
{
  TIMER_VERBOSE("inverse(p4d,p4d,inv)");
  const Geometry& geo = geo_reform(src.geo);
  sol.init(geo);
  FermionField4d ff_sol, ff_src;
  for (int j = 0; j < 4*NUM_COLOR; ++j) {
    set_fermion_field_from_propagator_col(ff_src, src, j);
    set_zero(ff_sol);
    inverse(ff_sol, ff_src, inv);
    set_propagator_col_from_fermion_field(sol, j, ff_sol);
  }
}

inline void set_fermion_field_point_src(FermionField4d& ff, const Coordinate& xg, const int cs, const Complex& value = 1.0)
  // ff need to be initialized
{
  TIMER("set_fermion_field_point_src");
  const Geometry& geo = ff.geo;
  set_zero(ff);
  const Coordinate xl = geo.coordinate_l_from_g(xg);
  if (geo.is_local(xl)) {
    ff.get_elem(xl)(cs) = value;
  }
}

inline void set_point_src(Propagator4d& prop, const Geometry& geo_input, const Coordinate& xg, const Complex& value = 1.0)
{
  TIMER_VERBOSE("set_point_src");
  const Geometry geo = geo_reform(geo_input);
  prop.init(geo);
  FermionField4d src;
  src.init(geo);
  for (int cs = 0; cs < 4 * NUM_COLOR; ++cs) {
    set_fermion_field_point_src(src, xg, cs, value);
    set_propagator_col_from_fermion_field(prop, cs, src);
  }
}

template <class Inverter>
inline void set_point_src_propagator(Propagator4d& prop, const Inverter& inv, const Coordinate& xg, const Complex& value = 1.0)
{
  TIMER_VERBOSE("set_point_src_propagator");
  const Geometry geo = geo_reform(inv.geo);
  prop.init(geo);
  FermionField4d sol, src;
  sol.init(geo);
  src.init(geo);
  for (int cs = 0; cs < 4 * NUM_COLOR; ++cs) {
    set_fermion_field_point_src(src, xg, cs, value);
    set_zero(sol);
    inverse(sol, src, inv);
    set_propagator_col_from_fermion_field(prop, cs, sol);
  }
}

inline CoordinateD lattice_mom_mult(const Coordinate& total_site)
{
  return 2*PI / CoordinateD(total_site);
}

inline CoordinateD lattice_mom_mult(const Geometry& geo)
{
  return lattice_mom_mult(geo.total_site());
}

inline void set_mom_src_fermion_field(FermionField4d& ff, const CoordinateD& lmom, const int cs)
  // ff need to be initialized beforehand
{
  const Geometry& geo = ff.geo;
  const CoordinateD mom = lmom * lattice_mom_mult(geo);
  set_zero(ff);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    double phase = 0.0;
    for (int i = 0; i < DIMN; ++i) {
      phase += mom[i] * xg[i];
    }
    ff.get_elem(xl)(cs) = std::polar(1.0, phase);
  }
}

inline void set_tslice_mom_src_fermion_field(FermionField4d& ff, const int tslice, const CoordinateD& lmom, const int cs)
  // ff need to be initialized beforehand
{
  qassert(lmom[3] == 0);
  const Geometry& geo = ff.geo;
  const CoordinateD mom = lmom * lattice_mom_mult(geo);
  set_zero(ff);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    if (xg[3] == tslice) {
      double phase = 0.0;
      for (int i = 0; i < DIMN; ++i) {
        phase += mom[i] * xg[i];
      }
      ff.get_elem(xl)(cs) = std::polar(1.0, phase);
    }
  }
}

template <class Inverter>
inline void set_mom_src_propagator(Propagator4d& prop, const CoordinateD& lmom, Inverter& inverter)
{
  TIMER_VERBOSE("set_mom_src_propagator");
  const Geometry& geo = geo_remult(inverter.geo);
  prop.init(geo);
  qassert(prop.geo == geo);
  FermionField4d src, sol;
  src.init(geo);
  sol.init(geo);
  for (int cs = 0; cs < 4*NUM_COLOR; ++cs) {
    set_mom_src_fermion_field(src, lmom, cs);
    set_zero(sol);
    inverse(sol, src, inverter);
    set_propagator_col_from_fermion_field(prop, cs, sol);
  }
}

template <class Inverter>
inline void set_tslice_mom_src_propagator(Propagator4d& prop, const int tslice, const CoordinateD& lmom, Inverter& inverter)
{
  TIMER_VERBOSE("set_tslice_mom_src_propagator");
  const Geometry& geo = geo_remult(inverter.geo);
  prop.init(geo);
  qassert(prop.geo == geo);
  FermionField4d src, sol;
  src.init(geo);
  sol.init(geo);
  for (int cs = 0; cs < 4*NUM_COLOR; ++cs) {
    set_tslice_mom_src_fermion_field(src, tslice, lmom, cs);
    set_zero(sol);
    inverse(sol, src, inverter);
    set_propagator_col_from_fermion_field(prop, cs, sol);
  }
}

inline void smear_propagator(Propagator4d& prop, const GaugeField& gf1,
    const double coef, const int step, const CoordinateD& mom = CoordinateD(), const bool smear_in_time_dir = false)
  // gf1 is left_expanded and refreshed
  // prop is of normal size
{
  TIMER_VERBOSE("smear_propagator");
  if (0 == step) {
    return;
  }
  const Geometry& geo = prop.geo;
  const Geometry geo1 = smear_in_time_dir
    ? geo_resize(geo, 1)
    : geo_resize(geo, Coordinate(1,1,1,0), Coordinate(1,1,1,0));
  const int n_avg = smear_in_time_dir ? 8 : 6;
  const int dir_limit = smear_in_time_dir ? 4 : 3;
  std::array<Complex,8> mom_factors;
  for (int i = 0; i < 8; ++i) {
    const int dir = i - 4;
    const double phase = dir >= 0 ? mom[dir] : -mom[-dir-1];
    mom_factors[i] = std::polar(coef/n_avg, -phase);
  }
  Propagator4d prop1;
  prop1.init(geo1);
  for (int i = 0; i < step; ++i) {
    prop1 = prop;
    refresh_expanded_1(prop1);
#pragma omp parallel for
    for (long index = 0; index < geo.local_volume(); ++index) {
      const Coordinate xl = geo.coordinate_from_index(index);
      WilsonMatrix& wm = prop.get_elem(xl);
      wm *= 1-coef;
      for (int dir = -dir_limit; dir < dir_limit; ++dir) {
        const Coordinate xl1 = coordinate_shifts(xl, dir);
        const ColorMatrix link = dir >= 0
          ? gf1.get_elem(xl, dir)
          : (ColorMatrix)matrix_adjoint(gf1.get_elem(coordinate_shifts(xl, dir), -dir-1));
        wm += mom_factors[dir+4] * link * prop1.get_elem(xl1);
      }
    }
  }
}

QLAT_END_NAMESPACE
