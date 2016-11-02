#pragma once

#include <qlat/field.h>
#include <qlat/field-fft.h>
#include <qlat/qed.h>

#include <Eigen/Eigen>

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
    Coordinate xl; geo.coordinate_from_index(xl, index);
    Vector<ColorMatrix> v = gf.get_elems(xl);
    for (int m = 0; m < geo.multiplicity; ++m) {
      unitarize(v[m]);
    }
  }
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
    Coordinate xl; geo.coordinate_from_index(xl, index);
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
