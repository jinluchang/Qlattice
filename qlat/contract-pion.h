#pragma once

#include <qlat/qcd.h>
#include <qlat/selected-field.h>
#include <qlat/selected-points.h>

namespace qlat
{  //

typedef SelectedField<WilsonMatrix> SelProp;

typedef SelectedPoints<WilsonMatrix> PselProp;

inline LatData contract_pion(const Propagator4d& prop, const int tslice_src)
{
  TIMER_VERBOSE("contract_pion(prop)");
  const Geometry& geo = prop.geo;
  const Coordinate total_site = geo.total_site();
  LatData ld;
  ld.info.push_back(lat_dim_number("tsep", 0, total_site[3] - 1));
  ld.info.push_back(lat_dim_re_im());
  lat_data_alloc(ld);
  set_zero(ld);
  Vector<Complex> ldv = lat_data_complex_get(ld, make_array<int>());
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const Complex val = qnorm(prop.get_elem(xl));
    const int tsep = mod(xg[3] - tslice_src, total_site[3]);
    ldv[tsep] += val;
  }
  glb_sum_lat_data(ld);
  return ld;
}

inline LatData contract_pion(const SelProp& prop, const FieldSelection& fsel,
                             const int tslice_src)
{
  TIMER_VERBOSE("contract_pion(s_prop,fsel)");
  const Geometry& geo = prop.geo;
  const Coordinate total_site = geo.total_site();
  LatData ld;
  ld.info.push_back(lat_dim_number("tsep", 0, total_site[3] - 1));
  ld.info.push_back(lat_dim_re_im());
  lat_data_alloc(ld);
  set_zero(ld);
  Vector<Complex> ldv = lat_data_complex_get(ld, make_array<int>());
  for (long idx = 0; idx < fsel.n_elems; ++idx) {
    const long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const Complex val = qnorm(prop.get_elem(idx));
    const int tsep = mod(xg[3] - tslice_src, total_site[3]);
    ldv[tsep] += val;
  }
  glb_sum_lat_data(ld);
  ld *= 1.0 / fsel.prob;
  return ld;
}

inline LatData contract_pion(const PselProp& prop, const PointSelection& psel,
                             const Geometry& geo, const int tslice_src)
{
  TIMER_VERBOSE("contract_pion(ps_prop,psel)");
  const long n_points = prop.n_points;
  qassert(n_points == (long)psel.size());
  const Coordinate total_site = geo.total_site();
  LatData ld;
  ld.info.push_back(lat_dim_number("tsep", 0, total_site[3] - 1));
  ld.info.push_back(lat_dim_re_im());
  lat_data_alloc(ld);
  set_zero(ld);
  Vector<Complex> ldv = lat_data_complex_get(ld, make_array<int>());
  for (long idx = 0; idx < n_points; ++idx) {
    const Coordinate& xg = psel[idx];
    const Complex val = qnorm(prop.get_elem(idx));
    const int tsep = mod(xg[3] - tslice_src, total_site[3]);
    ldv[tsep] += val;
  }
  ld *= (double)product(total_site) / (double)n_points;
  return ld;
}

}  // namespace qlat
