#pragma once

#include <qlat/hmc.h>
#include <qlat/wilson-flow.h>

namespace qlat
{  //

inline std::vector<double> get_gm_force_magnitudes(
    const GaugeMomentum& gm_force, const int n_elems)
// return the l1, l2, l4, l8, l16, ... norm of the gm_force magnitudes
// n_elems == mag_vec.size();
{
  TIMER("get_gm_force_magnitudes");
  qassert(n_elems >= 2);
  const Geometry geo = geo_reform(gm_force.geo, n_elems);
  Field<double> fd;
  fd.init(geo);
  set_zero(fd);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<ColorMatrix> gm_force_v = gm_force.get_elems_const(xl);
    qassert(gm_force_v.size() == 4);
    Vector<double> fdv = fd.get_elems(xl);
    for (int mu = 0; mu < 4; ++mu) {
      // multiply by additional small factor, will be compensated below
      const double l2 =
          sqr(1.0 / 15.0) * 2.0 * neg_half_tr_square(gm_force_v[mu]);
      const double l1 = std::sqrt(l2);
      fdv[0] += l1;
      fdv[1] += l2;
      double ln = l2;
      for (int m = 2; m < n_elems; ++m) {
        ln = sqr(ln);
        fdv[m] += ln;
      }
    }
  }
  std::vector<double> mag_vec(n_elems, 0.0);
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Vector<double> fdv = fd.get_elems_const(index);
    for (int m = 0; m < n_elems; ++m) {
      mag_vec[m] += fdv[m];
    }
  }
  glb_sum(get_data(mag_vec));
  for (int m = 0; m < n_elems; ++m) {
    mag_vec[m] *= (1.0 / 4.0) / (double)geo.total_volume();
    for (int i = 0; i < m; ++i) {
      mag_vec[m] = std::sqrt(mag_vec[m]);
    }
    // The compensate the additional factor introduced above
    mag_vec[m] *= 15.0;
  }
  return mag_vec;
}

inline std::string show_gm_force_magnitudes(const std::vector<double>& mag_vec)
{
  std::ostringstream out;
  out << "# norm_idx gm_force_magnitude" << std::endl;
  for (int i = 0; i < (int)mag_vec.size(); ++i) {
    out << ssprintf("%5d %24.17E", i, mag_vec[i]) << std::endl;
  }
  return out.str();
}

inline void display_gm_force_magnitudes(const GaugeMomentum& gm_force,
                                        const int n_elems)
// e.g. n_elems = 10
{
  TIMER_VERBOSE("display_gm_force_magnitudes");
  const std::vector<double> mag_vec =
      get_gm_force_magnitudes(gm_force, n_elems);
  display_info(show_gm_force_magnitudes(mag_vec));
}

inline std::vector<double> get_gauge_field_infos(const GaugeField& gf)
// Values:
// 0: avg(plaq_action)
// 1: avg(plaq_action^2)
// 2: tot(topo_density)
// 3: avg(topo_density^2)
{
  TIMER("get_gauge_field_infos");
  const int info_vec_size = 4;
  const Geometry geo = geo_reform(gf.geo);
  CloverLeafField clf1, clf2, clf3, clf4, clf5;
  gf_clover_leaf_field_5(clf1, clf2, clf3, clf4, clf5, gf);
  FieldM<double, 1> paf;
  clf_plaq_action_field(paf, clf1);
  FieldM<double, 1> topf;
  clf_topology_field_5(topf, clf1, clf2, clf3, clf4, clf5);
  std::vector<double> info_vec(info_vec_size, 0.0);
  for (long index = 0; index < geo.local_volume(); ++index) {
    const double pa = paf.get_elem(index);
    const double top = topf.get_elem(index);
    info_vec[0] += pa;
    info_vec[1] += sqr(pa);
    info_vec[2] += top;
    info_vec[3] += sqr(top);
  }
  glb_sum(get_data(info_vec));
  info_vec[0] *= 1.0 / (double)geo.total_volume();
  info_vec[1] *= 1.0 / (double)geo.total_volume();
  info_vec[3] *= 1.0 / (double)geo.total_volume();
  return info_vec;
}

inline std::string show_gauge_field_info_table(
    const std::vector<std::vector<double> >& dt)
{
  std::ostringstream out;
  out << "# smear_step avg(plaq_action) avg(plaq_action^2) tot(topo_density) "
         "avg(topo_density^2)"
      << std::endl;
  for (int i = 0; i < (int)dt.size(); ++i) {
    const std::vector<double>& v = dt[i];
    qassert(v.size() == 4);
    out << ssprintf("%5d %24.17E %24.17E %24.17E %24.17E", i, v[0], v[1], v[2],
                    v[3])
        << std::endl;
  }
  return out.str();
}

inline std::vector<std::vector<double> >
get_gauge_field_info_table_with_ape_smear(const GaugeField& gf,
                                          const double alpha, const int steps)
// DataTable.size() == steps + 1
//
{
  TIMER("get_gauge_field_info_table_with_ape_smear");
  std::vector<std::vector<double> > dt;
  dt.push_back(get_gauge_field_infos(gf));
  GaugeField gf1;
  gf1 = gf;
  for (int i = 0; i < steps; ++i) {
    gf_ape_smear(gf1, gf1, alpha);
    dt.push_back(get_gauge_field_infos(gf1));
  }
  return dt;
}

inline void display_gauge_field_info_table_with_ape_smear(const std::string& fn,
                                                          const GaugeField& gf,
                                                          const double alpha,
                                                          const int steps)
// e.g. alpha = 0.5; steps = 50;
{
  TIMER_VERBOSE("display_gauge_field_info_table_with_ape_smear(fn)");
  const std::vector<std::vector<double> > dt =
      get_gauge_field_info_table_with_ape_smear(gf, alpha, steps);
  const std::string str = show_gauge_field_info_table(dt);
  display_info(str);
  if (fn != "") {
    qtouch_info(fn, str);
  }
}

inline std::vector<std::vector<double> >
get_gauge_field_info_table_with_wilson_flow(const GaugeField& gf,
                                            const double flow_time,
                                            const int flow_steps,
                                            const int steps,
                                            const double c1 = 0.0)
// DataTable.size() == steps + 1
//
{
  TIMER("get_gauge_field_info_table_with_wilson_flow");
  std::vector<std::vector<double> > dt;
  dt.push_back(get_gauge_field_infos(gf));
  GaugeField gf1;
  gf1 = gf;
  for (int i = 0; i < steps; ++i) {
    gf_wilson_flow(gf1, flow_time, flow_steps, c1);
    dt.push_back(get_gauge_field_infos(gf1));
  }
  return dt;
}

inline void display_gauge_field_info_table_with_wilson_flow(
    const std::string& fn, const GaugeField& gf, const double flow_time,
    const int flow_steps, const int steps, const double c1 = 0.0)
// e.g. alpha = 0.5; steps = 50;
{
  TIMER_VERBOSE("display_gauge_field_info_table_with_wilson_flow(fn)");
  const std::vector<std::vector<double> > dt =
      get_gauge_field_info_table_with_wilson_flow(gf, flow_time, flow_steps,
                                                  steps, c1);
  const std::string str = show_gauge_field_info_table(dt);
  display_info(str);
  if (fn != "") {
    qtouch_info(fn, str);
  }
}

}  // namespace qlat
