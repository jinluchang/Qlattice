#pragma once

#include <qlat/hmc.h>
#include <qlat/wilson-flow.h>

namespace qlat
{  //

inline double get_field_max(const FieldM<double, 1>& fd)
{
  TIMER("get_field_max");
  const Geometry& geo = fd.geo();
  qassert(fd.geo().is_only_local());
  double m = fd.get_elem(0);
  for (long index = 1; index < geo.local_volume(); ++index) {
    m = std::max(m, fd.get_elem(index));
  }
  std::vector<double> ms(get_num_node(), 0.0);
  ms[get_id_node()] = m;
  glb_sum(get_data(ms));
  for (int i = 0; i < (int)ms.size(); ++i) {
    m = std::max(m, ms[i]);
  }
  return m;
}

inline std::vector<double> get_gm_force_magnitudes(
    const GaugeMomentum& gm_force, const int n_elems)
// return the l1, l2, ..., linf norm of the gm_force magnitudes
// n_elems == mag_vec.size();
// n_elems >= 3
{
  TIMER("get_gm_force_magnitudes");
  qassert(n_elems >= 2);
  const Geometry geo = geo_reform(gm_force.geo(), n_elems - 1);
  Field<double> fd;
  fd.init(geo);
  set_zero(fd);
  FieldM<double, 1> fd_max;
  fd_max.init(geo);
  set_zero(fd_max);
  qacc_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<ColorMatrix> gm_force_v = gm_force.get_elems_const(xl);
    qassert(gm_force_v.size() == 4);
    Vector<double> fdv = fd.get_elems(index);
    double linf = 0.0;
    for (int mu = 0; mu < 4; ++mu) {
      // multiply by additional small factor, will be compensated below
      const double l2 =
          sqr(1.0 / 15.0) * 2.0 * neg_half_tr_square(gm_force_v[mu]);
      const double l1 = std::sqrt(l2);
      fdv[0] += l1;
      fdv[1] += l2;
      double ln = l2;
      for (int m = 2; m < n_elems - 1; ++m) {
        ln = sqr(ln);
        fdv[m] += ln;
      }
      linf = std::max(linf, l1);
    }
    fd_max.get_elem(index) = 15.0 * linf;
  });
  std::vector<double> mag_vec(n_elems, 0.0);
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Vector<double> fdv = fd.get_elems_const(index);
    for (int m = 0; m < n_elems - 1; ++m) {
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
  mag_vec[n_elems - 1] = get_field_max(fd_max);
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

API inline std::vector<std::vector<double> >& get_gm_force_magnitudes_list()
{
  static std::vector<std::vector<double> > gm_force_magnitudes_list;
  return gm_force_magnitudes_list;
}

inline void display_gm_force_magnitudes(const GaugeMomentum& gm_force,
                                        const int n_elems)
// e.g. n_elems = 10
// Need to call: clear(get_gm_force_magnitudes_list()) to free memory
{
  TIMER_VERBOSE("display_gm_force_magnitudes");
  const std::vector<double> mag_vec =
      get_gm_force_magnitudes(gm_force, n_elems);
  display_info(show_gm_force_magnitudes(mag_vec));
  get_gm_force_magnitudes_list().push_back(mag_vec);
}

inline void save_gm_force_magnitudes_list(const std::string& fn = "")
{
  TIMER_VERBOSE("save_gm_force_magnitudes_list");
  std::vector<std::vector<double> >& db = get_gm_force_magnitudes_list();
  if (db.size() > 0 and fn != "") {
    const long idx_size = db.size();
    const long ln_size = db[0].size();
    LatData ld;
    ld.info.push_back(lat_dim_number("idx", 0, idx_size - 1));
    ld.info.push_back(lat_dim_number("ln", 0, ln_size - 1));
    lat_data_alloc(ld);
    for (long i = 0; i < idx_size; ++i) {
      for (long j = 0; j < ln_size; ++j) {
        lat_data_get(ld, make_array<long>(i, j))[0] = db[i][j];
      }
    }
    lat_data_save_info(fn, ld);
  }
  clear(get_gm_force_magnitudes_list());
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
  const Geometry geo = geo_reform(gf.geo());
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

inline std::string show_gauge_field_info_line(const int i,
                                              const std::vector<double>& v)
{
  qassert(v.size() == 4);
  return ssprintf("%5d %24.17E %24.17E %24.17E %24.17E", i, v[0], v[1], v[2],
                  v[3]);
}

inline LatData convert_gauge_field_info_table(
    const std::vector<std::vector<double> >& dt)
{
  TIMER("convert_gauge_field_info_table");
  LatData ld;
  ld.info.push_back(lat_dim_number("smear_step", 0, dt.size() - 1));
  ld.info.push_back(lat_dim_string(
      "type",
      make_array<std::string>("avg(plaq_action)", "avg(plaq_action^2)",
                              "tot(topo_density)", "avg(topo_density^2)")));
  lat_data_alloc(ld);
  for (int i = 0; i < (int)dt.size(); ++i) {
    const std::vector<double>& v = dt[i];
    Vector<double> ldv = lat_data_get(ld, make_array<int>(i));
    for (int j = 0; j < (int)v.size(); ++j) {
      ldv[j] = v[j];
    }
  }
  return ld;
}

inline LatData convert_energy_list(
    const std::vector<double>& energy_density_list, const double wilson_flow_step)
{
  TIMER("convert_energy_list");
  LatData ld;
  ld.info.push_back(lat_dim_number("flow_steps", 1, energy_density_list.size()));
  ld.info.push_back(lat_dim_string(
      "name", make_array<std::string>("flow_time", "energy_density")));
  lat_data_alloc(ld);
  for (long i = 0; i < (long)energy_density_list.size(); ++i) {
    Vector<double> ldv = lat_data_get(ld, make_array<long>(i));
    ldv[0] = (i + 1) * wilson_flow_step;
    ldv[1] = energy_density_list[i];
  }
  return ld;
}

inline LatData get_gauge_field_info_table_with_ape_smear(const GaugeField& gf,
                                                         const double alpha,
                                                         const int ape_steps,
                                                         const int steps)
// DataTable.size() == steps + 1
//
{
  TIMER("get_gauge_field_info_table_with_ape_smear");
  std::vector<std::vector<double> > dt;
  std::vector<double> v = get_gauge_field_infos(gf);
  displayln_info(show_gauge_field_info_line(0, v));
  dt.push_back(v);
  GaugeField gf1;
  gf1 = gf;
  for (int i = 0; i < steps; ++i) {
    gf_ape_smear(gf1, gf1, alpha, ape_steps);
    v = get_gauge_field_infos(gf1);
    displayln_info(show_gauge_field_info_line(i + 1, v));
    dt.push_back(v);
  }
  return convert_gauge_field_info_table(dt);
}

inline void display_gauge_field_info_table_with_ape_smear(const std::string& fn,
                                                          const GaugeField& gf,
                                                          const double alpha,
                                                          const int ape_steps,
                                                          const int steps)
// e.g. alpha = 0.5; steps = 50;
{
  TIMER_VERBOSE("display_gauge_field_info_table_with_ape_smear");
  const LatData ld =
      get_gauge_field_info_table_with_ape_smear(gf, alpha, ape_steps, steps);
  display_info(show(ld));
  if (fn != "") {
    lat_data_save_info(fn, ld);
  }
}

inline void get_gauge_field_info_table_with_wilson_flow(
    LatData& ld_gf_info, LatData& ld_wilson_flow_energy, const GaugeField& gf,
    const double flow_time, const int flow_steps, const int steps,
    const double c1 = 0.0)
// DataTable.size() == steps + 1
//
{
  TIMER("get_gauge_field_info_table_with_wilson_flow");
  std::vector<std::vector<double> > dt;
  std::vector<double> energy_density_list;
  std::vector<double> v = get_gauge_field_infos(gf);
  displayln_info(show_gauge_field_info_line(0, v));
  dt.push_back(v);
  GaugeField gf1;
  gf1 = gf;
  double existing_flow_time = 0.0;
  for (int i = 0; i < steps; ++i) {
    const std::vector<double> edl =
        gf_wilson_flow(gf1, existing_flow_time, flow_time, flow_steps, c1);
    existing_flow_time += flow_time;
    v = get_gauge_field_infos(gf1);
    displayln_info(show_gauge_field_info_line(i + 1, v));
    dt.push_back(v);
    vector_append(energy_density_list, edl);
  }
  ld_gf_info = convert_gauge_field_info_table(dt);
  ld_wilson_flow_energy =
      convert_energy_list(energy_density_list, flow_time / (double)flow_steps);
}

inline void display_gauge_field_info_table_with_wilson_flow(
    const std::string& fn_gf_info, const std::string& fn_wilson_flow_energy,
    const GaugeField& gf, const double flow_time, const int flow_steps,
    const int steps, const double c1 = 0.0)
// e.g. alpha = 0.5; steps = 50;
{
  TIMER_VERBOSE("display_gauge_field_info_table_with_wilson_flow");
  LatData ld_gf_info, ld_wilson_flow_energy;
  get_gauge_field_info_table_with_wilson_flow(
      ld_gf_info, ld_wilson_flow_energy, gf, flow_time, flow_steps, steps, c1);
  display_info(show(ld_wilson_flow_energy));
  display_info(show(ld_gf_info));
  if (fn_gf_info != "") {
    lat_data_save_info(fn_gf_info, ld_gf_info);
  }
  if (fn_wilson_flow_energy != "") {
    lat_data_save_info(fn_wilson_flow_energy, ld_wilson_flow_energy);
  }
}

}  // namespace qlat
