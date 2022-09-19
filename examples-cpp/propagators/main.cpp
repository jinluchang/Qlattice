#include "qlat-setup.h"

namespace qlat
{  //

typedef InverterDomainWall Inverter;

// -----------------------------------------------------------------------------------

inline void set_sparse_parameters(FieldSelection& fsel,
                                  const std::string& job_tag, const int traj)
{
  TIMER_VERBOSE("set_sparse_parameters");
  const std::string job_path = get_job_path(job_tag, traj);
  const std::string f_rank_path = job_path + "/f-rank.field";
  const Coordinate total_site = get_total_site(job_tag);
  const long spatial_vol = total_site[0] * total_site[1] * total_site[2];
  const long n_per_tslice = spatial_vol / 16;
  const RngState rs_sel = RngState("field-sel").split(job_tag).split(traj);
  set_field_selection(fsel, total_site, n_per_tslice, rs_sel);
  if (not does_file_exist_sync_node(f_rank_path)) {
    write_field_selection(fsel, f_rank_path);
  }
  FieldSelection fsel_load;
  read_field_selection(fsel_load, f_rank_path, n_per_tslice);
  FieldM<int64_t, 1> f_rank;
  f_rank = fsel_load.f_rank;
  f_rank -= fsel.f_rank;
  qassert(qnorm(f_rank) == 0.0);
}

inline bool compute_traj_do(const std::string& job_tag, const int traj)
{
  TIMER_VERBOSE("compute_traj_do");
  displayln_info(fname +
                 ssprintf(": job_tag='%s' ; traj=%d", job_tag.c_str(), traj));
  const std::string job_path = get_job_path(job_tag, traj);
  const Coordinate total_site = get_total_site(job_tag);
  Geometry geo;
  geo.init(total_site, 1);
  //
  FieldSelection fsel;
  set_sparse_parameters(fsel, job_tag, traj);
  //
  GaugeField gf;
  load_configuration(gf, job_tag, traj);
  gf_show_info(gf);
  qassert(is_matching_geo(gf.geo(), geo));
  //
  const FermionAction fa = get_fermion_actions(job_tag)[0];
  LowModes lm;
  // load_or_compute_low_modes(lm, get_low_modes_path(job_tag, traj), gf, fa,
  // get_lanc_arg(job_tag));
  Inverter inv;
  setup_inverter(inv, gf, fa, lm);
  //
  const Coordinate xg_point_src = Coordinate(0,0,0,0);
  //
  Propagator4d prop;
  set_point_src_propagator(prop, inv, xg_point_src);
  //
  // save prop (big endian)
  write_field_float_from_double(prop, job_path + "/psrc-prop-0.field");
  dist_write_field_float_from_double(prop, Coordinate(2,2,2,4), job_path + "/psrc-prop-0");
  write_selected_field_float_from_double(prop, job_path + "/psrc-prop-0.sfield", fsel);
  //
  // pion contraction
  const LatData ld = contract_pion(prop, xg_point_src[3]);
  if (0 == get_id_node()) {
    ld.save(job_path + "/pion-corr.lat");
  }
  qtouch_info(job_path + "/pion-corr.txt", show(ld));
  if (0 == get_id_node()) {
    print(ld);
  }
  //
  displayln_info(ssprintf("prop qnorm = %24.17E", qnorm(prop)));
  //
  displayln_info(ssprintf("prop[0] qnorm = %24.17E", qnorm(prop.get_elem(Coordinate()))));
  //
  qtouch_info(job_path + "/checkpoint.txt");
  return false;
}

inline bool compute_traj(const std::string& job_tag, const int traj)
{
  TIMER_VERBOSE("compute_traj");
  qmkdir_info(get_job_path(job_tag));
  qmkdir_info(get_job_path(job_tag) + "/logs");
  switch_monitor_file_info(get_job_path(job_tag) +
                           ssprintf("/logs/%010ld.txt", get_log_idx()));
  setup(job_tag);
  displayln_info(fname + ssprintf(": Checking '%s'.",
                                  get_job_path(job_tag, traj).c_str()));
  if (does_file_exist_sync_node(get_job_path(job_tag, traj) +
                                "/checkpoint.txt")) {
    displayln_info(fname + ssprintf(": Finished '%s'.",
                                    get_job_path(job_tag, traj).c_str()));
    return false;
  }
  if (get_config_fn(job_tag, traj) == "") {
    displayln_info(fname + ssprintf(": No config '%s'.",
                                    get_job_path(job_tag, traj).c_str()));
    return false;
  }
  // if (does_file_exist_sync_node(get_low_modes_path(job_tag, traj) +
  //                               "/checkpoint")) {
  //   displayln_info(fname + ssprintf(": No low modes '%s'.",
  //                                   get_job_path(job_tag, traj).c_str()));
  //   return false;
  // }
  if (obtain_lock(get_job_path(job_tag, traj) + "-lock")) {
    displayln_info(fname + ssprintf(": Start computing '%s'.",
                                    get_job_path(job_tag, traj).c_str()));
    const std::string job_path = get_job_path(job_tag, traj);
    qmkdir_info(job_path);
    qmkdir_sync_node(job_path + "/logs");
    switch_monitor_file_info(job_path +
                             ssprintf("/logs/%010ld.txt", get_log_idx()));
    const bool is_failed = compute_traj_do(job_tag, traj);
    Timer::display();
    switch_monitor_file_info(get_job_path(job_tag) +
                             ssprintf("/logs/%010ld.txt", get_log_idx()));
    release_lock();
    return is_failed;
  } else {
    displayln_info(fname + ssprintf(": Cannot obtain lock '%s'.",
                                    get_job_path(job_tag, traj).c_str()));
    return true;
  }
}

inline bool compute(const std::string& job_tag)
{
  TIMER_VERBOSE("compute");
  bool is_failed = false;
  const std::vector<int> trajs = get_trajs(job_tag);
  for (int i = 0; i < (int)trajs.size(); ++i) {
    const int traj = trajs[i];
    is_failed = compute_traj(job_tag, traj) or is_failed;
  }
  return is_failed;
}

// -----------------------------------------------------------------------------------

}  // namespace qlat

int main(int argc, char* argv[])
{
  using namespace qlat;
  std::vector<std::string> job_tags;
  job_tags.push_back("test-4nt8");
  job_tags.push_back("free-4nt8");
  std::vector<Coordinate> size_node_list;
  size_node_list.push_back(Coordinate(2,2,2,4));
  size_node_list.push_back(Coordinate(1,2,2,2));
  size_node_list.push_back(Coordinate(1,1,1,4));
  begin(&argc, &argv, size_node_list);
  display_geometry_node();
  setup_log_idx();
  setup();
  for (int k = 0; k < (int)job_tags.size(); ++k) {
    const std::string& job_tag = job_tags[k];
    if (not compute(job_tag)) {
      Timer::display();
      displayln_info(ssprintf("program finished successfully for '%s'.", job_tag.c_str()));
    }
  }
  displayln_info("CHECK: finished successfully.");
  Timer::display();
  end();
  return 0;
}
