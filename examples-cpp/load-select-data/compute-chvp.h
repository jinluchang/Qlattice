#pragma once

#include "compute-utils.h"
#include "data-load.h"

namespace qlat
{  //

inline std::string get_chvp_path(const std::string& job_tag, const int traj)
{
  return ssprintf("analysis/field-chvp/%s/results=%d", job_tag.c_str(), traj);
}

inline void compute_chvp_type(const std::string& job_tag, const int traj,
                              const std::vector<int>& type1_list,
                              const std::vector<int>& type2_list)
{
  check_sigterm();
  check_time_limit();
  Timer::autodisplay();
  const int num_type = type1_list.size();
  qassert(num_type == (int)type2_list.size());
  const std::string path = get_chvp_path(job_tag, traj);
  std::vector<std::string> fn_list(num_type);
  for (int i = 0; i < num_type; ++i) {
    const int type1 = type1_list[i];
    const int type2 = type2_list[i];
    fn_list[i] = path + ssprintf("/chvp-%d-%d.field", type1, type2);
  }
  bool is_complete = true;
  for (int i = 0; i < num_type; ++i) {
    if (not is_d_field(fn_list[i])) {
      is_complete = false;
      break;
    }
  }
  if (is_complete) {
    return;
  }
  TIMER_VERBOSE("compute_chvp_type");
  const PointSelection& psel = get_point_selection(job_tag, traj);
  const long n_points = psel.size();
  const FieldSelection& fsel = get_field_selection(job_tag, traj);
  const int tsep = tsep_op_wall_src(job_tag);
  std::map<std::string, FieldM<Complex, 8 * 8> > cache;
  std::map<std::string, long> counts;
  for (long n = 0; n < n_points; ++n) {
    const long xg_y_psel_idx = n;
    const Coordinate& xg_y = psel[xg_y_psel_idx];
    Timer::autodisplay();
    TIMER_VERBOSE("compute_chvp_type-iter");
    const ShiftShufflePlan ssp = make_shift_shuffle_plan(fsel, -xg_y);
    for (int i = 0; i < num_type; ++i) {
      const int type1 = type1_list[i];
      const int type2 = type2_list[i];
      if (get_point_src_info(job_tag, traj, xg_y, type1).size() == 0) {
        continue;
      }
      if (get_point_src_info(job_tag, traj, xg_y, type2).size() == 0) {
        continue;
      }
      const std::string tag = ssprintf("chvp-%d-%d", type1, type2);
      displayln_info(fname + ssprintf(":n=%ld types=%d-%d", n, type1, type2));
      SelectedField<Complex> s_chvp;
      contract_chvp_ama(s_chvp, job_tag, traj, xg_y, type1, type2);
      acc_field(cache[tag], 1.0, s_chvp, ssp);
      counts[tag] += 1;
    }
  }
  for (int i = 0; i < num_type; ++i) {
    const int type1 = type1_list[i];
    const int type2 = type2_list[i];
    const std::string tag = ssprintf("chvp-%d-%d", type1, type2);
    FieldM<Complex, 8 * 8>& chvp = cache[tag];
    const double coef = 1.0 / (double)counts[tag];
    chvp *= coef;
    // avg decay and fission to decay
    FieldM<Complex, 8 * 8> tmp;
    tmp = chvp;
    reflect_field(tmp);
    field_permute_mu_nu(tmp);
    field_conjugate_mu_nu(tmp);
    field_complex_conjugate(tmp);
    chvp += tmp;
    chvp *= 0.5;
  }
  for (int i = 0; i < num_type; ++i) {
    const int type1 = type1_list[i];
    const int type2 = type2_list[i];
    const std::string tag = ssprintf("chvp-%d-%d", type1, type2);
    const FieldM<Complex, 8 * 8>& chvp = cache[tag];
    qassert(fn_list[i] != "");
    write_field_float_from_double(chvp, fn_list[i]);
  }
}

inline void compute_chvp(const std::string& job_tag, const int traj)
{
  check_sigterm();
  Timer::autodisplay();
  TIMER_VERBOSE("compute_chvp");
  const std::string path = get_chvp_path(job_tag, traj);
  const std::string path_checkpoint = path + "/checkpoint.txt";
  if (does_file_exist_sync_node(path_checkpoint)) {
    return;
  }
  if (not(check_prop_psrc(job_tag, traj, 0) and
          check_prop_psrc(job_tag, traj, 1))) {
    return;
  }
  check_sigterm();
  check_time_limit();
  if (not obtain_lock(ssprintf("lock-chvp-%s-%d", job_tag.c_str(), traj))) {
    return;
  }
  setup(job_tag, traj);
  qmkdir_info("analysis/field-chvp");
  qmkdir_info(ssprintf("analysis/field-chvp/%s", job_tag.c_str()));
  qmkdir_info(path);
  std::vector<int> type1_list, type2_list;
  type1_list.push_back(0);
  type2_list.push_back(0);
  type1_list.push_back(0);
  type2_list.push_back(1);
  type1_list.push_back(1);
  type2_list.push_back(1);
  qassert(type1_list.size() == type2_list.size());
  compute_chvp_type(job_tag, traj, type1_list, type2_list);
  qtouch_info(path_checkpoint);
  release_lock();
}

inline void compute_chvp_light(const std::string& job_tag, const int traj)
{
  check_sigterm();
  Timer::autodisplay();
  TIMER_VERBOSE("compute_chvp_light");
  const std::string path = get_chvp_path(job_tag, traj);
  const std::string path_checkpoint = path + "/checkpoint.txt";
  if (does_file_exist_sync_node(path_checkpoint)) {
    return;
  }
  if (does_file_exist_sync_node(path + "/chvp-0-0.field")) {
    return;
  }
  if (not check_prop_psrc(job_tag, traj, 0)) {
    return;
  }
  check_sigterm();
  check_time_limit();
  if (not obtain_lock(ssprintf("lock-chvp-%s-%d", job_tag.c_str(), traj))) {
    return;
  }
  setup(job_tag, traj);
  qmkdir_info("analysis/field-chvp");
  qmkdir_info(ssprintf("analysis/field-chvp/%s", job_tag.c_str()));
  qmkdir_info(path);
  std::vector<int> type1_list, type2_list;
  type1_list.push_back(0);
  type2_list.push_back(0);
  qassert(type1_list.size() == type2_list.size());
  compute_chvp_type(job_tag, traj, type1_list, type2_list);
  release_lock();
}

}  // namespace qlat
