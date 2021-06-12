#pragma once

#include "compute-utils.h"
#include "data-load.h"

namespace qlat
{  //

inline std::string get_meson_vv_meson_path(const std::string& job_tag, const int traj)
{
  return ssprintf("analysis/field-meson-vv-meson/%s/results=%d", job_tag.c_str(),
                  traj);
}

inline void compute_meson_vv_meson_type(const std::string& job_tag,
                                        const int traj,
                                        const std::vector<int>& type1_list,
                                        const std::vector<int>& type2_list,
                                        const std::vector<int>& type3_list,
                                        const int type4)
{
  check_sigterm();
  check_time_limit();
  Timer::autodisplay();
  const int num_type = type1_list.size();
  qassert(num_type == (int)type2_list.size());
  qassert(num_type == (int)type3_list.size());
  const std::string path = get_meson_vv_meson_path(job_tag, traj);
  std::vector<std::string> fn_forward_list(num_type);
  for (int i = 0; i < num_type; ++i) {
    const int type1 = type1_list[i];
    const int type2 = type2_list[i];
    const int type3 = type3_list[i];
    if (type1 <= type2) {
      fn_forward_list[i] = path + ssprintf("/forward-%d-%d-%d-%d.field", type1,
                                           type2, type3, type4);
    }
  }
  bool is_complete = true;
  for (int i = 0; i < num_type; ++i) {
    if (fn_forward_list[i] != "") {
      if (not is_d_field(fn_forward_list[i])) {
        is_complete = false;
        break;
      }
    }
  }
  if (is_complete) {
    return;
  }
  TIMER_VERBOSE("compute_meson_vv_meson_type");
  const PointSelection& psel = get_point_selection(job_tag, traj);
  const long n_points = psel.size();
  const FieldSelection& fsel = get_field_selection(job_tag, traj);
  const int tsep = tsep_op_wall_src(job_tag);
  std::map<std::string, FieldM<Complex, 8 * 8> > cache;
  long iter = 0;
  for (long n = 0; n < n_points; ++n) {
    const long xg_y_psel_idx = n;
    const Coordinate& xg_y = psel[xg_y_psel_idx];
    if (get_point_src_info(job_tag, traj, xg_y, type4).size() == 0) {
      continue;
    }
    Timer::autodisplay();
    TIMER_VERBOSE("compute_meson_vv_meson_type-iter");
    iter += 1;
    const SelProp& prop4_x_y = get_prop_psrc_ama(job_tag, traj, xg_y, type4);
    const ShiftShufflePlan ssp = make_shift_shuffle_plan(fsel, -xg_y);
    for (int i = 0; i < num_type; ++i) {
      const int type1 = type1_list[i];
      const int type2 = type2_list[i];
      const int type3 = type3_list[i];
      FieldM<Complex, 8 * 8>& forward =
          cache[ssprintf("forward-%d-%d-%d-%d", type1, type2, type3, type4)];
      FieldM<Complex, 8 * 8>& backward =
          cache[ssprintf("backward-%d-%d-%d-%d", type1, type2, type3, type4)];
      displayln_info(fname + ssprintf(":n=%ld iter=%ld types=%d-%d-%d-%d", n,
                                      iter, type1, type2, type3, type4));
      const WallSrcProps& wsp1 = get_wall_src_props(job_tag, traj, type1);
      const WallSrcProps& wsp2 = get_wall_src_props(job_tag, traj, type2);
      const WallSrcProps& wsp3 = get_wall_src_props(job_tag, traj, type3);
      contract_meson_vv_meson_acc(forward, backward, wsp1, wsp2, wsp3,
                                  prop4_x_y, xg_y, xg_y_psel_idx, tsep, psel,
                                  fsel, ssp);
    }
  }
  const long n_iter = iter;
  for (int i = 0; i < num_type; ++i) {
    const int type1 = type1_list[i];
    const int type2 = type2_list[i];
    const int type3 = type3_list[i];
    FieldM<Complex, 8 * 8>& forward =
        cache[ssprintf("forward-%d-%d-%d-%d", type1, type2, type3, type4)];
    FieldM<Complex, 8 * 8>& backward =
        cache[ssprintf("backward-%d-%d-%d-%d", type1, type2, type3, type4)];
    const double coef = 1.0 / (double)n_iter;
    forward *= coef;
    backward *= coef;
    // avg forward and backward to forward
    reflect_field(backward);
    forward += backward;
    forward *= 0.5;
  }
  for (int i = 0; i < num_type; ++i) {
    const int type1 = type1_list[i];
    const int type2 = type2_list[i];
    const int type3 = type3_list[i];
    if (type1 <= type2) {
      FieldM<Complex, 8 * 8> avg;
      FieldM<Complex, 8 * 8>& f1 =
          cache[ssprintf("forward-%d-%d-%d-%d", type1, type2, type3, type4)];
      FieldM<Complex, 8 * 8>& f2 =
          cache[ssprintf("forward-%d-%d-%d-%d", type2, type1, type3, type4)];
      qassert(is_initialized(f1));
      qassert(is_initialized(f2));
      avg = f2;
      field_permute_mu_nu(avg);
      field_conjugate_mu_nu(avg);
      field_complex_conjugate(avg);
      avg += f1;
      avg *= 0.5;
      qassert(fn_forward_list[i] != "");
      write_field_float_from_double(avg, fn_forward_list[i]);
    }
  }
}

inline void compute_meson_vv_meson(const std::string& job_tag, const int traj)
{
  check_sigterm();
  Timer::autodisplay();
  TIMER_VERBOSE("compute_meson_vv_meson");
  const std::string path = get_meson_vv_meson_path(job_tag, traj);
  const std::string path_checkpoint = path + "/checkpoint.txt";
  if (does_file_exist_sync_node(path_checkpoint)) {
    return;
  }
  if (not(check_wall_src_props(job_tag, traj, 0) and
          check_wall_src_props(job_tag, traj, 1) and
          check_prop_psrc(job_tag, traj, 0) and
          check_prop_psrc(job_tag, traj, 1))) {
    return;
  }
  check_sigterm();
  check_time_limit();
  if (not obtain_lock(
          ssprintf("lock-meson-vv-meson-%s-%d", job_tag.c_str(), traj))) {
    return;
  }
  setup(job_tag, traj);
  qmkdir_info("analysis/field-meson-vv-meson");
  qmkdir_info(ssprintf("analysis/field-meson-vv-meson/%s", job_tag.c_str()));
  qmkdir_info(path);
  std::vector<int> type1_list, type2_list, type3_list;
  type1_list.push_back(0);
  type2_list.push_back(0);
  type3_list.push_back(0);
  type1_list.push_back(0);
  type2_list.push_back(0);
  type3_list.push_back(1);
  type1_list.push_back(1);
  type2_list.push_back(1);
  type3_list.push_back(0);
  type1_list.push_back(1);
  type2_list.push_back(0);
  type3_list.push_back(0);
  type1_list.push_back(0);
  type2_list.push_back(1);
  type3_list.push_back(0);
  qassert(type1_list.size() == type2_list.size());
  qassert(type1_list.size() == type3_list.size());
  compute_meson_vv_meson_type(job_tag, traj, type1_list, type2_list, type3_list,
                              0);
  compute_meson_vv_meson_type(job_tag, traj, type1_list, type2_list, type3_list,
                              1);
  qtouch_info(path_checkpoint);
  release_lock();
}

inline void compute_meson_vv_meson_light(const std::string& job_tag,
                                         const int traj)
{
  check_sigterm();
  Timer::autodisplay();
  TIMER_VERBOSE("compute_meson_vv_meson_light");
  const std::string path = get_meson_vv_meson_path(job_tag, traj);
  const std::string path_checkpoint = path + "/checkpoint.txt";
  if (does_file_exist_sync_node(path_checkpoint)) {
    return;
  }
  if (does_file_exist_sync_node(path + "/forward-0-0-0-0.field")) {
    return;
  }
  if (not(check_wall_src_props(job_tag, traj, 0) and
          check_prop_psrc(job_tag, traj, 0))) {
    return;
  }
  check_sigterm();
  check_time_limit();
  if (not obtain_lock(
          ssprintf("lock-meson-vv-meson-%s-%d", job_tag.c_str(), traj))) {
    return;
  }
  setup(job_tag, traj);
  qmkdir_info("analysis/field-meson-vv-meson");
  qmkdir_info(ssprintf("analysis/field-meson-vv-meson/%s", job_tag.c_str()));
  qmkdir_info(path);
  std::vector<int> type1_list, type2_list, type3_list;
  type1_list.push_back(0);
  type2_list.push_back(0);
  type3_list.push_back(0);
  qassert(type1_list.size() == type2_list.size());
  qassert(type1_list.size() == type3_list.size());
  compute_meson_vv_meson_type(job_tag, traj, type1_list, type2_list, type3_list,
                              0);
  release_lock();
}

}  // namespace qlat
