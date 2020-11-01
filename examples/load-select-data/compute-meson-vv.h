#pragma once

#include "compute-utils.h"
#include "data-load.h"

namespace qlat
{  //

inline std::string get_meson_vv_path(const std::string& job_tag, const int traj)
{
  return ssprintf("analysis/field-meson-vv/%s/results=%d", job_tag.c_str(),
                  traj);
}

inline void compute_meson_vv_type(const std::string& job_tag, const int traj,
                                  const std::vector<int>& type1_list,
                                  const std::vector<int>& type2_list,
                                  const int type3)
{
  check_sigint();
  check_time_limit();
  Timer::autodisplay();
  const int num_type = type1_list.size();
  qassert(num_type == (int)type2_list.size());
  const std::string path = get_meson_vv_path(job_tag, traj);
  std::vector<std::string> fn_decay_list(num_type);
  std::vector<std::string> fn_fission_list(num_type);
  for (int i = 0; i < num_type; ++i) {
    const int type1 = type1_list[i];
    const int type2 = type2_list[i];
    fn_decay_list[i] =
        path + ssprintf("/decay-%d-%d-%d.field", type1, type2, type3);
    fn_fission_list[i] =
        path + ssprintf("/fission-%d-%d-%d.field", type1, type2, type3);
  }
  bool is_complete = true;
  for (int i = 0; i < num_type; ++i) {
    if (not is_d_field(fn_decay_list[i])) {
      is_complete = false;
      break;
    }
    if (not is_d_field(fn_fission_list[i])) {
      is_complete = false;
      break;
    }
  }
  if (is_complete) {
    return;
  }
  TIMER_VERBOSE("compute_meson_vv_type");
  const PointSelection& psel = get_point_selection(job_tag, traj);
  const long n_points = psel.size();
  const FieldSelection& fsel = get_field_selection(job_tag, traj);
  const Geometry& geo = fsel.f_rank.geo;
  const int tsep = tsep_op_wall_src(job_tag);
  std::vector<FieldM<Complex, 8 * 8> > meson_vv_decay_list(num_type);
  std::vector<FieldM<Complex, 8 * 8> > meson_vv_fission_list(num_type);
  for (int i = 0; i < num_type; ++i) {
    meson_vv_decay_list[i].init(geo);
    meson_vv_fission_list[i].init(geo);
    set_zero(meson_vv_decay_list[i]);
    set_zero(meson_vv_fission_list[i]);
  }
  long iter = 0;
  for (long n = 0; n < n_points; ++n) {
    const long xg_y_psel_idx = n;
    const Coordinate& xg_y = psel[xg_y_psel_idx];
    if (get_point_src_info(job_tag, traj, xg_y, type3).size() == 0) {
      continue;
    }
    TIMER_VERBOSE("compute_meson_vv_type-iter");
    iter += 1;
    const SelProp& prop3_x_y = get_prop_psrc_ama(job_tag, traj, xg_y, type3);
    const ShiftShufflePlan ssp = make_shift_shuffle_plan(fsel, -xg_y);
    const ShiftShufflePlan ssp_reflect =
        make_shift_shuffle_plan(fsel, -xg_y, true);
    for (int i = 0; i < num_type; ++i) {
      const int type1 = type1_list[i];
      const int type2 = type2_list[i];
      displayln_info(fname + ssprintf(":n=%ld iter=%ld types=%d-%d-%d", n, iter,
                                      type1, type2, type3));
      const WallSrcProps& wsp1 = get_wall_src_props(job_tag, traj, type1);
      const WallSrcProps& wsp2 = get_wall_src_props(job_tag, traj, type2);
      contract_meson_vv_acc(meson_vv_decay_list[i], meson_vv_fission_list[i],
                            wsp1, wsp2, prop3_x_y, xg_y, xg_y_psel_idx, tsep,
                            psel, fsel, ssp, ssp_reflect);
    }
  }
  const long n_iter = iter;
  for (int i = 0; i < num_type; ++i) {
    const double coef = 1.0 / (double)n_iter;
    meson_vv_decay_list[i] *= coef;
    meson_vv_fission_list[i] *= coef;
  }
  for (int i = 0; i < num_type; ++i) {
    write_field_float_from_double(meson_vv_decay_list[i], fn_decay_list[i]);
    write_field_float_from_double(meson_vv_fission_list[i], fn_fission_list[i]);
  }
}

inline void compute_meson_vv(const std::string& job_tag, const int traj)
{
  check_sigint();
  Timer::autodisplay();
  TIMER_VERBOSE("compute_meson_vv");
  const std::string path = get_meson_vv_path(job_tag, traj);
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
  check_sigint();
  check_time_limit();
  if (not obtain_lock(ssprintf("lock-meson-vv-%s-%d", job_tag.c_str(), traj))) {
    return;
  }
  setup(job_tag, traj);
  qmkdir_info("analysis/field-meson-vv");
  qmkdir_info(ssprintf("analysis/field-meson-vv/%s", job_tag.c_str()));
  qmkdir_info(path);
  std::vector<int> type1_list, type2_list;
  type1_list.push_back(0);
  type2_list.push_back(0);
  type1_list.push_back(0);
  type2_list.push_back(1);
  type1_list.push_back(1);
  type2_list.push_back(1);
  qassert(type1_list.size() == type2_list.size());
  compute_meson_vv_type(job_tag, traj, type1_list, type2_list, 0);
  compute_meson_vv_type(job_tag, traj, type1_list, type2_list, 1);
  qtouch_info(path_checkpoint);
  release_lock();
}

}  // namespace qlat
