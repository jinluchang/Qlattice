#pragma once

#include "data-load.h"

namespace qlat
{  //

inline std::string get_psel_fsel_distribution_path(const std::string& job_tag,
                                                   const int traj)
{
  return ssprintf("analysis/field-psel-fsel-distribution/%s/results=%d",
                  job_tag.c_str(), traj);
}

inline void compute_psel_fsel_distribution_type(const std::string& job_tag,
                                                const int traj, const int type)
{
  check_sigterm();
  check_time_limit();
  Timer::autodisplay();
  const std::string path = get_psel_fsel_distribution_path(job_tag, traj);
  const std::string fn_pos = path + ssprintf("/pos-%d.field", type);
  const std::string fn_neg = path + ssprintf("/neg-%d.field", type);
  const std::string fn_avg = path + ssprintf("/avg-%d.field", type);
  bool is_complete = true;
  if (not is_d_field(fn_pos)) {
    is_complete = false;
  }
  if (not is_d_field(fn_neg)) {
    is_complete = false;
  }
  if (not is_d_field(fn_avg)) {
    is_complete = false;
  }
  if (is_complete) {
    return;
  }
  TIMER_VERBOSE("compute_psel_fsel_distribution_type");
  const PointSelection& psel = get_point_selection(job_tag, traj);
  const long n_points = psel.size();
  const FieldSelection& fsel = get_field_selection(job_tag, traj);
  FieldM<Complex, 1> pos;
  long iter = 0;
  for (long n = 0; n < n_points; ++n) {
    const long xg_y_psel_idx = n;
    const Coordinate& xg_y = psel[xg_y_psel_idx];
    if (get_point_src_info(job_tag, traj, xg_y, type).size() == 0) {
      continue;
    }
    Timer::autodisplay();
    TIMER_VERBOSE("compute_psel_fsel_distribution_type-iter");
    iter += 1;
    const ShiftShufflePlan ssp = make_shift_shuffle_plan(fsel, -xg_y);
    displayln_info(fname + ssprintf(":n=%ld iter=%ld type=%d", n, iter, type));
    contract_psel_fsel_distribution_acc(pos, xg_y, fsel, ssp);
  }
  const long n_iter = iter;
  const double coef = 1.0 / (double)n_iter;
  pos *= coef;
  FieldM<Complex, 1> neg;
  neg = pos;
  reflect_field(neg);
  FieldM<Complex, 1> avg;
  avg += pos;
  avg += neg;
  avg *= 0.5;
  write_field_double(pos, fn_pos);
  write_field_double(neg, fn_neg);
  write_field_double(avg, fn_avg);
}

inline void compute_psel_fsel_distribution(const std::string& job_tag, const int traj)
{
  check_sigterm();
  Timer::autodisplay();
  TIMER_VERBOSE("compute_psel_fsel_distribution");
  const std::string path = get_psel_fsel_distribution_path(job_tag, traj);
  const std::string path_checkpoint = path + "/checkpoint.txt";
  if (does_file_exist_sync_node(path_checkpoint)) {
    return;
  }
  if (not(check_sparse_parameters(job_tag, traj))) {
    return;
  }
  check_sigterm();
  check_time_limit();
  if (not obtain_lock(ssprintf("lock-psel-fsel-distribution-%s-%d", job_tag.c_str(), traj))) {
    return;
  }
  setup(job_tag, traj);
  qmkdir_info("analysis/field-psel-fsel-distribution");
  qmkdir_info(ssprintf("analysis/field-psel-fsel-distribution/%s", job_tag.c_str()));
  qmkdir_info(path);
  compute_psel_fsel_distribution_type(job_tag, traj, 0);
  compute_psel_fsel_distribution_type(job_tag, traj, 1);
  qtouch_info(path_checkpoint);
  release_lock();
}

}  // namespace qlat
