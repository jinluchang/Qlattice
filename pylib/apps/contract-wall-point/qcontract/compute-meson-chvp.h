#pragma once

#include "compute-meson-snk-src.h"
#include "compute-utils.h"
#include "data-load.h"

namespace qlat
{  //

inline std::string get_meson_chvp_path(const std::string& job_tag,
                                       const int traj)
{
  return ssprintf("analysis/field-meson-chvp/%s/results=%d", job_tag.c_str(),
                  traj);
}

inline std::string get_meson_snk_src_shift_weight_path(
    const std::string& job_tag, const int traj)
{
  return ssprintf("analysis/lat-meson-snk-src-shift-weight/%s/results=%d",
                  job_tag.c_str(), traj);
}

inline void compute_meson_chvp_type(const std::string& job_tag, const int traj,
                                    const std::vector<int>& type1_list,
                                    const std::vector<int>& type2_list,
                                    const std::vector<int>& type3_list,
                                    const std::vector<int>& type4_list)
{
  check_sigterm();
  check_time_limit();
  Timer::autodisplay();
  const int num_type_12 = type1_list.size();
  qassert(num_type_12 == (int)type2_list.size());
  const int num_type_34 = type3_list.size();
  qassert(num_type_34 == (int)type4_list.size());
  const int num_type = num_type_12 * num_type_34;
  const std::string path = get_meson_chvp_path(job_tag, traj);
  const std::string path_mss_s_w =
      get_meson_snk_src_shift_weight_path(job_tag, traj);
  const std::string path_mss = get_meson_snk_src(job_tag, traj);
  std::vector<std::string> fn_mss_list(num_type_12);
  std::vector<std::string> fn_mss_s_w_list(num_type);
  std::vector<std::string> fn_list(num_type);
  for (int i = 0; i < num_type_12; ++i) {
    const int type1 = type1_list[i];
    const int type2 = type2_list[i];
    fn_mss_list[i] =
        path_mss + ssprintf("/meson-snk-src-%d-%d.lat", type1, type2);
    if (type1 <= type2) {
      for (int j = 0; j < num_type_34; ++j) {
        const int type3 = type3_list[j];
        const int type4 = type4_list[j];
        std::string& fn = fn_list[i * num_type_34 + j];
        fn = path +
             ssprintf("/mchvp-%d-%d-%d-%d.field", type1, type2, type3, type4);
        std::string& fn_mss_s_w = fn_mss_s_w_list[i * num_type_34 + j];
        fn_mss_s_w = path_mss_s_w + ssprintf("/meson-snk-src-%d-%d-%d-%d.field",
                                             type1, type2, type3, type4);
      }
    }
  }
  bool is_complete = true;
  for (int i = 0; i < num_type; ++i) {
    if (fn_list[i] != "") {
      if (not is_d_field(fn_list[i])) {
        is_complete = false;
        break;
      }
      if (not does_file_exist_sync_node(fn_mss_s_w_list[i])) {
        is_complete = false;
        break;
      }
    }
  }
  if (is_complete) {
    return;
  }
  for (int i = 0; i < num_type_12; ++i) {
    if (not does_file_exist_sync_node(fn_mss_list[i])) {
      displayln_info(ssprintf("Missing '%s'", fn_mss_list[i].c_str()));
      return;
    }
  }
  TIMER_VERBOSE("compute_meson_chvp_type");
  const PointSelection& psel = get_point_selection(job_tag, traj);
  const long n_points = psel.size();
  const FieldSelection& fsel = get_field_selection(job_tag, traj);
  const Geometry& geo = fsel.f_rank.geo();
  const Coordinate total_site = geo.total_site();
  const int tsep = tsep_op_wall_src(job_tag);
  std::map<std::string, FieldM<Complex, 8 * 8> > cache;
  std::map<std::string, long> counts;
  std::map<std::string, LatData> ld_mss_shift_weight;
  std::vector<LatData> ld_mss_list(num_type_12);
  for (int i = 0; i < num_type_12; ++i) {
    ld_mss_list[i] = lat_data_load_info(fn_mss_list[i]);
  }
  for (long tslice = 0; tslice < total_site[3]; ++tslice) {
    Timer::autodisplay();
    TIMER_VERBOSE("compute_meson_chvp_type-tslice");
    std::vector<FieldM<Complex, 8 * 8> > chvp_list(num_type_34);
    std::vector<long> chvp_count_list(num_type_34, 0);
    for (long n = 0; n < n_points; ++n) {
      const long xg_y_psel_idx = n;
      const Coordinate& xg_y = psel[xg_y_psel_idx];
      if (xg_y[3] != tslice) {
        continue;
      }
      const ShiftShufflePlan ssp = make_shift_shuffle_plan(fsel, -xg_y);
      for (int j = 0; j < num_type_34; ++j) {
        const int type3 = type3_list[j];
        const int type4 = type4_list[j];
        if (get_point_src_info(job_tag, traj, xg_y, type3).size() == 0) {
          continue;
        }
        if (get_point_src_info(job_tag, traj, xg_y, type4).size() == 0) {
          continue;
        }
        displayln_info(fname +
                       ssprintf(": tslice=%ld n=%ld j=%d chvp types=%d-%d",
                                tslice, n, j, type3, type4));
        SelectedField<Complex> s_chvp;
        contract_chvp_ama(s_chvp, job_tag, traj, xg_y, type3, type4);
        acc_field(chvp_list[j], 1.0, s_chvp, ssp);
        chvp_count_list[j] += 1;
      }
    }
    for (int i = 0; i < num_type_12; ++i) {
      const int type1 = type1_list[i];
      const int type2 = type2_list[i];
      const LatData ld_mss_shifted =
          meson_snk_src_shift(ld_mss_list[i], -tslice);
      for (int j = 0; j < num_type_34; ++j) {
        if (chvp_count_list[j] == 0) {
          continue;
        }
        qassert(is_initialized(chvp_list[j]));
        const int type3 = type3_list[j];
        const int type4 = type4_list[j];
        const std::string tag =
            ssprintf("mchvp-%d-%d-%d-%d", type1, type2, type3, type4);
        contract_meson_chvp_acc(cache[tag], ld_mss_shifted, chvp_list[j], tsep);
        counts[tag] += chvp_count_list[j];
        ld_mss_shift_weight[tag] += chvp_count_list[j] * ld_mss_shifted;
      }
    }
  }
  for (int i = 0; i < num_type_12; ++i) {
    const int type1 = type1_list[i];
    const int type2 = type2_list[i];
    for (int j = 0; j < num_type_34; ++j) {
      const int type3 = type3_list[j];
      const int type4 = type4_list[j];
      const std::string tag =
          ssprintf("mchvp-%d-%d-%d-%d", type1, type2, type3, type4);
      const double coef = 1.0 / (double)counts[tag];
      cache[tag] *= coef;
      ld_mss_shift_weight[tag] *= coef;
    }
  }
  for (int i = 0; i < num_type_12; ++i) {
    const int type1 = type1_list[i];
    const int type2 = type2_list[i];
    if (type1 <= type2) {
      for (int j = 0; j < num_type_34; ++j) {
        const int type3 = type3_list[j];
        const int type4 = type4_list[j];
        const std::string tag12 =
            ssprintf("mchvp-%d-%d-%d-%d", type1, type2, type3, type4);
        const std::string tag21 =
            ssprintf("mchvp-%d-%d-%d-%d", type2, type1, type3, type4);
        const std::string& fn_mss_s_w = fn_mss_s_w_list[i * num_type_34 + j];
        const LatData& ld_mss_s_w = ld_mss_shift_weight[tag12];
        lat_data_save_info(fn_mss_s_w, ld_mss_s_w);
        FieldM<Complex, 8 * 8>& mchvp_1_2_3_4 = cache[tag12];
        FieldM<Complex, 8 * 8>& mchvp_2_1_3_4 = cache[tag21];
        FieldM<Complex, 8 * 8> tmp;
        tmp = mchvp_2_1_3_4;
        reflect_field(tmp);
        tmp += mchvp_1_2_3_4;
        tmp *= 0.5;
        mchvp_1_2_3_4 = tmp;
        field_permute_mu_nu(tmp);
        field_conjugate_mu_nu(tmp);
        field_complex_conjugate(tmp);
        mchvp_1_2_3_4 += tmp;
        mchvp_1_2_3_4 *= 0.5;
        const std::string& fn = fn_list[i * num_type_34 + j];
        qassert(fn != "");
        write_field_float_from_double(mchvp_1_2_3_4, fn);
      }
    }
  }
}

inline void compute_meson_chvp(const std::string& job_tag, const int traj)
{
  check_sigterm();
  Timer::autodisplay();
  TIMER_VERBOSE("compute_meson_chvp");
  const std::string path = get_meson_chvp_path(job_tag, traj);
  const std::string path_mss =
      get_meson_snk_src_shift_weight_path(job_tag, traj);
  const std::string path_checkpoint = path + "/checkpoint.txt";
  const std::string path_mss_checkpoint = path_mss + "/checkpoint.txt";
  if (does_file_exist_sync_node(path_checkpoint) and
      does_file_exist_sync_node(path_mss_checkpoint)) {
    return;
  }
  if (not does_file_exist_sync_node(get_meson_snk_src(job_tag, traj) +
                                    "/checkpoint.txt")) {
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
          ssprintf("lock-meson-chvp-%s-%d", job_tag.c_str(), traj))) {
    return;
  }
  setup(job_tag, traj);
  qmkdir_info("analysis/field-meson-chvp");
  qmkdir_info(ssprintf("analysis/field-meson-chvp/%s", job_tag.c_str()));
  qmkdir_info(path);
  qmkdir_info("analysis/lat-meson-snk-src-shift-weight");
  qmkdir_info(
      ssprintf("analysis/lat-meson-snk-src-shift-weight/%s", job_tag.c_str()));
  qmkdir_info(path_mss);
  std::vector<int> type1_list, type2_list;
  type1_list.push_back(0);
  type2_list.push_back(0);
  type1_list.push_back(0);
  type2_list.push_back(1);
  type1_list.push_back(1);
  type2_list.push_back(0);
  std::vector<int> type3_list, type4_list;
  type3_list.push_back(0);
  type4_list.push_back(0);
  type3_list.push_back(0);
  type4_list.push_back(1);
  type3_list.push_back(1);
  type4_list.push_back(1);
  qassert(type1_list.size() == type2_list.size());
  qassert(type3_list.size() == type4_list.size());
  compute_meson_chvp_type(job_tag, traj, type1_list, type2_list, type3_list,
                          type4_list);
  qtouch_info(path_checkpoint);
  qtouch_info(path_mss_checkpoint);
  release_lock();
}

inline void compute_meson_chvp_light(const std::string& job_tag, const int traj)
{
  check_sigterm();
  Timer::autodisplay();
  TIMER_VERBOSE("compute_meson_chvp_light");
  const std::string path = get_meson_chvp_path(job_tag, traj);
  const std::string path_mss =
      get_meson_snk_src_shift_weight_path(job_tag, traj);
  const std::string path_checkpoint = path + "/checkpoint.txt";
  const std::string path_mss_checkpoint = path_mss + "/checkpoint.txt";
  if (does_file_exist_sync_node(path_checkpoint) and
      does_file_exist_sync_node(path_mss_checkpoint)) {
    return;
  }
  if (does_file_exist_sync_node(path + "/mchvp-0-0-0-0.field") and
      does_file_exist_sync_node(path + "/mchvp-0-0-0-1.field") and
      does_file_exist_sync_node(path + "/mchvp-0-0-1-1.field") and
      does_file_exist_sync_node(path_mss + "/meson-snk-src-0-0-0-0.field") and
      does_file_exist_sync_node(path_mss + "/meson-snk-src-0-0-0-1.field") and
      does_file_exist_sync_node(path_mss + "/meson-snk-src-0-0-1-1.field")) {
    return;
  }
  if (not does_file_exist_sync_node(get_meson_snk_src(job_tag, traj) +
                                    "/meson-snk-src-0-0.lat")) {
    return;
  }
  if (not(check_wall_src_props(job_tag, traj, 0) and
          check_prop_psrc(job_tag, traj, 0) and
          check_prop_psrc(job_tag, traj, 1))) {
    return;
  }
  check_sigterm();
  check_time_limit();
  if (not obtain_lock(
          ssprintf("lock-meson-chvp-%s-%d", job_tag.c_str(), traj))) {
    return;
  }
  setup(job_tag, traj);
  qmkdir_info("analysis/field-meson-chvp");
  qmkdir_info(ssprintf("analysis/field-meson-chvp/%s", job_tag.c_str()));
  qmkdir_info(path);
  qmkdir_info("analysis/lat-meson-snk-src-shift-weight");
  qmkdir_info(
      ssprintf("analysis/lat-meson-snk-src-shift-weight/%s", job_tag.c_str()));
  qmkdir_info(path_mss);
  std::vector<int> type1_list, type2_list;
  type1_list.push_back(0);
  type2_list.push_back(0);
  std::vector<int> type3_list, type4_list;
  type3_list.push_back(0);
  type4_list.push_back(0);
  type3_list.push_back(0);
  type4_list.push_back(1);
  type3_list.push_back(1);
  type4_list.push_back(1);
  qassert(type1_list.size() == type2_list.size());
  qassert(type3_list.size() == type4_list.size());
  compute_meson_chvp_type(job_tag, traj, type1_list, type2_list, type3_list,
                          type4_list);
  release_lock();
}

}  // namespace qlat
