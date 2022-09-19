#include "configs.h"
#include "paths-old.h"
#include "paths.h"
#include "psrc-distribution.h"
#include "psrc-sample.h"
#include "qlat-setup.h"

namespace qlat
{  //

inline void collect_pis(const std::string& job_tag)
{
  const std::string old_job_tag = get_old_job_tag(job_tag);
  const std::string checkpoint =
      ssprintf("data/point-src-info/%s/checkpoint.txt", job_tag.c_str());
  if (does_file_exist_sync_node(checkpoint)) {
    return;
  }
  if (not obtain_lock(ssprintf("lock-point-src-info-%s", job_tag.c_str()))) {
    return;
  }
  {
    TIMER_VERBOSE("collect_pis");
    qmkdir_info(ssprintf("data/point-src-info"));
    qmkdir_info(ssprintf("data/point-src-info/%s", job_tag.c_str()));
    const std::vector<int> trajs = get_todo_trajs(old_job_tag);
    for (int i = 0; i < (int)trajs.size(); ++i) {
      const int traj = trajs[i];
      check_sigint();
      const std::string fn = get_point_src_info_path(job_tag, traj);
      const std::string fn_old = get_pis_path(old_job_tag, traj);
      if (fn_old != "") {
        const std::vector<PointInfo> pis = load_lbl_pis_info(fn_old);
        if (not does_file_exist_sync_node(fn)) {
          save_lbl_pis_info(pis, fn);
        }
        const std::vector<PointInfo> pis_load = load_lbl_pis_info(fn);
        qassert(pis_load == pis);
      }
    }
    qtouch_info(checkpoint);
    release_lock();
  }
  Timer::display();
}

inline void collect_pcs(const std::string& job_tag)
{
  const std::string old_job_tag = get_old_job_tag(job_tag);
  const std::string checkpoint =
      ssprintf("data/point-selection/%s/checkpoint.txt", job_tag.c_str());
  if (does_file_exist_sync_node(checkpoint)) {
    return;
  }
  if (not obtain_lock(ssprintf("lock-point-selection-%s", job_tag.c_str()))) {
    return;
  }
  {
    TIMER_VERBOSE("collect_pcs");
    qmkdir_info(ssprintf("data/point-selection"));
    qmkdir_info(ssprintf("data/point-selection/%s", job_tag.c_str()));
    const std::vector<int> trajs = get_todo_trajs(old_job_tag);
    for (int i = 0; i < (int)trajs.size(); ++i) {
      const int traj = trajs[i];
      check_sigint();
      const std::string fn = get_point_selection_path(job_tag, traj);
      const std::string fn_old = get_pis_path(old_job_tag, traj);
      if (fn_old != "") {
        const std::vector<PointInfo> pis = load_lbl_pis_info(fn_old);
        const std::vector<Coordinate> pcs = coordinates_from_point_infos(pis);
        if (not does_file_exist_sync_node(fn)) {
          save_lbl_pcs_info(pcs, fn);
        }
        const std::vector<Coordinate> pcs_load = load_lbl_pcs_info(fn);
        qassert(pcs_load == pcs);
      }
    }
    qtouch_info(checkpoint);
    release_lock();
  }
  Timer::display();
}

inline void collect_simulated_pis(const std::string& job_tag)
{
  const std::string old_job_tag = get_old_job_tag(job_tag);
  if (does_file_exist_sync_node(ssprintf(
          "data/point-src-info-simulated/%s.tar.gz", job_tag.c_str()))) {
    return;
  }
  const std::string checkpoint = ssprintf(
      "data/point-src-info-simulated/%s/checkpoint.txt", job_tag.c_str());
  if (does_file_exist_sync_node(checkpoint)) {
    return;
  }
  if (not obtain_lock(
          ssprintf("lock-point-src-info-simulated-%s", job_tag.c_str()))) {
    return;
  }
  {
    TIMER_VERBOSE("collect_simulated_pis");
    qmkdir_info(ssprintf("data/point-src-info-simulated"));
    qmkdir_info(ssprintf("data/point-src-info-simulated/%s", job_tag.c_str()));
    const int num_simulated_traj = 16384;
    for (int traj = 0; traj < num_simulated_traj; ++traj) {
      check_sigint();
      const std::string fn = get_point_src_info_simulated_path(job_tag, traj);
      const std::string fn_old = get_simulated_pis_path(old_job_tag) +
                                 ssprintf("/traj=%d ; pis.txt", traj);
      qassert(fn_old != "");
      const std::vector<PointInfo> pis = load_lbl_pis_info(fn_old);
      if (not does_file_exist_sync_node(fn)) {
        save_lbl_pis_info(pis, fn);
      }
      const std::vector<PointInfo> pis_load = load_lbl_pis_info(fn);
      qassert(pis_load == pis);
    }
    qtouch_info(checkpoint);
    release_lock();
  }
  Timer::display();
}

inline void collect_point_distribution(const std::string& job_tag)
{
  const std::string old_job_tag = get_old_job_tag(job_tag);
  const std::string checkpoint =
      ssprintf("data/point-distribution/%s/checkpoint.txt", job_tag.c_str());
  if (does_file_exist_sync_node(checkpoint)) {
    return;
  }
  if (not obtain_lock(
          ssprintf("lock-point-distribution-%s", job_tag.c_str()))) {
    return;
  }
  {
    TIMER_VERBOSE("collect_simulated_pis");
    qmkdir_info(ssprintf("data/point-distribution"));
    qmkdir_info(ssprintf("data/point-distribution/%s", job_tag.c_str()));
    {
      const std::string fn = get_point_distribution_path(job_tag);
      const std::string fn_old = get_pis_distribution_path(old_job_tag);
      qassert(fn_old != "");
      const Coordinate total_site = get_total_site(job_tag);
      PointDistribution pd;
      load_pd(pd, total_site, fn_old);
      if (not does_file_exist_sync_node(fn)) {
        save_pd(pd, total_site, fn);
      }
      PointDistribution pd_load;
      load_pd(pd_load, total_site, fn);
      qassert(pd_load == pd);
    }
    qtouch_info(checkpoint);
    release_lock();
  }
  Timer::display();
}

inline void collect_field_selection(const std::string& job_tag, const int traj)
{
  check_sigint();
  check_time_limit();
  const std::string old_job_tag = get_old_job_tag(job_tag);
  const std::string fn = get_field_selection_path(job_tag, traj);
  const std::string fn_old = get_f_rank_path(old_job_tag, traj);
  if (does_file_exist_sync_node(fn)) {
    return;
  }
  if (not does_file_exist_sync_node(fn_old)) {
    return;
  }
  if (not obtain_lock(
          ssprintf("lock-field-selection-%s-%d", job_tag.c_str(), traj))) {
    return;
  }
  {
    setup(job_tag, traj);
    TIMER_VERBOSE("collect_field_selection");
    qmkdir_info(ssprintf("data/field-selection"));
    qmkdir_info(ssprintf("data/field-selection/%s", job_tag.c_str()));
    const Coordinate total_site = get_total_site(job_tag);
    const long spatial_vol = total_site[0] * total_site[1] * total_site[2];
    const long n_per_tslice = spatial_vol / 16;
    const RngState rs_sel =
        RngState("field-sel").split(old_job_tag).split(traj);
    FieldSelection fsel;
    set_field_selection(fsel, total_site, n_per_tslice, rs_sel);
    FieldSelection fsel_load;
    read_field_selection(fsel_load, fn_old, n_per_tslice);
    FieldM<int64_t, 1> f_rank;
    f_rank = fsel_load.f_rank;
    f_rank -= fsel.f_rank;
    qassert(qnorm(f_rank) == 0.0);
    write_field_selection(fsel, fn);
    FieldSelection fsel_load2;
    read_field_selection(fsel_load2, fn, n_per_tslice);
    FieldM<int64_t, 1> f_rank2;
    f_rank2 = fsel_load2.f_rank;
    f_rank2 -= fsel.f_rank;
    qassert(qnorm(f_rank2) == 0.0);
    release_lock();
  }
  Timer::autodisplay();
}

inline void collect_gauge_transform(const std::string& job_tag, const int traj)
{
  check_sigint();
  check_time_limit();
  const std::string old_job_tag = get_old_job_tag(job_tag);
  const std::string fn = get_gauge_transform_path(job_tag, traj);
  const std::string fn_old = get_old_gauge_transform_path(old_job_tag, traj);
  if (does_file_exist_sync_node(fn)) {
    return;
  }
  if (not does_file_exist_sync_node(fn_old)) {
    return;
  }
  if (not obtain_lock(
          ssprintf("lock-gauge-transform-%s-%d", job_tag.c_str(), traj))) {
    return;
  }
  {
    setup(job_tag, traj);
    TIMER_VERBOSE("collect_gauge_transform");
    qmkdir_info(ssprintf("data/gauge-transform"));
    qmkdir_info(ssprintf("data/gauge-transform/%s", job_tag.c_str()));
    GaugeTransform gt;
    const long total_bytes = read_field_double(gt, fn_old);
    qassert(total_bytes > 0);
    qassert(is_initialized(gt));
    write_field_double(gt, fn);
    release_lock();
  }
  Timer::autodisplay();
}

inline std::vector<PointInfo> get_point_src_info(const std::string& job_tag,
                                                 const int traj)
{
  TIMER_VERBOSE("get_point_src_info");
  return load_lbl_pis_info(get_point_src_info_path(job_tag, traj));
}

inline void set_sparse_parameters(std::vector<Coordinate>& psel,
                                  FieldSelection& fsel,
                                  const std::string& job_tag, const int traj)
{
  TIMER_VERBOSE("set_sparse_parameters");
  const std::string fn_point_selection =
      get_point_selection_path(job_tag, traj);
  const std::string fn_field_selection =
      get_field_selection_path(job_tag, traj);
  qassert(fn_point_selection != "");
  qassert(fn_field_selection != "");
  psel = load_lbl_pcs_info(fn_point_selection);
  const Coordinate total_site = get_total_site(job_tag);
  const long spatial_vol = total_site[0] * total_site[1] * total_site[2];
  const long n_per_tslice = spatial_vol / 16;
  read_field_selection(fsel, fn_field_selection, n_per_tslice);
}

inline bool check_prop_psrc(const std::string& job_tag, const int traj,
                            const int type)
{
  TIMER_VERBOSE("check_prop_psrc");
  const std::string old_job_tag = get_old_job_tag(job_tag);
  if (is_old_collection_ensemble(job_tag)) {
    const std::string& path = get_old_psrc_prop_path(old_job_tag, traj);
    return does_file_exist_sync_node(path + "/checkpoint.txt");
  } else if (is_final_run_ensemble(job_tag)) {
    const std::string& path = get_old_psrc_prop_path(old_job_tag, traj);
    if (type == 0) {
      return does_file_exist_sync_node(path + "/sparse_point_src_props_light");
    } else if (type == 1) {
      if (job_tag == "48I") {
        if (not does_file_exist_sync_node(
                path + "/sparse_point_src_props_strange_new2")) {
          return false;
        }
      }
      return does_file_exist_sync_node(path +
                                       "/sparse_point_src_props_strange");
    } else {
      qassert(false);
    }
  } else {
    qassert(false);
  }
  return false;
}

inline void collect_prop_psrc_light(const std::string& job_tag, const int traj)
{
  check_sigint();
  check_time_limit();
  const int type = 0;  // light quark
  const std::string old_job_tag = get_old_job_tag(job_tag);
  const std::string path = get_prop_psrc_light_path(job_tag, traj);
  const std::string psel_path = get_psel_prop_psrc_light_path(job_tag, traj);
  if (does_file_exist_sync_node(path) and
      does_file_exist_sync_node(psel_path)) {
    return;
  }
  if (not check_prop_psrc(job_tag, traj, type)) {
    return;
  }
  if (not obtain_lock(
          ssprintf("lock-prop-psrc-light-%s-%d", job_tag.c_str(), traj))) {
    return;
  }
  {
    setup(job_tag, traj);
    TIMER_VERBOSE("collect_prop_psrc_light");
    qmkdir_info(ssprintf("data/prop-psrc-light"));
    qmkdir_info(ssprintf("data/prop-psrc-light/%s", job_tag.c_str()));
    qmkdir_info(ssprintf("data/psel-prop-psrc-light"));
    qmkdir_info(ssprintf("data/psel-prop-psrc-light/%s", job_tag.c_str()));
    qmkdir_info(psel_path);
    const std::vector<PointInfo> pis = get_point_src_info(job_tag, traj);
    std::vector<Coordinate> psel;
    FieldSelection fsel;
    set_sparse_parameters(psel, fsel, job_tag, traj);
    const Coordinate new_size_node = Coordinate(1, 1, 1, 16);
    ShuffledBitSet sbs = mk_shuffled_bitset(fsel.f_rank, psel, new_size_node);
    ShuffledFieldsWriter sfw(path + ".acc", new_size_node, true);
    if (is_old_collection_ensemble(job_tag)) {
      for (long i = 0; i < (long)pis.size(); ++i) {
        const PointInfo& pi = pis[i];
        if (pi.type == type) {
          const Coordinate xg = pi.xg;
          const int accuracy = pi.accuracy;
          const std::string psrc_tag = get_psrc_tag(xg, type, accuracy);
          const std::string fn_sp = psel_path + "/" + psrc_tag + ".lat";
          if (does_file_exist_sync_node(fn_sp)) {
            continue;
          }
          if (check_status()) {
            sfw.init();
            qquit(ssprintf("collect_prop_psrc_light %ld", i));
          }
          TIMER_VERBOSE("collect_prop_iter");
          const std::string old_path =
              get_old_psrc_prop_path(old_job_tag, traj, xg, type, accuracy);
          SelectedPoints<WilsonMatrix> sp;
          load_selected_points_complex(sp, old_path + ".lat");
          Propagator4d prop;
          const long total_bytes = read_selected_field_double_from_float(
              prop, old_path + ".sfield", fsel);
          qassert(total_bytes > 0);
          set_field_selected(prop, sp, psel);
          write_float_from_double(sfw, psrc_tag, prop, sbs);
          sync_node();
          save_selected_points_complex(sp, fn_sp);
          Timer::autodisplay();
        }
      }
    } else if (is_final_run_ensemble(job_tag)) {
      const std::string old_path = get_old_psrc_prop_path(old_job_tag, traj);
      ShuffledFieldsReader sfr(old_path + "/sparse_point_src_props_light");
      ShuffledFieldsReader sfr_psel(old_path +
                                    "/sparse_psrc_point_src_props_light");
      for (long i = 0; i < (long)pis.size(); ++i) {
        const PointInfo& pi = pis[i];
        if (pi.type == type) {
          const Coordinate xg = pi.xg;
          const int accuracy = pi.accuracy;
          const std::string psrc_tag = get_psrc_tag(xg, type, accuracy);
          const std::string fn_sp = psel_path + "/" + psrc_tag + ".lat";
          if (does_file_exist_sync_node(fn_sp)) {
            continue;
          }
          if (check_status()) {
            sfw.init();
            qquit(ssprintf("collect_prop_psrc_light %ld", i));
          }
          TIMER_VERBOSE("collect_prop_iter");
          const std::string tag = ssprintf(
              "point_src_prop ; job_tag=%s ; traj=%d ; xg=%s ; type=%d ; "
              "acc=%d",
              old_job_tag.c_str(), traj, show(xg).c_str(), type, accuracy);
          SelectedPoints<WilsonMatrix> sp;
          {
            Propagator4d prop_psel;
            const long total_bytes_psel =
                read_double_from_float(sfr_psel, tag, prop_psel);
            qassert(total_bytes_psel > 0);
            set_selected_points(sp, prop_psel, psel);
          }
          Propagator4d prop;
          const long total_bytes = read_double_from_float(sfr, tag, prop);
          qassert(total_bytes > 0);
          set_field_selected(prop, sp, psel);
          write_float_from_double(sfw, psrc_tag, prop, sbs);
          sync_node();
          save_selected_points_complex(sp, fn_sp);
          Timer::autodisplay();
        }
      }
    } else {
      qassert(false);
    }
    sfw.init();
    qrename_info(path + ".acc", path);
    release_lock();
  }
  Timer::autodisplay();
}

inline void collect_prop_psrc_strange(const std::string& job_tag,
                                      const int traj)
{
  check_sigint();
  check_time_limit();
  const int type = 1;  // strange quark
  const std::string old_job_tag = get_old_job_tag(job_tag);
  const std::string path = get_prop_psrc_strange_path(job_tag, traj);
  const std::string psel_path = get_psel_prop_psrc_strange_path(job_tag, traj);
  if (does_file_exist_sync_node(path) and
      does_file_exist_sync_node(psel_path)) {
    return;
  }
  if (not check_prop_psrc(job_tag, traj, type)) {
    return;
  }
  if (not obtain_lock(
          ssprintf("lock-prop-psrc-strange-%s-%d", job_tag.c_str(), traj))) {
    return;
  }
  {
    setup(job_tag, traj);
    TIMER_VERBOSE("collect_prop_psrc_strange");
    qmkdir_info(ssprintf("data/prop-psrc-strange"));
    qmkdir_info(ssprintf("data/prop-psrc-strange/%s", job_tag.c_str()));
    qmkdir_info(ssprintf("data/psel-prop-psrc-strange"));
    qmkdir_info(ssprintf("data/psel-prop-psrc-strange/%s", job_tag.c_str()));
    qmkdir_info(psel_path);
    const std::vector<PointInfo> pis = get_point_src_info(job_tag, traj);
    std::vector<Coordinate> psel;
    FieldSelection fsel;
    set_sparse_parameters(psel, fsel, job_tag, traj);
    const Coordinate new_size_node = Coordinate(1, 1, 1, 16);
    ShuffledBitSet sbs = mk_shuffled_bitset(fsel.f_rank, psel, new_size_node);
    ShuffledFieldsWriter sfw(path + ".acc", new_size_node, true);
    if (is_old_collection_ensemble(job_tag)) {
      for (long i = 0; i < (long)pis.size(); ++i) {
        const PointInfo& pi = pis[i];
        if (pi.type == type) {
          const Coordinate xg = pi.xg;
          const int accuracy = pi.accuracy;
          const std::string psrc_tag = get_psrc_tag(xg, type, accuracy);
          const std::string fn_sp = psel_path + "/" + psrc_tag + ".lat";
          if (does_file_exist_sync_node(fn_sp)) {
            continue;
          }
          if (check_status()) {
            sfw.init();
            qquit(ssprintf("collect_prop_psrc_strange %ld", i));
          }
          TIMER_VERBOSE("collect_prop_iter");
          const std::string old_path =
              get_old_psrc_prop_path(old_job_tag, traj, xg, type, accuracy);
          SelectedPoints<WilsonMatrix> sp;
          load_selected_points_complex(sp, old_path + ".lat");
          Propagator4d prop;
          const long total_bytes = read_selected_field_double_from_float(
              prop, old_path + ".sfield", fsel);
          qassert(total_bytes > 0);
          set_field_selected(prop, sp, psel);
          write_float_from_double(sfw, psrc_tag, prop, sbs);
          sync_node();
          save_selected_points_complex(sp, fn_sp);
          Timer::autodisplay();
        }
      }
    } else if (is_final_run_ensemble(job_tag)) {
      const std::string old_path = get_old_psrc_prop_path(old_job_tag, traj);
      ShuffledFieldsReader sfr(old_path + "/sparse_point_src_props_strange");
      ShuffledFieldsReader sfr_psel(old_path +
                                    "/sparse_psrc_point_src_props_strange");
      ShuffledFieldsReader sfr_new, sfr_new_psel;
      if (job_tag == "48I") {
        sfr_new.init(old_path + "/sparse_point_src_props_strange_new2");
        sfr_new_psel.init(old_path +
                          "/sparse_psrc_point_src_props_strange_new2");
      }
      for (long i = 0; i < (long)pis.size(); ++i) {
        const PointInfo& pi = pis[i];
        if (pi.type == type) {
          const Coordinate xg = pi.xg;
          const int accuracy = pi.accuracy;
          const std::string psrc_tag = get_psrc_tag(xg, type, accuracy);
          const std::string fn_sp = psel_path + "/" + psrc_tag + ".lat";
          if (does_file_exist_sync_node(fn_sp)) {
            continue;
          }
          if (check_status()) {
            sfw.init();
            qquit(ssprintf("collect_prop_psrc_strange %ld", i));
          }
          TIMER_VERBOSE("collect_prop_iter");
          const std::string tag = ssprintf(
              "point_src_prop ; job_tag=%s ; traj=%d ; xg=%s ; type=%d ; "
              "acc=%d",
              old_job_tag.c_str(), traj, show(xg).c_str(), type, accuracy);
          SelectedPoints<WilsonMatrix> sp;
          {
            Propagator4d prop_psel;
            long total_bytes_psel = 0;
            if (job_tag == "48I") {
              total_bytes_psel =
                  read_double_from_float(sfr_new_psel, tag, prop_psel);
            }
            if (total_bytes_psel == 0) {
              total_bytes_psel =
                  read_double_from_float(sfr_psel, tag, prop_psel);
            }
            qassert(total_bytes_psel > 0);
            set_selected_points(sp, prop_psel, psel);
          }
          Propagator4d prop;
          long total_bytes = 0;
          if (job_tag == "48I") {
            total_bytes = read_double_from_float(sfr_new, tag, prop);
          }
          if (total_bytes == 0) {
            total_bytes = read_double_from_float(sfr, tag, prop);
          }
          qassert(total_bytes > 0);
          set_field_selected(prop, sp, psel);
          write_float_from_double(sfw, psrc_tag, prop, sbs);
          sync_node();
          save_selected_points_complex(sp, fn_sp);
          Timer::autodisplay();
        }
      }
    } else {
      qassert(false);
    }
    sfw.init();
    qrename_info(path + ".acc", path);
    release_lock();
  }
  Timer::autodisplay();
}

inline void collect_prop_psrc_exact(const std::string& job_tag, const int traj)
{
  check_sigint();
  check_time_limit();
  const int accuracy = 2;  // exact prop
  const std::string old_job_tag = get_old_job_tag(job_tag);
  const std::string path = get_prop_psrc_exact_path(job_tag, traj);
  const std::string psel_path = get_psel_prop_psrc_exact_path(job_tag, traj);
  if (does_file_exist_sync_node(path) and
      does_file_exist_sync_node(psel_path)) {
    return;
  }
  if (not check_prop_psrc(job_tag, traj, 0)) {
    return;
  }
  if (not check_prop_psrc(job_tag, traj, 1)) {
    return;
  }
  if (not obtain_lock(
          ssprintf("lock-prop-psrc-exact-%s-%d", job_tag.c_str(), traj))) {
    return;
  }
  {
    setup(job_tag, traj);
    TIMER_VERBOSE("collect_prop_psrc_exact");
    qmkdir_info(ssprintf("data/prop-psrc-exact"));
    qmkdir_info(ssprintf("data/prop-psrc-exact/%s", job_tag.c_str()));
    qmkdir_info(ssprintf("data/psel-prop-psrc-exact"));
    qmkdir_info(ssprintf("data/psel-prop-psrc-exact/%s", job_tag.c_str()));
    qmkdir_info(psel_path);
    const std::vector<PointInfo> pis = get_point_src_info(job_tag, traj);
    std::vector<Coordinate> psel;
    FieldSelection fsel;
    set_sparse_parameters(psel, fsel, job_tag, traj);
    const Coordinate new_size_node = Coordinate(1, 1, 1, 8);
    ShuffledBitSet sbs = mk_shuffled_bitset(fsel.f_rank, psel, new_size_node);
    ShuffledFieldsWriter sfw(path + ".acc", new_size_node, true);
    if (is_old_collection_ensemble(job_tag)) {
      for (long i = 0; i < (long)pis.size(); ++i) {
        const PointInfo& pi = pis[i];
        if (pi.accuracy == accuracy) {
          const Coordinate xg = pi.xg;
          const int type = pi.type;
          const std::string psrc_tag = get_psrc_tag(xg, type, accuracy);
          const std::string fn_sp = psel_path + "/" + psrc_tag + ".lat";
          if (does_file_exist_sync_node(fn_sp)) {
            continue;
          }
          if (check_status()) {
            sfw.init();
            qquit(ssprintf("collect_prop_psrc_exact %ld", i));
          }
          TIMER_VERBOSE("collect_prop_iter");
          const std::string old_path =
              get_old_psrc_prop_path(old_job_tag, traj, xg, type, accuracy);
          SelectedPoints<WilsonMatrix> sp;
          load_selected_points_complex(sp, old_path + ".lat");
          Propagator4d prop;
          const long total_bytes = read_selected_field_double_from_float(
              prop, old_path + ".sfield", fsel);
          qassert(total_bytes > 0);
          set_field_selected(prop, sp, psel);
          write_float_from_double(sfw, psrc_tag, prop, sbs);
          sync_node();
          save_selected_points_complex(sp, fn_sp);
          Timer::autodisplay();
        }
      }
    } else if (is_final_run_ensemble(job_tag)) {
      const std::string old_path = get_old_psrc_prop_path(old_job_tag, traj);
      ShuffledFieldsReader sfr(old_path + "/sparse_point_src_props_light");
      ShuffledFieldsReader sfr_s(old_path + "/sparse_point_src_props_strange");
      ShuffledFieldsReader sfr_psel(old_path +
                                    "/sparse_psrc_point_src_props_light");
      ShuffledFieldsReader sfr_s_psel(old_path +
                                      "/sparse_psrc_point_src_props_strange");
      for (long i = 0; i < (long)pis.size(); ++i) {
        const PointInfo& pi = pis[i];
        if (pi.accuracy == accuracy) {
          const Coordinate xg = pi.xg;
          const int type = pi.type;
          const std::string psrc_tag = get_psrc_tag(xg, type, accuracy);
          const std::string fn_sp = psel_path + "/" + psrc_tag + ".lat";
          if (does_file_exist_sync_node(fn_sp)) {
            continue;
          }
          if (check_status()) {
            sfw.init();
            qquit(ssprintf("collect_prop_psrc_exact %ld", i));
          }
          TIMER_VERBOSE("collect_prop_iter");
          const std::string tag = ssprintf(
              "point_src_prop ; job_tag=%s ; traj=%d ; xg=%s ; type=%d ; "
              "acc=%d",
              old_job_tag.c_str(), traj, show(xg).c_str(), type, accuracy);
          SelectedPoints<WilsonMatrix> sp;
          {
            Propagator4d prop_psel;
            long total_bytes_psel = 0;
            total_bytes_psel = read_double_from_float(sfr_psel, tag, prop_psel);
            if (total_bytes_psel == 0) {
              total_bytes_psel =
                  read_double_from_float(sfr_s_psel, tag, prop_psel);
            }
            qassert(total_bytes_psel > 0);
            set_selected_points(sp, prop_psel, psel);
          }
          Propagator4d prop;
          long total_bytes = 0;
          total_bytes = read_double_from_float(sfr, tag, prop);
          if (total_bytes == 0) {
            total_bytes = read_double_from_float(sfr_s, tag, prop);
          }
          qassert(total_bytes > 0);
          set_field_selected(prop, sp, psel);
          write_float_from_double(sfw, psrc_tag, prop, sbs);
          sync_node();
          save_selected_points_complex(sp, fn_sp);
          Timer::autodisplay();
        }
      }
    } else {
      qassert(false);
    }
    sfw.init();
    qrename_info(path + ".acc", path);
    release_lock();
  }
  Timer::autodisplay();
}

inline bool check_prop_wsrc(const std::string& job_tag, const int traj,
                            const int type)
{
  TIMER_VERBOSE("check_prop_wsrc");
  const std::string old_job_tag = get_old_job_tag(job_tag);
  if (is_old_collection_ensemble(job_tag)) {
    const std::string& path = get_old_wsrc_prop_path(old_job_tag, traj);
    if (type == 0) {
      return does_file_exist_sync_node(path + "/checkpoint.txt");
    } else if (type == 1) {
      const Coordinate total_site = get_total_site(job_tag);
      if (not does_file_exist_sync_node(path + "/checkpoint.txt")) {
        return false;
      }
      for (int i = 0; i < total_site[3]; ++i) {
        if (not does_file_exist_sync_node(
                path + ssprintf("/tslice=%d ; type=1 ; accuracy=1.lat", i))) {
          return false;
        }
      }
      return true;
    }
  } else if (is_final_run_ensemble(job_tag)) {
    const std::string& path = get_old_wsrc_prop_path(old_job_tag, traj);
    if (type == 0) {
      return does_file_exist_sync_node(path + "/sparse_wall_src_props_light");
    } else if (type == 1) {
      return does_file_exist_sync_node(path + "/sparse_wall_src_props_strange");
    } else {
      qassert(false);
    }
  } else {
    qassert(false);
  }
  return false;
}

inline void collect_prop_wsrc_light(const std::string& job_tag, const int traj)
{
  check_sigint();
  check_time_limit();
  const int type = 0;  // light quark
  const std::string old_job_tag = get_old_job_tag(job_tag);
  const std::string path = get_prop_wsrc_light_path(job_tag, traj);
  const std::string psel_path = get_psel_prop_wsrc_light_path(job_tag, traj);
  if (does_file_exist_sync_node(path) and
      does_file_exist_sync_node(psel_path)) {
    return;
  }
  if (not check_prop_wsrc(job_tag, traj, type)) {
    return;
  }
  if (not obtain_lock(
          ssprintf("lock-prop-wsrc-light-%s-%d", job_tag.c_str(), traj))) {
    return;
  }
  {
    setup(job_tag, traj);
    TIMER_VERBOSE("collect_prop_wsrc_light");
    qmkdir_info(ssprintf("data/prop-wsrc-light"));
    qmkdir_info(ssprintf("data/prop-wsrc-light/%s", job_tag.c_str()));
    qmkdir_info(ssprintf("data/psel-prop-wsrc-light"));
    qmkdir_info(ssprintf("data/psel-prop-wsrc-light/%s", job_tag.c_str()));
    qmkdir_info(psel_path);
    const std::vector<PointInfo> pis = get_point_src_info(job_tag, traj);
    std::vector<Coordinate> psel;
    FieldSelection fsel;
    set_sparse_parameters(psel, fsel, job_tag, traj);
    const Coordinate new_size_node = Coordinate(1, 1, 1, 8);
    ShuffledBitSet sbs = mk_shuffled_bitset(fsel.f_rank, psel, new_size_node);
    ShuffledFieldsWriter sfw(path + ".acc", new_size_node, true);
    const Coordinate total_site = get_total_site(job_tag);
    if (is_old_collection_ensemble(job_tag)) {
      for (int accuracy = 2; accuracy >= 1; accuracy -= 1) {
        for (int tslice = 0; tslice < total_site[3]; ++tslice) {
          const std::string wsrc_tag = get_wsrc_tag(tslice, type, accuracy);
          const std::string fn_sp = psel_path + "/" + wsrc_tag + ".lat";
          if (does_file_exist_sync_node(fn_sp)) {
            continue;
          }
          if (check_status()) {
            sfw.init();
            qquit(ssprintf("collect_prop_wsrc_light tslice=%d accuracy=%d",
                           tslice, accuracy));
          }
          const std::string old_path =
              get_old_wsrc_prop_path(old_job_tag, traj, tslice, type, accuracy);
          if (accuracy == 2 and
              (not does_file_exist_sync_node(old_path + ".lat"))) {
            continue;
          }
          TIMER_VERBOSE("collect_prop_iter");
          qassert(does_file_exist_sync_node(old_path + ".lat"));
          qassert(does_file_exist_sync_node(old_path + ".sfield"));
          SelectedPoints<WilsonMatrix> sp;
          load_selected_points_complex(sp, old_path + ".lat");
          Propagator4d prop;
          const long total_bytes = read_selected_field_double_from_float(
              prop, old_path + ".sfield", fsel);
          qassert(total_bytes > 0);
          set_field_selected(prop, sp, psel);
          write_float_from_double(sfw, wsrc_tag, prop, sbs);
          sync_node();
          save_selected_points_complex(sp, fn_sp);
          Timer::autodisplay();
        }
      }
    } else if (is_final_run_ensemble(job_tag)) {
      const std::string old_path = get_old_wsrc_prop_path(old_job_tag, traj);
      ShuffledFieldsReader sfr(old_path + "/sparse_wall_src_props_light");
      ShuffledFieldsReader sfr_psel(old_path +
                                    "/sparse_psrc_wall_src_props_light");
      for (int accuracy = 2; accuracy >= 1; accuracy -= 1) {
        for (int tslice = 0; tslice < total_site[3]; ++tslice) {
          const std::string wsrc_tag = get_wsrc_tag(tslice, type, accuracy);
          const std::string fn_sp = psel_path + "/" + wsrc_tag + ".lat";
          if (does_file_exist_sync_node(fn_sp)) {
            continue;
          }
          if (check_status()) {
            sfw.init();
            qquit(ssprintf("collect_prop_wsrc_light tslice=%d accuracy=%d",
                           tslice, accuracy));
          }
          const std::string tag = ssprintf(
              "wall_src_prop ; job_tag=%s ; traj=%d ; tslice=%d ; type=%d ; "
              "acc=%d",
              old_job_tag.c_str(), traj, tslice, type, accuracy);
          SelectedPoints<WilsonMatrix> sp;
          Propagator4d prop_psel;
          const long total_bytes_psel =
              read_double_from_float(sfr_psel, tag, prop_psel);
          if (accuracy == 2 and total_bytes_psel == 0) {
            continue;
          }
          TIMER_VERBOSE("collect_prop_iter");
          qassert(total_bytes_psel > 0);
          set_selected_points(sp, prop_psel, psel);
          prop_psel.init();
          Propagator4d prop;
          const long total_bytes = read_double_from_float(sfr, tag, prop);
          qassert(total_bytes > 0);
          set_field_selected(prop, sp, psel);
          write_float_from_double(sfw, wsrc_tag, prop, sbs);
          sync_node();
          save_selected_points_complex(sp, fn_sp);
          Timer::autodisplay();
        }
      }
    } else {
      qassert(false);
    }
    sfw.init();
    qrename_info(path + ".acc", path);
    release_lock();
  }
  Timer::autodisplay();
}

inline void collect_prop_wsrc_strange(const std::string& job_tag,
                                      const int traj)
{
  check_sigint();
  check_time_limit();
  const int type = 1;  // strange quark
  const std::string old_job_tag = get_old_job_tag(job_tag);
  const std::string path = get_prop_wsrc_strange_path(job_tag, traj);
  const std::string psel_path = get_psel_prop_wsrc_strange_path(job_tag, traj);
  if (does_file_exist_sync_node(path) and
      does_file_exist_sync_node(psel_path)) {
    return;
  }
  if (not check_prop_wsrc(job_tag, traj, type)) {
    return;
  }
  if (not obtain_lock(
          ssprintf("lock-prop-wsrc-strange-%s-%d", job_tag.c_str(), traj))) {
    return;
  }
  {
    setup(job_tag, traj);
    TIMER_VERBOSE("collect_prop_wsrc_strange");
    qmkdir_info(ssprintf("data/prop-wsrc-strange"));
    qmkdir_info(ssprintf("data/prop-wsrc-strange/%s", job_tag.c_str()));
    qmkdir_info(ssprintf("data/psel-prop-wsrc-strange"));
    qmkdir_info(ssprintf("data/psel-prop-wsrc-strange/%s", job_tag.c_str()));
    qmkdir_info(psel_path);
    const std::vector<PointInfo> pis = get_point_src_info(job_tag, traj);
    std::vector<Coordinate> psel;
    FieldSelection fsel;
    set_sparse_parameters(psel, fsel, job_tag, traj);
    const Coordinate new_size_node = Coordinate(1, 1, 1, 8);
    ShuffledBitSet sbs = mk_shuffled_bitset(fsel.f_rank, psel, new_size_node);
    ShuffledFieldsWriter sfw(path + ".acc", new_size_node, true);
    const Coordinate total_site = get_total_site(job_tag);
    if (is_old_collection_ensemble(job_tag)) {
      for (int accuracy = 2; accuracy >= 1; accuracy -= 1) {
        for (int tslice = 0; tslice < total_site[3]; ++tslice) {
          const std::string wsrc_tag = get_wsrc_tag(tslice, type, accuracy);
          const std::string fn_sp = psel_path + "/" + wsrc_tag + ".lat";
          if (does_file_exist_sync_node(fn_sp)) {
            continue;
          }
          if (check_status()) {
            sfw.init();
            qquit(ssprintf("collect_prop_wsrc_strange tslice=%d accuracy=%d",
                           tslice, accuracy));
          }
          const std::string old_path =
              get_old_wsrc_prop_path(old_job_tag, traj, tslice, type, accuracy);
          if (accuracy == 2 and
              (not does_file_exist_sync_node(old_path + ".lat"))) {
            continue;
          }
          TIMER_VERBOSE("collect_prop_iter");
          qassert(does_file_exist_sync_node(old_path + ".lat"));
          qassert(does_file_exist_sync_node(old_path + ".sfield"));
          SelectedPoints<WilsonMatrix> sp;
          load_selected_points_complex(sp, old_path + ".lat");
          Propagator4d prop;
          const long total_bytes = read_selected_field_double_from_float(
              prop, old_path + ".sfield", fsel);
          qassert(total_bytes > 0);
          set_field_selected(prop, sp, psel);
          write_float_from_double(sfw, wsrc_tag, prop, sbs);
          sync_node();
          save_selected_points_complex(sp, fn_sp);
          Timer::autodisplay();
        }
      }
    } else if (is_final_run_ensemble(job_tag)) {
      const std::string old_path = get_old_wsrc_prop_path(old_job_tag, traj);
      ShuffledFieldsReader sfr(old_path + "/sparse_wall_src_props_strange");
      ShuffledFieldsReader sfr_psel(old_path +
                                    "/sparse_psrc_wall_src_props_strange");
      for (int accuracy = 2; accuracy >= 1; accuracy -= 1) {
        for (int tslice = 0; tslice < total_site[3]; ++tslice) {
          const std::string wsrc_tag = get_wsrc_tag(tslice, type, accuracy);
          const std::string fn_sp = psel_path + "/" + wsrc_tag + ".lat";
          if (does_file_exist_sync_node(fn_sp)) {
            continue;
          }
          if (check_status()) {
            sfw.init();
            qquit(ssprintf("collect_prop_wsrc_strange tslice=%d accuracy=%d",
                           tslice, accuracy));
          }
          const std::string tag = ssprintf(
              "wall_src_prop ; job_tag=%s ; traj=%d ; tslice=%d ; type=%d ; "
              "acc=%d",
              old_job_tag.c_str(), traj, tslice, type, accuracy);
          SelectedPoints<WilsonMatrix> sp;
          Propagator4d prop_psel;
          const long total_bytes_psel =
              read_double_from_float(sfr_psel, tag, prop_psel);
          if (accuracy == 2 and total_bytes_psel == 0) {
            continue;
          }
          TIMER_VERBOSE("collect_prop_iter");
          qassert(total_bytes_psel > 0);
          set_selected_points(sp, prop_psel, psel);
          prop_psel.init();
          Propagator4d prop;
          const long total_bytes = read_double_from_float(sfr, tag, prop);
          qassert(total_bytes > 0);
          set_field_selected(prop, sp, psel);
          write_float_from_double(sfw, wsrc_tag, prop, sbs);
          sync_node();
          save_selected_points_complex(sp, fn_sp);
          Timer::autodisplay();
        }
      }
    } else {
      qassert(false);
    }
    sfw.init();
    qrename_info(path + ".acc", path);
    release_lock();
  }
  Timer::autodisplay();
}

inline void compute_traj(const std::string& job_tag, const int traj)
{
  {
    setup(job_tag);
    TIMER_VERBOSE("compute_traj");
    // ADJUST ME
    collect_field_selection(job_tag, traj);
    collect_gauge_transform(job_tag, traj);
    collect_prop_wsrc_strange(job_tag, traj);
    collect_prop_wsrc_light(job_tag, traj);
    collect_prop_psrc_strange(job_tag, traj);
    collect_prop_psrc_light(job_tag, traj);
    collect_prop_psrc_exact(job_tag, traj);
    //
  }
  if (get_total_time() > 60.0) {
    Timer::display();
  }
  Timer::reset();
}

inline void compute(const std::string& job_tag)
{
  TIMER_VERBOSE("compute");
  {
    setup(job_tag);
    // ADJUST ME
    collect_point_distribution(job_tag);
    collect_simulated_pis(job_tag);
    collect_pis(job_tag);
    collect_pcs(job_tag);
    //
    if (get_total_time() > 60.0) {
      Timer::display();
    }
    Timer::reset();
  }
  const std::string old_job_tag = get_old_job_tag(job_tag);
  const std::vector<int> trajs = get_todo_trajs(old_job_tag);
  for (int i = 0; i < (int)trajs.size(); ++i) {
    const int traj = trajs[i];
    compute_traj(job_tag, traj);
  }
}

}  // namespace qlat

int main(int argc, char* argv[])
{
  using namespace qlat;
  std::vector<Coordinate> size_node_list;
  size_node_list.push_back(Coordinate(1, 1, 1, 1));
  size_node_list.push_back(Coordinate(1, 1, 1, 2));
  size_node_list.push_back(Coordinate(1, 1, 1, 4));
  size_node_list.push_back(Coordinate(1, 1, 1, 8));
  size_node_list.push_back(Coordinate(1, 1, 1, 16));
  size_node_list.push_back(Coordinate(1, 1, 2, 16));
  size_node_list.push_back(Coordinate(1, 2, 2, 16));
  size_node_list.push_back(Coordinate(2, 2, 2, 16));
  begin(&argc, &argv, size_node_list);
  setup_log_idx();
  setup();
  qmkdir_info(ssprintf("data"));
  std::vector<std::string> job_tags;
  // SADJUST ME
  job_tags.push_back("24D");
  job_tags.push_back("24DH");
  job_tags.push_back("32Dfine");
  job_tags.push_back("32D");
  job_tags.push_back("48I");
  job_tags.push_back("64I");
  //
  for (int k = 0; k < (int)job_tags.size(); ++k) {
    const std::string& job_tag = job_tags[k];
    compute(job_tag);
  }
  end();
  return 0;
}
