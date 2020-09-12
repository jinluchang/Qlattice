#pragma once

#include <qlat/qlat.h>

#include "qlat-setup.h"

namespace qlat
{  //

inline std::string get_new_job_tag(const std::string& old_job_tag)
{
  if (old_job_tag == "24D-0.00107") {
    return "24D";
  } else if (old_job_tag == "32D-0.00107") {
    return "32D";
  } else if (old_job_tag == "32Dfine-0.0001") {
    return "32Dfine";
  } else if (old_job_tag == "24D-0.0174") {
    return "24DH";
  } else if (old_job_tag == "48I-0.00078") {
    return "48I";
  } else if (old_job_tag == "64I-0.000678") {
    return "64I";
  }
  return old_job_tag;
}

inline std::string get_old_job_tag(const std::string& job_tag)
{
  if (job_tag == "24D") {
    return "24D-0.00107";
  } else if (job_tag == "32D") {
    return "32D-0.00107";
  } else if (job_tag == "32Dfine") {
    return "32Dfine-0.0001";
  } else if (job_tag == "24DH") {
    return "24D-0.0174";
  } else if (job_tag == "48I") {
    return "48I-0.00078";
  } else if (job_tag == "64I") {
    return "64I-0.000678";
  }
  return job_tag;
}

inline bool is_old_collection_ensemble(const std::string& job_tag)
{
  if (job_tag == "24D" or job_tag == "32D" or
      job_tag == "32Dfine" or job_tag == "24DH") {
    return true;
  } else {
    return false;
  }
}

inline bool is_final_run_ensemble(const std::string& job_tag)
{
  if (job_tag == "48I" or job_tag == "64I") {
    return true;
  } else {
    return false;
  }
}

inline std::string get_previous_data_path()
{
  return "/home/ljin/application/Public/Muon-GM2-cc/jobs/final-run/"
         "previous-data/";
}

inline std::string get_pis_distribution_path(const std::string& old_job_tag)
{
  std::string ret;
  if (old_job_tag == "64I-0.000678") {
    ret =
        "/home/ljin/application/Public/Muon-GM2-cc/jobs/final-run/64I/run/"
        "results/64I-0.000678/point-distribution.txt";
  } else if (old_job_tag == "48I-0.00078") {
    ret =
        "/home/ljin/application/Public/Muon-GM2-cc/jobs/final-run/48I/run/"
        "results/48I-0.00078/point-distribution.txt";
  } else if (old_job_tag == "32Dfine-0.0001") {
    ret =
        "/home/ljin/application/Public/Muon-GM2-cc/jobs/32Dfine/discon-2/"
        "results/points-distribution-summary.txt";
  } else if (old_job_tag == "24D-0.00107") {
    ret =
        "/home/ljin/application/Public/Muon-GM2-cc/jobs/24D/discon-1/results/"
        "points-distribution-summary.txt";
  } else if (old_job_tag == "24D-0.0174") {
    ret =
        "/home/ljin/application/Public/Muon-GM2-cc/jobs/24D-0.0174/discon-2/"
        "results/points-distribution-summary.txt";
  } else if (old_job_tag == "32D-0.00107") {
    ret =
        "/home/ljin/application/Public/Muon-GM2-cc/jobs/32D/discon-1/results/"
        "points-distribution-summary.txt";
  } else if (old_job_tag == "16I-0.01") {
    qassert(false);
  } else {
    qassert(false);
  }
  if (ret != "" and (not does_file_exist_sync_node(ret))) {
    ret = "";
  }
  return ret;
}

inline std::string get_simulated_pis_path(const std::string& old_job_tag)
{
  std::string ret;
  if (old_job_tag == "64I-0.000678") {
    ret =
        "/home/ljin/application/Public/Muon-GM2-cc/jobs/final-run/64I/run/"
        "results/64I-0.000678/point-distribution.dir";
  } else if (old_job_tag == "48I-0.00078") {
    ret =
        "/home/ljin/application/Public/Muon-GM2-cc/jobs/final-run/48I/run/"
        "results/48I-0.00078/point-distribution.dir";
  } else if (old_job_tag == "32Dfine-0.0001") {
    ret =
        "/home/ljin/application/Public/Muon-GM2-cc/jobs/32Dfine/discon-2/"
        "results/points-distribution";
  } else if (old_job_tag == "24D-0.00107") {
    ret =
        "/home/ljin/application/Public/Muon-GM2-cc/jobs/24D/discon-1/results/"
        "points-distribution";
  } else if (old_job_tag == "24D-0.0174") {
    ret =
        "/home/ljin/application/Public/Muon-GM2-cc/jobs/24D-0.0174/discon-2/"
        "results/points-distribution";
  } else if (old_job_tag == "32D-0.00107") {
    ret =
        "/home/ljin/application/Public/Muon-GM2-cc/jobs/32D/discon-1/results/"
        "points-distribution";
  } else if (old_job_tag == "16I-0.01") {
    qassert(false);
  } else {
    qassert(false);
  }
  if (ret != "" and (not does_file_exist_sync_node(ret))) {
    ret = "";
  }
  return ret;
}

inline std::string get_pis_path(const std::string& old_job_tag, const int traj)
{
  std::string ret;
  if (old_job_tag == "64I-0.000678") {
    ret = "/home/ljin/application/Public/Muon-GM2-cc/jobs/final-run/64I/run/" +
          get_job_path(old_job_tag, traj) + "/pis.txt";
  } else if (old_job_tag == "48I-0.00078") {
    ret = "/home/ljin/application/Public/Muon-GM2-cc/jobs/final-run/48I/run/" +
          get_job_path(old_job_tag, traj) + "/pis.txt";
  } else if (old_job_tag == "32Dfine-0.0001") {
    ret =
        get_previous_data_path() + get_job_path(old_job_tag, traj) + "/pis.txt";
  } else if (old_job_tag == "24D-0.00107") {
    ret =
        get_previous_data_path() + get_job_path(old_job_tag, traj) + "/pis.txt";
  } else if (old_job_tag == "24D-0.0174") {
    ret =
        get_previous_data_path() + get_job_path(old_job_tag, traj) + "/pis.txt";
  } else if (old_job_tag == "32D-0.00107") {
    ret =
        get_previous_data_path() + get_job_path(old_job_tag, traj) + "/pis.txt";
  } else if (old_job_tag == "16I-0.01") {
    qassert(false);
  } else {
    qassert(false);
  }
  if (ret != "" and (not does_file_exist_sync_node(ret))) {
    ret = "";
  }
  return ret;
}

inline std::string get_f_rank_path(const std::string& old_job_tag,
                                   const int traj)
{
  std::string ret;
  if (old_job_tag == "64I-0.000678") {
    ret = "/home/ljin/application/Public/Muon-GM2-cc/jobs/final-run/64I/run/" +
          get_job_path(old_job_tag, traj) + "/huge-data-keep/f-rank";
    return ret;
  } else if (old_job_tag == "48I-0.00078") {
    ret = "/home/ljin/application/Public/Muon-GM2-cc/jobs/final-run/48I/run/" +
          get_job_path(old_job_tag, traj) + "/huge-data-keep/f-rank";
    if (not does_file_exist_sync_node(ret)) {
      ret = ssprintf(
          "/home/ljin/application/Public/Muon-GM2-cc/jobs/final-run/"
          "previous-data/results/48I-0.00078/results=%d/f-rank.field",
          traj);
    }
    return ret;
  }
  const std::string job_path =
      get_previous_data_path() + get_job_path(old_job_tag, traj);
  return job_path + "/f-rank.field";
}

inline std::string get_old_gauge_transform_path(const std::string& old_job_tag,
                                                const int traj)
{
  std::string ret;
  if (old_job_tag == "64I-0.000678") {
    ret = "/home/ljin/application/Public/Muon-GM2-cc/jobs/final-run/64I/run/" +
          get_job_path(old_job_tag, traj) + "/huge-data-keep/gauge-transform";
    return ret;
  } else if (old_job_tag == "48I-0.00078") {
    ret = "/home/ljin/application/Public/Muon-GM2-cc/jobs/final-run/48I/run/" +
          get_job_path(old_job_tag, traj) + "/huge-data-keep/gauge-transform";
    if (not does_file_exist_sync_node(ret)) {
      ret = ssprintf(
          "/home/ljin/application/Public/Muon-GM2-cc/jobs/final-run/"
          "previous-data/results/48I-0.00078/results=%d/gauge-transform.field",
          traj);
    }
    return ret;
  }
  const std::string job_path =
      get_previous_data_path() + get_job_path(old_job_tag, traj);
  return job_path + "/gauge-transform.field";
}

inline std::string get_old_psrc_prop_path(const std::string& old_job_tag,
                                          const int traj)
{
  const std::string job_tag = get_new_job_tag(old_job_tag);
  if (is_old_collection_ensemble(job_tag)) {
    const std::string job_path =
        get_previous_data_path() + get_job_path(old_job_tag, traj);
    const std::string psrc_prop_path = job_path + "/prop-point-src";
    return psrc_prop_path;
  } else if (is_final_run_ensemble(job_tag)) {
    return ssprintf(
        "/home/ljin/application/Public/Muon-GM2-cc/jobs/final-run/%s/run/"
        "results/%s/results=%d/huge-data",
        job_tag.c_str(), old_job_tag.c_str(), traj);
  } else {
    qassert(false);
    return "";
  }
}

inline std::string get_old_psrc_prop_path(const std::string& old_job_tag,
                                          const int traj, const Coordinate& xg,
                                          const int type, const int accuracy)
{
  const std::string psrc_prop_path = get_old_psrc_prop_path(old_job_tag, traj);
  const std::string prop_path =
      psrc_prop_path + ssprintf("/xg=(%d,%d,%d,%d) ; type=%d ; accuracy=%d",
                                xg[0], xg[1], xg[2], xg[3], type, accuracy);
  return prop_path;
}

inline std::string get_old_wsrc_prop_path(const std::string& old_job_tag,
                                          const int traj)
{
  const std::string job_tag = get_new_job_tag(old_job_tag);
  if (is_old_collection_ensemble(job_tag)) {
    const std::string job_path =
        get_previous_data_path() + get_job_path(old_job_tag, traj);
    const std::string wsrc_prop_path = job_path + "/prop-wall-src";
    return wsrc_prop_path;
  } else if (is_final_run_ensemble(job_tag)) {
    return ssprintf(
        "/home/ljin/application/Public/Muon-GM2-cc/jobs/final-run/%s/run/"
        "results/%s/results=%d/huge-data",
        job_tag.c_str(), old_job_tag.c_str(), traj);
  } else {
    qassert(false);
    return "";
  }
}

inline std::string get_old_wsrc_prop_path(const std::string& old_job_tag,
                                          const int traj, const int tslice,
                                          const int type, const int accuracy)
// accuracy should be either 1 or 2 for wsrc prop
{
  const std::string wsrc_prop_path = get_old_wsrc_prop_path(old_job_tag, traj);
  const std::string prop_path =
      wsrc_prop_path +
      ssprintf("/tslice=%d ; type=%d ; accuracy=%d", tslice, type, accuracy);
  return prop_path;
}

}  // namespace qlat
