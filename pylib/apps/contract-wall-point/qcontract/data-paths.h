#pragma once

#include <qlat/qlat.h>

namespace qlat
{  //

inline std::string& get_data_path()
{
  static std::string path = "data";
  return path;
}

inline std::string get_point_src_info_path(const std::string& job_tag,
                                           const int traj)
{
  return ssprintf("%s/point-src-info/%s/traj=%d.txt", get_data_path().c_str(),
                  job_tag.c_str(), traj);
}

inline std::string get_point_selection_path(const std::string& job_tag,
                                            const int traj)
{
  return ssprintf("%s/point-selection/%s/traj=%d.txt", get_data_path().c_str(),
                  job_tag.c_str(), traj);
}

inline std::string get_point_src_info_simulated_path(const std::string& job_tag,
                                                     const int traj)
{
  return ssprintf("%s/point-src-info-simulated/%s/traj=%d.txt",
                  get_data_path().c_str(), job_tag.c_str(), traj);
}

inline std::string get_point_distribution_path(const std::string& job_tag)
{
  return ssprintf("%s/point-distribution/%s/point-distribution.txt",
                  get_data_path().c_str(), job_tag.c_str());
}

inline std::string get_field_selection_path(const std::string& job_tag,
                                            const int traj)
{
  return ssprintf("%s/field-selection/%s/traj=%d.field",
                  get_data_path().c_str(), job_tag.c_str(), traj);
}

inline std::string get_gauge_transform_path(const std::string& job_tag,
                                            const int traj)
{
  return ssprintf("%s/gauge-transform/%s/traj=%d.field",
                  get_data_path().c_str(), job_tag.c_str(), traj);
}

inline std::string get_prop_psrc_light_path(const std::string& job_tag,
                                            const int traj)
{
  return ssprintf("%s/prop-psrc-light/%s/traj=%d", get_data_path().c_str(),
                  job_tag.c_str(), traj);
}

inline std::string get_psel_prop_psrc_light_path(const std::string& job_tag,
                                                 const int traj)
{
  return ssprintf("%s/psel-prop-psrc-light/%s/traj=%d", get_data_path().c_str(),
                  job_tag.c_str(), traj);
}

inline std::string get_prop_psrc_strange_path(const std::string& job_tag,
                                              const int traj)
{
  return ssprintf("%s/prop-psrc-strange/%s/traj=%d", get_data_path().c_str(),
                  job_tag.c_str(), traj);
}

inline std::string get_psel_prop_psrc_strange_path(const std::string& job_tag,
                                                   const int traj)
{
  return ssprintf("%s/psel-prop-psrc-strange/%s/traj=%d",
                  get_data_path().c_str(), job_tag.c_str(), traj);
}

inline std::string get_prop_psrc_exact_path(const std::string& job_tag,
                                            const int traj)
{
  return ssprintf("%s/prop-psrc-exact/%s/traj=%d", get_data_path().c_str(),
                  job_tag.c_str(), traj);
}

inline std::string get_psel_prop_psrc_exact_path(const std::string& job_tag,
                                                 const int traj)
{
  return ssprintf("%s/psel-prop-psrc-exact/%s/traj=%d", get_data_path().c_str(),
                  job_tag.c_str(), traj);
}

inline std::string get_psrc_tag(const Coordinate& xg, const int type,
                                const int accuracy)
{
  return ssprintf("xg=(%d,%d,%d,%d) ; type=%d ; accuracy=%d", xg[0], xg[1],
                  xg[2], xg[3], type, accuracy);
}

inline std::string get_prop_wsrc_light_path(const std::string& job_tag,
                                            const int traj)
{
  return ssprintf("%s/prop-wsrc-light/%s/traj=%d", get_data_path().c_str(),
                  job_tag.c_str(), traj);
}

inline std::string get_psel_prop_wsrc_light_path(const std::string& job_tag,
                                                 const int traj)
{
  return ssprintf("%s/psel-prop-wsrc-light/%s/traj=%d", get_data_path().c_str(),
                  job_tag.c_str(), traj);
}

inline std::string get_prop_wsrc_strange_path(const std::string& job_tag,
                                              const int traj)
{
  return ssprintf("%s/prop-wsrc-strange/%s/traj=%d", get_data_path().c_str(),
                  job_tag.c_str(), traj);
}

inline std::string get_psel_prop_wsrc_strange_path(const std::string& job_tag,
                                                   const int traj)
{
  return ssprintf("%s/psel-prop-wsrc-strange/%s/traj=%d",
                  get_data_path().c_str(), job_tag.c_str(), traj);
}

inline std::string get_wsrc_tag(const int tslice, const int type,
                                const int accuracy)
{
  return ssprintf("tslice=%d ; type=%d ; accuracy=%d", tslice, type, accuracy);
}

inline std::string get_prop_psrc_path(const std::string& job_tag,
                                      const int traj, const int type)
{
  if (type == 0) {
    return get_prop_psrc_light_path(job_tag, traj);
  } else if (type == 1) {
    return get_prop_psrc_strange_path(job_tag, traj);
  } else {
    qassert(false);
    return "";
  }
}

inline std::string get_psel_prop_psrc_path(const std::string& job_tag,
                                           const int traj, const int type)
{
  if (type == 0) {
    return get_psel_prop_psrc_light_path(job_tag, traj);
  } else if (type == 1) {
    return get_psel_prop_psrc_strange_path(job_tag, traj);
  } else {
    qassert(false);
    return "";
  }
}

inline std::string get_prop_wsrc_path(const std::string& job_tag,
                                      const int traj, const int type)
{
  if (type == 0) {
    return get_prop_wsrc_light_path(job_tag, traj);
  } else if (type == 1) {
    return get_prop_wsrc_strange_path(job_tag, traj);
  } else {
    qassert(false);
    return "";
  }
}

inline std::string get_psel_prop_wsrc_path(const std::string& job_tag,
                                           const int traj, const int type)
{
  if (type == 0) {
    return get_psel_prop_wsrc_light_path(job_tag, traj);
  } else if (type == 1) {
    return get_psel_prop_wsrc_strange_path(job_tag, traj);
  } else {
    qassert(false);
    return "";
  }
}

}  // namespace qlat
