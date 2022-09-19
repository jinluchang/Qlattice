#pragma once

#include "data-load-base.h"

namespace qlat
{  //

struct WallInfo {
  int tslice;
  int type;      // index of a array; 0 is light, 1 is strange
  int accuracy;  // index of a array; 0 is sloppy, higher is more accurate
  //
  WallInfo() { init(); }
  WallInfo(const int tslice_, const int type_, const int accuracy_)
  {
    init(tslice_, type_, accuracy_);
  }
  //
  void init() { memset(this, 0, sizeof(WallInfo)); }
  void init(const int tslice_, const int type_, const int accuracy_)
  {
    tslice = tslice_;
    type = type_;
    accuracy = accuracy_;
  }
};

inline bool operator==(const WallInfo& x, const WallInfo& y)
{
  return 0 == memcmp(&x, &y, sizeof(WallInfo));
}

inline std::vector<WallInfo> load_wall_src_info(const std::string& path)
{
  TIMER_VERBOSE("load_wall_src_info");
  std::vector<WallInfo> wis;
  if (get_id_node() == 0) {
    const std::vector<std::string> lines = qgetlines(path);
    for (int i = 0; i < (int)lines.size(); ++i) {
      if (lines[i] == "") {
        continue;
      }
      const std::vector<std::string> strs = split_line_with_spaces(lines[i]);
      qassert(strs.size() == 4);
      qassert(read_long(strs[0]) == i);
      const WallInfo wi(read_long(strs[1]), read_long(strs[2]),
                        read_long(strs[3]));
      wis.push_back(wi);
    }
  }
  bcast(wis);
  return wis;
}

inline void save_wall_src_info(const std::vector<WallInfo>& wis,
                               const std::string& path)
{
  TIMER_VERBOSE("save_wall_src_info");
  if (get_id_node() == 0) {
    std::ostringstream out;
    for (int i = 0; i < (int)wis.size(); ++i) {
      const WallInfo& wi = wis[i];
      out << ssprintf("%5d  %3d  %2d %2d\n", i, wi.tslice, wi.type,
                      wi.accuracy);
    }
    qtouch(path, out.str());
  }
  if (true) {
    const std::vector<WallInfo> wis_load = load_wall_src_info(path);
    qassert(wis_load == wis);
  }
}

inline std::string get_wall_src_info_path(const std::string& job_tag,
                                          const int traj, const int type)
{
  if (type == 0) {
    return ssprintf("data/wall-src-info-light/%s/traj=%d.txt", job_tag.c_str(),
                    traj);
  } else if (type == 1) {
    return ssprintf("data/wall-src-info-strange/%s/traj=%d.txt",
                    job_tag.c_str(), traj);
  } else {
    qassert(false);
    return "";
  }
}

inline bool check_wall_src_info(const std::string& job_tag, const int traj,
                                const int type)
{
  TIMER_VERBOSE("check_wall_src_info");
  return get_does_file_exist(get_wall_src_info_path(job_tag, traj, type));
}

typedef Cache<std::string, std::vector<WallInfo> > WallSrcInfoCache;

inline WallSrcInfoCache& get_wall_src_info_cache()
{
  static WallSrcInfoCache cache("WallSrcInfoCache", 8, 2);
  return cache;
}

inline const std::vector<WallInfo>& get_wall_src_info(
    const std::string& job_tag, const int traj, const int type)
{
  const std::string key = ssprintf("%s,%d,%d,wis", job_tag.c_str(), traj, type);
  WallSrcInfoCache& cache = get_wall_src_info_cache();
  if (not cache.has(key)) {
    TIMER_VERBOSE("get_wall_src_info");
    const std::string fn = get_wall_src_info_path(job_tag, traj, type);
    qassert(get_does_file_exist(fn));
    cache[key] = load_wall_src_info(fn);
  }
  return cache[key];
}

inline void compute_wall_src_info(const std::string& job_tag, const int traj,
                                  const int type)
{
  if (check_wall_src_info(job_tag, traj, type)) {
    return;
  }
  if (check_prop_wsrc(job_tag, traj, type)) {
    check_sigterm();
    check_time_limit();
    if (not obtain_lock(ssprintf("lock-wall-src-info-%s-%d-%d", job_tag.c_str(),
                                 traj, type))) {
      return;
    }
    setup(job_tag, traj);
    TIMER_VERBOSE("compute_wall_src_info");
    std::string path;
    if (type == 0) {
      path = "data/wall-src-info-light";
    } else if (type == 1) {
      path = "data/wall-src-info-strange";
    }
    qmkdir_info(path);
    qmkdir_info(path + "/" + job_tag);
    std::vector<WallInfo> wis;
    const Coordinate total_site = get_total_site(job_tag);
    for (int tslice = 0; tslice < total_site[3]; ++tslice) {
      qassert(get_does_prop_wsrc_exist(job_tag, traj, tslice, type, 1));
      wis.push_back(WallInfo(tslice, type, 1));
      if (get_does_prop_wsrc_exist(job_tag, traj, tslice, type, 2)) {
        wis.push_back(WallInfo(tslice, type, 2));
      }
    }
    save_wall_src_info(wis, get_wall_src_info_path(job_tag, traj, type));
    release_lock();
  }
}

}  // namespace qlat
