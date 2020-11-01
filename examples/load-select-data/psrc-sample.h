#pragma once

#include <qlat/qlat.h>

#include "qlat-setup.h"

namespace qlat
{  //

struct PointInfo {
  Coordinate xg;
  int type;      // index of a array; 0 is light, 1 is strange
  int accuracy;  // index of a array; 0 is sloppy, higher is more accurate
  //
  PointInfo() { init(); }
  PointInfo(const Coordinate& xg_, const int type_, const int accuracy_)
  {
    init(xg_, type_, accuracy_);
  }
  //
  void init() { memset(this, 0, sizeof(PointInfo)); }
  void init(const Coordinate& xg_, const int type_, const int accuracy_)
  {
    init();
    xg = xg_;
    type = type_;
    accuracy = accuracy_;
  }
};

inline RngState mk_rs_psrc_locations(const std::string& job_tag, const int traj)
{
  int results_id = traj;
  if (job_tag == "48I-0.00078") {
    results_id = traj / 5;
  } else if (job_tag == "64I-0.000678") {
    results_id = traj / 10;
  } else if (job_tag == "16I-0.01") {
    results_id = traj / 100;
  }
  int id = results_id * 8;
  RngState rs = RngState().split("rng_psrc_short_dis_dir").split(id);
  return rs;
}

inline int random_int(qlat::RngState& rs, const uint64_t limit)
// 0 <= randomInt < limit
{
  return rand_gen(rs) % limit;
}

inline void setPointCoordinatesWithReferenceLimit(
    RngState& rng, std::vector<Coordinate>& pcs, const int npoints,
    const Coordinate& total_site, const Coordinate& xgref,
    // refLimitSq default to be very large (basically no limit)
    // if a limit is imposed, it had better be smaller then L/4 to avoid around
    // boundary effect
    const int refLimitSq = sqr(1024))
{
  for (int i = 0; i < npoints; ++i) {
    Coordinate c;
    while (true) {
      for (int mu = 0; mu < 4; ++mu) {
        c[mu] = random_int(rng, total_site[mu]);
      }
      if (sqr(smod(c - xgref, total_site)) > refLimitSq) {
        continue;
      }
      bool duplicate = false;
      for (int k = 0; k < (int)pcs.size(); ++k) {
        if (c == pcs[k]) {
          duplicate = true;
          break;
        }
      }
      if (duplicate) {
        continue;
      }
      break;
    }
    pcs.push_back(c);
  }
}

inline void setPointCoordinatesSpheres(RngState& rng,
                                       std::vector<Coordinate>& pcs,
                                       const int ngroups,
                                       const int npoints,  // per group
                                       const Coordinate& total_site,
                                       const int sizeSphereSq)
{
  for (int gi = 0; gi < ngroups; ++gi) {
    Coordinate xgref;
    for (int mu = 0; mu < 4; ++mu) {
      xgref[mu] = random_int(rng, total_site[mu]);
    }
    setPointCoordinatesWithReferenceLimit(rng, pcs, npoints, total_site, xgref,
                                          sizeSphereSq);
  }
}

inline void pointInfoFromCoordinates(std::vector<PointInfo>& pis,
                                     const std::vector<Coordinate>& pcs,
                                     const int type, const int accuracy)
{
  PointInfo pi;
  pi.type = type;
  pi.accuracy = accuracy;
  for (int i = 0; i < (int)pcs.size(); ++i) {
    pi.xg = pcs[i];
    pis.push_back(pi);
  }
}

inline void sampleFromCoordinates(RngState& rng, std::vector<Coordinate>& spcs,
                                  const int size,
                                  const std::vector<Coordinate>& pcs)
{
  qassert(size <= (int)pcs.size());
  std::vector<bool> htable(pcs.size(), false);
  for (int i = 0; i < size; ++i) {
    bool filled = false;
    while (!filled) {
      int k = random_int(rng, pcs.size());
      if (!htable[k]) {
        spcs.push_back(pcs[k]);
        htable[k] = true;
        filled = true;
      }
    }
  }
}

inline bool operator==(const PointInfo& x, const PointInfo& y)
{
  return 0 == memcmp(&x, &y, sizeof(PointInfo));
}

inline bool pointInfoComp(const PointInfo& x, const PointInfo& y)
{
  for (int i = 0; i < 4; ++i) {
    if (x.xg[i] < y.xg[i]) {
      return true;
    } else if (x.xg[i] > y.xg[i]) {
      return false;
    }
  }
  if (x.type < y.type) {
    return true;
  } else if (x.type > y.type) {
    return false;
  } else {
    return x.accuracy < y.accuracy;
  }
}

inline std::vector<PointInfo> mk_point_infos(const std::string& job_tag,
                                             const int traj)
// interface function
{
  TIMER_VERBOSE("mk_point_infos");
  RngState rng = mk_rs_psrc_locations(job_tag, traj);
  const Coordinate total_site = get_total_site(job_tag);
  std::vector<PointInfo> pis;
  std::vector<Coordinate> pcs, spcs, sspcs;
  if (job_tag == "48I-0.00078" or job_tag == "64I-0.000678") {
    setPointCoordinatesSpheres(rng, pcs, 256, 4, total_site, sqr(6));
    pointInfoFromCoordinates(pis, pcs, 0, 0);
    spcs.clear();
    sampleFromCoordinates(rng, spcs, pcs.size() / 2, pcs);
    pointInfoFromCoordinates(pis, spcs, 1, 0);
    sspcs.clear();
    sampleFromCoordinates(rng, sspcs, pcs.size() / 32, spcs);
    pointInfoFromCoordinates(pis, sspcs, 0, 1);
    pointInfoFromCoordinates(pis, sspcs, 1, 1);
    spcs.clear();
    sampleFromCoordinates(rng, spcs, pcs.size() / 128, sspcs);
    pointInfoFromCoordinates(pis, spcs, 0, 2);
    pointInfoFromCoordinates(pis, spcs, 1, 2);
  } else if (job_tag == "16I-0.01") {
    setPointCoordinatesSpheres(rng, pcs, 16, 4, total_site, sqr(4));
    pointInfoFromCoordinates(pis, pcs, 0, 0);
    spcs.clear();
    sampleFromCoordinates(rng, spcs, pcs.size() / 2, pcs);
    pointInfoFromCoordinates(pis, spcs, 1, 0);
    sspcs.clear();
    sampleFromCoordinates(rng, sspcs, pcs.size() / 8, spcs);
    pointInfoFromCoordinates(pis, sspcs, 0, 1);
    pointInfoFromCoordinates(pis, sspcs, 1, 1);
  } else {
    setPointCoordinatesSpheres(rng, pcs, 16, 4, total_site, sqr(4));
    pointInfoFromCoordinates(pis, pcs, 0, 0);
    spcs.clear();
    sampleFromCoordinates(rng, spcs, pcs.size() / 2, pcs);
    pointInfoFromCoordinates(pis, spcs, 1, 0);
    sspcs.clear();
    sampleFromCoordinates(rng, sspcs, pcs.size() / 8, spcs);
    pointInfoFromCoordinates(pis, sspcs, 0, 1);
    pointInfoFromCoordinates(pis, sspcs, 1, 1);
  }
  sort(pis.begin(), pis.end(), pointInfoComp);
  return pis;
}

inline void save_lbl_pis(const std::vector<PointInfo>& pis,
                         const std::string& path)
// interface function
{
  TIMER_VERBOSE("save_lbl_pis");
  FILE* fp = qopen(path, "w");
  fprintf(fp, "%d\n", (int)pis.size());
  for (int i = 0; i < (int)pis.size(); ++i) {
    const PointInfo& pi = pis[i];
    const Coordinate& c = pi.xg;
    fprintf(fp, "%5d    %3d %3d %3d %3d    %2d %2d\n", i, c[0], c[1], c[2],
            c[3], pi.type, pi.accuracy);
  }
  qclose(fp);
}

inline void save_lbl_pis_info(const std::vector<PointInfo>& pis,
                              const std::string& path)
// interface function
{
  TIMER_VERBOSE("save_lbl_pis_info");
  if (0 == get_id_node()) {
    save_lbl_pis(pis, path);
  }
}

inline std::vector<PointInfo> load_lbl_pis(const std::string& path)
{
  TIMER_VERBOSE("load_lbl_pis");
  const std::vector<std::string> lines = qgetlines(path);
  qassert(lines.size() > 0);
  const long len = read_long(lines[0]);
  qassert(len + 1 <= (long)lines.size());
  std::vector<PointInfo> pis;
  for (int k = 1; k < len + 1; ++k) {
    const std::vector<std::string> strs = split_line_with_spaces(lines[k]);
    if (strs.size() >= 7) {
      qassert(k - 1 == read_long(strs[0]));
      const Coordinate xg(read_long(strs[1]), read_long(strs[2]),
                          read_long(strs[3]), read_long(strs[4]));
      const int type = read_long(strs[5]);
      const int acc = read_long(strs[6]);
      const PointInfo pi(xg, type, acc);
      pis.push_back(pi);
    } else {
      displayln(fname + ssprintf(": line is '%s'.", lines[k].c_str()));
      qassert(false);
    }
  }
  return pis;
}

inline std::vector<PointInfo> load_lbl_pis_info(const std::string& path)
{
  TIMER_VERBOSE("load_lbl_pis_info");
  std::vector<PointInfo> pis;
  if (0 == get_id_node()) {
    pis = load_lbl_pis(path);
  }
  bcast(pis);
  return pis;
}

inline std::vector<Coordinate> coordinates_from_point_infos(const std::vector<PointInfo>& pis)
// interface function
{
  TIMER("coordinates_from_point_infos");
  std::vector<Coordinate> ret;
  if (pis.size() == 0) {
    return ret;
  }
  ret.push_back(pis[0].xg);
  for (int i = 1; i < (int)pis.size(); ++i) {
    if (pis[i].xg != ret.back()) {
      ret.push_back(pis[i].xg);
    }
  }
  return ret;
}

inline std::vector<Coordinate> coordinates_from_point_infos_strange(
    const std::vector<PointInfo>& pis)
// interface function
{
  TIMER("coordinates_from_point_infos");
  std::vector<Coordinate> ret;
  for (int i = 0; i < (int)pis.size(); ++i) {
    if (pis[i].type == 1) {
      if (ret.size() == 0 or pis[i].xg != ret.back()) {
        ret.push_back(pis[i].xg);
      }
    }
  }
  return ret;
}

inline void save_lbl_pcs(const std::vector<Coordinate>& pcs,
                         const std::string& path)
// interface function
{
  TIMER_VERBOSE("save_lbl_pcs");
  FILE* fp = qopen(path, "w");
  fprintf(fp, "%d\n", (int)pcs.size());
  for (int i = 0; i < (int)pcs.size(); ++i) {
    const Coordinate& c = pcs[i];
    fprintf(fp, "%5d    %3d %3d %3d %3d\n", i, c[0], c[1], c[2], c[3]);
  }
  qclose(fp);
}

inline void save_lbl_pcs_info(const std::vector<Coordinate>& pcs,
                              const std::string& path)
// interface function
{
  TIMER_VERBOSE("save_lbl_pcs_info");
  if (0 == get_id_node()) {
    save_lbl_pcs(pcs, path);
  }
}

inline std::vector<Coordinate> load_lbl_pcs(const std::string& path)
{
  TIMER_VERBOSE("load_lbl_pcs");
  const std::vector<std::string> lines = qgetlines(path);
  qassert(lines.size() > 0);
  const long len = read_long(lines[0]);
  qassert(len + 1 <= (long)lines.size());
  std::vector<Coordinate> pcs;
  for (int k = 1; k < len + 1; ++k) {
    const std::vector<std::string> strs = split_line_with_spaces(lines[k]);
    if (strs.size() >= 5) {
      qassert(k - 1 == read_long(strs[0]));
      const Coordinate xg(read_long(strs[1]), read_long(strs[2]),
                          read_long(strs[3]), read_long(strs[4]));
      pcs.push_back(xg);
    } else {
      displayln(fname + ssprintf(": line is '%s'.", lines[k].c_str()));
      qassert(false);
    }
  }
  return pcs;
}

inline std::vector<Coordinate> load_lbl_pcs_info(const std::string& path)
{
  TIMER_VERBOSE("load_lbl_pcs_info");
  std::vector<Coordinate> pcs;
  if (0 == get_id_node()) {
    pcs = load_lbl_pcs(path);
  }
  bcast(pcs);
  return pcs;
}

typedef std::vector<std::vector<int> > TypeAccuracyTable;

inline TypeAccuracyTable mk_type_accuracy_table(
    const std::vector<PointInfo>& pis)
// tat[type][acc] = counts of this type and acc
{
  TIMER_VERBOSE("mk_type_accuracy_table");
  TypeAccuracyTable tat;
  for (int i = 0; i < (int)pis.size(); ++i) {
    const PointInfo& pi = pis[i];
    const int type = pi.type;
    const int acc = pi.accuracy;
    if (type >= (int)tat.size()) {
      tat.resize(type + 1);
    }
    if (acc >= (int)tat[type].size()) {
      tat[type].resize(acc + 1, 0);
    }
    tat[type][acc] += 1;
  }
  qassert(tat.size() == 2);
  qassert(tat[0].size() == tat[1].size());
  return tat;
}

inline double get_type_accuracy_weight(const TypeAccuracyTable& tat,
                                       const int type, const int accuracy)
{
  qassert(type < (int)tat.size());
  qassert(accuracy < (int)tat[type].size());
  return (double)tat[0][0] / (double)tat[type][accuracy];
}

inline double get_accuracy_weight(const TypeAccuracyTable& tat, const int type,
                                  const int accuracy)
{
  qassert(type < (int)tat.size());
  qassert(accuracy < (int)tat[type].size());
  return (double)tat[type][0] / (double)tat[type][accuracy];
}

}  // namespace qlat
