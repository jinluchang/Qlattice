#pragma once

#include "psrc-sample.h"

namespace qlat
{  //

typedef std::vector<double> PointDistribution;

inline Coordinate normalize_coordinate(const Coordinate& c)
{
  const int x = std::abs(c[0]);
  const int y = std::abs(c[1]);
  const int z = std::abs(c[2]);
  const int t = std::abs(c[3]);
  if (x >= y and y >= z) {
    return Coordinate(x, y, z, t);
  } else if (x >= z and z >= y) {
    return Coordinate(x, z, y, t);
  } else if (y >= z and z >= x) {
    return Coordinate(y, z, x, t);
  } else if (y >= x and x >= z) {
    return Coordinate(y, x, z, t);
  } else if (z >= x and x >= y) {
    return Coordinate(z, x, y, t);
  } else if (z >= y and y >= x) {
    return Coordinate(z, y, x, t);
  } else {
    qassert(false);
    return c;
  }
}

inline long index_from_relative_coordinate(const Coordinate& rel,
                                           const Coordinate& total_site)
{
  const Coordinate limit = total_site / 2 + Coordinate(1, 1, 1, 1);
  const Coordinate c = normalize_coordinate(rel);
  return index_from_coordinate(c, limit);
}

inline Coordinate normalized_relative_coordinate_from_index(
    const long index, const Coordinate& total_site)
{
  const Coordinate limit = total_site / 2 + Coordinate(1, 1, 1, 1);
  return coordinate_from_index(index, limit);
}

inline void init_pd(PointDistribution& pd, const Coordinate& total_site)
{
  const Coordinate limit = total_site / 2 + Coordinate(1, 1, 1, 1);
  clear(pd);
  pd.resize(product(limit), 0.0);
}

inline void normalize_pd(PointDistribution& pd, const Coordinate& total_site)
// normalize probability (under the condition that distance is non-zero)
// NOTE: call "renormalize_pd" if you already have a prob distribution. This
// function ONLY works for counts as input
{
  TIMER_VERBOSE("normalize_pd");
  double count = 0;
  for (long index = 1; index < (long)pd.size(); ++index) {
    count += pd[index];
  }
  PointDistribution pdm;
  init_pd(pdm, total_site);
  for (long index = 0; index < product(total_site); ++index) {
    const Coordinate rel = relative_coordinate(
        coordinate_from_index(index, total_site), total_site);
    const long idx = index_from_relative_coordinate(rel, total_site);
    qassert(0 <= idx and idx < (long)pdm.size());
    pdm[idx] += 1;
  }
  qassert(pd.size() == pdm.size());
#pragma omp parallel for
  for (long index = 1; index < (long)pd.size(); ++index) {
    if (pdm[index] > 0) {
      pd[index] /= count * pdm[index];
    }
  }
  pd[0] = 1.0;
}

inline void renormalize_pd(PointDistribution& pd, const Coordinate& total_site)
// renormalize probability (under the condition that distance is non-zero)
{
  TIMER_VERBOSE("renormalize_pd");
  PointDistribution pdm;
  init_pd(pdm, total_site);
  for (long index = 0; index < product(total_site); ++index) {
    const Coordinate rel = relative_coordinate(
        coordinate_from_index(index, total_site), total_site);
    const long idx = index_from_relative_coordinate(rel, total_site);
    qassert(0 <= idx and idx < (long)pdm.size());
    pdm[idx] += 1;
  }
  qassert(pd.size() == pdm.size());
  double prob_sum = 0;
  for (long index = 1; index < (long)pd.size(); ++index) {
    prob_sum += pd[index] * pdm[index];
  }
#pragma omp parallel for
  for (long index = 1; index < (long)pd.size(); ++index) {
    if (pdm[index] > 0) {
      pd[index] /= prob_sum;
    }
  }
  pd[0] = 1.0;
}

inline void save_pd(const PointDistribution& pd, const Coordinate& total_site,
                    const std::string& path)
{
  TIMER_VERBOSE("save_pd");
  if (get_id_node() == 0) {
    FILE* fp = qopen(path + ".partial", "w");
    qassert(fp != NULL);
    for (long index = 0; index < (long)pd.size(); ++index) {
      const Coordinate rel =
          normalized_relative_coordinate_from_index(index, total_site);
      if (rel == normalize_coordinate(rel)) {
        const double v =
            pd[index];  // random select two different points, the probability
                        // that the relative coordinate is equal to ``rel''
        displayln(ssprintf("%4d %4d %4d %4d   %24.17E", rel[0], rel[1], rel[2],
                           rel[3], v),
                  fp);
      } else {
        qassert(pd[index] == 0.0);
      }
    }
    qclose(fp);
    qrename(path + ".partial", path);
  }
}

inline void load_pd(PointDistribution& pd, const Coordinate& total_site,
                    const std::string& path)
{
  TIMER_VERBOSE("load_pd");
  init_pd(pd, total_site);
  if (get_id_node() == 0) {
    FILE* fp = qopen(path, "r");
    qassert(fp != NULL);
    for (long index = 0; index < (long)pd.size(); ++index) {
      const Coordinate rel =
          normalized_relative_coordinate_from_index(index, total_site);
      if (rel == normalize_coordinate(rel)) {
        Coordinate read_rel;
        double& v = pd[index];
        fscanf(fp, "%d %d %d %d %lf", &read_rel[0], &read_rel[1], &read_rel[2],
               &read_rel[3], &v);
        qassert(rel == read_rel);
      } else {
        qassert(pd[index] == 0.0);
      }
    }
    qclose(fp);
  }
  bcast(get_data(pd));
}

inline double weight_from_pd(const PointDistribution& pd,
                             const Coordinate& total_site, const Coordinate& x,
                             const Coordinate& y)
// interface function
{
  return 1.0 / pd[index_from_relative_coordinate(
                   relative_coordinate(x - y, total_site), total_site)];
}

inline void acc_pd_with_pis(PointDistribution& pd,
                            const std::vector<PointInfo>& pis,
                            const Coordinate& total_site)
{
  TIMER_VERBOSE("acc_pd_with_pis");
  for (int i = 0; i < (int)pis.size(); ++i) {
    const PointInfo& pi1 = pis[i];
    if (pi1.type != 0 or pi1.accuracy != 0) {
      continue;
    }
    for (int j = 0; j < (int)pis.size(); ++j) {
      const PointInfo& pi2 = pis[j];
      if (pi2.type != 0 or pi2.accuracy != 0) {
        continue;
      }
      const Coordinate rel = relative_coordinate(pi1.xg - pi2.xg, total_site);
      pd[index_from_relative_coordinate(rel, total_site)] += 1;
    }
  }
}

inline PointDistribution compute_point_distribution(const std::string& job_tag)
// interface function
{
  const Coordinate total_site = get_total_site(job_tag);
  const std::string path = get_job_path(job_tag) + "/point-distribution";
  PointDistribution pd;
  init_pd(pd, total_site);
  if (does_file_exist_sync_node(path + ".txt")) {
    load_pd(pd, total_site, path + ".txt");
    return pd;
  }
  const std::string lock = path + "-lock";
  if (!obtain_lock(lock)) {
    qassert(false);
  }
  qmkdir_sync_node(path + ".dir");
  const int n_sims = 16 * 1024;
  for (int traj = 0; traj < n_sims; ++traj) {
    if (traj % get_num_node() == get_id_node()) {
      const std::vector<PointInfo> pis = mk_point_infos(job_tag, traj);
      acc_pd_with_pis(pd, pis, total_site);
      save_lbl_pis(pis, path + ".dir/traj=" + show(traj) + " ; pis.txt");
    }
  }
  glb_sum(get_data(pd));
  normalize_pd(pd, total_site);
  save_pd(pd, total_site, path + ".txt");
  sync_node();
  release_lock();
  Timer::display();
  return pd;
}

}  // namespace qlat
