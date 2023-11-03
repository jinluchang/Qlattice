#define QLAT_INSTANTIATE_SELECTED_POINTS

#include <qlat/selected-points.h>

namespace qlat
{  //

PointsSelection mk_tslice_point_selection(const int t_size, const int t_dir)
{
  PointsSelection psel;
  psel.resize(t_size);
  qassert(0 <= t_dir and t_dir < 4);
  const Coordinate xg_all = Coordinate(-1, -1, -1, -1);
  qthread_for(idx, t_size, {
    psel[idx] = xg_all;
    psel[idx][t_dir] = idx;
  });
  return psel;
}

PointsSelection mk_tslice_point_selection(const Coordinate& total_site,
                                          const int t_dir)
{
  return mk_tslice_point_selection(total_site[t_dir], t_dir);
}

PointsSelection mk_random_point_selection(const Coordinate& total_site,
                                         const Long num, const RngState& rs,
                                         const Long pool_factor)
// same rs for all node for uniform result
{
  TIMER_VERBOSE("mk_random_point_selection");
  if (num == 0) {
    PointsSelection psel;
    return psel;
  }
  qassert(num > 0);
  PointsSelection psel_pool(pool_factor * num);
#pragma omp parallel for
  for (Long i = 0; i < (Long)psel_pool.size(); ++i) {
    RngState rsi = rs.split(i);
    Coordinate xg;
    for (int m = 0; m < 4; ++m) {
      xg[m] = modl(rand_gen(rsi), total_site[m]);
    }
    psel_pool[i] = xg;
  }
  PointsSelection psel(num, Coordinate(-1, -1, -1, -1));
  Long idx = 0;
  for (Long i = 0; i < (Long)psel.size(); ++i) {
    while (idx < (Long)psel_pool.size()) {
      const Coordinate xg = psel_pool[idx];
      idx += 1;
      bool is_repeat = false;
      for (Long j = 0; j < i; ++j) {
        if (xg == psel[j]) {
          is_repeat = true;
          break;
        }
      }
      if (not is_repeat) {
        psel[i] = xg;
        break;
      }
    }
  }
  if (psel.back() != Coordinate(-1, -1, -1, -1)) {
    return psel;
  } else {
    displayln_info(
        fname +
        ssprintf(": pool_factor=%d is too small, rerun with larger factor.",
                 pool_factor));
    return mk_random_point_selection(total_site, num, rs, pool_factor + 2);
  }
}

void save_point_selection(const PointsSelection& psel, const std::string& path)
{
  TIMER_VERBOSE("save_point_selection");
  QFile qfile = qfopen(path + ".partial", "w");
  qfprintf(qfile, "%ld\n", (Long)psel.size());
  for (Long i = 0; i < (Long)psel.size(); ++i) {
    const Coordinate& c = psel[i];
    qfprintf(qfile, "%5ld    %3d %3d %3d %3d\n", i, c[0], c[1], c[2], c[3]);
  }
  qfclose(qfile);
  qrename(path + ".partial", path);
}

void save_point_selection_info(const PointsSelection& psel,
                               const std::string& path)
{
  TIMER_VERBOSE("save_point_selection_info");
  if (0 == get_id_node()) {
    save_point_selection(psel, path);
  }
}

PointsSelection load_point_selection(const std::string& path)
{
  TIMER_VERBOSE("load_point_selection");
  const std::vector<std::string> lines = qgetlines(path);
  qassert(lines.size() > 0);
  const Long len = read_long(lines[0]);
  qassert(len + 1 <= (Long)lines.size());
  PointsSelection psel;
  for (Long k = 1; k < len + 1; ++k) {
    const std::vector<std::string> strs = split_line_with_spaces(lines[k]);
    if (strs.size() >= 5) {
      qassert(k - 1 == read_long(strs[0]));
      const Coordinate xg(read_long(strs[1]), read_long(strs[2]),
                          read_long(strs[3]), read_long(strs[4]));
      psel.push_back(xg);
    } else {
      displayln(fname + ssprintf(": line is '%s'.", lines[k].c_str()));
      qassert(false);
    }
  }
  return psel;
}

PointsSelection load_point_selection_info(const std::string& path)
{
  TIMER_VERBOSE("load_point_selection_info");
  PointsSelection psel;
  if (0 == get_id_node()) {
    psel = load_point_selection(path);
  }
  bcast(psel);
  return psel;
}

void set_sqrt_field(SelectedPoints<RealD>& sp, const SelectedPoints<RealD>& sp1)
{
  TIMER("set_sqrt_field(sp,sp1)");
  const Long n_points = sp1.n_points;
  const Int multiplicity = sp1.multiplicity;
  sp.init(n_points, multiplicity);
  qthread_for(idx, n_points, {
    const Vector<RealD> spv1 = sp1.get_elems_const(idx);
    Vector<RealD> spv = sp.get_elems(idx);
    for (Int m = 0; m < multiplicity; ++m) {
      spv[m] = std::sqrt(spv1[m]);
    }
  });
}

}  // namespace qlat
