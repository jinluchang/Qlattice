#define QLAT_INSTANTIATE_SELECTED_POINTS

#include <qlat/selected-points.h>

namespace qlat
{  //

PointSelection mk_tslice_point_selection(const int t_size, const int t_dir)
{
  PointSelection psel;
  psel.resize(t_size);
  qassert(0 <= t_dir and t_dir < 4);
  qthread_for(idx, t_size, {
    psel[idx] = Coordinate();
    psel[idx][t_dir] = idx;
  });
  return psel;
}

PointSelection mk_tslice_point_selection(const Coordinate& total_site,
                                         const int t_dir)
{
  return mk_tslice_point_selection(total_site[t_dir], t_dir);
}

PointSelection mk_random_point_selection(const Coordinate& total_site,
                                         const long num, const RngState& rs,
                                         const long pool_factor)
// same rs for all node for uniform result
{
  TIMER_VERBOSE("mk_random_point_selection");
  if (num == 0) {
    PointSelection psel;
    return psel;
  }
  qassert(num > 0);
  PointSelection psel_pool(pool_factor * num);
#pragma omp parallel for
  for (long i = 0; i < (long)psel_pool.size(); ++i) {
    RngState rsi = rs.split(i);
    Coordinate xg;
    for (int m = 0; m < 4; ++m) {
      xg[m] = modl(rand_gen(rsi), total_site[m]);
    }
    psel_pool[i] = xg;
  }
  PointSelection psel(num, Coordinate(-1, -1, -1, -1));
  long idx = 0;
  for (long i = 0; i < (long)psel.size(); ++i) {
    while (idx < (long)psel_pool.size()) {
      const Coordinate xg = psel_pool[idx];
      idx += 1;
      bool is_repeat = false;
      for (long j = 0; j < i; ++j) {
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

void save_point_selection(const PointSelection& psel, const std::string& path)
{
  TIMER_VERBOSE("save_point_selection");
  QFile qfile = qfopen(path + ".partial", "w");
  qfprintf(qfile, "%ld\n", (long)psel.size());
  for (long i = 0; i < (long)psel.size(); ++i) {
    const Coordinate& c = psel[i];
    qfprintf(qfile, "%5ld    %3d %3d %3d %3d\n", i, c[0], c[1], c[2], c[3]);
  }
  qfclose(qfile);
  qrename(path + ".partial", path);
}

void save_point_selection_info(const PointSelection& psel,
                               const std::string& path)
{
  TIMER_VERBOSE("save_point_selection_info");
  if (0 == get_id_node()) {
    save_point_selection(psel, path);
  }
}

PointSelection load_point_selection(const std::string& path)
{
  TIMER_VERBOSE("load_point_selection");
  const std::vector<std::string> lines = qgetlines(path);
  qassert(lines.size() > 0);
  const long len = read_long(lines[0]);
  qassert(len + 1 <= (long)lines.size());
  PointSelection psel;
  for (long k = 1; k < len + 1; ++k) {
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

PointSelection load_point_selection_info(const std::string& path)
{
  TIMER_VERBOSE("load_point_selection_info");
  PointSelection psel;
  if (0 == get_id_node()) {
    psel = load_point_selection(path);
  }
  bcast(psel);
  return psel;
}

}  // namespace qlat
