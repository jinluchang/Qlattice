#define QLAT_INSTANTIATE_SELECTED_POINTS

#include <qlat/selected-points.h>

namespace qlat
{  //

PointsSelection mk_tslice_points_selection(const Coordinate& total_site,
                                           const Int t_dir)
// psel.size() == total_site[t_dir]
{
  qassert(0 <= t_dir and t_dir < 4);
  const Int t_size = total_site[t_dir];
  PointsSelection psel;
  psel.init(total_site, t_size);
  const Coordinate xg_all = Coordinate(-1, -1, -1, -1);
  qthread_for(idx, t_size, {
    psel[idx] = xg_all;
    psel[idx][t_dir] = idx;
  });
  return psel;
}

PointsSelection mk_random_points_selection(const Coordinate& total_site,
                                           const Long num, const RngState& rs,
                                           const Long pool_factor)
// same rs for all node for uniform result
{
  TIMER_VERBOSE("mk_random_points_selection");
  if (num == 0) {
    PointsSelection psel;
    return psel;
  }
  qassert(num > 0);
  PointsSelection psel_pool(total_site, pool_factor * num);
#pragma omp parallel for
  for (Long i = 0; i < (Long)psel_pool.size(); ++i) {
    RngState rsi = rs.split(i);
    Coordinate xg;
    for (Int m = 0; m < 4; ++m) {
      xg[m] = modl(rand_gen(rsi), total_site[m]);
    }
    psel_pool[i] = xg;
  }
  PointsSelection psel(total_site, num);
  qthread_for(i, num, { psel[i] = Coordinate(-1, -1, -1, -1); });
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
  if (psel[num - 1] != Coordinate(-1, -1, -1, -1)) {
    return psel;
  } else {
    displayln_info(
        fname +
        ssprintf(": pool_factor=%d is too small, rerun with larger factor.",
                 pool_factor));
    return mk_random_points_selection(total_site, num, rs, pool_factor + 2);
  }
}

void set_psel_full(PointsSelection& psel, const Geometry& geo)
{
  TIMER("set_psel_full");
  const Coordinate total_site = geo.total_site();
  const Long n_points = geo.local_volume();
  psel.init(total_site, n_points, PointsDistType::Full);
  qacc_for(index, n_points, {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    psel[index] = xg;
  });
}

void lat_data_from_points_selection(LatDataInt& ld, const PointsSelection& psel)
{
  TIMER("lat_data_from_points_selection");
  Long n_points = psel.size();
  LatDim dim;
  dim.name = "idx";
  dim.size = n_points + 1;
  dim.indices.resize(1);
  dim.indices[0] = "total_site";
  ld.init();
  ld.info.push_back(dim);
  ld.info.push_back(lat_dim_number("mu", 0, 3));
  lat_data_alloc(ld);
  Vector<Int> v = lat_data_get(ld, make_array<Int>(0));
  for (Int mu = 0; mu < 4; ++mu) {
    v[mu] = psel.total_site[mu];
  }
  for (Long idx = 0; idx < (Long)psel.size(); ++idx) {
    const Coordinate& c = psel[idx];
    Vector<Int> v = lat_data_get(ld, make_array<Int>(idx + 1));
    for (Int mu = 0; mu < 4; ++mu) {
      v[mu] = c[mu];
    }
  }
}

void points_selection_from_lat_data(PointsSelection& psel, const LatDataInt& ld,
                                    const PointsDistType points_dist_type)
{
  TIMER("points_selection_from_lat_data");
  qassert(ld.info.size() == 2);
  qassert(ld.info[0].name == "idx");
  qassert(ld.info[1].name == "mu");
  qassert(ld.info[0].size >= 1);
  qassert(ld.info[1].size == 4);
  qassert(ld.info[0].indices[0] == "total_site");
  const Long n_points = ld.info[0].size - 1;
  Coordinate total_site;
  const Vector<Int> v = lat_data_get_const(ld, make_array<Int>(0));
  for (Int mu = 0; mu < 4; ++mu) {
    total_site[mu] = v[mu];
  }
  psel.init();
  psel.init(total_site, n_points, points_dist_type);
  for (Long idx = 0; idx < psel.size(); ++idx) {
    Coordinate& c = psel[idx];
    const Vector<Int> v = lat_data_get_const(ld, make_array<Int>(idx + 1));
    for (Int mu = 0; mu < 4; ++mu) {
      c[mu] = v[mu];
    }
  }
}

static void save_points_selection_txt(const PointsSelection& psel,
                                      const std::string& path)
{
  TIMER_VERBOSE("save_points_selection_txt");
  qassert(ends_with(path, ".txt"));
  qassert(psel.points_dist_type == PointsDistType::Global);
  QFile qfile = qfopen(path + ".partial", "w");
  qfprintf(qfile, "%ld\n", (Long)psel.size());
  for (Long i = 0; i < (Long)psel.size(); ++i) {
    const Coordinate& c = psel[i];
    qfprintf(qfile, "%5ld    %3d %3d %3d %3d\n", i, c[0], c[1], c[2], c[3]);
  }
  qfclose(qfile);
  qrename(path + ".partial", path);
}

void save_points_selection(const PointsSelection& psel, const std::string& path)
// IMPORTANT: first point saved in lati is the total_site
// path has to end with ".lati"
{
  TIMER_VERBOSE("save_points_selection");
  qassert(psel.points_dist_type == PointsDistType::Global);
  if (ends_with(path, ".txt")) {
    qwarn(fname + ssprintf(": path='%s' with old format. Need to set "
                           "total_site manually when reading.",
                           path.c_str()));
    save_points_selection_txt(psel, path);
    return;
  }
  qassert(ends_with(path, ".lati"));
  // save_points_selection_txt(psel, remove_suffix(path, ".lati") + ".txt");
  LatDataInt ld;
  lat_data_from_points_selection(ld, psel);
  ld.save(path);
}

void save_points_selection_info(const PointsSelection& psel,
                                const std::string& path)
{
  TIMER_VERBOSE("save_points_selection_info");
  qassert(psel.points_dist_type == PointsDistType::Global);
  if (0 == get_id_node()) {
    save_points_selection(psel, path);
  }
  SYNC_NODE();
}

static PointsSelection load_points_selection_txt(const std::string& path)
{
  TIMER_VERBOSE("load_points_selection");
  qassert(ends_with(path, ".txt"));
  qwarn(
      fname +
      ssprintf(": path='%s' with old format. Need to set total_site manually.",
               path.c_str()));
  const std::vector<std::string> lines = qgetlines(path);
  qassert(lines.size() > 0);
  const Long len = read_long(lines[0]);
  qassert(len + 1 <= (Long)lines.size());
  PointsSelection psel(Coordinate(), len);
  for (Long idx = 0; idx < len; ++idx) {
    const Long k = idx + 1;
    const std::vector<std::string> strs = split_line_with_spaces(lines[k]);
    if (strs.size() >= 5) {
      qassert(idx == read_long(strs[0]));
      const Coordinate xg(read_long(strs[1]), read_long(strs[2]),
                          read_long(strs[3]), read_long(strs[4]));
      psel[idx] = xg;
    } else {
      displayln(fname + ssprintf(": line is '%s'.", lines[k].c_str()));
      qassert(false);
    }
  }
  return psel;
}

PointsSelection load_points_selection(const std::string& path)
{
  TIMER_VERBOSE("load_points_selection");
  if (ends_with(path, ".txt")) {
    return load_points_selection_txt(path);
  }
  qassert(ends_with(path, ".lati"));
  LatDataInt ld;
  ld.load(path);
  PointsSelection psel;
  points_selection_from_lat_data(psel, ld);
  return psel;
}

PointsSelection load_points_selection_info(const std::string& path)
{
  TIMER_VERBOSE("load_points_selection_info");
  PointsSelection psel;
  if (0 == get_id_node()) {
    psel = load_points_selection(path);
  }
  bcast(psel);
  return psel;
}

crc32_t crc32_par(const PointsSelection& psel)
{
  PointsSelection psel_ec;
  psel_ec = psel;
  to_from_big_endian(get_data(psel_ec.xgs), true);
  return crc32_par(get_data(psel_ec.xgs));
}

void set_sqrt_field(SelectedPoints<RealD>& sp, const SelectedPoints<RealD>& sp1)
{
  TIMER("set_sqrt_field(sp,sp1)");
  const Long n_points = sp1.n_points;
  const Int multiplicity = sp1.multiplicity;
  sp.init(n_points, multiplicity, sp1.points_dist_type);
  qthread_for(idx, n_points, {
    const Vector<RealD> spv1 = sp1.get_elems_const(idx);
    Vector<RealD> spv = sp.get_elems(idx);
    for (Int m = 0; m < multiplicity; ++m) {
      spv[m] = std::sqrt(spv1[m]);
    }
  });
}

}  // namespace qlat
