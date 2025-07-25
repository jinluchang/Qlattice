#define QLAT_INSTANTIATE_CORE

#include <qlat/core.h>
#include <qlat/selected-field.h>

namespace qlat
{  //

std::string show(const qlat::GeometryNode& geon)
{
  std::string s;
  s += ssprintf("{ initialized = %s\n", show(geon.initialized).c_str());
  s += ssprintf(", num_node    = %d\n", geon.num_node);
  s += ssprintf(", id_node     = %d\n", geon.id_node);
  s += ssprintf(", size_node   = %s\n", show(geon.size_node).c_str());
  s += ssprintf(", coor_node   = %s }", show(geon.coor_node).c_str());
  return s;
}

void Geometry::init(const Coordinate& total_site)
{
  if (!initialized) {
    init();
    geon = get_geometry_node();
    for (int i = 0; i < DIMN; ++i) {
      qassert(0 == total_site[i] % geon.size_node[i]);
      node_site[i] = total_site[i] / geon.size_node[i];
    }
    reset_node_site_expanded();
    initialized = true;
  }
}

std::string show(const qlat::Geometry& geo)
{
  std::string s;
  s += ssprintf("{ initialized   = %s\n", show(geo.initialized).c_str());
  s += ssprintf(", eo            = %s\n", show(geo.eo).c_str());
  s += ssprintf(", geon          =\n%s\n", show(geo.geon).c_str());
  s += ssprintf(", node_site     = %s\n", show(geo.node_site).c_str());
  s += ssprintf(", expan_left    = %s\n", show(geo.expansion_left).c_str());
  s += ssprintf(", expan_right   = %s\n", show(geo.expansion_right).c_str());
  s += ssprintf(", node_site_exp = %s\n", show(geo.node_site_expanded).c_str());
  s += ssprintf(", is_only_local = %s }", show(geo.is_only_local).c_str());
  return s;
}

// ----------------

std::string show(const PointsDistType points_dist_type)
{
  if (points_dist_type == PointsDistType::Global) {
    return "g";
  } else if (points_dist_type == PointsDistType::Full) {
    return "f";
  } else if (points_dist_type == PointsDistType::Local) {
    return "l";
  } else if (points_dist_type == PointsDistType::Random) {
    return "r";
  } else if (points_dist_type == PointsDistType::Other) {
    return "o";
  } else {
    qassert(false);
    return "";
  }
}

PointsDistType read_points_dist_type(const std::string& points_dist_type_str)
{
  if (points_dist_type_str == "g") {
    return PointsDistType::Global;
  } else if (points_dist_type_str == "f") {
    return PointsDistType::Full;
  } else if (points_dist_type_str == "l") {
    return PointsDistType::Local;
  } else if (points_dist_type_str == "r") {
    return PointsDistType::Random;
  } else if (points_dist_type_str == "o") {
    return PointsDistType::Other;
  } else {
    qassert(false);
    return PointsDistType::Global;
  }
}

void PointsSelection::init()
{
  initialized = false;
  points_dist_type = PointsDistType::Global;
  xgs.init();
}

void PointsSelection::init(const Coordinate& total_site_, const Long n_points_,
                           const PointsDistType points_dist_type_)
{
  initialized = true;
  points_dist_type = points_dist_type_;
  total_site = total_site_;
  xgs.clear();
  xgs.resize(n_points_);
}

void PointsSelection::init(const Coordinate& total_site_,
                           const std::vector<Coordinate>& xgs_)
{
  initialized = true;
  points_dist_type = PointsDistType::Global;
  total_site = total_site_;
  xgs = xgs_;
}

void PointsSelection::init(const Coordinate& total_site_,
                           const vector<Coordinate>& xgs_)
{
  initialized = true;
  points_dist_type = PointsDistType::Global;
  total_site = total_site_;
  xgs = xgs_;
}

void PointsSelection::init(const Coordinate& total_site_,
                           const SelectedPoints<Coordinate>& spx)
{
  qassert(spx.multiplicity == 1);
  initialized = spx.initialized;
  points_dist_type = spx.points_dist_type;
  total_site = total_site_;
  xgs = spx.points;
}

void PointsSelection::resize(const Long n_points) { xgs.resize(n_points); }

SelectedPoints<Coordinate> PointsSelection::view_sp() const
{
  TIMER("PointsSelection::view_sp");
  SelectedPoints<Coordinate> f;
  f.initialized = initialized;
  f.points_dist_type = points_dist_type;
  f.multiplicity = 1;
  f.n_points = xgs.size();
  f.points.set_view(xgs);
  return f;
}

void PointsSelection::push_back_slow(const Coordinate& xg)
{
  const Long n_points = xgs.size();
  xgs.resize(n_points + 1);
  xgs[n_points] = xg;
}

void qswap(PointsSelection& f1, PointsSelection& f2)
{
  std::swap(f1.initialized, f2.initialized);
  std::swap(f1.points_dist_type, f2.points_dist_type);
  std::swap(f1.total_site, f2.total_site);
  qswap(f1.xgs, f2.xgs);
}

bool operator==(const PointsSelection& psel1, const PointsSelection& psel2)
{
  if (psel1.initialized != psel2.initialized) {
    return false;
  }
  if (psel1.points_dist_type != psel2.points_dist_type) {
    return false;
  }
  if (psel1.total_site != psel2.total_site) {
    return false;
  }
  if (psel1.xgs.size() != psel2.xgs.size()) {
    return false;
  }
  for (Long i = 0; i < psel1.xgs.size(); ++i) {
    if (psel1.xgs[i] != psel2.xgs[i]) {
      return false;
    }
  }
  return true;
}

bool operator!=(const PointsSelection& psel1, const PointsSelection& psel2)
{
  return not(psel1 == psel2);
}

// ----------------

void FieldSelection::init()
{
  f_rank.init();
  f_local_idx.init();
  n_elems = 0;
  ranks.init();
  indices.init();
}

void set_psel_from_fsel(PointsSelection& psel, const FieldSelection& fsel)
// psel.points_dist_type will be PointsDistType::Local
{
  TIMER("set_psel_from_fsel(psel,fsel)");
  qassert(fsel.f_rank.initialized);
  const Geometry& geo = fsel.f_rank.geo();
  const Coordinate total_site = geo.total_site();
  const Long n_points = fsel.n_elems;
  psel.init(total_site, n_points);
  psel.points_dist_type = PointsDistType::Local;
  qthread_for(idx, n_points, {
    const Long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    psel[idx] = xg;
  });
}

void set_fsel_from_psel(FieldSelection& fsel, const PointsSelection& psel,
                        const Geometry& geo, const Long rank_psel)
// psel.points_dist_type can be either
// PointsDistType::Local
// or PointsDistType::Full
// or PointsDistType::Global
{
  TIMER("set_fsel_from_psel(fsel,psel,geo,rank_psel)");
  qassert(psel.initialized);
  qassert((psel.points_dist_type == PointsDistType::Local) or
          (psel.points_dist_type == PointsDistType::Full) or
          (psel.points_dist_type == PointsDistType::Global));
  qassert(psel.total_site == geo.total_site());
  fsel.init();
  mk_field_selection(fsel.f_rank, geo, -1);
  qthread_for(idx, psel.size(), {
    const Coordinate xg = psel[idx];
    const Coordinate xl = geo.coordinate_l_from_g(xg);
    if (geo.is_local(xl)) {
      fsel.f_rank.get_elem(xl) = rank_psel;
    } else {
      qassert(psel.points_dist_type == PointsDistType::Global);
    }
  });
  update_field_selection(fsel);
  if ((psel.points_dist_type == PointsDistType::Local) or
      (psel.points_dist_type == PointsDistType::Full)) {
    qassert(fsel.n_elems == psel.size());
  }
}

}  // namespace qlat
