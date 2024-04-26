#define QLAT_INSTANTIATE_CORE

#include <qlat/core.h>

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

void Geometry::init(const Coordinate& total_site, const int multiplicity_)
{
  if (!initialized) {
    init();
    geon = get_geometry_node();
    multiplicity = multiplicity_;
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
  s += ssprintf(", multiplicity  = %s\n", show(geo.multiplicity).c_str());
  s += ssprintf(", node_site     = %s\n", show(geo.node_site).c_str());
  s += ssprintf(", expan_left    = %s\n", show(geo.expansion_left).c_str());
  s += ssprintf(", expan_right   = %s\n", show(geo.expansion_right).c_str());
  s += ssprintf(", node_site_exp = %s\n", show(geo.node_site_expanded).c_str());
  s += ssprintf(", is_only_local = %s }", show(geo.is_only_local).c_str());
  return s;
}

void PointsSelection::init()
{
  initialized = false;
  distributed = false;
  xgs.init();
}

void PointsSelection::init(const Long n_points)
{
  initialized = true;
  distributed = false;
  xgs.clear();
  xgs.resize(n_points);
}

void PointsSelection::init(const Long n_points, const Coordinate& xg_init)
{
  initialized = true;
  distributed = false;
  xgs.clear();
  xgs.resize(n_points, xg_init);
}

void PointsSelection::init(const std::vector<Coordinate>& xgs_)
{
  initialized = true;
  distributed = false;
  xgs = xgs_;
}

void PointsSelection::init(const vector<Coordinate>& xgs_)
{
  initialized = true;
  distributed = false;
  xgs = xgs_;
}

void PointsSelection::init(const SelectedPoints<Coordinate>& spx)
{
  qassert(spx.multiplicity == 1);
  initialized = spx.initialized;
  distributed = spx.distributed;
  xgs = spx.points;
}

PointsSelection& PointsSelection::operator=(const std::vector<Coordinate>& xgs_)
{
  init(xgs_);
  return *this;
}

PointsSelection& PointsSelection::operator=(const vector<Coordinate>& xgs_)
{
  init(xgs_);
  return *this;
}

PointsSelection& PointsSelection::operator=(
    const SelectedPoints<Coordinate>& spx)
{
  init(spx);
  return *this;
}

void PointsSelection::resize(const Long n_points) { xgs.resize(n_points); }

void PointsSelection::resize(const Long n_points, const Coordinate& xg_init)
{
  xgs.resize(n_points, xg_init);
}

SelectedPoints<Coordinate> PointsSelection::view_sp()
{
  TIMER("PointsSelection::view_sp");
  SelectedPoints<Coordinate> f;
  f.initialized = initialized;
  f.distributed = distributed;
  f.multiplicity = 1;
  f.n_points = xgs.size();
  f.points = xgs.view();
  return f;
}

void PointsSelection::push_back_slow(const Coordinate& xg)
{
  const Long n_points = xgs.size();
  xgs.resize(n_points + 1);
  xgs[n_points] = xg;
}

bool operator==(const PointsSelection& psel1, const PointsSelection& psel2)
{
  if (psel1.initialized != psel2.initialized) {
    return false;
  }
  if (psel1.distributed != psel2.distributed) {
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

void FieldSelection::init()
{
  f_rank.init();
  f_local_idx.init();
  n_elems = 0;
  ranks.init();
  indices.init();
}

}  // namespace qlat
