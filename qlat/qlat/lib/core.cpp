#define QLAT_INSTANTIATE_CORE

#include <qlat/core.h>
#include <qlat/selected-field.h>

namespace qlat
{  //

qacc_no_inline GeometryNode::GeometryNode() { init(); }

qacc_no_inline GeometryNode::GeometryNode(const Coordinate& coor_node_,
                                          const Coordinate& size_node_)
{
  init(coor_node_, size_node_);
}

qacc_no_inline GeometryNode::GeometryNode(const int id_node_, const Coordinate& size_node_)
{
  init(id_node_, size_node_);
}

qacc_no_inline void GeometryNode::init()
{
  initialized = false;
  num_node = 0;
  id_node = 0;
  size_node.init();
  coor_node.init();
}

qacc_no_inline void GeometryNode::init(const Coordinate& coor_node_,
                                       const Coordinate& size_node_)
{
  const int id_node_ = index_from_coordinate(coor_node_, size_node_);
  const int num_node_ = product(size_node_);
  initialized = true;
  num_node = num_node_;
  id_node = id_node_;
  size_node = size_node_;
  coor_node = coor_node_;
}

qacc_no_inline void GeometryNode::init(const Int id_node_,
                                       const Coordinate& size_node_)
{
  const Coordinate coor_node_ = coordinate_from_index(id_node_, size_node_);
  init(coor_node_, size_node_);
}

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

qacc_no_inline bool operator==(const GeometryNode& geon1,
                               const GeometryNode& geon2)
{
  return geon1.initialized == geon2.initialized &&
         geon1.num_node == geon2.num_node && geon1.id_node == geon2.id_node &&
         geon1.size_node == geon2.size_node &&
         geon1.coor_node == geon2.coor_node;
}

qacc_no_inline bool operator!=(const GeometryNode& geon1,
                               const GeometryNode& geon2)
{
  return !(geon1 == geon2);
}

qacc_no_inline Geometry::Geometry() { init(); }

Geometry::Geometry(const Coordinate& total_site) { init(total_site); }

qacc_no_inline void Geometry::init()
{
  initialized = false;
  geon.init();
  eo = 0;
  node_site.init();
  expansion_left.init();
  expansion_right.init();
  node_site_expanded.init();
  is_only_local = false;
}

qacc_no_inline void Geometry::init(const GeometryNode& geon_,
                                   const Coordinate& node_site_)
{
  init();
  geon = geon_;
  node_site = node_site_;
  reset_node_site_expanded();
  initialized = true;
}

qacc_no_inline void Geometry::init(const Coordinate& coor_node_,
                                   const Coordinate& size_node_,
                                   const Coordinate& node_site_)
{
  GeometryNode geon;
  geon.init(coor_node_, size_node_);
  init(geon, node_site_);
}

qacc_no_inline void Geometry::init(const Int id_node_,
                                   const Coordinate& size_node_,
                                   const Coordinate& node_site_)
{
  GeometryNode geon;
  geon.init(id_node_, size_node_);
  init(geon, node_site_);
}

void Geometry::init(const Coordinate& total_site)
{
  const GeometryNode& geon_ = get_geometry_node();
  Coordinate node_site_;
  for (int i = 0; i < DIMN; ++i) {
    qassert(0 == total_site[i] % geon_.size_node[i]);
    node_site_[i] = total_site[i] / geon_.size_node[i];
  }
  init(geon_, node_site_);
}

qacc_no_inline void Geometry::reset_node_site_expanded()
{
  is_only_local = true;
  for (int i = 0; i < DIMN; ++i) {
    node_site_expanded[i] =
        expansion_left[i] + node_site[i] + expansion_right[i];
    if (expansion_left[i] != 0 or expansion_right[i] != 0) {
      is_only_local = false;
    }
  }
}

qacc_no_inline void Geometry::resize(const Coordinate& expansion_left_,
                                     const Coordinate& expansion_right_)
{
  expansion_left = expansion_left_;
  expansion_right = expansion_right_;
  for (int i = 0; i < DIMN; ++i) {
    if (geon.size_node[i] == 1) {
      expansion_left[i] = 0;
      expansion_right[i] = 0;
    }
  }
  reset_node_site_expanded();
}

qacc_no_inline void Geometry::resize(const int thick)
{
  const Coordinate expansion(thick, thick, thick, thick);
  resize(expansion, expansion);
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

qacc_no_inline bool operator==(const Geometry& geo1, const Geometry& geo2)
{
  return geo1.initialized == geo2.initialized && geo1.eo == geo2.eo &&
         geo1.geon == geo2.geon && geo1.node_site == geo2.node_site &&
         geo1.expansion_left == geo2.expansion_left &&
         geo1.expansion_right == geo2.expansion_right &&
         geo1.node_site_expanded == geo2.node_site_expanded &&
         geo1.is_only_local == geo2.is_only_local;
}

qacc_no_inline bool operator!=(const Geometry& geo1, const Geometry& geo2)
{
  return !(geo1 == geo2);
}

qacc_no_inline Geometry geo_resize(const Geometry& geo_, const int thick)
{
  Geometry geo = geo_;
  geo.resize(thick);
  return geo;
}

qacc_no_inline Geometry geo_resize(const Geometry& geo_,
                                   const Coordinate& expansion_left_,
                                   const Coordinate& expansion_right_)
{
  Geometry geo = geo_;
  geo.resize(expansion_left_, expansion_right_);
  return geo;
}

qacc_no_inline Geometry geo_eo(const Geometry& geo_, const int eo_)
// 0:regular; 1:odd; 2:even
{
  Geometry geo = geo_;
  geo.eo = eo_;
  return geo;
}

qacc_no_inline bool is_matching_geo(const Geometry& geo1, const Geometry& geo2)
{
  return geo1.initialized == geo2.initialized && geo1.geon == geo2.geon &&
         geo1.node_site == geo2.node_site;
}

qacc_no_inline bool is_matching_geo_included(const Geometry& geo1,
                                             const Geometry& geo2)
// return if geo1 is included in geo2
{
  bool include = is_matching_geo(geo1, geo2);
  for (int i = 0; i < 4; i++) {
    if (geo2.expansion_left[i] < geo1.expansion_left[i]) {
      include = false;
    }
  }
  for (int i = 0; i < 4; i++) {
    if (geo2.expansion_right[i] < geo1.expansion_right[i]) {
      include = false;
    }
  }
  return include;
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

void PointsSelection::set_view(const PointsSelection& psel)
{
  TIMER("PointsSelection::set_view");
  initialized = psel.initialized;
  points_dist_type = psel.points_dist_type;
  total_site = psel.total_site;
  xgs.set_view(psel.xgs);
}

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

void set_field_m(Field<Char>& f, const Field<Char>& f1, const Int m,
                 const Int m1, const Int sizeof_m)
// `m` and `m1` have NOT be multiplied by `sizeof_m` yet.
{
  TIMER_FLOPS("set_field_m(f,f1,m,m1,sizeof_m)");
  qassert(f.initialized);
  qassert(f1.initialized);
  const Geometry geo = f.get_geo();
  qassert(geo.is_only_local);
  qassert(geo == f1.get_geo());
  const Int multiplicity = f.multiplicity;
  const Int multiplicity1 = f1.multiplicity;
  qassert(multiplicity % sizeof_m == 0);
  qassert(multiplicity1 % sizeof_m == 0);
  const Int m_c = m * sizeof_m;
  const Int m1_c = m1 * sizeof_m;
  const Long local_volume = geo.local_volume();
  qacc_for(index, local_volume, {
    const Vector<Char> v1 = f1.get_elems_const(index);
    const Vector<Char> v1s(v1.p + m1_c, sizeof_m);
    Vector<Char> v = f.get_elems(index);
    Vector<Char> vs(v.p + m_c, sizeof_m);
    assign(vs, v1s);
  });
  timer.flops += local_volume * sizeof_m;
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

qacc_no_inline Geometry FieldSelection::get_geo() const
{
  return f_rank.geo.get();
}

void FieldSelection::set_mem_type(const MemType mem_type) const
{
  f_rank.set_mem_type(mem_type);
  f_local_idx.set_mem_type(mem_type);
  ranks.set_mem_type(mem_type);
  indices.set_mem_type(mem_type);
}

void FieldSelection::set_view(const FieldSelection& fsel)
{
  f_rank.set_view(fsel.f_rank);
  n_elems = fsel.n_elems;
  f_local_idx.set_view(fsel.f_local_idx);
  ranks.set_view(fsel.ranks);
  indices.set_view(fsel.indices);
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

void set_geo_from_psel(Geometry& geo, const PointsSelection& psel)
{
  TIMER("set_geo_from_psel");
  qassert(psel.points_dist_type == PointsDistType::Full);
  psel.set_mem_type(MemType::Cpu);
  const Coordinate total_site = psel.total_site;
  const Long n_points = psel.size();
  qassert(n_points > 0);
  Coordinate left = psel.xgs[0];
  Coordinate right = psel.xgs[n_points - 1];
  qfor(idx, n_points, {
    const Coordinate& xg = psel.xgs[idx];
    for (Int m = 0; m < 4; ++m) {
      if (xg[m] < left[m]) {
        left[m] = xg[m];
      }
      if (xg[m] > right[m]) {
        right[m] = xg[m];
      }
    }
  });
  Coordinate node_site, size_node, coor_node;
  node_site = right - left;
  for (Int m = 0; m < 4; ++m) {
    node_site[m] = right[m] - left[m] + 1;
    size_node[m] = total_site[m] / node_site[m];
    coor_node[m] = left[m] / node_site[m];
    qassert(size_node[m] * node_site[m] == total_site[m]);
    qassert(coor_node[m] * node_site[m] == left[m]);
  }
  geo.init(coor_node, size_node, node_site);
  qassert(geo.local_volume() == n_points);
  psel.set_mem_type(get_default_mem_type());
}

}  // namespace qlat
