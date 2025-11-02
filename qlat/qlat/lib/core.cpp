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

Geometry::Geometry(const Coordinate& total_site) { init(total_site); }

void Geometry::init(const Coordinate& total_site)
{
  const GeometryNode& geon_ = get_geometry_node();
  Coordinate node_site_;
  for (Int i = 0; i < DIMN; ++i) {
    Qassert(0 == total_site[i] % geon_.size_node[i]);
    node_site_[i] = total_site[i] / geon_.size_node[i];
  }
  init(geon_, node_site_);
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
    Qassert(false);
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
    Qassert(false);
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
  Qassert(spx.multiplicity == 1);
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
  Qassert(f.initialized);
  Qassert(f1.initialized);
  const Geometry geo = f.get_geo();
  Qassert(geo.is_only_local);
  Qassert(geo == f1.get_geo());
  const Int multiplicity = f.multiplicity;
  const Int multiplicity1 = f1.multiplicity;
  Qassert(multiplicity % sizeof_m == 0);
  Qassert(multiplicity1 % sizeof_m == 0);
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
  Qassert(fsel.f_rank.initialized);
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
  Qassert(psel.initialized);
  Qassert((psel.points_dist_type == PointsDistType::Local) or
          (psel.points_dist_type == PointsDistType::Full) or
          (psel.points_dist_type == PointsDistType::Global));
  Qassert(psel.total_site == geo.total_site());
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
    Qassert(fsel.n_elems == psel.size());
  }
}

void set_geo_from_psel(Geometry& geo, const PointsSelection& psel)
{
  TIMER("set_geo_from_psel");
  Qassert(psel.points_dist_type == PointsDistType::Full);
  psel.set_mem_type(MemType::Cpu);
  const Coordinate total_site = psel.total_site;
  const Long n_points = psel.size();
  Qassert(n_points > 0);
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
    Qassert(size_node[m] * node_site[m] == total_site[m]);
    Qassert(coor_node[m] * node_site[m] == left[m]);
  }
  geo.init(coor_node, size_node, node_site);
  Qassert(geo.local_volume() == n_points);
  psel.set_mem_type(get_default_mem_type());
}

}  // namespace qlat
