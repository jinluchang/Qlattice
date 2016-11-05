#pragma once

#include <qlat/config.h>
#include <qlat/utils.h>
#include <qlat/mpi.h>

QLAT_START_NAMESPACE

struct Geometry
{
  bool initialized;
  //
  GeometryNode geon;
  //
  int multiplicity;
  // number of elements on each lattice site
  //
  Coordinate node_site; // size of the coordinate on local node.
  Coordinate expansion_left;
  Coordinate expansion_right;
  //
  Coordinate node_site_expanded;
  // node_site_expanded[i] = expansion_left[i] + node_site[i] + expansion_right[i]
  //
  void reset_node_site_expanded()
  {
    for (int i = 0; i < DIM; ++i) {
      node_site_expanded[i] = expansion_left[i] + node_site[i] + expansion_right[i];
    }
  }
  //
  void init()
  {
    memset(this, 0, sizeof(Geometry));
  }
  void init(const Coordinate total_site, const int multiplicity_)
  {
    init();
    geon = get_geometry_node();
    multiplicity = multiplicity_;
    for (int i = 0; i < DIM; ++i) {
      qassert(0 == total_site[i] % geon.size_node[i]);
      node_site[i] = total_site[i] / geon.size_node[i];
    }
    reset_node_site_expanded();
    initialized = true;
  }
  void init(const GeometryNode& geon_,
      const int multiplicity_,
      const Coordinate& node_site_)
  {
    init();
    geon = geon_;
    multiplicity = multiplicity_;
    node_site = node_site_;
    reset_node_site_expanded();
    initialized = true;
  }
  //
  void remult(const int multiplicity_) {
    multiplicity = multiplicity_;
    reset_node_site_expanded();
  }
  //
  void resize(const Coordinate& expansion_left_, const Coordinate& expansion_right_)
  {
#ifdef USE_MULTI_NODE
    expansion_left = expansion_left_;
    expansion_right = expansion_right_;
#endif
    reset_node_site_expanded();
  }
  void resize(const int thick)
  {
    const Coordinate expansion(thick, thick, thick, thick);
    resize(expansion, expansion);
  }
  //
  Geometry()
  {
    init();
  }
  //
  long offset_from_coordinate(const Coordinate& x) const
  {
    Coordinate xe = x;
    xe = xe + expansion_left;
    return qlat::index_from_coordinate(xe, node_site_expanded) * multiplicity;
  }
  //
  Coordinate coordinate_from_offset(const long offset) const
    // 0 <= offset < local_volume_expanded() * multiplicity
  {
    Coordinate x = qlat::coordinate_from_index(offset/multiplicity, node_site_expanded);
    x = x - expansion_left;
    return x;
  }
  //
  long recordFromCoordinate(const Coordinate& x) const
  {
    Coordinate xe = x;
    xe = xe + expansion_left;
    return qlat::index_from_coordinate(xe, node_site_expanded);
  }
  //
  Coordinate coordinateFromRecord(long record) const
    // 0 <= offset < local_volume_expanded() * multiplicity
  {
    Coordinate x = qlat::coordinate_from_index(record, node_site_expanded);
    x = x - expansion_left;
    return x;
  }
  //
  long index_from_coordinate(const Coordinate& x) const
    // 0 <= index < local_volume()
  {
    return qlat::index_from_coordinate(x, node_site);
  }
  //
  Coordinate coordinate_from_index(const long index) const
    // get local coordinate from index
    // 0 <= index < local_volume()
  {
    return qlat::coordinate_from_index(index, node_site);
  }
	//
  long offset_from_index(const long index) const
  {
    return offset_from_coordinate(coordinate_from_index(index));
  }
  //
  long g_index_from_g_coordinate(const Coordinate& xg) const
  {
    return qlat::index_from_coordinate(xg, total_site());
  }
  //
  bool is_on_node(const Coordinate& x) const
  {
    return -expansion_left[0] <= x[0] && x[0] < node_site[0] + expansion_right[0]
      && -expansion_left[1] <= x[1] && x[1] < node_site[1] + expansion_right[1]
      && -expansion_left[2] <= x[2] && x[2] < node_site[2] + expansion_right[2]
      && -expansion_left[3] <= x[3] && x[3] < node_site[3] + expansion_right[3];
  }
  //
  bool is_local(const Coordinate& x) const
  {
    bool b = true;
    for (int mu = 0; mu < DIM; mu++) {
      b = b && 0 <= x[mu] && x[mu] < node_site[mu];
    }
    return b;
  }
  //
  bool is_only_local() const
  {
    bool b = true;
    for (int i = 0; i < 4; i++) {
      b = b && expansion_left[i] == 0 && expansion_right[i] == 0;
    }
    return b;
  }
  //
  long local_volume() const
  {
    return node_site[0] * node_site[1] * node_site[2] * node_site[3];
  }
  //
  long local_volume_expanded() const
  {
    return node_site_expanded[0] * node_site_expanded[1] * node_site_expanded[2] * node_site_expanded[3];
  }
  //
  Coordinate total_site() const
  {
    return node_site * geon.size_node;
  }
  //
  Coordinate global_size() const
  {
    warn("use total_site()");
    return total_site();
  }
  //
  long total_volume() const
  {
    return local_volume() * geon.num_node;
  }
  //
  Coordinate coordinate_g_from_l(const Coordinate& xl) const
  {
    Coordinate xg;
    for (int mu = 0; mu < 4; mu++) {
      xg[mu] = xl[mu] + geon.coor_node[mu] * node_site[mu];
    }
    return xg;
  }
  //
  Coordinate coordinate_l_from_g(const Coordinate& xg) const
  {
    Coordinate xl;
    for (int mu = 0; mu < 4; mu++) {
      xl[mu] = xg[mu] - geon.coor_node[mu] * node_site[mu];
    }
    return xl;
  }
};

inline bool operator==(const Geometry& geo1, const Geometry& geo2)
{
  return geo1.initialized == geo2.initialized
    && geo1.geon == geo2.geon
    && geo1.multiplicity == geo2.multiplicity
    && geo1.node_site == geo2.node_site
    && geo1.expansion_left == geo2.expansion_left
    && geo1.expansion_right == geo2.expansion_right
    && geo1.node_site_expanded == geo2.node_site_expanded;
}

inline bool operator!=(const Geometry& geo1, const Geometry& geo2)
{
  return !(geo1 == geo2);
}

inline Geometry geo_resize(const Geometry& geo_, const int thick = 0)
{
  Geometry geo = geo_;
  geo.resize(thick);
  return geo;
}

inline Geometry geo_remult(const Geometry& geo_, const int multiplicity_ = 1)
{
  Geometry geo = geo_;
  geo.remult(multiplicity_);
  return geo;
}

inline Geometry geo_reform(const Geometry& geo_, const int multiplicity_ = 1, const int thick = 0)
{
  Geometry geo = geo_;
  geo.remult(multiplicity_);
  geo.resize(thick);
  return geo;
}

std::string show(const Geometry& geo)
{
  std::string s;
  s += ssprintf("{ initialized  = %s\n", ::show(geo.initialized).c_str());
  s += ssprintf(", geon         =\n%s\n", show(geo.geon).c_str());
  s += ssprintf(", node_site    = %s\n", show(geo.node_site).c_str());
  s += ssprintf(", expanLeft    = %s\n", show(geo.expansion_left).c_str());
  s += ssprintf(", expanRight   = %s\n", show(geo.expansion_right).c_str());
  s += ssprintf(", node_siteExp = %s }", show(geo.node_site_expanded).c_str());
  return s;
}

void swap(Geometry& geo1, Geometry& geo2)
{
  Geometry geo = geo1;
  geo1 = geo2;
  geo2 = geo;
}

inline bool is_matching_geo(const Geometry& geo1, const Geometry& geo2)
{
  return geo1.initialized == geo2.initialized
    && geo1.geon == geo2.geon
    && geo1.node_site == geo2.node_site;
}

inline bool is_matching_geo_mult(const Geometry& geo1, const Geometry& geo2)
{
  return is_matching_geo(geo1, geo2) && geo1.multiplicity == geo2.multiplicity;
}

inline bool is_initialized(const Geometry& geo)
{
  return geo.initialized;
}

QLAT_END_NAMESPACE
