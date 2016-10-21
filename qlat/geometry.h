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
    initialized = false;
  }
  void init(const GeometryNode& geon_,
      const int multiplicity_,
      const Coordinate& node_site_,
      const Coordinate& expansion_left_,
      const Coordinate& expansion_right_)
  {
    init();
    geon = geon_;
    multiplicity = multiplicity_;
    node_site = node_site_;
#ifdef USE_MULTI_NODE
    expansion_left = expansion_left_;
    expansion_right = expansion_right_;
#endif
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
  void init(const Coordinate total_site, const int multiplicity_)
  {
    init();
    geon = get_geometry_node();
    multiplicity = multiplicity_;
    for (int i = 0; i < DIM; ++i) {
      assert(0 == total_site[i] % geon.size_node[i]);
      node_site[i] = total_site[i] / geon.size_node[i];
    }
    reset_node_site_expanded();
    initialized = true;
  }
  void init(const Geometry& geo_)
  {
    std::memcpy(this, &geo_, sizeof(Geometry));
  }
  void init_remult(const Geometry& geo_, const int multiplicity_ = 1)
  {
    init(geo_);
    remult(multiplicity_);
  }
  void init_resize(const Geometry& geo_, const int thick = 0)
  {
    init(geo_);
    resize(thick);
  }
  void init_reform(const Geometry& geo_, const int multiplicity_ = 1, const int thick = 0)
  {
    init(geo_);
    remult(multiplicity_);
    resize(thick);
  }
  void init(const Geometry& geo_, const int multiplicity_)
  {
    init(geo_);
    multiplicity = multiplicity_;
    reset_node_site_expanded();
  }
  void init(const Geometry& geo_, const int multiplicity_, const int thick)
  {
    init(geo_);
    multiplicity = multiplicity_;
#ifdef USE_MULTI_NODE
    const Coordinate expansion(thick, thick, thick, thick);
    expansion_left = expansion;
    expansion_right = expansion;
#endif
    reset_node_site_expanded();
  }
  //
  void copyOnlyLocal(const Geometry& geo_){
    this->init(geo_.geon, geo_.multiplicity, geo_.node_site);
    // only local
  }
  void copyButExpand(const Geometry& geo_, int thick){
    this->init(geo_); resize(thick);
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
  Geometry(const Geometry& geo_)
  {
    init(geo_);
  }
  //
  const Geometry& operator=(const Geometry& geo_)
  {
    init(geo_);
    assert(0 == memcmp(this, &geo_, sizeof(Geometry)));
    return *this;
  }
  //
  long offset_from_coordinate(const Coordinate& x) const
  {
    Coordinate xe = x;
    xe = xe + expansion_left;
    return qlat::index_from_coordinate(xe, node_site_expanded) * multiplicity;
  }
  //
  void coordinate_from_offset(Coordinate& x, int& m, const long offset) const
    // 0 <= offset < local_volume_expanded() * multiplicity
  {
    qlat::coordinate_from_index(x, offset/multiplicity, node_site_expanded);
    x = x - expansion_left;
    m = offset % multiplicity;
  }
  //
  long recordFromCoordinate(const Coordinate& x) const
  {
    Coordinate xe = x;
    xe = xe + expansion_left;
    return qlat::index_from_coordinate(xe, node_site_expanded);
  }
  //
  void coordinateFromRecord(Coordinate& x, long record) const
    // 0 <= offset < local_volume_expanded() * multiplicity
  {
    qlat::coordinate_from_index(x, record, node_site_expanded);
    x = x - expansion_left;
  }
  //
  long index_from_coordinate(const Coordinate& x) const
    // 0 <= index < local_volume()
  {
    return qlat::index_from_coordinate(x, node_site);
  }
  //
  void coordinate_from_index(Coordinate& x, const long index) const
    // get local coordinate from index
    // 0 <= index < local_volume()
  {
    qlat::coordinate_from_index(x, index, node_site);
  }
	//
  long offset_from_index(const long index) const
    // jtu
  {
    Coordinate coor;
    coordinate_from_index(coor, index);
    return offset_from_coordinate(coor);
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
    bool is_local_ = true;
    for(int mu = 0; mu < DIM; mu++) 
      is_local_ = is_local_ && 0 <= x[mu] && x[mu] < node_site[mu];
    return is_local_;
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
  int total_site(int mu) const
  {
    return node_site[mu] * geon.size_node[mu];
  }
  //
  Coordinate global_size() const
  {
    Coordinate ret;
    for(int i = 0; i < DIM; i++){
      ret[i] = total_site(i);
    }
    return ret;
  }
  //
  long total_volume() const
  {
    return local_volume() * geon.num_node;
  }
  //
  void coordinate_g_from_l(Coordinate& xg, const Coordinate& xl) const
  {
    for (int mu = 0; mu < 4; mu++) {
      xg[mu] = xl[mu] + geon.coor_node[mu] * node_site[mu];
    }
  }
  //
  void coordinate_l_from_g(Coordinate& xl, const Coordinate& xg) const
  {
    for (int mu = 0; mu < 4; mu++) {
      xl[mu] = xg[mu] - geon.coor_node[mu] * node_site[mu];
    }
  }
};

std::string show(const Geometry& geo)
{
  std::string s;
  s += ssprintf("{ initialized = %s\n", ::show(geo.initialized).c_str());
  s += ssprintf(", geon        =\n%s\n" , show(geo.geon).c_str());
  s += ssprintf(", node_site    = %s\n", show(geo.node_site).c_str());
  s += ssprintf(", expanLeft   = %s\n", show(geo.expansion_left).c_str());
  s += ssprintf(", expanRight  = %s\n", show(geo.expansion_right).c_str());
  s += ssprintf(", node_siteExp = %s }", show(geo.node_site_expanded).c_str());
  return s;
}

inline bool operator==(const Geometry& geo1, const Geometry& geo2)
{
  return 0 == memcmp(&geo1, &geo2, sizeof(Geometry));
}

inline bool operator!=(const Geometry& geo1, const Geometry& geo2)
{
  return !(geo1 == geo2);
}

inline bool is_matching_geo(const Geometry& geo1, const Geometry& geo2)
{
  bool b = true;
  b = b && geo1.initialized == geo2.initialized;
  b = b && geo1.geon == geo2.geon;
  b = b && geo1.multiplicity == geo2.multiplicity;
  for (int mu = 0; mu < 4; ++mu) {
    b = b && geo1.node_site[mu] == geo2.node_site[mu];
  }
  return b;
}

inline bool is_initialized(const Geometry& geo)
{
  return geo.initialized;
}

QLAT_END_NAMESPACE
