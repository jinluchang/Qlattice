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
  int eo; // 0:full; 1:odd ; 2:even
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
    for (int i = 0; i < DIMN; ++i) {
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
  void init(const GeometryNode& geon_,
      const int multiplicity_,
      const Coordinate& node_site_)
  {
    if (!initialized) {
      init();
      geon = geon_;
      multiplicity = multiplicity_;
      node_site = node_site_;
      reset_node_site_expanded();
      initialized = true;
    }
  }
  //
  void remult(const int multiplicity_) {
    multiplicity = multiplicity_;
    reset_node_site_expanded();
  }
  //
  void resize(const Coordinate& expansion_left_, const Coordinate& expansion_right_)
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
  Geometry(const Coordinate& total_site, const int multiplicity_)
  {
    init();
    init(total_site, multiplicity_);
  }
  //
  Coordinate mirror(const Coordinate& x) const
  {
    Coordinate ret = x;
    for (int mu = 0; mu < DIMN; ++mu) {
      if (geon.size_node[mu] == 1) {
        ret[mu] = mod(x[mu], node_site[mu]);
      }
    }
    return ret;
  }
  //
  long offset_from_coordinate(const Coordinate& x) const
  {
    Coordinate xe = mirror(x);
    if (eo == 0) {
      xe = xe + expansion_left;
      return qlat::index_from_coordinate(xe, node_site_expanded) * multiplicity;
    } else {
      qassert(eo == 1 or eo == 2);
      qassert(node_site % 2 == Coordinate());
      qassert((x[0] + x[1] + x[2] + x[3] + 16*1024*1024) % 2 == 2 - eo);
      xe = xe + expansion_left;
      return qlat::index_from_coordinate(xe, node_site_expanded)/2 * multiplicity;
    }
  }
  //
  Coordinate coordinate_from_offset(const long offset) const
    // 0 <= offset < local_volume_expanded() * multiplicity
  {
    Coordinate x;
    if (eo == 0) {
      x = qlat::coordinate_from_index(offset/multiplicity, node_site_expanded);
      x = x - expansion_left;
    } else {
      qassert(eo == 1 or eo == 2);
      qassert(node_site % 2 == Coordinate());
      x = qlat::coordinate_from_index(offset/multiplicity * 2, node_site_expanded);
      x = x - expansion_left;
      if ((x[0] + x[1] + x[2] + x[3] + 16*1024*1024) % 2 != 2 - eo) {
        x = qlat::coordinate_from_index(offset/multiplicity * 2 + 1, node_site_expanded);
        x = x - expansion_left;
      }
    }
    return x;
  }
  //
  long index_from_coordinate(const Coordinate& x) const
    // 0 <= index < local_volume()
  {
    const Coordinate xm = mirror(x);
    if (eo == 0) {
      return qlat::index_from_coordinate(xm, node_site);
    } else {
      qassert(eo == 1 or eo == 2);
      qassert(node_site % 2 == Coordinate());
      qassert((x[0] + x[1] + x[2] + x[3] + 16*1024*1024) % 2 == 2 - eo);
      return qlat::index_from_coordinate(xm, node_site) / 2;
    }
  }
  //
  Coordinate coordinate_from_index(const long index) const
    // get local coordinate from index
    // 0 <= index < local_volume()
  {
    if (eo == 0) {
      return qlat::coordinate_from_index(index, node_site);
    } else {
      qassert(eo == 1 or eo == 2);
      qassert(node_site % 2 == Coordinate());
      Coordinate x = qlat::coordinate_from_index(index * 2, node_site);
      if ((x[0] + x[1] + x[2] + x[3] + 16*1024*1024) % 2 != 2 - eo) {
        x = qlat::coordinate_from_index(index * 2 + 1, node_site);
      }
      return x;
    }
  }
	//
  long offset_from_index(const long index) const
  {
    return offset_from_coordinate(coordinate_from_index(index));
  }
  //
  long g_index_from_g_coordinate(const Coordinate& xg) const
  {
    const Coordinate ts = total_site();
    return qlat::index_from_coordinate(mod(xg, ts), ts);
  }
  //
  bool is_on_node(const Coordinate& x) const
  {
    for (int mu = 0; mu < DIMN; mu++) {
      if (not (-expansion_left[mu] <= x[mu] and x[mu] < node_site[mu] + expansion_right[mu] or geon.size_node[mu] == 1)) {
        return false;
      }
    }
    return eo == 0 or (x[0] + x[1] + x[2] + x[3] + 16*1024*1024) % 2 == 2 - eo;
  }
  //
  bool is_local(const Coordinate& x) const
  {
    for (int mu = 0; mu < DIMN; mu++) {
      if (not (0 <= x[mu] and x[mu] < node_site[mu] or geon.size_node[mu] == 1)) {
        return false;
      }
    }
    return eo == 0 or (x[0] + x[1] + x[2] + x[3] + 16*1024*1024) % 2 == 2 - eo;
  }
  //
  bool is_only_local() const
  {
    for (int i = 0; i < 4; i++) {
      if (expansion_left[i] != 0 or expansion_right[i] != 0) {
        return false;
      }
    }
    return true;
  }
  //
  long local_volume() const
  {
    if (eo == 0) {
      return node_site[0] * node_site[1] * node_site[2] * node_site[3];
    } else {
      qassert(eo == 1 or eo == 2);
      qassert(node_site[0] % 2 == 0);
      return node_site[0] * node_site[1] * node_site[2] * node_site[3] / 2;
    }
  }
  //
  long local_volume_expanded() const
  {
    if (eo == 0) {
      return node_site_expanded[0] * node_site_expanded[1] * node_site_expanded[2] * node_site_expanded[3];
    } else {
      qassert(eo == 1 or eo == 2);
      qassert(node_site[0] % 2 == 0);
      qassert(node_site_expanded[0] % 2 == 0);
      return node_site_expanded[0] * node_site_expanded[1] * node_site_expanded[2] * node_site_expanded[3] / 2;
    }
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
  //
  ///////////////////////////////////////////////////////////////////
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
};

inline bool operator==(const Geometry& geo1, const Geometry& geo2)
{
  return geo1.initialized == geo2.initialized
    && geo1.eo == geo2.eo
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

inline Geometry geo_resize(const Geometry& geo_, const Coordinate& expansion_left_, const Coordinate& expansion_right_)
{
  Geometry geo = geo_;
  geo.resize(expansion_left_, expansion_right_);
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

inline Geometry geo_eo(const Geometry& geo_, const int eo = 0)
  // 0:regular; 1:odd; 2:even
{
  Geometry geo = geo_;
  geo.eo = eo;
  return geo;
}

inline bool is_matching_geo(const Geometry& geo1, const Geometry& geo2)
{
  return geo1.initialized == geo2.initialized
    && geo1.geon == geo2.geon
    && geo1.node_site == geo2.node_site;
}

inline bool is_matching_geo_mult(const Geometry& geo1, const Geometry& geo2)
{
  return is_matching_geo(geo1, geo2)
    && geo1.eo == geo2.eo
    && geo1.multiplicity == geo2.multiplicity;
}

inline bool is_initialized(const Geometry& geo)
{
  return geo.initialized;
}

QLAT_END_NAMESPACE

namespace qshow {

inline std::string show(const qlat::Geometry& geo)
{
  std::string s;
  s += ssprintf("{ initialized   = %s\n", show(geo.initialized).c_str());
  s += ssprintf(", eo            = %s\n", show(geo.eo).c_str());
  s += ssprintf(", geon          =\n%s\n", show(geo.geon).c_str());
  s += ssprintf(", multiplicity  = %s\n", show(geo.multiplicity).c_str());
  s += ssprintf(", node_site     = %s\n", show(geo.node_site).c_str());
  s += ssprintf(", expan_left    = %s\n", show(geo.expansion_left).c_str());
  s += ssprintf(", expan_right   = %s\n", show(geo.expansion_right).c_str());
  s += ssprintf(", node_site_exp = %s }", show(geo.node_site_expanded).c_str());
  return s;
}

}

#ifndef USE_NAMESPACE
using namespace qshow;
#endif
