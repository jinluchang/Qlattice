#pragma once

#include <qlat/config.h>
#include <qlat/env.h>

#include <qlat-utils/core.h>
#include <qlat-utils/mat.h>
#include <qlat-utils/matrix-hmc.h>
#include <qlat-utils/utils-vec.h>
#include <qlat-utils/coordinate.h>
#include <qlat-utils/mat-vec.h>

namespace qlat
{  //

struct API GeometryNode {
  bool initialized;
  // About node geometry.
  Int num_node;
  // num_node = size_node[0] * size_node[1] * size_node[2] * size_node[3]
  Int id_node;
  // id_node = get_id_node()
  // 0 <= id_node < num_node
  Coordinate size_node;
  Coordinate coor_node;
  // 0 <= coor_node[i] < size_node[i]
  //
  qacc void init() { memset((void*)this, 0, sizeof(GeometryNode)); }
  qacc void init(const int id_node_, const Coordinate& size_node_)
  {
    initialized = true;
    num_node = product(size_node_);
    id_node = id_node_;
    size_node = size_node_;
    coor_node = coordinate_from_index(id_node_, size_node_);
  }
  //
  qacc GeometryNode() { init(); }
  qacc GeometryNode(const int id_node_, const Coordinate& size_node_)
  {
    init(id_node_, size_node_);
  }
};

std::string show(const GeometryNode& geon);

qacc bool operator==(const GeometryNode& geon1, const GeometryNode& geon2)
{
  return geon1.initialized == geon2.initialized &&
         geon1.num_node == geon2.num_node && geon1.id_node == geon2.id_node &&
         geon1.size_node == geon2.size_node &&
         geon1.coor_node == geon2.coor_node;
}

qacc bool operator!=(const GeometryNode& geon1, const GeometryNode& geon2)
{
  return !(geon1 == geon2);
}

API inline GeometryNode& get_geometry_node_internal()
{
  static GeometryNode geon;
  return geon;
}

inline const GeometryNode& get_geometry_node()
{
  return get_geometry_node_internal();
}

// --------------------

struct API Geometry {
  bool initialized;
  //
  GeometryNode geon;
  //
  Int eo;  // 0:full; 1:odd ; 2:even
  //
  Coordinate node_site;  // size of the coordinate on local node.
  Coordinate expansion_left;
  Coordinate expansion_right;
  //
  Coordinate node_site_expanded;
  // node_site_expanded[i] = expansion_left[i] + node_site[i] +
  // expansion_right[i]
  //
  bool is_only_local;
  //
  qacc void reset_node_site_expanded()
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
  //
  qacc void init() { memset((void*)this, 0, sizeof(Geometry)); }
  qacc void init(const GeometryNode& geon_, const Coordinate& node_site_)
  {
    if (!initialized) {
      init();
      geon = geon_;
      node_site = node_site_;
      reset_node_site_expanded();
      initialized = true;
    }
  }
  void init(const Coordinate& total_site);
  //
  qacc void resize(const Coordinate& expansion_left_,
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
  qacc void resize(const int thick)
  {
    const Coordinate expansion(thick, thick, thick, thick);
    resize(expansion, expansion);
  }
  //
  qacc Geometry() { init(); }
  //
  Geometry(const Coordinate& total_site)
  {
    init();
    init(total_site);
  }
  //
  qacc Coordinate mirror(const Coordinate& xl) const
  // avoid communicate in direction mu when geon.size_node[mu] == 1
  {
    Coordinate ret = xl;
    for (int mu = 0; mu < DIMN; ++mu) {
      if (geon.size_node[mu] == 1) {
        ret[mu] = mod(xl[mu], node_site[mu]);
      }
    }
    return ret;
  }
  //
  qacc Long offset_from_coordinate(const Coordinate& xl, const Int multiplicity) const
  {
    Coordinate xe = mirror(xl);
    if (eo == 0) {
      xe = xe + expansion_left;
      return qlat::index_from_coordinate(xe, node_site_expanded) * multiplicity;
    } else {
      qassert(eo == 1 or eo == 2);
      qassert(node_site % 2 == Coordinate());
      qassert(eo_from_coordinate(xl) == eo);
      xe = xe + expansion_left;
      return qlat::index_from_coordinate(xe, node_site_expanded) / 2 *
             multiplicity;
    }
  }
  //
  qacc Coordinate coordinate_from_offset(const Long offset, const Int multiplicity) const
  // 0 <= offset < local_volume_expanded() * multiplicity
  {
    Coordinate xl;
    if (eo == 0) {
      xl = qlat::coordinate_from_index(offset / multiplicity,
                                      node_site_expanded);
      xl = xl - expansion_left;
    } else {
      qassert(eo == 1 or eo == 2);
      qassert(node_site % 2 == Coordinate());
      xl = qlat::coordinate_from_index(offset / multiplicity * 2,
                                      node_site_expanded);
      xl = xl - expansion_left;
      if (eo_from_coordinate(xl) != eo) {
        xl = qlat::coordinate_from_index(offset / multiplicity * 2 + 1,
                                        node_site_expanded);
        xl = xl - expansion_left;
      }
    }
    return xl;
  }
  //
  qacc Long index_from_coordinate(const Coordinate& xl) const
  // 0 <= index < local_volume()
  {
    const Coordinate xm = mirror(xl);
    if (eo == 0) {
      return qlat::index_from_coordinate(xm, node_site);
    } else {
      qassert(eo == 1 or eo == 2);
      qassert(node_site % 2 == Coordinate());
      qassert(eo_from_coordinate(xl) == eo);
      return qlat::index_from_coordinate(xm, node_site) / 2;
    }
  }
  //
  qacc Coordinate coordinate_from_index(const Long index) const
  // get local coordinate from index
  // 0 <= index < local_volume()
  {
    if (eo == 0) {
      return qlat::coordinate_from_index(index, node_site);
    } else {
      qassert(eo == 1 or eo == 2);
      qassert(node_site % 2 == Coordinate());
      Coordinate xl = qlat::coordinate_from_index(index * 2, node_site);
      if (eo_from_coordinate(xl) != eo) {
        xl = qlat::coordinate_from_index(index * 2 + 1, node_site);
      }
      return xl;
    }
  }
  //
  qacc Long offset_from_index(const Long index, Int multiplicity) const
  {
    return offset_from_coordinate(coordinate_from_index(index), multiplicity);
  }
  //
  qacc Long g_index_from_g_coordinate(const Coordinate& xg) const
  {
    const Coordinate ts = total_site();
    return qlat::index_from_coordinate(mod(xg, ts), ts);
  }
  //
  qacc bool is_on_node(const Coordinate& xl) const
  {
    for (int mu = 0; mu < DIMN; mu++) {
      if (not((-expansion_left[mu] <= xl[mu] and
               xl[mu] < node_site[mu] + expansion_right[mu]) or
              geon.size_node[mu] == 1)) {
        return false;
      }
    }
    return eo == 0 or eo_from_coordinate(xl) == eo;
  }
  //
  qacc bool is_local(const Coordinate& xl) const
  {
    for (int mu = 0; mu < DIMN; mu++) {
      if (not((0 <= xl[mu] and xl[mu] < node_site[mu]) or
              geon.size_node[mu] == 1)) {
        return false;
      }
    }
    return eo == 0 or eo_from_coordinate(xl) == eo;
  }
  //
  qacc const Coordinate& local_site() const { return node_site; }
  //
  qacc Long local_volume() const
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
  qacc Long local_volume_expanded() const
  {
    if (eo == 0) {
      return node_site_expanded[0] * node_site_expanded[1] *
             node_site_expanded[2] * node_site_expanded[3];
    } else {
      qassert(eo == 1 or eo == 2);
      qassert(node_site[0] % 2 == 0);
      qassert(node_site_expanded[0] % 2 == 0);
      return node_site_expanded[0] * node_site_expanded[1] *
             node_site_expanded[2] * node_site_expanded[3] / 2;
    }
  }
  //
  qacc Coordinate total_site() const { return node_site * geon.size_node; }
  //
  qacc Long total_volume() const { return local_volume() * geon.num_node; }
  //
  qacc Coordinate coordinate_g_from_l(const Coordinate& xl) const
  {
    const Coordinate xg = xl + geon.coor_node * node_site;
    return xg;
  }
  //
  qacc Coordinate coordinate_l_from_g(const Coordinate& xg) const
  {
    const Coordinate xl = xg - geon.coor_node * node_site;
    return xl;
  }
};

std::string show(const qlat::Geometry& geo);

qacc bool operator==(const Geometry& geo1, const Geometry& geo2)
{
  return geo1.initialized == geo2.initialized && geo1.eo == geo2.eo &&
         geo1.geon == geo2.geon && geo1.node_site == geo2.node_site &&
         geo1.expansion_left == geo2.expansion_left &&
         geo1.expansion_right == geo2.expansion_right &&
         geo1.node_site_expanded == geo2.node_site_expanded &&
         geo1.is_only_local == geo2.is_only_local;
}

qacc bool operator!=(const Geometry& geo1, const Geometry& geo2)
{
  return !(geo1 == geo2);
}

qacc Geometry geo_resize(const Geometry& geo_, const int thick = 0)
{
  Geometry geo = geo_;
  geo.resize(thick);
  return geo;
}

qacc Geometry geo_resize(const Geometry& geo_,
                         const Coordinate& expansion_left_,
                         const Coordinate& expansion_right_)
{
  Geometry geo = geo_;
  geo.resize(expansion_left_, expansion_right_);
  return geo;
}

qacc Geometry geo_eo(const Geometry& geo_, const int eo_ = 0)
// 0:regular; 1:odd; 2:even
{
  Geometry geo = geo_;
  geo.eo = eo_;
  return geo;
}

qacc bool is_matching_geo(const Geometry& geo1, const Geometry& geo2)
{
  return geo1.initialized == geo2.initialized && geo1.geon == geo2.geon &&
         geo1.node_site == geo2.node_site;
}

qacc bool is_matching_geo_included(const Geometry& geo1, const Geometry& geo2)
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

// --------------------

template <class M>
struct API SelectedPoints;

enum struct PointsDistType {
  Global,
  Local,
  Random, // shuffle based on coordinate
};

std::string show(const PointsDistType points_dist_type);

PointsDistType read_points_dist_type(const std::string& points_dist_type_str);

struct API PointsSelection {
  bool initialized;
  PointsDistType points_dist_type;  // default PointsDistType::Global (all node has the same data)
  Coordinate total_site;
  vector_acc<Coordinate> xgs;
  //
  void init();
  void init(const Coordinate& total_site_, const Long n_points_,
            const PointsDistType points_dist_type_ = PointsDistType::Global);
  void init(const Coordinate& total_site_, const std::vector<Coordinate>& xgs_);
  void init(const Coordinate& total_site_, const vector<Coordinate>& xgs_);
  void init(const Coordinate& total_site_, const SelectedPoints<Coordinate>& spx);
  //
  PointsSelection() { init(); }
  PointsSelection(const PointsSelection&) = default;
  PointsSelection(PointsSelection&&) noexcept = default;
  PointsSelection(
      const Coordinate& total_site_, const Long n_points_,
      const PointsDistType points_dist_type_ = PointsDistType::Global)
  {
    init(total_site_, n_points_, points_dist_type_);
  }
  PointsSelection(const Coordinate& total_site_,
                  const std::vector<Coordinate>& xgs_)
  {
    init(total_site_, xgs_);
  }
  //
  PointsSelection& operator=(const PointsSelection& psel) = default;
  PointsSelection& operator=(PointsSelection&& psel) noexcept = default;
  //
  qacc Long size() const { return xgs.size(); }
  qacc const Coordinate* data() const { return xgs.data(); }
  qacc Coordinate* data() { return xgs.data(); }
  //
  qacc const Coordinate& operator[](const Long i) const { return xgs[i]; }
  qacc Coordinate& operator[](const Long i) { return xgs[i]; }
  //
  void resize(const Long size);
  void resize(const Long size, const Coordinate& xg_init);
  //
  SelectedPoints<Coordinate> view_sp() const;
  //
  void push_back_slow(const Coordinate& xg);  // Try to avoid. Very inefficient.
};

bool operator==(const PointsSelection& psel1, const PointsSelection& psel2);

bool operator!=(const PointsSelection& psel1, const PointsSelection& psel2);

// --------------------

template <class M>
struct API SelectedPoints {
  // Avoid copy constructor when possible
  // (it is likely not be what you think it is)
  //
  bool initialized;
  PointsDistType points_dist_type;  // default PointsDistType::Global (all node has the same data)
  Int multiplicity;
  Long n_points;
  vector_acc<M> points;  // global quantity, same on each node
  // points.size() == n_points * multiplicity if initialized = true
  //
  void init();
  void init(const Long n_points_, const int multiplicity_,
            const PointsDistType points_dist_type_);
  void init(const PointsSelection& psel, const int multiplicity_);
  //
  void init_zero(const Long n_points_, const int multiplicity_,
                 const PointsDistType points_dist_type_);
  void init_zero(const PointsSelection& psel, const int multiplicity);
  //
  SelectedPoints() { init(); }
  SelectedPoints(const SelectedPoints<M>&) = default;
  SelectedPoints(SelectedPoints<M>&&) noexcept = default;
  //
  SelectedPoints<M>& operator=(const SelectedPoints<M>&) = default;
  SelectedPoints<M>& operator=(SelectedPoints<M>&&) noexcept = default;
  //
  void set_view(const SelectedPoints<M>& sp)
  {
    TIMER("SelectedPoints::set_view");
    initialized = sp.initialized;
    points_dist_type = sp.points_dist_type;
    multiplicity = sp.multiplicity;
    n_points = sp.n_points;
    points.set_view(sp.points);
  }
  //
  template <class N>
  void set_view_cast(const SelectedPoints<N>& sp)
  {
    TIMER("SelectedPoints::set_view_cast");
    const Int total_size = sp.multiplicity * sizeof(N);
    initialized = sp.initialized;
    points_dist_type = sp.points_dist_type;
    multiplicity = total_size / sizeof(M);
    qassert(multiplicity * (Int)sizeof(M) == total_size);
    n_points = sp.n_points;
    points.set_view_cast(sp.points);
  }
  //
  SelectedPoints<M> view() const
  {
    TIMER("SelectedPoints::view");
    SelectedPoints<M> sp;
    sp.set_view(*this);
    return sp;
  }
  //
  SelectedPoints<Char> view_as_char() const
  {
    TIMER("SelectedPoints::view_as_char");
    SelectedPoints<Char> sp;
    sp.set_view_cast(*this);
    return sp;
  }
  //
  qacc bool is_view() { return points.is_copy; }
  //
  qacc M& get_elem(const Long& idx)
  {
    qassert(1 == multiplicity);
    return points[idx];
  }
  qacc const M& get_elem(const Long& idx) const
  {
    qassert(1 == multiplicity);
    return points[idx];
  }
  //
  qacc M& get_elem(const Long& idx, const int m)
  {
    qassert(0 <= m and m < multiplicity);
    return points[idx * multiplicity + m];
  }
  qacc const M& get_elem(const Long& idx, const int m) const
  {
    qassert(0 <= m and m < multiplicity);
    return points[idx * multiplicity + m];
  }
  //
  qacc Vector<M> get_elems(const Long idx)
  {
    return Vector<M>(&points[idx * multiplicity], multiplicity);
  }
  qacc Vector<M> get_elems_const(const Long idx) const
  // Be cautious about the const property
  // 改不改靠自觉
  {
    return Vector<M>(&points[idx * multiplicity], multiplicity);
  }
};

template <class M>
void SelectedPoints<M>::init()
{
  initialized = false;
  points_dist_type = PointsDistType::Global;
  multiplicity = 0;
  n_points = 0;
  points.init();
}

template <class M>
void SelectedPoints<M>::init(const Long n_points_, const int multiplicity_,
                             const PointsDistType points_dist_type_)
{
  if (initialized) {
    qassert(points_dist_type_ == points_dist_type);
    qassert(multiplicity_ == multiplicity);
    qassert(n_points_ == n_points);
    qassert((Long)points.size() == n_points * multiplicity);
  } else {
    TIMER("SelectedPoints::init(np,mult,dist)")
    init();
    initialized = true;
    points_dist_type = points_dist_type_;
    multiplicity = multiplicity_;
    n_points = n_points_;
    points.resize(n_points * multiplicity);
    if (1 == get_field_init()) {
      set_zero(*this);
    } else if (2 == get_field_init()) {
      set_u_rand_float(get_data(points), RngState(show(get_time())));
    } else {
      qassert(0 == get_field_init());
    }
  }
}

template <class M>
void SelectedPoints<M>::init(const PointsSelection& psel,
                             const int multiplicity)
{
  init(psel.size(), multiplicity, psel.points_dist_type);
}

template <class M>
void SelectedPoints<M>::init_zero(const Long n_points_, const int multiplicity_,
                                  const PointsDistType points_dist_type_)
{
  if (initialized) {
    qassert(points_dist_type_ == points_dist_type);
    qassert(multiplicity_ == multiplicity);
    qassert(n_points_ == n_points);
    qassert((Long)points.size() == n_points * multiplicity);
  } else {
    TIMER("SelectedPoints::init_zero(np,mult,dist)")
    init();
    initialized = true;
    points_dist_type = points_dist_type_;
    multiplicity = multiplicity_;
    n_points = n_points_;
    points.resize(n_points * multiplicity);
    set_zero(*this);
  }
}

template <class M>
void SelectedPoints<M>::init_zero(const PointsSelection& psel,
                                  const int multiplicity)
{
  init_zero(psel.size(), multiplicity, psel.points_dist_type);
}

template <class M>
Vector<M> get_data(const SelectedPoints<M>& sp)
{
  return get_data(sp.points);
}

template <class M>
void set_zero(SelectedPoints<M>& sp)
{
  TIMER("set_zero(SelectedPoints)");
  set_zero(get_data(sp));
}

template <class M>
void qswap(SelectedPoints<M>& f1, SelectedPoints<M>& f2)
{
  std::swap(f1.initialized, f2.initialized);
  std::swap(f1.n_points, f2.n_points);
  std::swap(f1.multiplicity, f2.multiplicity);
  qswap(f1.points, f2.points);
}

// --------------------

template <class M>
struct API Field {
  // Avoid copy constructor when possible
  // (it is likely not be what you think it is)
  //
  bool initialized;
  Int multiplicity;
  box_acc<Geometry> geo;
  vector_acc<M> field;
  //
  void init();
  void init(const Geometry& geo_, const int multiplicity_);
  void init(const Field<M>& f);
  //
  void init_zero(const Geometry& geo_, const int multiplicity_);
  //
  Field() { init(); }
  Field(const Field<M>&) = default;
  Field(Field<M>&&) noexcept = default;
  //
  Field<M>& operator=(const Field<M>& f);
  Field<M>& operator=(Field<M>&&) noexcept = default;
  //
  void set_view(const Field<M>& f)
  {
    TIMER("Field::set_view");
    initialized = f.initialized;
    geo.set_view(f.geo);
    multiplicity = f.multiplicity;
    field.set_view(f.field);
  }
  //
  template <class N>
  void set_view_cast(const Field<N>& f)
  {
    TIMER("Field::set_view_cast");
    const Int total_size = f.multiplicity * sizeof(N);
    initialized = f.initialized;
    multiplicity = total_size / sizeof(M);
    qassert(multiplicity * (Int)sizeof(M) == total_size);
    geo.set(f.geo());
    field.set_view_cast(f.field);
  }
  //
  SelectedPoints<M> view_sp() const
  {
    TIMER("Field::view_sp");
    qassert(geo().is_only_local);
    SelectedPoints<M> f;
    f.initialized = initialized;
    f.points_dist_type = PointsDistType::Local;
    f.multiplicity = multiplicity;
    f.n_points = geo().local_volume();
    f.points.set_view(field);
    return f;
  }
  //
  qacc const Geometry& get_geo() const { return geo(); }
  //
  qacc M& get_elem_offset(const Long offset)
  {
    qassert(0 <= offset && offset < (Long)field.size());
    return field[offset];
  }
  qacc const M& get_elem_offset(const Long offset) const
  {
    qassert(0 <= offset && offset < (Long)field.size());
    return field[offset];
  }
  //
  qacc M& get_elem(const Coordinate& x, const int m)
  {
    const Geometry& geo_v = geo();
    qassert(geo_v.is_on_node(x));
    qassert(0 <= m && m < multiplicity);
    const Long offset = geo_v.offset_from_coordinate(x, multiplicity) + m;
    return get_elem_offset(offset);
  }
  qacc const M& get_elem(const Coordinate& x, const int m) const
  {
    const Geometry& geo_v = geo();
    qassert(geo_v.is_on_node(x));
    qassert(0 <= m && m < multiplicity);
    const Long offset = geo_v.offset_from_coordinate(x, multiplicity) + m;
    return get_elem_offset(offset);
  }
  //
  qacc M& get_elem(const Coordinate& x)
  {
    qassert(1 == multiplicity);
    return get_elem(x, 0);
  }
  qacc const M& get_elem(const Coordinate& x) const
  {
    qassert(1 == multiplicity);
    return get_elem(x, 0);
  }
  //
  qacc Vector<M> get_elems(const Coordinate& x)
  {
    const Geometry& geo_v = geo();
    qassert(geo_v.is_on_node(x));
    const Long offset = geo_v.offset_from_coordinate(x, multiplicity);
    return Vector<M>(&field[offset], multiplicity);
  }
  qacc Vector<M> get_elems_const(const Coordinate& x) const
  // Be cautious about the const property
  // 改不改靠自觉
  {
    const Geometry& geo_v = geo();
    if (not geo_v.is_on_node(x)) {
#ifndef QLAT_IN_ACC
      displayln_c_stdout("Field::get_elems_const: x=" + show(x) +
                         "\ngeo=" + show(geo_v));
#endif
      qassert(false);
    }
    const Long offset = geo_v.offset_from_coordinate(x, multiplicity);
    return Vector<M>(&field[offset], multiplicity);
  }
  //
  qacc M& get_elem(const Long index, const int m)
  {
    const Geometry& geo_v = geo();
    qassert(geo_v.is_only_local);
    qassert(0 <= m && m < multiplicity);
    return get_elem_offset(index * multiplicity + m);
  }
  qacc const M& get_elem(const Long index, const int m) const
  {
    const Geometry& geo_v = geo();
    qassert(geo_v.is_only_local);
    qassert(0 <= m && m < multiplicity);
    return get_elem_offset(index * multiplicity + m);
  }
  //
  qacc M& get_elem(const Long index)
  {
    const Geometry& geo_v = geo();
    qassert(geo_v.is_only_local);
    if (1 != multiplicity) {
      qerr(ssprintf("Field::get_elem: mult=%d", multiplicity));
    }
    return get_elem_offset(index);
  }
  qacc const M& get_elem(const Long index) const
  {
    const Geometry& geo_v = geo();
    qassert(geo_v.is_only_local);
    qassert(1 == multiplicity);
    return get_elem_offset(index);
  }
  //
  qacc Vector<M> get_elems(const Long index)
  {
    const Geometry& geo_v = geo();
    qassert(geo_v.is_only_local);
    return Vector<M>(&field[index * multiplicity], multiplicity);
  }
  qacc Vector<M> get_elems_const(const Long index) const
  // Be cautious about the const property
  // 改不改靠自觉
  {
    const Geometry& geo_v = geo();
    qassert(geo_v.is_only_local);
    return Vector<M>(&field[index * multiplicity], multiplicity);
  }
};

template <class M, Int multiplicity>
struct API FieldM : Field<M> {
  void init() { Field<M>::init(); }
  void init(const Geometry& geo_) { Field<M>::init(geo_, multiplicity); }
  void init(const Geometry& geo_, const int multiplicity_)
  {
    qassert(multiplicity == multiplicity_);
    Field<M>::init(geo_, multiplicity);
  }
  void init(const Field<M>& f)
  {
    qassert(multiplicity == f.multiplicity);
    Field<M>::init(f);
  }
  //
  FieldM<M, multiplicity>() { init(); }
  FieldM(const FieldM<M, multiplicity>&) = default;
  FieldM(FieldM<M, multiplicity>&&) noexcept = default;
  FieldM<M, multiplicity>& operator=(FieldM<M, multiplicity>&&) noexcept =
      default;
  //
  FieldM<M, multiplicity>& operator=(const FieldM<M, multiplicity>& f)
  {
    qassert(f.multiplicity == multiplicity);
    Field<M>::operator=(f);
    return *this;
  }
};

template <class M>
void Field<M>::init()
{
  initialized = false;
  geo.init();
  field.init();
}

template <class M>
void Field<M>::init(const Geometry& geo_, const int multiplicity_)
// only initialize if uninitialized
// if initialized already, then check for matching geo (including
// multiplicity)
// can have different geo expansion
{
  if (initialized) {
    Geometry geo_new = geo_;
    if (not is_matching_geo_included(geo_new, geo())) {
      displayln("old geo = " + show(geo()));
      displayln("new geo = " + show(geo_new));
      qassert(false);
    }
    if (multiplicity != multiplicity_) {
      qerr(ssprintf("Field::init: mult=%d ; mult_=%d", multiplicity,
                    multiplicity_));
    }
  } else {
    TIMER("Field::init(geo,mult)");
    init();
    initialized = true;
    geo.set(geo_);
    multiplicity = multiplicity_;
    field.resize(geo().local_volume_expanded() * multiplicity);
    if (1 == get_field_init()) {
      set_zero(*this);
    } else if (2 == get_field_init()) {
      set_u_rand_float(*this, RngState(show(get_time())));
    } else {
      qassert(0 == get_field_init());
    }
  }
}

template <class M>
void Field<M>::init(const Field<M>& f)
// initialize to be identical to f if uninitilized
// otherwise use assignment operator
{
  if (initialized) {
    (*this) = f;
  } else {
    TIMER("Field::init(f)");
    initialized = f.initialized;
    multiplicity = f.multiplicity;
    geo = f.geo;
    field = f.field;
  }
}

template <class M>
void Field<M>::init_zero(const Geometry& geo_, const Int multiplicity_)
// only initialize and zero the field if uninitialized
// if initialized already, then check for matching geo (including
// multiplicity)
// can have different geo expansion (actual field needs to be larger)
// if check failed, the program crash
{
  if (initialized) {
    Geometry geo_new = geo_;
    if (not is_matching_geo_included(geo_new, geo())) {
      displayln("old geo = " + show(geo()));
      displayln("new geo = " + show(geo_new));
      qassert(false);
    }
    qassert(multiplicity == multiplicity_);
  } else {
    TIMER("Field::init_zero(geo,mult)");
    init();
    initialized = true;
    multiplicity = multiplicity_;
    geo.set(geo_);
    field.resize(geo().local_volume_expanded() * multiplicity);
    set_zero(*this);
  }
}

template <class M>
Field<M>& Field<M>::operator=(const Field<M>& f)
// skip if same object
// otherwise:
// 1. assert f is initialized
// 2. init with geo_resize(f.geo())
// 3. copy content
{
  if (this == &f) {
    return *this;
  }
  TIMER_FLOPS("Field::operator=");
  qassert(f.initialized);
  init(geo_resize(f.geo()), f.multiplicity);
  Field<M>& f0 = *this;
  qacc_for(index, geo().local_volume(), {
    const Geometry& geo_v = f0.geo();
    const Coordinate xl = geo_v.coordinate_from_index(index);
    const Vector<M> v = f.get_elems_const(xl);
    Vector<M> v0 = f0.get_elems(xl);
    for (int m = 0; m < f0.multiplicity; ++m) {
      v0[m] = v[m];
    }
  });
  timer.flops += get_data(f0).data_size();
  return *this;
}

template <class M>
void set_zero(Field<M>& f)
{
  TIMER("set_zero(Field)");
  set_zero(f.field);
}

template <class M>
void set_unit(Field<M>& f, const ComplexD& coef = 1.0)
{
  TIMER("set_unit(Field)");
  for (Long offset = 0; offset < f.field.size(); ++offset) {
    set_unit(f.get_elem_offset(offset), coef);
  }
}

template <class M>
qacc Vector<M> get_data(const Field<M>& f)
{
  return get_data(f.field);
}

template <class M>
void qswap(Field<M>& f1, Field<M>& f2)
{
  std::swap(f1.initialized, f2.initialized);
  std::swap(f1.multiplicity, f2.multiplicity);
  qswap(f1.geo, f2.geo);
  qswap(f1.field, f2.field);
}

// --------------------

template <class T = Real>
struct API GaugeFieldT : FieldM<ColorMatrixT<T>, 4> {
};

template <class T = Real>
struct API GaugeTransformT : FieldM<ColorMatrixT<T>, 1> {
};

template <class T = Real>
struct API Propagator4dT : FieldM<WilsonMatrixT<T>, 1> {
};

template <class T = Real>
struct API FermionField4dT : FieldM<WilsonVectorT<T>, 1> {
};

template <class T = Real>
struct API FermionField5dT : Field<WilsonVectorT<T> > {
};

using GaugeField = GaugeFieldT<>;

using GaugeTransform = GaugeTransformT<>;

using Propagator4d = Propagator4dT<>;

using FermionField4d = FermionField4dT<>;

using FermionField5d = FermionField5dT<>;

// --------------------

inline GaugeField& gf_from_field(Field<ColorMatrix>& f)
{
  qassert(f.multiplicity == 4);
  return (GaugeField&)f;
}

// --------------------

using FieldRank = FieldM<int64_t, 1>;

using FieldIndex = FieldM<Long, 1>;

struct API FieldSelection {
  FieldRank f_rank;  // rank when the points being selected (-1 if not selected)
  //
  // Update the following info with `void update_field_selection(FieldSelection& fsel)`.
  //
  Long n_elems;  // num points of this node
  //
  FieldIndex f_local_idx;  // idx of points on this node (-1 if not selected)
  //
  vector_acc<int64_t> ranks;  // rank of the selected points. `ranks.size() == n_elems`
  vector_acc<Long> indices;   // local indices of selected points. `indices.size() == n_elems`
  //
  void init();
  //
  FieldSelection() { init(); }
  //
  qacc const Geometry& get_geo() const { return f_rank.geo(); }
};

void set_psel_from_fsel(PointsSelection& psel, const FieldSelection& fsel);

// --------------------

template <class M>
struct API SelectedField {
  // Avoid copy constructor when possible
  // (it is likely not be what you think it is)
  //
  bool initialized;
  Long n_elems;
  Int multiplicity;
  box_acc<Geometry> geo;
  vector_acc<M> field;  // field.size() == n_elems * multiplicity
  //
  void init();
  void init(const Geometry& geo_, const Long n_elems_, const Int multiplicity_);
  void init(const FieldSelection& fsel, const Int multiplicity_);
  //
  void init_zero(const Geometry& geo_, const Long n_elems_, const Int multiplicity_);
  void init_zero(const FieldSelection& fsel, const Int multiplicity_);
  //
  SelectedField() { init(); }
  SelectedField(const SelectedField<M>&) = default;
  SelectedField(SelectedField<M>&&) noexcept = default;
  //
  SelectedField<M>& operator=(const SelectedField<M>&) = default;
  SelectedField<M>& operator=(SelectedField<M>&&) noexcept = default;
  //
  void set_view(const SelectedField<M>& sf)
  {
    TIMER("SelectedField::set_view");
    SelectedField<M> f;
    initialized = sf.initialized;
    n_elems = sf.n_elems;
    multiplicity = f.multiplicity;
    geo.set_view(sf.geo);
    field.set_view(sf.field);
  }
  //
  template <class N>
  void set_view_cast(const SelectedField<N>& sf)
  {
    TIMER("SelectedField::set_view_cast");
    const Int total_size = sf.multiplicity * sizeof(N);
    initialized = sf.initialized;
    n_elems = sf.n_elems;
    multiplicity = total_size / sizeof(M);
    qassert(multiplicity * (Int)sizeof(M) == total_size);
    geo.set(sf.geo());
    field.set_view_cast(sf.field);
  }
  //
  SelectedPoints<M> view_sp() const
  {
    TIMER("SelectedField::view_sp");
    SelectedPoints<M> f;
    f.initialized = initialized;
    f.points_dist_type = PointsDistType::Local;
    f.multiplicity = multiplicity;
    f.n_points = n_elems;
    f.points.set_view(field);
    return f;
  }
  //
  qacc const Geometry& get_geo() const { return geo(); }
  //
  qacc M& get_elem(const Long idx)
  {
    qassert(1 == multiplicity);
    return field[idx];
  }
  qacc const M& get_elem(const Long idx) const
  {
    qassert(1 == multiplicity);
    return field[idx];
  }
  qacc M& get_elem(const Long idx, const int m)
  {
    qassert(0 <= m and m < multiplicity);
    return field[idx * multiplicity + m];
  }
  qacc const M& get_elem(const Long idx, const int m) const
  {
    qassert(0 <= m and m < multiplicity);
    return field[idx * multiplicity + m];
  }
  //
  qacc Vector<M> get_elems(const Long idx)
  {
    return Vector<M>(&field[idx * multiplicity], multiplicity);
  }
  qacc Vector<M> get_elems_const(const Long idx) const
  // Be cautious about the const property
  // 改不改靠自觉
  {
    return Vector<M>(&field[idx * multiplicity], multiplicity);
  }
};

template <class M>
void SelectedField<M>::init()
{
  initialized = false;
  geo.init();
  field.init();
}

template <class M>
void SelectedField<M>::init(const Geometry& geo_, const Long n_elems_,
                            const Int multiplicity_)
{
  if (initialized) {
    qassert(geo() == geo_);
    qassert(n_elems == n_elems_);
    if (multiplicity != multiplicity_) {
      qerr(ssprintf("SelectedField::init: mult=%d ; mult_=%d", multiplicity,
                    multiplicity_));
    }
    qassert((Long)field.size() == n_elems * multiplicity);
  } else {
    TIMER("SelectedField::init(geo,n_elems,mult)")
    qassert(multiplicity_ != 0);
    init();
    initialized = true;
    geo.set(geo_);
    n_elems = n_elems_;
    multiplicity = multiplicity_;
    field.resize(n_elems * multiplicity);
    if (1 == get_field_init()) {
      set_zero(*this);
    } else if (2 == get_field_init()) {
      set_u_rand_float(get_data(field), RngState(show(get_time())));
    } else {
      qassert(0 == get_field_init());
    }
  }
}

template <class M>
void SelectedField<M>::init(const FieldSelection& fsel, const Int multiplicity_)
{
  init(fsel.f_rank.geo(), fsel.n_elems, multiplicity_);
}

template <class M>
void SelectedField<M>::init_zero(const Geometry& geo_, const Long n_elems_,
                                 const Int multiplicity_)
{
  if (initialized) {
    qassert(geo() == geo_);
    qassert(n_elems == n_elems_);
    qassert(multiplicity == multiplicity_);
    qassert((Long)field.size() == n_elems * multiplicity);
  } else {
    TIMER("SelectedField::init_zero(geo,n_elems,mult)");
    qassert(multiplicity_ != 0);
    init();
    initialized = true;
    geo.set(geo_);
    n_elems = n_elems_;
    multiplicity = multiplicity_;
    field.resize(n_elems * multiplicity);
    set_zero(*this);
  }
}

template <class M>
void SelectedField<M>::init_zero(const FieldSelection& fsel,
                                 const Int multiplicity_)
{
  init_zero(fsel.f_rank.geo(), fsel.n_elems, multiplicity_);
}

template <class M>
Vector<M> get_data(const SelectedField<M>& sf)
{
  return get_data(sf.field);
}

template <class M>
void set_zero(SelectedField<M>& sf)
{
  TIMER("set_zero(SelectedField)");
  set_zero(get_data(sf));
}

template <class M>
void qswap(SelectedField<M>& f1, SelectedField<M>& f2)
{
  std::swap(f1.initialized, f2.initialized);
  std::swap(f1.n_elems, f2.n_elems);
  std::swap(f1.multiplicity, f2.multiplicity);
  qswap(f1.geo, f2.geo);
  qswap(f1.field, f2.field);
}

// --------------------

using Prop = Propagator4d;

using SelProp = SelectedField<WilsonMatrix>;

using PselProp = SelectedPoints<WilsonMatrix>;

// --------------------

template <class M>
void set_u_rand_float(Field<M>& f, const RngState& rs, const double upper = 1.0,
                      const double lower = -1.0)
{
  TIMER("set_u_rand_float");
  const Geometry& geo = f.geo();
  qthread_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const Long gindex = geo.g_index_from_g_coordinate(xg);
    RngState rsi = rs.newtype(gindex);
    Vector<M> v = f.get_elems(xl);
    Vector<float> dv((float*)v.data(), v.data_size() / sizeof(float));
    for (int m = 0; m < dv.size(); ++m) {
      dv[m] = u_rand_gen(rsi, upper, lower);
    }
  });
}

template <class M>
void set_u_rand_double(Field<M>& f, const RngState& rs,
                       const double upper = 1.0, const double lower = -1.0)
{
  TIMER("set_u_rand_double");
  qthread_for(index, f.geo().local_volume(), {
    const Geometry& geo = f.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const Long gindex = geo.g_index_from_g_coordinate(xg);
    RngState rsi = rs.newtype(gindex);
    Vector<M> v = f.get_elems(xl);
    Vector<double> dv((double*)v.data(), v.data_size() / sizeof(double));
    for (int m = 0; m < dv.size(); ++m) {
      dv[m] = u_rand_gen(rsi, upper, lower);
    }
  });
}

template <class M>
void set_g_rand_double(Field<M>& f, const RngState& rs,
                       const double center = 0.0, const double sigma = 1.0)
{
  TIMER("set_g_rand_double");
  const Geometry& geo = f.geo();
  qthread_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const Long gindex = geo.g_index_from_g_coordinate(xg);
    RngState rsi = rs.newtype(gindex);
    Vector<M> v = f.get_elems(xl);
    Vector<double> dv((double*)v.data(), v.data_size() / sizeof(double));
    for (int m = 0; m < dv.size(); ++m) {
      dv[m] = g_rand_gen(rsi, center, sigma);
    }
  });
}

// --------------------

#define QLAT_CALL_WITH_TYPES(FUNC) \
  FUNC(ColorMatrix);               \
  FUNC(WilsonMatrix);              \
  FUNC(NonRelWilsonMatrix);        \
  FUNC(IsospinMatrix);             \
  FUNC(SpinMatrix);                \
  FUNC(WilsonVector);              \
  FUNC(ComplexD);                  \
  FUNC(ComplexF);                  \
  FUNC(RealD);                     \
  FUNC(RealF);                     \
  FUNC(Long);                      \
  FUNC(Int);                       \
  FUNC(Char)

#define QLAT_CALL_WITH_TYPES_1(FUNC, TYPENAME) \
  FUNC(ColorMatrix, TYPENAME);                 \
  FUNC(WilsonMatrix, TYPENAME);                \
  FUNC(NonRelWilsonMatrix, TYPENAME);          \
  FUNC(IsospinMatrix, TYPENAME);               \
  FUNC(SpinMatrix, TYPENAME);                  \
  FUNC(WilsonVector, TYPENAME);                \
  FUNC(ComplexD, TYPENAME);                    \
  FUNC(ComplexF, TYPENAME);                    \
  FUNC(RealD, TYPENAME);                       \
  FUNC(RealF, TYPENAME);                       \
  FUNC(Long, TYPENAME);                        \
  FUNC(Int, TYPENAME);                         \
  FUNC(Char, TYPENAME)

#define QLAT_CALL_WITH_TYPES_2(FUNC)                \
  QLAT_CALL_WITH_TYPES_1(FUNC, ColorMatrix);        \
  QLAT_CALL_WITH_TYPES_1(FUNC, WilsonMatrix);       \
  QLAT_CALL_WITH_TYPES_1(FUNC, NonRelWilsonMatrix); \
  QLAT_CALL_WITH_TYPES_1(FUNC, IsospinMatrix);      \
  QLAT_CALL_WITH_TYPES_1(FUNC, SpinMatrix);         \
  QLAT_CALL_WITH_TYPES_1(FUNC, WilsonVector);       \
  QLAT_CALL_WITH_TYPES_1(FUNC, ComplexD);           \
  QLAT_CALL_WITH_TYPES_1(FUNC, ComplexF);           \
  QLAT_CALL_WITH_TYPES_1(FUNC, RealD);              \
  QLAT_CALL_WITH_TYPES_1(FUNC, RealF);              \
  QLAT_CALL_WITH_TYPES_1(FUNC, Long);               \
  QLAT_CALL_WITH_TYPES_1(FUNC, Int);                \
  QLAT_CALL_WITH_TYPES_1(FUNC, Char)

#ifdef QLAT_INSTANTIATE_CORE
#define QLAT_EXTERN
#else
#define QLAT_EXTERN extern
#endif

#define QLAT_EXTERN_TEMPLATE(TYPENAME)                                        \
                                                                              \
  QLAT_EXTERN template struct Field<TYPENAME>;                                \
                                                                              \
  QLAT_EXTERN template struct SelectedField<TYPENAME>;                        \
                                                                              \
  QLAT_EXTERN template struct SelectedPoints<TYPENAME>;                       \
                                                                              \
  QLAT_EXTERN template void set_zero<TYPENAME>(Field<TYPENAME> & f);          \
                                                                              \
  QLAT_EXTERN template void set_zero<TYPENAME>(SelectedField<TYPENAME> & f);  \
                                                                              \
  QLAT_EXTERN template void set_zero<TYPENAME>(SelectedPoints<TYPENAME> & f); \
                                                                              \
  QLAT_EXTERN template void set_unit<TYPENAME>(Field<TYPENAME> & f,           \
                                               const ComplexD& coef);         \
                                                                              \
  QLAT_EXTERN template void qswap<TYPENAME>(Field<TYPENAME> & f1,             \
                                            Field<TYPENAME> & f2);            \
                                                                              \
  QLAT_EXTERN template void qswap<TYPENAME>(SelectedPoints<TYPENAME> & f1,    \
                                            SelectedPoints<TYPENAME> & f2);   \
                                                                              \
  QLAT_EXTERN template void qswap<TYPENAME>(SelectedField<TYPENAME> & f1,     \
                                            SelectedField<TYPENAME> & f2);    \
                                                                              \
  QLAT_EXTERN template void set_u_rand_double<TYPENAME>(                      \
      Field<TYPENAME> & f, const RngState& rs, const double upper,            \
      const double lower);                                                    \
                                                                              \
  QLAT_EXTERN template void set_u_rand_float<TYPENAME>(                       \
      Field<TYPENAME> & f, const RngState& rs, const double upper,            \
      const double lower);                                                    \
                                                                              \
  QLAT_EXTERN template void set_g_rand_double<TYPENAME>(                      \
      Field<TYPENAME> & f, const RngState& rs, const double center,           \
      const double sigma)

QLAT_CALL_WITH_TYPES(QLAT_EXTERN_TEMPLATE);
#undef QLAT_EXTERN_TEMPLATE

#define QLAT_EXTERN_CLASS                              \
                                                       \
  QLAT_EXTERN template struct FieldM<ColorMatrix, 4>;  \
                                                       \
  QLAT_EXTERN template struct FieldM<ColorMatrix, 1>;  \
                                                       \
  QLAT_EXTERN template struct FieldM<WilsonMatrix, 1>; \
                                                       \
  QLAT_EXTERN template struct FieldM<WilsonVector, 1>

QLAT_EXTERN_CLASS;
#undef QLAT_EXTERN_CLASS

#undef QLAT_EXTERN

// --------------------

}  // namespace qlat
