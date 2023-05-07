#pragma once

#include <qlat/config.h>
#include <qlat/env.h>

#include <qlat-utils/core.h>
#include <qlat-utils/mat.h>
#include <qlat-utils/matrix-hmc.h>
#include <qlat-utils/qutils-vec.h>
#include <qlat-utils/coordinate.h>
#include <qlat-utils/mat-vec.h>

namespace qlat
{  //

struct API GeometryNode {
  bool initialized;
  // About node geometry.
  int num_node;
  // num_node = size_node[0] * size_node[1] * size_node[2] * size_node[3]
  int id_node;
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
  int eo;  // 0:full; 1:odd ; 2:even
  //
  int multiplicity;
  // number of elements on each lattice site
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
  qacc void init(const GeometryNode& geon_, const Coordinate& node_site_,
                 const int multiplicity_)
  {
    if (!initialized) {
      init();
      geon = geon_;
      node_site = node_site_;
      multiplicity = multiplicity_;
      reset_node_site_expanded();
      initialized = true;
    }
  }
  void init(const Coordinate& total_site, const int multiplicity_);
  //
  qacc void remult(const int multiplicity_) { multiplicity = multiplicity_; }
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
  Geometry(const Coordinate& total_site, const int multiplicity_)
  {
    init();
    init(total_site, multiplicity_);
  }
  //
  qacc Coordinate mirror(const Coordinate& x) const
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
  qacc long offset_from_coordinate(const Coordinate& x) const
  {
    Coordinate xe = mirror(x);
    if (eo == 0) {
      xe = xe + expansion_left;
      return qlat::index_from_coordinate(xe, node_site_expanded) * multiplicity;
    } else {
      qassert(eo == 1 or eo == 2);
      qassert(node_site % 2 == Coordinate());
      qassert(eo_from_coordinate(x) == eo);
      xe = xe + expansion_left;
      return qlat::index_from_coordinate(xe, node_site_expanded) / 2 *
             multiplicity;
    }
  }
  //
  qacc Coordinate coordinate_from_offset(const long offset) const
  // 0 <= offset < local_volume_expanded() * multiplicity
  {
    Coordinate x;
    if (eo == 0) {
      x = qlat::coordinate_from_index(offset / multiplicity,
                                      node_site_expanded);
      x = x - expansion_left;
    } else {
      qassert(eo == 1 or eo == 2);
      qassert(node_site % 2 == Coordinate());
      x = qlat::coordinate_from_index(offset / multiplicity * 2,
                                      node_site_expanded);
      x = x - expansion_left;
      if (eo_from_coordinate(x) != eo) {
        x = qlat::coordinate_from_index(offset / multiplicity * 2 + 1,
                                        node_site_expanded);
        x = x - expansion_left;
      }
    }
    return x;
  }
  //
  qacc long index_from_coordinate(const Coordinate& x) const
  // 0 <= index < local_volume()
  {
    const Coordinate xm = mirror(x);
    if (eo == 0) {
      return qlat::index_from_coordinate(xm, node_site);
    } else {
      qassert(eo == 1 or eo == 2);
      qassert(node_site % 2 == Coordinate());
      qassert(eo_from_coordinate(x) == eo);
      return qlat::index_from_coordinate(xm, node_site) / 2;
    }
  }
  //
  qacc Coordinate coordinate_from_index(const long index) const
  // get local coordinate from index
  // 0 <= index < local_volume()
  {
    if (eo == 0) {
      return qlat::coordinate_from_index(index, node_site);
    } else {
      qassert(eo == 1 or eo == 2);
      qassert(node_site % 2 == Coordinate());
      Coordinate x = qlat::coordinate_from_index(index * 2, node_site);
      if (eo_from_coordinate(x) != eo) {
        x = qlat::coordinate_from_index(index * 2 + 1, node_site);
      }
      return x;
    }
  }
  //
  qacc long offset_from_index(const long index) const
  {
    return offset_from_coordinate(coordinate_from_index(index));
  }
  //
  qacc long g_index_from_g_coordinate(const Coordinate& xg) const
  {
    const Coordinate ts = total_site();
    return qlat::index_from_coordinate(mod(xg, ts), ts);
  }
  //
  qacc bool is_on_node(const Coordinate& x) const
  {
    for (int mu = 0; mu < DIMN; mu++) {
      if (not((-expansion_left[mu] <= x[mu] and
               x[mu] < node_site[mu] + expansion_right[mu]) or
              geon.size_node[mu] == 1)) {
        return false;
      }
    }
    return eo == 0 or eo_from_coordinate(x) == eo;
  }
  //
  qacc bool is_local(const Coordinate& x) const
  {
    for (int mu = 0; mu < DIMN; mu++) {
      if (not((0 <= x[mu] and x[mu] < node_site[mu]) or
              geon.size_node[mu] == 1)) {
        return false;
      }
    }
    return eo == 0 or eo_from_coordinate(x) == eo;
  }
  //
  qacc Coordinate local_site() const { return node_site; }
  //
  qacc long local_volume() const
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
  qacc long local_volume_expanded() const
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
  qacc long total_volume() const { return local_volume() * geon.num_node; }
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
  //
  ///////////////////////////////////////////////////////////////////
  //
  // Coordinate global_size() const
  // {
  //   warn("use total_site()");
  //   return total_site();
  // }
  // //
  // long recordFromCoordinate(const Coordinate& x) const
  // {
  //   Coordinate xe = x;
  //   xe = xe + expansion_left;
  //   return qlat::index_from_coordinate(xe, node_site_expanded);
  // }
  // //
  // Coordinate coordinateFromRecord(long record) const
  // // 0 <= offset < local_volume_expanded() * multiplicity
  // {
  //   Coordinate x = qlat::coordinate_from_index(record, node_site_expanded);
  //   x = x - expansion_left;
  //   return x;
  // }
};

std::string show(const qlat::Geometry& geo);

qacc bool operator==(const Geometry& geo1, const Geometry& geo2)
{
  return geo1.initialized == geo2.initialized && geo1.eo == geo2.eo &&
         geo1.geon == geo2.geon && geo1.multiplicity == geo2.multiplicity &&
         geo1.node_site == geo2.node_site &&
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

qacc Geometry geo_remult(const Geometry& geo_, const int multiplicity_ = 1)
{
  Geometry geo = geo_;
  geo.remult(multiplicity_);
  return geo;
}

qacc Geometry geo_reform(const Geometry& geo_, const int multiplicity_ = 1,
                         const int thick = 0)
// do not change eo
{
  Geometry geo = geo_;
  geo.remult(multiplicity_);
  geo.resize(thick);
  return geo;
}

qacc Geometry geo_reform(const Geometry& geo_, const int multiplicity_,
                         const Coordinate& expansion_left_,
                         const Coordinate& expansion_right_)
// do not change eo
{
  Geometry geo = geo_;
  geo.remult(multiplicity_);
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

qacc bool is_matching_geo_mult(const Geometry& geo1, const Geometry& geo2)
{
  return is_matching_geo(geo1, geo2) && geo1.eo == geo2.eo &&
         geo1.multiplicity == geo2.multiplicity;
}

qacc bool is_matching_geo_included(const Geometry& geo1, const Geometry& geo2)
{
  bool include = is_matching_geo_mult(geo1, geo2);
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
struct API Field {
  // Avoid copy constructor when possible
  // (it is likely not be what you think it is)
  //
  bool initialized;
  box_acc<Geometry> geo;
  vector_acc<M> field;
  //
  void init();
  void init(const Geometry& geo_, const int multiplicity_ = 0);
  void init(const Field<M>& f);
  //
  Field() { init(); }
  Field(const Field<M>&) = default;
  Field(Field<M>&&) noexcept = default;
  //
  const Field<M>& operator=(const Field<M>& f);
  //
  qacc const Geometry& get_geo() const { return geo(); }
  //
  qacc M& get_elem_offset(const long offset)
  {
    qassert(0 <= offset && offset < (long)field.size());
    return field[offset];
  }
  qacc const M& get_elem_offset(const long offset) const
  {
    qassert(0 <= offset && offset < (long)field.size());
    return field[offset];
  }
  //
  qacc M& get_elem(const Coordinate& x, const int m)
  {
    const Geometry& geo_v = geo();
    qassert(geo_v.is_on_node(x));
    qassert(0 <= m && m < geo_v.multiplicity);
    const long offset = geo_v.offset_from_coordinate(x) + m;
    return get_elem_offset(offset);
  }
  qacc const M& get_elem(const Coordinate& x, const int m) const
  {
    const Geometry& geo_v = geo();
    qassert(geo_v.is_on_node(x));
    qassert(0 <= m && m < geo_v.multiplicity);
    const long offset = geo_v.offset_from_coordinate(x) + m;
    return get_elem_offset(offset);
  }
  //
  qacc M& get_elem(const Coordinate& x)
  {
    qassert(1 == geo().multiplicity);
    return get_elem(x, 0);
  }
  qacc const M& get_elem(const Coordinate& x) const
  {
    qassert(1 == geo().multiplicity);
    return get_elem(x, 0);
  }
  //
  qacc Vector<M> get_elems(const Coordinate& x)
  {
    const Geometry& geo_v = geo();
    qassert(geo_v.is_on_node(x));
    const long offset = geo_v.offset_from_coordinate(x);
    return Vector<M>(&field[offset], geo_v.multiplicity);
  }
  qacc Vector<M> get_elems_const(const Coordinate& x) const
  // Be cautious about the const property
  // 改不改靠自觉
  {
    const Geometry& geo_v = geo();
    if (not geo_v.is_on_node(x)) {
#ifndef QLAT_IN_ACC
      displayln("Field::get_elems_const: x=" + show(x) +
                "\ngeo=" + show(geo_v));
#endif
      qassert(false);
    }
    const long offset = geo_v.offset_from_coordinate(x);
    return Vector<M>(&field[offset], geo_v.multiplicity);
  }
  //
  qacc M& get_elem(const long index, const int m)
  {
    const Geometry& geo_v = geo();
    qassert(geo_v.is_only_local);
    qassert(0 <= m && m < geo_v.multiplicity);
    return get_elem_offset(index * geo_v.multiplicity + m);
  }
  qacc const M& get_elem(const long index, const int m) const
  {
    const Geometry& geo_v = geo();
    qassert(geo_v.is_only_local);
    qassert(0 <= m && m < geo_v.multiplicity);
    return get_elem_offset(index * geo_v.multiplicity + m);
  }
  //
  qacc M& get_elem(const long index)
  {
    const Geometry& geo_v = geo();
    qassert(geo_v.is_only_local);
    qassert(1 == geo_v.multiplicity);
    return get_elem_offset(index);
  }
  qacc const M& get_elem(const long index) const
  {
    const Geometry& geo_v = geo();
    qassert(geo_v.is_only_local);
    qassert(1 == geo_v.multiplicity);
    return get_elem_offset(index);
  }
  //
  qacc Vector<M> get_elems(const long index)
  {
    const Geometry& geo_v = geo();
    qassert(geo_v.is_only_local);
    return Vector<M>(&field[index * geo_v.multiplicity], geo_v.multiplicity);
  }
  qacc Vector<M> get_elems_const(const long index) const
  // Be cautious about the const property
  // 改不改靠自觉
  {
    const Geometry& geo_v = geo();
    qassert(geo_v.is_only_local);
    return Vector<M>(&field[index * geo_v.multiplicity], geo_v.multiplicity);
  }
};

template <class M, int multiplicity>
struct API FieldM : Field<M> {
  using Field<M>::init;
  void init(const Geometry& geo_) { Field<M>::init(geo_, multiplicity); }
  void init(const Geometry& geo_, const int multiplicity_)
  {
    qassert(multiplicity == multiplicity_);
    Field<M>::init(geo_, multiplicity);
  }
  void init(const Field<M>& f)
  {
    qassert(multiplicity == f.geo().multiplicity);
    Field<M>::init(f);
  }
  //
  FieldM<M, multiplicity>() { init(); }
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
// only initialize if uninitilized
// if initialized already, then check for matching geo (including
// multiplicity)
// can have different geo expansion
{
  if (!initialized) {
    TIMER("Field::init(geo,mult)");
    init();
    if (multiplicity_ == 0) {
      geo.set(geo_);
    } else {
      geo.set(geo_remult(geo_, multiplicity_));
    }
    field.resize(geo().local_volume_expanded() * geo().multiplicity);
    if (1 == get_field_init()) {
      set_zero(*this);
    } else if (2 == get_field_init()) {
      set_u_rand_float(*this, RngState(show(get_time())));
    } else {
      qassert(0 == get_field_init());
    }
    initialized = true;
  } else {
    Geometry geo_new = geo_;
    if (multiplicity_ != 0) {
      geo_new.remult(multiplicity_);
    }
    if (not is_matching_geo_included(geo_new, geo())) {
      displayln("old geo = " + show(geo()));
      displayln("new geo = " + show(geo_new));
      qassert(false);
    }
  }
}

template <class M>
void Field<M>::init(const Field<M>& f)
// initialize to be identical to f if uninitilized
// otherwise use assignment operator
{
  if (!initialized) {
    TIMER("Field::init(f)");
    initialized = f.initialized;
    geo = f.geo;
    field = f.field;
  } else {
    (*this) = f;
  }
}

template <class M>
const Field<M>& Field<M>::operator=(const Field<M>& f)
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
  init(geo_resize(f.geo()));
  const int multiplicity = geo().multiplicity;
  Field<M>& f0 = *this;
  qacc_for(index, geo().local_volume(), {
    const Geometry& geo_v = f0.geo();
    const Coordinate xl = geo_v.coordinate_from_index(index);
    const Vector<M> v = f.get_elems_const(xl);
    Vector<M> v0 = f0.get_elems(xl);
    for (int m = 0; m < multiplicity; ++m) {
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
void set_unit(Field<M>& f, const Complex& coef = 1.0)
{
  TIMER("set_unit(Field)");
  for (long offset = 0; offset < f.field.size(); ++offset) {
    set_unit(f.get_elem_offset(offset), coef);
  }
}

template <class M>
qacc Vector<M> get_data(const Field<M>& f)
{
  return get_data(f.field);
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

using PointSelection = std::vector<Coordinate>;

// --------------------

template <class M>
struct API SelectedPoints {
  // Avoid copy constructor when possible
  // (it is likely not be what you think it is)
  //
  bool initialized;
  int multiplicity;
  long n_points;
  vector_acc<M> points;  // global quantity, same on each node
  // points.size() == n_points * multiplicity if initialized = true
  //
  void init();
  void init(const long n_points_, const int multiplicity_);
  void init(const PointSelection& psel, const int multiplicity);
  //
  SelectedPoints() { init(); }
  //
  qacc M& get_elem(const long& idx)
  {
    qassert(1 == multiplicity);
    return points[idx];
  }
  qacc const M& get_elem(const long& idx) const
  {
    qassert(1 == multiplicity);
    return points[idx];
  }
  //
  qacc M& get_elem(const long& idx, const int m)
  {
    qassert(0 <= m and m < multiplicity);
    return points[idx * multiplicity + m];
  }
  qacc const M& get_elem(const long& idx, const int m) const
  {
    qassert(0 <= m and m < multiplicity);
    return points[idx * multiplicity + m];
  }
  //
  qacc Vector<M> get_elems(const long idx)
  {
    return Vector<M>(&points[idx * multiplicity], multiplicity);
  }
  qacc Vector<M> get_elems_const(const long idx) const
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
  multiplicity = 0;
  n_points = 0;
  points.init();
}

template <class M>
void SelectedPoints<M>::init(const long n_points_, const int multiplicity_)
{
  if (initialized) {
    qassert(multiplicity_ == multiplicity);
    qassert(n_points_ == n_points);
    qassert((long)points.size() == n_points * multiplicity);
  } else {
    init();
    initialized = true;
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
void SelectedPoints<M>::init(const PointSelection& psel, const int multiplicity)
{
  init(psel.size(), multiplicity);
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

// --------------------

struct API FieldSelection {
  FieldM<int64_t, 1>
      f_rank;  // rank when the points being selected (-1 if not selected)
  //
  long n_per_tslice;  // num points per time slice (not enforced and should work
                      // properly if not true)
  double prob;        // (double)n_per_tslice / (double)spatial_vol
  //
  FieldM<long, 1>
      f_local_idx;  // idx of points on this node (-1 if not selected)
  long n_elems;     // num points of this node
  //
  vector_acc<int64_t> ranks;  // rank of the selected points
  vector_acc<long> indices;   // local indices of selected points
  //
  void init();
  //
  FieldSelection() { init(); }
  //
  qacc const Geometry& get_geo() const { return f_rank.geo(); }
};

// --------------------

template <class M>
struct API SelectedField {
  // Avoid copy constructor when possible
  // (it is likely not be what you think it is)
  //
  bool initialized;
  long n_elems;
  box_acc<Geometry> geo;
  vector_acc<M> field;  // field.size() == n_elems * multiplicity
  //
  void init();
  void init(const Geometry& geo_, const long n_elems_, const int multiplicity);
  void init(const FieldSelection& fsel, const int multiplicity);
  //
  SelectedField() { init(); }
  //
  qacc const Geometry& get_geo() const { return geo(); }
  //
  qacc M& get_elem(const long idx)
  {
    qassert(1 == geo().multiplicity);
    return field[idx];
  }
  qacc const M& get_elem(const long idx) const
  {
    qassert(1 == geo().multiplicity);
    return field[idx];
  }
  qacc M& get_elem(const long idx, const int m)
  {
    const int multiplicity = geo().multiplicity;
    qassert(0 <= m and m < multiplicity);
    return field[idx * multiplicity + m];
  }
  qacc const M& get_elem(const long idx, const int m) const
  {
    const int multiplicity = geo().multiplicity;
    qassert(0 <= m and m < multiplicity);
    return field[idx * multiplicity + m];
  }
  //
  qacc Vector<M> get_elems(const long idx)
  {
    const int multiplicity = geo().multiplicity;
    return Vector<M>(&field[idx * multiplicity], multiplicity);
  }
  qacc Vector<M> get_elems_const(const long idx) const
  // Be cautious about the const property
  // 改不改靠自觉
  {
    const int multiplicity = geo().multiplicity;
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
void SelectedField<M>::init(const Geometry& geo_, const long n_elems_,
                            const int multiplicity)
{
  if (initialized) {
    qassert(geo() == geo_remult(geo_, multiplicity));
    qassert(n_elems == n_elems_);
    qassert((long)field.size() == n_elems * multiplicity);
  } else {
    init();
    initialized = true;
    geo.set(geo_remult(geo_, multiplicity));
    n_elems = n_elems_;
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
void SelectedField<M>::init(const FieldSelection& fsel, const int multiplicity)
{
  init(fsel.f_rank.geo(), fsel.n_elems, multiplicity);
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
    const long gindex = geo.g_index_from_g_coordinate(xg);
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
    const long gindex = geo.g_index_from_g_coordinate(xg);
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
    const long gindex = geo.g_index_from_g_coordinate(xg);
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
  FUNC(double);                    \
  FUNC(float);                     \
  FUNC(int64_t);                   \
  FUNC(char);                      \
  FUNC(int8_t)

#define QLAT_CALL_WITH_TYPES_1(FUNC, TYPENAME) \
  FUNC(ColorMatrix, TYPENAME);                 \
  FUNC(WilsonMatrix, TYPENAME);                \
  FUNC(NonRelWilsonMatrix, TYPENAME);          \
  FUNC(IsospinMatrix, TYPENAME);               \
  FUNC(SpinMatrix, TYPENAME);                  \
  FUNC(WilsonVector, TYPENAME);                \
  FUNC(ComplexD, TYPENAME);                    \
  FUNC(ComplexF, TYPENAME);                    \
  FUNC(double, TYPENAME);                      \
  FUNC(float, TYPENAME);                       \
  FUNC(int64_t, TYPENAME);                     \
  FUNC(char, TYPENAME);                        \
  FUNC(int8_t, TYPENAME)

#define QLAT_CALL_WITH_TYPES_2(FUNC)                \
  QLAT_CALL_WITH_TYPES_1(FUNC, ColorMatrix);        \
  QLAT_CALL_WITH_TYPES_1(FUNC, WilsonMatrix);       \
  QLAT_CALL_WITH_TYPES_1(FUNC, NonRelWilsonMatrix); \
  QLAT_CALL_WITH_TYPES_1(FUNC, IsospinMatrix);      \
  QLAT_CALL_WITH_TYPES_1(FUNC, SpinMatrix);         \
  QLAT_CALL_WITH_TYPES_1(FUNC, WilsonVector);       \
  QLAT_CALL_WITH_TYPES_1(FUNC, ComplexD);           \
  QLAT_CALL_WITH_TYPES_1(FUNC, ComplexF);           \
  QLAT_CALL_WITH_TYPES_1(FUNC, double);             \
  QLAT_CALL_WITH_TYPES_1(FUNC, float);              \
  QLAT_CALL_WITH_TYPES_1(FUNC, int64_t);            \
  QLAT_CALL_WITH_TYPES_1(FUNC, char);               \
  QLAT_CALL_WITH_TYPES_1(FUNC, int8_t)

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
                                               const Complex& coef);          \
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

#define QLAT_EXTERN_CLASS                              \
                                                       \
  QLAT_EXTERN template struct FieldM<ColorMatrix, 4>;  \
                                                       \
  QLAT_EXTERN template struct FieldM<ColorMatrix, 1>;  \
                                                       \
  QLAT_EXTERN template struct FieldM<WilsonMatrix, 1>; \
                                                       \
  QLAT_EXTERN template struct FieldM<WilsonVector, 1>

QLAT_CALL_WITH_TYPES(QLAT_EXTERN_TEMPLATE);

QLAT_EXTERN_CLASS;

#undef QLAT_EXTERN
#undef QLAT_EXTERN_TEMPLATE
#undef QLAT_EXTERN_CLASS

// --------------------

}  // namespace qlat
