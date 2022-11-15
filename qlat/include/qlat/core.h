#pragma once

#include <qlat-utils/timer.h>
#include <qlat-utils/complex.h>
#include <qlat-utils/array.h>
#include <qlat-utils/qutils-vec.h>

#include <qlat/config.h>

#if defined NO_ALIGNED_ALLOC
#define QLAT_ALIGNED_BYTES 1
// #define EIGEN_MALLOC_ALREADY_ALIGNED 0
#define EIGEN_MAX_ALIGN_BYTES 0
#define EIGEN_MAX_STATIC_ALIGN_BYTES 0
#else
#define QLAT_ALIGNED_BYTES 16 // should divide all Eigen matrix sizes (which can convert with GPT).
// #ifndef EIGEN_MALLOC_ALREADY_ALIGNED
// #define EIGEN_MALLOC_ALREADY_ALIGNED 1
// #endif
#ifndef EIGEN_MAX_ALIGN_BYTES
#define EIGEN_MAX_ALIGN_BYTES QLAT_ALIGNED_BYTES
#endif
#ifndef EIGEN_MAX_STATIC_ALIGN_BYTES
#define EIGEN_MAX_STATIC_ALIGN_BYTES QLAT_ALIGNED_BYTES
#endif
#endif

// #define ALIGN alignas(QLAT_ALIGNED_BYTES)
#define ALIGN __attribute__((aligned(QLAT_ALIGNED_BYTES)))

#ifndef USE_SINGLE_NODE
#define USE_MULTI_NODE
#endif

namespace qlat
{  //

const int DIMN = 4;

const int NUM_COLOR = 3;

typedef Complex ComplexT;  // default Complex type

inline void warn(const std::string& str = "")
{
  if (str != "") {
    displayln_info(ssprintf("WARNING: %s", str.c_str()));
  }
}

// --------------------

struct API Coordinate : public array<int, DIMN> {
  qacc Coordinate() { array<int, DIMN>::fill(0); }
  qacc Coordinate(int first, int second, int third, int fourth)
  {
    int* p = data();
    p[0] = first;
    p[1] = second;
    p[2] = third;
    p[3] = fourth;
  }
};

inline std::string show(const qlat::Coordinate& x)
{
  return ssprintf("%dx%dx%dx%d", x[0], x[1], x[2], x[3]);
}

qacc long product(const Coordinate& coor)
{
  long ret = 1;
  for (int i = 0; i < (int)coor.size(); i++) {
    ret *= coor[i];
  }
  return ret;
}

qacc Coordinate coordinate_from_index(long index, const Coordinate& size)
{
  Coordinate x;
  x[0] = index % size[0];
  index /= size[0];
  x[1] = index % size[1];
  index /= size[1];
  x[2] = index % size[2];
  index /= size[2];
  x[3] = index % size[3];
  return x;
}

qacc long index_from_coordinate(const Coordinate& x, const Coordinate& size)
{
  return ((((long)x[3] * (long)size[2]) + (long)x[2]) * (long)size[1] +
          (long)x[1]) *
             (long)size[0] +
         (long)x[0];
}

qacc int eo_from_coordinate(const Coordinate& xl)
{
  return 2 - (xl[0] + xl[1] + xl[2] + xl[3] + 16 * 1024 * 1024) % 2;
}

qacc long sqr(const qlat::Coordinate& xg)
{
  return sqr((long)xg[0]) + sqr((long)xg[1]) + sqr((long)xg[2]) +
         sqr((long)xg[3]);
}

qacc bool operator==(const Coordinate &c1, const Coordinate &c2)
{
  return c1[0] == c2[0] and c1[1] == c2[1] and c1[2] == c2[2] and
         c1[3] == c2[3];
}

qacc bool operator!=(const Coordinate &c1, const Coordinate &c2)
{
  return !(c1 == c2);
}

qacc Coordinate operator+(const Coordinate &coor1, const Coordinate &coor2)
{
  return Coordinate(coor1[0] + coor2[0], coor1[1] + coor2[1],
                    coor1[2] + coor2[2], coor1[3] + coor2[3]);
}

qacc Coordinate operator-(const Coordinate &coor1, const Coordinate &coor2)
{
  return Coordinate(coor1[0] - coor2[0], coor1[1] - coor2[1],
                    coor1[2] - coor2[2], coor1[3] - coor2[3]);
}

qacc Coordinate operator-(const Coordinate &coor)
{
  return Coordinate(-coor[0], -coor[1], -coor[2], -coor[3]);
}

qacc Coordinate operator*(const int integer, const Coordinate &coor)
{
  return Coordinate(integer * coor[0], integer * coor[1], integer * coor[2],
                    integer * coor[3]);
}

qacc Coordinate operator*(const Coordinate &coor, const int integer)
{
  return integer * coor;
}

qacc Coordinate operator*(const Coordinate &coor1, const Coordinate &coor2)
{
  return Coordinate(coor1[0] * coor2[0], coor1[1] * coor2[1],
                    coor1[2] * coor2[2], coor1[3] * coor2[3]);
}

qacc Coordinate operator%(const Coordinate &coor1, const Coordinate &coor2)
{
  return Coordinate(coor1[0] % coor2[0], coor1[1] % coor2[1],
                    coor1[2] % coor2[2], coor1[3] % coor2[3]);
}

qacc Coordinate operator%(const Coordinate &coor, const int integer)
{
  return Coordinate(coor[0] % integer, coor[1] % integer, coor[2] % integer,
                    coor[3] % integer);
}

qacc Coordinate mod(const Coordinate& x, const Coordinate& size)
{
  Coordinate ret;
  ret[0] = mod(x[0], size[0]);
  ret[1] = mod(x[1], size[1]);
  ret[2] = mod(x[2], size[2]);
  ret[3] = mod(x[3], size[3]);
  return ret;
}

qacc Coordinate smod(const Coordinate& x, const Coordinate& size)
{
  Coordinate ret;
  ret[0] = smod(x[0], size[0]);
  ret[1] = smod(x[1], size[1]);
  ret[2] = smod(x[2], size[2]);
  ret[3] = smod(x[3], size[3]);
  return ret;
}

qacc Coordinate middle_mod(const Coordinate& x, const Coordinate& y,
                           const Coordinate& size)
{
  Coordinate ret;
  ret[0] = middle_mod(x[0], y[0], size[0]);
  ret[1] = middle_mod(x[1], y[1], size[1]);
  ret[2] = middle_mod(x[2], y[2], size[2]);
  ret[3] = middle_mod(x[3], y[3], size[3]);
  return ret;
}

// --------------------

struct API CoordinateD : public array<double, DIMN> {
  qacc CoordinateD() { memset(this, 0, sizeof(CoordinateD)); }
  qacc CoordinateD(const array<double, DIMN>& arr)
  {
    CoordinateD& c = *this;
    c = arr;
    qassert(false == qisnan(c));
  }
  qacc CoordinateD(const double x0, const double x1, const double x2,
                   const double x3)
  {
    qassert(DIMN == 4);
    CoordinateD& c = *this;
    c[0] = x0;
    c[1] = x1;
    c[2] = x2;
    c[3] = x3;
    qassert(false == qisnan(c));
  }
  qacc CoordinateD(const Coordinate& x)
  {
    CoordinateD& c = *this;
    for (int i = 0; i < DIMN; ++i) {
      c[i] = x[i];
    }
  }
};

qacc CoordinateD mod(const CoordinateD& x, const CoordinateD& size)
{
  CoordinateD ret;
  ret[0] = mod(x[0], size[0]);
  ret[1] = mod(x[1], size[1]);
  ret[2] = mod(x[2], size[2]);
  ret[3] = mod(x[3], size[3]);
  return ret;
}

qacc CoordinateD smod(const CoordinateD& x, const CoordinateD& size)
{
  CoordinateD ret;
  ret[0] = smod(x[0], size[0]);
  ret[1] = smod(x[1], size[1]);
  ret[2] = smod(x[2], size[2]);
  ret[3] = smod(x[3], size[3]);
  return ret;
}

qacc CoordinateD smod_sym(const CoordinateD& x, const CoordinateD& size)
{
  CoordinateD ret;
  ret[0] = smod_sym(x[0], size[0]);
  ret[1] = smod_sym(x[1], size[1]);
  ret[2] = smod_sym(x[2], size[2]);
  ret[3] = smod_sym(x[3], size[3]);
  return ret;
}

qacc CoordinateD middle_mod(const CoordinateD& x, const CoordinateD& y,
                            const CoordinateD& size)
{
  CoordinateD ret;
  ret[0] = middle_mod(x[0], y[0], size[0]);
  ret[1] = middle_mod(x[1], y[1], size[1]);
  ret[2] = middle_mod(x[2], y[2], size[2]);
  ret[3] = middle_mod(x[3], y[3], size[3]);
  return ret;
}

// --------------------

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
  Coordinate global_size() const
  {
    warn("use total_site()");
    return total_site();
  }
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

// --------------------




// --------------------



}  // namespace qlat
