#pragma once

#include <lqps/config.h>
#include <lqps/utils.h>
#include <lqps/mpi.h>

LQPS_START_NAMESPACE

struct Geometry
{
  bool initialized;
  //
  GeometryNode geon;
  //
  int multiplicity;
  // number of elements on each lattice site
  //
  Coordinate nodeSite;
  Coordinate expansionLeft;
  Coordinate expansionRight;
  //
  Coordinate nodeSiteExpended;
  // nodeSiteExpended[i] = expansionLeft[i] + nodeSite[i] + expansionRight[i]
  //
  void resetDirOffset()
  {
    for (int i = 0; i < DIM; ++i) {
      nodeSiteExpended[i] = expansionLeft[i] + nodeSite[i] + expansionRight[i];
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
      const Coordinate& nodeSite_,
      const Coordinate& expansionLeft_,
      const Coordinate& expansionRight_)
  {
    init();
    geon = geon_;
    multiplicity = multiplicity_;
    nodeSite = nodeSite_;
#ifdef USE_MULTI_NODE
    expansionLeft = expansionLeft_;
    expansionRight = expansionRight_;
#endif
    resetDirOffset();
    initialized = true;
  }
  void init(const GeometryNode& geon_,
      const int multiplicity_,
      const Coordinate& nodeSite_)
  {
    init();
    geon = geon_;
    multiplicity = multiplicity_;
    nodeSite = nodeSite_;
    resetDirOffset();
    initialized = true;
  }
  void init(const Coordinate totalSite, const int multiplicity_)
  {
    init();
    geon = getGeometryNode();
    multiplicity = multiplicity_;
    for (int i = 0; i < DIM; ++i) {
      assert(0 == totalSite[i] % geon.sizeNode[i]);
      nodeSite[i] = totalSite[i] / geon.sizeNode[i];
    }
    resetDirOffset();
    initialized = true;
  }
  void init(const Geometry& geo_)
  {
    std::memcpy(this, &geo_, sizeof(Geometry));
  }
  void init(const Geometry& geo_, const int multiplicity_)
  {
    init(geo_);
    multiplicity = multiplicity_;
    resetDirOffset();
  }
  void init(const Geometry& geo_, const int multiplicity_, const int thick)
  {
    init(geo_);
    multiplicity = multiplicity_;
#ifdef USE_MULTI_NODE
    const Coordinate expansion({ thick, thick, thick, thick });
    expansionLeft = expansion;
    expansionRight = expansion;
#endif
    resetDirOffset();
  }
  //
  void resize(const Coordinate& expansionLeft_, const Coordinate& expansionRight_)
  {
#ifdef USE_MULTI_NODE
    expansionLeft = expansionLeft_;
    expansionRight = expansionRight_;
#endif
    resetDirOffset();
  }
  void resize(const int thick)
  {
    const Coordinate expansion({ thick, thick, thick, thick });
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
  long offsetFromCoordinate(const Coordinate& x) const
  {
    Coordinate xe = x;
    shiftCoordinateAdd(xe, expansionLeft);
    return lqps::indexFromCoordinate(x, nodeSiteExpended) * multiplicity;
  }
  //
  void coordinateFromOffset(Coordinate& x, long offset) const
  {
    lqps::coordinateFromIndex(x, offset/multiplicity, nodeSiteExpended);
    shiftCoordinateSub(x, expansionLeft);
  }
  // 0 <= offset < localVolumeExpanded() * multiplicity
  //
  long indexFromCoordinate(const Coordinate& x) const
  {
    return lqps::indexFromCoordinate(x, nodeSite);
  }
  // 0 <= index < localVolume()
  //
  void coordinateFromIndex(Coordinate& x, const long index) const
  {
    lqps::coordinateFromIndex(x, index, nodeSite);
  }
  // get local coordinate from index
  // 0 <= index < localVolume()
  //
  bool isOnNode(const Coordinate& x) const
  {
    return -expansionLeft[0] <= x[0] && x[0] < nodeSite[0] + expansionRight[0]
    && -expansionLeft[1] <= x[1] && x[1] < nodeSite[1] + expansionRight[1]
    && -expansionLeft[2] <= x[2] && x[2] < nodeSite[2] + expansionRight[2]
    && -expansionLeft[3] <= x[3] && x[3] < nodeSite[3] + expansionRight[3];
  }
  //
  bool isLocal(const Coordinate& x) const
  {
    return 0 <= x[0] && x[0] < nodeSite[0]
    && 0 <= x[1] && x[1] < nodeSite[1]
    && 0 <= x[2] && x[2] < nodeSite[2]
    && 0 <= x[3] && x[3] < nodeSite[3];
  }
  //
  bool isOnlyLocal() const
  {
    bool b = true;
    for (int i = 0; i < 4; i++) {
      b = b && expansionLeft[i] == 0 && expansionRight[i] == 0;
    }
    return b;
  }
  //
  long localVolume() const
  {
    return nodeSite[0] * nodeSite[1] * nodeSite[2] * nodeSite[3];
  }
  //
  long localVolumeExpanded() const
  {
    return nodeSiteExpended[0] * nodeSiteExpended[1] * nodeSiteExpended[2] * nodeSiteExpended[3];
  }
  //
  int totalSite(int mu) const
  {
    return nodeSite[mu] * geon.sizeNode[mu];
  }
  //
  long totalVolume() const
  {
    return localVolume() * geon.numNode;
  }
  //
  void coordinateGfL(Coordinate& xg, const Coordinate& xl) const
  {
    for (int mu = 0; mu < 4; mu++) {
      xg[mu] = xl[mu] + geon.coorNode[mu] * nodeSite[mu];
    }
  }
  //
  void coordinateLfG(Coordinate& xl, const Coordinate& xg) const
  {
    for (int mu = 0; mu < 4; mu++) {
      xl[mu] = xg[mu] - geon.coorNode[mu] * nodeSite[mu];
    }
  }
};

std::string show(const Geometry& geo)
{
  std::string s;
  ssprintf(s, "{ initialized = %s\n", show(geo.initialized).c_str());
  ssprintf(s, ", geon        =\n%s\n" , show(geo.geon).c_str());
  ssprintf(s, ", nodeSite    = %s\n", show(geo.nodeSite).c_str());
  ssprintf(s, ", expanLeft   = %s\n", show(geo.expansionLeft).c_str());
  ssprintf(s, ", expanRight  = %s\n", show(geo.expansionRight).c_str());
  ssprintf(s, ", nodeSiteExp = %s }", show(geo.nodeSiteExpended).c_str());
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

inline bool isMatchingGeo(const Geometry& geo1, const Geometry& geo2)
{
  bool b = true;
  b = b && geo1.initialized == geo2.initialized;
  b = b && geo1.geon == geo2.geon;
  b = b && geo1.multiplicity == geo2.multiplicity;
  for (int mu = 0; mu < 4; ++mu) {
    b = b && geo1.nodeSite[mu] == geo2.nodeSite[mu];
  }
  return b;
}

inline bool isInitialized(const Geometry& geo)
{
  return geo.initialized;
}

LQPS_END_NAMESPACE
