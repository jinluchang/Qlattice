#pragma once

#include <lqps/config.h>
#include <lqps/utils.h>
#include <lqps/mpi.h>
#include <lqps/geometry.h>

#include <vector>

LQPS_START_NAMESPACE

template <class M>
struct Field
{
  bool initialized;
  Geometry geo;
  std::vector<M> field;
  //
  void init()
  {
    initialized = false;
    geo.init();
    field.clear();
  }
  void init(const Geometry& geo_)
  {
    init();
    geo = geo_;
    field.resize(geo.localVolumeExpanded() * geo.multiplicity);
    setZero(*this);
    initialized = true;
  }
  void init(const Geometry& geo_, const int multiplicity_)
  {
    init();
    geo.init(geo_, multiplicity_);
    field.resize(geo.localVolumeExpanded() * geo.multiplicity);
    setZero(*this);
    initialized = true;
  }
  void init(const Field& f)
  {
    init();
    geo.init(f.geo);
    field = f.field;
    initialized = true;
  }
  //
  Field()
  {
    init();
  }
  Field(const Geometry& geo_)
  {
    init(geo_);
  }
  Field(const Field& f)
  {
    assert(false);
  }
  //
  const Field& operator=(const Field& f)
  {
#pragma omp parallel for
    for (long index = 0; index < geo.localVolume(); ++index) {
      Coordinate xl; geo.coordinateFromIndex(xl, index);
      for (int m = 0; m < geo.multiplicity; ++m) {
        this->getElem(xl,m) = f.getElem(xl,m);
      }
    }
    return *this;
  }
  //
  M& getElem(const long offset)
  {
    assert(0 <= offset && offset < field.size());
    return field[offset];
  }
  const M& getElem(const long offset) const
  {
    assert(0 <= offset && offset < field.size());
    return field[offset];
  }
  M& getElem(const Coordinate& x, const int m)
  {
    assert(geo.isOnNode(x));
    assert(0 <= m && m < geo.multiplicity);
    long offset = geo.offsetFromCoordinate(x) + m;
    return field[offset];
  }
  const M& getElem(const Coordinate& x, const int m) const
  {
    assert(geo.isOnNode(x));
    assert(0 <= m && m < geo.multiplicity);
    long offset = geo.offsetFromCoordinate(x) + m;
    return field[offset];
  }
};

template <class M>
bool isInitialized(const Field<M>& f)
{
  return f.initialized;
}

template <class M>
void setZero(Field<M>& f)
{
  setZero(f.field);
}

template <class M>
Vector<M> getData(const Field<M>& f)
{
  return getData(f.field);
}

template <class M>
void swap(Field<M>& f1, Field<M>& f2)
{
  assert(isInitialized(f1));
  assert(isInitialized(f1));
  assert(f1.geo == f2.geo);
  swap(f1.field, f2.field);
}

LQPS_END_NAMESPACE
