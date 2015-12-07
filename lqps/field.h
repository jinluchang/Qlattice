#pragma once

#include <lqps/config.h>
#include <lqps/utils.h>
#include <lqps/mpi.h>
#include <lqps/geometry.h>

#include <omp.h>

#include <vector>
#include <cassert>

LQPS_START_NAMESPACE

template <class M>
struct Field
{
  bool initialized;
  Geometry geo;
  std::vector<M> field;
  //
  virtual const char* cname()
  {
    return "Field";
  }
  //
  virtual void init()
  {
    initialized = false;
    geo.init();
    field.clear();
  }
  virtual void init(const Geometry& geo_)
  {
    init();
    geo = geo_;
    field.resize(geo.localVolumeExpanded() * geo.multiplicity);
    setZero(*this);
    initialized = true;
  }
  virtual void init(const Geometry& geo_, const int multiplicity_)
  {
    init();
    geo.init(geo_, multiplicity_);
    field.resize(geo.localVolumeExpanded() * geo.multiplicity);
    setZero(*this);
    initialized = true;
  }
  virtual void init(const Field& f)
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
  Field(const Field& f)
  {
    assert(false);
  }
  //
  const Field& operator=(const Field& f)
  {
    assert(isMatchingGeo(geo, f.geo));
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
  //
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
  //
  M& getElem(const Coordinate& x)
  {
    assert(1 == geo.multiplicity);
    return getElem(x,0);
  }
  const M& getElem(const Coordinate& x) const
  {
    assert(1 == geo.multiplicity);
    return getElem(x,0);
  }
  //
  Vector<M> getElems(const Coordinate& x) const
  {
    assert(geo.isOnNode(x));
    long offset = geo.offsetFromCoordinate(x);
    return Vector<M>((M*)&field[offset], geo.multiplicity);
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

template<class M>
const Field<M>& operator+=(Field<M>& f, const Field<M>& f1)
{
  TIMER("fieldOperator");
  assert(isMatchingGeo(f.geo, f1.geo));
  const Geometry& geo = f.geo;
#pragma omp parallel for
  for (long index = 0; index < geo.localVolume(); ++index) {
    Coordinate x; geo.coordinateFromIndex(x, index);
    for (int m = 0; m < geo.multiplicity; ++m) {
      f.getElem(x,m) += f1.getElem(x,m);
    }
  }
  return f;
}

template<class M>
const Field<M>& operator-=(Field<M>& f, const Field<M>& f1)
{
  TIMER("fieldOperator");
  assert(isMatchingGeo(f.geo, f1.geo));
  const Geometry& geo = f.geo;
#pragma omp parallel for
  for (long index = 0; index < geo.localVolume(); index++) {
    Coordinate x; geo.coordinateFromIndex(x, index);
    for (int m = 0; m < geo.multiplicity; m++) {
      f.getElem(x,m) -= f1.getElem(x,m);
    }
  }
  return f;
}

template<class M>
const Field<M>& operator*=(Field<M>& f, const double factor)
{
  TIMER("fieldOperator");
  const Geometry& geo = f.geo;
#pragma omp parallel for
  for (long index = 0; index < geo.localVolume(); index++) {
    Coordinate x; geo.coordinateFromIndex(x, index);
    for (int m = 0; m < geo.multiplicity; m++) {
      f.getElem(x,m) *= factor;
    }
  }
  return f;
}

template<class M>
const Field<M>& operator*=(Field<M>& f, const Complex factor)
{
  TIMER("fieldOperator");
  const Geometry& geo = f.geo;
#pragma omp parallel for
  for (long index = 0; index < geo.localVolume(); index++) {
    Coordinate x; geo.coordinateFromIndex(x, index);
    for (int m = 0; m < geo.multiplicity; m++) {
      f.getElem(x,m) *= factor;
    }
  }
  return f;
}

template<class M>
double norm(const Field<M>& f)
{
  const Geometry& geo = f.geo;
  double sum = 0.0;
#pragma omp parallel
  {
    double psum = 0.0;
#pragma omp for nowait
    for (long index = 0; index < geo.localVolume(); ++index) {
      Coordinate x; geo.coordinateFromIndex(x, index);
      Vector<M> fx = f.getElems(x);
      for (int m = 0; m < geo.multiplicity; ++m) {
        psum += norm(fx[m]);
      }
    }
    for (int i = 0; i < omp_get_num_threads(); ++i) {
#pragma omp barrier
      if (omp_get_thread_num() == i) {
          sum += psum;
      }
    }
  }
  sumVector(Vector<double>(&sum, 1));
  return sum;
}

template <class M, int multiplicity>
struct FieldM : Field<M>
{
    using Field<M>::init;
    virtual void init(const Geometry& geo_)
    {
      Field<M>::init(geo_, multiplicity);
    }
    virtual void init(const Geometry& geo_, const int multiplicity_)
    {
      assert(multiplicity == multiplicity_);
      Field<M>::init(geo_, multiplicity);
    }
    virtual void init(const Field<M>& f)
    {
      assert(multiplicity == f.geo.multiplicity);
      Field<M>::init(f);
    }
    //
    FieldM()
    {
      init();
    }
    FieldM(const FieldM<M,multiplicity>& f)
    {
      assert(false);
    }
};

LQPS_END_NAMESPACE
