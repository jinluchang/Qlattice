#pragma once

#include <qlat/config.h>
#include <qlat/utils.h>
#include <qlat/mpi.h>
#include <qlat/geometry.h>

#include <omp.h>

#include <vector>
#include <ctime>
#include <fstream>

QLAT_START_NAMESPACE

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
    field.resize(geo.local_volume_expanded() * geo.multiplicity);
    set_zero(*this);
    initialized = true;
  }
  virtual void init(const Geometry& geo_, const int multiplicity_)
  {
    init();
    geo = geo_remult(geo_, multiplicity_);
    field.resize(geo.local_volume_expanded() * geo.multiplicity);
    set_zero(*this);
    initialized = true;
  }
  virtual void init(const Field& f)
  {
    init();
    geo = f.geo;
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
    qassert(false);
  }
  //
  const Field& operator=(const Field& f)
  {
    qassert(is_matching_geo_mult(geo, f.geo));
#pragma omp parallel for
    for (long index = 0; index < geo.local_volume(); ++index) {
      Coordinate xl = geo.coordinate_from_index(index);
      Vector<M> v = this->get_elems(xl);
      const Vector<M> v_ = f.get_elems_const(xl);
      for (int m = 0; m < geo.multiplicity; ++m) {
        v[m] = v_[m];
      }
    }
    return *this;
  }
  //
  M& get_elem(const long offset)
  {
    qassert(0 <= offset && offset < field.size());
    return field[offset];
  }
  const M& get_elem(const long offset) const
  {
    qassert(0 <= offset && offset < field.size());
    return field[offset];
  }
  //
  M& get_elem(const Coordinate& x, const int m)
  {
    qassert(geo.is_on_node(x));
    qassert(0 <= m && m < geo.multiplicity);
    long offset = geo.offset_from_coordinate(x) + m;
    return field[offset];
  }
  const M& get_elem(const Coordinate& x, const int m) const
  {
    qassert(geo.is_on_node(x));
    qassert(0 <= m && m < geo.multiplicity);
    long offset = geo.offset_from_coordinate(x) + m;
    return field[offset];
  }
  //
  M& get_elem(const Coordinate& x)
  {
    qassert(1 == geo.multiplicity);
    return get_elem(x,0);
  }
  const M& get_elem(const Coordinate& x) const
  {
    qassert(1 == geo.multiplicity);
    return get_elem(x,0);
  }
  //
  Vector<M> get_elems(const Coordinate& x)
  {
    qassert(geo.is_on_node(x));
    long offset = geo.offset_from_coordinate(x);
    return Vector<M>(&field[offset], geo.multiplicity);
  }
  Vector<M> get_elems_const(const Coordinate& x) const
    // Be cautious about the const property
    // 改不改靠自觉
  {
    qassert(geo.is_on_node(x));
    long offset = geo.offset_from_coordinate(x);
    return Vector<M>(&field[offset], geo.multiplicity);
  }
};

template <class M>
bool is_initialized(const Field<M>& f)
{
  return f.initialized;
}

template <class M>
Vector<M> get_data(const Field<M>& f)
{
  return get_data(f.field);
}

template <class M>
void swap(Field<M>& f1, Field<M>& f2)
{
  qassert(is_initialized(f1));
  qassert(is_initialized(f1));
  swap(f1.geo, f2.geo);
  swap(f1.field, f2.field);
}

template<class M>
const Field<M>& operator+=(Field<M>& f, const Field<M>& f1)
{
  TIMER("fieldOperator");
  qassert(is_matching_geo_mult(f.geo, f1.geo));
  const Geometry& geo = f.geo;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    Coordinate x = geo.coordinate_from_index(index);
    for (int m = 0; m < geo.multiplicity; ++m) {
      f.get_elem(x,m) += f1.get_elem(x,m);
    }
  }
  return f;
}

template<class M>
const Field<M>& operator-=(Field<M>& f, const Field<M>& f1)
{
  TIMER("fieldOperator");
  qassert(is_matching_geo_mult(f.geo, f1.geo));
  const Geometry& geo = f.geo;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); index++) {
    Coordinate x = geo.coordinate_from_index(index);
    for (int m = 0; m < geo.multiplicity; m++) {
      f.get_elem(x,m) -= f1.get_elem(x,m);
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
  for (long index = 0; index < geo.local_volume(); index++) {
    Coordinate x = geo.coordinate_from_index(index);
    for (int m = 0; m < geo.multiplicity; m++) {
      f.get_elem(x,m) *= factor;
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
  for (long index = 0; index < geo.local_volume(); index++) {
    Coordinate x = geo.coordinate_from_index(index);
    for (int m = 0; m < geo.multiplicity; m++) {
      f.get_elem(x,m) *= factor;
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
    for (long index = 0; index < geo.local_volume(); ++index) {
      Coordinate x = geo.coordinate_from_index(index);
      Vector<M> fx = f.get_elems(x);
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
  glb_sum(sum);
  return sum;
}

template <class M, int multiplicity>
struct FieldM : Field<M>
{
  virtual const char* cname()
  {
    return "FieldM";
  }
  //
  using Field<M>::init;
  virtual void init(const Geometry& geo_)
  {
    Field<M>::init(geo_, multiplicity);
  }
  virtual void init(const Geometry& geo_, const int multiplicity_)
  {
    qassert(multiplicity == multiplicity_);
    Field<M>::init(geo_, multiplicity);
  }
  virtual void init(const Field<M>& f)
  {
    qassert(multiplicity == f.geo.multiplicity);
    Field<M>::init(f);
  }
  //
  FieldM()
  {
    init();
  }
  FieldM(const FieldM<M,multiplicity>& f)
  {
    qassert(false);
  }
};

template <class M>
long get_data_size(const Field<M>& f)
  // NOT including the expended parts, only local volume data size
  // only size on one node
{
  return f.geo.local_volume() * f.geo.multiplicity * sizeof(M);
}

QLAT_END_NAMESPACE
