#pragma once

#include <qlat/config.h>
#include <qlat/geometry.h>
#include <qlat/mpi.h>
#include <qlat/utils.h>
// #include <qlat/field-utils.h>

#include <ctime>
#include <fstream>
#include <vector>

QLAT_START_NAMESPACE

template <class M>
struct Field {
  bool initialized;
  Geometry geo;
  std::vector<M> field;
  //
  virtual const std::string& cname()
  {
    static const std::string s = "Field";
    return s;
  }
  //
  virtual void init()
  {
    TIMER("Field::init()");
    initialized = false;
    geo.init();
    clear(field);
  }
  virtual void init(const Geometry& geo_)
  {
    if (!initialized) {
      TIMER("Field::init(geo)");
      init();
      geo = geo_;
      field.resize(geo.local_volume_expanded() * geo.multiplicity);
      set_zero(*this);
      initialized = true;
    } else {
      qassert(is_matching_geo_mult(geo_, geo));
    }
  }
  virtual void init(const Geometry& geo_, const int multiplicity_)
  {
    if (!initialized) {
      TIMER("Field::init(geo,mult)");
      init();
      geo = geo_remult(geo_, multiplicity_);
      field.resize(geo.local_volume_expanded() * geo.multiplicity);
      set_zero(*this);
      initialized = true;
    } else {
      if (not is_matching_geo_mult(geo_remult(geo_, multiplicity_), geo)) {
        displayln("old geo = " + show(geo));
        displayln("new geo = " + show(geo_remult(geo_, multiplicity_)));
        qassert(false);
      }
    }
  }
  virtual void init(const Field<M>& f)
  {
    if (!initialized) {
      TIMER("Field::init(f)");
      init();
      geo = f.geo;
      field = f.field;
      initialized = true;
    } else {
      (*this) = f;
    }
  }
  //
  Field<M>() { init(); }
  Field<M>(const Field<M>& f)
  {
    qassert(false == f.initialized);
    init();
  }
  //
  const Field<M>& operator=(const Field<M>& f)
  {
    TIMER("Field::operator=");
    if (this == &f) {
      return *this;
    }
    init(geo_resize(f.geo));
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
    qassert(0 <= offset && offset < (long)field.size());
    return field[offset];
  }
  const M& get_elem(const long offset) const
  {
    qassert(0 <= offset && offset < (long)field.size());
    return field[offset];
  }
  //
  M& get_elem(const Coordinate& x, const int m)
  {
    qassert(geo.is_on_node(x));
    qassert(0 <= m && m < geo.multiplicity);
    const long offset = geo.offset_from_coordinate(x) + m;
    return get_elem(offset);
  }
  const M& get_elem(const Coordinate& x, const int m) const
  {
    qassert(geo.is_on_node(x));
    qassert(0 <= m && m < geo.multiplicity);
    const long offset = geo.offset_from_coordinate(x) + m;
    return get_elem(offset);
  }
  //
  M& get_elem(const Coordinate& x)
  {
    qassert(1 == geo.multiplicity);
    return get_elem(x, 0);
  }
  const M& get_elem(const Coordinate& x) const
  {
    qassert(1 == geo.multiplicity);
    return get_elem(x, 0);
  }
  //
  Vector<M> get_elems(const Coordinate& x)
  {
    qassert(geo.is_on_node(x));
    const long offset = geo.offset_from_coordinate(x);
    return Vector<M>(&field[offset], geo.multiplicity);
  }
  Vector<M> get_elems_const(const Coordinate& x) const
  // Be cautious about the const property
  // 改不改靠自觉
  {
    if (not geo.is_on_node(x)) {
      displayln("Field::get_elems_const: x=" + show(x) + "\ngeo=" + show(geo));
      qassert(geo.is_on_node(x));
    }
    const long offset = geo.offset_from_coordinate(x);
    return Vector<M>(&field[offset], geo.multiplicity);
  }
  //
  Vector<M> get_elems(const long index)
  // qassert(geo.is_only_local())
  {
    return Vector<M>(&field[index * geo.multiplicity], geo.multiplicity);
  }
  Vector<M> get_elems_const(const long index) const
  // Be cautious about the const property
  // 改不改靠自觉
  // qassert(geo.is_only_local())
  {
    return Vector<M>(&field[index * geo.multiplicity], geo.multiplicity);
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
const Field<M>& operator+=(Field<M>& f, const Field<M>& f1)
{
  TIMER("field_operator+=");
  if (not is_initialized(f)) {
    f = f1;
    return f;
  }
  qassert(is_matching_geo_mult(f.geo, f1.geo));
  const Geometry& geo = f.geo;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    Coordinate x = geo.coordinate_from_index(index);
    for (int m = 0; m < geo.multiplicity; ++m) {
      f.get_elem(x, m) += f1.get_elem(x, m);
    }
  }
  return f;
}

template <class M>
const Field<M>& operator-=(Field<M>& f, const Field<M>& f1)
{
  TIMER("field_operator-=");
  if (not is_initialized(f)) {
    f.init(f1.geo);
    set_zero(f);
    f -= f1;
    return f;
  }
  qassert(is_matching_geo_mult(f.geo, f1.geo));
  const Geometry& geo = f.geo;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); index++) {
    Coordinate x = geo.coordinate_from_index(index);
    for (int m = 0; m < geo.multiplicity; m++) {
      f.get_elem(x, m) -= f1.get_elem(x, m);
    }
  }
  return f;
}

template <class M>
const Field<M>& operator*=(Field<M>& f, const double factor)
{
  TIMER("field_operator*=(F,D)");
  const Geometry& geo = f.geo;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); index++) {
    const Coordinate x = geo.coordinate_from_index(index);
    for (int m = 0; m < geo.multiplicity; m++) {
      f.get_elem(x, m) *= factor;
    }
  }
  return f;
}

template <class M>
const Field<M>& operator*=(Field<M>& f, const Complex factor)
{
  TIMER("field_operator*=(F,C)");
  const Geometry& geo = f.geo;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); index++) {
    Coordinate x = geo.coordinate_from_index(index);
    for (int m = 0; m < geo.multiplicity; m++) {
      f.get_elem(x, m) *= factor;
    }
  }
  return f;
}

template <class M>
double qnorm(const Field<M>& f)
{
  const Geometry& geo = f.geo;
  double sum = 0.0;
#pragma omp parallel
  {
    double psum = 0.0;
#pragma omp for nowait
    for (long index = 0; index < geo.local_volume(); ++index) {
      const Coordinate x = geo.coordinate_from_index(index);
      const Vector<M> fx = f.get_elems_const(x);
      for (int m = 0; m < geo.multiplicity; ++m) {
        psum += qnorm(fx[m]);
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

template <class M>
double qnorm_double(const Field<M>& f1, const Field<M>& f2)
{
  const Geometry& geo = f1.geo;
  qassert(geo.is_only_local());
  qassert(geo == f2.geo);
  double sum = qnorm_double(get_data(f1), get_data(f2));
  glb_sum(sum);
  return sum;
}

template <class M, int multiplicity>
struct FieldM : Field<M> {
  virtual const std::string& cname()
  {
    static const std::string s = "FieldM";
    return s;
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
  FieldM<M, multiplicity>() { init(); }
  FieldM<M, multiplicity>(const Field<M>& f)
  {
    qassert(false == f.initialized);
    init();
  }
  FieldM<M, multiplicity>(const FieldM<M, multiplicity>& f)
  {
    qassert(false == f.initialized);
    init();
  }
};

template <class M>
long get_data_size(const Field<M>& f)
// NOT including the expended parts, only local volume data size
// only size on one node
{
  return f.geo.local_volume() * f.geo.multiplicity * sizeof(M);
}

template <class M>
void qswap(Field<M>& f1, Field<M>& f2)
{
  std::swap(f1.initialized, f2.initialized);
  std::swap(f1.geo, f2.geo);
  std::swap(f1.field, f2.field);
}

QLAT_END_NAMESPACE
