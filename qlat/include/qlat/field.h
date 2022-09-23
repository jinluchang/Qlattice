#pragma once

#include <qlat/config.h>
#include <qlat/geometry.h>
#include <qlat/mpi.h>

#include <qlat-utils/vector.h>

#include <ctime>
#include <fstream>
#include <vector>

namespace qlat
{  //

inline int& get_field_init();

template <class M>
struct Field {
  // Avoid copy constructor when possible
  // (it is likely not be what you think it is)
  //
  bool initialized;
  box_acc<Geometry> geo;
  vector_acc<M> field;
  //
  void init()
  {
    initialized = false;
    geo.init();
    field.init();
  }
  void init(const Geometry& geo_)
  // only initialize if uninitilized
  // if initialized already, then check for matching geo (including
  // multiplicity)
  // can have different geo expansion
  {
    if (!initialized) {
      TIMER("Field::init(geo)");
      init();
      geo.set(geo_);
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
      if (not is_matching_geo_included(geo_, geo())) {
        displayln("old geo = " + show(geo()));
        displayln("new geo = " + show(geo_));
        qassert(false);
      }
    }
  }
  void init(const Geometry& geo_, const int multiplicity_)
  // only initialize if uninitilized
  // if initialized already, then check for matching geo (including
  // multiplicity)
  // can have different geo expansion
  {
    if (!initialized) {
      TIMER("Field::init(geo,mult)");
      init();
      geo.set(geo_remult(geo_, multiplicity_));
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
      if (not is_matching_geo_included(geo_remult(geo_, multiplicity_), geo())) {
        displayln("old geo = " + show(geo()));
        displayln("new geo = " + show(geo_remult(geo_, multiplicity_)));
        qassert(false);
      }
    }
  }
  void init(const Field<M>& f)
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
  //
  Field() { init(); }
  Field(const Field<M>&) = default;
  Field(Field<M>&&) noexcept = default;
  //
  const Field<M>& operator=(const Field<M>& f)
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
    const Geometry& geo_v = geo();
    const int multiplicity = geo_v.multiplicity;
    Field<M>& f0 = *this;
    qacc_for(index, geo_v.local_volume(), {
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
  //
  qacc M& get_elem(const long offset)
  {
    qassert(0 <= offset && offset < (long)field.size());
    return field[offset];
  }
  qacc const M& get_elem(const long offset) const
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
    return get_elem(offset);
  }
  qacc const M& get_elem(const Coordinate& x, const int m) const
  {
    const Geometry& geo_v = geo();
    qassert(geo_v.is_on_node(x));
    qassert(0 <= m && m < geo_v.multiplicity);
    const long offset = geo_v.offset_from_coordinate(x) + m;
    return get_elem(offset);
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
      displayln("Field::get_elems_const: x=" + show(x) + "\ngeo=" + show(geo_v));
#endif
      qassert(false);
    }
    const long offset = geo_v.offset_from_coordinate(x);
    return Vector<M>(&field[offset], geo_v.multiplicity);
  }
  //
  qacc Vector<M> get_elems(const long index)
  // qassert(geo().is_only_local())
  {
    const Geometry& geo_v = geo();
    return Vector<M>(&field[index * geo_v.multiplicity], geo_v.multiplicity);
  }
  qacc Vector<M> get_elems_const(const long index) const
  // Be cautious about the const property
  // 改不改靠自觉
  // qassert(geo().is_only_local())
  {
    const Geometry& geo_v = geo();
    return Vector<M>(&field[index * geo_v.multiplicity], geo_v.multiplicity);
  }
};

inline int get_field_init_from_env()
{
  std::string tag = get_env_default("q_field_init", "fast");
  if (tag == "fast") {
    displayln_info("set q_field_init=fast.");
    return 0;  // do not do anything
  } else if (tag == "zero") {
    displayln_info("set q_field_init=zero.");
    return 1;  // set_zero
  } else if (tag == "random") {
    displayln_info("set q_field_init=random.");
    return 2;  // set rand
  } else {
    qassert(false);
    return -1;
  }
}

API inline int& get_field_init()
// qlat parameter
{
  static int t = get_field_init_from_env();
  return t;
}

template <class M, class N>
Field<M>& qcast(Field<N>& x)
// IMPORTANT: will modify the multiplicity of x, need to cast back after finish.
{
  if (x.initialized) {
    const int size = x.geo().multiplicity * sizeof(N);
    x.geo().multiplicity = size / sizeof(M);
    qassert(x.geo().multiplicity * sizeof(M) == size);
  }
  return (Field<M>&)x;
}

template <class M, class N>
const Field<M>& qcast_const(const Field<N>& x)
// IMPORTANT: will modify the multiplicity of x, need to cast back after finish.
{
  return qcast((Field<N>&)x);
}

template <class M>
bool is_initialized(const Field<M>& f)
{
  return f.initialized;
}

template <class M>
qacc Vector<M> get_data(const Field<M>& f)
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
  qassert(is_matching_geo_mult(f.geo(), f1.geo()));
  const Geometry& geo = f.geo();
  qacc_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<M> v1 = f1.get_elems_const(xl);
    Vector<M> v = f.get_elems(xl);
    for (int m = 0; m < geo.multiplicity; ++m) {
      v[m] += v1[m];
    }
  });
  return f;
}

template <class M>
const Field<M>& operator-=(Field<M>& f, const Field<M>& f1)
{
  TIMER("field_operator-=");
  if (not is_initialized(f)) {
    f.init(f1.geo());
    set_zero(f);
    f -= f1;
    return f;
  }
  qassert(is_matching_geo_mult(f.geo(), f1.geo()));
  const Geometry& geo = f.geo();
  qacc_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<M> v1 = f1.get_elems_const(xl);
    Vector<M> v = f.get_elems(xl);
    for (int m = 0; m < geo.multiplicity; ++m) {
      v[m] -= v1[m];
    }
  });
  return f;
}

template <class M>
const Field<M>& operator*=(Field<M>& f, const double factor)
{
  TIMER("field_operator*=(F,D)");
  const Geometry& geo = f.geo();
  qacc_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<M> v = f.get_elems(xl);
    for (int m = 0; m < geo.multiplicity; ++m) {
      v[m] *= factor;
    }
  });
  return f;
}

template <class M>
const Field<M>& operator*=(Field<M>& f, const Complex& factor)
{
  TIMER("field_operator*=(F,C)");
  const Geometry& geo = f.geo();
  qacc_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<M> v = f.get_elems(xl);
    for (int m = 0; m < geo.multiplicity; ++m) {
      v[m] *= factor;
    }
  });
  return f;
}

template <class M>
double qnorm(const Field<M>& f)
{
  const Geometry& geo = f.geo();
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
  const Geometry& geo = f1.geo();
  qassert(geo.is_only_local());
  qassert(geo == f2.geo());
  double sum = qnorm_double(get_data(f1), get_data(f2));
  glb_sum(sum);
  return sum;
}

template <class M, int multiplicity>
struct FieldM : Field<M> {
  using Field<M>::init;
  void init(const Geometry& geo_)
  {
    Field<M>::init(geo_, multiplicity);
  }
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
qacc long get_data_size(const Field<M>& f)
// NOT including the expended parts, only local volume data size
// only size on one node
{
  return f.geo().local_volume() * f.geo().multiplicity * sizeof(M);
}

template <class M>
void qswap(Field<M>& f1, Field<M>& f2)
{
  std::swap(f1.initialized, f2.initialized);
  qswap(f1.geo, f2.geo);
  qswap(f1.field, f2.field);
}

}  // namespace qlat
