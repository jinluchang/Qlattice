#pragma once

#include <qlat/setup.h>
#include <qlat/geometry.h>
#include <qlat/mpi.h>

#include <qlat-utils/vector.h>

#include <ctime>
#include <fstream>
#include <vector>

namespace qlat
{  //

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
  qassert(geo.is_only_local);
  qassert(geo == f2.geo());
  double sum = qnorm_double(get_data(f1), get_data(f2));
  glb_sum(sum);
  return sum;
}

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
