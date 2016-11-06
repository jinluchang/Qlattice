#pragma once

#include <qlat/config.h>
#include <qlat/field.h>

QLAT_START_NAMESPACE

template <class M>
void set_zero(Field<M>& f)
{
  set_zero(f.field);
}

template <class M>
void set_unit(Field<M>& f)
{
  for (size_t offset = 0; offset < f.field.size(); ++offset) {
    set_unit(f.get_elem(offset));
  }
}

template<class M>
std::vector<M> field_sum(const Field<M>& f)
{
  const Geometry& geo = f.geo;
  const int multiplicity = geo.multiplicity;
  std::vector<M> vec(multiplicity);
  set_zero(vec);
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<M> fvec = f.get_elems_const(xl);
    for (int m = 0; m < multiplicity; ++m) {
      vec[m] += fvec[m];
    }
  }
  return vec;
}

template<class M>
std::vector<M> field_glb_sum_double(const Field<M>& f)
{
  std::vector<M> vec = field_sum(f);
  glb_sum_double_vec(Vector<M>(vec));
  return vec;
}

template<class M>
std::vector<M> field_glb_sum_long(const Field<M>& f)
{
  std::vector<M> vec = field_sum(f);
  glb_sum_long_vec(Vector<M>(vec));
  return vec;
}

QLAT_END_NAMESPACE
