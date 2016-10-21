#pragma once

#include <qlat/config.h>
#include <qlat/field.h>

QLAT_START_NAMESPACE

template<class M>
void field_sum(Vector<M> vec, const Field<M>& f)
{
  const int multiplicity = f.geo.multiplicity;
  assert(vec.size() == multiplicity);
  const Geometry& geo = f.geo;
  set_zero(vec);
  for (long index = 0; index < geo.local_volume(); ++index) {
    Coordinate x; geo.coordinate_from_index(x, index);
    const Vector<M> fvec = f.get_elems_const(x);
    for (int m = 0; m < multiplicity; ++m) {
      vec[m] += fvec[m];
    }
  }
}

template<class M>
inline void field_glb_sum_double(Vector<M> vec, const Field<M>& f)
{
  field_sum(vec, f);
  glb_sum_double(vec);
}

template<class M>
inline void field_glb_sum_long(Vector<M> vec, const Field<M>& f)
{
  field_sum(vec, f);
  glb_sum_long(vec);
}

QLAT_END_NAMESPACE
