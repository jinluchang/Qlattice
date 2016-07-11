#pragma once

#include <qlat/config.h>
#include <qlat/field.h>

QLAT_START_NAMESPACE

template<class M>
void fieldSum(Vector<M> vec, const Field<M>& f)
{
  const int multiplicity = f.geo.multiplicity;
  assert(vec.size() == multiplicity);
  const Geometry& geo = f.geo;
  setZero(vec);
  for (long index = 0; index < geo.localVolume(); ++index) {
    Coordinate x; geo.coordinateFromIndex(x, index);
    const Vector<M> fvec = f.getElemsConst(x);
    for (int m = 0; m < multiplicity; ++m) {
      vec[m] += fvec[m];
    }
  }
}

template<class M>
inline void fieldGlbSumDouble(Vector<M> vec, const Field<M>& f)
{
  fieldSum(vec, f);
  glbSumDouble(vec);
}

template<class M>
inline void fieldGlbSumLong(Vector<M> vec, const Field<M>& f)
{
  fieldSum(vec, f);
  glbSumLong(vec);
}

QLAT_END_NAMESPACE
