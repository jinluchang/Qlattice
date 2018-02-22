#pragma once

#include <qlat/config.h>
#include <qlat/field.h>
#include <qlat/matrix.h>
#include <qlat/coordinate-d.h>

QLAT_START_NAMESPACE

// field.h and field-utils.h are including each other. Need this forward declaration.
// template <class M>
// struct Field;
// End of forward declaration.

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

template <class M>
std::vector<M> field_project_mom(const Field<M>& f, const CoordinateD& mom)
  // mom is in lattice unit (1/a)
  // project to component with momentum 'mom'
  // use glb_sum_double_vec to perform glb_sum
{
  TIMER("field_project_mom");
  const Geometry& geo = f.geo;
  std::vector<M> ret(geo.multiplicity);
  set_zero(ret);
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    double phase = 0;
    for (int k = 0; k < DIMN; ++k) {
      phase += mom[k] * xg[k];
    }
    const Vector<M> v = f.get_elems_const(xl);
    for (int m = 0; m < geo.multiplicity; ++m) {
      ret[m] += std::polar(1.0, -phase) * v[m];
    }
  }
  glb_sum_double_vec(get_data(ret));
  return ret;
}

template <class M>
std::vector<M> field_get_elems(const Field<M>& f, const Coordinate& xg)
{
  const Geometry& geo = f.geo;
  const Coordinate xl = geo.coordinate_l_from_g(xg);
  std::vector<M> ret(geo.multiplicity);
  if (geo.is_local(xl)) {
    assign(ret, f.get_elems_const(xl));
  } else {
    set_zero(ret);
  }
  glb_sum_byte_vec(get_data(ret));
  return ret;
}

template <class M>
M field_get_elem(const Field<M>& f, const Coordinate& xg, const int m)
{
  return field_get_elems(f, xg)[m];
}

template <class M>
M field_get_elem(const Field<M>& f, const Coordinate& xg)
{
  const Geometry& geo = f.geo;
  qassert(geo.multiplicity == 1);
  return field_get_elem(f, xg, 0);
}

QLAT_END_NAMESPACE
