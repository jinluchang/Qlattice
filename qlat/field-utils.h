#pragma once

#include <qlat/config.h>
#include <qlat/coordinate-d.h>
#include <qlat/field.h>
#include <qlat/matrix.h>

QLAT_START_NAMESPACE

// field.h and field-utils.h are including each other. Need this forward
// declaration. template <class M> struct Field; End of forward declaration.

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

template <class M>
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

template <class M>
std::vector<M> field_glb_sum_double(const Field<M>& f)
{
  std::vector<M> vec = field_sum(f);
  glb_sum_double_vec(Vector<M>(vec));
  return vec;
}

template <class M>
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

template <class M>
void field_shift_dir(Field<M>& f, const Field<M>& f1, const int dir,
                     const int shift)
// shift f1 in 'dir' direction for 'shift' steps
{
  TIMER("field_shift_dir");
  qassert(0 <= dir and dir < 4);
  const Geometry geo = geo_resize(f1.geo);
  const Coordinate total_site = geo.total_site();
  f.init(geo);
  qassert(is_matching_geo_mult(f.geo, f1.geo));
  Coordinate nvec;
  nvec[dir] = 1;
  Field<M> tmp, tmp1;
  tmp.init(geo);
  tmp1.init(geo);
  tmp1 = f1;
  for (int i = 0; i < geo.geon.size_node[dir]; ++i) {
#pragma omp parallel for
    for (long index = 0; index < geo.local_volume(); ++index) {
      const Coordinate xl = geo.coordinate_from_index(index);
      const Coordinate xg = geo.coordinate_g_from_l(xl);
      const Coordinate xg1 =
          mod(xg - (shift + i * geo.node_site[dir]) * nvec, total_site);
      const Coordinate xl1 = geo.coordinate_l_from_g(xg1);
      if (geo.is_local(xl1)) {
        assign(f.get_elems(xl), tmp1.get_elems_const(xl1));
      }
    }
    if (i < geo.geon.size_node[dir] - 1) {
      get_data_plus_mu(get_data(tmp), get_data(tmp1), dir);
      std::swap(tmp1, tmp);
    }
  }
}

template <class M>
void field_shift(Field<M>& f, const Field<M>& f1, const Coordinate& shift)
// shift f1 with 'shift'
{
  TIMER("field_shift");
  Field<M> tmp, tmp1;
  field_shift_dir(tmp, f1, 0, shift[0]);
  field_shift_dir(tmp1, tmp, 1, shift[1]);
  field_shift_dir(tmp, tmp1, 2, shift[2]);
  field_shift_dir(f, tmp, 3, shift[3]);
}

template <class M>
void set_field_u_rand_double(Field<M>& f, const RngState& rs,
                             const double upper = 1.0,
                             const double lower = -1.0)
{
  TIMER("set_field_u_rand_double");
  const Geometry& geo = f.geo;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    RngState rsi(rs, show(xg));
    Vector<WilsonVector> v = f.get_elems(xl);
    Vector<double> dv((double*)v.data(), v.data_size() / sizeof(double));
    for (int m = 0; m < dv.size(); ++m) {
      dv[m] = u_rand_gen(rsi, 1.0, -1.0);
    }
  }
}

QLAT_END_NAMESPACE
