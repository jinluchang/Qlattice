#include "lib.h"

namespace qlat
{  //

template <class M>
PyObject* mk_field_ctype(PyObject* p_geo, const int multiplicity)
{
  Field<M>* p_field = new Field<M>();
  Field<M>& f = *p_field;
  if (p_geo != NULL) {
    const Geometry& geo = py_convert_type<Geometry>(p_geo);
    if (multiplicity == 0) {
      f.init(geo);
    } else {
      pqassert(multiplicity > 0);
      f.init(geo, multiplicity);
    }
  }
  return py_convert((void*)p_field);
}

template <class M>
PyObject* free_field_ctype(PyObject* p_field)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  delete &f;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_field_ctype(PyObject* p_field_new, PyObject* p_field)
{
  Field<M>& f_new = py_convert_type_field<M>(p_field_new);
  const Field<M>& f = py_convert_type_field<M>(p_field);
  f_new = f;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_add_field_ctype(PyObject* p_field_new, PyObject* p_field)
{
  Field<M>& f_new = py_convert_type_field<M>(p_field_new);
  const Field<M>& f = py_convert_type_field<M>(p_field);
  f_new += f;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_sub_field_ctype(PyObject* p_field_new, PyObject* p_field)
{
  Field<M>& f_new = py_convert_type_field<M>(p_field_new);
  const Field<M>& f = py_convert_type_field<M>(p_field);
  f_new -= f;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_mul_field_ctype(PyObject* p_field, const double& factor)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  f *= factor;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_mul_field_ctype(PyObject* p_field, const Complex& factor)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  f *= factor;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_mul_field_ctype(PyObject* p_field, const FieldM<Complex, 1>& f_factor)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  f *= f_factor;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_zero_field_ctype(PyObject* p_field)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  set_zero(f);
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_unit_field_ctype(PyObject* p_field, const Complex& coef)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  set_unit(f, coef);
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_u_rand_double_field_ctype(PyObject* p_field, const RngState& rs,
                                        const double upper, const double lower)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  set_u_rand_double(f, rs, upper, lower);
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_g_rand_double_field_ctype(PyObject* p_field, const RngState& rs,
                                        const double center, const double sigma)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  set_g_rand_double(f, rs, center, sigma);
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_checkers_double_field_ctype(PyObject* p_field)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  set_checkers_double(f);
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_complex_from_double_field_ctype(PyObject* p_field, PyObject* p_sf)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  Field<double>& sf = py_convert_type_field<double>(p_sf);
  set_complex_from_double(f, sf);
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_double_from_complex_field_ctype(PyObject* p_field, PyObject* p_cf)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  Field<Complex>& cf = py_convert_type_field<Complex>(p_cf);
  set_double_from_complex(f, cf);
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_abs_from_complex_field_ctype(PyObject* p_field, PyObject* p_cf)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  Field<Complex>& cf = py_convert_type_field<Complex>(p_cf);
  set_abs_from_complex(f, cf);
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_ratio_double_field_ctype(PyObject* p_field, PyObject* p_sf1, PyObject* p_sf2)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  Field<double>& sf1 = py_convert_type_field<double>(p_sf1);
  Field<double>& sf2 = py_convert_type_field<double>(p_sf2);
  set_ratio_double(f, sf1, sf2);
  Py_RETURN_NONE;
}

template <class M>
PyObject* less_than_double_field_ctype(PyObject* p_sf1, PyObject* p_sf2, PyObject* p_mask)
{
  Field<M>& sf1 = py_convert_type_field<M>(p_sf1);
  Field<double>& sf2 = py_convert_type_field<double>(p_sf2);
  Field<double>& mask = py_convert_type_field<double>(p_mask);
  less_than_double(sf1, sf2, mask);
  Py_RETURN_NONE;
}

template <class M>
PyObject* invert_double_field_ctype(PyObject* p_sf)
{
  Field<M>& sf = py_convert_type_field<M>(p_sf);
  invert_double(sf);
  Py_RETURN_NONE;
}

template <class M>
PyObject* multiply_double_field_ctype(PyObject* p_sf, PyObject* p_factor)
{
  Field<M>& sf = py_convert_type_field<M>(p_sf);
  Field<double>& factor = py_convert_type_field<double>(p_factor);
  multiply_double(sf, factor);
  Py_RETURN_NONE;
}

template <class M>
PyObject* get_total_site_field_ctype(PyObject* p_field)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  const Coordinate ret = f.geo().total_site();
  return py_convert(ret);
}

template <class M>
PyObject* get_multiplicity_field_ctype(PyObject* p_field)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  const long ret = f.geo().multiplicity;
  return py_convert(ret);
}

template <class M>
PyObject* get_sizeof_m_field_ctype(PyObject* p_field)
{
  (void)p_field;
  return py_convert(sizeof(M));
}

template <class M>
PyObject* set_geo_field_ctype(Geometry& geo, PyObject* p_field)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  geo = f.geo();
  Py_RETURN_NONE;
}

template <class M>
PyObject* qnorm_field_ctype(PyObject* p_field)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  const double ret = qnorm(f);
  return py_convert(ret);
}

template <class M>
PyObject* crc32_field_ctype(PyObject* p_field)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  const crc32_t ret = field_crc32(f);
  return py_convert((long)ret);
}

template <class M>
PyObject* get_mview_field_ctype(PyObject* p_field)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  Vector<M> fv = get_data(f);
  PyObject* p_mview = py_convert_mview(fv);
  pqassert(p_field != NULL);
  Py_INCREF(p_field);
  pqassert(!((PyMemoryViewObject*)p_mview)->mbuf->master.obj);
  ((PyMemoryViewObject*)p_mview)->mbuf->master.obj = p_field;
  return p_mview;
}

} // namespace qlat

EXPORT(mk_field, {
  using namespace qlat;
  PyObject* p_ctype = NULL;
  PyObject* p_geo = NULL;
  int multiplicity = 0;
  if (!PyArg_ParseTuple(args, "O|Oi", &p_ctype, &p_geo, &multiplicity)) {
    return NULL;
  }
  const std::string ctype = py_convert_data<std::string>(p_ctype);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, mk_field_ctype, ctype, p_geo, multiplicity);
  return p_ret;
})

EXPORT(free_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, free_field_ctype, ctype, p_field);
  return p_ret;
})

EXPORT(set_field, {
  using namespace qlat;
  PyObject* p_field_new = NULL;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field_new, &p_field)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  pqassert(py_get_ctype(p_field_new) == ctype);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_field_ctype, ctype, p_field_new, p_field);
  return p_ret;
})

EXPORT(set_add_field, {
  using namespace qlat;
  PyObject* p_field_new = NULL;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field_new, &p_field)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  pqassert(py_get_ctype(p_field_new) == ctype);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_add_field_ctype, ctype, p_field_new, p_field);
  return p_ret;
})

EXPORT(set_sub_field, {
  using namespace qlat;
  PyObject* p_field_new = NULL;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field_new, &p_field)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  pqassert(py_get_ctype(p_field_new) == ctype);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_sub_field_ctype, ctype, p_field_new, p_field);
  return p_ret;
})

EXPORT(set_mul_double_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  double factor = 0.0;
  if (!PyArg_ParseTuple(args, "Od", &p_field, &factor)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_mul_field_ctype, ctype, p_field, factor);
  return p_ret;
})

EXPORT(set_mul_complex_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  Complex factor = 0.0;
  if (!PyArg_ParseTuple(args, "OD", &p_field, &factor)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_mul_field_ctype, ctype, p_field, factor);
  return p_ret;
})

EXPORT(set_mul_cfield_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_cfield = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field, &p_cfield)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  FieldM<Complex, 1>& f_factor = py_convert_type_field<Complex, 1>(p_cfield);
  pqassert(f_factor.geo().multiplicity == 1);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_mul_field_ctype, ctype, p_field, f_factor);
  return p_ret;
})

EXPORT(set_zero_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_zero_field_ctype, ctype, p_field);
  return p_ret;
})

EXPORT(set_unit_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  Complex coef = 1.0;
  if (!PyArg_ParseTuple(args, "O|D", &p_field, &coef)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_unit_field_ctype, ctype, p_field, coef);
  return p_ret;
})

EXPORT(set_u_rand_double_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_rng = NULL;
  double upper = 1.0;
  double lower = 0.0;
  if (!PyArg_ParseTuple(args, "OO|dd", &p_field, &p_rng, &upper, &lower)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  const RngState& rng = py_convert_type<RngState>(p_rng);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_u_rand_double_field_ctype, ctype, p_field, rng,
                 upper, lower);
  return p_ret;
})

EXPORT(set_g_rand_double_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_rng = NULL;
  double center = 0.0;
  double sigma = 1.0;
  if (!PyArg_ParseTuple(args, "OO|dd", &p_field, &p_rng, &center, &sigma)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  const RngState& rng = py_convert_type<RngState>(p_rng);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_g_rand_double_field_ctype, ctype, p_field, rng,
                 center, sigma);
  return p_ret;
})

EXPORT(set_checkers_double_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_checkers_double_field_ctype, ctype, p_field);
  return p_ret;
})

EXPORT(set_complex_from_double_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_sf = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field, &p_sf)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_complex_from_double_field_ctype, ctype, p_field, p_sf);
  return p_ret;
})

EXPORT(set_double_from_complex_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_cf = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field, &p_cf)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_double_from_complex_field_ctype, ctype, p_field, p_cf);
  return p_ret;
})

EXPORT(set_abs_from_complex_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_cf = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field, &p_cf)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_abs_from_complex_field_ctype, ctype, p_field, p_cf);
  return p_ret;
})

EXPORT(set_ratio_double_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_sf1 = NULL;
  PyObject* p_sf2 = NULL;
  if (!PyArg_ParseTuple(args, "OOO", &p_field, &p_sf1, &p_sf2)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_ratio_double_field_ctype, ctype, p_field, p_sf1, p_sf2);
  return p_ret;
})

EXPORT(less_than_double_field, {
  using namespace qlat;
  PyObject* p_sf1 = NULL;
  PyObject* p_sf2 = NULL;
  PyObject* p_mask = NULL;
  if (!PyArg_ParseTuple(args, "OOO", &p_sf1, &p_sf2, &p_mask)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_sf1);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, less_than_double_field_ctype, ctype, p_sf1, p_sf2, p_mask);
  return p_ret;
})

EXPORT(invert_double_field, {
  using namespace qlat;
  PyObject* p_sf = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_sf)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_sf);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, invert_double_field_ctype, ctype, p_sf);
  return p_ret;
})

EXPORT(multiply_double_field, {
  using namespace qlat;
  PyObject* p_sf = NULL;
  PyObject* p_factor = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_sf, &p_factor)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_sf);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, multiply_double_field_ctype, ctype, p_sf, p_factor);
  return p_ret;
})

EXPORT(get_total_site_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, get_total_site_field_ctype, ctype, p_field);
  return p_ret;
})

EXPORT(get_multiplicity_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, get_multiplicity_field_ctype, ctype, p_field);
  return p_ret;
})

EXPORT(get_sizeof_m_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, get_sizeof_m_field_ctype, ctype, p_field);
  return p_ret;
})

EXPORT(set_geo_field, {
  using namespace qlat;
  PyObject* p_geo = NULL;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_geo, &p_field)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  Geometry& geo = py_convert_type<Geometry>(p_geo);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_geo_field_ctype, ctype, geo, p_field);
  return p_ret;
})

EXPORT(qnorm_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, qnorm_field_ctype, ctype, p_field);
  return p_ret;
})

EXPORT(crc32_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, crc32_field_ctype, ctype, p_field);
  return p_ret;
})

EXPORT(get_mview_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, get_mview_field_ctype, ctype, p_field);
  return p_ret;
})
