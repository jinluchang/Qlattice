#include "lib.h"

namespace qlat
{  //

template <class M>
PyObject* mk_field_ctype(PyObject* p_geo, const Int multiplicity)
{
  Field<M>* p_field = new Field<M>();
  Field<M>& f = *p_field;
  if (p_geo != NULL) {
    const Geometry& geo = py_convert_type<Geometry>(p_geo);
    if (multiplicity == 0) {
      f.init(geo);
    } else {
      qassert(multiplicity > 0);
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
PyObject* set_mul_field_ctype(PyObject* p_field, const RealD& factor)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  f *= factor;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_mul_field_ctype(PyObject* p_field, const ComplexD& factor)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  f *= factor;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_mul_field_ctype(PyObject* p_field, const Field<ComplexD>& f_factor)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  f *= f_factor;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_mul_field_ctype(PyObject* p_field, const Field<RealD>& f_factor)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  f *= f_factor;
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
  const Long ret = f.multiplicity;
  return py_convert(ret);
}

template <class M>
PyObject* get_sizeof_m_field_ctype(PyObject* p_field)
{
  (void)p_field;
  const Long size = sizeof(M);
  return py_convert(size);
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
  const RealD ret = qnorm(f);
  return py_convert(ret);
}

template <class M>
PyObject* crc32_field_ctype(PyObject* p_field)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  const crc32_t ret = field_crc32(f);
  return py_convert((Long)ret);
}

template <class M>
PyObject* get_mview_field_ctype(PyObject* p_field)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  Vector<M> fv = get_data(f);
  PyObject* p_mview = py_convert_mview(fv);
  qassert(p_field != NULL);
  Py_INCREF(p_field);
  qassert(!((PyMemoryViewObject*)p_mview)->mbuf->master.obj);
  ((PyMemoryViewObject*)p_mview)->mbuf->master.obj = p_field;
  return p_mview;
}

} // namespace qlat

EXPORT(set_add_field, {
  using namespace qlat;
  PyObject* p_field_new = NULL;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field_new, &p_field)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  qassert(py_get_ctype(p_field_new) == ctype);
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
  qassert(py_get_ctype(p_field_new) == ctype);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_sub_field_ctype, ctype, p_field_new, p_field);
  return p_ret;
})

EXPORT(set_mul_double_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  RealD factor = 0.0;
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
  ComplexD factor = 0.0;
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
  const std::string ctype_c = py_get_ctype(p_cfield);
  PyObject* p_ret = NULL;
  if (ctype_c == "ComplexD") {
    Field<ComplexD>& f_factor = py_convert_type_field<ComplexD>(p_cfield);
    FIELD_DISPATCH(p_ret, set_mul_field_ctype, ctype, p_field, f_factor);
  } else if (ctype_c == "RealD") {
    Field<RealD>& f_factor = py_convert_type_field<RealD>(p_cfield);
    FIELD_DISPATCH(p_ret, set_mul_field_ctype, ctype, p_field, f_factor);
  } else {
    qassert(false);
  }
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
