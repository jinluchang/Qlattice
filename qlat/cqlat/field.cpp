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
PyObject* set_mul_field_ctype(PyObject* p_field,
                              const Field<ComplexD>& f_factor)
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

}  // namespace qlat

EXPORT(set_mul_complex_field, {  // tested: cqlat-fields
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

EXPORT(get_mview_field, {  // tested: cqlat-fields
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
