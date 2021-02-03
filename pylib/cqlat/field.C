#include "lib.h"

namespace qlat
{  //

template <class M>
PyObject* mk_field_ctype(PyObject* p_geo, const int multiplicity)
{
  Field<M>* pf = new Field<M>();
  Field<M>& f = *pf;
  if (p_geo != NULL) {
    const Geometry& geo = py_convert_type<Geometry>(p_geo);
    if (multiplicity == 0) {
      f.init(geo);
    } else {
      f.init(geo, multiplicity);
    }
  }
  return py_convert((void*)pf);
}

template <class M>
PyObject* free_field_ctype(PyField& pf)
{
  Field<M>& f = *(Field<M>*)pf.cdata;
  delete &f;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_zero_field_ctype(PyField& pf)
{
  Field<M>& f = *(Field<M>*)pf.cdata;
  set_zero(f);
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_unit_field_ctype(PyField& pf, const Complex& coef = 1.0)
{
  Field<M>& f = *(Field<M>*)pf.cdata;
  set_unit(f, coef);
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_geo_field_ctype(Geometry& geo, void* pfield)
{
  Field<M>* pf = (Field<M>*)pfield;
  Field<M>& f = *pf;
  geo = f.geo();
  Py_RETURN_NONE;
}

template <class M>
PyObject* get_mview_field_ctype(void* pfield, PyObject* p_field)
{
  Field<M>* pf = (Field<M>*)pfield;
  Field<M>& f = *pf;
  Vector<M> fv = get_data(f);
  PyObject* p_mview = py_convert(fv);
  pqassert(p_field != NULL);
  Py_INCREF(p_field);
  pqassert(!((PyMemoryViewObject*)p_mview)->mbuf->master.obj);
  ((PyMemoryViewObject*)p_mview)->mbuf->master.obj = p_field;
  return p_mview;
}

}  // namespace qlat

EXPORT(mk_field, {
  using namespace qlat;
  PyObject* p_ctype = NULL;
  PyObject* p_geo = NULL;
  int multiplicity = 0;
  if (!PyArg_ParseTuple(args, "O|Oi", &p_ctype, &p_geo, &multiplicity)) {
    return NULL;
  }
  std::string ctype;
  py_convert(ctype, p_ctype);
  PyObject* pfield;
  FIELD_DISPATCH(pfield, mk_field_ctype, ctype, p_geo, multiplicity);
  return pfield;
});

EXPORT(free_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, free_field_ctype, pf.ctype, pf);
  return p_ret;
});

EXPORT(set_zero_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_zero_field_ctype, pf.ctype, pf);
  return p_ret;
});

EXPORT(set_unit_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  Complex coef = 1.0;
  if (!PyArg_ParseTuple(args, "O|D", &p_field, &coef)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_unit_field_ctype, pf.ctype, pf, coef);
  return p_ret;
});

EXPORT(set_geo_field, {
  using namespace qlat;
  Geometry* pgeo = NULL;
  PyObject* p_ctype = NULL;
  void* pfield = NULL;
  if (!PyArg_ParseTuple(args, "lOl", &pgeo, &p_ctype, &pfield)) {
    return NULL;
  }
  Geometry& geo = *pgeo;
  std::string ctype;
  py_convert(ctype, p_ctype);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_geo_field_ctype, ctype, geo, pfield);
  return p_ret;
});

EXPORT(get_mview_field, {
  using namespace qlat;
  PyObject* p_ctype = NULL;
  void* pfield = NULL;
  PyObject* p_field= NULL;
  if (!PyArg_ParseTuple(args, "OlO", &p_ctype, &pfield, &p_field)) {
    return NULL;
  }
  std::string ctype;
  py_convert(ctype, p_ctype);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, get_mview_field_ctype, ctype, pfield, p_field);
  return p_ret;
});
