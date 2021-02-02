#include "convert.h"
#include "dispatch.h"
#include "exceptions.h"

namespace qlat
{  //

template <class M>
PyObject* mk_field_ctype(Geometry* pgeo, const int multiplicity)
{
  Field<M>* pf = new Field<M>();
  Field<M>& f = *pf;
  if (pgeo != NULL) {
    Geometry& geo = *pgeo;
    if (multiplicity == 0) {
      f.init(geo);
    } else {
      f.init(geo, multiplicity);
    }
  }
  return py_convert((void*)pf);
}

template <class M>
PyObject* free_field_ctype(void* pfield)
{
  Field<M>* pf = (Field<M>*)pfield;
  delete pf;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_zero_field_ctype(void* pfield)
{
  Field<M>* pf = (Field<M>*)pfield;
  Field<M>& f = *pf;
  set_zero(f);
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_unit_field_ctype(void* pfield, const Complex& coef = 1.0)
{
  Field<M>* pf = (Field<M>*)pfield;
  Field<M>& f = *pf;
  set_unit(f, coef);
  Py_RETURN_NONE;
}

template <class M>
PyObject* get_geo_field_ctype(void* pfield)
{
  Field<M>* pf = (Field<M>*)pfield;
  Field<M>& f = *pf;
  Geometry* pgeo = new Geometry(f.geo());
  return py_convert((void*)pgeo);
}

template <class M>
PyObject* get_mview_field_ctype(void* pfield)
{
  Field<M>* pf = (Field<M>*)pfield;
  Field<M>& f = *pf;
  Vector<M> fv = get_data(f);
  return py_convert(fv);
}

}  // namespace qlat

EXPORT(mk_field, {
  using namespace qlat;
  PyObject* p_ctype = NULL;
  Geometry* pgeo = NULL;
  int multiplicity = 0;
  if (!PyArg_ParseTuple(args, "O|li", &p_ctype, &pgeo, &multiplicity)) {
    return NULL;
  }
  std::string ctype;
  py_convert(ctype, p_ctype);
  PyObject* pfield;
  FIELD_DISPATCH(pfield, mk_field_ctype, ctype, pgeo, multiplicity);
  return pfield;
});

EXPORT(free_field, {
  using namespace qlat;
  PyObject* p_ctype = NULL;
  void* pfield = NULL;
  if (!PyArg_ParseTuple(args, "Ol", &p_ctype, &pfield)) {
    return NULL;
  }
  std::string ctype;
  py_convert(ctype, p_ctype);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, free_field_ctype, ctype, pfield);
  return p_ret;
});

EXPORT(set_zero_field, {
  using namespace qlat;
  PyObject* p_ctype = NULL;
  void* pfield = NULL;
  if (!PyArg_ParseTuple(args, "Ol", &p_ctype, &pfield)) {
    return NULL;
  }
  std::string ctype;
  py_convert(ctype, p_ctype);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_zero_field_ctype, ctype, pfield);
  return p_ret;
});

EXPORT(set_unit_field, {
  using namespace qlat;
  PyObject* p_ctype = NULL;
  void* pfield = NULL;
  Complex coef = 1.0;
  if (!PyArg_ParseTuple(args, "Ol|D", &p_ctype, &pfield, &coef)) {
    return NULL;
  }
  std::string ctype;
  py_convert(ctype, p_ctype);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_unit_field_ctype, ctype, pfield, coef);
  return p_ret;
});

EXPORT(get_geo_field, {
  using namespace qlat;
  PyObject* p_ctype = NULL;
  void* pfield = NULL;
  if (!PyArg_ParseTuple(args, "Ol", &p_ctype, &pfield)) {
    return NULL;
  }
  std::string ctype;
  py_convert(ctype, p_ctype);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, get_geo_field_ctype, ctype, pfield);
  return p_ret;
});

EXPORT(get_mview_field, {
  using namespace qlat;
  PyObject* p_ctype = NULL;
  void* pfield = NULL;
  if (!PyArg_ParseTuple(args, "Ol", &p_ctype, &pfield)) {
    return NULL;
  }
  std::string ctype;
  py_convert(ctype, p_ctype);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, get_mview_field_ctype, ctype, pfield);
  return p_ret;
});
