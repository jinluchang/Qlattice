#include "convert.h"
#include "dispatch.h"
#include "exceptions.h"

namespace qlat
{  //

template <class M>
PyObject* mk_field_ctype(Geometry* pgeo, const int multiplicity)
{
  Field<M>* p_f = new Field<M>();
  Field<M>& f = *p_f;
  if (pgeo != NULL) {
    Geometry& geo = *pgeo;
    if (multiplicity == 0) {
      f.init(geo);
    } else {
      f.init(geo, multiplicity);
    }
  }
  return PyLong_FromVoidPtr(p_f);
}

template <class M>
PyObject* free_field_ctype(void* p_field)
{
  Field<M>* p_f = (Field<M>*)p_field;
  delete p_f;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_zero_field_ctype(void* p_field)
{
  Field<M>* p_f = (Field<M>*)p_field;
  Field<M>& f = *p_f;
  set_zero(f);
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_unit_field_ctype(void* p_field, const Complex& coef = 1.0)
{
  Field<M>* p_f = (Field<M>*)p_field;
  Field<M>& f = *p_f;
  set_unit(f, coef);
  Py_RETURN_NONE;
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
  PyObject* p_field;
  field_dispatch(p_field, mk_field, ctype, pgeo, multiplicity);
  return p_field;
});

EXPORT(free_field, {
  using namespace qlat;
  PyObject* p_ctype = NULL;
  void* p_field = NULL;
  if (!PyArg_ParseTuple(args, "Ol", &p_ctype, &p_field)) {
    return NULL;
  }
  std::string ctype;
  py_convert(ctype, p_ctype);
  PyObject* p_ret = NULL;
  field_dispatch(p_ret, free_field, ctype, p_field);
  return p_ret;
});

EXPORT(set_zero_field, {
  using namespace qlat;
  PyObject* p_ctype = NULL;
  void* p_field = NULL;
  if (!PyArg_ParseTuple(args, "Ol", &p_ctype, &p_field)) {
    return NULL;
  }
  std::string ctype;
  py_convert(ctype, p_ctype);
  PyObject* p_ret = NULL;
  field_dispatch(p_ret, set_zero_field, ctype, p_field);
  return p_ret;
});

EXPORT(set_unit_field, {
  using namespace qlat;
  PyObject* p_ctype = NULL;
  void* p_field = NULL;
  Complex coef = 1.0;
  if (!PyArg_ParseTuple(args, "Ol|D", &p_ctype, &p_field, &coef)) {
    return NULL;
  }
  std::string ctype;
  py_convert(ctype, p_ctype);
  PyObject* p_ret = NULL;
  field_dispatch(p_ret, set_unit_field, ctype, p_field, coef);
  return p_ret;
});
