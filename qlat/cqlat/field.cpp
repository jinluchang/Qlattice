#include "lib.h"

namespace qlat
{  //

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
