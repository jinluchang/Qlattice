#include "lib.h"

EXPORT(get_elem_wm_prop, {
  using namespace qlat;
  PyObject* p_wm = NULL;
  PyObject* p_obj = NULL;
  PyObject* p_xg = NULL;
  long m = 0;
  if (!PyArg_ParseTuple(args, "OOO|l", &p_wm, &p_obj, &p_xg, &m)) {
    return NULL;
  }
  WilsonMatrix& wm = py_convert_type<WilsonMatrix>(p_wm);
  const Field<WilsonMatrix>& obj = py_convert_type_field<WilsonMatrix>(p_obj);
  const Coordinate xg = py_convert_data<Coordinate>(p_xg);
  wm = obj.get_elem(xg, m);
  Py_RETURN_NONE;
})

EXPORT(get_elem_wm_sprop, {
  using namespace qlat;
  PyObject* p_wm = NULL;
  PyObject* p_obj = NULL;
  long idx = 0;
  long m = 0;
  if (!PyArg_ParseTuple(args, "OOl|l", &p_wm, &p_obj, &idx, &m)) {
    return NULL;
  }
  WilsonMatrix& wm = py_convert_type<WilsonMatrix>(p_wm);
  const SelectedField<WilsonMatrix>& obj =
      py_convert_type_sfield<WilsonMatrix>(p_obj);
  wm = obj.get_elem(idx, m);
  Py_RETURN_NONE;
})

EXPORT(get_elem_wm_psprop, {
  using namespace qlat;
  PyObject* p_wm = NULL;
  PyObject* p_obj = NULL;
  long idx = 0;
  long m = 0;
  if (!PyArg_ParseTuple(args, "OOl|l", &p_wm, &p_obj, &idx, &m)) {
    return NULL;
  }
  WilsonMatrix& wm = py_convert_type<WilsonMatrix>(p_wm);
  const SelectedPoints<WilsonMatrix>& obj =
      py_convert_type_spoints<WilsonMatrix>(p_obj);
  wm = obj.get_elem(idx, m);
  Py_RETURN_NONE;
})
