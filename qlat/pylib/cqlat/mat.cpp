#include "lib.h"

EXPORT(mk_wilson_matrix, {
  using namespace qlat;
  WilsonMatrix* p = new WilsonMatrix();
  return py_convert((void*)p);
})

EXPORT(free_wilson_matrix, {
  using namespace qlat;
  return free_obj<WilsonMatrix>(args);
})

EXPORT(set_wilson_matrix, {
  using namespace qlat;
  return set_obj<WilsonMatrix>(args);
})

EXPORT(set_zero_wilson_matrix, {
  using namespace qlat;
  PyObject* p_obj = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_obj)) {
    return NULL;
  }
  WilsonMatrix& obj = py_convert_type<WilsonMatrix>(p_obj);
  set_zero(obj);
  Py_RETURN_NONE;
})

// -----------------------------------------

EXPORT(mk_spin_matrix, {
  using namespace qlat;
  SpinMatrix* p = new SpinMatrix();
  return py_convert((void*)p);
})

EXPORT(free_spin_matrix, {
  using namespace qlat;
  return free_obj<SpinMatrix>(args);
})

EXPORT(set_spin_matrix, {
  using namespace qlat;
  return set_obj<SpinMatrix>(args);
})

EXPORT(set_zero_spin_matrix, {
  using namespace qlat;
  PyObject* p_obj = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_obj)) {
    return NULL;
  }
  SpinMatrix& obj = py_convert_type<SpinMatrix>(p_obj);
  set_zero(obj);
  Py_RETURN_NONE;
})

// -----------------------------------------

EXPORT(set_g5_herm_wilson_matrix, {
  using namespace qlat;
  PyObject* p_obj = NULL;
  PyObject* p_obj1 = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_obj, &p_obj1)) {
    return NULL;
  }
  WilsonMatrix& obj = py_convert_type<WilsonMatrix>(p_obj);
  const WilsonMatrix& obj1 = py_convert_type<WilsonMatrix>(p_obj1);
  const box_acc<SpinMatrixConstants>& smc = get_spin_matrix_constants();
  const SpinMatrix& gamma5 = smc().gamma5;
  obj = gamma5 * (WilsonMatrix)matrix_adjoint(obj1) * gamma5;
  Py_RETURN_NONE;
})

EXPORT(set_wm_mul_wm_wm, {
  using namespace qlat;
  PyObject* p_obj = NULL;
  PyObject* p_obj1 = NULL;
  PyObject* p_obj2 = NULL;
  if (!PyArg_ParseTuple(args, "OOO", &p_obj, &p_obj1, &p_obj2)) {
    return NULL;
  }
  WilsonMatrix& obj = py_convert_type<WilsonMatrix>(p_obj);
  const WilsonMatrix& obj1 = py_convert_type<WilsonMatrix>(p_obj1);
  const WilsonMatrix& obj2 = py_convert_type<WilsonMatrix>(p_obj2);
  obj = obj1 * obj2;
  Py_RETURN_NONE;
})

EXPORT(set_wm_mul_wm_sm, {
  using namespace qlat;
  PyObject* p_obj = NULL;
  PyObject* p_obj1 = NULL;
  PyObject* p_obj2 = NULL;
  if (!PyArg_ParseTuple(args, "OOO", &p_obj, &p_obj1, &p_obj2)) {
    return NULL;
  }
  WilsonMatrix& obj = py_convert_type<WilsonMatrix>(p_obj);
  const WilsonMatrix& obj1 = py_convert_type<WilsonMatrix>(p_obj1);
  const SpinMatrix& obj2 = py_convert_type<SpinMatrix>(p_obj2);
  obj = obj1 * obj2;
  Py_RETURN_NONE;
})

EXPORT(set_wm_mul_sm_wm, {
  using namespace qlat;
  PyObject* p_obj = NULL;
  PyObject* p_obj1 = NULL;
  PyObject* p_obj2 = NULL;
  if (!PyArg_ParseTuple(args, "OOO", &p_obj, &p_obj1, &p_obj2)) {
    return NULL;
  }
  WilsonMatrix& obj = py_convert_type<WilsonMatrix>(p_obj);
  const SpinMatrix& obj1 = py_convert_type<SpinMatrix>(p_obj1);
  const WilsonMatrix& obj2 = py_convert_type<WilsonMatrix>(p_obj2);
  obj = obj1 * obj2;
  Py_RETURN_NONE;
})

EXPORT(set_sm_mul_sm_sm, {
  using namespace qlat;
  PyObject* p_obj = NULL;
  PyObject* p_obj1 = NULL;
  PyObject* p_obj2 = NULL;
  if (!PyArg_ParseTuple(args, "OOO", &p_obj, &p_obj1, &p_obj2)) {
    return NULL;
  }
  SpinMatrix& obj = py_convert_type<SpinMatrix>(p_obj);
  const SpinMatrix& obj1 = py_convert_type<SpinMatrix>(p_obj1);
  const SpinMatrix& obj2 = py_convert_type<SpinMatrix>(p_obj2);
  obj = obj1 * obj2;
  Py_RETURN_NONE;
})

EXPORT(set_wm_mul_a_wm, {
  using namespace qlat;
  PyObject* p_obj = NULL;
  Complex a = 0;
  PyObject* p_obj2 = NULL;
  if (!PyArg_ParseTuple(args, "ODO", &p_obj, &a, &p_obj2)) {
    return NULL;
  }
  WilsonMatrix& obj = py_convert_type<WilsonMatrix>(p_obj);
  const WilsonMatrix& obj2 = py_convert_type<WilsonMatrix>(p_obj2);
  obj = a * obj2;
  Py_RETURN_NONE;
})

EXPORT(set_sm_mul_a_sm, {
  using namespace qlat;
  PyObject* p_obj = NULL;
  Complex a = 0;
  PyObject* p_obj2 = NULL;
  if (!PyArg_ParseTuple(args, "ODO", &p_obj, &a, &p_obj2)) {
    return NULL;
  }
  SpinMatrix& obj = py_convert_type<SpinMatrix>(p_obj);
  const SpinMatrix& obj2 = py_convert_type<SpinMatrix>(p_obj2);
  obj = a * obj2;
  Py_RETURN_NONE;
})

EXPORT(trace_wm, {
  using namespace qlat;
  PyObject* p_obj = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_obj)) {
    return NULL;
  }
  const WilsonMatrix& obj = py_convert_type<WilsonMatrix>(p_obj);
  return py_convert(matrix_trace(obj));
})

EXPORT(trace_wm_wm, {
  using namespace qlat;
  PyObject* p_obj = NULL;
  PyObject* p_obj1 = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_obj, &p_obj1)) {
    return NULL;
  }
  const WilsonMatrix& obj = py_convert_type<WilsonMatrix>(p_obj);
  const WilsonMatrix& obj1 = py_convert_type<WilsonMatrix>(p_obj1);
  return py_convert(matrix_trace(obj, obj1));
})

EXPORT(trace_sm_wm, {
  using namespace qlat;
  PyObject* p_obj = NULL;
  PyObject* p_obj1 = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_obj, &p_obj1)) {
    return NULL;
  }
  const SpinMatrix& obj = py_convert_type<SpinMatrix>(p_obj);
  const WilsonMatrix& obj1 = py_convert_type<WilsonMatrix>(p_obj1);
  return py_convert(matrix_trace(obj, obj1));
})

// -----------------------------------------

EXPORT(get_elem_wm_prop, {
  using namespace qlat;
  PyObject* p_obj = NULL;
  PyObject* p_xg = NULL;
  long m = 0;
  if (!PyArg_ParseTuple(args, "OO|l", &p_obj, &p_xg, &m)) {
    return NULL;
  }
  const Field<WilsonMatrix>& obj = py_convert_type_field<WilsonMatrix>(p_obj);
  const Coordinate xg = py_convert_data<Coordinate>(p_xg);
  const WilsonMatrix& wm = obj.get_elem(xg, m);
  const long cdata = (long)&wm;
  return py_convert(cdata);
})

EXPORT(get_elem_wm_sprop, {
  using namespace qlat;
  PyObject* p_obj = NULL;
  long idx = 0;
  long m = 0;
  if (!PyArg_ParseTuple(args, "Ol|l", &p_obj, &idx, &m)) {
    return NULL;
  }
  const SelectedField<WilsonMatrix>& obj =
      py_convert_type_sfield<WilsonMatrix>(p_obj);
  const WilsonMatrix& wm = obj.get_elem(idx, m);
  const long cdata = (long)&wm;
  return py_convert(cdata);
})

EXPORT(get_elem_wm_psprop, {
  using namespace qlat;
  PyObject* p_obj = NULL;
  long idx = 0;
  long m = 0;
  if (!PyArg_ParseTuple(args, "Ol|l", &p_obj, &idx, &m)) {
    return NULL;
  }
  const SelectedPoints<WilsonMatrix>& obj =
      py_convert_type_spoints<WilsonMatrix>(p_obj);
  const WilsonMatrix& wm = obj.get_elem(idx, m);
  const long cdata = (long)&wm;
  return py_convert(cdata);
})

// -----------------------------------------

EXPORT(set_state_wm, {
  using namespace qlat;
  PyObject* p_obj = NULL;
  PyObject* p_obj1 = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_obj, &p_obj1)) {
    return NULL;
  }
  WilsonMatrix& obj = py_convert_type<WilsonMatrix>(p_obj);
  const char* buffer = PyBytes_AsString(p_obj1);
  const long len = PyBytes_Size(p_obj1);
  pqassert(len == sizeof(WilsonMatrix));
  std::memcpy((void*)&obj, (void*)buffer, len);
  Py_RETURN_NONE;
})

EXPORT(get_state_wm, {
  using namespace qlat;
  PyObject* p_obj = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_obj)) {
    return NULL;
  }
  WilsonMatrix& obj = py_convert_type<WilsonMatrix>(p_obj);
  return PyBytes_FromStringAndSize((const char*)&obj, sizeof(WilsonMatrix));
})

EXPORT(set_state_sm, {
  using namespace qlat;
  PyObject* p_obj = NULL;
  PyObject* p_obj1 = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_obj, &p_obj1)) {
    return NULL;
  }
  SpinMatrix& obj = py_convert_type<SpinMatrix>(p_obj);
  const char* buffer = PyBytes_AsString(p_obj1);
  const long len = PyBytes_Size(p_obj1);
  pqassert(len == sizeof(SpinMatrix));
  std::memcpy((void*)&obj, (void*)buffer, len);
  Py_RETURN_NONE;
})

EXPORT(get_state_sm, {
  using namespace qlat;
  PyObject* p_obj = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_obj)) {
    return NULL;
  }
  SpinMatrix& obj = py_convert_type<SpinMatrix>(p_obj);
  return PyBytes_FromStringAndSize((const char*)&obj, sizeof(SpinMatrix));
})
