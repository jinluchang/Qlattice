#include "lib.h"

namespace qlat
{  //

template <class M>
PyObject* set_add_sfield_ctype(PyObject* pf_new, PyObject* pf)
{
  SelectedField<M>& f_new = py_convert_type_sfield<M>(pf_new);
  SelectedField<M>& f = py_convert_type_sfield<M>(pf);
  f_new += f;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_mul_sfield_ctype(PyObject* pf, const RealD& factor)
{
  SelectedField<M>& f = py_convert_type_sfield<M>(pf);
  f *= factor;
  Py_RETURN_NONE;
}

template <class M>
PyObject* acc_field_sfield_ctype(PyObject* p_field, PyObject* p_sfield,
                                 const FieldSelection& fsel)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  const SelectedField<M>& sf = py_convert_type_sfield<M>(p_sfield);
  acc_field(f, sf, fsel);
  Py_RETURN_NONE;
}

template <class M>
PyObject* acc_field_spfield_ctype(PyObject* p_field, PyObject* p_spfield,
                                  const Geometry& geo,
                                  const PointsSelection& psel)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  const SelectedPoints<M>& sp = py_convert_type_spoints<M>(p_spfield);
  acc_field(f, sp, geo, psel);
  Py_RETURN_NONE;
}

template <class M>
PyObject* glb_sum_tslice_long_sfield_ctype(PyObject* p_spfield,
                                           PyObject* p_field,
                                           const FieldSelection& fsel,
                                           const Int t_dir)
{
  SelectedPoints<M>& sp = py_convert_type_spoints<M>(p_spfield);
  const SelectedField<M>& f = py_convert_type_sfield<M>(p_field);
  field_glb_sum_tslice(sp, f, fsel, t_dir);
  Py_RETURN_NONE;
}

}  // namespace qlat

EXPORT(set_add_sfield, {  // tested: cqlat-selected-field
  using namespace qlat;
  PyObject* p_field_new = NULL;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field_new, &p_field)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  qassert(py_get_ctype(p_field_new) == ctype);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_add_sfield_ctype, ctype, p_field_new, p_field);
  return p_ret;
})

EXPORT(set_mul_double_sfield, {  // tested: cqlat-selected-field
  using namespace qlat;
  PyObject* p_field = NULL;
  RealD factor = 0.0;
  if (!PyArg_ParseTuple(args, "Od", &p_field, &factor)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_mul_sfield_ctype, ctype, p_field, factor);
  return p_ret;
})

EXPORT(acc_field_sfield, {  // tested: cqlat-selected-field
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_sfield = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field, &p_sfield)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_sfield);
  qassert(py_get_ctype(p_sfield) == ctype);
  QLAT_PUSH_DIAGNOSTIC_DISABLE_DANGLING_REF;
  const FieldSelection& fsel =
      py_convert_type<FieldSelection>(p_sfield, "fsel");
  QLAT_DIAGNOSTIC_POP;
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, acc_field_sfield_ctype, ctype, p_field, p_sfield, fsel);
  return p_ret;
})

EXPORT(acc_field_spfield, {  // tested: cqlat-selected-field
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_spfield = NULL;
  PyObject* p_geo = NULL;
  PyObject* p_psel = NULL;
  if (!PyArg_ParseTuple(args, "OOOO", &p_field, &p_spfield, &p_geo, &p_psel)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_spfield);
  qassert(py_get_ctype(p_field) == ctype);
  const Geometry& geo = py_convert_type<Geometry>(p_geo);
  const PointsSelection& psel = py_convert_type<PointsSelection>(p_psel);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, acc_field_spfield_ctype, ctype, p_field, p_spfield,
                 geo, psel);
  return p_ret;
})

EXPORT(glb_sum_tslice_long_sfield, {  // tested: cqlat-selected-field
  using namespace qlat;
  PyObject* p_spfield = NULL;
  PyObject* p_field = NULL;
  Int t_dir = 3;
  if (!PyArg_ParseTuple(args, "OO|i", &p_spfield, &p_field, &t_dir)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  qassert(py_get_ctype(p_spfield) == ctype);
  QLAT_PUSH_DIAGNOSTIC_DISABLE_DANGLING_REF;
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_field, "fsel");
  QLAT_DIAGNOSTIC_POP;
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, glb_sum_tslice_long_sfield_ctype, ctype, p_spfield,
                 p_field, fsel, t_dir);
  return p_ret;
})
