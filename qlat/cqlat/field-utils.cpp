#include "lib.h"

namespace qlat
{  //

template <class M>
PyObject* refresh_expanded_field_ctype(PyObject* p_field, PyObject* p_comm_plan)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  if (NULL == p_comm_plan) {
    refresh_expanded(f);
  } else {
    const CommPlan& cp = py_convert_type<CommPlan>(p_comm_plan);
    refresh_expanded(f, cp);
  }
  Py_RETURN_NONE;
}

template <class M>
PyObject* refresh_expanded_1_field_ctype(PyObject* p_field)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  refresh_expanded_1(f);
  Py_RETURN_NONE;
}

template <class M>
PyObject* merge_fields_ms_ctype(PyObject* p_field,
                                const std::vector<PyObject*>& p_f_vec,
                                const std::vector<int>& m_vec)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  const std::string ctype = py_get_ctype(p_field);
  const Int multiplicity = p_f_vec.size();
  std::vector<ConstHandle<Field<M>>> vec(multiplicity);
  for (Int m = 0; m < multiplicity; ++m) {
    qassert(py_get_ctype(p_f_vec[m]) == ctype);
    const Field<M>& fm = py_convert_type_field<M>(p_f_vec[m]);
    vec[m].init(fm);
  }
  merge_fields_ms(f, vec, m_vec);
  Py_RETURN_NONE;
}

}  // namespace qlat

EXPORT(make_field_expand_comm_plan, {  // tested: cqlat-fields
  using namespace qlat;
  PyObject* p_comm_plan = NULL;
  PyObject* p_comm_marks = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_comm_plan, &p_comm_marks)) {
    return NULL;
  }
  CommPlan& cp = py_convert_type<CommPlan>(p_comm_plan);
  const CommMarks& marks = py_convert_type<CommMarks>(p_comm_marks);
  cp = make_comm_plan(marks);
  Py_RETURN_NONE;
})

EXPORT(set_marks_field_all, {  // tested: cqlat-fields
  using namespace qlat;
  PyObject* p_comm_marks = NULL;
  PyObject* p_geo = NULL;
  Int multiplicity = 0;
  PyObject* p_tag = NULL;
  if (!PyArg_ParseTuple(args, "OOiO", &p_comm_marks, &p_geo, &multiplicity,
                        &p_tag)) {
    return NULL;
  }
  CommMarks& marks = py_convert_type<CommMarks>(p_comm_marks);
  const Geometry& geo = py_convert_type<Geometry>(p_geo);
  std::string tag = py_convert_data<std::string>(p_tag);
  set_marks_field_all(marks, geo, multiplicity, tag);
  Py_RETURN_NONE;
})

EXPORT(refresh_expanded_field, {  // tested: cqlat-fields
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_comm_plan = NULL;
  if (!PyArg_ParseTuple(args, "O|O", &p_field, &p_comm_plan)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, refresh_expanded_field_ctype, ctype, p_field,
                 p_comm_plan);
  return p_ret;
})

EXPORT(refresh_expanded_1_field, {  // tested: cqlat-fields
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, refresh_expanded_1_field_ctype, ctype, p_field);
  return p_ret;
})

EXPORT(merge_fields_ms_field, {  // tested: cqlat-fields
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_field_vec = NULL;
  PyObject* p_m_vec = NULL;
  if (!PyArg_ParseTuple(args, "OOO", &p_field, &p_field_vec, &p_m_vec)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  const std::vector<PyObject*> p_f_vec =
      py_convert_data<std::vector<PyObject*>>(p_field_vec);
  const std::vector<int> m_vec = py_convert_data<std::vector<int>>(p_m_vec);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, merge_fields_ms_ctype, ctype, p_field, p_f_vec, m_vec);
  return p_ret;
})
