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

template <class M, class N>
PyObject* assign_as_field_ctype(Field<N>& f, PyObject* p_field1)
{
  const Field<M>& f1 = py_convert_type_field<M>(p_field1);
  assign(f, f1);
  Py_RETURN_NONE;
}

template <class M, class N>
PyObject* assign_from_field_ctype(PyObject* p_field, const Field<N>& f1)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  assign(f, f1);
  Py_RETURN_NONE;
}

template <class M>
PyObject* get_elems_field_ctype(PyObject* p_field, const Coordinate& xg)
{
  const Field<M>& f = py_convert_type_field<M>(p_field);
  return py_convert(field_get_elems(f, xg));
}

template <class M>
PyObject* get_elem_field_ctype(PyObject* p_field, const Coordinate& xg,
                               const Int m)
{
  const Field<M>& f = py_convert_type_field<M>(p_field);
  if (m >= 0) {
    return py_convert(field_get_elem(f, xg, m));
  } else {
    return py_convert(field_get_elem(f, xg));
  }
}

template <class M>
PyObject* set_elems_field_ctype(PyObject* p_field, const Coordinate& xg,
                                PyObject* p_val)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  const Int multiplicity = f.multiplicity;
  qassert((Long)PyBytes_Size(p_val) == (Long)multiplicity * (Long)sizeof(M));
  const Vector<M> val((M*)PyBytes_AsString(p_val), multiplicity);
  field_set_elems(f, xg, val);
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_elem_field_ctype(PyObject* p_field, const Coordinate& xg,
                               const Int m, PyObject* p_val)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  qassert(PyBytes_Size(p_val) == sizeof(M));
  const M& val = *(M*)PyBytes_AsString(p_val);
  if (m >= 0) {
    field_set_elem(f, xg, m, val);
  } else {
    field_set_elem(f, xg, val);
  }
  Py_RETURN_NONE;
}

template <class M>
PyObject* get_elems_field_ctype(PyObject* p_field, const Long index)
{
  const Field<M>& f = py_convert_type_field<M>(p_field);
  return py_convert(f.get_elems_const(index));
}

template <class M>
PyObject* get_elem_field_ctype(PyObject* p_field, const Long index, const Int m)
{
  const Field<M>& f = py_convert_type_field<M>(p_field);
  if (m >= 0) {
    return py_convert(f.get_elem(index, m));
  } else {
    return py_convert(f.get_elem(index));
  }
}

template <class M>
PyObject* set_elems_field_ctype(PyObject* p_field, const Long index,
                                PyObject* p_val)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  const Int multiplicity = f.multiplicity;
  qassert((Long)PyBytes_Size(p_val) == (Long)multiplicity * (Long)sizeof(M));
  const Vector<M> val((M*)PyBytes_AsString(p_val), multiplicity);
  assign(f.get_elems(index), val);
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_elem_field_ctype(PyObject* p_field, const Long index, const Int m,
                               PyObject* p_val)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  qassert(PyBytes_Size(p_val) == sizeof(M));
  const M& val = *(M*)PyBytes_AsString(p_val);
  if (m >= 0) {
    f.get_elem(index, m) = val;
  } else {
    f.get_elem(index) = val;
  }
  Py_RETURN_NONE;
}

template <class M>
PyObject* fft_fields_ctype(const std::vector<PyObject*>& p_field_vec,
                           const std::vector<int>& fft_dirs,
                           const std::vector<bool>& fft_is_forwards,
                           Int mode_fft = 1)
{
  const Long n_field = p_field_vec.size();
  std::vector<Handle<Field<M>>> vec(n_field);
  for (Long i = 0; i < n_field; ++i) {
    vec[i].init(py_convert_type_field<M>(p_field_vec[i]));
  }
  fft_complex_fields(vec, fft_dirs, fft_is_forwards, mode_fft);
  Py_RETURN_NONE;
}

template <class M>
PyObject* field_shift_field_ctype(PyObject* p_field_new, PyObject* p_field,
                                  const Coordinate& shift)
{
  Field<M>& f_new = py_convert_type_field<M>(p_field_new);
  const Field<M>& f = py_convert_type_field<M>(p_field);
  field_shift(f_new, f, shift);
  Py_RETURN_NONE;
}

template <class M>
PyObject* reflect_field_ctype(PyObject* p_field)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  reflect_field(f);
  Py_RETURN_NONE;
}

template <class M>
PyObject* split_fields_field_ctype(const std::vector<PyObject*>& p_f_vec,
                                   PyObject* p_field)
{
  const std::string ctype = py_get_ctype(p_field);
  const Field<M>& f = py_convert_type_field<M>(p_field);
  const Int nf = p_f_vec.size();
  std::vector<Handle<Field<M>>> vec(nf);
  for (Int i = 0; i < nf; ++i) {
    qassert(py_get_ctype(p_f_vec[i]) == ctype);
    Field<M>& fi = py_convert_type_field<M>(p_f_vec[i]);
    vec[i].init(fi);
  }
  split_fields(vec, f);
  Py_RETURN_NONE;
}

template <class M>
PyObject* merge_fields_field_ctype(PyObject* p_field,
                                   const std::vector<PyObject*>& p_f_vec)
{
  const std::string ctype = py_get_ctype(p_field);
  Field<M>& f = py_convert_type_field<M>(p_field);
  const Int nf = p_f_vec.size();
  std::vector<ConstHandle<Field<M>>> vec(nf);
  for (Int i = 0; i < nf; ++i) {
    qassert(py_get_ctype(p_f_vec[i]) == ctype);
    const Field<M>& fi = py_convert_type_field<M>(p_f_vec[i]);
    vec[i].init(fi);
  }
  merge_fields(f, vec);
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

template <class M>
PyObject* qnorm_field_field_ctype(FieldM<RealD, 1>& f, PyObject* p_field1)
{
  const Field<M>& f1 = py_convert_type_field<M>(p_field1);
  qnorm_field(f, f1);
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

EXPORT(reflect_field, {  // tested: cqlat-fields
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, reflect_field_ctype, ctype, p_field);
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

EXPORT(set_sqrt_field, {  // tested: cqlat-fields
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_field1 = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field, &p_field1)) {
    return NULL;
  }
  Field<RealD>& f = py_convert_type_field<RealD>(p_field);
  const Field<RealD>& f1 = py_convert_type_field<RealD>(p_field1);
  const Geometry geo = geo_resize(f1.geo());
  qassert(geo.is_only_local);
  f.init();
  f.init(geo, f1.multiplicity);
  qacc_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<RealD> f1v = f1.get_elems_const(xl);
    Vector<RealD> fv = f.get_elems(xl);
    for (Int m = 0; m < f1.multiplicity; ++m) {
      fv[m] = std::sqrt(f1v[m]);
    }
  });
  Py_RETURN_NONE;
})
