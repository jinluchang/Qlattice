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
PyObject* assign_as_complex_field_ctype(Field<Complex>& f, PyObject* p_field1)
{
  const Field<M>& f1 = py_convert_type_field<M>(p_field1);
  assign(f, f1);
  Py_RETURN_NONE;
}

template <class M>
PyObject* assign_from_complex_field_ctype(PyObject* p_field, const Field<Complex>& f1)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  assign(f, f1);
  Py_RETURN_NONE;
}

template <class M>
PyObject* get_elems_field_ctype(PyField& pf, const Coordinate& xg)
{
  const Field<M>& f = *(Field<M>*)pf.cdata;
  return py_convert(field_get_elems(f, xg));
}

template <class M>
PyObject* get_elem_field_ctype(PyField& pf, const Coordinate& xg, const int m)
{
  const Field<M>& f = *(Field<M>*)pf.cdata;
  if (m >= 0) {
    return py_convert(field_get_elem(f, xg, m));
  } else {
    return py_convert(field_get_elem(f, xg));
  }
}

template <class M>
PyObject* glb_sum_double_field_ctype(PyField& pf)
{
  const Field<M>& f = *(Field<M>*)pf.cdata;
  return py_convert(field_glb_sum_double(f));
}

template <class M>
PyObject* glb_sum_long_field_ctype(PyField& pf)
{
  const Field<M>& f = *(Field<M>*)pf.cdata;
  return py_convert(field_glb_sum_long(f));
}

template <class M>
PyObject* glb_sum_tslice_double_field_ctype(PyField& pspf, PyField& pf)
{
  SelectedPoints<M>& sp = *(SelectedPoints<M>*)pspf.cdata;
  const Field<M>& f = *(Field<M>*)pf.cdata;
  field_glb_sum_tslice_double(sp, f);
  Py_RETURN_NONE;
}

template <class M>
PyObject* glb_sum_tslice_long_field_ctype(PyField& pspf, PyField& pf)
{
  SelectedPoints<M>& sp = *(SelectedPoints<M>*)pspf.cdata;
  const Field<M>& f = *(Field<M>*)pf.cdata;
  field_glb_sum_tslice_long(sp, f);
  Py_RETURN_NONE;
}

template <class M>
PyObject* fft_dir_complex_field_ctype(PyField& pf1, PyField& pf, const int fft_dir, const bool is_forward)
{
  Field<M>& f1 = *(Field<M>*)pf1.cdata;
  const Field<M>& f = *(Field<M>*)pf.cdata;
  fft_complex_field_dir(f1, f, fft_dir, is_forward);
  Py_RETURN_NONE;
}

template <class M>
PyObject* split_fields_field_ctype(std::vector<PyField>& pf_vec, PyField& pf)
{
  Field<M>& f = *(Field<M>*)pf.cdata;
  const int nf = pf_vec.size();
  std::vector<Handle<Field<M> > > vec(nf);
  for (int i = 0; i < nf; ++i) {
    pqassert(pf_vec[i].ctype == pf.ctype);
    vec[i].init(*(Field<M>*)pf_vec[i].cdata);
  }
  split_fields(vec, f);
  Py_RETURN_NONE;
}

template <class M>
PyObject* merge_fields_field_ctype(PyField& pf,
                                   const std::vector<PyField>& pf_vec)
{
  Field<M>& f = *(Field<M>*)pf.cdata;
  const int nf = pf_vec.size();
  std::vector<ConstHandle<Field<M> > > vec(nf);
  for (int i = 0; i < nf; ++i) {
    pqassert(pf_vec[i].ctype == pf.ctype);
    vec[i].init(*(Field<M>*)pf_vec[i].cdata);
  }
  merge_fields(f, vec);
  Py_RETURN_NONE;
}

template <class M>
PyObject* merge_fields_ms_ctype(PyObject* p_field,
                                      const std::vector<PyObject*>& p_f_vec,
                                      const std::vector<int> m_vec)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  const int multiplicity = p_f_vec.size();
  std::vector<ConstHandle<Field<M> > > vec(multiplicity);
  for (int m = 0; m < multiplicity; ++m) {
    const std::string ctype = py_get_ctype(p_f_vec[m]);
    pqassert(check_ctype_name<M>(ctype));
    vec[m].init(py_convert_type_field<M>(p_f_vec[m]));
  }
  merge_fields_ms(f, vec, m_vec);
  Py_RETURN_NONE;
}

}  // namespace qlat

EXPORT(mk_field_expand_comm_plan, {
  using namespace qlat;
  CommPlan* pcp = new CommPlan();
  return py_convert((void*)pcp);
});

EXPORT(free_field_expand_comm_plan, {
  using namespace qlat;
  return free_obj<CommPlan>(args);
});

EXPORT(set_field_expand_comm_plan, {
  using namespace qlat;
  return set_obj<CommPlan>(args);
});

EXPORT(make_field_expand_comm_plan, {
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
});

EXPORT(set_marks_field_all, {
  using namespace qlat;
  PyObject* p_comm_marks = NULL;
  PyObject* p_geo = NULL;
  PyObject* p_tag = NULL;
  if (!PyArg_ParseTuple(args, "OOO", &p_comm_marks, &p_geo, &p_tag)) {
    return NULL;
  }
  CommMarks& marks = py_convert_type<CommMarks>(p_comm_marks);
  const Geometry& geo = py_convert_type<Geometry>(p_geo);
  std::string tag = py_convert_data<std::string>(p_tag);
  set_marks_field_all(marks, geo, tag);
  Py_RETURN_NONE;
});

EXPORT(refresh_expanded_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_comm_plan = NULL;
  if (!PyArg_ParseTuple(args, "O|O", &p_field, &p_comm_plan)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, refresh_expanded_field_ctype, ctype, p_field, p_comm_plan);
  return p_ret;
});

EXPORT(refresh_expanded_1_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, refresh_expanded_1_field_ctype, ctype, p_field);
  return p_ret;
});

EXPORT(set_phase_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_lmom = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field, &p_lmom)) {
    return NULL;
  }
  const CoordinateD lmom = py_convert_data<CoordinateD>(p_lmom);
  FieldM<Complex, 1>& f = py_convert_type_field<Complex, 1>(p_field);
  set_phase_field(f, lmom);
  Py_RETURN_NONE;
});

EXPORT(assign_as_complex_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_field1 = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field, &p_field1)) {
    return NULL;
  }
  Field<Complex>& f = py_convert_type_field<Complex>(p_field);
  const std::string ctype = py_get_ctype(p_field1);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, assign_as_complex_field_ctype, ctype, f, p_field1);
  return p_ret;
});

EXPORT(assign_from_complex_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_field1 = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field, &p_field1)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  const Field<Complex>& f1 = py_convert_type_field<Complex>(p_field1);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, assign_from_complex_field_ctype, ctype, p_field, f1);
  return p_ret;
});

EXPORT(get_elems_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_xg = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field, &p_xg)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  Coordinate xg;
  py_convert(xg, p_xg);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, get_elems_field_ctype, pf.ctype, pf, xg);
  return p_ret;
});

EXPORT(get_elem_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_xg = NULL;
  long m = -1;
  if (!PyArg_ParseTuple(args, "OO|l", &p_field, &p_xg, &m)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  Coordinate xg;
  py_convert(xg, p_xg);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, get_elem_field_ctype, pf.ctype, pf, xg, m);
  return p_ret;
});

EXPORT(glb_sum_double_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, glb_sum_double_field_ctype, pf.ctype, pf);
  return p_ret;
});

EXPORT(glb_sum_long_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, glb_sum_long_field_ctype, pf.ctype, pf);
  return p_ret;
});

EXPORT(glb_sum_tslice_double_field, {
  using namespace qlat;
  PyObject* p_spfield = NULL;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_spfield, &p_field)) {
    return NULL;
  }
  PyField pspf = py_convert_field(p_spfield);
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, glb_sum_tslice_double_field_ctype, pf.ctype, pspf, pf);
  return p_ret;
});

EXPORT(glb_sum_tslice_long_field, {
  using namespace qlat;
  PyObject* p_spfield = NULL;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_spfield, &p_field)) {
    return NULL;
  }
  PyField pspf = py_convert_field(p_spfield);
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, glb_sum_tslice_long_field_ctype, pf.ctype, pspf, pf);
  return p_ret;
});

EXPORT(fft_dir_complex_field, {
  // forward compute
  // field1(k) <- \sum_{x} exp( - ii * 2 pi * k * x ) field(x)
  // backwards compute
  // field1(x) <- \sum_{k} exp( + ii * 2 pi * k * x ) field(k)
  using namespace qlat;
  PyObject* p_field1 = NULL;
  PyObject* p_field = NULL;
  int fft_dir = 0;
  bool is_forward = true;
  if (!PyArg_ParseTuple(args, "OOib", &p_field1, &p_field, &fft_dir, &is_forward)) {
    return NULL;
  }
  PyField pf1 = py_convert_field(p_field1);
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, fft_dir_complex_field_ctype, pf.ctype, pf1, pf, fft_dir, is_forward);
  return p_ret;
});

EXPORT(split_fields_field, {
  using namespace qlat;
  PyObject* p_field_vec = NULL;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field_vec, &p_field)) {
    return NULL;
  }
  std::vector<PyField> pf_vec;
  py_convert(pf_vec, p_field_vec);
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, split_fields_field_ctype, pf.ctype, pf_vec, pf);
  return p_ret;
});

EXPORT(merge_fields_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_field_vec = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field, &p_field_vec)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  std::vector<PyField> pf_vec;
  py_convert(pf_vec, p_field_vec);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, merge_fields_field_ctype, pf.ctype, pf, pf_vec);
  return p_ret;
}
);

EXPORT(merge_fields_ms_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_field_vec = NULL;
  PyObject* p_m_vec = NULL;
  if (!PyArg_ParseTuple(args, "OOO", &p_field, &p_field_vec, &p_m_vec)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  const std::vector<PyObject*> p_f_vec =
      py_convert_data<std::vector<PyObject*> >(p_field_vec);
  const std::vector<int> m_vec = py_convert_data<std::vector<int> >(p_m_vec);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, merge_fields_ms_ctype, ctype, p_field, p_f_vec, m_vec);
  return p_ret;
});
