#include "lib.h"
#include <qlat/vector_utils/utils_FFT_GPU.h>

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
PyObject* set_elems_field_ctype(PyObject* p_field, const Coordinate& xg,
                                PyObject* p_val)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  const int multiplicity = f.geo().multiplicity;
  pqassert((long)PyBytes_Size(p_val) == (long)multiplicity * (long)sizeof(M));
  const Vector<M> val((M*)PyBytes_AsString(p_val), multiplicity);
  field_set_elems(f, xg, val);
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_elem_field_ctype(PyObject* p_field, const Coordinate& xg,
                               const int m, PyObject* p_val)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  pqassert(PyBytes_Size(p_val) == sizeof(M));
  const M& val = *(M*)PyBytes_AsString(p_val);
  field_set_elem(f, xg, m, val);
  Py_RETURN_NONE;
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
PyObject* glb_sum_tslice_double_field_ctype(PyObject* p_spfield,
                                            PyObject* p_field, const int t_dir)
{
  SelectedPoints<M>& sp = py_convert_type_spoints<M>(p_spfield);
  const Field<M>& f = py_convert_type_field<M>(p_field);
  field_glb_sum_tslice_double(sp, f, t_dir);
  Py_RETURN_NONE;
}

template <class M>
PyObject* glb_sum_tslice_long_field_ctype(PyObject* p_spfield,
                                          PyObject* p_field, const int t_dir)
{
  SelectedPoints<M>& sp = py_convert_type_spoints<M>(p_spfield);
  const Field<M>& f = py_convert_type_field<M>(p_field);
  field_glb_sum_tslice_long(sp, f, t_dir);
  Py_RETURN_NONE;
}

template <class M>
PyObject* fft_fields_ctype(const std::vector<PyObject*> p_field_vec,
                           const std::vector<int> fft_dirs,
                           const std::vector<bool> fft_is_forwards, int mode_fft = 1)
{
  const long n_field = p_field_vec.size();
  std::vector<Handle<Field<M> > > vec(n_field);
  for (long i = 0; i < n_field; ++i) {
    vec[i].init(py_convert_type_field<M>(p_field_vec[i]));
  }
  int use_plan = 0;
  bool fft_direction = false;
  bool ft4D = false;
  if (mode_fft == 1) {
    ////check 3D
    if (fft_dirs.size() == 3 and fft_dirs[0] == 0 and fft_dirs[1] == 1 and
        fft_dirs[2] == 2) {
      if (fft_is_forwards[0] == fft_is_forwards[1] and
          fft_is_forwards[0] == fft_is_forwards[2]) {
        use_plan = 1;
        fft_direction = fft_is_forwards[0];
        ft4D = false;
      }
    }
    ////check 4D
    if (fft_dirs.size() == 4 and fft_dirs[0] == 0 and fft_dirs[1] == 1 and
        fft_dirs[2] == 2 and fft_dirs[3] == 3) {
      if (fft_is_forwards[0] == fft_is_forwards[1] and
          fft_is_forwards[0] == fft_is_forwards[2])
        if (fft_is_forwards[0] == fft_is_forwards[3]) {
          use_plan = 1;
          fft_direction = fft_is_forwards[0];
          ft4D = true;
        }
    }
  }
  if (use_plan == 0) {
    for (long i = 0; i < n_field; ++i) {
      Field<M> ft;
      for (long k = 0; k < (long)fft_dirs.size(); ++k) {
        ft = vec[i]();
        const int fft_dir = fft_dirs[k];
        const bool is_forward = fft_is_forwards[k];
        fft_complex_field_dir(vec[i](), ft, fft_dir, is_forward);
      }
    }
  }
  if (use_plan == 1) {
    fft_fieldM(vec, fft_direction, ft4D);
  }
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
  const std::string ctype = py_get_ctype(p_field);
  const int multiplicity = p_f_vec.size();
  std::vector<ConstHandle<Field<M> > > vec(multiplicity);
  for (int m = 0; m < multiplicity; ++m) {
    pqassert(ctype == py_get_ctype(p_f_vec[m]));
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
})

EXPORT(free_field_expand_comm_plan, {
  using namespace qlat;
  return free_obj<CommPlan>(args);
})

EXPORT(set_field_expand_comm_plan, {
  using namespace qlat;
  return set_obj<CommPlan>(args);
})

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
})

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
})

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
})

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
})

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
})

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
})

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
})

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
})

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
})

EXPORT(set_elems_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_xg = NULL;
  PyObject* p_val = NULL;
  if (!PyArg_ParseTuple(args, "OOO", &p_field, &p_xg, &p_val)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  const Coordinate xg = py_convert_data<Coordinate>(p_xg);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_elems_field_ctype, ctype, p_field, xg, p_val);
  return p_ret;
})

EXPORT(set_elem_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_xg = NULL;
  long m = 0;
  PyObject* p_val = NULL;
  if (!PyArg_ParseTuple(args, "OOlO", &p_field, &p_xg, &m, &p_val)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  const Coordinate xg = py_convert_data<Coordinate>(p_xg);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_elem_field_ctype, ctype, p_field, xg, m, p_val);
  return p_ret;
})

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
})

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
})

EXPORT(glb_sum_tslice_double_field, {
  using namespace qlat;
  PyObject* p_spfield = NULL;
  PyObject* p_field = NULL;
  int t_dir = 3;
  if (!PyArg_ParseTuple(args, "OO|i", &p_spfield, &p_field, &t_dir)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, glb_sum_tslice_double_field_ctype, ctype, p_spfield,
                 p_field, t_dir);
  return p_ret;
})

EXPORT(glb_sum_tslice_long_field, {
  using namespace qlat;
  PyObject* p_spfield = NULL;
  PyObject* p_field = NULL;
  int t_dir = 3;
  if (!PyArg_ParseTuple(args, "OO|i", &p_spfield, &p_field, &t_dir)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, glb_sum_tslice_long_field_ctype, ctype, p_spfield,
                 p_field, t_dir);
  return p_ret;
})

EXPORT(fft_fields, {
  // forward compute
  // field(k) <- \sum_{x} exp( - ii * 2 pi * k * x ) field(x)
  // backwards compute
  // field(x) <- \sum_{k} exp( + ii * 2 pi * k * x ) field(k)
  using namespace qlat;
  PyObject* p_fields = NULL;
  PyObject* p_fft_dirs = NULL;
  PyObject* p_fft_is_forwards = NULL;
  int mode_fft = 1;
  if (!PyArg_ParseTuple(args, "OOO|i", &p_fields, &p_fft_dirs,
                        &p_fft_is_forwards, &mode_fft)) {
    return NULL;
  }
  const std::vector<PyObject*> p_field_vec =
      py_convert_data<std::vector<PyObject*> >(p_fields);
  pqassert(p_field_vec.size() >= 1);
  const std::string ctype = py_get_ctype(p_field_vec[0]);
  for (long i = 0; i < (long)p_field_vec.size(); ++i) {
    pqassert(ctype == py_get_ctype(p_field_vec[i]));
  }
  const std::vector<int> fft_dirs =
      py_convert_data<std::vector<int> >(p_fft_dirs);
  const std::vector<bool> fft_is_forwards =
      py_convert_data<std::vector<bool> >(p_fft_is_forwards);
  pqassert(fft_dirs.size() == fft_is_forwards.size());
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, fft_fields_ctype, ctype, p_field_vec, fft_dirs,
                 fft_is_forwards, mode_fft);
  return p_ret;
})

EXPORT(field_shift_field, {
  using namespace qlat;
  PyObject* p_field_new = NULL;
  PyObject* p_field = NULL;
  PyObject* p_shift = NULL;
  if (!PyArg_ParseTuple(args, "OOO", &p_field_new, &p_field, &p_shift)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  pqassert(py_get_ctype(p_field_new) == py_get_ctype(p_field));
  const Coordinate shift = py_convert_data<Coordinate>(p_shift);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, field_shift_field_ctype, ctype, p_field_new, p_field,
                 shift);
  return p_ret;
})

EXPORT(reflect_field, {
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
})

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
})

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
})
