#include "lib.h"

namespace qlat
{  //

template <class M>
PyObject* mk_spfield_ctype(const long n_points, const int multiplicity)
{
  SelectedPoints<M>* pspfield = new SelectedPoints<M>();
  if (n_points >= 0) {
    SelectedPoints<M>& spfield = *pspfield;
    spfield.init(n_points, multiplicity);
  }
  return py_convert((void*)pspfield);
}

template <class M>
PyObject* free_spfield_ctype(PyField& pf)
{
  SelectedPoints<M>& f = *(SelectedPoints<M>*)pf.cdata;
  delete &f;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_spfield_ctype(PyField& pf_new, PyField& pf)
{
  SelectedPoints<M>& f_new = *(SelectedPoints<M>*)pf_new.cdata;
  SelectedPoints<M>& f = *(SelectedPoints<M>*)pf.cdata;
  f_new = f;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_spfield_field_ctype(PyField& pspf, PyField& pf,
                                 const PointSelection& psel)
{
  SelectedPoints<M>& sp = *(SelectedPoints<M>*)pspf.cdata;
  const Field<M>& f = *(Field<M>*)pf.cdata;
  set_selected_points(sp, f, psel);
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_spfield_sfield_ctype(PyField& pspf, PyField& psf,
                                  const PointSelection& psel,
                                  const FieldSelection& fsel)
{
  SelectedPoints<M>& sp = *(SelectedPoints<M>*)pspf.cdata;
  const SelectedField<M>& sf = *(SelectedField<M>*)psf.cdata;
  set_selected_points(sp, sf, psel, fsel);
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_add_spfield_ctype(PyField& pf_new, PyField& pf)
{
  SelectedPoints<M>& f_new = *(SelectedPoints<M>*)pf_new.cdata;
  SelectedPoints<M>& f = *(SelectedPoints<M>*)pf.cdata;
  f_new += f;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_sub_spfield_ctype(PyField& pf_new, PyField& pf)
{
  SelectedPoints<M>& f_new = *(SelectedPoints<M>*)pf_new.cdata;
  SelectedPoints<M>& f = *(SelectedPoints<M>*)pf.cdata;
  f_new -= f;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_mul_spfield_ctype(PyField& pf, const Complex& factor)
{
  SelectedPoints<M>& f = *(SelectedPoints<M>*)pf.cdata;
  f *= factor;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_mul_spfield_ctype(PyField& pf, const double& factor)
{
  SelectedPoints<M>& f = *(SelectedPoints<M>*)pf.cdata;
  f *= factor;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_zero_spfield_ctype(PyField& pf)
{
  SelectedPoints<M>& f = *(SelectedPoints<M>*)pf.cdata;
  set_zero(f);
  Py_RETURN_NONE;
}

template <class M>
PyObject* get_n_points_spfield_ctype(PyField& pf)
{
  SelectedPoints<M>& spf = *(SelectedPoints<M>*)pf.cdata;
  const long ret = spf.n_points;
  return py_convert(ret);
}

template <class M>
PyObject* get_multiplicity_spfield_ctype(PyField& pf)
{
  SelectedPoints<M>& spf = *(SelectedPoints<M>*)pf.cdata;
  const long ret = spf.multiplicity;
  return py_convert(ret);
}

template <class M>
PyObject* qnorm_spfield_ctype(PyField& pf)
{
  SelectedPoints<M>& f = *(SelectedPoints<M>*)pf.cdata;
  const double ret = qnorm(f);
  return py_convert(ret);
}

}  // namespace qlat

EXPORT(mk_spfield, {
  using namespace qlat;
  PyObject* p_ctype = NULL;
  long n_points = -1;
  int multiplicity = 0;
  if (!PyArg_ParseTuple(args, "O|li", &p_ctype, &n_points, &multiplicity)) {
    return NULL;
  }
  std::string ctype;
  py_convert(ctype, p_ctype);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, mk_spfield_ctype, ctype, n_points, multiplicity);
  return p_ret;
});

EXPORT(free_spfield, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, free_spfield_ctype, pf.ctype, pf);
  return p_ret;
});

EXPORT(set_spfield, {
  using namespace qlat;
  PyObject* p_spfield_new = NULL;
  PyObject* p_spfield = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_spfield_new, &p_spfield)) {
    return NULL;
  }
  PyField pf_new = py_convert_field(p_spfield_new);
  PyField pf = py_convert_field(p_spfield);
  pqassert(pf_new.ctype == pf.ctype);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_spfield_ctype, pf.ctype, pf_new, pf);
  return p_ret;
});

EXPORT(set_spfield_field, {
  using namespace qlat;
  PyObject* p_spfield = NULL;
  PyObject* p_field = NULL;
  PyObject* p_psel = NULL;
  if (!PyArg_ParseTuple(args, "OOO", &p_spfield, &p_field, &p_psel)) {
    return NULL;
  }
  PyField pspf = py_convert_field(p_spfield);
  PyField pf = py_convert_field(p_field);
  const PointSelection& psel = py_convert_type<PointSelection>(p_psel);
  pqassert(pspf.ctype == pf.ctype);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_spfield_field_ctype, pf.ctype, pspf, pf, psel);
  return p_ret;
});

EXPORT(set_spfield_sfield, {
  using namespace qlat;
  PyObject* p_spfield = NULL;
  PyObject* p_sfield = NULL;
  PyObject* p_psel = NULL;
  PyObject* p_fsel = NULL;
  if (!PyArg_ParseTuple(args, "OOOO", &p_spfield, &p_sfield, &p_psel, &p_fsel)) {
    return NULL;
  }
  PyField pspf = py_convert_field(p_spfield);
  PyField psf = py_convert_field(p_sfield);
  const PointSelection& psel = py_convert_type<PointSelection>(p_psel);
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_fsel);
  pqassert(pspf.ctype == psf.ctype);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_spfield_sfield_ctype, psf.ctype, pspf, psf, psel, fsel);
  return p_ret;
});

EXPORT(set_add_spfield, {
  using namespace qlat;
  PyObject* p_field_new = NULL;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field_new, &p_field)) {
    return NULL;
  }
  PyField pf_new = py_convert_field(p_field_new);
  PyField pf = py_convert_field(p_field);
  pqassert(pf_new.ctype == pf.ctype);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_add_spfield_ctype, pf.ctype, pf_new, pf);
  return p_ret;
});

EXPORT(set_sub_spfield, {
  using namespace qlat;
  PyObject* p_field_new = NULL;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field_new, &p_field)) {
    return NULL;
  }
  PyField pf_new = py_convert_field(p_field_new);
  PyField pf = py_convert_field(p_field);
  pqassert(pf_new.ctype == pf.ctype);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_sub_spfield_ctype, pf.ctype, pf_new, pf);
  return p_ret;
});

EXPORT(set_mul_double_spfield, {
  using namespace qlat;
  PyObject* p_field = NULL;
  double factor = 0.0;
  if (!PyArg_ParseTuple(args, "Od", &p_field, &factor)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_mul_spfield_ctype, pf.ctype, pf, factor);
  return p_ret;
});

EXPORT(set_zero_spfield, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_zero_spfield_ctype, pf.ctype, pf);
  return p_ret;
});

EXPORT(get_n_points_spfield, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, get_n_points_spfield_ctype, pf.ctype, pf);
  return p_ret;
});

EXPORT(get_multiplicity_spfield, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, get_multiplicity_spfield_ctype, pf.ctype, pf);
  return p_ret;
});

EXPORT(qnorm_spfield, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, qnorm_spfield_ctype, pf.ctype, pf);
  return p_ret;
});
