#include "lib.h"

namespace qlat
{  //

template <class M>
PyObject* mk_sfield_ctype(PyObject* p_geo, const long n_elems, const int multiplicity)
{
  SelectedField<M>* pf = new SelectedField<M>();
  SelectedField<M>& f = *pf;
  if (p_geo != NULL) {
    const Geometry& geo = py_convert_type<Geometry>(p_geo);
    pqassert(n_elems > 0);
    pqassert(multiplicity > 0);
    f.init(geo, n_elems, multiplicity);
  }
  return py_convert((void*)pf);
}

template <class M>
PyObject* mk_sfield_fsel_ctype(PyObject* p_fsel, const int multiplicity)
{
  SelectedField<M>* pf = new SelectedField<M>();
  SelectedField<M>& f = *pf;
  pqassert(p_fsel != NULL);
  pqassert(multiplicity > 0);
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_fsel);
  f.init(fsel, multiplicity);
  return py_convert((void*)pf);
}

template <class M>
PyObject* free_sfield_ctype(PyField& pf)
{
  SelectedField<M>& f = *(SelectedField<M>*)pf.cdata;
  delete &f;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_sfield_ctype(PyField& pf_new, PyField& pf)
{
  SelectedField<M>& f_new = *(SelectedField<M>*)pf_new.cdata;
  SelectedField<M>& f = *(SelectedField<M>*)pf.cdata;
  f_new = f;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_sfield_field_ctype(PyField& psf, PyField& pf,
                                 const FieldSelection& fsel)
{
  SelectedField<M>& sf = *(SelectedField<M>*)psf.cdata;
  const Field<M>& f = *(Field<M>*)pf.cdata;
  set_selected_field(sf, f, fsel);
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_sfield_sfield_ctype(PyField& psf, PyField& psf0,
                                  const FieldSelection& fsel,
                                  const FieldSelection& fsel0)
{
  SelectedField<M>& sf = *(SelectedField<M>*)psf.cdata;
  const SelectedField<M>& sf0 = *(SelectedField<M>*)psf0.cdata;
  set_selected_field(sf, sf0, fsel, fsel0);
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_add_sfield_ctype(PyField& pf_new, PyField& pf)
{
  SelectedField<M>& f_new = *(SelectedField<M>*)pf_new.cdata;
  SelectedField<M>& f = *(SelectedField<M>*)pf.cdata;
  f_new += f;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_sub_sfield_ctype(PyField& pf_new, PyField& pf)
{
  SelectedField<M>& f_new = *(SelectedField<M>*)pf_new.cdata;
  SelectedField<M>& f = *(SelectedField<M>*)pf.cdata;
  f_new -= f;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_mul_sfield_ctype(PyField& pf, const Complex& factor)
{
  SelectedField<M>& f = *(SelectedField<M>*)pf.cdata;
  f *= factor;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_mul_sfield_ctype(PyField& pf, const double& factor)
{
  SelectedField<M>& f = *(SelectedField<M>*)pf.cdata;
  f *= factor;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_zero_sfield_ctype(PyField& pf)
{
  SelectedField<M>& f = *(SelectedField<M>*)pf.cdata;
  set_zero(f);
  Py_RETURN_NONE;
}

template <class M>
PyObject* get_n_elems_sfield_ctype(PyField& pf)
{
  SelectedField<M>& sf = *(SelectedField<M>*)pf.cdata;
  const long ret = sf.n_elems;
  return py_convert(ret);
}

template <class M>
PyObject* get_total_site_sfield_ctype(PyField& pf)
{
  SelectedField<M>& f = *(SelectedField<M>*)pf.cdata;
  const Coordinate ret = f.geo().total_site();
  return py_convert(ret);
}

template <class M>
PyObject* get_multiplicity_sfield_ctype(PyField& pf)
{
  SelectedField<M>& sf = *(SelectedField<M>*)pf.cdata;
  const long ret = sf.geo().multiplicity;
  return py_convert(ret);
}

template <class M>
PyObject* set_geo_sfield_ctype(Geometry& geo, PyField& pf)
{
  const SelectedField<M>& sf = *(SelectedField<M>*)pf.cdata;
  geo = sf.geo();
  Py_RETURN_NONE;
}

template <class M>
PyObject* qnorm_sfield_ctype(PyField& pf)
{
  SelectedField<M>& f = *(SelectedField<M>*)pf.cdata;
  const double ret = qnorm(f);
  return py_convert(ret);
}

}  // namespace qlat

EXPORT(mk_sfield, {
  using namespace qlat;
  PyObject* p_ctype = NULL;
  PyObject* p_geo = NULL;
  long n_elems = 0;
  int multiplicity = 0;
  if (!PyArg_ParseTuple(args, "O|Oli", &p_ctype, &p_geo, &n_elems, &multiplicity)) {
    return NULL;
  }
  std::string ctype;
  py_convert(ctype, p_ctype);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, mk_sfield_ctype, ctype, p_geo, n_elems, multiplicity);
  return p_ret;
});

EXPORT(mk_sfield_fsel, {
  using namespace qlat;
  PyObject* p_ctype = NULL;
  PyObject* p_fsel = NULL;
  int multiplicity = 0;
  if (!PyArg_ParseTuple(args, "OOi", &p_ctype, &p_fsel, &multiplicity)) {
    return NULL;
  }
  std::string ctype;
  py_convert(ctype, p_ctype);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, mk_sfield_fsel_ctype, ctype, p_fsel, multiplicity);
  return p_ret;
});

EXPORT(free_sfield, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, free_sfield_ctype, pf.ctype, pf);
  return p_ret;
});

EXPORT(set_sfield, {
  using namespace qlat;
  PyObject* p_sfield_new = NULL;
  PyObject* p_sfield = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_sfield_new, &p_sfield)) {
    return NULL;
  }
  PyField pf_new = py_convert_field(p_sfield_new);
  PyField pf = py_convert_field(p_sfield);
  pqassert(pf_new.ctype == pf.ctype);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_sfield_ctype, pf.ctype, pf_new, pf);
  return p_ret;
});

EXPORT(set_sfield_field, {
  using namespace qlat;
  PyObject* p_sfield = NULL;
  PyObject* p_field = NULL;
  PyObject* p_fsel = NULL;
  if (!PyArg_ParseTuple(args, "OOO", &p_sfield, &p_field, &p_fsel)) {
    return NULL;
  }
  PyField psf = py_convert_field(p_sfield);
  PyField pf = py_convert_field(p_field);
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_fsel);
  pqassert(psf.ctype == pf.ctype);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_sfield_field_ctype, pf.ctype, psf, pf, fsel);
  return p_ret;
});

EXPORT(set_sfield_sfield, {
  using namespace qlat;
  PyObject* p_sfield = NULL;
  PyObject* p_sfield0 = NULL;
  PyObject* p_fsel = NULL;
  PyObject* p_fsel0 = NULL;
  if (!PyArg_ParseTuple(args, "OOOO", &p_sfield, &p_sfield0, &p_fsel, &p_fsel0)) {
    return NULL;
  }
  PyField psf = py_convert_field(p_sfield);
  PyField psf0 = py_convert_field(p_sfield0);
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_fsel);
  const FieldSelection& fsel0 = py_convert_type<FieldSelection>(p_fsel0);
  pqassert(psf.ctype == psf0.ctype);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_sfield_sfield_ctype, psf0.ctype, psf, psf0, fsel, fsel0);
  return p_ret;
});

EXPORT(set_add_sfield, {
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
  FIELD_DISPATCH(p_ret, set_add_sfield_ctype, pf.ctype, pf_new, pf);
  return p_ret;
});

EXPORT(set_sub_sfield, {
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
  FIELD_DISPATCH(p_ret, set_sub_sfield_ctype, pf.ctype, pf_new, pf);
  return p_ret;
});

EXPORT(set_mul_double_sfield, {
  using namespace qlat;
  PyObject* p_field = NULL;
  double factor = 0.0;
  if (!PyArg_ParseTuple(args, "Od", &p_field, &factor)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_mul_sfield_ctype, pf.ctype, pf, factor);
  return p_ret;
});

EXPORT(set_zero_sfield, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_zero_sfield_ctype, pf.ctype, pf);
  return p_ret;
});

EXPORT(get_n_elems_sfield, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, get_n_elems_sfield_ctype, pf.ctype, pf);
  return p_ret;
});

EXPORT(get_total_site_sfield, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, get_total_site_sfield_ctype, pf.ctype, pf);
  return p_ret;
});

EXPORT(get_multiplicity_sfield, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, get_multiplicity_sfield_ctype, pf.ctype, pf);
  return p_ret;
});

EXPORT(set_geo_sfield, {
  using namespace qlat;
  PyObject* p_geo = NULL;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_geo, &p_field)) {
    return NULL;
  }
  Geometry& geo = py_convert_type<Geometry>(p_geo);
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_geo_sfield_ctype, pf.ctype, geo, pf);
  return p_ret;
});

EXPORT(qnorm_sfield, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, qnorm_sfield_ctype, pf.ctype, pf);
  return p_ret;
});
