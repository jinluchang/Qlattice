#include "lib.h"

namespace qlat
{  //

template <class M>
PyObject* mk_sfield_ctype(int dummy)
{
  (void)dummy;
  SelectedField<M>* pf = new SelectedField<M>();
  return py_convert((void*)pf);
}

template <class M>
PyObject* mk_sfield_fsel_ctype(const FieldSelection& fsel, const int multiplicity)
{
  SelectedField<M>* pf = new SelectedField<M>();
  SelectedField<M>& f = *pf;
  pqassert(multiplicity > 0);
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
PyObject* set_field_sfield_ctype(PyField& pf, PyField& psf,
                                 const FieldSelection& fsel)
{
  Field<M>& f = *(Field<M>*)pf.cdata;
  const SelectedField<M>& sf = *(SelectedField<M>*)psf.cdata;
  set_field_selected(f, sf, fsel);
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
PyObject* acc_field_sfield_ctype(PyObject* p_field, PyObject* p_sfield,
                                 const FieldSelection& fsel)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  const SelectedField<M>& sf = py_convert_type_sfield<M>(p_sfield);
  acc_field(f, sf, fsel);
  Py_RETURN_NONE;
}

template <class M>
PyObject* field_shift_sfield_ctype(PyObject* p_sfield_new, PyObject* p_sfield,
                                   const Coordinate& shift, const bool is_reflect)
{
  FieldSelection& fsel_new = py_convert_type<FieldSelection>(p_sfield_new, "fsel");
  SelectedField<M>& sf_new = py_convert_type_sfield<M>(p_sfield_new);
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_sfield, "fsel");
  const SelectedField<M>& sf = py_convert_type_sfield<M>(p_sfield);
  field_shift(sf_new, fsel_new, sf, fsel, shift, is_reflect);
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

template <class M>
PyObject* get_elems_sfield_ctype(PyObject* p_sfield, const long idx)
{
  const SelectedField<M>& f = py_convert_type_sfield<M>(p_sfield);
  return py_convert(f.get_elems_const(idx));
}

template <class M>
PyObject* get_elem_sfield_ctype(PyObject* p_sfield, const long idx, const int m)
{
  const SelectedField<M>& f = py_convert_type_sfield<M>(p_sfield);
  if (m >= 0) {
    return py_convert(f.get_elem(idx, m));
  } else {
    return py_convert(f.get_elem(idx));
  }
}

template <class M>
PyObject* set_elems_sfield_ctype(PyObject* p_field, const long idx,
                                 PyObject* p_val)
{
  SelectedField<M>& f = py_convert_type_sfield<M>(p_field);
  const int multiplicity = f.geo().multiplicity;
  pqassert((long)PyBytes_Size(p_val) == (long)multiplicity * (long)sizeof(M));
  const Vector<M> val((M*)PyBytes_AsString(p_val), multiplicity);
  assign(f.get_elems(idx), val);
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_elem_sfield_ctype(PyObject* p_field, const long idx,
                                const int m, PyObject* p_val)
{
  SelectedField<M>& f = py_convert_type_sfield<M>(p_field);
  pqassert(PyBytes_Size(p_val) == sizeof(M));
  const M& val = *(M*)PyBytes_AsString(p_val);
  f.get_elem(idx, m) = val;
  Py_RETURN_NONE;
}

template <class M>
PyObject* glb_sum_tslice_double_sfield_ctype(PyObject* p_spfield,
                                             PyObject* p_field,
                                             const FieldSelection& fsel,
                                             const int t_dir)
{
  SelectedPoints<M>& sp = py_convert_type_spoints<M>(p_spfield);
  const SelectedField<M>& f = py_convert_type_sfield<M>(p_field);
  field_glb_sum_tslice_double(sp, f, fsel, t_dir);
  Py_RETURN_NONE;
}

template <class M>
PyObject* glb_sum_tslice_long_sfield_ctype(PyObject* p_spfield,
                                           PyObject* p_field,
                                           const FieldSelection& fsel,
                                           const int t_dir)
{
  SelectedPoints<M>& sp = py_convert_type_spoints<M>(p_spfield);
  const SelectedField<M>& f = py_convert_type_sfield<M>(p_field);
  field_glb_sum_tslice_long(sp, f, fsel, t_dir);
  Py_RETURN_NONE;
}

template <class M>
PyObject* save_sfield_ctype(PyField& pf, const std::string& path,
                            const FieldSelection& fsel)
{
  const SelectedField<M>& f = *(SelectedField<M>*)pf.cdata;
  const long ret = write_selected_field(f, path, fsel);
  return py_convert(ret);
}

template <class M>
PyObject* load_sfield_ctype(PyField& pf, const std::string& path,
                            const FieldSelection& fsel)
{
  SelectedField<M>& f = *(SelectedField<M>*)pf.cdata;
  const long ret = read_selected_field(f, path, fsel);
  return py_convert(ret);
}

template <class M>
PyObject* convert_float_from_double_sfield_ctype(PyField& pf_new, PyField& pf)
{
  pqassert(pf_new.ctype == "float");
  SelectedField<float>& f_new = *(SelectedField<float>*)pf_new.cdata;
  const SelectedField<M>& f = *(SelectedField<M>*)pf.cdata;
  convert_field_float_from_double(f_new, f);
  Py_RETURN_NONE;
}

template <class M>
PyObject* convert_double_from_float_sfield_ctype(PyField& pf_new, PyField& pf)
{
  pqassert(pf.ctype == "float");
  const SelectedField<float>& f = *(SelectedField<float>*)pf.cdata;
  SelectedField<M>& f_new = *(SelectedField<M>*)pf_new.cdata;
  convert_field_double_from_float(f_new, f);
  Py_RETURN_NONE;
}

template <class M>
PyObject* to_from_endianness_sfield_ctype(PyField& pf,
                                          const std::string& endianness_tag)
{
  SelectedField<M>& f = *(SelectedField<M>*)pf.cdata;
  if ("big_32" == endianness_tag) {
    to_from_big_endian_32(get_data(f));
  } else if ("big_64" == endianness_tag) {
    to_from_big_endian_64(get_data(f));
  } else if ("little_32" == endianness_tag) {
    to_from_little_endian_32(get_data(f));
  } else if ("little_64" == endianness_tag) {
    to_from_little_endian_64(get_data(f));
  } else {
    pqassert(false);
  }
  Py_RETURN_NONE;
}

}  // namespace qlat

EXPORT(mk_sfield, {
  using namespace qlat;
  PyObject* p_ctype = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_ctype)) {
    return NULL;
  }
  std::string ctype;
  py_convert(ctype, p_ctype);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, mk_sfield_ctype, ctype, 0);
  return p_ret;
})

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
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_fsel);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, mk_sfield_fsel_ctype, ctype, fsel, multiplicity);
  return p_ret;
})

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
})

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
})

EXPORT(set_sfield_field, {
  using namespace qlat;
  PyObject* p_sfield = NULL;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_sfield, &p_field)) {
    return NULL;
  }
  PyField psf = py_convert_field(p_sfield);
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_sfield, "fsel");
  PyField pf = py_convert_field(p_field);
  pqassert(psf.ctype == pf.ctype);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_sfield_field_ctype, pf.ctype, psf, pf, fsel);
  return p_ret;
})

EXPORT(set_field_sfield, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_sfield = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field, &p_sfield)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  PyField psf = py_convert_field(p_sfield);
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_sfield, "fsel");
  pqassert(psf.ctype == pf.ctype);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_field_sfield_ctype, pf.ctype, pf, psf, fsel);
  return p_ret;
})

EXPORT(set_sfield_sfield, {
  using namespace qlat;
  PyObject* p_sfield = NULL;
  PyObject* p_sfield0 = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_sfield, &p_sfield0)) {
    return NULL;
  }
  PyField psf = py_convert_field(p_sfield);
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_sfield, "fsel");
  PyField psf0 = py_convert_field(p_sfield0);
  const FieldSelection& fsel0 = py_convert_type<FieldSelection>(p_sfield0, "fsel");
  pqassert(psf.ctype == psf0.ctype);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_sfield_sfield_ctype, psf0.ctype, psf, psf0, fsel, fsel0);
  return p_ret;
})

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
})

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
})

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
})

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
})

EXPORT(acc_field_sfield, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_sfield = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field, &p_sfield)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_sfield);
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_sfield, "fsel");
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, acc_field_sfield_ctype, ctype, p_field, p_sfield, fsel);
  return p_ret;
})

EXPORT(field_shift_sfield, {
  using namespace qlat;
  PyObject* p_sfield_new = NULL;
  PyObject* p_sfield = NULL;
  PyObject* p_shift = NULL;
  bool is_reflect = false;
  if (!PyArg_ParseTuple(args, "OOO|b", &p_sfield_new, &p_sfield, &p_shift,
                        &is_reflect)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_sfield);
  pqassert(py_get_ctype(p_sfield_new) == py_get_ctype(p_sfield));
  const Coordinate shift = py_convert_data<Coordinate>(p_shift);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, field_shift_sfield_ctype, ctype, p_sfield_new, p_sfield,
                 shift, is_reflect);
  return p_ret;
})

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
})

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
})

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
})

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
})

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
})

EXPORT(get_elems_sfield, {
  using namespace qlat;
  PyObject* p_field = NULL;
  long idx = -1;
  if (!PyArg_ParseTuple(args, "Ol", &p_field, &idx)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, get_elems_sfield_ctype, ctype, p_field, idx);
  return p_ret;
})

EXPORT(get_elem_sfield, {
  using namespace qlat;
  PyObject* p_field = NULL;
  long idx = -1;
  long m = -1;
  if (!PyArg_ParseTuple(args, "Ol|l", &p_field, &idx, &m)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, get_elem_sfield_ctype, ctype, p_field, idx, m);
  return p_ret;
})

EXPORT(set_elems_sfield, {
  using namespace qlat;
  PyObject* p_field = NULL;
  long idx = -1;
  PyObject* p_val = NULL;
  if (!PyArg_ParseTuple(args, "OlO", &p_field, &idx, &p_val)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_elems_sfield_ctype, ctype, p_field, idx, p_val);
  return p_ret;
})

EXPORT(set_elem_sfield, {
  using namespace qlat;
  PyObject* p_field = NULL;
  long idx = -1;
  long m = -1;
  PyObject* p_val = NULL;
  if (!PyArg_ParseTuple(args, "OllO", &p_field, &idx, &m, &p_val)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_elem_sfield_ctype, ctype, p_field, idx, m, p_val);
  return p_ret;
})

EXPORT(glb_sum_tslice_double_sfield, {
  using namespace qlat;
  PyObject* p_spfield = NULL;
  PyObject* p_field = NULL;
  int t_dir = 3;
  if (!PyArg_ParseTuple(args, "OO|i", &p_spfield, &p_field, &t_dir)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_field, "fsel");
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, glb_sum_tslice_double_sfield_ctype, ctype, p_spfield,
                 p_field, fsel, t_dir);
  return p_ret;
})

EXPORT(glb_sum_tslice_long_sfield, {
  using namespace qlat;
  PyObject* p_spfield = NULL;
  PyObject* p_field = NULL;
  int t_dir = 3;
  if (!PyArg_ParseTuple(args, "OO|i", &p_spfield, &p_field, &t_dir)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_field, "fsel");
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, glb_sum_tslice_long_sfield_ctype, ctype, p_spfield,
                 p_field, fsel, t_dir);
  return p_ret;
})

EXPORT(save_sfield, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field, &p_path)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_field, "fsel");
  std::string path;
  py_convert(path, p_path);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, save_sfield_ctype, pf.ctype, pf, path, fsel);
  return p_ret;
})

EXPORT(load_sfield, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field, &p_path)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_field, "fsel");
  std::string path;
  py_convert(path, p_path);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, load_sfield_ctype, pf.ctype, pf, path, fsel);
  return p_ret;
})

EXPORT(convert_float_from_double_sfield, {
  using namespace qlat;
  PyObject* p_field_new = NULL;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field_new, &p_field)) {
    return NULL;
  }
  PyField pf_new = py_convert_field(p_field_new);
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, convert_float_from_double_sfield_ctype, pf.ctype,
                 pf_new, pf);
  return p_ret;
})

EXPORT(convert_double_from_float_sfield, {
  using namespace qlat;
  PyObject* p_field_new = NULL;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field_new, &p_field)) {
    return NULL;
  }
  PyField pf_new = py_convert_field(p_field_new);
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, convert_double_from_float_sfield_ctype, pf_new.ctype,
                 pf_new, pf);
  return p_ret;
})

EXPORT(to_from_endianness_sfield, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_endianness_tag = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field, &p_endianness_tag)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  std::string endianness_tag;
  py_convert(endianness_tag, p_endianness_tag);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, to_from_endianness_sfield_ctype, pf.ctype, pf,
                 endianness_tag);
  return p_ret;
})
