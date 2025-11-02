#include "lib.h"

namespace qlat
{  //

template <class M>
PyObject* mk_sfield_ctype(Int dummy)
{
  (void)dummy;
  SelectedField<M>* pf = new SelectedField<M>();
  return py_convert((void*)pf);
}

template <class M>
PyObject* mk_sfield_fsel_ctype(const FieldSelection& fsel, const Int multiplicity)
{
  SelectedField<M>* pf = new SelectedField<M>();
  SelectedField<M>& f = *pf;
  qassert(multiplicity > 0);
  f.init(fsel, multiplicity);
  return py_convert((void*)pf);
}

template <class M>
PyObject* free_sfield_ctype(PyObject* p_field)
{
  SelectedField<M>& f = py_convert_type_sfield<M>(p_field);
  delete &f;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_sfield_ctype(PyObject* pf_new, PyObject* pf)
{
  SelectedField<M>& f_new = py_convert_type_sfield<M>(pf_new);
  SelectedField<M>& f = py_convert_type_sfield<M>(pf);
  f_new = f;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_sfield_field_ctype(PyObject* psf, PyObject* pf,
                                 const FieldSelection& fsel)
{
  SelectedField<M>& sf = py_convert_type_sfield<M>(psf);
  const Field<M>& f = py_convert_type_field<M>(pf);
  set_selected_field(sf, f, fsel);
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_sfield_sfield_ctype(PyObject* psf, PyObject* psf0,
                                  const FieldSelection& fsel,
                                  const FieldSelection& fsel0)
{
  SelectedField<M>& sf = py_convert_type_sfield<M>(psf);
  const SelectedField<M>& sf0 = py_convert_type_sfield<M>(psf0);
  set_selected_field(sf, sf0, fsel, fsel0);
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_field_sfield_ctype(PyObject* pf, PyObject* psf,
                                 const FieldSelection& fsel)
{
  Field<M>& f = py_convert_type_field<M>(pf);
  const SelectedField<M>& sf = py_convert_type_sfield<M>(psf);
  set_field_selected(f, sf, fsel);
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_add_sfield_ctype(PyObject* pf_new, PyObject* pf)
{
  SelectedField<M>& f_new = py_convert_type_sfield<M>(pf_new);
  SelectedField<M>& f = py_convert_type_sfield<M>(pf);
  f_new += f;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_sub_sfield_ctype(PyObject* pf_new, PyObject* pf)
{
  SelectedField<M>& f_new = py_convert_type_sfield<M>(pf_new);
  SelectedField<M>& f = py_convert_type_sfield<M>(pf);
  f_new -= f;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_mul_sfield_ctype(PyObject* pf, const ComplexD& factor)
{
  SelectedField<M>& f = py_convert_type_sfield<M>(pf);
  f *= factor;
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
PyObject* set_zero_sfield_ctype(PyObject* pf)
{
  SelectedField<M>& f = py_convert_type_sfield<M>(pf);
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
  QLAT_PUSH_DIAGNOSTIC_DISABLE_DANGLING_REF;
  FieldSelection& fsel_new = py_convert_type<FieldSelection>(p_sfield_new, "fsel");
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_sfield, "fsel");
  QLAT_DIAGNOSTIC_POP;
  SelectedField<M>& sf_new = py_convert_type_sfield<M>(p_sfield_new);
  const SelectedField<M>& sf = py_convert_type_sfield<M>(p_sfield);
  field_shift(sf_new, fsel_new, sf, fsel, shift, is_reflect);
  Py_RETURN_NONE;
}

template <class M>
PyObject* get_n_elems_sfield_ctype(PyObject* pf)
{
  SelectedField<M>& sf = py_convert_type_sfield<M>(pf);
  const Long ret = sf.n_elems;
  return py_convert(ret);
}

template <class M>
PyObject* get_total_site_sfield_ctype(PyObject* pf)
{
  SelectedField<M>& f = py_convert_type_sfield<M>(pf);
  const Coordinate ret = f.geo().total_site();
  return py_convert(ret);
}

template <class M>
PyObject* get_multiplicity_sfield_ctype(PyObject* pf)
{
  SelectedField<M>& sf = py_convert_type_sfield<M>(pf);
  const Long ret = sf.multiplicity;
  return py_convert(ret);
}

template <class M>
PyObject* set_geo_sfield_ctype(Geometry& geo, PyObject* pf)
{
  const SelectedField<M>& sf = py_convert_type_sfield<M>(pf);
  geo = sf.geo();
  Py_RETURN_NONE;
}

template <class M>
PyObject* qnorm_sfield_ctype(PyObject* pf)
{
  SelectedField<M>& f = py_convert_type_sfield<M>(pf);
  const RealD ret = qnorm(f);
  return py_convert(ret);
}

template <class M>
PyObject* qnorm_field_sfield_ctype(SelectedField<RealD>& f, PyObject* p_field1)
{
  const SelectedField<M>& f1 = py_convert_type_sfield<M>(p_field1);
  qnorm_field(f, f1);
  Py_RETURN_NONE;
}

template <class M>
PyObject* get_elems_sfield_ctype(PyObject* p_sfield, const Long idx)
{
  const SelectedField<M>& f = py_convert_type_sfield<M>(p_sfield);
  return py_convert(f.get_elems_const(idx));
}

template <class M>
PyObject* get_elem_sfield_ctype(PyObject* p_sfield, const Long idx, const Int m)
{
  const SelectedField<M>& f = py_convert_type_sfield<M>(p_sfield);
  if (m >= 0) {
    return py_convert(f.get_elem(idx, m));
  } else {
    return py_convert(f.get_elem(idx));
  }
}

template <class M>
PyObject* set_elems_sfield_ctype(PyObject* p_field, const Long idx,
                                 PyObject* p_val)
{
  SelectedField<M>& f = py_convert_type_sfield<M>(p_field);
  const Int multiplicity = f.multiplicity;
  qassert((Long)PyBytes_Size(p_val) == (Long)multiplicity * (Long)sizeof(M));
  const Vector<M> val((M*)PyBytes_AsString(p_val), multiplicity);
  assign(f.get_elems(idx), val);
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_elem_sfield_ctype(PyObject* p_field, const Long idx,
                                const Int m, PyObject* p_val)
{
  SelectedField<M>& f = py_convert_type_sfield<M>(p_field);
  qassert(PyBytes_Size(p_val) == sizeof(M));
  const M& val = *(M*)PyBytes_AsString(p_val);
  f.get_elem(idx, m) = val;
  Py_RETURN_NONE;
}

template <class M>
PyObject* glb_sum_tslice_double_sfield_ctype(PyObject* p_spfield,
                                             PyObject* p_field,
                                             const FieldSelection& fsel,
                                             const Int t_dir)
{
  SelectedPoints<M>& sp = py_convert_type_spoints<M>(p_spfield);
  const SelectedField<M>& f = py_convert_type_sfield<M>(p_field);
  field_glb_sum_tslice(sp, f, fsel, t_dir);
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

template <class M>
PyObject* save_sfield_ctype(PyObject* pf, const std::string& path,
                            const FieldSelection& fsel)
{
  const SelectedField<M>& f = py_convert_type_sfield<M>(pf);
  const Long ret = write_selected_field(f, path, fsel);
  return py_convert(ret);
}

template <class M>
PyObject* load_sfield_ctype(PyObject* pf, const std::string& path,
                            const FieldSelection& fsel)
{
  SelectedField<M>& f = py_convert_type_sfield<M>(pf);
  const Long ret = read_selected_field(f, path, fsel);
  return py_convert(ret);
}

template <class M>
PyObject* convert_float_from_double_sfield_ctype(PyObject* pf_new, PyObject* pf)
{
  SelectedField<RealF>& f_new = py_convert_type_sfield<RealF>(pf_new);
  const SelectedField<M>& f = py_convert_type_sfield<M>(pf);
  convert_field_float_from_double(f_new, f);
  Py_RETURN_NONE;
}

template <class M>
PyObject* convert_double_from_float_sfield_ctype(PyObject* pf_new, PyObject* pf)
{
  const SelectedField<RealF>& f = py_convert_type_sfield<RealF>(pf);
  SelectedField<M>& f_new = py_convert_type_sfield<M>(pf_new);
  convert_field_double_from_float(f_new, f);
  Py_RETURN_NONE;
}

template <class M>
PyObject* to_from_endianness_sfield_ctype(PyObject* pf,
                                          const std::string& endianness_tag)
{
  SelectedField<M>& f = py_convert_type_sfield<M>(pf);
  if ("big_32" == endianness_tag) {
    to_from_big_endian(get_data(f));
  } else if ("big_64" == endianness_tag) {
    to_from_big_endian(get_data(f));
  } else if ("little_32" == endianness_tag) {
    to_from_little_endian(get_data(f));
  } else if ("little_64" == endianness_tag) {
    to_from_little_endian(get_data(f));
  } else {
    qassert(false);
  }
  Py_RETURN_NONE;
}

}  // namespace qlat


EXPORT(set_add_sfield, {
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

EXPORT(set_sub_sfield, {
  using namespace qlat;
  PyObject* p_field_new = NULL;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field_new, &p_field)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  qassert(py_get_ctype(p_field_new) == ctype);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_sub_sfield_ctype, ctype, p_field_new, p_field);
  return p_ret;
})

EXPORT(set_mul_double_sfield, {
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

EXPORT(set_zero_sfield, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_zero_sfield_ctype, ctype, p_field);
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
  qassert(py_get_ctype(p_sfield) == ctype);
  QLAT_PUSH_DIAGNOSTIC_DISABLE_DANGLING_REF;
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_sfield, "fsel");
  QLAT_DIAGNOSTIC_POP;
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
  qassert(py_get_ctype(p_sfield_new) == ctype);
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
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, get_n_elems_sfield_ctype, ctype, p_field);
  return p_ret;
})

EXPORT(get_total_site_sfield, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, get_total_site_sfield_ctype, ctype, p_field);
  return p_ret;
})

EXPORT(get_multiplicity_sfield, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, get_multiplicity_sfield_ctype, ctype, p_field);
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
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_geo_sfield_ctype, ctype, geo, p_field);
  return p_ret;
})

EXPORT(qnorm_sfield, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, qnorm_sfield_ctype, ctype, p_field);
  return p_ret;
})

EXPORT(qnorm_field_sfield, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_field1 = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field, &p_field1)) {
    return NULL;
  }
  SelectedField<RealD>& f = py_convert_type_sfield<RealD>(p_field);
  const std::string ctype = py_get_ctype(p_field1);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, qnorm_field_sfield_ctype, ctype, f, p_field1);
  return p_ret;
})

EXPORT(set_sqrt_double_sfield, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_field1 = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field, &p_field1)) {
    return NULL;
  }
  SelectedField<RealD>& f = py_convert_type_sfield<RealD>(p_field);
  const SelectedField<RealD>& f1 = py_convert_type_sfield<RealD>(p_field1);
  const Geometry& geo = f1.geo();
  qassert(geo.is_only_local);
  f.init();
  f.init(geo, f1.n_elems, f1.multiplicity);
  qacc_for(idx, f.n_elems, {
    const Vector<RealD> f1v = f1.get_elems_const(idx);
    Vector<RealD> fv = f.get_elems(idx);
    for (Int m = 0; m < f.multiplicity; ++m) {
      fv[m] = std::sqrt(f1v[m]);
    }
  });
  Py_RETURN_NONE;
})

EXPORT(get_elems_sfield, {
  using namespace qlat;
  PyObject* p_field = NULL;
  Long idx = -1;
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
  Long idx = -1;
  Long m = -1;
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
  Long idx = -1;
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
  Long idx = -1;
  Long m = -1;
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
  FIELD_DISPATCH(p_ret, glb_sum_tslice_double_sfield_ctype, ctype, p_spfield,
                 p_field, fsel, t_dir);
  return p_ret;
})

EXPORT(glb_sum_tslice_long_sfield, {
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

EXPORT(save_sfield, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field, &p_path)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  QLAT_PUSH_DIAGNOSTIC_DISABLE_DANGLING_REF;
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_field, "fsel");
  QLAT_DIAGNOSTIC_POP;
  const std::string path = py_convert_data<std::string>(p_path);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, save_sfield_ctype, ctype, p_field, path, fsel);
  return p_ret;
})

EXPORT(load_sfield, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field, &p_path)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  QLAT_PUSH_DIAGNOSTIC_DISABLE_DANGLING_REF;
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_field, "fsel");
  QLAT_DIAGNOSTIC_POP;
  const std::string path = py_convert_data<std::string>(p_path);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, load_sfield_ctype, ctype, p_field, path, fsel);
  return p_ret;
})

EXPORT(convert_float_from_double_sfield, {
  using namespace qlat;
  PyObject* p_field_new = NULL;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field_new, &p_field)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, convert_float_from_double_sfield_ctype, ctype,
                 p_field_new, p_field);
  return p_ret;
})

EXPORT(convert_double_from_float_sfield, {
  using namespace qlat;
  PyObject* p_field_new = NULL;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field_new, &p_field)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field_new);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, convert_double_from_float_sfield_ctype, ctype,
                 p_field_new, p_field);
  return p_ret;
})

EXPORT(to_from_endianness_sfield, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_endianness_tag = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field, &p_endianness_tag)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  const std::string endianness_tag =
      py_convert_data<std::string>(p_endianness_tag);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, to_from_endianness_sfield_ctype, ctype, p_field,
                 endianness_tag);
  return p_ret;
})
