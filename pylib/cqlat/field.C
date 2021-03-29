#include "lib.h"

namespace qlat
{  //

template <class M>
PyObject* mk_field_ctype(PyObject* p_geo, const int multiplicity)
{
  Field<M>* pf = new Field<M>();
  Field<M>& f = *pf;
  if (p_geo != NULL) {
    const Geometry& geo = py_convert_type<Geometry>(p_geo);
    if (multiplicity == 0) {
      f.init(geo);
    } else {
      pqassert(multiplicity > 0);
      f.init(geo, multiplicity);
    }
  }
  return py_convert((void*)pf);
}

template <class M>
PyObject* free_field_ctype(PyField& pf)
{
  Field<M>& f = *(Field<M>*)pf.cdata;
  delete &f;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_field_ctype(PyField& pf_new, PyField& pf)
{
  Field<M>& f_new = *(Field<M>*)pf_new.cdata;
  const Field<M>& f = *(Field<M>*)pf.cdata;
  f_new = f;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_add_field_ctype(PyField& pf_new, PyField& pf)
{
  Field<M>& f_new = *(Field<M>*)pf_new.cdata;
  Field<M>& f = *(Field<M>*)pf.cdata;
  f_new += f;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_sub_field_ctype(PyField& pf_new, PyField& pf)
{
  Field<M>& f_new = *(Field<M>*)pf_new.cdata;
  Field<M>& f = *(Field<M>*)pf.cdata;
  f_new -= f;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_mul_field_ctype(PyField& pf, const Complex& factor)
{
  Field<M>& f = *(Field<M>*)pf.cdata;
  f *= factor;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_mul_field_ctype(PyField& pf, const double& factor)
{
  Field<M>& f = *(Field<M>*)pf.cdata;
  f *= factor;
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_zero_field_ctype(PyField& pf)
{
  Field<M>& f = *(Field<M>*)pf.cdata;
  set_zero(f);
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_unit_field_ctype(PyField& pf, const Complex& coef)
{
  Field<M>& f = *(Field<M>*)pf.cdata;
  set_unit(f, coef);
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_u_rand_double_field_ctype(PyField& pf, const RngState& rs,
                                        const double upper, const double lower)
{
  Field<M>& f = *(Field<M>*)pf.cdata;
  set_u_rand_double(f, rs, upper, lower);
  Py_RETURN_NONE;
}

template <class M>
PyObject* get_total_site_field_ctype(PyField& pf)
{
  Field<M>& f = *(Field<M>*)pf.cdata;
  const Coordinate ret = f.geo().total_site();
  return py_convert(ret);
}

template <class M>
PyObject* get_multiplicity_field_ctype(PyField& pf)
{
  Field<M>& f = *(Field<M>*)pf.cdata;
  const long ret = f.geo().multiplicity;
  return py_convert(ret);
}

template <class M>
PyObject* set_geo_field_ctype(Geometry& geo, PyField& pf)
{
  const Field<M>& f = *(Field<M>*)pf.cdata;
  geo = f.geo();
  Py_RETURN_NONE;
}

template <class M>
PyObject* qnorm_field_ctype(PyField& pf)
{
  Field<M>& f = *(Field<M>*)pf.cdata;
  const double ret = qnorm(f);
  return py_convert(ret);
}

template <class M>
PyObject* crc32_field_ctype(PyField& pf)
{
  const Field<M>& f = *(Field<M>*)pf.cdata;
  const crc32_t ret = field_crc32(f);
  return py_convert((long)ret);
}

template <class M>
PyObject* get_mview_field_ctype(PyField& pf, PyObject* p_field)
{
  Field<M>& f = *(Field<M>*)pf.cdata;
  Vector<M> fv = get_data(f);
  PyObject* p_mview = py_convert_mview(fv);
  pqassert(p_field != NULL);
  Py_INCREF(p_field);
  pqassert(!((PyMemoryViewObject*)p_mview)->mbuf->master.obj);
  ((PyMemoryViewObject*)p_mview)->mbuf->master.obj = p_field;
  return p_mview;
}

}  // namespace qlat

EXPORT(mk_field, {
  using namespace qlat;
  PyObject* p_ctype = NULL;
  PyObject* p_geo = NULL;
  int multiplicity = 0;
  if (!PyArg_ParseTuple(args, "O|Oi", &p_ctype, &p_geo, &multiplicity)) {
    return NULL;
  }
  std::string ctype;
  py_convert(ctype, p_ctype);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, mk_field_ctype, ctype, p_geo, multiplicity);
  return p_ret;
});

EXPORT(free_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, free_field_ctype, pf.ctype, pf);
  return p_ret;
});

EXPORT(set_field, {
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
  FIELD_DISPATCH(p_ret, set_field_ctype, pf.ctype, pf_new, pf);
  return p_ret;
});

EXPORT(set_add_field, {
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
  FIELD_DISPATCH(p_ret, set_add_field_ctype, pf.ctype, pf_new, pf);
  return p_ret;
});

EXPORT(set_sub_field, {
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
  FIELD_DISPATCH(p_ret, set_sub_field_ctype, pf.ctype, pf_new, pf);
  return p_ret;
});

EXPORT(set_mul_double_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  double factor = 0.0;
  if (!PyArg_ParseTuple(args, "Od", &p_field, &factor)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_mul_field_ctype, pf.ctype, pf, factor);
  return p_ret;
});

EXPORT(set_zero_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_zero_field_ctype, pf.ctype, pf);
  return p_ret;
});

EXPORT(set_unit_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  Complex coef = 1.0;
  if (!PyArg_ParseTuple(args, "O|D", &p_field, &coef)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_unit_field_ctype, pf.ctype, pf, coef);
  return p_ret;
});

EXPORT(set_u_rand_double_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_rng = NULL;
  double upper = 1.0;
  double lower = 0.0;
  if (!PyArg_ParseTuple(args, "OO|dd", &p_field, &p_rng, &upper, &lower)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  const RngState& rng = py_convert_type<RngState>(p_rng);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_u_rand_double_field_ctype, pf.ctype, pf, rng, upper,
                 lower);
  return p_ret;
});

EXPORT(get_total_site_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, get_total_site_field_ctype, pf.ctype, pf);
  return p_ret;
});

EXPORT(get_multiplicity_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, get_multiplicity_field_ctype, pf.ctype, pf);
  return p_ret;
});

EXPORT(set_geo_field, {
  using namespace qlat;
  PyObject* p_geo = NULL;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_geo, &p_field)) {
    return NULL;
  }
  Geometry& geo = py_convert_type<Geometry>(p_geo);
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_geo_field_ctype, pf.ctype, geo, pf);
  return p_ret;
});

EXPORT(qnorm_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, qnorm_field_ctype, pf.ctype, pf);
  return p_ret;
});

EXPORT(crc32_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, crc32_field_ctype, pf.ctype, pf);
  return p_ret;
});

EXPORT(get_mview_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, get_mview_field_ctype, pf.ctype, pf, p_field);
  return p_ret;
});


