#include "exceptions.h"
#include "convert.h"

namespace qlat
{  //

template <class M>
PyObject* mk_field_type(Geometry* p_geo, const int multiplicty)
{
  Field<M>* p_f = new Field<M>();
  Field<M>& f = *p_f;
  if (p_geo != NULL) {
    Geometry& geo = *p_geo;
    if (multiplicty == 0) {
      f.init(geo);
    } else {
      f.init(geo, multiplicty);
    }
  }
  return PyLong_FromVoidPtr(p_f);
}

template <class M>
void free_field_type(void* p_field)
{
  Field<M>* p_f = (Field<M>*)p_field;
  delete p_f;
}

}  // namespace qlat

EXPORT(mk_field, {
  using namespace qlat;
  PyObject* p_type = NULL;
  Geometry* p_geo = NULL;
  int multiplicty = 0;
  if (!PyArg_ParseTuple(args, "O|li", &p_type, &p_geo, &multiplicty)) {
    return NULL;
  }
  std::string type;
  py_convert(type, p_type);
  PyObject* p_field;
  if ("ColorMatrix" == type) {
    p_field = mk_field_type<ColorMatrix>(p_geo, multiplicty);
  } else if ("WilsonMatrix" == type) {
    p_field = mk_field_type<WilsonMatrix>(p_geo, multiplicty);
  } else if ("WilsonVector" == type) {
    p_field = mk_field_type<WilsonVector>(p_geo, multiplicty);
  } else if ("double" == type) {
    p_field = mk_field_type<double>(p_geo, multiplicty);
  } else if ("Complex" == type) {
    p_field = mk_field_type<Complex>(p_geo, multiplicty);
  } else {
    pqerr("mk_field type='%s' does not exist.", type.c_str());
  }
  return p_field;
});

EXPORT(free_field, {
  using namespace qlat;
  PyObject* p_type = NULL;
  void* p_field = NULL;
  if (!PyArg_ParseTuple(args, "Ol", &p_type, &p_field)) {
    return NULL;
  }
  std::string type;
  py_convert(type, p_type);
  if ("ColorMatrix" == type) {
    free_field_type<ColorMatrix>(p_field);
  } else if ("WilsonMatrix" == type) {
    free_field_type<WilsonMatrix>(p_field);
  } else if ("WilsonVector" == type) {
    free_field_type<WilsonVector>(p_field);
  } else if ("double" == type) {
    free_field_type<double>(p_field);
  } else if ("Complex" == type) {
    free_field_type<Complex>(p_field);
  } else {
    pqerr("mk_field type='%s' does not exist.", type.c_str());
  }
  return PyLong_FromLong(0);
});
