#include "lib.h"

namespace qlat
{  //

template <class M>
PyObject* write_sfw_field_ctype(ShuffledFieldsWriter& sfw, const std::string& fn, PyField& pf)
{
  const Field<M>& f = *(Field<M>*)pf.cdata;
  const long ret = write(sfw, fn, f);
  return py_convert(ret);
}

template <class M>
PyObject* read_sfr_field_ctype(ShuffledFieldsReader& sfr, const std::string& fn, PyField& pf)
{
  Field<M>& f = *(Field<M>*)pf.cdata;
  const long ret = read(sfr, fn, f);
  return py_convert(ret);
}

template <class M>
PyObject* write_sfw_sfield_ctype(ShuffledFieldsWriter& sfw,
                                 const std::string& fn, PyField& pf,
                                 const ShuffledBitSet& sbs)
{
  const SelectedField<M>& f = *(SelectedField<M>*)pf.cdata;
  const long ret = write(sfw, fn, f, sbs);
  return py_convert(ret);
}

template <class M>
PyObject* read_sfr_sfield_ctype(ShuffledFieldsReader& sfr,
                                const std::string& fn,
                                const ShuffledBitSet& sbs, PyField& pf)
{
  SelectedField<M>& f = *(SelectedField<M>*)pf.cdata;
  const long ret = read(sfr, fn, sbs, f);
  return py_convert(ret);
}

}  // namespace qlat

EXPORT(mk_sfw, {
  using namespace qlat;
  PyObject* p_path = NULL;
  PyObject* p_new_size_node = NULL;
  bool is_append = false;
  if (!PyArg_ParseTuple(args, "OO|b", &p_path, &p_new_size_node, &is_append)) {
    return NULL;
  }
  std::string path;
  py_convert(path, p_path);
  Coordinate new_size_node;
  py_convert(new_size_node, p_new_size_node);
  ShuffledFieldsWriter* psfw = new ShuffledFieldsWriter(path, new_size_node, is_append);
  return py_convert((void*)psfw);
});

EXPORT(free_sfw, {
  using namespace qlat;
  return free_obj<ShuffledFieldsWriter>(args);
});

EXPORT(get_new_size_node_sfw, {
  using namespace qlat;
  PyObject* p_obj = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_obj)) {
    return NULL;
  }
  const ShuffledFieldsWriter& obj = py_convert_type<ShuffledFieldsWriter>(p_obj);
  return py_convert(obj.new_size_node);
});

EXPORT(mk_sfr, {
  using namespace qlat;
  PyObject* p_path = NULL;
  PyObject* p_new_size_node = NULL;
  if (!PyArg_ParseTuple(args, "O|O", &p_path, &p_new_size_node)) {
    return NULL;
  }
  std::string path;
  py_convert(path, p_path);
  Coordinate new_size_node;
  if (p_new_size_node != NULL) {
    py_convert(new_size_node, p_new_size_node);
  }
  ShuffledFieldsReader* psfr = new ShuffledFieldsReader(path, new_size_node);
  return py_convert((void*)psfr);
});

EXPORT(free_sfr, {
  using namespace qlat;
  return free_obj<ShuffledFieldsReader>(args);
});

EXPORT(get_new_size_node_sfr, {
  using namespace qlat;
  PyObject* p_obj = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_obj)) {
    return NULL;
  }
  const ShuffledFieldsReader& obj = py_convert_type<ShuffledFieldsReader>(p_obj);
  return py_convert(obj.new_size_node);
});

EXPORT(mk_sbs, {
  using namespace qlat;
  PyObject* p_fsel = NULL;
  PyObject* p_new_size_node = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_fsel, &p_new_size_node)) {
    return NULL;
  }
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_fsel);
  Coordinate new_size_node;
  py_convert(new_size_node, p_new_size_node);
  ShuffledBitSet* psbs = new ShuffledBitSet();
  *psbs = mk_shuffled_bitset(fsel, new_size_node);
  return py_convert((void*)psbs);
});

EXPORT(free_sbs, {
  using namespace qlat;
  return free_obj<ShuffledBitSet>(args);
});

EXPORT(write_sfw_field, {
  using namespace qlat;
  PyObject* p_sfw = NULL;
  PyObject* p_fn = NULL;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "OOO", &p_sfw, &p_fn, &p_field)) {
    return NULL;
  }
  ShuffledFieldsWriter& sfw = py_convert_type<ShuffledFieldsWriter>(p_sfw);
  PyField pf = py_convert_field(p_field);
  std::string fn;
  py_convert(fn, p_fn);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, write_sfw_field_ctype, pf.ctype, sfw, fn, pf);
  return p_ret;
});

EXPORT(read_sfr_field, {
  using namespace qlat;
  PyObject* p_sfr = NULL;
  PyObject* p_fn = NULL;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "OOO", &p_sfr, &p_fn, &p_field)) {
    return NULL;
  }
  ShuffledFieldsReader& sfr = py_convert_type<ShuffledFieldsReader>(p_sfr);
  PyField pf = py_convert_field(p_field);
  std::string fn;
  py_convert(fn, p_fn);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, read_sfr_field_ctype, pf.ctype, sfr, fn, pf);
  return p_ret;
});

EXPORT(write_sfw_sfield, {
  using namespace qlat;
  PyObject* p_sfw = NULL;
  PyObject* p_fn = NULL;
  PyObject* p_field = NULL;
  PyObject* p_sbs = NULL;
  if (!PyArg_ParseTuple(args, "OOOO", &p_sfw, &p_fn, &p_field, &p_sbs)) {
    return NULL;
  }
  ShuffledFieldsWriter& sfw = py_convert_type<ShuffledFieldsWriter>(p_sfw);
  std::string fn;
  py_convert(fn, p_fn);
  PyField pf = py_convert_field(p_field);
  const ShuffledBitSet& sbs = py_convert_type<ShuffledBitSet>(p_sbs);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, write_sfw_sfield_ctype, pf.ctype, sfw, fn, pf, sbs);
  return p_ret;
});

EXPORT(read_sfr_sfield, {
  using namespace qlat;
  PyObject* p_sfr = NULL;
  PyObject* p_fn = NULL;
  PyObject* p_sbs = NULL;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "OOOO", &p_sfr, &p_fn, &p_sbs, &p_field)) {
    return NULL;
  }
  ShuffledFieldsReader& sfr = py_convert_type<ShuffledFieldsReader>(p_sfr);
  std::string fn;
  py_convert(fn, p_fn);
  const ShuffledBitSet& sbs = py_convert_type<ShuffledBitSet>(p_sbs);
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, read_sfr_sfield_ctype, pf.ctype, sfr, fn, sbs, pf);
  return p_ret;
});

EXPORT(list_sfr, {
  using namespace qlat;
  PyObject* p_sfr = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_sfr)) {
    return NULL;
  }
  ShuffledFieldsReader& sfr = py_convert_type<ShuffledFieldsReader>(p_sfr);
  return py_convert(list_fields(sfr));
});
