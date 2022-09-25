#include "lib.h"

#include <qlat-utils/lat-io.h>

EXPORT(mk_lat_data, {
  using namespace qlat;
  LatData* pld = new LatData();
  return py_convert((void*)pld);
})

EXPORT(free_lat_data, {
  using namespace qlat;
  return free_obj<LatData>(args);
})

EXPORT(set_lat_data, {
  using namespace qlat;
  return set_obj<LatData>(args);
})

EXPORT(set_zero_lat_data, {
  using namespace qlat;
  PyObject* p_ld = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_ld)) {
    return NULL;
  }
  LatData& ld = py_convert_type<LatData>(p_ld);
  set_zero(ld);
  Py_RETURN_NONE;
})

EXPORT(show_lat_data, {
  using namespace qlat;
  PyObject* p_ld = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_ld)) {
    return NULL;
  }
  LatData& ld = py_convert_type<LatData>(p_ld);
  return py_convert(show(ld));
})

EXPORT(qnorm_lat_data, {
  using namespace qlat;
  PyObject* p_ld = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_ld)) {
    return NULL;
  }
  LatData& ld = py_convert_type<LatData>(p_ld);
  return py_convert(qnorm(ld));
})

EXPORT(load_lat_data, {
  using namespace qlat;
  PyObject* p_ld = NULL;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_ld, &p_path)) {
    return NULL;
  }
  LatData& ld = py_convert_type<LatData>(p_ld);
  std::string path;
  py_convert(path, p_path);
  ld.load(path);
  Py_RETURN_NONE;
})

EXPORT(save_lat_data, {
  using namespace qlat;
  PyObject* p_ld = NULL;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_ld, &p_path)) {
    return NULL;
  }
  const LatData& ld = py_convert_type<LatData>(p_ld);
  std::string path;
  py_convert(path, p_path);
  ld.save(path);
  Py_RETURN_NONE;
})

EXPORT(is_matching_lat_data, {
  using namespace qlat;
  PyObject* p_ld = NULL;
  PyObject* p_ld1 = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_ld, &p_ld1)) {
    return NULL;
  }
  LatData& ld = py_convert_type<LatData>(p_ld);
  LatData& ld1 = py_convert_type<LatData>(p_ld1);
  return py_convert(is_matching(ld, ld1));
})

EXPORT(is_complex_lat_data, {
  using namespace qlat;
  PyObject* p_ld = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_ld)) {
    return NULL;
  }
  LatData& ld = py_convert_type<LatData>(p_ld);
  const bool is_complex_ld = is_lat_info_complex(ld.info);
  return py_convert(is_complex_ld);
})

EXPORT(get_ndim_lat_data, {
  using namespace qlat;
  PyObject* p_ld = NULL;
  bool is_always_double = false;
  if (!PyArg_ParseTuple(args, "O|b", &p_ld, &is_always_double)) {
    return NULL;
  }
  LatData& ld = py_convert_type<LatData>(p_ld);
  long ndim = ld.info.size();
  if (not is_always_double) {
    const bool is_complex_ld = is_lat_info_complex(ld.info);
    if (is_complex_ld) {
      ndim -= 1;
    }
  }
  return py_convert(ndim);
})

EXPORT(get_dim_sizes_lat_data, {
  using namespace qlat;
  PyObject* p_ld = NULL;
  bool is_always_double = false;
  if (!PyArg_ParseTuple(args, "O|b", &p_ld, &is_always_double)) {
    return NULL;
  }
  LatData& ld = py_convert_type<LatData>(p_ld);
  long ndim = ld.info.size();
  if (not is_always_double) {
    const bool is_complex_ld = is_lat_info_complex(ld.info);
    if (is_complex_ld) {
      ndim -= 1;
    }
  }
  std::vector<long> dim_sizes(ndim);
  for (long i = 0; i < ndim; ++i) {
    dim_sizes[i] = ld.info[i].size;
  }
  return py_convert(dim_sizes);
})

EXPORT(get_dim_name_lat_data, {
  using namespace qlat;
  PyObject* p_ld = NULL;
  long dim = 0;
  if (!PyArg_ParseTuple(args, "Ol", &p_ld, &dim)) {
    return NULL;
  }
  LatData& ld = py_convert_type<LatData>(p_ld);
  pqassert(0 <= dim and dim < (long)ld.info.size());
  return py_convert(ld.info[dim].name);
})

EXPORT(get_dim_size_lat_data, {
  using namespace qlat;
  PyObject* p_ld = NULL;
  long dim = 0;
  if (!PyArg_ParseTuple(args, "Ol", &p_ld, &dim)) {
    return NULL;
  }
  LatData& ld = py_convert_type<LatData>(p_ld);
  pqassert(0 <= dim and dim < (long)ld.info.size());
  return py_convert(ld.info[dim].size);
})

EXPORT(get_dim_indices_lat_data, {
  using namespace qlat;
  PyObject* p_ld = NULL;
  long dim = 0;
  if (!PyArg_ParseTuple(args, "Ol", &p_ld, &dim)) {
    return NULL;
  }
  LatData& ld = py_convert_type<LatData>(p_ld);
  pqassert(0 <= dim and dim < (long)ld.info.size());
  return py_convert(ld.info[dim].indices);
})

EXPORT(set_dim_sizes_lat_data, {
  using namespace qlat;
  PyObject* p_ld = NULL;
  PyObject* p_dim_sizes = NULL;
  bool is_complex = true;
  if (!PyArg_ParseTuple(args, "OO|b", &p_ld, &p_dim_sizes, &is_complex)) {
    return NULL;
  }
  LatData& ld = py_convert_type<LatData>(p_ld);
  std::vector<long> dim_sizes;
  py_convert(dim_sizes, p_dim_sizes);
  clear(ld.info);
  if (is_complex) {
    ld.info.resize(dim_sizes.size() + 1);
    ld.info.back() = lat_dim_re_im();
  } else {
    ld.info.resize(dim_sizes.size());
  }
  for (long i = 0; i < (long)dim_sizes.size(); ++i) {
    ld.info[i].size = dim_sizes[i];
  }
  lat_data_alloc(ld);
  Py_RETURN_NONE;
})

EXPORT(set_dim_name_lat_data, {
  using namespace qlat;
  PyObject* p_ld = NULL;
  long dim = 0;
  PyObject* p_name = NULL;
  PyObject* p_indices = NULL;
  if (!PyArg_ParseTuple(args, "OlO|O", &p_ld, &dim, &p_name, &p_indices)) {
    return NULL;
  }
  LatData& ld = py_convert_type<LatData>(p_ld);
  long ndim = ld.info.size();
  pqassert(0 <= dim and dim < ndim);
  py_convert(ld.info[dim].name, p_name);
  if (NULL == p_indices) {
    clear(ld.info[dim].indices);
  } else {
    py_convert(ld.info[dim].indices, p_indices);
  }
  Py_RETURN_NONE;
})

EXPORT(peek_lat_data, {
  using namespace qlat;
  PyObject* p_ld = NULL;
  PyObject* p_idx = NULL;
  bool is_always_double = false;
  if (!PyArg_ParseTuple(args, "O|Ob", &p_ld, &p_idx, &is_always_double)) {
    return NULL;
  }
  LatData& ld = py_convert_type<LatData>(p_ld);
  std::vector<long> idx;
  if (p_idx != NULL) {
    py_convert(idx, p_idx);
  }
  const bool is_complex_ld = (not is_always_double) and is_lat_info_complex(ld.info);
  if (is_complex_ld) {
    return py_convert(lat_data_cget_const(ld, idx));
  } else {
    return py_convert(lat_data_get_const(ld, idx));
  }
})

EXPORT(poke_lat_data, {
  using namespace qlat;
  PyObject* p_ld = NULL;
  PyObject* p_idx = NULL;
  PyObject* p_val = NULL;
  bool is_always_double = false;
  if (!PyArg_ParseTuple(args, "OOO|b", &p_ld, &p_idx, &p_val,
                        &is_always_double)) {
    return NULL;
  }
  LatData& ld = py_convert_type<LatData>(p_ld);
  std::vector<long> idx;
  py_convert(idx, p_idx);
  const bool is_complex_ld = (not is_always_double) and is_lat_info_complex(ld.info);
  if (is_complex_ld) {
    py_convert(lat_data_cget(ld, idx), p_val);
  } else {
    py_convert(lat_data_get(ld, idx), p_val);
  }
  Py_RETURN_NONE;
})

EXPORT(set_add_lat_data, {
  using namespace qlat;
  PyObject* p_ld_new = NULL;
  PyObject* p_ld = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_ld_new, &p_ld)) {
    return NULL;
  }
  LatData& ld_new = py_convert_type<LatData>(p_ld_new);
  const LatData& ld = py_convert_type<LatData>(p_ld);
  ld_new += ld;
  Py_RETURN_NONE;
})

EXPORT(set_sub_lat_data, {
  using namespace qlat;
  PyObject* p_ld_new = NULL;
  PyObject* p_ld = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_ld_new, &p_ld)) {
    return NULL;
  }
  LatData& ld_new = py_convert_type<LatData>(p_ld_new);
  const LatData& ld = py_convert_type<LatData>(p_ld);
  ld_new -= ld;
  Py_RETURN_NONE;
})

EXPORT(set_mul_double_lat_data, {
  using namespace qlat;
  PyObject* p_ld_new = NULL;
  double factor = 0.0;
  if (!PyArg_ParseTuple(args, "Od", &p_ld_new, &factor)) {
    return NULL;
  }
  LatData& ld_new = py_convert_type<LatData>(p_ld_new);
  ld_new *= factor;
  Py_RETURN_NONE;
})
