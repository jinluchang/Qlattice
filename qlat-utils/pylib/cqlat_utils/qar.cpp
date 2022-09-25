#include "lib.h"

#include <qlat-utils/qar-cache.h>

EXPORT(does_file_exist_qar, {
  using namespace qlat;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_path)) {
    return NULL;
  }
  const std::string path = py_convert_data<std::string>(p_path);
  const bool ret = does_file_exist_qar(path);
  return py_convert(ret);
})

EXPORT(does_file_or_directory_exist_qar, {
  using namespace qlat;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_path)) {
    return NULL;
  }
  const std::string path = py_convert_data<std::string>(p_path);
  const bool ret = does_file_or_directory_exist_qar(path);
  return py_convert(ret);
})

EXPORT(get_qar_multi_vol_max_size, {
  using namespace qlat;
  PyObject* p_size = NULL;
  if (!PyArg_ParseTuple(args, "|O", &p_size)) {
    return NULL;
  }
  if (NULL != p_size) {
    const long size = py_convert_data<long>(p_size);
    get_qar_multi_vol_max_size() = size;
  }
  return py_convert(get_qar_multi_vol_max_size());
})

EXPORT(qar_create, {
  using namespace qlat;
  PyObject* p_path_qar = NULL;
  PyObject* p_path_folder = NULL;
  bool is_remove_folder_after = false;
  if (!PyArg_ParseTuple(args, "OO|b", &p_path_qar, &p_path_folder, &is_remove_folder_after)) {
    return NULL;
  }
  const std::string path_qar = py_convert_data<std::string>(p_path_qar);
  const std::string path_folder = py_convert_data<std::string>(p_path_folder);
  const int ret = qar_create(path_qar, path_folder, is_remove_folder_after);
  return py_convert(ret);
})

EXPORT(qar_extract, {
  using namespace qlat;
  PyObject* p_path_qar = NULL;
  PyObject* p_path_folder = NULL;
  bool is_remove_qar_after = false;
  if (!PyArg_ParseTuple(args, "OO|b", &p_path_qar, &p_path_folder,
                        &is_remove_qar_after)) {
    return NULL;
  }
  const std::string path_qar = py_convert_data<std::string>(p_path_qar);
  const std::string path_folder = py_convert_data<std::string>(p_path_folder);
  const int ret = qar_extract(path_qar, path_folder, is_remove_qar_after);
  return py_convert(ret);
})

EXPORT(qcopy_file, {
  using namespace qlat;
  PyObject* p_path_src = NULL;
  PyObject* p_path_dst = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_path_src, &p_path_dst)) {
    return NULL;
  }
  const std::string path_src = py_convert_data<std::string>(p_path_src);
  const std::string path_dst = py_convert_data<std::string>(p_path_dst);
  const int ret = qcopy_file(path_src, path_dst);
  return py_convert(ret);
})

EXPORT(list_qar, {
  using namespace qlat;
  PyObject* p_path_qar = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_path_qar)) {
    return NULL;
  }
  const std::string path_qar = py_convert_data<std::string>(p_path_qar);
  return py_convert(list_qar(path_qar));
})

EXPORT(qcat, {
  using namespace qlat;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_path)) {
    return NULL;
  }
  const std::string path = py_convert_data<std::string>(p_path);
  const std::string ret = qcat(path);
  return py_convert(ret);
})

EXPORT(qcat_bytes, {
  using namespace qlat;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_path)) {
    return NULL;
  }
  const std::string path = py_convert_data<std::string>(p_path);
  const std::string ret = qcat(path);
  return py_convert(get_data(ret));
})
