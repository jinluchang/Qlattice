#include "lib.h"

EXPORT(qmkdir_sync_node, {
  using namespace qlat;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_path)) {
    return NULL;
  }
  const std::string path = py_convert_data<std::string>(p_path);
  const int ret = qmkdir_sync_node(path);
  return py_convert(ret);
})

EXPORT(obtain_lock, {
  using namespace qlat;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_path)) {
    return NULL;
  }
  const std::string path = py_convert_data<std::string>(p_path);
  const bool ret = obtain_lock(path);
  return py_convert(ret);
})

EXPORT(release_lock, {
  using namespace qlat;
  release_lock();
  Py_RETURN_NONE;
})

EXPORT(does_file_exist_sync_node, {
  using namespace qlat;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_path)) {
    return NULL;
  }
  const std::string path = py_convert_data<std::string>(p_path);
  const bool ret = does_file_exist_sync_node(path);
  return py_convert(ret);
})

EXPORT(does_regular_file_exist_qar_sync_node, {
  using namespace qlat;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_path)) {
    return NULL;
  }
  const std::string path = py_convert_data<std::string>(p_path);
  const bool ret = does_regular_file_exist_qar_sync_node(path);
  return py_convert(ret);
})

EXPORT(does_file_exist_qar_sync_node, {
  using namespace qlat;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_path)) {
    return NULL;
  }
  const std::string path = py_convert_data<std::string>(p_path);
  const bool ret = does_file_exist_qar_sync_node(path);
  return py_convert(ret);
})

EXPORT(is_directory_sync_node, {
  using namespace qlat;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_path)) {
    return NULL;
  }
  const std::string path = py_convert_data<std::string>(p_path);
  const bool ret = is_directory_sync_node(path);
  return py_convert(ret);
})

EXPORT(is_regular_file_sync_node, {
  using namespace qlat;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_path)) {
    return NULL;
  }
  const std::string path = py_convert_data<std::string>(p_path);
  const bool ret = is_regular_file_sync_node(path);
  return py_convert(ret);
})

EXPORT(qls_sync_node, {
  using namespace qlat;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_path)) {
    return NULL;
  }
  const std::string path = py_convert_data<std::string>(p_path);
  const std::vector<std::string> ret = qls_sync_node(path);
  return py_convert(ret);
})

EXPORT(qls_all_sync_node, {
  using namespace qlat;
  PyObject* p_path = NULL;
  bool is_folder_before_files = false;
  if (!PyArg_ParseTuple(args, "O|b", &p_path, &is_folder_before_files)) {
    return NULL;
  }
  const std::string path = py_convert_data<std::string>(p_path);
  const std::vector<std::string> ret =
      qls_all_sync_node(path, is_folder_before_files);
  return py_convert(ret);
})

EXPORT(qcat_sync_node, {
  using namespace qlat;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_path)) {
    return NULL;
  }
  const std::string path = py_convert_data<std::string>(p_path);
  const std::string ret = qcat_sync_node(path);
  return py_convert(ret);
})

EXPORT(qcat_bytes_sync_node, {
  using namespace qlat;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_path)) {
    return NULL;
  }
  const std::string path = py_convert_data<std::string>(p_path);
  const std::string ret = qcat_sync_node(path);
  return py_convert(get_data(ret));
})

EXPORT(qar_create_info, {
  using namespace qlat;
  PyObject* p_path_qar = NULL;
  PyObject* p_path_folder = NULL;
  bool is_remove_folder_after = false;
  if (!PyArg_ParseTuple(args, "OO|b", &p_path_qar, &p_path_folder, &is_remove_folder_after)) {
    return NULL;
  }
  const std::string path_qar = py_convert_data<std::string>(p_path_qar);
  const std::string path_folder = py_convert_data<std::string>(p_path_folder);
  const int ret = qar_create_info(path_qar, path_folder, is_remove_folder_after);
  return py_convert(ret);
})

EXPORT(qar_extract_info, {
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
  const int ret = qar_extract_info(path_qar, path_folder, is_remove_qar_after);
  return py_convert(ret);
})

EXPORT(qcopy_file_info, {
  using namespace qlat;
  PyObject* p_path_src = NULL;
  PyObject* p_path_dst = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_path_src, &p_path_dst)) {
    return NULL;
  }
  const std::string path_src = py_convert_data<std::string>(p_path_src);
  const std::string path_dst = py_convert_data<std::string>(p_path_dst);
  const int ret = qcopy_file_info(path_src, path_dst);
  return py_convert(ret);
})

EXPORT(qload_datatable_sync_node, {
  using namespace qlat;
  PyObject* p_path = NULL;
  bool is_par = false;
  if (!PyArg_ParseTuple(args, "O|b", &p_path, &is_par)) {
    return NULL;
  }
  const std::string path = py_convert_data<std::string>(p_path);
  const DataTable dt = qload_datatable_sync_node(path, is_par);
  return py_convert(dt);
})

EXPORT(check_time_limit, {
  using namespace qlat;
  PyObject* p_budget = NULL;
  if (!PyArg_ParseTuple(args, "|O", &p_budget)) {
    return NULL;
  }
  if (NULL == p_budget) {
    check_time_limit();
  } else {
    const double budget = py_convert_data<double>(p_budget);
    check_time_limit(budget);
  }
  Py_RETURN_NONE;
})

EXPORT(check_stop, {
  using namespace qlat;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "|O", &p_path)) {
    return NULL;
  }
  if (NULL == p_path) {
    check_stop();
  } else {
    const std::string path = py_convert_data<std::string>(p_path);
    check_stop(path);
  }
  Py_RETURN_NONE;
})

EXPORT(get_time_limit, {
  using namespace qlat;
  PyObject* p_time_limit = NULL;
  if (!PyArg_ParseTuple(args, "|O", &p_time_limit)) {
    return NULL;
  }
  if (NULL != p_time_limit) {
    const double time_limit = py_convert_data<double>(p_time_limit);
    get_time_limit() = time_limit;
  }
  return py_convert(get_time_limit());
})

EXPORT(get_default_budget, {
  using namespace qlat;
  PyObject* p_budget = NULL;
  if (!PyArg_ParseTuple(args, "|O", &p_budget)) {
    return NULL;
  }
  if (NULL != p_budget) {
    const double budget = py_convert_data<double>(p_budget);
    get_default_budget() = budget;
  }
  return py_convert(get_default_budget());
})

EXPORT(qquit, {
  using namespace qlat;
  PyObject* p_msg = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_msg)) {
    return NULL;
  }
  const std::string msg = py_convert_data<std::string>(p_msg);
  qquit(msg);
  Py_RETURN_NONE;
})
