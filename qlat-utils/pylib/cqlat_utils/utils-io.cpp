#include "lib.h"

#include <qlat-utils/qar-cache.h>

EXPORT(flush, {
  using namespace qlat;
  fflush(get_output_file());
  Py_RETURN_NONE;
})

EXPORT(qremove, {
  using namespace qlat;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_path)) {
    return NULL;
  }
  const std::string path = py_convert_data<std::string>(p_path);
  const int ret = qremove(path);
  return py_convert(ret);
})

EXPORT(qremove_all, {
  using namespace qlat;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_path)) {
    return NULL;
  }
  const std::string path = py_convert_data<std::string>(p_path);
  const int ret = qremove_all(path);
  return py_convert(ret);
})

EXPORT(qmkdir, {
  using namespace qlat;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_path)) {
    return NULL;
  }
  const std::string path = py_convert_data<std::string>(p_path);
  const int ret = qmkdir(path);
  return py_convert(ret);
})

EXPORT(qmkdir_info, {
  using namespace qlat;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_path)) {
    return NULL;
  }
  const std::string path = py_convert_data<std::string>(p_path);
  const int ret = qmkdir_info(path);
  return py_convert(ret);
})

EXPORT(does_file_exist, {
  using namespace qlat;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_path)) {
    return NULL;
  }
  const std::string path = py_convert_data<std::string>(p_path);
  const bool ret = does_file_exist(path);
  return py_convert(ret);
})

EXPORT(is_directory, {
  using namespace qlat;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_path)) {
    return NULL;
  }
  const std::string path = py_convert_data<std::string>(p_path);
  const bool ret = is_directory(path);
  return py_convert(ret);
})

EXPORT(is_regular_file, {
  using namespace qlat;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_path)) {
    return NULL;
  }
  const std::string path = py_convert_data<std::string>(p_path);
  const bool ret = is_regular_file(path);
  return py_convert(ret);
})

EXPORT(qtouch, {
  using namespace qlat;
  PyObject* p_path = NULL;
  PyObject* p_content = NULL;
  if (!PyArg_ParseTuple(args, "O|O", &p_path, &p_content)) {
    return NULL;
  }
  const std::string path = py_convert_data<std::string>(p_path);
  if (p_content == NULL) {
    const int ret = qtouch(path);
    return py_convert(ret);
  } else {
    std::string content;
    py_convert(content, p_content);
    const int ret = qtouch(path, content);
    return py_convert(ret);
  }
})

EXPORT(qtouch_info, {
  using namespace qlat;
  PyObject* p_path = NULL;
  PyObject* p_content = NULL;
  if (!PyArg_ParseTuple(args, "O|O", &p_path, &p_content)) {
    return NULL;
  }
  const std::string path = py_convert_data<std::string>(p_path);
  if (p_content == NULL) {
    const int ret = qtouch_info(path);
    return py_convert(ret);
  } else {
    std::string content;
    py_convert(content, p_content);
    const int ret = qtouch_info(path, content);
    return py_convert(ret);
  }
})

EXPORT(qappend, {
  using namespace qlat;
  PyObject* p_path = NULL;
  PyObject* p_content = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_path, &p_content)) {
    return NULL;
  }
  const std::string path = py_convert_data<std::string>(p_path);
  std::string content;
  py_convert(content, p_content);
  const int ret = qappend(path, content);
  return py_convert(ret);
})

EXPORT(qappend_info, {
  using namespace qlat;
  PyObject* p_path = NULL;
  PyObject* p_content = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_path, &p_content)) {
    return NULL;
  }
  const std::string path = py_convert_data<std::string>(p_path);
  std::string content;
  py_convert(content, p_content);
  const int ret = qappend_info(path, content);
  return py_convert(ret);
})

EXPORT(qrename, {
  using namespace qlat;
  PyObject* p_old_path = NULL;
  PyObject* p_new_path = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_old_path, &p_new_path)) {
    return NULL;
  }
  const std::string old_path = py_convert_data<std::string>(p_old_path);
  const std::string new_path = py_convert_data<std::string>(p_new_path);
  const int ret = qrename(old_path, new_path);
  return py_convert(ret);
})

EXPORT(qrename_info, {
  using namespace qlat;
  PyObject* p_old_path = NULL;
  PyObject* p_new_path = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_old_path, &p_new_path)) {
    return NULL;
  }
  const std::string old_path = py_convert_data<std::string>(p_old_path);
  const std::string new_path = py_convert_data<std::string>(p_new_path);
  const int ret = qrename_info(old_path, new_path);
  return py_convert(ret);
})

EXPORT(qls, {
  using namespace qlat;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_path)) {
    return NULL;
  }
  const std::string path = py_convert_data<std::string>(p_path);
  const std::vector<std::string> ret = qls(path);
  return py_convert(ret);
})

EXPORT(qls_all, {
  using namespace qlat;
  PyObject* p_path = NULL;
  bool is_folder_before_files = false;
  if (!PyArg_ParseTuple(args, "O|b", &p_path, &is_folder_before_files)) {
    return NULL;
  }
  const std::string path = py_convert_data<std::string>(p_path);
  const std::vector<std::string> ret = qls_all(path, is_folder_before_files);
  return py_convert(ret);
})

EXPORT(compute_crc32, {
  using namespace qlat;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_path)) {
    return NULL;
  }
  const std::string path = py_convert_data<std::string>(p_path);
  const long ret = compute_crc32(path);
  return py_convert(ret);
})

EXPORT(check_all_files_crc32_info, {
  using namespace qlat;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_path)) {
    return NULL;
  }
  const std::string path = py_convert_data<std::string>(p_path);
  check_all_files_crc32_info(path);
  Py_RETURN_NONE;
})

EXPORT(qload_datatable, {
  using namespace qlat;
  PyObject* p_path = NULL;
  bool is_par = false;
  if (!PyArg_ParseTuple(args, "O|O", &p_path, &is_par)) {
    return NULL;
  }
  const std::string path = py_convert_data<std::string>(p_path);
  const DataTable dt = qload_datatable(path, is_par);
  return py_convert(dt);
})
