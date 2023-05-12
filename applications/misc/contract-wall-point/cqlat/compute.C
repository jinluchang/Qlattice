#include "lib.h"

#include "qlat/qlat-setup.h"
#include "data-load.h"

EXPORT(setup, {
  using namespace qlat;
  PyObject* p_job_tag = NULL;
  if (!PyArg_ParseTuple(args, "|O", &p_job_tag)) {
    return NULL;
  }
  if (p_job_tag == NULL) {
    setup();
  } else {
    std::string job_tag;
    py_convert(job_tag, p_job_tag);
    setup(job_tag);
  }
  Py_RETURN_NONE;
});

EXPORT(setup_log_idx, {
  using namespace qlat;
  setup_log_idx();
  Py_RETURN_NONE;
});

EXPORT(clear_all_data_cache, {
  using namespace qlat;
  clear_all_data_cache();
  Py_RETURN_NONE;
});

EXPORT(set_data_path, {
  using namespace qlat;
  PyObject* p_data_path = NULL;
  if (!PyArg_ParseTuple(args, "|O", &p_data_path)) {
    return NULL;
  }
  if (p_data_path != NULL) {
    py_convert(get_data_path(), p_data_path);
  }
  return py_convert(get_data_path());
});
