#include "lib.h"

#include "compute-meson-vv.h"

EXPORT(compute_meson_vv, {
  using namespace qlat;
  PyObject* p_job_tag = NULL;
  int traj = 0;
  if (!PyArg_ParseTuple(args, "Oi", &p_job_tag, &traj)) {
    return NULL;
  }
  std::string job_tag;
  py_convert(job_tag, p_job_tag);
  compute_meson_vv(job_tag, traj);
  Py_RETURN_NONE;
});

