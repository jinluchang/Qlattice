#include "lib.h"

#include "compute-chvp.h"

EXPORT(compute_chvp, {
  using namespace qlat;
  PyObject* p_job_tag = NULL;
  int traj = 0;
  if (!PyArg_ParseTuple(args, "Oi", &p_job_tag, &traj)) {
    return NULL;
  }
  std::string job_tag;
  py_convert(job_tag, p_job_tag);
  compute_chvp(job_tag, traj);
  Py_RETURN_NONE;
});

