#include "lib.h"

#include "compute-three-point-func.h"

EXPORT(compute_three_point_func, {
  using namespace qlat;
  PyObject* p_job_tag = NULL;
  int traj = 0;
  if (!PyArg_ParseTuple(args, "Oi", &p_job_tag, &traj)) {
    return NULL;
  }
  std::string job_tag;
  py_convert(job_tag, p_job_tag);
  compute_three_point_func(job_tag, traj);
  Py_RETURN_NONE;
});

