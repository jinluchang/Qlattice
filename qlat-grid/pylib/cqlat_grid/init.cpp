#include <qlat-grid/qlat-grid.h>

#include "lib.h"

EXPORT(begin_with_grid, {
  using namespace qlat;
  PyObject* p_sargs = NULL;
  PyObject* p_size_node_list = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_sargs, &p_size_node_list)) {
    return NULL;
  }
  static std::vector<std::string> sargs =
      py_convert_data<std::vector<std::string> >(p_sargs);
  const std::vector<Coordinate> node_size_list =
      py_convert_data<std::vector<Coordinate> >(p_size_node_list);
  // make cargs
  static std::vector<const char*> cargs(sargs.size() + 1);
  for (long i = 0; i < (long)sargs.size(); ++i) {
    cargs[i] = sargs[i].c_str();
  }
  cargs.back() = NULL;
  //
  int argc = (int)sargs.size();
  char** argv = (char**)&cargs[0];
  //
  grid_begin(&argc, &argv, node_size_list);
  Py_RETURN_NONE;
})
