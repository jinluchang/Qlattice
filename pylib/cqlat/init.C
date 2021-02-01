#include "exceptions.h"
#include "convert.h"

EXPORT(begin, {
  using namespace qlat;
  PyObject* p_sargs = NULL;
  PyObject* p_node_size_list = NULL;
  if (!PyArg_ParseTuple(args, "|OO", &p_sargs, &p_node_size_list)) {
    return NULL;
  }
  std::vector<std::string> sargs;
  if (p_sargs != NULL) {
    py_convert(sargs, p_sargs);
  }
  std::vector<Coordinate> node_size_list;
  if (p_node_size_list != NULL) {
    py_convert(node_size_list, p_node_size_list);
  }
  // make cargs
  std::vector<const char*> cargs(sargs.size() + 1);
  for (long i = 0; i < (long)sargs.size(); ++i) {
    cargs[i] = sargs[i].c_str();
  }
  cargs.back() = NULL;
  //
  int argc = (int)sargs.size();
  char** argv = (char**)&cargs[0];
  //
  begin(&argc, &argv, node_size_list);
  //
  return PyLong_FromLong(0);
});

EXPORT(end, {
  using namespace qlat;
  end();
  return PyLong_FromLong(0);
});
