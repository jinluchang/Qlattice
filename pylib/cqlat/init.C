#include "exceptions.h"
#include "convert.h"

EXPORT(begin, {
  using namespace qlat;
  PyObject* p_v1 = NULL;
  PyObject* p_v2 = NULL;
  if (!PyArg_ParseTuple(args, "|OO", &p_v1, &p_v2)) {
    return NULL;
  }
  if (p_v1 != NULL and PyLong_Check(p_v1) and p_v2 != NULL) {
    // initialize with existing MPI
    int id_node = 0;
    Coordinate size_node;
    py_convert(id_node, p_v1);
    py_convert(size_node, p_v2);
    begin(id_node, size_node);
  } else {
    // initialize MPI by itself
    PyObject* p_sargs = p_v1;
    PyObject* p_node_size_list = p_v2;
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
  }
  //
  Py_RETURN_NONE;
});

EXPORT(end, {
  using namespace qlat;
  end();
  Py_RETURN_NONE;
});

EXPORT(get_id_node, {
  using namespace qlat;
  return py_convert(get_id_node());
});

EXPORT(get_num_node, {
  using namespace qlat;
  return py_convert(get_num_node());
});

EXPORT(get_coor_node, {
  using namespace qlat;
  return py_convert(get_coor_node());
});

EXPORT(get_size_node, {
  using namespace qlat;
  return py_convert(get_size_node());
});

