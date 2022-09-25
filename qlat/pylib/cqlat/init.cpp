#include "lib.h"

EXPORT(begin, {
  // id_node, size_node, color = 0
  // sys.argv, node_size_list = []
  using namespace qlat;
  PyObject* p_v1 = NULL;
  PyObject* p_v2 = NULL;
  PyObject* p_v3 = NULL;
  if (!PyArg_ParseTuple(args, "|OOi", &p_v1, &p_v2, &p_v3)) {
    return NULL;
  }
  if (p_v1 != NULL and PyLong_Check(p_v1) and p_v2 != NULL) {
    // initialize with existing MPI
    int id_node = 0;
    Coordinate size_node;
    int color = 0;
    py_convert(id_node, p_v1);
    py_convert(size_node, p_v2);
    if (p_v3 != NULL) {
      py_convert(color, p_v3);
    }
    begin(id_node, size_node, color);
  } else {
    // initialize MPI by itself
    PyObject* p_sargs = p_v1;
    PyObject* p_node_size_list = p_v2;
    static std::vector<std::string> sargs;
    if (p_sargs != NULL) {
      py_convert(sargs, p_sargs);
    }
    std::vector<Coordinate> node_size_list;
    if (p_node_size_list != NULL) {
      py_convert(node_size_list, p_node_size_list);
    }
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
    begin(&argc, &argv, node_size_list);
  }
  //
  qset_line_buf(get_output_file());
  //
  Py_RETURN_NONE;
})

EXPORT(end, {
  using namespace qlat;
  bool is_preserving_cache = false;
  if (!PyArg_ParseTuple(args, "|b", &is_preserving_cache)) {
    return NULL;
  }
  end(is_preserving_cache);
  Py_RETURN_NONE;
})
