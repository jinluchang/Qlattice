#include "lib.h"

EXPORT(mk_geo, {
  using namespace qlat;
  Geometry* pgeo = new Geometry();
  return py_convert((void*)pgeo);
})

EXPORT(free_geo, {
  using namespace qlat;
  return free_obj<Geometry>(args);
})

EXPORT(set_geo, {
  using namespace qlat;
  return set_obj<Geometry>(args);
})

EXPORT(set_geo_total_site, {
  using namespace qlat;
  PyObject* p_geo = NULL;
  PyObject* p_total_site = NULL;
  int multiplicity = 1;
  if (!PyArg_ParseTuple(args, "OO|i", &p_geo, &p_total_site, &multiplicity)) {
    return NULL;
  }
  Geometry& geo = py_convert_type<Geometry>(p_geo);
  const Coordinate total_site = py_convert_data<Coordinate>(p_total_site);
  geo.init(total_site, multiplicity);
  Py_RETURN_NONE;
})

EXPORT(set_geo_reform, {
  using namespace qlat;
  PyObject* p_geo = NULL;
  int multiplicity = 1;
  PyObject* p_expansion_left = NULL;
  PyObject* p_expansion_right = NULL;
  if (!PyArg_ParseTuple(args, "O|iOO", &p_geo, &multiplicity, &p_expansion_left,
                        &p_expansion_right)) {
    return NULL;
  }
  Geometry& geo = py_convert_type<Geometry>(p_geo);
  Coordinate expansion_left, expansion_right;
  if (NULL != p_expansion_left) {
    py_convert(expansion_left, p_expansion_left);
  }
  if (NULL != p_expansion_right) {
    py_convert(expansion_right, p_expansion_right);
  }
  geo.remult(multiplicity);
  geo.resize(expansion_left, expansion_right);
  Py_RETURN_NONE;
})

EXPORT(set_geo_eo, {
  using namespace qlat;
  PyObject* p_geo = NULL;
  int eo = 0;
  if (!PyArg_ParseTuple(args, "O|i", &p_geo, &eo)) {
    return NULL;
  }
  Geometry& geo = py_convert_type<Geometry>(p_geo);
  geo.eo = eo;
  Py_RETURN_NONE;
})

EXPORT(get_total_site_geo, {
  using namespace qlat;
  PyObject* p_geo = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_geo)) {
    return NULL;
  }
  const Geometry& geo = py_convert_type<Geometry>(p_geo);
  return py_convert(geo.total_site());
})

EXPORT(get_local_volume_geo, {
  using namespace qlat;
  PyObject* p_geo = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_geo)) {
    return NULL;
  }
  const Geometry& geo = py_convert_type<Geometry>(p_geo);
  return py_convert(geo.local_volume());
})

EXPORT(get_total_volume_geo, {
  using namespace qlat;
  PyObject* p_geo = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_geo)) {
    return NULL;
  }
  const Geometry& geo = py_convert_type<Geometry>(p_geo);
  return py_convert(geo.total_volume());
})

EXPORT(get_multiplicity_geo, {
  using namespace qlat;
  PyObject* p_geo = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_geo)) {
    return NULL;
  }
  const Geometry& geo = py_convert_type<Geometry>(p_geo);
  return py_convert(geo.multiplicity);
})

EXPORT(get_node_site_geo, {
  using namespace qlat;
  PyObject* p_geo = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_geo)) {
    return NULL;
  }
  const Geometry& geo = py_convert_type<Geometry>(p_geo);
  return py_convert(geo.node_site);
})

EXPORT(get_eo_geo, {
  using namespace qlat;
  PyObject* p_geo = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_geo)) {
    return NULL;
  }
  const Geometry& geo = py_convert_type<Geometry>(p_geo);
  return py_convert(geo.eo);
})

EXPORT(get_expansion_left_geo, {
  using namespace qlat;
  PyObject* p_geo = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_geo)) {
    return NULL;
  }
  const Geometry& geo = py_convert_type<Geometry>(p_geo);
  return py_convert(geo.expansion_left);
})

EXPORT(get_expansion_right_geo, {
  using namespace qlat;
  PyObject* p_geo = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_geo)) {
    return NULL;
  }
  const Geometry& geo = py_convert_type<Geometry>(p_geo);
  return py_convert(geo.expansion_right);
})

EXPORT(get_id_node_geo, {
  using namespace qlat;
  PyObject* p_geo = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_geo)) {
    return NULL;
  }
  const Geometry& geo = py_convert_type<Geometry>(p_geo);
  return py_convert(geo.geon.id_node);
})

EXPORT(get_num_node_geo, {
  using namespace qlat;
  PyObject* p_geo = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_geo)) {
    return NULL;
  }
  const Geometry& geo = py_convert_type<Geometry>(p_geo);
  return py_convert(geo.geon.num_node);
})

EXPORT(get_coor_node_geo, {
  using namespace qlat;
  PyObject* p_geo = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_geo)) {
    return NULL;
  }
  const Geometry& geo = py_convert_type<Geometry>(p_geo);
  return py_convert(geo.geon.coor_node);
})

EXPORT(get_size_node_geo, {
  using namespace qlat;
  PyObject* p_geo = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_geo)) {
    return NULL;
  }
  const Geometry& geo = py_convert_type<Geometry>(p_geo);
  return py_convert(geo.geon.size_node);
})

EXPORT(coordinate_g_from_l_geo, {
  using namespace qlat;
  PyObject* p_geo = NULL;
  PyObject* p_xl = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_geo, &p_xl)) {
    return NULL;
  }
  const Geometry& geo = py_convert_type<Geometry>(p_geo);
  Coordinate xl;
  py_convert(xl, p_xl);
  return py_convert(geo.coordinate_g_from_l(xl));
})

EXPORT(coordinate_l_from_g_geo, {
  using namespace qlat;
  PyObject* p_geo = NULL;
  PyObject* p_xg = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_geo, &p_xg)) {
    return NULL;
  }
  const Geometry& geo = py_convert_type<Geometry>(p_geo);
  Coordinate xg;
  py_convert(xg, p_xg);
  return py_convert(geo.coordinate_l_from_g(xg));
})

EXPORT(is_local_geo, {
  using namespace qlat;
  PyObject* p_geo = NULL;
  PyObject* p_xl = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_geo, &p_xl)) {
    return NULL;
  }
  const Geometry& geo = py_convert_type<Geometry>(p_geo);
  const Coordinate xl = py_convert_data<Coordinate>(p_xl);
  return py_convert(geo.is_local(xl));
})

EXPORT(is_local_xg_geo, {
  using namespace qlat;
  PyObject* p_geo = NULL;
  PyObject* p_xg = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_geo, &p_xg)) {
    return NULL;
  }
  const Geometry& geo = py_convert_type<Geometry>(p_geo);
  const Coordinate xg = py_convert_data<Coordinate>(p_xg);
  const Coordinate xl = geo.coordinate_l_from_g(xg);
  return py_convert(geo.is_local(xl));
})

EXPORT(get_xg_list, {
  // return xg for all local sites
  using namespace qlat;
  PyObject* p_geo = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_geo)) {
    return NULL;
  }
  const Geometry& geo = py_convert_type<Geometry>(p_geo);
  vector<Coordinate> xgs(geo.local_volume());
  qthread_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    xgs[index] = xg;
  })
  return py_convert(get_data(xgs));
})
