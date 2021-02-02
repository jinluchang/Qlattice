#include "exceptions.h"
#include "convert.h"

EXPORT(mk_geo, {
  using namespace qlat;
  PyObject* p_total_site = NULL;
  int multiplicity = 1;
  if (!PyArg_ParseTuple(args, "O|i", &p_total_site, &multiplicity)) {
    return NULL;
  }
  Coordinate total_site;
  py_convert(total_site, p_total_site);
  Geometry* pgeo = new Geometry(total_site, multiplicity);
  return py_convert((void*)pgeo);
});

EXPORT(free_geo, {
  using namespace qlat;
  Geometry* pgeo = NULL;
  if (!PyArg_ParseTuple(args, "l", &pgeo)) {
    return NULL;
  }
  pqassert(pgeo != NULL);
  delete pgeo;
  Py_RETURN_NONE;
});

EXPORT(geo_reform, {
  using namespace qlat;
  Geometry* pgeo = NULL;
  int multiplicity = 1;
  PyObject* p_expansion_left = NULL;
  PyObject* p_expansion_right = NULL;
  if (!PyArg_ParseTuple(args, "l|iOO", &pgeo, &multiplicity, &p_expansion_left, &p_expansion_right)) {
    return NULL;
  }
  pqassert(pgeo != NULL);
  const Geometry& geo = *pgeo;
  Coordinate expansion_left, expansion_right;
  if (NULL != p_expansion_left) {
    py_convert(expansion_left, p_expansion_left);
  }
  if (NULL != p_expansion_right) {
    py_convert(expansion_right, p_expansion_right);
  }
  Geometry* pgeo_new = new Geometry();
  *pgeo_new = geo_resize(geo, expansion_left, expansion_right);
  return py_convert((void*)pgeo_new);
});

EXPORT(geo_eo, {
  using namespace qlat;
  Geometry* pgeo = NULL;
  int eo = 0;
  if (!PyArg_ParseTuple(args, "l|i", &pgeo, &eo)) {
    return NULL;
  }
  pqassert(pgeo != NULL);
  const Geometry& geo = *pgeo;
  Geometry* pgeo_new = new Geometry();
  *pgeo_new = geo_eo(geo, eo);
  return py_convert((void*)pgeo_new);
});

EXPORT(get_total_site_geo, {
  using namespace qlat;
  Geometry* pgeo = NULL;
  if (!PyArg_ParseTuple(args, "l", &pgeo)) {
    return NULL;
  }
  pqassert(pgeo != NULL);
  const Geometry& geo = *pgeo;
  return py_convert(geo.total_site());
});

EXPORT(get_multiplicity_geo, {
  using namespace qlat;
  Geometry* pgeo = NULL;
  if (!PyArg_ParseTuple(args, "l", &pgeo)) {
    return NULL;
  }
  pqassert(pgeo != NULL);
  const Geometry& geo = *pgeo;
  return py_convert(geo.multiplicity);
});

EXPORT(get_node_site_geo, {
  using namespace qlat;
  Geometry* pgeo = NULL;
  if (!PyArg_ParseTuple(args, "l", &pgeo)) {
    return NULL;
  }
  pqassert(pgeo != NULL);
  const Geometry& geo = *pgeo;
  return py_convert(geo.node_site);
});

EXPORT(get_eo_geo, {
  using namespace qlat;
  Geometry* pgeo = NULL;
  if (!PyArg_ParseTuple(args, "l", &pgeo)) {
    return NULL;
  }
  pqassert(pgeo != NULL);
  const Geometry& geo = *pgeo;
  return py_convert(geo.eo);
});

EXPORT(get_expansion_left_geo, {
  using namespace qlat;
  Geometry* pgeo = NULL;
  if (!PyArg_ParseTuple(args, "l", &pgeo)) {
    return NULL;
  }
  pqassert(pgeo != NULL);
  const Geometry& geo = *pgeo;
  return py_convert(geo.expansion_left);
});

EXPORT(get_expansion_right_geo, {
  using namespace qlat;
  Geometry* pgeo = NULL;
  if (!PyArg_ParseTuple(args, "l", &pgeo)) {
    return NULL;
  }
  pqassert(pgeo != NULL);
  const Geometry& geo = *pgeo;
  return py_convert(geo.expansion_right);
});

EXPORT(get_id_node_geo, {
  using namespace qlat;
  Geometry* pgeo = NULL;
  if (!PyArg_ParseTuple(args, "l", &pgeo)) {
    return NULL;
  }
  pqassert(pgeo != NULL);
  const Geometry& geo = *pgeo;
  return py_convert(geo.geon.id_node);
});

EXPORT(get_num_node_geo, {
  using namespace qlat;
  Geometry* pgeo = NULL;
  if (!PyArg_ParseTuple(args, "l", &pgeo)) {
    return NULL;
  }
  pqassert(pgeo != NULL);
  const Geometry& geo = *pgeo;
  return py_convert(geo.geon.num_node);
});

EXPORT(get_coor_node_geo, {
  using namespace qlat;
  Geometry* pgeo = NULL;
  if (!PyArg_ParseTuple(args, "l", &pgeo)) {
    return NULL;
  }
  pqassert(pgeo != NULL);
  const Geometry& geo = *pgeo;
  return py_convert(geo.geon.coor_node);
});

EXPORT(get_size_node_geo, {
  using namespace qlat;
  Geometry* pgeo = NULL;
  if (!PyArg_ParseTuple(args, "l", &pgeo)) {
    return NULL;
  }
  pqassert(pgeo != NULL);
  const Geometry& geo = *pgeo;
  return py_convert(geo.geon.size_node);
});
