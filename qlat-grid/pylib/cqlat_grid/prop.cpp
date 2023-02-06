#include <qlat-grid/qlat-grid.h>

#include "lib.h"

EXPORT(save_prop_float, {
  using namespace qlat;
  TIMER_VERBOSE("save_prop_float");
  PyObject* p_prop = NULL;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_prop, &p_path)) {
    return NULL;
  }
  const Propagator4d& prop = py_convert_type<Propagator4d>(p_prop);
  const std::string path = py_convert_data<std::string>(p_path);
  const std::string path_dir = dirname(path);
  qmkdir_p_info(path_dir);
  const Geometry& geo = prop.geo();
  const Coordinate total_site = geo.total_site();
  const Coordinate size_node = geo.geon.size_node;
  Grid::GridCartesian UGrid(
      grid_convert(total_site),
      Grid::GridDefaultSimd(Grid::Nd, Grid::vComplexF::Nsimd()),
      grid_convert(size_node));
  Grid::LatticePropagatorF gprop(&UGrid);
  grid_convert(gprop, prop);
  Grid::emptyUserRecord record;
  Grid::ScidacWriter sw(UGrid.IsBoss());
  sw.open(path);
  sw.writeScidacFieldRecord(gprop, record);
  sw.close();
  Py_RETURN_NONE;
})

EXPORT(load_prop_float, {
  // p_prop need to have correct geometry
  using namespace qlat;
  TIMER_VERBOSE("load_prop_float");
  PyObject* p_prop = NULL;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_prop, &p_path)) {
    return NULL;
  }
  Propagator4d& prop = py_convert_type<Propagator4d>(p_prop);
  const std::string path = py_convert_data<std::string>(p_path);
  if (not does_file_exist_sync_node(path)) {
    // if file does not exist, clear the prop obj
    displayln_info(fname +
                   ssprintf(": file='%s' does not exist.", path.c_str()));
    prop.init();
    Py_RETURN_NONE;
  }
  const Geometry& geo = prop.geo();
  const Coordinate total_site = geo.total_site();
  const Coordinate size_node = geo.geon.size_node;
  Grid::GridCartesian UGrid(
      grid_convert(total_site),
      Grid::GridDefaultSimd(Grid::Nd, Grid::vComplexF::Nsimd()),
      grid_convert(size_node));
  Grid::LatticePropagatorF gprop(&UGrid);
  Grid::emptyUserRecord record;
  Grid::ScidacReader sr;
  sr.open(path);
  sr.readScidacFieldRecord(gprop, record);
  sr.close();
  grid_convert(prop, gprop);
  Py_RETURN_NONE;
})
