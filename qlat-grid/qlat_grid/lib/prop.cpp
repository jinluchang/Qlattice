#include <qlat-grid/qlat-grid.h>

namespace qlat
{  //

void save_grid_prop_float(const Field<WilsonMatrix>& prop, const std::string& path)
{
  TIMER_VERBOSE("save_grid_prop_float");
  const std::string path_dir = dirname(path);
  qmkdir_p_info(path_dir);
  const Geometry& geo = prop.geo();
  qassert(geo.multiplicity == 1);
  const Coordinate total_site = geo.total_site();
  const Coordinate size_node = geo.geon.size_node;
  Grid::GridCartesian UGrid(
      grid_convert(total_site),
      Grid::GridDefaultSimd(Grid::Nd, Grid::vComplexF::Nsimd()),
      grid_convert(size_node));
  Grid::LatticePropagatorF gprop(&UGrid);
  grid_convert(gprop, (const Propagator4d&)prop);
  Grid::emptyUserRecord record;
  Grid::ScidacWriter sw(UGrid.IsBoss());
  sw.open(path + ".partial");
  sw.writeScidacFieldRecord(gprop, record);
  sw.close();
  qrename(path + ".partial", path);
}

void load_grid_prop_float(Field<WilsonMatrix>& prop, const std::string& path)
{
  TIMER_VERBOSE("load_grid_prop_float");
  // p_prop need to have correct geometry
  if (not does_file_exist_sync_node(path)) {
    // if file does not exist, clear the prop obj
    displayln_info(fname +
                   ssprintf(": file='%s' does not exist.", path.c_str()));
    prop.init();
    return;
  }
  const Geometry& geo = prop.geo();
  qassert(geo.multiplicity == 1);
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
  grid_convert((Propagator4d&)prop, gprop);
}

void save_grid_prop_double(const Field<WilsonMatrix>& prop, const std::string& path)
{
  TIMER_VERBOSE("save_grid_prop_double");
  const std::string path_dir = dirname(path);
  qmkdir_p_info(path_dir);
  const Geometry& geo = prop.geo();
  qassert(geo.multiplicity == 1);
  const Coordinate total_site = geo.total_site();
  const Coordinate size_node = geo.geon.size_node;
  Grid::GridCartesian UGrid(
      grid_convert(total_site),
      Grid::GridDefaultSimd(Grid::Nd, Grid::vComplexD::Nsimd()),
      grid_convert(size_node));
  Grid::LatticePropagatorD gprop(&UGrid);
  grid_convert(gprop, (const Propagator4d&)prop);
  Grid::emptyUserRecord record;
  Grid::ScidacWriter sw(UGrid.IsBoss());
  sw.open(path + ".partial");
  sw.writeScidacFieldRecord(gprop, record);
  sw.close();
  qrename(path + ".partial", path);
}

void load_grid_prop_double(Field<WilsonMatrix>& prop, const std::string& path)
{
  TIMER_VERBOSE("load_grid_prop_double");
  // p_prop need to have correct geometry
  if (not does_file_exist_sync_node(path)) {
    // if file does not exist, clear the prop obj
    displayln_info(fname +
                   ssprintf(": file='%s' does not exist.", path.c_str()));
    prop.init();
    return;
  }
  const Geometry& geo = prop.geo();
  qassert(geo.multiplicity == 1);
  const Coordinate total_site = geo.total_site();
  const Coordinate size_node = geo.geon.size_node;
  Grid::GridCartesian UGrid(
      grid_convert(total_site),
      Grid::GridDefaultSimd(Grid::Nd, Grid::vComplexD::Nsimd()),
      grid_convert(size_node));
  Grid::LatticePropagatorD gprop(&UGrid);
  Grid::emptyUserRecord record;
  Grid::ScidacReader sr;
  sr.open(path);
  sr.readScidacFieldRecord(gprop, record);
  sr.close();
  grid_convert((Propagator4d&)prop, gprop);
}

}  // namespace qlat
