#include <qlat-cps/qlat-cps.h>

namespace qlat
{  //

void save_cps_prop_float(const Field<WilsonMatrix>& prop, const std::string& path)
{
  TIMER_VERBOSE("save_cps_prop_float");
  const std::string path_dir = dirname(path);
  qmkdir_p_info(path_dir);
  const Geometry& geo = prop.geo();
  qassert(geo.multiplicity == 1);
  const Coordinate total_site = geo.total_site();
  const Coordinate size_node = geo.geon.size_node;
  write_field(prop, path);
  qrename(path + ".partial", path);
}

void load_cps_prop_float(Field<WilsonMatrix>& prop, const std::string& path)
{
  TIMER_VERBOSE("load_cps_prop_float");
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
  read_field(prop, path);
}

void save_cps_prop_double(const Field<WilsonMatrix>& prop, const std::string& path)
{
  TIMER_VERBOSE("save_cps_prop_double");
  const std::string path_dir = dirname(path);
  qmkdir_p_info(path_dir);
  const Geometry& geo = prop.geo();
  qassert(geo.multiplicity == 1);
  const Coordinate total_site = geo.total_site();
  const Coordinate size_node = geo.geon.size_node;
  write_field(prop, path);
  qrename(path + ".partial", path);
}

void load_cps_prop_double(Field<WilsonMatrix>& prop, const std::string& path)
{
  TIMER_VERBOSE("load_cps_prop_double");
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
  read_field(prop, path);
}

}  // namespace qlat
