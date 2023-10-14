#include <qlat-cps/qlat-cps.h>
#include <util/qio_readPropagator.h>
#include <util/qio_writePropagator.h>

namespace qlat
{  //

void save_cps_prop_double(const Field<WilsonMatrix>& prop,
                          const std::string& path)
{
  TIMER_VERBOSE("save_cps_prop_double");
  const std::string path_dir = dirname(path);
  qmkdir_p_info(path_dir);
  const Geometry& geo = prop.geo();
  qassert(geo.multiplicity == 1);
  Field<WilsonMatrix> prop_src;
  prop_src.init(geo);
  Vector<WilsonMatrix> prop_vec = get_data(prop);
  Vector<WilsonMatrix> prop_src_vec = get_data(prop_src);
  const std::string path_partial = path + ".partial";
  char* outfile = (char*)path_partial.c_str();
  const cps::QIO_PROP_SOURCE_TYPES sType = cps::QIO_FULL_SOURCE;
  void* prop_data = prop_vec.data();
  void* prop_src_data = prop_src_vec.data();
  const int volFormat = QIO_SINGLEFILE;
  const cps::FP_FORMAT floatFormat = cps::FP_AUTOMATIC;
  cps::qio_writePropagator(outfile, sType, prop_data, prop_src_data, volFormat,
                           floatFormat);
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
  Field<WilsonMatrix> prop_src;
  prop_src.init(geo);
  Vector<WilsonMatrix> prop_vec = get_data(prop);
  Vector<WilsonMatrix> prop_src_vec = get_data(prop_src);
  char* infile = (char*)path.c_str();
  const cps::QIO_PROP_SOURCE_TYPES sType = cps::QIO_FULL_SOURCE;
  void* prop_data = prop_vec.data();
  void* prop_src_data = prop_src_vec.data();
  const int volFormat = QIO_UNKNOWN;
  cps::qio_readPropagator(infile, sType, prop_data, prop_src_data, volFormat);
}

}  // namespace qlat
