#define QLAT_INSTANTIATE_FIELD_IO

#include <qlat/field-io.h>

namespace qlat
{  //

bool is_dist_field(const std::string& path)
{
  TIMER("is_dist_field");
  Long nfile = 0;
  if (get_id_node() == 0) {
    if (does_file_exist(path + "/geo-info.txt") and
        does_file_exist(path + "/checkpoint")) {
      nfile = 1;
    }
  }
  bcast(get_data(nfile));
  return nfile > 0;
}

bool is_field(const std::string& path)
{
  TIMER("is_field");
  Long nfile = 0;
  if (get_id_node() == 0) {
    QFile qfile = qfopen(path, "r");
    if (not qfile.null()) {
      const std::string header = "BEGIN_FIELD_HEADER\n";
      std::vector<char> check_line(header.size(), 0);
      if (1 == qfread(check_line.data(), header.size(), 1, qfile)) {
        if (std::string(check_line.data(), check_line.size()) == header) {
          nfile = 1;
        }
      }
    }
    qfclose(qfile);
  }
  bcast(get_data(nfile));
  return nfile > 0;
}

bool is_d_field(const std::string& path)
{
  return is_dist_field(path) or is_field(path);
}

bool dist_repartition(const Coordinate& new_size_node, const std::string& path,
                      const std::string& new_path)
// interface_function
{
  bool is_failed = false;
  const std::string npath = remove_trailing_slashes(path);
  if (std::string(npath, npath.length() - 4, 4) == ".tmp") {
    return true;
  }
  if (does_file_exist_sync_node(npath + "-lock")) {
    return true;
  }
  const std::string new_npath = remove_trailing_slashes(new_path);
  const bool is_dir = is_directory_sync_node(path);
  if (is_dir and new_size_node != Coordinate(1, 1, 1, 1) and
      (new_npath == npath or new_npath == "")) {
    if (is_dist_field(npath)) {
      Geometry geo;
      Int multiplicity;
      Int sizeof_M;
      Coordinate size_node;
      dist_read_geo_info(geo, multiplicity, sizeof_M, size_node, npath);
      if (size_node == new_size_node and
          (new_path == "" or new_npath == npath)) {
        displayln_info(
            ssprintf("repartition: size_node=%s ; no need to repartition '%s'.",
                     show(size_node).c_str(), npath.c_str()));
        return true;
      }
    } else {
      displayln_info(
          ssprintf("repartition: WARNING: not a folder to partition: '%s'.",
                   npath.c_str()));
    }
  }
  if (not is_dir and new_size_node == Coordinate(1, 1, 1, 1) and
      (new_npath == npath or new_npath == "")) {
    displayln_info(
        ssprintf("repartition: size_node=%s ; no need to repartition '%s'.",
                 show(new_size_node).c_str(), npath.c_str()));
    return true;
  }
  if (not obtain_lock(npath + "-lock")) {
    return true;
  }
  TIMER_VERBOSE("dist_repartition");
  Field<RealF> f;
  read_field(f, npath);
  if (get_incorrect_field_read_sizeof_M() != 0) {
    get_force_field_write_sizeof_M() = get_incorrect_field_read_sizeof_M();
    get_incorrect_field_read_sizeof_M() = 0;
  }
  if (new_npath == npath or new_npath == "") {
    qassert(not does_file_exist_sync_node(npath + "-repartition-new.tmp"));
    qassert(not does_file_exist_sync_node(npath + "-repartition-old.tmp"));
    if (new_size_node == Coordinate(1, 1, 1, 1)) {
      write_field(f, npath + "-repartition-new.tmp");
    } else {
      dist_write_field(f, new_size_node, npath + "-repartition-new.tmp");
    }
    qrename_info(npath, npath + "-repartition-old.tmp");
    qrename_info(npath + "-repartition-new.tmp", npath);
    qremove_all_info(npath + "-repartition-old.tmp");
  } else {
    if (new_size_node == Coordinate(1, 1, 1, 1)) {
      write_field(f, new_npath);
    } else {
      dist_write_field(f, new_size_node, new_npath);
    }
  }
  release_lock();
  return is_failed;
}

}  // namespace qlat
