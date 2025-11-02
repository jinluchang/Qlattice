// vim: set ts=2 sw=2 expandtab:

#pragma once

#include <qlat-utils/qar.h>
#include <qlat/field-shuffle.h>
#include <qlat/utils-io.h>

namespace qlat
{  //

const Int DATA_READ_WRITE_NUMBER_OF_DIRECTORIES = 32;

API inline Int& dist_write_par_limit()
// qlat parameter
{
  static Int npar = get_env_long_default("q_write_par_limit", 16);
  return npar;
}

API inline Int& dist_read_par_limit()
// qlat parameter
{
  static Int npar = get_env_long_default("q_read_par_limit", 16);
  return npar;
}

API inline Int& get_force_field_write_sizeof_M()
{
  static Int sizeof_M = 0;
  return sizeof_M;
}

API inline Int& get_incorrect_field_read_sizeof_M()
{
  static Int sizeof_M = 0;
  return sizeof_M;
}

API inline bool& is_checksum_mismatch()
// qlat parameter
// if initially false, then will be set to true when checksum check fail
// if initially true, then the program will crash when checksum fail
// by default, the program will crash if checksums fail
{
  static bool b = true;
  return b;
}

inline Int dist_mkdir(const std::string& path, const Int num_node,
                      const mode_t mode = default_dir_mode())
{
  Int ret = 0;
  qmkdir_p(path, mode);
  for (Int i = 0; i < std::min(num_node, DATA_READ_WRITE_NUMBER_OF_DIRECTORIES);
       ++i) {
    const std::string dir = path + ssprintf("/%02d", i);
    ret += qmkdir(dir, mode);
  }
  return ret;
}

inline Int compute_dist_file_dir_id(const Int id_node, const Int num_node)
{
  const Int ndir = std::min(num_node, DATA_READ_WRITE_NUMBER_OF_DIRECTORIES);
  const Int dir_size = (num_node - 1) / ndir + 1;
  const Int id_dir = id_node / dir_size;
  qassert(0 <= id_dir && id_dir < ndir);
  return id_dir;
}

inline std::string dist_file_name(const std::string& path, const Int id_node,
                                  const Int num_node)
{
  return path + ssprintf("/%02d/%010d",
                         compute_dist_file_dir_id(id_node, num_node), id_node);
}

inline QFile dist_open(const std::string& path, const Int id_node,
                       const Int num_node, const std::string& fmode)
// fmode: "w" for write, "r" for read, and "a" for append
{
  const std::string fn = dist_file_name(path, id_node, num_node);
  QFile ret = qfopen(fn, fmode);
  if (ret.null()) {
    qwarn("dist_open: " + ssprintf("failed to open '%s'.", fn.c_str()));
  }
  return ret;
}

inline void dist_close(QFile& fp) { qfclose(fp); }

template <class M>
struct DistData {
  Int id_node;
  Vector<M> data;
};

template <class M>
Vector<M> get_data(const DistData<M>& dd)
{
  return dd.data;
}

template <class M>
std::vector<crc32_t> dist_crc32s(const std::vector<DistData<M>>& dds,
                                 const Int num_node)
{
  SYNC_NODE();
  TIMER_VERBOSE_FLOPS("dist_crc32s");
  Long total_bytes = 0;
  std::vector<crc32_t> ret(num_node, 0);
  for (Int k = 0; k < (int)dds.size(); ++k) {
    const DistData<M>& dd = dds[k];
    const Int id_node = dd.id_node;
    ret[id_node] = crc32_par(ret[id_node], dd.data);
    total_bytes += dd.data.data_size();
  }
  glb_sum(get_data_char(ret));
  glb_sum(total_bytes);
  timer.flops += total_bytes;
  return ret;
}

inline crc32_t dist_crc32(const std::vector<crc32_t>& crcs)
{
  std::vector<crc32_t> d = crcs;
  to_from_big_endian(get_data(d));
  return crc32(get_data(d));
}

template <class M>
crc32_t dist_crc32(const std::vector<DistData<M>>& dds, const Int num_node)
{
  return dist_crc32(dist_crc32s(dds, num_node));
}

template <class M>
crc32_t field_dist_crc32(const Field<M>& f)
{
  std::vector<DistData<M>> dds(1);
  dds[0].id_node = f.geo().geon.id_node;
  dds[0].data = get_data(f);
  return dist_crc32(dds, get_num_node());
}

inline void dist_write_geo_info(const Geometry& geo, const Int multiplicity,
                                const Int sizeof_M, const std::string& path,
                                const mode_t mode = default_dir_mode())
{
  (void)mode;
  TIMER("dist_write_geo_info");
  const Int id_node = geo.geon.id_node;
  qassert(geo.is_only_local);
  if (0 == id_node) {
    const std::string fn = path + "/geo-info.txt";
    QFile fp = qfopen(fn, "w");
    qwrite_data(ssprintf("node_file_size = %ld\n",
                         sizeof_M * multiplicity * geo.local_volume()),
                fp);
    qwrite_data(ssprintf("multiplicity = %d\n", multiplicity), fp);
    qwrite_data(ssprintf("sizeof(M) = %d\n", sizeof_M), fp);
    qwrite_data(ssprintf("geo.geon.num_node = %d\n", geo.geon.num_node), fp);
    qwrite_data(ssprintf("geo.geon.size_node[0] = %d\n", geo.geon.size_node[0]),
                fp);
    qwrite_data(ssprintf("geo.geon.size_node[1] = %d\n", geo.geon.size_node[1]),
                fp);
    qwrite_data(ssprintf("geo.geon.size_node[2] = %d\n", geo.geon.size_node[2]),
                fp);
    qwrite_data(ssprintf("geo.geon.size_node[3] = %d\n", geo.geon.size_node[3]),
                fp);
    qwrite_data(ssprintf("geo.local_volume() = %ld\n", geo.local_volume()), fp);
    qwrite_data(ssprintf("geo.node_site[0] = %d\n", geo.node_site[0]), fp);
    qwrite_data(ssprintf("geo.node_site[1] = %d\n", geo.node_site[1]), fp);
    qwrite_data(ssprintf("geo.node_site[2] = %d\n", geo.node_site[2]), fp);
    qwrite_data(ssprintf("geo.node_site[3] = %d\n", geo.node_site[3]), fp);
    qwrite_data(ssprintf("PI = %.20f\n", PI), fp);
    const char* pic = (const char*)&PI;
    qwrite_data(
        ssprintf("PI_double = %hhx %hhx %hhx %hhx %hhx %hhx %hhx %hhx\n",
                 pic[0], pic[1], pic[2], pic[3], pic[4], pic[5], pic[6],
                 pic[7]),
        fp);
    const RealF PIf = PI;
    const char* pifc = (const char*)&PIf;
    qwrite_data(ssprintf("PI_float = %hhx %hhx %hhx %hhx\n", pifc[0], pifc[1],
                         pifc[2], pifc[3]),
                fp);
    qfclose(fp);
  }
}

inline void dist_read_geo_info(Geometry& geo, Int& multiplicity, Int& sizeof_M,
                               Coordinate& new_size_node,
                               const std::string& path)
{
  TIMER("dist_read_geo_info");
  Coordinate size_node;
  Coordinate node_site;
  if (get_id_node() == 0) {
    const std::string fn = path + "/geo-info.txt";
    const std::vector<std::string> lines = qgetlines(fn);
    reads(multiplicity,
          info_get_prop(lines, "multiplicity = ", "geo.multiplicity = "));
    reads(sizeof_M, info_get_prop(lines, "sizeof(M) = "));
    for (Int i = 0; i < 4; ++i) {
      reads(size_node[i],
            info_get_prop(lines, ssprintf("geo.geon.size_node[%d] = ", i),
                          ssprintf("geo.sizeNode[%d] = ", i)));
      reads(node_site[i],
            info_get_prop(lines, ssprintf("geo.node_site[%d] = ", i),
                          ssprintf("geo.nodeSite[%d] = ", i)));
    }
    Long node_file_size;
    Int num_node;
    Long local_volume;
    reads(node_file_size,
          info_get_prop(lines, "node_file_size = ", "nodeFileSize = "));
    reads(num_node,
          info_get_prop(lines, "geo.geon.num_node = ", "geo.numNode = "));
    reads(
        local_volume,
        info_get_prop(lines, "geo.local_volume() = ", "geo.localVolume() = "));
    qassert(num_node == product(size_node));
    qassert(local_volume == product(node_site));
    qassert(node_file_size == local_volume * multiplicity * sizeof_M);
  }
  bcast(get_data(multiplicity));
  bcast(get_data(sizeof_M));
  bcast(get_data(size_node));
  bcast(get_data(node_site));
  geo.init();
  geo.init(size_node * node_site);
  new_size_node = size_node;
}

template <class M>
Long dist_write_dist_data(const std::vector<DistData<M>>& dds,
                          const Int num_node, const std::string& path)
// interface_function
{
  SYNC_NODE();
  TIMER_VERBOSE_FLOPS("dist_write_dist_data");
  Long total_bytes = 0;
  Long total_ops = 0;
  const Int n_cycle = std::max(1, num_node / dist_write_par_limit());
  std::vector<Long> id_counts(num_node, 0);
  for (Int cycle = 0; cycle < n_cycle; cycle++) {
    Long bytes = 0;
    Long ops = 0;
    for (size_t k = 0; k < dds.size(); ++k) {
      const Int id_node = dds[k].id_node;
      qassert(0 <= id_node && id_node < num_node);
      if (id_node % n_cycle == cycle) {
        if (id_counts[id_node] == 0) {
          QFile fp = dist_open(path, id_node, num_node, "w");
          qassert(not fp.null());
          for (size_t l = k; l < dds.size(); ++l) {
            const DistData<M>& dd = dds[l];
            if (id_node == dd.id_node) {
              bytes += qwrite_data(get_data(dd), fp);
              ops += 1;
              id_counts[id_node] += 1;
            }
          }
          dist_close(fp);
        }
      }
    }
    glb_sum(bytes);
    glb_sum(ops);
    total_bytes += bytes;
    total_ops += ops;
    displayln_info(fname + ssprintf(": cycle / n_cycle = %4d / %4d ; total_ops "
                                    "= %10ld ; total_bytes = %15ld",
                                    cycle + 1, n_cycle, total_ops, total_bytes));
  }
  std::vector<Long> id_exists(num_node, 0);
  for (size_t id = 0; id < id_exists.size(); ++id) {
    id_exists[id] = id_counts[id] > 0 ? 1 : 0;
  }
  glb_sum(get_data(id_exists));
  glb_sum(get_data(id_counts));
  for (size_t id = 0; id < id_counts.size(); ++id) {
    qassert(id_exists[id] ==
            1);  // every id_node exist on one node, and one node only
    qassert(id_counts[id] ==
            id_counts[0]);  // every id_node has the same number of fields
  }
  std::vector<crc32_t> crcs = dist_crc32s(dds, num_node);
  crc32_t crc = dist_crc32(crcs);
  if (get_id_node() == 0) {
    const std::string fn = path + "/checksums.txt";
    QFile fp = qfopen(fn, "w");
    qassert(not fp.null());
    qwrite_data(ssprintf("%08X\n", crc), fp);
    qwrite_data("\n", fp);
    for (size_t i = 0; i < crcs.size(); ++i) {
      qwrite_data(ssprintf("%08X\n", crcs[i]), fp);
    }
    qfclose(fp);
  }
  qtouch_info(path + "/checkpoint");
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
Long dist_write_fields(const std::vector<ConstHandle<Field<M>>>& fs,
                       const Int num_node, const std::string& path)
{
  for (Int k = 0; k < (int)fs.size(); ++k) {
    const Field<M>& f = fs[k]();
    const Int id_node = f.geo().geon.id_node;
    if (id_node == 0) {
      dist_mkdir(path + ".partial", id_node);
      if (get_force_field_write_sizeof_M() == 0) {
        dist_write_geo_info(f.geo(), f.multiplicity, sizeof(M), path + ".partial");
      } else {
        const Int sizeof_M = get_force_field_write_sizeof_M();
        qassert((f.multiplicity * sizeof(M)) % sizeof_M == 0);
        const Int multiplicity = (f.multiplicity * sizeof(M)) / sizeof_M;
        dist_write_geo_info(f.geo(), multiplicity, sizeof_M, path + ".partial");
        get_force_field_write_sizeof_M() = 0;
      }
      break;
    }
  }
  SYNC_NODE();
  std::vector<DistData<M>> dds(fs.size());
  for (size_t i = 0; i < dds.size(); ++i) {
    dds[i].id_node = fs[i]().geo().geon.id_node;
    dds[i].data = get_data(fs[i]());
  }
  const Long total_bytes =
      dist_write_dist_data(dds, num_node, path + ".partial");
  qrename_info(path + ".partial", path);
  return total_bytes;
}

template <class M>
Long dist_write_fields(const std::vector<Field<M>>& fs, const Int num_node,
                       const std::string& path)
{
  std::vector<ConstHandle<Field<M>>> fhs(fs.size());
  for (size_t i = 0; i < fs.size(); ++i) {
    fhs[i].init(fs[i]);
  }
  return dist_write_fields(fhs, num_node, path);
}

template <class M>
Long dist_write_field(const Field<M>& f, const std::string& path)
// interface_function
{
  TIMER_VERBOSE("dist_write_field");
  qassert(f.geo().is_only_local);
  std::vector<ConstHandle<Field<M>>> fs(1);
  fs[0].init(f);
  return dist_write_fields(fs, get_num_node(), path);
}

template <class M>
Long dist_read_dist_data(const std::vector<DistData<M>>& dds,
                         const Int num_node, const std::string& path)
// interface_function
{
  if (!does_file_exist_sync_node(path + "/checkpoint")) {
    return 0;
  }
  SYNC_NODE();
  TIMER_VERBOSE_FLOPS("dist_read_dist_data");
  Long total_bytes = 0;
  Long total_ops = 0;
  const Int n_cycle = std::max(1, num_node / dist_read_par_limit());
  std::vector<Long> id_counts(num_node, 0);
  for (Int cycle = 0; cycle < n_cycle; cycle++) {
    Long bytes = 0;
    Long ops = 0;
    for (size_t k = 0; k < dds.size(); ++k) {
      const Int id_node = dds[k].id_node;
      qassert(0 <= id_node && id_node < num_node);
      if (id_node % n_cycle == cycle) {
        if (id_counts[id_node] == 0) {
          QFile fp = dist_open(path, id_node, num_node, "r");
          qassert(not fp.null());
          for (size_t l = k; l < dds.size(); ++l) {
            const DistData<M>& dd = dds[l];
            if (id_node == dd.id_node) {
              bytes += qread_data(get_data(dd), fp);
              ops += 1;
              id_counts[id_node] += 1;
            }
          }
          dist_close(fp);
        }
      }
    }
    glb_sum(bytes);
    glb_sum(ops);
    total_bytes += bytes;
    total_ops += ops;
    displayln_info(fname + ssprintf(": cycle / n_cycle = %4d / %4d ; total_ops "
                                    "= %10ld ; total_bytes = %15ld",
                                    cycle + 1, n_cycle, total_ops, total_bytes));
  }
  std::vector<Long> id_exists(num_node, 0);
  for (size_t id = 0; id < id_exists.size(); ++id) {
    id_exists[id] = id_counts[id] > 0 ? 1 : 0;
  }
  glb_sum(get_data(id_exists));
  glb_sum(get_data(id_counts));
  for (size_t id = 0; id < id_counts.size(); ++id) {
    qassert(id_exists[id] ==
            1);  // every id_node exist on one node, and one node only
    qassert(id_counts[id] ==
            id_counts[0]);  // every id_node has the same number of fields
  }
  std::vector<crc32_t> crcs = dist_crc32s(dds, num_node);
  crc32_t crc = dist_crc32(crcs);
  const bool is_checking = is_checksum_mismatch();
  is_checksum_mismatch() = false;
  if (get_id_node() == 0) {
    const std::string fn = path + "/checksums.txt";
    if (not does_file_exist(fn)) {
      displayln(fname + ssprintf(": WARNING '%s' does not exist.", fn.c_str()));
    } else {
      const std::vector<std::string> lines = qgetlines(fn);
      for (size_t i = 0; i < crcs.size(); ++i) {
        if (read_crc32(lines[i + 2]) != crcs[i]) {
          displayln(
              fname +
              ssprintf(
                  ": WARNING: checksums of file ; i=%d ; checksum.txt=%08X ; "
                  "computed=%08X ; path=%s",
                  i, read_crc32(lines[i + 2]), crcs[i], path.c_str()));
          is_checksum_mismatch() = true;
        }
      }
      if (read_crc32(lines[0]) != crc) {
        displayln(
            fname +
            ssprintf(": WARNING: checksums of files ; checksum.txt=%08X ; "
                     "computed=%08X ; path=%s",
                     read_crc32(lines[0]), crc, path.c_str()));
        is_checksum_mismatch() = true;
      }
    }
  }
  if (is_checking) {
    qassert(is_checksum_mismatch() == false);
  }
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
Long dist_read_fields(std::vector<Field<M>>& fs, Geometry& geo, Int& multiplicity,
                      Coordinate& new_size_node, const std::string& path)
// will clear fs before read
{
  if (!does_file_exist_sync_node(path + "/checkpoint")) {
    return 0;
  }
  clear(fs);
  Int sizeof_M;
  dist_read_geo_info(geo, multiplicity, sizeof_M, new_size_node, path);
  get_incorrect_field_read_sizeof_M() = 0;
  if ((int)sizeof(M) != sizeof_M) {
    displayln_info(
        "dist_read_fields: WARNING: sizeof(M) do not match with data on disk.");
    get_incorrect_field_read_sizeof_M() = sizeof_M;
    if (multiplicity * sizeof_M % sizeof(M) != 0) {
      displayln_info(
          ssprintf("dist_read_fields: ERROR: multiplicity = %d ; sizeof_M "
                   "= %d ; sizeof(M) = %d",
                   multiplicity, sizeof_M, sizeof(M)));
      qassert(false);
    }
    multiplicity = multiplicity * sizeof_M / sizeof(M);
  }
  std::vector<Geometry> new_geos =
      make_dist_io_geos(geo.total_site(), new_size_node);
  fs.resize(new_geos.size());
  std::vector<DistData<M>> dds(fs.size());
  for (size_t i = 0; i < fs.size(); ++i) {
    fs[i].init(new_geos[i], multiplicity);
    dds[i].id_node = fs[i].geo().geon.id_node;
    dds[i].data = get_data(fs[i]);
  }
  return dist_read_dist_data(dds, product(new_size_node), path);
}

}  // namespace qlat
