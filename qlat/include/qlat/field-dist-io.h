// vim: set ts=2 sw=2 expandtab:

#pragma once

#include <qlat/config.h>
#include <qlat/field-shuffle.h>
#include <qlat/utils-io.h>

#include <qlat-utils/qar-cache.h>

namespace qlat
{  //

const int DATA_READ_WRITE_NUMBER_OF_DIRECTORIES = 32;

API inline int& dist_write_par_limit()
// qlat parameter
{
  static int npar = get_env_long_default("q_write_par_limit", 16);
  return npar;
}

API inline int& dist_read_par_limit()
// qlat parameter
{
  static int npar = get_env_long_default("q_read_par_limit", 16);
  return npar;
}

API inline int& get_force_field_write_sizeof_M()
{
  static int sizeof_M = 0;
  return sizeof_M;
}

API inline int& get_incorrect_field_read_sizeof_M()
{
  static int sizeof_M = 0;
  return sizeof_M;
}

API inline bool& is_checksum_missmatch()
// qlat parameter
{
  static bool b = false;
  return b;
}

inline int dist_mkdir(const std::string& path, const int num_node,
                      const mode_t mode = default_dir_mode())
{
  int ret = 0;
  qmkdir(path, mode);
  for (int i = 0; i < std::min(num_node, DATA_READ_WRITE_NUMBER_OF_DIRECTORIES);
       ++i) {
    const std::string dir = path + ssprintf("/%02d", i);
    ret += qmkdir(dir, mode);
  }
  return ret;
}

inline int compute_dist_file_dir_id(const int id_node, const int num_node)
{
  const int ndir = std::min(num_node, DATA_READ_WRITE_NUMBER_OF_DIRECTORIES);
  const int dir_size = (num_node - 1) / ndir + 1;
  const int id_dir = id_node / dir_size;
  qassert(0 <= id_dir && id_dir < ndir);
  return id_dir;
}

inline std::string dist_file_name(const std::string& path, const int id_node,
                                  const int num_node)
{
  return path + ssprintf("/%02d/%010d",
                         compute_dist_file_dir_id(id_node, num_node), id_node);
}

inline FILE* dist_open(const std::string& path, const int id_node,
                       const int num_node,
                       const std::string& fmode,  // "w" for write, "r" for
                                                  // read, and "a" for append
                       const mode_t mode = default_dir_mode())
{
  const std::string fn = dist_file_name(path, id_node, num_node);
  FILE* ret = qopen(fn, fmode);
  if (ret == NULL && fmode != "r") {
    check_dir(path, mode);
    check_dir(
        path + ssprintf("/%02d", compute_dist_file_dir_id(id_node, num_node)),
        mode);
    ret = qopen(fn, fmode);
  }
  if (ret == NULL) {
    qwarn("dist_open: " + ssprintf("failed to open '%s'.", fn.c_str()));
  }
  return ret;
}

inline int dist_close(FILE*& fp) { return qclose(fp); }

template <class M>
struct DistData {
  int id_node;
  Vector<M> data;
};

template <class M>
Vector<M> get_data(const DistData<M>& dd)
{
  return dd.data;
}

template <class M>
std::vector<crc32_t> dist_crc32s(const std::vector<DistData<M> >& dds,
                                 const int num_node)
{
  sync_node();
  TIMER_VERBOSE_FLOPS("dist_crc32s");
  long total_bytes = 0;
  std::vector<crc32_t> ret(num_node, 0);
  for (int k = 0; k < (int)dds.size(); ++k) {
    const DistData<M>& dd = dds[k];
    const int id_node = dd.id_node;
    ret[id_node] = crc32_par(ret[id_node], dd.data);
    total_bytes += dd.data.data_size();
  }
  glb_sum_byte_vec(get_data(ret));
  glb_sum(total_bytes);
  timer.flops += total_bytes;
  return ret;
}

inline crc32_t dist_crc32(const std::vector<crc32_t>& crcs)
{
  std::vector<crc32_t> d = crcs;
  to_from_big_endian_32(get_data(d));
  return crc32(get_data(d));
}

template <class M>
crc32_t dist_crc32(const std::vector<DistData<M> >& dds, const int num_node)
{
  return dist_crc32(dist_crc32s(dds, num_node));
}

template <class M>
crc32_t field_dist_crc32(const Field<M>& f)
{
  std::vector<DistData<M> > dds(1);
  dds[0].id_node = f.geo().geon.id_node;
  dds[0].data = get_data(f);
  return dist_crc32(dds, get_num_node());
}

inline void dist_write_geo_info(const Geometry& geo, const int sizeof_M,
                                const std::string& path,
                                const mode_t mode = default_dir_mode())
{
  TIMER("dist_write_geo_info");
  const int id_node = geo.geon.id_node;
  qassert(geo.is_only_local());
  if (0 == id_node) {
    check_dir(path, mode);
    const std::string fn = path + "/geo-info.txt";
    FILE* fp = qopen(fn, "w");
    displayln(ssprintf("node_file_size = %ld",
                       sizeof_M * geo.multiplicity * geo.local_volume()),
              fp);
    displayln(ssprintf("geo.multiplicity = %d", geo.multiplicity), fp);
    displayln(ssprintf("sizeof(M) = %d", sizeof_M), fp);
    displayln(ssprintf("geo.geon.num_node = %d", geo.geon.num_node), fp);
    displayln(ssprintf("geo.geon.size_node[0] = %d", geo.geon.size_node[0]),
              fp);
    displayln(ssprintf("geo.geon.size_node[1] = %d", geo.geon.size_node[1]),
              fp);
    displayln(ssprintf("geo.geon.size_node[2] = %d", geo.geon.size_node[2]),
              fp);
    displayln(ssprintf("geo.geon.size_node[3] = %d", geo.geon.size_node[3]),
              fp);
    displayln(ssprintf("geo.local_volume() = %ld", geo.local_volume()), fp);
    displayln(ssprintf("geo.node_site[0] = %d", geo.node_site[0]), fp);
    displayln(ssprintf("geo.node_site[1] = %d", geo.node_site[1]), fp);
    displayln(ssprintf("geo.node_site[2] = %d", geo.node_site[2]), fp);
    displayln(ssprintf("geo.node_site[3] = %d", geo.node_site[3]), fp);
    displayln(ssprintf("PI = %.20f", PI), fp);
    const char* pic = (const char*)&PI;
    displayln(
        ssprintf("PI_double = %hhx %hhx %hhx %hhx %hhx %hhx %hhx %hhx", pic[0],
                 pic[1], pic[2], pic[3], pic[4], pic[5], pic[6], pic[7]),
        fp);
    const float PIf = PI;
    const char* pifc = (const char*)&PIf;
    displayln(ssprintf("PI_float = %hhx %hhx %hhx %hhx", pifc[0], pifc[1],
                       pifc[2], pifc[3]),
              fp);
    qclose(fp);
  }
}

inline std::string info_get_prop(const std::vector<std::string>& lines,
                                 const std::string& prop)
{
  for (size_t i = 0; i < lines.size(); ++i) {
    if (lines[i].compare(0, prop.size(), prop) == 0) {
      return std::string(lines[i], prop.size());
    }
  }
  return std::string("");
}

inline std::string info_get_prop(const std::vector<std::string>& lines,
                                 const std::string& prop,
                                 const std::string& prop1)
{
  const std::string ret = info_get_prop(lines, prop);
  if (ret != std::string("")) {
    return ret;
  } else {
    return info_get_prop(lines, prop1);
  }
}

inline void dist_read_geo_info(Geometry& geo, int& sizeof_M,
                               Coordinate& new_size_node,
                               const std::string& path)
{
  TIMER("dist_read_geo_info");
  int multiplicity;
  Coordinate size_node;
  Coordinate node_site;
  if (get_id_node() == 0) {
    const std::string fn = path + "/geo-info.txt";
    const std::vector<std::string> lines = qgetlines(fn);
    reads(multiplicity, info_get_prop(lines, "geo.multiplicity = "));
    reads(sizeof_M, info_get_prop(lines, "sizeof(M) = "));
    for (int i = 0; i < 4; ++i) {
      reads(size_node[i],
            info_get_prop(lines, ssprintf("geo.geon.size_node[%d] = ", i),
                          ssprintf("geo.sizeNode[%d] = ", i)));
      reads(node_site[i],
            info_get_prop(lines, ssprintf("geo.node_site[%d] = ", i),
                          ssprintf("geo.nodeSite[%d] = ", i)));
    }
    long node_file_size;
    int num_node;
    long local_volume;
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
  geo.init(size_node * node_site, multiplicity);
  new_size_node = size_node;
}

template <class M>
long dist_write_dist_data(const std::vector<DistData<M> >& dds,
                          const int num_node, const std::string& path,
                          const mode_t mode = default_dir_mode())
// interface_function
{
  sync_node();
  TIMER_VERBOSE_FLOPS("dist_write_dist_data");
  long total_bytes = 0;
  long total_ops = 0;
  const int n_cycle = std::max(1, num_node / dist_write_par_limit());
  std::vector<long> id_counts(num_node, 0);
  for (int i = 0; i < n_cycle; i++) {
    long bytes = 0;
    long ops = 0;
    for (size_t k = 0; k < dds.size(); ++k) {
      const int id_node = dds[k].id_node;
      qassert(0 <= id_node && id_node < num_node);
      if (id_node % n_cycle == i) {
        if (id_counts[id_node] == 0) {
          FILE* fp = dist_open(path, id_node, num_node, "w", mode);
          qassert(fp != NULL);
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
                                    i + 1, n_cycle, total_ops, total_bytes));
  }
  std::vector<long> id_exists(num_node, 0);
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
    FILE* fp = qopen(fn, "w");
    qassert(fp != NULL);
    displayln(ssprintf("%08X", crc), fp);
    displayln("", fp);
    for (size_t i = 0; i < crcs.size(); ++i) {
      displayln(ssprintf("%08X", crcs[i]), fp);
    }
    qclose(fp);
  }
  qtouch_info(path + "/checkpoint");
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
long dist_write_fields(const std::vector<ConstHandle<Field<M> > >& fs,
                       const int num_node, const std::string& path,
                       const mode_t mode = default_dir_mode())
{
  for (int k = 0; k < (int)fs.size(); ++k) {
    const Field<M>& f = fs[k]();
    const int id_node = f.geo().geon.id_node;
    if (id_node == 0) {
      dist_mkdir(path + ".partial", id_node, mode);
      if (get_force_field_write_sizeof_M() == 0) {
        dist_write_geo_info(f.geo(), sizeof(M), path + ".partial", mode);
      } else {
        const int sizeof_M = get_force_field_write_sizeof_M();
        qassert((f.geo().multiplicity * sizeof(M)) % sizeof_M == 0);
        const int multiplicity = (f.geo().multiplicity * sizeof(M)) / sizeof_M;
        dist_write_geo_info(geo_remult(f.geo(), multiplicity), sizeof_M,
                            path + ".partial", mode);
        get_force_field_write_sizeof_M() = 0;
      }
      break;
    }
  }
  sync_node();
  std::vector<DistData<M> > dds(fs.size());
  for (size_t i = 0; i < dds.size(); ++i) {
    dds[i].id_node = fs[i]().geo().geon.id_node;
    dds[i].data = get_data(fs[i]());
  }
  const long total_bytes =
      dist_write_dist_data(dds, num_node, path + ".partial", mode);
  qrename_info(path + ".partial", path);
  return total_bytes;
}

template <class M>
long dist_write_fields(const std::vector<Field<M> >& fs, const int num_node,
                       const std::string& path,
                       const mode_t mode = default_dir_mode())
{
  std::vector<ConstHandle<Field<M> > > fhs(fs.size());
  for (size_t i = 0; i < fs.size(); ++i) {
    fhs[i].init(fs[i]);
  }
  return dist_write_fields(fhs, num_node, path, mode);
}

template <class M>
long dist_write_field(const Field<M>& f, const std::string& path,
                      const mode_t mode = default_dir_mode())
// interface_function
{
  TIMER_VERBOSE("dist_write_field");
  qassert(f.geo().is_only_local());
  std::vector<ConstHandle<Field<M> > > fs(1);
  fs[0].init(f);
  return dist_write_fields(fs, get_num_node(), path, mode);
}

template <class M>
long dist_read_dist_data(const std::vector<DistData<M> >& dds,
                         const int num_node, const std::string& path)
// interface_function
{
  if (!does_file_exist_sync_node(path + "/checkpoint")) {
    return 0;
  }
  sync_node();
  TIMER_VERBOSE_FLOPS("dist_read_dist_data");
  long total_bytes = 0;
  long total_ops = 0;
  const int n_cycle = std::max(1, num_node / dist_read_par_limit());
  std::vector<long> id_counts(num_node, 0);
  for (int i = 0; i < n_cycle; i++) {
    long bytes = 0;
    long ops = 0;
    for (size_t k = 0; k < dds.size(); ++k) {
      const int id_node = dds[k].id_node;
      qassert(0 <= id_node && id_node < num_node);
      if (id_node % n_cycle == i) {
        if (id_counts[id_node] == 0) {
          FILE* fp = dist_open(path, id_node, num_node, "r");
          qassert(fp != NULL);
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
                                    i + 1, n_cycle, total_ops, total_bytes));
  }
  std::vector<long> id_exists(num_node, 0);
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
  is_checksum_missmatch() = false;
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
          is_checksum_missmatch() = true;
        }
      }
      if (read_crc32(lines[0]) != crc) {
        displayln(
            fname +
            ssprintf(": WARNING: checksums of files ; checksum.txt=%08X ; "
                     "computed=%08X ; path=%s",
                     read_crc32(lines[0]), crc, path.c_str()));
        is_checksum_missmatch() = true;
      }
    }
  }
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
long dist_read_fields(std::vector<Field<M> >& fs, Geometry& geo,
                      Coordinate& new_size_node, const std::string& path)
// will clear fs before read
{
  if (!does_file_exist_sync_node(path + "/checkpoint")) {
    return 0;
  }
  clear(fs);
  int sizeof_M;
  dist_read_geo_info(geo, sizeof_M, new_size_node, path);
  get_incorrect_field_read_sizeof_M() = 0;
  if ((int)sizeof(M) != sizeof_M) {
    displayln_info(
        "dist_read_fields: WARNING: sizeof(M) do not match with data on disk.");
    get_incorrect_field_read_sizeof_M() = sizeof_M;
    if (geo.multiplicity * sizeof_M % sizeof(M) != 0) {
      displayln_info(
          ssprintf("dist_read_fields: ERROR: geo.multiplicity = %d ; sizeof_M "
                   "= %d ; sizeof(M) = %d",
                   geo.multiplicity, sizeof_M, sizeof(M)));
      qassert(false);
    }
    geo.remult(geo.multiplicity * sizeof_M / sizeof(M));
  }
  std::vector<Geometry> new_geos =
      make_dist_io_geos(geo.total_site(), geo.multiplicity, new_size_node);
  fs.resize(new_geos.size());
  std::vector<DistData<M> > dds(fs.size());
  for (size_t i = 0; i < fs.size(); ++i) {
    fs[i].init(new_geos[i]);
    dds[i].id_node = fs[i].geo().geon.id_node;
    dds[i].data = get_data(fs[i]);
  }
  return dist_read_dist_data(dds, product(new_size_node), path);
}

template <class M>
long dist_write_field(const Field<M>& f, const Coordinate new_size_node,
                      const std::string& path,
                      const mode_t mode = default_dir_mode())
// interface_function
{
  TIMER_VERBOSE_FLOPS("dist_write_field");
  displayln_info(fname + ssprintf(": fn='%s'.", path.c_str()));
  std::vector<Field<M> > fs;
  shuffle_field(fs, f, new_size_node);
  long total_bytes = dist_write_fields(fs, product(new_size_node), path, mode);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
long dist_read_field(Field<M>& f, const std::string& path)
// interface_function
{
  TIMER_VERBOSE_FLOPS("dist_read_field");
  displayln_info(fname + ssprintf(": fn='%s'.", path.c_str()));
  Geometry geo;
  std::vector<Field<M> > fs;
  Coordinate new_size_node;
  const long total_bytes = dist_read_fields(fs, geo, new_size_node, path);
  if (total_bytes == 0) {
    return 0;
  } else {
    f.init(geo);
    qassert(f.geo() == geo);
    shuffle_field_back(f, fs, new_size_node);
    timer.flops += total_bytes;
    return total_bytes;
  }
}

template <class M, class N>
void convert_field_float_from_double(Field<N>& ff, const Field<M>& f)
// interface_function
{
  TIMER("convert_field_float_from_double");
  qassert(f.geo().is_only_local());
  qassert(sizeof(M) % sizeof(double) == 0);
  qassert(sizeof(N) % sizeof(float) == 0);
  qassert(f.geo().multiplicity * sizeof(M) / 2 % sizeof(N) == 0);
  const int multiplicity = f.geo().multiplicity * sizeof(M) / 2 / sizeof(N);
  const Geometry geo = geo_remult(f.geo(), multiplicity);
  ff.init(geo);
  const Vector<M> fdata = get_data(f);
  const Vector<double> fd((double*)fdata.data(),
                          fdata.data_size() / sizeof(double));
  Vector<N> ffdata = get_data(ff);
  Vector<float> ffd((float*)ffdata.data(), ffdata.data_size() / sizeof(float));
  qassert(ffd.size() == fd.size());
  qacc_for(i, ffd.size(), {
    ffd[i] = fd[i];
  });
}

template <class M, class N>
void convert_field_double_from_float(Field<N>& ff, const Field<M>& f)
// interface_function
{
  TIMER("convert_field_double_from_float");
  qassert(f.geo().is_only_local());
  qassert(sizeof(M) % sizeof(float) == 0);
  qassert(sizeof(N) % sizeof(double) == 0);
  qassert(f.geo().multiplicity * sizeof(M) * 2 % sizeof(N) == 0);
  const int multiplicity = f.geo().multiplicity * sizeof(M) * 2 / sizeof(N);
  const Geometry geo = geo_remult(f.geo(), multiplicity);
  ff.init(geo);
  const Vector<M> fdata = get_data(f);
  const Vector<float> fd((float*)fdata.data(),
                         fdata.data_size() / sizeof(float));
  Vector<N> ffdata = get_data(ff);
  Vector<double> ffd((double*)ffdata.data(),
                     ffdata.data_size() / sizeof(double));
  qassert(ffd.size() == fd.size());
  qacc_for(i, ffd.size(), {
    ffd[i] = fd[i];
  });
}

template <class M>
long dist_write_field_float_from_double(const Field<M>& f,
                                        const std::string& path,
                                        const mode_t mode = default_dir_mode())
// interface_function
{
  TIMER_VERBOSE_FLOPS("dist_write_field_float_from_double");
  Field<float> ff;
  convert_field_float_from_double(ff, f);
  to_from_big_endian_32(get_data(ff));
  const long total_bytes = dist_write_field(ff, path, mode);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
long dist_write_field_float_from_double(const Field<M>& f,
                                        const Coordinate& new_size_node,
                                        const std::string& path,
                                        const mode_t mode = default_dir_mode())
// interface_function
{
  TIMER_VERBOSE_FLOPS("dist_write_field_float_from_double");
  Field<float> ff;
  convert_field_float_from_double(ff, f);
  to_from_big_endian_32(get_data(ff));
  const long total_bytes = dist_write_field(ff, new_size_node, path, mode);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
long dist_read_field_double_from_float(Field<M>& f, const std::string& path)
// interface_function
{
  TIMER_VERBOSE_FLOPS("dist_read_field_double_from_float");
  Field<float> ff;
  const long total_bytes = dist_read_field(ff, path);
  if (total_bytes == 0) {
    return 0;
  } else {
    to_from_big_endian_32(get_data(ff));
    convert_field_double_from_float(f, ff);
    timer.flops += total_bytes;
    return total_bytes;
  }
}

template <class M>
long dist_write_field_double(const Field<M>& f, const std::string& path,
                             const mode_t mode = default_dir_mode())
// interface_function
{
  TIMER_VERBOSE_FLOPS("dist_write_field_double");
  Field<M> ff;
  ff.init(f);
  to_from_big_endian_64(get_data(ff));
  const long total_bytes = dist_write_field(ff, path, mode);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
long dist_write_field_double(const Field<M>& f, const Coordinate& new_size_node,
                             const std::string& path,
                             const mode_t mode = default_dir_mode())
// interface_function
{
  TIMER_VERBOSE_FLOPS("dist_write_field_double");
  Field<M> ff;
  ff.init(f);
  to_from_big_endian_64(get_data(ff));
  const long total_bytes = dist_write_field(ff, new_size_node, path, mode);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
long dist_read_field_double(Field<M>& f, const std::string& path)
// interface_function
{
  TIMER_VERBOSE_FLOPS("dist_read_field_double");
  const long total_bytes = dist_read_field(f, path);
  if (total_bytes == 0) {
    return 0;
  } else {
    to_from_big_endian_64(get_data(f));
    timer.flops += total_bytes;
    return total_bytes;
  }
}

inline bool is_dist_field(const std::string& path)
{
  TIMER("is_dist_field");
  long nfile = 0;
  if (get_id_node() == 0) {
    if (does_file_exist(path + "/geo-info.txt") and
        does_file_exist(path + "/checkpoint")) {
      nfile = 1;
    }
  }
  bcast(get_data(nfile));
  return nfile > 0;
}

}  // namespace qlat
