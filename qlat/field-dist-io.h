// vim: set ts=2 sw=2 expandtab:

#pragma once

#include <qlat/config.h>
#include <qlat/utils.h>
#include <qlat/crc32.h>

QLAT_START_NAMESPACE

const int DATA_READ_WRITE_NUMBER_OF_DIRECTORIES = 32;

inline int& dist_write_par_limit() {
  static int npar = 3;
  return npar;
}

inline int& dist_read_par_limit() {
  static int npar = 3;
  return npar;
}

inline std::vector<Geometry> make_dist_io_geos(const Coordinate& total_site, const int multiplicity, const Coordinate& new_size_node)
{
  TIMER("make_dist_io_geos");
  const Coordinate new_node_site = total_site / new_size_node;
  const int new_num_node = product(new_size_node);
  std::vector<Geometry> ret;
  const int min_size_chunk = new_num_node / get_num_node();
  const int remain = new_num_node % get_num_node();
  const int size_chunk = get_id_node() < remain ? min_size_chunk + 1 : min_size_chunk;
  const int chunk_start = get_id_node() * min_size_chunk + (get_id_node() < remain ? get_id_node() : remain);
  const int chunk_end = std::min(new_num_node, chunk_start + size_chunk);
  for (int new_id_node = chunk_start; new_id_node < chunk_end; ++new_id_node) {
    GeometryNode geon;
    geon.initialized = true;
    geon.num_node = new_num_node;
    geon.id_node = new_id_node;
    geon.size_node = new_size_node;
    geon.coor_node = coordinate_from_index(new_id_node, new_size_node);
    Geometry new_geo;
    new_geo.init(geon, multiplicity, new_node_site);
    ret.push_back(new_geo);
  }
  return ret;
}

inline int dist_mkdir(const std::string& path, const int num_node, const mode_t mode = default_dir_mode())
{
  int ret = 0;
  qmkdir(path, mode);
  for (int i = 0; i < std::min(num_node, DATA_READ_WRITE_NUMBER_OF_DIRECTORIES); ++i) {
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

inline std::string dist_file_name(const std::string& path, const int id_node, const int num_node)
{
  return path + ssprintf("/%02d/%010d", compute_dist_file_dir_id(id_node, num_node), id_node);
}

inline FILE* dist_open(const std::string& path, const int id_node, const int num_node,
    const std::string& fmode, // "w" for write, "r" for read, and "a" for append
    const mode_t mode = default_dir_mode())
{
  const std::string fn = dist_file_name(path, id_node, num_node);
  FILE* ret = qopen(fn, fmode);
  if (ret == NULL && fmode != "r") {
    check_dir(path, mode);
    check_dir(path + ssprintf("/%02d", compute_dist_file_dir_id(id_node, num_node)), mode);
    ret = qopen(fn, fmode);
  }
  qassert(ret != NULL);
  return ret;
}

inline int dist_close(FILE*& fp)
{
  return qclose(fp);
}

template <class M>
long qwrite_data(const Vector<M>& v, FILE* fp)
{
  return sizeof(M) * std::fwrite((void*)v.p, sizeof(M), v.n, fp);
}

template <class M>
long qread_data(const Vector<M>& v, FILE* fp)
{
  return sizeof(M) * std::fread((void*)v.p, sizeof(M), v.n, fp);
}

template <class M>
struct DistData
{
  int id_node;
  Vector<M> data;
};

template <class M>
Vector<M> get_data(const DistData<M>& dd)
{
  return dd.data;
}

template <class M>
std::vector<crc32_t> dist_crc32s(const std::vector<DistData<M> >& dds, const int num_node)
{
  sync_node();
  TIMER_VERBOSE_FLOPS("dist_crc32s");
  long total_bytes = 0;
  std::vector<crc32_t> ret(num_node, 0);
  for (int k = 0; k < dds.size(); ++k) {
    const DistData<M>& dd = dds[k];
    const int id_node = dd.id_node;
    ret[id_node] = crc32(ret[id_node], dd.data);
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
crc32_t field_crc32(const Field<M>& f)
{
  std::vector<DistData<M> > dds(1);
  dds[0].id_node = f.geo.geon.id_node;
  dds[0].data = get_data(f);
  return dist_crc32(dds, get_num_node());
}

inline void dist_write_geo_info(const Geometry& geo, const size_t sizeof_M,
    const std::string& path, const mode_t mode = default_dir_mode())
{
  TIMER("dist_write_geo_info");
  const int id_node = geo.geon.id_node;
  const int num_node = geo.geon.num_node;
  qassert(geo.is_only_local());
  if (0 == id_node) {
    check_dir(path, mode);
    const std::string fn = path + "/geo-info.txt";
    FILE* fp = qopen(fn, "w");
    displayln(ssprintf("node_file_size = %ld", sizeof_M * geo.multiplicity * geo.local_volume()), fp);
    displayln(ssprintf("geo.multiplicity = %d", geo.multiplicity), fp);
    displayln(ssprintf("sizeof(M) = %d", sizeof_M), fp);
    displayln(ssprintf("geo.geon.num_node = %d", geo.geon.num_node), fp);
    displayln(ssprintf("geo.geon.size_node[0] = %d", geo.geon.size_node[0]), fp);
    displayln(ssprintf("geo.geon.size_node[1] = %d", geo.geon.size_node[1]), fp);
    displayln(ssprintf("geo.geon.size_node[2] = %d", geo.geon.size_node[2]), fp);
    displayln(ssprintf("geo.geon.size_node[3] = %d", geo.geon.size_node[3]), fp);
    displayln(ssprintf("geo.local_volume() = %ld", geo.local_volume()), fp);
    displayln(ssprintf("geo.node_site[0] = %d", geo.node_site[0]), fp);
    displayln(ssprintf("geo.node_site[1] = %d", geo.node_site[1]), fp);
    displayln(ssprintf("geo.node_site[2] = %d", geo.node_site[2]), fp);
    displayln(ssprintf("geo.node_site[3] = %d", geo.node_site[3]), fp);
    displayln(ssprintf("PI = %.20f", PI), fp);
    const char* pic = (const char*)&PI;
    displayln(ssprintf("PI_double = %hhx %hhx %hhx %hhx %hhx %hhx %hhx %hhx",
        pic[0], pic[1], pic[2], pic[3], pic[4], pic[5], pic[6], pic[7]), fp);
    const float PIf = PI;
    const char* pifc = (const char*)&PIf;
    displayln(ssprintf("PI_float = %hhx %hhx %hhx %hhx",
        pifc[0], pifc[1], pifc[2], pifc[3]), fp);
    qclose(fp);
  }
}

inline std::string geo_info_get_prop(const std::vector<std::string>& lines, const std::string& prop)
{
  for (size_t i = 0; i < lines.size(); ++i) {
    if (lines[i].compare(0, prop.size(), prop) == 0) {
      return std::string(lines[i], prop.size());
    }
  }
  return std::string("");
}

inline std::string geo_info_get_prop(const std::vector<std::string>& lines, const std::string& prop, const std::string& prop1)
{
  const std::string ret = geo_info_get_prop(lines, prop);
  if (ret != std::string("")) {
    return ret;
  } else {
    return geo_info_get_prop(lines, prop1);
  }
}

inline void dist_read_geo_info(Geometry& geo, size_t& sizeof_M, Coordinate& new_size_node,
    const std::string& path)
{
  TIMER("dist_read_geo_info");
  int multiplicity;
  Coordinate size_node;
  Coordinate node_site;
  if (get_id_node() == 0) {
    const std::string fn = path + "/geo-info.txt";
    FILE* fp = qopen(fn, "r");
    qassert(fp != NULL);
    const std::vector<std::string> lines = qgetlines(fp);
    qclose(fp);
    reads(multiplicity, geo_info_get_prop(lines, "geo.multiplicity = "));
    reads(sizeof_M, geo_info_get_prop(lines, "sizeof(M) = "));
    for (int i = 0; i < 4; ++i) {
      reads(size_node[i], geo_info_get_prop(lines,
            ssprintf("geo.geon.size_node[%d] = ", i), ssprintf("geo.sizeNode[%d] = ", i)));
      reads(node_site[i], geo_info_get_prop(lines,
            ssprintf("geo.node_site[%d] = ", i), ssprintf("geo.nodeSite[%d] = ", i)));
    }
    long node_file_size;
    int num_node;
    long local_volume;
    reads(node_file_size, geo_info_get_prop(lines, "node_file_size = ", "nodeFileSize = "));
    reads(num_node, geo_info_get_prop(lines, "geo.geon.num_node = ", "geo.numNode = "));
    reads(local_volume, geo_info_get_prop(lines, "geo.local_volume() = ", "geo.localVolume() = "));
    qassert(num_node == product(size_node));
    qassert(local_volume == product(node_site));
    qassert(node_file_size == local_volume * multiplicity * sizeof_M);
  }
  bcast(get_data(multiplicity));
  bcast(get_data((long)sizeof_M));
  bcast(get_data(size_node));
  bcast(get_data(node_site));
  geo.init();
  geo.init(size_node * node_site, multiplicity);
  new_size_node = size_node;
}

template <class M>
long dist_write_dist_data(const std::vector<DistData<M> >& dds, const int num_node,
    const std::string& path, const mode_t mode = default_dir_mode())
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
    displayln_info(ssprintf("%s: cycle / n_cycle = %4d / %4d ; total_ops = %10ld ; total_bytes = %15ld",
          fname, i + 1, n_cycle, total_ops, total_bytes));
  }
  std::vector<long> id_exists(num_node, 0);
  for (size_t id = 0; id < id_exists.size(); ++id) {
    id_exists[id] = id_counts[id] > 0 ? 1 : 0;
  }
  glb_sum(get_data(id_exists));
  glb_sum(get_data(id_counts));
  for (size_t id = 0; id < id_counts.size(); ++id) {
    qassert(id_exists[id] == 1); // every id_node exist on one node, and one node only
    qassert(id_counts[id] == id_counts[0]); // every id_node has the same number of fields
  }
  std::vector<crc32_t> crcs = dist_crc32s(dds, num_node);
  crc32_t crc = dist_crc32(crcs);
  if (get_id_node() == 0) {
    const std::string fn = path + "/checksums.txt";
    FILE* fp = qopen(fn, "w");
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
    const int num_node,
    const std::string& path, const mode_t mode = default_dir_mode())
{
  for (int k = 0; k < fs.size(); ++k) {
    const Field<M>& f = fs[k]();
    const int id_node = f.geo.geon.id_node;
    if (id_node == 0) {
      dist_mkdir(path, id_node, mode);
      dist_write_geo_info(f.geo, sizeof(M), path, mode);
      break;
    }
  }
  sync_node();
  std::vector<DistData<M> > dds(fs.size());
  for (size_t i = 0; i < dds.size(); ++i) {
    dds[i].id_node = fs[i]().geo.geon.id_node;
    dds[i].data = get_data(fs[i]());
  }
  return dist_write_dist_data(dds, num_node, path, mode);
}

template <class M>
long dist_write_fields(const std::vector<Field<M> >& fs,
    const int num_node,
    const std::string& path, const mode_t mode = default_dir_mode())
{
  std::vector<ConstHandle<Field<M> > > fhs(fs.size());
  for (size_t i = 0; i < fs.size(); ++i) {
    fhs[i].init(fs[i]);
  }
  return dist_write_fields(fhs, num_node, path, mode);
}

template <class M>
long dist_write_field(const Field<M>& f,
    const std::string& path, const mode_t mode = default_dir_mode())
  // interface_function
{
  TIMER_VERBOSE("dist_write_field");
  qassert(f.geo.is_only_local());
  std::vector<ConstHandle<Field<M> > > fs(1);
  fs[0].init(f);
  return dist_write_fields(fs, get_num_node(), path, mode);
}

template <class M>
long dist_read_dist_data(const std::vector<DistData<M> >& dds, const int num_node, const std::string& path)
  // interface_function
{
  sync_node();
  TIMER_VERBOSE_FLOPS("dist_read_dist_data");
  if (!does_file_exist_sync_node(path + "/checkpoint")) {
    return 0;
  }
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
    displayln_info(ssprintf("%s: cycle / n_cycle = %4d / %4d ; total_ops = %10ld ; total_bytes = %15ld",
          fname, i + 1, n_cycle, total_ops, total_bytes));
  }
  std::vector<long> id_exists(num_node, 0);
  for (size_t id = 0; id < id_exists.size(); ++id) {
    id_exists[id] = id_counts[id] > 0 ? 1 : 0;
  }
  glb_sum(get_data(id_exists));
  glb_sum(get_data(id_counts));
  for (size_t id = 0; id < id_counts.size(); ++id) {
    qassert(id_exists[id] == 1); // every id_node exist on one node, and one node only
    qassert(id_counts[id] == id_counts[0]); // every id_node has the same number of fields
  }
  std::vector<crc32_t> crcs = dist_crc32s(dds, num_node);
  crc32_t crc = dist_crc32(crcs);
  if (get_id_node() == 0) {
    const std::string fn = path + "/checksums.txt";
    FILE* fp = qopen(fn, "r");
    crc32_t crc_read;
    std::fscanf(fp, "%X\n", &crc_read);
    qassert(crc == crc_read);
    std::fscanf(fp, "\n");
    for (size_t i = 0; i < crcs.size(); ++i) {
      std::fscanf(fp, "%X\n", &crc_read);
      qassert(crcs[i] == crc_read);
    }
    qclose(fp);
  }
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
long dist_read_fields(std::vector<Field<M> >& fs, Geometry& geo, Coordinate& new_size_node, const std::string& path)
  // will clear fs before read
{
  if (!does_file_exist_sync_node(path + "/checkpoint")) {
    return 0;
  }
  fs.clear();
  size_t sizeof_M;
  dist_read_geo_info(geo, sizeof_M, new_size_node, path);
  std::vector<Geometry> new_geos = make_dist_io_geos(geo.total_site(), geo.multiplicity, new_size_node);
  fs.resize(new_geos.size());
  std::vector<DistData<M> > dds(fs.size());
  for (size_t i = 0; i < fs.size(); ++i) {
    fs[i].init(new_geos[i]);
    dds[i].id_node = fs[i].geo.geon.id_node;
    dds[i].data = get_data(fs[i]);
  }
  return dist_read_dist_data(dds, product(new_size_node), path);
}

inline int get_id_node_from_new_id_node(const int new_id_node, const int new_num_node, const int num_node)
{
  const int min_size_chunk = new_num_node / num_node;
  const int remain = new_num_node % num_node;
  const int limit = remain * min_size_chunk + remain;
  if (new_id_node <= limit) {
    return new_id_node / (min_size_chunk + 1);
  } else {
    return (new_id_node - limit) / min_size_chunk + remain;
  }
}

struct ShufflePlanKey
{
  Coordinate total_site;
  Coordinate new_size_node;
};

inline bool operator<(const ShufflePlanKey& x, const ShufflePlanKey& y)
{
  if (x.new_size_node < y.new_size_node) {
    return true;
  } else if (y.new_size_node < x.new_size_node) {
    return false;
  } else {
    return x.total_site < y.total_site;
  }
}

struct ShufflePlanRecvPackInfo
{
  int local_geos_idx;
  long field_idx;
  long buffer_idx;
  long size;
};

struct ShufflePlanSendPackInfo
{
  long field_idx;
  long buffer_idx;
  long size;
};

struct ShufflePlanMsgInfo
{
  int id_node;
  long idx;
  long size;
};

struct ShufflePlan
{
  long total_send_size; // total send buffer size
  std::vector<ShufflePlanMsgInfo> send_msg_infos; // corresponds to every sent msg
  std::vector<ShufflePlanSendPackInfo> send_pack_infos; // corresponds to how to create send buffer from local field
  long total_recv_size;
  std::vector<ShufflePlanMsgInfo> recv_msg_infos; // corresponds to every recv msg
  std::vector<ShufflePlanRecvPackInfo> recv_pack_infos; // corresponds to how to copy recv buffer to new local fields
};

inline ShufflePlan make_shuffle_plan(const ShufflePlanKey& spk)
{
  TIMER_VERBOSE("make_shuffle_plan");
  const Coordinate& total_site = spk.total_site;
  const Coordinate& new_size_node = spk.new_size_node;
  const Coordinate new_node_site = total_site / new_size_node;
  qassert(new_size_node * new_node_site == total_site);
  const int new_num_node = product(new_size_node);
  Geometry geo;
  geo.init(total_site, 1);
  const int num_node = geo.geon.num_node;
  std::vector<Geometry> new_geos = make_dist_io_geos(geo.total_site(), geo.multiplicity, new_size_node);
  ShufflePlan ret;
  // total_send_size
  ret.total_send_size = geo.local_volume();
  // send_id_node_size
  // send_new_id_node_size
  std::map<int,long> send_id_node_size;
  std::map<int,long> send_new_id_node_size;
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const Coordinate new_coor_node = xg / new_node_site;
    const int new_id_node = index_from_coordinate(new_coor_node, new_size_node);
    const int id_node = get_id_node_from_new_id_node(new_id_node, new_num_node, num_node);
    send_id_node_size[id_node] += 1;
    send_new_id_node_size[new_id_node] += 1;
  }
  // send_id_node_idx
  std::map<int,long> send_id_node_idx;
  {
    long count = 0;
    for (std::map<int,long>::const_iterator it = send_id_node_size.cbegin(); it != send_id_node_size.cend(); ++it) {
      const int id_node = it->first;
      const int node_size = it->second;
      ShufflePlanMsgInfo mi;
      mi.id_node = id_node;
      mi.idx = count;
      mi.size = node_size;
      ret.send_msg_infos.push_back(mi);
      send_id_node_idx[id_node] = count;
      count += node_size;
    }
    qassert(count == geo.local_volume());
    qassert(count == ret.total_send_size);
  }
  // send_new_id_node_idx
  std::map<int,long> send_new_id_node_idx;
  {
    long count = 0;
    for (std::map<int,long>::const_iterator it = send_new_id_node_size.cbegin(); it != send_new_id_node_size.cend(); ++it) {
      const int new_id_node = it->first;
      const int node_size = it->second;
      send_new_id_node_idx[new_id_node] = count;
      count += node_size;
    }
    qassert(count == geo.local_volume());
    qassert(count == ret.total_send_size);
  }
  // send_pack_infos
  {
    long last_buffer_idx = -1;
    for (long index = 0; index < geo.local_volume(); ++index) {
      const Coordinate xl = geo.coordinate_from_index(index);
      const Coordinate xg = geo.coordinate_g_from_l(xl);
      const Coordinate new_coor_node = xg / new_node_site;
      const int new_id_node = index_from_coordinate(new_coor_node, new_size_node);
      const long buffer_idx = send_new_id_node_idx[new_id_node];
      // const int id_node = get_id_node_from_new_id_node(new_id_node, new_num_node, num_node);
      if (buffer_idx != last_buffer_idx) {
        last_buffer_idx = buffer_idx;
        ShufflePlanSendPackInfo pi;
        pi.field_idx = index;
        pi.buffer_idx = buffer_idx;
        pi.size = 1;
        ret.send_pack_infos.push_back(pi);
      } else {
        ShufflePlanSendPackInfo& pi = ret.send_pack_infos.back();
        pi.size += 1;
      }
      send_new_id_node_idx[new_id_node] += 1;
      last_buffer_idx += 1;
    }
  }
  // total_recv_size
  ret.total_recv_size = 0;
  for (size_t i = 0; i < new_geos.size(); ++i) {
    const Geometry& new_geo = new_geos[i];
    ret.total_recv_size += new_geo.local_volume();
  }
  // recv_id_node_size
  std::map<int,long> recv_id_node_size;
  for (size_t i = 0; i < new_geos.size(); ++i) {
    const Geometry& new_geo = new_geos[i];
    for (long index = 0; index < new_geo.local_volume(); ++index) {
      const Coordinate xl = new_geo.coordinate_from_index(index);
      const Coordinate xg = new_geo.coordinate_g_from_l(xl);
      const Coordinate coor_node = xg / geo.node_site;
      const long id_node = index_from_coordinate(coor_node, geo.geon.size_node);
      recv_id_node_size[id_node] += 1;
    }
  }
  // recv_id_node_idx
  std::map<int,long> recv_id_node_idx;
  {
    long count = 0;
    for (std::map<int,long>::const_iterator it = recv_id_node_size.cbegin(); it != recv_id_node_size.cend(); ++it) {
      const int id_node = it->first;
      const int node_size = it->second;
      ShufflePlanMsgInfo mi;
      mi.id_node = id_node;
      mi.idx = count;
      mi.size = node_size;
      ret.recv_msg_infos.push_back(mi);
      recv_id_node_idx[id_node] = count;
      count += node_size;
    }
    qassert(count == ret.total_recv_size);
  }
  // recv_pack_infos
  {
    for (size_t i = 0; i < new_geos.size(); ++i) {
      const Geometry& new_geo = new_geos[i];
      long last_buffer_idx = -1;
      for (long index = 0; index < new_geo.local_volume(); ++index) {
        const Coordinate xl = new_geo.coordinate_from_index(index);
        const Coordinate xg = new_geo.coordinate_g_from_l(xl);
        const Coordinate coor_node = xg / geo.node_site;
        const int id_node = index_from_coordinate(coor_node, geo.geon.size_node);
        const long buffer_idx = recv_id_node_idx[id_node];
        if (buffer_idx != last_buffer_idx) {
          last_buffer_idx = buffer_idx;
          ShufflePlanRecvPackInfo pi;
          pi.local_geos_idx = i;
          pi.field_idx = index;
          pi.buffer_idx = buffer_idx;
          pi.size = 1;
          ret.recv_pack_infos.push_back(pi);
        } else {
          ShufflePlanRecvPackInfo& pi = ret.recv_pack_infos.back();
          pi.size += 1;
        }
        recv_id_node_idx[id_node] += 1;
        last_buffer_idx += 1;
      }
    }
  }
  displayln_info(ssprintf("%s: send_pack_infos.size() = %10ld", fname, ret.send_pack_infos.size()));
  displayln_info(ssprintf("%s: recv_pack_infos.size() = %10ld", fname, ret.send_pack_infos.size()));
  displayln_info(ssprintf("%s: send_msg_infos.size()  = %10ld", fname, ret.send_msg_infos.size()));
  displayln_info(ssprintf("%s: recv_msg_infos.size()  = %10ld", fname, ret.recv_msg_infos.size()));
  return ret;
}

inline Cache<ShufflePlanKey,ShufflePlan>& get_shuffle_plan_cache()
{
  static Cache<ShufflePlanKey,ShufflePlan> cache(16);
  return cache;
}

inline const ShufflePlan& get_shuffle_plan(const ShufflePlanKey& spk)
{
  if (!get_shuffle_plan_cache().has(spk)) {
    get_shuffle_plan_cache()[spk] = make_shuffle_plan(spk);
  }
  return get_shuffle_plan_cache()[spk];
}

inline const ShufflePlan& get_shuffle_plan(const Coordinate& total_site, const Coordinate& new_size_node)
{
  ShufflePlanKey spk;
  spk.total_site = total_site;
  spk.new_size_node = new_size_node;
  return get_shuffle_plan(spk);
}

template <class M>
void shuffle_field(std::vector<Field<M> >& fs, const Field<M>& f, const Coordinate& new_size_node)
{
  fs.clear();
  const Geometry& geo = f.geo;
  if (geo.geon.size_node == new_size_node) {
    fs.resize(1);
    fs[0] = f;
  }
  TIMER_VERBOSE_FLOPS("shuffle_field");
  const ShufflePlan& sp = get_shuffle_plan(geo.total_site(), new_size_node);
  const long total_bytes = sp.total_send_size * geo.multiplicity * sizeof(M) * get_num_node();
  timer.flops += total_bytes;
  std::vector<M> send_buffer(sp.total_send_size * geo.multiplicity);
#pragma omp parallel for
  for (size_t i = 0; i < sp.send_pack_infos.size(); ++i) {
    const ShufflePlanSendPackInfo& pi = sp.send_pack_infos[i];
    memcpy(
        &send_buffer[pi.buffer_idx * geo.multiplicity],
        f.get_elems_const(pi.field_idx).data(),
        pi.size * geo.multiplicity * sizeof(M));
  }
  std::vector<M> recv_buffer(sp.total_recv_size * geo.multiplicity);
  {
    sync_node();
    TIMER_VERBOSE_FLOPS("shuffle_field-comm");
    timer.flops += total_bytes;
    std::vector<MPI_Request> send_reqs(sp.send_msg_infos.size());
    std::vector<MPI_Request> recv_reqs(sp.recv_msg_infos.size());
    {
      TIMER_VERBOSE("shuffle_field-comm-init");
      const int mpi_tag = 4;
      for (size_t i = 0; i < sp.send_msg_infos.size(); ++i) {
        const ShufflePlanMsgInfo& mi = sp.send_msg_infos[i];
        MPI_Isend(&send_buffer[mi.idx * geo.multiplicity], mi.size * geo.multiplicity * sizeof(M), MPI_BYTE, mi.id_node,
            mpi_tag, get_comm(), &send_reqs[i]);
      }
      for (size_t i = 0; i < sp.recv_msg_infos.size(); ++i) {
        const ShufflePlanMsgInfo& mi = sp.recv_msg_infos[i];
        MPI_Irecv(&recv_buffer[mi.idx * geo.multiplicity], mi.size * geo.multiplicity * sizeof(M), MPI_BYTE, mi.id_node,
            mpi_tag, get_comm(), &recv_reqs[i]);
      }
    }
    for (size_t i = 0; i < recv_reqs.size(); ++i) {
      MPI_Wait(&recv_reqs[i], MPI_STATUS_IGNORE);
    }
    for (size_t i = 0; i < send_reqs.size(); ++i) {
      MPI_Wait(&send_reqs[i], MPI_STATUS_IGNORE);
    }
    sync_node();
  }
  send_buffer.clear();
  const std::vector<Geometry> new_geos = make_dist_io_geos(geo.total_site(), geo.multiplicity, new_size_node);
  fs.resize(new_geos.size());
  for (size_t i = 0; i < fs.size(); ++i) {
    fs[i].init(new_geos[i]);
  }
#pragma omp parallel for
  for (size_t i = 0; i < sp.recv_pack_infos.size(); ++i) {
    const ShufflePlanRecvPackInfo& pi = sp.recv_pack_infos[i];
    qassert(0 <= pi.local_geos_idx && pi.local_geos_idx < fs.size());
    memcpy(
        fs[pi.local_geos_idx].get_elems(pi.field_idx).data(),
        &recv_buffer[pi.buffer_idx * geo.multiplicity],
        pi.size * geo.multiplicity * sizeof(M));
  }
}

template <class M>
long dist_write_field(const Field<M>& f,
    const Coordinate new_size_node,
    const std::string& path, const mode_t mode = default_dir_mode())
  // interface_function
{
  TIMER_VERBOSE_FLOPS("dist_write_field");
  std::vector<Field<M> > fs;
  shuffle_field(fs, f, new_size_node);
  long total_bytes = dist_write_fields(fs, product(new_size_node), path, mode);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
void shuffle_field_back(Field<M>& f, const std::vector<Field<M> >& fs, const Coordinate& new_size_node)
{
  const Geometry& geo = f.geo;
  if (geo.geon.size_node == new_size_node) {
    qassert(fs.size() == 1);
    f = fs[0];
  }
  TIMER_VERBOSE_FLOPS("shuffle_field_back");
  const ShufflePlan& sp = get_shuffle_plan(geo.total_site(), new_size_node);
  const long total_bytes = sp.total_send_size * geo.multiplicity * sizeof(M) * get_num_node();
  timer.flops += total_bytes;
  std::vector<M> recv_buffer(sp.total_recv_size * geo.multiplicity);
#pragma omp parallel for
  for (size_t i = 0; i < sp.recv_pack_infos.size(); ++i) {
    const ShufflePlanRecvPackInfo& pi = sp.recv_pack_infos[i];
    memcpy(
        &recv_buffer[pi.buffer_idx * geo.multiplicity],
        fs[pi.local_geos_idx].get_elems_const(pi.field_idx).data(),
        pi.size * geo.multiplicity * sizeof(M));
  }
  std::vector<M> send_buffer(sp.total_send_size * geo.multiplicity);
  {
    sync_node();
    TIMER_VERBOSE_FLOPS("shuffle_field-comm");
    timer.flops += total_bytes;
    std::vector<MPI_Request> send_reqs(sp.send_msg_infos.size());
    std::vector<MPI_Request> recv_reqs(sp.recv_msg_infos.size());
    {
      TIMER("shuffle_field-comm-init");
      const int mpi_tag = 5;
      for (size_t i = 0; i < sp.recv_msg_infos.size(); ++i) {
        const ShufflePlanMsgInfo& mi = sp.recv_msg_infos[i];
        MPI_Isend(&recv_buffer[mi.idx * geo.multiplicity], mi.size * geo.multiplicity * sizeof(M), MPI_BYTE, mi.id_node,
            mpi_tag, get_comm(), &recv_reqs[i]);
      }
      for (size_t i = 0; i < sp.send_msg_infos.size(); ++i) {
        const ShufflePlanMsgInfo& mi = sp.send_msg_infos[i];
        MPI_Irecv(&send_buffer[mi.idx * geo.multiplicity], mi.size * geo.multiplicity * sizeof(M), MPI_BYTE, mi.id_node,
            mpi_tag, get_comm(), &send_reqs[i]);
      }
    }
    for (size_t i = 0; i < recv_reqs.size(); ++i) {
      MPI_Wait(&recv_reqs[i], MPI_STATUS_IGNORE);
    }
    for (size_t i = 0; i < send_reqs.size(); ++i) {
      MPI_Wait(&send_reqs[i], MPI_STATUS_IGNORE);
    }
    sync_node();
  }
  recv_buffer.clear();
#pragma omp parallel for
  for (size_t i = 0; i < sp.send_pack_infos.size(); ++i) {
    const ShufflePlanSendPackInfo& pi = sp.send_pack_infos[i];
    memcpy(
        f.get_elems(pi.field_idx).data(),
        &send_buffer[pi.buffer_idx * geo.multiplicity],
        pi.size * geo.multiplicity * sizeof(M));
  }
}

template <class M>
long dist_read_field(Field<M>& f, const std::string& path)
  // interface_function
{
  TIMER_VERBOSE_FLOPS("dist_read_field");
  Geometry geo;
  std::vector<Field<M> > fs;
  Coordinate new_size_node;
  const long total_bytes = dist_read_fields(fs, geo, new_size_node, path);
  f.init(geo);
  qassert(f.geo == geo);
  shuffle_field_back(f, fs, new_size_node);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M, class N>
void convert_field_float_from_double(Field<N>& ff, const Field<M>&f)
  // interface_function
{
  TIMER("convert_field_float_from_double");
  qassert(f.geo.is_only_local());
  qassert(sizeof(M) % sizeof(double) == 0);
  qassert(sizeof(N) % sizeof(float) == 0);
  qassert(f.geo.multiplicity * sizeof(M) / 2 % sizeof(N) == 0);
  const int multiplicity = f.geo.multiplicity * sizeof(M) / 2 / sizeof(N);
  const Geometry geo = geo_remult(f.geo, multiplicity);
  ff.init(geo);
  const Vector<M> fdata = get_data(f);
  const Vector<double> fd((double*)fdata.data(), fdata.data_size() / sizeof(double));
  Vector<N> ffdata = get_data(ff);
  Vector<float> ffd((float*)ffdata.data(), ffdata.data_size() / sizeof(float));
  qassert(ffd.size() == fd.size());
  for (long i = 0; i < ffd.size(); ++i) {
    ffd[i] = fd[i];
  }
}

template <class M, class N>
void convert_field_double_from_float(Field<N>& ff, const Field<M>&f)
  // interface_function
{
  TIMER("convert_field_double_from_float");
  qassert(f.geo.is_only_local());
  qassert(sizeof(M) % sizeof(float) == 0);
  qassert(sizeof(N) % sizeof(double) == 0);
  qassert(f.geo.multiplicity * sizeof(M) * 2 % sizeof(N) == 0);
  const int multiplicity = f.geo.multiplicity * sizeof(M) * 2 / sizeof(N);
  const Geometry geo = geo_remult(f.geo, multiplicity);
  ff.init(geo);
  const Vector<M> fdata = get_data(f);
  const Vector<float> fd((float*)fdata.data(), fdata.data_size() / sizeof(float));
  Vector<N> ffdata = get_data(ff);
  Vector<double> ffd((double*)ffdata.data(), ffdata.data_size() / sizeof(double));
  qassert(ffd.size() == fd.size());
  for (long i = 0; i < ffd.size(); ++i) {
    ffd[i] = fd[i];
  }
}

template <class M>
long dist_write_field_float_from_double(const Field<M>& f,
    const std::string& path, const mode_t mode = default_dir_mode())
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
    const std::string& path, const mode_t mode = default_dir_mode())
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
  to_from_big_endian_32(get_data(ff));
  convert_field_double_from_float(f, ff);
  timer.flops += total_bytes;
  return total_bytes;
}

QLAT_END_NAMESPACE
