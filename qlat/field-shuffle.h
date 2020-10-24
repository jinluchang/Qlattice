// vim: set ts=2 sw=2 expandtab:

#pragma once

#include <qlat/config.h>
#include <qlat/utils.h>

QLAT_START_NAMESPACE

struct ShufflePlanRecvPackInfo {
  int local_geos_idx;  // idx of the field that the data belong to
  long field_idx;      // idx of the data in the field
  long buffer_idx;     // idx of the data in the buffer
  long size;           // number of data site for this pack of data
};

struct ShufflePlanSendPackInfo {
  long field_idx;   // idx of the data in the field
  long buffer_idx;  // idx of the data in the buffer
  long size;        // number of data site for this pack of data
};

struct ShufflePlanMsgInfo {
  int id_node;  // index for the target send/recv node
  long idx;     // idx of the starting location in send/recv buffer
  long size;    // number of data site for this msg
};

struct ShufflePlan {
  Geometry geo_send;                // geo of the send field
  std::vector<Geometry> geos_recv;  // geos of the recv fields
  long total_send_size;             // total send buffer size
  std::vector<ShufflePlanMsgInfo>
      send_msg_infos;  // corresponds to every sent msg
  std::vector<ShufflePlanSendPackInfo>
      send_pack_infos;  // corresponds to how to create send buffer from local
                        // field
  long total_recv_size;
  std::vector<ShufflePlanMsgInfo>
      recv_msg_infos;  // corresponds to every recv msg
  std::vector<ShufflePlanRecvPackInfo>
      recv_pack_infos;  // corresponds to how to copy recv buffer to new local
                        // fields
};

template <class M>
void shuffle_field(std::vector<Field<M> >& fs, const Field<M>& f,
                   const ShufflePlan& sp)
{
  TIMER_FLOPS("shuffle_field(fs,f,sp)");
  const Geometry& geo = f.geo;
  qassert(sp.geo_send == geo_reform(geo, 1, 0));
  clear(fs);
  sync_node();
  const long total_bytes =
      sp.total_send_size * geo.multiplicity * sizeof(M) * get_num_node();
  timer.flops += total_bytes;
  std::vector<M> send_buffer(sp.total_send_size * geo.multiplicity);
#pragma omp parallel for
  for (size_t i = 0; i < sp.send_pack_infos.size(); ++i) {
    const ShufflePlanSendPackInfo& pi = sp.send_pack_infos[i];
    memcpy(&send_buffer[pi.buffer_idx * geo.multiplicity],
           f.get_elems_const(pi.field_idx).data(),
           pi.size * geo.multiplicity * sizeof(M));
  }
  std::vector<M> recv_buffer(sp.total_recv_size * geo.multiplicity);
  {
    sync_node();
    TIMER_FLOPS("shuffle_field-comm");
    timer.flops += total_bytes;
    std::vector<MPI_Request> send_reqs(sp.send_msg_infos.size());
    std::vector<MPI_Request> recv_reqs(sp.recv_msg_infos.size());
    {
      TIMER("shuffle_field-comm-init");
      const int mpi_tag = 4;
      for (size_t i = 0; i < sp.send_msg_infos.size(); ++i) {
        const ShufflePlanMsgInfo& mi = sp.send_msg_infos[i];
        MPI_Isend(&send_buffer[mi.idx * geo.multiplicity],
                  mi.size * geo.multiplicity * sizeof(M), MPI_BYTE, mi.id_node,
                  mpi_tag, get_comm(), &send_reqs[i]);
      }
      for (size_t i = 0; i < sp.recv_msg_infos.size(); ++i) {
        const ShufflePlanMsgInfo& mi = sp.recv_msg_infos[i];
        MPI_Irecv(&recv_buffer[mi.idx * geo.multiplicity],
                  mi.size * geo.multiplicity * sizeof(M), MPI_BYTE, mi.id_node,
                  mpi_tag, get_comm(), &recv_reqs[i]);
      }
    }
    MPI_Waitall(recv_reqs.size(), recv_reqs.data(), MPI_STATUS_IGNORE);
    MPI_Waitall(send_reqs.size(), send_reqs.data(), MPI_STATUS_IGNORE);
    sync_node();
  }
  clear(send_buffer);
  std::vector<Geometry> geos_recv = sp.geos_recv;
  fs.resize(geos_recv.size());
  for (size_t i = 0; i < fs.size(); ++i) {
    geos_recv[i].remult(geo.multiplicity);
    fs[i].init(geos_recv[i]);
  }
#pragma omp parallel for
  for (size_t i = 0; i < sp.recv_pack_infos.size(); ++i) {
    const ShufflePlanRecvPackInfo& pi = sp.recv_pack_infos[i];
    qassert(0 <= pi.local_geos_idx && pi.local_geos_idx < (int)fs.size());
    std::memcpy(fs[pi.local_geos_idx].get_elems(pi.field_idx).data(),
                &recv_buffer[pi.buffer_idx * geo.multiplicity],
                pi.size * geo.multiplicity * sizeof(M));
  }
}

template <class M>
void shuffle_field_back(Field<M>& f, const std::vector<Field<M> >& fs,
                        const ShufflePlan& sp)
{
  TIMER_FLOPS("shuffle_field_back(f,fs,sp)");
  const Geometry& geo = f.geo;
  sync_node();
  const long total_bytes =
      sp.total_send_size * geo.multiplicity * sizeof(M) * get_num_node();
  timer.flops += total_bytes;
  std::vector<M> recv_buffer(sp.total_recv_size * geo.multiplicity);
#pragma omp parallel for
  for (size_t i = 0; i < sp.recv_pack_infos.size(); ++i) {
    const ShufflePlanRecvPackInfo& pi = sp.recv_pack_infos[i];
    memcpy(&recv_buffer[pi.buffer_idx * geo.multiplicity],
           fs[pi.local_geos_idx].get_elems_const(pi.field_idx).data(),
           pi.size * geo.multiplicity * sizeof(M));
  }
  std::vector<M> send_buffer(sp.total_send_size * geo.multiplicity);
  {
    sync_node();
    TIMER_FLOPS("shuffle_field-comm");
    timer.flops += total_bytes;
    std::vector<MPI_Request> send_reqs(sp.send_msg_infos.size());
    std::vector<MPI_Request> recv_reqs(sp.recv_msg_infos.size());
    {
      TIMER("shuffle_field-comm-init");
      const int mpi_tag = 5;
      for (size_t i = 0; i < sp.recv_msg_infos.size(); ++i) {
        const ShufflePlanMsgInfo& mi = sp.recv_msg_infos[i];
        MPI_Isend(&recv_buffer[mi.idx * geo.multiplicity],
                  mi.size * geo.multiplicity * sizeof(M), MPI_BYTE, mi.id_node,
                  mpi_tag, get_comm(), &recv_reqs[i]);
      }
      for (size_t i = 0; i < sp.send_msg_infos.size(); ++i) {
        const ShufflePlanMsgInfo& mi = sp.send_msg_infos[i];
        MPI_Irecv(&send_buffer[mi.idx * geo.multiplicity],
                  mi.size * geo.multiplicity * sizeof(M), MPI_BYTE, mi.id_node,
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
  clear(recv_buffer);
#pragma omp parallel for
  for (size_t i = 0; i < sp.send_pack_infos.size(); ++i) {
    const ShufflePlanSendPackInfo& pi = sp.send_pack_infos[i];
    memcpy(f.get_elems(pi.field_idx).data(),
           &send_buffer[pi.buffer_idx * geo.multiplicity],
           pi.size * geo.multiplicity * sizeof(M));
  }
}

// field dist write shuffle

inline long& get_shuffle_max_msg_size()
// qlat parameter
{
  static long size = 16 * 1024;
  return size;
}

inline long& get_shuffle_max_pack_size()
// qlat parameter
{
  static long size = 1024;
  return size;
}

inline std::vector<GeometryNode> make_dist_io_geons(
    const Coordinate& new_size_node)
{
  TIMER("make_dist_io_geons");
  const int new_num_node = product(new_size_node);
  const int num_node = get_num_node();
  std::vector<GeometryNode> ret;
  const int min_size_chunk = new_num_node / num_node;
  const int remain = new_num_node % num_node;
  const int id_node_in_shuffle =
      get_num_node() == new_num_node
          ? get_id_node()
          : get_id_node_in_shuffle(get_id_node(), new_num_node, num_node);
  const int size_chunk =
      id_node_in_shuffle < remain ? min_size_chunk + 1 : min_size_chunk;
  const int chunk_start =
      id_node_in_shuffle * min_size_chunk +
      (id_node_in_shuffle < remain ? id_node_in_shuffle : remain);
  const int chunk_end = std::min(new_num_node, chunk_start + size_chunk);
  for (int new_id_node = chunk_start; new_id_node < chunk_end; ++new_id_node) {
    GeometryNode geon;
    geon.initialized = true;
    geon.num_node = new_num_node;
    geon.id_node = new_id_node;
    geon.size_node = new_size_node;
    geon.coor_node = coordinate_from_index(new_id_node, new_size_node);
    ret.push_back(geon);
  }
  return ret;
}

inline std::vector<Geometry> make_dist_io_geos(const Coordinate& total_site,
                                               const int multiplicity,
                                               const Coordinate& new_size_node)
{
  TIMER("make_dist_io_geos");
  const std::vector<GeometryNode> geons = make_dist_io_geons(new_size_node);
  std::vector<Geometry> ret;
  const Coordinate new_node_site = total_site / new_size_node;
  for (int i = 0; i < (int)geons.size(); ++i) {
    Geometry geo_recv;
    geo_recv.init(geons[i], new_node_site, multiplicity);
    ret.push_back(geo_recv);
  }
  return ret;
}

inline int get_id_node_in_shuffle_from_new_id_node(const int new_id_node,
                                                   const int new_num_node,
                                                   const int num_node)
{
  const int min_size_chunk = new_num_node / num_node;
  const int remain = new_num_node % num_node;
  const int limit = remain * min_size_chunk + remain;
  int id_node_in_shuffle;
  if (new_id_node <= limit) {
    id_node_in_shuffle = new_id_node / (min_size_chunk + 1);
  } else {
    id_node_in_shuffle = (new_id_node - limit) / min_size_chunk + remain;
  }
  return id_node_in_shuffle;
}

inline int get_id_node_from_new_id_node(const int new_id_node,
                                        const int new_num_node,
                                        const int num_node)
{
  return get_id_node_from_id_node_in_shuffle(
      get_id_node_in_shuffle_from_new_id_node(new_id_node, new_num_node,
                                              num_node),
      new_num_node, num_node);
}

struct ShufflePlanKey {
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

inline ShufflePlan make_shuffle_plan(const ShufflePlanKey& spk)
// assume the order is preserved in transfer.
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
  std::vector<Geometry> geos_recv =
      make_dist_io_geos(geo.total_site(), geo.multiplicity, new_size_node);
  // to be shared
  ShufflePlan ret;
  // geo_send
  ret.geo_send = geo;
  // geos_recv
  ret.geos_recv = geos_recv;
  // total_send_size
  ret.total_send_size = geo.local_volume();
  // send_id_node_size
  // send_new_id_node_size
  std::map<int, long> send_id_node_size;
  std::map<int, long> send_new_id_node_size;
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    // custom
    const Coordinate new_coor_node = xg / new_node_site;
    const int new_id_node = index_from_coordinate(new_coor_node, new_size_node);
    // custom
    const int id_node_in_shuffle = get_id_node_in_shuffle_from_new_id_node(
        new_id_node, new_num_node, num_node);
    //
    send_id_node_size[id_node_in_shuffle] += 1;
    send_new_id_node_size[new_id_node] += 1;
  }
  {
    long count = 0;
    for (std::map<int, long>::const_iterator it = send_id_node_size.begin();
         it != send_id_node_size.end(); ++it) {
      const int id_node_in_shuffle = it->first;
      const long node_size = it->second;
      long node_size_remain = node_size;
      while (node_size_remain > 0) {
        ShufflePlanMsgInfo mi;
        mi.id_node = get_id_node_from_id_node_in_shuffle(
            id_node_in_shuffle, new_num_node, num_node);
        mi.idx = count;
        mi.size = std::min(node_size_remain, get_shuffle_max_msg_size());
        ret.send_msg_infos.push_back(mi);
        node_size_remain -= mi.size;
        count += mi.size;
      }
    }
    qassert(count == geo.local_volume());
    qassert(count == ret.total_send_size);
  }
  // send_new_id_node_idx
  std::map<int, long> send_new_id_node_idx;
  {
    long count = 0;
    for (std::map<int, long>::const_iterator it = send_new_id_node_size.begin();
         it != send_new_id_node_size.end(); ++it) {
      const int new_id_node = it->first;
      const long node_size = it->second;
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
      // custom
      const Coordinate new_coor_node = xg / new_node_site;
      const int new_id_node =
          index_from_coordinate(new_coor_node, new_size_node);
      //
      const long buffer_idx = send_new_id_node_idx[new_id_node];
      if (buffer_idx == last_buffer_idx and
          ret.send_pack_infos.back().size < get_shuffle_max_pack_size()) {
        ret.send_pack_infos.back().size += 1;
      } else {
        last_buffer_idx = buffer_idx;
        ShufflePlanSendPackInfo pi;
        pi.field_idx = index;
        pi.buffer_idx = buffer_idx;
        pi.size = 1;
        ret.send_pack_infos.push_back(pi);
      }
      send_new_id_node_idx[new_id_node] += 1;
      last_buffer_idx += 1;
    }
  }
  // total_recv_size
  ret.total_recv_size = 0;
  for (size_t i = 0; i < ret.geos_recv.size(); ++i) {
    const Geometry& geos_recv = ret.geos_recv[i];
    ret.total_recv_size += geos_recv.local_volume();
  }
  // recv_id_node_size
  std::map<int, long> recv_id_node_size;
  for (size_t i = 0; i < ret.geos_recv.size(); ++i) {
    const Geometry& geo_recv = ret.geos_recv[i];
    for (long index = 0; index < geo_recv.local_volume(); ++index) {
      const Coordinate xl = geo_recv.coordinate_from_index(index);
      const Coordinate xg = geo_recv.coordinate_g_from_l(xl);
      const Coordinate coor_node = xg / geo.node_site;
      // custom
      const long id_node = index_from_coordinate(coor_node, geo.geon.size_node);
      //
      recv_id_node_size[id_node] += 1;
    }
  }
  {
    long count = 0;
    for (std::map<int, long>::const_iterator it = recv_id_node_size.begin();
         it != recv_id_node_size.end(); ++it) {
      const int id_node = it->first;
      const long node_size = it->second;
      long node_size_remain = node_size;
      while (node_size_remain > 0) {
        ShufflePlanMsgInfo mi;
        mi.id_node = id_node;
        mi.idx = count;
        mi.size = std::min(node_size_remain, get_shuffle_max_msg_size());
        ret.recv_msg_infos.push_back(mi);
        node_size_remain -= mi.size;
        count += mi.size;
      }
    }
    qassert(count == ret.total_recv_size);
  }
  // recv_id_node_idx
  std::map<int, long> recv_id_node_idx;
  {
    long count = 0;
    for (std::map<int, long>::const_iterator it = recv_id_node_size.begin();
         it != recv_id_node_size.end(); ++it) {
      const int id_node = it->first;
      const long node_size = it->second;
      recv_id_node_idx[id_node] = count;
      count += node_size;
    }
    qassert(count == ret.total_recv_size);
  }
  // recv_pack_infos
  {
    for (size_t i = 0; i < ret.geos_recv.size(); ++i) {
      const Geometry& geo_recv = ret.geos_recv[i];
      long last_buffer_idx = -1;
      for (long index = 0; index < geo_recv.local_volume(); ++index) {
        const Coordinate xl = geo_recv.coordinate_from_index(index);
        const Coordinate xg = geo_recv.coordinate_g_from_l(xl);
        // custom
        const Coordinate coor_node = xg / geo.node_site;
        const int id_node =
            index_from_coordinate(coor_node, geo.geon.size_node);
        //
        const long buffer_idx = recv_id_node_idx[id_node];
        if (buffer_idx == last_buffer_idx and
            ret.recv_pack_infos.back().size < get_shuffle_max_pack_size()) {
          ret.recv_pack_infos.back().size += 1;
        } else {
          last_buffer_idx = buffer_idx;
          ShufflePlanRecvPackInfo pi;
          pi.local_geos_idx = i;
          pi.field_idx = index;
          pi.buffer_idx = buffer_idx;
          pi.size = 1;
          ret.recv_pack_infos.push_back(pi);
        }
        recv_id_node_idx[id_node] += 1;
        last_buffer_idx += 1;
      }
    }
  }
  long num_send_packs = ret.send_pack_infos.size();
  long num_recv_packs = ret.recv_pack_infos.size();
  long num_send_msgs = ret.send_msg_infos.size();
  long num_recv_msgs = ret.recv_msg_infos.size();
  displayln_info(ssprintf("%s::%s: num_send_packs = %10ld", cname().c_str(),
                          fname, num_send_packs));
  displayln_info(ssprintf("%s::%s: num_recv_packs = %10ld", cname().c_str(),
                          fname, num_recv_packs));
  displayln_info(ssprintf("%s::%s: num_send_msgs  = %10ld", cname().c_str(),
                          fname, num_send_msgs));
  displayln_info(ssprintf("%s::%s: num_recv_msgs  = %10ld", cname().c_str(),
                          fname, num_recv_msgs));
  glb_sum(num_send_packs);
  glb_sum(num_recv_packs);
  glb_sum(num_send_msgs);
  glb_sum(num_recv_msgs);
  displayln_info(ssprintf("%s::%s: total num_send_packs = %10ld",
                          cname().c_str(), fname, num_send_packs));
  displayln_info(ssprintf("%s::%s: total num_recv_packs = %10ld",
                          cname().c_str(), fname, num_recv_packs));
  displayln_info(ssprintf("%s::%s: total num_send_msgs  = %10ld",
                          cname().c_str(), fname, num_send_msgs));
  displayln_info(ssprintf("%s::%s: total num_recv_msgs  = %10ld",
                          cname().c_str(), fname, num_recv_msgs));
  return ret;
}

inline Cache<ShufflePlanKey, ShufflePlan>& get_shuffle_plan_cache()
{
  static Cache<ShufflePlanKey, ShufflePlan> cache("ShufflePlanCache", 16);
  return cache;
}

inline const ShufflePlan& get_shuffle_plan(const ShufflePlanKey& spk)
{
  if (!get_shuffle_plan_cache().has(spk)) {
    get_shuffle_plan_cache()[spk] = make_shuffle_plan(spk);
  }
  return get_shuffle_plan_cache()[spk];
}

inline const ShufflePlan& get_shuffle_plan(const Coordinate& total_site,
                                           const Coordinate& new_size_node)
{
  ShufflePlanKey spk;
  spk.total_site = total_site;
  spk.new_size_node = new_size_node;
  return get_shuffle_plan(spk);
}

template <class M>
void shuffle_field(std::vector<Field<M> >& fs, const Field<M>& f,
                   const Coordinate& new_size_node)
{
  clear(fs);
  const Geometry& geo = f.geo;
  if (geo.geon.size_node == new_size_node) {
    fs.resize(1);
    fs[0] = f;
    return;
  }
  sync_node();
  TIMER_VERBOSE_FLOPS("shuffle_field");
  displayln_info(
      fname +
      ssprintf(
          ": %s -> %s (total_site: %s ; site_size: %d ; total_size: %.3lf GB)",
          show(geo.geon.size_node).c_str(), show(new_size_node).c_str(),
          show(geo.total_site()).c_str(), geo.multiplicity * (int)sizeof(M),
          (double)(product(geo.total_site()) * geo.multiplicity * sizeof(M) *
                   std::pow(0.5, 30))));
  const ShufflePlan& sp = get_shuffle_plan(geo.total_site(), new_size_node);
  shuffle_field(fs, f, sp);
  const long total_bytes =
      geo.local_volume() * geo.multiplicity * sizeof(M) * get_num_node();
  timer.flops += total_bytes;
}

template <class M>
void shuffle_field_back(Field<M>& f, const std::vector<Field<M> >& fs,
                        const Coordinate& new_size_node)
{
  const Geometry& geo = f.geo;
  if (geo.geon.size_node == new_size_node) {
    qassert(fs.size() == 1);
    f = fs[0];
    return;
  }
  sync_node();
  TIMER_VERBOSE_FLOPS("shuffle_field_back");
  displayln_info(
      fname +
      ssprintf(
          ": %s -> %s (total_site: %s ; site_size: %d ; total_size: %.3lf GB)",
          show(new_size_node).c_str(), show(geo.geon.size_node).c_str(),
          show(geo.total_site()).c_str(), geo.multiplicity * (int)sizeof(M),
          (double)(product(geo.total_site()) * geo.multiplicity * sizeof(M) *
                   std::pow(0.5, 30))));
  const ShufflePlan& sp = get_shuffle_plan(geo.total_site(), new_size_node);
  shuffle_field_back(f, fs, sp);
  const long total_bytes =
      geo.local_volume() * geo.multiplicity * sizeof(M) * get_num_node();
  timer.flops += total_bytes;
}

// fft shuffle

inline long index_from_coordinate_perp_dir(const Coordinate& coor,
                                           const Coordinate& size,
                                           const int dir)
{
  std::array<int, 3> dirs;
  int cur_dir = 0;
  for (int i = 0; i < 4; ++i) {
    if (i != dir) {
      dirs[cur_dir] = i;
      cur_dir += 1;
    }
  }
  return coor[dirs[0]] +
         size[dirs[0]] * (coor[dirs[1]] + size[dirs[1]] * coor[dirs[2]]);
}

inline Coordinate coordinate_from_index_perp_dir(long index,
                                                 const Coordinate& size,
                                                 const int dir,
                                                 const int coor_dir)
{
  std::array<int, 3> dirs;
  int cur_dir = 0;
  for (int i = 0; i < 4; ++i) {
    if (i != dir) {
      dirs[cur_dir] = i;
      cur_dir += 1;
    }
  }
  Coordinate x;
  x[dirs[0]] = index % size[dirs[0]];
  index /= size[dirs[0]];
  x[dirs[1]] = index % size[dirs[1]];
  index /= size[dirs[1]];
  x[dirs[2]] = index % size[dirs[2]];
  x[dir] = coor_dir;
  return x;
}

inline long get_id_node_fft(const Coordinate& xl, const Coordinate& node_site,
                            const Coordinate& coor_node,
                            const Coordinate& size_node, const int dir)
{
  const long vol_perp_dir = product(node_site) / node_site[dir];
  const long idx = index_from_coordinate_perp_dir(xl, node_site, dir);
  Coordinate new_coor_node = coor_node;
  new_coor_node[dir] = find_worker(idx, vol_perp_dir, size_node[dir]);
  qassert(new_coor_node[dir] < size_node[dir]);
  return index_from_coordinate(new_coor_node, size_node);
}

inline ShufflePlan make_shuffle_plan_fft(const Coordinate& total_site,
                                         const int dir)
{
  TIMER_VERBOSE("make_shuffle_plan_fft");
  Geometry geo;
  geo.init(total_site, 1);
  // vol_perp_dir -> vpd
  const long vol_perp_dir =
      product(geo.node_site) / geo.node_site[dir];  // local volume perp to dir
  const Coordinate size_node = geo.geon.size_node;
  const Coordinate coor_node = geo.geon.coor_node;
  const long num_node_dir = size_node[dir];
  const long id_node_dir = coor_node[dir];
  long vpd_start;
  long vpd_size;  // share of the ``vol_perp_dir'' on this node
  split_work(vpd_start, vpd_size, vol_perp_dir, num_node_dir, id_node_dir);
  std::vector<Geometry> geos_recv(total_site[dir]);
  const Coordinate new_size_node(1, 1, total_site[dir], geo.geon.num_node);
  for (int i = 0; i < total_site[dir]; ++i) {
    const Coordinate new_coor_node(0, 0, i, geo.geon.id_node);
    const GeometryNode geon(index_from_coordinate(new_coor_node, new_size_node),
                            new_size_node);
    geos_recv[i].init(geon, Coordinate(vpd_size, 1, 1, 1), 1);
  }
  // to be shared
  ShufflePlan ret;
  // geo_send
  ret.geo_send = geo;
  // geos_recv
  ret.geos_recv = geos_recv;
  // total_send_size
  ret.total_send_size = geo.local_volume();
  // send_id_node_size
  // send_new_id_node_size
  std::map<int, long> send_id_node_size;
  std::map<int, long> send_new_id_node_size;
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    // custom
    const int id_node =
        get_id_node_fft(xl, geo.node_site, coor_node, size_node, dir);
    // custom
    const Coordinate new_coor_node(1, 1, xg[dir], id_node);
    const int new_id_node = index_from_coordinate(new_coor_node, new_size_node);
    //
    send_id_node_size[id_node] += 1;
    send_new_id_node_size[new_id_node] += 1;
  }
  {
    long count = 0;
    for (std::map<int, long>::const_iterator it = send_id_node_size.begin();
         it != send_id_node_size.end(); ++it) {
      const int id_node = it->first;
      const long node_size = it->second;
      long node_size_remain = node_size;
      while (node_size_remain > 0) {
        ShufflePlanMsgInfo mi;
        mi.id_node = id_node;
        mi.idx = count;
        mi.size = std::min(node_size_remain, get_shuffle_max_msg_size());
        ret.send_msg_infos.push_back(mi);
        node_size_remain -= mi.size;
        count += mi.size;
      }
    }
    qassert(count == geo.local_volume());
    qassert(count == ret.total_send_size);
  }
  // send_new_id_node_idx
  std::map<int, long> send_new_id_node_idx;
  {
    long count = 0;
    for (std::map<int, long>::const_iterator it = send_new_id_node_size.begin();
         it != send_new_id_node_size.end(); ++it) {
      const int new_id_node = it->first;
      const long node_size = it->second;
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
      // custom
      const int id_node =
          get_id_node_fft(xl, geo.node_site, coor_node, size_node, dir);
      // custom
      const Coordinate new_coor_node(1, 1, xg[dir], id_node);
      const int new_id_node =
          index_from_coordinate(new_coor_node, new_size_node);
      //
      const long buffer_idx = send_new_id_node_idx[new_id_node];
      if (buffer_idx == last_buffer_idx and
          ret.send_pack_infos.back().size < get_shuffle_max_pack_size()) {
        ret.send_pack_infos.back().size += 1;
      } else {
        last_buffer_idx = buffer_idx;
        ShufflePlanSendPackInfo pi;
        pi.field_idx = index;
        pi.buffer_idx = buffer_idx;
        pi.size = 1;
        ret.send_pack_infos.push_back(pi);
      }
      send_new_id_node_idx[new_id_node] += 1;
      last_buffer_idx += 1;
    }
  }
  // total_recv_size
  ret.total_recv_size = 0;
  for (size_t i = 0; i < ret.geos_recv.size(); ++i) {
    const Geometry& geo_recv = ret.geos_recv[i];
    ret.total_recv_size += geo_recv.local_volume();
  }
  // recv_id_node_size
  std::map<int, long> recv_id_node_size;
  for (size_t i = 0; i < ret.geos_recv.size(); ++i) {
    const Geometry& geo_recv = ret.geos_recv[i];
    for (long index = 0; index < geo_recv.local_volume(); ++index) {
      const Coordinate xl = geo_recv.coordinate_from_index(index);
      const Coordinate xg = geo_recv.coordinate_g_from_l(xl);
      // custom
      Coordinate coor_node_orig = coordinate_from_index(xg[3], size_node);
      coor_node_orig[dir] = xg[2] / geo.node_site[dir];
      const long id_node = index_from_coordinate(coor_node_orig, size_node);
      //
      recv_id_node_size[id_node] += 1;
    }
  }
  {
    long count = 0;
    for (std::map<int, long>::const_iterator it = recv_id_node_size.begin();
         it != recv_id_node_size.end(); ++it) {
      const int id_node = it->first;
      const long node_size = it->second;
      long node_size_remain = node_size;
      while (node_size_remain > 0) {
        ShufflePlanMsgInfo mi;
        mi.id_node = id_node;
        mi.idx = count;
        mi.size = std::min(node_size_remain, get_shuffle_max_msg_size());
        ret.recv_msg_infos.push_back(mi);
        node_size_remain -= mi.size;
        count += mi.size;
      }
    }
    if (count != ret.total_recv_size) {
      displayln(fname + ssprintf(": count = %ld", count));
      displayln(fname +
                ssprintf(": ret.total_recv_size = %ld", ret.total_recv_size));
    }
    qassert(count == ret.total_recv_size);
  }
  // recv_id_node_idx
  std::map<int, long> recv_id_node_idx;
  {
    long count = 0;
    for (std::map<int, long>::const_iterator it = recv_id_node_size.begin();
         it != recv_id_node_size.end(); ++it) {
      const int id_node = it->first;
      const long node_size = it->second;
      recv_id_node_idx[id_node] = count;
      count += node_size;
    }
    qassert(count == ret.total_recv_size);
  }
  // recv_pack_infos
  {
    for (size_t i = 0; i < ret.geos_recv.size(); ++i) {
      const Geometry& geo_recv = ret.geos_recv[i];
      long last_buffer_idx = -1;
      for (long index = 0; index < geo_recv.local_volume(); ++index) {
        const Coordinate xl = geo_recv.coordinate_from_index(index);
        const Coordinate xg = geo_recv.coordinate_g_from_l(xl);
        // custom
        Coordinate coor_node_orig = coordinate_from_index(xg[3], size_node);
        coor_node_orig[dir] = xg[2] / geo.node_site[dir];
        const long id_node = index_from_coordinate(coor_node_orig, size_node);
        //
        const long buffer_idx = recv_id_node_idx[id_node];
        if (buffer_idx == last_buffer_idx and
            ret.recv_pack_infos.back().size < get_shuffle_max_pack_size()) {
          ret.recv_pack_infos.back().size += 1;
        } else {
          last_buffer_idx = buffer_idx;
          ShufflePlanRecvPackInfo pi;
          pi.local_geos_idx = i;
          pi.field_idx = index;
          pi.buffer_idx = buffer_idx;
          pi.size = 1;
          ret.recv_pack_infos.push_back(pi);
        }
        recv_id_node_idx[id_node] += 1;
        last_buffer_idx += 1;
      }
    }
  }
  long num_send_packs = ret.send_pack_infos.size();
  long num_recv_packs = ret.recv_pack_infos.size();
  long num_send_msgs = ret.send_msg_infos.size();
  long num_recv_msgs = ret.recv_msg_infos.size();
  displayln_info(ssprintf("%s::%s: num_send_packs = %10ld", cname().c_str(),
                          fname, num_send_packs));
  displayln_info(ssprintf("%s::%s: num_recv_packs = %10ld", cname().c_str(),
                          fname, num_recv_packs));
  displayln_info(ssprintf("%s::%s: num_send_msgs  = %10ld", cname().c_str(),
                          fname, num_send_msgs));
  displayln_info(ssprintf("%s::%s: num_recv_msgs  = %10ld", cname().c_str(),
                          fname, num_recv_msgs));
  glb_sum(num_send_packs);
  glb_sum(num_recv_packs);
  glb_sum(num_send_msgs);
  glb_sum(num_recv_msgs);
  displayln_info(ssprintf("%s::%s: total num_send_packs = %10ld",
                          cname().c_str(), fname, num_send_packs));
  displayln_info(ssprintf("%s::%s: total num_recv_packs = %10ld",
                          cname().c_str(), fname, num_recv_packs));
  displayln_info(ssprintf("%s::%s: total num_send_msgs  = %10ld",
                          cname().c_str(), fname, num_send_msgs));
  displayln_info(ssprintf("%s::%s: total num_recv_msgs  = %10ld",
                          cname().c_str(), fname, num_recv_msgs));
  return ret;

  // reflection shuffle
}

QLAT_END_NAMESPACE
