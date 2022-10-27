// vim: set ts=2 sw=2 expandtab:

#pragma once

#include <qlat-utils/cache.h>
#include <qlat/setup.h>
#include <qlat/selected-field.h>


namespace qlat
{  //

struct API ShufflePlanMsgInfo {
  int id_node;  // index for the target send/recv node
  long idx;     // idx of the starting location in send/recv buffer
  long size;    // number of data site for this msg
};

struct API ShuffleCommPlan {
  long global_comm_size;  // global comm data size for Timer flops
  long total_send_size;   // total send buffer size
  std::vector<ShufflePlanMsgInfo>
      send_msg_infos;    // corresponds to every sent msg
  long total_recv_size;  // total recv buffer size
  std::vector<ShufflePlanMsgInfo>
      recv_msg_infos;  // corresponds to every recv msg
  //
  ShuffleCommPlan()
  {
    global_comm_size = 0;
    total_send_size = 0;
    total_recv_size = 0;
  }
};

struct API ShufflePlanRecvPackInfo {
  int local_geos_idx;  // idx of the field that the data belong to
  long field_idx;      // idx of the data in the field
  long buffer_idx;     // idx of the data in the buffer
  long size;           // number of data site for this pack of data
};

struct API ShufflePlanSendPackInfo {
  long field_idx;   // idx of the data in the field
  long buffer_idx;  // idx of the data in the buffer
  long size;        // number of data site for this pack of data
};

struct API ShufflePlan {
  bool is_no_shuffle;
  Coordinate new_size_node;
  Geometry geo_send;                // geo of the send field
  std::vector<Geometry> geos_recv;  // geos of the recv fields
  long n_elems_send;  // n_elems for only send field, useful in sparse field
  std::vector<long>
      n_elems_recv;  // n_elems for recv fields, only useful in sparse field
  std::vector<ShufflePlanSendPackInfo>
      send_pack_infos;  // corresponds to how to create send buffer from local
                        // field
  std::vector<ShufflePlanRecvPackInfo>
      recv_pack_infos;  // corresponds to how to copy recv buffer to new local
                        // fields
  ShuffleCommPlan scp;  // plan for comm
  //
  ShufflePlan()
  {
    is_no_shuffle = false;
    n_elems_send = 0;
  }
};

inline bool is_no_shuffle(const ShufflePlan& sp) { return sp.is_no_shuffle; }

inline bool is_initialized(const ShufflePlan& sp)
{
  return is_initialized(sp.geo_send) or is_no_shuffle(sp);
}

template <class M>
void shuffle_field_comm(Vector<M> recv_buffer, const Vector<M> send_buffer,
                        const ShuffleCommPlan& scp, const int multiplicity)
{
  sync_node();
  TIMER_FLOPS("shuffle_field_comm(recv,send,scp,multiplicity)");
  const long total_bytes = scp.global_comm_size * multiplicity * sizeof(M);
  timer.flops += total_bytes;
  qassert(send_buffer.size() == scp.total_send_size * multiplicity);
  qassert(recv_buffer.size() == scp.total_recv_size * multiplicity);
  std::vector<MPI_Request> reqs;
  {
    TIMER("shuffle_field_comm-init");
    const int mpi_tag = 4;
    for (size_t i = 0; i < scp.recv_msg_infos.size(); ++i) {
      const ShufflePlanMsgInfo& mi = scp.recv_msg_infos[i];
      mpi_irecv(&recv_buffer[mi.idx * multiplicity],
                mi.size * multiplicity * sizeof(M), MPI_BYTE, mi.id_node,
                mpi_tag, get_comm(), reqs);
    }
    for (size_t i = 0; i < scp.send_msg_infos.size(); ++i) {
      const ShufflePlanMsgInfo& mi = scp.send_msg_infos[i];
      mpi_isend(&send_buffer[mi.idx * multiplicity],
                mi.size * multiplicity * sizeof(M), MPI_BYTE, mi.id_node,
                mpi_tag, get_comm(), reqs);
    }
  }
  mpi_waitall(reqs);
  sync_node();
}

template <class M>
void shuffle_field_comm_back(Vector<M> send_buffer, const Vector<M> recv_buffer,
                             const ShuffleCommPlan& scp, const int multiplicity)
// name is reversed
{
  sync_node();
  TIMER_FLOPS("shuffle_field_comm_back(send,recv,scp,multiplicity)");
  const long total_bytes = scp.global_comm_size * multiplicity * sizeof(M);
  timer.flops += total_bytes;
  std::vector<MPI_Request> reqs;
  {
    TIMER("shuffle_field_comm_back-init");
    const int mpi_tag = 5;
    for (size_t i = 0; i < scp.send_msg_infos.size(); ++i) {
      const ShufflePlanMsgInfo& mi = scp.send_msg_infos[i];
      mpi_irecv(&send_buffer[mi.idx * multiplicity],
                mi.size * multiplicity * sizeof(M), MPI_BYTE, mi.id_node,
                mpi_tag, get_comm(), reqs);
    }
    for (size_t i = 0; i < scp.recv_msg_infos.size(); ++i) {
      const ShufflePlanMsgInfo& mi = scp.recv_msg_infos[i];
      mpi_isend(&recv_buffer[mi.idx * multiplicity],
                mi.size * multiplicity * sizeof(M), MPI_BYTE, mi.id_node,
                mpi_tag, get_comm(), reqs);
    }
  }
  mpi_waitall(reqs);
  sync_node();
}

template <class M>
void shuffle_field_pack_send(
    Vector<M> send_buffer, const Vector<M> fdata,
    const std::vector<ShufflePlanSendPackInfo>& send_pack_infos,
    const int multiplicity)
{
  TIMER_FLOPS("shuffle_field_pack_send(send_buffer,fdata)");
#pragma omp parallel for
  for (size_t i = 0; i < send_pack_infos.size(); ++i) {
    const ShufflePlanSendPackInfo& pi = send_pack_infos[i];
    memcpy(&send_buffer[pi.buffer_idx * multiplicity],
           &fdata[pi.field_idx * multiplicity],
           pi.size * multiplicity * sizeof(M));
  }
}

template <class M>
void shuffle_field_unpack_send(
    Vector<M> fdata, const Vector<M> send_buffer,
    const std::vector<ShufflePlanSendPackInfo>& send_pack_infos,
    const int multiplicity)
{
  TIMER_FLOPS("shuffle_field_unpack_send(fdata,send_buffer)");
#pragma omp parallel for
  for (size_t i = 0; i < send_pack_infos.size(); ++i) {
    const ShufflePlanSendPackInfo& pi = send_pack_infos[i];
    memcpy(&fdata[pi.field_idx * multiplicity],
           &send_buffer[pi.buffer_idx * multiplicity],
           pi.size * multiplicity * sizeof(M));
  }
}

template <class M>
void shuffle_field_unpack_recv(
    vector<Vector<M> >& fsdata, const Vector<M> recv_buffer,
    const std::vector<ShufflePlanRecvPackInfo>& recv_pack_infos,
    const int multiplicity)
{
  TIMER_FLOPS("shuffle_field_unpack_recv(fsdata,recv_buffer)");
#pragma omp parallel for
  for (size_t i = 0; i < recv_pack_infos.size(); ++i) {
    const ShufflePlanRecvPackInfo& pi = recv_pack_infos[i];
    qassert(0 <= pi.local_geos_idx && pi.local_geos_idx < (int)fsdata.size());
    std::memcpy(&fsdata[pi.local_geos_idx][pi.field_idx * multiplicity],
                &recv_buffer[pi.buffer_idx * multiplicity],
                pi.size * multiplicity * sizeof(M));
  }
}

template <class M>
void shuffle_field_pack_recv(
    Vector<M> recv_buffer, const vector<Vector<M> >& fsdata,
    const std::vector<ShufflePlanRecvPackInfo>& recv_pack_infos,
    const int multiplicity)
{
  TIMER_FLOPS("shuffle_field_pack_recv(recv_buffer,fsdata)");
#pragma omp parallel for
  for (size_t i = 0; i < recv_pack_infos.size(); ++i) {
    const ShufflePlanRecvPackInfo& pi = recv_pack_infos[i];
    qassert(0 <= pi.local_geos_idx && pi.local_geos_idx < (int)fsdata.size());
    qassert(fsdata[pi.local_geos_idx].size() >=
            (pi.field_idx + pi.size) * multiplicity);
    qassert(recv_buffer.size() >= (pi.buffer_idx + pi.size) * multiplicity);
    std::memcpy(&recv_buffer[pi.buffer_idx * multiplicity],
                &fsdata[pi.local_geos_idx][pi.field_idx * multiplicity],
                pi.size * multiplicity * sizeof(M));
  }
}

template <class M>
void shuffle_field(std::vector<Field<M> >& fs, const Field<M>& f,
                   const ShufflePlan& sp)
{
  qassert(is_initialized(sp));
  if (is_no_shuffle(sp)) {
    clear(fs);
    fs.resize(1);
    fs[0] = f;
    return;
  }
  sync_node();
  TIMER_VERBOSE_FLOPS("shuffle_field(fs,f,sp)");
  const Geometry& geo = f.geo();
  if (sp.new_size_node != Coordinate()) {
    displayln_info(
        0,
        fname + ssprintf(": %s -> %s (total_site: %s ; site_size: %d ; "
                         "total_size: %.3lf GB)",
                         show(geo.geon.size_node).c_str(),
                         show(sp.new_size_node).c_str(),
                         show(geo.total_site()).c_str(),
                         geo.multiplicity * (int)sizeof(M),
                         (double)(sp.scp.global_comm_size * geo.multiplicity *
                                  sizeof(M) * std::pow(0.5, 30))));
  }
  qassert(sp.geo_send == geo_reform(geo, 1, 0));
  clear(fs);
  const long total_bytes =
      sp.scp.global_comm_size * geo.multiplicity * sizeof(M);
  timer.flops += total_bytes;
  vector<M> send_buffer(sp.scp.total_send_size * geo.multiplicity);
  shuffle_field_pack_send(get_data(send_buffer), get_data(f),
                          sp.send_pack_infos, geo.multiplicity);
  vector<M> recv_buffer(sp.scp.total_recv_size * geo.multiplicity);
  shuffle_field_comm(get_data(recv_buffer), get_data(send_buffer), sp.scp,
                     geo.multiplicity);
  clear(send_buffer);
  vector<Geometry> geos_recv = sp.geos_recv;
  fs.resize(geos_recv.size());
  for (size_t i = 0; i < fs.size(); ++i) {
    geos_recv[i].remult(geo.multiplicity);
    fs[i].init(geos_recv[i]);
  }
  vector<Vector<M> > fsdata(fs.size());
  for (size_t i = 0; i < fs.size(); ++i) {
    fsdata[i] = get_data(fs[i]);
  }
  shuffle_field_unpack_recv(fsdata, get_data(recv_buffer), sp.recv_pack_infos,
                            geo.multiplicity);
  sync_node();
}

template <class M>
void shuffle_field_back(Field<M>& f, const std::vector<Field<M> >& fs,
                        const ShufflePlan& sp)
// f needs to have correct size
{
  qassert(is_initialized(sp));
  if (is_no_shuffle(sp)) {
    qassert(fs.size() == 1);
    qassert(is_initialized(fs[0]));
    f.init();
    f = fs[0];
    return;
  }
  sync_node();
  TIMER_VERBOSE_FLOPS("shuffle_field_back(f,fs,sp)");
  qassert(is_initialized(f));
  const Geometry& geo = f.geo();
  if (sp.new_size_node != Coordinate()) {
    displayln_info(
        0,
        fname + ssprintf(": %s -> %s (total_site: %s ; site_size: %d ; "
                         "total_size: %.3lf GB)",
                         show(sp.new_size_node).c_str(),
                         show(geo.geon.size_node).c_str(),
                         show(geo.total_site()).c_str(),
                         geo.multiplicity * (int)sizeof(M),
                         (double)(sp.scp.global_comm_size * geo.multiplicity *
                                  sizeof(M) * std::pow(0.5, 30))));
  }
  const long total_bytes =
      sp.scp.global_comm_size * geo.multiplicity * sizeof(M);
  timer.flops += total_bytes;
  vector<Vector<M> > fsdata(fs.size());
  for (size_t i = 0; i < fs.size(); ++i) {
    fsdata[i] = get_data(fs[i]);
  }
  vector<M> recv_buffer(sp.scp.total_recv_size * geo.multiplicity);
  shuffle_field_pack_recv(get_data(recv_buffer), fsdata, sp.recv_pack_infos,
                          geo.multiplicity);
  vector<M> send_buffer(sp.scp.total_send_size * geo.multiplicity);
  shuffle_field_comm_back(get_data(send_buffer), get_data(recv_buffer), sp.scp,
                          geo.multiplicity);
  clear(recv_buffer);
  shuffle_field_unpack_send(get_data(f), get_data(send_buffer),
                            sp.send_pack_infos, geo.multiplicity);
  sync_node();
}

template <class M>
void shuffle_field(std::vector<SelectedField<M> >& fs,
                   const SelectedField<M>& f, const ShufflePlan& sp)
{
  qassert(is_initialized(sp));
  if (is_no_shuffle(sp)) {
    clear(fs);
    fs.resize(1);
    fs[0] = f;
    return;
  }
  sync_node();
  TIMER_VERBOSE_FLOPS("shuffle_field(sel_fs,sel_f,sp)");
  const Geometry& geo = f.geo();
  displayln_info(
      0, fname + ssprintf(": %s -> %s (total_site: %s ; site_size: %d ; "
                          "total_size: %.3lf GB)",
                          show(geo.geon.size_node).c_str(),
                          show(sp.new_size_node).c_str(),
                          show(geo.total_site()).c_str(),
                          geo.multiplicity * (int)sizeof(M),
                          (double)(sp.scp.global_comm_size * geo.multiplicity *
                                   sizeof(M) * std::pow(0.5, 30))));
  qassert(sp.geo_send == geo_reform(geo, 1, 0));
  clear(fs);
  const long total_bytes =
      sp.scp.global_comm_size * geo.multiplicity * sizeof(M);
  timer.flops += total_bytes;
  qassert(sp.n_elems_send * geo.multiplicity == (long)f.field.size());
  vector<M> send_buffer(sp.scp.total_send_size * geo.multiplicity);
  shuffle_field_pack_send(get_data(send_buffer), get_data(f),
                          sp.send_pack_infos, geo.multiplicity);
  vector<M> recv_buffer(sp.scp.total_recv_size * geo.multiplicity);
  shuffle_field_comm(get_data(recv_buffer), get_data(send_buffer), sp.scp,
                     geo.multiplicity);
  clear(send_buffer);
  vector<Geometry> geos_recv = sp.geos_recv;
  fs.resize(geos_recv.size());
  for (size_t i = 0; i < fs.size(); ++i) {
    geos_recv[i].remult(geo.multiplicity);
    fs[i].init(geos_recv[i], sp.n_elems_recv[i], geo.multiplicity);
  }
  vector<Vector<M> > fsdata(fs.size());
  for (size_t i = 0; i < fs.size(); ++i) {
    fsdata[i] = get_data(fs[i]);
  }
  shuffle_field_unpack_recv(fsdata, get_data(recv_buffer), sp.recv_pack_infos,
                            geo.multiplicity);
  sync_node();
}

template <class M>
void shuffle_field_back(SelectedField<M>& f,
                        const std::vector<SelectedField<M> >& fs,
                        const ShufflePlan& sp)
// f needs to have correct size
{
  qassert(is_initialized(sp));
  if (is_no_shuffle(sp)) {
    qassert(fs.size() == 1);
    qassert(is_initialized(fs[0]));
    f.init();
    f = fs[0];
    return;
  }
  sync_node();
  TIMER_VERBOSE_FLOPS("shuffle_field_back(sel_f,sel_fs,sp)");
  qassert(is_initialized(f));
  const Geometry& geo = f.geo();
  displayln_info(
      0, fname + ssprintf(": %s -> %s (total_site: %s ; site_size: %d ; "
                          "total_size: %.3lf GB)",
                          show(sp.new_size_node).c_str(),
                          show(geo.geon.size_node).c_str(),
                          show(geo.total_site()).c_str(),
                          geo.multiplicity * (int)sizeof(M),
                          (double)(sp.scp.global_comm_size * geo.multiplicity *
                                   sizeof(M) * std::pow(0.5, 30))));
  const long total_bytes =
      sp.scp.global_comm_size * geo.multiplicity * sizeof(M);
  timer.flops += total_bytes;
  vector<Vector<M> > fsdata(fs.size());
  for (size_t i = 0; i < fs.size(); ++i) {
    fsdata[i] = get_data(fs[i]);
  }
  vector<M> recv_buffer(sp.scp.total_recv_size * geo.multiplicity);
  shuffle_field_pack_recv(get_data(recv_buffer), fsdata, sp.recv_pack_infos,
                          geo.multiplicity);
  vector<M> send_buffer(sp.scp.total_send_size * geo.multiplicity);
  shuffle_field_comm_back(get_data(send_buffer), get_data(recv_buffer), sp.scp,
                          geo.multiplicity);
  clear(recv_buffer);
  shuffle_field_unpack_send(get_data(f), get_data(send_buffer),
                            sp.send_pack_infos, geo.multiplicity);
  sync_node();
}

template <class M>
void shuffle_field(SelectedField<M>& sf, const SelectedField<M>& sf0,
                   const ShufflePlan& sp)
{
  qassert(is_initialized(sp));
  if (is_no_shuffle(sp)) {
    sf = sf0;
    return;
  }
  TIMER_VERBOSE_FLOPS("shuffle_field(sf,sf0,sp)");
  qassert(sp.geos_recv.size() == 1);
  const Geometry& geo = sf0.geo();
  timer.flops += sp.scp.global_comm_size * geo.multiplicity * sizeof(M);
  std::vector<SelectedField<M> > sfs;
  shuffle_field(sfs, sf0, sp);
  qassert(sfs.size() == 1);
  qswap(sf, sfs[0]);
}

template <class M>
void shuffle_field_back(SelectedField<M>& sf, const SelectedField<M>& sf0,
                        const ShufflePlan& sp)
// sf can be empty
{
  qassert(is_initialized(sp));
  if (is_no_shuffle(sp)) {
    sf = sf0;
    return;
  }
  TIMER_VERBOSE_FLOPS("shuffle_field_back(sf,sf0,sp)");
  qassert(sp.geos_recv.size() == 1);
  const Geometry& geo = sf0.geo();
  timer.flops += sp.scp.global_comm_size * geo.multiplicity * sizeof(M);
  std::vector<SelectedField<M> > sfs(1);
  sfs[0] = sf0;
  sf.init(geo, sp.n_elems_send, geo.multiplicity);
  shuffle_field_back(sf, sfs, sp);
}

// making shuffle plan generic

API inline long& get_shuffle_max_msg_size()
// qlat parameter
{
  static long size = 16 * 1024;
  return size;
}

API inline long& get_shuffle_max_pack_size()
// qlat parameter
{
  static long size = 32;
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

template <class Func>
ShufflePlan make_shuffle_plan_generic(std::vector<FieldSelection>& fsels,
                                      const FieldSelection& fsel,
                                      const Coordinate& new_size_node,
                                      const Func& func)
// gindex_s = func(gindex)
// return sp
// shuffle_field(fs, f0, sp)
// fs{gindex_s} == f0{gindex}
{
  TIMER_VERBOSE("make_shuffle_plan_generic");
  fsels.clear();
  ShufflePlan sp;
  sp.new_size_node = new_size_node;
  sp.is_no_shuffle = false;
  sp.geo_send = fsel.f_rank.geo();
  sp.geos_recv = make_dist_io_geos(sp.geo_send.total_site(),
                                   sp.geo_send.multiplicity, sp.new_size_node);
  fsels.resize(sp.geos_recv.size());
  sp.n_elems_send = fsel.n_elems;
  const Coordinate total_site = sp.geo_send.total_site();
  const Coordinate new_node_site = total_site / sp.new_size_node;
  qassert(sp.new_size_node * new_node_site == total_site);
  const int new_num_node = product(sp.new_size_node);
  const int num_node = sp.geo_send.geon.num_node;
  // global index for each selected site
  SelectedField<long> sf_rank;
  sf_rank.init(fsel, 1);
  // shuffled global index for each selected site
  SelectedField<long> sf_gindex_s;
  sf_gindex_s.init(fsel, 1);
  // target node id_node (real id_node) for each selected site
  SelectedField<int> sf_id_node_send;
  sf_id_node_send.init(fsel, 1);
#pragma omp parallel for
  for (long idx = 0; idx < fsel.n_elems; ++idx) {
    const long index = fsel.indices[idx];
    const Coordinate xl = sp.geo_send.coordinate_from_index(index);
    const Coordinate xg = sp.geo_send.coordinate_g_from_l(xl);
    const long gindex = index_from_coordinate(xg, total_site);
    const long gindex_s = func(gindex);
    const Coordinate xg_s = coordinate_from_index(gindex_s, total_site);
    const Coordinate new_coor_node = xg_s / new_node_site;
    const int new_id_node =
        index_from_coordinate(new_coor_node, sp.new_size_node);
    const int id_node =
        get_id_node_from_new_id_node(new_id_node, new_num_node, num_node);
    sf_rank.get_elem(idx) = fsel.f_rank.get_elem(xl);
    sf_gindex_s.get_elem(idx) = gindex_s;
    sf_id_node_send.get_elem(idx) = id_node;
  }
  std::map<int, long> send_id_node_size, send_id_node_buffer_idx;
  for (long idx = 0; idx < fsel.n_elems; ++idx) {
    send_id_node_size[sf_id_node_send.get_elem(idx)] += 1;
  }
  // send_msg_infos
  {
    long count = 0;
    for (auto it = send_id_node_size.cbegin(); it != send_id_node_size.cend();
         ++it) {
      const int id_node = it->first;
      const long node_size = it->second;
      send_id_node_buffer_idx[id_node] = count;
      long node_size_remain = node_size;
      while (node_size_remain > 0) {
        ShufflePlanMsgInfo mi;
        mi.id_node = id_node;
        mi.idx = count;
        mi.size = std::min(node_size_remain, get_shuffle_max_msg_size());
        sp.scp.send_msg_infos.push_back(mi);
        node_size_remain -= mi.size;
        count += mi.size;
      }
    }
    sp.scp.total_send_size = count;
    qassert(count == fsel.n_elems);
  }
  qassert(sp.scp.total_send_size == sp.n_elems_send);
  sp.scp.global_comm_size = sp.scp.total_send_size;
  glb_sum(sp.scp.global_comm_size);
  // send_pack_infos
  {
    long last_buffer_idx = -1;
    for (long idx = 0; idx < fsel.n_elems; ++idx) {
      const int id_node = sf_id_node_send.get_elem(idx);
      const long buffer_idx = send_id_node_buffer_idx[id_node];
      if (buffer_idx == last_buffer_idx and
          sp.send_pack_infos.back().size < get_shuffle_max_pack_size()) {
        sp.send_pack_infos.back().size += 1;
      } else {
        last_buffer_idx = buffer_idx;
        ShufflePlanSendPackInfo pi;
        pi.field_idx = idx;
        pi.buffer_idx = buffer_idx;
        pi.size = 1;
        sp.send_pack_infos.push_back(pi);
      }
      send_id_node_buffer_idx[id_node] += 1;
      last_buffer_idx += 1;
    }
  }
  // communicate to determine recv msg size from each node
  std::map<int, long> recv_id_node_size;
  {
    vector<long> send_size(num_node, 0), recv_size(num_node, -1);
    for (auto it = send_id_node_size.cbegin(); it != send_id_node_size.cend();
         ++it) {
      const int id_node = it->first;
      const long node_size = it->second;
      qassert(0 <= id_node and id_node < num_node);
      qassert(0 <= node_size);
      send_size[id_node] = node_size;
    }
    long neg_count = 0;
    do {
      std::vector<MPI_Request> reqs;
      const int mpi_tag = 12;
      for (int i = 0; i < num_node; ++i) {
        mpi_irecv(&recv_size[i], 1, MPI_LONG, i, mpi_tag, get_comm(), reqs);
      }
      for (int i = 0; i < num_node; ++i) {
        mpi_isend(&send_size[i], 1, MPI_LONG, i, mpi_tag, get_comm(), reqs);
      }
      mpi_waitall(reqs);
      neg_count = 0;
      for (int i = 0; i < num_node; ++i) {
        if (recv_size[i] < 0) {
          qwarn(fname +
                ssprintf(": id_node=%d i=%d recv_size[i]=%ld neg_count=%ld",
                         get_id_node(), i, recv_size[i], neg_count));
          neg_count += 1;
        } else if (recv_size[i] > 0) {
          recv_id_node_size[i] = recv_size[i];
        }
      }
      glb_sum(neg_count);
      if (neg_count > 0) {
        if (get_id_node() == 0) {
          qwarn(fname +
                ssprintf(
                    ": WARNING: total neg_count=%ld but should be 0. Retrying.",
                    neg_count));
        }
      }
    } while (neg_count > 0);
  }
  // recv_msg_infos
  {
    long count = 0;
    for (auto it = recv_id_node_size.cbegin(); it != recv_id_node_size.cend();
         ++it) {
      const int id_node = it->first;
      const long node_size = it->second;
      long node_size_remain = node_size;
      while (node_size_remain > 0) {
        ShufflePlanMsgInfo mi;
        mi.id_node = id_node;
        mi.idx = count;
        mi.size = std::min(node_size_remain, get_shuffle_max_msg_size());
        sp.scp.recv_msg_infos.push_back(mi);
        node_size_remain -= mi.size;
        count += mi.size;
      }
    }
    sp.scp.total_recv_size = count;
  }
  // exec shuffle for sf_gindex_s
  vector<long> send_buffer(sp.scp.total_send_size);
  shuffle_field_pack_send(get_data(send_buffer), get_data(sf_gindex_s),
                          sp.send_pack_infos, 1);
  vector<long> recv_buffer(sp.scp.total_recv_size);
  shuffle_field_comm(get_data(recv_buffer), get_data(send_buffer), sp.scp, 1);
  shuffle_field_pack_send(get_data(send_buffer), get_data(sf_rank),
                          sp.send_pack_infos, 1);
  vector<long> recv_buffer_rank(sp.scp.total_recv_size);
  shuffle_field_comm(get_data(recv_buffer_rank), get_data(send_buffer), sp.scp,
                     1);
  clear(send_buffer);
  // recv_pack_infos
  {
    vector<int> local_geos_idx_from_new_id_node(new_num_node, -1);
    for (int i = 0; i < (int)fsels.size(); ++i) {
      const int local_geos_idx = i;
      const Geometry& geo_recv = sp.geos_recv[i];
      local_geos_idx_from_new_id_node[geo_recv.geon.id_node] = local_geos_idx;
    }
    std::vector<FieldM<int64_t, 1> > f_ranks(fsels.size());
    for (int i = 0; i < (int)f_ranks.size(); ++i) {
      f_ranks[i].init(sp.geos_recv[i], 1);
#pragma omp parallel for
      for (long index = 0; index < f_ranks[i].geo().local_volume(); ++index) {
        int64_t& rank = f_ranks[i].get_elem(index);
        rank = -1;
      }
    }
#pragma omp parallel for
    for (long buffer_idx = 0; buffer_idx < (long)recv_buffer.size();
         ++buffer_idx) {
      const long gindex_s = recv_buffer[buffer_idx];
      const Coordinate xg_s = coordinate_from_index(gindex_s, total_site);
      const Coordinate new_coor_node = xg_s / new_node_site;
      const int new_id_node =
          index_from_coordinate(new_coor_node, sp.new_size_node);
      const int local_geos_idx = local_geos_idx_from_new_id_node[new_id_node];
      const Geometry& geo_recv = sp.geos_recv[local_geos_idx];
      const Coordinate xl_s = geo_recv.coordinate_l_from_g(xg_s);
      f_ranks[local_geos_idx].get_elem(xl_s) = recv_buffer_rank[buffer_idx];
    }
    for (int i = 0; i < (int)fsels.size(); ++i) {
      set_field_selection(fsels[i], f_ranks[i], fsel.n_per_tslice);
    }
    long last_local_geos_idx = -1;
    vector<long> last_field_idx(sp.geos_recv.size(), 0);
    for (long buffer_idx = 0; buffer_idx < (long)recv_buffer.size();
         ++buffer_idx) {
      const long gindex_s = recv_buffer[buffer_idx];
      const Coordinate xg_s = coordinate_from_index(gindex_s, total_site);
      const Coordinate new_coor_node = xg_s / new_node_site;
      const int new_id_node =
          index_from_coordinate(new_coor_node, sp.new_size_node);
      const int local_geos_idx = local_geos_idx_from_new_id_node[new_id_node];
      qassert(local_geos_idx >= 0);
      qassert(
          get_id_node_from_new_id_node(new_id_node, new_num_node, num_node) ==
          sp.geo_send.geon.id_node);
      const Geometry& geo_recv = sp.geos_recv[local_geos_idx];
      qassert(geo_recv.geon.id_node == new_id_node);
      bool is_continue = local_geos_idx == last_local_geos_idx;
      long field_idx;
      const Coordinate xl_s = geo_recv.coordinate_l_from_g(xg_s);
      field_idx = fsels[local_geos_idx].f_local_idx.get_elem(xl_s);
      is_continue = is_continue and field_idx == last_field_idx[local_geos_idx];
      last_field_idx[local_geos_idx] = field_idx;
      last_field_idx[local_geos_idx] += 1;
      if (is_continue and
          sp.recv_pack_infos.back().size < get_shuffle_max_pack_size()) {
        sp.recv_pack_infos.back().size += 1;
      } else {
        last_local_geos_idx = local_geos_idx;
        ShufflePlanRecvPackInfo pi;
        pi.local_geos_idx = local_geos_idx;
        pi.field_idx = field_idx;
        pi.buffer_idx = buffer_idx;
        pi.size = 1;
        sp.recv_pack_infos.push_back(pi);
      }
    }
    sp.n_elems_recv.resize(fsels.size());
    for (int i = 0; i < (int)fsels.size(); ++i) {
      sp.n_elems_recv[i] = fsels[i].n_elems;
    }
  }
  long num_send_packs = sp.send_pack_infos.size();
  long num_recv_packs = sp.recv_pack_infos.size();
  long num_send_msgs = sp.scp.send_msg_infos.size();
  long num_recv_msgs = sp.scp.recv_msg_infos.size();
  displayln_info(0,
                 fname + ssprintf(": num_send_packs = %10ld", num_send_packs));
  displayln_info(0,
                 fname + ssprintf(": num_recv_packs = %10ld", num_recv_packs));
  displayln_info(0,
                 fname + ssprintf(": num_send_msgs  = %10ld", num_send_msgs));
  displayln_info(0,
                 fname + ssprintf(": num_recv_msgs  = %10ld", num_recv_msgs));
  glb_sum(num_send_packs);
  glb_sum(num_recv_packs);
  glb_sum(num_send_msgs);
  glb_sum(num_recv_msgs);
  displayln_info(
      0, fname + ssprintf(": total num_send_packs = %10ld", num_send_packs));
  displayln_info(
      0, fname + ssprintf(": total num_recv_packs = %10ld", num_recv_packs));
  displayln_info(
      0, fname + ssprintf(": total num_send_msgs  = %10ld", num_send_msgs));
  displayln_info(
      0, fname + ssprintf(": total num_recv_msgs  = %10ld", num_recv_msgs));
  displayln_info(0, fname + ssprintf(": global_comm_size = %10ld",
                                     sp.scp.global_comm_size));
  return sp;
}

inline ShufflePlan make_shuffle_plan(std::vector<FieldSelection>& fsels,
                                     const FieldSelection& fsel,
                                     const Coordinate& new_size_node)
{
  TIMER("make_shuffle_plan");
  if (new_size_node == fsel.f_rank.geo().geon.size_node) {
    fsels.clear();
    fsels.resize(1);
    fsels[0] = fsel;
    ShufflePlan sp;
    sp.is_no_shuffle = true;
    sp.new_size_node = new_size_node;
    sp.scp.global_comm_size = fsel.f_rank.geo().total_volume();
    return sp;
  }
  return make_shuffle_plan_generic(fsels, fsel, new_size_node, identity<long>);
}

// field dist write shuffle

struct API ShufflePlanKey {
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
{
  FieldSelection fsel;
  set_field_selection(fsel, spk.total_site);
  std::vector<FieldSelection> fsels;
  return make_shuffle_plan(fsels, fsel, spk.new_size_node);
}

API inline Cache<ShufflePlanKey, ShufflePlan>& get_shuffle_plan_cache()
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
  sync_node();
  TIMER_FLOPS("shuffle_field");
  const Geometry& geo = f.geo();
  clear(fs);
  const ShufflePlan& sp = get_shuffle_plan(geo.total_site(), new_size_node);
  shuffle_field(fs, f, sp);
  const long total_bytes =
      sp.scp.global_comm_size * geo.multiplicity * sizeof(M);
  timer.flops += total_bytes;
}

template <class M>
void shuffle_field_back(Field<M>& f, const std::vector<Field<M> >& fs,
                        const Coordinate& new_size_node)
// f needs to have correct size
{
  sync_node();
  TIMER_FLOPS("shuffle_field_back");
  const ShufflePlan& sp = get_shuffle_plan(f.geo().total_site(), new_size_node);
  shuffle_field_back(f, fs, sp);
  const long total_bytes =
      sp.scp.global_comm_size * f.geo().multiplicity * sizeof(M);
  timer.flops += total_bytes;
}

// sparse field dist write shuffle

template <class M>
void shuffle_field(std::vector<SelectedField<M> >& fs,
                   std::vector<FieldSelection>& fsels,
                   const SelectedField<M>& f, const Coordinate& new_size_node,
                   const FieldSelection& fsel)
{
  sync_node();
  TIMER_FLOPS("shuffle_field(fs,f,nsn,fsel)");
  const Geometry& geo = f.geo();
  const ShufflePlan sp = make_shuffle_plan(fsels, fsel, new_size_node);
  shuffle_field(fs, f, sp);
  const long total_bytes =
      sp.scp.global_comm_size * geo.multiplicity * sizeof(M);
  timer.flops += total_bytes;
}

template <class M>
void shuffle_field_back(SelectedField<M>& f,
                        const std::vector<SelectedField<M> >& fs,
                        const Coordinate& new_size_node,
                        const FieldSelection& fsel)
// f needs to have correct size
{
  sync_node();
  TIMER_FLOPS("shuffle_field_back");
  std::vector<FieldSelection> fsels;
  const ShufflePlan sp = make_shuffle_plan(fsels, fsel, new_size_node);
  shuffle_field_back(f, fs, sp);
  const long total_bytes =
      sp.scp.global_comm_size * f.geo().multiplicity * sizeof(M);
  timer.flops += total_bytes;
}

// fft shuffle (can not yet use generic code)

inline long index_from_coordinate_perp_dir(const Coordinate& coor,
                                           const Coordinate& size,
                                           const int dir)
{
  array<int, 3> dirs;
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
  array<int, 3> dirs;
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
  ret.scp.total_send_size = geo.local_volume();
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
        ret.scp.send_msg_infos.push_back(mi);
        node_size_remain -= mi.size;
        count += mi.size;
      }
    }
    qassert(count == geo.local_volume());
    qassert(count == ret.scp.total_send_size);
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
    qassert(count == ret.scp.total_send_size);
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
  ret.scp.total_recv_size = 0;
  for (size_t i = 0; i < ret.geos_recv.size(); ++i) {
    const Geometry& geo_recv = ret.geos_recv[i];
    ret.scp.total_recv_size += geo_recv.local_volume();
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
        ret.scp.recv_msg_infos.push_back(mi);
        node_size_remain -= mi.size;
        count += mi.size;
      }
    }
    if (count != ret.scp.total_recv_size) {
      qwarn(fname + ssprintf(": count = %ld", count));
      qwarn(fname + ssprintf(": ret.scp.total_recv_size = %ld",
                             ret.scp.total_recv_size));
    }
    qassert(count == ret.scp.total_recv_size);
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
    qassert(count == ret.scp.total_recv_size);
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
  long num_send_msgs = ret.scp.send_msg_infos.size();
  long num_recv_msgs = ret.scp.recv_msg_infos.size();
  displayln_info(0,
                 fname + ssprintf(": num_send_packs = %10ld", num_send_packs));
  displayln_info(0,
                 fname + ssprintf(": num_recv_packs = %10ld", num_recv_packs));
  displayln_info(0,
                 fname + ssprintf(": num_send_msgs  = %10ld", num_send_msgs));
  displayln_info(0,
                 fname + ssprintf(": num_recv_msgs  = %10ld", num_recv_msgs));
  glb_sum(num_send_packs);
  glb_sum(num_recv_packs);
  glb_sum(num_send_msgs);
  glb_sum(num_recv_msgs);
  displayln_info(
      0, fname + ssprintf(": total num_send_packs = %10ld", num_send_packs));
  displayln_info(
      0, fname + ssprintf(": total num_recv_packs = %10ld", num_recv_packs));
  displayln_info(
      0, fname + ssprintf(": total num_send_msgs  = %10ld", num_send_msgs));
  displayln_info(
      0, fname + ssprintf(": total num_recv_msgs  = %10ld", num_recv_msgs));
  ret.scp.global_comm_size = ret.scp.total_send_size;
  glb_sum(ret.scp.global_comm_size);
  displayln_info(0, fname + ssprintf(": global_comm_size = %10ld",
                                     ret.scp.global_comm_size));
  return ret;
}

// reflection shuffle

struct API ShuffleReflectGIndexMap {
  Coordinate total_site;
  long operator()(const long gindex) const
  {
    const Coordinate xg = coordinate_from_index(gindex, total_site);
    const Coordinate xg_reflect = mod(-xg, total_site);
    return index_from_coordinate(xg_reflect, total_site);
  }
};

template <class M>
void reflect_field(Field<M>& f)
{
  TIMER_VERBOSE_FLOPS("reflect_field");
  timer.flops += get_data_size(f) * f.geo().geon.num_node;
  const Geometry& geo = f.geo();
  const Coordinate total_site = geo.total_site();
  ShuffleReflectGIndexMap func;
  func.total_site = total_site;
  FieldSelection fsel;
  set_field_selection(fsel, total_site);
  std::vector<FieldSelection> fsels;
  ShufflePlan sp =
      make_shuffle_plan_generic(fsels, fsel, geo.geon.size_node, func);
  qassert(fsels.size() == 1);
  std::vector<Field<M> > fs;
  shuffle_field(fs, f, sp);
  qassert(fs.size() == 1);
  qswap(f, fs[0]);
}

// shift shuffle

struct API ShuffleShiftGIndexMap {
  Coordinate total_site;
  Coordinate shift;
  bool is_reflect;  // shift and then reflect
  long operator()(const long gindex) const
  {
    const Coordinate xg = coordinate_from_index(gindex, total_site);
    const Coordinate xg_shift = is_reflect ? mod(-(xg + shift), total_site)
                                           : mod(xg + shift, total_site);
    return index_from_coordinate(xg_shift, total_site);
  }
};

template <class M>
void field_shift_shuffle(Field<M>& f, const Field<M>& f0,
                         const Coordinate& shift, const bool is_reflect = false)
{
  TIMER_VERBOSE_FLOPS("field_shift_shuffle");
  timer.flops += get_data_size(f0) * f0.geo().geon.num_node;
  const Geometry& geo = f0.geo();
  const Coordinate total_site = geo.total_site();
  ShuffleShiftGIndexMap func;
  func.total_site = total_site;
  func.shift = shift;
  func.is_reflect = is_reflect;
  FieldSelection fsel;
  set_field_selection(fsel, total_site);
  std::vector<FieldSelection> fsels;
  ShufflePlan sp =
      make_shuffle_plan_generic(fsels, fsel, geo.geon.size_node, func);
  qassert(fsels.size() == 1);
  std::vector<Field<M> > fs;
  shuffle_field(fs, f0, sp);
  qassert(fs.size() == 1);
  qswap(f, fs[0]);
}

struct API ShiftShufflePlan {
  Coordinate shift;
  bool is_reflect;      // shift and then reflect
  FieldSelection fsel;  // fsel after the shift
  ShufflePlan sp;
};

inline ShufflePlan make_shuffle_plan_shift(FieldSelection& fsel_shift,
                                           const FieldSelection& fsel,
                                           const Coordinate& shift,
                                           const bool is_reflect = false)
{
  TIMER_VERBOSE("make_shuffle_plan_shift");
  const Geometry& geo = fsel.f_rank.geo();
  const Coordinate& new_size_node = geo.geon.size_node;
  ShuffleShiftGIndexMap func;
  func.total_site = geo.total_site();
  func.shift = shift;
  func.is_reflect = is_reflect;
  std::vector<FieldSelection> fsels;
  const ShufflePlan sp =
      make_shuffle_plan_generic(fsels, fsel, new_size_node, func);
  qassert(fsels.size() == 1);
  fsel_shift = fsels[0];
  return sp;
}

inline ShiftShufflePlan make_shift_shuffle_plan(const FieldSelection& fsel,
                                                const Coordinate& shift,
                                                const bool is_reflect = false)
{
  TIMER_VERBOSE("make_shift_shuffle_plan");
  ShiftShufflePlan ssp;
  ssp.shift = shift;
  ssp.is_reflect = is_reflect;
  ssp.sp = make_shuffle_plan_shift(ssp.fsel, fsel, shift, is_reflect);
  return ssp;
}

template <class M>
void field_shift(SelectedField<M>& sf, const SelectedField<M>& sf0,
                 const ShiftShufflePlan& ssp)
// sf{gindex_s} = sf0{gindex}
// xg_s <=> gindex_s
// xg <=> gindex
// xg_s = mod(xg + shift, total_site)
{
  shuffle_field(sf, sf0, ssp.sp);
}

template <class M>
void field_shift(SelectedField<M>& sf, FieldSelection& fsel,
                 const SelectedField<M>& sf0, const FieldSelection& fsel0,
                 const Coordinate& shift, const bool is_reflect = false)
// sf{gindex_s} = sf0{gindex}
// xg_s <=> gindex_s
// xg <=> gindex
// xg_s = mod(xg + shift, total_site)
{
  TIMER_VERBOSE_FLOPS("field_shift(sf,fsel,sf0,fsel0,shift)");
  const ShiftShufflePlan ssp =
      make_shift_shuffle_plan(fsel0, shift, is_reflect);
  const ShufflePlan& sp = ssp.sp;
  const Geometry& geo = sf0.geo();
  timer.flops += sp.scp.global_comm_size * geo.multiplicity * sizeof(M);
  fsel = ssp.fsel;
  shuffle_field(sf, sf0, ssp.sp);
}

// old code

inline ShufflePlan make_shuffle_plan_nofsel_nofunc(const ShufflePlanKey& spk)
// assume the order is preserved in transfer.
// obsolete
{
  TIMER_VERBOSE("make_shuffle_plan_nofsel_nofunc");
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
  ShufflePlan sp;
  // geo_send
  sp.geo_send = geo;
  // geos_recv
  sp.geos_recv = geos_recv;
  // total_send_size
  sp.scp.total_send_size = geo.local_volume();
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
        sp.scp.send_msg_infos.push_back(mi);
        node_size_remain -= mi.size;
        count += mi.size;
      }
    }
    qassert(count == geo.local_volume());
    qassert(count == sp.scp.total_send_size);
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
    qassert(count == sp.scp.total_send_size);
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
          sp.send_pack_infos.back().size < get_shuffle_max_pack_size()) {
        sp.send_pack_infos.back().size += 1;
      } else {
        last_buffer_idx = buffer_idx;
        ShufflePlanSendPackInfo pi;
        pi.field_idx = index;
        pi.buffer_idx = buffer_idx;
        pi.size = 1;
        sp.send_pack_infos.push_back(pi);
      }
      send_new_id_node_idx[new_id_node] += 1;
      last_buffer_idx += 1;
    }
  }
  // total_recv_size
  sp.scp.total_recv_size = 0;
  for (size_t i = 0; i < sp.geos_recv.size(); ++i) {
    const Geometry& geos_recv = sp.geos_recv[i];
    sp.scp.total_recv_size += geos_recv.local_volume();
  }
  // recv_id_node_size
  std::map<int, long> recv_id_node_size;
  for (size_t i = 0; i < sp.geos_recv.size(); ++i) {
    const Geometry& geo_recv = sp.geos_recv[i];
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
        sp.scp.recv_msg_infos.push_back(mi);
        node_size_remain -= mi.size;
        count += mi.size;
      }
    }
    qassert(count == sp.scp.total_recv_size);
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
    qassert(count == sp.scp.total_recv_size);
  }
  // recv_pack_infos
  {
    for (size_t i = 0; i < sp.geos_recv.size(); ++i) {
      const Geometry& geo_recv = sp.geos_recv[i];
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
            sp.recv_pack_infos.back().size < get_shuffle_max_pack_size()) {
          sp.recv_pack_infos.back().size += 1;
        } else {
          last_buffer_idx = buffer_idx;
          ShufflePlanRecvPackInfo pi;
          pi.local_geos_idx = i;
          pi.field_idx = index;
          pi.buffer_idx = buffer_idx;
          pi.size = 1;
          sp.recv_pack_infos.push_back(pi);
        }
        recv_id_node_idx[id_node] += 1;
        last_buffer_idx += 1;
      }
    }
  }
  long num_send_packs = sp.send_pack_infos.size();
  long num_recv_packs = sp.recv_pack_infos.size();
  long num_send_msgs = sp.scp.send_msg_infos.size();
  long num_recv_msgs = sp.scp.recv_msg_infos.size();
  displayln_info(0,
                 fname + ssprintf(": num_send_packs = %10ld", num_send_packs));
  displayln_info(0,
                 fname + ssprintf(": num_recv_packs = %10ld", num_recv_packs));
  displayln_info(0,
                 fname + ssprintf(": num_send_msgs  = %10ld", num_send_msgs));
  displayln_info(0,
                 fname + ssprintf(": num_recv_msgs  = %10ld", num_recv_msgs));
  glb_sum(num_send_packs);
  glb_sum(num_recv_packs);
  glb_sum(num_send_msgs);
  glb_sum(num_recv_msgs);
  displayln_info(
      0, fname + ssprintf(": total num_send_packs = %10ld", num_send_packs));
  displayln_info(
      0, fname + ssprintf(": total num_recv_packs = %10ld", num_recv_packs));
  displayln_info(
      0, fname + ssprintf(": total num_send_msgs  = %10ld", num_send_msgs));
  displayln_info(
      0, fname + ssprintf(": total num_recv_msgs  = %10ld", num_recv_msgs));
  sp.scp.global_comm_size = sp.scp.total_send_size;
  glb_sum(sp.scp.global_comm_size);
  displayln_info(0, fname + ssprintf(": global_comm_size = %10ld",
                                     sp.scp.global_comm_size));
  return sp;
}

}  // namespace qlat
