#pragma once

#include <qlat/setup.h>
#include <qlat/utils-coordinate.h>

#include <qlat-utils/cache.h>
#include <qlat-utils/lat-io.h>

#ifndef QLAT_NO_MALLOPT
#include <malloc.h>
#endif

#ifdef USE_MULTI_NODE
#include <mpi.h>
#endif

namespace qlat
{  //

struct API Q_Comm {
  MPI_Comm comm;
  Coordinate size_node;
  RngState sync_node_rs;
  //
  Q_Comm() {}
  Q_Comm(MPI_Comm comm_, const Coordinate& size_node_,
         const RngState& sync_node_rs_)
  {
    comm = comm_;
    size_node = size_node_;
    sync_node_rs = sync_node_rs_;
  }
};

API inline std::vector<Q_Comm>& get_comm_list()
{
  static std::vector<Q_Comm> comm_list;
  return comm_list;
}

API inline MPI_Comm& get_comm_internal()
{
  static MPI_Comm comm;
  return comm;
}

inline MPI_Comm get_comm() { return get_comm_internal(); }

qacc bool is_initialized(const GeometryNode& geon) { return geon.initialized; }

qacc void init(GeometryNode& geon) { geon.init(); }

inline int id_node_from_coor_node(const Coordinate& coor_node)
{
  return index_from_coordinate(coor_node, get_geometry_node().size_node);
}

inline Coordinate coor_node_from_id_node(int id_node)
{
  return coordinate_from_index(id_node, get_geometry_node().size_node);
}

inline const Coordinate& get_size_node()
{
  return get_geometry_node().size_node;
}

inline const Coordinate& get_coor_node()
{
  return get_geometry_node().coor_node;
}

struct API GeometryNodeNeighbor {
  int dest[2][DIMN];
  // dest[dir][mu]
  // dir = 0, 1 for Plus dir or Minus dir
  // 0 <= mu < 4 for different directions
  //
  void init()
  {
    const Coordinate& coor_node = get_geometry_node().coor_node;
    const Coordinate& size_node = get_geometry_node().size_node;
    for (int mu = 0; mu < DIMN; ++mu) {
      Coordinate coor;
      coor = coor_node;
      ++coor[mu];
      regularize_coordinate(coor, size_node);
      dest[0][mu] = id_node_from_coor_node(coor);
      coor = coor_node;
      --coor[mu];
      regularize_coordinate(coor, size_node);
      dest[1][mu] = id_node_from_coor_node(coor);
    }
  }
  //
  GeometryNodeNeighbor() {}
  GeometryNodeNeighbor(bool initialize)
  {
    if (initialize) {
      init();
    }
  }
};

API inline const GeometryNodeNeighbor& get_geometry_node_neighbor()
{
  static GeometryNodeNeighbor geonb(true);
  return geonb;
}

inline int mpi_send(const void* buf, long count, MPI_Datatype datatype,
                    int dest, int tag, MPI_Comm comm)
{
  const long int_max = INT_MAX;
  if (count <= int_max) {
    return MPI_Send(buf, count, datatype, dest, tag, comm);
  } else {
    int type_size = 0;
    MPI_Type_size(datatype, &type_size);
    uint8_t* cbuf = (uint8_t*)buf;
    while (count > int_max) {
      MPI_Send(cbuf, int_max, datatype, dest, tag, comm);
      cbuf += (long)int_max * type_size;
      count -= int_max;
    }
    return MPI_Send(cbuf, count, datatype, dest, tag, comm);
  }
}

inline int mpi_recv(void* buf, long count, MPI_Datatype datatype, int source,
                    int tag, MPI_Comm comm, MPI_Status* status)
{
  const long int_max = INT_MAX;
  if (count <= int_max) {
    return MPI_Recv(buf, count, datatype, source, tag, comm, status);
  } else {
    int type_size = 0;
    MPI_Type_size(datatype, &type_size);
    uint8_t* cbuf = (uint8_t*)buf;
    while (count > int_max) {
      MPI_Recv(cbuf, int_max, datatype, source, tag, comm, status);
      cbuf += (long)int_max * type_size;
      count -= int_max;
    }
    return MPI_Recv(cbuf, count, datatype, source, tag, comm, status);
  }
}

inline int mpi_isend(const void* buf, long count, MPI_Datatype datatype,
                     int dest, int tag, MPI_Comm comm,
                     std::vector<MPI_Request>& requests)
{
  const long int_max = INT_MAX;
  if (count <= int_max) {
    MPI_Request r;
    int ret = MPI_Isend(buf, count, datatype, dest, tag, comm, &r);
    requests.push_back(r);
    qassert(ret == MPI_SUCCESS);
    return MPI_SUCCESS;
  } else {
    int type_size = 0;
    MPI_Type_size(datatype, &type_size);
    uint8_t* cbuf = (uint8_t*)buf;
    while (count > int_max) {
      mpi_isend(cbuf, int_max, datatype, dest, tag, comm, requests);
      cbuf += (long)int_max * type_size;
      count -= int_max;
    }
    return mpi_isend(cbuf, count, datatype, dest, tag, comm, requests);
  }
}

inline int mpi_irecv(void* buf, long count, MPI_Datatype datatype, int source,
                     int tag, MPI_Comm comm, std::vector<MPI_Request>& requests)
{
  const long int_max = INT_MAX;
  if (count <= int_max) {
    MPI_Request r;
    int ret = MPI_Irecv(buf, count, datatype, source, tag, comm, &r);
    requests.push_back(r);
    qassert(ret == MPI_SUCCESS);
    return MPI_SUCCESS;
  } else {
    int type_size = 0;
    MPI_Type_size(datatype, &type_size);
    uint8_t* cbuf = (uint8_t*)buf;
    while (count > int_max) {
      mpi_irecv(cbuf, int_max, datatype, source, tag, comm, requests);
      cbuf += (long)int_max * type_size;
      count -= int_max;
    }
    return mpi_irecv(cbuf, count, datatype, source, tag, comm, requests);
  }
}

inline int mpi_waitall(std::vector<MPI_Request>& requests)
{
  TIMER("mpi_waitall");
  if (requests.size() > 0) {
    int ret = MPI_Waitall(requests.size(), requests.data(), MPI_STATUS_IGNORE);
    qassert(ret == MPI_SUCCESS);
  }
  return MPI_SUCCESS;
}

template <class M>
inline std::vector<char> pad_flag_data(const int64_t flag, const M& data)
{
  std::vector<char> fdata(8 + sizeof(M), (char)0);
  std::memcpy(&fdata[0], &flag, sizeof(int64_t));
  std::memcpy(&fdata[8], &data, sizeof(M));
  return fdata;
}

template <class M>
inline void extract_flag_data(int64_t& flag, M& data,
                              const std::vector<char>& fdata)
{
  flag = fdata[0];
  std::memcpy(&flag, &fdata[0], sizeof(int64_t));
  std::memcpy(&data, &fdata[8], sizeof(M));
}

template <class N>
inline int receive_job(int64_t& flag, N& data, const int root = 0)
{
  const int mpi_tag = 3;
  const int count = sizeof(int64_t) + sizeof(N);
  std::vector<char> fdata(count, (char)0);
  int ret = mpi_recv(fdata.data(), count, MPI_BYTE, root, mpi_tag, get_comm(),
                     MPI_STATUS_IGNORE);
  extract_flag_data(flag, data, fdata);
  return ret;
}

template <class M>
inline int send_result(const int64_t flag, const M& data, const int root = 0)
{
  const int mpi_tag = 2;
  std::vector<char> fdata = pad_flag_data(flag, data);
  return mpi_send(fdata.data(), fdata.size(), MPI_BYTE, root, mpi_tag,
                  get_comm());
}

template <class N>
inline int send_job(const int64_t flag, const N& data, const int dest)
{
  const int mpi_tag = 3;
  std::vector<char> fdata = pad_flag_data(flag, data);
  return mpi_send(fdata.data(), fdata.size(), MPI_BYTE, dest, mpi_tag,
                  get_comm());
}

template <class M>
inline int receive_result(int& source, int64_t& flag, M& result)
{
  const int mpi_tag = 2;
  const int count = sizeof(int64_t) + sizeof(M);
  std::vector<char> fdata(count, (char)0);
  MPI_Status status;
  const int ret = mpi_recv(fdata.data(), fdata.size(), MPI_BYTE, MPI_ANY_SOURCE,
                           mpi_tag, get_comm(), &status);
  source = status.MPI_SOURCE;
  extract_flag_data(flag, result, fdata);
  return ret;
}

template <class M>
int get_data_dir(Vector<M> recv, const Vector<M>& send, const int dir)
// dir = 0, 1 for Plus dir or Minus dir
{
  TIMER_FLOPS("get_data_dir");
  const int mpi_tag = 0;
  qassert(recv.size() == send.size());
  const long size = recv.size() * sizeof(M);
  timer.flops += size;
#ifdef USE_MULTI_NODE
  const int self_ID = get_id_node();
  const int idf = (self_ID + 1 - 2 * dir + get_num_node()) % get_num_node();
  const int idt = (self_ID - 1 + 2 * dir + get_num_node()) % get_num_node();
  ;
  MPI_Request req;
  MPI_Isend((void*)send.data(), size, MPI_BYTE, idt, mpi_tag, get_comm(), &req);
  const int ret = MPI_Recv(recv.data(), size, MPI_BYTE, idf, mpi_tag,
                           get_comm(), MPI_STATUS_IGNORE);
  MPI_Wait(&req, MPI_STATUS_IGNORE);
  return ret;
#else
  memcpy(recv.data(), send.data(), size);
  return 0;
#endif
}

template <class M>
int get_data_dir_mu(Vector<M> recv, const Vector<M>& send, const int dir,
                    const int mu)
// dir = 0, 1 for Plus dir or Minus dir
// 0 <= mu < 4 for different directions
{
  TIMER_FLOPS("get_data_dir_mu");
  const int mpi_tag = 1;
  qassert(recv.size() == send.size());
  const long size = recv.size() * sizeof(M);
  timer.flops += size;
#ifdef USE_MULTI_NODE
  const GeometryNodeNeighbor& geonb = get_geometry_node_neighbor();
  const int idf = geonb.dest[dir][mu];
  const int idt = geonb.dest[1 - dir][mu];
  MPI_Request req;
  MPI_Isend((void*)send.data(), size, MPI_BYTE, idt, mpi_tag, get_comm(), &req);
  const int ret = mpi_recv(recv.data(), size, MPI_BYTE, idf, mpi_tag,
                           get_comm(), MPI_STATUS_IGNORE);
  MPI_Wait(&req, MPI_STATUS_IGNORE);
  return ret;
#else
  memcpy(recv.data(), send.data(), size);
  return 0;
#endif
}

template <class M>
int get_data_plus_mu(Vector<M> recv, const Vector<M>& send, const int mu)
{
  return get_data_dir_mu(recv, send, 0, mu);
}

template <class M>
int get_data_minus_mu(Vector<M> recv, const Vector<M>& send, const int mu)
{
  return get_data_dir_mu(recv, send, 1, mu);
}

inline int glb_sum(Vector<double> recv, const Vector<double>& send)
{
  qassert(recv.size() == send.size());
#ifdef USE_MULTI_NODE
  return MPI_Allreduce((double*)send.data(), recv.data(), recv.size(),
                       MPI_DOUBLE, MPI_SUM, get_comm());
#else
  memmove(recv.data(), send.data(), recv.size() * sizeof(double));
  return 0;
#endif
}

inline int glb_sum(Vector<float> recv, const Vector<float>& send)
{
  qassert(recv.size() == send.size());
#ifdef USE_MULTI_NODE
  return MPI_Allreduce((float*)send.data(), recv.data(), recv.size(), MPI_FLOAT,
                       MPI_SUM, get_comm());
#else
  memmove(recv.data(), send.data(), recv.size() * sizeof(float));
  return 0;
#endif
}

inline int glb_sum(Vector<Complex> recv, const Vector<Complex>& send)
{
  return glb_sum(Vector<double>((double*)recv.data(), recv.size() * 2),
                 Vector<double>((double*)send.data(), send.size() * 2));
}

inline int glb_sum(Vector<ComplexF> recv, const Vector<ComplexF>& send)
{
  return glb_sum(Vector<float>((float*)recv.data(), recv.size() * 2),
                 Vector<float>((float*)send.data(), send.size() * 2));
}

inline int glb_sum(Vector<long> recv, const Vector<long>& send)
{
  qassert(recv.size() == send.size());
#ifdef USE_MULTI_NODE
  return MPI_Allreduce((long*)send.data(), recv.data(), recv.size(), MPI_LONG,
                       MPI_SUM, get_comm());
#else
  memmove(recv.data(), send.data(), recv.size() * sizeof(long));
  return 0;
#endif
}

inline int glb_sum(Vector<char> recv, const Vector<char>& send)
{
  qassert(recv.size() == send.size());
#ifdef USE_MULTI_NODE
  return MPI_Allreduce((char*)send.data(), (char*)recv.data(), recv.size(),
                       MPI_BYTE, MPI_BXOR, get_comm());
#else
  memmove(recv.data(), send.data(), recv.data_size());
  return 0;
#endif
}

inline int glb_sum(Vector<double> vec)
{
  if (1 == get_num_node()) {
    return 0;
  }
  std::vector<double> tmp(vec.size());
  assign(tmp, vec);
  return glb_sum(vec, tmp);
}

inline int glb_sum(Vector<float> vec)
{
  if (1 == get_num_node()) {
    return 0;
  }
  std::vector<float> tmp(vec.size());
  assign(tmp, vec);
  return glb_sum(vec, tmp);
}

inline int glb_sum(Vector<Complex> vec)
{
  if (1 == get_num_node()) {
    return 0;
  }
  std::vector<Complex> tmp(vec.size());
  assign(tmp, vec);
  return glb_sum(vec, tmp);
}

inline int glb_sum(Vector<ComplexF> vec)
{
  if (1 == get_num_node()) {
    return 0;
  }
  std::vector<ComplexF> tmp(vec.size());
  assign(tmp, vec);
  return glb_sum(vec, tmp);
}

inline int glb_sum(Vector<long> vec)
{
  if (1 == get_num_node()) {
    return 0;
  }
  std::vector<long> tmp(vec.size());
  assign(tmp, vec);
  return glb_sum(vec, tmp);
}

inline int glb_sum(Vector<char> vec)
{
  if (1 == get_num_node()) {
    return 0;
  }
  std::vector<char> tmp(vec.size());
  assign(tmp, vec);
  return glb_sum(vec, tmp);
}

inline int glb_sum(double& x) { return glb_sum(Vector<double>(x)); }

inline int glb_sum(float& x) { return glb_sum(Vector<float>(x)); }

inline int glb_sum(long& x) { return glb_sum(Vector<long>(x)); }

inline int glb_sum(Complex& c)
{
  return glb_sum(Vector<double>((double*)&c, 2));
}

inline int glb_sum(ComplexF& c)
{
  return glb_sum(Vector<float>((float*)&c, 2));
}

template <class M>
int glb_sum_double_vec(Vector<M> x)
{
  return glb_sum(
      Vector<double>((double*)x.data(), x.data_size() / sizeof(double)));
}

template <class M>
int glb_sum_float_vec(Vector<M> x)
{
  return glb_sum(
      Vector<float>((float*)x.data(), x.data_size() / sizeof(float)));
}

template <class M>
inline int glb_sum_long_vec(Vector<M> x)
{
  return glb_sum(Vector<long>((long*)x.data(), x.data_size() / sizeof(long)));
}

template <class M>
inline int glb_sum_byte_vec(Vector<M> x)
{
  return glb_sum(Vector<char>((char*)x.data(), x.data_size()));
}

template <class M>
inline int glb_sum_double(M& x)
{
  return glb_sum(Vector<double>((double*)&x, sizeof(M) / sizeof(double)));
}

template <class M>
inline int glb_sum_float(M& x)
{
  return glb_sum(Vector<float>((float*)&x, sizeof(M) / sizeof(float)));
}

template <class M>
inline int glb_sum_long(M& x)
{
  return glb_sum(Vector<long>((long*)&x, sizeof(M) / sizeof(long)));
}

template <class M>
inline int glb_sum_byte(M& x)
{
  return glb_sum(Vector<char>((char*)&x, sizeof(M)));
}

inline int glb_sum_lat_data(LatData& ld)
{
  return glb_sum_double_vec(get_data(ld.res));
}

inline int glb_sum(LatData& ld) { return glb_sum_lat_data(ld); }

template <class M>
void all_gather(Vector<M> recv, const Vector<M>& send)
{
  qassert(recv.size() == send.size() * get_num_node());
#ifdef USE_MULTI_NODE
  MPI_Allgather((void*)send.data(), send.data_size(), MPI_BYTE,
                (void*)recv.data(), send.data_size(), MPI_BYTE, get_comm());
#else
  const long sendsize = send.size() * sizeof(M);
  memmove(recv, send, sendsize);
#endif
}

template <class M>
void bcast(Vector<M> recv, const int root = 0)
{
#ifdef USE_MULTI_NODE
  if (1 == get_num_node()) {
    return;
  }
  MPI_Bcast((void*)recv.data(), recv.data_size(), MPI_BYTE, root, get_comm());
#endif
}

inline void bcast(int& x, const int root = 0)
{
  bcast(get_data_one_elem(x), root);
}

inline void bcast(long& x, const int root = 0)
{
  bcast(get_data_one_elem(x), root);
}

inline void bcast(uint32_t& x, const int root = 0)
{
  bcast(get_data_one_elem(x), root);
}

inline void bcast(double& x, const int root = 0)
{
  bcast(get_data_one_elem(x), root);
}

inline void bcast(Coordinate& x, const int root = 0)
{
  bcast(get_data_one_elem(x), root);
}

template <class M>
void bcast(std::vector<M>& recv, const int root = 0)
{
#ifdef USE_MULTI_NODE
  if (1 == get_num_node()) {
    return;
  }
  long size = recv.size();
  bcast(get_data(size), root);
  recv.resize(size);
  bcast(get_data(recv), root);
#endif
}

inline void bcast(std::string& recv, const int root = 0)
{
#ifdef USE_MULTI_NODE
  if (1 == get_num_node()) {
    return;
  }
  long size = recv.size();
  bcast(get_data(size), root);
  recv.resize(size);
  bcast(get_data(recv), root);
#endif
}

inline void bcast(std::vector<std::string>& recv, const int root = 0)
{
#ifdef USE_MULTI_NODE
  if (1 == get_num_node()) {
    return;
  }
  long size = recv.size();
  bcast(get_data(size), root);
  recv.resize(size);
  for (long i = 0; i < size; ++i) {
    bcast(recv[i], root);
  }
#endif
}

inline void bcast(LatData& ld, const int root = 0)
{
#ifdef USE_MULTI_NODE
  if (1 == get_num_node()) {
    return;
  }
  std::string info_str;
  if (get_id_node() == root) {
    info_str = show(ld.info);
  }
  bcast(info_str, root);
  if (get_id_node() != root) {
    ld.info = read_lat_info(info_str);
    lat_data_alloc(ld);
  }
  bcast(get_data(ld.res), root);
#endif
}

template <class M>
void bcast(std::vector<std::vector<M> >& datatable, const int root = 0)
{
#ifdef USE_MULTI_NODE
  if (1 == get_num_node()) {
    return;
  }
  long nrow, total_size;
  std::vector<long> row_sizes;
  std::vector<M> data;
  if (get_id_node() == root) {
    row_sizes = vector_map_size(datatable);
    data = vector_concat(datatable);
    nrow = row_sizes.size();
    total_size = data.size();
  }
  bcast(get_data(nrow), root);
  bcast(get_data(total_size), root);
  if (get_id_node() != root) {
    row_sizes.resize(nrow);
    data.resize(total_size);
  }
  bcast(get_data(row_sizes), root);
  bcast(get_data(data), root);
  if (get_id_node() != root) {
    datatable = vector_split(data, row_sizes);
  }
#endif
}

inline void sync_node()
{
  RngState& rs = get_comm_list().back().sync_node_rs;
  const long v = rand_gen(rs) % (1024 * 1024);
  long s = v;
  glb_sum(s);
  qassert(s == v * get_num_node());
}

inline std::vector<int> mk_id_node_list_for_shuffle_rs(const RngState& rs)
{
  TIMER_VERBOSE("mk_id_node_list_for_shuffle_rs");
  const int num_node = get_num_node();
  std::vector<int> list(num_node);
  for (int i = 0; i < num_node; ++i) {
    list[i] = i;
  }
  random_permute(list, rs);
  for (int i = 0; i < num_node; ++i) {
    if (0 == list[i]) {
      list[i] = list[0];
      list[0] = 0;
      break;
    }
  }
  return list;
}

inline std::vector<int> mk_id_node_list_for_shuffle_step_size(const int step_size_)
{
  TIMER_VERBOSE("mk_id_node_list_for_shuffle_step_size");
  const int num_node = get_num_node();
  const int step_size =
      (step_size_ < num_node and num_node % step_size_ == 0) ? step_size_ : 1;
  std::vector<int> list(num_node);
  for (int i = 0; i < num_node; ++i) {
    const int id_node_in_shuffle = i;
    const int id_node =
        mod(i * step_size, num_node) + (i * step_size / num_node);
    list[id_node_in_shuffle] = id_node;
  }
  return list;
}

inline std::vector<int> mk_id_node_list_for_shuffle_node()
// return list
// list[id_node_in_shuffle] = id_node
{
  TIMER_VERBOSE("mk_id_node_list_for_shuffle_node");
  // global id
  int globalRank;
  MPI_Comm_rank(get_comm(), &globalRank);
  qassert(globalRank == get_id_node());
  // node local comm
  MPI_Comm nodeComm;
  MPI_Comm_split_type(get_comm(), MPI_COMM_TYPE_SHARED, globalRank,
                      MPI_INFO_NULL, &nodeComm);
  // id within the node
  int localRank;
  MPI_Comm_rank(nodeComm, &localRank);
  if (0 == get_id_node()) {
    qassert(localRank == 0);
  }
  // number of process in this node
  int localSize;
  MPI_Comm_size(nodeComm, &localSize);
  // comm across node (each node select one process with the same local rank)
  MPI_Comm masterComm;
  MPI_Comm_split(get_comm(), localRank, globalRank, &masterComm);
  // id across node
  int masterRank;
  MPI_Comm_rank(masterComm, &masterRank);
  // size of each master comm
  int masterSize;
  MPI_Comm_size(masterComm, &masterSize);
  // calculate number of node
  long num_of_node = masterSize;
  MPI_Bcast(&num_of_node, 1, MPI_LONG, 0, nodeComm);
  // calculate id of node (master rank of the 0 local rank process)
  long id_of_node = masterRank;
  MPI_Bcast(&id_of_node, 1, MPI_LONG, 0, nodeComm);
  qassert(id_of_node < num_of_node);
  // calculate number of processes for each node
  std::vector<long> num_process_for_each_node(num_of_node, 0);
  num_process_for_each_node[id_of_node] = 1;
  glb_sum(get_data(num_process_for_each_node));
  qassert(num_process_for_each_node[id_of_node] == localSize);
  // calculate the number of master comm (the maximum in num_process_for_each_node)
  long num_of_master_comm = 0;
  for (long i = 0; i < (long)num_process_for_each_node.size(); ++i) {
    if (num_process_for_each_node[i] > num_of_master_comm) {
      num_of_master_comm = num_process_for_each_node[i];
    }
  }
  // calculate the id of the master comm (same as local rank)
  long id_of_master_comm = localRank;
  qassert(id_of_master_comm < num_of_master_comm);
  // calculate number of processes for each masterComm
  std::vector<long> num_process_for_each_master_comm(num_of_master_comm, 0);
  num_process_for_each_master_comm[id_of_master_comm] = 1;
  glb_sum(get_data(num_process_for_each_master_comm));
  qassert(num_process_for_each_master_comm[id_of_master_comm] == masterSize);
  // calculate id_node_in_shuffle
  long id_node_in_shuffle = masterRank;
  for (long i = 0; i < id_of_master_comm; ++i) {
    id_node_in_shuffle += num_process_for_each_master_comm[i];
  }
  // calculate the list of id_node for each id_node_in_shuffle
  std::vector<long> list_long(get_num_node(), 0);
  list_long[id_node_in_shuffle] = get_id_node();
  glb_sum(get_data(list_long));
  std::vector<int> list(get_num_node(), 0);
  for (long i = 0; i < get_num_node(); ++i) {
    list[i] = list_long[i];
  }
  // checking
  qassert(list[0] == 0);
  for (long i = 0; i < get_num_node(); ++i) {
    qassert(0 <= list[i]);
    qassert(list[i] < get_num_node());
    for (long j = 0; j < i; ++j) {
      qassert(list[i] != list[j]);
    }
  }
  return list;
}

inline std::vector<int> mk_id_node_list_for_shuffle()
// use env variable "q_mk_id_node_in_shuffle_seed"
// if env variable start with "seed_", then the rest will be used as seed for random assignment
// else env variable will be viewed as int for step_size
{
  TIMER_VERBOSE("mk_id_node_list_for_shuffle")
  const std::string seed = get_env_default("q_mk_id_node_in_shuffle_seed", "");
  const std::string seed_prefix = "seed_";
  if (seed == "") {
    return mk_id_node_list_for_shuffle_node();
  } else if (seed.compare(0, seed_prefix.size(), seed_prefix) == 0) {
    RngState rs(seed.substr(seed_prefix.size()));
    return mk_id_node_list_for_shuffle_rs(rs);
  } else {
    const long step_size = read_long(seed);
    return mk_id_node_list_for_shuffle_step_size(step_size);
  }
}

API inline std::vector<int>& get_id_node_list_for_shuffle()
// qlat parameter
// initialized in begin_comm with mk_id_node_list_for_shuffle()
// return list
// list[id_node_in_shuffle] = id_node
{
  static std::vector<int> list;
  return list;
}

inline std::vector<int> mk_id_node_in_shuffle_list()
// return list_new
// list_new[id_node] = id_node_in_shuffle
{
  TIMER_VERBOSE("mk_id_node_in_shuffle_list")
  const std::vector<int>& list = get_id_node_list_for_shuffle();
  const int num_node = list.size();
  qassert(num_node == get_num_node());
  std::vector<int> list_new(num_node, 0);
  for (int i = 0; i < num_node; ++i) {
    const int id_node_in_shuffle = i;
    const int id_node = list[i];
    qassert(0 <= id_node_in_shuffle);
    qassert(id_node_in_shuffle < num_node);
    qassert(0 <= id_node);
    qassert(id_node < num_node);
    list_new[id_node] = id_node_in_shuffle;
  }
  return list_new;
}

API inline std::vector<int>& get_id_node_in_shuffle_list()
// qlat parameter
// initialized in begin_comm with mk_id_node_in_shuffle_list()
// return list
// list[id_node] = id_node_in_shuffle
{
  static std::vector<int> list;
  return list;
}

inline int get_id_node_in_shuffle(const int id_node, const int new_num_node,
                                  const int num_node)
// not called very often
{
  qassert(0 <= id_node);
  qassert(id_node < num_node);
  if (new_num_node == num_node) {
    return id_node;
  } else {
    const std::vector<int>& list = get_id_node_in_shuffle_list();
    qassert((long)list.size() == num_node);
    qassert(list[0] == 0);
    return list[id_node];
  }
}

inline int get_id_node_from_id_node_in_shuffle(const int id_node_in_shuffle,
                                               const int new_num_node,
                                               const int num_node)
{
  qassert(0 <= id_node_in_shuffle);
  qassert(id_node_in_shuffle < num_node);
  if (new_num_node == num_node) {
    return id_node_in_shuffle;
  } else {
    const std::vector<int>& list = get_id_node_list_for_shuffle();
    qassert((long)list.size() == num_node);
    qassert(list[0] == 0);
    return list[id_node_in_shuffle];
  }
}

inline std::string get_hostname()
{
  char name[MPI_MAX_PROCESSOR_NAME];
  int len;
  MPI_Get_processor_name(name, &len);
  return std::string(name, len);
}

inline void display_geometry_node()
{
  TIMER("display_geometry_node");
  const GeometryNode& geon = get_geometry_node();
  for (int i = 0; i < geon.num_node; ++i) {
    if (i == geon.id_node) {
      displayln(std::string(fname) + " : " +
                ssprintf("id_node = %5d ; coor_node = %s ; id_node_in_shuffle "
                         "= %5d ; hostname = %s",
                         geon.id_node, show(geon.coor_node).c_str(),
                         get_id_node_in_shuffle(), get_hostname().c_str()));
      fflush(get_output_file());
    }
    sync_node();
  }
  fflush(get_output_file());
  sync_node();
}

inline Coordinate plan_size_node(const int num_node)
{
  // assuming MPI is initialized ...
  int dims[] = {0, 0, 0, 0};
  MPI_Dims_create(num_node, DIMN, dims);
  return Coordinate(dims[3], dims[2], dims[1], dims[0]);
}

inline bool is_MPI_initialized()
{
  int b;
  MPI_Initialized(&b);
  return b;
}

inline int init_mpi(int* argc, char** argv[])
{
  if (!is_MPI_initialized()) MPI_Init(argc, argv);
  int num_node;
  MPI_Comm_size(MPI_COMM_WORLD, &num_node);
  int id_node;
  MPI_Comm_rank(MPI_COMM_WORLD, &id_node);
  if (0 == id_node) {
    displayln("qlat::begin(): " +
              ssprintf("MPI Initialized. num_node = %d", num_node));
  }
  return num_node;
}

inline void set_global_geon(const Coordinate& size_node)
{
  int num_node;
  MPI_Comm_size(get_comm(), &num_node);
  qassert(num_node == product(size_node));
  int id_node;
  MPI_Comm_rank(get_comm(), &id_node);
  GeometryNode& geon = get_geometry_node_internal();
  geon.init(id_node, size_node);
  get_id_node_internal() = geon.id_node;
  get_num_node_internal() = geon.num_node;
  qassert(geon.num_node == num_node);
}

inline void begin_comm(const MPI_Comm comm, const Coordinate& size_node)
// begin Qlat with existing comm (assuming MPI already initialized)
{
  get_comm_list().push_back(
      Q_Comm(comm, size_node, RngState("sync_node:" + show(size_node))));
  get_comm_internal() = get_comm_list().back().comm;
  set_global_geon(get_comm_list().back().size_node);
  sync_node();
  displayln_info(ssprintf(
      "qlat::begin_comm(comm,size_node): get_comm_list().push_back()"));
  displayln_info(
      ssprintf("qlat::begin_comm(comm,size_node): get_comm_list().size() = %d",
               (int)get_comm_list().size()));
  const GeometryNode& geon = get_geometry_node();
  if (get_env("OMP_NUM_THREADS") == "") {
    const long num_threads = get_env_long_default("q_num_threads", 2);
    omp_set_num_threads(num_threads);
  }
  displayln_info("qlat::begin(): q_num_threads = " +
                 show(omp_get_max_threads()));
#ifndef QLAT_NO_MALLOPT
  std::string q_malloc_mmap_threshold =
      get_env_default("q_malloc_mmap_threshold", "");
  if (q_malloc_mmap_threshold != "") {
    mallopt(M_MMAP_THRESHOLD, read_long(q_malloc_mmap_threshold));
  }
#endif
  displayln_info("qlat::begin(): GeometryNode =\n" + show(geon));
  fflush(get_output_file());
  displayln_info(ssprintf("Timer::get_timer_database().size() = %ld",
                          Timer::get_timer_database().size()));
  displayln_info(ssprintf("Timer::get_timer_stack().size() = %ld",
                          Timer::get_timer_stack().size()));
  get_id_node_list_for_shuffle() = mk_id_node_list_for_shuffle();
  get_id_node_in_shuffle_list() = mk_id_node_in_shuffle_list();
  get_id_node_in_shuffle_internal() =
      get_id_node_in_shuffle(get_id_node(), 0, get_num_node());
  // display_geometry_node();
  // install_qhandle_sig();
  clear_all_caches();
  sync_node();
}

inline void begin(const int id_node, const Coordinate& size_node,
                  const int color = 0)
// begin Qlat with existing id_node mapping (assuming MPI already initialized)
{
  if (get_comm_list().empty()) {
    get_comm_list().push_back(
        Q_Comm(MPI_COMM_WORLD, Coordinate(), RngState("sync_node")));
  }
  qassert(0 <= id_node and id_node < product(size_node));
  qassert(0 <= color);
  MPI_Comm comm;
  const int ret =
      MPI_Comm_split(get_comm_list().back().comm, color, id_node, &comm);
  qassert(ret == MPI_SUCCESS);
  begin_comm(comm, size_node);
  qassert(get_id_node() == id_node);
}

inline void begin(int* argc, char** argv[], const Coordinate& size_node)
// not recommended
{
  const int num_node = init_mpi(argc, argv);
  qassert(num_node == product(size_node));
  begin_comm(MPI_COMM_WORLD, size_node);
}

inline void begin(
    int* argc, char** argv[],
    const std::vector<Coordinate>& size_node_list = std::vector<Coordinate>())
// begin Qlat and initialize a new comm
{
  const int num_node = init_mpi(argc, argv);
  Coordinate size_node;
  for (int i = 0; i < (int)size_node_list.size(); ++i) {
    size_node = size_node_list[i];
    if (num_node == product(size_node)) {
      break;
    }
  }
  if (num_node != product(size_node)) {
    size_node = plan_size_node(num_node);
  }
  begin_comm(MPI_COMM_WORLD, size_node);
}

inline void end(const bool is_preserving_cache = false)
{
  if (get_comm_list().empty()) {
    displayln_info(ssprintf("qlat::end(): get_comm_list().empty() = true."));
    if (not is_preserving_cache) {
      clear_all_caches();
    }
  } else {
    qassert(get_comm_list().back().comm == get_comm());
    if (get_comm() == MPI_COMM_WORLD) {
      if (not is_preserving_cache) {
        clear_all_caches();
      }
      sync_node();
      displayln_info(ssprintf("qlat::end(): get_comm_list().pop_back()"));
      get_comm_list().pop_back();
      displayln_info(ssprintf("qlat::end(): get_comm_list().size() = %d.",
                              (int)get_comm_list().size()));
      qassert(get_comm_list().size() == 0);
      displayln_info("qlat::end(): Finalize MPI.");
      if (is_MPI_initialized()) MPI_Finalize();
      displayln_info("qlat::end(): MPI Finalized.");
    } else {
      if (not is_preserving_cache) {
        clear_all_caches();
      }
      sync_node();
      displayln_info(ssprintf("qlat::end(): get_comm_list().pop_back()"));
      MPI_Comm comm = get_comm();
      get_comm_list().pop_back();
      displayln_info(ssprintf("qlat::end(): get_comm_list().size() = %d.",
                              (int)get_comm_list().size()));
      get_comm_internal() = get_comm_list().back().comm;
      if (get_comm_list().back().size_node != Coordinate()) {
        displayln_info(
            ssprintf("qlat::end(): Switch to old comm and setup global geon."));
        set_global_geon(get_comm_list().back().size_node);
      } else {
        displayln_info(ssprintf("qlat::end(): Switch to old comm (foreign)."));
      }
      displayln_info(ssprintf("qlat::end(): MPI_Comm_free ended comm."));
      MPI_Comm_free(&comm);
      sync_node();
    }
  }
}

}  // namespace qlat
