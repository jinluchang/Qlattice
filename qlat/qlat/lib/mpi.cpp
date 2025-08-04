#include <qlat-utils/core.h>
#include <qlat/mpi.h>

namespace qlat
{

int mpi_send(const void* buf, Long count, MPI_Datatype datatype, int dest,
             int tag, MPI_Comm comm)
{
  const Long int_max = INT_MAX;
  if (count <= int_max) {
    return MPI_Send(buf, count, datatype, dest, tag, comm);
  } else {
    Int type_size = 0;
    MPI_Type_size(datatype, &type_size);
    Char* cbuf = (Char*)buf;
    while (count > int_max) {
      MPI_Send(cbuf, int_max, datatype, dest, tag, comm);
      cbuf += int_max * type_size;
      count -= int_max;
    }
    return MPI_Send(cbuf, count, datatype, dest, tag, comm);
  }
}

int mpi_recv(void* buf, Long count, MPI_Datatype datatype, int source, int tag,
             MPI_Comm comm, MPI_Status* status)
{
  const Long int_max = INT_MAX;
  if (count <= int_max) {
    return MPI_Recv(buf, count, datatype, source, tag, comm, status);
  } else {
    qassert(status == MPI_STATUS_IGNORE);
    Int type_size = 0;
    MPI_Type_size(datatype, &type_size);
    Char* cbuf = (Char*)buf;
    while (count > int_max) {
      MPI_Recv(cbuf, int_max, datatype, source, tag, comm, MPI_STATUS_IGNORE);
      cbuf += int_max * type_size;
      count -= int_max;
    }
    return MPI_Recv(cbuf, count, datatype, source, tag, comm,
                    MPI_STATUS_IGNORE);
  }
}

int mpi_isend(const void* buf, Long count, MPI_Datatype datatype, int dest,
              int tag, MPI_Comm comm, std::vector<MPI_Request>& requests)
{
  const Long int_max = INT_MAX;
  if (count <= int_max) {
    MPI_Request r;
    int ret = MPI_Isend(buf, count, datatype, dest, tag, comm, &r);
    requests.push_back(r);
    qassert(ret == MPI_SUCCESS);
    return MPI_SUCCESS;
  } else {
    Int type_size = 0;
    MPI_Type_size(datatype, &type_size);
    Char* cbuf = (Char*)buf;
    while (count > int_max) {
      mpi_isend(cbuf, int_max, datatype, dest, tag, comm, requests);
      cbuf += int_max * type_size;
      count -= int_max;
    }
    return mpi_isend(cbuf, count, datatype, dest, tag, comm, requests);
  }
}

int mpi_irecv(void* buf, Long count, MPI_Datatype datatype, int source, int tag,
              MPI_Comm comm, std::vector<MPI_Request>& requests)
{
  const Long int_max = INT_MAX;
  if (count <= int_max) {
    MPI_Request r;
    int ret = MPI_Irecv(buf, count, datatype, source, tag, comm, &r);
    requests.push_back(r);
    qassert(ret == MPI_SUCCESS);
    return MPI_SUCCESS;
  } else {
    Int type_size = 0;
    MPI_Type_size(datatype, &type_size);
    Char* cbuf = (Char*)buf;
    while (count > int_max) {
      mpi_irecv(cbuf, int_max, datatype, source, tag, comm, requests);
      cbuf += int_max * type_size;
      count -= int_max;
    }
    return mpi_irecv(cbuf, count, datatype, source, tag, comm, requests);
  }
}

int mpi_waitall(std::vector<MPI_Request>& requests)
{
  TIMER("mpi_waitall");
  if (requests.size() > 0) {
    int ret = MPI_Waitall(requests.size(), requests.data(), MPI_STATUS_IGNORE);
    qassert(ret == MPI_SUCCESS);
  }
  requests.resize(0);
  return MPI_SUCCESS;
}

static const std::vector<Int>& get_random_rank_order(const Int size)
{
  static std::vector<Int> order;
  if ((Int)order.size() == size) {
    return order;
  }
  order.resize(size);
  for (Int i = 0; i < size; ++i) {
    order[i] = i;
  }
  random_permute(order, RngState("get_random_rank_order"));
  return order;
}

static int mpi_alltoallv_custom(const void* sendbuf, const Long* sendcounts,
                                const Long* sdispls, MPI_Datatype sendtype,
                                void* recvbuf, const Long* recvcounts,
                                const Long* rdispls, MPI_Datatype recvtype,
                                MPI_Comm comm)
// Send and Recv data as MPI_BYTE of corresponding type size.
{
  TIMER_FLOPS("mpi_alltoallv_custom");
  //
  static const Int q_mpi_alltoallv_max_parallel_transfer =
      get_env_long_default("q_mpi_alltoallv_max_parallel_transfer", 128);
  //
  Int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  //
  const std::vector<Int>& order = get_random_rank_order(size);
  Int rank_idx = -1;
  for (Int i = 0; i < size; ++i) {
    if (order[i] == rank) {
      rank_idx = i;
      break;
    }
  }
  qassert(rank_idx >= 0);
  qassert(order[rank_idx] == rank);
  //
  const Int mpi_tag = 13;
  // 计算数据类型大小
  Int sendtype_size_i, recvtype_size_i;
  MPI_Type_size(sendtype, &sendtype_size_i);
  MPI_Type_size(recvtype, &recvtype_size_i);
  const Long sendtype_size = sendtype_size_i;
  const Long recvtype_size = recvtype_size_i;
  //
  // order[rank_idx] == rank
  // where
  // rank_idx is the virtual rank in this all toall transfer
  // rank is the actual MPI rank
  //
  std::vector<MPI_Request> requests;
  for (Int i = 0; i < size; i += q_mpi_alltoallv_max_parallel_transfer) {
    // 非阻塞接收阶段
    for (Int j = 0; j < q_mpi_alltoallv_max_parallel_transfer; j++) {
      if (i + j >= size) {
        continue;
      }
      const Int rank_from = order[(rank_idx + size - i - j) % size];
      if (recvcounts[rank_from] > 0) {
        char* recv_ptr = (char*)recvbuf + rdispls[rank_from] * recvtype_size;
        mpi_irecv(recv_ptr, recvcounts[rank_from] * recvtype_size, MPI_BYTE,
                  rank_from, mpi_tag, comm, requests);
        timer.flops += recvcounts[rank_from] * recvtype_size / 2;
      }
    }
    // 非阻塞发送阶段
    for (Int j = 0; j < q_mpi_alltoallv_max_parallel_transfer; j++) {
      if (i + j >= size) {
        continue;
      }
      const Int rank_to = order[(rank_idx + i + j) % size];
      if (sendcounts[rank_to] > 0) {
        const char* send_ptr =
            (const char*)sendbuf + sdispls[rank_to] * sendtype_size;
        mpi_isend(send_ptr, sendcounts[rank_to] * sendtype_size, MPI_BYTE,
                  rank_to, mpi_tag, comm, requests);
        timer.flops += sendcounts[rank_to] * sendtype_size / 2;
      }
    }
    // 等待所有通信完成
    mpi_waitall(requests);
  }
  return 0;
}

static int mpi_alltoallv_native(const void* sendbuf, const Long* sendcounts,
                                 const Long* sdispls, MPI_Datatype sendtype,
                                 void* recvbuf, const Long* recvcounts,
                                 const Long* rdispls, MPI_Datatype recvtype,
                                 MPI_Comm comm)
{
  TIMER_FLOPS("mpi_alltoallv_native");
  Int sendtype_size_i, recvtype_size_i;
  MPI_Type_size(sendtype, &sendtype_size_i);
  MPI_Type_size(recvtype, &recvtype_size_i);
  const Long sendtype_size = sendtype_size_i;
  const Long recvtype_size = recvtype_size_i;
  Int size;
  MPI_Comm_size(comm, &size);
  vector<Int> sendcounts_i(size);
  vector<Int> sdispls_i(size);
  vector<Int> recvcounts_i(size);
  vector<Int> rdispls_i(size);
  for (Int k = 0; k < size; ++k) {
    qassert(sendcounts[k] < INT_MAX);
    qassert(sdispls[k] < INT_MAX);
    qassert(recvcounts[k] < INT_MAX);
    qassert(rdispls[k] < INT_MAX);
    sendcounts_i[k] = sendcounts[k];
    sdispls_i[k] = sdispls[k];
    recvcounts_i[k] = recvcounts[k];
    rdispls_i[k] = rdispls[k];
    timer.flops += sendcounts[k] * sendtype_size / 2;
    timer.flops += recvcounts[k] * recvtype_size / 2;
  }
  return MPI_Alltoallv(sendbuf, sendcounts_i.data(), sdispls_i.data(), sendtype,
                       recvbuf, recvcounts_i.data(), rdispls_i.data(), recvtype,
                       comm);
}

int mpi_alltoallv(const void* sendbuf, const Long* sendcounts,
                  const Long* sdispls, MPI_Datatype sendtype, void* recvbuf,
                  const Long* recvcounts, const Long* rdispls,
                  MPI_Datatype recvtype, MPI_Comm comm)
{
  static const std::string q_env_type =
      get_env_default("q_mpi_alltoallv_type", "custom");
  if (q_env_type == "custom") {
    return mpi_alltoallv_custom(sendbuf, sendcounts, sdispls, sendtype, recvbuf,
                                recvcounts, rdispls, recvtype, comm);
  } else if (q_env_type == "native") {
    return mpi_alltoallv_native(sendbuf, sendcounts, sdispls, sendtype, recvbuf,
                                recvcounts, rdispls, recvtype, comm);
  } else {
    qerr(ssprintf("mpi_alltoallv: q_mpi_alltoallv_type='%s'.",
                  q_env_type.c_str()));
  }
  return -1;
}

static int mpi_bcast_native(void* buffer, const Long count,
                            MPI_Datatype datatype, const int root,
                            MPI_Comm comm)
{
  TIMER_FLOPS("mpi_bcast_native");
  Int type_size = 0;
  MPI_Type_size(datatype, &type_size);
  timer.flops += count * type_size;
  qassert(count < INT_MAX);
  const Int count_i = count;
  return MPI_Bcast(buffer, count_i, datatype, root, comm);
}

static int mpi_bcast_custom(void* buffer, const Long count,
                            MPI_Datatype datatype, const int root,
                            MPI_Comm comm)
{
  TIMER_FLOPS("mpi_bcast_custom");
  Int type_size = 0;
  MPI_Type_size(datatype, &type_size);
  timer.flops += count * type_size;
  Int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  const int mpi_tag = 13;
  if (rank == root) {
    // Root process sends to all other processes
    for (int i = 0; i < size; i++) {
      if (i != root) {
        mpi_send(buffer, count, datatype, i, mpi_tag, comm);
      }
    }
  } else {
    // Other processes receive from root
    mpi_recv(buffer, count, datatype, root, mpi_tag, comm);
  }
  return 0;
}

int mpi_bcast(void* buffer, const Long count, MPI_Datatype datatype,
              const int root, MPI_Comm comm)
{
  static const std::string q_env_type =
      get_env_default("q_mpi_bcast", "custom");
  if (q_env_type == "custom") {
    return mpi_bcast_custom(buffer, count, datatype, root, comm);
  } else if (q_env_type == "native") {
    return mpi_bcast_native(buffer, count, datatype, root, comm);
  } else {
    qerr(ssprintf("mpi_bcast: q_mpi_bcast='%s'.", q_env_type.c_str()));
  }
  return -1;
}

static int mpi_allreduce_custom(const void* sendbuf, void* recvbuf,
                                const Long count, MPI_Datatype datatype,
                                MPI_Op op, MPI_Comm comm)
{
  TIMER_FLOPS("mpi_allreduce_custom");
  Int type_size;
  MPI_Type_size(datatype, &type_size);
  timer.flops += count * type_size;
  //
  Int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  // First copy sendbuf to recvbuf
  //
  std::memcpy(recvbuf, sendbuf, count * type_size);
  // Tree-based reduction to root (rank 0)
  int mask = 1;
  while (mask < size) {
    int partner = rank ^ mask;
    if (partner < size) {
      vector<Char> temp_vec(count * type_size, MemType::Comm);
      void* tempbuf = temp_vec.data();
      if (rank < partner) {
        mpi_recv(tempbuf, count, datatype, partner, 0, comm);
        // Perform sum operation based on datatype
        if (datatype == MPI_INT and op == MPI_SUM) {
          int* rbuf = (int*)recvbuf;
          int* tbuf = (int*)tempbuf;
          for (Long i = 0; i < count; i++) {
            rbuf[i] += tbuf[i];
          }
        } else if (datatype == MPI_LONG and op == MPI_SUM) {
          long* rbuf = (long*)recvbuf;
          long* tbuf = (long*)tempbuf;
          for (Long i = 0; i < count; i++) {
            rbuf[i] += tbuf[i];
          }
        } else if (datatype == MPI_FLOAT and op == MPI_SUM) {
          float* rbuf = (float*)recvbuf;
          float* tbuf = (float*)tempbuf;
          for (Long i = 0; i < count; i++) {
            rbuf[i] += tbuf[i];
          }
        } else if (datatype == MPI_DOUBLE and op == MPI_SUM) {
          double* rbuf = (double*)recvbuf;
          double* tbuf = (double*)tempbuf;
          for (Long i = 0; i < count; i++) {
            rbuf[i] += tbuf[i];
          }
        } else if (datatype == MPI_CHAR and op == MPI_BXOR) {
          char* rbuf = (char*)recvbuf;
          char* tbuf = (char*)tempbuf;
          for (Long i = 0; i < count; i++) {
            rbuf[i] ^= tbuf[i];
          }
        }
      } else {
        mpi_send(recvbuf, count, datatype, partner, 0, comm);
        break;
      }
    }
    mask <<= 1;
  }
  // Broadcast the result from root to all processes
  return mpi_bcast(recvbuf, count, datatype, 0, comm);
}

static int mpi_allreduce_native(void* sendbuf, void* recvbuf, const Long count,
                                MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  TIMER_FLOPS("mpi_allreduce_native");
  Int type_size;
  MPI_Type_size(datatype, &type_size);
  timer.flops += count * type_size;
  //
  const Long int_max = INT_MAX;
  qassert(count < int_max);
  return MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
}

int mpi_allreduce(void* sendbuf, void* recvbuf, const Long count,
                  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  static const std::string q_env_type =
      get_env_default("q_mpi_allreduce", "custom");
  if (q_env_type == "custom") {
    return mpi_allreduce_custom(sendbuf, recvbuf, count, datatype, op, comm);
  } else if (q_env_type == "native") {
    return mpi_allreduce_native(sendbuf, recvbuf, count, datatype, op, comm);
  } else {
    qerr(ssprintf("mpi_allreduce: q_mpi_allreduce='%s'.", q_env_type.c_str()));
  }
  return -1;
}

int glb_sum(Vector<RealD> recv, const Vector<RealD>& send)
{
  TIMER("glb_sum(RealD)");
  qassert(recv.size() == send.size());
  return mpi_allreduce((RealD*)send.data(), recv.data(), recv.size(),
                       MPI_DOUBLE, MPI_SUM, get_comm());
}

int glb_sum(Vector<RealF> recv, const Vector<RealF>& send)
{
  TIMER("glb_sum(RealF)");
  qassert(recv.size() == send.size());
  return mpi_allreduce((RealF*)send.data(), recv.data(), recv.size(), MPI_FLOAT,
                       MPI_SUM, get_comm());
}

int glb_sum(Vector<Long> recv, const Vector<Long>& send)
{
  TIMER("glb_sum(Long)");
  qassert(recv.size() == send.size());
  return mpi_allreduce((Long*)send.data(), recv.data(), recv.size(),
                       MPI_INT64_T, MPI_SUM, get_comm());
}

int glb_sum(Vector<Int> recv, const Vector<Int>& send)
{
  TIMER("glb_sum(Int)");
  qassert(recv.size() == send.size());
  return mpi_allreduce((Int*)send.data(), recv.data(), recv.size(), MPI_INT32_T,
                       MPI_SUM, get_comm());
}

int glb_sum(Vector<Char> recv, const Vector<Char>& send)
// not SUM but BXOR
{
  TIMER("glb_sum(Char)");
  qassert(recv.size() == send.size());
  return mpi_allreduce((char*)send.data(), (char*)recv.data(), recv.size(),
                       MPI_BYTE, MPI_BXOR, get_comm());
}

int glb_sum(Vector<char> recv, const Vector<char>& send)
// not SUM but BXOR
{
  TIMER("glb_sum(char)");
  qassert(recv.size() == send.size());
  return mpi_allreduce((char*)send.data(), (char*)recv.data(), recv.size(),
                       MPI_BYTE, MPI_BXOR, get_comm());
}

bool glb_any(const bool b)
{
  Long ret = 0;
  if (b) {
    ret = 1;
  }
  glb_sum(ret);
  return ret > 0;
}

bool glb_all(const bool b)
{
  Long ret = 0;
  if (not b) {
    ret = 1;
  }
  glb_sum(ret);
  return ret == 0;
}

int bcast(Vector<Char> recv, const int root)
{
  TIMER("bcast(Vector<Char>)");
  if (1 == get_num_node()) {
    return 0;
  }
  return mpi_bcast((void*)recv.data(), recv.size(), MPI_BYTE, root, get_comm());
}

int bcast(std::string& recv, const int root)
{
  TIMER("bcast(std::string&)");
  if (1 == get_num_node()) {
    return 0;
  }
  int ret = 0;
  Long size = recv.size();
  ret += bcast(get_data(size), root);
  recv.resize(size);
  ret += bcast(get_data(recv), root);
  return ret;
}

int bcast(std::vector<std::string>& recv, const int root)
{
  TIMER("bcast(std::vector<std::string>&)");
  if (1 == get_num_node()) {
    return 0;
  }
  int ret = 0;
  Long size = recv.size();
  ret += bcast(get_data(size), root);
  recv.resize(size);
  for (Long i = 0; i < size; ++i) {
    ret += bcast(recv[i], root);
  }
  return ret;
}

int bcast(PointsSelection& psel, const int root)
{
  TIMER("bcast(PointsSelection&)");
  if (1 == get_num_node()) {
    return 0;
  }
  if (get_id_node() == root) {
    qassert(psel.initialized);
    qassert(psel.points_dist_type == PointsDistType::Global);
  } else {
    psel.initialized = true;
    psel.points_dist_type = PointsDistType::Global;
  }
  int ret = 0;
  ret += bcast(psel.total_site, root);
  ret += bcast(psel.xgs, root);
  return ret;
}

Int bcast_any(Vector<Char> xx, const bool b)
// bcast to all nodes from any node if `b == true`.
// `glb_all(b)` should be `true`.
// The sizes of `xx` should be the same even when `b == false`.
// The value of `xx` when `b == true`, should be the same.
// If the condition is not met, this function will return `-1`.
{
  Int code = 0;
  const Int num_node = get_num_node();
  Char bc = 0;
  if (b) {
    bc = 1;
  }
  vector<Char> all_b(num_node, MemType::Comm);
  code = all_gather(all_b, bc);
  if (code != 0) {
    return code;
  }
  bool b_any = false;
  for (Int i = 0; i < num_node; ++i) {
    if (all_b[i] == 1) {
      b_any = true;
      break;
    }
  }
  if (not b_any) {
    return -1;
  }
  const Long size = xx.size();
  vector<Long> all_size(num_node, MemType::Comm);
  code = all_gather(all_size, size);
  if (code != 0) {
    return code;
  }
  for (Int i = 0; i < num_node; ++i) {
    if (all_size[i] != size) {
      return -2;
    }
  }
  vector<Char> all_v(size * num_node, MemType::Comm);
  code = all_gather(all_v, xx);
  if (code != 0) {
    return code;
  }
  b_any = false;
  for (Int i = 0; i < num_node; ++i) {
    if (all_b[i] == 1) {
      const Vector<Char> v(all_v.data() + i * size, size);
      if (b_any == false) {
        assign(xx, v);
        b_any = true;
      } else {
        if (xx != v) {
          return -3;
        }
      }
    }
  }
  return 0;
}

Int all_gather(Vector<Char> recv, const Vector<Char> send)
{
  qassert(recv.size() == send.size() * get_num_node());
  return MPI_Allgather((void*)send.data(), send.data_size(), MPI_BYTE,
                       (void*)recv.data(), send.data_size(), MPI_BYTE,
                       get_comm());
}

std::vector<Int> mk_id_node_list_for_shuffle_rs(const RngState& rs)
// return list
// list[id_node_in_shuffle] = id_node
// list[0] = 0
{
  TIMER_VERBOSE("mk_id_node_list_for_shuffle_rs");
  const int num_node = get_num_node();
  std::vector<Int> list(num_node);
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

std::vector<Int> mk_id_node_list_for_shuffle_step_size(const int step_size_)
// return list
// list[id_node_in_shuffle] = id_node
// list[0] = 0
{
  TIMER_VERBOSE("mk_id_node_list_for_shuffle_step_size");
  const int num_node = get_num_node();
  const int step_size =
      (step_size_ < num_node and num_node % step_size_ == 0) ? step_size_ : 1;
  std::vector<Int> list(num_node);
  for (int i = 0; i < num_node; ++i) {
    const int id_node_in_shuffle = i;
    const int id_node =
        mod(i * step_size, num_node) + (i * step_size / num_node);
    list[id_node_in_shuffle] = id_node;
  }
  return list;
}

std::vector<Int> mk_id_node_list_for_shuffle_node()
// return list
// list[id_node_in_shuffle] = id_node
// list[0] = 0
// Assign `id_node_in_shuffle` mainly for IO when one physical node may run
// multiple MPI processes. Physical node is identified as `id_of_node`. MPI
// processe is identified as `id_node`. `id_node_in_shuffle` is new id for each
// MPI processes, such that it iterate through physical nodes first.
// Illustration:
// id_of_node & id_node & id_node_in_shuffle
// 0 & 0 & 0
// 0 & 1 & 3
// 1 & 2 & 1
// 1 & 3 & 4
// 2 & 4 & 2
// 2 & 5 & 5
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
  Long num_of_node = masterSize;
  mpi_bcast(&num_of_node, 1, MPI_LONG, 0, nodeComm);
  // calculate id of node (master rank of the 0 local rank process)
  Long id_of_node = masterRank;
  mpi_bcast(&id_of_node, 1, MPI_LONG, 0, nodeComm);
  qassert(id_of_node < num_of_node);
  // calculate number of processes for each node
  std::vector<Long> num_process_for_each_node(num_of_node, 0);
  num_process_for_each_node[id_of_node] = 1;
  glb_sum(get_data(num_process_for_each_node));
  qassert(num_process_for_each_node[id_of_node] == localSize);
  // calculate the number of master comm (the maximum in
  // num_process_for_each_node)
  Long num_of_master_comm = 0;
  for (Long i = 0; i < (Long)num_process_for_each_node.size(); ++i) {
    if (num_process_for_each_node[i] > num_of_master_comm) {
      num_of_master_comm = num_process_for_each_node[i];
    }
  }
  // calculate the id of the master comm (same as local rank)
  Long id_of_master_comm = localRank;
  qassert(id_of_master_comm < num_of_master_comm);
  // calculate number of processes for each masterComm
  std::vector<Long> num_process_for_each_master_comm(num_of_master_comm, 0);
  num_process_for_each_master_comm[id_of_master_comm] = 1;
  glb_sum(get_data(num_process_for_each_master_comm));
  qassert(num_process_for_each_master_comm[id_of_master_comm] == masterSize);
  // calculate id_node_in_shuffle
  Long id_node_in_shuffle = masterRank;
  for (Long i = 0; i < id_of_master_comm; ++i) {
    id_node_in_shuffle += num_process_for_each_master_comm[i];
  }
  // calculate the list of id_node for each id_node_in_shuffle
  std::vector<Long> list_long(get_num_node(), 0);
  list_long[id_node_in_shuffle] = get_id_node();
  glb_sum(get_data(list_long));
  std::vector<Int> list(get_num_node(), 0);
  for (Long i = 0; i < get_num_node(); ++i) {
    list[i] = list_long[i];
  }
  // checking
  qassert(list[0] == 0);
  for (Long i = 0; i < get_num_node(); ++i) {
    qassert(0 <= list[i]);
    qassert(list[i] < get_num_node());
    for (Long j = 0; j < i; ++j) {
      qassert(list[i] != list[j]);
    }
  }
  return list;
}

std::vector<Int> mk_id_node_list_for_shuffle()
// Use env variable "q_mk_id_node_in_shuffle_seed".
// If the env variable is empty, then assignment based on physical node. (Should
// be a good idea.) If env variable start with "seed_", then the rest will be
// used as seed for random assignment. Otherwise, env variable will be viewed as
// int for step_size.
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
    const Long step_size = read_long(seed);
    return mk_id_node_list_for_shuffle_step_size(step_size);
  }
}

std::vector<Int> mk_id_node_in_shuffle_list()
// return list_new
// list_new[id_node] = id_node_in_shuffle
{
  TIMER_VERBOSE("mk_id_node_in_shuffle_list")
  const std::vector<Int>& list = get_id_node_list_for_shuffle();
  const int num_node = list.size();
  qassert(num_node == get_num_node());
  std::vector<Int> list_new(num_node, 0);
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

int get_id_node_in_shuffle(const int id_node, const int new_num_node,
                           const int num_node)
// not called very often
{
  qassert(0 <= id_node);
  qassert(id_node < num_node);
  if (new_num_node == num_node) {
    return id_node;
  } else {
    const std::vector<Int>& list = get_id_node_in_shuffle_list();
    qassert((Long)list.size() == num_node);
    qassert(list[0] == 0);
    return list[id_node];
  }
}

int get_id_node_from_id_node_in_shuffle(const int id_node_in_shuffle,
                                        const int new_num_node,
                                        const int num_node)
{
  qassert(0 <= id_node_in_shuffle);
  qassert(id_node_in_shuffle < num_node);
  if (new_num_node == num_node) {
    return id_node_in_shuffle;
  } else {
    const std::vector<Int>& list = get_id_node_list_for_shuffle();
    qassert((Long)list.size() == num_node);
    qassert(list[0] == 0);
    return list[id_node_in_shuffle];
  }
}

void set_node_rank_size(int& localRank, int& localSize)
{
  // global id
  int globalRank;
  MPI_Comm_rank(get_comm(), &globalRank);
  qassert(globalRank == get_id_node());
  // node local comm
  MPI_Comm nodeComm;
  MPI_Comm_split_type(get_comm(), MPI_COMM_TYPE_SHARED, globalRank,
                      MPI_INFO_NULL, &nodeComm);
  // id within the node
  // int localRank;
  MPI_Comm_rank(nodeComm, &localRank);
  if (0 == get_id_node()) {
    qassert(localRank == 0);
  }
  // number of process in this node
  // int localSize;
  MPI_Comm_size(nodeComm, &localSize);
}

std::string get_hostname()
{
  char name[MPI_MAX_PROCESSOR_NAME];
  int len;
  MPI_Get_processor_name(name, &len);
  return std::string(name, len);
}

void display_geometry_node()
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
      flush();
    }
    sync_node();
  }
  flush();
  sync_node();
}

Coordinate plan_size_node(const int num_node)
{
  // assuming MPI is initialized ...
  int dims[] = {0, 0, 0, 0};
  MPI_Dims_create(num_node, DIMN, dims);
  return Coordinate(dims[3], dims[2], dims[1], dims[0]);
}

bool is_MPI_initialized()
{
  int b;
  MPI_Initialized(&b);
  return b;
}

int init_mpi(int* argc, char** argv[])
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

void set_global_geon(const Coordinate& size_node)
{
  int num_node;
  MPI_Comm_size(get_comm(), &num_node);
  qassert(num_node == product(size_node));
  int id_node;
  MPI_Comm_rank(get_comm(), &id_node);
  GeometryNode& geon = get_geometry_node_internal();
  geon.init(id_node, size_node);
  qassert(geon.num_node == num_node);
  get_id_node_internal() = geon.id_node;
  get_num_node_internal() = geon.num_node;
  get_sync_node_rs_ptr() = &(get_comm_list().back().sync_node_rs);
  get_glb_sum_long_vec_ptr() = glb_sum_long_vec_mpi;
  get_glb_sum_int_vec_ptr() = glb_sum_int_vec_mpi;
  get_glb_sum_real_d_vec_ptr() = glb_sum_real_d_vec_mpi;
  get_glb_sum_real_f_vec_ptr() = glb_sum_real_f_vec_mpi;
  get_glb_sum_byte_vec_ptr() = glb_sum_byte_vec_mpi;
  get_bcast_byte_vec_ptr() = bcast_byte_vec_mpi;
}

void set_cuda_device()
{
#ifdef __NVCC__
  TIMER_VERBOSE("set_cuda_device");
  int local_rank = 0;
  int local_size = 0;
  set_node_rank_size(local_rank, local_size);
  int num_devices = 0;
  qacc_GetDeviceCount(&num_devices);
  if (num_devices > 0) {
    displayln_info(fname +
                   ssprintf(": num_devices=%d (local_rank=%d local_size=%d)",
                            num_devices, local_rank, local_size));
    qacc_SetDevice(local_rank % num_devices);
  }
#endif
}

void display_qlat_banner()
{
  displayln_info(
      "======================================================================");
  displayln_info("");
  displayln_info("                              Qlattice");
  displayln_info("");
  displayln_info("                   Copyright (C) 2023 Luchang Jin");
  displayln_info("");
  displayln_info("                          " + get_qlat_version());
  displayln_info("");
  displayln_info(
      "This program is free software; you can redistribute it and/or modify\n"
      "it under the terms of the GNU General Public License as published by\n"
      "the Free Software Foundation; either version 3 of the License, or\n"
      "(at your option) any later version.\n"
      "\n"
      "This program is distributed in the hope that it will be useful,\n"
      "but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
      "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
      "GNU General Public License for more details.");
  displayln_info("");
  displayln_info(
      "======================================================================");
}

void initialize_qlat_comm()
{
  display_qlat_banner();
  if (get_env("OMP_NUM_THREADS") == "") {
    const Long num_threads = get_env_long_default("q_num_threads", 2);
    omp_set_num_threads(num_threads);
  }
  displayln_info("qlat::begin(): q_num_threads = " +
                 show(omp_get_max_threads()));
  const GeometryNode& geon = get_geometry_node();
  displayln_info("qlat::begin(): GeometryNode =\n" + show(geon));
  displayln_info(ssprintf(
      "qlat::begin_comm(comm,size_node): get_comm_list().push_back()"));
  displayln_info(
      ssprintf("qlat::begin_comm(comm,size_node): get_comm_list().size() = %d",
               (int)get_comm_list().size()));
  // Do not set cuda device
  // Rely on the environment variable
  // Can use the bind-gpu-qlat.sh scripts
  // set_cuda_device();
  qset_line_buf(stdout);
  displayln_info(ssprintf("Timer::get_timer_database().size() = %ld",
                          Timer::get_timer_database().size()));
  displayln_info(ssprintf("Timer::get_timer_stack().size() = %ld",
                          Timer::get_timer_stack().size()));
#ifndef QLAT_NO_MALLOPT
  std::string q_malloc_mmap_threshold =
      get_env_default("q_malloc_mmap_threshold", "");
  if (q_malloc_mmap_threshold != "") {
    mallopt(M_MMAP_THRESHOLD, read_long(q_malloc_mmap_threshold));
  }
#endif
  flush();
  sync_node();
}

Long& mpi_level_count()
{
  static Long c = 0;
  return c;
}

void begin_comm(const MPI_Comm comm, const Coordinate& size_node)
// begin Qlat with existing comm (assuming MPI already initialized)
{
  mpi_level_count() += 1;
  get_comm_list().push_back(
      Q_Comm(comm, size_node, RngState("sync_node:" + show(size_node))));
  get_comm_internal() = get_comm_list().back().comm;
  set_global_geon(get_comm_list().back().size_node);
  sync_node();
  if (mpi_level_count() == 1) {
    initialize_qlat_comm();
  }
  get_id_node_list_for_shuffle() = mk_id_node_list_for_shuffle();
  get_id_node_in_shuffle_list() = mk_id_node_in_shuffle_list();
  get_id_node_in_shuffle_internal() =
      get_id_node_in_shuffle(get_id_node(), 0, get_num_node());
  // display_geometry_node();
  // install_qhandle_sig();
  clear_all_caches();
  sync_node();
}

void begin(const int id_node, const Coordinate& size_node, const int color)
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

void begin(int* argc, char** argv[], const Coordinate& size_node)
// not recommended
{
  const int num_node = init_mpi(argc, argv);
  qassert(num_node == product(size_node));
  begin_comm(MPI_COMM_WORLD, size_node);
}

void begin(int* argc, char** argv[],
           const std::vector<Coordinate>& size_node_list)
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

void end(const bool is_preserving_cache)
{
  if (get_comm_list().empty()) {
    displayln_info(ssprintf("qlat::end(): get_comm_list().empty() = true."));
    if (not is_preserving_cache) {
      clear_all_caches();
    }
  } else {
    mpi_level_count() -= 1;
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
