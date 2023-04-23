#pragma once

#include <qlat-utils/cache.h>
#include <qlat-utils/lat-io.h>
#include <qlat/core.h>
#include <qlat/utils-coordinate.h>

#ifndef QLAT_NO_MALLOPT
#include <malloc.h>
#endif

#include <mpi.h>

namespace qlat
{  //

API inline std::vector<int>& get_id_node_list_for_shuffle()
// qlat parameter
// initialized in begin_comm with mk_id_node_list_for_shuffle()
// return list
// list[id_node_in_shuffle] = id_node
{
  static std::vector<int> list;
  return list;
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

// ----------------------------------

int mpi_send(const void* buf, long count, MPI_Datatype datatype, int dest,
             int tag, MPI_Comm comm);

int mpi_recv(void* buf, long count, MPI_Datatype datatype, int source, int tag,
             MPI_Comm comm, MPI_Status* status);

int mpi_isend(const void* buf, long count, MPI_Datatype datatype, int dest,
              int tag, MPI_Comm comm, std::vector<MPI_Request>& requests);

int mpi_irecv(void* buf, long count, MPI_Datatype datatype, int source, int tag,
              MPI_Comm comm, std::vector<MPI_Request>& requests);

int mpi_waitall(std::vector<MPI_Request>& requests);

int glb_sum(Vector<double> recv, const Vector<double>& send);

int glb_sum(Vector<float> recv, const Vector<float>& send);

int glb_sum(Vector<Complex> recv, const Vector<Complex>& send);

int glb_sum(Vector<ComplexF> recv, const Vector<ComplexF>& send);

int glb_sum(Vector<long> recv, const Vector<long>& send);

int glb_sum(Vector<char> recv, const Vector<char>& send);

int glb_sum(Vector<double> vec);

int glb_sum(Vector<float> vec);

int glb_sum(Vector<Complex> vec);

int glb_sum(Vector<ComplexF> vec);

int glb_sum(Vector<long> vec);

int glb_sum(Vector<char> vec);

int glb_sum(double& x);

int glb_sum(float& x);

int glb_sum(long& x);

int glb_sum(Complex& c);

int glb_sum(ComplexF& c);

int glb_sum_lat_data(LatData& ld);

int glb_sum(LatData& ld);

void bcast(int& x, const int root = 0);

void bcast(long& x, const int root = 0);

void bcast(uint32_t& x, const int root = 0);

void bcast(double& x, const int root = 0);

void bcast(Coordinate& x, const int root = 0);

void bcast(std::string& recv, const int root = 0);

void bcast(std::vector<std::string>& recv, const int root = 0);

void bcast(LatData& ld, const int root = 0);

void sync_node();

std::vector<int> mk_id_node_list_for_shuffle_rs(const RngState& rs);

std::vector<int> mk_id_node_list_for_shuffle_step_size(const int step_size_);

std::vector<int> mk_id_node_list_for_shuffle_node();

std::vector<int> mk_id_node_list_for_shuffle();

std::vector<int> mk_id_node_in_shuffle_list();

int get_id_node_in_shuffle(const int id_node, const int new_num_node,
                           const int num_node);

int get_id_node_from_id_node_in_shuffle(const int id_node_in_shuffle,
                                        const int new_num_node,
                                        const int num_node);

void set_node_rank_size(int& node_rank, int& node_size);

std::string get_hostname();

void display_geometry_node();

Coordinate plan_size_node(const int num_node);

bool is_MPI_initialized();

int init_mpi(int* argc, char** argv[]);

void set_global_geon(const Coordinate& size_node);

void set_cuda_device();

void display_qlat_banner();

void initialize_qlat_comm();

void begin_comm(const MPI_Comm comm, const Coordinate& size_node);

void begin(const int id_node, const Coordinate& size_node, const int color = 0);

void begin(int* argc, char** argv[], const Coordinate& size_node);

void begin(
    int* argc, char** argv[],
    const std::vector<Coordinate>& size_node_list = std::vector<Coordinate>());

void end(const bool is_preserving_cache = false);

// ----------------------------------

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
  const int self_ID = get_id_node();
  const int idf = (self_ID + 1 - 2 * dir + get_num_node()) % get_num_node();
  const int idt = (self_ID - 1 + 2 * dir + get_num_node()) % get_num_node();
  //
  MPI_Request req;
  MPI_Isend((void*)send.data(), size, MPI_BYTE, idt, mpi_tag, get_comm(), &req);
  const int ret = MPI_Recv(recv.data(), size, MPI_BYTE, idf, mpi_tag,
                           get_comm(), MPI_STATUS_IGNORE);
  MPI_Wait(&req, MPI_STATUS_IGNORE);
  return ret;
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
  const GeometryNodeNeighbor& geonb = get_geometry_node_neighbor();
  const int idf = geonb.dest[dir][mu];
  const int idt = geonb.dest[1 - dir][mu];
  MPI_Request req;
  MPI_Isend((void*)send.data(), size, MPI_BYTE, idt, mpi_tag, get_comm(), &req);
  const int ret = mpi_recv(recv.data(), size, MPI_BYTE, idf, mpi_tag,
                           get_comm(), MPI_STATUS_IGNORE);
  MPI_Wait(&req, MPI_STATUS_IGNORE);
  return ret;
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
int glb_sum_long_vec(Vector<M> x)
{
  return glb_sum(Vector<long>((long*)x.data(), x.data_size() / sizeof(long)));
}

template <class M>
int glb_sum_byte_vec(Vector<M> x)
{
  return glb_sum(Vector<char>((char*)x.data(), x.data_size()));
}

template <class M>
int glb_sum_double(M& x)
{
  return glb_sum(Vector<double>((double*)&x, sizeof(M) / sizeof(double)));
}

template <class M>
int glb_sum_float(M& x)
{
  return glb_sum(Vector<float>((float*)&x, sizeof(M) / sizeof(float)));
}

template <class M>
int glb_sum_long(M& x)
{
  return glb_sum(Vector<long>((long*)&x, sizeof(M) / sizeof(long)));
}

template <class M>
int glb_sum_byte(M& x)
{
  return glb_sum(Vector<char>((char*)&x, sizeof(M)));
}

template <class M>
void all_gather(Vector<M> recv, const Vector<M>& send)
{
  qassert(recv.size() == send.size() * get_num_node());
  MPI_Allgather((void*)send.data(), send.data_size(), MPI_BYTE,
                (void*)recv.data(), send.data_size(), MPI_BYTE, get_comm());
}

template <class M>
void bcast(Vector<M> recv, const int root = 0)
{
  if (1 == get_num_node()) {
    return;
  }
  MPI_Bcast((void*)recv.data(), recv.data_size(), MPI_BYTE, root, get_comm());
}

template <class M>
void bcast(std::vector<M>& recv, const int root = 0)
{
  if (1 == get_num_node()) {
    return;
  }
  long size = recv.size();
  bcast(get_data(size), root);
  recv.resize(size);
  bcast(get_data(recv), root);
}

template <class M>
void bcast(std::vector<std::vector<M> >& datatable, const int root = 0)
{
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
}

}  // namespace qlat
