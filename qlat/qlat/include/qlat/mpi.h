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

API inline std::vector<Int>& get_id_node_list_for_shuffle()
// qlat parameter
// initialized in begin_comm with mk_id_node_list_for_shuffle()
// return list
// list[id_node_in_shuffle] = id_node
{
  static std::vector<Int> list;
  return list;
}

API inline std::vector<Int>& get_id_node_in_shuffle_list()
// qlat parameter
// initialized in begin_comm with mk_id_node_in_shuffle_list()
// return list
// list[id_node] = id_node_in_shuffle
{
  static std::vector<Int> list;
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

// -------------------

struct API MpiDataType {
  MPI_Datatype mpi_dtype;
  //
  MpiDataType() { mpi_dtype = MPI_DATATYPE_NULL; }
  ~MpiDataType()
  {
    if (mpi_dtype != MPI_DATATYPE_NULL) {
      MPI_Type_free(&mpi_dtype);
      mpi_dtype = MPI_DATATYPE_NULL;
    }
  }
  //
  void set_contiguous(const Int size)
  {
    const Int mpi_ret = MPI_Type_contiguous(size, MPI_BYTE, &mpi_dtype);
    qassert(mpi_ret == 0);
    MPI_Type_commit(&mpi_dtype);
  }
};

API inline Cache<std::string, MpiDataType>& get_mpi_data_type_cache()
{
  static Cache<std::string, MpiDataType> cache("MpiDataTypeCache", 16);
  return cache;
}

inline const MpiDataType& get_mpi_data_type_contiguous(const Int& size)
{
  const std::string key = ssprintf("set_contiguous: %d", size);
  if (!get_mpi_data_type_cache().has(key)) {
    MpiDataType& mpi_dtype = get_mpi_data_type_cache()[key];
    mpi_dtype.set_contiguous(size);
  }
  return get_mpi_data_type_cache()[key];
}

// ----------------------------------

int mpi_send(const void* buf, Long count, MPI_Datatype datatype, int dest,
             int tag, MPI_Comm comm);

int mpi_recv(void* buf, Long count, MPI_Datatype datatype, int source, int tag,
             MPI_Comm comm, MPI_Status* status = MPI_STATUS_IGNORE);

int mpi_isend(const void* buf, Long count, MPI_Datatype datatype, int dest,
              int tag, MPI_Comm comm, std::vector<MPI_Request>& requests);

int mpi_irecv(void* buf, Long count, MPI_Datatype datatype, int source, int tag,
              MPI_Comm comm, std::vector<MPI_Request>& requests);

int mpi_waitall(std::vector<MPI_Request>& requests);

int mpi_alltoallv(const void* sendbuf, const Long* sendcounts,
                  const Long* sdispls, MPI_Datatype sendtype, void* recvbuf,
                  const Long* recvcounts, const Long* rdispls,
                  MPI_Datatype recvtype, MPI_Comm comm);

int mpi_bcast(void* buffer, const Long count, MPI_Datatype datatype,
              const int root, MPI_Comm comm);

int mpi_allreduce(void* sendbuf, void* recvbuf, const Long count,
                  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);

int glb_sum(Vector<RealD> recv, const Vector<RealD>& send);

int glb_sum(Vector<RealF> recv, const Vector<RealF>& send);

int glb_sum(Vector<Long> recv, const Vector<Long>& send);

int glb_sum(Vector<Int> recv, const Vector<Int>& send);

int glb_sum(Vector<Char> recv, const Vector<Char>& send);

int glb_sum(Vector<char> recv, const Vector<char>& send);

template <class T, QLAT_ENABLE_IF(is_get_data_type<T>())>
int glb_sum(T& xx)
{
  TIMER("glb_sum(T&xx)");
  using E = typename IsDataVectorType<T>::ElementaryType;
  Vector<E> vec = get_data_in_elementary_type(xx);
  vector<E> tmp_vec(vec.size(), MemType::Comm);
  vector<E> tmp2_vec(vec.size(), MemType::Comm);
  assign(tmp_vec, vec);
  const int ret = glb_sum(get_data(tmp2_vec), get_data(tmp_vec));
  assign(vec, tmp2_vec);
  return ret;
}

template <class M, QLAT_ENABLE_IF(is_data_value_type<M>())>
int glb_sum(Vector<M> xx)
// so that xx don't have to be a reference
{
  return glb_sum<Vector<M>>(xx);
}

template <class M, QLAT_ENABLE_IF(is_data_value_type<M>())>
M f_glb_sum(const M& xx)
// so that xx don't have to be a reference
{
  M v = xx;
  const int code = glb_sum(v);
  qassert(code == 0);
  return v;
}

bool glb_any(const bool b);

bool glb_all(const bool b);

int bcast(Vector<Char> recv, const int root = 0);

template <class T, QLAT_ENABLE_IF(is_get_data_type<T>())>
int bcast(T& xx, const int root = 0)
{
  Vector<Char> vec = get_data_char(xx);
  return bcast(vec, root);
}

template <class M, QLAT_ENABLE_IF(is_data_value_type<M>())>
int bcast(Vector<M> xx, const int root = 0)
// so that xx don't have to be a reference
{
  return bcast<Vector<M>>(xx, root);
}

template <class M, QLAT_ENABLE_IF(is_data_value_type<M>())>
int bcast(std::vector<M>& recv, const int root = 0)
{
  if (1 == get_num_node()) {
    return 0;
  }
  int ret = 0;
  Long size = recv.size();
  ret += bcast<Long>(size, root);
  recv.resize(size);
  ret += bcast(get_data(recv), root);
  return ret;
}

template <class M, QLAT_ENABLE_IF(is_data_value_type<M>())>
int bcast(vector<M>& recv, const int root = 0)
{
  if (1 == get_num_node()) {
    return 0;
  }
  int ret = 0;
  Long size = recv.size();
  ret += bcast<Long>(size, root);
  recv.resize(size);
  vector<M> buffer(size, MemType::Comm);
  if (get_id_node() == root) {
    buffer = recv;
  }
  ret += bcast(get_data(buffer), root);
  recv = buffer;
  return ret;
}

template <class M, QLAT_ENABLE_IF(is_data_value_type<M>())>
int bcast(std::vector<std::vector<M>>& datatable, const int root = 0)
{
  if (1 == get_num_node()) {
    return 0;
  }
  int ret = 0;
  Long nrow, total_size;
  vector<Long> row_sizes;
  vector<M> data;
  if (get_id_node() == root) {
    row_sizes = vector_map_size(datatable);
    data = vector_concat(datatable);
    nrow = row_sizes.size();
    total_size = data.size();
  }
  ret += bcast(nrow, root);
  ret += bcast(total_size, root);
  if (get_id_node() != root) {
    row_sizes.resize(nrow);
    data.resize(total_size);
  }
  ret += bcast(row_sizes, root);
  ret += bcast(data, root);
  if (get_id_node() != root) {
    datatable = vector_split(data, row_sizes);
  }
  return ret;
}

int bcast(std::string& recv, const int root = 0);

int bcast(std::vector<std::string>& recv, const int root = 0);

int bcast(PointsSelection& psel, const int root = 0);

Int bcast_any(Vector<Char> xx, const bool b);

template <class T, QLAT_ENABLE_IF(is_get_data_type<T>())>
Int bcast_any(T& xx, const bool b)
// bcast to all nodes from any node if `b == true`.
// `glb_any(b)` should be `true`, otherwise will return `-1`.
{
  return bcast_any(get_data_char(xx), b);
}

template <class M, QLAT_ENABLE_IF(is_data_value_type<M>())>
M f_bcast_any(const M& xx, const bool b)
{
  M v = xx;
  const int code = bcast_any(v, b);
  qassert(code == 0);
  return v;
}

std::vector<Int> mk_id_node_list_for_shuffle_rs(const RngState& rs);

std::vector<Int> mk_id_node_list_for_shuffle_step_size(const int step_size_);

std::vector<Int> mk_id_node_list_for_shuffle_node();

std::vector<Int> mk_id_node_list_for_shuffle();

std::vector<Int> mk_id_node_in_shuffle_list();

// `id_node` is the usual ID for MPI processes.
// `id_node_in_shuffle` is a shuffled ID for MPI processes, which is used for
// purposes like IO. This is useful for example when one physical node runs
// multiple MPI processes, but usually IO happens for MPI processes with the
// first few IDs. To prevent IO only uses the all MPI processes within the first
// few nodes, we can use `id_node_in_shuffle` for parallel IO.

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

API Long& mpi_level_count();

void begin_comm(const MPI_Comm comm, const Coordinate& size_node);

void begin(const int id_node, const Coordinate& size_node, const int color = 0);

void begin(int* argc, char** argv[], const Coordinate& size_node);

void begin(
    int* argc, char** argv[],
    const std::vector<Coordinate>& size_node_list = std::vector<Coordinate>());

void end(const bool is_preserving_cache = false);

inline void begin_once(const int id_node, const Coordinate& size_node,
                       const int color = 0)
// Only actually call begin if we haven't call it before.
// Can be called any times.
{
  if (mpi_level_count() == 0) {
    begin(id_node, size_node, color);
  }
}

// ----------------------------------

template <class M>
std::vector<char> pad_flag_data(const int64_t flag, const M& data)
{
  std::vector<char> fdata(8 + sizeof(M), (char)0);
  std::memcpy(&fdata[0], &flag, sizeof(int64_t));
  std::memcpy(&fdata[8], &data, sizeof(M));
  return fdata;
}

template <class M>
void extract_flag_data(int64_t& flag, M& data, const std::vector<char>& fdata)
{
  flag = fdata[0];
  std::memcpy(&flag, &fdata[0], sizeof(int64_t));
  std::memcpy(&data, &fdata[8], sizeof(M));
}

template <class N>
int receive_job(int64_t& flag, N& data, const int root = 0)
{
  const int mpi_tag = 3;
  const int count = sizeof(int64_t) + sizeof(N);
  std::vector<char> fdata(count, (char)0);
  int ret = mpi_recv(fdata.data(), count, MPI_BYTE, root, mpi_tag, get_comm());
  extract_flag_data(flag, data, fdata);
  return ret;
}

template <class M>
int send_result(const int64_t flag, const M& data, const int root = 0)
{
  const int mpi_tag = 2;
  std::vector<char> fdata = pad_flag_data(flag, data);
  return mpi_send(fdata.data(), fdata.size(), MPI_BYTE, root, mpi_tag,
                  get_comm());
}

template <class N>
int send_job(const int64_t flag, const N& data, const int dest)
{
  const int mpi_tag = 3;
  std::vector<char> fdata = pad_flag_data(flag, data);
  return mpi_send(fdata.data(), fdata.size(), MPI_BYTE, dest, mpi_tag,
                  get_comm());
}

template <class M>
int receive_result(int& source, int64_t& flag, M& result)
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
  const Long size = recv.size() * sizeof(M);
  timer.flops += size;
  const int self_ID = get_id_node();
  const int idf = (self_ID + 1 - 2 * dir + get_num_node()) % get_num_node();
  const int idt = (self_ID - 1 + 2 * dir + get_num_node()) % get_num_node();
  //
  MPI_Request req;
  MPI_Isend((void*)send.data(), size, MPI_BYTE, idt, mpi_tag, get_comm(), &req);
  const int ret =
      mpi_recv(recv.data(), size, MPI_BYTE, idf, mpi_tag, get_comm());
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
  const Long size = recv.size() * sizeof(M);
  timer.flops += size;
  const GeometryNodeNeighbor& geonb = get_geometry_node_neighbor();
  const int idf = geonb.dest[dir][mu];
  const int idt = geonb.dest[1 - dir][mu];
  MPI_Request req;
  MPI_Isend((void*)send.data(), size, MPI_BYTE, idt, mpi_tag, get_comm(), &req);
  const int ret =
      mpi_recv(recv.data(), size, MPI_BYTE, idf, mpi_tag, get_comm());
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
  if (sizeof(Long) == sizeof(int64_t)) {
    return glb_sum_int64_vec(x);
  } else if (sizeof(Long) == sizeof(int32_t)) {
    return glb_sum_int32_vec(x);
  } else {
    qassert(false);
    return 0;
  }
}

template <class M>
int glb_sum_int64_vec(Vector<M> x)
{
  return glb_sum(
      Vector<int64_t>((int64_t*)x.data(), x.data_size() / sizeof(int64_t)));
}

template <class M>
int glb_sum_int32_vec(Vector<M> x)
{
  return glb_sum(
      Vector<int32_t>((int32_t*)x.data(), x.data_size() / sizeof(int32_t)));
}

template <class M>
int glb_sum_byte_vec(Vector<M> x)
{
  return glb_sum(Vector<char>((char*)x.data(), x.data_size()));
}

template <class M>
int glb_sum_vec(Vector<M> x)
{
  if (is_composed_of_real_d<M>()) {
    return glb_sum_double_vec(x);
  } else if (is_composed_of_long<M>()) {
    return glb_sum_int64_vec(x);
  } else if (is_composed_of_real_f<M>()) {
    return glb_sum_float_vec(x);
  } else {
    qerr(ssprintf("glb_sum_vec get_type_name(M)='%s'",
                  get_type_name<M>().c_str()));
    return 0;
  }
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
int glb_sum_byte(M& x)
{
  return glb_sum(Vector<char>((char*)&x, sizeof(M)));
}

Int all_gather(Vector<Char> recv, const Vector<Char> send);

template <class T1, class T2,
          class E1 = typename IsDataVectorType<T1>::ElementaryType,
          class E2 = typename IsDataVectorType<T2>::ElementaryType,
          QLAT_ENABLE_IF(is_data_vector_type<T1>() and
                         is_data_vector_type<T2>() and (is_same<E1, E2>()))>
Int all_gather(T1& recv, const T2& send)
{
  Vector<Char> vec_recv = get_data_char(recv);
  const Vector<Char> vec_send = get_data_char(send);
  return all_gather(vec_recv, vec_send);
}

template <class M1, class T2,
          class E1 = typename IsDataValueType<M1>::ElementaryType,
          class E2 = typename IsDataVectorType<T2>::ElementaryType,
          QLAT_ENABLE_IF(is_data_value_type<M1>() and
                         is_data_vector_type<T2>() and (is_same<E1, E2>()))>
Int all_gather(Vector<M1> recv, const T2& send)
{
  Vector<Char> vec_recv = get_data_char(recv);
  const Vector<Char> vec_send = get_data_char(send);
  return all_gather(vec_recv, vec_send);
}

// ----------------------------------

inline RngState& get_sync_node_rs_mpi()
{
  return get_comm_list().back().sync_node_rs;
}

inline Int glb_sum_long_vec_mpi(void* ptr, const Long size)
{
  const Long n = size / sizeof(Long);
  Vector<Long> data((Long*)ptr, n);
  qassert(data.data_size() == size);
  return glb_sum(data);
}

inline Int glb_sum_int_vec_mpi(void* ptr, const Long size)
{
  const Long n = size / sizeof(Int);
  Vector<Int> data((Int*)ptr, n);
  qassert(data.data_size() == size);
  return glb_sum(data);
}

inline Int glb_sum_real_d_vec_mpi(void* ptr, const Long size)
{
  const Long n = size / sizeof(RealD);
  Vector<RealD> data((RealD*)ptr, n);
  qassert(data.data_size() == size);
  return glb_sum(data);
}

inline Int glb_sum_real_f_vec_mpi(void* ptr, const Long size)
{
  const Long n = size / sizeof(RealF);
  Vector<RealF> data((RealF*)ptr, n);
  qassert(data.data_size() == size);
  return glb_sum(data);
}

inline Int glb_sum_byte_vec_mpi(void* ptr, const Long size)
{
  Vector<Char> data((Char*)ptr, size);
  return glb_sum(data);
}

inline Int bcast_byte_vec_mpi(void* ptr, const Long size, const Int root)
{
  Vector<Char> data((Char*)ptr, size);
  return bcast(data, root);
}

}  // namespace qlat
