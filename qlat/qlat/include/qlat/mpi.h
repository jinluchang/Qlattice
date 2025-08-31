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

template <class M>
struct IsMpiDataType {
  static constexpr bool value = false;
  static constexpr bool is_complex = false;
  static const std::string get_type_name() { return "unknown_type"; }
  static MPI_Datatype get_mpi_datatype() { return MPI_BYTE; }
  using ElementaryType = M;
};

template <>
struct IsMpiDataType<int8_t> {
  static constexpr bool value = true;
  static constexpr bool is_complex = false;
  static const std::string get_type_name() { return "Int8t"; }
  static MPI_Datatype get_mpi_datatype() { return MPI_INT8_T; }
  using ElementaryType = int8_t;
};

template <>
struct IsMpiDataType<int32_t> {
  static constexpr bool value = true;
  static constexpr bool is_complex = false;
  static const std::string get_type_name() { return "Int32t"; }
  static MPI_Datatype get_mpi_datatype() { return MPI_INT32_T; }
  using ElementaryType = int32_t;
};

template <>
struct IsMpiDataType<int64_t> {
  static constexpr bool value = true;
  static constexpr bool is_complex = false;
  static const std::string get_type_name() { return "Int64t"; }
  static MPI_Datatype get_mpi_datatype() { return MPI_INT64_T; }
  using ElementaryType = int64_t;
};

template <>
struct IsMpiDataType<RealF> {
  static constexpr bool value = true;
  static constexpr bool is_complex = false;
  static const std::string get_type_name() { return "RealF"; }
  static MPI_Datatype get_mpi_datatype() { return MPI_FLOAT; }
  using ElementaryType = RealF;
};

template <>
struct IsMpiDataType<RealD> {
  static constexpr bool value = true;
  static constexpr bool is_complex = false;
  static const std::string get_type_name() { return "RealD"; }
  static MPI_Datatype get_mpi_datatype() { return MPI_DOUBLE; }
  using ElementaryType = RealD;
};

template <class M>
qacc constexpr bool is_mpi_datatype()
// Char, Int, Long, RealF, RealD
{
  return IsMpiDataType<M>::value;
}

template <class M, QLAT_ENABLE_IF(is_mpi_datatype<M>())>
MPI_Datatype get_mpi_datatype()
{
  return IsMpiDataType<M>::get_mpi_datatype();
}

// ------------------------------

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

inline Int id_node_from_coor_node(const Coordinate& coor_node)
{
  return index_from_coordinate(coor_node, get_geometry_node().size_node);
}

inline Coordinate coor_node_from_id_node(Int id_node)
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
  Int dest[2][DIMN];
  // dest[dir][mu]
  // dir = 0, 1 for Plus dir or Minus dir
  // 0 <= mu < 4 for different directions
  //
  void init()
  {
    const Coordinate& coor_node = get_geometry_node().coor_node;
    const Coordinate& size_node = get_geometry_node().size_node;
    for (Int mu = 0; mu < DIMN; ++mu) {
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
    Qassert(mpi_ret == 0);
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

Int mpi_send(const void* buf, Long count, MPI_Datatype datatype, Int dest,
             Int tag, MPI_Comm comm);

Int mpi_recv(void* buf, Long count, MPI_Datatype datatype, Int source, Int tag,
             MPI_Comm comm, MPI_Status* status = MPI_STATUS_IGNORE);

Int mpi_isend(const void* buf, Long count, MPI_Datatype datatype, Int dest,
              Int tag, MPI_Comm comm, std::vector<MPI_Request>& requests);

Int mpi_irecv(void* buf, Long count, MPI_Datatype datatype, Int source, Int tag,
              MPI_Comm comm, std::vector<MPI_Request>& requests);

Int mpi_waitall(std::vector<MPI_Request>& requests);

Int mpi_alltoallv(const void* sendbuf, const Long* sendcounts,
                  const Long* sdispls, MPI_Datatype sendtype, void* recvbuf,
                  const Long* recvcounts, const Long* rdispls,
                  MPI_Datatype recvtype, MPI_Comm comm);

Int mpi_bcast(void* buffer, const Long count, MPI_Datatype datatype,
              const Int root, MPI_Comm comm);

Int mpi_allreduce(const void* sendbuf, void* recvbuf, const Long count,
                  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);

template <class M, QLAT_ENABLE_IF(is_mpi_datatype<M>())>
Int mpi_allreduce(const Vector<M>& sendbuf, Vector<M> recvbuf, MPI_Op op,
                  MPI_Comm comm)
{
  Qassert(sendbuf.size() == recvbuf.size());
  MPI_Datatype datatype = get_mpi_datatype<M>();
  return mpi_allreduce(sendbuf.data(), recvbuf.data(), sendbuf.size(), datatype,
                       op, comm);
}

template <class M, QLAT_ENABLE_IF(is_mpi_datatype<M>())>
Int glb_sum(Vector<M> recvbuf, const Vector<M>& sendbuf)
{
  if (is_same<M, Char>()) {
    return mpi_allreduce(sendbuf, recvbuf, MPI_BXOR, get_comm());
  }
  return mpi_allreduce(sendbuf, recvbuf, MPI_SUM, get_comm());
}

template <class M, QLAT_ENABLE_IF(is_mpi_datatype<M>())>
Int glb_max(Vector<M> recvbuf, const Vector<M>& sendbuf)
{
  return mpi_allreduce(sendbuf, recvbuf, MPI_MAX, get_comm());
}

template <class M, QLAT_ENABLE_IF(is_mpi_datatype<M>())>
Int glb_min(Vector<M> recvbuf, const Vector<M>& sendbuf)
{
  return mpi_allreduce(sendbuf, recvbuf, MPI_MIN, get_comm());
}

template <class T, QLAT_ENABLE_IF(is_get_data_type<T>())>
Int glb_sum(T& xx)
{
  TIMER("glb_sum(T&xx)");
  using E = typename IsDataVectorType<T>::ElementaryType;
  Vector<E> vec = get_data_in_elementary_type(xx);
  vector<E> tmp_vec(vec.size(), MemType::Comm);
  vector<E> tmp2_vec(vec.size(), MemType::Comm);
  assign(tmp_vec, vec);
  const Int ret = glb_sum(get_data(tmp2_vec), get_data(tmp_vec));
  assign(vec, tmp2_vec);
  return ret;
}

template <class M, QLAT_ENABLE_IF(is_data_value_type<M>())>
Int glb_sum(Vector<M> xx)
// so that xx don't have to be a reference
{
  return glb_sum<Vector<M>>(xx);
}

template <class M, QLAT_ENABLE_IF(is_data_value_type<M>())>
M f_glb_sum(const M& xx)
// so that xx don't have to be a reference
{
  M v = xx;
  const Int code = glb_sum(v);
  Qassert(code == 0);
  return v;
}

template <class T, QLAT_ENABLE_IF(is_get_data_type<T>())>
Int glb_max(T& xx)
// Both send and recv buffer will be copied to and from memory of type
// MemType::Comm.
{
  TIMER("glb_max(T&xx)");
  using E = typename IsDataVectorType<T>::ElementaryType;
  Vector<E> vec = get_data_in_elementary_type(xx);
  vector<E> tmp_vec(vec.size(), MemType::Comm);
  vector<E> tmp2_vec(vec.size(), MemType::Comm);
  assign(tmp_vec, vec);
  const Int ret = glb_max(get_data(tmp2_vec), get_data(tmp_vec));
  assign(vec, tmp2_vec);
  return ret;
}

template <class M, QLAT_ENABLE_IF(is_data_value_type<M>())>
Int glb_max(Vector<M> xx)
// so that xx don't have to be a reference
{
  return glb_max<Vector<M>>(xx);
}

template <class M, QLAT_ENABLE_IF(is_data_value_type<M>())>
M f_glb_max(const M& xx)
// so that xx don't have to be a reference
{
  M v = xx;
  const Int code = glb_max(v);
  Qassert(code == 0);
  return v;
}

template <class T, QLAT_ENABLE_IF(is_get_data_type<T>())>
Int glb_min(T& xx)
{
  TIMER("glb_min(T&xx)");
  using E = typename IsDataVectorType<T>::ElementaryType;
  Vector<E> vec = get_data_in_elementary_type(xx);
  vector<E> tmp_vec(vec.size(), MemType::Comm);
  vector<E> tmp2_vec(vec.size(), MemType::Comm);
  assign(tmp_vec, vec);
  const Int ret = glb_min(get_data(tmp2_vec), get_data(tmp_vec));
  assign(vec, tmp2_vec);
  return ret;
}

template <class M, QLAT_ENABLE_IF(is_data_value_type<M>())>
Int glb_min(Vector<M> xx)
// so that xx don't have to be a reference
{
  return glb_min<Vector<M>>(xx);
}

template <class M, QLAT_ENABLE_IF(is_data_value_type<M>())>
M f_glb_min(const M& xx)
// so that xx don't have to be a reference
{
  M v = xx;
  const Int code = glb_min(v);
  Qassert(code == 0);
  return v;
}

bool glb_any(const bool b);

bool glb_all(const bool b);

Int bcast(Vector<Char> recv, const Int root = 0);

template <class T, QLAT_ENABLE_IF(is_get_data_type<T>())>
Int bcast(T& xx, const Int root = 0)
{
  Vector<Char> vec = get_data_char(xx);
  return bcast(vec, root);
}

template <class M, QLAT_ENABLE_IF(is_data_value_type<M>())>
Int bcast(Vector<M> xx, const Int root = 0)
// so that xx don't have to be a reference
{
  return bcast<Vector<M>>(xx, root);
}

template <class M, QLAT_ENABLE_IF(is_data_value_type<M>())>
Int bcast(std::vector<M>& recv, const Int root = 0)
{
  if (1 == get_num_node()) {
    return 0;
  }
  Int ret = 0;
  Long size = recv.size();
  ret += bcast<Long>(size, root);
  recv.resize(size);
  ret += bcast(get_data(recv), root);
  return ret;
}

template <class M, QLAT_ENABLE_IF(is_data_value_type<M>())>
Int bcast(vector<M>& recv, const Int root = 0)
{
  if (1 == get_num_node()) {
    return 0;
  }
  Int ret = 0;
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
Int bcast(std::vector<std::vector<M>>& datatable, const Int root = 0)
{
  if (1 == get_num_node()) {
    return 0;
  }
  Int ret = 0;
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

Int bcast(std::string& recv, const Int root = 0);

Int bcast(std::vector<std::string>& recv, const Int root = 0);

Int bcast(PointsSelection& psel, const Int root = 0);

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
  const Int code = bcast_any(v, b);
  Qassert(code == 0);
  return v;
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

// -----------------------------------

std::vector<Int> mk_id_node_list_for_shuffle_rs(const RngState& rs);

std::vector<Int> mk_id_node_list_for_shuffle_step_size(const Int step_size_);

std::vector<Int> mk_id_node_list_for_shuffle_node();

std::vector<Int> mk_id_node_list_for_shuffle();

std::vector<Int> mk_id_node_in_shuffle_list();

// `id_node` is the usual ID for MPI processes.
// `id_node_in_shuffle` is a shuffled ID for MPI processes, which is used for
// purposes like IO. This is useful for example when one physical node runs
// multiple MPI processes, but usually IO happens for MPI processes with the
// first few IDs. To prevent IO only uses the all MPI processes within the first
// few nodes, we can use `id_node_in_shuffle` for parallel IO.

Int get_id_node_in_shuffle(const Int id_node, const Int new_num_node,
                           const Int num_node);

Int get_id_node_from_id_node_in_shuffle(const Int id_node_in_shuffle,
                                        const Int new_num_node,
                                        const Int num_node);

void set_node_rank_size(Int& node_rank, Int& node_size);

std::string get_hostname();

void display_geometry_node();

Coordinate plan_size_node(const Int num_node);

bool is_MPI_initialized();

Int init_mpi(Int* argc, char** argv[]);

void set_global_geon(const Coordinate& size_node);

void set_cuda_device();

void display_qlat_banner();

void initialize_qlat_comm();

API Long& mpi_level_count();

void begin_comm(const MPI_Comm comm, const Coordinate& size_node);

void begin(const Int id_node, const Coordinate& size_node, const Int color = 0);

void begin(Int* argc, char** argv[], const Coordinate& size_node);

void begin(
    Int* argc, char** argv[],
    const std::vector<Coordinate>& size_node_list = std::vector<Coordinate>());

void begin_once(const Int id_node, const Coordinate& size_node,
                const Int color = 0);

void end(const bool is_preserving_cache = false);

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
Int receive_job(int64_t& flag, N& data, const Int root = 0)
{
  const Int mpi_tag = 3;
  const Int count = sizeof(int64_t) + sizeof(N);
  std::vector<char> fdata(count, (char)0);
  Int ret = mpi_recv(fdata.data(), count, MPI_BYTE, root, mpi_tag, get_comm());
  extract_flag_data(flag, data, fdata);
  return ret;
}

template <class M>
Int send_result(const int64_t flag, const M& data, const Int root = 0)
{
  const Int mpi_tag = 2;
  std::vector<char> fdata = pad_flag_data(flag, data);
  return mpi_send(fdata.data(), fdata.size(), MPI_BYTE, root, mpi_tag,
                  get_comm());
}

template <class N>
Int send_job(const int64_t flag, const N& data, const Int dest)
{
  const Int mpi_tag = 3;
  std::vector<char> fdata = pad_flag_data(flag, data);
  return mpi_send(fdata.data(), fdata.size(), MPI_BYTE, dest, mpi_tag,
                  get_comm());
}

template <class M>
Int receive_result(Int& source, int64_t& flag, M& result)
{
  const Int mpi_tag = 2;
  const Int count = sizeof(int64_t) + sizeof(M);
  std::vector<char> fdata(count, (char)0);
  MPI_Status status;
  const Int ret = mpi_recv(fdata.data(), fdata.size(), MPI_BYTE, MPI_ANY_SOURCE,
                           mpi_tag, get_comm(), &status);
  source = status.MPI_SOURCE;
  extract_flag_data(flag, result, fdata);
  return ret;
}

template <class M>
Int get_data_dir(Vector<M> recv, const Vector<M>& send, const Int dir)
// dir = 0, 1 for Plus dir or Minus dir
{
  TIMER_FLOPS("get_data_dir");
  const Int mpi_tag = 0;
  Qassert(recv.size() == send.size());
  const Long size = recv.size() * sizeof(M);
  timer.flops += size;
  const Int self_ID = get_id_node();
  const Int idf = (self_ID + 1 - 2 * dir + get_num_node()) % get_num_node();
  const Int idt = (self_ID - 1 + 2 * dir + get_num_node()) % get_num_node();
  //
  MPI_Request req;
  MPI_Isend((void*)send.data(), size, MPI_BYTE, idt, mpi_tag, get_comm(), &req);
  const Int ret =
      mpi_recv(recv.data(), size, MPI_BYTE, idf, mpi_tag, get_comm());
  MPI_Wait(&req, MPI_STATUS_IGNORE);
  return ret;
}

template <class M>
Int get_data_dir_mu(Vector<M> recv, const Vector<M>& send, const Int dir,
                    const Int mu)
// dir = 0, 1 for Plus dir or Minus dir
// 0 <= mu < 4 for different directions
{
  TIMER_FLOPS("get_data_dir_mu");
  const Int mpi_tag = 1;
  Qassert(recv.size() == send.size());
  const Long size = recv.size() * sizeof(M);
  timer.flops += size;
  const GeometryNodeNeighbor& geonb = get_geometry_node_neighbor();
  const Int idf = geonb.dest[dir][mu];
  const Int idt = geonb.dest[1 - dir][mu];
  MPI_Request req;
  MPI_Isend((void*)send.data(), size, MPI_BYTE, idt, mpi_tag, get_comm(), &req);
  const Int ret =
      mpi_recv(recv.data(), size, MPI_BYTE, idf, mpi_tag, get_comm());
  MPI_Wait(&req, MPI_STATUS_IGNORE);
  return ret;
}

template <class M>
Int get_data_plus_mu(Vector<M> recv, const Vector<M>& send, const Int mu)
{
  return get_data_dir_mu(recv, send, 0, mu);
}

template <class M>
Int get_data_minus_mu(Vector<M> recv, const Vector<M>& send, const Int mu)
{
  return get_data_dir_mu(recv, send, 1, mu);
}

// ----------------------------------

RngState& get_sync_node_rs_mpi();

Int glb_sum_long_vec_mpi(void* ptr, const Long size);

Int glb_sum_int_vec_mpi(void* ptr, const Long size);

Int glb_sum_real_d_vec_mpi(void* ptr, const Long size);

Int glb_sum_real_f_vec_mpi(void* ptr, const Long size);

Int glb_sum_byte_vec_mpi(void* ptr, const Long size);

Int bcast_byte_vec_mpi(void* ptr, const Long size, const Int root);

}  // namespace qlat
