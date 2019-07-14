#pragma once

#include <qlat/config.h>
#include <qlat/utils-coordinate.h>
#include <qlat/utils.h>

#ifdef USE_MULTI_NODE
#include <mpi.h>
#endif

QLAT_START_NAMESPACE

inline MPI_Comm& get_comm()
{
  static MPI_Comm comm;
  return comm;
}

struct GeometryNode {
  bool initialized;
  // About node geometry.
  int num_node;
  // num_node = size_node[0] * size_node[1] * size_node[2] * size_node[3]
  int id_node;
  // id_node = get_id_node()
  // 0 <= id_node < num_node
  Coordinate size_node;
  Coordinate coor_node;
  // 0 <= coor_node[i] < size_node[i]
  //
  inline void init() { memset(this, 0, sizeof(GeometryNode)); }
  inline void init(const int id_node_, const Coordinate& size_node_)
  {
    initialized = true;
    num_node = product(size_node_);
    id_node = id_node_;
    size_node = size_node_;
    coor_node = coordinate_from_index(id_node_, size_node_);
  }
  //
  GeometryNode() { init(); }
  GeometryNode(const int id_node_, const Coordinate& size_node_)
  {
    init(id_node_, size_node_);
  }
};

inline bool is_initialized(const GeometryNode& geon)
{
  return geon.initialized;
}

inline void init(GeometryNode& geon) { geon.init(); }

inline GeometryNode& get_geometry_node_internal()
{
  static GeometryNode geon;
  return geon;
}

inline const GeometryNode& get_geometry_node()
{
  return get_geometry_node_internal();
}

inline bool operator==(const GeometryNode& geon1, const GeometryNode& geon2)
{
  return geon1.initialized == geon2.initialized &&
         geon1.num_node == geon2.num_node && geon1.id_node == geon2.id_node &&
         geon1.size_node == geon2.size_node &&
         geon1.coor_node == geon2.coor_node;
}

inline bool operator!=(const GeometryNode& geon1, const GeometryNode& geon2)
{
  return !(geon1 == geon2);
}

inline int id_node_from_coor_node(const Coordinate& coor_node)
{
  return index_from_coordinate(coor_node, get_geometry_node().size_node);
}

inline Coordinate coor_node_from_id_node(int id_node)
{
  return coordinate_from_index(id_node, get_geometry_node().size_node);
}

inline int get_num_node() { return get_geometry_node().num_node; }

inline int get_id_node() { return get_geometry_node().id_node; }

inline const Coordinate& get_size_node()
{
  return get_geometry_node().size_node;
}

inline const Coordinate& get_coor_node()
{
  return get_geometry_node().coor_node;
}

struct GeometryNodeNeighbor {
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

inline const GeometryNodeNeighbor& get_geometry_node_neighbor()
{
  static GeometryNodeNeighbor geonb(true);
  return geonb;
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
  int ret = MPI_Recv(fdata.data(), count, MPI_BYTE, root, mpi_tag, get_comm(),
                     MPI_STATUS_IGNORE);
  extract_flag_data(flag, data, fdata);
  return ret;
}

template <class M>
inline int send_result(const int64_t flag, const M& data, const int root = 0)
{
  const int mpi_tag = 2;
  std::vector<char> fdata = pad_flag_data(flag, data);
  return MPI_Send(fdata.data(), fdata.size(), MPI_BYTE, root, mpi_tag,
                  get_comm());
}

template <class N>
inline int send_job(const int64_t flag, const N& data, const int dest)
{
  const int mpi_tag = 3;
  std::vector<char> fdata = pad_flag_data(flag, data);
  return MPI_Send(fdata.data(), fdata.size(), MPI_BYTE, dest, mpi_tag,
                  get_comm());
}

template <class M>
inline int receive_result(int& source, int64_t& flag, M& result)
{
  const int mpi_tag = 2;
  const int count = sizeof(int64_t) + sizeof(M);
  std::vector<char> fdata(count, (char)0);
  MPI_Status status;
  const int ret = MPI_Recv(fdata.data(), fdata.size(), MPI_BYTE, MPI_ANY_SOURCE,
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
  std::vector<double> tmp(vec.size());
  assign(tmp, vec);
  return glb_sum(vec, tmp);
}

inline int glb_sum(Vector<float> vec)
{
  std::vector<float> tmp(vec.size());
  assign(tmp, vec);
  return glb_sum(vec, tmp);
}

inline int glb_sum(Vector<Complex> vec)
{
  std::vector<Complex> tmp(vec.size());
  assign(tmp, vec);
  return glb_sum(vec, tmp);
}

inline int glb_sum(Vector<ComplexF> vec)
{
  std::vector<ComplexF> tmp(vec.size());
  assign(tmp, vec);
  return glb_sum(vec, tmp);
}

inline int glb_sum(Vector<long> vec)
{
  std::vector<long> tmp(vec.size());
  assign(tmp, vec);
  return glb_sum(vec, tmp);
}

inline int glb_sum(Vector<char> vec)
{
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

template <class M>
void all_gather(Vector<M> recv, const Vector<M>& send)
{
  qassert(recv.size() == send.size() * get_num_node());
  const long sendsize = send.size() * sizeof(M);
#ifdef USE_MULTI_NODE
  MPI_Allgather((void*)send.data(), send.data_size(), MPI_BYTE,
                (void*)recv.data(), send.data_size(), MPI_BYTE, get_comm());
#else
  memmove(recv, send, sendsize);
#endif
}

template <class M>
inline void bcast(Vector<M> recv, const int root = 0)
{
#ifdef USE_MULTI_NODE
  MPI_Bcast((void*)recv.data(), recv.data_size(), MPI_BYTE, root, get_comm());
#endif
}

template <class M>
inline void concat_vector(std::vector<long>& idx, std::vector<M>& data,
                          const std::vector<std::vector<M> >& datatable)
{
  idx.resize(datatable.size());
  size_t total_size = 0;
  for (size_t i = 0; i < datatable.size(); ++i) {
    const std::vector<M>& row = datatable[i];
    idx[i] = row.size();
    total_size += row.size();
  }
  data.resize(total_size);
  size_t count = 0;
  for (size_t i = 0; i < datatable.size(); ++i) {
    const std::vector<M>& row = datatable[i];
    for (size_t j = 0; j < row.size(); ++j) {
      data[count] = row[j];
      count += 1;
    }
  }
}

template <class M>
inline void split_vector(std::vector<std::vector<M> >& datatable,
                         const std::vector<long>& idx,
                         const std::vector<M>& data)
{
  clear(datatable);
  datatable.resize(idx.size());
  size_t count = 0;
  for (size_t i = 0; i < datatable.size(); ++i) {
    std::vector<M>& row = datatable[i];
    row.resize(idx[i]);
    for (size_t j = 0; j < row.size(); ++j) {
      row[j] = data[count];
      count += 1;
    }
  }
}

template <class M>
inline void bcast(std::vector<std::vector<M> >& datatable, const int root = 0)
{
#ifdef USE_MULTI_NODE
  long nrow, total_size;
  std::vector<long> idx;
  std::vector<M> data;
  if (get_id_node() == root) {
    concat_vector(idx, data, datatable);
    nrow = idx.size();
    total_size = data.size();
  }
  bcast(get_data(nrow), root);
  bcast(get_data(total_size), root);
  if (get_id_node() != root) {
    idx.resize(nrow);
    data.resize(total_size);
  }
  bcast(get_data(idx), root);
  bcast(get_data(data), root);
  if (get_id_node() != root) {
    split_vector(datatable, idx, data);
  }
#endif
}

inline void sync_node()
{
  long v = 1;
  glb_sum(Vector<long>(&v, 1));
}

inline void display_geometry_node()
{
  TIMER("display_geo_node");
  const GeometryNode& geon = get_geometry_node();
  for (int i = 0; i < geon.num_node; ++i) {
    if (i == geon.id_node) {
      displayln(std::string(fname) + " : " +
                ssprintf("id_node = %5d ; coor_node = %s", geon.id_node,
                         show(geon.coor_node).c_str()));
      fflush(get_output_file());
    }
    sync_node();
  }
  fflush(get_output_file());
  sync_node();
}

inline std::string show(const qlat::GeometryNode& geon)
{
  std::string s;
  s += ssprintf("{ initialized = %s\n", show(geon.initialized).c_str());
  s += ssprintf(", num_node    = %d\n", geon.num_node);
  s += ssprintf(", id_node     = %d\n", geon.id_node);
  s += ssprintf(", size_node   = %s\n", show(geon.size_node).c_str());
  s += ssprintf(", coor_node   = %s }", show(geon.coor_node).c_str());
  return s;
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
  displayln_info("qlat::begin(): " +
                 ssprintf("MPI Initialized. num_node = %d", num_node));
  return num_node;
}

inline void begin_comm(const MPI_Comm comm, const Coordinate& size_node)
// begin Qlat with existing comm (assuming MPI already initialized)
{
  get_comm() = comm;
  int id_node;
  MPI_Comm_rank(get_comm(), &id_node);
  GeometryNode& geon = get_geometry_node_internal();
  geon.init(id_node, size_node);
  sync_node();
  displayln_info("qlat::begin(): OMP_NUM_THREADS = " +
                 show(omp_get_max_threads()));
  displayln_info("qlat::begin(): GeometryNode =\n" + show(geon));
  fflush(get_output_file());
  sync_node();
  display_geometry_node();
}

inline void begin(const int id_node, const Coordinate& size_node)
// begin Qlat with existing id_node maping (assuming MPI already initialized)
{
  MPI_Comm comm;
  MPI_Comm_split(MPI_COMM_WORLD, 0, id_node, &comm);
  begin_comm(comm, size_node);
}

inline void begin(int* argc, char** argv[], const Coordinate& size_node)
// begin Qlat and initialize a new comm
{
  const int num_node = init_mpi(argc, argv);
  qassert(num_node == product(size_node));
  begin_comm(MPI_COMM_WORLD, size_node);
}

inline void begin(int* argc, char** argv[], const std::vector<Coordinate>& size_node_list)
{
  const int num_node = init_mpi(argc, argv);
  for (int i = 0; i < (int)size_node_list.size(); ++i) {
    const Coordinate& size_node = size_node_list[i];
    if (num_node == product(size_node)) {
      begin_comm(MPI_COMM_WORLD, size_node);
      return;
    }
  }
  qassert(false);
}

inline void begin(int* argc, char** argv[])
// begin Qlat and initialize a new comm with default topology
{
  const int num_node = init_mpi(argc, argv);
  begin_comm(MPI_COMM_WORLD, plan_size_node(num_node));
}

inline void end()
{
  if (is_MPI_initialized()) MPI_Finalize();
  displayln_info("qlat::end(): MPI Finalized.");
}

QLAT_END_NAMESPACE
