#pragma once

#include <lqps/config.h>
#include <lqps/utils.h>

#include <mpi.h>
#include <timer.h>

#include <array>

LQPS_START_NAMESPACE

inline MPI_Comm& getLqpsComm()
{
  static MPI_Comm comm;
  return comm;
}

inline MPI_Comm*& getPtrComm()
{
  static MPI_Comm* p_comm = NULL;
  return p_comm;
}

inline MPI_Comm& getComm()
{
  return *getPtrComm();
}

struct GeometryNode
{
  bool initialized;
  // About node geometry.
  int numNode;
  // numNode = sizeNode[0] * sizeNode[1] * sizeNode[2] * sizeNode[3]
  int idNode;
  // idNode = getIdNode()
  // 0 <= idNode < numNode
  Coordinate sizeNode;
  Coordinate coorNode;
  // 0 <= coorNode[i] < sizeNode[i]
  //
  void init();
  //
  GeometryNode()
  {
    initialized = false;
    memset(this, 0, sizeof(GeometryNode));
  }
  GeometryNode(const bool initialize)
  {
    memset(this, 0, sizeof(GeometryNode));
    if (initialize) {
      init();
    }
  }
  //
  const GeometryNode& operator=(const GeometryNode& geon)
  {
    memcpy(this, &geon, sizeof(GeometryNode));
    return *this;
  }
};

void GeometryNode::init()
{
  if (initialized) {
    return;
  }
#ifdef USE_MULTI_NODE
  MPI_Comm_size(getComm(), &numNode);
  MPI_Comm_rank(getComm(), &idNode);
  int ndims;
  MPI_Cartdim_get(getComm(), &ndims);
  assert(DIM == ndims);
  Coordinate periods;
  MPI_Cart_get(getComm(), DIM, sizeNode.data(), periods.data(), coorNode.data());
  for (int i = 0; i < DIM; ++i) {
    assert(0 != periods[i]);
  }
  assert(indexFromCoordinate(coorNode, sizeNode) == idNode);
  assert(sizeNode[0] * sizeNode[1] * sizeNode[2] * sizeNode[3] == numNode);
#else
  numNode = 1;
  idNode = 0;
  for (int i = 0; i < DIM; ++i) {
    sizeNode[i] = 1;
    coorNode[i] = 1;
  }
#endif
  initialized = true;
}

inline bool isInitialized(const GeometryNode& geon)
{
  return geon.initialized;
}

inline void init(GeometryNode& geon)
{
  geon.init();
}

inline const GeometryNode& getGeometryNode()
{
  static GeometryNode geon(true);
  return geon;
}

std::string show(const GeometryNode& geon) {
  std::string s;
  ssprintf(s, "{ initialized = %s\n", show(geon.initialized).c_str());
  ssprintf(s, ", numNode     = %d\n", geon.numNode);
  ssprintf(s, ", idNode      = %d\n", geon.idNode);
  ssprintf(s, ", sizeNode    = %s\n", show(geon.sizeNode).c_str());
  ssprintf(s, ", coorNode    = %s }", show(geon.coorNode).c_str());
  return s;
}

inline bool operator==(const GeometryNode& geon1, const GeometryNode& geon2)
{
  return 0 == memcmp(&geon1, &geon2, sizeof(GeometryNode));
}

inline int getNumNode()
{
  return getGeometryNode().numNode;
}

inline int getIdNode()
{
  return getGeometryNode().idNode;
}

inline const Coordinate& getSizeNode()
{
  return getGeometryNode().sizeNode;
}

inline const Coordinate& getCoorNode()
{
  return getGeometryNode().coorNode;
}

struct GeometryNodeNeighbor
{
  int dest[2][DIM];
  // dest[dir][mu]
  // dir = 0, 1 for Plus dir or Minus dir
  // 0 <= mu < 4 for different directions
  //
  void init()
  {
    const Coordinate& coorNode = getGeometryNode().coorNode;
    const Coordinate& sizeNode = getGeometryNode().sizeNode;
    for (int mu = 0; mu < DIM; ++mu) {
      Coordinate coor;
      coor = coorNode;
      ++coor[mu];
      regularizeCoordinate(coor, sizeNode);
      dest[0][mu] = indexFromCoordinate(coor, sizeNode);
      coor = coorNode;
      --coor[mu];
      regularizeCoordinate(coor, sizeNode);
      dest[1][mu] = indexFromCoordinate(coor, sizeNode);
    }
  }
  //
  GeometryNodeNeighbor()
  {
  }
  GeometryNodeNeighbor(bool initialize)
  {
    if (initialize) {
      init();
    }
  }
};

inline const GeometryNodeNeighbor& getGeometryNodeNeighbor()
{
  static GeometryNodeNeighbor geonb(true);
  return geonb;
}

template <class M>
int getDataDirMu(Vector<M> recv, const Vector<M>& send, const int dir, const int mu)
  // dir = 0, 1 for Plus dir or Minus dir
  // 0 <= mu < 4 for different directions
{
  TIMER_FLOPS("getDataDirMu");
  assert(recv.size() == send.size());
  const long size = recv.size()*sizeof(M);
  timer.flops += size;
#ifdef USE_MULTI_NODE
  const GeometryNodeNeighbor& geonb = getGeometryNodeNeighbor();
  const int idf = geonb.dest[dir][mu];
  const int idt = geonb.dest[1-dir][mu];
  MPI_Request req;
  MPI_Isend((void*)send.data(), size, MPI_BYTE, idt, 0, getComm(), &req);
  const int ret = MPI_Recv(recv.data(), size, MPI_BYTE, idf, 0, getComm(), MPI_STATUS_IGNORE);
  MPI_Wait(&req, MPI_STATUS_IGNORE);
  return ret;
#else
  memcpy(recv.data(), send.data(), size);
  return 0;
#endif
}

template <class M>
int getDataPlusMu(Vector<M> recv, const Vector<M>& send, const int mu)
{
  return getDataDirMu(recv, send, 0, mu);
}

template <class M>
int getDataMinusMu(Vector<M> recv, const Vector<M>& send, const int mu)
{
  return getDataDirMu(recv, send, 1, mu);
}

inline int sumVector(Vector<long> recv, const Vector<long>& send)
{
  assert(recv.size() == send.size());
#ifdef USE_MULTI_NODE
  return MPI_Allreduce((long*)send.data(), recv.data(), recv.size(), MPI_LONG, MPI_SUM, getComm());
#else
  memmove(recv.data(), send.data(), recv.size()* sizeof(long));
  return 0;
#endif
}

inline int sumVector(Vector<double> recv, const Vector<double>& send)
{
  assert(recv.size() == send.size());
#ifdef USE_MULTI_NODE
  return MPI_Allreduce((double*)send.data(), recv.data(), recv.size(), MPI_DOUBLE, MPI_SUM, getComm());
#else
  memmove(recv.data(), send.data(), recv.size()* sizeof(double));
  return 0;
#endif
}

template <class M>
int sumVector(Vector<M> vec)
{
  std::vector<M> tmp(vec.size());
  assign(tmp, vec);
  return sumVector(vec, tmp);
}

template <class M>
void allGather(Vector<M> recv, const Vector<M>& send)
{
  assert(recv.size() == send.size() * getNumNode());
  const long sendsize = send.size() * sizeof(M);
#ifdef USE_MULTI_NODE
  MPI_Allgather((void*)send, sendsize, MPI_BYTE, recv, sendsize, MPI_BYTE, getComm());
#else
  memmove(recv, send, sendsize);
#endif
  allGather(recv.data(), send.data(), send.size() * sizeof(M));
}

inline void syncNode()
{
  long v;
  sumVector(Vector<long>(&v,1));
}

inline void begin(int* argc, char** argv[], const Coordinate dims)
  // can also initialize by set getPtrComm() to point to some externally constructed COMM object
  // eg. getPtrComm() = &QMP_COMM_WORLD;
{
  MPI_Init(argc, argv);
  int numNode;
  MPI_Comm_size(MPI_COMM_WORLD, &numNode);
  DisplayInfo(cname, "begin", "MPI Initialized. NumNode=%d\n", numNode);
  getPtrComm() = &getLqpsComm();
  const Coordinate periods({ 1, 1, 1, 1 });
  MPI_Cart_create(MPI_COMM_WORLD, DIM, (int*)dims.data(), (int*)periods.data(), 1, &getLqpsComm());
  DisplayInfo(cname, "begin", "MPI Cart created. GeometryNode=\n%s\n", show(getGeometryNode()).c_str());
}

inline void end()
{
  MPI_Finalize();
  DisplayInfo(cname, "end", "MPI Finalized.\n");
}

LQPS_END_NAMESPACE
