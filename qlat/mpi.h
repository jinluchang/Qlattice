#pragma once

#include <qlat/config.h>
#include <qlat/utils.h>

#include <mpi.h>
#include <timer.h>

#include <array>

QLAT_START_NAMESPACE

inline MPI_Comm& getQlatComm()
{
  static MPI_Comm comm;
  return comm;
}

inline MPI_Comm*& getPtrComm()
{
  static MPI_Comm* p_comm = &getQlatComm();
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
  inline void init();
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

inline int idNodeFromCoorNode(const Coordinate& coor)
{
  int rank;
  MPI_Cart_rank(getComm(), (int*)coor.data(), &rank);
  return rank;
}

inline void GeometryNode::init()
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
  assert(sizeNode[0] * sizeNode[1] * sizeNode[2] * sizeNode[3] == numNode);
  Display(cname, "GeometryNode::init()", "idNode = %5d ; coorNode = %s\n", idNode, show(coorNode).c_str());
  assert(idNodeFromCoorNode(coorNode) == idNode);
#else
  numNode = 1;
  idNode = 0;
  for (int i = 0; i < DIM; ++i) {
    sizeNode[i] = 1;
    coorNode[i] = 0;
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
      dest[0][mu] = idNodeFromCoorNode(coor);
      coor = coorNode;
      --coor[mu];
      regularizeCoordinate(coor, sizeNode);
      dest[1][mu] = idNodeFromCoorNode(coor);
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

inline int glbSum(Vector<double> recv, const Vector<double>& send)
{
  assert(recv.size() == send.size());
#ifdef USE_MULTI_NODE
  return MPI_Allreduce((double*)send.data(), recv.data(), recv.size(), MPI_DOUBLE, MPI_SUM, getComm());
#else
  memmove(recv.data(), send.data(), recv.size()* sizeof(double));
  return 0;
#endif
}

inline int glbSum(Vector<long> recv, const Vector<long>& send)
{
  assert(recv.size() == send.size());
#ifdef USE_MULTI_NODE
  return MPI_Allreduce((long*)send.data(), recv.data(), recv.size(), MPI_LONG, MPI_SUM, getComm());
#else
  memmove(recv.data(), send.data(), recv.size()* sizeof(long));
  return 0;
#endif
}

inline int glbSum(Vector<double> vec)
{
  std::vector<double> tmp(vec.size());
  assign(tmp, vec);
  return glbSum(vec, tmp);
}

inline int glbSum(Vector<long> vec)
{
  std::vector<long> tmp(vec.size());
  assign(tmp, vec);
  return glbSum(vec, tmp);
}

inline int glbSum(double& x)
{
  glbSum(Vector<double>(x));
}

inline int glbSum(long& x)
{
  glbSum(Vector<long>(x));
}

template <class M>
inline int glbSumDouble(M& x)
{
  glbSum(Vector<double>(&x, sizeof(M)/sizeof(double)));
}

template <class M>
inline int glbSumLong(M& x)
{
  glbSum(Vector<long>(&x, sizeof(M)/sizeof(long)));
}

template <class M>
inline int glbSumDouble(Vector<M>& x)
{
  glbSum(Vector<double>(x.data(), x.size()*sizeof(M)/sizeof(double)));
}

template <class M>
inline int glbSumLong(Vector<M>& x)
{
  glbSum(Vector<long>(x.data(), x.size()*sizeof(M)/sizeof(long)));
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
  glbSum(Vector<long>(&v,1));
}

inline Coordinate planSizeNode(const int numNode)
{
  if (numNode == 1) {
    return Coordinate(1, 1, 1, 1);
  } else if (numNode == 16) {
    return Coordinate(1, 2, 2, 4);
  } else {
    return Coordinate(0, 0, 0, 0);
  }
}

inline int beginMpi(int* argc, char** argv[])
{
  MPI_Init(argc, argv);
  int numNode;
  MPI_Comm_size(MPI_COMM_WORLD, &numNode);
  DisplayInfo(cname, "begin", "MPI Initialized. NumNode = %d\n", numNode);
  return numNode;
}

inline void beginCart(const MPI_Comm& mpiComm, const Coordinate& sizeNode)
{
  const Coordinate periods(1, 1, 1, 1);
  MPI_Cart_create(mpiComm, DIM, (int*)sizeNode.data(), (int*)periods.data(), 1, &getComm());
  const GeometryNode& geon = getGeometryNode();
  syncNode();
  DisplayInfo(cname, "begin", "MPI Cart created. GeometryNode =\n%s\n", show(geon).c_str());
  syncNode();
}

inline void begin(int* argc, char** argv[], const Coordinate& sizeNode)
{
  beginMpi(argc, argv);
  beginCart(MPI_COMM_WORLD, sizeNode);
}

inline void begin(int* argc, char** argv[])
{
  int numNode = beginMpi(argc, argv);
  beginCart(MPI_COMM_WORLD, planSizeNode(numNode));
}

inline void end()
{
  MPI_Finalize();
  DisplayInfo(cname, "end", "MPI Finalized.\n");
}

QLAT_END_NAMESPACE
