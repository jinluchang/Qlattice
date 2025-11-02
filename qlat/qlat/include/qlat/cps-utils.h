#pragma once

CPS_START_NAMESPACE

const std::complex<RealD> ii(0.0, 1.0);

template <class M>
bool isInitialized(const GridComm<M>& gc)
{
  return isInitialized(gc.getGeometry());
}

bool isInitialized(const Geometry& geo) { return geo.initialized; }

template <class T>
inline void setZeroN(T xs[], const Int xn)
{
  for (Int i = 0; i < xn; i++) {
    setZero(xs[i]);
  }
}

template <class M, Int N>
struct Array {
  M data[N];
  //
  const M& operator[](Int n) const
  {
    assert(0 <= n && n < N);
    return data[n];
  }
  M& operator[](Int n)
  {
    assert(0 <= n && n < N);
    return data[n];
  }
  //
  const Array& operator=(const Array& v)
  {
    memcpy(this, &v, N * sizeof(M));
    return *this;
  }
  //
  Int size() const { return N; }
};

template <class M, Int N>
const Array<M, N>& operator+=(Array<M, N>& v, const Array<M, N>& v1)
{
  for (Int i = 0; i < N; ++i) {
    v.data[i] += v1.data[i];
  }
  return v;
}

template <class M, Int N>
const Array<M, N>& operator-=(Array<M, N>& v, const Array<M, N>& v1)
{
  for (Int i = 0; i < N; ++i) {
    v.data[i] -= v1.data[i];
  }
  return v;
}

template <class M, Int N>
const Array<M, N>& operator*=(Array<M, N>& v, RealD factor)
{
  for (Int i = 0; i < N; ++i) {
    v.data[i] *= factor;
  }
  return v;
}

template <class M, Int N>
const Array<M, N>& operator*=(Array<M, N>& v, ComplexD factor)
{
  for (Int i = 0; i < N; ++i) {
    v.data[i] *= factor;
  }
  return v;
}

template <class M, Int N>
Array<M, N> operator+(const Array<M, N>& v1, const Array<M, N>& v2)
{
  Array<M, N> v;
  for (Int i = 0; i < N; ++i) {
    v.data[i] = v1.data[i] + v2.data[i];
  }
  return v;
}

template <class M, Int N>
Array<M, N> operator-(const Array<M, N>& v1, const Array<M, N>& v2)
{
  Array<M, N> v;
  for (Int i = 0; i < N; ++i) {
    v.data[i] = v1.data[i] - v2.data[i];
  }
  return v;
}

template <class M, Int N>
Array<M, N> operator*(const Array<M, N>& v1, RealD factor)
{
  Array<M, N> v;
  for (Int i = 0; i < N; ++i) {
    v.data[i] = v1.data[i] * factor;
  }
  return v;
}

template <class M, Int N>
Array<M, N> operator*(const Array<M, N>& v1, ComplexD factor)
{
  Array<M, N> v;
  for (Int i = 0; i < N; ++i) {
    v.data[i] = v1.data[i] * factor;
  }
  return v;
}

template <class M, Int N>
Array<M, N> operator*(RealD factor, const Array<M, N>& v1)
{
  Array<M, N> v;
  for (Int i = 0; i < N; ++i) {
    v.data[i] = v1.data[i] * factor;
  }
  return v;
}

template <class M, Int N>
Array<M, N> operator*(ComplexD factor, const Array<M, N>& v1)
{
  Array<M, N> v;
  for (Int i = 0; i < N; ++i) {
    v.data[i] = v1.data[i] * factor;
  }
  return v;
}

template <class M, Int N>
void setZero(Array<M, N>& v)
{
  memset(v.data, 0, N * sizeof(M));
}

template <class M>
void setZero(std::vector<M>& v)
{
  memset(v.data(), 0, v.size() * sizeof(M));
}

void print(const ComplexD& v) { printf("%.16E:+%.16E", v.real(), v.imag()); }

template <class M, Int N>
void print(const Array<M, N>& v)
{
  printf("Array:");
  for (Int i = 0; i < N; ++i) {
    printf(" ");
    print(v[i]);
  }
  printf("\n");
}

struct SpinorMatrix : public Array<ComplexD, 2 * 2> {
};

inline bool notnan(const RealD& x) { return !isnan(x); }

inline bool notnan(const ComplexD& x)
{
  return notnan(x.real()) && notnan(x.imag());
}

template <class M>
inline bool notnanDoubles(const M& x)
{
  const RealD* v = (const RealD*)&x;
  for (Int i = 0; i < sizeof(M) / sizeof(RealD); ++i) {
    if (isnan(v[i])) {
      return false;
    }
  }
  return true;
}

inline bool notnan(const SpinorMatrix& x) { return notnanDoubles(x); }

template <class M>
inline bool notnan(const GridComm<M>& gc)
{
  const Geometry& geo = gc.getGeometry();
  for (Long index = 0; index < geo.localVolume(); ++index) {
    Int xl[4];
    geo.coordinateFromIndex(xl, index);
    for (Int m = 0; m < geo.multiplicity; ++m) {
      if (!notnan(gc.getElemV(xl, m))) {
        return false;
      }
    }
  }
  return true;
}

struct Coordinate {
 public:
  Int x[4];
  //
  Coordinate() { memset(x, 0, 4 * sizeof(int)); }
  Coordinate(Int x_, Int y_, Int z_, Int t_) { init(x_, y_, z_, t_); }
  Coordinate(const Int x_[4]) { init(x_); }
  void init(Int x_, Int y_, Int z_, Int t_)
  {
    x[0] = x_;
    x[1] = y_;
    x[2] = z_;
    x[3] = t_;
  }
  void init(const Int x_[4]) { memcpy(x, x_, 4 * sizeof(int)); }
  //
  Int& operator[](Int mu)
  {
    assert(0 <= mu && mu < 4);
    return x[mu];
  }
  const Int& operator[](Int mu) const
  {
    assert(0 <= mu && mu < 4);
    return x[mu];
  }
};

inline Coordinate operator+(const Coordinate& c1, const Coordinate& c2)
{
  Coordinate c;
  for (Int mu = 0; mu < 4; mu++) {
    c[mu] = c1[mu] + c2[mu];
  }
  return c;
}

inline Coordinate operator-(const Coordinate& c1, const Coordinate& c2)
{
  Coordinate c;
  for (Int mu = 0; mu < 4; mu++) {
    c[mu] = c1[mu] - c2[mu];
  }
  return c;
}

inline Coordinate operator%(const Coordinate& c1, const Coordinate& c2)
{
  Coordinate c;
  for (Int mu = 0; mu < 4; mu++) {
    c[mu] = c1[mu] % c2[mu];
  }
  return c;
}

inline bool operator==(const Coordinate& x, const Coordinate& y)
{
  return 0 == memcmp(x.x, y.x, sizeof(Coordinate));
}

inline bool operator!=(const Coordinate& x, const Coordinate& y)
{
  return 0 != memcmp(x.x, y.x, sizeof(Coordinate));
}

inline Int mod(const Int x, const Int len)
{
  assert(0 < len);
  const Int m = x % len;
  if (0 <= m) {
    return m;
  } else {
    return m + len;
  }
}

inline Int signMod(const Int x, const Int len)
{
  assert(0 < len);
  const Int m = mod(x, len);
  if (m * 2 < len) {
    return m;
  } else {
    return m - len;
  }
}

inline Int middleMod(const Int x, const Int y, const Int len)
{
  assert(0 < len);
  const Int xm = mod(x, len);
  const Int ym = mod(y, len);
  if (xm <= ym) {
    const Int r = signMod(ym - xm, len);
    return mod(xm + r / 2, len);
  } else {
    const Int r = signMod(xm - ym, len);
    return mod(ym + r / 2, len);
  }
}

inline RealD modf(const RealD x, const Int len)
{
  assert(0 < len);
  const RealD m = std::fmod(x, len);
  if (0 <= m) {
    return m;
  } else {
    return m + len;
  }
}

inline RealD signModf(const RealD x, const Int len)
{
  assert(0 < len);
  const RealD m = modf(x, len);
  if (m * 2 < len) {
    return m;
  } else {
    return m - len;
  }
}

inline RealD middleModf(const RealD x, const RealD y, const Int len)
{
  assert(0 < len);
  const RealD xm = modf(x, len);
  const RealD ym = modf(y, len);
  if (xm <= ym) {
    const RealD r = signModf(ym - xm, len);
    return modf(xm + r / 2, len);
  } else {
    const RealD r = signModf(xm - ym, len);
    return modf(ym + r / 2, len);
  }
}

inline void regularizeCoordinate(Int x[4], const Int size[4])
{
  x[0] = mod(x[0], size[0]);
  x[1] = mod(x[1], size[1]);
  x[2] = mod(x[2], size[2]);
  x[3] = mod(x[3], size[3]);
}

inline void regularizeCoordinateG(Int xg[4], const Geometry& geo)
{
  for (Int mu = 0; mu < 4; mu++) {
    xg[mu] = mod(xg[mu], geo.totalSite(mu));
  }
}

inline void relativeCoordinateG(Int xgrel[4], const Geometry& geo,
                                const Int xg[4], const Int xgref[4])
{
  for (Int mu = 0; mu < 4; mu++) {
    xgrel[mu] = signMod(xg[mu] - xgref[mu], geo.totalSite(mu));
  }
}

inline void absolueCoordinateG(Int xg[4], const Geometry& geo,
                               const Int xgrel[4], const Int xgref[4])
{
  for (Int mu = 0; mu < 4; mu++) {
    xg[mu] = mod(xgrel[mu] + xgref[mu], geo.totalSite(mu));
  }
}

inline void middleCoordinateG(Int xmg[4], const Geometry& geo, const Int x1g[4],
                              const Int x2g[4])
{
  for (Int mu = 0; mu < 4; mu++) {
    xmg[mu] = middleMod(x1g[mu], x2g[mu], geo.totalSite(mu));
  }
}

void coordinateRegularize(Coordinate& c, const Geometry& geo)
{
  regularizeCoordinateG(c.x, geo);
}

inline Coordinate mkCoordinateFromRefRel(const Geometry& geo,
                                         const Coordinate& c1,
                                         const Coordinate& c2)
{
  Coordinate c;
  absolueCoordinateG(c.x, geo, c1.x, c2.x);
  return c;
}

inline Long distance2RelativeCoordinateG(const Int xg[4])
{
  return sqr((Long)xg[0]) + sqr((Long)xg[1]) + sqr((Long)xg[2]) +
         sqr((Long)xg[3]);
}

inline RealD distanceRelativeCoordinateG(const Int xg[4])
{
  return sqrt(distance2RelativeCoordinateG(xg));
}

inline Long distance2TwoCoordinateG(const Int xg1[4], const Int xg2[4],
                                    const Geometry& geo)
{
  Int xgref[4];
  relativeCoordinateG(xgref, geo, xg1, xg2);
  return distance2RelativeCoordinateG(xgref);
}

inline Long distance2TwoCoordinateG(const Coordinate& xg1,
                                    const Coordinate& xg2, const Geometry& geo)
{
  Coordinate xgref;
  relativeCoordinateG(xgref.x, geo, xg1.x, xg2.x);
  return distance2RelativeCoordinateG(xgref.x);
}

inline void coordinateFromIndex(Int x[4], Long index, const Int size[4])
{
  x[0] = index % size[0];
  index /= size[0];
  x[1] = index % size[1];
  index /= size[1];
  x[2] = index % size[2];
  index /= size[2];
  x[3] = index % size[3];
}

inline Long indexFromCoordinate(const Int x[4], const Int size[4])
{
  return (((x[3] * size[2]) + x[2]) * size[1] + x[1]) * size[0] + x[0];
}

inline void sumDoubleArray(RealD* vs, const Long n_elem)
{
#ifdef USE_QMP
  QMP_sum_double_array(vs, n_elem);
#endif
}

inline Int sumArray(Long* recv, const Long* send, const Long n_elem)
{
#ifdef USE_QMP
  return mpi_allreduce((Long*)send, recv, n_elem, MPI_LONG, MPI_SUM,
                       QMP_COMM_WORLD);
#else
  memmove(recv, send, n_elem * sizeof(Long));
  return 0;
#endif
}

inline Int sumArray(RealD* recv, const RealD* send, const Long n_elem)
{
#ifdef USE_QMP
  return mpi_allreduce((RealD*)send, recv, n_elem, MPI_DOUBLE, MPI_SUM,
                       QMP_COMM_WORLD);
#else
  memmove(recv, send, n_elem * sizeof(RealD));
  return 0;
#endif
}

template <class M>
Int sumArray(M* vs, const Long n_elem)
{
  // M can be double or long
  Int status = 0;
#ifdef USE_QMP
  M tmp[n_elem];
  status = sumArray(tmp, vs, n_elem);
  memcpy(vs, tmp, n_elem * sizeof(M));
#endif
  return status;
}

inline void syncNode()
{
  Long v;
  sumArray(&v, 1);
}

inline void allGather(void* recvbuf, const void* sendbuf, const Long sendsize)
{
#ifdef USE_QMP
  MPI_Allgather((void*)sendbuf, sendsize, MPI_BYTE, recvbuf, sendsize, MPI_BYTE,
                QMP_COMM_WORLD);
#else
  assert(recvsize == sendsize);
  memmove(recvbuf, sendbuf, recvsize);
#endif
}

inline void bcast(void* recvbuf, const Long size, const Int root = 0)
{
#ifdef USE_QMP
  MPI_Bcast(recvbuf, size, MPI_BYTE, root, QMP_COMM_WORLD);
#else
#endif
}

struct GeometryNode {
  // About qmp geometry.
  Int numNode;
  // numNode = sizeNode[0] * sizeNode[1] * sizeNode[2] * sizeNode[3]
  Int idNode;
  // idNode = UniqueID()
  // 0 <= idNode < numNode
  Int sizeNode[4];
  Int coorNode[4];
  // 0 <= coorNode[i] < sizeNode[i]
  //
  void init()
  {
#ifdef USE_QMP
    MPI_Comm_size(QMP_COMM_WORLD, &numNode);
    MPI_Comm_rank(QMP_COMM_WORLD, &idNode);
    sizeNode[0] = SizeX();
    sizeNode[1] = SizeY();
    sizeNode[2] = SizeZ();
    sizeNode[3] = SizeT();
    coorNode[0] = CoorX();
    coorNode[1] = CoorY();
    coorNode[2] = CoorZ();
    coorNode[3] = CoorT();
    assert(indexFromCoordinate(coorNode, sizeNode) == idNode);
    assert(sizeNode[0] * sizeNode[1] * sizeNode[2] * sizeNode[3] == numNode);
#else
    numNode = 1;
    idNode = 0;
    for (Int i = 0; i < 4; i++) {
      sizeNode[i] = 1;
      coorNode[i] = 1;
    }
#endif
  }
  //
  GeometryNode(const bool initialize = false)
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

API inline const GeometryNode& getGeometryNode()
{
  static GeometryNode geon(true);
  return geon;
}

inline Int getNumNode() { return getGeometryNode().numNode; }

inline Int getIdNode() { return getGeometryNode().idNode; }

struct GeometryNodeNeighbor {
  Int dest[2][4];
  // dest[dir][mu]
  // dir = 0, 1 for Plus dir or Minus dir
  // 0 <= mu < 4 for different directions
  //
  void init()
  {
    const Int* coorNode = getGeometryNode().coorNode;
    const Int* sizeNode = getGeometryNode().sizeNode;
    for (Int mu = 0; mu < 4; mu++) {
      Int coor[4];
      memcpy(coor, coorNode, 4 * sizeof(int));
      coor[mu]++;
      regularizeCoordinate(coor, sizeNode);
      dest[0][mu] = indexFromCoordinate(coor, sizeNode);
      memcpy(coor, coorNode, 4 * sizeof(int));
      coor[mu]--;
      regularizeCoordinate(coor, sizeNode);
      dest[1][mu] = indexFromCoordinate(coor, sizeNode);
    }
  }
  //
  GeometryNodeNeighbor() { init(); }
};

API inline Int getDataDirMu(void* recv, void* send, const Long size, const Int dir,
                        const Int mu)
{
  // dir = 0, 1 for Plus dir or Minus dir
  // 0 <= mu < 4 for different directions
  // Old implementation is
  // getData((char*)recv, (char*)send, size, mu, 0 == dir ? 1 : -1);
  TIMER_FLOPS("getDataDirMu");
  timer.flops += size;
#ifdef USE_QMP
  static GeometryNodeNeighbor geonb;
  const Int idf = geonb.dest[dir][mu];
  const Int idt = geonb.dest[1 - dir][mu];
  MPI_Request req;
  MPI_Isend(send, size, MPI_BYTE, idt, 0, QMP_COMM_WORLD, &req);
  const Int ret = mpi_recv(recv, size, MPI_BYTE, idf, 0, QMP_COMM_WORLD);
  MPI_Wait(&req, MPI_STATUS_IGNORE);
  return ret;
#else
  memcpy(recv, send, size);
  return 0;
#endif
}

inline Int getDataPlusMu(void* recv, void* send, const Long size, const Int mu)
{
  return getDataDirMu(recv, send, size, 0, mu);
}

inline Int getDataMinusMu(void* recv, void* send, const Long size, const Int mu)
{
  return getDataDirMu(recv, send, size, 1, mu);
}

inline std::string emptyString() { return ""; }

template <class T>
std::string show(const T& x)
{
  std::ostringstream out;
  out << x;
  return out.str();
}

inline bool doesFileExist(const std::string& fn)
{
  struct stat sb;
  return 0 == stat(fn.c_str(), &sb);
}

inline bool DoesFileExist(const std::string& fn)
{
  Long nfile = 0;
  if (0 == getIdNode()) {
    if (doesFileExist(fn)) {
      nfile = 1;
    }
  }
  sumArray(&nfile, 1);
  return 0 != nfile;
}

API inline mode_t& defaultDirMode()
{
  static mode_t mode = 0775;
  return mode;
}

inline Int checkdir(const std::string& path,
                    const mode_t mode = defaultDirMode())
{
  TIMER("checkdir");
  Int ret = 0;
  while (!doesFileExist(path)) {
    ret = mkdir(path.c_str(), mode);
    sleep(0.001);
  }
  return ret;
}

inline Int makedir(const std::string& path,
                   const mode_t mode = defaultDirMode())
{
  TIMER("makedir");
  mkdir(path.c_str(), mode);
  return checkdir(path, mode);
}

inline Int Makedir(const std::string& path,
                   const mode_t mode = defaultDirMode())
{
  TIMER("Makedir");
  if (0 == getIdNode()) {
    makedir(path, mode);
  }
  return checkdir(path, mode);
}

inline Int Mkdir(const std::string& path, const mode_t mode = defaultDirMode())
{
  TIMER("Mkdir");
  Long ret = 0;
  if (0 == getIdNode()) {
    ret = mkdir(path.c_str(), mode);
  }
  sumArray(&ret, 1);
  return ret;
}

inline Int Rmdir(const std::string& path)
{
  TIMER("Rmdir");
  Long ret = 0;
  if (0 == getIdNode()) {
    ret = rmdir(path.c_str());
  }
  sumArray(&ret, 1);
  return ret;
}

API inline std::string& jobLock()
{
  static std::string lock;
  return lock;
}

inline Int sleep(const RealD seconds)
{
  return usleep((useconds_t)(seconds * 1.0e6));
}

inline void saveStringAllNode(const std::string& filename,
                              const std::string& content, bool isAppend = false)
{
  FILE* file = fopen(filename.c_str(), isAppend ? "a" : "w");
  while (NULL == file) {
    usleep(10000);
    file = fopen(filename.c_str(), isAppend ? "a" : "w");
  }
  fprintf(file, "%s", content.c_str());
  fclose(file);
}

inline void saveString(const std::string& filename, const std::string& content,
                       bool isAppend = false)
{
  if (0 == getIdNode()) {
    saveStringAllNode(filename, content, isAppend);
  }
  syncNode();
}

inline std::string showWilsonMatrix(const WilsonMatrix& mat,
                                    bool human_readable = false)
{
  Float* p = (Float*)&mat;
  const Int sizewm = sizeof(WilsonMatrix) / sizeof(WilsonVector);
  const Int sizewv = sizeof(WilsonVector) / sizeof(Float);
  std::ostringstream out;
  if (human_readable) {
    out.precision(1);
    out << std::fixed;
  } else {
    out.precision(16);
    out << std::scientific;
  }
  for (Int i = 0; i < sizewm; i++) {
    out << p[i * sizewv];
    for (Int j = 1; j < sizewv; j++) {
      out << " " << p[i * sizewv + j];
    }
    out << std::endl;
  }
  return out.str();
}

inline std::string showSpinMatrix(const WilsonMatrix& mat, const Int color1 = 0,
                                  const Int color2 = 0,
                                  bool human_readable = false)
{
  Float* p = (Float*)&mat;
  const Int sizewm = sizeof(WilsonMatrix) / sizeof(WilsonVector);
  const Int sizewv = sizeof(WilsonVector) / sizeof(Float);
  std::ostringstream out;
  if (human_readable) {
    out.precision(1);
    out << std::fixed;
  } else {
    out.precision(16);
    out << std::scientific;
  }
  for (Int i = color1; i < sizewm; i += 3) {
    for (Int j = color2 * 2; j < sizewv; j += 6) {
      out << std::setw(8) << p[i * sizewv + j] << " " << std::setw(8)
          << p[i * sizewv + j + 1] << " ";
    }
    out << std::endl;
  }
  return out.str();
}

inline std::string showSpinMatrix(const SpinMatrix& mat,
                                  bool human_readable = false)
{
  Float* p = (Float*)&mat;
  std::ostringstream out;
  if (human_readable) {
    out.precision(1);
    out << std::fixed;
  } else {
    out.precision(16);
    out << std::scientific;
  }
  for (Int i = 0; i < 4; i++) {
    for (Int j = 0; j < 4; j++) {
      out << std::setw(8) << mat(i, j).real() << " " << std::setw(8)
          << mat(i, j).imag() << " ";
    }
    out << std::endl;
  }
  return out.str();
}

inline void setZero(std::complex<RealD>& x) { x = 0.0; }

inline void setZero(Matrix& m) { memset(&m, 0, sizeof(m)); }

inline void setZero(WilsonMatrix& m) { memset(&m, 0, sizeof(m)); }

inline void setZero(SpinMatrix& m) { memset(&m, 0, sizeof(m)); }

inline void setUnit(ComplexD& x) { x = 1; }

inline void setUnit(WilsonMatrix& m)
{
  setZero(m);
  for (Int s = 0; s < 4; s++) {
    for (Int c = 0; c < 3; c++) {
      m(s, c, s, c) = 1.0;
    }
  }
}

inline void setUnit(SpinMatrix& m)
{
  setZero(m);
  for (Int s = 0; s < 4; s++) {
    m(s, s) = 1.0;
  }
}

inline void multiplyPlus(SpinMatrix& m, const SpinMatrix& m1,
                         const SpinMatrix& m2)
{
  const Int dim = 4;
  for (Int i = 0; i < dim; i++) {
    for (Int j = 0; j < dim; j++) {
      for (Int k = 0; k < dim; k++) {
        m(i, j) += m1(i, k) * m2(k, j);
      }
    }
  }
}

inline SpinMatrix gammaInit(const Int dir)
{
  static const char* fname = "gammaInit(sm)";
  VRB.Result(cname, fname, "dir=%d\n", dir);
  SpinMatrix m;
  setZero(m);
  switch (dir) {
    case 0:
      m(0, 3) += ii;
      m(1, 2) += ii;
      m(2, 1) -= ii;
      m(3, 0) -= ii;
      break;
    case 1:
      m(0, 3) -= 1;
      m(1, 2) += 1;
      m(2, 1) += 1;
      m(3, 0) -= 1;
      break;
    case 2:
      m(0, 2) += ii;
      m(1, 3) -= ii;
      m(2, 0) -= ii;
      m(3, 1) += ii;
      break;
    case 3:
      m(0, 2) += 1;
      m(1, 3) += 1;
      m(2, 0) += 1;
      m(3, 1) += 1;
      break;
    case -5:
      m(0, 0) += 1;
      m(1, 1) += 1;
      m(2, 2) -= 1;
      m(3, 3) -= 1;
      break;
    default:
      ERR.General(cname, fname, "dir=%d\n", dir);
      break;
  };
  return m;
}

API inline const SpinMatrix& gamma(const Int dir)
{
  static const char* fname = "gamma(sm)";
  static SpinMatrix gamma0 = gammaInit(0);
  static SpinMatrix gamma1 = gammaInit(1);
  static SpinMatrix gamma2 = gammaInit(2);
  static SpinMatrix gamma3 = gammaInit(3);
  static SpinMatrix gamma5 = gammaInit(-5);
  switch (dir) {
    case 0:
      return gamma0;
    case 1:
      return gamma1;
    case 2:
      return gamma2;
    case 3:
      return gamma3;
    case -5:
      return gamma5;
    default:
      ERR.General(cname, fname, "dir=%d\n", dir);
      break;
  }
  return gamma5;
}

inline SpinMatrix sigmaInit(const Int dir)
{
  TIMER_VERBOSE("sigmaInit");
  switch (dir) {
    case 0:
      return (gamma(1) * gamma(2) - gamma(2) * gamma(1)) * (1.0 / (2.0 * ii));
    case 1:
      return (gamma(2) * gamma(0) - gamma(0) * gamma(2)) * (1.0 / (2.0 * ii));
    case 2:
      return (gamma(0) * gamma(1) - gamma(1) * gamma(0)) * (1.0 / (2.0 * ii));
    default:
      assert(false);
  }
}

API inline const SpinMatrix& sigma(const Int dir)
{
  static SpinMatrix sigma0 = sigmaInit(0);
  static SpinMatrix sigma1 = sigmaInit(1);
  static SpinMatrix sigma2 = sigmaInit(2);
  switch (dir) {
    case 0:
      return sigma0;
    case 1:
      return sigma1;
    case 2:
      return sigma2;
    default:
      assert(false);
  }
}

inline void gammaLeftPlus(SpinMatrix& m, const SpinMatrix& m1, const Int dir)
{
  multiplyPlus(m, gamma(dir), m1);
}

inline void gammaRightPlus(SpinMatrix& m, const SpinMatrix& m1, const Int dir)
{
  multiplyPlus(m, m1, gamma(dir));
}

inline void gammaLeft(WilsonMatrix& m, const Int dir) { m.gl(dir); }

inline void gammaLeft(SpinMatrix& m, const Int dir)
{
  SpinMatrix m1 = m;
  setZero(m);
  gammaLeftPlus(m, m1, dir);
}

inline ComplexD trace(const WilsonMatrix& m)
{
  ComplexD t = 0;
  for (Int s = 0; s < 4; s++) {
    for (Int c = 0; c < 3; c++) {
      t += m(s, c, s, c);
    }
  }
  return t;
}

inline ComplexD trace(const SpinMatrix& m)
{
  ComplexD t = 0;
  for (Int s = 0; s < 4; s++) {
    t += m(s, s);
  }
  return t;
}

CPS_END_NAMESPACE
