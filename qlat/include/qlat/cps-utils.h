#pragma once

CPS_START_NAMESPACE

const std::complex<double> ii(0.0, 1.0);

template <class M>
bool isInitialized(const GridComm<M>& gc)
{
  return isInitialized(gc.getGeometry());
}

bool isInitialized(const Geometry& geo) { return geo.initialized; }

template <class T>
inline void setZeroN(T xs[], const int xn)
{
  for (int i = 0; i < xn; i++) {
    setZero(xs[i]);
  }
}

template <class M, int N>
struct Array {
  M data[N];
  //
  const M& operator[](int n) const
  {
    assert(0 <= n && n < N);
    return data[n];
  }
  M& operator[](int n)
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
  int size() const { return N; }
};

template <class M, int N>
const Array<M, N>& operator+=(Array<M, N>& v, const Array<M, N>& v1)
{
  for (int i = 0; i < N; ++i) {
    v.data[i] += v1.data[i];
  }
  return v;
}

template <class M, int N>
const Array<M, N>& operator-=(Array<M, N>& v, const Array<M, N>& v1)
{
  for (int i = 0; i < N; ++i) {
    v.data[i] -= v1.data[i];
  }
  return v;
}

template <class M, int N>
const Array<M, N>& operator*=(Array<M, N>& v, double factor)
{
  for (int i = 0; i < N; ++i) {
    v.data[i] *= factor;
  }
  return v;
}

template <class M, int N>
const Array<M, N>& operator*=(Array<M, N>& v, Complex factor)
{
  for (int i = 0; i < N; ++i) {
    v.data[i] *= factor;
  }
  return v;
}

template <class M, int N>
Array<M, N> operator+(const Array<M, N>& v1, const Array<M, N>& v2)
{
  Array<M, N> v;
  for (int i = 0; i < N; ++i) {
    v.data[i] = v1.data[i] + v2.data[i];
  }
  return v;
}

template <class M, int N>
Array<M, N> operator-(const Array<M, N>& v1, const Array<M, N>& v2)
{
  Array<M, N> v;
  for (int i = 0; i < N; ++i) {
    v.data[i] = v1.data[i] - v2.data[i];
  }
  return v;
}

template <class M, int N>
Array<M, N> operator*(const Array<M, N>& v1, double factor)
{
  Array<M, N> v;
  for (int i = 0; i < N; ++i) {
    v.data[i] = v1.data[i] * factor;
  }
  return v;
}

template <class M, int N>
Array<M, N> operator*(const Array<M, N>& v1, Complex factor)
{
  Array<M, N> v;
  for (int i = 0; i < N; ++i) {
    v.data[i] = v1.data[i] * factor;
  }
  return v;
}

template <class M, int N>
Array<M, N> operator*(double factor, const Array<M, N>& v1)
{
  Array<M, N> v;
  for (int i = 0; i < N; ++i) {
    v.data[i] = v1.data[i] * factor;
  }
  return v;
}

template <class M, int N>
Array<M, N> operator*(Complex factor, const Array<M, N>& v1)
{
  Array<M, N> v;
  for (int i = 0; i < N; ++i) {
    v.data[i] = v1.data[i] * factor;
  }
  return v;
}

template <class M, int N>
void setZero(Array<M, N>& v)
{
  memset(v.data, 0, N * sizeof(M));
}

template <class M>
void setZero(std::vector<M>& v)
{
  memset(v.data(), 0, v.size() * sizeof(M));
}

void print(const Complex& v) { printf("%.16E:+%.16E", v.real(), v.imag()); }

template <class M, int N>
void print(const Array<M, N>& v)
{
  printf("Array:");
  for (int i = 0; i < N; ++i) {
    printf(" ");
    print(v[i]);
  }
  printf("\n");
}

struct SpinorMatrix : public Array<Complex, 2 * 2> {
};

inline bool notnan(const double& x) { return !isnan(x); }

inline bool notnan(const Complex& x)
{
  return notnan(x.real()) && notnan(x.imag());
}

template <class M>
inline bool notnanDoubles(const M& x)
{
  const double* v = (const double*)&x;
  for (int i = 0; i < sizeof(M) / sizeof(double); ++i) {
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
  for (long index = 0; index < geo.localVolume(); ++index) {
    int xl[4];
    geo.coordinateFromIndex(xl, index);
    for (int m = 0; m < geo.multiplicity; ++m) {
      if (!notnan(gc.getElemV(xl, m))) {
        return false;
      }
    }
  }
  return true;
}

class Coordinate
{
 public:
  int x[4];
  //
  Coordinate() { memset(x, 0, 4 * sizeof(int)); }
  Coordinate(int x_, int y_, int z_, int t_) { init(x_, y_, z_, t_); }
  Coordinate(const int x_[4]) { init(x_); }
  void init(int x_, int y_, int z_, int t_)
  {
    x[0] = x_;
    x[1] = y_;
    x[2] = z_;
    x[3] = t_;
  }
  void init(const int x_[4]) { memcpy(x, x_, 4 * sizeof(int)); }
  //
  int& operator[](int mu)
  {
    assert(0 <= mu && mu < 4);
    return x[mu];
  }
  const int& operator[](int mu) const
  {
    assert(0 <= mu && mu < 4);
    return x[mu];
  }
};

inline Coordinate operator+(const Coordinate& c1, const Coordinate& c2)
{
  Coordinate c;
  for (int mu = 0; mu < 4; mu++) {
    c[mu] = c1[mu] + c2[mu];
  }
  return c;
}

inline Coordinate operator-(const Coordinate& c1, const Coordinate& c2)
{
  Coordinate c;
  for (int mu = 0; mu < 4; mu++) {
    c[mu] = c1[mu] - c2[mu];
  }
  return c;
}

inline Coordinate operator%(const Coordinate& c1, const Coordinate& c2)
{
  Coordinate c;
  for (int mu = 0; mu < 4; mu++) {
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

inline int mod(const int x, const int len)
{
  assert(0 < len);
  const int m = x % len;
  if (0 <= m) {
    return m;
  } else {
    return m + len;
  }
}

inline int signMod(const int x, const int len)
{
  assert(0 < len);
  const int m = mod(x, len);
  if (m * 2 < len) {
    return m;
  } else {
    return m - len;
  }
}

inline int middleMod(const int x, const int y, const int len)
{
  assert(0 < len);
  const int xm = mod(x, len);
  const int ym = mod(y, len);
  if (xm <= ym) {
    const int r = signMod(ym - xm, len);
    return mod(xm + r / 2, len);
  } else {
    const int r = signMod(xm - ym, len);
    return mod(ym + r / 2, len);
  }
}

inline double modf(const double x, const int len)
{
  assert(0 < len);
  const double m = std::fmod(x, len);
  if (0 <= m) {
    return m;
  } else {
    return m + len;
  }
}

inline double signModf(const double x, const int len)
{
  assert(0 < len);
  const double m = modf(x, len);
  if (m * 2 < len) {
    return m;
  } else {
    return m - len;
  }
}

inline double middleModf(const double x, const double y, const int len)
{
  assert(0 < len);
  const double xm = modf(x, len);
  const double ym = modf(y, len);
  if (xm <= ym) {
    const double r = signModf(ym - xm, len);
    return modf(xm + r / 2, len);
  } else {
    const double r = signModf(xm - ym, len);
    return modf(ym + r / 2, len);
  }
}

inline void regularizeCoordinate(int x[4], const int size[4])
{
  x[0] = mod(x[0], size[0]);
  x[1] = mod(x[1], size[1]);
  x[2] = mod(x[2], size[2]);
  x[3] = mod(x[3], size[3]);
}

inline void regularizeCoordinateG(int xg[4], const Geometry& geo)
{
  for (int mu = 0; mu < 4; mu++) {
    xg[mu] = mod(xg[mu], geo.totalSite(mu));
  }
}

inline void relativeCoordinateG(int xgrel[4], const Geometry& geo,
                                const int xg[4], const int xgref[4])
{
  for (int mu = 0; mu < 4; mu++) {
    xgrel[mu] = signMod(xg[mu] - xgref[mu], geo.totalSite(mu));
  }
}

inline void absolueCoordinateG(int xg[4], const Geometry& geo,
                               const int xgrel[4], const int xgref[4])
{
  for (int mu = 0; mu < 4; mu++) {
    xg[mu] = mod(xgrel[mu] + xgref[mu], geo.totalSite(mu));
  }
}

inline void middleCoordinateG(int xmg[4], const Geometry& geo, const int x1g[4],
                              const int x2g[4])
{
  for (int mu = 0; mu < 4; mu++) {
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

inline long distance2RelativeCoordinateG(const int xg[4])
{
  return sqr((long)xg[0]) + sqr((long)xg[1]) + sqr((long)xg[2]) +
         sqr((long)xg[3]);
}

inline double distanceRelativeCoordinateG(const int xg[4])
{
  return sqrt(distance2RelativeCoordinateG(xg));
}

inline long distance2TwoCoordinateG(const int xg1[4], const int xg2[4],
                                    const Geometry& geo)
{
  int xgref[4];
  relativeCoordinateG(xgref, geo, xg1, xg2);
  return distance2RelativeCoordinateG(xgref);
}

inline long distance2TwoCoordinateG(const Coordinate& xg1,
                                    const Coordinate& xg2, const Geometry& geo)
{
  Coordinate xgref;
  relativeCoordinateG(xgref.x, geo, xg1.x, xg2.x);
  return distance2RelativeCoordinateG(xgref.x);
}

inline void coordinateFromIndex(int x[4], long index, const int size[4])
{
  x[0] = index % size[0];
  index /= size[0];
  x[1] = index % size[1];
  index /= size[1];
  x[2] = index % size[2];
  index /= size[2];
  x[3] = index % size[3];
}

inline long indexFromCoordinate(const int x[4], const int size[4])
{
  return (((x[3] * size[2]) + x[2]) * size[1] + x[1]) * size[0] + x[0];
}

inline void sumDoubleArray(double* vs, const long n_elem)
{
#ifdef USE_QMP
  QMP_sum_double_array(vs, n_elem);
#endif
}

inline int sumArray(long* recv, const long* send, const long n_elem)
{
#ifdef USE_QMP
  return MPI_Allreduce((long*)send, recv, n_elem, MPI_LONG, MPI_SUM,
                       QMP_COMM_WORLD);
#else
  memmove(recv, send, n_elem * sizeof(long));
  return 0;
#endif
}

inline int sumArray(double* recv, const double* send, const long n_elem)
{
#ifdef USE_QMP
  return MPI_Allreduce((double*)send, recv, n_elem, MPI_DOUBLE, MPI_SUM,
                       QMP_COMM_WORLD);
#else
  memmove(recv, send, n_elem * sizeof(double));
  return 0;
#endif
}

template <class M>
int sumArray(M* vs, const long n_elem)
{
  // M can be double or long
  int status = 0;
#ifdef USE_QMP
  M tmp[n_elem];
  status = sumArray(tmp, vs, n_elem);
  memcpy(vs, tmp, n_elem * sizeof(M));
#endif
  return status;
}

inline void syncNode()
{
  long v;
  sumArray(&v, 1);
}

inline void allGather(void* recvbuf, const void* sendbuf, const long sendsize)
{
#ifdef USE_QMP
  MPI_Allgather((void*)sendbuf, sendsize, MPI_BYTE, recvbuf, sendsize, MPI_BYTE,
                QMP_COMM_WORLD);
#else
  assert(recvsize == sendsize);
  memmove(recvbuf, sendbuf, recvsize);
#endif
}

inline void bcast(void* recvbuf, const long size, const int root = 0)
{
#ifdef USE_QMP
  MPI_Bcast(recvbuf, size, MPI_BYTE, root, QMP_COMM_WORLD);
#else
#endif
}

struct GeometryNode {
  // About qmp geometry.
  int numNode;
  // numNode = sizeNode[0] * sizeNode[1] * sizeNode[2] * sizeNode[3]
  int idNode;
  // idNode = UniqueID()
  // 0 <= idNode < numNode
  int sizeNode[4];
  int coorNode[4];
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
    for (int i = 0; i < 4; i++) {
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

inline int getNumNode() { return getGeometryNode().numNode; }

inline int getIdNode() { return getGeometryNode().idNode; }

struct GeometryNodeNeighbor {
  int dest[2][4];
  // dest[dir][mu]
  // dir = 0, 1 for Plus dir or Minus dir
  // 0 <= mu < 4 for different directions
  //
  void init()
  {
    const int* coorNode = getGeometryNode().coorNode;
    const int* sizeNode = getGeometryNode().sizeNode;
    for (int mu = 0; mu < 4; mu++) {
      int coor[4];
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

API inline int getDataDirMu(void* recv, void* send, const long size, const int dir,
                        const int mu)
{
  // dir = 0, 1 for Plus dir or Minus dir
  // 0 <= mu < 4 for different directions
  // Old implementation is
  // getData((char*)recv, (char*)send, size, mu, 0 == dir ? 1 : -1);
  TIMER_FLOPS("getDataDirMu");
  timer.flops += size;
#ifdef USE_QMP
  static GeometryNodeNeighbor geonb;
  const int idf = geonb.dest[dir][mu];
  const int idt = geonb.dest[1 - dir][mu];
  MPI_Request req;
  MPI_Isend(send, size, MPI_BYTE, idt, 0, QMP_COMM_WORLD, &req);
  const int ret =
      MPI_Recv(recv, size, MPI_BYTE, idf, 0, QMP_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Wait(&req, MPI_STATUS_IGNORE);
  return ret;
#else
  memcpy(recv, send, size);
  return 0;
#endif
}

inline int getDataPlusMu(void* recv, void* send, const long size, const int mu)
{
  return getDataDirMu(recv, send, size, 0, mu);
}

inline int getDataMinusMu(void* recv, void* send, const long size, const int mu)
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
  long nfile = 0;
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

inline int checkdir(const std::string& path,
                    const mode_t mode = defaultDirMode())
{
  TIMER("checkdir");
  int ret = 0;
  while (!doesFileExist(path)) {
    ret = mkdir(path.c_str(), mode);
    sleep(0.001);
  }
  return ret;
}

inline int makedir(const std::string& path,
                   const mode_t mode = defaultDirMode())
{
  TIMER("makedir");
  mkdir(path.c_str(), mode);
  return checkdir(path, mode);
}

inline int Makedir(const std::string& path,
                   const mode_t mode = defaultDirMode())
{
  TIMER("Makedir");
  if (0 == getIdNode()) {
    makedir(path, mode);
  }
  return checkdir(path, mode);
}

inline int Mkdir(const std::string& path, const mode_t mode = defaultDirMode())
{
  TIMER("Mkdir");
  long ret = 0;
  if (0 == getIdNode()) {
    ret = mkdir(path.c_str(), mode);
  }
  sumArray(&ret, 1);
  return ret;
}

inline int Rmdir(const std::string& path)
{
  TIMER("Rmdir");
  long ret = 0;
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

inline int sleep(const double seconds)
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
  const int sizewm = sizeof(WilsonMatrix) / sizeof(WilsonVector);
  const int sizewv = sizeof(WilsonVector) / sizeof(Float);
  std::ostringstream out;
  if (human_readable) {
    out.precision(1);
    out << std::fixed;
  } else {
    out.precision(16);
    out << std::scientific;
  }
  for (int i = 0; i < sizewm; i++) {
    out << p[i * sizewv];
    for (int j = 1; j < sizewv; j++) {
      out << " " << p[i * sizewv + j];
    }
    out << std::endl;
  }
  return out.str();
}

inline std::string showSpinMatrix(const WilsonMatrix& mat, const int color1 = 0,
                                  const int color2 = 0,
                                  bool human_readable = false)
{
  Float* p = (Float*)&mat;
  const int sizewm = sizeof(WilsonMatrix) / sizeof(WilsonVector);
  const int sizewv = sizeof(WilsonVector) / sizeof(Float);
  std::ostringstream out;
  if (human_readable) {
    out.precision(1);
    out << std::fixed;
  } else {
    out.precision(16);
    out << std::scientific;
  }
  for (int i = color1; i < sizewm; i += 3) {
    for (int j = color2 * 2; j < sizewv; j += 6) {
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
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      out << std::setw(8) << mat(i, j).real() << " " << std::setw(8)
          << mat(i, j).imag() << " ";
    }
    out << std::endl;
  }
  return out.str();
}

inline void setZero(std::complex<double>& x) { x = 0.0; }

inline void setZero(Matrix& m) { memset(&m, 0, sizeof(m)); }

inline void setZero(WilsonMatrix& m) { memset(&m, 0, sizeof(m)); }

inline void setZero(SpinMatrix& m) { memset(&m, 0, sizeof(m)); }

inline void setUnit(Complex& x) { x = 1; }

inline void setUnit(WilsonMatrix& m)
{
  setZero(m);
  for (int s = 0; s < 4; s++) {
    for (int c = 0; c < 3; c++) {
      m(s, c, s, c) = 1.0;
    }
  }
}

inline void setUnit(SpinMatrix& m)
{
  setZero(m);
  for (int s = 0; s < 4; s++) {
    m(s, s) = 1.0;
  }
}

inline void multiplyPlus(SpinMatrix& m, const SpinMatrix& m1,
                         const SpinMatrix& m2)
{
  const int dim = 4;
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      for (int k = 0; k < dim; k++) {
        m(i, j) += m1(i, k) * m2(k, j);
      }
    }
  }
}

inline SpinMatrix gammaInit(const int dir)
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

API inline const SpinMatrix& gamma(const int dir)
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

inline SpinMatrix sigmaInit(const int dir)
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

API inline const SpinMatrix& sigma(const int dir)
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

inline void gammaLeftPlus(SpinMatrix& m, const SpinMatrix& m1, const int dir)
{
  multiplyPlus(m, gamma(dir), m1);
}

inline void gammaRightPlus(SpinMatrix& m, const SpinMatrix& m1, const int dir)
{
  multiplyPlus(m, m1, gamma(dir));
}

inline void gammaLeft(WilsonMatrix& m, const int dir) { m.gl(dir); }

inline void gammaLeft(SpinMatrix& m, const int dir)
{
  SpinMatrix m1 = m;
  setZero(m);
  gammaLeftPlus(m, m1, dir);
}

inline Complex trace(const WilsonMatrix& m)
{
  Complex t = 0;
  for (int s = 0; s < 4; s++) {
    for (int c = 0; c < 3; c++) {
      t += m(s, c, s, c);
    }
  }
  return t;
}

inline Complex trace(const SpinMatrix& m)
{
  Complex t = 0;
  for (int s = 0; s < 4; s++) {
    t += m(s, s);
  }
  return t;
}

CPS_END_NAMESPACE
