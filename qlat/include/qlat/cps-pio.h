#pragma once

CPS_START_NAMESPACE

inline std::string showIntWP(const int n, const int w, const char pad)
{
  std::ostringstream out;
  out << std::setfill(pad) << std::setw(w) << n;
  return out.str();
}

struct DataSizeNmemb {
  void* data;
  long size;
  long nmemb;
  //
  void init()
  {
    data = NULL;
    size = 0;
    nmemb = 0;
  }
  void init(void* data_, const long size_, const long nmemb_)
  {
    data = data_;
    size = size_;
    nmemb = nmemb_;
  }
  //
  DataSizeNmemb() { init(); }
  DataSizeNmemb(void* data_, const long size_, const long nmemb_)
  {
    init(data_, size_, nmemb_);
  }
};

inline qlat::crc32_t dataCRC32(qlat::crc32_t* chksums,
                               const std::vector<DataSizeNmemb>& dsns)
{
  // Length of chksums should be getNumNode()
  TIMER_FLOPS("dataCRC32");
  for (int k = 0; k < dsns.size(); k++) {
    timer.flops += dsns[k].size * dsns[k].nmemb;
  }
  qlat::crc32_t chksums1[getNumNode()];
  if (NULL == chksums) {
    chksums = chksums1;
  }
  memset(chksums, 0, getNumNode() * sizeof(qlat::crc32_t));
  qlat::crc32_t chksum = 0;
  for (int k = 0; k < dsns.size(); k++) {
    chksum = qlat::crc32(chksum, dsns[k].data, dsns[k].size * dsns[k].nmemb);
  }
  allGather(chksums, &chksum, sizeof(qlat::crc32_t));
  return qlat::crc32(chksums, getNumNode() * sizeof(qlat::crc32_t));
}

inline qlat::crc32_t dataCRC32(qlat::crc32_t* chksums, const void* data,
                               const long size, const long nmemb)
{
  // Length of chksums should be getNumNode()
  DataSizeNmemb dsn((void*)data, size, nmemb);
  std::vector<DataSizeNmemb> dsns;
  dsns.push_back(dsn);
  return dataCRC32(chksums, dsns);
}

template <class M>
qlat::crc32_t gcCRC32(qlat::crc32_t* chksums, const GridComm<M>& gc)
{
  // Length of chksums should be getNumNode()
  TIMER_FLOPS("gcCRC32");
  timer.flops += gc.getAllocSize();
  const Geometry& geo = gc.getGeometry();
  assert(geo.isOnlyLocal());
  return dataCRC32(chksums, gc.getField(), sizeof(M) * geo.multiplicity,
                   geo.localVolume());
}

template <class M>
void gcWriteInfo(const GridComm<M>& gc, const std::string& path)
{
  TIMER("gcWriteInfo");
  const Geometry& geo = gc.getGeometry();
  assert(geo.isOnlyLocal());
  if (0 == getIdNode()) {
    const std::string filename = path + "/geo-info.txt";
    FILE* file = fopen(filename.c_str(), "w");
    fprintf(file, "nodeFileSize = %ld\n", gc.getAllocSize());
    fprintf(file, "geo.multiplicity = %d\n", geo.multiplicity);
    fprintf(file, "sizeof(M) = %d\n", sizeof(M));
    fprintf(file, "geo.numNode = %d\n", geo.numNode);
    fprintf(file, "geo.sizeNode[0] = %d\n", geo.sizeNode[0]);
    fprintf(file, "geo.sizeNode[1] = %d\n", geo.sizeNode[1]);
    fprintf(file, "geo.sizeNode[2] = %d\n", geo.sizeNode[2]);
    fprintf(file, "geo.sizeNode[3] = %d\n", geo.sizeNode[3]);
    fprintf(file, "geo.localVolume() = %ld\n", geo.localVolume());
    fprintf(file, "geo.nodeSite[0] = %d\n", geo.nodeSite[0]);
    fprintf(file, "geo.nodeSite[1] = %d\n", geo.nodeSite[1]);
    fprintf(file, "geo.nodeSite[2] = %d\n", geo.nodeSite[2]);
    fprintf(file, "geo.nodeSite[3] = %d\n", geo.nodeSite[3]);
    fprintf(file, "PI = %.20f\n", PI);
    const char* pic = (const char*)&PI;
    fprintf(file, "PI_double = %hhx %hhx %hhx %hhx %hhx %hhx %hhx %hhx\n",
            pic[0], pic[1], pic[2], pic[3], pic[4], pic[5], pic[6], pic[7]);
    const float PIf = PI;
    const char* pifc = (const char*)&PIf;
    fprintf(file, "PI_float = %hhx %hhx %hhx %hhx\n", pifc[0], pifc[1], pifc[2],
            pifc[3]);
    fclose(file);
  }
}

inline void dataWriteInfo(const std::vector<DataSizeNmemb>& dsns,
                          const std::string& path)
{
  TIMER_FLOPS("dataWriteInfo");
  for (int k = 0; k < dsns.size(); k++) {
    timer.flops += dsns[k].size * dsns[k].nmemb;
  }
  qlat::crc32_t chksums[getNumNode()];
  qlat::crc32_t chksum = dataCRC32(chksums, dsns);
  if (0 == getIdNode()) {
    const std::string filename = path + "/checksums.txt";
    FILE* file = fopen(filename.c_str(), "w");
    fprintf(file, "%X\n", chksum);
    fprintf(file, "\n");
    for (int i = 0; i < getNumNode(); i++) {
      fprintf(file, "%X\n", chksums[i]);
    }
    fclose(file);
  }
}

inline void dataReadInfo(const std::vector<DataSizeNmemb>& dsns,
                         const std::string& path)
{
  TIMER_FLOPS("dataReadInfo");
  for (int k = 0; k < dsns.size(); k++) {
    timer.flops += dsns[k].size * dsns[k].nmemb;
  }
  qlat::crc32_t chksums[getNumNode()];
  qlat::crc32_t chksum = dataCRC32(chksums, dsns);
  if (0 == getIdNode()) {
    qlat::crc32_t chksum_read;
    const std::string filename = path + "/checksums.txt";
    FILE* file = fopen(filename.c_str(), "r");
    qassert(file != NULL);
    fscanf(file, "%X\n", &chksum_read);
    fscanf(file, "\n");
    for (int i = 0; i < getNumNode(); i++) {
      fscanf(file, "%X\n", &chksum_read);
      qassert(chksums[i] == chksum_read);
    }
    if (chksum != chksum_read) {
      qlat::displayln(
          fname +
          qlat::ssprintf(": WARNING mismatch total chksum=%08X chksum_read=%08X",
                   chksum, chksum_read));
    }
    fclose(file);
  }
}

const int DATA_READ_WRITE_FILENAME_WIDTH = 10;
const int DATA_READ_WRITE_NUMBER_OF_DIRECTORIES = 32;

API inline int& dataWriteParNumber()
{
  static int npar = 3;
  return npar;
}

API inline int& dataReadParNumber()
{
  static int npar = 3;
  return npar;
}

inline long dataWriteParNode(const std::vector<DataSizeNmemb>& dsns,
                             const std::string& path,
                             const mode_t mode = defaultDirMode())
{
  TIMER_VERBOSE_FLOPS("dataWriteParNode");
  for (int k = 0; k < dsns.size(); k++) {
    timer.flops += dsns[k].size * dsns[k].nmemb * getNumNode();
  }
  const int size_ndir =
      std::min(DATA_READ_WRITE_NUMBER_OF_DIRECTORIES, getNumNode());
  Makedir(path.c_str(), mode);
  if (getIdNode() < DATA_READ_WRITE_NUMBER_OF_DIRECTORIES) {
    makedir((path + "/" + showIntWP(getIdNode(), 2, '0')).c_str(), mode);
  }
  const int dir_size = (getNumNode() - 1) / size_ndir + 1;
  const int idDir = getIdNode() / dir_size;
  assert(0 <= idDir && idDir < size_ndir);
  const std::string pathId = path + "/" + showIntWP(idDir, 2, '0');
  checkdir(pathId.c_str(), mode);
  const std::string filename =
      pathId + "/" +
      showIntWP(getIdNode(), DATA_READ_WRITE_FILENAME_WIDTH, '0');
  const int n_cycle = std::max(1, getNumNode() / dataWriteParNumber());
  long total_bytes = 0;
  {
    TIMER_FLOPS("dataWriteParNode-fwrite");
    for (int k = 0; k < dsns.size(); k++) {
      timer.flops += dsns[k].size * dsns[k].nmemb * getNumNode();
    }
    for (int i = 0; i < n_cycle; i++) {
      TIMER_VERBOSE_FLOPS("dataWriteParNode-fwrite-flops");
      long bytes = 0;
      if (getIdNode() % n_cycle == i) {
        FILE* file = fopen(filename.c_str(), "w");
        while (NULL == file) {
          makedir(pathId.c_str(), mode);
          file = fopen(filename.c_str(), "w");
        }
        for (int k = 0; k < dsns.size(); k++) {
          bytes += dsns[k].size *
                   fwrite(dsns[k].data, dsns[k].size, dsns[k].nmemb, file);
        }
        fclose(file);
      }
      sumArray(&bytes, 1);
      total_bytes += bytes;
      qlat::DisplayInfo(cname, fname.c_str(),
                  "cycle / n_cycle = %d / %d ; total_bytes = %ld\n", i + 1,
                  n_cycle, total_bytes);
      timer.flops += bytes;
    }
  }
  dataWriteInfo(dsns, path);
  saveString(path + "/checkpoint", "");
  return total_bytes;
}

inline long dataWriteParNode(const void* data, const long size,
                             const long nmemb, const std::string& path,
                             const mode_t mode = defaultDirMode())
{
  DataSizeNmemb dsn((void*)data, size, nmemb);
  std::vector<DataSizeNmemb> dsns;
  dsns.push_back(dsn);
  return dataWriteParNode(dsns, path, mode);
}

template <class M>
long gcWriteParNode(const GridComm<M>& gc, const std::string& path,
                    const mode_t mode = defaultDirMode())
{
  TIMER_VERBOSE_FLOPS("gcWriteParNode");
  const Geometry& geo = gc.getGeometry();
  assert(geo.isOnlyLocal());
  long total_bytes =
      dataWriteParNode(gc.getField(), sizeof(M) * geo.multiplicity,
                       geo.localVolume(), path, mode);
  gcWriteInfo(gc, path);
  timer.flops += total_bytes;
  return total_bytes;
}

inline long dataReadParNode(const std::vector<DataSizeNmemb>& dsns,
                            const std::string& path)
{
  if (!DoesFileExist(path + "/checkpoint")) {
    qlat::DisplayInfo(cname, "dataReadParNode", "'%s' do not exist.\n", path.c_str());
    return 0;
  }
  TIMER_VERBOSE_FLOPS("dataReadParNode");
  for (int k = 0; k < dsns.size(); k++) {
    timer.flops += dsns[k].size * dsns[k].nmemb * getNumNode();
  }
  const int size_ndir =
      std::min(DATA_READ_WRITE_NUMBER_OF_DIRECTORIES, getNumNode());
  const int dir_size = (getNumNode() - 1) / size_ndir + 1;
  const int idDir = getIdNode() / dir_size;
  assert(0 <= idDir && idDir < size_ndir);
  const std::string pathId = path + "/" + showIntWP(idDir, 2, '0');
  const std::string filename =
      pathId + "/" +
      showIntWP(getIdNode(), DATA_READ_WRITE_FILENAME_WIDTH, '0');
  const int n_cycle = std::max(1, getNumNode() / dataReadParNumber());
  long total_bytes = 0;
  {
    TIMER_FLOPS("dataReadParNode-fread");
    for (int k = 0; k < dsns.size(); k++) {
      timer.flops += dsns[k].size * dsns[k].nmemb * getNumNode();
    }
    for (int i = 0; i < n_cycle; i++) {
      TIMER_VERBOSE_FLOPS("dataReadParNode-fread-flops");
      long bytes = 0;
      if (getIdNode() % n_cycle == i) {
        FILE* file = fopen(filename.c_str(), "r");
        assert(NULL != file);
        for (int k = 0; k < dsns.size(); k++) {
          bytes += dsns[k].size *
                   fread(dsns[k].data, dsns[k].size, dsns[k].nmemb, file);
        }
        fclose(file);
      }
      sumArray(&bytes, 1);
      total_bytes += bytes;
      qlat::DisplayInfo(cname, fname.c_str(),
                  "cycle / n_cycle = %d / %d ; total_bytes = %ld\n", i + 1,
                  n_cycle, total_bytes);
      timer.flops += bytes;
    }
  }
  dataReadInfo(dsns, path);
  return total_bytes;
}

inline long dataReadParNode(void* data, const long size, const long nmemb,
                            const std::string& path)
{
  DataSizeNmemb dsn((void*)data, size, nmemb);
  std::vector<DataSizeNmemb> dsns;
  dsns.push_back(dsn);
  return dataReadParNode(dsns, path);
}

template <class M>
long gcReadParNode(GridComm<M>& gc, const std::string& path)
{
  TIMER_VERBOSE_FLOPS("gcReadParNode");
  const Geometry& geo = gc.getGeometry();
  assert(geo.isOnlyLocal());
  long total_bytes = dataReadParNode(
      gc.getField(), sizeof(M) * geo.multiplicity, geo.localVolume(), path);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M, class N>
void gcFloatFromDouble(GridComm<N>& gcf, const GridComm<M>& gcd)
{
  // gcf can be uninitialized
  TIMER_FLOPS("gcFloatFromDouble");
  timer.flops += gcd.getAllocSize() / sizeof(double);
  const Geometry& geod = gcd.getGeometry();
  assert(geod.isOnlyLocal());
  assert(0 == sizeof(M) % sizeof(double));
  assert(0 == sizeof(N) % sizeof(float));
  assert(0 == (geod.multiplicity * sizeof(M) / sizeof(double)) %
                  (sizeof(N) / sizeof(float)));
  Geometry geof;
  geof.init(geod, 0,
            geod.multiplicity * sizeof(M) / sizeof(double) /
                (sizeof(N) / sizeof(float)));
  gcf.init(geof);
  assert(gcf.getGeometry() == geof);
  const double* gcdf = (const double*)gcd.getField();
  float* gcff = (float*)gcf.getField();
  for (long i = 0; i < gcf.getAllocSize() / sizeof(float); i++) {
    gcff[i] = gcdf[i];
  }
}

template <class M, class N>
void gcDoubleFromFloat(GridComm<M>& gcd, const GridComm<N>& gcf)
{
  // gcd can be uninitialized
  TIMER_FLOPS("gcDoubleFromFloat");
  timer.flops += gcf.getAllocSize() / sizeof(float);
  const Geometry& geof = gcf.getGeometry();
  assert(0 == sizeof(N) % sizeof(float));
  assert(0 == sizeof(M) % sizeof(double));
  assert(0 == (geof.multiplicity * sizeof(N) / sizeof(float)) %
                  (sizeof(M) / sizeof(double)));
  Geometry geod;
  geod.init(geof, 0,
            geof.multiplicity * sizeof(N) / sizeof(float) /
                (sizeof(M) / sizeof(double)));
  gcd.init(geod);
  assert(gcd.getGeometry() == geod);
  const float* gcff = (const float*)gcf.getField();
  double* gcdf = (double*)gcd.getField();
  for (long i = 0; i < gcd.getAllocSize() / sizeof(double); i++) {
    gcdf[i] = gcff[i];
  }
}

template <class M>
long gcReadParNodeDoubleFromFloat(GridComm<M>& gc, const std::string& path)
{
  GridComm<float> gcf;
  gcFloatFromDouble(gcf, gc);
  long total_bytes = gcReadParNode(gcf, path);
  if (0 != total_bytes) {
    gcDoubleFromFloat(gc, gcf);
  }
  return total_bytes;
}

template <class M>
long gcWriteParNodeFloatFromDouble(const GridComm<M>& gc,
                                   const std::string& path,
                                   const mode_t mode = defaultDirMode())
{
  GridComm<float> gcf;
  gcFloatFromDouble(gcf, gc);
  long total_bytes = gcWriteParNode(gcf, path, mode);
  return total_bytes;
}

CPS_END_NAMESPACE
