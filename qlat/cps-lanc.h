#pragma once

CPS_START_NAMESPACE

long lanczosWriteParNode(const Lanczos& lanc, const std::string& path)
{
  TIMER_VERBOSE_FLOPS("lanczosWriteParNode");
  timer.flops += lanc.size * lanc.vec_size * getNumNode();
  Makedir(path);
  if (0 == getIdNode()) {
    const std::string filename = path + "/eigen-values.txt";
    FILE* file = fopen(filename.c_str(), "w");
    fprintf(file, "%d\n", lanc.size);
    for (int i = 0; i < lanc.size; i++) {
      double val = lanc.getVal(i);
      fprintf(file, "%.20lE\n", val);
    }
    fclose(file);
  }
  std::vector<DataSizeNmemb> dsns;
  dsns.resize(lanc.size);
  for (int k = 0; k < lanc.size; k++) {
    dsns[k].data = (void*)lanc.getVec(k);
    dsns[k].size = lanc.vec_size;
    dsns[k].nmemb = 1;
  }
  qlat::DisplayInfo(cname, fname.c_str(), "Writing %d vectors.\n", lanc.size);
  long total_bytes = dataWriteParNode(dsns, path);
  saveString(path + "/checkpoint", "");
  return total_bytes;
}

long lanczosReadParNode(Lanczos& lanc, const std::string& path)
{
  if (!DoesFileExist(path + "/checkpoint")) {
    qlat::DisplayInfo(cname, "lanczosReadParNode", "'%s' do not exist.\n",
                path.c_str());
    return 0;
  }
  TIMER_VERBOSE_FLOPS("lanczosReadParNode");
  long nvec = 0;
  if (0 == getIdNode()) {
    const std::string filename = path + "/eigen-values.txt";
    FILE* file = fopen(filename.c_str(), "r");
    fscanf(file, "%ld\n", &nvec);
    fclose(file);
  }
  sumArray(&nvec, 1);
  lanc.free_evecs();
  lanc.resize(nvec);
  timer.flops += lanc.size * lanc.vec_size * getNumNode();
  double vals[lanc.size];
  memset(vals, 0, sizeof(vals));
  qlat::DisplayInfo(cname, fname.c_str(), "Reading %d eigen-values.\n", lanc.size);
  if (0 == getIdNode()) {
    const std::string filename = path + "/eigen-values.txt";
    FILE* file = fopen(filename.c_str(), "r");
    fscanf(file, "%ld\n", &nvec);
    for (int i = 0; i < lanc.size; i++) {
      fscanf(file, "%lE\n", &vals[i]);
      std::cout << sqrt(vals[i]) << std::endl;
    }
    fclose(file);
  }
  sumArray(vals, lanc.size);
  for (int k = 0; k < lanc.size; k++) {
    lanc.getVal(k) = vals[k];
  }
  std::vector<DataSizeNmemb> dsns;
  dsns.resize(lanc.size);
  {
    TIMER_VERBOSE_FLOPS("lanczosReadParNode-alloc");
    timer.flops += lanc.size * lanc.vec_size;
    for (int k = 0; k < lanc.size; k++) {
      lanc.alloc(k);
    }
  }
  for (int k = 0; k < lanc.size; k++) {
    dsns[k].data = lanc.getVec(k);
    dsns[k].size = lanc.vec_size;
    dsns[k].nmemb = 1;
  }
  qlat::DisplayInfo(cname, fname.c_str(), "Reading %d vectors.\n", lanc.size);
  long total_bytes = dataReadParNode(dsns, path);
  return total_bytes;
}

void lanczosFill(Lanczos& lanc, const int size)
{
  TIMER_VERBOSE_FLOPS("lanczosFill");
  timer.flops += size * lanc.vec_size;
  lanc.resize(size);
  for (int k = 0; k < size; k++) {
    lanc.alloc(k);
    lanc.getVal(k) = k;
    float* vec = lanc.getVec(k);
    for (int i = 0; i < lanc.vec_size / sizeof(float); i++) {
      vec[i] = i;
    }
  }
}

CPS_END_NAMESPACE
