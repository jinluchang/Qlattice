#include <qlat/qcd.h>
#include <sys/sysinfo.h>
#include "utils_test_unified.h"

int main(int argc, char* argv[])
{
  using namespace qlat;

  std::vector<Coordinate> size_node_list;
  size_node_list.push_back(Coordinate(1, 1, 1,  1));
  size_node_list.push_back(Coordinate(1, 1, 1,  2));
  size_node_list.push_back(Coordinate(1, 1, 1,  4));
  size_node_list.push_back(Coordinate(1, 1, 1,  8));
  size_node_list.push_back(Coordinate(1, 1, 1, 16));
  size_node_list.push_back(Coordinate(1, 1, 1, 32));
  size_node_list.push_back(Coordinate(1, 1, 1, 64));
  size_node_list.push_back(Coordinate(4, 4, 8, 16));
  size_node_list.push_back(Coordinate(4, 8, 8, 16));
  begin(&argc, &argv, size_node_list);

  #ifdef QLAT_USE_ACC
  int num_node;MPI_Comm_size(get_comm(), &num_node);
  int id_node;MPI_Comm_rank(get_comm(), &id_node);

  int num_gpus = 0;
  cudaGetDeviceCount(&num_gpus);
  ////cudaDeviceReset();
  cudaSetDevice(id_node % num_gpus);
  int gpu_id = -1;
  cudaGetDevice(&gpu_id);
  printf("CPU node %d (of %d) uses CUDA device %d\n", id_node, num_node, gpu_id);
  fflush(stdout);
  MPI_Barrier(get_comm());
  #endif

  omp_set_num_threads(omp_get_max_threads());

  //qlat::displayln_info(qlat::ssprintf(
  //    "===nthreads %8d %8d, max %8d \n",qlat::qacc_num_threads(),omp_get_num_threads(),omp_get_max_threads()
  //  ));

  //Coordinate total_site = Coordinate(nx, ny, nz, nt);
  //Geometry geo;
  //geo.init(total_site, 1); 

  //test_unified ei(1);
  //////ei.initiallize_mass(20);
  //ei.prop_L();

  qlat::vector<qlat::Complex > a;
  qlat::vector<qlat::Complex > b;
  qlat::vector<qlat::Complex > c;
  int m = 1024;
  int n = 1024;
  int w = 1024*64;
  a.resize(m*w);b.resize(n*w);c.resize(m*n);
  qacc_for(coff, long(a.size()),{a[coff] = std::cos(coff*0.7);});
  qacc_for(coff, long(b.size()),{b[coff] = qlat::Complex(std::cos(coff*1.7), 0.5/(coff+1));});

  {TIMER("Loop A");
  qacc_for(coff, long(m*n),{
    long mi = coff/n;
    long ni = coff%n;
    qlat::Complex buf = 0.0;
    for(int wi=0;wi<w;wi++){buf += a[mi*w + wi] * b[ni*w + wi];}
    c[coff] += buf;
  });}

  qlat::displayln_info(qlat::ssprintf("a0 qacc for call"));


  m = 2024;
  n = 2023;
  w = 2024*12;

  a.resize(m*w);b.resize(n*w);c.resize(m*n);
  qacc_for(coff, long(a.size()),{a[coff] = std::cos(coff*0.7);});
  qacc_for(coff, long(b.size()),{b[coff] = qlat::Complex(std::cos(coff*1.7), 0.5/(coff+1));});

  {TIMER("Loop B");
  qacc_for(coff, long(m*n),{
    long mi = coff/n;
    long ni = coff%n;
    qlat::Complex buf = 0.0;
    for(int wi=0;wi<w;wi++){buf += a[mi*w + wi] * b[ni*w + wi];}
    c[coff] += buf;
  });}

  qlat::displayln_info(qlat::ssprintf("a1 qacc for call"));


  /////alpha.resize(0);
  /////alpha.resize(200);
  /////qlat::displayln_info(qlat::ssprintf("a1 qacc for call"));



  qlat::Timer::display();

  qlat::end();
  return 0;
}

