#include "utils_io_vec.h"
#include "utils_reduce_vec.h"

////#define Complexq qlat::ComplexF

int main(int argc, char* argv[])
{
  using namespace qlat;

  /////namespace qcd
  //init_machine_thread(argc, argv,false);
  //timer walltime;walltime.start("over all");

  std::vector<Coordinate> size_node_list;
  size_node_list.push_back(Coordinate(1, 1, 1, 1));
  size_node_list.push_back(Coordinate(1, 1, 1, 2));
  size_node_list.push_back(Coordinate(1, 1, 1, 4));
  size_node_list.push_back(Coordinate(1, 1, 1, 8));
  size_node_list.push_back(Coordinate(1, 1, 1, 16));
  size_node_list.push_back(Coordinate(1, 1, 2, 16));
  size_node_list.push_back(Coordinate(1, 2, 2, 16));
  size_node_list.push_back(Coordinate(2, 2, 2, 16));
  size_node_list.push_back(Coordinate(2, 2, 4, 16));
  size_node_list.push_back(Coordinate(2, 4, 4, 16));
  size_node_list.push_back(Coordinate(4, 4, 4, 16));
  begin(&argc, &argv, size_node_list);

  //fft_desc_basic fd();

  //Coordinate size_node = Coordinate(fd.mx, fd.my, fd.mz, fd.mt);
  //begin(fd.rank, size_node);
  //begin(MPI_COMM_WORLD, size_node);

  int nx,ny,nz,nt;
  nx = 24;
  ny = 24;
  nz = 24;
  nt = 64;

  //int n_vec = 1;

  Coordinate total_site = Coordinate(nx, ny, nz, nt);
  Geometry geo;
  geo.init(total_site, 1); 

  std::vector<int > nv(4);
  std::vector<int > Nv(4);
  ///int nx,ny,nz,nt;
  int Nx,Ny,Nz,Nt;

  for(int i=0;i<4;i++){Nv[i]=geo.node_site[i];nv[i] = geo.node_site[i] * geo.geon.size_node[i];}
  nx = nv[0];ny = nv[1];nz = nv[2];nt = nv[3];
  Nx = Nv[0];Ny = Nv[1];Nz = Nv[2];Nt = Nv[3];

  //int vol   =  nx*ny*nz*nt;
  size_t noden =  Nx*Ny*Nz*Nt;

  qlat::vector<Complexq > p;p.resize(12*noden);
  for(size_t isp=0;isp<noden*12;isp++)
  {
    p[isp] = Complexq(std::cos((isp+3)*1.0),0.0);
  }

  double res0 = 0.0;
  double res1 = 0.0;
  for(size_t isp=0;isp<noden*12;isp++)
  {
    res0 += (p[isp].real()*p[isp].real() + p[isp].imag()*p[isp].imag());
  }

  #pragma omp parallel for reduction(+: res1)
  for(size_t isp=0;isp<noden*12;isp++)
  {
    res1 += (p[isp].real()*p[isp].real() + p[isp].imag()*p[isp].imag());
  }

  print0("res0 %.6e, res1 %.6e ,diff %.6e \n",res0,res1,(res0-res1)/res0);



  fflush_MPI();
  qlat::Timer::display();



  qlat::end();
  return 0;
}

