#include <sys/sysinfo.h>
#include "io_gwu.h"
#include "general_funs.h"
////#include "utils_smear_src.h"
#include "check_fun.h"
#include "utils_FFT_redistribute.h"
////#include "utils_low_rho.h"

int main(int argc, char* argv[])
{
  using namespace qlat;

  std::vector<Coordinate> size_node_list;
  size_node_list.push_back(Coordinate(1, 1, 1,  1));
  size_node_list.push_back(Coordinate(1, 1, 1,  2));
  size_node_list.push_back(Coordinate(1, 1, 3,  1));
  size_node_list.push_back(Coordinate(1, 1, 1,  4));
  //size_node_list.push_back(Coordinate(1, 1, 3,  2));
  size_node_list.push_back(Coordinate(1, 2, 3,  1));
  //size_node_list.push_back(Coordinate(1, 1, 1,  8));
  size_node_list.push_back(Coordinate(1, 2, 4,  1));
  size_node_list.push_back(Coordinate(1, 1, 1, 12));
  size_node_list.push_back(Coordinate(1, 1, 1, 16));
  //size_node_list.push_back(Coordinate(1, 1, 1, 24));
  size_node_list.push_back(Coordinate(1, 1, 6,  4));
  size_node_list.push_back(Coordinate(1, 1, 1, 32));
  //size_node_list.push_back(Coordinate(1, 8, 4, 2 ));
  size_node_list.push_back(Coordinate(2, 4, 4, 2 ));
  //size_node_list.push_back(Coordinate(2, 8, 4, 1 ));
  //size_node_list.push_back(Coordinate(1, 1, 4, 16));
  //size_node_list.push_back(Coordinate(1, 1, 1, 64));
  //size_node_list.push_back(Coordinate(1, 1, 1, 48));
  //size_node_list.push_back(Coordinate(1, 1, 1, 96));
  //size_node_list.push_back(Coordinate(1, 1, 1,128));
  size_node_list.push_back(Coordinate(4, 4, 8, 16));
  size_node_list.push_back(Coordinate(4, 8, 8, 16));
  //size_node_list.push_back(Coordinate(1, 2, 2, 16));
  //size_node_list.push_back(Coordinate(1, 1, 2, 16));
  //size_node_list.push_back(Coordinate(1, 2, 2, 16));
  //size_node_list.push_back(Coordinate(2, 2, 2, 16));
  //size_node_list.push_back(Coordinate(2, 2, 4, 16));
  //size_node_list.push_back(Coordinate(2, 4, 4, 16));
  //size_node_list.push_back(Coordinate(4, 4, 4, 16));

  //begin_thread(&argc, &argv, size_node_list);
  begin(&argc, &argv, size_node_list);
  set_GPU();
  ///set_GPU_threads();

  //fft_desc_basic fd();

  //Coordinate size_node = Coordinate(fd.mx, fd.my, fd.mz, fd.mt);
  //begin(fd.rank, size_node);
  //begin(MPI_COMM_WORLD, size_node);

  inputpara in;
  in.load_para("input.txt");

  int nx,ny,nz,nt;
  nx = in.nx;
  ny = in.ny;
  nz = in.nz;
  nt = in.nt;

  int icfg  = in.icfg;
  int ionum = in.ionum;

  ////int vini  = 0;
  int n_vec = in.nvec;

  omp_set_num_threads(omp_get_max_threads());
  print0("===nthreads %8d %8d, max %8d \n",qlat::qacc_num_threads(),omp_get_num_threads(),omp_get_max_threads());

  Coordinate total_site = Coordinate(nx, ny, nz, nt);
  Geometry geo;
  geo.init(total_site, 1); 
  fflush_MPI();

  fft_desc_basic fd(geo);
  //FFT_single fft_large(fd, 1);
  FFT_single fft_large(fd);


  /////Need to understand
  std::vector<int > secT;secT.resize(fd.Nmpi);
  for(int i=0;i<secT.size();i++){secT[i] = fd.Nt;}


  int Nt = fd.Nt;

  int nvol = in.nx*in.ny*in.nz;
  qlat::vector<qlat::Complex > sendbuf,recvbuf;
  qlat::vector<qlat::Complex > databuf;
  int biva = 2*(fd.vol*fd.Nt)/(fd.Nvol);int civ = 6;
  int Nvec = biva*fd.Nvol/(fd.vol*fd.Nt);
  int b0 = Nvec; 
  int c0 = 2*civ;

  sendbuf.resize(b0*civ*nvol*fd.Nt);
  recvbuf.resize(b0*civ*nvol*fd.Nt);
  databuf.resize(b0*civ*nvol*fd.Nt);
  for(long bi=0;bi<biva;bi++)
  for(long long vi=0;vi<fd.Nvol;vi++)
  for(long ci=0;ci<civ;ci++)
  {
    long long off = (bi*fd.Nvol+vi)*civ+ci;
    //////Need a function in fd to get the coordinates?
    int ti = fd.Pos0[fd.rank][3] +  vi/(fd.Nz*fd.Ny*fd.Nx);
    int zi = fd.Pos0[fd.rank][2] + (vi%(fd.Nz*fd.Ny*fd.Nx))/(fd.Ny*fd.Nx);
    int yi = fd.Pos0[fd.rank][1] + (vi%(fd.Ny*fd.Nx))/fd.Nx;
    int xi = fd.Pos0[fd.rank][0] + (vi%fd.Nx);
    databuf[off] = ((ti*900+zi)*900+yi)*900+xi + std::cos(ci) + qlat::Complex(0.0,std::cos(bi));
  }

  ///for(int i=0;i<biva*12*12/civ;i++){memcpy(&sendbuf[i*Nt*N0*N1*N2*civ + 0],&srcE_civ[i](0),sizeof(Ftype)*2*Nt*N0*N1*N2*civ);}
  memcpy(&sendbuf[0], &databuf[0], 2*sizeof(double)*sendbuf.size());

  /////Data will be modified for sendbuf and recvbuf, results on sendbuf
  {
  TIMER("Single FFT");
  fft_large.reorder((double*) &sendbuf[0],(double*) &recvbuf[0],b0,c0 ,  0);
  } 

  double diff = 0.0;
  for(long ni=0;ni<b0;ni++)
  for(int Nti=0;Nti<fd.Nt;Nti++)
  for(long long vi=0;vi<nvol;vi++)
  for(long ci=0;ci<civ;ci++)
  {
    long long off = ((ni*fd.Nt+Nti)*nvol+vi)*civ+ci;
    //////Need a function in fd to get the coordinates?
    int ti = fd.Pos0[fd.rank][3] + Nti;
    /////Order by ti 
    int zi = (vi%(fd.nz*fd.ny*fd.nx))/(fd.ny*fd.nx);
    int yi = (vi%(fd.ny*fd.nx))/fd.nx;
    int xi = (vi%fd.nx);
    ///////Need a function for this offset from coordinate

    //int bi = ni*fd.mz*fd.my*fd.mx + fd.get_mi_curr();
    /////Becarefule about the order
    int bi = fd.get_mi_curr()*b0 + ni;

    qlat::Complex tem = sendbuf[off] - (((ti*900+zi)*900+yi)*900+xi + std::cos(ci) + qlat::Complex(0.0,std::cos(bi)));
    diff += qnorm(tem);
  }
  sum_all_size(&diff, 1);
  print0("Diff rotate 0 %.3e .\n",diff);


  fft_large.reorder((double*) &sendbuf[0],(double*) &recvbuf[0],b0,c0 ,100);

  diff = 0.0;
  for(size_t i = 0;i< sendbuf.size();i++)
  {
    diff += qnorm(sendbuf[i] - databuf[i]);
  }
  sum_all_size(&diff, 1);
  print0("Diff rotate 1 %.3e .\n",diff);


  fflush_MPI();
  qlat::Timer::display();

  qlat::end();
  return 0;
}

