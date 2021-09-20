#include <sys/sysinfo.h>
#include "utils_Matrix_prod.h"
#include "general_funs.h"
#include "check_fun.h"
//#include "io_vec.h"
////#include "utils_low_rho.h"
//#include "utils_eigensys.h"
//#include "utils_construction.h"

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

  //begin_thread(&argc, &argv, size_node_list);
  begin(&argc, &argv, size_node_list);
  set_GPU();

  inputpara in;
  in.load_para("input.txt");

  ////int icfg  = in.icfg;
  ////int ionum = in.ionum;

  //////int vini  = 0;
  /////int n_vec = in.nvec;

  omp_set_num_threads(omp_get_max_threads());
  print0("===nthreads %8d %8d, max %8d \n",qlat::qacc_num_threads(),omp_get_num_threads(),omp_get_max_threads());
  fflush_MPI();

  int ncpy = 10;
  int nsrc = 12;
  int nvec = 10;
  int modeGPU = 0;
  int IConj = 0;
  bool Conj = false;

  if(in.paraI != "None"){
    std::vector<std::string > Li = stringtolist(in.paraI);
    print0("Li %s, size %d \n", in.paraI.c_str(),int(Li.size()) );
    fflush_MPI();
    ncpy    = stringtonum(Li[0]);
    nsrc    = stringtonum(Li[1]);
    nvec    = stringtonum(Li[2]);
    modeGPU = stringtonum(Li[3]);
    IConj   = stringtonum(Li[4]);
  }
  if(IConj == 1){Conj = true;}

  LInt L = ncpy;
  LInt m = nsrc;
  LInt n = nvec;
  LInt w = in.nx*in.ny*in.nz*in.nt*6 + nvec;


  EigenV a( L*m*w);
  EigenV b( L*w*n);
  EigenV c0(L*m*n);
  EigenV c1(L*m*n);

  print0("start memory allocation! \n");
  fflush_MPI();

  //c0.resize(L*m*n);
  //cudaDeviceSynchronize();
  //print0("END 0 memory allocation! \n");
  //fflush_MPI();

  //cudaDeviceSynchronize();
  //print0("END 1 memory allocation! \n");
  //fflush_MPI();

  //c1.resize(L*m*n);
  //cudaDeviceSynchronize();
  //cudaDeviceSynchronize();

  //a.resize( L*m*w);
  //cudaDeviceSynchronize();
  //cudaDeviceSynchronize();
  //b.resize( L*w*n);
  //cudaDeviceSynchronize();

  //cudaDeviceSynchronize();
  //print0("END 0 memory allocation! \n");fflush_MPI();

  zeroE(c0,1);zeroE(c1,1);
  ran_EigenM( a,1);ran_EigenM( b,1);

  //cudaDeviceSynchronize();
  //print0("END 1 memory allocation! \n");fflush_MPI();

  ////for(int i=0;i< in.debuga;i++){matrix_prod_cpu(&a[0],&b[0],&c0[0], m,n,w,L);}
  for(int i=0;i< in.debuga;i++){matrix_prod_cpu(&a[0],&b[0],&c0[0], m,n,w,L, Conj);fflush_MPI();}
  print0("END CPU Multi! \n");
  fflush_MPI();

  ////cudaDeviceSynchronize();
  for(int i=0;i< in.debuga;i++){matrix_prod_gpu(&a[0],&b[0],&c1[0], m,n,w,L, Conj, true,false, modeGPU);fflush_MPI();}
    

  diff_EigenM(c0,c1, "Matrix prod ");

  ////matrix_prod_cpu(&a[0],&b[0],&c0[0], m,n,w,L);
  ////for(int i=0;i< in.debuga;i++){
  ////  zeroE(c1);matrix_prod_gpu(&a[0],&b[0],&c1[0], m,n,w,L, modeGPU);
  ////}


  //for(int i=0;i< in.debuga;i++){matrix_prod_cpu(&a[0],&b[0],&c1[0], m,n,w,L);}
  //print0("END CPU Multi! \n");
  //fflush_MPI();



  fflush_MPI();
  qlat::Timer::display();

  qlat::end();
  return 0;
}

