#include <sys/sysinfo.h>
#include "utils_Matrix_prod.h"
#include "general_funs.h"
#include "check_fun.h"
////#include "utils_low_rho.h"
//#include "utils_eigensys.h"
//#include "utils_construction.h"

int main(int argc, char* argv[])
{
  using namespace qlat;

  inputpara in;begin_Lat(&argc, &argv, in);

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
  LInt w = in.nx*in.ny*in.nz*in.nt;

  qlat::vector_gpu<Complexq > a ;a.resize( L * m * w);
  qlat::vector_gpu<Complexq > b ;b.resize( L * w * m);
  qlat::vector_gpu<Complexq > c0;c0.resize(L * m * n);
  qlat::vector_gpu<Complexq > c1;c1.resize(L * m * n);


  print0("start memory allocation! \n");
  fflush_MPI();


  c0.set_zero();c1.set_zero();
  random_Ty(a.data(), a.size(), 1, in.seed + 1);
  random_Ty(b.data(), b.size(), 1, in.seed + 2);

  ////cudaDeviceSynchronize();
  //for(int i=0;i< in.debuga;i++){matrix_prod_gpu(a.data(), b.data(), c1.data(), m,n,w,L, Conj, true,false, modeGPU);fflush_MPI();}
  for(int i=0;i< in.debuga;i++){matrix_prod(a.data(), b.data(), c1.data(), m,n,w,L, Conj);fflush_MPI();}
  fflush_MPI();
    
  //for(int i=0;i< in.debuga;i++){matrix_prod_cpu(&a[0],&b[0],&c0[0], m,n,w,L, Conj);fflush_MPI();}
  //print0("END CPU Multi! \n");
  //diff_EigenM(c0,c1, "Matrix prod ");

  ////for(int i=0;i< in.debuga;i++){matrix_prod_cpu(&a[0],&b[0],&c0[0], m,n,w,L);}
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

