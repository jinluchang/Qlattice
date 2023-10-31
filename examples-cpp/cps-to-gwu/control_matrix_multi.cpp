#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/sysinfo.h>
#include "utils_Matrix_prod.h"
#include "general_funs.h"
#include "check_fun.h"
////#include "utils_low_rho.h"
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
  int test_pointers = 0;
  in.find_para(std::string("ncpy" ), ncpy);
  in.find_para(std::string("nsrc" ), nsrc);
  in.find_para(std::string("nvec" ), nvec);
  in.find_para(std::string("modeGPU" ), modeGPU);
  in.find_para(std::string("IConj" ), IConj);
  in.find_para(std::string("test_pointers" ), test_pointers);

  bool Conj = false;
  if(IConj == 1){Conj = true;}

  LInt L = ncpy;
  LInt m = nsrc;
  LInt n = nvec;
  LInt w = in.nx*in.ny*in.nz*in.nt;

  qlat::vector_gpu<Complexq > a ;a.resize( L * m * w);
  qlat::vector_gpu<Complexq > b ;b.resize( L * w * n);
  qlat::vector_gpu<Complexq > c0;c0.resize(L * m * n);
  qlat::vector_gpu<Complexq > c1;c1.resize(L * m * n);

  print0("end memory allocation! \n");
  fflush_MPI();

  c0.set_zero();c1.set_zero();
  random_Ty(a.data(), a.size(), 1, in.seed + 1);
  random_Ty(b.data(), b.size(), 1, in.seed + 2);

  ////cudaDeviceSynchronize();
  //for(int i=0;i< in.debuga;i++){matrix_prod_gpu(a.data(), b.data(), c1.data(), m,n,w,L, Conj, true,false, modeGPU);fflush_MPI();}
  for(int i=0;i< in.debuga;i++){
    matrix_prod(a.data(), b.data(), c1.data(), m,n,w,L, Conj, false, modeGPU);fflush_MPI();
  }
  Complexq res = c1.norm2();
  print0("===result %.8e %.8e \n", res.real(), res.imag());

  if(test_pointers >= 1)
  {
    c1.set_zero();
    qlat::vector_acc<Complexq* > aP;aP.resize(L);Complexq* A = a.data();
    qlat::vector_acc<Complexq* > bP;bP.resize(L);Complexq* B = b.data();
    qlat::vector_acc<Complexq* > cP;cP.resize(L);Complexq* C = c1.data();
    for(size_t i=0;i < L;i++){
      aP[i] = &A[i* m*w];
      bP[i] = &B[i* w*n];
      cP[i] = &C[i* m*n];
    }
    for(int i=0;i< in.debuga;i++){
      matrix_prodP(aP.data(), bP.data(), cP.data(), m,n,w, L, Conj, false, modeGPU);fflush_MPI();
    }
    Complexq ra = c1.norm2();
    print0("===result %.8e %.8e, %.8e %.8e \n", ra.real(), ra.imag(), (ra-res).real(), (ra-res).imag());
  }
  if(test_pointers >= 2)
  {
    c1.set_zero();
    qlat::vector_acc<Complexq* > aP;aP.resize(L);Complexq* A = a.data();
    qlat::vector_acc<Complexq* > bP;bP.resize(L);Complexq* B = b.data();
    qlat::vector_acc<Complexq* > cP;cP.resize(L);Complexq* C = c1.data();
    qGPU_for(i, long(L), true, {
      aP[i] = &A[i* m*w];
      bP[i] = &B[i* w*n];
      cP[i] = &C[i* m*n];
    });
    for(int i=0;i< in.debuga;i++){
      matrix_prodP(aP.data(), bP.data(), cP.data(), m,n,w,L, Conj, false, modeGPU);fflush_MPI();
    }
    Complexq ra = c1.norm2();
    print0("===result %.8e %.8e, %.8e %.8e \n", ra.real(), ra.imag(), (ra-res).real(), (ra-res).imag());
  }

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

  return end_Lat();
}

