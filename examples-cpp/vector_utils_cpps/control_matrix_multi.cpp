#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/sysinfo.h>
#include "utils_Matrix_prod.h"
#include "general_funs.h"
#include "utils_check_fun.h"
////#include "utils_low_rho.h"
//#include "utils_construction.h"

#define ComplexM qlat::ComplexT<double>

int main(int argc, char* argv[])
{
  using namespace qlat;

  inputpara in;begin_Lat(&argc, &argv, in);

  ////int icfg  = in.icfg;
  ////int ionum = in.ionum;

  //////int vini  = 0;
  /////int n_vec = in.nvec;

  omp_set_num_threads(omp_get_max_threads());
  qmessage("===nthreads %8d %8d, max %8d \n",qlat::qacc_num_threads(),omp_get_num_threads(),omp_get_max_threads());
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

  qlat::vector_gpu<ComplexM > a ;a.resize( L * m * w);
  qlat::vector_gpu<ComplexM > b ;b.resize( L * w * n);
  qlat::vector_gpu<ComplexM > c0;c0.resize(L * m * n);
  qlat::vector_gpu<ComplexM > c1;c1.resize(L * m * n);

  qmessage("end memory allocation! \n");
  fflush_MPI();

  c0.set_zero();c1.set_zero();
  random_Ty(a.data(), a.size(), 1, in.seed + 1);
  random_Ty(b.data(), b.size(), 1, in.seed + 2);

  ////qlat_GPU_DeviceSynchronize();
  //for(int i=0;i< in.debuga;i++){matrix_prod_gpu(a.data(), b.data(), c1.data(), m,n,w,L, Conj, true,false, modeGPU);fflush_MPI();}
  for(int i=0;i< in.debuga;i++){
    matrix_prod(a.data(), b.data(), c1.data(), m,n,w,L, Conj, false, modeGPU);fflush_MPI();
  }
  ComplexM res = c1.norm2();
  qmessage("===result %.12e %.12e \n", res.real(), res.imag());

  if(test_pointers >= 1)
  {
    c1.set_zero();
    qlat::vector<ComplexM* > aP;aP.resize(L);ComplexM* A = a.data();
    qlat::vector<ComplexM* > bP;bP.resize(L);ComplexM* B = b.data();
    qlat::vector<ComplexM* > cP;cP.resize(L);ComplexM* C = c1.data();
    for(size_t i=0;i < L;i++){
      aP[i] = &A[i* m*w];
      bP[i] = &B[i* w*n];
      cP[i] = &C[i* m*n];
    }
    for(int i=0;i< in.debuga;i++){
      matrix_prodP(aP.data(), bP.data(), cP.data(), m,n,w, L, Conj, false, modeGPU);fflush_MPI();
    }
    ComplexM ra = c1.norm2();
    qmessage("===result %.12e %.12e, %.12e %.12e \n", ra.real(), ra.imag(), (ra-res).real(), (ra-res).imag());
  }
  if(test_pointers >= 2)
  {
    c1.set_zero();
    qlat::vector<ComplexM* > aP;aP.resize(L);ComplexM* A = a.data();
    qlat::vector<ComplexM* > bP;bP.resize(L);ComplexM* B = b.data();
    qlat::vector<ComplexM* > cP;cP.resize(L);ComplexM* C = c1.data();
    qGPU_for(i, Long(L), true, {
      aP[i] = &A[i* m*w];
      bP[i] = &B[i* w*n];
      cP[i] = &C[i* m*n];
    });
    for(int i=0;i< in.debuga;i++){
      matrix_prodP(aP.data(), bP.data(), cP.data(), m,n,w,L, Conj, false, modeGPU);fflush_MPI();
    }
    ComplexM ra = c1.norm2();
    qmessage("===result %.12e %.12e, %.12e %.12e \n", ra.real(), ra.imag(), (ra-res).real(), (ra-res).imag());
  }

  return end_Lat();
}

