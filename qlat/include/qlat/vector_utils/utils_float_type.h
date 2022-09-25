// utils_float_type.h
// Gen Wang
// Jun. 2021

#ifndef UTILS_FLOAT_TYPE_H
#define UTILS_FLOAT_TYPE_H
#pragma once

#include <string.h>
#include <sys/resource.h>
#include <mpi.h>
#include <time.h>
#include <typeinfo>
#include <iterator>
#include <sys/sysinfo.h>

#include <qlat/qcd.h>

namespace qlat{

/////__CLEAR_SMEAR_MEM__ FLAG for smear memory clear

#ifdef __ENABLE_DOUBLE__
#define Enablefloat 0
#else
#define Enablefloat 1
#endif

//#define __NO_GPU_DIRECT__
//#ifdef  __HIP_PLATFORM_HCC__
//#define __NO_GPU_DIRECT__
//#endif

#ifdef __DEBUG_VECUTILS__
#define PRINT_TIMER 10
#else
#define PRINT_TIMER  0
#endif

////#include <cuda_runtime.h>
//
////#include <Eigen/Dense>
////#include "fftw3.h"
////#include "fftw3-mpi.h"

#define LInt unsigned long

#define large_vuse Elarge_vector
#if Enablefloat==0
#define Complexq qlat::Complex
#define Ftype double
//////needed for contraction change to small power of 2 if shared memory too small
#define BFACG_SHARED 4
#endif

#if Enablefloat==1
#define Complexq qlat::ComplexF
#define Ftype float
#define BFACG_SHARED 8
#endif

#define Evector qlat::vector_acc<Complexq >
#define EigenV   qlat::vector_acc<Complexq >

/////May have some errors for python, mem buf, mem leak
////#define EigenM qlat::vector<Evector >
#define EigenM  std::vector<Evector >
////

#define qnoi qlat::FieldM<Complexq, 1>
//////dim 12*12 --> Nt --> Nxyz
#define qprop  qlat::FieldM<Complexq, 12*12>

#define qnoiT qlat::FieldM<Ty, 1>
//////dim 12*12 --> Nt --> Nxyz
#define qpropT  qlat::FieldM<Ty, 12*12>

/////Fields for staggered fermion
#define colorFD qlat::FieldM<Complex , 3>
#define colorFF qlat::FieldM<Complexq, 3>
#define colorFT qlat::FieldM<Ty, 3>

/////read input length of each line
#define LINE_LIMIT 3000

#define PI 3.1415926535898

inline void print_NONE(const char *filename)
{
  (void)filename;
  return;
}

#define TIMERA(name) print_NONE(name);
#define TIMERB(name) print_NONE(name);
#define TIMERC(name) print_NONE(name);
#define TIMERD(name) print_NONE(name);
#define TIMERE(name) print_NONE(name);
#define TIMERZ(name) TIMER(name);

#if PRINT_TIMER>0
#undef TIMERT
#define TIMERT(name) TIMER((std::string("T_T ") + std::string(name)).c_str());
#undef TIMERZ
#define TIMERZ(name) TIMER(name);
#endif

#if PRINT_TIMER>1
#undef TIMERE
#define TIMERE(name) TIMER((std::string("E_T ") + std::string(name)).c_str());
#endif

#if PRINT_TIMER>2
#undef TIMERD
#define TIMERD(name) TIMER((std::string("D_T ") + std::string(name)).c_str());
#endif

#if PRINT_TIMER>3
#undef  TIMERC
#define TIMERC(name) TIMER((std::string("C_T ") + std::string(name)).c_str());
#endif

#if PRINT_TIMER>4
#undef  TIMERB
#define TIMERB(name) TIMER((std::string("B_T ") + std::string(name)).c_str());
#endif

#if PRINT_TIMER>5
#undef  TIMERA
#define TIMERA(name) TIMER((std::string("A_T ") + std::string(name)).c_str());
#endif


#define ckpoint abort_r("check point \n");

#ifdef QLAT_USE_ACC
// *************** FOR ERROR CHECKING *******************
#ifndef CUDA_RT_CALL
#define CUDA_RT_CALL( call )                                                                                           \
    {                                                                                                                  \
        auto status = static_cast<cudaError_t>( call );                                                                \
        if ( status != cudaSuccess )                                                                                   \
            fprintf( stderr,                                                                                           \
                     "ERROR: CUDA RT call \"%s\" in line %d of file %s failed "                                        \
                     "with "                                                                                           \
                     "%s (%d).\n",                                                                                     \
                     #call,                                                                                            \
                     __LINE__,                                                                                         \
                     __FILE__,                                                                                         \
                     cudaGetErrorString( status ),                                                                     \
                     status );                                                                                         \
    }
#endif  // CUDA_RT_CALL

#ifndef CUFFT_CALL
#define CUFFT_CALL( call )                                                                                             \
    {                                                                                                                  \
        auto status = static_cast<cufftResult>( call );                                                                \
        if ( status != CUFFT_SUCCESS )                                                                                 \
            fprintf( stderr,                                                                                           \
                     "ERROR: CUFFT call \"%s\" in line %d of file %s failed "                                          \
                     "with "                                                                                           \
                     "code (%d).\n",                                                                                   \
                     #call,                                                                                            \
                     __LINE__,                                                                                         \
                     __FILE__,                                                                                         \
                     status );                                                                                         \
    }
#endif  // CUFFT_CALL
// *************** FOR ERROR CHECKING *******************
#endif

#ifdef QLAT_USE_ACC
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      qlat::displayln(qlat::ssprintf("cuda error: %s %s %d\n", cudaGetErrorString(code), file, line ));
      qassert(false);
   }
}
#endif

inline void gpuFree(void* res)
{
  if(res!=NULL){
    #ifdef QLAT_USE_ACC
    gpuErrchk(cudaFree(res));
    #else
    //delete [] res;
    free(res);
    #endif
    //res = NULL;
  }
}

inline void* aligned_alloc_no_acc(const size_t size)
{
  const size_t alignment = get_alignment();
#if defined NO_ALIGNED_ALLOC
  return malloc(size);
#else
  return aligned_alloc(alignment, size);
#endif
}


#ifdef QLAT_USE_ACC
#define gpuMalloc(bres,bsize, Ty) {gpuErrchk(cudaMalloc(&bres, bsize*sizeof(Ty)));}
#else
#define gpuMalloc(bres,bsize, Ty) {bres = (Ty *)aligned_alloc_no_acc(bsize*sizeof(Ty));}
#endif

inline void free_buf(void* buf, bool GPU){
  if(buf != NULL){if(GPU){gpuFree(buf);}else{free(buf);}}
  buf = NULL;
}

//template<typename Ty>
//inline void alloc_buf(Ty* buf, size_t size, bool GPU){
//  free_buf(buf, GPU);
//  buf = (Ty*) malloc(size * sizeof(Ty));
//  //if(GPU){gpuMalloc(buf, size, Ty);}
//  //else{buf = (Ty*) malloc(size * sizeof(Ty));}
//}


#ifndef QLAT_USE_ACC
#ifdef __GNUC__
///////Multiply of different types of complex

//inline std::complex<double> operator*(const std::complex<float> &a, const double &b) {
//    return std::complex<double>(a.real()*b, a.imag()*b);
//}
//
//inline std::complex<double> operator/(const std::complex<float> &a, const double &b) {
//    return std::complex<double>(a.real()/b, a.imag()/b);
//}

inline std::complex<double> operator*(const std::complex<float> &a, const std::complex<double > &b) {
    return std::complex<double>(a.real() * b.real() - a.imag()*b.imag(), a.imag()*b.real() + a.real()*b.imag());
}

inline std::complex<double> operator*(const std::complex<double > &b, const std::complex<float> &a) {
    return std::complex<double>(a.real() * b.real() - a.imag()*b.imag(), a.imag()*b.real() + a.real()*b.imag());
}

inline std::complex<double> operator-(const std::complex<double > &a, const std::complex<float> &b) {
    return std::complex<double>(a.real() - b.real() , a.imag() - b.imag());
}
inline std::complex<double> operator-(const std::complex<float > &a, const std::complex<double> &b) {
    return std::complex<double>(a.real() - b.real() , a.imag() - b.imag());
}
#endif
#endif

template<typename Ty>
void zero_Ty(Ty* a, long size,int GPU=0, bool dummy=true)
{
  (void)GPU;
  (void)dummy;
  #ifdef QLAT_USE_ACC
  if(GPU == 1){
    cudaMemsetAsync(a, 0, size*sizeof(Ty));
    if(dummy)qacc_barrier(dummy);
    return ;
  }
  #endif

  #pragma omp parallel for
  for(long isp=0;isp<size;isp++){  a[isp] = 0;}
}

#define print0 if(qlat::get_id_node() == 0) printf

inline unsigned int get_node_rank_funs0()
{
  int rank;
  //MPI_Comm_rank(get_comm(), &rank);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
}

inline void abort_r(std::string stmp=std::string(""))
{
  if(stmp!=std::string(""))print0("%s\n",stmp.c_str());
  //MPI_Barrier(get_comm());
  MPI_Barrier(MPI_COMM_WORLD);
  fflush(stdout);
  ////MPI_Finalize();
  qlat::end();
  //qassert(false);
  abort();
}



}

#endif
