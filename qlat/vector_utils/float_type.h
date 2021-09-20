// float_type.h
// Gen Wang
// Jun. 2021

#ifndef FLOAT_TYPE_H
#define FLOAT_TYPE_H
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

#define Enablefloat 0

#define PRINT_TIMER 10
//#define PRINT_TIMER 100

#define DEBUGM 0

#define __NO_MEMCACHE_LOG__

#if DEBUGM==1
#undef __NO_MEMCACHE_LOG__
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
#endif

#if Enablefloat==1
#define Complexq qlat::ComplexF
#define Ftype float
#endif

////#define Evector qlat::vector<std::complex< Ftype > >
////#define EA Eigen::Map<Eigen::Array<std::complex<Ftype >,Eigen::Dynamic,1 > >
////#define EM Eigen::Map< Eigen::Matrix< std::complex<Ftype >, Eigen::Dynamic, Eigen::Dynamic ,Eigen::RowMajor> > 

#define Evector qlat::vector_acc<Complexq >
#define EigenV  qlat::vector_acc<Complexq >
#define EA Eigen::Map<Eigen::Array<Complexq ,Eigen::Dynamic,1 > >
#define EM Eigen::Map< Eigen::Matrix<Complexq , Eigen::Dynamic, Eigen::Dynamic ,Eigen::RowMajor> > 
#define EMC Eigen::Map< Eigen::Matrix<Complexq , Eigen::Dynamic, Eigen::Dynamic ,Eigen::ColMajor> > 
////#define EM Eigen::Map< Eigen::Matrix<std::complex<Ftype > , Eigen::Dynamic, Eigen::Dynamic ,Eigen::RowMajor> > 
////#define EMC Eigen::Map< Eigen::Matrix<std::complex<Ftype > , Eigen::Dynamic, Eigen::Dynamic ,Eigen::ColMajor> > 

/////May have some errors for python, mem buf, mem leak
////#define EigenM qlat::vector<Evector >
#define EigenM std::vector<Evector >

#define EG  Eigen::Matrix<Complexq , Eigen::Dynamic, Eigen::Dynamic ,Eigen::RowMajor>  
#define EGC Eigen::Matrix<Complexq , Eigen::Dynamic, Eigen::Dynamic ,Eigen::ColMajor> 


void print_NONE(const char *filename){return ;}

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
    res = NULL;
    #endif
  }
}

#ifdef QLAT_USE_ACC
#define gpuMalloc(bres,bsize, Ty) {gpuErrchk(cudaMalloc(&bres, bsize*sizeof(Ty)));}
#else
#define gpuMalloc(bres,bsize, Ty) {bres = (Ty *)malloc(bsize*sizeof(Ty));}
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
#endif
#endif


}

#endif
