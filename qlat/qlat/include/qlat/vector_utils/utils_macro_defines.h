// utils_macro_defines.h
// Gen Wang
// Sep. 2025

#ifndef UTILS_MACRO_DEFINES_H
#define UTILS_MACRO_DEFINES_H
#pragma once

#include <string.h>
#include <sys/resource.h>
#include <mpi.h>
#include <time.h>
#include <typeinfo>
#include <iterator>
#include <cstdarg>
#include <qlat-utils/mat-vec.h>
#include <qlat-utils/eigen.h>
#include <qlat/qcd.h>

namespace qlat
{

/////__CLEAR_SMEAR_MEM__ FLAG for smear memory clear

#ifdef __ENABLE_DOUBLE__
#define Enablefloat 0
#else
#define Enablefloat 1
#endif

#define MAX_VECTOR_GPU_BUF 100
#define QLAT_FILE_IO_SIZE  30

// AMD machine MPI have tag issues....
#define QLAT_VECTOR_UTILS_MPI_TAG 8712

////q_io_vec_ionum
////q_file_io_each_size

//#define __NO_GPU_DIRECT__
//#ifdef  __HIP_PLATFORM_HCC__
//#define __NO_GPU_DIRECT__
//#endif

#ifdef __DEBUG_VECUTILS__
#define PRINT_TIMER 10
#else
#define PRINT_TIMER  0
#endif

#define LInt unsigned long


#define large_vuse Elarge_vector
#if Enablefloat==0
#define Complexq qlat::ComplexT<double >
#define Ftype double
//////needed for contraction change to small power of 2 if shared memory too small
#endif

#if Enablefloat==1
#define Complexq qlat::ComplexT<float >
#define Ftype float
#endif

#ifdef __QLAT_BARYON_SHARED_SMALL__
#define BFACG_SHARED 1
#else
#define BFACG_SHARED 4
#endif

#define Evector qlat::vector<Complexq >
#define EigenV   qlat::vector<Complexq >

/////May have some errors for python, mem buf, mem leak
////#define EigenM qlat::vector<Evector >
#define EigenM  std::vector<Evector >
////

#define qnoi qlat::FieldM<Complexq, 1>
#define qnoiT qlat::FieldM<Ty, 1>

//////dim 12*12 --> Nt --> Nxyz
/////only memory size and geo are used, donot use others
#define qprop   qlat::FieldM<Complexq, 12*12>
#define qpropT  qlat::FieldM<Ty, 12*12>

/////Fields for staggered fermion
#define colorFD qlat::FieldM<ComplexD , 3>
#define colorFF qlat::FieldM<Complexq, 3>
#define colorFT qlat::FieldM<Ty, 3>

/////read input length of each line
#define LINE_LIMIT 3000

//#define QLAT_PI_LOCAL 3.1415926535898
#define QLAT_PI_LOCAL 3.14159265358979323846

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

//// *************** FOR ERROR CHECKING *******************
#ifdef QLAT_USE_ACC
#define qacc_Errchk(ans) { qacc_Err((ans), __FILE__, __LINE__); }
#endif

#define qGPU_for(iter, num, GPU, ...) \
  Qassert(int(GPU) != -2); \
  if(bool(GPU) == true){qacc_for(iter, num, {__VA_ARGS__});} \
  else{qthread_for(iter, num, {__VA_ARGS__});}

#define qGPU_forNB(iter, num, GPU, ...) \
  Qassert(int(GPU) != -2); \
  if(bool(GPU) == true){qacc_forNB(iter, num, {__VA_ARGS__});} \
  else{qthread_for(iter, num, {__VA_ARGS__});}

/*
*
//#define qGPU_for2d(iter1, num1, iter2, num2, GPU, ...) \
//  Qassert(int(GPU) != -2); \
//  if(bool(GPU) == true){qacc_for2d(iter1, num1, iter2, num2, {__VA_ARGS__});} \
//  else{qthread_for2d(iter1, num1, iter2, num2, {__VA_ARGS__});}
//
//#define qGPU_for2dNB(iter1, num1, iter2, num2, GPU, ...) \
//  Qassert(int(GPU) != -2); \
//  if(bool(GPU) == true){qacc_for2dNB(iter1, num1, iter2, num2, {__VA_ARGS__});} \
//  else{qthread_for2d(iter1, num1, iter2, num2, {__VA_ARGS__});}
*
*/

/*
  checkpoint with environment settings
  display_mem_type();maybe used for debug
*/
#define Qcheck(num)                                                   \
  {                                                                   \
    static Int do_check =                                             \
      qlat::get_env_long_default(std::string("qlat_checkpoint"), -1); \
    if(num == do_check){Qassert(false);}                              \
  }

}

#endif
