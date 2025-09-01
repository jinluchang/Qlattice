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
#include <qlat-utils/mat-vec.h>
#include <qlat-utils/eigen.h>
#include <qlat/qcd.h>
#include <cstdarg>

#ifndef QLAT_NO_SYSINFO
#include <sys/sysinfo.h>
#endif


namespace qlat{

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

//enum QBOOL{
//  QFALSE = 0,
//  QTRUE  = 1
//};

enum QBOOL{
  QFALSE,
  QTRUE,
};

enum QMEM{
  QMSYNC = -1,
  QMCPU  =  0,
  QMGPU  =  1
};

enum QMEM_ORDER{
  QLAT_DEFAULT,// default t,z,y,x : multiplicity
  QLAT_OUTTER  //  multiplicity : t, z, y, x
};

inline QMEM get_type_mem(QMEM g){
  return g;
}

inline QMEM get_type_mem(int g){
  if(g == -1){return QMSYNC;}
  if(g ==  0){return QMCPU ;}
  if(g ==  1){return QMGPU ;}
  return QMCPU;
}

inline int check_GPU_same(QMEM a, QMEM b)
{
  if(a == QMSYNC){return 1;}
  if(b == QMSYNC){return 1;}
  if(a == b){return 1;}
  return 0;
}

inline int check_GPU_same(QMEM a, bool b)
{
  if(a == QMSYNC){return 1;}
  if(a == QMGPU and b == 1){return 1;}
  if(a == QMCPU and b == 0){return 1;}
  return 0;
}

inline int check_GPU_multi(QMEM a, QMEM b)
{
  if(a == QMSYNC){return b;}
  if(b == QMSYNC){return a;}
  if(a == b){return a;}
  return -2;
}

inline int check_GPU_multi(QMEM a, bool b)
{
  if(a == QMSYNC){return b;}
  if(a == QMGPU and b == 1){return 1;}
  if(a == QMCPU and b == 0){return 0;}
  return -2;
}

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

inline void gpuFree(void* res)
{
  if(res!=NULL){
    #ifdef QLAT_USE_ACC
    qacc_Errchk(qacc_Free(res));
    #else
    //delete [] res;
    free(res);
    #endif
    //res = NULL;
  }
}

inline void* aligned_alloc_no_acc(const size_t min_size)
{
  const size_t alignment = get_alignment(MemType::Cpu);
  const size_t n_elem = 1 + (min_size - 1) / alignment;
  const size_t size = n_elem * alignment;
#if defined NO_ALIGNED_ALLOC
  return malloc(size);
#else
  return aligned_alloc(alignment, size);
#endif
}


inline void qmessage(const char* fmt, ...)
{
  if(qlat::get_id_node() == 0){
    va_list args;
    va_start(args, fmt);

    char* cstr;
    int ret = vasprintf(&cstr, fmt, args);
    if (ret < 0) {
      assert(false);
    }
    const std::string str = std::string(cstr);
    printf("%s", str.c_str());
    std::free(cstr);
    va_end(args);
  }
}

////#define print0   if(qlat::get_id_node() == 0) printf

inline void abort_r(std::string stmp=std::string(""))
{
  if(stmp!=std::string(""))qmessage("%s\n",stmp.c_str());
  //MPI_Barrier(get_comm());
  MPI_Barrier(MPI_COMM_WORLD);
  fflush(stdout);
  ////MPI_Finalize();
  qlat::end();
  abort();
}

#ifdef QLAT_USE_ACC
#define gpuMalloc(bres, bsize, Ty, GPU) { \
  if(int(GPU) == -1){qacc_Errchk(qacc_MallocManaged(&bres, bsize*sizeof(Ty)));} \
  if(int(GPU) ==  0){bres = (Ty*) aligned_alloc_no_acc(bsize * sizeof(Ty));} \
  if(int(GPU) ==  1){qacc_Errchk(qacc_Malloc(&bres, bsize*sizeof(Ty)));} }
#else
#define gpuMalloc(bres,bsize, Ty, GPU) {bres = (Ty *)aligned_alloc_no_acc(bsize*sizeof(Ty));}
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

inline void free_buf(void* buf, const int GPU){
  if(buf != NULL){if(GPU){gpuFree(buf);}else{free(buf);}}
  buf = NULL;
}

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

template<typename T>
struct qlat_is_pointer { static const bool value = false; };

template<typename T>
struct qlat_is_pointer<T*> { static const bool value = true; };

template<typename Ty>
void zero_Ty(Ty* a, size_t size,int GPU=0, QBOOL dummy=QTRUE)
{
  TIMERA("zero_Ty")
  if(qlat_is_pointer<Ty>::value){
    ////printf("type is pointers! \n");
    //void* pr = (void*) a;
    //qGPU_for(isp, size, GPU, {pr[isp] = NULL;});
    return ;
  } ///donot set zero to pointers!!!

  (void)GPU;
  (void)dummy;
  #ifdef QLAT_USE_ACC
  if(GPU == 1 or GPU == -1){
    qacc_ErrCheck(qacc_MemsetAsync(a, 0, size*sizeof(Ty)));
    if(dummy==QTRUE){qacc_barrier(dummy);}
    return ;
  }
  #endif

  ////#pragma omp parallel for
  ////for(size_t isp=0;isp<size;isp++){  a[isp] = 0;}
  const Long Nv = omp_get_max_threads();
  const Long dN = (size + Nv - 1) / Nv;
  #pragma omp parallel for
  for(Long isp=0;isp<Nv;isp++)
  {
    const Long off  = isp * dN;
    const Long step = (off + dN) <= Long(size) ? dN : (size - off);
    if(step > 0 and off < Long(size)){
      Ty* tmp = &a[off + 0];
      memset((void*) tmp, 0, sizeof(Ty) * step);
    }
  }
}

template<typename Ty>
void clear_qv(qlat::vector<Ty > &G, QBOOL dummy=QTRUE)
{
  zero_Ty(G.data(), G.size(), 1 , dummy);
}

template<typename Ty>
inline crc32_t quick_checksum(Ty* buf, size_t Nsize, const Long Nsum =11, const Long hits = 13)
{
  TIMERA("quick_checksum");
  crc32_t sum = 0;
  if(hits < 0){sum = crc32_par((void*) buf, Nsize * sizeof(Ty));}
  if(hits > 0){
    Long Nuse = Nsum;
    if(Long(Nsize / hits) < Nuse){Nuse = Nsize / hits;}
    for(Long hi=0;hi < hits; hi++)
    {
      sum += crc32_par((void*) &buf[hi * Nuse], Nuse * sizeof(Ty));
    }
  }
  return sum;
}


inline unsigned int get_node_rank_funs0()
{
  int rank;
  //MPI_Comm_rank(get_comm(), &rank);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
}

}

#endif
