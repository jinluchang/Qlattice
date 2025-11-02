// utils_float_type.h
// Gen Wang
// Jun. 2021

#ifndef UTILS_FLOAT_TYPE_H
#define UTILS_FLOAT_TYPE_H
#pragma once

#include "utils_geo_cache.h"

namespace qlat{

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

inline QMEM get_type_mem(Int g){
  if(g == -1){return QMSYNC;}
  if(g ==  0){return QMCPU ;}
  if(g ==  1){return QMGPU ;}
  return QMCPU;
}

inline Int check_GPU_same(QMEM a, QMEM b)
{
  if(a == QMSYNC){return 1;}
  if(b == QMSYNC){return 1;}
  if(a == b){return 1;}
  return 0;
}

inline Int check_GPU_same(QMEM a, bool b)
{
  if(a == QMSYNC){return 1;}
  if(a == QMGPU and b == 1){return 1;}
  if(a == QMCPU and b == 0){return 1;}
  return 0;
}

inline Int check_GPU_multi(QMEM a, QMEM b)
{
  if(a == QMSYNC){return b;}
  if(b == QMSYNC){return a;}
  if(a == b){return a;}
  return -2;
}

inline Int check_GPU_multi(QMEM a, bool b)
{
  if(a == QMSYNC){return b;}
  if(a == QMGPU and b == 1){return 1;}
  if(a == QMCPU and b == 0){return 0;}
  return -2;
}

inline void qmessage(const char* fmt, ...)
{
  if(qlat::get_id_node() == 0){
    va_list args;
    va_start(args, fmt);

    char* cstr;
    Int ret = vasprintf(&cstr, fmt, args);
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

inline void free_buf(void* buf, const Int GPU){
  if(buf != NULL){if(GPU){gpuFree(buf);}else{free(buf);}}
  buf = NULL;
}

#ifdef QLAT_USE_ACC
#define gpuMalloc(bres, bsize, Ty, GPU) { \
  if(int(GPU) == -1){qacc_Errchk(qacc_MallocManaged(&bres, bsize*sizeof(Ty)));} \
  if(int(GPU) ==  0){bres = (Ty*) aligned_alloc_no_acc(bsize * sizeof(Ty));} \
  if(int(GPU) ==  1){qacc_Errchk(qacc_Malloc(&bres, bsize*sizeof(Ty)));} }
#else
#define gpuMalloc(bres,bsize, Ty, GPU) {bres = (Ty *)aligned_alloc_no_acc(bsize*sizeof(Ty));}
#endif

#ifndef QLAT_USE_ACC
#ifdef __GNUC__
///////Multiply of different types of complex

//inline std::complex<RealD> operator*(const std::complex<float> &a, const double &b) {
//    return std::complex<RealD>(a.real()*b, a.imag()*b);
//}
//
//inline std::complex<RealD> operator/(const std::complex<float> &a, const double &b) {
//    return std::complex<RealD>(a.real()/b, a.imag()/b);
//}

inline std::complex<RealD> operator*(const std::complex<float> &a, const std::complex<RealD > &b) {
    return std::complex<RealD>(a.real() * b.real() - a.imag()*b.imag(), a.imag()*b.real() + a.real()*b.imag());
}

inline std::complex<RealD> operator*(const std::complex<RealD > &b, const std::complex<float> &a) {
    return std::complex<RealD>(a.real() * b.real() - a.imag()*b.imag(), a.imag()*b.real() + a.real()*b.imag());
}

inline std::complex<RealD> operator-(const std::complex<RealD > &a, const std::complex<float> &b) {
    return std::complex<RealD>(a.real() - b.real() , a.imag() - b.imag());
}
inline std::complex<RealD> operator-(const std::complex<float > &a, const std::complex<RealD> &b) {
    return std::complex<RealD>(a.real() - b.real() , a.imag() - b.imag());
}
#endif
#endif

template<typename T>
struct qlat_is_pointer { static const bool value = false; };

template<typename T>
struct qlat_is_pointer<T*> { static const bool value = true; };

template<typename Ty>
void zero_Ty(Ty* a, size_t size,Int GPU=0, QBOOL dummy=QTRUE)
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
  Int rank;
  ///MPI_Comm_rank(get_comm(), &rank);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
}

}

#endif
