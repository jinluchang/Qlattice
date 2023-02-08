// utils_COPY_data.h
// Gen Wang
// Jul. 2021

#ifndef UTILS_COPY_DATA_H
#define UTILS_COPY_DATA_H
#pragma once

#include <qlat/qcd.h>
#include "utils_float_type.h"

#define QLAT_COPY_LIMIT 1e-20

namespace qlat{

#ifdef QLAT_USE_ACC
template <typename T, typename TInt, typename TI0, typename TI1>
__global__ void cpy_data_from_index_global(T* Pres, T* Psrc, const TInt* map_res, const TInt* map_src, const TI1 bfac)
{
  TI0  index =  blockIdx.y*gridDim.x + blockIdx.x;
  unsigned int  tid    =  threadIdx.x;
  unsigned int   nt    = blockDim.x;

  float* r = (float*) &Pres[map_res[index]*bfac];
  float* s = (float*) &Psrc[map_src[index]*bfac];
  TI0 Mend = (bfac*sizeof(T))/sizeof(float);
  TI0 off = tid;
  while(off < Mend){r[off] = s[off];off += nt;}
}

#endif

//////Faster when T* is GPU memory without unified
template <typename T, typename TInt, typename TI0, typename TI1>
void cpy_data_from_index(T* Pres, T* Psrc, const TInt* map_res, const TInt* map_src, const TI0 Nvol, const TI1 bfac, int GPU=1, bool dummy=true)
{
  TIMERB("copy data form index");
  (void)dummy;
  (void)GPU;

  #ifdef QLAT_USE_ACC
  if(GPU == 1){
  long Mend = (bfac*sizeof(T));
  qassert(Mend%sizeof(float) == 0);
  long nt = 32;if(Mend < 12){nt = 12;}
  dim3 dimBlock(   nt,   1, 1);
  dim3 dimGrid(  Nvol,  1, 1);
  cpy_data_from_index_global<T,TInt,TI0,TI1  ><<< dimGrid, dimBlock >>>(Pres, Psrc, &map_res[0], &map_src[0], bfac);
  if(dummy)qacc_barrier(dummy);
  /////qacc_for(iv, Nvol, {for(int bi=0;bi<bfac;bi++)Pres[map_res[iv]*bfac+bi] = Psrc[map_src[iv]*bfac+bi];});
  return ;}
  #endif

  #pragma omp parallel for
  for(TI0 iv=0;iv<Nvol;iv++)
  {
    //memcpy(&Pres[map_res[iv]*bfac], &Psrc[map_src[iv]*bfac], bfac*sizeof(T));
    T* res = &Pres[map_res[iv]*bfac];
    T* src = &Psrc[map_src[iv]*bfac];
    for(TI1 j=0;j<bfac;j++){res[j] = src[j];}
  }

}


template <typename T>
void cpy_data_from_index(qlat::vector<T >& res, qlat::vector<T >& src, const qlat::vector<long >& map_res, const qlat::vector<long >& map_src, const long bfac, int GPU=1, bool dummy=true)
{
  //qassert(map_res.size() ==     src.size());
  //qassert(map_res.size() ==     res.size());
  qassert(map_res.size() <= map_src.size());
  /////qlat vector correct data pointer
  T* s1 = (T*) qlat::get_data(res).data();
  T* s0 = (T*) qlat::get_data(src).data();
  long* m1 = (long*) qlat::get_data(map_res).data();
  long* m0 = (long*) qlat::get_data(map_src).data();
  cpy_data_from_index(s1, s0, m1, m0, map_res.size(), bfac, GPU, dummy);
}

#ifdef QLAT_USE_ACC
template <typename T0, typename T1, typename TInt, int bfac, int ADD_FAC, typename Tadd>
__global__ void cpy_data_thread_global(T0* Pres, const T1* Psrc,  const TInt Nvol, const Tadd ADD)
{
  TInt off = blockIdx.x*blockDim.x*bfac + threadIdx.x;
  for(int i=0;i<bfac;i++)
  {
    if(ADD_FAC== 0)if(off < Nvol){Pres[off]  = Psrc[off];}
    if(ADD_FAC== 1)if(off < Nvol){Pres[off] += ADD*Psrc[off];}
    off += blockDim.x;
  }
}
#endif

//////Copy data thread, cannot give to zero with ADD = 0
template <typename T0, typename T1, typename TInt, typename Tadd>
void CPY_data_thread_basic(T0* Pres, const T1* Psrc, const TInt Nvol, int GPU, bool dummy, const Tadd ADD)
{
  (void)dummy;
  if(GPU != 0 and GPU != 1){qassert(false);}
  bool do_copy = true;
  if(qlat::qnorm(ADD) <  QLAT_COPY_LIMIT){do_copy = true ;}
  if(qlat::qnorm(ADD) >= QLAT_COPY_LIMIT){do_copy = false;}

  #ifdef QLAT_USE_ACC
  if(GPU == 1){
  /////qacc_forNB(i, Nvol, {Pres[i]=Psrc[i];});

  const int Threads = 64;const int Biva = (4*16+sizeof(T0)-1)/sizeof(T0);
  long Nb = (Nvol + Threads*Biva -1)/(Threads*Biva);
  dim3 dimBlock(    Threads,    1, 1);
  dim3 dimGrid(     Nb,    1, 1);

  if( do_copy)cpy_data_thread_global<T0, T1, TInt , Biva, 0, Tadd><<< dimGrid, dimBlock >>>(Pres, Psrc, Nvol, ADD);
  if(!do_copy)cpy_data_thread_global<T0, T1, TInt , Biva, 1, Tadd><<< dimGrid, dimBlock >>>(Pres, Psrc, Nvol, ADD);

  if(dummy)qacc_barrier(dummy);
  return ;}
  #endif

  //////===from device to device, mode 0
  if(do_copy){
    #pragma omp parallel for
    for(TInt i=0;i<Nvol;i++)
    {
      Pres[i] = Psrc[i];
    }
  }else{
    //////print0("value add %.3e %.3e", ADD.real(), ADD.imag());
    #pragma omp parallel for
    for(TInt i=0;i<Nvol;i++)
    {
      ////Pres[i] = Pres[i] + (ADD*Psrc[i]);
      Pres[i] += ADD*Psrc[i];
      //Pres[i] = 0.0;
    }
  }

}

//////Copy data thread
//0--> host to host, 1 device to device
//2--> ===from host to device, 3 ===from device to host
template <typename T0, typename T1,  typename TInt, typename Tadd>
void cpy_data_threadT(T0* Pres, const T1* Psrc, const TInt Nvol, int GPU, bool dummy, const Tadd ADD)
{
  if((void*) Pres == (void*) Psrc and qlat::qnorm(ADD) <  QLAT_COPY_LIMIT){return ;}////return if points are the same
  TIMERA("cpy_data_thread");
  ////if(Pres == Psrc){return ;}
  ////0--> host to host, 1 device to device
  if(GPU == 0 or GPU == 1){CPY_data_thread_basic(Pres, Psrc, Nvol, GPU, dummy, ADD);return ;}
  /////if(GPU == 2 or GPU == 3){if(ADD != 0){qassert(false);}}

  #ifdef QLAT_USE_ACC
  //////===from host to device
  if(GPU ==  2){
  /////qassert(sizeof(T0) == sizeof(T1));
  if(sizeof(T0) == sizeof(T1) and qlat::qnorm(ADD) <  QLAT_COPY_LIMIT){
    gpuErrchk(cudaMemcpyAsync(Pres, Psrc , Nvol*sizeof(T0), cudaMemcpyHostToDevice));
    if(dummy)qacc_barrier(dummy);
  }else{
    qlat::vector_acc< T0 > buf;buf.resize(Nvol);T0* s0 = (T0*) qlat::get_data(buf).data();
    /////host to host
    CPY_data_thread_basic(s0, Psrc, Nvol, 0, false, 0.0);
    /////devic to device
    CPY_data_thread_basic(Pres, s0, Nvol, 1, dummy, ADD);
  }
  return ;}

  //////===from device to host
  if(GPU ==  3){
  ////qassert(sizeof(T0) == sizeof(T1));
  if(sizeof(T0) == sizeof(T1) and qlat::qnorm(ADD) <  QLAT_COPY_LIMIT){
    gpuErrchk(cudaMemcpyAsync(Pres, Psrc , Nvol*sizeof(T0), cudaMemcpyDeviceToHost));
    if(dummy)qacc_barrier(dummy);
  }else{
    qlat::vector_acc< T0 > buf;buf.resize(Nvol);T0* s0 = (T0*) qlat::get_data(buf).data();
    /////device to device
    CPY_data_thread_basic(s0, Psrc, Nvol, 1,  true, 0.0);
    /////host to host
    CPY_data_thread_basic(Pres, s0, Nvol, 0, false, ADD);
  }
  return ;}

  #else
  CPY_data_thread_basic(Pres, Psrc, Nvol, 0, false, ADD);
  #endif
}

template <typename T0, typename T1,  typename TInt>
void cpy_data_thread(T0* Pres, const T1* Psrc, const TInt Nvol, int GPU=1, bool dummy=true, const double ADD = 0)
{
  cpy_data_threadT<T0, T1, TInt, double>(Pres, Psrc, Nvol, GPU, dummy, ADD);
}

template <typename T0, typename T1,  typename TInt>
void cpy_data_threadC(T0* Pres, const T1* Psrc, const TInt Nvol, int GPU=1, bool dummy=true, const T1 ADD = 0)
{
  cpy_data_threadT<T0, T1, TInt, T1>(Pres, Psrc, Nvol, GPU, dummy, ADD);
}

template <typename Ty>
void touch_GPU(Ty* Mres, long long Msize,long long offM = 0,long long size = -1, int mode = 1)
{
  (void)Mres;
  (void)mode;
  if(offM <= 0)return;
  long long Total = size;

  if(offM >= Msize or offM <= 0)return;
  if(Total == -1){Total = Msize - offM;}
  if(Total + offM > Msize){Total = Msize - offM;}

  #ifdef QLAT_USE_ACC
  int gpu_id = -1;
  cudaGetDevice(&gpu_id);
  cudaMemAdvise(&Mres[offM], Total*sizeof(Ty), cudaMemAdviseSetReadMostly, gpu_id);
  if(mode == 1){
  cudaMemPrefetchAsync(&Mres[offM], Total*sizeof(Ty), gpu_id, cudaStreamLegacy);}
  #endif
}

template <typename Ty>
void untouch_GPU(Ty* Mres , long long Msize)
{
  (void)Mres;
  (void)Msize;
  #ifdef QLAT_USE_ACC
  size_t Total = Msize;
  int gpu_id = -1;
  cudaGetDevice(&gpu_id);
  cudaMemAdvise(&Mres[0], Total*sizeof(Ty), cudaMemAdviseUnsetReadMostly, gpu_id);
  #endif
}

////Vol * N0 for res <== Vol * N1 for src
////N0 >= N1
template<typename Ty>
void copy_buffers_vecs(Ty *res, Ty *src,long N0, long N1, long Ncopy, size_t Vol, int GPU = 1)
{
  int GPU_set = GPU;
  #ifndef QLAT_USE_ACC
  GPU_set = 0;
  #endif
  qassert(Ncopy <= N0 );
  qassert(Ncopy <= N1 );
  if(N0 == N1 and Ncopy == N0){cpy_data_thread(res, src, Vol * N0, GPU_set, false);}
  else{
    #pragma omp parallel for
    for(size_t isp=0;isp<Vol;isp++){
      if(GPU_set==0){memcpy(&res[isp*N0+0], &src[isp*N1+0], Ncopy*sizeof(Ty));}
      #ifdef QLAT_USE_ACC
      if(GPU_set==1){cudaMemcpyAsync(&res[isp*N0+0], &src[isp*N1+0], Ncopy*sizeof(Ty), cudaMemcpyDeviceToDevice);}
      if(GPU_set==2){cudaMemcpyAsync(&res[isp*N0+0], &src[isp*N1+0], Ncopy*sizeof(Ty), cudaMemcpyHostToDevice);}
      if(GPU_set==3){cudaMemcpyAsync(&res[isp*N0+0], &src[isp*N1+0], Ncopy*sizeof(Ty), cudaMemcpyDeviceToHost);}
      #endif
    }
  }
  qacc_barrier(dummy);
}


}

#undef QLAT_COPY_LIMIT

#endif
