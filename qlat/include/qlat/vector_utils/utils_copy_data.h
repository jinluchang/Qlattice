// utils_copy_data.h
// Gen Wang
// Jul. 2021

#ifndef UTILS_COPY_DATA_H
#define UTILS_COPY_DATA_H
#pragma once

#include <qlat/qcd.h>
#include "utils_float_type.h"

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
template <typename T0, typename T1, typename TInt, int bfac, int ADD_FAC>
__global__ void cpy_data_thread_global(T0* Pres, const T1* Psrc,  const TInt Nvol, const double ADD)
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
template <typename T0, typename T1, typename TInt>
void CPY_data_thread_basic(T0* Pres, const T1* Psrc, const TInt Nvol, int GPU=1, bool dummy=true, const double ADD = 0)
{
  (void)dummy;
  if(GPU != 0 and GPU != 1){qassert(false);}
  bool do_copy = true;
  if(qlat::qnorm(ADD) <  1e-13){do_copy = true ;}
  if(qlat::qnorm(ADD) >= 1e-13){do_copy = false;}

  #ifdef QLAT_USE_ACC
  if(GPU == 1){
  /////qacc_forNB(i, Nvol, {Pres[i]=Psrc[i];});

  const int Threads = 64;const int Biva = (4*16+sizeof(T0)-1)/sizeof(T0);
  long Nb = (Nvol + Threads*Biva -1)/(Threads*Biva);
  dim3 dimBlock(    Threads,    1, 1);
  dim3 dimGrid(     Nb,    1, 1);

  if( do_copy)cpy_data_thread_global<T0, T1, TInt , Biva, 0><<< dimGrid, dimBlock >>>(Pres, Psrc, Nvol, ADD);
  if(!do_copy)cpy_data_thread_global<T0, T1, TInt , Biva, 1><<< dimGrid, dimBlock >>>(Pres, Psrc, Nvol, ADD);

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
template <typename T0, typename T1,  typename TInt>
void cpy_data_thread(T0* Pres, const T1* Psrc, const TInt Nvol, int GPU=1, bool dummy=true, const double ADD = 0)
{
  if((void*) Pres == (void*) Psrc){return ;}////return if points are the same
  TIMERA("cpy_data_thread");
  ////if(Pres == Psrc){return ;}
  ////0--> host to host, 1 device to device
  if(GPU == 0 or GPU == 1){CPY_data_thread_basic(Pres, Psrc, Nvol, GPU, dummy, ADD);return ;}
  /////if(GPU == 2 or GPU == 3){if(ADD != 0){qassert(false);}}

  #ifdef QLAT_USE_ACC
  //////===from host to device
  if(GPU ==  2){
  /////qassert(sizeof(T0) == sizeof(T1));
  if(sizeof(T0) == sizeof(T1) and qlat::qnorm(ADD) <  1e-13){
    gpuErrchk(cudaMemcpyAsync(Pres, Psrc , Nvol*sizeof(T0), cudaMemcpyHostToDevice));
    if(dummy)qacc_barrier(dummy);
  }else{
    qlat::vector_acc< T0 > buf;buf.resize(Nvol);T0* s0 = (T0*) qlat::get_data(buf).data();
    /////host to host
    CPY_data_thread_basic(s0, Psrc, Nvol, 0, false);
    /////devic to device
    CPY_data_thread_basic(Pres, s0, Nvol, 1, dummy, ADD);
  }
  return ;}

  //////===from device to host
  if(GPU ==  3){
  ////qassert(sizeof(T0) == sizeof(T1));
  if(sizeof(T0) == sizeof(T1) and qlat::qnorm(ADD) <  1e-13){
    gpuErrchk(cudaMemcpyAsync(Pres, Psrc , Nvol*sizeof(T0), cudaMemcpyDeviceToHost));
    if(dummy)qacc_barrier(dummy);
  }else{
    qlat::vector_acc< T0 > buf;buf.resize(Nvol);T0* s0 = (T0*) qlat::get_data(buf).data();
    /////device to device
    CPY_data_thread_basic(s0, Psrc, Nvol, 1, true);
    /////host to host
    CPY_data_thread_basic(Pres, s0, Nvol, 0, false, ADD);
  }
  return ;}

  #else
  CPY_data_thread_basic(Pres, Psrc, Nvol, 0, false, ADD);
  #endif

}

template <typename Ty>
void touch_GPU(Ty* Mres, long long Msize,long long offM = 0,long long size = -1, int mode = 1)
{
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
  #ifdef QLAT_USE_ACC
  size_t Total = Msize;
  int gpu_id = -1;
  cudaGetDevice(&gpu_id);
  cudaMemAdvise(&Mres[0], Total*sizeof(Ty), cudaMemAdviseUnsetReadMostly, gpu_id);
  #endif
}



}



#endif
