// utils_copy_data.h
// Gen Wang
// Jul. 2021

#ifndef UTILS_COPY_DATA_H
#define UTILS_COPY_DATA_H
#pragma once

#include <qlat/qcd.h>

namespace qlat{

#ifdef QLAT_USE_ACC
template <class T, class TInt, class TI0, class TI1>
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

template <class T, class TInt, class TI0, class TI1>
void cpy_data_from_index(T* Pres, T* Psrc, const TInt* map_res, const TInt* map_src, const TI0 Nvol, const TI1 bfac, int cpu=0, bool dummy=true)
{
  TIMER("copy data form index");

  #ifdef QLAT_USE_ACC
  if(cpu == 0){
  long Mend = (bfac*sizeof(T));
  qassert(Mend%sizeof(float) == 0);
  long nt = 32;if(Mend < 12){nt = 12;}
  dim3 dimBlock(   nt,   1, 1);
  dim3 dimGrid(  Nvol,  1, 1);
  cpy_data_from_index_global<T,TInt,TI0,TI1  ><<< dimGrid, dimBlock >>>(Pres, Psrc, &map_res[0], &map_src[0], bfac);
  if(dummy)qacc_barrier(dummy);
  return ;}
  #endif

  #pragma omp parallel for
  for(TI0 iv=0;iv<Nvol;iv++)
  {
    memcpy(&Pres[map_res[iv]*bfac], &Psrc[map_src[iv]*bfac], bfac*sizeof(T));
    //T* res = &Pres[map_res[iv]*bfac];
    //T* src = &Psrc[map_src[iv]*bfac];
    //for(TI1 j=0;j<bfac;j++){res[j] = src[j];}
  }

}


template <class T>
void cpy_data_from_index(qlat::vector<T >& res, qlat::vector<T >& src, const qlat::vector<long >& map_res, const qlat::vector<long >& map_src, const long bfac, int cpu=0, bool dummy=true)
{
  //qassert(map_res.size() ==     src.size());
  //qassert(map_res.size() ==     res.size());
  qassert(map_res.size() <= map_src.size());
  cpy_data_from_index(&res[0], &src[0], &map_res[0], &map_src[0], res.size(), bfac, cpu, dummy);
}

}



#endif
