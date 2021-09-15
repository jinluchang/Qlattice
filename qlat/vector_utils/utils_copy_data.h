// utils_copy_data.h
// Gen Wang
// Jul. 2021

#ifndef UTILS_COPY_DATA_H
#define UTILS_COPY_DATA_H
#pragma once

#include <qlat/qcd.h>
#include "float_type.h"

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
void cpy_data_from_index(T* Pres, T* Psrc, const TInt* map_res, const TInt* map_src, const TI0 Nvol, const TI1 bfac, int cpu=0, bool dummy=true)
{
  TIMERB("copy data form index");

  #ifdef QLAT_USE_ACC
  if(cpu == 0){
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
    memcpy(&Pres[map_res[iv]*bfac], &Psrc[map_src[iv]*bfac], bfac*sizeof(T));
    //T* res = &Pres[map_res[iv]*bfac];
    //T* src = &Psrc[map_src[iv]*bfac];
    //for(TI1 j=0;j<bfac;j++){res[j] = src[j];}
  }

}


template <typename T>
void cpy_data_from_index(qlat::vector<T >& res, qlat::vector<T >& src, const qlat::vector<long >& map_res, const qlat::vector<long >& map_src, const long bfac, int cpu=0, bool dummy=true)
{
  //qassert(map_res.size() ==     src.size());
  //qassert(map_res.size() ==     res.size());
  qassert(map_res.size() <= map_src.size());
  cpy_data_from_index(&res[0], &src[0], &map_res[0], &map_src[0], res.size(), bfac, cpu, dummy);
}

//////Copy data thread
template <typename T0, typename T1, typename TInt>
void cpy_data_thread(T0* Pres, T1* Psrc, const TInt Nvol, int cpu=0, bool dummy=true)
{
  TIMERB("cpy_data_thread");

  #ifdef QLAT_USE_ACC
  if(cpu == 0){
  qacc_forNB(i, Nvol, {Pres[i]=Psrc[i];});
  if(dummy)qacc_barrier(dummy);
  return ;}
  #endif

  #pragma omp parallel for
  for(TInt i=0;i<Nvol;i++)
  {
    Pres[i]=Psrc[i];
  }
}

////flag = 1 --> biva * sizeF * civ * size_inner --> biva * civ * sizeF * size_inner
#ifdef QLAT_USE_ACC
template <typename Ty, bool flag>
__global__ void move_index_global(Ty* src, Ty* res, int size_inner)
{
  //  dim3 dimBlock(  civ ,    1, 1);
  //  dim3 dimGrid(  sizeF, biva, 1);
  size_t   sizeF = gridDim.x;
  size_t   si    = blockIdx.x;

  //unsigned int   biva  = gridDim.y;
  unsigned int   bi    = blockIdx.y;

  unsigned int   civ   = blockDim.x;
  unsigned int   ci    = threadIdx.x;

  size_t off_0 = ((bi*sizeF + si)*civ + ci)*size_inner;
  size_t off_1 = ((bi*civ+ci)*sizeF   + si)*size_inner;
  Ty* s0 = NULL; 
  Ty* s1 = NULL; 
  if( flag){s0 = &src[off_0];s1 = &res[off_1];}
  if(!flag){s0 = &src[off_1];s1 = &res[off_0];}
  for(int i=0;i<size_inner;i++)s1[i] = s0[i];

}
#endif


////TODO change into Ty*
struct move_index
{
  bool GPU;

  void* buf;
  size_t buf_size;
  qlat::vector<char* > psrc;

  move_index(bool GPU_set=false){
    #ifndef QLAT_USE_ACC
    GPU = false;
    #else
    GPU = GPU_set;
    #endif
    buf = NULL;
    buf_size = 0;
  }

  void set_mem(int civ, size_t Bsize)
  {
    TIMERA("move_index set_mem");
    if(buf_size != Bsize){
      free_mem();
      if(GPU){gpuMalloc(buf, Bsize, char);}
      else{buf = (void *)malloc(Bsize);}
      buf_size = Bsize;
    }
    psrc.resize(civ);
  }


  ////flag = 1 --> biva * sizeF * civ * size_inner --> biva * civ * sizeF * size_inner
  template <typename Ty >
  void dojob(Ty* src,Ty* res,int biva,int civ,long sizeF0,int flag, int size_inner)
  {
  size_t sizeF = sizeF0;

  TIMERB("reorder index");
  if(size_inner < 1){qlat::displayln_info(qlat::ssprintf("size_inner too small %d !\n", size_inner));
    MPI_Barrier(get_comm());fflush(stdout);qassert(false);
  }

  size_t bufN = biva*civ*size_inner*sizeof(Ty)*sizeF;

  #ifdef QLAT_USE_ACC
  if(GPU){
  
  Ty* tem;
  if(src == res){set_mem(civ, bufN);tem = (Ty*)buf;}else{tem = res;}

  if(src == res)if(bufN % sizeof(qlat::ComplexF) != 0){
    qlat::displayln_info(qlat::ssprintf("size not divided by 16, too small. \n"));qassert(false);}

  {
  dim3 dimBlock(  civ ,    1, 1);
  dim3 dimGrid(  sizeF, biva, 1);
  if(flag==0)move_index_global<Ty, false ><<< dimGrid, dimBlock >>>(src, tem, size_inner);
  if(flag==1)move_index_global<Ty, true  ><<< dimGrid, dimBlock >>>(src, tem, size_inner);
  }
  qacc_barrier(dummy);

  if(src == res){
  long Nvol = long(bufN/sizeof(qlat::ComplexF));
  qlat::ComplexF* s0 = (qlat::ComplexF*) tem;
  qlat::ComplexF* s1 = (qlat::ComplexF*) res;
  qacc_for(i, Nvol, {s1[i] = s0[i];});
  }

  return ;}
  #endif

  set_mem(civ, bufN);
  Ty* tem;tem = (Ty*)buf;
  if(flag == 1){
    ///memcpy((Ty*)&tem[0],(Ty*)&src[0],biva*sizeof(Ty)*sizeF*civ*size_inner);
    cpy_data_thread((Ty*)&tem[0], (Ty*)&res[0], biva*sizeF*civ*size_inner, 1);
  }
  for(size_t bi=0;bi<size_t(biva);bi++)
  {
    for(int ci=0;ci<civ;ci++)
    {
      if(flag==0)psrc[ci] = (char*) &src[(bi*civ*sizeF + ci*sizeF + 0)*size_inner];
      if(flag==1)psrc[ci] = (char*) &res[(bi*civ*sizeF + ci*sizeF + 0)*size_inner];
    }

    #pragma omp parallel for
    for(size_t si=0;si<sizeF;si++)
    for(int    ci=0;ci<civ;ci++)
    {
      Ty* s0 = (Ty*) &tem[(bi*sizeF*civ + si*civ + ci)*size_inner];
      Ty* s1 = (Ty*) &(psrc[ci][si*sizeof(Ty)*size_inner]);
      if(flag==0){
        //for(int i0=0;i0<size_inner;i0++){s0[i0] = s1[i0];}
        memcpy(s0, s1,sizeof(Ty)*size_inner);
      }
      if(flag==1){
        //for(int i0=0;i0<size_inner;i0++){s1[i0] = s0[i0];}
        memcpy(s1, s0,sizeof(Ty)*size_inner);
      }
    }
  }
  if(flag == 0){
    ////memcpy((Ty*)&res[0],(Ty*)&tem[0],biva*sizeof(Ty)*sizeF*civ*size_inner);
    cpy_data_thread((Ty*)&res[0], (Ty*)&tem[0], biva*sizeF*civ*size_inner, 1);
  }


  }


  void free_mem(){
    if(buf != NULL){
      if(GPU){gpuFree(buf);}else{free(buf);}
      buf = NULL;
    }
    buf_size = 0;
    psrc.resize(0);
  }

  ~move_index(){
    free_mem();
  }

};


}



#endif
