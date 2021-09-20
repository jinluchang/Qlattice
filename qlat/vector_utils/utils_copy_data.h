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
    //memcpy(&Pres[map_res[iv]*bfac], &Psrc[map_src[iv]*bfac], bfac*sizeof(T));
    T* res = &Pres[map_res[iv]*bfac];
    T* src = &Psrc[map_src[iv]*bfac];
    for(TI1 j=0;j<bfac;j++){res[j] = src[j];}
  }

}


template <typename T>
void cpy_data_from_index(qlat::vector<T >& res, qlat::vector<T >& src, const qlat::vector<long >& map_res, const qlat::vector<long >& map_src, const long bfac, int cpu=0, bool dummy=true)
{
  //qassert(map_res.size() ==     src.size());
  //qassert(map_res.size() ==     res.size());
  qassert(map_res.size() <= map_src.size());
  /////qlat vector correct data pointer
  T* s1 = (T*) qlat::get_data(res).data();
  T* s0 = (T*) qlat::get_data(src).data();
  long* m1 = (long*) qlat::get_data(map_res).data();
  long* m0 = (long*) qlat::get_data(map_src).data();
  cpy_data_from_index(s1, s0, m1, m1, res.size(), bfac, cpu, dummy);
}

#ifdef QLAT_USE_ACC
template <typename T0, typename T1, typename TInt, int bfac>
__global__ void cpy_data_thread_global(T0* Pres, T1* Psrc,  const TInt Nvol)
{
  TInt off = blockIdx.x*blockDim.x*bfac + threadIdx.x;
  for(int i=0;i<bfac;i++)
  {
    if(off < Nvol){Pres[off] = Psrc[off];}
    off += blockDim.x;
  }
}
#endif

//////Copy data thread
template <typename T0, typename T1, typename TInt>
void cpy_data_thread(T0* Pres, T1* Psrc, const TInt Nvol, int cpu=0, bool dummy=true)
{
  TIMERA("cpy_data_thread");

  #ifdef QLAT_USE_ACC
  if(cpu == 0){
  //qacc_forNB(i, Nvol, {Pres[i]=Psrc[i];});

  const int Threads = 64;const int Biva = (4*16+sizeof(T0)-1)/sizeof(T0);
  long Nb = (Nvol + Threads*Biva -1)/(Threads*Biva);
  dim3 dimBlock(    Threads,    1, 1);
  dim3 dimGrid(     Nb,    1, 1);
  cpy_data_thread_global<T0, T1, TInt , Biva><<< dimGrid, dimBlock >>>(Pres, Psrc, Nvol);

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
template <typename Ty, bool flag, int Threads, int Biva>
__global__ void move_index_global(Ty* src, Ty* res, long sizeF, int civ, int inner)
{
  __shared__ Ty buf[Threads*Biva];

  int    tid = threadIdx.x;
  long s0    = blockIdx.x*blockDim.x;

  int Total = Threads*civ*inner;
  if(s0 + Threads > sizeF){Total = (sizeF - s0) * civ*inner;}

  int nB    = (Total + Threads-1)/Threads;
  int nC    = (Total + Biva*Threads-1)/(Biva*Threads);
    
  int ci, si, i0;
  long z0 = 0;long off = 0;long off1 = 0;
  for(int ni=0;ni < nC; ni++)
  {
    if(z0 >= Total){break;}
    if(flag){
    off = z0 + tid;
    for(int xi=0;xi<Biva;xi++)
    {
      if(off < Total){buf[xi*Threads + tid] = src[s0*civ*inner + off];off += Threads;}
    }
    __syncthreads();
    }

    off = tid;
    for(int xi=0;xi<nB;xi++)
    {
      ci = off/(Threads*inner);
      si = (off/inner)%Threads;
      i0 = off%inner;

      off1 = (si*civ + ci)*inner + i0 - z0;
      if(off1 >= 0)
      if((off1 < Threads*Biva) and (off1 < (Total - z0)) )
      {
        if( flag){res[(ci*sizeF+s0+si)*inner + i0] = buf[off1];}
        if(!flag){buf[off1] = src[(ci*sizeF+s0+si)*inner + i0];}
      }
      off += Threads;
    }
    __syncthreads();

    if(!flag){
    off = z0 + tid;
    for(int xi=0;xi<Biva;xi++)
    {
      if(off < Total){res[s0*civ*inner + off] = buf[xi*Threads + tid];off += Threads;}
    }
    __syncthreads();
    }

    z0 += Threads*Biva;
  }

}
#endif


////TODO change into Ty*
struct move_index
{
  bool GPU;

  void* buf;
  size_t buf_size;
  //qlat::vector<char* > pciv;

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
    //////psrc.resize(civ);
  }


  ////flag = 1 --> biva * sizeF * civ * size_inner --> biva * civ * sizeF * size_inner
  template <typename Ty >
  void dojob(Ty* src,Ty* res,int biva,int civ,long sizeF,int flag, int size_inner)
  {
  if(biva == 0 or civ == 0 or sizeF == 0 or size_inner == 0){return ;}
  /////size_t sizeF = sizeF0;

  ////size_t bufN = biva*civ*size_inner*sizeof(Ty)*sizeF;
  size_t Off = civ*sizeF*size_inner;
  #if PRINT_TIMER>5
  TIMER_FLOPS("reorder index");
  timer.flops += biva*Off*sizeof(Ty); 
  #endif

  ////TIMERB("reorder index");
  if(size_inner < 1){qlat::displayln_info(qlat::ssprintf("size_inner too small %d !\n", size_inner));
    MPI_Barrier(get_comm());fflush(stdout);qassert(false);
  }

  if(src == res){set_mem(civ, Off*sizeof(Ty));}
  //pciv.resize(civ);
  Ty* s0;Ty *s1;
  //#ifdef QLAT_USE_ACC
  //if(GPU)
  if(src == res)if((Off*sizeof(Ty)) % sizeof(qlat::ComplexF) != 0){
    qlat::displayln_info(qlat::ssprintf("size not divided by 16, too small. \n"));qassert(false);}
  ///#endif
 
  for(int bi=0;bi<biva;bi++){
    s0 = &src[bi*Off];
    if(src == res){s1 = (Ty*)buf;}else{s1 = (Ty*) &res[bi*Off];}
    #ifdef QLAT_USE_ACC
    if(GPU){

      {
      const int Threads = 32;const int Biva =  (16*16+sizeof(Ty)-1)/sizeof(Ty);
      long Nb = (sizeF + Threads -1)/Threads;
      dim3 dimBlock(    Threads,    1, 1);
      dim3 dimGrid(     Nb,    1, 1);
      if(flag==0)move_index_global<Ty, false , Threads, Biva><<< dimGrid, dimBlock >>>(s0, s1, sizeF, civ, size_inner);
      if(flag==1)move_index_global<Ty, true  , Threads, Biva><<< dimGrid, dimBlock >>>(s0, s1, sizeF, civ, size_inner);
      qacc_barrier(dummy);
      }

      if(src == res){
      long Nvol = long(Off*sizeof(Ty)/sizeof(qlat::ComplexF));
      cpy_data_thread((qlat::ComplexF*) &res[bi*Off], (qlat::ComplexF*) s1, Nvol, 0);
      //cpy_data_thread((Ty*) &res[bi*Off], (Ty*) s1, Off, 0);
      }

    continue ;}
    #endif

    #pragma omp parallel for
    for(long   si=0;si<sizeF;si++)
    for(int    ci=0;ci<civ;ci++)
    {
      Ty* p0;Ty* p1;
      if(flag == 1){
        p0 = (Ty*) &s0[(si*civ   + ci)*size_inner];
        p1 = (Ty*) &s1[(ci*sizeF + si)*size_inner];
      }
      if(flag == 0){
        p0 = (Ty*) &s0[(ci*sizeF + si)*size_inner];
        p1 = (Ty*) &s1[(si*civ   + ci)*size_inner];
      }
      memcpy(p1, p0, sizeof(Ty)*size_inner);
    }

    if(src == res){
      long Nvol = long(Off*sizeof(Ty)/sizeof(qlat::ComplexF));
      cpy_data_thread((qlat::ComplexF*) &res[bi*Off], (qlat::ComplexF*) s1, Nvol, 1);
    }

  }


  }


  void free_mem(){
    free_buf(buf, GPU);
    ////if(buf != NULL){
    ////  if(GPU){gpuFree(buf);}else{free(buf);}
    ////  buf = NULL;
    ////}
    buf_size = 0;
    //psrc.resize(0);
  }

  ~move_index(){
    free_mem();
  }

};


}



#endif
