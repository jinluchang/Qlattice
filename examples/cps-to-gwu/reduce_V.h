#ifndef reduce_V_h
#define reduce_V_h
#pragma once

/////Reference https://developer.download.nvidia.com/assets/cuda/files/reduction.pdf

#include <qlat/qlat.h>

#ifdef QLAT_USE_ACC
template <unsigned int blockSize, typename Ty>
__global__ void reduce6(const Ty *g_idata, Ty *g_odata, unsigned long n,unsigned int divide){
  extern __shared__ Ty sdata[];
  unsigned int tid = threadIdx.x;
  unsigned int iv  = blockIdx.y;
  sdata[tid] = 0;
  unsigned long i  = iv*n + (blockIdx.x*(blockSize) + tid)*divide;
  unsigned long i_end = i + divide;if(i_end > (iv+1)*n){i_end = (iv+1)*n;}
  while (i < i_end) { sdata[tid] += g_idata[i]; i += 1; }
  __syncthreads();

  if (blockSize >=1024) { if (tid < 512) { sdata[tid] += sdata[tid + 512]; } __syncthreads(); }
  if (blockSize >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
  if (blockSize >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
  if (blockSize >= 128) { if (tid <  64) { sdata[tid] += sdata[tid +  64]; } __syncthreads(); }

  unsigned int start = 32;if(blockSize<=32){start=blockSize/2;}
  for (unsigned int s=start; s>0; s>>=1) {
  if (tid < s) {
  sdata[tid] += sdata[tid + s];
  }
  __syncthreads();
  }
  if (tid == 0) g_odata[iv*gridDim.x + blockIdx.x] = sdata[0];

}
#endif

namespace qlat{

template<typename Ty>
inline void reduce_T_global6(const Ty* src,Ty* res,const long n, const int nv,long nt, long blockS)
{
  ///const long blockS = (n + nt - 1)/(nt);
  dim3 dimGrid(blockS, nv, 1);
  dim3 dimBlock(nt, 1, 1);

  long threads = nt;
  long smemSize = nt*sizeof(Ty);
  unsigned int divide = (n+nt*blockS-1)/(nt*blockS);
  /////reduce6 <<<cu_blocks,cu_threads, nt*sizeof(Ty) >>>(src,res,n);

  #ifdef QLAT_USE_ACC
  switch (threads)
  {
    case 1024:reduce6<1024,Ty><<< dimGrid, dimBlock, smemSize >>>(src, res, n, divide); break;
    case 512: reduce6< 512,Ty><<< dimGrid, dimBlock, smemSize >>>(src, res, n, divide); break;
    case 256: reduce6< 256,Ty><<< dimGrid, dimBlock, smemSize >>>(src, res, n, divide); break;
    case 128: reduce6< 128,Ty><<< dimGrid, dimBlock, smemSize >>>(src, res, n, divide); break;
    case 64:  reduce6<  64,Ty><<< dimGrid, dimBlock, smemSize >>>(src, res, n, divide); break;
    case 32:  reduce6<  32,Ty><<< dimGrid, dimBlock, smemSize >>>(src, res, n, divide); break;
    case 16:  reduce6<  16,Ty><<< dimGrid, dimBlock, smemSize >>>(src, res, n, divide); break;
    case 8:   reduce6<   8,Ty><<< dimGrid, dimBlock, smemSize >>>(src, res, n, divide); break;
    case 4:   reduce6<   4,Ty><<< dimGrid, dimBlock, smemSize >>>(src, res, n, divide); break;
    case 2:   reduce6<   2,Ty><<< dimGrid, dimBlock, smemSize >>>(src, res, n, divide); break;
    case 1:   reduce6<   1,Ty><<< dimGrid, dimBlock, smemSize >>>(src, res, n, divide); break;
  }
  #endif
  qacc_barrier(dummy);
}

template<typename Ty>
inline void reduce_cpu(const Ty *src,Ty &res,const unsigned long n){
  ////#pragma omp parallel for reduction(+: tem)
  int Nv = omp_get_num_threads();
  if(Nv == 1){
  for(unsigned long index=0;index<n;index++){
    res += src[index];
  }}
  else{
    std::vector<Ty > buf;buf.resize(Nv);
    for(int iv=0;iv<Nv;iv++){buf[iv]=0.0;}
    #pragma omp parallel
    for(unsigned long index=0;index<n;index++){
      buf[omp_get_thread_num()] += src[index];
    }
    for(int iv=0;iv<Nv;iv++){res += buf[iv];}
  }
}

template<typename Ty>
inline void reduce_gpu2d_6(const Ty* src,Ty* res,long n, int nv=1,
    int thread_pow2 = 8,int divide=128,int fac=16)
{
  #ifdef QLAT_USE_ACC
  long nthreads = qlat::qacc_num_threads();
  nthreads = nextPowerOf2(nthreads);

  ////long nfac = 1<<thread_pow2;
  long nfac = thread_pow2;
  long cutN = nthreads*fac;
  long nt  = nthreads*nfac;if(nt > 1024){nt = 1024;nfac=nt/nthreads;}

  long ntL = nt*divide;
  long Ny0 = (n  + ntL - 1)/(ntL);
  long Ny1 = (Ny0 + ntL - 1)/(ntL);
  Ty *psrc;Ty *pres;Ty *tem;
  qlat::vector<Ty > buf0,buf1;
  buf0.resize(nv*Ny0);buf1.resize(nv*Ny1);
  pres = &buf0[0];

  if(n <= cutN){
    //for(int i=0;i<nv;i++)reduce_cpu(&src[i*n],res[i],n);return;
    reduce_T_global6(&src[0],&pres[0], n, nv, nt, 1);
    for(int i=0;i<nv;i++)res[i] += pres[i];return;
  }

  
  reduce_T_global6(src,pres, n, nv, nt, Ny0);
  ////Nres0 = Ny;
  psrc = &buf0[0];pres=&buf1[0];

  ////for(int i=0;i<nv;i++)reduce_cpu(&psrc[i*Ny0],res[i],Ny0);return;

  for(int si=0;si<1000;si++){
    if(Ny0 <= cutN){
      for(int i=0;i<nv;i++)reduce_cpu(&psrc[i*Ny0],res[i],Ny0);return;
      //reduce_T_global6(psrc,pres, Ny0, nv, nt, 1);
      //for(int i=0;i<nv;i++)res[i] += pres[i];return;
    }
    Ny1 = (Ny0 + ntL - 1)/(ntL);
    reduce_T_global6(psrc,pres, Ny0, nv, nt, Ny1);
    /////Switch psrc, pres
    tem = pres;pres = psrc;psrc = tem;
    Ny0 = Ny1;
  }
  reduce_T_global6(psrc,pres, Ny0, nv, nt, 1);
  for(int i=0;i<nv;i++)res[i] += pres[i];return;
  /////for(int i=0;i<nv;i++)reduce_cpu(&psrc[i*Ny0],res[i],Ny0);return;
  #endif
  
  #ifndef QLAT_USE_ACC
  for(int i=0;i<nv;i++)reduce_cpu(&src[i*n],res[i],n);return;
  #endif
}

template<typename Ty>
inline unsigned long reduce_T(const Ty *src,Ty *res,const unsigned long n,const int nv,const unsigned long Lx){
  unsigned long Ny = n/Lx;
  unsigned long Nx = ((n+Ny-1)/Ny);
  qacc_for(index, Ny, {
    for(int iv=0;iv<nv;iv++){
      res[iv*Ny+index] = 0;
      unsigned long tid = iv*n + index*Nx;
      unsigned long end = tid + Nx;
      if(end>(iv+1)*n){end=(iv+1)*n;}
      for(unsigned long i=tid;i<end;i++){
        res[iv*Ny+index] += src[i];
      }
    }
  });

  return Ny;
}

template<typename Ty>
inline void reduce_gpu(const Ty *src,Ty *res,const long n,const int nv=1,
  const int Ld=128,const int Ld0=8,const int fac=16)
{
  #ifdef QLAT_USE_ACC
  const int cutN = qlat::qacc_num_threads()*fac;
  unsigned long Ny = n/Ld;
  if(n <= cutN){for(int i=0;i<nv;i++)reduce_cpu(&src[i*n],res[i],n);return;}

  Ty *psrc;Ty *pres;Ty *tem;
  long Nres;

  qlat::vector<Ty > buf0,buf1;
  buf0.resize(nv*Ny);
  pres = &buf0[0];
  Nres = reduce_T(src,pres, n, nv, Ld);
  psrc = &buf0[0];
  buf1.resize(nv*Nres);pres = &buf1[0];

  for(int si=0;si<1000;si++){
    if(Nres <= cutN){for(int i=0;i<nv;i++)reduce_cpu(&psrc[i*Nres],res[i],Nres);return;}
    Nres = reduce_T(psrc,pres, Nres, nv, Ld0);
    /////Switch psrc, pres
    tem = pres;pres = psrc;psrc = tem;
  }
  for(int i=0;i<nv;i++)reduce_cpu(&psrc[i*Nres],res[i],Nres);return;
  #endif

  #ifndef QLAT_USE_ACC
  for(int i=0;i<nv;i++)reduce_cpu(&src[i*n],res[i],n);return;
  #endif
}

template<typename Ty>
void reduce_vec(const Ty* src,Ty* res,long n, int nv=1)
{
  int thread_pow2 = 8;
  int divide = 128;
  int fac = 16;

  reduce_gpu2d_6(src,res,n,nv, thread_pow2,divide, fac);
}


//template<typename Ty>
//qacc void reduce_gpu(const qlat::Vector<Ty > src,Ty &res)
//{
//  const long n = src.size()*sizeof(M);
//  reduce_gpu(src.p,res,n);
//}



}


#endif

