// utils_reduce_vec.h
// Gen Wang
// Apr. 2021

#ifndef UTILS_REDUCE_VEC_H
#define UTILS_REDUCE_VEC_H
#pragma once

/////Reference https://developer.download.nvidia.com/assets/cuda/files/reduction.pdf

#include <qlat/qcd.h>

#ifdef QLAT_USE_ACC
template <unsigned int blockSize, typename Ty>
__global__ void reduce6(const Ty *g_idata, Ty *g_odata, unsigned long n,unsigned int divide){
  ////extern __shared__ Ty sdata[];
  __shared__ Ty sdata[blockSize];
  unsigned int tid = threadIdx.x;
  unsigned int iv  = blockIdx.y;
  sdata[tid] = 0;
  //unsigned long i  = iv*n + (blockIdx.x*(blockSize) + tid)*divide;
  //unsigned long i_end = i + divide;if(i_end > (iv+1)*n){i_end = (iv+1)*n;}
  //while (i < i_end) { sdata[tid] += g_idata[i]; i += 1; }
  unsigned long i  = iv*n + blockIdx.x*blockSize + tid;
  unsigned long i_end = (iv+1)*n; unsigned long off_divide = gridDim.x*blockSize;
  while (i < i_end) { sdata[tid] += g_idata[i]; i += off_divide; }
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

/////if not power of 2, return the first larger number with power of 2
inline unsigned long nextPowerOf2(unsigned long n)
{
    unsigned long count = 0;
    // First n in the below condition
    // is for the case where n is 0
    if (n && !(n & (n - 1)))
        return n;
    while( n != 0)
    {
        n >>= 1;
        count += 1;
    }
    return 1 << count;
}

#ifdef QLAT_USE_ACC
template<typename Ty>
inline void reduce_T_global6(const Ty* src,Ty* res,const long n, const int nv,long nt, long blockS)
{
  ///const long blockS = (n + nt - 1)/(nt);
  dim3 dimGrid(blockS, nv, 1);
  dim3 dimBlock(nt, 1, 1);

  long threads = nt;
  ////long smemSize = nt*sizeof(Ty);
  unsigned int divide = (n+nt*blockS-1)/(nt*blockS);
  /////reduce6 <<<cu_blocks,cu_threads, nt*sizeof(Ty) >>>(src,res,n);

  switch (threads)
  {
    case 1024:reduce6<1024,Ty><<< dimGrid, dimBlock >>>(src, res, n, divide); break;
    case 512: reduce6< 512,Ty><<< dimGrid, dimBlock >>>(src, res, n, divide); break;
    case 256: reduce6< 256,Ty><<< dimGrid, dimBlock >>>(src, res, n, divide); break;
    case 128: reduce6< 128,Ty><<< dimGrid, dimBlock >>>(src, res, n, divide); break;
    case 64:  reduce6<  64,Ty><<< dimGrid, dimBlock >>>(src, res, n, divide); break;
    case 32:  reduce6<  32,Ty><<< dimGrid, dimBlock >>>(src, res, n, divide); break;
    case 16:  reduce6<  16,Ty><<< dimGrid, dimBlock >>>(src, res, n, divide); break;
    case 8:   reduce6<   8,Ty><<< dimGrid, dimBlock >>>(src, res, n, divide); break;
    case 4:   reduce6<   4,Ty><<< dimGrid, dimBlock >>>(src, res, n, divide); break;
    case 2:   reduce6<   2,Ty><<< dimGrid, dimBlock >>>(src, res, n, divide); break;
    case 1:   reduce6<   1,Ty><<< dimGrid, dimBlock >>>(src, res, n, divide); break;
  }
  qacc_barrier(dummy);
}
#endif

template<typename Ty>
void reduce_cpu(const Ty *src,Ty &res,const long n){
  //#pragma omp parallel for reduction(+: res)
  //for(unsigned long index=0;index<n;index++){
  //  res += src[index];
  //}
  //int Nv = omp_get_num_threads();
  int Nv = omp_get_max_threads();
  //if(n < 3*Nv)Nv=1;
  if(n < 10*Nv)Nv=1;
  //int Nv = 1;
  if(Nv == 1){
  for(long index=0;index<n;index++){
    res += src[index];
  }}
  else{
    ////print0("====Reduce omp \n");
    omp_set_num_threads(omp_get_max_threads());
    std::vector<Ty > buf;buf.resize(Nv);
    for(int iv=0;iv<Nv;iv++){buf[iv]=0.0;}
    size_t bsize = (n + Nv-1)/Nv;
    #pragma omp parallel for
    for(int ompi=0;ompi<Nv;ompi++)
    {
      int temi = omp_get_thread_num();
      Ty &pres = buf[temi];
      const Ty *psrc = &src[temi*bsize];
      size_t bsize_end = bsize;
      if(temi*bsize + bsize > (LInt) n){bsize_end = n - temi*bsize;}
      for(size_t isp=0;isp<bsize_end;isp++){pres += psrc[isp];}
    }

    //#pragma omp parallel for
    //for(unsigned long index=0;index<n;index++){
    //  buf[omp_get_thread_num()] += src[index];
    //}
    for(int iv=0;iv<Nv;iv++){res += buf[iv];}
  }
}


template<typename Ty>
inline void reduce_gpu2d_6(const Ty* src,Ty* res,long n, int nv=1,
    int thread_pow2 = 8,int divide=128,int fac=16)
{
  (void)thread_pow2;
  (void)divide;
  (void)fac;
  #ifdef QLAT_USE_ACC
  //long nthreads = qlat::qacc_num_threads();
  long nthreads = 32;
  nthreads = nextPowerOf2(nthreads);

  ////long nfac = 1<<thread_pow2;
  long nfac = thread_pow2;
  long cutN = nthreads*fac;
  long nt  = nthreads*nfac;if(nt > 1024){nt = 1024;nfac=nt/nthreads;}

  long ntL = nt*divide;
  long Ny0 = (n  + ntL - 1)/(ntL);
  long Ny1 = (Ny0 + ntL - 1)/(ntL);
  Ty *psrc;Ty *pres;Ty *tem;
  qlat::vector_acc<Ty > buf0,buf1;
  buf0.resize(nv*Ny0);buf1.resize(nv*Ny1);
  pres = &buf0[0];

  if(n <= cutN){
    //for(int i=0;i<nv;i++)reduce_cpu(&src[i*n],res[i],n);return;
    reduce_T_global6(&src[0],&pres[0], n, nv, nt, 1);
    #pragma omp parallel for
    for(int i=0;i<nv;i++)res[i] += pres[i];
    return;
  }

  /////for(int iv=0;iv<nv;iv++)reduce_cpu(&src[iv*n],res[iv],n);return;

  reduce_T_global6(src,pres, n, nv, nt, Ny0);
  ////Nres0 = Ny;
  psrc = &buf0[0];pres=&buf1[0];

  //for(int i=0;i<nv;i++)reduce_cpu(&psrc[i*Ny0],res[i],Ny0);return;

  for(int si=0;si<1000;si++){
    if(Ny0 <= cutN){
      #pragma omp parallel for
      for(int i=0;i<nv;i++)reduce_cpu(&psrc[i*Ny0],res[i],Ny0);
      return;
      //reduce_T_global6(psrc,pres, Ny0, nv, nt, 1);
      //#pragma omp parallel for
      //for(int i=0;i<nv;i++)res[i] += pres[i];return;
    }
    Ny1 = (Ny0 + ntL - 1)/(ntL);
    reduce_T_global6(psrc,pres, Ny0, nv, nt, Ny1);
    /////Switch psrc, pres
    tem = pres;pres = psrc;psrc = tem;
    Ny0 = Ny1;
  }
  reduce_T_global6(psrc,pres, Ny0, nv, nt, 1);
  #pragma omp parallel for
  for(int i=0;i<nv;i++)res[i] += pres[i];
  return;
  /////for(int i=0;i<nv;i++)reduce_cpu(&psrc[i*Ny0],res[i],Ny0);return;
  #endif
  
  #ifndef QLAT_USE_ACC
  for(int i=0;i<nv;i++)reduce_cpu(&src[i*n],res[i],n);
  return;
  #endif
}

template<typename Ty>
inline unsigned long reduce_T(const Ty *src,Ty *res,const unsigned long n,const int nv,const unsigned long Lx){
  unsigned long Ny = n/Lx;
  unsigned long Nx = ((n+Ny-1)/Ny);
  qacc_for(index, (long)Ny, {
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
void reduce_cpu(const Ty* src,Ty* res,long n, int nv)
{
  for(long i=0;i<nv;i++){reduce_cpu(&src[i*n], res[i],n);}
}

template<typename Ty>
inline void reduce_gpu(const Ty *src,Ty *res,const long n,const int nv=1,
  const int Ld=128,const int Ld0=8,const int fac=16)
{
  (void)Ld;
  (void)Ld0;
  (void)fac;
  #ifdef QLAT_USE_ACC
  //const int cutN = qlat::qacc_num_threads()*fac;
  const int cutN = 32*fac;
  unsigned long Ny = n/Ld;
  if(n <= cutN){for(int i=0;i<nv;i++){reduce_cpu(&src[i*n],res[i],n);}return;}

  Ty *psrc;Ty *pres;Ty *tem;
  long Nres;

  qlat::vector<Ty > buf0,buf1;
  buf0.resize(nv*Ny);
  pres = &buf0[0];
  Nres = reduce_T(src,pres, n, nv, Ld);
  psrc = &buf0[0];
  buf1.resize(nv*Nres);pres = &buf1[0];

  for(int si=0;si<1000;si++){
    if(Nres <= cutN){for(int i=0;i<nv;i++){reduce_cpu(&psrc[i*Nres],res[i],Nres);}return;}
    Nres = reduce_T(psrc,pres, Nres, nv, Ld0);
    /////Switch psrc, pres
    tem = pres;pres = psrc;psrc = tem;
  }
  for(int i=0;i<nv;i++)reduce_cpu(&psrc[i*Nres],res[i],Nres);
  return;
  #endif

  #ifndef QLAT_USE_ACC
  for(int i=0;i<nv;i++)reduce_cpu(&src[i*n],res[i],n);
  return;
  #endif
}


template<typename Ty>
void reduce_vec(const Ty* src, Ty* res,long n, int nv=1)
{
  #ifndef QLAT_USE_ACC
  reduce_cpu(src, res, n , nv);
  return ;
  #else
  int thread_pow2 = 1;
  int divide = 256;
  int fac    = 4;

  //if(nt_use == 1)nt_use = omp_get_num_threads();
  //unsigned long nv = nt_use*Aoper;
  //#ifdef QLAT_USE_ACC
  //cudaDeviceProp prop;
  //cudaGetDeviceProperties(&prop, 0);
  //unsigned int nthreads = omp_get_num_threads();
  //unsigned long cores = prop.multiProcessorCount;
  //unsigned long maxthreads = prop.maxThreadsPerMultiProcessor;
  //unsigned long Fullthreads = maxthreads*cores;
  //unsigned long maxblock = 8*(Fullthreads + nv-1)/nv;
  //print0("====cores %8d, maxthreads %8d, maxblock %8d \n",cores,maxthreads,maxblock);
  //if(blockS_use > maxblock)blockS_use = maxblock;
  //#endif
  reduce_gpu2d_6(src, res, n, nv, thread_pow2,divide, fac);
  #endif

}


}


#endif

