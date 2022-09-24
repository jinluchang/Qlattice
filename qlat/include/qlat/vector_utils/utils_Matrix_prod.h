// utils_GPU_Matrix_prod.h
// Gen Wang
// Jun. 2021
#ifndef UTILS_MATRIX_PROD_H
#define UTILS_MATRIX_PROD_H

#pragma once
#include "utils_float_type.h"

#define EML  Eigen::Map< Eigen::Matrix<Ty , Eigen::Dynamic, Eigen::Dynamic ,Eigen::RowMajor> >
#define EMLC Eigen::Map< Eigen::Matrix<Ty , Eigen::Dynamic, Eigen::Dynamic ,Eigen::ColMajor> >

//////  https://docs.nvidia.com/cuda/pdf/CUDA_C_Best_Practices_Guide.pdf

namespace qlat{
////need Conjugate in first element

#ifdef QLAT_USE_ACC
//////shared memory 16KB or 48KB per thread
template <unsigned int blockSize, unsigned int sm, unsigned int sn, unsigned int ch,unsigned int Conj, unsigned int trans, typename Ty >
__global__ void matrix_prod_global2(Ty* a, Ty* b, Ty* c, const long m, const long n, const long w)
{
  __shared__ Ty as[sm+1][blockSize+1];
  __shared__ Ty bs[blockSize+1][sn+1];

  ////////======c, mxn, a, mxw, b, nxw
  unsigned int  tid=  threadIdx.y*blockDim.x + threadIdx.x;

  unsigned long ni =  blockIdx.x * blockDim.x + threadIdx.x;
  unsigned long mi =  blockIdx.y * blockDim.y + threadIdx.y;
  unsigned long li =  blockIdx.z;

  unsigned long ci =  (li*m + mi)*n + ni;
  Ty* ap;
  Ty* bp;

  Ty buf;buf = 0.0;

  //unsigned long wini = 0*blockSize;
  //while(wini < w){
  //  ap = &a[(li*m + blockIdx.y * blockDim.y + 0)*w + wini];
  //  bp = &b[(li*n + blockIdx.x * blockDim.x + 0)*w + wini];

  //  for(int i=0;i<sm;i++){as[i][tid] = ap[i*w + tid ];}

  //  for(int i=0;i<sn;i++){bs[tid][i] = bp[i*w + tid ];}

  //  __syncthreads();

  //  for(int i = 0; i < blockSize; i++)
  //  {
  //    buf += as[threadIdx.y][i] * bs[i][threadIdx.x];
  //  }
  //  wini += blockSize;
  //}
  //c[ci] += buf;

  unsigned long wini = 0*blockSize;
  unsigned int  wcut =   blockSize;
  unsigned int smcut = sm;
  unsigned int sncut = sn;

  if(ch==1)if(blockIdx.y * blockDim.y + sm > m){smcut = m - blockIdx.y * blockDim.y;}
  if(ch==1)if(blockIdx.x * blockDim.x + sn > n){sncut = n - blockIdx.x * blockDim.x;}

  while(wini < w){
    if(ch == 1)if(wini + blockSize > w){wcut = w - wini;}

    ap = &a[(li*m + blockIdx.y * blockDim.y + 0)*w + wini];
    if(trans == 0){bp = &b[(li*n + blockIdx.x * blockDim.x + 0)*w + wini];}
    if(trans == 1){bp = &b[(li*w + wini)*n + blockIdx.x * blockDim.x + 0];}

    if(tid < wcut){
      for(int i=0;i<smcut;i++){as[i][tid] = ap[i*w + tid ];}
      if(trans == 0)for(int i=0;i<sncut;i++){bs[tid][i] = bp[i*w + tid ];}
      if(trans == 1)for(int i=0;i<sncut;i++){bs[tid][i] = bp[tid*n + i];}
    }

    __syncthreads();

    if(ch == 1){if(threadIdx.y < smcut and threadIdx.x < sncut)
         for(int i = 0; i < wcut; i++){
          if(Conj==0)buf += as[threadIdx.y][i] * bs[i][threadIdx.x];
          if(Conj==1)buf += qlat::qconj(as[threadIdx.y][i]) * bs[i][threadIdx.x];
    }}
    else{for(int i = 0; i < wcut; i++){
          if(Conj==0)buf += as[threadIdx.y][i] * bs[i][threadIdx.x];
          if(Conj==1)buf += qlat::qconj(as[threadIdx.y][i]) * bs[i][threadIdx.x];
    }}

    wini += blockSize;
  }

  if(ch == 1){if(mi < m and ni < n){ c[ci] += buf;}}
  else{                              c[ci] += buf;}


}

template<typename Ty>
void matrix_prod_gpu2(Ty* a, Ty* b, Ty* c, const long m, const long n, const long w, const long L=1, bool Conj=true, bool trans=false)
{
  int nt = 32;int sm = 8;int sn = 4;

  //int ntLm[] = {32,32,16,16,8,8,4,4,2,2,1};
  //int ntLn[] = {32,16,16, 8,8,4,4,2,2,1,1};

  //if(sizeof(Ty) == 8){nt = 128; sm=16;sn=8;}
  //if(sizeof(Ty) == 8){nt = 64;sm=8;sn=8;}
  //int sm=1;int sn=1;

  //for(int i=0;i< 11;i++){
  //  if(nt == ntLm[i]*ntLn[i]){sm=ntLm[i];sn=ntLn[i];break;}
  //}

  int cnt = nt;
  int mc = m/sm;
  int nc = n/sn;

  if(m%sm != 0 or m < sm){cnt = nt+1;mc += 1;}
  if(n%sn != 0 or n < sn){cnt = nt+1;nc += 1;}
  if(w%nt != 0 or w < nt){cnt = nt+1;}
  if(Conj )cnt += 1000;
  if(trans)cnt += 7000;

  //if(m%sm != 0)  qlat::displayln_info(qlat::ssprintf("m Dimension not dividable! "));
  //if(n%sn != 0)  qlat::displayln_info(qlat::ssprintf("n Dimension not dividable! "));
  //if(w%nt != 0)  qlat::displayln_info(qlat::ssprintf("w Dimension not dividable! "));

  //if(m%sm != 0 or n%sn != 0 or w%nt != 0){
  //  qlat::displayln_info(qlat::ssprintf("Dimension not dividable by 4 ! "));
  //  qassert(false);
  //}

  ////////memory L --> m --> n
  dim3 dimGrid( nc, mc, L); 
  dim3 dimBlock(sn, sm, 1); 

  bool jobdo = true;
  switch (cnt)
  {
    //case 1024:matrix_prod_global2<1024,32,32,0, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    //case 512: matrix_prod_global2< 512,32,16,0, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    //case 256: matrix_prod_global2< 256,16,16,0, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    //case 4:   matrix_prod_global2<   4, 2, 2,0, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    //case 2:   matrix_prod_global2<   2, 2, 1,0, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    //case 1:   matrix_prod_global2<   1, 1, 1,0, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    ////#if Enablefloat == 1
    //if(sizeof(Ty) == 8){
    //case 128: matrix_prod_global2< 128,16, 8,0, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    //case 64:  matrix_prod_global2<  64, 8, 8,0, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    //case 129: matrix_prod_global2< 128,16, 8,1, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    //case 65:  matrix_prod_global2<  64, 8, 8,1, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    //}
    ////#endif

    case   32:  matrix_prod_global2<  32, 8, 4,0,0,0, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    case   16:  matrix_prod_global2<  16, 4, 4,0,0,0, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    case    8:  matrix_prod_global2<   8, 4, 2,0,0,0, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;

    case   33:  matrix_prod_global2<  32, 8, 4,1,0,0, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    case   17:  matrix_prod_global2<  16, 4, 4,1,0,0, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    case    9:  matrix_prod_global2<   8, 4, 2,1,0,0, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;

    case 1032:  matrix_prod_global2<  32, 8, 4,0,1,0, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    case 1016:  matrix_prod_global2<  16, 4, 4,0,1,0, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    case 1008:  matrix_prod_global2<   8, 4, 2,0,1,0, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;

    case 1033:  matrix_prod_global2<  32, 8, 4,1,1,0, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    case 1017:  matrix_prod_global2<  16, 4, 4,1,1,0, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    case 1009:  matrix_prod_global2<   8, 4, 2,1,1,0, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;


    case 7032:  matrix_prod_global2<  32, 8, 4,0,0,1, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    case 7016:  matrix_prod_global2<  16, 4, 4,0,0,1, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    case 7008:  matrix_prod_global2<   8, 4, 2,0,0,1, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;

    case 7033:  matrix_prod_global2<  32, 8, 4,1,0,1, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    case 7017:  matrix_prod_global2<  16, 4, 4,1,0,1, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    case 7009:  matrix_prod_global2<   8, 4, 2,1,0,1, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;

    case 8032:  matrix_prod_global2<  32, 8, 4,0,1,1, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    case 8016:  matrix_prod_global2<  16, 4, 4,0,1,1, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    case 8008:  matrix_prod_global2<   8, 4, 2,0,1,1, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;

    case 8033:  matrix_prod_global2<  32, 8, 4,1,1,1, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    case 8017:  matrix_prod_global2<  16, 4, 4,1,1,1, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    case 8009:  matrix_prod_global2<   8, 4, 2,1,1,1, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;

    jobdo = false;
  }
  qassert(jobdo);

}


template <unsigned int sn, unsigned int ch,unsigned int Conj, typename Ty >
__global__ void matrix_prod_global1(Ty* a, Ty* b, Ty* c, const long m, const long n, const long w)
{
  __shared__ Ty as[sn+1][sn+1];
  __shared__ Ty bs[sn+1][sn+1];

  ////////======c, mxn, a, mxw, b, nxw
  /////unsigned int  tid= threadIdx.y*blockDim.x + threadIdx.x;

  unsigned long ni =  blockIdx.x * blockDim.x + threadIdx.x;
  unsigned long mi =  blockIdx.y * blockDim.y + threadIdx.y;
  unsigned long li =  blockIdx.z;

  unsigned long ci =  (li*m + mi)*n + ni;

  unsigned long wini = 0*sn;
  unsigned int  wcut =   sn;
  Ty buf;buf = 0.0;

  while(wini < w){
    if(ch == 1){
      if(wini + sn > w){wcut = w - wini;}
      if(mi < m and threadIdx.x < wcut){
            as[threadIdx.y][threadIdx.x] = a[li*m*w + (mi)*w + wini + threadIdx.x];
      }else{as[threadIdx.y][threadIdx.x]=0.0;}
      if(ni < n and threadIdx.y < wcut){
            bs[threadIdx.y][threadIdx.x] = b[li*n*w + (ni)*w + wini + threadIdx.y];
      }else{bs[threadIdx.y][threadIdx.x]=0.0;}
    }else{
      as[threadIdx.y][threadIdx.x] = a[li*m*w + (mi)*w + wini + threadIdx.x];
      bs[threadIdx.y][threadIdx.x] = b[li*n*w + (ni)*w + wini + threadIdx.y];
    }
    __syncthreads();

    for(int i = 0; i < sn; i++)
    {
      //////buf += as[threadIdx.y][i] * bs[i][threadIdx.x];
      if(Conj==0){buf += as[threadIdx.y][i] * bs[i][threadIdx.x];}
      if(Conj==1){buf += qlat::qconj(as[threadIdx.y][i]) * bs[i][threadIdx.x];}
    }
    wini += sn;
  }
  if(ch == 1){if(mi < m and ni < n){ c[ci] += buf;}}
  else{                              c[ci] += buf;}


}

template<typename Ty>
void matrix_prod_gpu1(Ty* a, Ty* b, Ty* c, const long m, const long n, const long w, const long L=1, bool Conj=true)
{
  int sn = 4;
  int nt = sn*sn;

  ///if(m%sn != 0 or n%sn != 0 or w%sn != 0){
  ///  qlat::displayln_info(qlat::ssprintf("Dimension not dividable by 4 ! "));
  ///  qassert(false);
  ///}

  int nc = n/sn ;
  int mc = m/sn ;

  int csn = sn;
  if(m%sn != 0 or m<sn){csn = sn+2000;mc += 1;}
  if(n%sn != 0 or n<sn){csn = sn+2000;nc += 1;}
  if(w%nt != 0 or w<nt){csn = sn+2000;}
  if(Conj)csn += 3000;
  //if(m%sn != 0)  qlat::displayln_info(qlat::ssprintf("m Dimension not dividable! "));
  //if(n%sn != 0)  qlat::displayln_info(qlat::ssprintf("n Dimension not dividable! "));
  //if(w%nt != 0)  qlat::displayln_info(qlat::ssprintf("w Dimension not dividable! "));


  ////////memory L --> m --> n
  dim3 dimGrid( nc, mc, L); 
  dim3 dimBlock(sn, sn, 1); 

  bool jobdo = true;
  switch (csn)
  {
    //case   32:  matrix_prod_global1< 32,0,0, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    //case   16:  matrix_prod_global1< 16,0,0, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    case    8:  matrix_prod_global1<  8,0,0, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    case    4:  matrix_prod_global1<  4,0,0, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    case    2:  matrix_prod_global1<  2,0,0, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    case    1:  matrix_prod_global1<  1,0,0, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    //case 2032:  matrix_prod_global1< 32,1,0, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    //case 2016:  matrix_prod_global1< 16,1,0, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    case 2008:  matrix_prod_global1<  8,1,0, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    case 2004:  matrix_prod_global1<  4,1,0, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    case 2002:  matrix_prod_global1<  2,1,0, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    case 2001:  matrix_prod_global1<  1,1,0, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;

    //case 3032:  matrix_prod_global1< 32,0,1, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    //case 3016:  matrix_prod_global1< 16,0,1, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    case 3008:  matrix_prod_global1<  8,0,1, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    case 3004:  matrix_prod_global1<  4,0,1, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    case 3002:  matrix_prod_global1<  2,0,1, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    case 3001:  matrix_prod_global1<  1,0,1, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    //case 5032:  matrix_prod_global1< 32,1,1, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    //case 5016:  matrix_prod_global1< 16,1,1, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    case 5008:  matrix_prod_global1<  8,1,1, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    case 5004:  matrix_prod_global1<  4,1,1, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    case 5002:  matrix_prod_global1<  2,1,1, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    case 5001:  matrix_prod_global1<  1,1,1, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;

    jobdo = false;
  }
  qassert(jobdo);

}


template <unsigned int Conj,typename Ty>
__global__ void matrix_prod_global0(Ty* a, Ty* b, Ty* c, const long m, const long n, const long w)
{
  //////======c, mxn, a, mxw, b, nxw
  unsigned long ni =  blockIdx.x * blockDim.x + threadIdx.x;
  unsigned long mi =  blockIdx.y * blockDim.y + threadIdx.y;
  unsigned long li =  blockIdx.z;

  unsigned long ci =  (li*m + mi)*n + ni;
  Ty* ap;
  Ty* bp;
  ap = &a[li*m*w + mi*w];
  bp = &b[li*n*w + ni*w];

  if(ni < n and mi < m){
    Ty buf;buf = 0.0;
    unsigned long wi = 0;
    while(wi < w){
      if(Conj==0){buf +=             ap[wi] * bp[wi];}
      if(Conj==1){buf += qlat::qconj(ap[wi]) * bp[wi];}
      wi += 1;
    }
    c[ci] += buf;
  }

}


template<typename Ty>
void matrix_prod_gpu0(Ty* a, Ty* b, Ty* c, const long m, const long n, const long w, const long L=1, bool Conj=true)
{
  //case 1024:reduce6<1024,Ty><<< dimGrid, dimBlock >>>(src, res, n, divide); break;
  int ntLm[] = {32,32,16,16,8,8,4,4,2,2,1};
  int ntLn[] = {32,16,16, 8,8,4,4,2,2,1,1};

  int nt = 64;
  int sm=1;int sn=1;
  for(int i=0;i< 11;i++){
    if(nt == ntLm[i]*ntLn[i]){sm=ntLm[i];sn=ntLn[i];break;}
  }

  int nc = n/sn + 1;
  int mc = m/sm + 1;
  ////////memory L --> m --> n
  dim3 dimGrid( nc, mc, L); 
  dim3 dimBlock(sn, sm, 1); 
  int snt = 0;if(Conj){snt += 1;}

  bool jobdo = true;
  switch (snt)
  {
    case 0:matrix_prod_global0< 0, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    case 1:matrix_prod_global0< 1, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    jobdo = false;
  }
  qassert(jobdo);

}
#endif

template<typename Ty>
void matrix_prod_cpu(Ty* a, Ty* b, Ty* c, const long m, const long n, const long w, const long L=1, bool Conj=true, bool trans=false)
{
  TIMER_FLOPS("==Matrix Multi CPU");
  size_t offA = m*w;
  size_t offB = n*w;
  size_t offC = m*n;

  //#pragma omp parallel for
  //for(long li=0;li<L;li++)
  //{
  //  EML  A(&a[li*offA], m, w);
  //  EMLC B(&b[li*offB], w, n);
  //  EML  C(&c[li*offC], m, n);
  //  C += A * B;
  //}

  //Eigen::initParallel();
  //Eigen::setNbThreads(0);

  //for(long li=0;li<L;li++)
  //{
  //  EM  A( &a[li*offA], m, w);
  //  EMC B( &b[li*offB], w, n);
  //  EM  C( &c[li*offC], m, n);
  //    if(Conj){C += A.conjugate() * B;}else{C += A * B;}
  //}

  int Nv = omp_get_max_threads();
  #ifdef QLAT_USE_ACC
  Nv = 1;
  #endif
  ////Nv = 1;
  if(m*L < Nv){Nv = 1;}

  if(Nv == 1){
  for(long li=0;li<L;li++)
  {
    EML  A( &a[li*offA], m, w);
    EMLC B0( &b[li*offB], w, n);
    EML  B1( &b[li*offB], w, n);
    EML  C( &c[li*offC], m, n);
    if(trans){if(Conj){C += A.conjugate() * B1;}else{C += A * B1;}}
    else{     if(Conj){C += A.conjugate() * B0;}else{C += A * B0;}}
  }}else{
    //Eigen::setNbThreads(1);
    //int Nm = get_threads_GPU(m,Nv);
    //int Nfac = m/Nm;
    int Nfac = Nv;
    int Nm = (m+Nv-1)/Nv;
    #pragma omp parallel for
    for(long off=0;off<L*Nfac;off++)
    {
      int li   = off/Nfac;
      int mi   = off%Nfac;
      int mcut = Nm;if(mi*Nm + mcut > m){mcut = m - mi*Nm;if(mcut <=0){continue;}}

      /////if(mi*Nm + Nm > m){mcut = m - mi*Nm;}
      EML  A(&a[li*offA+(mi*Nm)*w], mcut, w);
      /////EMLC B(&b[li*offB], w, n);
      EMLC B0(&b[li*offB], w, n);
      EML  B1(&b[li*offB], w, n);
      EML  C(&c[li*offC+(mi*Nm)*n], mcut, n);
      if(trans){if(Conj){C += A.conjugate() * B1;}else{C += A * B1;}}
      else{     if(Conj){C += A.conjugate() * B0;}else{C += A * B0;}}
    }
    //Eigen::setNbThreads(0);
  }


  long long vGb = L*m*n*w;
  int Fcount0   = 6 + 2; 
  timer.flops  += vGb*Fcount0;

  //double Gsize = (m*n + m*w + n*w)*sizeof(Complexq)/(1024.0*1024*1024);
  /////qlat::displayln_info(qlat::ssprintf("Total memory size %.3e GB \n",Gsize));
}

template<typename Ty>
void matrix_prod_gpu(Ty* a, Ty* b, Ty* c, const long m, const long n, const long w, const long L=1, bool Conj=true, bool dummy = true, bool trans=false, int modeGPU = 2)
{
  #ifdef QLAT_USE_ACC

  TIMER_FLOPS("==Matrix Multi GPU");
  long long vGb = L*m*n*w;
  int Fcount0   = 6 + 2; 
  timer.flops += vGb*Fcount0;

  /////int modeGPU = 1;
  if(modeGPU == 0){matrix_prod_gpu0(a, b, c, m, n, w, L, Conj);}
  if(modeGPU == 1){matrix_prod_gpu1(a, b, c, m, n, w, L, Conj);}
  ///////Trans only works for modeGPU 2 
  if(modeGPU == 2){matrix_prod_gpu2(a, b, c, m, n, w, L, Conj, trans);}

  if(dummy)qacc_barrier(dummy);

  #else
  ////matrix_prod_gpu(bool Conj=true, bool dummy = true, bool trans=false, int modeGPU = 2)
  matrix_prod_cpu(a,b,c , m,n,w,L, Conj, trans);
  #endif
}


////// c = a b; c dim m x n, a dim m x w, b dim w x n, done it l times
template<typename Ty>
void matrix_prod(Ty* a, Ty* b, Ty* c, const long m, const long n, const long w, const long L=1, bool Conj=true)
{

  #ifdef QLAT_USE_ACC
  matrix_prod_gpu(a,b,c , m,n,w,L, Conj);
  #else
  matrix_prod_cpu(a,b,c , m,n,w,L, Conj);
  #endif

  /////double Gsize = (m*n + m*w + n*w)*sizeof(Ty)/1024.0*1024*1024;
  /////qlat::displayln_info(qlat::ssprintf("Total memory size, GFlop %.3e \n",Gsize));

}


}

#undef EML
#undef EMLC

#endif
