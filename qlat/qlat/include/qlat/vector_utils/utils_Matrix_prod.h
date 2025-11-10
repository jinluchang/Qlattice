// utils_GPU_Matrix_prod.h
// Gen Wang
// Jun. 2021
#ifndef UTILS_MATRIX_PROD_H
#define UTILS_MATRIX_PROD_H

#pragma once
#include "utils_vector_GPU.h"

#define EML  Eigen::Map< Eigen::Matrix<Ty , Eigen::Dynamic, Eigen::Dynamic ,Eigen::RowMajor> >
#define EMLC Eigen::Map< Eigen::Matrix<Ty , Eigen::Dynamic, Eigen::Dynamic ,Eigen::ColMajor> >

//////  https://docs.nvidia.com/cuda/pdf/CUDA_C_Best_Practices_Guide.pdf

namespace qlat{
////need Conjugate in first element

#ifdef QLAT_USE_ACC
//////shared memory 16KB or 48KB per thread
template <unsigned int blockSize, unsigned int sm, unsigned int sn, unsigned int ch,unsigned int Conj, unsigned int trans, typename Ty >
__global__ void matrix_prod_global2(Ty** a, Ty** b, Ty** c, const Long m, const Long n, const Long w)
{
  __shared__ Ty as[sm+1][blockSize+1];
  __shared__ Ty bs[blockSize+1][sn+1];

  ////////======c, mxn, a, mxw, b, nxw
  unsigned int  tid=  threadIdx.y*blockDim.x + threadIdx.x;

  unsigned long ni =  blockIdx.x * blockDim.x + threadIdx.x;
  unsigned long mi =  blockIdx.y * blockDim.y + threadIdx.y;
  unsigned long li =  blockIdx.z;

  unsigned long ci =  (mi)*n + ni;
  Ty* ap;
  Ty* bp;

  Ty buf;buf = 0.0;

  //unsigned long wini = 0*blockSize;
  //while(wini < w){
  //  ap = &a[(li*m + blockIdx.y * blockDim.y + 0)*w + wini];
  //  bp = &b[(li*n + blockIdx.x * blockDim.x + 0)*w + wini];

  //  for(Int i=0;i<sm;i++){as[i][tid] = ap[i*w + tid ];}

  //  for(Int i=0;i<sn;i++){bs[tid][i] = bp[i*w + tid ];}

  //  __syncthreads();

  //  for(Int i = 0; i < blockSize; i++)
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
  Ty* al = a[li];
  Ty* bl = b[li];
  Ty* cl = c[li];

  if(ch==1)if(blockIdx.y * blockDim.y + sm > m){smcut = m - blockIdx.y * blockDim.y;}
  if(ch==1)if(blockIdx.x * blockDim.x + sn > n){sncut = n - blockIdx.x * blockDim.x;}

  while(wini < w){
    if(ch == 1)if(wini + blockSize > w){wcut = w - wini;}

    ap = &al[(blockIdx.y * blockDim.y + 0)*w + wini];
    if(trans == 0){bp = &bl[(blockIdx.x * blockDim.x + 0)*w + wini];}
    if(trans == 1){bp = &bl[(wini)*n + blockIdx.x * blockDim.x + 0];}

    if(tid < wcut){
      for(Int i=0;i<smcut;i++){as[i][tid] = ap[i*w + tid ];}
      if(trans == 0)for(Int i=0;i<sncut;i++){bs[tid][i] = bp[i*w + tid ];}
      if(trans == 1)for(Int i=0;i<sncut;i++){bs[tid][i] = bp[tid*n + i];}
    }

    __syncthreads();

    if(ch == 1){if(threadIdx.y < smcut and threadIdx.x < sncut)
         for(Int i = 0; i < wcut; i++){
          if(Conj==0)buf += as[threadIdx.y][i] * bs[i][threadIdx.x];
          if(Conj==1)buf += qlat::qconj(as[threadIdx.y][i]) * bs[i][threadIdx.x];
    }}
    else{for(Int i = 0; i < wcut; i++){
          if(Conj==0)buf += as[threadIdx.y][i] * bs[i][threadIdx.x];
          if(Conj==1)buf += qlat::qconj(as[threadIdx.y][i]) * bs[i][threadIdx.x];
    }}

    wini += blockSize;
  }

  if(ch == 1){if(mi < m and ni < n){ cl[ci] += buf;}}
  else{                              cl[ci] += buf;}
}

template<typename Ty>
void matrix_prod_gpu2(Ty** a, Ty** b, Ty** c, const Long m, const Long n, const Long w, const Long L=1, bool Conj=true, bool trans=false)
{
  Int nt = 32;int sm = 8;int sn = 4;

  //int ntLm[] = {32,32,16,16,8,8,4,4,2,2,1};
  //int ntLn[] = {32,16,16, 8,8,4,4,2,2,1,1};

  //if(sizeof(Ty) == 8){nt = 128; sm=16;sn=8;}
  //if(sizeof(Ty) == 8){nt = 64;sm=8;sn=8;}
  //int sm=1;int sn=1;

  //for(Int i=0;i< 11;i++){
  //  if(nt == ntLm[i]*ntLn[i]){sm=ntLm[i];sn=ntLn[i];break;}
  //}

  Int cnt = nt;
  Int mc = m/sm;
  Int nc = n/sn;

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
  //  Qassert(false);
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
  Qassert(jobdo);

}


template <unsigned int sn, unsigned int ch,unsigned int Conj, typename Ty >
__global__ void matrix_prod_global1(Ty** a, Ty** b, Ty** c, const Long m, const Long n, const Long w)
{
  __shared__ Ty as[sn+1][sn+1];
  __shared__ Ty bs[sn+1][sn+1];

  ////////======c, mxn, a, mxw, b, nxw
  /////unsigned int  tid= threadIdx.y*blockDim.x + threadIdx.x;

  unsigned long ni =  blockIdx.x * blockDim.x + threadIdx.x;
  unsigned long mi =  blockIdx.y * blockDim.y + threadIdx.y;
  unsigned long li =  blockIdx.z;

  //unsigned long ci =  (li*m + mi)*n + ni;
  unsigned long ci =  (mi)*n + ni;

  unsigned long wini = 0*sn;
  unsigned int  wcut =   sn;
  Ty buf;buf = 0.0;

  while(wini < w){
    if(ch == 1){
      if(wini + sn > w){wcut = w - wini;}
      if(mi < m and threadIdx.x < wcut){
            as[threadIdx.y][threadIdx.x] = a[li][ (mi)*w + wini + threadIdx.x];
      }else{as[threadIdx.y][threadIdx.x]=0.0;}
      if(ni < n and threadIdx.y < wcut){
            bs[threadIdx.y][threadIdx.x] = b[li][ (ni)*w + wini + threadIdx.y];
      }else{bs[threadIdx.y][threadIdx.x]=0.0;}
    }else{
      as[threadIdx.y][threadIdx.x] = a[li][ (mi)*w + wini + threadIdx.x];
      bs[threadIdx.y][threadIdx.x] = b[li][ (ni)*w + wini + threadIdx.y];
    }
    __syncthreads();

    for(Int i = 0; i < sn; i++)
    {
      //////buf += as[threadIdx.y][i] * bs[i][threadIdx.x];
      if(Conj==0){buf += as[threadIdx.y][i] * bs[i][threadIdx.x];}
      if(Conj==1){buf += qlat::qconj(as[threadIdx.y][i]) * bs[i][threadIdx.x];}
    }
    wini += sn;
  }
  if(ch == 1){if(mi < m and ni < n){ c[li][ci] += buf;}}
  else{                              c[li][ci] += buf;}


}

template<typename Ty>
void matrix_prod_gpu1(Ty** a, Ty** b, Ty** c, const Long m, const Long n, const Long w, const Long L=1, bool Conj=true)
{
  Int sn = 4;
  Int nt = sn*sn;

  ///if(m%sn != 0 or n%sn != 0 or w%sn != 0){
  ///  qlat::displayln_info(qlat::ssprintf("Dimension not dividable by 4 ! "));
  ///  Qassert(false);
  ///}

  Int nc = n/sn ;
  Int mc = m/sn ;

  Int csn = sn;
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
  Qassert(jobdo);

}


template <unsigned int Conj,typename Ty>
__global__ void matrix_prod_global0(Ty** a, Ty** b, Ty** c, const Long m, const Long n, const Long w)
{
  //////======c, mxn, a, mxw, b, nxw
  unsigned long ni =  blockIdx.x * blockDim.x + threadIdx.x;
  unsigned long mi =  blockIdx.y * blockDim.y + threadIdx.y;
  unsigned long li =  blockIdx.z;

  unsigned long ci =  ( mi)*n + ni;
  Ty* ap;
  Ty* bp;
  //ap = &a[li*m*w + mi*w];
  //bp = &b[li*n*w + ni*w];

  ap = &a[li][ mi*w];
  bp = &b[li][ ni*w];

  if(ni < n and mi < m){
    Ty buf;buf = 0.0;
    unsigned long wi = 0;
    while(wi < w){
      if(Conj==0){buf +=             ap[wi] * bp[wi];}
      if(Conj==1){buf += qlat::qconj(ap[wi]) * bp[wi];}
      wi += 1;
    }
    c[li][ci] += buf;
  }

}

template<typename Ty>
void matrix_prod_gpu0(Ty** a, Ty** b, Ty** c, const Long m, const Long n, const Long w, const Long L=1, bool Conj=true)
{
  //case 1024:reduce6<1024,Ty><<< dimGrid, dimBlock >>>(src, res, n, divide); break;
  Int ntLm[] = {32,32,16,16,8,8,4,4,2,2,1};
  Int ntLn[] = {32,16,16, 8,8,4,4,2,2,1,1};

  Int nt = 64;
  Int sm=1;int sn=1;
  for(Int i=0;i< 11;i++){
    if(nt == ntLm[i]*ntLn[i]){sm=ntLm[i];sn=ntLn[i];break;}
  }

  Int nc = n/sn + 1;
  Int mc = m/sm + 1;
  ////////memory L --> m --> n
  dim3 dimGrid( nc, mc, L); 
  dim3 dimBlock(sn, sm, 1); 
  Int snt = 0;if(Conj){snt += 1;}

  bool jobdo = true;
  switch (snt)
  {
    case 0:matrix_prod_global0< 0, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    case 1:matrix_prod_global0< 1, Ty ><<< dimGrid, dimBlock >>>(a, b, c, m, n, w); break;
    jobdo = false;
  }
  Qassert(jobdo);

}
#endif

template<typename Ty>
void matrix_prod_cpu(Ty** a, Ty** b, Ty** c, const Long m, const Long n, const Long w, const Long L=1, bool Conj=true, bool trans=false)
{
  TIMER_FLOPS("==Matrix Multi CPU");
  //size_t offA = m*w;
  //size_t offB = n*w;
  //size_t offC = m*n;

  //#pragma omp parallel for
  //for(Long li=0;li<L;li++)
  //{
  //  EML  A(&a[li*offA], m, w);
  //  EMLC B(&b[li*offB], w, n);
  //  EML  C(&c[li*offC], m, n);
  //  C += A * B;
  //}

  //Eigen::initParallel();
  //Eigen::setNbThreads(0);

  //for(Long li=0;li<L;li++)
  //{
  //  EM  A( &a[li*offA], m, w);
  //  EMC B( &b[li*offB], w, n);
  //  EM  C( &c[li*offC], m, n);
  //    if(Conj){C += A.conjugate() * B;}else{C += A * B;}
  //}

  Int Nv = omp_get_max_threads();
  /////#ifdef QLAT_USE_ACC
  /////Nv = 1;
  /////#endif
  ////////Nv = 1;
  if(m*L < Nv){Nv = 1;}

  if(Nv == 1){
  for(Long li=0;li<L;li++)
  {
    EML  A(  &a[li][0], m, w);
    EMLC B0( &b[li][0], w, n);
    EML  B1( &b[li][0], w, n);
    EML  C(  &c[li][0], m, n);
    if(trans){if(Conj){C += A.conjugate() * B1;}else{C += A * B1;}}
    else{     if(Conj){C += A.conjugate() * B0;}else{C += A * B0;}}
  }}else{
    //Eigen::setNbThreads(1);
    //int Nm = get_threads_GPU(m,Nv);
    //int Nfac = m/Nm;
    Int Nfac = Nv;
    Int Nm = (m+Nv-1)/Nv;
    #pragma omp parallel for
    for(Long off=0;off<L*Nfac;off++)
    {
      Int li   = off/Nfac;
      Int mi   = off%Nfac;
      Int mcut = Nm;if(mi*Nm + mcut > m){mcut = m - mi*Nm;if(mcut <=0){continue;}}

      /////if(mi*Nm + Nm > m){mcut = m - mi*Nm;}
      EML  A( &a[li][(mi*Nm)*w], mcut, w);
      /////EMLC B(&b[li*offB], w, n);
      EMLC B0(&b[li][0], w, n);
      EML  B1(&b[li][0], w, n);
      EML  C( &c[li][(mi*Nm)*n], mcut, n);
      if(trans){if(Conj){C += A.conjugate() * B1;}else{C += A * B1;}}
      else{     if(Conj){C += A.conjugate() * B0;}else{C += A * B0;}}
    }
    //Eigen::setNbThreads(0);
  }


  long long vGb = L*m*n*w;
  Int Fcount0   = 6 + 2; 
  timer.flops  += vGb*Fcount0;

  //double Gsize = (m*n + m*w + n*w)*sizeof(Complexq)/(1024.0*1024*1024);
  /////qlat::displayln_info(qlat::ssprintf("Total memory size %.3e GB \n",Gsize));
}

template<typename Ty>
void matrix_prod_gpu(Ty** a, Ty** b, Ty** c, const Long m, const Long n, const Long w, const Long L=1, bool Conj=true, QBOOL dummy = QTRUE, bool trans=false, Int modeGPU = 2)
{
  (void)dummy;
  (void)modeGPU;
  #ifdef QLAT_USE_ACC

  TIMER_FLOPS("==Matrix Multi GPU");
  long long vGb = L*m*n*w;
  Int Fcount0   = 6 + 2; 
  timer.flops += vGb*Fcount0;

  /////int modeGPU = 1;
  if(modeGPU == 0){matrix_prod_gpu0(a, b, c, m, n, w, L, Conj);}
  if(modeGPU == 1){matrix_prod_gpu1(a, b, c, m, n, w, L, Conj);}
  ///////Trans only works for modeGPU 2 
  if(modeGPU == 2){matrix_prod_gpu2(a, b, c, m, n, w, L, Conj, trans);}

  if(dummy == QTRUE)qacc_barrier(dummy);

  #else
  ////matrix_prod_gpu(bool Conj=true, bool dummy = true, bool trans=false, Int modeGPU = 2)
  matrix_prod_cpu(a,b,c , m,n,w,L, Conj, trans);
  #endif
}

/* 
  C = A B; C dim m x n, A dim m x w, B dim w x n, do it l times
  mem A --> [L][mi * w + wi] if Conj, a[i] -> a^*[i]
  mem B --> [L][ni * w + wi] if tran [wi * n + ni]
  mem C --> [L][mi * n + ni]
*/
template<typename Ty>
void matrix_prodPT(Ty** a, Ty** b, Ty** c, const Long m, const Long n, const Long w, const Long L=1, bool Conj=true, bool trans=false, bool GPU = true, QBOOL dummy = QTRUE)
{
  if( GPU){matrix_prod_gpu(a,b,c , m,n,w,L, Conj, dummy, trans);}
  if(!GPU){matrix_prod_cpu(a,b,c , m,n,w,L, Conj, trans);}
}

template<typename Ty>
void matrix_prodT(Ty* A, Ty* B, Ty* C, const Long m, const Long n, const Long w, const Long L=1, bool Conj=true, bool trans=false, bool GPU = true, QBOOL dummy = QTRUE)
{
  const size_t offA = m*w;
  const size_t offB = n*w;
  const size_t offC = m*n;
  qlat::vector_gpu<Ty* > aP;aP.resize(L, GPU);Ty** a = (Ty**) aP.data();
  qlat::vector_gpu<Ty* > bP;bP.resize(L, GPU);Ty** b = (Ty**) bP.data();
  qlat::vector_gpu<Ty* > cP;cP.resize(L, GPU);Ty** c = (Ty**) cP.data();
  qGPU_for(i, L, GPU, {
    a[i] = &A[i*offA + 0];
    b[i] = &B[i*offB + 0];
    c[i] = &C[i*offC + 0];
  });
  if( GPU){matrix_prod_gpu(a, b, c, m,n,w,L, Conj, dummy, trans);}
  if(!GPU){matrix_prod_cpu(a, b, c, m,n,w,L, Conj, trans);}
  if(dummy == QFALSE){abort_r("buffers of pointers undefined if job not complete!");}
  if(dummy == QTRUE ){qacc_barrier(dummy);}
}

void matrix_prodP(qlat::ComplexT<double>** a, qlat::ComplexT<double>** b, qlat::ComplexT<double>** c, const Long m, const Long n, const Long w, const Long L=1, bool Conj=true, bool trans=false, bool GPU = true, QBOOL dummy = QTRUE);
void matrix_prodP(qlat::ComplexT<float >** a, qlat::ComplexT<float >** b, qlat::ComplexT<float >** c, const Long m, const Long n, const Long w, const Long L=1, bool Conj=true, bool trans=false, bool GPU = true, QBOOL dummy = QTRUE);

void matrix_prod(qlat::ComplexT<double>* A, qlat::ComplexT<double>* B, qlat::ComplexT<double>* C, const Long m, const Long n, const Long w, const Long L=1, bool Conj=true, bool trans=false, bool GPU = true, QBOOL dummy = QTRUE);
void matrix_prod(qlat::ComplexT<float >* A, qlat::ComplexT<float >* B, qlat::ComplexT<float >* C, const Long m, const Long n, const Long w, const Long L=1, bool Conj=true, bool trans=false, bool GPU = true, QBOOL dummy = QTRUE);


//#ifdef QLAT_INSTANTIATE_MATRIX_PROD
//#define QLAT_EXTERN
//#else
//#define QLAT_EXTERN extern
//#endif
//
//QLAT_EXTERN template void prop_smear_qlat_convension<Real, Real>(
//    Propagator4d&, const GaugeField&, const double, const int,
//    const CoordinateD&, const bool, const int, const int, const int);
//
//void matrix_prod<double>(qlat::ComplexT<double>* , qlat::ComplexT<double>* , qlat::ComplexT<double>* ,
//  const , const , const , const , bool , bool , bool , QBOOL );
//
//
//
//#undef QLAT_EXTERN


}

#undef EML
#undef EMLC

#endif
