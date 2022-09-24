#ifndef utils_low_rho_h
#define utils_low_rho_h
#pragma once


//#include <qlat/reduce_V.h>
///#include "cach_reduce.h"
#include "utils_reduce_vec.h"
#include "utils_gammas.h"
#include "utils_Matrix_prod.h"
#include "utils_fft_desc.h"

#define EigenVq Eigen::Matrix< Complexq, 1, Eigen::Dynamic ,Eigen::RowMajor>
//////#define EigenAq Eigen::Array< Complexq, Eigen::Dynamic , 1>
#define Aoper  16

namespace qlat{

#ifdef QLAT_USE_ACC
//__device__ __constant__  signed char  Gmap0C[32];
//__device__ __constant__  signed char  Gmap1C[32];
//__device__ __constant__  signed char  Gmap2C[32];
//__device__ __constant__  signed char  Gmap3C[32];

__global__ void multiplyNab_global(const Complexq* Nab, Ftype *Mres,const Ftype *Mvalues,const int nt,const int nmass,const unsigned long bufN0, signed char* GmapM)
{
  __shared__ signed char G0[32];
  __shared__ signed char G1[32];
  __shared__ signed char G2[32];
  __shared__ signed char G3[32];
  extern __shared__ Complexq NMv_values[];
  /////__shared__ Complexq Nabv_multi[16*3];

  unsigned int tid = threadIdx.x;
  long bi = blockIdx.x;
  //long bN = gridDim.x;
  ///unsigned int bi = ji/Aoper;
  ///unsigned int ipr  = ji%Aoper;

  ////Load gammas
  int off = tid;
  //          while(off < 32){G0[off] = Gmap0C[off];off += blockDim.x;}
  //off = tid;while(off < 32){G1[off] = Gmap1C[off];off += blockDim.x;}
  //off = tid;while(off < 32){G2[off] = Gmap2C[off];off += blockDim.x;}
  //off = tid;while(off < 32){G3[off] = Gmap3C[off];off += blockDim.x;}
            while(off < 32){G0[off] = GmapM[0*32+off];off += blockDim.x;}
  off = tid;while(off < 32){G1[off] = GmapM[1*32+off];off += blockDim.x;}
  off = tid;while(off < 32){G2[off] = GmapM[2*32+off];off += blockDim.x;}
  off = tid;while(off < 32){G3[off] = GmapM[3*32+off];off += blockDim.x;}

  //off = tid;unsigned long offAB = bi*nt*16;
  //while(off < nt*16){NMv_values[off] = Nab[offAB + off];off += blockDim.x;}
  ////off = tid;
  unsigned long offAB = bufN0*16;/////off = bi*16;
  /////if(tid < 16){for(int ti=0;ti<nt;ti++)NMv_values[ti*16 + tid] = Nab[ti*offAB + off + tid];}
  off = tid;
  while(off < nt*16){
    int ti = off/16;
    int op = off%16;
    NMv_values[ti*16 + op] = Nab[ti*offAB + bi*16 + op];
    off += blockDim.x;
  }
  ///off = tid;int offT = nt*16;
  ///while(off<nmass*2){NMv_values[offT + off] = Mvalues[off];off += blockDim.x;}
  __syncthreads();

  const Complexq *Av,*Bv;
  Complexq v0[16];
  Complexq v1[16];
  ////for(int ipr=0;ipr<16*nt;ipr++){v1[ipr] = 0.0;}

  off = tid;offAB = bi*nmass*2;
  Ftype *buf = (Ftype*) &NMv_values[nt*16+0];
  Ftype *src = (Ftype*) &Mvalues[offAB];
  while(off<nmass*2){buf[off] = src[off];off += blockDim.x;}
  __syncthreads();

  unsigned int toff = tid;
  unsigned long shiftM = bufN0*nt*nt;
  /////unsigned int t0 = tid;
  for(unsigned int t0=0;t0<nt;t0++)
  {
    ////for(int ipr=0;ipr<16;ipr++){v0[ipr] = 0.0;v1[ipr] = 0.0;}
    unsigned int t1 = (t0+toff)%nt;
    ///Bv = &Avs[0];
    Bv = &NMv_values[t0*16 + 0];
    Av = &NMv_values[t1*16 + 0];

    for(int ipr=0;ipr<16;ipr++){
      v0[ipr] = Av[G0[ipr*2+0]]*qlat::qconj(Bv[G1[ipr*2+0]])
       *Ftype(G0[ipr*2+1])*Ftype(         G1[ipr*2+1]);
      v1[ipr] = Av[G2[ipr*2+0]]*qlat::qconj(Bv[G3[ipr*2+0]])
       *Ftype(G2[ipr*2+1])*Ftype(         G3[ipr*2+1]);
    }

    offAB = (((0*nmass+0)*bufN0 + bi)*nt+t0)*nt + toff;
    for(int ipr=0;ipr<16;ipr++)
    {
      Ftype v00 = v0[ipr].real();
      Ftype v10 = v1[ipr].real();
      for(int mi=0;mi<nmass;mi++)
      {
        ////long jobN = Aoper*nt;
        Mres[offAB] += (buf[mi*2+0]*v00 + buf[mi*2+1]*v10);
        offAB += shiftM;
        ///mi += blockDim.x;
      }
    }
  }

}

__global__ void prodab_global(const Complexq *a,const Complexq *b, Complexq *fd,const int Nvol,const int Nsum)
{
  Complexq as[12];
  Complexq bs[12];
  //extern __shared__ Complexq resab[];
  __shared__ Complexq resab[32*16];
  unsigned int tid = threadIdx.x;
  unsigned long isp = blockIdx.x*blockDim.x + threadIdx.x;

  if(isp < Nvol){
    for(int dc=0;dc<12;dc++){as[dc] = a[dc*Nvol+isp];}
    for(int dc=0;dc<12;dc++){bs[dc] = b[dc*Nvol+isp];}

    {
    for(int bi=0;bi<4;bi++)
    {
      int iv = bi*4 + 0;
      for(int ai=0;ai<4;ai++)
      {
        //#ifndef __HIP_PLATFORM_HCC__
        //Eigen::Map<const EigenVq > aM(&as[ai*3+0],3);
        //Eigen::Map<const EigenVq > bM(&bs[bi*3+0],3);
        //resab[tid*16 + iv] =  bM.dot(aM);
        //#else

        resab[tid*16 + iv] = 0; 
        for(int doti=0;doti<3;doti++){resab[tid*16 + iv] += qlat::qconj(bs[bi*3+doti]) * as[ai*3+doti];}

        //#endif
        iv += 1;
      }
    }
    }

  }else{
    for(int iv=0;iv<16;iv++){resab[tid*16 + iv] = 0.0;}
  }
  __syncthreads();

  ///resab --> isp --> 16 --> reduce by a factor of 4/8/16
  ////Assume thread number 32
  if(tid<16){for(int is= 1;is<16;is++){resab[ 0*16 +tid   ] += resab[is*16+tid   ];}}
  else{      for(int is=17;is<32;is++){resab[16*16 +tid-16] += resab[is*16+tid-16];}}
  __syncthreads();

  if(tid < 16){resab[ 0*16 + tid] += resab[16*16 + tid];}
  __syncthreads();

  const int it  = isp/Nsum;
  ////blockIdx.x*blockDim.x
  const unsigned long offv = Nsum/32;
  const unsigned long off0 = it*16*offv + blockIdx.x%offv;
  Complexq *f0 = &fd[off0];

  ////if(tid==0){for(int iv=0;iv<16;iv++)f0[iv*offv] = resab[iv];}
  if(tid < 16){f0[tid*offv] = resab[tid];}

}
#endif


void prodab(Complexq* a0,Complexq* b0, const qlat::Geometry &geo, Complexq *fM, int mode_reduce=1)
{
  ////const qlat::Geometry &geo = a0.geo();
  ////const Coordinate vg = geo.total_site();
  ////int nt = vg[3];

  unsigned long Nvol = geo.local_volume();
  int Nt = geo.node_site[3];
  long Nsum = Nvol/Nt;
  /////qassert(Nsum%32 == 0);

  Complexq* a = a0;
  Complexq* b = b0;
  //Complexq* a = (Complexq* ) &(a0.get_elem(0));
  //Complexq* b = (Complexq* ) &(b0.get_elem(0));
  Complexq *f0 = &fM[0];

  //////GPU version and function return
  if(mode_reduce == 1){
    #ifdef QLAT_USE_ACC
    const int nthreads = 32;
    size_t bN = (Nvol+nthreads-1)/nthreads;
    /////#pragma acc host_data use_device (a,b,f0)
    ////size_t bSize = nthreads*16*sizeof(Complexq);
    ////prodab_global<<< bN, nthreads, bSize >>>(a,b,f0, Nvol, Nsum);
    prodab_global<<< bN, nthreads >>>(a,b,f0, Nvol, Nsum);
    qacc_barrier(dummy);
    return ;
    #endif
  }

  if(mode_reduce == 0){
    /////ELSE do with CPU
    long m = 4;
    long n = 4;
    long w = 3*Nsum;
    #pragma omp parallel for
    for(long i=0;i<Nt*4*4;i++){f0[i] = 0;}
    //#pragma omp parallel for
    //for(int it=0;it<Nt;it++){
    //  matrix_prod_cpu(&b[it*m*w], &a[it*n*w], &f0[it*m*n] , m,n, w, 1, true, true);
    //}
    matrix_prod_cpu(b,a, f0 , m, n, w, Nt, true, false);
    
    return ;
  }
}

inline void reducefM(qlat::vector_acc<Complexq > &fd,Complexq* NabL, long bufN, std::vector<ga_M > &gL,const Geometry &geo,const int nvec,const Ftype facvol, unsigned long bufi, int mode_reduce=1)
{
  /////const qlat::Geometry &geo = a0.geo();
  ////const Coordinate vg = geo.total_site();
  ////int nt = vg[3];

  unsigned long Nvol = geo.local_volume();
  int Nt = geo.node_site[3];
  qlat::vector_acc<Complexq > reduce_sum;reduce_sum.resize((nvec*Nt)*16);

  if(mode_reduce == 1)
  {
    long Nsum = Nvol/Nt;
    if(Nsum%32 !=0){print0("Assumed Nsum divide 32 == 0, %8d \n", int(Nsum%32));qassert(false);}
    Nsum = Nsum/32;

    set_zero(reduce_sum);

    int bSum = 256/4;
    int cutN  = 32;
    reduce_gpu2d_6(fd.data(),&reduce_sum[0],Nsum,nvec*Nt*16,  1, bSum, cutN);
  }
  if(mode_reduce == 0)
  {
    #pragma omp parallel for
    for(int i=0;i<reduce_sum.size();i++){reduce_sum[i] = fd[i];}
  }
  

  ///TODO correct the reduce_gamma to GPU
  //qlat::vector_acc<Complexq > NabL_tem;NabL_tem.resize(nvec*Nt*16);
  //qthread_for(op0, nvec*Nt*16, {
  //  int ivec = op0/(Nt*16);int op = op0%(Nt*16);
  //  int it = op/16; int gi = op%16;
  //  NabL_tem[ivec*Nt*16 + it*16 + gi] = reduce_gamma(&reduce_sum[ivec*Nt*16 + it*16+0], gL[gi])/facvol;
  //});

  //qacc_for(op0, nvec*Nt*16, {
  //  int ivec = op0/(Nt*16);int op = op0%(Nt*16);
  //  int it = op/16; int gi = op%16;
  //  NabL[ivec*Nt*bufN*16 + it*bufN*16 + bufi*16 + gi] += NabL_tem[ivec*Nt*16 + it*16 + gi];
  //});

  qlat::vector_acc<Complexq* > gP; qlat::vector_acc<int* > iP;get_g_pointer(gL, gP, iP);
  /////qacc_for(i0, 1,{
  /////  NabL[0] = gP[0][0];
  /////});

  qacc_for(op0, nvec*Nt*16, {
    int ivec = op0/(Nt*16);int op = op0%(Nt*16);
    int it = op/16; int gi = op%16;
    NabL[ivec*Nt*bufN*16 + it*bufN*16 + bufi*16 + gi] += reduce_gamma(&reduce_sum[ivec*Nt*16 + it*16+0], gP[gi], iP[gi])/facvol;
  });
  /////============================


  //qacc_for(op0, nvec*Nt*16, {
  //  NabL[op0] = NabL_tem[op0];
  //});

  //if(mode_nt == 1)
  //{
  //qacc_for(op0, nvec*Nt*16, {
  //  int ivec = op0/(Nt*16);int op = op0%(Nt*16);
  //  int it = op/16; int gi = op%16;
  //  ////Nab[ivec*Nt*16 + it*16 + gi] += reduce_gamma(&reduce_sum[ivec*Nt*16 + it*16+0],gL[gi])/facvol;
  //  NabL[ivec*Nt*bufN*16 + it*bufN*16 + bufi*16 + gi] += reduce_gamma(&reduce_sum[ivec*Nt*16 + it*16+0],gL[gi])/facvol;
  //});
  /////unsigned long bufN = NabL.size()/(nt*16);
  ////#pragma omp parallel for
  ////for(int op0=0;op0<nvec*Nt*16;op0++){
  ////  int ivec = op0/(Nt*16);int op = op0%(Nt*16);
  ////  int it = op/16; int gi = op%16;
  ////  ////Nab[ivec*Nt*16 + it*16 + gi] += reduce_gamma(&reduce_sum[ivec*Nt*16 + it*16+0],gL[gi])/facvol;
  ////  NabL[ivec*Nt*bufN*16 + it*bufN*16 + bufi*16 + gi] += reduce_gamma(&reduce_sum[ivec*Nt*16 + it*16+0],gL[gi])/facvol;
  ////}
  //}

  //if(mode_nt == 0){
  //Coordinate xl = geo.coordinate_from_index(0);xl[3] = 0;
  //Coordinate xg = geo.coordinate_g_from_l(xl);
  //int tini = xg[3];
  //int  Nmpi   = qlat::get_num_node();
  ///////unsigned long bufN = NabL.size()/(Nmpi*nt*16);
  //qacc_for(op0, nvec*Nt*16, {
  //  int ivec = op0/(Nt*16);int op = op0%(Nt*16);
  //  int it = op/16; int gi = op%16;
  //  ////Nab[ivec*Nt*16 + it*16 + gi] += reduce_gamma(&reduce_sum[ivec*Nt*16 + it*16+0],gL[gi])/facvol;
  //  NabL[ivec*nt*bufN*16 + (it+tini)*bufN*16 + bufi*16 + gi] += reduce_gamma(&reduce_sum[ivec*Nt*16 + it*16+0],gL[gi])/facvol;
  //});

  //#pragma omp parallel for
  //for(int op0=0;op0<nvec*Nt*16;op0++){
  //  int ivec = op0/(Nt*16);int op = op0%(Nt*16);
  //  int it = op/16; int gi = op%16;
  //  ////Nab[ivec*Nt*16 + it*16 + gi] += reduce_gamma(&reduce_sum[ivec*Nt*16 + it*16+0],gL[gi])/facvol;
  //  NabL[ivec*nt*bufN*16 + (it+tini)*bufN*16 + bufi*16 + gi] += reduce_gamma(&reduce_sum[ivec*Nt*16 + it*16+0],gL[gi])/facvol;
  //}
  //}

}

inline void multiplyNab_Global(const Complexq* Nab, qlat::vector_acc<Ftype > &Mres,std::vector<int > avL, std::vector<int > bvL,const qlat::vector_acc<Complexq > &values, qlat::vector_acc<signed char> &GmapM,const int &nmass,const int &nt,const int nzero,const unsigned long bufN0, int mode_reduce=1)
{
  unsigned long bufN = avL.size();
  if(bufN == 0)return;

  /////int nmass = values.size()/n_vec;
  /////int nt = Nab.size()/(Aoper);

  /////Set up Mvalues
  qlat::vector_acc<Ftype > MvaluesV;MvaluesV.resize(bufN*nmass*2);
  Ftype* Mvalues = MvaluesV.data();
  Ftype* MresP   = Mres.data();
  #pragma omp parallel for
  for(unsigned long bmi=0;bmi<bufN*nmass;bmi++){
    unsigned long bi = bmi/nmass;
    int mi = bmi%nmass;
    int av = avL[bi];
    int bv = bvL[bi];

    int caseab = 2;
    Ftype fac_ab = 2.0;
    if(av==bv)fac_ab = 1.0;
    if(av >= nzero and bv >= nzero ){caseab = 2;}
    if(av >= nzero and bv <  nzero ){caseab = 1;}
    if(av <  nzero and bv <  nzero ){caseab = 0;}

    Complexq la = values[av*nmass + mi];
    Complexq lb = values[bv*nmass + mi];

    unsigned long offM = bi*nmass*2 + mi*2;
    if(caseab == 2){
      Mvalues[offM+0] = 2*fac_ab*(la*lb).real();
      Mvalues[offM+1] = 2*fac_ab*(la*qlat::qconj(lb)).real();
    }
    if(caseab == 1)
    {
      Mvalues[offM+0] = 2*fac_ab*(la*lb).real();
      Mvalues[offM+1] = 0.0;
    }

    if(caseab == 0)
    {
      Mvalues[offM+0] = fac_ab*(la*lb).real();
      Mvalues[offM+1] = 0.0;
    }
  }

  ////nt distributed
  //int rank = qlat::get_id_node();
  //long offNab = rank*nt*bufN0*16;
  //if(mode_nt == 1){offNab = 0;}
  long offNab = 0;

  ////GPU prod
  if(mode_reduce == 1){
    #ifdef QLAT_USE_ACC
    int rank = qlat::get_id_node();

    int  nthreads = nt;
    long nB = bufN;
    //long largeB = nt*16;if(nmass > largeB){largeB = nmass;}
    //int sizeB = largeB*sizeof(Complexq);
    long largeB = nt*16;largeB += nmass;
    int sizeB = largeB*sizeof(Complexq);
    /////Sum over time
    /////if(nt < 16){print0("time too short for production. \n");qassert(false);}
    multiplyNab_global<<< nB, nthreads, sizeB >>>(&Nab[offNab],MresP,&Mvalues[0],nt,nmass,bufN0, GmapM.data());
    /////multiplyNab_global<<< nB, nthreads, sizeB >>>(Nab, 0,&Mres[0],&Mvalues[0],nt,nmass,bufN0, &GmapM[0]);
    qacc_barrier(dummy);
    return ;
    #endif
  }

  ///CPU version
  if(mode_reduce == 0){
    signed char G0[32];
    signed char G1[32];
    signed char G2[32];
    signed char G3[32];

    /////unsigned int tid = threadIdx.x;
    /////long bi = blockIdx.x;
    for(unsigned int off=0;off<32;off++){G0[off] = GmapM[0*32+off];}
    for(unsigned int off=0;off<32;off++){G1[off] = GmapM[1*32+off];}
    for(unsigned int off=0;off<32;off++){G2[off] = GmapM[2*32+off];}
    for(unsigned int off=0;off<32;off++){G3[off] = GmapM[3*32+off];}

    /////NabL[ivec*nt*bufN*16 + it*bufN*16 + bufi*16 + gi]
    /////Each rank do one multiply
    ////long offNab = 0;
    #pragma omp parallel for
    for(LInt bi=0; bi<bufN; bi++)
    {
      std::vector<Complexq > Mv; Mv.resize(nt*16);
      /////Copy original data
      for(int ti=0;ti<nt;ti++)
      for(unsigned int op=0;op<16;op++)
      {
        Mv[ti*16 + op] = Nab[offNab + ti*(bufN0*16) + bi*16 + op];
      }

      /////Copy masses
      unsigned long offAB = bi*nmass*2;
      std::vector<Ftype > buf;buf.resize(nmass*2);
      Ftype *src = (Ftype*) &Mvalues[offAB];
      for(int off=0;off<nmass*2;off++){buf[off] = src[off];}

      /////Buffer for A,B values
      const Complexq *Av,*Bv;
      /////Buffer for prods
      Complexq v0[16];Complexq v1[16];

      /////unsigned long shiftM = bufN0*nt*nt;
      ////Do calculations
      for(int t0=0;t0<nt;t0++)
      for(int t1=0;t1<nt;t1++)
      {
        unsigned int to = (t0+t1)%nt;
        Bv = &Mv[t0*16 + 0];
        Av = &Mv[to*16 + 0];
        for(int ipr=0;ipr<16;ipr++){
          v0[ipr] = Av[G0[ipr*2+0]]*qlat::qconj(Bv[G1[ipr*2+0]])
           *Ftype(G0[ipr*2+1])*Ftype(         G1[ipr*2+1]);
          v1[ipr] = Av[G2[ipr*2+0]]*qlat::qconj(Bv[G3[ipr*2+0]])
           *Ftype(G2[ipr*2+1])*Ftype(         G3[ipr*2+1]);
        }

        //offAB = (((0*nmass+0)*bufN0 + bi)*nt+t0)*nt + t1;
        Ftype *res = (Ftype*) &MresP[((bi*nt+t0)*nt+t1)*16*nmass];

        for(int ipr=0;ipr<16;ipr++)
        {
          Ftype v00 = v0[ipr].real();
          Ftype v10 = v1[ipr].real();
          for(int mi=0;mi<nmass;mi++)
          {
            res[ipr*nmass + mi] += (buf[mi*2+0]*v00 + buf[mi*2+1]*v10);
            //Mres[offAB] += (buf[mi*2+0]*v00 + buf[mi*2+1]*v10);
            //offAB += shiftM;
          }
        }
      }
    }
    return ;
  }

 
}


inline std::vector<unsigned long > get_loop_cut(int Nx,int Ny, int Nycut, int Nxcut){
  std::vector<unsigned long >jobL;
  jobL.resize(Ny*Nx);
  int Nx_bound = (Nx+Nxcut-1)/Nxcut;
  int Ny_bound = (Ny+Nycut-1)/Nycut;
  long count = 0;
  for(int lyi=0;lyi<Ny_bound;lyi++)
  for(int lxi=0;lxi<Nx_bound;lxi++)
  for(int ayi=0;ayi<Nycut;ayi++)
  for(int axi=0;axi<Nxcut;axi++)
  {
    long yi = lyi*Nycut + ayi;
    long xi = lxi*Nxcut + axi;
    if(xi < Nx and yi < Ny){
      jobL[count] = yi*Nx + xi;count += 1;
    }
  }
  return jobL;
}

inline void get_map_gammaL(std::vector<ga_M > &g0,std::vector<ga_M > &gL,qlat::vector_acc<signed char > &Gmap){
  Gmap.resize(32);
  for(int i=0;i<16;i++){
    unsigned long r0;unsigned long r1;
    unsigned long a0;unsigned long a1;
    signed char findi =-1;
    signed char sign = 1;
    g0[i].check_sum(r0,r1);
    for(int j=0;j<16;j++){
      gL[j].check_sum(a0,a1);
      if(r0==a0){
        if(findi != -1){print0("WRONG!!!!\n");}
        findi = j;
      }
      if(r0==a1){
        if(findi != -1){print0("WRONG!!!!\n");}
        findi = j;
        sign  = -1;
      }
    }
    if(findi == -1){print0("WRONG!!!! %5d \n",findi);}
    Gmap[i*2+0] = findi;
    Gmap[i*2+1] = sign;
  }
  
}

struct Nab_distribute{

  int mxyz;
  int NabL_size;
  int bufN;
  int rank;
  ////Complexq* NabN;
  //qlat::vector<Complexq > NabN;
  qlat::vector_gpu<Complexq > NabN;
  MPI_Comm xyz_comm;
  MPI_Comm t_comm;
  std::vector<int > rank_map;
  std::vector<int > send,recv,spls,rpls;

  Nab_distribute(const fft_desc_basic &fd, int bufN_or){
    TIMER("Create Nab_distribute");

    NabL_size = 16*fd.Nt*fd.Nmpi;
    mxyz = fd.mx*fd.my*fd.mz;
    rank = fd.rank;

    int color_xyz = fd.init;
    MPI_Comm_split(get_comm(), color_xyz, fd.rank, &xyz_comm);

      color_xyz = (fd.iniz * fd.ny + fd.iniy)*fd.nx + fd.inix;
    int rank_t    = fd.init;
    MPI_Comm_split(get_comm(), color_xyz, rank_t, &t_comm);
    rank_map.resize(fd.mt);

    {
      int rank;
      MPI_Comm_rank(t_comm, &rank);
      rank_map[rank] = fd.rank;
    }
    sum_all_size(&rank_map[0], rank_map.size() , 0, &t_comm);

    ////NabN = NULL;
    set_bufN(fd, bufN_or);
  
  }

  void set_bufN(const fft_desc_basic &fd, int bufN_or)
  {
    bufN = bufN_or;
    long size_c = fd.Nt*bufN*16 * sizeof(Complexq);
    send.resize(fd.mt);
    recv.resize(fd.mt);
    spls.resize(fd.mt);
    rpls.resize(fd.mt);

    for(int ri=0;ri<fd.mt;ri++)
    {
      send[ri] = size_c;
      spls[ri] = size_c*rank_map[ri];

      recv[ri] = size_c;
      rpls[ri] = size_c*ri;
    }

    NabN.resize(bufN*NabL_size);
    /////gpuFree(NabN);gpuMalloc(NabN, bufN*NabL_size, Complexq);
    //NabN.resize(bufN*NabL_size);
  }

  void communicate(Complexq* NabL)
  {
    TIMER("Reduce Nab");
    if(mxyz!=1){sum_all_size(&NabL[0], NabN.data(), bufN*NabL_size , 1, &xyz_comm);}
    else{
      cpy_data_thread(NabN.data(), &NabL[0], bufN*NabL_size, 1);
      ///#ifdef QLAT_USE_ACC
      ///cudaMemcpy(&NabN[0], &NabL[0], bufN*NabL_size*sizeof(Complexq), cudaMemcpyDeviceToDevice);
      ///#else
      ///memcpy(&NabN[0], &NabL[0], bufN*NabL_size*sizeof(Complexq));
      ///#endif
    }

    //MPI_Alltoallv(NabN.data(),(int*) &send[0],(int*) &spls[0], MPI_CHAR,
    //              &NabL[0]   ,(int*) &recv[0],(int*) &rpls[0], MPI_CHAR, t_comm);

    int GPU = 1;
    MPI_Alltoallv_mode(NabN.data(),(int*) &send[0],(int*) &spls[0],
                       &NabL[0]   ,(int*) &recv[0],(int*) &rpls[0], t_comm, 1, GPU);
  }

  //~Nab_distribute(){
  //  NabN.resize(0);
  //  /////gpuFree(NabN);NabN=NULL;
  //}

};

//void sum_NabL(Complexq* NabL, Nab_distribute& dis)
//{
//  TIMER("Reduce Nab");
//  dis()
//  ///qlat::vector<int > nv,Nv,mv;
//  ///geo_to_nv(geo, nv,Nv,mv);
//
//  /////rank  --> Nt/nt --> bufi --> 16
//  //if(mode_nt == 1){Redistribute_all_Nt(reinterpret_cast<Ftype* > (&NabL[0]), 2*sumsize,geo, 1);}
//
//  //if(mode_nt == 0)
//  //{
//  //  sum_all_size(reinterpret_cast<Ftype* > (&NabL[0]), 2*sumsize , 1);
//  //  //long bufN = NabL.size()/(fd.nt*fd.Nmpi  *16);
//  //  //MPI_Comm vec_comm;
//  //  //int tini = fd->init;
//  //  //int color_xyz = tini;
//  //  //MPI_Comm_split(get_comm() ,color_xyz, fd.rank,&vec_comm);
//
//  //  //MPI_Datatype curr = MPI_DOUBLE;unsigned int M_size = sizeof(double);
//  //  //M_size = get_MPI_type<Ftype >(curr );qassert(M_size <= FLOATIND+3);
//
//  //  //Ftype* src = reinterpret_cast<Ftype* > (&NabL[0]);
//  //  //long size_sum = fd.Nt*fd.Nmpi*16 * 2;
//  //  //long src_off = (tini/fd.Nt)*size_sum
//  //  //MPI_Allreduce(&src[src_off],&src[src_off], size_sum, curr, MPI_SUM, vec_comm);
//
//  //  //MPI_Gather( void* send_data,
//  //  //    int send_count,
//  //  //    MPI_Datatype send_datatype,
//  //  //    void* recv_data,
//  //  //    int recv_count,
//  //  //    MPI_Datatype recv_datatype,
//  //  //    int root,
//  //  //    get_comm());
//  //}
//
//}

inline void get_low_rho(std::vector<qlat::FieldM<Complexq, 12>  > &eigen,const qlat::vector_acc<Complexq > &values,const int &nzero,qlat::vector_acc<Ftype > &Mres,const qlat::Geometry &geo, int GPUFM=1)
{
  ////Input must be chiral vectors, eigen_chi, n_vec --> chi --> d/2 --> t,y,z,x --> c --> complex
  ////values --> n_vec --> massi
  ////
  const int n_vec = eigen.size();
  const int nmass = values.size()/n_vec;
  const Coordinate vg = geo.total_site();
  const int nt = vg[3];

  fft_desc_basic fd(geo);

  Ftype facvol = std::sqrt(vg[0]*vg[1]*vg[2]);

  ////0 for nt not on MPIs, 1 for nt on MPIs 
  //int mode_nt = 1;
  ////0 use CPU to reduce, 1 use GPU to reduce
  int mode_reduce = 1;
  //if(GPUFM == 0){mode_nt = 0;mode_reduce = 0;}
  //if(GPUFM == 1){mode_nt = 1;mode_reduce = 1;}
  //if(GPUFM == 2){mode_nt = 0;mode_reduce = 1;}
  if(GPUFM == 0){mode_reduce = 0;}
  if(GPUFM == 1){mode_reduce = 1;}
  if(GPUFM == 2){mode_reduce = 1;}


  /////Get map list
  unsigned short Nt = geo.node_site[3];

  std::vector<ga_M > gL;gL.resize(Aoper);
  qlat::vector_acc<signed char> GmapM;GmapM.resize(32*4);

  {
  TIMER("Copy gammas");
  ga_matrices_cps   ga_cps;
  std::vector<ga_M > g0;g0.resize(Aoper);
  std::vector<ga_M > g05;g05.resize(Aoper);
  std::vector<ga_M > g1;g1.resize(Aoper);
  std::vector<ga_M > g15;g15.resize(Aoper);
  //////0 , 1, 2, 3, 4, 5
  //////1-2, 1-3, 1-4, 1-5
  //////2-3, 2-4, 2-5
  //////3-4, 3-5
  //////4-5
  {int o=0;
  for(int i=0;i<6;i++){gL[o] = ga_cps.ga[0][i];o+=1;}
  for(int i=2;i<6;i++){gL[o] = ga_cps.ga[1][i];o+=1;}
  for(int i=3;i<6;i++){gL[o] = ga_cps.ga[2][i];o+=1;}
  for(int i=4;i<6;i++){gL[o] = ga_cps.ga[3][i];o+=1;}
  for(int i=5;i<6;i++){gL[o] = ga_cps.ga[4][i];o+=1;}}

  ga_M g5;
  g5 = ga_cps.ga[0][5];
  ////match gammas with twopt functions
  for(int o=0;o<16;o++){
    ////gL[o] = ga_cps.ga[0][5] * gL[o];
    gL[o] = gL[o] * g5;
    /////gL[o] = g5 * gL[o];
  }
  ////GL

  for(int i=0;i<Aoper;i++){
     g0[i] = gL[i];
    g05[i] = (g5*gL[i])*g5;
     g1[i] = gL[i]*g5;
    g15[i] = g5*gL[i];
  }

  ///std::vector<std::vector<int > > Gmap;Gmap.resize(4);
  ///for(int gi=0;gi<4;gi++){Gmap[gi].resize(32);}
  //qlat::vector<signed char> Gmap;///g0
  qlat::vector_acc<signed char> Gmap0;///g0
  qlat::vector_acc<signed char> Gmap1;///g05
  qlat::vector_acc<signed char> Gmap2;///g1
  qlat::vector_acc<signed char> Gmap3;///g15
  get_map_gammaL(g0 ,gL, Gmap0);for(int i=0;i<32;i++){GmapM[0*32+i] = Gmap0[i];}
  get_map_gammaL(g05,gL, Gmap1);for(int i=0;i<32;i++){GmapM[1*32+i] = Gmap1[i];}
  get_map_gammaL(g1 ,gL, Gmap2);for(int i=0;i<32;i++){GmapM[2*32+i] = Gmap2[i];}
  get_map_gammaL(g15,gL, Gmap3);for(int i=0;i<32;i++){GmapM[3*32+i] = Gmap3[i];}


  //#ifdef QLAT_USE_ACC
  //cudaMemcpyToSymbol(Gmap0C, &Gmap0[0],32*sizeof(signed char),0 , cudaMemcpyHostToDevice);
  //cudaMemcpyToSymbol(Gmap1C, &Gmap1[0],32*sizeof(signed char),0 , cudaMemcpyHostToDevice);
  //cudaMemcpyToSymbol(Gmap2C, &Gmap2[0],32*sizeof(signed char),0 , cudaMemcpyHostToDevice);
  //cudaMemcpyToSymbol(Gmap3C, &Gmap3[0],32*sizeof(signed char),0 , cudaMemcpyHostToDevice);
  //#endif

  }

  /////for(int i=0;i<32;i++){GmapM[0*32+i] = Gmap0[i];}
  /////for(int i=0;i<32;i++){GmapM[1*32+i] = Gmap1[i];}
  /////for(int i=0;i<32;i++){GmapM[2*32+i] = Gmap2[i];}
  /////for(int i=0;i<32;i++){GmapM[3*32+i] = Gmap3[i];}

  /////int Ncut = n_vec;
  int Ncut = (n_vec-1)/2 + 1;

  int  Nmpi   = qlat::get_num_node();
  long npoints = eigen[0].geo().local_volume()*12;
  ////double vGb_vec = npoints*Nmpi*2.0/(1024.0*1024*1024);
  double vGb_vec = npoints*2.0/(1024.0*1024*1024);

  int meas = 4;int Fcount = 3 + 1;////((complex multi 6 + plus 2)/2)
  double vGb     = vGb_vec*meas*Fcount;
  print0("==total Eigen %.3e Gb \n", vGb_vec*Nmpi*(sizeof(Complexq)/2.0)*n_vec);

  ////double length = (geo.local_volume()*pow(0.5,30))*12*sizeof(Complexq);
  int modeCopy = 0;

  ///int N_bound = (n_vec+Ncut-1)/Ncut;
  ///Ncut = n_vec;
  Complexq* a0p;Complexq* b0p;

  qlat::vector_acc<Complexq > prodFMV;
  if(mode_reduce == 0)prodFMV.resize(Nmpi*Nt*16);
  if(mode_reduce == 1)prodFMV.resize(Nmpi*geo.local_volume()*16/32);
  Complexq* prodFM = prodFMV.data();

  long NabL_size  = 16*Nt*Nmpi;
  long MresL_size = nmass*16*nt*nt;

  //NabL_size = 16*nt;
  //if(mode_nt==1){NabL_size = 16*nt;}
  //if(mode_nt==0){NabL_size = 16*nt*Nmpi;}

  print0("===job start 0 \n");
  fflush_MPI();

  int Ncutbuf = Ncut;
  /////Buffer for Nab products
  unsigned long bufN = 1;
  /////Get Ncut
  #ifdef QLAT_USE_ACC
  int bufa0 = -1;int bufa1 = -1;
  int bufb0 = -1;int bufb1 = -1;
  int facbufN = 1;

  size_t freeM = 0;size_t totalM = 0;double extra = 0.2;double fac_extra=1.5;
  #ifdef __HIP_PLATFORM_HCC__
  extra = 0.1;fac_extra = 1.1; 
  #endif

  modeCopy = 1;

  long total_vol = Nmpi*geo.local_volume();
  bufN = 2*((total_vol*12)/(nmass*16*nt*facbufN));
  double membufN = (MresL_size + NabL_size)*sizeof(Complexq)*pow(0.5,30);
  /////CPU memory limit, only bufN 1 or a small number

  ////if(modeCopy == 1){extra = 0.0;facbufN = 5;}
  cudaMemGetInfo(&freeM,&totalM);
  double freeD = freeM*pow(0.5,30);double totalD = totalM*pow(0.5,30);
  if(membufN * bufN > extra*totalD){bufN = int(extra*totalD/(membufN));if(bufN == 0)bufN = 1;}
  freeD = freeD - membufN * bufN;

  int Nfull = (freeD/(vGb_vec*sizeof(Complexq)/2.0));
  if(n_vec < Nfull/(fac_extra)){
    //Ncut = (n_vec-1)/2 + 1;
    Ncut = n_vec;
    Ncutbuf = Ncut;
  }
  else{Ncut = Nfull/(2.5);Ncutbuf = 2*Ncut;}///(or 2.5)

  ////==propagate setups
  {if(fd.rank != 0){Ncutbuf = 0;Ncut = 0;}
  sum_all_size(&Ncutbuf, 1);sum_all_size(&Ncut   , 1);}
  print0("==rank %d, n_vec %8d, Ncut %5d/%5d , Fac %.3e , free %.3e GB, total %.3e GB \n",
      qlat::get_id_node(), n_vec,Ncut,Nfull,n_vec*1.0/Ncut,freeD, totalD);
  ////==propagate setups

  #ifdef __HIP_PLATFORM_HCC__
  if(n_vec != Ncut ){abort_r("HIPCC has propblem with memory move!\n");}
  #endif

  #endif

  unsigned long bufi = 0;
  print0("===rank %d, bufN %lu, mem MresL %.3e, NabL %.3e \n", qlat::get_id_node(), 
      bufN, bufN*MresL_size*sizeof(Complexq)*pow(0.5,30), bufN*NabL_size*sizeof(Complexq)*pow(0.5,30));
  print0("===rank %d, bufE %d, mem %.3e \n", qlat::get_id_node(), Ncutbuf, Ncutbuf*npoints*sizeof(Complexq)*pow(0.5,30));

  print0("===job start 1 \n");
  fflush_MPI();

  std::vector<Complexq* > bufE;
  if(modeCopy == 1){
    TIMER("CUDA mem allocate");
    bufE.resize(Ncutbuf);
    for(unsigned long iv=0;iv<bufE.size();iv++){gpuMalloc(bufE[iv], npoints, Complexq);}
    //for(unsigned long iv=0;iv<bufE.size();iv++){bufE[iv].init(eigen[0].geo());}
    //for(unsigned long iv=0;iv<bufE.size();iv++){cudaMalloc(&bufE[iv],  npoints*sizeof(Complexq));}
  }

  print0("===job start 2 \n");
  fflush_MPI();

  //Mres.resize(nmass*16*nt*nt);set_zero(Mres);
  qlat::vector_acc<Ftype > MresL;
  {TIMER("CUDA mem allocate");MresL.resize(bufN*MresL_size);set_zero(MresL);}
  //qlat::vector<Complexq > Nab;Nab.resize(16*nt);
  //set_zero(Nab);

  Nab_distribute dis(fd, bufN);

  qlat::vector_gpu<Complexq > NabV;NabV.resize( bufN*NabL_size);
  NabV.set_zero();

  Complexq* NabL = NabV.data();
  ///{TIMER("CUDA mem allocate");gpuMalloc(NabL, bufN*NabL_size, Complexq);}
  //gpuErrchk(cudaMalloc(&NabL, bufN*NabL_size*sizeof(Complexq)))
  //qlat::vector<Complexq > NabL;NabL.resize(bufN*NabL_size);
  //NabL.resize(bufN*16*nt);
  //set_zero(NabL);
  //#ifdef QLAT_USE_ACC
  //cudaMemset(NabL, 0, bufN*NabL_size*sizeof(Complexq));
  //#else
  //qacc_for(i, bufN*NabL_size, {NabL[i] = 0.0;});
  //#endif

  print0("===job start 3 \n");
  fflush_MPI();


  std::vector<unsigned long > jobL = get_loop_cut(n_vec,n_vec,Ncut,Ncut);
  int countrun = 0;int totrun =  0;timeval tm0,tm1,tm2;
  gettimeofday(&tm0, NULL);gettimeofday(&tm1, NULL);gettimeofday(&tm2, NULL);
  int eachrun  = 0;
  for(LInt jobi=0;jobi<jobL.size();jobi++){
    int avi = jobL[jobi];
    int av = avi/n_vec;
    int bv = avi%n_vec;
    if(bv > av){continue;}totrun +=1;
  }

  print0("===job start 4 \n");
  fflush_MPI();


  ////qlat::vector<Complexq > NabS;NabS.resize(16*nt);
  /////Rederive chiral forms
  std::vector<int > avL,bvL;
  avL.resize(0);bvL.resize(0);
  std::vector<int > avL_local,bvL_local;
  avL_local.resize(0);bvL_local.resize(0);

  print0("===job start n \n");
  fflush_MPI();

  ////Buffer index for av,bv
  //////#pragma omp parallel for
  for(LInt jobi=0;jobi<jobL.size();jobi++)
  {
    TIMER("Kernel jobs");
    int avi = jobL[jobi];
    int av = avi/n_vec;
    int bv = avi%n_vec;
    if(bv > av)continue;

    avL.push_back(av);
    bvL.push_back(bv);

    if(modeCopy == 0){
      a0p = (Complexq* ) qlat::get_data(eigen[av]).data();
      b0p = (Complexq* ) qlat::get_data(eigen[bv]).data();
    }

    //Buffer for bv
    #ifdef QLAT_USE_ACC
    if(modeCopy == 1)
    {
    TIMER("Copy memory to Device");

    if(bv >= bufb0 and bv < bufb1){b0p = bufE[bv%Ncut];}else{
      if(bv % Ncut == 0){
        for(int iv=0;iv<Ncut;iv++){
          if(bv + iv < n_vec)cudaMemcpyAsync(bufE[iv], qlat::get_data(eigen[bv+iv]).data(),
            npoints*sizeof(Complexq), cudaMemcpyHostToDevice);
        }
        qacc_barrier(dummy);
        bufb0 = (bv/Ncut)*Ncut;bufb1 = bufb0 + Ncut;
      }
      b0p = bufE[bv%Ncut];
    }

    if(av >= bufb0 and av < bufb1){a0p = bufE[av%Ncut];}else{
    if(av >= bufa0 and av < bufa1){a0p = bufE[Ncut+av%Ncut];}else{
      if(av % Ncut == 0){
        for(int iv=0;iv<Ncut;iv++){
          if(av + iv < n_vec)cudaMemcpyAsync(bufE[Ncut + iv], qlat::get_data(eigen[av+iv]).data(),
            npoints*sizeof(Complexq), cudaMemcpyHostToDevice);
        }
        qacc_barrier(dummy);
        bufa0 = (av/Ncut)*Ncut;bufa1 = bufa0 + Ncut;
      }
      a0p = bufE[Ncut + av%Ncut];
    }
    }
    }
    #endif

    ////Vector reduce to Nab
    long off_FM = 0;
    if(mode_reduce == 0)off_FM=((countrun%Nmpi)*Nt*16);
    if(mode_reduce == 1)off_FM=((countrun%Nmpi)*geo.local_volume()*16)/32;
    {
      TIMER("Prod core a b");
      prodab(a0p, b0p, geo,&prodFM[off_FM], mode_reduce);
      //{TIMER("Touch 1");touchmem(*a0p,fd);}
      //{TIMER("Touch 2");touchmem(*b0p,fd);}
    }
    countrun += 1;

    if((countrun%Nmpi == 0) or countrun == totrun)
    {
      int nvec = Nmpi;if(countrun%Nmpi != 0)nvec = countrun%Nmpi;
      {TIMER("Reduce prodFM");reducefM(prodFMV, NabL, bufN,gL,geo,nvec,facvol, bufi, mode_reduce);}
      ////bufi += nvec;
      bufi += 1;

      int rank = qlat::get_id_node();
      if(rank < nvec){
        avL_local.push_back(avL[rank]);bvL_local.push_back(bvL[rank]);
      }

      if(avL_local.size() == bufN or countrun == totrun )
      {
        //////NabL[ivec*Nt*bufN*16 + it*bufN*16 + bufi*16 + gi]
        //////nvec --> it --> bufi --> 16
        //sum_NabL(NabL, bufN*NabL_size, bufN, geo, fd, mode_nt);

        dis.communicate(NabL);

        {TIMER("Sum at zero node.");
        multiplyNab_Global(NabL,MresL,avL_local,bvL_local,values, GmapM, nmass, nt, nzero, bufN, mode_reduce);}

        avL_local.resize(0);bvL_local.resize(0);

        NabV.set_zero();
        /////qacc_for(i, bufN*NabL_size, {NabL[i] = 0.0;});
        /////set_zero(NabL);
        //{
        //TIMER("Set zero Nab")
        //#ifdef QLAT_USE_ACC
        //cudaMemset(NabL, 0, bufN*NabL_size*sizeof(Complexq));
        //#else
        //qacc_for(i, bufN*NabL_size, {NabL[i] = 0.0;});
        //#endif
        //}
        bufi = 0;
      }
      avL.resize(0);bvL.resize(0);
      ////set_zero(Nab);
    }

    double perc = countrun*1.0/totrun;
    eachrun  += 1;

    if(jobi%(n_vec)==0){
      gettimeofday(&tm1, NULL);
      double time0 = tm1.tv_sec - tm0.tv_sec;
      time0 += (tm1.tv_usec - tm0.tv_usec)/1000000.0;

      double time1 = tm1.tv_sec - tm2.tv_sec;
      time1 += (tm1.tv_usec - tm2.tv_usec)/1000000.0;

      double flops_pers = vGb*countrun/(1.0*time0);
      double flops_pers_round = vGb*eachrun/(1.0*time1);eachrun=0;gettimeofday(&tm2, NULL);
      print0("==jobi %10d, ai %5d , bi %5d , per %.3f, use %.3e sec, %.3f Gflops, %.3f Gflops/r . \n",
        int(jobi),av,bv, perc, time0,flops_pers,flops_pers_round);
    }
  }

  {
    gettimeofday(&tm1, NULL);
    double time0 = tm1.tv_sec - tm0.tv_sec;
    time0 += (tm1.tv_usec - tm0.tv_usec)/1000000.0;

    double flops_pers = vGb*totrun/(1.0*time0);
    print0("==Total use %.3e sec, average %.3f Gflops . \n", time0,flops_pers);
  }

  #ifdef QLAT_USE_ACC
  qacc_barrier(dummy);
  #endif

  ///ipr --> mi --> bi -- > t0
  Mres.resize(nmass*16*nt*nt);set_zero(Mres);
  /////long Msum = Mres.size();
  {
  TIMER("Copy final result");
  long Nsize = nmass*16*nt;
  qacc_for(isp, Nsize, {
    int mi  =  isp/(16*nt); 
    int ipr = (isp%(16*nt))/(nt);
    int t0  = (isp)%(nt);
    Ftype *res =  &Mres[((mi*16 + ipr)*nt+t0)*nt + 0];

    if(mode_reduce == 1)
    for(unsigned long bi=0;bi<bufN;bi++){
      Ftype *src =  &MresL[(((ipr*nmass+mi)*bufN+bi)*nt+t0)*nt + 0];
      ///#pragma omp parallel for
      for(int ti=0;ti<nt;ti++){res[ti] += src[ti];}
    }

    if(mode_reduce == 0)
    for(unsigned long bi=0;bi<bufN;bi++){
      Ftype *src =  &MresL[((bi*nt+t0)*nt +0)*16*nmass + ipr*nmass + mi];
      ///#pragma omp parallel for
      for(int ti=0;ti<nt;ti++){res[ti] += src[ti*16*nmass];}
    }
  
  });

  //for(int mi=0;mi<nmass;mi++)
  //for(int ipr=0;ipr<16;ipr++)
  //for(unsigned int t0=0;t0<nt;t0++)
  //{
  //  Ftype *res =  &Mres[((mi*16 + ipr)*nt+t0)*nt + 0];
  //  if(mode_reduce == 1)
  //  for(unsigned long bi=0;bi<bufN;bi++){
  //    Ftype *src =  &MresL[(((ipr*nmass+mi)*bufN+bi)*nt+t0)*nt + 0];
  //    #pragma omp parallel for
  //    for(int ti=0;ti<nt;ti++){res[ti] += src[ti];}
  //  }

  //  if(mode_reduce == 0)
  //  for(unsigned long bi=0;bi<bufN;bi++){
  //    Ftype *src =  &MresL[((bi*nt+t0)*nt +0)*16*nmass + ipr*nmass + mi];
  //    #pragma omp parallel for
  //    for(int ti=0;ti<nt;ti++){res[ti] += src[ti*16*nmass];}
  //  }

  //}
  }

  {TIMER("Final sum Mres");sum_all_size(Mres.data(),Mres.size());}

  {TIMER("Free memory");
  if(modeCopy == 1){
    for(unsigned long iv=0;iv<bufE.size();iv++){gpuFree(bufE[iv]);bufE[iv] = NULL;}
    bufE.resize(0);
  }
  ////gpuFree(NabL);NabL = NULL;
  }


}



}

#endif

