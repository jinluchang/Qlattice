// utils_corr_meson.h
// Gen Wang
// Jul. 2021

#ifndef UTILS_CORR_MESON_H
#define UTILS_CORR_MESON_H

#pragma once

#include "utils_float_type.h"
#include "utils_gammas.h"
#include "utils_fft_desc.h"
#include "utils_reduce_vec.h"
#include "utils_grid_src.h"
#include "utils_io_vec.h"
#include "utils_corr_prop.h"

namespace qlat{

////ga1 sink gammas, ga2 src gammas
template <typename Td>
void meson_vectorE(std::vector<Propagator4dT<Td > > &pV1, std::vector<Propagator4dT<Td > > &pV2, ga_M &ga1,ga_M &ga2,
        qlat::vector_acc<qlat::ComplexT<Td > > &res, qlat::fft_desc_basic &fd,int clear=1){
  TIMER("Meson_vectorE");
  qassert(fd.order_ch == 0);
  ///////check_prop_size(prop1);check_prop_size(prop2);
  int  NTt  = fd.Nv[3];
  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  int  nmass = pV1.size();
  if(nmass == 0){res.resize(0);return;}

  if(clear == 1){ini_resE(res,nmass,fd);}

  if(res.size()%NTt !=0 or res.size()==0){print0("Size of res wrong. \n");qassert(false);}
  qassert(pV1.size() == pV2.size());

  for(int mi=0;mi<nmass;mi++)
  {
  Propagator4dT<Td >& pL1 = pV1[mi];
  Propagator4dT<Td >& pL2 = pV2[mi];

  qacc_for(isp, long(pV1[0].geo().local_volume()),{ 
    int ti = isp/Nxyz;
    int xi = isp%Nxyz;
      qlat::ComplexT<Td > pres;pres = 0.0;
      const qlat::WilsonMatrix& p1 =  pL1.get_elem_offset(isp);
      const qlat::WilsonMatrix& p2 =  pL2.get_elem_offset(isp);

      for(int d1=0;d1<4;d1++)
      for(int c1=0;c1<3;c1++)
      for(int d2=0;d2<4;d2++)
      {
      const qlat::ComplexT<Td > g_tem = ga2.g[d2]*ga1.g[d1];
      for(int c2=0;c2<3;c2++)
      {
        pres += g_tem * 
          p1(ga1.ind[d1]*3+c1,d2*3+c2) * qlat::qconj(p2(d1*3+c1,ga2.ind[d2]*3+c2)) ;
      }
      }
      res[(mi*NTt + ti)*Nxyz + xi%Nxyz] += pres;
  });
  }

}

template <typename Ty >
void meson_vectorE(std::vector<qpropT >& prop1, std::vector<qpropT >& prop2, ga_M &ga1,ga_M &ga2,
        qlat::vector_acc<Ty > &res, int clear=1, int invmode=1, const Ty factor = Ty(1.0, 0.0)){
  TIMER("Meson_vectorE");
  const qlat::Geometry &geo = prop1[0].geo();
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);

  int  NTt  = fd.Nv[3];
  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  int  nmass = prop1.size();  ////(12*12*NTt)
  if(nmass == 0){res.resize(0);return;}
  if(clear == 1){ini_resE(res, nmass, fd);}
  if(res.size()%NTt != 0 or res.size() == 0){print0("Size of res wrong. \n");qassert(false);}

  qassert(prop1.size() == prop2.size());
  qlat::vector_acc<Ty* > p1 = EigenM_to_pointers(prop1);
  qlat::vector_acc<Ty* > p2 = EigenM_to_pointers(prop2);

  for(int d2=0;d2<4;d2++)
  for(int c2=0;c2<3;c2++)
  for(int d1=0;d1<4;d1++)
  for(int c1=0;c1<3;c1++)
  {
  //#pragma omp parallel for
  for(int ji=0;ji<nmass*NTt;ji++)
  {
    int massi = ji/NTt;
    int ti    = ji%NTt;

    int off1 = (d2*3+c2)*12+ga1.ind[d1]*3+c1;
    int off2 = (ga2.ind[d2]*3+c2)*12+d1*3+c1;

    const Ty g_tem = ga2.g[d2]*ga1.g[d1];

    Ty* tp1 = &p1[massi][(off1*NTt+ti) * Nxyz];
    Ty* tp2 = &p2[massi][(off2*NTt+ti) * Nxyz];

    Ty* tr0 = &((res.data())[(massi*NTt + ti)*Nxyz]);

    #if USEQACC==1
    if(invmode == 1){qacc_forNB(i, long(Nxyz),{ tr0[i] += factor * (tp1[i]*qlat::qconj(tp2[i]) * g_tem);});}
    if(invmode == 0){qacc_forNB(i, long(Nxyz),{ tr0[i] += factor * (tp1[i]*           (tp2[i]) * g_tem);});}
    #else
    EAy vp1(tp1,Nxyz);
    EAy vp2(tp2,Nxyz);
    EAy vr0(tr0,Nxyz);
    if(invmode == 1)vr0 += factor * (vp1*vp2.conjugate() * g_tem);
    if(invmode == 0)vr0 += factor * (vp1*vp2             * g_tem);
    #endif
  }
  qacc_barrier(dummy);
  }
}

/////merge to each gamma with sortted cases
//template <typename Ty >
//void meson_vectorE(std::vector<qpropT >& prop1, std::vector<qpropT >& prop2,
//        qlat::vector_acc<Ty > &res, int clear=1, int invmode=1){
//  TIMER("Meson_vectorE");
//  const qlat::Geometry &geo = prop1[0].geo();
//  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);
//
//  int  NTt  = fd.Nv[3];
//  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
//  int  nmass = prop1.size();  ////(12*12*NTt)
//  if(nmass == 0){res.resize(0);return;}
//  if(clear == 1){ini_resE(res, nmass*16, fd);}
//  if(res.size()%NTt != 0 or res.size() == 0){print0("Size of res wrong. \n");qassert(false);}
//
//  qassert(prop1.size() == prop2.size());
//  qlat::vector_acc<Ty* > p1 = EigenM_to_pointers(prop1);
//  qlat::vector_acc<Ty* > p2 = EigenM_to_pointers(prop2);
//
//  for(int ds=0;ds<4;ds++)
//  for(int d2=0;d2<4;d2++)
//  for(int d1=0;d1<4;d1++)
//  for(int c2=0;c2<3;c2++)
//  for(int c1=0;c1<3;c1++)
//  {
//  //#pragma omp parallel for
//  for(int ji=0;ji<nmass*NTt;ji++)
//  {
//    int massi = ji/NTt;
//    int ti    = ji%NTt;
//    const int offdi = d2*4+d1;
//
//    int off1 = (ds*3+c2)*12+d2*3+c1;
//    int off2 = (ds*3+c2)*12+d1*3+c1;
//
//    Ty* tp1 = &p1[massi][(off1*NTt+ti) * Nxyz];
//    Ty* tp2 = &p2[massi][(off2*NTt+ti) * Nxyz];
//
//    Ty* tr0 = &((res.data())[((massi*16+offdi)*NTt + ti)*Nxyz]);
//
//    #if USEQACC==1
//    if(invmode == 1){qacc_forNB(i, long(Nxyz),{ tr0[i] += (tp1[i]*qlat::qconj(tp2[i]));});}
//    if(invmode == 0){qacc_forNB(i, long(Nxyz),{ tr0[i] += (tp1[i]*           (tp2[i]));});}
//    #else
//    EAy vp1(tp1,Nxyz);
//    EAy vp2(tp2,Nxyz);
//    EAy vr0(tr0,Nxyz);
//    if(invmode == 1)vr0 += (vp1*vp2.conjugate());
//    if(invmode == 0)vr0 += (vp1*vp2            );
//    #endif
//  }
//  qacc_barrier(dummy);
//  }
//}

//template <typename Ta >
//void meson_vectorE(EigenMTa &prop1, EigenMTa &prop2, ga_M &ga1,ga_M &ga2,
//        EigenVTa &res, qlat::fft_desc_basic &fd,int clear=1, int invmode=1){
//  TIMER("Meson_vectorE");
//  check_prop_size(prop1);check_prop_size(prop2);
//  ///////check_prop_size(prop1);check_prop_size(prop2);
//  int  NTt  = fd.Nv[3];
//  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
//  int  nmass = prop1.size()/(12*12*NTt);
//  if(nmass == 0){res.resize(0);return;}
//
//  if(clear == 1){ini_resE(res,nmass,fd);}
//
//  if(res.size()%NTt !=0 or res.size()==0){print0("Size of res wrong. \n");qassert(false);}
//  qassert(prop1.size() == prop2.size());
//
//  for(int d2=0;d2<4;d2++)
//  for(int c2=0;c2<3;c2++)
//  for(int d1=0;d1<4;d1++)
//  for(int c1=0;c1<3;c1++)
//  {
//  //#pragma omp parallel for
//  for(int ji=0;ji<nmass*NTt;ji++)
//  {
//    int massi = ji/NTt;
//    int ti    = ji%NTt;
//
//    int off1 = massi*12*12 + (d2*3+c2)*12+ga1.ind[d1]*3+c1;
//    int off2 = massi*12*12 + (ga2.ind[d2]*3+c2)*12+d1*3+c1;
//
//    Ta g_tem = ga2.g[d2]*ga1.g[d1];
//
//    Ta* tp1 = prop1[off1*NTt+ti].data();
//    Ta* tp2 = prop2[off2*NTt+ti].data();
//    Ta* tr0 = &((res.data())[(massi*NTt + ti)*Nxyz]);
//
//    #if USEQACC==1
//    if(invmode == 1){qacc_forNB(i, long(Nxyz),{ tr0[i] += (tp1[i]*qlat::qconj(tp2[i]) * g_tem);});}
//    if(invmode == 0){qacc_forNB(i, long(Nxyz),{ tr0[i] += (tp1[i]*           (tp2[i]) * g_tem);});}
//    #else
//    EAa vp1(tp1,Nxyz);
//    EAa vp2(tp2,Nxyz);
//    EAa vr0(tr0,Nxyz);
//    if(invmode == 1)vr0 += (vp1*vp2.conjugate() * g_tem);
//    if(invmode == 0)vr0 += (vp1*vp2             * g_tem);
//    #endif
//  }
//  qacc_barrier(dummy);
//  }
//
//}

#ifdef QLAT_USE_ACC
template <typename Ty, int invmode, int bfac, int Blocks>
__global__ void meson_vectorEV_global(Ty** p1, Ty** p2, Ty* resP, 
  char** gPP, unsigned char** oPP, const int* ivP,
  const int nmass, const int NTt, const long Nxyz, const int Ngv){
  const unsigned long gi =  blockIdx.x;
  const unsigned int tid = threadIdx.y*blockDim.x+ threadIdx.x;
  const long Ntotal = nmass * NTt * Nxyz;
  const long Nbfac  = Ntotal/bfac;
  __shared__ Ty P1[bfac*12*12];
  __shared__ Ty P2[bfac*12*12];

  if(gi*bfac < Ntotal){

  int bi0= 0;int dc = 0;
  int ji = 0;int massi = 0;int ti = 0;

  int jobN = bfac*12*12;
  unsigned int off = tid;
  while(off < jobN){
    bi0= off/(12*12);
    dc = off%(12*12);
    ji    = (bi0*Nbfac + gi)/Nxyz;
    massi = ji/NTt;
    ti    = ji%NTt;
    long ixyz = (bi0*Nbfac + gi)%Nxyz;
    P1[dc*bfac + bi0] = p1[(massi*12*12 + dc)*NTt + ti][ixyz];
    P2[dc*bfac + bi0] = p2[(massi*12*12 + dc)*NTt + ti][ixyz];
    off += Blocks;
  }
  __syncthreads();

  int ini = 0;
  int dv = 0;
  unsigned int MAX = 0;

  const int bfacC = bfac;
  const int Nth   = Blocks/bfac;
  const unsigned int Each  =  4*bfacC;
  const unsigned int GROUP = (Blocks/bfac)*Each;
  unsigned char* s0 = NULL;
           char* s1 = NULL;

  const int bi =  threadIdx.y;
  const int ai =  threadIdx.x;

  const int ba = (threadIdx.y/bfacC)*bfacC + 0;
  const int aa = (threadIdx.y%bfacC)*blockDim.x + threadIdx.x;

  __shared__ Ty buf[bfacC*Blocks];
  __shared__ unsigned char pos[3*GROUP];
  __shared__ char g0[2*GROUP];

  for(int iv=0;iv<Ngv;iv++)
  {
    for(int bz=0;bz<bfacC;bz++){buf[bz*Blocks + tid] = 0;}
    MAX = ivP[iv];
    jobN = (MAX + GROUP - 1 )/GROUP;
    ini = 0; dv = GROUP;
    for(int ji=0;ji<jobN;ji++){
      ////if(ini >= MAX){continue;}
      if(ini + dv >= MAX){dv = MAX - ini;}
      s0 = &(oPP[iv][ini*2]);
      s1 = &(gPP[iv][ini*2]);

      off = tid;
      while(off < dv*2){pos[off] = s0[off];off += Blocks;}
      off = tid;
      while(off < dv*2){g0[off]  = s1[off];off += Blocks;}
      __syncthreads();

      off = aa;
      while(off < dv){
        const Ty* t1 = &P1[(pos[off*2+0])*bfac + ba];
        const Ty* t2 = &P2[(pos[off*2+1])*bfac + ba];
        Ty gtem = Ty(g0[off*2+0], g0[off*2+1]);
        if(bfacC == 1){
          if(invmode == 0){ buf[aa*bfac+ba] += (t1[0]*           (t2[0]) * gtem); }
          if(invmode == 1){ buf[aa*bfac+ba] += (t1[0]*qlat::qconj(t2[0]) * gtem); }
        }else{
          Ty* b0 = &buf[aa*bfac+ba];
          for(int bz=0;bz<bfacC;bz++){
            if(invmode == 0){b0[bz] += (t1[bz]*           (t2[bz]) * gtem); }
            if(invmode == 1){b0[bz] += (t1[bz]*qlat::qconj(t2[bz]) * gtem); }
          }
        }
        off += Nth*bfacC;
      }
      __syncthreads();

      ini += dv;
    }

    for(int atem=1;atem<bfacC;atem++){buf[(0*Nth+ai)*bfac+bi] += buf[(atem*Nth+ai)*bfac+bi];} __syncthreads();

    if(Nth >=256){if(ai <128){buf[ai*bfac + bi] += buf[(ai+128)*bfac + bi];}__syncthreads();}
    if(Nth >=128){if(ai < 64){buf[ai*bfac + bi] += buf[(ai+ 64)*bfac + bi];}__syncthreads();}
    if(Nth >= 64){if(ai < 32){buf[ai*bfac + bi] += buf[(ai+ 32)*bfac + bi];}__syncthreads();}
    if(Nth >= 32){if(ai < 16){buf[ai*bfac + bi] += buf[(ai+ 16)*bfac + bi];}__syncthreads();}
    if(Nth >= 16){if(ai <  8){buf[ai*bfac + bi] += buf[(ai+  8)*bfac + bi];}__syncthreads();}
    if(Nth >=  8){if(ai <  4){buf[ai*bfac + bi] += buf[(ai+  4)*bfac + bi];}__syncthreads();}
    if(Nth >=  4){if(ai <  2){buf[ai*bfac + bi] += buf[(ai+  2)*bfac + bi];}__syncthreads();}

    if(ai == 0){
      resP[iv*Ntotal + bi*Nbfac + gi] += (buf[bi] + buf[bfac+bi]);
      //if(clear == 0){resP[iv*Ntotal + bi*Nbfac + gi] += (buf[bi] + buf[bfac+bi]);}
      //if(clear == 1){resP[iv*Ntotal + bi*Nbfac + gi]  = (buf[bi] + buf[bfac+bi]);}
    }
    __syncthreads();

  }

  }

}
#endif

template <typename Ty, int invmode, int bfac>
void meson_vectorEV_kernel(Ty** p1, Ty** p2, Ty* resP, 
  char** gPP, unsigned char** oPP, const int* ivP,
  const int nmass, const int NTt, const long Nxyz, const int Ngv){
  long Ntotal  = nmass*NTt*Nxyz;
  if(Ntotal % bfac != 0){abort_r("Please correct your bfac! \n");}
  long Nbfac = Ntotal/bfac;
  #if USEGLOBAL==1
  const int nt =  8;
  const int Blocks = nt*bfac;
  dim3 dimBlock(    nt, bfac, 1);
  dim3 dimGrid(  Nbfac,  1, 1);
  meson_vectorEV_global<Ty, invmode, bfac, Blocks><<<dimGrid, dimBlock>>>(p1, 
        p2, resP, gPP, oPP, ivP, nmass, NTt, Nxyz, Ngv);
  qacc_barrier(dummy);
  #else
  if((nmass*NTt) % bfac != 0){abort_r("Please correct your bfac! \n");}
  qacc_for(gi, Nbfac ,
  {
    Ty buf[bfac+1];
    Ty P1[bfac*12*12+1];
    Ty P2[bfac*12*12+1];

    long ixyz = gi%Nxyz;
    int ji    = (gi/Nxyz)*bfac + 0;
    int massi = ji/NTt;
    int ti    = ji%NTt;
    const long offR0 = (massi*NTt + ti)*Nxyz + ixyz;

    for(int bi=0;bi<bfac;bi++)
    {
      massi = (ji+bi)/NTt;
      ti    = (ji+bi)%NTt;

      for(int dc=0;dc<12*12;dc++){
        P1[dc*bfac + bi] = p1[(massi*12*12 + dc)*NTt + ti][ixyz];
        P2[dc*bfac + bi] = p2[(massi*12*12 + dc)*NTt + ti][ixyz];
      }
    }

    for(int iv=0;iv<Ngv;iv++){
      for(int bi=0;bi<bfac;bi++){buf[bi] = 0;}
      for(int off=0;off<ivP[iv];off++)
      {
        const Ty* t1 = &P1[(oPP[iv][off*2+0])*bfac];
        const Ty* t2 = &P2[(oPP[iv][off*2+1])*bfac];
        const Ty gtem = Ty(gPP[iv][off*2+0], gPP[iv][off*2+1]);
        for(int bi=0;bi<bfac;bi++)
        {
          if(invmode == 1){ buf[bi] += (t1[bi]*qlat::qconj(t2[bi]) * gtem); }
          if(invmode == 0){ buf[bi] += (t1[bi]*           (t2[bi]) * gtem); }
        }
      }

      long offR = iv * Ntotal;
      Ty* r0 = &resP[offR + offR0];
      for(int bi=0;bi<bfac; bi++){
        r0[bi*Nxyz] += buf[bi];
        //if(clear == 0){r0[bi*Nxyz] += buf[bi];}
        //if(clear == 1){r0[bi*Nxyz]  = buf[bi];}
      }
    }
  });
  #endif

}

template <typename Ty>
void meson_vectorEV(Ty** p1, Ty** p2, Ty* resP,  int nmass, 
    std::vector<ga_M > &ga1V, std::vector<ga_M > &ga2V,
    qlat::fft_desc_basic &fd, int clear=1, int invmode=1){
  TIMER("meson_vectorEV");
  ///////check_prop_size(prop1);check_prop_size(prop2);
  int  NTt  = fd.Nv[3];
  long Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  qassert(ga1V.size() == ga2V.size());
  int Ngv = ga1V.size();
  if(clear == 1){zero_Ty(resP, Ngv*nmass*NTt*Nxyz , 1);}

  qlat::vector_acc<Ty > gMap;
  qlat::vector_acc<unsigned char > IMap;
  gMap.resize(Ngv*4*2);IMap.resize(Ngv*4*2);
  for(int iv=0;iv<Ngv;iv++){
    for(int i=0;i<4;i++){
      int j = iv*4 + i;
      gMap[0*Ngv*4+j] = ga1V[iv].g[i];
      gMap[1*Ngv*4+j] = ga2V[iv].g[i];
      IMap[0*Ngv*4+j] = ga1V[iv].ind[i];
      IMap[1*Ngv*4+j] = ga2V[iv].ind[i];
    }
  }

  Ty* gC_P = gMap.data();
  unsigned char*      gI_P = IMap.data();

  #if USEKERNEL==1
  std::vector<std::vector<char > > giEL;giEL.resize(Ngv);
  std::vector<std::vector<unsigned char   > > oiL ;oiL.resize(Ngv );

  ////reformulate index
  for(int iv=0;iv<Ngv;iv++){
  oiL[iv].resize(0);
  giEL[iv].resize(0);
  const int j1 = 0*Ngv*4 + iv*4 ;
  const int j2 = 1*Ngv*4 + iv*4 ;
  const Ty* gC1 = &(gC_P[j1]);
  const Ty* gC2 = &(gC_P[j2]);
  const unsigned char* gI1 = &(gI_P[j1]);
  const unsigned char* gI2 = &(gI_P[j2]);
  for(int d2=0;d2<4;d2++)
  for(int c2=0;c2<3;c2++)
  for(int d1=0;d1<4;d1++)
  for(int c1=0;c1<3;c1++)
  {
    const char off1 = (d2*3+c2)*12+gI1[d1]*3+c1;
    const char off2 = (gI2[d2]*3+c2)*12+d1*3+c1;
    const Ty g_tem = gC2[d2]*gC1[d1];
    const double norm = qlat::qnorm(g_tem);
    if(norm < 1e-20)continue;

    oiL[iv].push_back(off1);
    oiL[iv].push_back(off2);
    giEL[iv].push_back(char(g_tem.real()));
    giEL[iv].push_back(char(g_tem.imag()));
  }
  }

  std::vector<qlat::vector_gpu<char > > giEG;giEG.resize(Ngv);
  for(int iv=0;iv<Ngv;iv++){giEG[iv].copy_from(giEL[iv]);}
  qlat::vector_acc<char* > gP = EigenM_to_pointers(giEG);
  char** gPP = gP.data();

  std::vector<qlat::vector_gpu<unsigned char   > > oiG ; oiG.resize(Ngv);
  for(int iv=0;iv<Ngv;iv++){oiG[iv].copy_from(oiL[iv]);}
  qlat::vector_acc<unsigned char* > oP = EigenM_to_pointers(oiG);
  unsigned char** oPP = oP.data();

  qlat::vector_acc<int > iv_size;iv_size.resize(Ngv);
  for(int iv=0;iv<Ngv;iv++){iv_size[iv] = giEL[iv].size()/2;}
  int*  ivP = iv_size.data();

  //////long Ntotal  = nmass*NTt*Nxyz;
  //int mode = 0;mode = invmode*2 + clear;
  //const int BFACG = BFACG_SHARED;


  #if USEGLOBAL==1
  const int BFACG = BFACG_SHARED;
  #else
  const int BFACG = 1;
  #endif

  if(invmode==0)meson_vectorEV_kernel<Ty,0, BFACG>(p1, p2, resP, gPP, oPP, ivP, nmass, NTt, Nxyz, Ngv);
  if(invmode==1)meson_vectorEV_kernel<Ty,1, BFACG>(p1, p2, resP, gPP, oPP, ivP, nmass, NTt, Nxyz, Ngv);
  //if(mode==2)meson_vectorEV_kernel<Ty,1,0, BFACG>(p1, p2, resP, gPP, oPP, ivP, nmass, NTt, Nxyz, Ngv);
  //if(mode==3)meson_vectorEV_kernel<Ty,1,1, BFACG>(p1, p2, resP, gPP, oPP, ivP, nmass, NTt, Nxyz, Ngv);
  qacc_barrier(dummy);
  #endif

  #if USEKERNEL==0
  for(int iv=0;iv<Ngv;iv++){
    int j1 = 0*Ngv*4 + iv*4 ;
    int j2 = 1*Ngv*4 + iv*4 ;
    Ty* gC1 = &(gC_P[j1]);
    Ty* gC2 = &(gC_P[j2]);
    unsigned char* gI1 = &(gI_P[j1]);
    unsigned char* gI2 = &(gI_P[j2]);
    long offR = iv*nmass*NTt * Nxyz;
    for(int d2=0;d2<4;d2++)
    for(int c2=0;c2<3;c2++)
    for(int d1=0;d1<4;d1++)
    for(int c1=0;c1<3;c1++)
    {
    #pragma omp parallel for
    for(int ji=0;ji<nmass*NTt;ji++)
    {
      int massi = ji/NTt;
      int ti    = ji%NTt;

      int off1 = massi*12*12 + (d2*3+c2)*12+gI1[d1]*3+c1;
      int off2 = massi*12*12 + (gI2[d2]*3+c2)*12+d1*3+c1;

      Ty g_tem = gC2[d2]*gC1[d1];

      Ty* tp1 = p1[off1*NTt+ti];
      Ty* tp2 = p2[off2*NTt+ti];
      Ty* tr0 = &(resP[offR + (massi*NTt + ti)*Nxyz]);

      #if USEQACC==1
      if(invmode == 1){qacc_forNB(i, long(Nxyz),{ tr0[i] += (tp1[i]*qlat::qconj(tp2[i]) * g_tem);});}
      if(invmode == 0){qacc_forNB(i, long(Nxyz),{ tr0[i] += (tp1[i]*           (tp2[i]) * g_tem);});}
      #else
      EAy vp1(tp1,Nxyz);
      EAy vp2(tp2,Nxyz);
      EAy vr0(tr0,Nxyz);
      if(invmode == 1)vr0 += (vp1*vp2.conjugate() * g_tem);
      if(invmode == 0)vr0 += (vp1*vp2             * g_tem);
      #endif

    }
    qacc_barrier(dummy);
    }
  }
  #endif

}

//template <typename Ta>
//void meson_vectorEV(EigenMTa &prop1, EigenMTa &prop2, EigenVTa &res, std::vector<ga_M > &ga1V, std::vector<ga_M > &ga2V,
//        qlat::fft_desc_basic &fd, int clear=1, int invmode=1){
//  check_prop_size(prop1);check_prop_size(prop2);
//  int  nmass = prop1.size()/(12*12*fd.Nv[3]);
//  if(nmass == 0){res.resize(0);return;}
//  int Ngv = ga1V.size();
//  long resL = Ngv * nmass * fd.Nv[0]*fd.Nv[1]*fd.Nv[2] * fd.Nv[3];
//  if(clear == 1){if(res.size()!= resL){res.resize(resL);}}
//
//  if(res.size() != resL){print0("Size of res wrong. \n");qassert(false);}
//  qassert(prop1.size() == prop2.size());
//
//  qlat::vector_acc<Ta* > prop1P = EigenM_to_pointers(prop1);
//  qlat::vector_acc<Ta* > prop2P = EigenM_to_pointers(prop2);
//
//  Ta** p1 = prop1P.data();
//  Ta** p2 = prop2P.data();
//  Ta* resP = res.data();
//  meson_vectorEV(p1, p2, resP, nmass, ga1V, ga2V, fd, clear, invmode);
//}

template <typename Ty >
void meson_vectorEV(EigenTy& prop1, EigenTy& prop2, qlat::vector_gpu<Ty > &res
  ,std::vector<ga_M > &ga1V, std::vector<ga_M > &ga2V,
  qlat::fft_desc_basic &fd, int clear=1, int invmode=1)
{
  check_prop_size(prop1, fd);check_prop_size(prop2, fd);
  int  nmass = prop1.size();
  if(nmass == 0){res.resize(0);return;}
  int Ngv = ga1V.size();
  const unsigned long resL = Ngv * nmass * fd.Nv[0]*fd.Nv[1]*fd.Nv[2] * fd.Nv[3];
  const long Nxyz= fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  if(clear == 1){if(res.size()!= resL){res.resize(resL);}}

  if(res.size() != resL){print0("Size of res wrong. \n");qassert(false);}
  qassert(prop1.size() == prop2.size());
  for(int mi=0;mi<nmass;mi++)
  {
    qassert(prop1[mi].size() == 12 * 12 * fd.Nvol);
    qassert(prop2[mi].size() == 12 * 12 * fd.Nvol);
  }

  qlat::vector_acc<Ty* > prop1P = EigenM_to_pointers(prop1, Nxyz);
  qlat::vector_acc<Ty* > prop2P = EigenM_to_pointers(prop2, Nxyz);

  Ty** p1 = prop1P.data();
  Ty** p2 = prop2P.data();
  Ty* resP = res.data();
  meson_vectorEV(p1, p2, resP, nmass, ga1V, ga2V, fd, clear, invmode);
}


template<typename Td, typename Ta>
void meson_corrE(std::vector<Propagator4dT<Td > > &pV1, std::vector<Propagator4dT<Td >> &pV2,  ga_M &ga1, ga_M &ga2,
  qlat::vector_acc<Ta > &res, qlat::fft_desc_basic &fd,int clear=1,const Coordinate& mom = Coordinate()){
  ///int NTt  = fd.Nv[3];
  ///LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  int nmass = pV1.size();
  ///int nt = fd.nt;

  qlat::vector_acc<Ta > resE;
  ini_resE(resE,nmass,fd);

  meson_vectorE(pV1,pV2,ga1,ga2,resE,fd,1);

  vec_corrE(resE,res,fd, clear, mom);
}

template<typename Ty>
void meson_corrE(std::vector<qpropT > &prop1, std::vector<qpropT > &prop2,  ga_M &ga1, ga_M &ga2,
  qlat::vector_acc<Ty >& res, int clear=1,const Coordinate& mom = Coordinate()){
  qlat::vector_acc<Ty > resE;

  const qlat::Geometry &geo = prop1[0].geo();
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);

  meson_vectorE(prop1,prop2,ga1,ga2,resE,1);
  vec_corrE(resE, res, fd, clear, mom);
}

template<typename Td>
void print_pion(std::vector<Propagator4dT<Td > > &prop1, std::vector<Propagator4dT<Td > > &prop2, const std::string& tag=std::string(""), double factor = 1.0){
  ga_matrices_cps   ga_cps;
  const qlat::Geometry &geo = prop1[0].geo();
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);

  qlat::vector_acc<qlat::ComplexT<Td > > resC;

  meson_corrE(prop1, prop2, ga_cps.ga[0][0], ga_cps.ga[0][0], resC, fd);

  int nv = resC.size()/fd.nt;
  for(int iv=0;iv<nv;iv++)
  for(int t=0;t<fd.nt;t++)
  {
    qlat::ComplexT<Td > v = resC[iv*fd.nt + t] * qlat::ComplexT<Td >(factor, 0.0);
    print0("%s iv %d, t %d, v %.6e %.6e \n", tag.c_str(), iv, t, v.real(), v.imag());
  }
}


//template <typename Ta>
//void meson_corr_write(Propagator4dT<Ta > &propVa, Propagator4dT<Ta > &propVb, int pos, std::vector<double > &write, int offw, const Geometry &geo, int a=0, int b=0, int c=0 , int d=0){
//  print_mem_info();
//  fft_desc_basic fd(geo);
//  //qlat::vector<int > nv, Nv, mv;
//  //geo_to_nv(geo, nv, Nv, mv);
//  int nt = fd.nt;
//
//  ///char output[500];
//  ///sprintf(output,   out_n.c_str());
//  ///print0("output %s \n", output);
//
//  EigenMTa propa,propb;
//  copy_prop4d_to_propE(propa, propVa, fd);
//  copy_prop4d_to_propE(propb, propVb, fd);
//  ////copy_propE(propVa, propa, fd );
//  ////copy_propE(propVb, propb, fd );
//
//  ///Coordinate xg1;
//  ///xg1[0] = pos/10000000;xg1[1] = (pos%10000000)/100000;xg1[2] = (pos%100000)/1000;xg1[3] = pos%1000;
//  int t0 = pos%1000;
//
//  EigenVTa res;ga_matrices_cps   ga_cps;
//  meson_corrE(propa, propb, ga_cps.ga[a][b],ga_cps.ga[c][d],  res, fd);
//  ///std::vector<double > write;write.resize(2*nt);
//  for(int ti=0;ti<nt;ti++)
//  {
//    double v0 = res[ti].real();
//    double v1 = res[ti].imag();
//    write[offw + ((ti- t0 +nt)%nt)*2+0]= v0;
//    write[offw + ((ti- t0 +nt)%nt)*2+1]= v1;
//  }
//  ////write_data(write,output);
//
//}

//inline void meson_corr_write(std::string prop_a, std::string prop_b, std::string src_n, std::string out_n, const Geometry &geo, int a=0, int b=0, int c=0 , int d=0){
//  print_mem_info();
//
//  qlat::vector_acc<int > nv, Nv, mv;
//  geo_to_nv(geo, nv, Nv, mv);
//  int nt = nv[3];
//
//  qlat::FieldM<Complexq, 1> noi;
//  noi.init(geo);
//  Propagator4d propVa;propVa.init(geo);
//  Propagator4d propVb;propVb.init(geo);
//
//  char prop_na[500],prop_nb[500],noi_name[500];
//  char output[500];
//  sprintf(prop_na, "%s",prop_a.c_str() );
//  sprintf(prop_nb, "%s",prop_b.c_str() );
//
//  sprintf(noi_name ,"%s",src_n.c_str()  );
//  sprintf(output,   "%s",out_n.c_str());
//
//  print0("Noise %s \n",noi_name);
//  print0("Prop  %s %s \n",prop_na, prop_nb);
//  print0("output %s \n", output);
//
//  qlat::set_zero(noi);
//  load_gwu_noi(noi_name,noi);
//  load_gwu_prop(prop_na, propVa);
//  if(prop_a == prop_b){propVb = propVa;}
//  else{load_gwu_prop(prop_nb, propVb);}
//  
//  ////std::vector<qlat::vector_acc<Complexq > > propa,propb;
//  std::vector<qprop > propa, propb;
//  propa.resize(1);propa[0].init(geo);
//  propb.resize(1);propb[0].init(geo);
//  prop4d_to_qprop(propa[0], propVa);
//  prop4d_to_qprop(propb[0], propVb);
//
//  Coordinate pos;Coordinate off_L;
//  check_noise_pos(noi, pos, off_L);
//
//  ////Coordinate xg1;
//  ////xg1[0] = pos/10000000;xg1[1] = (pos%10000000)/100000;xg1[2] = (pos%100000)/1000;xg1[3] = pos%1000;
//
//  qlat::vector_acc<qlat::Complex > res;ga_matrices_cps   ga_cps;
//  meson_corrE(propa, propb, ga_cps.ga[a][b],ga_cps.ga[c][d],  res);
//  std::vector<double > write;write.resize(2*nt);
//  for(unsigned int ti=0;ti<write.size()/2;ti++){
//    double v0 = res[ti].real();
//    double v1 = res[ti].imag();
//    write[((ti-pos[3]+nt)%nt)*2+0]= v0;
//    write[((ti-pos[3]+nt)%nt)*2+1]= v1;
//  }
//
//  write_data(write,output);
//
//}

//template <typename Ta>
//void print_meson(Propagator4dT<Ta > &propVa, Propagator4dT<Ta > &propVb, std::string tag=std::string(""), int a=0, int b=0, int c=0 , int d=0){
//  const qlat::Geometry &geo = propVa.geo();
//  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);
//
//  int nt = fd.nt;
//
//  std::vector<qprop > propa, propb;
//  propa.resize(1);propa[0].init(geo);
//  propb.resize(1);propb[0].init(geo);
//  prop4d_to_qprop(propa[0], propVa);
//  prop4d_to_qprop(propb[0], propVb);
//
//  ////copy_propE(propVa, propa, fd );
//  ////copy_propE(propVb, propb, fd );
//
//  qlat::vector_acc<qlat::Complex > res;ga_matrices_cps   ga_cps;
//  meson_corrE(propa, propb, ga_cps.ga[a][b],ga_cps.ga[c][d],  res);
//  for(int ti=0;ti<nt;ti++)
//  {
//    double v0 = res[ti].real();
//    double v1 = res[ti].imag();
//    print0("%s ti %5d , v  %.8e   %.8e \n", tag.c_str(), ti, v0, v1);
//  }
//}

//template<typename Ta>
//void print_pion(qlat::FieldM<Ta, 12*12 >& propM, const std::string& tag=std::string(""), double factor = 1.0){
//  const Geometry& geo = propM.geo();
//  fft_desc_basic fd(geo);
//
//  Propagator4dT<Ta > prop4d;prop4d.init(geo);
//  std::vector<qlat::vector_acc<Ta > > propE;
//
//  copy_noise_to_prop(propM, prop4d, 1);
//  copy_prop4d_to_propE(propE, prop4d, fd);
//  ////copy_propE(prop4d, propE, fd);
//
//  ga_matrices_cps   ga_cps;
//  EigenVTa res;EigenVTa corr;
//  meson_vectorE(propE, propE, ga_cps.ga[0][0], ga_cps.ga[0][0],res, fd);
//
//  vec_corrE(res, corr, fd, 1 );
//
//  int nv = corr.size()/fd.nt;
//  for(int iv=0;iv<nv;iv++)
//  for(int t=0;t<fd.nt;t++)
//  {
//    Ta v = corr[iv*fd.nt + t] * Ta(factor, 0.0);
//    print0("%s iv %d, t %d, v %.6e %.6e \n", tag.c_str(), iv, t, v.real(), v.imag());
//  }
//}

}

#endif

