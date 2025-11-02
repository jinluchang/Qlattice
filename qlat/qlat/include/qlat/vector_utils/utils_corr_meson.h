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
#include "utils_field_gpu.h"

namespace qlat{

/*
  ga1 sink gammas, ga2 src gammas
  invmode 1 : P1_{ab} P2_{ab}^*
  invmode 0 : P1_{ab} P2_{ab}
*/
template <typename Td>
void meson_vectorE(std::vector<Propagator4dT<Td >* > &pV1, std::vector<Propagator4dT<Td >* > &pV2, ga_M &ga1,ga_M &ga2,
        qlat::ComplexT<Td >* res, qlat::fft_desc_basic &fd, Int invmode=1){
  TIMER("Meson_vectorE");
  Qassert(fd.order_ch == 0);
  ///////check_prop_size(prop1);check_prop_size(prop2);
  Int  NTt  = fd.Nv[3];
  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  Int  nmass = pV1.size();

  //if(nmass == 0){res.resize(0);return;}
  //if(clear == 1){ini_resE(res,nmass,fd);}
  //if(res.size()%NTt !=0 or res.size()==0){qmessage("Size of res wrong. \n");Qassert(false);}
  Qassert(pV1.size() == pV2.size());

  for(Int mi=0;mi<nmass;mi++)
  {
  Propagator4dT<Td >& pL1 = *pV1[mi];
  Propagator4dT<Td >& pL2 = *pV2[mi];

  qacc_for(isp, Long(pV1[0]->geo().local_volume()),{ 
    Int ti = isp/Nxyz;
    Int xi = isp%Nxyz;
      qlat::ComplexT<Td > pres;pres = 0.0;
      const qlat::WilsonMatrixT<Td>& p1 =  pL1.get_elem_offset(isp);
      const qlat::WilsonMatrixT<Td>& p2 =  pL2.get_elem_offset(isp);

      for(Int d1=0;d1<4;d1++)
      for(Int c1=0;c1<3;c1++)
      for(Int d2=0;d2<4;d2++)
      {
      const qlat::ComplexT<Td > g_tem = ga2.g[d2]*ga1.g[d1];
      for(Int c2=0;c2<3;c2++)
      {
        if(invmode == 1){
          pres += g_tem * 
            p1(ga1.ind[d1]*3+c1,d2*3+c2) * qlat::qconj(p2(d1*3+c1,ga2.ind[d2]*3+c2)) ;
        }
        if(invmode == 0){
          pres += g_tem * 
            p1(ga1.ind[d1]*3+c1,d2*3+c2) *            (p2(d1*3+c1,ga2.ind[d2]*3+c2)) ;
        }

      }
      }
      res[(mi*NTt + ti)*Nxyz + xi%Nxyz] += pres;
  });
  }

}

template <typename Td>
void meson_vectorE(std::vector<Propagator4dT<Td > > &pV1, std::vector<Propagator4dT<Td > > &pV2, ga_M &ga1,ga_M &ga2,
        qlat::vector<qlat::ComplexT<Td > > &res, qlat::fft_desc_basic &fd,Int clear=1, Int invmode=1){
  std::vector<Propagator4dT<Td >* > p1;
  std::vector<Propagator4dT<Td >* > p2;
  Qassert(pV1.size() == pV2.size());
  const Int Nprop = pV1.size();
  p1.resize(Nprop);
  p2.resize(Nprop);
  for(Int i=0;i<Nprop;i++){
    p1[i] = &pV1[i];
    p2[i] = &pV2[i];
  }

  const Int  NTt  = fd.Nv[3];
  // const LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  const Int  nmass = pV1.size();
  if(nmass == 0){res.resize(0);return;}
  if(clear == 1){ini_resE(res,nmass,fd);}
  if(res.size()%NTt !=0 or res.size()==0){qmessage("Size of res wrong. \n");Qassert(false);}

  meson_vectorE(p1, p2, ga1, ga2, res.data(), fd, invmode);
}

template <typename Td>
void meson_vectorE(Propagator4dT<Td > &pV1, Propagator4dT<Td > &pV2, ga_M &ga1,ga_M &ga2,
        qlat::ComplexT<Td >* res, Int clear=1, Int invmode=1){
  Qassert(pV1.initialized and pV2.initialized);
  std::vector<Propagator4dT<Td >* > p1;
  std::vector<Propagator4dT<Td >* > p2;
  const Int Nprop = 1;
  p1.resize(Nprop);
  p2.resize(Nprop);
  p1[0] = &pV1;
  p2[0] = &pV2;
  const qlat::Geometry& geo = pV1.geo();
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);
  if(clear == 1){
    const Int  NTt  = fd.Nv[3];
    const LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
    zero_Ty(res, NTt * Nxyz, 1);
  }

  meson_vectorE(p1, p2, ga1, ga2, res, fd, invmode);
}

template <typename Ty >
void meson_vectorE(std::vector<qpropT >& prop1, std::vector<qpropT >& prop2, ga_M &ga1,ga_M &ga2,
        qlat::vector<Ty > &res, Int clear=1, Int invmode=1, const Ty factor = Ty(1.0, 0.0)){
  TIMER("Meson_vectorE");
  const qlat::Geometry& geo = prop1[0].geo();
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);

  Int  NTt  = fd.Nv[3];
  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  Int  nmass = prop1.size();  ////(12*12*NTt)
  if(nmass == 0){res.resize(0);return;}
  if(clear == 1){ini_resE(res, nmass, fd);}
  if(res.size()%NTt != 0 or res.size() == 0){qmessage("Size of res wrong. \n");Qassert(false);}

  Qassert(prop1.size() == prop2.size());
  qlat::vector<Ty* > p1 = FieldM_to_pointers(prop1);
  qlat::vector<Ty* > p2 = FieldM_to_pointers(prop2);

  for(Int d2=0;d2<4;d2++)
  for(Int c2=0;c2<3;c2++)
  for(Int d1=0;d1<4;d1++)
  for(Int c1=0;c1<3;c1++)
  {
  //#pragma omp parallel for
  for(Int ji=0;ji<nmass*NTt;ji++)
  {
    Int massi = ji/NTt;
    Int ti    = ji%NTt;

    Int off1 = (d2*3+c2)*12+ga1.ind[d1]*3+c1;
    Int off2 = (ga2.ind[d2]*3+c2)*12+d1*3+c1;

    const Ty g_tem = ga2.g[d2]*ga1.g[d1];

    Ty* tp1 = &p1[massi][(off1*NTt+ti) * Nxyz];
    Ty* tp2 = &p2[massi][(off2*NTt+ti) * Nxyz];

    Ty* tr0 = &((res.data())[(massi*NTt + ti)*Nxyz]);

    #if USEQACC==1
    if(invmode == 1){qacc_forNB(i, Long(Nxyz),{ tr0[i] += factor * (tp1[i]*qlat::qconj(tp2[i]) * g_tem);});}
    if(invmode == 0){qacc_forNB(i, Long(Nxyz),{ tr0[i] += factor * (tp1[i]*           (tp2[i]) * g_tem);});}
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

#ifdef QLAT_USE_ACC
template <typename Ty, Int invmode, Int bfac, Int Blocks>
__global__ void meson_vectorEV_global(Ty** p1, Ty** p2, Ty** resP, 
  int8_t** gPP, uint8_t** oPP, const Int* ivP,
  const Int nmass, const Int NTt, const Long Nxyz, const Int Ngv)
{
  const unsigned long gi =  blockIdx.x;
  const unsigned int tid = threadIdx.y*blockDim.x+ threadIdx.x;
  const Long Ntotal = nmass * NTt * Nxyz;
  const Long Nbfac  = Ntotal/bfac;
  __shared__ Ty P1[bfac*12*12];
  __shared__ Ty P2[bfac*12*12];

  if(gi*bfac < Ntotal){

  Long bi0= 0;int dc = 0;
  Int  ji = 0;int massi = 0;int ti = 0;

  Int jobN = bfac*12*12;
  unsigned int off = tid;
  while(off < jobN){
    bi0= off/(12*12);
    dc = off%(12*12);
    ji    = (bi0*Nbfac + gi)/Nxyz;
    massi = ji/NTt;
    ti    = ji%NTt;
    Long ixyz = (bi0*Nbfac + gi)%Nxyz;
    P1[dc*bfac + bi0] = p1[(massi*12*12 + dc)*NTt + ti][ixyz];
    P2[dc*bfac + bi0] = p2[(massi*12*12 + dc)*NTt + ti][ixyz];
    off += Blocks;
  }
  __syncthreads();


  Int ini = 0;
  Int dv = 0;
  unsigned int MAX = 0;

  const Int bfacC = bfac;
  const Int Nth   = Blocks/bfac;
  const unsigned int Each  =  4*bfacC;
  const unsigned int GROUP = (Blocks/bfac)*Each;
  uint8_t* s0 = NULL;
    int8_t* s1 = NULL;

  const Int bi =  threadIdx.y;
  const Int ai =  threadIdx.x;

  const Int ba = (threadIdx.y/bfacC)*bfacC + 0;
  const Int aa = (threadIdx.y%bfacC)*blockDim.x + threadIdx.x;

  ///Long off v = iv*Ntotal + bi*Nbfac + gi;
  Long idx_res = bi*Nbfac + gi;
  ji    = idx_res / Nxyz;
  massi = ji/NTt;
  ti    = ji%NTt;
  Long off_res = ti * Nxyz + (idx_res) % Nxyz;

  __shared__ Ty buf[bfacC*Blocks];
  __shared__ uint8_t pos[3*GROUP];
  __shared__  signed  char g0[2*GROUP];

  for(Int iv=0;iv<Ngv;iv++)
  {
    for(Int bz=0;bz<bfacC;bz++){buf[bz*Blocks + tid] = 0;}
    MAX = ivP[iv];
    jobN = (MAX + GROUP - 1 )/GROUP;
    ini = 0; dv = GROUP;
    for(Int ji=0;ji<jobN;ji++){
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
          for(Int bz=0;bz<bfacC;bz++){
            if(invmode == 0){b0[bz] += (t1[bz]*           (t2[bz]) * gtem); }
            if(invmode == 1){b0[bz] += (t1[bz]*qlat::qconj(t2[bz]) * gtem); }
          }
        }
        off += Nth*bfacC;
      }
      __syncthreads();

      ini += dv;
    }

    for(Int atem=1;atem<bfacC;atem++){buf[(0*Nth+ai)*bfac+bi] += buf[(atem*Nth+ai)*bfac+bi];} __syncthreads();

    if(Nth >=256){if(ai <128){buf[ai*bfac + bi] += buf[(ai+128)*bfac + bi];}__syncthreads();}
    if(Nth >=128){if(ai < 64){buf[ai*bfac + bi] += buf[(ai+ 64)*bfac + bi];}__syncthreads();}
    if(Nth >= 64){if(ai < 32){buf[ai*bfac + bi] += buf[(ai+ 32)*bfac + bi];}__syncthreads();}
    if(Nth >= 32){if(ai < 16){buf[ai*bfac + bi] += buf[(ai+ 16)*bfac + bi];}__syncthreads();}
    if(Nth >= 16){if(ai <  8){buf[ai*bfac + bi] += buf[(ai+  8)*bfac + bi];}__syncthreads();}
    if(Nth >=  8){if(ai <  4){buf[ai*bfac + bi] += buf[(ai+  4)*bfac + bi];}__syncthreads();}
    if(Nth >=  4){if(ai <  2){buf[ai*bfac + bi] += buf[(ai+  2)*bfac + bi];}__syncthreads();}

    //iv ops
    if(ai == 0){
      resP[iv*nmass + massi][off_res] += (buf[bi] + buf[bfac+bi]);
      //if(clear == 0){resP[iv*Ntotal + bi*Nbfac + gi] += (buf[bi] + buf[bfac+bi]);}
      //if(clear == 1){resP[iv*Ntotal + bi*Nbfac + gi]  = (buf[bi] + buf[bfac+bi]);}
    }
    __syncthreads();

  }

  }
}
#endif

/*
  vector pointers
  p1, p2 [nvec, ti][Nxyz]
  resP  [nops, nvec ][ti x Nxyz]
*/
template <typename Ty, Int invmode, Int bfac>
void meson_vectorEV_kernel(Ty** p1, Ty** p2, Ty** resP, 
  int8_t** gPP, uint8_t** oPP, const Int* ivP,
  const Int nmass, const Int NTt, const Long Nxyz, const Int Ngv)
{
  Long Ntotal  = nmass*NTt*Nxyz;
  if(Ntotal % bfac != 0){abort_r("Please correct your bfac! \n");}
  Long Nbfac = Ntotal/bfac;
  #if USEGLOBAL==1
  const Int nt =  8;
  const Int Blocks = nt*bfac;
  dim3 dimBlock(    nt, bfac, 1);
  dim3 dimGrid(  Nbfac,  1, 1);
  meson_vectorEV_global<Ty, invmode, bfac, Blocks><<<dimGrid, dimBlock>>>(p1, 
        p2, resP, gPP, oPP, ivP, nmass, NTt, Nxyz, Ngv);
  qacc_barrier(dummy);
  #else
  if((nmass*NTt) % bfac != 0){abort_r("Please correct your bfac! \n");}
  // could try spatial group 4 vectorizations
  qacc_for(gi, Nbfac ,
  {
    Ty buf[bfac+1];
    Ty P1[bfac*12*12+1];
    Ty P2[bfac*12*12+1];

    Long ixyz = gi%Nxyz;
    Int ji    = (gi/Nxyz)*bfac + 0;
    Int massi = ji/NTt;
    Int ti    = ji%NTt;
    //const Long offR0 = (ti)*Nxyz + ixyz;
    const Int j0    = (gi/Nxyz)*bfac + 0;

    for(Int bi=0;bi<bfac;bi++)
    {
      massi = (ji+bi)/NTt;
      ti    = (ji+bi)%NTt;

      for(Int dc=0;dc<12*12;dc++){
        P1[dc*bfac + bi] = p1[(massi*12*12 + dc)*NTt + ti][ixyz];
        P2[dc*bfac + bi] = p2[(massi*12*12 + dc)*NTt + ti][ixyz];
      }
    }

    for(Int iv=0;iv<Ngv;iv++){
      for(Int bi=0;bi<bfac;bi++){buf[bi] = 0;}
      for(Int off=0;off<ivP[iv];off++)
      {
        const Ty* t1 = &P1[(oPP[iv][off*2+0])*bfac];
        const Ty* t2 = &P2[(oPP[iv][off*2+1])*bfac];
        const Ty gtem = Ty(gPP[iv][off*2+0], gPP[iv][off*2+1]);
        for(Int bi=0;bi<bfac;bi++)
        {
          if(invmode == 1){ buf[bi] += (t1[bi]*qlat::qconj(t2[bi]) * gtem); }
          if(invmode == 0){ buf[bi] += (t1[bi]*           (t2[bi]) * gtem); }
        }
      }

      //Long offR = iv * Ntotal;
      // may need more optimization ?
      for(Int bi=0;bi<bfac; bi++){
        const Int ji    = j0 + bi;
        const Int massi = ji/NTt;
        const Int ti    = ji%NTt;
        //const Long offR0 = (ti)*Nxyz + ixyz;

        resP[iv*nmass + massi][ti*Nxyz + ixyz] += buf[bi];
        //r0[bi*Nxyz] += buf[bi];
      }
    }
  });
  #endif
}

template <typename Ty>
void meson_vectorEV(Ty** p1, Ty** p2, Ty** resP,  Int nmass, 
    std::vector<ga_M > &ga1V, std::vector<ga_M > &ga2V,
    qlat::fft_desc_basic &fd, Int clear=1, Int invmode=1){
  TIMER("meson_vectorEV");
  ///////check_prop_size(prop1);check_prop_size(prop2);
  Int  NTt  = fd.Nv[3];
  Long Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  Qassert(ga1V.size() == ga2V.size());
  Int Ngv = ga1V.size();
  if(clear == 1){
    for(Long iv=0;iv<Ngv*nmass;iv++){
      zero_Ty(resP[iv], NTt*Nxyz , 1);
    }
  }

  qlat::vector<Ty > gMap;
  qlat::vector<uint8_t > IMap;
  gMap.resize(Ngv*4*2);IMap.resize(Ngv*4*2);
  for(Int iv=0;iv<Ngv;iv++){
    for(Int i=0;i<4;i++){
      Int j = iv*4 + i;
      gMap[0*Ngv*4+j] = ga1V[iv].g[i];
      gMap[1*Ngv*4+j] = ga2V[iv].g[i];
      IMap[0*Ngv*4+j] = ga1V[iv].ind[i];
      IMap[1*Ngv*4+j] = ga2V[iv].ind[i];
    }
  }

  Ty* gC_P = gMap.data();
  uint8_t*      gI_P = IMap.data();

  #if USEKERNEL==1
  std::vector<std::vector<int8_t > > giEL;giEL.resize(Ngv);
  std::vector<std::vector<uint8_t   > > oiL ;oiL.resize(Ngv );

  ////reformulate index
  for(Int iv=0;iv<Ngv;iv++){
  oiL[iv].resize(0);
  giEL[iv].resize(0);
  const Int j1 = 0*Ngv*4 + iv*4 ;
  const Int j2 = 1*Ngv*4 + iv*4 ;
  const Ty* gC1 = &(gC_P[j1]);
  const Ty* gC2 = &(gC_P[j2]);
  const uint8_t* gI1 = &(gI_P[j1]);
  const uint8_t* gI2 = &(gI_P[j2]);
  for(Int d2=0;d2<4;d2++)
  for(Int c2=0;c2<3;c2++)
  for(Int d1=0;d1<4;d1++)
  for(Int c1=0;c1<3;c1++)
  {
    const int8_t off1 = (d2*3+c2)*12+gI1[d1]*3+c1;
    const int8_t off2 = (gI2[d2]*3+c2)*12+d1*3+c1;
    const Ty g_tem = gC2[d2]*gC1[d1];
    const RealD norm = qlat::qnorm(g_tem);
    if(norm < 1e-20)continue;

    oiL[iv].push_back(off1);
    oiL[iv].push_back(off2);
    giEL[iv].push_back(int8_t(g_tem.real()));
    giEL[iv].push_back(int8_t(g_tem.imag()));
  }
  }

  std::vector<qlat::vector_gpu<int8_t > > giEG;giEG.resize(Ngv);
  for(Int iv=0;iv<Ngv;iv++){giEG[iv].copy_from(giEL[iv]);}
  qlat::vector<int8_t* > gP = EigenM_to_pointers(giEG);
  int8_t** gPP = gP.data();

  std::vector<qlat::vector_gpu<uint8_t   > > oiG ; oiG.resize(Ngv);
  for(Int iv=0;iv<Ngv;iv++){oiG[iv].copy_from(oiL[iv]);}
  qlat::vector<uint8_t* > oP = EigenM_to_pointers(oiG);
  uint8_t** oPP = oP.data();

  qlat::vector<Int > iv_size;iv_size.resize(Ngv);
  for(Int iv=0;iv<Ngv;iv++){iv_size[iv] = giEL[iv].size()/2;}
  Int*  ivP = iv_size.data();

  //////Long Ntotal  = nmass*NTt*Nxyz;
  //int mode = 0;mode = invmode*2 + clear;
  //const Int BFACG = BFACG_SHARED;


  #if USEGLOBAL==1
  const Int BFACG = BFACG_SHARED;
  #else
  const Int BFACG = 1;
  #endif
  if(invmode==0)meson_vectorEV_kernel<Ty,0, BFACG>(p1, p2, resP, gPP, oPP, ivP, nmass, NTt, Nxyz, Ngv);
  if(invmode==1)meson_vectorEV_kernel<Ty,1, BFACG>(p1, p2, resP, gPP, oPP, ivP, nmass, NTt, Nxyz, Ngv);
  //if(mode==2)meson_vectorEV_kernel<Ty,1,0, BFACG>(p1, p2, resP, gPP, oPP, ivP, nmass, NTt, Nxyz, Ngv);
  //if(mode==3)meson_vectorEV_kernel<Ty,1,1, BFACG>(p1, p2, resP, gPP, oPP, ivP, nmass, NTt, Nxyz, Ngv);
  qacc_barrier(dummy);
  #endif

  #if USEKERNEL==0
  for(Int iv=0;iv<Ngv;iv++){
    Int j1 = 0*Ngv*4 + iv*4 ;
    Int j2 = 1*Ngv*4 + iv*4 ;
    Ty* gC1 = &(gC_P[j1]);
    Ty* gC2 = &(gC_P[j2]);
    uint8_t* gI1 = &(gI_P[j1]);
    uint8_t* gI2 = &(gI_P[j2]);
    for(Int d2=0;d2<4;d2++)
    for(Int c2=0;c2<3;c2++)
    for(Int d1=0;d1<4;d1++)
    for(Int c1=0;c1<3;c1++)
    {
    #pragma omp parallel for
    for(Int ji=0;ji<nmass*NTt;ji++)
    {
      const Int massi = ji/NTt;
      const Int ti    = ji%NTt;
      const Long offv = iv*nmass + massi;

      const Int off1 = massi*12*12 + (d2*3+c2)*12+gI1[d1]*3+c1;
      const Int off2 = massi*12*12 + (gI2[d2]*3+c2)*12+d1*3+c1;

      Ty g_tem = gC2[d2]*gC1[d1];

      Ty* tp1 = p1[off1*NTt+ti];
      Ty* tp2 = p2[off2*NTt+ti];
      Ty* tr0 = &(resP[offv][ti*Nxyz]);

      #if USEQACC==1
      if(invmode == 1){qacc_forNB(i, Long(Nxyz),{ tr0[i] += (tp1[i]*qlat::qconj(tp2[i]) * g_tem);});}
      if(invmode == 0){qacc_forNB(i, Long(Nxyz),{ tr0[i] += (tp1[i]*           (tp2[i]) * g_tem);});}
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
//        qlat::fft_desc_basic &fd, Int clear=1, Int invmode=1){
//  check_prop_size(prop1);check_prop_size(prop2);
//  Int  nmass = prop1.size()/(12*12*fd.Nv[3]);
//  if(nmass == 0){res.resize(0);return;}
//  Int Ngv = ga1V.size();
//  Long resL = Ngv * nmass * fd.Nv[0]*fd.Nv[1]*fd.Nv[2] * fd.Nv[3];
//  if(clear == 1){if(res.size()!= resL){res.resize(resL);}}
//
//  if(res.size() != resL){qmessage("Size of res wrong. \n");Qassert(false);}
//  Qassert(prop1.size() == prop2.size());
//
//  qlat::vector<Ta* > prop1P = EigenM_to_pointers(prop1);
//  qlat::vector<Ta* > prop2P = EigenM_to_pointers(prop2);
//
//  Ta** p1 = prop1P.data();
//  Ta** p2 = prop2P.data();
//  Ta* resP = res.data();
//  meson_vectorEV(p1, p2, resP, nmass, ga1V, ga2V, fd, clear, invmode);
//}

template <typename Ty >
void meson_vectorEV(EigenTy& prop1, EigenTy& prop2, qlat::vector_gpu<Ty > &res
  ,std::vector<ga_M > &ga1V, std::vector<ga_M > &ga2V,
  qlat::fft_desc_basic &fd, Int clear=1, Int invmode=1)
{
  check_prop_size(prop1, fd);check_prop_size(prop2, fd);
  Int  nmass = prop1.size();
  if(nmass == 0){res.resize(0);return;}
  Int Ngv = ga1V.size();
  const Long resL = Ngv * nmass * fd.Nv[0]*fd.Nv[1]*fd.Nv[2] * fd.Nv[3];
  const Long Nxyz= fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  if(clear == 1){if(Long(res.size())!= resL){res.resize(resL);}}

  if(Long(res.size()) != resL){qmessage("Size of res wrong. \n");Qassert(false);}
  Qassert(prop1.size() == prop2.size());
  for(Int mi=0;mi<nmass;mi++)
  {
    Qassert(prop1[mi].size() == 12 * 12 * fd.Nvol);
    Qassert(prop2[mi].size() == 12 * 12 * fd.Nvol);
  }

  qlat::vector<Ty* > prop1P = EigenM_to_pointers(prop1, Nxyz);
  qlat::vector<Ty* > prop2P = EigenM_to_pointers(prop2, Nxyz);

  Ty** p1 = prop1P.data();
  Ty** p2 = prop2P.data();
  qlat::vector<Ty* > resvP;resvP.resize(Ngv * nmass);
  for(Int iv = 0;iv<Ngv*nmass;iv++){
    resvP[iv] = &res[iv* fd.Nv[3] * Nxyz];
  }
  //Ty* resP = res.data();
  meson_vectorEV(p1, p2, resvP.data(), nmass, ga1V, ga2V, fd, clear, invmode);
}

template <typename Ty >
void meson_vectorEV(std::vector<qlat::FieldG<Ty > >& prop1, std::vector<qlat::FieldG<Ty > >& prop2,
  std::vector<qlat::FieldG<Ty > > &res, std::vector<ga_M > &ga1V, std::vector<ga_M > &ga2V, Int clear=1, Int invmode=1)
{
  const Int nmass = prop1.size();
  if(nmass == 0){res.resize(0);return;}
  Qassert(prop1[0].initialized);
  const Geometry& geo = prop1[0].geo();
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);

  const Int Ngv = ga1V.size();
  const Long Nres = Ngv * nmass;
  const Long Nxyz= fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  if(clear == 1){
    if(Long(res.size()) != Nres){
      res.resize(0);res.resize(Nres);
    }
    for(Long vi=0;vi<Nres;vi++){
      //  print_mem_info(ssprintf("vi %5d ", int(vi)));
      if(!res[vi].initialized){
        res[vi].init(geo, 1, QMGPU, QLAT_OUTTER);
      }
      set_zero(res[vi]);
    }
  }
  if(Long(res.size()) != Nres or !res[0].initialized){
    qmessage("Size of res wrong. \n");Qassert(false);
  }
  Qassert(prop1.size() == prop2.size());

  for(Int mi=0;mi<nmass;mi++)
  {
    Qassert(prop1[mi].multiplicity == 144 and prop2[mi].multiplicity == 144);
    Qassert(prop1[mi].mem_order == QLAT_OUTTER and prop2[mi].mem_order == QLAT_OUTTER);
  }

  qlat::vector<Ty* > prop1P = FieldG_to_pointers(prop1, Nxyz);
  qlat::vector<Ty* > prop2P = FieldG_to_pointers(prop2, Nxyz);
  qlat::vector<Ty* > resP   = FieldG_to_pointers(res  , fd.Nv[3] * Nxyz);

  Ty** p1 = prop1P.data();
  Ty** p2 = prop2P.data();
  //  Ty* resP = res.data();
  meson_vectorEV(p1, p2, resP.data(), nmass, ga1V, ga2V, fd, clear, invmode);
}

template<typename Td, typename Ta>
void meson_corrE(std::vector<Propagator4dT<Td > > &pV1, std::vector<Propagator4dT<Td >> &pV2,  ga_M &ga1, ga_M &ga2,
  qlat::vector<Ta > &res, qlat::fft_desc_basic &fd,Int clear=1,const Coordinate& mom = Coordinate(), const Int invmode=1, const Int tini = 0){
  ///int NTt  = fd.Nv[3];
  ///LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  Int nmass = pV1.size();
  ///int nt = fd.nt;

  qlat::vector<Ta > resE;
  ini_resE(resE,nmass,fd);

  meson_vectorE(pV1,pV2,ga1,ga2,resE,fd, 1, invmode);

  vec_corrE(resE,res,fd, clear, mom, qlat::ComplexT<Td>(1.0, 0.0), tini);
}

template<typename Td, typename Ta>
void meson_corrE(Propagator4dT<Td > &p1, Propagator4dT<Td > &p2,  ga_M &ga1, ga_M &ga2,
  qlat::vector<Ta > &res, Int clear=1,const Coordinate& mom = Coordinate(), Int invmode=1, const Int tini = 0){
  const qlat::Geometry& geo = p1.geo();

  std::vector<Propagator4dT<Td > > P1;
  std::vector<Propagator4dT<Td > > P2;
  P1.resize(1);P2.resize(1);

  P1[0].init(geo);P2[0].init(geo);

  P1[0] = p1;P2[0] = p2;
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);
  meson_corrE(P1, P2, ga1, ga2, res, fd, clear, mom, invmode, tini);
}

template<typename Td>
void meson_corrE(Propagator4dT<Td > &p1, Propagator4dT<Td > &p2, const Int ga, const Int gb,
  const std::string& filename, const Coordinate& mom = Coordinate(), Int invmode=1, const Int tini = 0,
  const std::string& info = std::string("NONE"), const Int shift_end = 1)
{
  ga_matrices_cps ga_cps;
  std::vector<ga_M > gL;gL.resize(16);
  {int o=0;
  for(Int i=0;i<6;i++){gL[o] = ga_cps.ga[0][i];o+=1;}
  for(Int i=2;i<6;i++){gL[o] = ga_cps.ga[1][i];o+=1;}
  for(Int i=3;i<6;i++){gL[o] = ga_cps.ga[2][i];o+=1;}
  for(Int i=4;i<6;i++){gL[o] = ga_cps.ga[3][i];o+=1;}
  for(Int i=5;i<6;i++){gL[o] = ga_cps.ga[4][i];o+=1;}}

  qlat::vector<qlat::ComplexT<RealD > > Res;
  meson_corrE(p1, p2, gL[ga], gL[gb], Res, 1, mom, invmode, tini);
  const size_t sizen = get_file_size_MPI(filename, true);
  corr_dat<RealD > corr(std::string(""));
  Int nt = Res.size();
  if(sizen > 0){
    corr.read_dat(filename, 1);
    if(shift_end == 1){
      corr.shift_end();
    }
  }
  else{
    std::string ktem = ssprintf("%d  %d  %d", 1, nt, 2);
    std::string dtem = ssprintf("nsrc nt complex");
    corr.create_dat(ktem, dtem);
  }
  if(info != std::string("NONE") and info.size() != 0){
    corr.INFOA.push_back(info);
  }

  corr.write_corr(Res.data(), Res.size());
  corr.write_dat(filename);
}

template<typename Ty>
void meson_corrE(std::vector<qpropT > &prop1, std::vector<qpropT > &prop2,  ga_M &ga1, ga_M &ga2,
  qlat::vector<Ty >& res, Int clear=1,const Coordinate& mom = Coordinate(), Int invmode=1){
  qlat::vector<Ty > resE;

  const qlat::Geometry& geo = prop1[0].geo();
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);

  meson_vectorE(prop1,prop2,ga1,ga2,resE,1, invmode);
  vec_corrE(resE, res, fd, clear, mom);
}

template<typename Td>
void print_pion(std::vector<Propagator4dT<Td > > &prop1, std::vector<Propagator4dT<Td > > &prop2, const std::string& tag=std::string(""), RealD factor = 1.0){
  ga_matrices_cps   ga_cps;
  const qlat::Geometry& geo = prop1[0].geo();
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);

  qlat::vector<qlat::ComplexT<Td > > resC;

  meson_corrE(prop1, prop2, ga_cps.ga[0][0], ga_cps.ga[0][0], resC, fd);

  Int nv = resC.size()/fd.nt;
  for(Int iv=0;iv<nv;iv++)
  for(Int t=0;t<fd.nt;t++)
  {
    qlat::ComplexT<Td > v = resC[iv*fd.nt + t] * qlat::ComplexT<Td >(factor, 0.0);
    qmessage("%s iv %2d, t %3d, v %.15e %.15e \n", tag.c_str(), iv, t, v.real(), v.imag());
  }
}

template <typename Ty >
void pion_corr_simple(qlat::FieldG<Ty >& prop1, qlat::FieldG<Ty >& prop2, vector<Ty>& corr, 
  const Coordinate& mom = Coordinate(), const Ty src_phase = Ty(1.0, 0.0), const Int tini = 0)
{
  ga_matrices_cps   ga_cps;
  std::vector<ga_M> gI;gI.resize(1);gI[0] = ga_cps.ga[0][0];
  std::vector<qlat::FieldG<Ty > > p1;
  std::vector<qlat::FieldG<Ty > > p2;
  std::vector<qlat::FieldG<Ty > > res;
  p1.resize(1);p2.resize(1);
  p1[0].set_pointer(prop1);
  p2[0].set_pointer(prop2);
  meson_vectorEV(p1, p2, res, gI, gI, 1, 1);

  fft_desc_basic& fd = get_fft_desc_basic_plan(prop1.geo());
  const Int clear = 1;
  vec_corrE(res[0].field_gpu, corr, fd, clear, mom, src_phase, tini);
}

template<typename Td>
void print_pion(Propagator4dT<Td > &p1, Propagator4dT<Td >&p2, const std::string& tag=std::string(""), RealD factor = 1.0){
  std::vector<Propagator4dT<Td > > prop1;
  std::vector<Propagator4dT<Td > > prop2;
  const qlat::Geometry& geo = p1.geo();
  prop1.resize(1);prop2.resize(1);
  prop1[0].init(geo);prop2[0].init(geo);
  prop1[0] = p1;prop2[0] = p2;
  print_pion(prop1, prop2, tag, factor);
}

}

#endif

