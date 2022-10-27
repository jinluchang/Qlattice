// utils_corr_baryon.h
// Gen Wang
// Jul. 2021

#ifndef UTILS_CORR_BARYON_H
#define UTILS_CORR_BARYON_H

#pragma once

#include "utils_float_type.h"
#include "utils_gammas.h"
#include "utils_fft_desc.h"
#include "utils_reduce_vec.h"
#include "utils_grid_src.h"
#include "utils_io_vec.h"
#include "utils_corr_prop.h"

namespace qlat{

#ifdef QLAT_USE_ACC
template<typename Ty, int bfac, int Blocks>
__global__ void baryon_vectorEV_global(Ty** p1, Ty** p2, Ty** p3, Ty* resP,
  char** gPP, unsigned char** oPP, const int* ivP,
  const int nmass, const int NTt, const long Nxyz, const int Ngv)
{
  //unsigned long gi =  blockIdx.x*blockDim.x + threadIdx.x;
  const unsigned long gi =  blockIdx.x;
  //const unsigned int tid = threadIdx.x;
  const unsigned int tid = threadIdx.y*blockDim.x+ threadIdx.x;
  /////const unsigned int Tid = blockDim.x;
  const long Ntotal = nmass * NTt * Nxyz;
  const long Nbfac  = Ntotal/bfac;
  __shared__ Ty P1[bfac*12*12];
  __shared__ Ty P2[bfac*12*12];
  __shared__ Ty P3[bfac*12*12];

  if(gi*bfac < Ntotal){

  //__shared__ Ty buf[bfac*16+1];
  ////const long offR0 = gi;
  int bi0= 0;int dc = 0;
  int ji = 0;int massi = 0;int ti = 0;
  ////const long ixyz = (bi0*Nbfac + gi)%Nxyz;

  int jobN = bfac*12*12;
  unsigned int off = tid;
  while(off < jobN){
    bi0= off/(12*12);
    dc = off%(12*12);
    ji    = (bi0*Nbfac + gi)/Nxyz;
    massi = ji/NTt;
    ti    = ji%NTt;
    long ixyz = (bi0*Nbfac + gi)%Nxyz;
    //massi = (ji + bi)/NTt;
    //ti    = (ji + bi)%NTt;
    P1[dc*bfac + bi0] = p1[(massi*12*12 + dc)*NTt + ti][ixyz];
    P2[dc*bfac + bi0] = p2[(massi*12*12 + dc)*NTt + ti][ixyz];
    P3[dc*bfac + bi0] = p3[(massi*12*12 + dc)*NTt + ti][ixyz];
    off += Blocks;
  }
  __syncthreads();

  int ini = 0;
  int dv = 0;
  unsigned int MAX = 0;

  const int bfacC = bfac;
  const int Nth   = Blocks/bfac;
  const unsigned int Each  = 16*bfacC;
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
    jobN = (MAX + GROUP -1 )/GROUP;
    ini = 0; dv = GROUP;
    for(int ji=0;ji<jobN;ji++){
      if(ini + dv >= MAX){dv = MAX - ini;}
      s0 = &(oPP[iv][ini*3]);
      s1 = &(gPP[iv][ini*2]);

      off = tid;
      while(off < dv*3){pos[off] = s0[off];off += Blocks;}
      off = tid;
      while(off < dv*2){g0[off]  = s1[off];off += Blocks;}
      __syncthreads();

      off = aa;
      while(off < dv){
        const Ty* t1 = &P1[(pos[off*3+0])*bfac + ba];
        const Ty* t2 = &P2[(pos[off*3+1])*bfac + ba];
        const Ty* t3 = &P3[(pos[off*3+2])*bfac + ba];
        Ty gtem = Ty(g0[off*2+0], g0[off*2+1]);
        if(bfacC == 1){
          buf[aa*bfac+ba] += (t1[0] * t2[0] * t3[0]*gtem);
        }else{
          Ty* b0 = &buf[aa*bfac+ba];
          for(int z=0;z<bfacC;z++){b0[z] += (t1[z] * t2[z] * t3[z]*gtem);}
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
      //if(clear == 0){resP[iv*Ntotal + bi*Nbfac + gi] += (buf[bi] + buf[bfac+bi]);}
      //if(clear == 1){resP[iv*Ntotal + bi*Nbfac + gi]  = (buf[bi] + buf[bfac+bi]);}
      resP[iv*Ntotal + bi*Nbfac + gi] += (buf[bi] + buf[bfac+bi]); 
    }
    __syncthreads();

  }

  }


}
#endif

////baryon on GPU, 
////USEGLOBAL then use global functions, else use qacc
template<typename Ty, int bfac>
void baryon_vectorEV_kernel(Ty** p1, Ty** p2, Ty** p3, Ty* resP,
  char** gPP, unsigned char** oPP, const int* ivP,
  const int nmass, const int NTt, const long Nxyz, const int Ngv)
{
  long Ntotal  = nmass*NTt*Nxyz;
  ////print0("nmass %d, NTt %d, Nxyz %d \n", int(nmass), int(NTt), int(Nxyz));
  if(Ntotal % bfac != 0){abort_r("Please correct your bfac! \n");}
  long Nbfac = Ntotal/bfac;

  #if USEGLOBAL==1
  const int nt = 16;
  const int Blocks = nt*bfac;
  dim3 dimBlock(    nt, bfac, 1);
  dim3 dimGrid(  Nbfac,  1, 1);
  baryon_vectorEV_global<Ty, bfac, Blocks><<<dimGrid, dimBlock>>>(p1, p2, p3, resP, gPP, oPP, ivP, nmass, NTt, Nxyz, Ngv);
  qacc_barrier(dummy);
  #else
  if((nmass*NTt) % bfac != 0){abort_r("Please correct your bfac! \n");}
  qacc_for(gi, Nbfac ,
  {
    Ty buf[bfac+1];
    Ty P1[bfac*12*12+1];
    Ty P2[bfac*12*12+1];
    Ty P3[bfac*12*12+1];

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
        P3[dc*bfac + bi] = p3[(massi*12*12 + dc)*NTt + ti][ixyz];
      }
    }

    for(int iv=0;iv<Ngv;iv++)
    {
      for(int bi=0;bi<bfac;bi++){buf[bi] = 0;}

      for(int off=0;off<ivP[iv];off++)
      {
        const Ty* t1 = &P1[(oPP[iv][off*3+0])*bfac];
        const Ty* t2 = &P2[(oPP[iv][off*3+1])*bfac];
        const Ty* t3 = &P3[(oPP[iv][off*3+2])*bfac];
        const Ty gtem = Ty(gPP[iv][off*2+0], gPP[iv][off*2+1]);
        for(int bi=0;bi<bfac;bi++)
        {
          buf[bi] += (t1[bi] * t2[bi] * t3[bi] * gtem);
        }
      }

      long offR = iv * Ntotal;
      Ty* r0 = &resP[offR + offR0];
      for(int bi=0;bi<bfac; bi++){
        r0[bi*Nxyz]  += buf[bi];
        //if(clear == 0){r0[bi*Nxyz] += buf[bi];}
        //if(clear == 1){r0[bi*Nxyz]  = buf[bi];}
      }
    }

  });
  #endif

}

////default gpu use kernel, cpu use c++ eigen
/////A source gamma, B sink Gamma, G projections with fermion sign, mL shape of diagram
template <typename Ty>
void baryon_vectorEV(Ty** p1, Ty** p2, Ty** p3, Ty* resP, int nmass,
  ga_M &A, ga_M &B, qlat::vector_acc<Ty > &GV, 
  qlat::vector_acc<int > &mLV, fft_desc_basic &fd, int clear=1)
{
  TIMER("Proton_vectorEV");
  int NTt  = fd.Nv[3];
  long Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  int Ngv = GV.size()/16;
  qassert(GV.size()  == 16*Ngv);
  qassert(mLV.size() == 3*Ngv);

  if(clear == 1){zero_Ty(resP, Ngv*nmass*NTt*Nxyz , 1);}

  qlat::vector_acc<Ty > epslV;epslV.resize(9);
  for(int i=0;i<3;i++){epslV[i*3+i]=0;epslV[i*3 + (i+1)%3]=1;epslV[i*3 + (i+2)%3]=-1;}
  qlat::vector_acc<Ty > gMap;
  qlat::vector_acc<int > IMap;
  gMap.resize(4*2);IMap.resize(4*2);
  for(int i=0;i<4;i++){
    /////int j = + i;
    gMap[0*4+i] = A.g[i];
    gMap[1*4+i] = B.g[i];
    IMap[0*4+i] = A.ind[i];
    IMap[1*4+i] = B.ind[i];
  }

  const Ty* epsl = epslV.data();
  const Ty* gCA = &((gMap.data())[0*4]);
  const Ty* gCB = &((gMap.data())[1*4]);
  const int* gIA = &((IMap.data())[0*4]);
  const int* gIB = &((IMap.data())[1*4]);
  const Ty* GVP = GV.data();
  const int*  mLP     = mLV.data();

  /////contraction Kernel
  #if USEKERNEL==1
  ////long Ntotal  = nmass*NTt*Nxyz;
  /////const int Loff = 3*3*3*3*4*4*4*4;
  std::vector<std::vector<char > > giEL;giEL.resize(Ngv);//giEL.resize(  Ngv*Loff);
  std::vector<std::vector<unsigned char   > > oiL ;oiL.resize(Ngv );//oiL.resize(3*Ngv*Loff);
  int bmL[3];
  int nmL[3];
  for(int iv=0;iv<Ngv;iv++)
  {
    oiL[iv].resize(0);
    giEL[iv].resize(0);

    const Ty* G  = &GVP[iv*16];
    const int*      mL = &mLP[iv*3];

    for(int a1=0;a1<3;a1++)
    for(int ia=1;ia<3;ia++)
    for(int b1=0;b1<3;b1++)
    for(int ib=1;ib<3;ib++)
    {
      int b2=(b1+ib)%3,b3=(b2+ib)%3;
      int a2=(a1+ia)%3,a3=(a2+ia)%3;
      for(int m2=0;m2<4;m2++)
      for(int m1=0;m1<4;m1++)
      for(int n2=0;n2<4;n2++)
      for(int n1=0;n1<4;n1++)
      {
        const Ty Gtem =  G[m1*4+n1];
        const double norm = qlat::qnorm(Gtem);
        if(norm < 1e-20)continue;

        const int m3 = gIA[m2];
        const int n3 = gIB[n2];
        const Ty giE = epsl[a1*3 + a2]*epsl[b1*3 + b2]*gCA[m2]*gCB[n2]*G[m1*4+n1];
        nmL[0] = n1;nmL[1] = n2;nmL[2] = n3;
        bmL[0] = b1;bmL[1] = b2;bmL[2] = b3;
        const int nm1 = nmL[mL[0]];
        const int nm2 = nmL[mL[1]];
        const int nm3 = nmL[mL[2]];

        const int bm1 = bmL[mL[0]];
        const int bm2 = bmL[mL[1]];
        const int bm3 = bmL[mL[2]];

        const int o1 = (m1*3+a1)*12+(nm1*3+bm1);
        const int o2 = (m2*3+a2)*12+(nm2*3+bm2);
        const int o3 = (m3*3+a3)*12+(nm3*3+bm3);

        ////buf += (P1[o1] * P2[o2] *P3[o3] * giE);
        oiL[iv].push_back(o1);
        oiL[iv].push_back(o2);
        oiL[iv].push_back(o3);
        giEL[iv].push_back(char(giE.real()));
        giEL[iv].push_back(char(giE.imag()));
        //giEL[iv].push_back(giE);

      }
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
  ///int maxNv = iv_size[0];
  ///for(int iv=0;iv<Ngv;iv++){if(iv_size[iv] > maxNv){maxNv = iv_size[iv];}}

  //{
  //long total = 0;
  //for(int iv=0;iv<Ngv;iv++){total += iv_size[iv];}
  //print0("==Ngv %d, total %d \n", int(Ngv), int(total));
  //}

  ////int mode = 0;mode = clear;
  const int BFACG = BFACG_SHARED;
  baryon_vectorEV_kernel<Ty, BFACG>(p1, p2, p3, resP, gPP, oPP, ivP, nmass, NTt, Nxyz, Ngv);

  #endif

  #if USEKERNEL==0
  for(int iv=0;iv<Ngv;iv++)
  {
    long offR = iv*nmass*NTt * Nxyz;
    const Ty* G  = &(GVP[iv*16 + 0]);
    const int*      mL = &(mLP[iv*3 + 0]);
    int bmL[3];
    int nmL[3];

    for(int a1=0;a1<3;a1++)
    for(int ia=1;ia<3;ia++)
    for(int b1=0;b1<3;b1++)
    for(int ib=1;ib<3;ib++)
    {
      int b2=(b1+ib)%3,b3=(b2+ib)%3;
      int a2=(a1+ia)%3,a3=(a2+ia)%3;
      for(int m2=0;m2<4;m2++)
      for(int m1=0;m1<4;m1++)
      for(int n2=0;n2<4;n2++)
      for(int n1=0;n1<4;n1++)
      {
        Ty Gtem =  G[m1*4+n1];
        double norm = qlat::qnorm(Gtem);
        if(norm < 1e-20)continue;

        int m3 = gIA[m2];
        int n3 = gIB[n2];
        Ty giE = epsl[a1*3 + a2]*epsl[b1*3 + b2]*gCA[m2]*gCB[n2]*G[m1*4+n1];
        nmL[0] = n1;nmL[1] = n2;nmL[2] = n3;
        bmL[0] = b1;bmL[1] = b2;bmL[2] = b3;
        int nm1 = nmL[mL[0]];
        int nm2 = nmL[mL[1]];
        int nm3 = nmL[mL[2]];

        int bm1 = bmL[mL[0]];
        int bm2 = bmL[mL[1]];
        int bm3 = bmL[mL[2]];

        #pragma omp parallel for
        for(int ji=0;ji<nmass*NTt;ji++)
        {
          int massi = ji/NTt;
          int ti    = ji%NTt;

          int o1 = massi*12*12 + (m1*3+a1)*12+(nm1*3+bm1);
          int o2 = massi*12*12 + (m2*3+a2)*12+(nm2*3+bm2);
          int o3 = massi*12*12 + (m3*3+a3)*12+(nm3*3+bm3);

          Ty* tp1 = p1[o1*NTt+ti];
          Ty* tp2 = p2[o2*NTt+ti];
          Ty* tp3 = p3[o3*NTt+ti];
          Ty* tr0 = &(resP[offR + (massi*NTt + ti)*Nxyz]);

          #if USEQACC==1
          qacc_forNB(i, long(Nxyz),{ tr0[i] += (tp1[i]*tp2[i]*tp3[i] * giE); });
          #else
          EAy vp1(tp1,Nxyz);
          EAy vp2(tp2,Nxyz);
          EAy vp3(tp3,Nxyz);
          EAy vr0(tr0,Nxyz);
          vr0 += (vp1*vp2*vp3 * giE);
          #endif

        }
        qacc_barrier(dummy);
      }
    }
  }
  #endif

}

////container 
template <typename Ta>
void baryon_vectorEV(EigenMTa &prop1, EigenMTa &prop2, EigenMTa &prop3, EigenVTa &res,
  ga_M &A, ga_M &B, qlat::vector_acc<Ta > &GV, qlat::vector_acc<int > &mLV,
  fft_desc_basic &fd,int clear=1){
  int NTt  = fd.Nv[3];
  long Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];

  check_prop_size(prop1);check_prop_size(prop2);check_prop_size(prop3);
  int nmass = prop1.size()/(12*12*NTt);
  qassert(prop1.size() == prop2.size());
  qassert(prop1.size() == prop3.size());
  int Ngv = GV.size()/16;
  long resL = Ngv * nmass*NTt * Nxyz;
  if(clear == 1){if(res.size()!= resL){res.resize(resL); } }
  if(res.size() != resL){print0("Size of res wrong. \n");qassert(false);}

  qlat::vector_acc<Ta* > prop1P = EigenM_to_pointers(prop1);
  qlat::vector_acc<Ta* > prop2P = EigenM_to_pointers(prop2);
  qlat::vector_acc<Ta* > prop3P = EigenM_to_pointers(prop3);
  Ta** p1 = prop1P.data();
  Ta** p2 = prop2P.data();
  Ta** p3 = prop3P.data();
  Ta* resP = res.data();

  baryon_vectorEV(p1, p2, p3, resP, nmass, A, B, GV, mLV, fd, clear);
}

////simple gpu version
///////Proton contractions
template <typename Ta>
void proton_vectorE(EigenMTa &prop1, EigenMTa &prop2, EigenMTa &prop3,
  const ga_M &ga2,const int ind2, const ga_M &ga1, const int ind1,
      EigenVTa &res, fft_desc_basic &fd,int clear=1){
  TIMER("Proton_vectorE");
  ////ga2/ind2 for source, gam1/ind1 for sink
  ////"[]|+N" type diagram
  check_prop_size(prop1);check_prop_size(prop2);check_prop_size(prop3);
  int NTt  = fd.Nv[3];
  LInt Nxyz = prop1[0].size();
  int nmass = prop1.size()/(12*12*NTt);
  qassert(prop1.size() == prop2.size());
  qassert(prop1.size() == prop3.size());
  if(clear == 1){ini_resE(res,nmass,fd);}
  if(clear == 0){qassert(res.size() == long(nmass*NTt * Nxyz));}
    
  ////Prop format, src d-4, c-3, sink d-4, c-3, Nt, EigenVTa<Nxyz>
  if(res.size()%NTt !=0 or res.size()==0){print0("Size of res wrong. \n");qassert(false);}

  Ta epsl[3][3];
  for(int i=0;i<3;i++)for(int j=0;j<3;j++){epsl[i][j] = 0;}
  for(int i=0;i<3;i++){epsl[i][i]=0;epsl[i][(i+1)%3]=Ta(1,0);epsl[i][(i+2)%3]=Ta(-1,0);}

  for(int d2=0;d2<4;d2++)
  for(int c21=0;c21<3;c21++)
  for(int ib=1;ib<3;ib++)
  {
    int c22=(c21+ib)%3,c23=(c22+ib)%3;
    for(int d1=0;d1<4;d1++)
    for(int c11=0;c11<3;c11++)
    for(int ia=1;ia<3;ia++)
    {
      int c12=(c11+ia)%3,c13=(c12+ia)%3;
      Ta giE = epsl[c11][c12]*epsl[c21][c22]*ga1.g[d1]*ga2.g[d2];

      #pragma omp parallel for
      for(int ji=0;ji<nmass*NTt;ji++)
      {
        int massi = ji/NTt;
        int ti    = ji%NTt;

        int m1 = massi*12*12 + (ind2*3+c21)*12+ind1*3+c11;
        int m2 = massi*12*12 + (ga2.ind[d2]*3+c22)*12+d1*3+c12;
        int m3 = massi*12*12 + (d2*3+c23)*12+ga1.ind[d1]*3+c13;

        int n1 = massi*12*12 + (ind2*3+c21)*12+ga1.ind[d1]*3+c11;
        int n2 = massi*12*12 + (ga2.ind[d2]*3+c22)*12+d1*3+c12;
        int n3 = massi*12*12 + (d2*3+c23)*12+ind1*3+c13;

        Ta* tp1 = prop1[m1*NTt+ti].data();
        Ta* tp2 = prop2[m2*NTt+ti].data();
        Ta* tp3 = prop3[m3*NTt+ti].data();
        Ta* tn1 = prop1[n1*NTt+ti].data();
        Ta* tn2 = prop2[n2*NTt+ti].data();
        Ta* tn3 = prop3[n3*NTt+ti].data();
        Ta* tr0 = &(res.data()[(massi*NTt + ti)*Nxyz]);

        #if USEQACC==1
        qacc_forNB(i, long(Nxyz),{ tr0[i] -= ((tp1[i]*tp2[i]*tp3[i] + tn1[i]*tn2[i]*tn3[i])*giE); });
        #else
        EAa vp1(tp1,Nxyz);
        EAa vp2(tp2,Nxyz);
        EAa vp3(tp3,Nxyz);
        EAa vn1(tn1,Nxyz);
        EAa vn2(tn2,Nxyz);
        EAa vn3(tn3,Nxyz);
        EAa vr0(tr0,Nxyz);
        vr0 -= ((vp1*vp2*vp3 + vn1*vn2*vn3)*giE);
        #endif

      }
      qacc_barrier(dummy);
    }
  }

}

///container for 
template <typename Ta>
void proton_vectorE(EigenMTa &prop1, EigenMTa &prop2, EigenMTa &prop3,
        EigenVTa &res, fft_desc_basic &fd, const ga_M &ga1,int t0,int dT,int clear=1,int oppo=0){
  TIMER("Proton_vectorE");
  int NTt  = fd.Nv[3];
  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  int nmass = prop1.size()/(12*12*NTt);
  qassert(prop1.size() == prop2.size());
  qassert(prop1.size() == prop3.size());

  if(clear == 1){ini_resE(res,nmass,fd);}

  //int nv = res.size();int Nsize = res[0].size();
  EigenVTa resE0;resE0.resize(res.size());
  EigenVTa resE1;resE1.resize(res.size());
  ////qlat::set_zero(resE0);qlat::set_zero(resE1);

  proton_vectorE(prop1,prop2,prop3,ga1,0,ga1,0,resE0,fd,1);
  proton_vectorE(prop1,prop2,prop3,ga1,1,ga1,1,resE0,fd,0);
  proton_vectorE(prop1,prop2,prop3,ga1,2,ga1,2,resE1,fd,1);
  proton_vectorE(prop1,prop2,prop3,ga1,3,ga1,3,resE1,fd,0);

  std::vector<int > map_sec = get_sec_map(dT,fd.nt);
  //////int Nt = fd.Nt;

  /////int t0 = 0;
  int nt = fd.nt;
  ///for(int massi=0;massi<nmass;massi++)
  ///for(int ti = 0;ti<Nt;ti++)
  #pragma omp parallel for
  for(int ji=0;ji<nmass*NTt;ji++)
  {
    int massi = ji/NTt;
    int ti    = ji%NTt;
    int t = ti + fd.Pos0[fd.rank][3];
    Ta* tr0 = &(res.data()[(massi*NTt+ti)*Nxyz]);
    Ta* tv0 = &(resE0.data()[(massi*NTt+ti)*Nxyz]);
    Ta* tv1 = &(resE1.data()[(massi*NTt+ti)*Nxyz]);

    #if USEQACC==0
    EAa r0(tr0,Nxyz);
    EAa v0(tv0,Nxyz);
    EAa v1(tv1,Nxyz);
    #endif

    if(map_sec[(t-t0+nt)%nt]%2==0)
    {
      #if USEQACC==1
      if(oppo==0)qacc_forNB(i, long(Nxyz), { tr0[i] += tv0[i]; });
      if(oppo==1)qacc_forNB(i, long(Nxyz), { tr0[i] += tv1[i]; });
      #else
      if(oppo==0){r0 += v0;}
      if(oppo==1){r0 += v1;}
      #endif

    }
    if(map_sec[(t-t0+nt)%nt]%2==1)
    {
      #if USEQACC==1
      if(oppo==0)qacc_forNB(i, long(Nxyz), { tr0[i] += tv1[i]; });
      if(oppo==1)qacc_forNB(i, long(Nxyz), { tr0[i] += tv0[i]; });
      #else
      if(oppo==0){r0 += v1;}
      if(oppo==1){r0 += v0;}
      #endif

    }
  }
  qacc_barrier(dummy);
}

////container
template <typename Ta>
void proton_corrE(EigenMTa &prop1, EigenMTa &prop2, EigenMTa &prop3,
  const ga_M &ga2,const int ind2, const ga_M &ga1,const int ind1,
  EigenVTa &res, fft_desc_basic &fd,int clear=1,const Coordinate& mom = Coordinate()){
  EigenVTa resE;
  proton_vectorE(prop1,prop2,prop3,ga2,ind2,ga1,ind1,resE,fd,1);
  vec_corrE(resE,res,fd,clear,mom);
}

////container
template <typename Ta>
void proton_corrE(EigenMTa &prop1, EigenMTa &prop2, EigenMTa &prop3,
 EigenVTa &res, fft_desc_basic &fd, const ga_M &ga1,const int t0,const int dT,int clear=1,const Coordinate& mom = Coordinate()){
  EigenVTa resE;
  proton_vectorE(prop1,prop2,prop3,resE,fd, ga1, t0,dT,1);

  vec_corrE(resE,res,fd,clear,mom);
}


////simple gpu version
/////A source gamma, B sink Gamma, G projections with fermion sign, mL shape of diagram
template <typename Ta>
void baryon_vectorE(EigenMTa &prop1, EigenMTa &prop2, EigenMTa &prop3,
  ga_M &A, ga_M &B, qlat::vector_acc<Ta > &G, qlat::vector_acc<int > &mL,
        EigenVTa &res, fft_desc_basic &fd,int clear=1){
  TIMER("Proton_vectorE");
  check_prop_size(prop1);check_prop_size(prop2);check_prop_size(prop3);
  int NTt  = fd.Nv[3];
  LInt Nxyz = prop1[0].size();
  int nmass = prop1.size()/(12*12*NTt);
  qassert(prop1.size() == prop2.size());
  qassert(prop1.size() == prop3.size());
  qassert(G.size()  == 16);
  qassert(mL.size() == 3);
  if(clear == 1){ini_resE(res,nmass,fd);}

  //if(res.size()%NTt !=0 or res.size()==0){print0("Size of res wrong. \n");qassert(false);}
  if(res.size()==0){print0("Size of res wrong. \n");qassert(false);}

  Ta epsl[3][3];
  for(int i=0;i<3;i++)for(int j=0;j<3;j++){epsl[i][j] = 0;}
  for(int i=0;i<3;i++){epsl[i][i]=0;epsl[i][(i+1)%3]=1;epsl[i][(i+2)%3]=-1;}

  //mL = {};
  //std::vector<int > mL;mL.resize(3);
  //mL[0] = 0;mL[1] = 1;mL[2] = 2;
  std::vector<int > nmL;nmL.resize(3);
  std::vector<int > bmL;bmL.resize(3);

  {
    for(int a1=0;a1<3;a1++)
    for(int ia=1;ia<3;ia++)
    for(int b1=0;b1<3;b1++)
    for(int ib=1;ib<3;ib++)
    {
      int b2=(b1+ib)%3,b3=(b2+ib)%3;
      int a2=(a1+ia)%3,a3=(a2+ia)%3;
      for(int m2=0;m2<4;m2++)
      for(int m1=0;m1<4;m1++)
      for(int n2=0;n2<4;n2++)
      for(int n1=0;n1<4;n1++)
      {
        Ta Gv =  G[m1*4+n1];
        double norm = std::sqrt(Gv.real()*Gv.real() + Gv.imag()*Gv.imag());
        if(norm < 1e-20)continue;

        int m3 = A.ind[m2];
        int n3 = B.ind[n2];
        Ta giE = epsl[a1][a2]*epsl[b1][b2]*A.g[m2]*B.g[n2]*G[m1*4+n1];
        nmL[0] = n1;nmL[1] = n2;nmL[2] = n3;
        bmL[0] = b1;bmL[1] = b2;bmL[2] = b3;
        int nm1 = nmL[mL[0]];
        int nm2 = nmL[mL[1]];
        int nm3 = nmL[mL[2]];

        int bm1 = bmL[mL[0]];
        int bm2 = bmL[mL[1]];
        int bm3 = bmL[mL[2]];

        #pragma omp parallel for
        for(int ji=0;ji<nmass*NTt;ji++)
        {
          int massi = ji/NTt;
          int ti    = ji%NTt;

          int o1 = massi*12*12 + (m1*3+a1)*12+(nm1*3+bm1);
          int o2 = massi*12*12 + (m2*3+a2)*12+(nm2*3+bm2);
          int o3 = massi*12*12 + (m3*3+a3)*12+(nm3*3+bm3);

          Ta* tp1 = prop1[o1*NTt+ti].data();
          Ta* tp2 = prop2[o2*NTt+ti].data();
          Ta* tp3 = prop3[o3*NTt+ti].data();
          Ta* tr0 = &(res.data()[(massi*NTt + ti)*Nxyz]);

          #if USEQACC==1
          qacc_forNB(i, long(Nxyz),{ tr0[i] += (tp1[i]*tp2[i]*tp3[i] * giE); });
          #else
          EAa vp1(tp1,Nxyz);
          EAa vp2(tp2,Nxyz);
          EAa vp3(tp3,Nxyz);
          EAa vr0(tr0,Nxyz);
          vr0 += (vp1*vp2*vp3 * giE);
          #endif

        }
        qacc_barrier(dummy);
      }
    }
  }

}

////3pt insertion version
////A source gamma, B sink Gamma, G projections with fermion sign, mL shape of diagram
template <typename Ta>
void baryon_vectorE(EigenMTa &prop1, EigenMTa &prop2, EigenMTa &prop3,
  ga_M &A, ga_M &B, qlat::vector_acc<Ta > &G, qlat::vector_acc<int > &mL, int insertion,
        EigenMTa &resP, fft_desc_basic &fd,int clear=1)
{
  TIMER("Baryon_vectorE_insertion");
  check_prop_size(prop1);check_prop_size(prop2);check_prop_size(prop3);
  int NTt  = fd.Nv[3];
  LInt Nxyz = prop1[0].size();
  int nmass = prop1.size()/(12*12*NTt);
  qassert(prop1.size() == prop2.size());
  qassert(prop1.size() == prop3.size());
  qassert(G.size()  == 16);
  qassert(mL.size() == 3);
  qassert(fd.order_ch == 0);

  if(clear==1){ini_propE(resP,nmass,fd);}
  /////check_prop_size(resP);

  Ta epsl[3][3];
  for(int i=0;i<3;i++)for(int j=0;j<3;j++){epsl[i][j] = 0;}
  for(int i=0;i<3;i++){epsl[i][i]=0;epsl[i][(i+1)%3]=1;epsl[i][(i+2)%3]=-1;}

  std::vector<int > nmL;nmL.resize(3);
  std::vector<int > bmL;bmL.resize(3);

  {
    for(int a1=0;a1<3;a1++)
    for(int ia=1;ia<3;ia++)
    for(int b1=0;b1<3;b1++)
    for(int ib=1;ib<3;ib++)
    {
      int b2=(b1+ib)%3,b3=(b2+ib)%3;
      int a2=(a1+ia)%3,a3=(a2+ia)%3;
      for(int m2=0;m2<4;m2++)
      for(int m1=0;m1<4;m1++)
      for(int n2=0;n2<4;n2++)
      for(int n1=0;n1<4;n1++)
      {
        Ta Gv =  G[m1*4+n1];
        double norm = std::sqrt(Gv.real()*Gv.real() + Gv.imag()*Gv.imag());
        if(norm < 1e-20)continue;

        int m3 = A.ind[m2];
        int n3 = B.ind[n2];
        Ta giE = epsl[a1][a2]*epsl[b1][b2]*A.g[m2]*B.g[n2]*G[m1*4+n1];
        nmL[0] = n1;nmL[1] = n2;nmL[2] = n3;
        bmL[0] = b1;bmL[1] = b2;bmL[2] = b3;
        int nm1 = nmL[mL[0]];
        int nm2 = nmL[mL[1]];
        int nm3 = nmL[mL[2]];

        int bm1 = bmL[mL[0]];
        int bm2 = bmL[mL[1]];
        int bm3 = bmL[mL[2]];

        #pragma omp parallel for
        for(int ji=0;ji<nmass*NTt;ji++)
        {
          int massi = ji/NTt;
          int ti    = ji%NTt;

          int o1 = massi*12*12 + (m1*3+a1)*12+(nm1*3+bm1);
          int o2 = massi*12*12 + (m2*3+a2)*12+(nm2*3+bm2);
          int o3 = massi*12*12 + (m3*3+a3)*12+(nm3*3+bm3);

          int r0 = 0;
          if(insertion == 0){r0 = massi*12*12 + (m1*3+a1)*12+(nm1*3+bm1);}
          if(insertion == 1){r0 = massi*12*12 + (m2*3+a2)*12+(nm2*3+bm2);}
          if(insertion == 2){r0 = massi*12*12 + (m3*3+a3)*12+(nm3*3+bm3);}

          Ta* tp1 = prop1[o1*NTt+ti].data();
          Ta* tp2 = prop2[o2*NTt+ti].data();
          Ta* tp3 = prop3[o3*NTt+ti].data();
          Ta* tr0 = resP[r0*NTt+ti].data();

          #if USEQACC==1
          if(insertion == 0)qacc_forNB(i, long(Nxyz),{ tr0[i] += (tp2[i]*tp3[i] * giE); });
          if(insertion == 1)qacc_forNB(i, long(Nxyz),{ tr0[i] += (tp1[i]*tp3[i] * giE); });
          if(insertion == 2)qacc_forNB(i, long(Nxyz),{ tr0[i] += (tp1[i]*tp2[i] * giE); });
          #else
          EAa vp1(tp1,Nxyz);
          EAa vp2(tp2,Nxyz);
          EAa vp3(tp3,Nxyz);
          EAa vr0(tr0,Nxyz);
          if(insertion == 0){vr0 += vp2*vp3 * giE;}
          if(insertion == 1){vr0 += vp1*vp3 * giE;}
          if(insertion == 2){vr0 += vp1*vp2 * giE;}
          #endif
        }
        qacc_barrier(dummy);
      }
    }
  }

}

template <typename Ta>
void baryon_corrE(EigenMTa &prop1, EigenMTa &prop2, EigenMTa &prop3,
   ga_M &ga2,int ind2,ga_M &ga1,int ind1,
  EigenVTa &res, fft_desc_basic &fd,int clear=1,const Coordinate& mom = Coordinate())
{
  int NTt  = fd.Nv[3];
  ////LInt Nxyz = prop1[0].size();
  int nmass = prop1.size()/(12*12*NTt);
  ////int nt = fd.nt;

  EigenVTa resE;
  ini_resE(resE,nmass,fd);

  qlat::vector_acc<Ta > G;G.resize(16);
  qlat::vector_acc<int > mL;mL.resize(3);

  clear_qv(G);G[ind2*4 + ind1] = +1.0;
  mL[0] = 0;mL[1] = 1;mL[2] = 2;
  baryon_vectorE(prop1,prop2,prop3, ga2,ga1, G, mL, resE, fd, 1);
  clear_qv(G);G[ind2*4 + ind1] = -1.0;
  mL[0] = 1;mL[1] = 0;mL[2] = 2;
  baryon_vectorE(prop1,prop2,prop3, ga2,ga1, G, mL, resE, fd, 0);

  ////proton_vectorE(prop1,prop2,prop3,ga2,ind2,ga1,ind1,resE,fd,1);

  vec_corrE(resE,res,fd,clear,mom);
}

template <typename Ta>
void Omega_corrE(EigenMTa &prop1, EigenMTa &prop2, EigenMTa &prop3,
   ga_M &ga2,int ind2,ga_M &ga1,int ind1,
  EigenVTa &res, fft_desc_basic &fd,int clear=1,const Coordinate& mom = Coordinate())
{
  int NTt  = fd.Nv[3];
  ///LInt Nxyz = prop1[0].size();
  int nmass = prop1.size()/(12*12*NTt);
  ///int nt = fd.nt;

  EigenVTa resE;
  ini_resE(resE,nmass,fd);

  qlat::vector_acc<Ta > G;G.resize(16);
  qlat::vector_acc<int > mL;mL.resize(3);

  std::vector<int > dia;dia.resize(6);
  std::vector<int > sn ;sn.resize(6);
  dia[0] = 9012;sn[0] =  1;
  dia[1] = 9102;sn[1] = -1;
  dia[2] = 9021;sn[2] = -1;
  dia[3] = 9201;sn[3] =  1;
  dia[4] = 9210;sn[4] = -1;
  dia[5] = 9120;sn[5] =  1;

  for(int di=0;di<6;di++)
  {
    clear_qv(G);G[ind2*4 + ind1] = sn[di];
    mL[0] = (dia[di]/100)%10;mL[1] =  (dia[di]%100)/10;mL[2] = dia[di]%10;
    baryon_vectorE(prop1,prop2,prop3, ga2,ga1, G, mL, resE, fd, 0);
  }
     
  ////clear_qv(G);G[ind2*4 + ind1] = +1.0;
  ////mL[0] = 0;mL[1] = 1;mL[2] = 2;
  ////baryon_vectorE(prop1,prop2,prop3, ga2,ga1, G, mL, resE, fd, 0);

  ////clear_qv(G);G[ind2*4 + ind1] = -1.0;
  ////mL[0] = 1;mL[1] = 0;mL[2] = 2;
  ////baryon_vectorE(prop1,prop2,prop3, ga2,ga1, G, mL, resE, fd, 0);

  ////clear_qv(G);G[ind2*4 + ind1] = -1.0;
  ////mL[0] = 0;mL[1] = 2;mL[2] = 1;
  ////baryon_vectorE(prop1,prop2,prop3, ga2,ga1, G, mL, resE, fd, 0);

  ////clear_qv(G);G[ind2*4 + ind1] = +1.0;
  ////mL[0] = 2;mL[1] = 0;mL[2] = 1;
  ////baryon_vectorE(prop1,prop2,prop3, ga2,ga1, G, mL, resE, fd, 0);

  ////clear_qv(G);G[ind2*4 + ind1] = -1.0;
  ////mL[0] = 2;mL[1] = 1;mL[2] = 0;
  ////baryon_vectorE(prop1,prop2,prop3, ga2,ga1, G, mL, resE, fd, 0);

  ////clear_qv(G);G[ind2*4 + ind1] = +1.0;
  ////mL[0] = 1;mL[1] = 2;mL[2] = 0;
  ////baryon_vectorE(prop1,prop2,prop3, ga2,ga1, G, mL, resE, fd, 0);


  ////proton_vectorE(prop1,prop2,prop3,ga2,ind2,ga1,ind1,resE,fd,1);

  vec_corrE(resE,res,fd,clear,mom);
}

}

#endif

