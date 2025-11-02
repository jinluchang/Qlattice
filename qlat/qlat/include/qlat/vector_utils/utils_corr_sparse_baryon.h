// utils_corr_sparse_baryon.h
// Gen Wang
// Jun. 2025

#ifndef UTILS_CORR_SPARSE_BARYON_H
#define UTILS_CORR_SPARSE_BARYON_H

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
template<typename Ty, Int bfac, Int Blocks>
__global__ void baryon_sparse_global(Ty** p1, Ty** p2, Ty** p3, Ty** resP,
  signed char** gPP, unsigned char** oPP, const Int* ivP,
  const Int Nsrc, const Long Nvol, const Int Ngv)
{
  //unsigned long gi =  blockIdx.x*blockDim.x + threadIdx.x;
  const unsigned long gi =  blockIdx.x;
  //const unsigned int tid = threadIdx.x;
  const unsigned int tid = threadIdx.y*blockDim.x+ threadIdx.x;
  /////const unsigned int Tid = blockDim.x;
  const Long Ntotal = Nsrc * Nvol;
  const Long Nbfac  = Ntotal/bfac;

  // compute loop related
  const Int bfacC = bfac;
  const Int Nth   = Blocks/bfac;
  const unsigned int Each  = 16*bfacC;
  const unsigned int GROUP = (Blocks/bfac)*Each;

  __shared__ Ty P1[bfac*12*12];
  __shared__ Ty P2[bfac*12*12];
  __shared__ Ty P3[bfac*12*12];
  __shared__ Ty buf[bfacC*Blocks];
  __shared__ unsigned char pos[3*GROUP];
  __shared__ signed   char g0[2*GROUP];

  // gi * bfac == Ntotal
  const unsigned int jobN = bfac*12*12;
  unsigned int off = tid;
  {
  unsigned int bi, dc;
  Long ga  , isp, o1;
  Int  si  ;

  // Blocks : number of threads within each block == nt * bfac
  while(off < jobN){
    bi  = off/(12*12);
    dc  = off%(12*12);
    ga  = bi * Nbfac + gi;
    si  = ga / Nvol;
    isp = ga % Nvol;

    o1 = dc * Nvol;

    P1[dc*bfac + bi] = p1[si][o1 + isp];
    P2[dc*bfac + bi] = p2[si][o1 + isp]; 
    P3[dc*bfac + bi] = p3[si][o1 + isp]; 
    off += Blocks;
  }
  }
  __syncthreads();

  // create loops within each Ngv
  Int ini = 0;
  Int dv = 0;
  unsigned int MAX = 0;

  const Int bi   =  threadIdx.y;
  const Int ai   =  threadIdx.x;
  const Long ga  = bi * Nbfac + gi;
  const Int  si  = ga / Nvol;
  const Long isp = ga % Nvol;

  const Int ba = (threadIdx.y/bfacC)*bfacC + 0;
  const Int aa = (threadIdx.y%bfacC)*blockDim.x + threadIdx.x;

  for(Int iv=0;iv<Ngv;iv++)
  {
    for(Int bz=0;bz<bfacC;bz++){buf[bz*Blocks + tid] = 0;}
    MAX = ivP[iv];
    Int jobN = (MAX + GROUP - 1 )/GROUP;
    ini = 0; dv = GROUP;
    for(Int ji=0;ji<jobN;ji++){
      if(ini + dv >= MAX){dv = MAX - ini;}
      if(dv <= 0){continue;}
      const unsigned char* s0 = &(oPP[iv][ini*3]);
      const   signed char* s1 = &(gPP[iv][ini*2]);

      off = tid;
      while(off < dv*3){pos[off] = s0[off];off += Blocks;}
      off = tid;
      while(off < dv*2){g0[off]  = s1[off];off += Blocks;}
      __syncthreads();

      off = aa;
      while(off < dv){
        // previous productions without int() may have issues ?
        const Ty* t1 = &P1[int(pos[off*3+0])*bfac + ba];
        const Ty* t2 = &P2[int(pos[off*3+1])*bfac + ba];
        const Ty* t3 = &P3[int(pos[off*3+2])*bfac + ba];
        const Ty gtem = Ty(g0[off*2+0], g0[off*2+1]);
        Ty* b0 = &buf[aa*bfac+ba];
        for(Int z=0;z<bfacC;z++){
          b0[z] += (t1[z] * t2[z] * t3[z] * gtem);
        }
        off += Nth*bfacC;
      }
      __syncthreads();

      ini += dv;
    }

    /*
      reduce within each Ngv
        Sum outside Nth first
    */
    for(Int atem=1;atem<bfacC;atem++){
      buf[(0*Nth+ai)*bfac+bi] += buf[(atem*Nth+ai)*bfac+bi];
    }
    __syncthreads();

    if(Nth >=256){if(ai <128){buf[ai*bfac + bi] += buf[(ai+128)*bfac + bi];}__syncthreads();}
    if(Nth >=128){if(ai < 64){buf[ai*bfac + bi] += buf[(ai+ 64)*bfac + bi];}__syncthreads();}
    if(Nth >= 64){if(ai < 32){buf[ai*bfac + bi] += buf[(ai+ 32)*bfac + bi];}__syncthreads();}
    if(Nth >= 32){if(ai < 16){buf[ai*bfac + bi] += buf[(ai+ 16)*bfac + bi];}__syncthreads();}
    if(Nth >= 16){if(ai <  8){buf[ai*bfac + bi] += buf[(ai+  8)*bfac + bi];}__syncthreads();}
    if(Nth >=  8){if(ai <  4){buf[ai*bfac + bi] += buf[(ai+  4)*bfac + bi];}__syncthreads();}
    if(Nth >=  4){if(ai <  2){buf[ai*bfac + bi] += buf[(ai+  2)*bfac + bi];}__syncthreads();}

    if(ai == 0){
      resP[iv*Nsrc + si][isp] += (buf[bi] + buf[bfac+bi]); 
    }
    __syncthreads();
  }

}
#endif

#ifdef QLAT_USE_ACC
/*
  Asssumed Ngv is large enough
*/
template<typename Ty, Int bfac, Int insertion>
__global__ void baryon_sparse_global_insertion(Ty** p1, Ty** p2, Ty** p3, Ty** resP,
  signed char** gPP, unsigned char** oPP, const Int* ivP,
  const Int Nsrc, const Long Nvol, const Int Ngv)
{
  //unsigned long gi =  blockIdx.x*blockDim.x + threadIdx.x;
  const unsigned long gi =  blockIdx.x;
  //const unsigned int tid = threadIdx.x;
  const unsigned int tid = threadIdx.y*blockDim.x+ threadIdx.x;
  const unsigned int Blocks = blockDim.y * blockDim.x;
  /////const unsigned int Tid = blockDim.x;
  const Long Ntotal = Nsrc * Nvol;
  const Long Nbfac  = Ntotal/bfac;

  __shared__ Ty P1[bfac*12*12];
  __shared__ Ty P2[bfac*12*12];
  __shared__ Ty P3[bfac*12*12];

  // gi * bfac == Ntotal
  const unsigned int jobN = bfac*12*12;
  unsigned int off = tid;
  {
  unsigned int bi, dc;
  Long ga  , isp, o1;
  Int  si  ;

  // Blocks : number of threads within each block == nt * bfac
  while(off < jobN){
    bi  = off/(12*12);
    dc  = off%(12*12);
    ga  = bi * Nbfac + gi;
    si  = ga / Nvol;
    isp = ga % Nvol;

    o1 = dc * Nvol;

    P1[dc*bfac + bi] = p1[si][o1 + isp];
    P2[dc*bfac + bi] = p2[si][o1 + isp]; 
    P3[dc*bfac + bi] = p3[si][o1 + isp]; 
    off += Blocks;
  }
  }
  __syncthreads();

  const Int GROUP = 16;
  unsigned char pos[3*GROUP];
  signed   char g0[2*GROUP];
  Int ini, dv;
  const Int bsize   = insertion < 0 ? bfac : bfac * 144 * 3;
  //const Int Ngv_tot = insertion < 0 ? Ngv  : Ngv  * 3;
  Ty buf[bsize];

  off = tid;
  while(off < Ngv)
  {
    //const Int iv  = insertion < 0 ? off : off % Ngv;
    //const Int ins = off / Ngv;
    //  const Int iv = tid / 3;
    const Int iv = off;
    const Int MAX = ivP[iv];

    for(Int bi=0;bi<bsize;bi++){buf[bi] = 0.0;}

    const Int jobN = (MAX + GROUP - 1 )/GROUP;
    ini = 0; dv = GROUP;
    for(Int ji=0;ji<jobN;ji++){
      if(ini + dv >= MAX){dv = MAX - ini;}
      if(dv <= 0){continue;}
      const unsigned char* s0 = &(oPP[iv][ini*3]);
      const   signed char* s1 = &(gPP[iv][ini*2]);

      // load the index
      for(Int k=0;k<dv*3;k++){pos[k] = s0[k];}
      for(Int k=0;k<dv*2;k++){g0[k]  = s1[k];}

      for(Int idv=0;idv<dv;idv++){
        // get prop buf pos
        if(insertion < 0){
          const Ty* t1 = &P1[int(pos[idv*3+0]) * bfac + 0];
          const Ty* t2 = &P2[int(pos[idv*3+1]) * bfac + 0];
          const Ty* t3 = &P3[int(pos[idv*3+2]) * bfac + 0];
          const Ty gtem = Ty(g0[idv*2+0], g0[idv*2+1]);
          for(Int z=0;z<bfac;z++){
            buf[z] += (t1[z] * t2[z] * t3[z] * gtem);
          }
        }

        if(insertion >= 0){
          const Int o1 = pos[idv*3+0];
          const Int o2 = pos[idv*3+1];
          const Int o3 = pos[idv*3+2];

          const Ty* t1 = &P1[o1*bfac];
          const Ty* t2 = &P2[o2*bfac];
          const Ty* t3 = &P3[o3*bfac];
          const Ty gtem = Ty(g0[idv*2+0], g0[idv*2+1]);

          Ty* b1 = &buf[0*144*bfac + o1*bfac + 0];
          Ty* b2 = &buf[1*144*bfac + o2*bfac + 0];
          Ty* b3 = &buf[2*144*bfac + o3*bfac + 0];

          for(Int bi=0;bi<bfac;bi++){b1[bi] += (t2[bi] * t3[bi] * gtem);}
          for(Int bi=0;bi<bfac;bi++){b2[bi] += (t1[bi] * t3[bi] * gtem);}
          for(Int bi=0;bi<bfac;bi++){b3[bi] += (t1[bi] * t2[bi] * gtem);}

          //if(ins == 0){
          //  Ty* b1 = &buf[o1*bfac + 0];
          //  for(Int bi=0;bi<bfac;bi++){b1[bi] += (t2[bi] * t3[bi] * gtem);}
          //}
          //if(ins == 1){
          //  Ty* b2 = &buf[o2*bfac + 0];
          //  for(Int bi=0;bi<bfac;bi++){b2[bi] += (t1[bi] * t3[bi] * gtem);}
          //}
          //if(ins == 2){
          //  Ty* b3 = &buf[o3*bfac + 0];
          //  for(Int bi=0;bi<bfac;bi++){b3[bi] += (t1[bi] * t2[bi] * gtem);}
          //}
        }
      }
      ini += dv;
    }

    if(insertion <  0){
      for(Int bi=0;bi<bfac;bi++){
        const Long ga  = bi * Nbfac + gi;
        const Int  si  = ga / Nvol;
        const Long isp = ga % Nvol;
        resP[iv*Nsrc + si][isp] += buf[bi]; 
      }
    }

    if(insertion >= 0){
      for(Int bi=0;bi<bfac;bi++){
        const Long ga  = bi * Nbfac + gi;
        const Int  si  = ga / Nvol;
        const Long isp = ga % Nvol;
        for(Int ins=0;ins<3;ins++)
        {
          Ty* rp = resP[ins*Ngv*Nsrc + iv*Nsrc + si];
          for(Int dc=0;dc<144;dc++)
          {
            rp[dc * Nvol + isp] += buf[(ins*144 + dc) * bfac + bi]; 
            ////  rp[dc * Nvol + isp] += buf[dc * bfac + bi]; 
          }
        }
      }
    }
    off += Blocks;
  }
}
#endif

/*
  u0 (u1 G d) : (u0 (u1 G d))
  prop1 == (u0 u0)
  prop2 == (u1 u1) prop1 and prop2 will swap 
  prop3 == (d  d )
  Gres[0, ...] diagram (u0 u0) (u1 u1) (d d)
  Gres[1, ...] diagram (u0 u1) (u1 u0) (d d)

  A source gamma, B sink Gamma
  G projections with fermion sign : defines the sign of terms
  mL shape of diagram : defines the fermion connections
  mem orders Nsrc, 12 * 12, Nvol
*/
template<typename Ty, Int bfac, Int insertion>
void baryon_sparse_kernel(Ty** p1, Ty** p2, Ty** p3, Ty** resP,
  signed char** gPP, unsigned char** oPP, const Int* ivP,
  const Int Nsrc, const Long Nvol, const Int Ngv)
{
  //Qassert(bfac < 255);// int8_t limit and cache sizes ?
  const Long Ntotal  = Nsrc * Nvol;
  if(Ntotal % bfac != 0){
    qmessage("Please correct your bfac! Nsrc %5d, Nvol %5d, bfac %d. \n", int(Nsrc), int(Nvol), bfac );
    abort_r("");
  }
  const Long Nbfac = Ntotal/bfac;
  Qassert(insertion == -1 or insertion == 1);
  //Qassert(insertion == -1 or insertion == 0 or insertion == 1 or insertion == 2);
  //qmessage("sites %5d %5d \n", int(Nbfac), int(bfac));

  bool job_done = false;

  #ifdef QLAT_USE_ACC
  /*
    baryon_sparse_global : if Ngv not large enough
    baryon_sparse_global_insertion : Ngv large same performance
  */
  if(insertion < 0){
    const Int nt = 16;
    const Int Blocks = nt*bfac;
    dim3 dimBlock(    nt, bfac, 1);
    dim3 dimGrid(  Nbfac,  1, 1);
    baryon_sparse_global<Ty, bfac, Blocks><<<dimGrid, dimBlock>>>(p1, p2, p3, resP, gPP, oPP, ivP, Nsrc, Nvol, Ngv);
    job_done = true;
  }
  // slower than naive qacc ........ due to write out not in thread sequence
  //if(insertion >= 0){
  //  //const Int Ngv_tot = insertion < 0 ? Ngv : Ngv * 3 ;
  //  //const Int nt = ( Ngv_tot + bfac -1 ) / bfac;
  //  const Int nt = ( Ngv + bfac -1 ) / bfac;
  //  dim3 dimBlock(    nt, bfac, 1);
  //  dim3 dimGrid(  Nbfac,  1, 1);
  //  baryon_sparse_global_insertion<Ty, bfac, insertion><<<dimGrid, dimBlock>>>(p1, p2, p3, resP, gPP, oPP, ivP, Nsrc, Nvol, Ngv);
  //}
  qacc_barrier(dummy);
  #endif
  if(!job_done){
  // each gi and bi will control one tzyx
  Qassert(Nvol % bfac == 0);
  qacc_for(gi, Nbfac, {
    const Long Nbuf = insertion < 0 ? bfac : bfac * 3 * 144;
    //const Long Nbuf = bfac;// no buf needed when doing insertion
    Ty buf[Nbuf+1];
    Ty P1[ bfac*12*12];
    Ty P2[ bfac*12*12];
    Ty P3[ bfac*12*12];

    const Long si   = (gi * bfac + 0 ) / Nvol;
    const Long isp  = (gi * bfac + 0 ) % Nvol;
    {
      Ty* f1 = &p1[si][           isp + 0];
      Ty* f2 = &p2[si][           isp + 0];
      Ty* f3 = &p3[si][           isp + 0];
      for(Int dc=0;dc<12*12;dc++){
        const Long o0 = dc * bfac;
        const Long o1 = dc * Nvol;
        for(Int bi=0;bi<bfac;bi++)P1[o0 + bi] = f1[o1 + bi];
        for(Int bi=0;bi<bfac;bi++)P2[o0 + bi] = f2[o1 + bi];
        for(Int bi=0;bi<bfac;bi++)P3[o0 + bi] = f3[o1 + bi];
      }
    }
    for(Int iv=0;iv<Ngv;iv++)
    {
      if(insertion < 0){
        for(Int bi=0;bi<Nbuf;bi++){buf[bi] = 0;}
        for(Int off=0;off<ivP[iv];off++)
        {
          const Ty* t1 = &P1[int(oPP[iv][off*3+0])*bfac];
          const Ty* t2 = &P2[int(oPP[iv][off*3+1])*bfac];
          const Ty* t3 = &P3[int(oPP[iv][off*3+2])*bfac];
          const Ty gtem = Ty(gPP[iv][off*2+0], gPP[iv][off*2+1]);
          for(Int bi=0;bi<bfac;bi++)
          {
            buf[bi] += (t1[bi] * t2[bi] * t3[bi] * gtem);
          }
        }

        Ty* r0 = &resP[iv*Nsrc + si][isp + 0];
        for(Int bi=0;bi<bfac; bi++){
          r0[bi] += buf[bi];
        }
      }

      if(insertion >= 0){
        for(Int bi=0;bi<Nbuf;bi++){buf[bi] = 0;}
        for(Int off=0;off<ivP[iv];off++)
        {
          const Int o1 = int(oPP[iv][off*3+0]);
          const Int o2 = int(oPP[iv][off*3+1]);
          const Int o3 = int(oPP[iv][off*3+2]);

          const Ty* t1 = &P1[o1*bfac];
          const Ty* t2 = &P2[o2*bfac];
          const Ty* t3 = &P3[o3*bfac];
          const Ty gtem = Ty(gPP[iv][off*2+0], gPP[iv][off*2+1]);
          Ty* b1 = &buf[0*144*bfac + o1*bfac + 0];
          Ty* b2 = &buf[1*144*bfac + o2*bfac + 0];
          Ty* b3 = &buf[2*144*bfac + o3*bfac + 0];
          for(Int bi=0;bi<bfac;bi++){b1[bi] += (t2[bi] * t3[bi] * gtem);}
          for(Int bi=0;bi<bfac;bi++){b2[bi] += (t1[bi] * t3[bi] * gtem);}
          for(Int bi=0;bi<bfac;bi++){b3[bi] += (t1[bi] * t2[bi] * gtem);}
        }

        for(Int ins=0;ins<3;ins++)
        {
          Ty* rB = resP[ins*Ngv*Nsrc + (iv*Nsrc + si)];
          for(Int dc=0;dc<144;dc++)
          {
            const Ty* br = &buf[ins*144*bfac + dc*bfac + 0];
                  Ty* rr = &rB[dc*Nvol + isp + 0];
            for(Int bi=0;bi<bfac;bi++){
              rr[bi] += br[bi];
            }
          }
        }
      }
    }
  });
  }
}

/////A source gamma, B sink Gamma, G projections with fermion sign, mL shape of diagram
/*
  p1, p2, p3 [Nsrc][Nvol]
*/
template <typename Ty, Int ins>
void baryon_sparseP(Ty** p1, Ty** p2, Ty** p3, Ty** resP, const Int Nsrc,
  ga_M &A, ga_M &B, qlat::vector<Ty > &GV, 
  qlat::vector<Int > &mLV, const Long Nvol, bool clear=true)
{
  TIMER_FLOPS("baryon_sparse");
  // double or single complex for buffers
  Qassert(sizeof(Ty) == 16 or sizeof(Ty) == 8);
  const Int Ngv = GV.size()/16;
  Qassert(GV.size()  == 16*Ngv);
  Qassert(mLV.size() == 3*Ngv);
  const Int Nins = ins < 0 ? 1 : 3;
  const Long Nres = Nins * Ngv * Nsrc;

  /*
    clear results
    will add to it if not clear
  */
  if(clear){
    for(Long si=0;si<Nres;si++){
      if(ins <  0){zero_Ty(resP[si], Nvol     , 1, QFALSE);}// clear vectors
      if(ins >= 0){zero_Ty(resP[si], Nvol*144 , 1, QFALSE);}// clear props
    }
    qacc_barrier(dummy);
  }

  /*
    prepare local contraction orders
  */
  signed char** gPP = NULL;
  unsigned char** oPP = NULL;
  Int*  ivP = NULL;
  std::vector<qlat::vector_gpu<signed char > > giEG;
  std::vector<qlat::vector_gpu<unsigned char   > > oiG;
  qlat::vector<Int > iv_size;

  giEG.resize(Ngv);
  oiG.resize(Ngv);
  iv_size.resize(Ngv);
  const Ty* GVP = GV.data();
  const Int*  mLP     = mLV.data();

  qlat::vector<Ty > epslV;epslV.resize(9);
  qlat::vector<Ty > gMap;
  qlat::vector<Int > IMap;

  for(Int i=0;i<3;i++){epslV[i*3+i]=0;epslV[i*3 + (i+1)%3]=1;epslV[i*3 + (i+2)%3]=-1;}
  gMap.resize(4*2);IMap.resize(4*2);
  for(Int i=0;i<4;i++){
    /////int j = + i;
    gMap[0*4+i] = A.g[i];
    gMap[1*4+i] = B.g[i];
    IMap[0*4+i] = A.ind[i];
    IMap[1*4+i] = B.ind[i];
  }

  const Ty* epsl = epslV.data();
  const Ty* gCA = &((gMap.data())[0*4]);
  const Ty* gCB = &((gMap.data())[1*4]);
  const Int* gIA = &((IMap.data())[0*4]);
  const Int* gIB = &((IMap.data())[1*4]);

  Int count_flops  = 0;
  {
    TIMER("Baryon color spin local orders");
    std::vector<std::vector<signed char > > giEL;giEL.resize(Ngv);//giEL.resize(  Ngv*Loff);
    std::vector<std::vector<unsigned char   > > oiL ;oiL.resize(Ngv );//oiL.resize(3*Ngv*Loff);
    Int bmL[3];
    Int nmL[3];
    for(Int iv=0;iv<Ngv;iv++)
    {
      oiL[iv].resize(0);
      giEL[iv].resize(0);
      const Int*      mL = &mLP[iv*3];

      const Ty* G  = &GVP[iv*16];

      for(Int a1=0;a1<3;a1++)
      for(Int ia=1;ia<3;ia++)
      for(Int b1=0;b1<3;b1++)
      for(Int ib=1;ib<3;ib++)
      {
        Int b2=(b1+ib)%3,b3=(b2+ib)%3;
        Int a2=(a1+ia)%3,a3=(a2+ia)%3;
        for(Int m2=0;m2<4;m2++)
        for(Int m1=0;m1<4;m1++)
        for(Int n2=0;n2<4;n2++)
        for(Int n1=0;n1<4;n1++)
        {
          const Ty Gtem =  G[m1*4+n1];
          const double norm = qlat::qnorm(Gtem);
          if(norm < 1e-20)continue;

          const Int m3 = gIA[m2];
          const Int n3 = gIB[n2];
          const Ty giE = epsl[a1*3 + a2]*epsl[b1*3 + b2]*gCA[m2]*gCB[n2]*G[m1*4+n1];
          nmL[0] = n1;nmL[1] = n2;nmL[2] = n3;
          bmL[0] = b1;bmL[1] = b2;bmL[2] = b3;
          const Int nm1 = nmL[mL[0]];
          const Int nm2 = nmL[mL[1]];
          const Int nm3 = nmL[mL[2]];

          const Int bm1 = bmL[mL[0]];
          const Int bm2 = bmL[mL[1]];
          const Int bm3 = bmL[mL[2]];

          const Int o1 = (m1*3+a1)*12+(nm1*3+bm1);
          const Int o2 = (m2*3+a2)*12+(nm2*3+bm2);
          const Int o3 = (m3*3+a3)*12+(nm3*3+bm3);

          ////buf += (P1[o1] * P2[o2] *P3[o3] * giE);
          oiL[iv].push_back(o1);
          oiL[iv].push_back(o2);
          oiL[iv].push_back(o3);
          giEL[iv].push_back((signed char)(giE.real()));
          giEL[iv].push_back((signed char)(giE.imag()));
          //giEL[iv].push_back(giE);
        }
      }
    }

    for(Int iv=0;iv<Ngv;iv++){giEG[iv].copy_from(giEL[iv]);}
    qlat::vector<signed char* > gP = EigenM_to_pointers(giEG);
    gPP = gP.data();

    for(Int iv=0;iv<Ngv;iv++){oiG[iv].copy_from(oiL[iv]);}
    qlat::vector<unsigned char* > oP = EigenM_to_pointers(oiG);
    oPP = oP.data();

    for(Int iv=0;iv<Ngv;iv++){
      iv_size[iv] = giEL[iv].size()/2;
      if(ins <  0){count_flops += (3 * 6 + 2) * iv_size[iv];}
      if(ins >= 0){count_flops += (3 * 2 * 6 + 3 * 2) * iv_size[iv];}
    }
    ivP = iv_size.data();
  }
  timer.flops  += Nsrc * Nvol * count_flops    ;

  {
  TIMER_FLOPS("baryon_sparse kernel");
  timer.flops  += Nsrc * Nvol * count_flops;
  const Int BFACG_DEFAULT = sizeof(Ty) == 16 ? BFACG_SHARED : BFACG_SHARED * 2;
  bool get = false;
  #define baryon_macros(ba) if(Nvol % ba == 0 and get == false){get = true; \
    baryon_sparse_kernel<Ty, ba, ins>(p1, p2, p3, resP, gPP, oPP, ivP, Nsrc, Nvol, Ngv);}

  #ifdef QLAT_USE_ACC
  // if it's not OK, then abort
  baryon_macros(BFACG_DEFAULT);
  #else
  const Int Bfac = BFACG_DEFAULT * 2;
  //// on CPU buffer size is not an issue and NTt not large enough ?
  baryon_macros(Bfac);
  baryon_macros(16);
  baryon_macros(10);
  baryon_macros(8 );
  baryon_macros(4 );
  baryon_macros(1 );
  #endif

  #undef baryon_macros
  Qassert(get == true);
  }
}

template <typename Ty>
void baryon_sparseI(Ty** p1, Ty** p2, Ty** p3, Ty** resP, const Int Nsrc,
  ga_M &A, ga_M &B, qlat::vector<Ty > &GV, 
  qlat::vector<Int > &mLV, const Long Nvol, const bool clear=1, const Int insertion = -1)
{
  if(insertion == -1)baryon_sparseP<Ty, -1>(p1, p2, p3, resP, Nsrc, A, B, GV, mLV, Nvol, clear);
  if(insertion ==  1)baryon_sparseP<Ty,  1>(p1, p2, p3, resP, Nsrc, A, B, GV, mLV, Nvol, clear);
}

template <typename Ty>
void baryon_sparse(std::vector<FieldG<Ty> >& prop1, std::vector<FieldG<Ty> >& prop2, std::vector<FieldG<Ty> >& prop3,
  std::vector<FieldG<Ty> >& res, ga_M &A, ga_M &B, qlat::vector<Ty > &G, qlat::vector<Int > &mL, const bool clear=true, const Int insertion = -1){
  const unsigned int Nsrc = prop1.size();
  Qassert(Nsrc != 0);
  Qassert(prop2.size() == Nsrc and prop3.size() == Nsrc);
  const Int Ngv = G.size()/16;
  const Long Ndc = 12 * 12;
  for(unsigned int si=0;si<Nsrc;si++){
    Qassert(prop1[si].initialized and prop2[si].initialized and prop3[si].initialized);
    Qassert(prop1[si].multiplicity == Ndc and prop2[si].multiplicity == Ndc and prop3[si].multiplicity == Ndc);
    Qassert(prop1[si].mem_order == QLAT_OUTTER and prop2[si].mem_order == QLAT_OUTTER and prop3[si].mem_order == QLAT_OUTTER);
  }
  // will results in 48 props.....
  const Geometry& geo = prop1[0].geo();
  const Int Nins = insertion < 0 ? 1 : 3;
  const Long Nres = Nins * Ngv * Nsrc;
  
  if(Long(res.size()) != Nres){res.resize(0);res.resize(Nres);}
  for(Long ri=0;ri < Nres;ri++)
  {
    if(!res[ri].initialized){
      if(insertion < 0){
        res[ri].init(geo, 1, QMGPU, QLAT_OUTTER);
      }else{
        res[ri].init(geo, 144, QMGPU, QLAT_OUTTER);
      }
      set_zero(res[ri]);
    }
    if(insertion < 0){Qassert(res[ri].multiplicity == 1  );}
    else{             Qassert(res[ri].multiplicity == 144);}
  }
  const Long Nvol = geo.local_volume();
  vector<Ty* > p1 = FieldG_to_pointers(prop1 );
  vector<Ty* > p2 = FieldG_to_pointers(prop2 );
  vector<Ty* > p3 = FieldG_to_pointers(prop3 );
  vector<Ty* > rP = FieldG_to_pointers(res   );
  baryon_sparseI(p1.data(), p2.data(), p3.data(), rP.data(), Nsrc, A, B, G, mL, Nvol, clear, insertion);
}

template <typename Ty>
void baryon_sparse(FieldG<Ty>& prop1, FieldG<Ty>& prop2, FieldG<Ty>& prop3,
  std::vector<FieldG<Ty> >& res, ga_M &A, ga_M &B, qlat::vector<Ty > &G, qlat::vector<Int > &mL, const bool clear=true, const Int insertion = -1){
  const unsigned int Nsrc = 1;
  const Int Ngv = G.size()/16;
  const Long Ndc = 12 * 12;
  Qassert(prop1.initialized and prop2.initialized and prop3.initialized);
  Qassert(prop1.multiplicity == Ndc and prop2.multiplicity == Ndc and prop3.multiplicity == Ndc);
  Qassert(prop1.mem_order == QLAT_OUTTER and prop2.mem_order == QLAT_OUTTER and prop3.mem_order == QLAT_OUTTER);
  const Geometry& geo = prop1.geo();
  const Int Nins = insertion < 0 ? 1 : 3;
  const Long Nres = Nins * Ngv * Nsrc;

  if(res.size() != Nres){res.resize(0);res.resize(Nres);}
  for(Long ri=0;ri < Nres;ri++)
  {
    if(!res[ri].initialized){
      if(insertion < 0){
        res[ri].init(geo, 1, QMGPU, QLAT_OUTTER);
      }else{
        res[ri].init(geo, 144, QMGPU, QLAT_OUTTER);
      }
      set_zero(res[ri]);
    }
    if(insertion < 0){Qassert(res[ri].multiplicity == 1  );}
    else{             Qassert(res[ri].multiplicity == 144);}
  }

  const Long Nvol = geo.local_volume();
  vector<Ty* > p1;p1.resize(1);
  vector<Ty* > p2;p2.resize(1);
  vector<Ty* > p3;p3.resize(1);
  vector<Ty* > rP = FieldG_to_pointers(res   );

  p1[0] = (Ty*) get_data(prop1 ).data();
  p2[0] = (Ty*) get_data(prop2 ).data();
  p3[0] = (Ty*) get_data(prop3 ).data();
  baryon_sparseI(p1.data(), p2.data(), p3.data(), rP.data(), Nsrc, A, B, G, mL, Nvol, clear, insertion);
}

/*
  input 
    props, sections for positive and negative projections
  output 
    [Nsrc, ud(2), prj(4)]
*/
template <typename Ty>
void baryon_sparse_insertion(std::vector<FieldG<Ty> >& prop1, std::vector<FieldG<Ty> >& prop2, std::vector<FieldG<Ty> >& prop3,
  std::vector<FieldG<Ty> >& res, const vector<Int >& map_sec, std::vector<FieldG<Ty> >& bufV){

  TIMER("baryon_sparse_insertion");
  const Int Nsrc = prop1.size();
  const Geometry& geo = prop1[0].geo();
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);

  // results ud, 4 poly
  if(Long(res.size()) != Nsrc * 2 * 4){res.resize(Nsrc * 2 * 4);}
  for(unsigned int i=0;i<res.size();i++){
    if(!res[i].initialized)res[i].init(geo, 144, QMGPU, QLAT_OUTTER);
    Qassert(res[i].mem_order == QLAT_OUTTER and res[i].multiplicity == 144);
    // clean results
    set_zero(res[i]);
  }
  vector<Int > maps;maps.resize(map_sec.size());
  for(unsigned int i=0;i<maps.size();i++){
    maps[i] = map_sec[i];
  }

  ga_matrices_cps   ga_cps;
  std::vector<ga_M > gP(8);
  gP[0] = ga_cps.ga[0][0];gP[1] = ga_cps.ga[0][4];
  
  for(Int i=1;i<4;i++){
    gP[i*2 + 0] =                   ga_cps.ga[i][5]; //  need reverse or not
    gP[i*2 + 1] = ga_cps.ga[0][4] * ga_cps.ga[i][5];
  }
  
  std::vector<Int > udL = {0,0,  1};
  std::vector<Ty > factor(4);
  factor[0] = 1.0/2.0;
  //for(Int i=1;i<4;i++){factor[i] = Ty(0.0, 1.0/2.0);}
  for(Int i=1;i<4;i++){factor[i] = Ty(1.0/2.0, 0.0);}  // complex sign convention here?

  const Int Nt = fd.Nt;
  const Long Nvol = geo.local_volume();
  const Long Nxyz = geo.local_volume() / Nt;
  //const Long Ntot = geo.local_volume() * 144 / Nt;

  const Int Ngroup = 2;
  Qassert(16 % Ngroup == 0);
  vector<Ty > Gp ;Gp.resize(Ngroup*16);
  vector<Ty > Gm ;Gm.resize(Ngroup*16);
  vector<Int      > mp;mp.resize(Ngroup*3);
  vector<Int      > mm;mm.resize(Ngroup*3);
  ga_M &ga2 = ga_cps.ga[1][3];
  ga_M &ga1 = ga_cps.ga[1][3];

  /*
    shrink the buffer usage
  */
  for(Int gi=0;gi<16/Ngroup;gi++){
    //std::vector<Int > curr_ops;
    //curr_ops.resize(Ngroup);
    //for(Int j=0;j<Ngroup;j++){
    //  curr_ops.push_back(gi*Ngroup + j );
    //}
    const Int op_min = gi*Ngroup + 0;
    const Int op_max = op_min    + Ngroup;

    clear_qv(Gp);clear_qv(mp);
    clear_qv(Gm);clear_qv(mm);
    Int count = 0;
    for(Int ind2=0;ind2<4;ind2++)
    for(Int ind1=0;ind1<4;ind1++)
    {
      const Int ioff = ind2 * 4 + ind1;
      if(ioff < op_min or ioff >= op_max){continue;}

      Gp[count*16 + ioff] = +1.0;
      mp[count*3  + 0   ] =   0 ;
      mp[count*3  + 1   ] =   1 ;
      mp[count*3  + 2   ] =   2 ;
    
      Gm[count*16 + ioff] = -1.0;
      mm[count*3  + 0   ] =   1 ;
      mm[count*3  + 1   ] =   0 ;
      mm[count*3  + 2   ] =   2 ;
      count += 1;
    }
    Qassert(count == Ngroup);

    // bufV sizes ins * Ngroup
    // for(unsigned int iv=0;iv<bufV.size();iv++){set_zero(bufV[iv]);}
    baryon_sparse(prop1, prop2, prop3, bufV, ga2, ga1, Gp, mp, true , 1);
    baryon_sparse(prop1, prop2, prop3, bufV, ga2, ga1, Gm, mm, false, 1);

    //// choose sections and baryon polarization projections
    {
    TIMER("Sum baryon_sparse_insertion");
    for(Int srci=0;srci<Nsrc;srci++)
    for(Int ins=0;ins<3;ins++)
    for(Int prj=0;prj<4;prj++)
    {
      const Int ud = udL[ins];
      for(Int si = 0; si < 2 ;si++)//// 1 + r4 sectors
      for(Int gd=0;gd<4;gd++)
      {
        //  const Int bind = gP[prj*2 + si].ind[gd]*4 + gd ; //   need reverse or not
        const Int bind = gd*4 + gP[prj*2 + si].ind[gd] - op_min; //   need reverse or not
        if(bind < 0 or bind >= Ngroup){continue;}
        //if(bind < op_min or bind >= op_max){continue;}
        //Qassert(bind > 0 )
    
        const Ty sign  = gP[prj*2 + si].g[gd]           * factor[prj];
        //  Ty* sP = (Ty*) qlat::get_data(bufV[(ins*16 + bind)*Nsrc + srci]).data();
        const Ty* sP = (Ty*) qlat::get_data(bufV[(ins*Ngroup + bind)*Nsrc + srci]).data();
        Ty* rP = (Ty*) qlat::get_data(res[  srci*2*4 + (ud*4 + prj)]).data();
        const Int tini = fd.init;
        //  if(seci == 1 and si == 1){sign_use = (-1.0) * sign;} //   reverse sign if seci == 1
        qacc_for(isp, Nvol, {
          const Long ti  = isp / Nxyz;
          //const Long idx = isp % Nxyz;
          const Int seci = maps[ti + tini]%2;
          const Ty sign_use = (seci == 1 and si == 1) ? (-1.0) * sign : sign;
          for(Int ops=0;ops<144;ops++){
            rP[ops*Nvol + isp] += sign_use * sP[ops*Nvol + isp];
          }
        });
      }
    }
    }
  }

}


}

#endif

