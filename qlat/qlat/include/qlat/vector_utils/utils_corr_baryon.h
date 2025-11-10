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

/*
  simple gpu version, under Yi-bo's implimentation
  Proton contractions, for checks purpose only
  u0 (u1 G d) : (u0 (u1 G d))
  prop1 == (u0 u0)
  prop2 == (u1 u1) prop1 and prop2 will swap 
  prop3 == (d  d )
  Gres[0, ...] diagram (u0 u0) (u1 u1) (d d)
  Gres[1, ...] diagram (u0 u1) (u1 u0) (d d)
*/
template <typename Ty>
void proton_vectorE_gwu(std::vector<qpropT >& prop1, std::vector<qpropT >& prop2, std::vector<qpropT >& prop3,
  qlat::vector<Ty > &res, const ga_M &ga2,const Int ind2, const ga_M &ga1, const Int ind1, Int clear=1){
  TIMER("Proton_vectorE");
  if(prop1.size() == 0){res.resize(0); return ;}
  const qlat::Geometry& geo = prop1[0].geo();
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);

  /* 
    ga2/ind2 for source, gam1/ind1 for sink
    "[]|+N" type diagram
    check_prop_size(prop1);check_prop_size(prop2);check_prop_size(prop3);
  */
  Int NTt  = fd.Nv[3];
  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  Int nmass = prop1.size();
  Qassert(prop1.size() == prop2.size());
  Qassert(prop1.size() == prop3.size());
  if(clear == 1){ini_resE(res, nmass, fd);}
  if(clear == 0){Qassert(res.size() == Long(nmass*NTt * Nxyz));}
    
  //  Prop format, src d-4, c-3, sink d-4, c-3, Nt, EigenVTa<Nxyz>
  if(res.size()%NTt !=0 or res.size()==0){qmessage("Size of res wrong. \n");Qassert(false);}

  qlat::vector<Ty* > p1 = FieldM_to_pointers(prop1);
  // swap definitions to match qlat, usually p2 = p3
  qlat::vector<Ty* > p3 = FieldM_to_pointers(prop2); // up  quark
  qlat::vector<Ty* > p2 = FieldM_to_pointers(prop3); //down quark

  Ty epsl[3][3];
  for(Int i=0;i<3;i++)for(Int j=0;j<3;j++){epsl[i][j] = 0;}
  for(Int i=0;i<3;i++){epsl[i][i]=0;epsl[i][(i+1)%3]=Ty(1,0);epsl[i][(i+2)%3]=Ty(-1,0);}

  for(Int d2=0;d2<4;d2++)
  for(Int c21=0;c21<3;c21++)
  for(Int ib=1;ib<3;ib++)
  {
    Int c22=(c21+ib)%3,c23=(c22+ib)%3;
    for(Int d1=0;d1<4;d1++)
    for(Int c11=0;c11<3;c11++)
    for(Int ia=1;ia<3;ia++)
    {
      Int c12=(c11+ia)%3,c13=(c12+ia)%3;
      Ty giE = epsl[c11][c12]*epsl[c21][c22]*ga1.g[d1]*ga2.g[d2];

      #pragma omp parallel for
      for(Int ji=0;ji<nmass*NTt;ji++)
      {
        Int massi = ji/NTt;
        Int ti    = ji%NTt;

        Int m1 = (ind2*3+c21)*12+ind1*3+c11;
        Int m2 = (ga2.ind[d2]*3+c22)*12+d1*3+c12;
        Int m3 = (d2*3+c23)*12+ga1.ind[d1]*3+c13;

        Int n1 = (ind2*3+c21)*12+ga1.ind[d1]*3+c11;
        Int n2 = (ga2.ind[d2]*3+c22)*12+d1*3+c12;
        Int n3 = (d2*3+c23)*12+ind1*3+c13;

        Ty* tp1 = &p1[massi][(m1*NTt+ti)*Nxyz];
        Ty* tp2 = &p2[massi][(m2*NTt+ti)*Nxyz];
        Ty* tp3 = &p3[massi][(m3*NTt+ti)*Nxyz];

        Ty* tn1 = &p1[massi][(n1*NTt+ti)*Nxyz];
        Ty* tn2 = &p2[massi][(n2*NTt+ti)*Nxyz];
        Ty* tn3 = &p3[massi][(n3*NTt+ti)*Nxyz];
        Ty* tr0 = &(res.data()[(massi*NTt + ti)*Nxyz]);

        #if USEQACC==1
        qacc_forNB(i, Long(Nxyz),{ tr0[i] -= ((tp1[i]*tp2[i]*tp3[i] + tn1[i]*tn2[i]*tn3[i])*giE); });
        #else
        EAy vp1(tp1,Nxyz);
        EAy vp2(tp2,Nxyz);
        EAy vp3(tp3,Nxyz);
        EAy vn1(tn1,Nxyz);
        EAy vn2(tn2,Nxyz);
        EAy vn3(tn3,Nxyz);
        EAy vr0(tr0,Nxyz);
        vr0 -= ((vp1*vp2*vp3 + vn1*vn2*vn3)*giE);
        #endif

      }
      qacc_barrier(dummy);
    }
  }

}

///proton sector corr with prop gwu convention
template <typename Ty>
void proton_vectorE_gwu(std::vector<qpropT >& prop1, std::vector<qpropT >& prop2, std::vector<qpropT >& prop3,
        qlat::vector<Ty > &res, const ga_M &ga1,Int t0,Int dT,Int clear=1,Int oppo=0){
  TIMER("Proton_vectorE");
  if(prop1.size() == 0){res.resize(0); return;}

  Int nmass = prop1.size();
  if(prop1.size() == 0){res.resize(0); return ;}
  const qlat::Geometry& geo = prop1[0].geo();
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);

  Int NTt  = fd.Nv[3];
  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];

  Qassert(prop1.size() == prop2.size());
  Qassert(prop1.size() == prop3.size());

  if(clear == 1){ini_resE(res,nmass,fd);}

  //int nv = res.size();int Nsize = res[0].size();
  qlat::vector<Ty >  resE0;resE0.resize(res.size());
  qlat::vector<Ty >  resE1;resE1.resize(res.size());
  ////qlat::set_zero(resE0);qlat::set_zero(resE1);

  proton_vectorE_gwu(prop1,prop2,prop3,resE0,ga1,0,ga1,0,1);
  proton_vectorE_gwu(prop1,prop2,prop3,resE0,ga1,1,ga1,1,0);
  proton_vectorE_gwu(prop1,prop2,prop3,resE1,ga1,2,ga1,2,1);
  proton_vectorE_gwu(prop1,prop2,prop3,resE1,ga1,3,ga1,3,0);

  std::vector<Int > map_sec = get_map_sec(dT,fd.nt);
  //////int Nt = fd.Nt;

  /////int t0 = 0;
  Int nt = fd.nt;
  ///for(Int massi=0;massi<nmass;massi++)
  ///for(Int ti = 0;ti<Nt;ti++)
  #pragma omp parallel for
  for(Int ji=0;ji<nmass*NTt;ji++)
  {
    Int massi = ji/NTt;
    Int ti    = ji%NTt;
    Int t = ti + fd.Pos0[fd.rank][3];
    Ty* tr0 = &(res.data()[(massi*NTt+ti)*Nxyz]);
    Ty* tv0 = &(resE0.data()[(massi*NTt+ti)*Nxyz]);
    Ty* tv1 = &(resE1.data()[(massi*NTt+ti)*Nxyz]);

    #if USEQACC==0
    EAy r0(tr0,Nxyz);
    EAy v0(tv0,Nxyz);
    EAy v1(tv1,Nxyz);
    #endif

    if(map_sec[(t-t0+nt)%nt]%2==0)
    {
      #if USEQACC==1
      if(oppo==0)qacc_forNB(i, Long(Nxyz), { tr0[i] += tv0[i]; });
      if(oppo==1)qacc_forNB(i, Long(Nxyz), { tr0[i] += tv1[i]; });
      #else
      if(oppo==0){r0 += v0;}
      if(oppo==1){r0 += v1;}
      #endif

    }
    if(map_sec[(t-t0+nt)%nt]%2==1)
    {
      #if USEQACC==1
      if(oppo==0)qacc_forNB(i, Long(Nxyz), { tr0[i] += tv1[i]; });
      if(oppo==1)qacc_forNB(i, Long(Nxyz), { tr0[i] += tv0[i]; });
      #else
      if(oppo==0){r0 += v1;}
      if(oppo==1){r0 += v0;}
      #endif

    }
  }
  qacc_barrier(dummy);
}

////container
template <typename Ty>
void proton_corrE(std::vector<qpropT >& prop1, std::vector<qpropT >& prop2, std::vector<qpropT >& prop3,
  const ga_M &ga2,const Int ind2, const ga_M &ga1,const Int ind1,
  qlat::vector<Ty > &res, Int clear=1,const Coordinate& mom = Coordinate()){
  qlat::vector<Ty > resE;
  proton_vectorE_gwu(prop1,prop2,prop3,ga2,ind2,ga1,ind1,resE,1);
  vec_corrE(resE,res,clear,mom);
}

//  container
template <typename Ty>
void proton_corrE(std::vector<qpropT >& prop1, std::vector<qpropT >& prop2, std::vector<qpropT >& prop3,
 qlat::vector<Ty > &res, const ga_M &ga1,const Int t0,const Int dT,Int clear=1,const Coordinate& mom = Coordinate()){
  qlat::vector<Ty >  resE;
  proton_vectorE_gwu(prop1,prop2,prop3,resE, ga1, t0,dT,1);

  vec_corrE(resE,res,clear,mom);
}

/*  
    simple gpu version
    A source gamma, B sink Gamma
    G projections with fermion sign : defines the sign of terms
    mL shape of diagram : defines the fermion connections
*/
template <typename Ty>
void baryon_vectorE(std::vector<qpropT >& prop1, std::vector<qpropT >& prop2, std::vector<qpropT >& prop3,
  qlat::vector<Ty > &res, ga_M &A, ga_M &B, qlat::vector<Ty > &G, qlat::vector<Int > &mL, Int clear=1){
  TIMER("Proton_vectorE");
  if(prop1.size() == 0){res.resize(0); return ;}
  ////check_prop_size(prop1);check_prop_size(prop2);check_prop_size(prop3);
  const qlat::Geometry& geo = prop1[0].geo();
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);

  Int NTt  = fd.Nv[3];
  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  Int nmass = prop1.size();
  Qassert(prop1.size() == prop2.size());
  Qassert(prop1.size() == prop3.size());
  Qassert(G.size()  == 16);
  Qassert(mL.size() == 3);
  Qassert(fd.order_ch == 0);
  if(clear == 1){ini_resE(res,nmass,fd);}

  ////if(res.size()%NTt !=0 or res.size()==0){qmessage("Size of res wrong. \n");Qassert(false);}
  if(res.size()==0){qmessage("Size of res wrong. \n");Qassert(false);}

  qlat::vector<Ty* > p1 = FieldM_to_pointers(prop1);
  qlat::vector<Ty* > p2 = FieldM_to_pointers(prop2);
  qlat::vector<Ty* > p3 = FieldM_to_pointers(prop3);

  Ty epsl[3][3];
  for(Int i=0;i<3;i++)for(Int j=0;j<3;j++){epsl[i][j] = 0;}
  for(Int i=0;i<3;i++){epsl[i][i]=0;epsl[i][(i+1)%3]=1;epsl[i][(i+2)%3]=-1;}

  /*
    mL = {};
    std::vector<Int > mL;mL.resize(3);
    mL[0] = 0;mL[1] = 1;mL[2] = 2;
  */
  std::vector<Int > nmL;nmL.resize(3);
  std::vector<Int > bmL;bmL.resize(3);

  {
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
        Ty Gv =  G[m1*4+n1];
        double norm = std::sqrt(Gv.real()*Gv.real() + Gv.imag()*Gv.imag());
        if(norm < 1e-20)continue;

        Int m3 = A.ind[m2];
        Int n3 = B.ind[n2];
        Ty giE = epsl[a1][a2]*epsl[b1][b2]*A.g[m2]*B.g[n2]*G[m1*4+n1];
        nmL[0] = n1;nmL[1] = n2;nmL[2] = n3;
        bmL[0] = b1;bmL[1] = b2;bmL[2] = b3;
        Int nm1 = nmL[mL[0]];
        Int nm2 = nmL[mL[1]];
        Int nm3 = nmL[mL[2]];

        Int bm1 = bmL[mL[0]];
        Int bm2 = bmL[mL[1]];
        Int bm3 = bmL[mL[2]];

        #pragma omp parallel for
        for(Int ji=0;ji<nmass*NTt;ji++)
        {
          Int massi = ji/NTt;
          Int ti    = ji%NTt;

          Int o1 = (m1*3+a1)*12+(nm1*3+bm1);
          Int o2 = (m2*3+a2)*12+(nm2*3+bm2);
          Int o3 = (m3*3+a3)*12+(nm3*3+bm3);

          Ty* tp1 = &p1[massi][(o1*NTt+ti)*Nxyz];
          Ty* tp2 = &p2[massi][(o2*NTt+ti)*Nxyz];
          Ty* tp3 = &p3[massi][(o3*NTt+ti)*Nxyz];
          Ty* tr0 = &(res.data()[(massi*NTt + ti)*Nxyz]);

          #if USEQACC==1
          qacc_forNB(i, Long(Nxyz),{ tr0[i] += (tp1[i]*tp2[i]*tp3[i] * giE); });
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
}

/*  
  3pt insertion version
  A source gamma, B sink Gamma, G projections with fermion sign, mL shape of diagram
*/
template <typename Ty>
void baryon_vectorE(std::vector<qpropT >& prop1, std::vector<qpropT >& prop2, std::vector<qpropT >& prop3,
  std::vector<qpropT >& resP, ga_M &A, ga_M &B, qlat::vector<Ty > &G, qlat::vector<Int > &mL, Int insertion,Int clear=1)
{
  TIMER("Baryon_vectorE_insertion");
  if(prop1.size() == 0){resP.resize(0); return ;}
  ////check_prop_size(prop1);check_prop_size(prop2);check_prop_size(prop3);
  const qlat::Geometry& geo = prop1[0].geo();
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);
  Qassert(insertion == 0 or insertion == 1 or insertion == 2);

  Int NTt  = fd.Nv[3];
  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  const Int nmass = prop1.size();
  Qassert(prop1.size() == prop2.size());
  Qassert(prop1.size() == prop3.size());
  Qassert(G.size()  == 16);
  Qassert(mL.size() == 3);
  Qassert(fd.order_ch == 0);

  if(clear==1){
    if(int(resP.size()) != nmass){resP.resize(nmass);}
    for(Int mi=0;mi<nmass;mi++){
      if(!resP[mi].initialized){resP[mi].init(geo);}
      else{qlat::set_zero(resP[mi]);}
    }
  }
  Qassert(int(resP.size()) == nmass);
  for(Int mi=0;mi<nmass;mi++){Qassert(resP[mi].initialized);}

  /////check_prop_size(resP);

  qlat::vector<Ty* > p1 = FieldM_to_pointers(prop1);
  qlat::vector<Ty* > p2 = FieldM_to_pointers(prop2);
  qlat::vector<Ty* > p3 = FieldM_to_pointers(prop3);
  qlat::vector<Ty* > rP = FieldM_to_pointers(resP);

  Ty epsl[3][3];
  for(Int i=0;i<3;i++)for(Int j=0;j<3;j++){epsl[i][j] = 0;}
  for(Int i=0;i<3;i++){epsl[i][i]=0;epsl[i][(i+1)%3]=1;epsl[i][(i+2)%3]=-1;}

  std::vector<Int > nmL;nmL.resize(3);
  std::vector<Int > bmL;bmL.resize(3);

  {
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
        Ty Gv =  G[m1*4+n1];
        double norm = std::sqrt(Gv.real()*Gv.real() + Gv.imag()*Gv.imag());
        if(norm < 1e-20)continue;

        Int m3 = A.ind[m2];
        Int n3 = B.ind[n2];
        Ty giE = epsl[a1][a2]*epsl[b1][b2]*A.g[m2]*B.g[n2]*G[m1*4+n1];
        nmL[0] = n1;nmL[1] = n2;nmL[2] = n3;
        bmL[0] = b1;bmL[1] = b2;bmL[2] = b3;
        Int nm1 = nmL[mL[0]];
        Int nm2 = nmL[mL[1]];
        Int nm3 = nmL[mL[2]];

        Int bm1 = bmL[mL[0]];
        Int bm2 = bmL[mL[1]];
        Int bm3 = bmL[mL[2]];

        #pragma omp parallel for
        for(Int ji=0;ji<nmass*NTt;ji++)
        {
          Int massi = ji/NTt;
          Int ti    = ji%NTt;

          Int o1 = (m1*3+a1)*12+(nm1*3+bm1);
          Int o2 = (m2*3+a2)*12+(nm2*3+bm2);
          Int o3 = (m3*3+a3)*12+(nm3*3+bm3);

          Int r0 = 0;
          if(insertion == 0){r0 = (m1*3+a1)*12+(nm1*3+bm1);}
          if(insertion == 1){r0 = (m2*3+a2)*12+(nm2*3+bm2);}
          if(insertion == 2){r0 = (m3*3+a3)*12+(nm3*3+bm3);}

          Ty* tp1 = &p1[massi][(o1*NTt+ti)*Nxyz];
          Ty* tp2 = &p2[massi][(o2*NTt+ti)*Nxyz];
          Ty* tp3 = &p3[massi][(o3*NTt+ti)*Nxyz];
          Ty* tr0 = &rP[massi][(r0*NTt+ti)*Nxyz];

          #if USEQACC==1
          if(insertion == 0)qacc_forNB(i, Long(Nxyz),{ tr0[i] += (tp2[i]*tp3[i] * giE); });
          if(insertion == 1)qacc_forNB(i, Long(Nxyz),{ tr0[i] += (tp1[i]*tp3[i] * giE); });
          if(insertion == 2)qacc_forNB(i, Long(Nxyz),{ tr0[i] += (tp1[i]*tp2[i] * giE); });
          #else
          EAy vp1(tp1,Nxyz);
          EAy vp2(tp2,Nxyz);
          EAy vp3(tp3,Nxyz);
          EAy vr0(tr0,Nxyz);
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

template <typename Ty>
void baryon_vectorE_ud_insert(std::vector<qpropT >& prop1, std::vector<qpropT >& prop2, std::vector<qpropT >& prop3,
  std::vector< std::vector<qpropT > >& resP, const Int ud, const std::vector<Int >& map_sec, std::vector< std::vector<qpropT > >& buf) 
{
  if(prop1.size() == 0){resP.resize(0); return ;}
  const qlat::Geometry& geo = prop1[0].geo();
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);

  const Long Nt = fd.Nt;
  const Long Nxyz = fd.Nx * fd.Ny * fd.Nz;
  ga_matrices_cps   ga_cps;
  ga_M& A = ga_cps.ga[1][3];
  ga_M& B = ga_cps.ga[1][3];

  ////get all 4 prj
  if(buf.size() != 16){buf.resize(16);}
  for(Int i=0;i<16;i++){
    init_qpropT(buf[i], prop1.size(), geo);
    clear_qpropT(buf[i]);
  }
  if(resP.size()!=4){resP.resize(4);}
  for(Int i=0;i<4;i++){
    init_qpropT(resP[i], prop1.size(), geo);
    {clear_qpropT(resP[i]);} ////must clear!
  }

  std::vector<Ty > factor(4);
  factor[0] = 1.0/2.0;
  ////for(Int i=1;i<4;i++){factor[i] = Ty(0.0, 1.0/2.0);}
  for(Int i=1;i<4;i++){factor[i] = Ty(1.0/2.0, 0.0);}  // complex sign convention here?

  std::vector<ga_M > gL(8);
  gL[0] = ga_cps.ga[0][0];gL[1] = ga_cps.ga[0][4];

  for(Int i=1;i<4;i++){
    gL[i*2 + 0] =                   ga_cps.ga[i][5]; //  need reverse or not
    gL[i*2 + 1] = ga_cps.ga[0][4] * ga_cps.ga[i][5];
  }

  std::vector<Int > udL = {0,0,  1};

  qlat::vector<Ty > G;G.resize(16);
  qlat::vector<Int > mL;mL.resize(3);

  //std::vector<Int > inL = {0,0,  1,1,  2,2,  3,3,  0,2,  1,3,  2,0,  3,1 };
  //for(Int ind=0;ind<8;ind++)
  //for(Int insertion=0;insertion<3;insertion++)
  //{   
  //  if(udL[insertion] == ud) 
  //  {   
  //    Int ind1 = inL[ind*2 + 0];
  //    Int ind2 = inL[ind*2 + 1];
  //    clear_qv(G);G[ind2*4 + ind1] = +1.0/2.0; //// (1+r4)/2, could be non-zero for all elements...
  //    mL[0] = 0;mL[1] = 1;mL[2] = 2;
  //    baryon_vectorE(prop1, prop2, prop3, resP[0], A, B, G, mL, insertion, 0); 
  //    ////baryon_vectorE(prop1,prop2,prop3, resE, ga2,ga1, G, mL, fd, 1); 
  //    clear_qv(G);G[ind2*4 + ind1] = -1.0/2.0;
  //    mL[0] = 1;mL[1] = 0;mL[2] = 2;
  //    baryon_vectorE(prop1, prop2, prop3, resP[0], A, B, G, mL, insertion, 0); //// no clear sum
  //  }   
  //}

  for(Int ind=0;ind<16;ind++)
  for(Int insertion=0;insertion<3;insertion++)
  {   
    if(udL[insertion] == ud) 
    {   
      Int ind1 = ind/4;
      Int ind2 = ind%4;
      qlat::vector<Ty > G;G.resize(16);
      qlat::vector<Int > mL;mL.resize(3);
      clear_qv(G);G[ind2*4 + ind1] = +1.0/1.0; //// (1+r4)/2, could be non-zero for all elements...
      mL[0] = 0;mL[1] = 1;mL[2] = 2;
      baryon_vectorE(prop1, prop2, prop3, buf[ind], A, B, G, mL, insertion, 0); 
      ////baryon_vectorE(prop1,prop2,prop3, resE, ga2,ga1, G, mL, fd, 1); 
      clear_qv(G);G[ind2*4 + ind1] = -1.0/1.0;
      mL[0] = 1;mL[1] = 0;mL[2] = 2;
      baryon_vectorE(prop1, prop2, prop3, buf[ind], A, B, G, mL, insertion, 0); //// no clear sum
    }   
  }

  for(Int prj=0;prj<4;prj++)
  {
    for(Int ti = 0; ti < Nt; ti++)
    {
      Int seci = map_sec[ti + fd.init]%2;
      for(Int si = 0; si < 2 ;si++)//// 1 + r4 sectors
      for(Int gd=0;gd<4;gd++)
      {
        Ty sign  = gL[prj*2 + si].g[gd]           * factor[prj];
        Int bind = gL[prj*2 + si].ind[gd]*4 + gd ; //// need reverse or not

        if(seci == 1 and si == 1){sign = (-1.0) * sign;} //// reverse sign if seci == 1
        for(Int da = 0; da < 12*12; da++)
        for(unsigned int vi = 0; vi < prop1.size();vi++)
        {
          Ty* src = (Ty*) qlat::get_data(buf[bind][vi]).data();
          Ty* res = (Ty*) qlat::get_data(resP[prj][vi]).data();
          const Long off_d = da*Nt*Nxyz + ti*Nxyz;
          cpy_data_threadC(&res[off_d], &src[off_d], Nxyz, 1, QFALSE, sign);
        }
        qacc_barrier(dummy);
      }
    }
  }

}



#ifdef QLAT_USE_ACC
template<typename Ty, Int bfac, Int Blocks>
__global__ void baryon_vectorEV_global(Ty** p1, Ty** p2, Ty** p3, Ty* resP,
  signed char** gPP, unsigned char** oPP, const Int* ivP,
  const Int nmass, const Int NTt, const Long Nxyz, const Int Ngv)
{
  //unsigned long gi =  blockIdx.x*blockDim.x + threadIdx.x;
  const unsigned long gi =  blockIdx.x;
  //const unsigned int tid = threadIdx.x;
  const unsigned int tid = threadIdx.y*blockDim.x+ threadIdx.x;
  /////const unsigned int Tid = blockDim.x;
  const Long Ntotal = nmass * NTt * Nxyz;
  const Long Nbfac  = Ntotal/bfac;

  const Int bfacC = bfac;
  const Int Nth   = Blocks/bfac;
  const unsigned int Each  = 16*bfacC;
  const unsigned int GROUP = (Blocks/bfac)*Each;

  __shared__ Ty P1[bfac*12*12];
  __shared__ Ty P2[bfac*12*12];
  __shared__ Ty P3[bfac*12*12];
  __shared__ Ty buf[bfacC*Blocks];
  __shared__ unsigned char pos[3*GROUP];
  __shared__ signed char g0[2*GROUP];

  //  if(gi*bfac < Ntotal)

  //__shared__ Ty buf[bfac*16+1];
  ////const Long offR0 = gi;
  Int bi0= 0;int dc = 0;
  Int ji = 0;int massi = 0;int ti = 0;
  ////const Long ixyz = (bi0*Nbfac + gi)%Nxyz;

  Int jobN = bfac*12*12;
  unsigned int off = tid;
  while(off < jobN){
    bi0= off/(12*12);
    dc = off%(12*12);
    ji    = (bi0*Nbfac + gi)/Nxyz;
    massi = ji/NTt;
    ti    = ji%NTt;
    Long ixyz = (bi0*Nbfac + gi)%Nxyz;
    //massi = (ji + bi)/NTt;
    //ti    = (ji + bi)%NTt;
    P1[dc*bfac + bi0] = p1[(massi*12*12 + dc)*NTt + ti][ixyz];
    P2[dc*bfac + bi0] = p2[(massi*12*12 + dc)*NTt + ti][ixyz];
    P3[dc*bfac + bi0] = p3[(massi*12*12 + dc)*NTt + ti][ixyz];
    off += Blocks;
  }
  __syncthreads();

  Int ini = 0;
  Int dv = 0;
  unsigned int MAX = 0;

  unsigned char* s0 = NULL;
    signed char* s1 = NULL;

  const Int bi =  threadIdx.y;
  const Int ai =  threadIdx.x;

  const Int ba = (threadIdx.y/bfacC)*bfacC + 0;
  const Int aa = (threadIdx.y%bfacC)*blockDim.x + threadIdx.x;

  for(Int iv=0;iv<Ngv;iv++)
  {
    for(Int bz=0;bz<bfacC;bz++){buf[bz*Blocks + tid] = 0;}
    MAX = ivP[iv];
    jobN = (MAX + GROUP -1 )/GROUP;
    ini = 0; dv = GROUP;
    for(Int ji=0;ji<jobN;ji++){
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
        const Ty* t1 = &P1[int(pos[off*3+0])*bfac + ba];
        const Ty* t2 = &P2[int(pos[off*3+1])*bfac + ba];
        const Ty* t3 = &P3[int(pos[off*3+2])*bfac + ba];
        Ty gtem = Ty(g0[off*2+0], g0[off*2+1]);
        if(bfacC == 1){
          buf[aa*bfac+ba] += (t1[0] * t2[0] * t3[0]*gtem);
        }else{
          Ty* b0 = &buf[aa*bfac+ba];
          for(Int z=0;z<bfacC;z++){b0[z] += (t1[z] * t2[z] * t3[z]*gtem);}
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

    if(ai == 0){
      //if(clear == 0){resP[iv*Ntotal + bi*Nbfac + gi] += (buf[bi] + buf[bfac+bi]);}
      //if(clear == 1){resP[iv*Ntotal + bi*Nbfac + gi]  = (buf[bi] + buf[bfac+bi]);}
      resP[iv*Ntotal + bi*Nbfac + gi] += (buf[bi] + buf[bfac+bi]); 
    }
    __syncthreads();

  }

}
#endif

////baryon on GPU, 
////USEGLOBAL then use global functions, else use qacc
template<typename Ty, Int bfac>
void baryon_vectorEV_kernel(Ty** p1, Ty** p2, Ty** p3, Ty* resP,
  signed char** gPP, unsigned char** oPP, const Int* ivP,
  const Int nmass, const Int NTt, const Long Nxyz, const Int Ngv)
{
  Long Ntotal  = nmass*NTt*Nxyz;
  if(Ntotal % bfac != 0){
    qmessage("Please correct your bfac! nmass %5d, NTt %5d, vol %ld, bfac %d. \n", int(nmass), int(NTt), long(Nxyz), bfac );
    Qassert(false);
  }
  Long Nbfac = Ntotal/bfac;

  #if USEGLOBAL==1
  const Int nt = 16;
  const Int Blocks = nt*bfac;
  dim3 dimBlock(    nt, bfac, 1);
  dim3 dimGrid(  Nbfac,  1, 1);
  baryon_vectorEV_global<Ty, bfac, Blocks><<<dimGrid, dimBlock>>>(p1, p2, p3, resP, gPP, oPP, ivP, nmass, NTt, Nxyz, Ngv);
  qacc_barrier(dummy);
  #else
  if((nmass*NTt) % bfac != 0){
    qmessage("Please correct your bfac! nmass %5d, NTt %5d, vol %ld, bfac %d. \n", int(nmass), int(NTt), long(Nxyz), bfac );
    Qassert(false);
  }
  qacc_for(gi, Nbfac ,
  {
    Ty buf[bfac+1];
    Ty P1[bfac*12*12+1];
    Ty P2[bfac*12*12+1];
    Ty P3[bfac*12*12+1];

    Long ixyz = gi%Nxyz;
    Int ji    = (gi/Nxyz)*bfac + 0;
    Int massi = ji/NTt;
    Int ti    = ji%NTt;
    const Long offR0 = (massi*NTt + ti)*Nxyz + ixyz;

    for(Int bi=0;bi<bfac;bi++)
    {
      massi = (ji+bi)/NTt;
      ti    = (ji+bi)%NTt;

      for(Int dc=0;dc<12*12;dc++){
        P1[dc*bfac + bi] = p1[(massi*12*12 + dc)*NTt + ti][ixyz];
        P2[dc*bfac + bi] = p2[(massi*12*12 + dc)*NTt + ti][ixyz];
        P3[dc*bfac + bi] = p3[(massi*12*12 + dc)*NTt + ti][ixyz];
      }
    }

    for(Int iv=0;iv<Ngv;iv++)
    {
      for(Int bi=0;bi<bfac;bi++){buf[bi] = 0;}

      for(Int off=0;off<ivP[iv];off++)
      {
        const Ty* t1 = &P1[(oPP[iv][off*3+0])*bfac];
        const Ty* t2 = &P2[(oPP[iv][off*3+1])*bfac];
        const Ty* t3 = &P3[(oPP[iv][off*3+2])*bfac];
        const Ty gtem = Ty(gPP[iv][off*2+0], gPP[iv][off*2+1]);
        for(Int bi=0;bi<bfac;bi++)
        {
          buf[bi] += (t1[bi] * t2[bi] * t3[bi] * gtem);
        }
      }

      Long offR = iv * Ntotal;
      Ty* r0 = &resP[offR + offR0];
      for(Int bi=0;bi<bfac; bi++){
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
void baryon_vectorEV(Ty** p1, Ty** p2, Ty** p3, Ty* resP, Int nmass,
  ga_M &A, ga_M &B, qlat::vector<Ty > &GV, 
  qlat::vector<Int > &mLV, fft_desc_basic &fd, Int clear=1)
{
  TIMER("Proton_vectorEV");
  Qassert(sizeof(Ty) == 16 or sizeof(Ty) == 8);
  Int NTt  = fd.Nv[3];
  Long Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  Int Ngv = GV.size()/16;
  Qassert(GV.size()  == 16*Ngv);
  Qassert(mLV.size() == 3*Ngv);

  if(clear == 1){zero_Ty(resP, Ngv*nmass*NTt*Nxyz , 1);}

  qlat::vector<Ty > epslV;epslV.resize(9);
  for(Int i=0;i<3;i++){epslV[i*3+i]=0;epslV[i*3 + (i+1)%3]=1;epslV[i*3 + (i+2)%3]=-1;}
  qlat::vector<Ty > gMap;
  qlat::vector<Int > IMap;
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
  const Ty* GVP = GV.data();
  const Int*  mLP     = mLV.data();

  /////contraction Kernel
  #if USEKERNEL==1
  ////Long Ntotal  = nmass*NTt*Nxyz;
  /////const Int Loff = 3*3*3*3*4*4*4*4;
  std::vector<std::vector<signed   char > > giEL;giEL.resize(Ngv);//giEL.resize(  Ngv*Loff);
  std::vector<std::vector<unsigned char > > oiL ;oiL.resize(Ngv );//oiL.resize(3*Ngv*Loff);
  Int bmL[3];
  Int nmL[3];
  Int count_flops  = 0;
  for(Int iv=0;iv<Ngv;iv++)
  {
    oiL[iv].resize(0);
    giEL[iv].resize(0);

    const Ty* G  = &GVP[iv*16];
    const Int*      mL = &mLP[iv*3];

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

  std::vector<qlat::vector_gpu<signed char > > giEG;giEG.resize(Ngv);
  for(Int iv=0;iv<Ngv;iv++){giEG[iv].copy_from(giEL[iv]);}
  qlat::vector<signed char* > gP = EigenM_to_pointers(giEG);
  signed char** gPP = gP.data();

  std::vector<qlat::vector_gpu<unsigned char   > > oiG ; oiG.resize(Ngv);
  for(Int iv=0;iv<Ngv;iv++){oiG[iv].copy_from(oiL[iv]);}
  qlat::vector<unsigned char* > oP = EigenM_to_pointers(oiG);
  unsigned char** oPP = oP.data();

  qlat::vector<Int > iv_size;iv_size.resize(Ngv);
  for(Int iv=0;iv<Ngv;iv++){
    iv_size[iv] = giEL[iv].size()/2;
    count_flops += (3 * 6 + 2) * iv_size[iv];
  }
  Int*  ivP = iv_size.data();
  ///int maxNv = iv_size[0];
  ///for(Int iv=0;iv<Ngv;iv++){if(iv_size[iv] > maxNv){maxNv = iv_size[iv];}}

  //{
  //Long total = 0;
  //for(Int iv=0;iv<Ngv;iv++){total += iv_size[iv];}
  //qmessage("==Ngv %d, total %d \n", int(Ngv), int(total));
  //}

  ////int mode = 0;mode = clear;
  {
  TIMER_FLOPS("baryon vectorEV kernel");
  timer.flops  += nmass*NTt *Nxyz * count_flops;
  const Int BFACG_DEFAULT = sizeof(Ty) == 16 ? BFACG_SHARED : BFACG_SHARED * 2;

  #ifdef QLAT_USE_ACC
  baryon_vectorEV_kernel<Ty, BFACG_DEFAULT>(p1, p2, p3, resP, gPP, oPP, ivP, nmass, NTt, Nxyz, Ngv);
  #else
  bool  get = false;
  const Long Ntot = nmass*NTt;
  const Int Bfac = BFACG_DEFAULT * 2;
  #define baryon_macros(ba) if(Ntot % ba == 0 and get == false){get = true; \
    baryon_vectorEV_kernel<Ty, ba>(p1, p2, p3, resP, gPP, oPP, ivP, nmass, NTt, Nxyz, Ngv);}
  baryon_macros(Bfac);
  baryon_macros(16);
  baryon_macros(10);
  baryon_macros(11);
  baryon_macros(8);
  baryon_macros(5);
  baryon_macros(3);
  baryon_macros(2);
  /////slowest mode ...
  baryon_macros(1);
  Qassert(get);
  #undef baryon_macros
  #endif
  }

  #endif

  #if USEKERNEL==0
  for(Int iv=0;iv<Ngv;iv++)
  {
    Long offR = iv*nmass*NTt * Nxyz;
    const Ty* G  = &(GVP[iv*16 + 0]);
    const Int*      mL = &(mLP[iv*3 + 0]);
    Int bmL[3];
    Int nmL[3];

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
        Ty Gtem =  G[m1*4+n1];
        double norm = qlat::qnorm(Gtem);
        if(norm < 1e-20)continue;

        Int m3 = gIA[m2];
        Int n3 = gIB[n2];
        Ty giE = epsl[a1*3 + a2]*epsl[b1*3 + b2]*gCA[m2]*gCB[n2]*G[m1*4+n1];
        nmL[0] = n1;nmL[1] = n2;nmL[2] = n3;
        bmL[0] = b1;bmL[1] = b2;bmL[2] = b3;
        Int nm1 = nmL[mL[0]];
        Int nm2 = nmL[mL[1]];
        Int nm3 = nmL[mL[2]];

        Int bm1 = bmL[mL[0]];
        Int bm2 = bmL[mL[1]];
        Int bm3 = bmL[mL[2]];

        #pragma omp parallel for
        for(Int ji=0;ji<nmass*NTt;ji++)
        {
          Int massi = ji/NTt;
          Int ti    = ji%NTt;

          Int o1 = massi*12*12 + (m1*3+a1)*12+(nm1*3+bm1);
          Int o2 = massi*12*12 + (m2*3+a2)*12+(nm2*3+bm2);
          Int o3 = massi*12*12 + (m3*3+a3)*12+(nm3*3+bm3);

          Ty* tp1 = p1[o1*NTt+ti];
          Ty* tp2 = p2[o2*NTt+ti];
          Ty* tp3 = p3[o3*NTt+ti];
          Ty* tr0 = &(resP[offR + (massi*NTt + ti)*Nxyz]);

          #if USEQACC==1
          qacc_forNB(i, Long(Nxyz),{ tr0[i] += (tp1[i]*tp2[i]*tp3[i] * giE); });
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
template <typename Ty>
void baryon_vectorEV(EigenTy& prop1, EigenTy& prop2, EigenTy& prop3,
  qlat::vector_gpu<Ty > &res, ga_M &A, ga_M &B, qlat::vector<Ty > &GV, qlat::vector<Int > &mLV,
  fft_desc_basic& fd, Int clear=1){

  if(prop1.size() == 0){res.resize(0); return ;}
  check_prop_size(prop1, fd);check_prop_size(prop2, fd);check_prop_size(prop3, fd);

  Int NTt  = fd.Nv[3];
  Long Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  ////check_prop_size(prop1);check_prop_size(prop2);check_prop_size(prop3);
  Int nmass = prop1.size();
  Qassert(prop1.size() == prop2.size());
  Qassert(prop1.size() == prop3.size());
  Int Ngv = GV.size()/16;
  const unsigned long resL = Ngv * nmass*NTt * Nxyz;
  if(clear == 1){if(res.size()!= resL){res.resize(resL); } }
  if(res.size() != resL){qmessage("Size of res wrong. \n");Qassert(false);}

  qlat::vector<Ty* > prop1P = EigenM_to_pointers(prop1, Nxyz);
  qlat::vector<Ty* > prop2P = EigenM_to_pointers(prop2, Nxyz);
  qlat::vector<Ty* > prop3P = EigenM_to_pointers(prop3, Nxyz);
  Ty** p1 = prop1P.data();
  Ty** p2 = prop2P.data();
  Ty** p3 = prop3P.data();
  Ty* resP = res.data();

  baryon_vectorEV(p1, p2, p3, resP, nmass, A, B, GV, mLV, fd, clear);
}



template <typename Ty>
void baryon_corrE(EigenTy& prop1, EigenTy& prop2, EigenTy& prop3,
  qlat::vector_gpu<Ty > &res,  ga_M &ga2,Int ind2,ga_M &ga1,Int ind1,
  fft_desc_basic& fd, Int clear=1,const Coordinate& mom = Coordinate())
{
  if(prop1.size() == 0){res.resize(0); return ;}
  //int NTt  = fd.Nv[3];
  ////LInt Nxyz = prop1[0].size();
  Int nmass = prop1.size();
  ////int nt = fd.nt;

  qlat::vector<Ty > resE;
  ini_resE(resE, nmass,fd);

  qlat::vector<Ty > G;G.resize(16);
  qlat::vector<Int > mL;mL.resize(3);

  clear_qv(G);G[ind2*4 + ind1] = +1.0;
  mL[0] = 0;mL[1] = 1;mL[2] = 2;
  baryon_vectorE(prop1,prop2,prop3, resE, ga2,ga1, G, mL, fd, 1);
  clear_qv(G);G[ind2*4 + ind1] = -1.0;
  mL[0] = 1;mL[1] = 0;mL[2] = 2;
  baryon_vectorE(prop1,prop2,prop3, resE, ga2,ga1, G, mL, fd, 0);

  vec_corrE(resE,res,fd,clear,mom);
}

template <typename Ty>
void Omega_corrE(EigenTy& prop1, EigenTy& prop2, EigenTy& prop3,
  qlat::vector_gpu<Ty > &res, ga_M &ga2,Int ind2,ga_M &ga1,Int ind1, Int clear=1,const Coordinate& mom = Coordinate())
{
  if(prop1.size() == 0){res.resize(0); return ;}
  const qlat::Geometry& geo = prop1[0].geo();
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);

  //int NTt  = fd.Nv[3];
  ///LInt Nxyz = prop1[0].size();
  Int nmass = prop1.size();
  ///int nt = fd.nt;

  qlat::vector<Ty > resE;
  ini_resE(resE,nmass,fd);

  qlat::vector<Ty > G;G.resize(16);
  qlat::vector<Int > mL;mL.resize(3);

  std::vector<Int > dia;dia.resize(6);
  std::vector<Int > sn ;sn.resize(6);
  dia[0] = 9012;sn[0] =  1;
  dia[1] = 9102;sn[1] = -1;
  dia[2] = 9021;sn[2] = -1;
  dia[3] = 9201;sn[3] =  1;
  dia[4] = 9210;sn[4] = -1;
  dia[5] = 9120;sn[5] =  1;

  for(Int di=0;di<6;di++)
  {
    clear_qv(G);G[ind2*4 + ind1] = sn[di];
    mL[0] = (dia[di]/100)%10;mL[1] =  (dia[di]%100)/10;mL[2] = dia[di]%10;
    baryon_vectorE(prop1,prop2,prop3, resE, ga2,ga1, G, mL, fd, 0);
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


  vec_corrE(resE,res,fd,clear,mom);
}

}

#endif

