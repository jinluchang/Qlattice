// utils_construction.h
// Gen Wang
// Jul. 2021

#ifndef UTILS_CONSTRUCTION_H
#define UTILS_CONSTRUCTION_H

#pragma once

#include "utils_float_type.h"
#include "utils_props_type.h"
#include "utils_gammas.h"
#include "utils_fft_desc.h"
#include "utils_reduce_vec.h"
#include "utils_grid_src.h"
#include "utils_io_vec.h"

#ifdef QLAT_USE_ACC
#define USEKERNEL 1
#define USEGLOBAL 1
#define USEQACC   1
#else
////#define USEKERNEL 0
////#define USEGLOBAL 0
////#define USEQACC   0
#define USEKERNEL 1
#define USEGLOBAL 0
#define USEQACC   0
#endif

#define EigenTy std::vector<qlat::vector_gpu<Ty > >

///#define EigenMTa std::vector<qlat::vector_acc<Ta > >
//#define EigenVTa qlat::vector_acc<Ta >
#define EAy   Eigen::Map<Eigen::Array<Ty ,Eigen::Dynamic,1 > >
//#define EAa   Eigen::Map<Eigen::Array<Ta ,Eigen::Dynamic,1 > >

#include "utils_corr_prop.h"
#include "utils_corr_meson.h"
#include "utils_corr_baryon.h"
#include "utils_corr_seq.h"

namespace qlat{


template<typename Ty>
void prop_to_vec(qlat::vector_acc<Ty* >& propP, qlat::vector_gpu<Ty >& resTa, fft_desc_basic& fd)
{
  TIMERB("Get corr vec");
  //check_prop_size(Eprop, fd);

  ga_matrices_cps   ga_cps;
  ////ga_matrices_PS   ga_cps;
  std::vector<ga_M > gL;gL.resize(16);
  {int o=0;
  for(int i=0;i<6;i++){gL[o] = ga_cps.ga[0][i];o+=1;}
  for(int i=2;i<6;i++){gL[o] = ga_cps.ga[1][i];o+=1;}
  for(int i=3;i<6;i++){gL[o] = ga_cps.ga[2][i];o+=1;}
  for(int i=4;i<6;i++){gL[o] = ga_cps.ga[3][i];o+=1;}
  for(int i=5;i<6;i++){gL[o] = ga_cps.ga[4][i];o+=1;}}

  ///////===new contractions
  const Long NTt = fd.Nv[3];
  const Long Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  //int nmass = Eprop.size();
  Qassert(propP.size() % ( 12 * 12 * fd.Nt) == 0);
  const int nmass = propP.size() / ( 12 * 12 * fd.Nt);

  //qlat::vector_gpu<Ty > resT0, resT1;////, resTa;
  //resT0.resize(16 * nmass * NTt * Nxyz);
  //resT1.resize(16 * nmass * NTt * Nxyz);
  resTa.resize(32 * nmass * NTt * Nxyz);//qlat::set_zero(resTa);////resTa.set_zero();

  ////gamma matrix follow current prec
  qlat::vector_acc<Ty > G ;G.resize( 2*16*16);
  qlat::vector_acc<int      > mL;mL.resize(2*16*3);
  ga_M &ga2 = ga_cps.ga[1][3];
  ga_M &ga1 = ga_cps.ga[1][3];

  clear_qv(G );clear_qv(mL);
  for(int ind2=0;ind2<4;ind2++)
  for(int ind1=0;ind1<4;ind1++)
  {
    int ioff = ind2*4 + ind1;
    G[ioff*16 + ind2*4 + ind1] = +1.0;
    mL[ioff*3 + 0] = 0;
    mL[ioff*3 + 1] = 1;
    mL[ioff*3 + 2] = 2;

    G[1*16*16 + ioff*16 + ind2*4 + ind1] = -1.0;
    mL[1*16*3 + ioff*3 + 0] = 1;
    mL[1*16*3 + ioff*3 + 1] = 0;
    mL[1*16*3 + ioff*3 + 2] = 2;
  }

  ////qlat::vector_acc<Ty* > propP = EigenM_to_pointers(Eprop, Nxyz);
  ////cps to PS

  Ty** p1 = propP.data();

  ////print0("Call!\n");

  Ty* ra = resTa.data();
  Ty* rb = &(resTa.data()[resTa.size()/2]);

  ////p1 must be pointters
  baryon_vectorEV(p1, p1, p1, ra, nmass, ga2,ga1, G, mL, fd, 1);
  /////add baryon two contractions
  cpy_data_thread(rb, ra, resTa.size()/2, 1, QTRUE,  1.0);
  
  meson_vectorEV( p1, p1, ra, nmass, gL, gL, fd, 1);

  //cpy_data_thread( (resTa.data()              ), resT0.data(), resT0.size(), 1, true);
  //vec_corrE(resTa, Eres, fd, clear);
}

template<typename Ty>
void prop_to_vec(std::vector<qlat::vector_gpu<Ty > >& Eprop, qlat::vector_gpu<Ty >& resTa, fft_desc_basic& fd)
{
  check_prop_size(Eprop, fd);

  const Long Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  ////Eprop, nmass --> 12 * 12 * Nvol
  qlat::vector_acc<Ty* > propP = EigenM_to_pointers(Eprop, Nxyz);
  prop_to_vec(propP, resTa, fd);
}

template<typename Ty>
void prop_to_vec(std::vector<qpropT >& Eprop, qlat::vector_gpu<Ty >& resTa, fft_desc_basic& fd)
{
  const Long Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  ////Eprop, nmass --> 12 * 12 * Nvol
  qlat::vector_acc<Ty* > propP = EigenM_to_pointers(Eprop, Nxyz);
  prop_to_vec(propP, resTa, fd);
}

template<typename Td>
void prop_to_vec(std::vector<qlat::Propagator4dT<Td >* >& prop4d, qlat::vector_gpu<qlat::ComplexT<Td> >& resTa, fft_desc_basic& fd)
{
  std::vector<qlat::FieldM<qlat::ComplexT<Td >, 12*12> > Eprop;
  const int nmass = prop4d.size();
  Eprop.resize(nmass);
  for(int mi=0;mi<nmass;mi++){
    prop4d_to_qprop(Eprop[mi], *prop4d[mi]);
  }
  prop_to_vec(Eprop, resTa, fd);
}

template<typename Td>
void prop_to_vec(std::vector<qlat::Propagator4dT<Td > >& prop4d, qlat::vector_gpu<qlat::ComplexT<Td> >& resTa, fft_desc_basic& fd)
{
  std::vector<qlat::Propagator4dT<Td >* > prop4dP;
  const int nmass = prop4d.size();
  prop4dP.resize(nmass);
  for(int mi=0;mi<nmass;mi++){
    prop4dP[mi] = &prop4d[mi];
  }
  prop_to_vec(prop4dP, resTa, fd);
}

template<typename Ty>
void prop_to_corr_mom0(std::vector<qlat::vector_gpu<Ty > >& Eprop, qlat::vector_acc<Ty >& Eres, 
  fft_desc_basic& fd, qlat::vector_gpu<Ty >& resTa, int clear = 1)
{
  prop_to_vec(Eprop, resTa, fd);  
  vec_corrE(resTa, Eres, fd, 0);
}

template<typename Td>
void prop_corrE(Propagator4dT<Td > &p1,
  qlat::vector_acc<qlat::ComplexT<Td > >& Eres, const Coordinate& mom = Coordinate(), const int tini = 0)
{
  const Geometry& geo = p1.geo();
  std::vector<Propagator4dT<Td >* > pL;pL.resize(1);
  pL[0] = &p1;

  qlat::vector_gpu<qlat::ComplexT<Td > > resTa;
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);
  prop_to_vec(pL, resTa, fd);
  vec_corrE(resTa, Eres, fd, 1, mom);
  if(tini != 0){
    shift_result_t(Eres, fd.nt, tini);
  }
}

template<typename Td>
void prop_corrE(Propagator4dT<Td > &p1,
  const std::string& filename, const Coordinate& mom = Coordinate(), const int tini = 0,
  const std::string& info = std::string("NONE"), const int shift_end = 1)
{
  const Geometry& geo = p1.geo();
  std::vector<Propagator4dT<Td >* > pL;pL.resize(1);
  pL[0] = &p1;

  qlat::vector_gpu<qlat::ComplexT<Td > > resTa;
  qlat::vector_acc<qlat::ComplexT<Td > > Eres;
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);
  prop_to_vec(pL, resTa, fd);
  vec_corrE(resTa, Eres, fd, 1, mom);
  if(tini != 0){
    shift_result_t(Eres, fd.nt, tini);
  }

  const size_t sizen = get_file_size_MPI(filename, true);
  corr_dat<double > corr(std::string(""));
  if(sizen > 0){
    corr.read_dat(filename, 1);
    if(shift_end == 1){
      corr.shift_end();
    }
  }
  else{
    //char key_T[1000], dimN[1000];
    //sprintf(key_T, "%d  %d %d  %d", 1, 32, fd.nt, 2); 
    //sprintf(dimN , "nsrc nop nt complex");
    std::string ktem = ssprintf("%d  %d %d  %d", 1, 32, fd.nt, 2);
    std::string dtem = ssprintf("nsrc nop nt complex");
    corr.create_dat(ktem, dtem);
  }
  if(info != std::string("NONE") and info.size() != 0){
    corr.INFOA.push_back(info);
  }

  corr.write_corr(Eres.data(), Eres.size());
  corr.write_dat(filename);
}

}

#undef  EigenTy
///#undef  EigenMTa
//#undef  EigenVTa
#undef  EAy
//#undef  EAa


#endif

