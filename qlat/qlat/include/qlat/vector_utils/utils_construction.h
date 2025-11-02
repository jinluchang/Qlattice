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

///#define EigenMTa std::vector<qlat::vector<Ta > >
//#define EigenVTa qlat::vector<Ta >
#define EAy   Eigen::Map<Eigen::Array<Ty ,Eigen::Dynamic,1 > >
//#define EAa   Eigen::Map<Eigen::Array<Ta ,Eigen::Dynamic,1 > >

#include "utils_corr_prop.h"
#include "utils_corr_meson.h"
#include "utils_corr_baryon.h"
#include "utils_corr_seq.h"

namespace qlat{


template<typename Ty>
void prop_to_vec(qlat::vector<Ty* >& propP, qlat::vector_gpu<Ty >& resTa, fft_desc_basic& fd)
{
  TIMERB("Get corr vec");
  //check_prop_size(Eprop, fd);

  ga_matrices_cps   ga_cps;
  ////ga_matrices_PS   ga_cps;
  std::vector<ga_M > gL;gL.resize(16);
  {int o=0;
  for(Int i=0;i<6;i++){gL[o] = ga_cps.ga[0][i];o+=1;}
  for(Int i=2;i<6;i++){gL[o] = ga_cps.ga[1][i];o+=1;}
  for(Int i=3;i<6;i++){gL[o] = ga_cps.ga[2][i];o+=1;}
  for(Int i=4;i<6;i++){gL[o] = ga_cps.ga[3][i];o+=1;}
  for(Int i=5;i<6;i++){gL[o] = ga_cps.ga[4][i];o+=1;}}

  ///////===new contractions
  const Long NTt = fd.Nv[3];
  const Long Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  //int nmass = Eprop.size();
  Qassert(propP.size() % ( 12 * 12 * fd.Nt) == 0);
  const Int nmass = propP.size() / ( 12 * 12 * fd.Nt);

  //qlat::vector_gpu<Ty > resT0, resT1;////, resTa;
  //resT0.resize(16 * nmass * NTt * Nxyz);
  //resT1.resize(16 * nmass * NTt * Nxyz);
  resTa.resize(32 * nmass * NTt * Nxyz);//qlat::set_zero(resTa);////resTa.set_zero();

  ////gamma matrix follow current prec
  qlat::vector<Ty > G ;G.resize( 2*16*16);
  qlat::vector<Int      > mL;mL.resize(2*16*3);
  ga_M &ga2 = ga_cps.ga[1][3];
  ga_M &ga1 = ga_cps.ga[1][3];

  clear_qv(G );clear_qv(mL);
  for(Int ind2=0;ind2<4;ind2++)
  for(Int ind1=0;ind1<4;ind1++)
  {
    Int ioff = ind2*4 + ind1;
    G[ioff*16 + ioff] = +1.0;
    mL[ioff*3 + 0] = 0;
    mL[ioff*3 + 1] = 1;
    mL[ioff*3 + 2] = 2;

    G[1*16*16 + ioff*16 + ioff] = -1.0;
    mL[1*16*3 + ioff*3 + 0] = 1;
    mL[1*16*3 + ioff*3 + 1] = 0;
    mL[1*16*3 + ioff*3 + 2] = 2;
  }

  ////qlat::vector<Ty* > propP = EigenM_to_pointers(Eprop, Nxyz);
  ////cps to PS

  Ty** p1 = propP.data();

  ////qmessage("Call!\n");

  Ty* ra = resTa.data();
  Ty* rb = &(resTa.data()[resTa.size()/2]);

  ////p1 must be pointters
  baryon_vectorEV(p1, p1, p1, ra, nmass, ga2,ga1, G, mL, fd, 1);
  /////add baryon two contractions
  cpy_data_thread(rb, ra, resTa.size()/2, 1, QTRUE,  1.0);
  
  qlat::vector<Ty* > resvP;resvP.resize(16 * nmass);
  for(Int iv = 0;iv<16*nmass;iv++){
    resvP[iv] = &ra[iv * NTt * Nxyz];
  }
  meson_vectorEV( p1, p1, resvP.data(), nmass, gL, gL, fd, 1);

  //meson_vectorEV( p1, p1, ra, nmass, gL, gL, fd, 1);

  //cpy_data_thread( (resTa.data()              ), resT0.data(), resT0.size(), 1, true);
  //vec_corrE(resTa, Eres, fd, clear);
}

template<typename Ty>
void prop_to_vec(std::vector<qlat::vector_gpu<Ty > >& Eprop, qlat::vector_gpu<Ty >& resTa, fft_desc_basic& fd)
{
  check_prop_size(Eprop, fd);

  const Long Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  ////Eprop, nmass --> 12 * 12 * Nvol
  qlat::vector<Ty* > propP = EigenM_to_pointers(Eprop, Nxyz);
  prop_to_vec(propP, resTa, fd);
}

template<typename Ty>
void prop_to_vec(FieldG<Ty >& prop, qlat::vector_gpu<Ty >& resTa)
{
  Qassert(prop.initialized and prop.multiplicity == 12 * 12 and prop.mem_order == QLAT_OUTTER);
  const Geometry& geo = prop.geo();
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);
  const Long Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  std::vector<FieldG<Ty>> buf;buf.resize(1);
  buf[0].set_pointer(prop);
  ////Eprop, nmass --> 12 * 12 * Nvol
  qlat::vector<Ty* > propP = FieldG_to_pointers(buf, Nxyz);
  prop_to_vec(propP, resTa, fd);
}

template<typename Ty>
void prop_to_vec(std::vector<qpropT >& Eprop, qlat::vector_gpu<Ty >& resTa, fft_desc_basic& fd)
{
  const Long Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  ////Eprop, nmass --> 12 * 12 * Nvol
  //qlat::vector<Ty* > propP = EigenM_to_pointers(Eprop, Nxyz);
  qlat::vector<Ty* > propP = FieldM_to_Tpointers(Eprop, Nxyz);
  prop_to_vec(propP, resTa, fd);
}

template<typename Td>
void prop_to_vec(std::vector<qlat::Propagator4dT<Td >* >& prop4d, qlat::vector_gpu<qlat::ComplexT<Td> >& resTa, fft_desc_basic& fd)
{
  std::vector<qlat::FieldM<qlat::ComplexT<Td >, 12*12> > Eprop;
  const Int nmass = prop4d.size();
  Eprop.resize(nmass);
  for(Int mi=0;mi<nmass;mi++){
    prop4d_to_qprop(Eprop[mi], *prop4d[mi]);
  }
  prop_to_vec(Eprop, resTa, fd);
}

template<typename Td>
void prop_to_vec(std::vector<qlat::Propagator4dT<Td > >& prop4d, qlat::vector_gpu<qlat::ComplexT<Td> >& resTa, fft_desc_basic& fd)
{
  std::vector<qlat::Propagator4dT<Td >* > prop4dP;
  const Int nmass = prop4d.size();
  prop4dP.resize(nmass);
  for(Int mi=0;mi<nmass;mi++){
    prop4dP[mi] = &prop4d[mi];
  }
  prop_to_vec(prop4dP, resTa, fd);
}

template<typename Ty>
void prop_to_corr_mom0(std::vector<qlat::vector_gpu<Ty > >& Eprop, qlat::vector<Ty >& Eres, 
  fft_desc_basic& fd, qlat::vector_gpu<Ty >& resTa, Int clear = 1)
{
  prop_to_vec(Eprop, resTa, fd);  
  vec_corrE(resTa, Eres, fd, 0, clear);
}

template<typename Td>
void prop_corrE(Propagator4dT<Td > &p1,
  qlat::vector<qlat::ComplexT<Td > >& Eres, const Coordinate& mom = Coordinate(), const Int tini = 0)
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
  const std::string& filename, const Coordinate& mom = Coordinate(), const Int tini = 0,
  const std::string& info = std::string("NONE"), const Int shift_end = 1)
{
  const Geometry& geo = p1.geo();
  std::vector<Propagator4dT<Td >* > pL;pL.resize(1);
  pL[0] = &p1;

  qlat::vector_gpu<qlat::ComplexT<Td > > resTa;
  qlat::vector<qlat::ComplexT<Td > > Eres;
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

