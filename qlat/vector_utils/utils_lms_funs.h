// utils_lms_funs.h
// Gen Wang
// Jul. 2021

#ifndef UTILS_LMS_FUNS_H
#define UTILS_LMS_FUNS_H

#pragma once

#include "float_type.h"
#include "gammas.h"
#include "utils_momentum.h"
#include "fft_desc.h"
#include "utils_shift_vecs.h"
#include "utils_eigensys.h"
#include "utils_construction.h"
#include "utils_FFT_GPU.h"
#include "utils_grid_src.h"

////#define LOCALUSEACC 1

namespace qlat{

/////nmass * operator * nt * 2
void prop_to_corr(EigenM& Eprop, EigenV& Eres, fft_desc_basic& fd, int clear = 1)
{
  TIMERB("Get corr zero mom");
  ga_matrices_cps   ga_cps;
  std::vector<ga_M > gL;gL.resize(16);
  {int o=0;
  for(int i=0;i<6;i++){gL[o] = ga_cps.ga[0][i];o+=1;}
  for(int i=2;i<6;i++){gL[o] = ga_cps.ga[1][i];o+=1;}
  for(int i=3;i<6;i++){gL[o] = ga_cps.ga[2][i];o+=1;}
  for(int i=4;i<6;i++){gL[o] = ga_cps.ga[3][i];o+=1;}
  for(int i=5;i<6;i++){gL[o] = ga_cps.ga[4][i];o+=1;}}
  check_prop_size(Eprop);

  ///////===new contractions
  long NTt = fd.Nv[3];
  long Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  int nmass = Eprop.size()/(12*12*NTt);

  qlat::vector_acc<Complexq > resT0, resT1, resTa;
  resT0.resize(16 * nmass * NTt * Nxyz);
  resT1.resize(16 * nmass * NTt * Nxyz);
  resTa.resize(32 * nmass * NTt * Nxyz);

  qlat::vector_acc<Complexq > G ;G.resize( 2*16*16);
  qlat::vector_acc<int      > mL;mL.resize(2*16*3);
  ga_M &ga2 = ga_cps.ga[4][1];
  ga_M &ga1 = ga_cps.ga[4][1];

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
  
  qlat::vector_acc<Complexq* > propP = EigenM_to_pointers(Eprop);
  Complexq** p1 = propP.data();

  ////print0("Call!\n");

  meson_vectorEV( p1, p1, resT0.data(), nmass, gL, gL, fd, 1);
  baryon_vectorEV(p1, p1, p1, resTa.data(), nmass, ga2,ga1, G, mL, fd, 1);

  ///cpy_data_thread(resT1.data(), resTa.data(), resT1.size(), 1, true);
  //cpy_data_thread(resT1.data(), &(resTa.data()[resT1.size()]), resT1.size(), 1, true,  1.0);
  //cpy_data_thread(resT1.data(), resTa.data(), resT1.size(), 1, true);

  //cpy_data_thread(resTa.data(), &(resTa.data()[resT1.size()]), resT1.size(), 1, true,  1.0);
  //cpy_data_thread(&(resTa.data()[resT1.size()]), resT0.data(), resT0.size(), 1, true);

  cpy_data_thread(&(resTa.data()[resT1.size()]), resTa.data(), resT1.size(), 1, true,  1.0);
  cpy_data_thread( (resTa.data()              ), resT0.data(), resT0.size(), 1, true);


  vec_corrE(resTa, Eres, fd, clear, 505050);

  //EigenV Eres0, Eres1;
  //vec_corrE(resT0, Eres0, fd, 1, 505050);
  //vec_corrE(resT1, Eres1, fd, 1, 505050);

  //shift_result_t(Eres0, fd.nt, tini);
  //shift_result_t(Eres1, fd.nt, tini);

  //res.write_corr((Ftype*) &Eres0[0], 2*Eres0.size());
  //res.write_corr((Ftype*) &Eres1[0], 2*Eres1.size());
  ///////===new contractions

  ///////===old contractions
  //EigenV Eres;
  //for(LInt i=0;i<gL.size();i++)
  //{
  //  meson_corrE(Eprop, Eprop, gL[i], gL[i], Eres, fd);
  //  res.write_corr((Ftype*) &Eres[0],2*Eres.size());
  //}
  //for(int ind2=0;ind2<4;ind2++)
  //for(int ind1=0;ind1<4;ind1++)
  //{
  //  baryon_corrE(Eprop, Eprop, Eprop, ga_cps.ga[4][1], ind2, ga_cps.ga[4][1], ind1, Eres, fd, 1, 505050);
  //  res.write_corr((Ftype*) &Eres[0],2*Eres.size());
  //}
  ///////===old contractions

  //int nt = fd.nt;
  //for(int op=0;op<16;op++)
  //{
  //  for(int mi=0;mi<nmass;mi++)
  //  {
  //    long off = (op*nmass + mi)*2*nt;
  //    res.write_corr((Ftype*) &Eres0[off], 2*nt);
  //  }
  //  for(int mi=0;mi<nmass;mi++)
  //  {
  //    long off = (op*nmass + mi)*2*nt;
  //    res.write_corr((Ftype*) &Eres1[off], 2*nt);
  //  }
  //}


}

/////propH , pi --> 12*12 --> Nt*Nxyz
template<typename Ty>
void point_corr(std::vector<qnoiT >& src, std::vector<qpropT >& propH,
    std::vector<double >& massL, eigen_ov& ei, fft_desc_basic& fd, corr_dat& res, int mode_sm = 0, int shift_t = 1, int lms = 0, int SRC_PROP_WITH_LOW = 0)
{
  if(propH.size() != src.size()*massL.size()){abort_r("prop size, mass size not match!\n");}
  int Ns = src.size();
  long Size_prop = fd.get_prop_size();

  int GPU = 1;bool rotate = false;
  int nmass = massL.size();

  std::vector<Coordinate > pos;pos.resize(Ns);
  std::vector<qlat::vector_acc<int > > off_L;off_L.resize(Ns);
  std::vector<qlat::vector_acc<int > > Ngrid;Ngrid.resize(Ns);
  for(int si=0;si<Ns;si++){check_noise_pos(src[si], pos[si], off_L[si]);}
  for(int si=0;si<Ns;si++){grid_list_pos(off_L[si], Ngrid[si]);}

  int tini = pos[0][3];
  if(shift_t == 1)for(int si=0;si<Ns;si++){if(tini != pos[si][3]){abort_r("Shift t for each src not coded currently!\n");}}
  for(int si=0;si<Ns;si++){if(Ngrid[0].size() != Ngrid[si].size()){
    abort_r("grid source not consistent, split the src production! \n");}}

  /////print0("Ns %d, Ngrid %d \n", Ns, int(Ngrid[0].size()) );

  /////need modification of positions
  ///for(int iv=0;iv<src.size();iv++){check_noise_pos(src[iv], pos[iv], off_L);}
  EigenV EresH;EresH.resize(32 * nmass * fd.nt);clear_qv(EresH );
  EigenV EresL;EresL.resize(32 * nmass * fd.nt);clear_qv(EresL );
  EigenV EresA;EresA.resize(32 * nmass * fd.nt);clear_qv(EresA );

  /////do corr
  EigenM Eprop;

  qlat::vector_gpu<Complexq > stmp, low_prop, high_prop;
  /////copy src vectors
  stmp.resize(Ns*Size_prop);

  copy_eigen_src_to_FeildM(high_prop, propH, ei.b_size, fd, 1, GPU, rotate);
  copy_eigen_prop_to_EigenM(high_prop.data(), Eprop, ei.b_size, nmass, fd, 0, GPU);
  low_prop.resize(high_prop.size());

  /////reduce the high prop
  if(SRC_PROP_WITH_LOW == 1){
    stmp.set_zero();
    low_prop.set_zero();
    std::vector<qlat::FieldM<Ty , 12*12> > src_prop;
    FieldM_src_to_FieldM_prop(src, src_prop, GPU);
    copy_eigen_src_to_FeildM(stmp, src_prop, ei.b_size, fd, 1, GPU, rotate);
    prop_L_device(ei, stmp.data(), low_prop.data(), 12*Ns, massL, mode_sm);
    high_prop -= low_prop;
  }
  /////reduce the high prop

  prop_to_corr(Eprop, EresH, fd, 0);  

  /////set memory for low_prop
  int Nlms = 1;
  if(lms == -1){Nlms = Ngrid[0].size();}
  if(lms == 0 ){Nlms = 1;}
  if(lms >  0 ){Nlms = lms;}

  for(int gi=0;gi<Nlms;gi++){
    stmp.set_zero();
    low_prop.set_zero();

    if(lms == 0)
    {
      std::vector<qlat::FieldM<Ty , 12*12> > src_prop;
      FieldM_src_to_FieldM_prop(src, src_prop, GPU);
      copy_eigen_src_to_FeildM(stmp, src_prop, ei.b_size, fd, 1, GPU, rotate);
    }

    if(lms != 0)
    for(int si=0;si<Ns;si++){
      Coordinate cur_pos = get_grid_off(Ngrid[si][gi], off_L[si], pos[si], fd);
      ////print0("position %d %d %d %d !\n", cur_pos[0], cur_pos[1], cur_pos[2], cur_pos[3]);
      write_grid_point_to_src(&stmp.data()[si*Size_prop], src[si], cur_pos, ei.b_size, fd);

      ////std::vector<qlat::FieldM<Ty , 12*12> > src_prop;
      ////copy_eigen_src_to_FeildM(stmp, src_prop, ei.b_size, fd, 0, GPU, true);
      ////for(LInj isp=0;isp<fd.noden;isp++){
      ////  for(long d=0;d<12*12;d++){
      ////    Ty v = src_prop[0].get_elems(isp)[d];
      ////    if(qlat::qnorm(v) > 1e-15) 
      ////    {
      ////      printf("isp %d, d %d, v %.3e %.3e. \n", int(isp), int(d), v.real(), v.imag());
      ////    }
      ////  }
      ////}
    }

    /////get low mode prop
    prop_L_device(ei, stmp.data(), low_prop.data(), 12*Ns, massL, mode_sm);

    //////zeroE(Eprop, 1);
    copy_eigen_prop_to_EigenM(low_prop.data(),  Eprop, ei.b_size, nmass, fd, 0, GPU);
    prop_to_corr(Eprop, EresL, fd, 0); 

    low_prop += high_prop;
    copy_eigen_prop_to_EigenM(low_prop.data(),  Eprop, ei.b_size, nmass, fd, 0, GPU);
    prop_to_corr(Eprop, EresA, fd, 0);  

  }

  if(shift_t == 1){
    shift_result_t(EresL, fd.nt, tini);
    shift_result_t(EresH, fd.nt, tini);
    shift_result_t(EresA, fd.nt, tini);
  }

  res.write_corr((Ftype*) &EresL[0], 2*EresL.size());
  res.write_corr((Ftype*) &EresH[0], 2*EresH.size());
  res.write_corr((Ftype*) &EresA[0], 2*EresA.size());

}

/////all to all prop naive do, test for all-to-all low eigen
void test_all_prop_corr(std::vector<double >& massL, eigen_ov& ei, fft_desc_basic& fd, corr_dat& res, int mode_sm = 0)
{
  /////if(propH.size() != src.size()*massL.size()){abort_r("prop size, mass size not match!\n");}
  Geometry geo;fd.get_geo(geo );

  int GPU = 1;
  int nmass = massL.size();

  std::vector<qnoi > src; src.resize(1);src[0].init(geo);
  EigenV Eres;Eres.resize(nmass*32*fd.nt);
  qlat::set_zero(Eres);

  ///std::vector<int > pos;pos.resize(src.size());qlat::vector_acc<int > off_L;
  ///check_noise_pos(src[0], pos[0], off_L);
  ///int tini = pos[0]%1000;
  qlat::vector_gpu<Complexq > stmp, low_prop;
  stmp.resize(fd.get_prop_size());
  low_prop.resize(nmass*fd.get_prop_size());
  /////do corr
  EigenM Eprop;

  int tini = 0;
  for(int zi=0;zi<fd.nz;zi++)
  for(int yi=0;yi<fd.ny;yi++)
  for(int xi=0;xi<fd.nx;xi++)
  {
    print0("===pos %d %d %d \n", zi, yi, xi);
    stmp.set_zero();
    low_prop.set_zero();

    LInt index0 = fd.index_g_from_local(0);
    LInt indexL = fd.index_g_from_g_coordinate(tini, zi, yi, xi);
    if(indexL >= index0 and indexL < index0 + fd.noden)
    {
      for(int di=0;di<12;di++){
        stmp[(di*12 + di)*fd.noden + indexL%fd.noden] = 1.0;
      }
    }

    /////get low mode prop
    prop_L_device(ei, stmp.data(), low_prop.data(), 12, massL, mode_sm);

    copy_eigen_prop_to_EigenM(low_prop.data(),  Eprop, ei.b_size, nmass, fd, 0, GPU);
    prop_to_corr(Eprop, Eres, fd,  0); 
    shift_result_t(Eres, fd.nt, tini);

    ////copy_eigen_prop_to_EigenM(high_prop.data(), Eprop, ei.b_size, nmass, fd, 0, GPU);
    ////prop_to_corr(Eprop, Eres, fd, tini);  res.write_corr((Ftype*) &Eres[0], 2*Eres.size());

    ////low_prop += high_prop;

    ////copy_eigen_prop_to_EigenM(low_prop.data(),  Eprop, ei.b_size, nmass, fd, 0, GPU);
    ////prop_to_corr(Eprop, Eres, fd, tini);  

  }
  for(int i=0;i<Eres.size();i++){Eres[i] = Eres[i]/(Ftype(fd.nx*fd.ny*fd.nz*1.0));}
  res.write_corr((Ftype*) &Eres[0], 2*Eres.size());
}


//struct lms_schedule{
//  fft_desc_basic* fdp;
//  eigen_ov*       eip;
//
//  fft_schedule   fft;
//  shift_vec      svec;
//
//  bool GPU;
//  int max_points;
//
//  lms_schedule(eigen_ov& ei,fft_desc_basic& fd, bool GPU_set=true):fft(fd, GPU_set),svec(fd, GPU_set)
//  {
//    fdp = &fd;
//    eip = &ei;
//    #ifndef QLAT_USE_ACC
//    GPU = false;
//    #else
//    GPU = GPU_set;
//    #endif
//
//    max_points = 1;
//
//  }
//
//  void set_paras(int maxP, int )
//  {
//    max_points = maxP;
//
//  }
//
//  void get_lms_corr(Complexq *src, Complexq *props, int Ns,std::vector<double> &mass);
//
//  ~lms_schedule(){
//    fdp = NULL;
//    eip = NULL;
//  }
//
//};
//
//void lms_schedule::get_lms_corr(Complexq *src, Complexq *props, int Ns,std::vector<double> &mass)
//{
//  prop_L_device(*eip, src, props, Ns, mass);
//}


}


#endif

