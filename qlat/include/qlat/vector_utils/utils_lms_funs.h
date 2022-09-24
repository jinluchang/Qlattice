// utils_lms_funs.h
// Gen Wang
// Jul. 2021

#ifndef UTILS_LMS_FUNS_H
#define UTILS_LMS_FUNS_H

#pragma once

#include "utils_float_type.h"
#include "utils_gammas.h"
#include "utils_momentum.h"
#include "utils_fft_desc.h"
#include "utils_shift_vecs.h"
#include "utils_eigensys.h"
#include "utils_construction.h"
#include "utils_FFT_GPU.h"
#include "utils_grid_src.h"

////#define LOCALUSEACC 1

namespace qlat{

struct lms_para{

  Coordinate ini_pos;
  Coordinate off_L;

  int SRC_PROP_WITH_LOW ;
  int lms;
  int combineT;

  int src_smear_iter;
  double src_smear_kappa;
  int sink_smear_iter;
  double sink_smear_kappa;

  int mode_eig_sm;
  int ionum;

  int do_all_low ;

  int mom_cut;
  int save_zero_corr;

  std::string name_mom_vecs;
  std::string name_zero_vecs;
  std::string name_zero;

  std::string INFO;
  std::vector<std::string > INFOA;

  /////initial with src noise

  void init(){
    for(int i=0;i<4;i++){ini_pos[i] = 0;off_L[i] = 1;}

    SRC_PROP_WITH_LOW = 0;
    lms = 0;
    combineT = 0;

    ///0  pt to pt, 1 pt to sm, 2 sm to pt, 3 sm to sm
    mode_eig_sm      = 0;
    src_smear_iter   = 0;
    sink_smear_iter  = 0;
    src_smear_kappa  = 0.0;
    sink_smear_kappa = 0.0;
    do_all_low       = 1;

    ionum            = 16;

    mom_cut          = 4;

    name_mom_vecs    = std::string("NONE");
    name_zero_vecs   = std::string("NONE");
    name_zero        = std::string("NONE");
    save_zero_corr   = 1;

    INFO = std::string("NONE");
    INFOA.resize(0);
  }

  void print(){
    print0("===Source Info %s \n", INFO.c_str());
    print0("  prop with low %3d, lms %3d, combineT %3d \n", SRC_PROP_WITH_LOW, lms, combineT);
    print0("  init %5d %5d %5d %5d, grid %5d %5d %5d %5d \n", 
      ini_pos[0], ini_pos[1], ini_pos[2], ini_pos[3],
      off_L[0], off_L[1], off_L[2], off_L[3] );

    int save_vecs = 0; if(name_mom_vecs != std::string("NONE")){save_vecs = 1;}
    print0("  mode sm %1d, src %5d %.3f, sink %5d %.3f , save_low %2d, mom_cut %5d , saveV %1d \n", 
      mode_eig_sm, src_smear_iter, src_smear_kappa, sink_smear_iter, sink_smear_kappa, do_all_low, mom_cut, save_vecs);
    //print0("  mom vecs %s \n" ,name_mom_vecs);
    //print0("  zero vecs %s \n",name_zero_vecs);
  }

};


template<typename Ty>
void pick_mom_data(qlat::vector_gpu<Ty >& res, qlat::vector_gpu<Ty >& src,
  const int nvec, qlat::vector_acc<long >& mapA, qlat::vector_acc<long >& mapB, const Geometry& geo)
{
  TIMER("save FFT");
  ////false to match to gwu code convention
  ////exp(+ i p y) need src exp(-i p x)
  fft_fieldM(src.data(), nvec, 1, geo, false);
  res.set_zero();

  const long Nvol = src.size()/nvec;
  const long Mvol = res.size()/nvec;
  qassert(mapA.size() <= Mvol);

  long* A = mapA.data();
  long* B = mapB.data();
  qacc_for(isp, mapA.size() ,{
    long i0 = A[isp];
    long i1 = B[isp];
    for(int iv=0;iv<nvec;iv++){res[iv*Mvol + i1] = src[iv*Nvol + i0];}
  });

  sum_all_size(res.data(), res.size(), 1);
}

template<typename Ty>
void prop_to_vec(std::vector<qlat::vector_acc<Ty > >& Eprop, qlat::vector_gpu<Ty >& resTa, fft_desc_basic& fd)
{
  TIMERB("Get corr vec");
  check_prop_size(Eprop);

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
  long NTt = fd.Nv[3];
  long Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  int nmass = Eprop.size()/(12*12*NTt);

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

  qlat::vector_acc<Ty* > propP = EigenM_to_pointers(Eprop);

  ////cps to PS
  //Geometry geo;fd.get_geo(geo);
  //std::vector< Propagator4dT<Ty > > prop;prop.resize(nmass);
  //for(unsigned int pi=0;pi<prop.size();pi++)prop[pi].init(geo);
  //std::vector<qlat::vector_acc<Ty > > Eprop_tmp;

  Ty** p1 = propP.data();

  ////print0("Call!\n");

  Ty* ra = resTa.data();
  Ty* rb = &(resTa.data()[resTa.size()/2]);

  baryon_vectorEV(p1, p1, p1, ra, nmass, ga2,ga1, G, mL, fd, 1);
  /////add baryon two contractions
  cpy_data_thread(rb, ra, resTa.size()/2, 1, true,  1.0);
  
  meson_vectorEV( p1, p1, ra, nmass, gL, gL, fd, 1);

  //cpy_data_thread( (resTa.data()              ), resT0.data(), resT0.size(), 1, true);
  //vec_corrE(resTa, Eres, fd, clear);
}

template<typename Ty>
void prop_to_corr_mom0(std::vector<qlat::vector_acc<Ty > >& Eprop, qlat::vector_acc<Ty >& Eres, 
  fft_desc_basic& fd, qlat::vector_gpu<Ty >& resTa, int clear = 1)
{
  prop_to_vec(Eprop, resTa, fd);  
  vec_corrE(resTa, Eres, fd, 0);
}

/////propH , pi --> 12*12 --> Nt*Nxyz
////input propH with or without sink smear
////      src without smear, mode_sm = 2 smtopt, mode_sm = 3  smtosm
template<typename Ty, typename Ta>
void point_corr(qnoiT& src, std::vector<qpropT >& propH,
    std::vector<double >& massL, eigen_ov& ei, fft_desc_basic& fd, corr_dat<Ta >& res, lms_para& srcI, int shift_t = 1)
{
  TIMER("point corr");
  print_time();
  print_mem_info("Before point_corr");
  if(propH.size() != massL.size()){abort_r("prop size, mass size not match!\n");}
  ////int Ns = src.size();
  const long Size_prop = fd.get_prop_size();
  const Geometry& geo = src.geo();

  const int GPU = 1;const bool rotate = false;
  const int nmass = massL.size();

  Coordinate Lat;for(int i=0;i<4;i++){Lat[i] = fd.nv[i];}
  Coordinate pos;Coordinate off_L;
  std::vector<PointSelection > Ngrid;
  check_noise_pos(src, pos, off_L);
  grid_list_posT(Ngrid, off_L, pos, srcI.combineT, Lat);
  for(unsigned int ic=0;ic<Ngrid.size();ic++)
  {
    for(unsigned int id=0;id<Ngrid[ic].size();id++){
      print0("%s ", qlat::show(Ngrid[ic][id]).c_str());
    }
    print0(" \n");
  }
  const int tini = pos[3];

  srcI.off_L   = off_L;
  srcI.ini_pos = pos;
  srcI.print();

  //////===container for fft vecs
  ////int off_FFT = 0;
  bool saveFFT = false;bool savezero = false;bool save_zero_corr = false;
  if(srcI.name_mom_vecs != std::string("NONE")){saveFFT = true;}
  if(srcI.name_zero_vecs != std::string("NONE")){savezero = true;}
  if(srcI.save_zero_corr == 1){save_zero_corr = true;}
  ////save fft vecs (all of them, no need to shrink mem now)

  /////need modification of positions
  qlat::vector_acc<Ty > EresH;qlat::vector_acc<Ty > EresL;qlat::vector_acc<Ty > EresA;
  if(save_zero_corr){
  EresH.resize(32 * nmass * fd.nt);clear_qv(EresH );
  EresL.resize(32 * nmass * fd.nt);clear_qv(EresL );
  EresA.resize(32 * nmass * fd.nt);clear_qv(EresA );}
  qlat::vector_gpu<Ty > resTa;const int nvecs = 32 * nmass;
  qlat::vector_gpu<Ty > resZero;

  /////do corr
  std::vector<qlat::vector_acc<Ty > > Eprop;

  qlat::vector_gpu<Ty > stmp, low_prop, high_prop;
  /////copy src vectors
  stmp.resize(Size_prop);

  copy_eigen_src_to_FieldM(high_prop, propH, ei.b_size, fd, 1, GPU, rotate);
  low_prop.resize(high_prop.size());

  std::vector<qlat::FieldM<Ty , 12*12> > src_prop;src_prop.resize(1);
  /////set memory for low_prop
  int Nlms = 1;
  if(srcI.lms == -1){Nlms = Ngrid.size();}
  if(srcI.lms == 0 ){Nlms = 1;}
  if(srcI.lms >  0 ){Nlms = srcI.lms;}

  const int mc = srcI.mom_cut*2 + 1;
  momentum_dat mdat(geo, srcI.mom_cut);
  std::vector<qlat::vector_gpu<Ty > > FFT_data;FFT_data.resize(1 + 1);
  //std::vector<qlat::vector_gpu<Ty > > FFT_data;FFT_data.resize(1 + Nlms);
  //qlat::vector_gpu<Ty > FFT_data;long Nfdata = long(32)*massL.size()*fd.nt*mc*mc*mc ;
  //qlat::vector_gpu<Ty > FFT_data_high;
  //qlat::vector_acc<long > mapA, mapB;

  ///////check production for this source
  ///////low mode ignored if check point enabled
  if(savezero){
    bool flag_do_job = false;
    if(get_file_size_MPI(srcI.name_zero_vecs.c_str()) == 0){flag_do_job = true;}
    if(flag_do_job == false){
      print0("Pass %s \n", srcI.name_zero_vecs.c_str());
      return ;
    }
  }

  //if(saveFFT){
  //  bool flag_do_job = false;
  //  TIMER("saveFFT");
  //  fft_fieldM(resTa.data(), 32*nmass, 1, geo, false);
  //  mdat.pick_mom_from_vecs(FFT_data[0], resTa);
  //  char namesm[500], nameQ[550], name_info[550], name_mom_tmp[100];
  //  sprintf(name_sm, "%s", srcI.name_mom_vecs.c_str());
  //  sprintf(nameQ, "%s.Gsrc", namesm );
  //  sprintf(name_info, "%s.GInfo", namesm );

  //  if(get_file_size_MPI(name_mom) == 0){flag_do_job = true;}
  //  if(get_file_size_MPI(srcI.name_zero_vecs.c_str()) == 0){flag_do_job = true;}
  //  for(int gi=0;  gi  < Nlms + 1; gi++){
  //    sprintf(name_mom_tmp, "%09d", gi);
  //    if(mdat.read(FFT_data, std::string(nameQ), std::string(namei)) == 0){flag_do_job = true;}
  //  }
  //  if(flag_do_job == false){return ;}
  //}

  char key_T[1000], dimN[1000];
  sprintf(key_T, "%d", 1);sprintf(dimN , "src");

  io_vec io_use(geo,srcI.ionum);
  std::string POS_LIST, POS_CUR;
  if(saveFFT){
    //get_mom_single_node(mapA, mapB, geo, srcI.mom_cut);
    //FFT_data.resize(Nfdata);FFT_data.set_zero();
    //FFT_data_high.resize(Nfdata);FFT_data_high.set_zero();
    //sprintf(key_T, "  %d   %d  %d  %d  %d %d %d %d", Nlms+1, 32, int(massL.size()), fd.nt, mc, mc, mc, 2);
    //sprintf(dimN , "grids+1 operator masses nt pz py px complex");
    sprintf(key_T, "%d  %d  %d %d %d %d", int(massL.size()), fd.nt, mc, mc, mc, 2);
    sprintf(dimN , "masses nt pz py px complex");
  }
  std::string ktem(key_T);
  std::string dtem(dimN);
  corr_dat<Ta > mom_res(ktem, dtem);
  mom_res.INFOA = srcI.INFOA;
  char  name_mom[500], name_mom_tmp[500];
  sprintf(name_mom, "%s.Gsrc", srcI.name_mom_vecs.c_str());
  //////===container for fft vecs

  ////===high mode contractions
  /////reduce the high prop
  if(srcI.SRC_PROP_WITH_LOW == 1){
    stmp.set_zero();
    low_prop.set_zero();
    FieldM_src_to_FieldM_prop(src, src_prop[0], GPU);
    copy_eigen_src_to_FieldM(stmp, src_prop, ei.b_size, fd, 1, GPU, rotate);
    prop_L_device(ei, stmp.data(), low_prop.data(), 12, massL, srcI.mode_eig_sm);
    high_prop -= low_prop;
  }
  /////reduce the high prop

  print0("Do high %s \n", POS_CUR.c_str());
  copy_eigen_prop_to_EigenM(high_prop.data(), Eprop, ei.b_size, nmass, fd, 0, GPU);
  prop_to_vec(Eprop, resTa, fd); 
  if(save_zero_corr){vec_corrE(resTa.data(), EresH, fd, nvecs, 0);}

  POS_CUR = std::string("");write_pos_to_string(POS_CUR, pos);POS_LIST += POS_CUR;
  if(savezero){
    resZero.resize(resTa.size());
    cpy_data_thread(resZero.data(), resTa.data(), resTa.size(), 1, true, -1.0*Nlms);
  }
  if(saveFFT){
    TIMER("saveFFT");
    fft_fieldM(resTa.data(), 32*nmass, 1, geo, false);
    mdat.pick_mom_from_vecs(FFT_data[0], resTa);
    sprintf(name_mom_tmp, "%09d", 0);mdat.write( FFT_data[0], std::string(name_mom), name_mom_tmp, true );

    //print0("norm0 %.8f \n", FFT_data[0].norm().real());

    //sprintf(key_T, "%d  %d  %d %d %d %d",  int(massL.size()), fd.nt, mc, mc, mc, 2);
    //sprintf(dimN , "masses nt pz py px complex");

    //std::string ktem(key_T);
    //std::string dtem(dimN);
    //corr_dat<Ta > mom_tmp(ktem, dtem);
    //mom_tmp.INFOA = srcI.INFOA;

    //long Nfdata = long(32)*massL.size()*fd.nt*mc*mc*mc ;
    //qlat::vector_gpu<Ty > FFT_data_high;FFT_data_high.resize(Nfdata);
    //pick_mom_data(FFT_data_high, resTa, 32*nmass, mdat.mapA, mdat.mapB, geo);
    //sprintf(name_mom, "%s.%05d.Gsrc", srcI.name_mom_vecs.c_str(),    0);
    //mom_tmp.INFO_LIST = POS_CUR;mom_tmp.shift_zero();
    //mom_tmp.write_corr(FFT_data_high.data(), FFT_data_high.size(), 3);
    //mom_tmp.write_dat(name_mom);
  }

  print_mem_info();

  for(int gi=0;gi<Nlms;gi++){
    POS_CUR = std::string("");write_pos_to_string(POS_CUR, Ngrid[gi][0]);POS_LIST += POS_CUR;
    print0("Do, Nlms %5d, gi %5d, %s \n", Nlms, gi, POS_CUR.c_str());
    fflush_MPI();

    stmp.set_zero();
    low_prop.set_zero();

    if(srcI.lms == 0)
    {
      FieldM_src_to_FieldM_prop(src, src_prop[0], GPU);
      copy_eigen_src_to_FieldM(stmp, src_prop, ei.b_size, fd, 1, GPU, rotate);
    }

    if(srcI.lms != 0){
      ////print0("gsize %d \n", int(Ngrid[gi].size()));
      write_grid_point_to_src(stmp.data(), src, Ngrid[gi], ei.b_size, fd);
    }

    /////get low mode prop
    prop_L_device(ei, stmp.data(), low_prop.data(), 12, massL, srcI.mode_eig_sm);

    //////format suitable for contractions
    ////low mode corr contractions
    if(srcI.do_all_low == 1){
      copy_eigen_prop_to_EigenM(low_prop.data(),  Eprop, ei.b_size, nmass, fd, 0, GPU);
      prop_to_vec(Eprop, resTa, fd);  
      if(save_zero_corr){vec_corrE(resTa.data(), EresL, fd, nvecs, 0);}
      //Ty norm0 = low_prop.norm();
      //Ty norm1 = resTa.norm();
      //print0("Check value %.3f %.3f, %.3f %.3f, %.3f %.3f \n", EresL[0].real(), EresL[0].imag(), 
      //  norm0.real(), norm0.imag(), norm1.real(), norm1.imag());
    }

    low_prop += high_prop;
    //////format suitable for contractions
    copy_eigen_prop_to_EigenM(low_prop.data(),  Eprop, ei.b_size, nmass, fd, 0, GPU);
    prop_to_vec(Eprop, resTa, fd);
    if(save_zero_corr){vec_corrE(resTa.data(), EresA, fd, nvecs, 0);}

    if(savezero){resZero += resTa;}
    if(saveFFT){
      TIMER("saveFFT");
      fft_fieldM(resTa.data(), 32*nmass, 1, geo, false);
      mdat.pick_mom_from_vecs(FFT_data[1], resTa);
      FFT_data[1] -= FFT_data[0];
      sprintf(name_mom_tmp, "%09d", gi + 1);mdat.write( FFT_data[1], std::string(name_mom), name_mom_tmp, false );
      ////print0("norm0 %.8f \n", FFT_data[gi+1].norm().real());
      //pick_mom_data(FFT_data, resTa, 32*nmass, mapA, mapB, geo);
      //sprintf(name_mom, "%s.%05d.Gsrc", srcI.name_mom_vecs.c_str(), gi+1);
      //mom_res.INFO_LIST = POS_CUR;mom_res.shift_zero();
      //////subtract high mode
      //mom_res.write_corr(FFT_data.data(), FFT_data.size(), 3);
      //mom_res.write_dat(name_mom);
    }
  }

  if(saveFFT){
    TIMER("saveFFT");
    //sprintf(name_mom, "%s.Gsrc", srcI.name_mom_vecs.c_str());
    //char name[500];bool clean = true;
    //for(unsigned int i=0;i<FFT_data.size();i++){
    //  sprintf(name, "%09d", i);if(i != 0){clean = false;}
    //  mdat.write( FFT_data[i], std::string(name_mom), name, clean );
    //}

    sprintf(name_mom, "%s.GInfo", srcI.name_mom_vecs.c_str());
    mom_res.INFO_LIST = POS_LIST;
    mom_res.write_dat(name_mom);
  }

  if(savezero){
    std::vector<qlat::FieldM<Ty, 1> > Vzero_data;
    Vzero_data.resize(32*nmass);
    for(unsigned int iv=0;iv<Vzero_data.size();iv++){
      Vzero_data[iv].init(geo);Ty* res = (Ty*) qlat::get_data(Vzero_data[iv]).data();
      cpy_data_thread(res, &resZero[iv*geo.local_volume()], geo.local_volume(), 1, true);
    }
    save_qlat_noises(srcI.name_zero_vecs.c_str(), Vzero_data, true, POS_LIST);
  }

  //if(saveFFT){
  //  mom_res.INFO_LIST = POS_LIST;
  //  mom_res.write_dat(srcI.name_mom_vecs.c_str());
  //  mom_res.print_info();
  //}

  if(save_zero_corr){
    if(shift_t == 1){
      shift_result_t(EresL, fd.nt, tini);
      shift_result_t(EresH, fd.nt, tini);
      shift_result_t(EresA, fd.nt, tini);
    }
    ////subtract high mode from EresA
    cpy_data_thread(EresA.data(), EresH.data(), EresA.size(), 1, true, -1.0*Nlms);
    res.write_corr((Ty*) EresL.data(), EresL.size());
    res.write_corr((Ty*) EresH.data(), EresH.size());
    res.write_corr((Ty*) EresA.data(), EresA.size());
  }

  print_mem_info("END point_corr");
}

/////all to all prop naive do, test for all-to-all low eigen
template<typename Ta>
void test_all_prop_corr(std::vector<double >& massL, eigen_ov& ei, fft_desc_basic& fd, corr_dat<Ta >& res, int mode_sm = 0)
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
  qlat::vector_gpu<Complexq > resTa;
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
    prop_to_corr_mom0(Eprop, Eres, fd, resTa,  0); 
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


}

#endif

