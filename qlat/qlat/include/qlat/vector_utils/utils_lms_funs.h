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
#include "utils_eigen_ov.h"
#include "utils_construction.h"
#include "utils_FFT_GPU.h"
#include "utils_grid_src.h"

////#define LOCALUSEACC 1

namespace qlat{

////buffers of lms included
template<typename Ty >
struct lms_para{

  Coordinate ini_pos;
  Coordinate off_L;

  Int SRC_PROP_WITH_LOW ;
  Int lms;
  Int combineT;

  Int src_smear_iter;
  double src_smear_kappa;
  Int sink_smear_iter;
  double sink_smear_kappa;

  Int mode_eig_sm;
  Int ckpoint;

  Int do_all_low ;

  Int mom_cut;
  Int save_zero_corr;
  Int save_full_vec;
  Int save_low_fft;


  // sequential related parameters and buffers
  Int Ngroup_seq;
  double sink_hfactor; 
  //eigen related
  eigen_ov* eig_P;
  double eig_mass;
  Int    eig_mode_sm;
  Int    eig_one_minus_halfD;
  std::vector<FieldG<Ty> > bufV ;
  std::vector<qlat::FieldG<Ty> > GresP;
  std::vector< std::vector<qlat::FieldG<Ty> > > shift_srcL;
  std::vector< std::vector<qlat::FieldG<Ty> > > srcL_pointers;
  std::vector< std::vector<qlat::FieldG<Ty> > > resL_pointers;
  std::vector<qlat::FieldG<Ty> > buf_res_props;// buffer for sink to current results, size Ngroup * Nlh
  std::vector<std::vector<FieldG<Ty > >> buf_shift;//sizes 2, Nsrc, 144 * Vol
  std::vector< std::vector<qlat::FieldG<Ty > > > curr_sbuf;
  std::vector< std::vector<qlat::FieldG<Ty > > > curr_s1L_buf;

  std::vector<FieldG<Ty> > buf_src ;// sizes of Ngroup, maybe larger than 4/8
  std::vector<FieldG<Ty> > buf_res ;// sizes of Ngroup, maybe larger than 4/8
  vector_gpu<Ty> buf_seq_prop;// buffer for seq_low_vecs
  //vector_gpu<Ty> buf_vec;// buffer for seq_high_vecs
  std::vector<vector_gpu<Ty>> sum_d  ;// buffer for seq_high_vecs
  std::vector<vector_gpu<Ty > > sink_phases;
  std::vector<std::vector<FieldG<Ty> >> curr_shift_buf;

  /* 
    will not be cleared once loaded
    noise and time slice for loaded props
    name_props   : name of noise and props
    //tag_props    : whether it's src or sink
    noise_props  : noises 
    prop_tposL   : tlist of the noise
    prop_npointsL : how many points in all time slice
  */
  std::vector<std::string > name_props ;
  std::vector<FieldG<Ty> > noise_props ;
  std::vector<std::vector<Int >    > prop_tposL ;
  std::vector<Long > prop_npointsL;
  std::vector<std::vector<Coordinate> > prop_coordL;
  std::vector<FieldG<Ty> > full_props ;
  std::vector<FieldG<Ty> > props_bufs ;

  // end

  FieldSelection fsel;
  ShuffledBitSet sbs ;

  std::string sf_tag;
  std::vector<qlat::SelectedField<Ty > > sfL;
  qlat::SelectedField<Ty > sf_single;
  Coordinate new_size_node;

  Int check_prop_norm;
  Int sparse_src;
  Int do_hadron_contra;

  std::string name_sparse_prop;
  std::string name_mom_vecs;
  std::string name_zero_vecs;
  std::string name_zero;

  std::string INFO;
  std::vector<std::string > INFOA;

  ////buffers
  qlat::vector<Ty > EresH;qlat::vector<Ty > EresL;qlat::vector<Ty > EresA;
  qlat::vector_gpu<Ty > stmp, low_prop, high_prop;
  std::vector<qlat::FieldM<Ty , 12*12> > src_prop;
  std::vector<qlat::vector_gpu<Ty > > FFT_data;
  qlat::vector_gpu<Ty > resTa;
  qlat::vector_gpu<Ty > resZero;
  std::vector<qlat::FieldM<Ty, 1> > Vzero_data;
  std::vector<qlat::vector_gpu<Ty > > Eprop;
  qlat::Field<Ty > tmp_prop;

  /////initial with src noise

  inline void free_buf(){
    EresH.resize(0);
    EresL.resize(0);
    EresA.resize(0);
    stmp.resize(0);
    low_prop.resize(0);
    high_prop.resize(0);
    src_prop.resize(0);
    FFT_data.resize(0);
    resTa.resize(0);
    resZero.resize(0);
    Vzero_data.resize(0);
    Eprop.resize(0);

    // sequential related parameters and buffers
    bufV.resize(0);
    GresP.resize(0);
    shift_srcL.resize(0);
    srcL_pointers.resize(0);
    resL_pointers.resize(0);
    buf_res_props.resize(0);// buffer for sink to current results, size Ngroup * Nlh
    buf_shift.resize(0);//sizes 2, Nsrc, 144 * Vol
    curr_sbuf.resize(0);
    curr_s1L_buf.resize(0);

    buf_src.resize(0);
    buf_res.resize(0);
    buf_seq_prop.resize(0);// buffer for seq_low_vecs
    //buf_vec.resize(0);// buffer for seq_high_vecs
    sum_d.resize(0) ;// buffer for seq_high_vecs
    sink_phases.resize(0);
    curr_shift_buf.resize(0);
    // end
  }

  inline void init(){
    for(Int i=0;i<4;i++){ini_pos[i] = 0;off_L[i] = 1;}

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
    sparse_src       = 0;
    save_low_fft     = 0;

    // sequential related parameters and buffers
    Ngroup_seq       = 8  ;
    sink_hfactor     = 1.0;

    //eigen related
    eig_P = NULL;
    eig_mass = 0.0;
    eig_mode_sm = 0;
    eig_one_minus_halfD = 1;

    // only save sparse prop
    do_hadron_contra = 1;

    mom_cut          = 4;
    ckpoint          = 1;

    name_sparse_prop = std::string("NONE");
    name_mom_vecs    = std::string("NONE");
    name_zero_vecs   = std::string("NONE");
    name_zero        = std::string("NONE");
    save_zero_corr   = 1;
    save_full_vec    = 0;////1 if all grid source vec need to be saved

    Int ionum = 16;
    std::string val = get_env(std::string("q_io_sparse_vec_ionum"));
    if(val == ""){ionum = 16;}else{
      Int tem = stringtonum(val);
      if(tem <= 8){ionum = 8;}
      if(tem > 8  and tem <= 16){ionum = 16;}
      if(tem > 16){ionum = 32;}
    }
    Qassert(ionum == 8 or ionum == 16 or ionum == 32);
  
    if(ionum == 8){new_size_node = Coordinate(1, 2, 2, 2);}
    if(ionum ==16){new_size_node = Coordinate(1, 2, 2, 4);}
    if(ionum ==32){new_size_node = Coordinate(2, 2, 2, 4);}
    
    INFO = std::string("NONE");
    INFOA.resize(0);
  }

  inline void set_eig(eigen_ov* eig, const double mass, const Int mode_sm, const Int one_minus = 1){
    Qassert(mass >= 0.0);
    Qassert(mode_sm >= 0);
    Qassert(one_minus == 1);
    eig_P = eig;
    eig_mass = mass;
    eig_mode_sm = mode_sm;
    eig_one_minus_halfD = one_minus;
  }

  inline void print(){
    qmessage("===Source Info %s \n", INFO.c_str());
    qmessage("  prop with low %3d, lms %3d, combineT %3d \n", SRC_PROP_WITH_LOW, lms, combineT);
    qmessage("  init %5d %5d %5d %5d, grid %5d %5d %5d %5d \n", 
      ini_pos[0], ini_pos[1], ini_pos[2], ini_pos[3],
      off_L[0], off_L[1], off_L[2], off_L[3] );

    Int save_vecs = 0; if(name_mom_vecs != std::string("NONE")){save_vecs = 1;}
    qmessage("eigen mode sm %1d, src %5d %.3f, sink %5d %.3f , save_low %2d, mom_cut %5d , saveV %1d \n", 
      mode_eig_sm, src_smear_iter, src_smear_kappa, sink_smear_iter, sink_smear_kappa, do_all_low, mom_cut, save_vecs);
    //qmessage("  mom vecs %s \n" ,name_mom_vecs);
    //qmessage("  zero vecs %s \n",name_zero_vecs);
  }

};

template<typename Ty>
void check_nan_GPU(qlat::vector_gpu<Ty >& resTa)
{
  TIMERA("check_nan_GPU");
  const Int GPU = resTa.GPU;
  const Long V  = resTa.size();
  const Ty* s   = resTa.data();
  const Long Bfac  = 64;
  const Long Neach = (V + Bfac - 1 ) / Bfac;
  qlat::vector<double > buf;buf.resize(Neach);
  qGPU_for(ieach, Neach, GPU, {
    buf[ieach] = 0;
    for(Long bi=0;bi<Bfac;bi++)
    {
      const Long isp = ieach * Bfac + bi;
      if(isp < V){
        if(qisnan(s[isp])){
          buf[ieach] = 1;
          break;
        }
      }
    }
  })
  for(Int ieach=0;ieach<Neach;ieach++)
  {
    if(buf[ieach] == 1)
    {
      qqwarn(ssprintf("WARNING: isnan corr"));
      break;
    }
  }
}

//template<typename Ty>
//void pick_mom_data(qlat::vector_gpu<Ty >& res, qlat::vector_gpu<Ty >& src,
//  const Int nvec, qlat::vector<Long >& mapA, qlat::vector<Long >& mapB, const Geometry& geo)
//{
//  TIMER("save FFT");
//  ////false to match to gwu code convention
//  ////exp(+ i p y) need src exp(-i p x)
//  fft_fieldM(src.data(), nvec, 1, geo, false);
//  res.set_zero();
//
//  const Long Nvol = src.size()/nvec;
//  const Long Mvol = res.size()/nvec;
//  Qassert(mapA.size() <= Mvol);
//
//  Long* A = mapA.data();
//  Long* B = mapB.data();
//  qacc_for(isp, mapA.size() ,{
//    Long i0 = A[isp];
//    Long i1 = B[isp];
//    for(Int iv=0;iv<nvec;iv++){res[iv*Mvol + i1] = src[iv*Nvol + i0];}
//  });
//
//  sum_all_size(res.data(), res.size(), 1);
//}


template<typename Ty>
void save_sink_vecs(qlat::vector_gpu<Ty >& src, const Geometry& geo, const Int nmass, lms_para<Ty >& srcI, 
  const std::string& GRID_INFO, const std::string tag = "")
{
  TIMER("lms savezero");
  std::vector<qlat::FieldM<Ty, 1> >& Vzero_data = srcI.Vzero_data;
  const size_t Vol = geo.local_volume();
  Qassert(src.size() == 32*nmass*Vol);
  const Long nvec = src.size()/Vol;
  //qmessage("=====vec %d \n", int(nvec));
  if(Long(Vzero_data.size() ) != nvec){
    Vzero_data.resize(0);Vzero_data.resize(nvec);
    for(Long iv=0;iv<nvec;iv++){
      if(!Vzero_data[iv].initialized){Vzero_data[iv].init(geo);}
    }
  }
  for(Long iv=0;iv<nvec;iv++){
    Ty* resP = (Ty*) qlat::get_data(Vzero_data[iv]).data();
    cpy_data_thread(resP, &src[iv*Vol], geo.local_volume(), 1, QTRUE);
  }

  std::string namew = srcI.name_zero_vecs;
  if(tag != std::string("")){
    namew = ssprintf("%s.%s", srcI.name_zero_vecs.c_str(), tag.c_str());
  }
  save_qlat_noises(namew.c_str(), Vzero_data, true, GRID_INFO);
}

template<typename Ty>
void Save_sparse_prop(std::vector<qlat::vector_gpu<Ty > >& src, lms_para<Ty >& srcI, const std::string& tag, const Int count, const bool save_files = false)
{
  if(srcI.name_sparse_prop == std::string("NONE")){return ;};
  TIMER("Save_sparse_prop");
  Qassert(srcI.fsel.n_elems > 0);
  const Int Nmass = src.size();

  const Long Ndc = 12 * 12;
  Qassert(int(srcI.sfL.size()) == Nmass);
  Qassert(srcI.sfL[0].multiplicity / Ndc > count );
  Qassert(srcI.sf_single.field.size() % Ndc == 0);
  const Long Ndata = srcI.sf_single.field.size() / (Ndc);
  Qassert(srcI.sfL[0].field.size() % Ndc * Ndata == 0);
  //Ty* p0 = (Ty*) qlat::get_data(srcI.sf_single.field).data();
  //Ty* pr = (Ty*) qlat::get_data(srcI.sf.field).data();

  for(Int mi=0;mi<Nmass;mi++)
  {
    ////if(mi != 0){append = true;}
    prop_gpu_to_qprop(srcI.tmp_prop, src[mi]);
    Qassert(srcI.tmp_prop.multiplicity == Ndc);
    set_selected_field(srcI.sf_single, srcI.tmp_prop, srcI.fsel);
    copy_sparse_fields(srcI.sfL[mi], srcI.sf_single, Ndc, count, 0);

    //copy back to src to test results ...
    //clear_fields(srcI.tmp_prop);
    //set_field_selected(srcI.tmp_prop, srcI.sf_single, srcI.fsel, true);
    //qprop_to_prop_gpu(src[mi], srcI.tmp_prop);
  }
  if(tag.size() > 0){
    srcI.sf_tag += tag + "_";
  }

  // save each mass of sparse prop 
  if(save_files){
    // save files
    bool append = false;
    for(Int mi=0;mi<Nmass;mi++){
      const std::string nameQ = ssprintf("%s.mass%02d", srcI.name_sparse_prop.c_str(), mi);
      //bool append = true;if(tag == "-1" or clean == true){append = false;}
      if(append == false)
      {
        if(0 == qlat::get_id_node()){
          qlat::qremove_all(nameQ);
        }
      }

      ShuffledFieldsWriter sfw(nameQ, srcI.new_size_node, append);
      //std::string tag_ = ssprintf("%s", tag.c_str());
      //std::string tag_ = ssprintf("%s.mass%02d", tag.c_str(), mi);
      const std::string tag_ =  srcI.sf_tag.substr(0, srcI.sf_tag.size() - 1);
      //  default single precision files
      qlat::write(sfw, tag_, srcI.sbs, srcI.sfL[mi]);
      sfw.close();
    }
    srcI.sf_tag = "";
  }
}

/////propH , pi --> 12*12 --> Nt*Nxyz
////input propH with or without sink smear
////      src without smear, mode_sm = 2 smtopt, mode_sm = 3  smtosm
template<typename Ty, typename Ta>
void point_corr(qnoiT& src, std::vector<qpropT >& propH,
    std::vector<double >& massL, eigen_ov& ei, fft_desc_basic& fd, corr_dat<Ta >& res, lms_para<Ty >& srcI, momentum_dat& mdat, Int shift_t = 1)
{
  TIMER("point corr");
  print_time();
  print_mem_info("Before point_corr");
  if(propH.size() != massL.size()){abort_r("prop size, mass size not match!\n");}
  ////int Ns = src.size();
  const Long Size_prop = fd.get_prop_size();
  const Geometry& geo = src.geo();

  const Int GPU = 1;const bool rotate = false;
  const Int nmass = massL.size();
  const size_t vol = size_t(fd.nx) * fd.ny * fd.nz * fd.nt;
  const Int do_hadron_contra = srcI.do_hadron_contra;
  //const size_t Vol = geo.local_volume();

  Coordinate Lat;for(Int i=0;i<4;i++){Lat[i] = fd.nv[i];}
  Coordinate pos;Coordinate off_L;
  std::vector<PointsSelection > Ngrid;
  {
  TIMER("check noise");
  if(srcI.sparse_src == 0){
    check_noise_pos(src, pos, off_L);
    grid_list_posT(Ngrid, off_L, pos, srcI.combineT, Lat);
  }else{
    // sparse src, no grid off at all
    for(Int i=0;i<4;i++){off_L[i] = 0;pos[i] = 0;}

    // will not work for combineT
    std::vector<Coordinate > grids;
    std::vector<Int > Zlist;
    //get_noise_pos(src, grids, Zlist, 3, 1);
    get_noise_pos(src, grids, Zlist);
    pos = grids[0];
    Ngrid.resize(grids.size());
    std::vector<Int > tlist;
    for(LInt gi=0;gi<grids.size();gi++){
      Ngrid[gi].resize(1);
      Ngrid[gi][0] = grids[gi];
      const Int tini = grids[gi][3];
      if(std_find(tlist, tini) < 0){tlist.push_back(tini);}
    }
    // sparsen one could not do  combineT
    if(srcI.combineT == 1){Qassert(tlist.size() == 1)}
  }
  }
  // grid info
  std::string GRID_INFO, tmp;
  {
    tmp = "";
    write_pos_to_string(tmp, Lat);
    GRID_INFO += ssprintf(" Lat %s", tmp.c_str());

    tmp = "";
    write_pos_to_string(tmp, pos);
    GRID_INFO += ssprintf(" pos %s", tmp.c_str());

    tmp = "";
    write_pos_to_string(tmp, off_L);
    GRID_INFO += ssprintf(" grid %s", tmp.c_str());

    GRID_INFO += ssprintf(" combineT %1d ", srcI.combineT);
  }

  // only do eigen props on desire src time slices
  std::vector<Int > tsrcL;tsrcL.resize(0);
  for(unsigned int ic=0;ic<Ngrid.size();ic++)
  {
    for(unsigned int id=0;id<Ngrid[ic].size();id++){
      qmessage("%s ", qlat::show(Ngrid[ic][id]).c_str());
      const Int tcur = Ngrid[ic][id][3];
      bool add_time = true;
      for(unsigned int si=0;si<tsrcL.size();si++){
        if(tsrcL[si] == tcur){add_time = false;}
      }
      if(add_time){tsrcL.push_back(tcur);}
    }
    qmessage(" \n");
  }
  //for(Int si=0;si<int(tsrcL.size());si++)
  //{
  //  qmessage("Low sink time si %3d : %6d \n", si, int(tsrcL[si]));
  //}
  const Int tini = pos[3];

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

  qlat::vector<Ty >& EresH = srcI.EresH;
  qlat::vector<Ty >& EresL = srcI.EresL;
  qlat::vector<Ty >& EresA = srcI.EresA;

  qlat::vector_gpu<Ty >& stmp = srcI.stmp;
  qlat::vector_gpu<Ty >& low_prop = srcI.low_prop;
  qlat::vector_gpu<Ty >& high_prop = srcI.high_prop;
  std::vector<qlat::FieldM<Ty , 12*12> >& src_prop = srcI.src_prop;
  std::vector<qlat::vector_gpu<Ty > >& FFT_data = srcI.FFT_data;

  qlat::vector_gpu<Ty >& resTa = srcI.resTa;
  qlat::vector_gpu<Ty >& resZero = srcI.resZero;

  std::vector<qlat::vector_gpu<Ty > >& Eprop = srcI.Eprop;

  ////FFT_data do not need clean 
  if(src_prop.size() != 1){src_prop.resize(1);}
  if(FFT_data.size() != 2){FFT_data.resize(2);}
  ///qlat::vector_gpu<Ty > stmp, low_prop, high_prop;
  ///std::vector<qlat::FieldM<Ty , 12*12> > src_prop;src_prop.resize(1);
  ///std::vector<qlat::vector_gpu<Ty > > FFT_data;FFT_data.resize(1 + 1);

  /////need modification of positions
  //qlat::vector<Ty > EresH;qlat::vector<Ty > EresL;qlat::vector<Ty > EresA;
  //qlat::vector_gpu<Ty > resTa;

  if(save_zero_corr){
    EresH.resize(32 * nmass * fd.nt); 
    EresL.resize(32 * nmass * fd.nt); 
    EresA.resize(32 * nmass * fd.nt); 
    clear_qv(EresH , QFALSE);
    clear_qv(EresL , QFALSE);
    clear_qv(EresA , QFALSE);
    qacc_barrier(dummy);
  }
  const Int nvecs = 32 * nmass;
  /////do corr

  /////copy src vectors
  stmp.resize(Size_prop);

  //clean sparse prop
  copy_FieldM_to_bsize_prop(high_prop, propH, ei.b_size, fd, GPU, rotate);
  low_prop.resize(high_prop.size());
  if(srcI.check_prop_norm){qmessage("===high norm 0 %ld , ", long(high_prop.size()));high_prop.print_norm2();}

  //{
  //qmessage("===check norm");
  //Ty* tmp = (Ty*) qlat::get_data(propH[0]).data();
  //qmessage("value %+.8e %.8e ", tmp[0].real(), tmp[1].real());
  //high_prop.print_norm2();
  //}

  /////set memory for low_prop
  Int Nlms = 1;
  if(srcI.lms == -1){Nlms = Ngrid.size();}
  if(srcI.lms == 0 ){Nlms = 1;}
  if(srcI.lms >  0 ){Nlms = srcI.lms;}

  const Int mc = srcI.mom_cut*2 + 1;
  ////momentum_dat mdat(geo, srcI.mom_cut);
  //std::vector<qlat::vector_gpu<Ty > > FFT_data;FFT_data.resize(1 + Nlms);
  //qlat::vector_gpu<Ty > FFT_data;Long Nfdata = Long(32)*massL.size()*fd.nt*mc*mc*mc ;
  //qlat::vector_gpu<Ty > FFT_data_high;
  //qlat::vector<Long > mapA, mapB;

  // initialize sparse parameters
  if(srcI.name_sparse_prop != std::string("NONE")){
    Qassert(srcI.fsel.n_elems > 0);
    bool ini = false;
    const Int Ndc = 12 * 12 ;
    const Int Nsparse = 1 + Nlms * (1 + srcI.do_all_low);
    if(!srcI.sf_single.initialized){ini = true;}
    if(srcI.sfL.size() == 0){ini = true;}
    else{
      if(long(srcI.sfL.size()) != nmass){ini = true;}
      for(Int mi=0;mi<nmass;mi++){
        if(srcI.sfL[mi].multiplicity != Nsparse*Ndc){ini = true;}
      }
    }

    if(ini){
      srcI.sfL.resize(0);
      srcI.sfL.resize(nmass);
      for(Int mi=0;mi<nmass;mi++){
        srcI.sfL[mi].init(srcI.fsel, Nsparse*Ndc);
      }
      srcI.sf_single.init(srcI.fsel, Ndc);
      srcI.sbs = mk_shuffled_bitset(srcI.fsel, srcI.new_size_node);
      srcI.tmp_prop.init(geo, Ndc);
    }
  }

  ///////check production for this source
  ///////low mode ignored if check point enabled
  if(savezero and srcI.ckpoint == 1){
    bool flag_do_job = false;
    ////single prec assummed
    if(get_file_size_MPI(srcI.name_zero_vecs.c_str()) < vol * 32 * nmass * 8){
      flag_do_job = true;
    }
    if(flag_do_job == false){
      qmessage("Pass %s \n", srcI.name_zero_vecs.c_str());
      return ;
    }
  }

  //const Long nZero = 1;////number of saved zeros
  //if(srcI.save_full_vec == 1){nZero = 1 + Nlms;}

  //char key_T[1000], dimN[1000];

  std::string POS_LIST, POS_CUR;
  //char  name_mom[500], name_mom_tmp[500];
  // ssprintf(name_mom, "%s.Gsrc", srcI.name_mom_vecs.c_str());
  //////===container for fft vecs
  std::string name_mom = ssprintf("%s.Gsrc", srcI.name_mom_vecs.c_str());
  std::string name_mom_tmp ;

  //qmessage("===check norm");low_prop.print_norm2();high_prop.print_norm2();
  ////===high mode contractions
  /////reduce the high prop
  if(srcI.SRC_PROP_WITH_LOW == 1){
    stmp.set_zero();
    low_prop.set_zero();
    FieldM_src_to_FieldM_prop(src, src_prop[0], GPU);
    copy_FieldM_to_bsize_prop(stmp, src_prop, ei.b_size, fd, GPU, rotate);
    if(srcI.check_prop_norm){qmessage("===src  norm 0 %ld , ", long(stmp.size())); stmp.print_norm2();}
    prop_L_device(ei, stmp.data(), low_prop.data(), 12, massL, srcI.mode_eig_sm, tsrcL);
    if(srcI.check_prop_norm){qmessage("===low  norm 0 %ld , ", long(low_prop.size())); low_prop.print_norm2();}
    high_prop -= low_prop;
  }

  if(srcI.check_prop_norm){qmessage("===high norm 1 %ld , ", long(high_prop.size())); high_prop.print_norm2();}
  //qmessage("===check norm");low_prop.print_norm2();high_prop.print_norm2();
  /////reduce the high prop

  qmessage("Do high %s \n", POS_CUR.c_str());
  copy_eigen_prop_to_EigenG(Eprop, high_prop.data(), ei.b_size, nmass, fd, GPU);
  Save_sparse_prop(Eprop, srcI, std::string("High"), 0, false);
  if(do_hadron_contra){
    prop_to_vec(Eprop, resTa, fd); 
    if(save_zero_corr){vec_corrE(resTa.data(), EresH, fd, nvecs, 0);}
    if(srcI.check_prop_norm){qmessage("===resTa norm 0 %ld , ", long(resTa.size()));resTa.print_norm2();}
    if(savezero){
      resZero.resize(resTa.size());resZero.set_zero();
      if(srcI.save_full_vec == 0){
        cpy_data_thread(resZero.data(), resTa.data(), resTa.size(), 1, QTRUE, -1.0*Nlms);
      }
      if(srcI.save_full_vec == 1){
        cpy_data_thread(resZero.data(), resTa.data(), resTa.size(), 1, QTRUE);
        save_sink_vecs(resZero, geo, nmass, srcI, GRID_INFO);
      }
    }
    if(saveFFT){
      TIMER("saveFFT");
      check_nan_GPU(resTa);
      fft_fieldM(resTa.data(), 32*nmass, 1, geo, false);
      mdat.pick_mom_from_vecs(FFT_data[0], resTa);
      //// ssprintf(name_mom_tmp, "%09d", 0);
      name_mom_tmp = ssprintf("%09d", 0);
      //name_mom_tmp = ssprintf("High");
      mdat.write( FFT_data[0], name_mom, name_mom_tmp, true );
    }
  }
  POS_CUR = std::string("");write_pos_to_string(POS_CUR, pos);POS_LIST += POS_CUR;

  print_mem_info();

  for(Int gi=0;gi<Nlms;gi++){
    POS_CUR = std::string("");write_pos_to_string(POS_CUR, Ngrid[gi][0]);POS_LIST += POS_CUR;
    qmessage("Do, Nlms %5d, gi %5d, %s \n", Nlms, gi, POS_CUR.c_str());
    fflush_MPI();

    stmp.set_zero();
    low_prop.set_zero();

    if(srcI.lms == 0)
    {
      FieldM_src_to_FieldM_prop(src, src_prop[0], GPU);
      copy_FieldM_to_bsize_prop(stmp, src_prop, ei.b_size, fd, GPU, rotate);
    }

    if(srcI.lms != 0){
      ////qmessage("gsize %d \n", int(Ngrid[gi].size()));
      write_grid_point_to_src(stmp.data(), src, Ngrid[gi], ei.b_size, fd);
    }

    /////get low mode prop
    prop_L_device(ei, stmp.data(), low_prop.data(), 12, massL, srcI.mode_eig_sm, tsrcL);

    //////format suitable for contractions
    ////low mode corr contractions
    if(srcI.do_all_low == 1){
      copy_eigen_prop_to_EigenG(Eprop, low_prop.data(), ei.b_size, nmass, fd, GPU);
      std::string sparse_tag = "";
      if(gi == Nlms - 1){sparse_tag = ssprintf("Low%04d", Nlms);}
      Save_sparse_prop(Eprop, srcI, sparse_tag, Nlms + gi + 1, false);
      if(do_hadron_contra){
        prop_to_vec(Eprop, resTa, fd);  
        if(save_zero_corr){vec_corrE(resTa.data(), EresL, fd, nvecs, 0);}
        if(saveFFT and srcI.save_low_fft == 1){
          TIMER("saveFFT");
          check_nan_GPU(resTa);
          fft_fieldM(resTa.data(), 32*nmass, 1, geo, false);
          mdat.pick_mom_from_vecs(FFT_data[1], resTa);
          name_mom_tmp = ssprintf("Low%09d", gi + 1);
          mdat.write( FFT_data[1], std::string(name_mom), name_mom_tmp, false );
        }
      }

      //Ty norm0 = low_prop.norm();
      //Ty norm1 = resTa.norm();
      //qmessage("Check value %.3f %.3f, %.3f %.3f, %.3f %.3f \n", EresL[0].real(), EresL[0].imag(), 
      //  norm0.real(), norm0.imag(), norm1.real(), norm1.imag());
    }
    //qmessage("===check norm");low_prop.print_norm2();high_prop.print_norm2();

    low_prop += high_prop;
    //qmessage("===check norm");low_prop.print_norm2();high_prop.print_norm2();
    if(srcI.check_prop_norm){qmessage("===low  norm 1 %ld , ", long(low_prop.size()));low_prop.print_norm2();}

    //////format suitable for contractions
    copy_eigen_prop_to_EigenG(Eprop, low_prop.data(), ei.b_size, nmass, fd, GPU);
    if(gi != Nlms - 1){
      Save_sparse_prop(Eprop, srcI, "", gi + 1, false);
    }else{
      Save_sparse_prop(Eprop, srcI, ssprintf("LH%04d", Nlms), gi + 1, true );
    }

    if(do_hadron_contra){
      prop_to_vec(Eprop, resTa, fd);
      if(save_zero_corr){vec_corrE(resTa.data(), EresA, fd, nvecs, 0);}
      if(srcI.check_prop_norm){qmessage("===resTa norm 0 %ld , ", long(resTa.size()));resTa.print_norm2();}

      if(savezero){
        if(srcI.save_full_vec == 0){
          cpy_data_thread(resZero.data(), resTa.data(), resTa.size(), 1, QTRUE, +1.0);
        }
        if(srcI.save_full_vec == 1){
          cpy_data_thread(resZero.data(), resTa.data(), resTa.size(), 1, QTRUE);
          std::string tag = ssprintf("LH%04d", gi + 1);
          save_sink_vecs(resZero, geo, nmass, srcI, GRID_INFO, tag);
        }
      }
      if(saveFFT){
        TIMER("saveFFT");
        check_nan_GPU(resTa);
        fft_fieldM(resTa.data(), 32*nmass, 1, geo, false);
        mdat.pick_mom_from_vecs(FFT_data[1], resTa);
        FFT_data[1] -= FFT_data[0];
        name_mom_tmp = ssprintf("%09d", gi + 1);
        //name_mom_tmp = ssprintf("LH%09d", gi);
        mdat.write( FFT_data[1], std::string(name_mom), name_mom_tmp, false );
      }
    }
  }

  ////  if(saveFFT)
  // save info
  {
    TIMER("save info");
    //if(saveFFT)
    std::string key_T = ssprintf("%d", 1);
    std::string dimN  = ssprintf("src");
    {
      key_T = ssprintf("%d  %d  %d %d %d %d", int(massL.size()), fd.nt, mc, mc, mc, 2);
      dimN  = ssprintf("masses nt pz py px complex");
    }
    std::string ktem(key_T);
    std::string dtem(dimN);
    corr_dat<Ta > mom_res(ktem, dtem);
    mom_res.INFOA = srcI.INFOA;

    std::string name_info = ssprintf("%s.GInfo", srcI.name_mom_vecs.c_str());
    mom_res.INFO_LIST = POS_LIST;
    mom_res.INFOA.push_back( GRID_INFO );
    readuce_input_Coordinate_info(mom_res.INFO_LIST, mom_res.INFOA);
    mom_res.write_dat(name_info);
  }

  if(do_hadron_contra){
    if(savezero and srcI.save_full_vec == 0){
      save_sink_vecs(resZero, geo, nmass, srcI, GRID_INFO);
      //TIMER("lms savezero");
      //////std::vector<qlat::FieldM<Ty, 1> > Vzero_data;
      //Qassert(resZero0.size() == 32*nmass*Vol);
      //const Long nvec = resZero0.size()/Vol;
      //qmessage("=====vec %d \n", int(nvec));
      //if(Long(Vzero_data.size() ) != nvec){
      //  Vzero_data.resize(0);Vzero_data.resize(nvec);
      //  for(Long iv=0;iv<nvec;iv++){
      //    if(!Vzero_data[iv].initialized){Vzero_data[iv].init(geo);}
      //  }
      //}
      //for(Long iv=0;iv<nvec;iv++){
      //  Ty* resP = (Ty*) qlat::get_data(Vzero_data[iv]).data();
      //  cpy_data_thread(resP, &resZero[iv*Vol], geo.local_volume(), 1, QTRUE);
      //}

      ////save_qlat_noises(srcI.name_zero_vecs.c_str(), Vzero_data, true, POS_LIST);
      //save_qlat_noises(srcI.name_zero_vecs.c_str(), Vzero_data, true, GRID_INFO);
    }

    if(save_zero_corr){
      if(shift_t == 1){
        shift_result_t(EresL, fd.nt, tini);
        shift_result_t(EresH, fd.nt, tini);
        shift_result_t(EresA, fd.nt, tini);
      }
      ////subtract high mode from EresA
      cpy_data_thread(EresA.data(), EresH.data(), EresA.size(), 1, QTRUE, -1.0*Nlms);
      res.write_corr((Ty*) EresL.data(), EresL.size());
      res.write_corr((Ty*) EresH.data(), EresH.size());
      res.write_corr((Ty*) EresA.data(), EresA.size());
    }
  }

  print_mem_info("END point_corr");
}

/////all to all prop naive do, test for all-to-all low eigen
template<typename Ty>
void test_all_prop_corr(std::vector<double >& massL, eigen_ov& ei, fft_desc_basic& fd, corr_dat<Ty >& res, Int mode_sm = 0)
{
  /////if(propH.size() != src.size()*massL.size()){abort_r("prop size, mass size not match!\n");}
  const Geometry& geo = fd.geo();

  Int GPU = 1;
  Int nmass = massL.size();

  std::vector<qnoi > src; src.resize(1);src[0].init(geo);
  EigenV Eres;Eres.resize(nmass*32*fd.nt);
  qlat::set_zero(Eres);

  ///std::vector<Int > pos;pos.resize(src.size());qlat::vector<Int > off_L;
  ///check_noise_pos(src[0], pos[0], off_L);
  ///int tini = pos[0]%1000;
  qlat::vector_gpu<Ty > stmp, low_prop;
  stmp.resize(fd.get_prop_size());
  low_prop.resize(nmass*fd.get_prop_size());
  qlat::vector_gpu<Ty > resTa;
  /////do corr
  std::vector<qlat::vector_gpu<Ty > >  Eprop;

  Int tini = 0;
  for(Int zi=0;zi<fd.nz;zi++)
  for(Int yi=0;yi<fd.ny;yi++)
  for(Int xi=0;xi<fd.nx;xi++)
  {
    qmessage("===pos %d %d %d \n", zi, yi, xi);
    stmp.set_zero();
    low_prop.set_zero();

    LInt index0 = fd.index_g_from_local(0);
    LInt indexL = fd.index_g_from_g_coordinate(tini, zi, yi, xi);
    if(indexL >= index0 and indexL < index0 + fd.noden)
    {
      for(Int di=0;di<12;di++){
        stmp[(di*12 + di)*fd.noden + indexL%fd.noden] = 1.0;
      }
    }

    /////get low mode prop
    prop_L_device(ei, stmp.data(), low_prop.data(), 12, massL, mode_sm);

    copy_eigen_prop_to_EigenG(Eprop, low_prop.data(),  ei.b_size, nmass, fd, GPU);
    prop_to_corr_mom0(Eprop, Eres, fd, resTa,  0); 
    shift_result_t(Eres, fd.nt, tini);
  }
  for(Int i=0;i<Eres.size();i++){Eres[i] = Eres[i]/(Ftype(fd.nx*fd.ny*fd.nz*1.0));}
  res.write_corr((Ftype*) &Eres[0], 2*Eres.size());
}


}

#endif

