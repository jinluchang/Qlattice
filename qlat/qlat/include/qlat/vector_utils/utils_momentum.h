// utils_momentum.h
// Gen Wang
// Jun. 2022

#ifndef _UTILS_MOMENTUM_H_
#define _UTILS_MOMENTUM_H_

#pragma once

#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <qlat/selected-field-io.h>
#include <qlat/fields-io.h>
#include "utils_field_gpu.h"
///////#include "../kentucky/utils_lat_exchanger.h"

namespace qlat{

inline void get_mom_single_nodeA(qlat::vector<Long >& mapA,
    const Geometry& geo, const Int mom_cut, const std::vector< Coordinate>& mom_off = std::vector< Coordinate>(0))
{
  TIMERA("get_mom_single_nodeA");
  qlat::vector<Int > nv,Nv,mv;
  geo_to_nv(geo, nv, Nv, mv);

  qlat::vector<Long > mapB;
  mapB.resize(geo.local_volume());Long* Bi = mapB.data();
  ///const Int mc = mom_cut*2 + 1;

  for(unsigned long momG=0;momG<mom_off.size();momG++){
    Qassert(mom_off[momG][3] == 0);////only spatial momentum pick
    for(Int i=0;i<3;i++){
      Qassert(mom_off[momG][i] >= 0 and mom_off[momG][i] < nv[i]);
    }
  }
  qlat::vector<long> momL;momL.resize(mom_off.size() * 4);
  for(unsigned long momG=0;momG<mom_off.size();momG++){
    for(Int i=0;i<4;i++){
      momL[momG*4 + i] = mom_off[momG][i];
    }
  }
  if(momL.size() == 0){momL.resize(4);for(Int i=0;i<4;i++){momL[i] = 0;}}
  Qassert(momL.size() % 4 == 0);
  const long Nmom_off = momL.size() / 4;

  ////mom_off
  qacc_for(isp,  geo.local_volume(),{
    Coordinate xl  = geo.coordinate_from_index(isp);
    Coordinate xg  = geo.coordinate_g_from_l(xl);
    Int Nmom_flag = 0;
    //Coordinate xo  = xg;
    ////xg = xg - mom_off;
    for(long momG=0;momG<Nmom_off;momG++)
    {
      bool flag = true;
      //xo = xg;
      Coordinate mom  = xg;
      ////any direction larger than mom_cut, than disable
      for(Int i=0;i<3;i++){
        mom[i] = (xg[i] - momL[momG*4 + i] + nv[i])%nv[i];
      } ////move to positive
      for(Int i=0;i<3;i++){
        if(mom[i] > nv[i]/2){
          mom[i] = (nv[i] - mom[i]);
          if(mom[i] > mom_cut)
          {
            //mom[i] = mom_cut*2 + 1 - mom[i];
            flag = false;
          }
        }else{
          if(mom[i] > mom_cut){flag = false;}
        }
      }
      if(flag == true){Nmom_flag += 1;}
    }
    //if(flag == true ){ Bi[isp] = ((mom[3]*mc + mom[2])*mc + mom[1])*mc+ mom[0]; }
    if(Nmom_flag != 0){Bi[isp] =  1;}
    if(Nmom_flag == 0){Bi[isp] = -1;}
  });

  std::vector<Long > A0;//std::vector<Long > B0;
  for(Long isp=0;isp < mapB.size();isp++)
  {
    if(mapB[isp] >= 0){
      A0.push_back(  isp );//B0.push_back(mapB[isp]);
    }
  }
  mapA.resize(A0.size());//mapB.resize(A0.size());
  ////Long* A = mapA.data(); Long* B  = mapB.data();
  ////Long* Av = A0.data();  Long* Bv = B0.data();
  ////copy back to mapB
  qthread_for(isp,  mapA.size(),{
    mapA[isp] = A0[isp];// mapB[isp] = B0[isp];
  });
}

inline void get_mom_single_nodeB(qlat::vector<Long >& mapA, qlat::vector<Long >& mapB,
    const Geometry& geo, const Int mom_cut, const Coordinate& mom_off = Coordinate(0, 0, 0, 0))
{
  TIMERA("get_mom_single_nodeB");
  qlat::vector<Int > nv,Nv,mv;
  geo_to_nv(geo, nv, Nv, mv);

  const Long* A   = mapA.data();
  const Long Mvol = mapA.size();
  const Int mc = mom_cut*2 + 1;

  mapB.resize(Mvol);Long* Bi = mapB.data();
  Qassert(mom_off[3] == 0);////only spatial momentum pick
  for(Int i=0;i<3;i++){
    Qassert(mom_off[i] >= 0 and mom_off[i] < nv[i]);
  }

  qacc_for(isp,  Mvol,{
    const Long ilocal = A[isp];
    const Coordinate xl   = geo.coordinate_from_index(ilocal);
    Coordinate xg   = geo.coordinate_g_from_l(xl);
    ////any direction not less than mom_cut, disable
    bool flag = true;
    for(Int i=0;i<3;i++){
      xg[i] = (xg[i] - mom_off[i] + nv[i])%nv[i];
    } ////move to positive
    Coordinate mom  = xg;
    for(Int i=0;i<3;i++){
      if(xg[i] > nv[i]/2){
        mom[i] = (nv[i] - xg[i]);
        if(mom[i] >  mom_cut){flag = false;}
        else{mom[i] = mom_cut*2 + 1 - mom[i];}
      }else{if(mom[i] > mom_cut){flag = false;}}
    }
    if(flag == true ){Bi[isp] = ((mom[3]*mc + mom[2])*mc + mom[1])*mc+ mom[0]; }
    if(flag == false){Bi[isp] = -1;}
  });

}

#define TWOPT_TYPE  qlat::ComplexF

/////write date in float prec 
struct momentum_dat{
  ////int nx,ny,nz,nt;
  Geometry geo;
  Int mom_cut;
  ////Long Nvol;
  ////int  nt;

  /////corr_dat<Ftype > mom_res;

  /////default single precision writing and readings
  qlat::SelectedField<TWOPT_TYPE  > sf;
  qlat::SelectedField<TWOPT_TYPE  > sf_1;

  ///qlat fn_list
  std::string file_name;
  std::vector<std::string > fn_list;

  ////PointSelection pconf;
  FieldSelection fsel;
  FieldSelection fsel_1;
  ////write_float_from_double(sfw, tag, sf, sbs);
  qlat::vector<Long > mapA;
  qlat::vector<Long > fsel_map;
  ShuffledBitSet sbs ;
  Coordinate new_size_node;

  qlat::vector<Int > nv;
  Coordinate  cur_pos;
  qlat::vector_gpu<Complexq > phases;
  Coordinate cur_shift;
  /////int mc;
  ///int Mvol;

  ////save default mom_off
  long mapB_size ;
  qlat::vector<Long > mapB_buf;
  Coordinate mom_off_buf;

  Int nvec_copy;

  template<typename Ty > 
  void pick_mom_from_vecs(qlat::vector_gpu<Ty >& resF, qlat::vector_gpu<Ty >& resV){
    TIMER("pick_mom_from_vecs");
    const Long Mvol = mapA.size();
    ////const Long Nvol_ = V;
    const Long Nvol = geo.local_volume();
    Int nvec = resV.size()/Nvol;
    resF.resize(Mvol * nvec);
    ////sum_value_mpi(nvec);
    nvec_copy = nvec;

    const Long* A = mapA.data();
    Ty* resFP = resF.data();
    Ty* resVP = resV.data();
    if(mapA.size() != 0)
    qacc_for(isp, mapA.size() ,{
      const Long i0 = A[isp];
      for(Int iv=0;iv<nvec;iv++){
        ////resF[isp*nvec + iv] = resV[i0*nvec + iv];
        resFP[isp*nvec + iv] = resVP[iv*Nvol + i0];
      }
    });
  }

  template<typename Ty > 
  void pick_mom_from_vecs(qlat::vector_gpu<Ty >& resF, std::vector<FieldG<Ty> >& resV){
    TIMER("pick_mom_from_vecs");
    if(resV.size() == 0){resF.resize(0);return;}
    Qassert(resV[0].initialized and geo == resV[0].geo());

    vector<Ty*> resVP;resVP.resize(resV.size());
    for(unsigned int si=0;si<resV.size();si++){
      Qassert(resV[si].initialized and resV[si].multiplicity == 1);
      resVP[si] = (Ty*) get_data(resV[si]).data();
    }
    const Long Mvol = mapA.size();
    ////const Long Nvol_ = V;
    //const Long Nvol = geo.local_volume();
    Int nvec = resV.size();
    resF.resize(Mvol * nvec);
    nvec_copy = nvec;

    const Long* A = mapA.data();
    Ty* resFP = resF.data();
    ////Ty* resVP = resV.data();
    if(mapA.size() != 0)
    qacc_for(isp, mapA.size() ,{
      const Long i0 = A[isp];
      for(Int iv=0;iv<nvec;iv++){
        ////resF[isp*nvec + iv] = resV[i0*nvec + iv];
        resFP[isp*nvec + iv] = resVP[iv][i0];
      }
    });
  }

  template<typename Ty, typename Ty1 > 
  Int copy_momF_to_sf(qlat::SelectedField<Ty1  >& sf_, qlat::vector_gpu<Ty >& srcF, Int dir = 0 ){
    TIMERA("copy_momF_to_sf");
    const Long Mvol = mapA.size();
    ///const Long Mvol_ = Mvol;
    ///const Long Nvol_ = V;
    Int nvec  = 0;
    if(dir == 0){
     if(Mvol != 0){nvec = srcF.size()/Mvol;}
    }
    if(dir == 1){
      Qassert(sf_.initialized);
      //Qassert(sf.n_elems % Mvol == 0)
      //nvec = sf.n_elems / Mvol;
      //nvec = sf.n_elems;
      Qassert(Mvol == sf.n_elems);
      nvec = sf_.multiplicity;
      srcF.resize(Mvol * nvec);
      //qmessage("nvec %d, Mvol %d, elems %d \n", int(nvec), int(Mvol), int(sf.n_elems));
    }
    sum_value_mpi(nvec);
    nvec_copy = nvec;
    if(!sf_.initialized or sf_.multiplicity != nvec){sf_.init(fsel, nvec);}

    Ty* src = srcF.data();
    const Long* F = fsel_map.data();
    if(Mvol != 0)
    qacc_for(isp, Mvol, {
      ////const Long i0 = A[isp];
      const Long si = F[isp];
      Ty1* sfP = (Ty1*) qlat::get_data(sf_.get_elems(si)).data();
      Ty* reP = (Ty*) &src[isp*nvec + 0];
      if(dir == 0){for(Int iv=0;iv<nvec;iv++){sfP[iv] = reP[iv];} }
      if(dir == 1){for(Int iv=0;iv<nvec;iv++){reP[iv] = sfP[iv];} }
    });
    return nvec;
  }

  template<typename Ty, typename Ty1 > 
  Int copy_sf_to_momF(qlat::vector_gpu<Ty >& srcF, qlat::SelectedField<Ty1 >& sf_ )
  {
    return copy_momF_to_sf(sf_, srcF, 1);
  }

  template<typename Ty > 
  void write(qlat::vector_gpu<Ty >& srcF, const std::string& nameQ, const std::string& tag_ = std::string("-1"), const bool clean = false){
    TIMERA("Qlat write mdat");
    //const Long Mvol =  mapA.size();
    //printf("Long %d %ld %ld %ld \n",  qlat::get_id_node(), srcF.size() ,mapA.size(), srcF.size()%mapA.size());
    //Qassert(srcF.size() % Mvol == 0);
    if(mapA.size() !=0 ){Qassert(srcF.size() % mapA.size() == 0);}

    //qmessage("===check norm");srcF.print_norm2();

    ////const ShuffledBitSet sbs = mk_shuffled_bitset(fsel, new_size_node);

    bool append = true;if(tag_ == "-1" or clean == true){append = false;}
    // clean if not exist
    if(!does_file_exist_qar_sync_node(nameQ + "/geon-info.txt")){
      append = false;
    }
    if(append == false)
    {
      if(0 == qlat::get_id_node()){
        qlat::qremove_all(nameQ);
      }
    }
    ShuffledFieldsWriter sfw(nameQ, new_size_node, append);

    const Int nread = copy_momF_to_sf(sf, srcF);
    Qassert( nread != 0);// must read or write something
    std::string tag = ssprintf("%s.momcut%05d", tag_.c_str(), mom_cut);
    //qlat::write(sfw, tag, sf, sbs);
    ////default single precision files
    qlat::write(sfw, tag, sbs, sf);
    sfw.close();
  }

  /////template<typename Ty > 
  //Long read(qlat::SelectedField<TWOPT_TYPE  >& sf_, const Int nvec, const std::string& nameQ, const std::string& tag_){
  //  if(!sf_.initialized or sf_.n_elems != nvec){sf_.init(fsel, nvec);}
  //  ShuffledFieldsReader sfr(nameQ);
  //  std::string tag = ssprintf("%s.momcut%05d", tag_.c_str(), mom_cut);
  //  Long total_bytes = qlat::read(sfr, tag, sf_, fsel);
  //  ////if(total_bytes ==  0){srcF.resize(0); return  total_bytes;}
  //  ////copy_sf_to_momF(srcF);
  //  return total_bytes;
  //}

  template<typename Ty > 
  void shift_t(qlat::vector_gpu<Ty >& s1, qlat::vector_gpu<Ty >& s0, const Coordinate& shift_){
    TIMERA("shift_t");
    ////Qassert(t0 < nv[3]);

    Coordinate shift(0, 0, 0, 0);
    for(Int i=0;i<4;i++){Qassert(shift[i] < nv[i]); shift[i] = -1 * shift_[i];}
    //shift[3] = (nv[3] - t0 + nv[3])%nv[3];
    //shift[3] = -t0;
    Int nread = copy_momF_to_sf(sf, s0);
    Qassert( nread != 0);// must read or write something

    qlat::field_shift(sf_1, fsel_1, sf, fsel, shift);
        nread = copy_sf_to_momF(s1, sf_1);
    Qassert( nread != 0);// must read or write something

    //Qassert(nvec_copy != 0);
    //Ty* src = s1.data();
    //const Long* A = mapA.data();
    //const Long  Mvol = mapA.size();
    //qacc_for(isp, Mvol, {
    //  Long i0 = A[isp];
    //  ////const Long i0 = A[isp];
    //  const Long si = fsel_1.f_local_idx.get_elem_offset(i0);
    //  TWOPT_TYPE* sfP = (TWOPT_TYPE*) qlat::get_data(sf_1.get_elems(si)).data();
    //  Ty* reP = (Ty*) &src[isp*nvec_copy + 0];
    //  ////if(dir == 0){for(Int iv=0;iv<nvec_copy;iv++){sfP[iv] = reP[iv];} }
    //  ////if(dir == 1){for(Int iv=0;iv<nvec_copy;iv++){reP[iv] = sfP[iv];} }
    //  for(Int iv=0;iv<nvec_copy;iv++){reP[iv] = sfP[iv];}
    //});
  }

  template<typename Ty > 
  void shift_t(qlat::vector_gpu<Ty >& s1, qlat::vector_gpu<Ty >& s0, const Int t0){
    shift_t(s1, s0, Coordinate(0,0,0, t0));
  }

  inline void get_fn_list(const std::string& nameQ)
  {
    if(file_name != nameQ){
      ShuffledFieldsReader sfr(nameQ);
      fn_list = qlat::list_fields(sfr);
      file_name = nameQ;
      sfr.close();
    }
  }

  inline bool check_fn(const std::string& nameQ, const std::string& tag){
    get_fn_list(nameQ);
    bool find = false;
    for(unsigned int fi=0;fi<fn_list.size();fi++)
    {
      if(fn_list[fi] == tag){
        find = true;
        break ;
      }
    }
    return find;
  }

  inline bool check_fn_momcut(const std::string& nameQ, const std::string& tag_){
    std::string tag = ssprintf("%s.momcut%05d", tag_.c_str(), mom_cut);
    return check_fn(nameQ, tag);
  }

  // change interface to return number of vectors read
  template<typename Ty > 
  Int read(qlat::vector_gpu<Ty >& srcF, const std::string& nameQ, const std::string& tag){
    TIMERA("Qlat read mdat");
    ShuffledFieldsReader sfr(nameQ);
    if(!check_fn(nameQ, tag)){
      qmessage("File %s , tag %s not found! \n", nameQ.c_str(), tag.c_str());MPI_Barrier(MPI_COMM_WORLD);
      fflush(stdout);qlat::end();abort();};
    Long total_bytes = qlat::read(sfr, tag, sf, fsel);
    sfr.close();
    if(total_bytes ==  0){srcF.resize(0); return  total_bytes;}
    const Int nread = copy_sf_to_momF(srcF, sf);
    Qassert( nread != 0);// must read or write something
    //return total_bytes;
    return nread;
  }

  template<typename Ty > 
  Int read_momcut(qlat::vector_gpu<Ty >& srcF, const std::string& nameQ, const std::string& tag_){
    std::string tag = ssprintf("%s.momcut%05d", tag_.c_str(), mom_cut);
    return read(srcF, nameQ, tag);
  }

  inline void update_mapB_mom_off(const Coordinate& mom_off)
  {
    if(mom_off_buf != mom_off or mapB_size == -1)
    {
      get_mom_single_nodeB(mapA, mapB_buf, geo, mom_cut, mom_off);
      mom_off_buf = mom_off;
      mapB_size = mapB_buf.size();
    }
  }

  /////mom_cut should be consistent with your production indicated in the saving file
  momentum_dat(const Geometry& geo_, const Int mom_cut_, const std::vector< Coordinate>& mom_off = std::vector< Coordinate>(0)){
    TIMERA("momentum_dat");
    geo = geo_;
    mom_cut = mom_cut_;

    qlat::vector<Int > Nv, mv;
    geo_to_nv(geo, nv, Nv, mv);
    ////nt = nv[3];
    ////Nvol = geo.local_volume();

    get_mom_single_nodeA(mapA, geo, mom_cut, mom_off);
    const Long Mvol = mapA.size();
    mapB_size = -1;

    ///Mvol = mapA.size();
    ////const Int mc = mom_cut*2 + 1;
    //Mvol = nt * mc*mc*mc;
    ////read_field_selection(fsel, nameQF, -1);

    ////TODO
    //const Int Nmpi = qlat::get_num_node();
    //const Int rank  = qlat::get_id_node();
    //std::vector<Int > pconf_vec;
    //pconf_vec.resize(Nmpi * Mvol * 4);
    //qthread_for(ai , Mvol , {
    //  const Long isp = mapA[ai];
    //  const Coordinate xl  = geo.coordinate_from_index(isp);
    //  const Coordinate xg  = geo.coordinate_g_from_l(xl);
    //  Int* res = &pconf_vec[(rank * Mvol + ai) * 4 + 0];
    //  for(unsigned int i = 0; i < 4 ; i++){res[i] = xg[i];}
    //});
    //sum_all_size(pconf_vec.data(), pconf_vec.size(), 0);

    //PointSelection pconf;
    //pconf.resize(Nmpi * Mvol);
    //qthread_for(ai , Nmpi * Mvol , {
    //  Int* src = &pconf_vec[ai * 4 + 0];
    //  for(unsigned int i = 0; i < 4 ; i++){pconf[ai][i] = src[i];}
    //});

    PointsSelection pconf;
    pconf.init(geo.total_site(), Mvol);
    if(Mvol != 0)
    qthread_for(ai , Mvol , {
      const Long isp = mapA[ai];
      const Coordinate xl  = geo.coordinate_from_index(isp);
      const Coordinate xg  = geo.coordinate_g_from_l(xl);
      pconf[ai] = xg;
    });

    ////TODO
    Coordinate total_site = Coordinate(nv[0], nv[1], nv[2], nv[3]);
    Long n_per_tslice =  0;qlat::RngState rs(321);
    {
      TIMERA("qlat::set_field_selection");
      qlat::set_field_selection(fsel, total_site, n_per_tslice, rs, pconf);
    }

    //fsel.init();fsel.f_rank.init(geo);
    //add_field_selection(fsel.f_rank, pconf);
    ////////update_field_selection(fsel);
    //update_field_selection(fsel, 0);

    fsel_map.resize(Mvol);
    if(Mvol != 0)
    qthread_for(isp, Mvol,{
      Long i0 = mapA[isp];
      //Long si = fsel.f_local_idx.get_elems_const(i0)[0];
      const Long si = fsel.f_local_idx.get_elem_offset(i0);
      //Qassert(si != -1);
      fsel_map[isp] = si;
    });

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

    sbs = mk_shuffled_bitset(fsel, new_size_node);

    cur_pos = Coordinate(0, 0, 0, 0);
    nvec_copy = 0;

    file_name = "NONE";
  }

  inline void get_mom_pos(const std::vector<Coordinate >& mom_pick, std::vector<bool >& Bsite, std::vector<Long >& Zsite)
  {
    const Long Mvol = mapA.size();
    Long* A = mapA.data();
    const Geometry& geo_ = geo;
    Bsite.resize(mom_pick.size());
    Zsite.resize(mom_pick.size());
    for(unsigned int bi=0;bi<Bsite.size();bi++){Bsite[bi] = false; Zsite[bi] = 0;}
    /////Bool is not write thread safe ...... be careful.....
    if(Mvol != 0)
    for(unsigned int momi=0;momi<mom_pick.size();momi++){
    qthread_for(isp, Mvol, {
      const Long ilocal = A[isp];
      const Coordinate xl   = geo_.coordinate_from_index(ilocal);
      const Coordinate mom  = geo_.coordinate_g_from_l(xl);
      if(mom == mom_pick[momi]){
        Bsite[momi] = true; Zsite[momi] = isp;
      }
    });}
  }

  /////calculate source phases with coordinate `shift`
  /////src phases e^{-i p ( x - y)}
  template<typename Ty > 
  void update_phases(const Coordinate& src , const Coordinate& shift = Coordinate(0,0,0,0), const Int sign = -1)
  {
    TIMERA("update_phases");
    const Long Mvol = mapA.size();
    Qassert(sizeof(Ty)%sizeof(Complexq) == 0);
    const Int fac = sizeof(Ty)/sizeof(Complexq);
    if(cur_pos == src and Long(phases.size()) == Mvol * fac and cur_shift == shift){return ;}
    cur_shift = shift;
    phases.resize(Mvol * fac);Ty* resP = (Ty*) phases.data();
    qlat::vector<Int >& Lat = nv;
    Long* A = mapA.data();
    const Geometry& geo_ = geo;
    if(Mvol != 0)
    qacc_for(isp, Mvol, {
      const Long ilocal = A[isp];
      const Coordinate xl  = geo_.coordinate_from_index(ilocal);
      const Coordinate mom  = geo_.coordinate_g_from_l(xl);
      double v0 = 0.0; 
      for(Int i=0;i<3;i++){v0 += (2.0* QLAT_PI_LOCAL * src[i] * ((mom[i] + shift[i] + Lat[i])%(Lat[i]))/Lat[i]);}
      resP[isp] = Ty(std::cos(v0), sign * std::sin(v0));
    });
    cur_pos = src;
  }

  /////src phases e^{-i p ( x - y)}
  template<typename Ty > 
  void apply_src_phases(qlat::vector_gpu<Ty >& vec, const Coordinate& src , const Coordinate& shift = Coordinate(0,0,0,0) , const Int sign = -1)
  {
    TIMERA("apply_src_phases");
    const Long Mvol = mapA.size();
    if(Mvol != 0){Qassert(vec.size()%Mvol == 0);}
    ////const Long nvec = src.size() / Mvol;
    Qassert(nvec_copy != 0);
    const Long nvec = nvec_copy;

    update_phases<Ty >(src, shift, sign);

    Ty* momP = (Ty*) phases.data();
    Ty* resV = vec.data();
    if(Mvol != 0)
    qacc_for(isp, Mvol, {
      Ty ph = momP[isp];
      Ty* tmp = &resV[isp*nvec + 0];
      for(Int iv=0;iv<nvec;iv++){
        //tmp[iv] = tmp[iv] * ph;
        tmp[iv] *= ph;
      }
    });
  }

};

#undef TWOPT_TYPE
  
template<typename Ty >
void fft_local_to_global(qlat::vector_gpu<Ty >& FG, qlat::vector_gpu<Ty >& FL, momentum_dat& mdat, const Coordinate& mom_off = Coordinate(0, 0, 0, 0))
{
  TIMERA("fft_local_to_global");
  //std::vector<Int > nv, Nv, mv;
  //geo_to_nv(mdat.geo, nv, Nv, mv);
  const Int    mc = mdat.mom_cut*2 + 1;

  Long nvec = 0;
  if(mdat.mapA.size() != 0){ nvec = FL.size()/mdat.mapA.size();}
  sum_value_mpi(nvec);

  const Long Mvol = mdat.nv[3]*mc*mc*mc;
  Qassert(nvec > 0);
  FG.resize(nvec * Mvol);FG.set_zero();

  mdat.update_mapB_mom_off(mom_off);
  const Long* PmapB = (Long*) qlat::get_data(mdat.mapB_buf).data();
  const Long NmapB  = mdat.mapB_buf.size();
  //qmessage("===check norm");FL.print_norm2();
  ////Qassert(NmapB <= Mvol);

  Ty* PFG = FG.data();
  Ty* PFL = FL.data();

  if(NmapB != 0)
  qacc_for(isp, NmapB, {
    const Long i1 = PmapB[isp];
    /////mapA may have more data than mapB due to mom_off
    if(i1 >= 0)
    {
      for(Long iv=0;iv<nvec;iv++){PFG[iv*Mvol + i1] = PFL[isp*nvec + iv];}
    }
  });

  sum_all_size(FG.data(), FG.size(), FG.GPU);
}

template<typename Ty, typename Ta>
void copy_sparse_fields(qlat::SelectedField<Ty >& res, qlat::SelectedField<Ta >& src, const Int Ndc, Int nr=0, const Int ns=0)
{
  Qassert(src.initialized and res.initialized);
  Qassert(src.field.size() % Ndc == 0);
  Qassert(res.field.size() % Ndc == 0);
  Qassert(src.multiplicity % Ndc == 0);
  Qassert(res.multiplicity % Ndc == 0);

  const Long Ndata = src.field.size() / src.multiplicity;
  ///// printf("=== %8d %8d \n", int(src.field.size() / src.multiplicity), int(res.field.size() / res.multiplicity));
  Qassert(Ndata == (res.field.size() / res.multiplicity));

  const Int Nr = res.field.size() / (Ndata * Ndc);
  const Int Ns = src.field.size() / (Ndata * Ndc);
  Qassert(Nr > nr and Ns > ns);

  Ty* pr = (Ty*) qlat::get_data(res.field).data();
  Ta* ps = (Ta*) qlat::get_data(src.field).data();

  qacc_for(idx, Ndata, {
    for(Int i=0;i<Ndc;i++)
    {
      pr[(idx*Nr + nr)*Ndc+i] = ps[(idx*Ns + ns)*Ndc + i];
    }
  })
}

}

#endif
