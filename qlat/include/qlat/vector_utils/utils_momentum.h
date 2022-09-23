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
///////#include "../kentucky/utils_lat_exchanger.h"

namespace qlat{

//template<typename Ty, int civ>
//void FieldM_to_vector_gpu(qlat::vector_gpu<Ty >& res, std::vector<qlat::FieldM<Ty , civ> >& src, int dir = 0){
//  qassert(src.size() > 0);qassert(src[0].initialized);
//  const Geometry geo = src[0].geo();
//  const long V = geo.local_volume();
//  const int nvec = src.size();
//  if(dir == 0){res.resize(V * nvec * civ);}
//  if(dir == 1){src.resize(nvec); for(unsigned int i=0;i<src.size();i++){src[i].init(geo);}}
//  Ty* rb = res.data();
//  for(int iv=0;iv<nvec;iv++){
//    Ty* ra = (Ty*) qlat::get_data(src[iv]).data();
//    if(dir ==  0)cpy_data_thread(&rb[iv * V*civ], ra, V*civ, 1, false);
//    if(dir ==  1)cpy_data_thread(ra, &rb[iv * V*civ], V*civ, 1, false);
//  }
//  qacc_barrier(dummy);
//}
//
//template<typename Ty, int civ>
//void vector_gpu_to_FieldM(std::vector<qlat::FieldM<Ty , civ> >& src, qlat::vector_gpu<Ty >& res){
//  FieldM_to_vector_gpu(res, src, 1);
//}

inline void get_mom_single_node(qlat::vector_acc<long >& mapA, qlat::vector_acc<long >& mapB,
    const Geometry& geo, const int mom_cut)
{
  qlat::vector_acc<int > nv,Nv,mv;
  geo_to_nv(geo, nv, Nv, mv);

  mapB.resize(geo.local_volume());long* Bi = mapB.data();
  const int mc = mom_cut*2 + 1;

  qacc_for(isp,  geo.local_volume(),{
    Coordinate xl  = geo.coordinate_from_index(isp);
    Coordinate xg  = geo.coordinate_g_from_l(xl);
    bool flag = true;
    Coordinate mom = xg;
    for(int i=0;i<3;i++){
      if(xg[i] > nv[i]/2){
        mom[i] = (nv[i] - xg[i]);
        if(mom[i] >  mom_cut){flag = false;}
        else{mom[i] = mom_cut*2 + 1 - mom[i];}
      }else{if(mom[i] > mom_cut){flag = false;}}
    }
    if(flag == true){ Bi[isp] = ((mom[3]*mc + mom[2])*mc + mom[1])*mc+ mom[0]; }
    if(flag == false){Bi[isp] = -1;}
  });

  std::vector<long > A0;std::vector<long > B0;
  for(long isp=0;isp < mapB.size();isp++)
  {
    if(mapB[isp] >= 0){
      A0.push_back(  isp );B0.push_back(mapB[isp]);
    }
  }
  mapA.resize(A0.size());mapB.resize(A0.size());
  ////long* A = mapA.data(); long* B  = mapB.data();
  ////long* Av = A0.data();  long* Bv = B0.data();
  ////copy back to mapB
  qthread_for(isp,  mapB.size(),{
    mapA[isp] = A0[isp]; mapB[isp] = B0[isp];
  });
}

#define TWOPT_TYPE  qlat::ComplexF

/////write date in float prec 
struct momentum_dat{
  ////int nx,ny,nz,nt;
  Geometry geo;
  int mom_cut;
  ////long Nvol;
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
  qlat::vector_acc<long > mapA, mapB;
  qlat::vector_acc<long > fsel_map;

  qlat::vector_acc<int > nv;
  Coordinate  cur_pos;
  qlat::vector_gpu<Complexq > phases;
  Coordinate cur_shift;
  /////int mc;
  ///int Mvol;

  int nvec_copy;

  template<typename Ty > 
  void pick_mom_from_vecs(qlat::vector_gpu<Ty >& resF, qlat::vector_gpu<Ty >& resV){
    TIMER("pick_mom_from_vecs");
    const long Mvol = mapA.size();
    ////const long Nvol_ = V;
    const long Nvol = geo.local_volume();
    int nvec = resV.size()/Nvol;
    resF.resize(Mvol * nvec);

    sum_value_mpi(nvec);
    nvec_copy = nvec;

    const long* A = mapA.data();
    Ty* resFP = resF.data();
    Ty* resVP = resV.data();
    if(mapA.size() != 0)
    qacc_for(isp, mapA.size() ,{
      const long i0 = A[isp];
      for(int iv=0;iv<nvec;iv++){
        ////resF[isp*nvec + iv] = resV[i0*nvec + iv];
        resFP[isp*nvec + iv] = resVP[iv*Nvol + i0];
      }
    });
  }

  template<typename Ty, typename Ty1 > 
  void copy_momF_to_sf(qlat::SelectedField<Ty1  >& sf_, qlat::vector_gpu<Ty >& srcF, int dir = 0 ){
    TIMERA("copy_momF_to_sf");
    const long Mvol = mapA.size();
    ///const long Mvol_ = Mvol;
    ///const long Nvol_ = V;
    int nvec  = 0;
    if(dir == 0){
     if(Mvol != 0){nvec = srcF.size()/Mvol;}
    }
    if(dir == 1){
      qassert(sf_.initialized);
      //qassert(sf.n_elems % Mvol == 0)
      //nvec = sf.n_elems / Mvol;
      //nvec = sf.n_elems;
      nvec = sf_.geo().multiplicity;
      srcF.resize(Mvol * nvec);
    }
    sum_value_mpi(nvec);
    nvec_copy = nvec;
    if(!sf_.initialized or sf_.geo().multiplicity != nvec){sf_.init(fsel, nvec);}

    Ty* src = srcF.data();
    const long* F = fsel_map.data();
    if(Mvol != 0)
    qacc_for(isp, Mvol, {
      ////const long i0 = A[isp];
      const long si = F[isp];
      Ty1* sfP = (Ty1*) qlat::get_data(sf_.get_elems(si)).data();
      Ty* reP = (Ty*) &src[isp*nvec + 0];
      if(dir == 0){for(int iv=0;iv<nvec;iv++){sfP[iv] = reP[iv];} }
      if(dir == 1){for(int iv=0;iv<nvec;iv++){reP[iv] = sfP[iv];} }
    });
  }

  template<typename Ty, typename Ty1 > 
  void copy_sf_to_momF(qlat::vector_gpu<Ty >& srcF, qlat::SelectedField<Ty1 >& sf_ )
  {copy_momF_to_sf(sf_, srcF, 1);}

  template<typename Ty > 
  void write(qlat::vector_gpu<Ty >& srcF, const std::string& nameQ, const std::string& tag_ = std::string("-1"), const bool clean = false){
    TIMERA("Qlat write mdat");
    //const long Mvol =  mapA.size();
    //printf("long %d %ld %ld %ld \n",  qlat::get_id_node(), srcF.size() ,mapA.size(), srcF.size()%mapA.size());
    //qassert(srcF.size() % Mvol == 0);
    if(mapA.size() !=0 ){qassert(srcF.size() % mapA.size() == 0);}

    Coordinate new_size_node = Coordinate(1, 1, 2, 4);
    const ShuffledBitSet sbs = mk_shuffled_bitset(fsel, new_size_node);

    bool append = true;if(tag_ == "-1" or clean == true){append = false;}
    ShuffledFieldsWriter sfw(nameQ, new_size_node, append);

    copy_momF_to_sf(sf, srcF);
    std::string tag = ssprintf("%s.momcut%05d", tag_.c_str(), mom_cut);
    qlat::write(sfw, tag, sf, sbs);
  }

  /////template<typename Ty > 
  //long read(qlat::SelectedField<TWOPT_TYPE  >& sf_, const int nvec, const std::string& nameQ, const std::string& tag_){
  //  if(!sf_.initialized or sf_.n_elems != nvec){sf_.init(fsel, nvec);}
  //  ShuffledFieldsReader sfr(nameQ);
  //  std::string tag = ssprintf("%s.momcut%05d", tag_.c_str(), mom_cut);
  //  long total_bytes = qlat::read(sfr, tag, sf_, fsel);
  //  ////if(total_bytes ==  0){srcF.resize(0); return  total_bytes;}
  //  ////copy_sf_to_momF(srcF);
  //  return total_bytes;
  //}

  template<typename Ty > 
  void shift_t(qlat::vector_gpu<Ty >& s1, qlat::vector_gpu<Ty >& s0, const Coordinate& shift_){
    TIMERA("shift_t");
    ////qassert(t0 < nv[3]);

    Coordinate shift(0, 0, 0, 0);
    for(int i=0;i<4;i++){qassert(shift[i] < nv[i]); shift[i] = -1 * shift_[i];}
    //shift[3] = (nv[3] - t0 + nv[3])%nv[3];
    //shift[3] = -t0;
    copy_momF_to_sf(sf, s0);
    qlat::field_shift(sf_1, fsel_1, sf, fsel, shift);
    copy_sf_to_momF(s1, sf_1);

    //qassert(nvec_copy != 0);
    //Ty* src = s1.data();
    //const long* A = mapA.data();
    //const long  Mvol = mapA.size();
    //qacc_for(isp, Mvol, {
    //  long i0 = A[isp];
    //  ////const long i0 = A[isp];
    //  const long si = fsel_1.f_local_idx.get_elem(i0);
    //  TWOPT_TYPE* sfP = (TWOPT_TYPE*) qlat::get_data(sf_1.get_elems(si)).data();
    //  Ty* reP = (Ty*) &src[isp*nvec_copy + 0];
    //  ////if(dir == 0){for(int iv=0;iv<nvec_copy;iv++){sfP[iv] = reP[iv];} }
    //  ////if(dir == 1){for(int iv=0;iv<nvec_copy;iv++){reP[iv] = sfP[iv];} }
    //  for(int iv=0;iv<nvec_copy;iv++){reP[iv] = sfP[iv];}
    //});
  }

  template<typename Ty > 
  void shift_t(qlat::vector_gpu<Ty >& s1, qlat::vector_gpu<Ty >& s0, const int t0){
    shift_t(s1, s0, Coordinate(0,0,0, t0));
  }

  inline void get_fn_list(const std::string& nameQ)
  {
    if(file_name != nameQ){
      ShuffledFieldsReader sfr(nameQ);
      fn_list = qlat::list_fields(sfr);
      file_name = nameQ;
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

  template<typename Ty > 
  long read(qlat::vector_gpu<Ty >& srcF, const std::string& nameQ, const std::string& tag){
    TIMERA("Qlat read mdat");
    ShuffledFieldsReader sfr(nameQ);
    if(!check_fn(nameQ, tag)){
      print0("File %s , tag %s not found! \n", nameQ.c_str(), tag.c_str());MPI_Barrier(MPI_COMM_WORLD);
      fflush(stdout);qlat::end();abort();};
    long total_bytes = qlat::read(sfr, tag, sf, fsel);
    if(total_bytes ==  0){srcF.resize(0); return  total_bytes;}
    copy_sf_to_momF(srcF, sf);
    return total_bytes;
  }

  template<typename Ty > 
  long read_momcut(qlat::vector_gpu<Ty >& srcF, const std::string& nameQ, const std::string& tag_){
    std::string tag = ssprintf("%s.momcut%05d", tag_.c_str(), mom_cut);
    return read(srcF, nameQ, tag);
  }

  /////mom_cut should be consistent with your production indicated in the saving file
  momentum_dat(const Geometry& geo_, const int mom_cut_){
    geo = geo_;
    mom_cut = mom_cut_;

    qlat::vector_acc<int > Nv, mv;
    geo_to_nv(geo, nv, Nv, mv);
    ////nt = nv[3];
    ////Nvol = geo.local_volume();

    get_mom_single_node(mapA, mapB, geo, mom_cut);
    const long Mvol = mapA.size();

    ///Mvol = mapA.size();
    ////const int mc = mom_cut*2 + 1;
    //Mvol = nt * mc*mc*mc;
    ////read_field_selection(fsel, nameQF, -1);

    ////TODO
    //const int Nmpi = qlat::get_num_node();
    //const int rank  = qlat::get_id_node();
    //std::vector<int > pconf_vec;
    //pconf_vec.resize(Nmpi * Mvol * 4);
    //qthread_for(ai , Mvol , {
    //  const long isp = mapA[ai];
    //  const Coordinate xl  = geo.coordinate_from_index(isp);
    //  const Coordinate xg  = geo.coordinate_g_from_l(xl);
    //  int* res = &pconf_vec[(rank * Mvol + ai) * 4 + 0];
    //  for(unsigned int i = 0; i < 4 ; i++){res[i] = xg[i];}
    //});
    //sum_all_size(pconf_vec.data(), pconf_vec.size(), 0);

    //PointSelection pconf;
    //pconf.resize(Nmpi * Mvol);
    //qthread_for(ai , Nmpi * Mvol , {
    //  int* src = &pconf_vec[ai * 4 + 0];
    //  for(unsigned int i = 0; i < 4 ; i++){pconf[ai][i] = src[i];}
    //});

    PointSelection pconf;pconf.resize(Mvol);
    if(Mvol != 0)
    qthread_for(ai , Mvol , {
      const long isp = mapA[ai];
      const Coordinate xl  = geo.coordinate_from_index(isp);
      const Coordinate xg  = geo.coordinate_g_from_l(xl);
      pconf[ai] = xg;
    });

    ////TODO
    Coordinate total_site = Coordinate(nv[0], nv[1], nv[2], nv[3]);
    long n_per_tslice =  0;qlat::RngState rs(321);
    qlat::set_field_selection(fsel, total_site, n_per_tslice, rs, pconf);

    //fsel.init();fsel.f_rank.init(geo);
    //add_field_selection(fsel.f_rank, pconf);
    ////////update_field_selection(fsel);
    //update_field_selection(fsel, 0);

    fsel_map.resize(Mvol);
    if(Mvol != 0)
    qthread_for(isp, Mvol,{
      long i0 = mapA[isp];
      //long si = fsel.f_local_idx.get_elems_const(i0)[0];
      const long si = fsel.f_local_idx.get_elem(i0);
      qassert(si != -1);
      fsel_map[isp] = si;
    });

    cur_pos = Coordinate(0, 0, 0, 0);
    nvec_copy = 0;

    file_name = "NONE";
  }

  inline void get_mom_pos(const std::vector<Coordinate >& mom_pick, std::vector<bool >& Bsite, std::vector<long >& Zsite)
  {
    const long Mvol = mapA.size();
    long* A = mapA.data();
    const Geometry& geo_ = geo;
    Bsite.resize(mom_pick.size());
    Zsite.resize(mom_pick.size());
    for(unsigned int bi=0;bi<Bsite.size();bi++){Bsite[bi] = false; Zsite[bi] = 0;}
    /////Bool is not write thread safe ...... be careful.....
    if(Mvol != 0)
    for(unsigned int momi=0;momi<mom_pick.size();momi++){
    qthread_for(isp, Mvol, {
      const long ilocal = A[isp];
      const Coordinate xl   = geo_.coordinate_from_index(ilocal);
      const Coordinate mom  = geo_.coordinate_g_from_l(xl);
      if(mom == mom_pick[momi]){
        Bsite[momi] = true; Zsite[momi] = isp;
      }
    });}
  }

  template<typename Ty > 
  void update_phases(const Coordinate& src , const Coordinate& shift = Coordinate(0,0,0,0), const int sign = -1)
  {
    TIMERA("update_phases");
    const long Mvol = mapA.size();
    qassert(sizeof(Ty)%sizeof(Complexq) == 0);
    const int fac = sizeof(Ty)/sizeof(Complexq);
    if(cur_pos == src and long(phases.size()) == Mvol * fac and cur_shift == shift){return ;}
    cur_shift = shift;
    phases.resize(Mvol * fac);Ty* resP = (Ty*) phases.data();
    qlat::vector_acc<int >& Lat = nv;
    long* A = mapA.data();
    const Geometry& geo_ = geo;
    if(Mvol != 0)
    qacc_for(isp, Mvol, {
      const long ilocal = A[isp];
      const Coordinate xl  = geo_.coordinate_from_index(ilocal);
      const Coordinate mom  = geo_.coordinate_g_from_l(xl);
      double v0 = 0.0; 
      for(int i=0;i<3;i++){v0 += (2.0*PI * src[i] * ((mom[i] + shift[i] + Lat[i])%(Lat[i]))/Lat[i]);}
      resP[isp] = Ty(std::cos(v0), sign * std::sin(v0));
    });
    cur_pos = src;
  }

  template<typename Ty > 
  void apply_src_phases(qlat::vector_gpu<Ty >& vec, const Coordinate& src , const Coordinate& shift = Coordinate(0,0,0,0) , const int sign = -1)
  {
    TIMERA("apply_src_phases");
    const long Mvol = mapA.size();
    if(Mvol != 0){qassert(vec.size()%Mvol == 0);}
    ////const long nvec = src.size() / Mvol;
    qassert(nvec_copy != 0);
    const long nvec = nvec_copy;

    update_phases<Ty >(src, shift, sign);

    Ty* momP = (Ty*) phases.data();
    Ty* resV = vec.data();
    if(Mvol != 0)
    qacc_for(isp, Mvol, {
      Ty ph = momP[isp];
      Ty* tmp = &resV[isp*nvec + 0];
      for(int iv=0;iv<nvec;iv++){
        //tmp[iv] = tmp[iv] * ph;
        tmp[iv] *= ph;
      }
    });
  }

};

#undef TWOPT_TYPE
  
template<typename Ty >
void fft_local_to_global(qlat::vector_gpu<Ty >& FG, qlat::vector_gpu<Ty >& FL, momentum_dat& mdat)
{
  TIMERA("fft_local_to_global");
  //std::vector<int > nv, Nv, mv;
  //geo_to_nv(mdat.geo, nv, Nv, mv);
  const int    mc = mdat.mom_cut*2 + 1;

  long nvec = 0;
  if(mdat.mapA.size() != 0){ nvec = FL.size()/mdat.mapA.size();}
  sum_value_mpi(nvec);

  const long Mvol = mdat.nv[3]*mc*mc*mc;
  qassert(mdat.mapA.size() <= Mvol);
  qassert(nvec > 0);
  FG.resize(nvec * Mvol);FG.set_zero();

  const long* PmapB = (long*) qlat::get_data(mdat.mapB).data();
  const long NmapB = mdat.mapB.size();

  Ty* PFG = FG.data();
  Ty* PFL = FL.data();

  if(NmapB != 0)
  qacc_for(isp, NmapB, {
    const long i1 = PmapB[isp];
    for(long iv=0;iv<nvec;iv++){PFG[iv*Mvol + i1] = PFL[isp*nvec + iv];}
  });

  sum_all_size(FG.data(), FG.size(), 1);
}



}


#endif
