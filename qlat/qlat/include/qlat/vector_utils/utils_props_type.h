// utils_props_type.h
// Gen Wang
// Mar. 2024

#ifndef UTILS_PROPS_TYPE_H
#define UTILS_PROPS_TYPE_H

#pragma once
#include "general_funs.h"
#include "utils_fft_desc.h"
#include "utils_field_operations.h"
#include "utils_field_gpu.h"

////////dim 12*12 --> Nt --> Nxyz
///////only memory size and geo are used, donot use others
//#define qprop   qlat::FieldM<Complexq, 12*12>
//#define qpropT  qlat::FieldM<Ty, 12*12>

////std::vector<qlat::FermionField4dT<Tc > >, gwu format, need to rotate bases

namespace qlat{

template<typename Td, typename Tc>
void prop4d_to_Fermion(std::vector<qlat::FermionField4dT<Tc > > &buf, Propagator4dT<Td>& prop, Int dir=1){

  if(dir==1){buf.resize(0);buf.resize(12);for(Int iv=0;iv<12;iv++){
    if(!buf[iv].initialized){buf[iv].init(prop.geo());}
  }}
  if(dir==0){Qassert(buf.size() == 12);if(!prop.initialized){prop.init(buf[0].geo());}}

  #pragma omp parallel for
  for (Long index = 0; index < prop.geo().local_volume(); ++index)
  {
    qlat::WilsonMatrixT<Td>& src =  prop.get_elem_offset(index);
    for(Int d0=0;d0<12;d0++)
    {
      ////v0(s*3 + c0, ga.ind[d0]*3 + c1)
      qlat::ComplexT<Tc>* res = (qlat::ComplexT<Tc>*)&(buf[d0].get_elem_offset(index));
      for(Int d1=0;d1<12;d1++)
      {
        if(dir==1){res[d1] = src(d1, d0);}
        if(dir==0){src(d1, d0) = res[d1];}
      }

    }
  }

}

template<typename Td, typename Tc>
void Fermion_to_prop4d(Propagator4dT<Td>& prop, std::vector<qlat::FermionField4dT<Tc > > &buf){
  Qassert(buf.size() == 12);
  prop4d_to_Fermion(buf, prop, 0);
}

/////V -- 12a x 12b   to   12b x 12a -- V
template<class Ty, typename Td>
void prop4d_to_qprop(qpropT& res, Propagator4dT<Td>& src, Int dir = 1){
  TIMERA("prop4d_to_qprop");
  if(dir == 1){
    Qassert(src.initialized);
    if(!res.initialized or res.geo() != src.geo()){
      res.init();
      res.init(src.geo());
    }
  }
  if(dir == 0){
    Qassert(res.initialized);
    if(!src.initialized or src.geo() != res.geo()){
      src.init();
      src.init(res.geo());
    }
  }

  Qassert(res.geo().is_only_local and src.geo().is_only_local);
  Long sizeF = src.geo().local_volume();

  move_index mv_civ;
  qlat::ComplexT<Td>* ps; Ty* pt;
  ps = (qlat::ComplexT<Td>* ) qlat::get_data(src).data();
  pt = (Ty*) qlat::get_data(res).data();

  ////V x 12 a x 12 b to 12b x 12a x V
  if(dir == 1){
    qthread_for(isp, Long(sizeF),{
      qlat::ComplexT<Td> buf[12*12];for(unsigned int i=0;i<12*12;i++){buf[i] = ps[isp*12*12 + i];}
      for(unsigned int d0=0;d0<12;d0++)
      for(unsigned int d1=0;d1<12;d1++)
      {
        pt[(isp*12+d0)*12+d1] = buf[d1*12+d0];
      }
    });
    mv_civ.move_civ_out(pt, pt, 1, sizeF, 12*12, 1, false);
  }
  ////12 a x 12 b x V to V x 12b x 12a
  if(dir == 0){
    mv_civ.move_civ_in(pt, pt, 1, 12*12, sizeF, 1, false);
    qthread_for(isp, Long(sizeF),{
      qlat::ComplexT<Td> buf[12*12];for(unsigned int i=0;i<12*12;i++){buf[i] = pt[isp*12*12 + i];}
      for(unsigned int d0=0;d0<12;d0++)
      for(unsigned int d1=0;d1<12;d1++)
      {
        ps[(isp*12+d0)*12+d1] = buf[d1*12+d0];
      }
    });
  }
}

template<typename Td, typename Ty>
void qprop_to_prop4d(Propagator4dT<Td>& res, qpropT& src){
  prop4d_to_qprop(src, res, 0);
}

/////V -- 12a x 12b   to   12b x 12a -- V
template<class Ty, typename Td>
void prop4d_to_fieldG(FieldG<Ty >& res, Propagator4dT<Td>& src, Int dir = 1){
  TIMERA("prop4d_to_fieldG");
  Qassert(src.initialized and res.initialized);
  if(dir == 0){Qassert(res.mem_order == QLAT_OUTTER);}

  qlat::ComplexT<Td>* srcP = (qlat::ComplexT<Td>* ) qlat::get_data(src).data();
  Ty* resP = (qlat::ComplexT<Td>* ) qlat::get_data(res).data();;

  Long sizeF = src.geo().local_volume();
  Qassert(res.field_gpu.GPU != QMCPU);
  Qassert(res.geo().is_only_local and src.geo().is_only_local);

  if(dir == 1){
    qacc_for(isp, Long(sizeF),{
      qlat::ComplexT<Td> buf[12*12];
      for(unsigned int i=0;i<12*12;i++){buf[i] = srcP[isp*12*12 + i];}
      for(unsigned int d0=0;d0<12;d0++)
      for(unsigned int d1=0;d1<12;d1++)
      {
        resP[(isp*12+d0)*12+d1] = buf[d1*12+d0];
      }
    });
    res.mem_order = QLAT_DEFAULT;
    switch_orders(res, QLAT_OUTTER);
    //move_index mv_civ;
    //mv_civ.move_civ_out(resP, resP, 1, sizeF, 12*12, 1, false);
    //res.mem_order = QLAT_OUTTER;
  }

  if(dir == 0){
    cpy_GPU(srcP, resP, sizeF*12*12, 1, 1);
    move_index mv_civ;
    mv_civ.move_civ_in(srcP, srcP, 1, 12*12, sizeF, 1);
    qacc_for(isp, Long(sizeF),{
      qlat::ComplexT<Td> buf[12*12];
      for(unsigned int i=0;i<12*12;i++){buf[i] = srcP[isp*12*12+i];}
      for(unsigned int d0=0;d0<12;d0++)
      for(unsigned int d1=0;d1<12;d1++)
      {
        srcP[(isp*12+d0)*12+d1] = buf[d1*12+d0];
      }
    });
  }
}

template<typename Td, typename Ty>
void fieldG_to_prop4d(Propagator4dT<Td>& res, FieldG<Ty >& src){
  prop4d_to_fieldG(src, res, 0);
}

// Field res Nvol --> 12 x 12
// Field src 12 x 12 --> Nvol
template<typename Td, typename Ty>
void prop_gpu_to_qprop(qlat::Field<Td>& res, qlat::vector_gpu<Ty >& src, Int dir = 1){
  const Long Ndc = 12 * 12;
  Qassert(res.initialized and res.multiplicity == Ndc);
  if(dir == 1){Qassert(Long(src.size()) == res.geo().local_volume() * Ndc);}
  if(dir == 0){src.resize(res.geo().local_volume() * Ndc);}
  const Long Ncopy = src.size();
  Td* rP = (Td*) qlat::get_data(res).data();
  Ty* sP = (Ty*) qlat::get_data(src).data();
  move_index mv_civ;
  if(dir == 1){
    cpy_GPU(rP, sP, Ncopy, 1, 1);
    mv_civ.move_civ_in(rP, rP, 1, Ndc, Ncopy/Ndc, 1, true);
  }
  if(dir == 0){
    cpy_GPU(sP, rP, Ncopy, 1, 1);
    mv_civ.move_civ_out(sP, sP, 1, Ncopy/Ndc, Ndc, 1, true);
  }
}

template<typename Td, typename Ty>
void qprop_to_prop_gpu(qlat::vector_gpu<Ty >& res, qlat::Field<Td>& src){
  prop_gpu_to_qprop(src, res, 0);
}

////
/*
  assumed civ == n*12 with n the source indices, 12 the sink indices 
  c_add : 0 clear res, 1 add, -1 subtract
*/
template <typename Ty, class Fieldy, Int c_add >
void copy_bsize_prop_to_FieldP(std::vector<Fieldy >& res, Ty* src, const LInt nV, LInt b_size, qlat::fft_desc_basic& fd, Int GPU = 1, bool rotate = false, Int dir = 0)
{
  TIMERA("copy_bisze_src_to_FieldM");
  const Int civ = 12 * 12;
  if(civ%12 != 0){abort_r("FieldM type not supported!\n");}
  //unsigned int nV = 0;
  Int cfac = civ/12;
  move_index mv_civ;

  //Ty --> double float ...
  Qassert(GetBasicDataType<Fieldy>::get_type_name() != std::string("unknown_type"));
  using Dy = typename GetBasicDataType<Fieldy>::ElementaryType;
  Qassert(IsBasicTypeReal<Dy>());
  if(c_add != 0){Qassert(rotate == false);}

  const Int  NTt  = fd.Nv[3];
  const LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  const LInt sizeF = NTt*Nxyz;
  const LInt total = 6*sizeF;
  if(total % b_size != 0){abort_r("bsize vec configurations wrong! \n");}

  if(dir == 0){
    //Long dsize = src.size();
    //if(dsize%(2*total) != 0){abort_r("src size wrong!\n");};
    //nV  = dsize/(2*total);
    if(nV%(cfac) != 0){abort_r("res civ wrong!\n");}
    unsigned int ntem = nV/cfac;

    bool do_ini = false;
    if(res.size() != ntem){do_ini = true;}
    if(do_ini == false){
      for(LInt iv=0;iv<res.size();iv++){
        if(!res[iv].initialized){
          do_ini = true;
        }
      }
    }
    // could not add if need initialization
    if(do_ini == true){Qassert(c_add == 0);}

    if(do_ini){
      const Geometry& geo = fd.geo();
      res.resize(0);res.resize(ntem);
      for(LInt iv=0;iv<res.size();iv++){res[iv].init(geo, civ);}
    }
  }
  if(dir == 1){
    //nV = res.size() * cfac;
    //src.resize(nV * 2*total);
    Qassert(nV == res.size() * cfac )
    Qassert(res[0].multiplicity == civ);
  }

  /////rotate FieldM, from Vol->civ to civ->Vol
  if(dir == 1 and rotate == true){
    for(LInt iv=0;iv<res.size();iv++){
      ComplexT<Dy >* s0 = (ComplexT<Dy >*) qlat::get_data(res[iv]).data();
      mv_civ.dojob(s0, s0, 1, civ, sizeF, 1, 1, GPU);
    }
  }

  const Long bfac = total/(b_size);
        Long each  = Nxyz; if(b_size < Nxyz){each = b_size;}
  const Long group = (2*total)/each;

  //Ty* psrc       = src.data();
  Ty* psrc       = src;

  ////buffers for result pointers
  qlat::vector<ComplexT<Dy >* > resP;resP.resize(0);resP.resize(res.size());
  for(unsigned int d0=0;d0<res.size();d0++){
    resP[d0] = (Ty*) qlat::get_data(res[d0]).data();
  }

  const Long Nfac = get_threads(32, each);Qassert(each % Nfac == 0);
  const Long Nsize = each / Nfac;
  const Long dsum = Long(nV) * group * Nsize;
  qGPU_for(isp, dsum, GPU, {
    const Long d0  =  isp / (group * Nsize);
    const Long gi  = (isp % (group * Nsize)) / Nsize;
    const Long offi = isp % Nsize;
    const Long mi = gi*each;

    ////index for res
    const Long d1 =  mi/(NTt*Nxyz);
    const Long ti = (mi/(Nxyz))%NTt;
    const Long vi =  mi%(Nxyz);
    const Long d0a = d0/cfac;
    const Long d0b = d0%cfac;

    ////index for src
    const Int chi = mi/(total);
    const Long xi = mi%(total);
    const Long bi = xi/b_size;
    const Long bj = xi%b_size;
    Ty* s0            = &psrc[(chi*bfac+bi)*nV*b_size  + d0*b_size + bj + offi * Nfac];
    ComplexT<Dy >* s1 = &resP[d0a][((d0b*12 + d1)*NTt+ti)*Nxyz + vi + offi * Nfac];
    if(dir == 0){
      for(Long i=0;i<Nfac;i++){
        if( c_add == 0){s1[i]  = s0[i];}
        if( c_add ==-1){s1[i] -= s0[i];}
        if( c_add == 1){s1[i] += s0[i];}
      }
    }
    if(dir == 1){
      for(Long i=0;i<Nfac;i++){
        if( c_add == 0){s0[i]  = s1[i];}
        if( c_add ==-1){s0[i] -= s1[i];}
        if( c_add == 1){s0[i] += s1[i];}
      }
    }
  })

  qacc_barrier(dummy);

  if(dir == 0 and rotate == true){
    for(LInt iv=0;iv<res.size();iv++){
      ComplexT<Dy >* s0 = (ComplexT<Dy >*) qlat::get_data(res[iv]).data();
      mv_civ.dojob(s0, s0, 1, civ, sizeF, 0, 1, GPU);
    }
  }
}

//  assumed civ == n*12 with n the source indices, 12 the sink indices 
template <typename Ty, class Fieldy >
void copy_bsize_prop_to_FieldM(std::vector<Fieldy >& res, qlat::vector_gpu<Ty >& src, LInt b_size, qlat::fft_desc_basic& fd, Int GPU = 1, bool rotate = false, Int dir = 0)
{
  const LInt total = 12 * fd.Nvol ;
  LInt nV = 0;
  if(dir == 0){
    const Long dsize = src.size();
    if(dsize%(total) != 0){abort_r("src size wrong!\n");};
    nV  = dsize/(total);
  }
  if(dir == 1){
    const Int cfac = 12;
    nV = res.size() * cfac;
    src.resize(nV * total);
  }
  Ty* srcP = src.data();
  copy_bsize_prop_to_FieldP<Ty, Fieldy, 0>(res, srcP, nV, b_size, fd, GPU, rotate, dir);
}

template <typename Ty, class Fieldy >
void copy_FieldM_to_bsize_prop(qlat::vector_gpu<Ty >& res, std::vector<Fieldy >& src, LInt b_size, qlat::fft_desc_basic& fd, Int GPU = 1, bool rotate = false)
{
  copy_bsize_prop_to_FieldM(src, res, b_size,fd, GPU, rotate, 1);
}

//template <typename Ty, class Fieldy >
//void copy_bsize_prop_to_FieldG(std::vector<Fieldy >& res, qlat::vector_gpu<Ty >& src, LInt b_size, const Geometry& geo, Int dir = 0)
//{
//  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);
//  copy_bsize_prop_to_FieldM(res, src, b_size, fd, dir );
//  for(LInt iv=0;iv<res.size();iv++){
//    res[iv].mem_order = QLAT_OUTTER;
//  }
//}
//
//template <typename Ty, class Fieldy >
//void copy_FieldG_to_bsize_prop(qlat::vector_gpu<Ty >& res, std::vector<Fieldy >& src, LInt b_size, const Geometry& geo)
//{
//  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);
//  for(LInt iv=0;iv<res.size();iv++){
//    Qassert(res[iv].mem_order == QLAT_OUTTER);
//  }
//  copy_bsize_prop_to_FieldM(src, res, b_size, fd, 1 );
//}

template<typename Ty>
void ini_propG(std::vector<qlat::vector_gpu<Ty > >& prop, const Long nmass, size_t Nsize, bool clear = true){
  if(Long(prop.size()) != nmass){prop.resize(nmass);}
  for(unsigned long i=0;i<prop.size();i++){
    if(prop[i].size() != Nsize){
      prop[i].resize(Nsize);
    }
    else{
      if(clear){prop[i].set_zero();}
    }
  }
}

////dir == 1, from src to EigenG res
////dir == 0  from EigenG res to src
////EigenG is the prop type needed for fast contractions
////resG, nmass, --> 12 x 12 --> Nvol
template <typename T, typename Ty>
void copy_eigen_prop_to_EigenG(std::vector<qlat::vector_gpu<Ty > >& resG, T* src, 
  LInt b_size, Int nmass, qlat::fft_desc_basic& fd, Int GPU = 1, Int dir = 1)
{
  TIMERA("copy_eigen_prop_to_EigenG");
  if(nmass == 0){resG.resize(0); return ;}
  if(dir == 1){
    ini_propG(resG, nmass, 12*12*size_t(fd.Nvol), false);
  }

  Int Ns    = nmass*12;
  Int  NTt  = fd.Nv[3];
  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  LInt total = 6*NTt*Nxyz;
  if(total % b_size != 0){abort_r("eigen system configurations wrong! \n");}
  LInt bfac = total/(b_size);
  LInt each  = Nxyz; if(b_size < Nxyz){each = b_size;}
  LInt group = (2*total)/each;

  ////buffers for result pointers
  qlat::vector<Ty* > resP;resP.resize(0);resP.resize(resG.size());
  for(unsigned int d0=0;d0<resG.size();d0++){
    resP[d0] = &resG[d0][0];
  }
  const Long Nfac = get_threads(32, each);Qassert(each % Nfac == 0);
  const Long Nsize = each / Nfac;
  const Long dsum = Long(Ns) * group * Nsize;
  qGPU_for(isp, dsum, GPU, {
    const Long d0  =  isp / (group * Nsize);
    const Long gi  = (isp % (group * Nsize)) / Nsize;
    const Long offi = isp % Nsize;

    const Int massi = d0/12;
    const Int d0i   = d0%12;
    const Long mi = gi*each;

    const Long d1 =  mi/(NTt*Nxyz);
    const Long ti = (mi/(Nxyz))%NTt;
    const Long vi =  mi%(Nxyz);

    ////index for src
    Int chi = mi/(total);
    const Long xi = mi%(total);
    const Long bi = xi/b_size;
    const Long bj = xi%b_size;
    T*  s0  = &src[(chi*bfac+bi)*Ns*b_size  + d0*b_size + bj + offi * Nfac];
    Ty* s1  = (Ty*) &resP[massi][((d0i*12 + d1)*NTt+ti)*Nxyz + vi + offi * Nfac];
 
    if(dir == 0){
      for(Long i=0;i<Nfac;i++){
        s0[i] = s1[i];
      }
    }
    if(dir == 1){
      for(Long i=0;i<Nfac;i++){
        s1[i] = s0[i];
      }
    }
  })

  //for(Int d0=0;d0<Ns;d0++)
  //for(LInt gi=0;gi<group;gi++)
  //{
  //  Int massi = d0/12;
  //  Int d0i   = d0%12;
  //  LInt mi = gi*each;

  //  ////index for res
  //  LInt d1 =  mi/(NTt*Nxyz);
  //  LInt ti = (mi/(Nxyz))%NTt;
  //  LInt vi =  mi%(Nxyz);

  //  ////index for src
  //  Int chi = mi/(total);
  //  LInt xi = mi%(total);
  //  Long bi = xi/b_size;
  //  Long bj = xi%b_size;

  //  T* s0   = &src[(chi*bfac+bi)*Ns*b_size  + d0*b_size + bj];
  //  Ty* s1  = (Ty*) &res[massi][((d0i*12 + d1)*NTt+ti)*Nxyz + vi];
  //  if(dir == 1){cpy_data_thread(s1, s0, each , GPU, QFALSE);}
  //  if(dir == 0){cpy_data_thread(s0, s1, each , GPU, QFALSE);}

  //}
  qacc_barrier(dummy);

}

template <typename T, typename Ty>
void copy_EigenG_to_eigen_prop(T* res, std::vector<qlat::vector_gpu<Ty > >& src,
  LInt b_size, Int nmass, qlat::fft_desc_basic& fd, Int GPU = 1)
{
  copy_eigen_prop_to_EigenG(src, res, b_size, nmass, fd, GPU, 0);
}

/////res in format src 12 * sink 12 --> Nt * Nxyz, diagonal sources
template <typename Ty >
void FieldM_src_to_FieldM_prop(qlat::FieldM<Ty , 1>& src, qlat::FieldM<Ty , 12*12>& res, Int GPU = true, bool dummy = true)
{
  qlat::Geometry& geo = src.geo();

  if(!res.initialized){res.init(geo);}
  Qassert(res.geo() == geo);

  //bool do_ini = true;
  //if(res.size() == src.size())if(res[src.size()-1].initialized){do_ini = false;}
  //if(do_ini){res.resize(nV);for(Int iv=0;iv<nV;iv++){res[iv].init(geo);}}
  clear_fields(res);//set all zero... Bug introduced...

  //std::vector<Int > nv, Nv, mv;
  //geo_to_nv(geo, nv, Nv,mv);
  Long Ncopy = geo.local_volume();

  Ty* s0 = NULL; Ty* s1 = NULL;Ty* st = NULL;
  ///for(Int iv=0;iv<nV;iv++)
  s0 = (Ty*) qlat::get_data(src).data();
  st = (Ty*) qlat::get_data(res).data();
  for(unsigned int d0=0;d0<12;d0++)
  {
    //////diagonal elements
    s1 = &st[(d0*12+d0)*Ncopy + 0];
    cpy_data_thread(s1, s0, Ncopy , GPU, QFALSE);
  }
  if (dummy) {
    qacc_barrier(dummy);
  }
}

template <typename Ty >
void FieldM_src_to_FieldM_prop(std::vector<qlat::FieldM<Ty , 1> >& src, std::vector<qlat::FieldM<Ty , 12*12> >& res, Int GPU = true)
{
  if(src.size() == 0){return ;}
  ////qlat::Geometry& geo = src[0].geo();
  Long nV = src.size();
  if(res.size() != src.size()){res.resize(nV);}
  for(Int iv=0;iv<nV;iv++)FieldM_src_to_FieldM_prop(src[iv], res[iv], GPU, false);
  qacc_barrier(dummy);

}

template <typename Ty >
void FieldG_src_to_FieldG_prop(qlat::FieldG<Ty>& src, qlat::FieldG<Ty>& res)
{
  Qassert(src.initialized and src.multiplicity == 1);
  qlat::Geometry& geo = src.geo();
  const Int Ndc = 12 * 12;

  if(!res.initialized){res.init(geo, Ndc, QMGPU, QLAT_OUTTER);}
  Qassert(res.geo() == geo);
  clear_fields(res);//set all zero... Bug introduced..., below only diagnal will be set

  Long Ncopy = geo.local_volume();

  Ty* s0 = NULL; Ty* s1 = NULL;Ty* st = NULL;
  ///for(Int iv=0;iv<nV;iv++)
  s0 = (Ty*) qlat::get_data(src).data();
  st = (Ty*) qlat::get_data(res).data();
  const Int GPU = 1;
  for(unsigned int d0=0;d0<12;d0++)
  {
    //  diagonal elements
    s1 = &st[(d0*12+d0)*Ncopy + 0];
    cpy_data_thread(s1, s0, Ncopy , GPU, QFALSE);
  }
  qacc_barrier(dummy);
}

template <class Ty>
void copy_color_prop(qlat::vector_gpu<Ty >& res, std::vector<colorFT >& src, Int dir = 1)
{
  Qassert(src.size() == 3);
  qlat::vector<Ty* > srcP;srcP.resize(3);
  for(Int ic=0;ic<3;ic++){
    Qassert(src[ic].initialized);
    srcP[ic] = (Ty*) qlat::get_data(src[ic]).data();
  }

  const qlat::Geometry& geo = src[0].geo();
  const Long V = geo.local_volume();
  Qassert(Long(res.size()) == V*9);
  Ty* resP = res.data();

  if(dir == 1){
  qacc_for(isp, V, {
    for(Int c0=0;c0<3;c0++)
    for(Int c1=0;c1<3;c1++)
    {
      resP[isp*9 + c1*3 + c0 ] = srcP[c0][isp*3 + c1];
    }
  });}

  if(dir == 0){
  qacc_for(isp, V, {
    for(Int c0=0;c0<3;c0++)
    for(Int c1=0;c1<3;c1++)
    {
      srcP[c0][isp*3 + c1] = resP[isp*9 + c1*3 + c0];
    }
  });}
}
template <class Ty>
void copy_to_color_prop(std::vector<colorFT >& res, qlat::vector_gpu<Ty >& src)
{
  copy_color_prop(src, res, 0);
}



}

#endif

