// utils_props_type.h
// Gen Wang
// Mar. 2024

#ifndef UTILS_PROPS_TYPE_H
#define UTILS_PROPS_TYPE_H

#pragma once
#include "general_funs.h"
#include "utils_fft_desc.h"

////////dim 12*12 --> Nt --> Nxyz
///////only memory size and geo are used, donot use others
//#define qprop   qlat::FieldM<Complexq, 12*12>
//#define qpropT  qlat::FieldM<Ty, 12*12>

////std::vector<qlat::FermionField4dT<Tc > >, gwu format, need to rotate bases

namespace qlat{

template<typename Td, typename Tc>
void prop4d_to_Fermion(std::vector<qlat::FermionField4dT<Tc > > &buf, Propagator4dT<Td>& prop, int dir=1){

  if(dir==1){buf.resize(0);buf.resize(12);for(int iv=0;iv<12;iv++){
    if(!buf[iv].initialized){buf[iv].init(prop.geo());}
  }}
  if(dir==0){Qassert(buf.size() == 12);if(!prop.initialized){prop.init(buf[0].geo());}}

  #pragma omp parallel for
  for (Long index = 0; index < prop.geo().local_volume(); ++index)
  {
    qlat::WilsonMatrixT<Td>& src =  prop.get_elem_offset(index);
    for(int d0=0;d0<12;d0++)
    {
      ////v0(s*3 + c0, ga.ind[d0]*3 + c1)
      qlat::ComplexT<Tc>* res = (qlat::ComplexT<Tc>*)&(buf[d0].get_elem_offset(index));
      for(int d1=0;d1<12;d1++)
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
void prop4d_to_qprop(qpropT& res, Propagator4dT<Td>& src, int dir = 1){
  TIMERA("prop4d_to_qprop");
  if(dir == 1){Qassert(src.initialized);res.init();res.init(src.geo());}
  if(dir == 0){Qassert(res.initialized);src.init();src.init(res.geo());}

  Long sizeF = src.geo().local_volume();

  move_index mv_civ;
  qlat::ComplexT<Td>* ps; Ty* pt;
  ps = (qlat::ComplexT<Td>* ) qlat::get_data(src).data();
  pt = (Ty*) qlat::get_data(res).data();

  ////V x 12 a x 12 b to 12b x 12a x V
  if(dir == 1){
    qthread_for(isp, Long(sizeF),{
      QLAT_ALIGN(QLAT_ALIGNED_BYTES) qlat::ComplexT<Td> buf[12*12];for(unsigned int i=0;i<12*12;i++){buf[i] = ps[isp*12*12 + i];}
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
      QLAT_ALIGN(QLAT_ALIGNED_BYTES) qlat::ComplexT<Td> buf[12*12];for(unsigned int i=0;i<12*12;i++){buf[i] = pt[isp*12*12 + i];}
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

////assumed civ == n*12 with n the source indices, 12 the sink indices 
template <typename Ty, int civ >
void copy_eigen_src_to_FieldM(qlat::vector_gpu<Ty >& src, std::vector<qlat::FieldM<Ty , civ> >& res, LInt b_size, qlat::fft_desc_basic& fd, int dir = 0, int GPU = 1, bool rotate = false)
{
  TIMERA("copy_eigen_src_to_FieldM");
  if(civ%12 != 0){abort_r("FieldM type not supported!\n");}
  unsigned int nV = 0;int cfac = civ/12;
  move_index mv_civ;

  int  NTt  = fd.Nv[3];
  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  LInt sizeF = NTt*Nxyz;
  LInt total = 6*sizeF;
  if(total % b_size != 0){abort_r("eigen system configurations wrong! \n");}

  if(dir == 0){
    Long dsize = src.size();
    if(dsize%(2*total) != 0){abort_r("src size wrong!\n");};
    nV  = dsize/(2*total);
    if(nV%(cfac) != 0){abort_r("res civ wrong!\n");}
    unsigned int ntem = nV/cfac;

    bool do_ini = true;if(res.size() == ntem)if(res[ntem-1].initialized){do_ini = false;}
    if(do_ini){
      //////print0("initial Fprop. \n");
      Geometry geo;fd.get_geo(geo);
      res.resize(0);res.resize(ntem);
      for(LInt iv=0;iv<res.size();iv++){res[iv].init(geo);}}
  }
  if(dir == 1){
    nV = res.size() * cfac;
    src.resize(nV * 2*total);
  }

  /////rotate FieldM, from Vol->civ to civ->Vol
  if(dir == 1 and rotate == true){
    for(LInt iv=0;iv<res.size();iv++){
      Ty* s0 = (Ty*) qlat::get_data(res[iv]).data();
      mv_civ.dojob(s0, s0, 1, civ, sizeF, 1, 1, GPU);
    }
  }

  const Long bfac = total/(b_size);
        Long each  = Nxyz; if(b_size < Nxyz){each = b_size;}
  const Long group = (2*total)/each;

  Ty* psrc       = src.data();

  ////buffers for result pointers
  qlat::vector_acc<Ty* > resP;resP.resize(0);resP.resize(res.size());
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
    const int chi = mi/(total);
    const Long xi = mi%(total);
    const Long bi = xi/b_size;
    const Long bj = xi%b_size;
    Ty* s0 = &psrc[(chi*bfac+bi)*nV*b_size  + d0*b_size + bj + offi * Nfac];
    Ty* s1 = &resP[d0a][((d0b*12 + d1)*NTt+ti)*Nxyz + vi + offi * Nfac];
    if(dir == 0){
      for(Long i=0;i<Nfac;i++){
        s1[i] = s0[i];
      }
    }
    if(dir == 1){
      for(Long i=0;i<Nfac;i++){
        s0[i] = s1[i];
      }
    }
  })

  //Ty* s0 = NULL; Ty* s1 = NULL;Ty* st = NULL;
  //for(Long d0=0;d0<nV;d0++)
  //for(Long gi=0;gi<group;gi++)
  //{
  //  LInt mi = gi*each;

  //  ////index for res
  //  LInt d1 =  mi/(NTt*Nxyz);
  //  LInt ti = (mi/(Nxyz))%NTt;
  //  LInt vi =  mi%(Nxyz);
  //  int d0a = d0/cfac;
  //  int d0b = d0%cfac;

  //  ////index for src
  //  int chi = mi/(total);
  //  LInt xi = mi%(total);
  //  Long bi = xi/b_size;
  //  Long bj = xi%b_size;

  //  s0 = &psrc[(chi*bfac+bi)*nV*b_size  + d0*b_size + bj];
  //  st = (Ty*) qlat::get_data(res[d0a]).data();
  //  s1 = &st[((d0b*12 + d1)*NTt+ti)*Nxyz + vi];

  //  if(dir == 0){cpy_data_thread(s1, s0, each , GPU, QFALSE);}
  //  if(dir == 1){cpy_data_thread(s0, s1, each , GPU, QFALSE);}
  //}

  //print0("nV %d, each %d, group %d, bfac %d, b_size %d \n", int(nV), int(each), int(group), int(bfac), int(b_size));
  //print0("===RES norm ");src.print_norm2();
  //for(int i=0;i<nV/cfac;i++){
  //  //Ty* r = (Ty*) qlat::get_data(res[i]).data();
  //  Ty* r = psrc;
  //  print0("==value %+.8e %+.8e \n", r[i * 17].real(), r[i * 17].imag());
  //} 

  qacc_barrier(dummy);

  if(dir == 0 and rotate == true){
    for(LInt iv=0;iv<res.size();iv++){
      Ty* s0 = (Ty*) qlat::get_data(res[iv]).data();
      mv_civ.dojob(s0, s0, 1, civ, sizeF, 0, 1, GPU);
    }
  }
}

template <typename Ty, int civ >
void copy_FieldM_to_eigen_src(std::vector<qlat::FieldM<Ty , civ> >& src, qlat::vector_gpu<Ty >& res, LInt b_size, qlat::fft_desc_basic& fd, int GPU = 1, bool rotate = false)
{
  copy_eigen_src_to_FieldM(res, src, b_size, 1,fd, GPU, rotate);
}

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
  LInt b_size, int nmass, qlat::fft_desc_basic& fd, int GPU = 1, int dir = 1)
{
  TIMERA("copy_eigen_prop_to_EigenG");
  if(nmass == 0){resG.resize(0); return ;}
  if(dir == 1){
    ini_propG(resG, nmass, 12*12*size_t(fd.Nvol), false);
  }

  int Ns    = nmass*12;
  int  NTt  = fd.Nv[3];
  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  LInt total = 6*NTt*Nxyz;
  if(total % b_size != 0){abort_r("eigen system configurations wrong! \n");}
  LInt bfac = total/(b_size);
  LInt each  = Nxyz; if(b_size < Nxyz){each = b_size;}
  LInt group = (2*total)/each;

  ////buffers for result pointers
  qlat::vector_acc<Ty* > resP;resP.resize(0);resP.resize(resG.size());
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

    const int massi = d0/12;
    const int d0i   = d0%12;
    const Long mi = gi*each;

    const Long d1 =  mi/(NTt*Nxyz);
    const Long ti = (mi/(Nxyz))%NTt;
    const Long vi =  mi%(Nxyz);

    ////index for src
    int chi = mi/(total);
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

  //for(int d0=0;d0<Ns;d0++)
  //for(LInt gi=0;gi<group;gi++)
  //{
  //  int massi = d0/12;
  //  int d0i   = d0%12;
  //  LInt mi = gi*each;

  //  ////index for res
  //  LInt d1 =  mi/(NTt*Nxyz);
  //  LInt ti = (mi/(Nxyz))%NTt;
  //  LInt vi =  mi%(Nxyz);

  //  ////index for src
  //  int chi = mi/(total);
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

/////res in format src 12 * sink 12 --> Nt * Nxyz, diagonal sources
template <typename Ty >
void FieldM_src_to_FieldM_prop(qlat::FieldM<Ty , 1>& src, qlat::FieldM<Ty , 12*12>& res, int GPU = true, bool dummy = true)
{
  qlat::Geometry& geo = src.geo();

  if(!res.initialized){res.init(geo);}

  //bool do_ini = true;
  //if(res.size() == src.size())if(res[src.size()-1].initialized){do_ini = false;}
  //if(do_ini){res.resize(nV);for(int iv=0;iv<nV;iv++){res[iv].init(geo);}}

  //std::vector<int > nv, Nv, mv;
  //geo_to_nv(geo, nv, Nv,mv);
  Long Ncopy = geo.local_volume();

  Ty* s0 = NULL; Ty* s1 = NULL;Ty* st = NULL;
  ///for(int iv=0;iv<nV;iv++)
  s0 = (Ty*) qlat::get_data(src).data();
  st = (Ty*) qlat::get_data(res).data();
  for(unsigned int d0=0;d0<12;d0++)
  {
    //////diagonal elements
    s1 = &st[(d0*12+d0)*Ncopy + 0];
    cpy_data_thread(s1, s0, Ncopy , GPU, QFALSE);
  }
  if(dummy)qacc_barrier(dummy);
}

template <typename Ty >
void FieldM_src_to_FieldM_prop(std::vector<qlat::FieldM<Ty , 1> >& src, std::vector<qlat::FieldM<Ty , 12*12> >& res, int GPU = true)
{
  if(src.size() == 0){return ;}
  ////qlat::Geometry& geo = src[0].geo();
  Long nV = src.size();
  if(res.size() != src.size()){res.resize(nV);}
  for(int iv=0;iv<nV;iv++)FieldM_src_to_FieldM_prop(src[iv], res[iv], GPU, false);
  qacc_barrier(dummy);

}



}

#endif

