// utils_corr_prop.h
// Gen Wang
// Oct. 2022

#ifndef UTILS_CORR_PROP_H
#define UTILS_CORR_PROP_H

#pragma once

#include "utils_float_type.h"
#include "utils_gammas.h"
#include "utils_fft_desc.h"
#include "utils_reduce_vec.h"
#include "utils_grid_src.h"
#include "utils_shift_vecs.h"
#include "utils_field_operations.h"
#include "utils_field_gpu.h"
#include "utils_sector_funs.h"

namespace qlat{

template<typename Ty>
void init_qpropT(std::vector<qpropT >& res, const unsigned int size, const Geometry& geo)
{
  if(res.size() != size){
    res.resize(size);
  }
  for(unsigned int i=0;i<res.size();i++){
    if(!res[i].initialized){res[i].init(geo);}
  }

}

template<typename Ty>
void clear_qpropT(std::vector<qpropT >& res)
{
  for(unsigned int i=0;i<res.size();i++){
    qlat::set_zero(res[i]);
  }
}

template<typename Td>
void prop4d_conj(Propagator4dT<Td >& prop, Int rotate = 1){
  TIMERA("prop4d_conj");
  ////Rowmajor (a,b), b is continues in memory
  qacc_for(isp, Long(prop.geo().local_volume()),{
    qlat::WilsonMatrixT<Td>& v0 =  prop.get_elem_offset(isp);
    qlat::WilsonMatrixT<Td>  v1 = v0;

    for(Int c0 = 0;c0< 3; c0++)
    for(Int d0 = 0;d0< 4; d0++)
    for(Int c1 = 0;c1< 3; c1++)
    for(Int d1 = 0;d1< 4; d1++)
    {
      if(rotate == 1){
        v0(d0*3 + c0, d1*3 + c0) = qlat::qconj( v1(d1*3 + c1, d0*3 + c0) ) ;
      }
      if(rotate == 0){
        v0(d0*3 + c0, d1*3 + c0) = qlat::qconj( v1(d0*3 + c0, d1*3 + c1) ) ;
      }
    }
  });
}

template<typename Ty, Int dir, bool conj>
void FieldM_gammaT(qlat::FieldM<Ty , 12>& vec, ga_M& ga){
  TIMERA("FieldM_gammaT");
  const Geometry& geo = vec.geo();
  //const long Nvol = geo.local_volume();
  Ty* src = (Ty*) qlat::get_data(vec).data();
  qacc_for(isp, Long(geo.local_volume()),{
    Ty buf[12];
    Ty* v = &src[isp* 12 + 0];
    for(Int i=0;i<12;i++){buf[i] = v[i];}

    for(Int d0 = 0;d0 < 4; ++d0)
    {
      /////Source multiply
      if(dir == 0)
      for(Int c0=0;c0<3;c0++){
        if(!conj){v[ga.ind[d0]*3 + c0] =             ga.g[d0] * buf[d0*3 + c0] ;}
        if( conj){v[ga.ind[d0]*3 + c0] = qlat::qconj(ga.g[d0] * buf[d0*3 + c0]);}
      }
      /////Sink   multiply
      if(dir == 1)
      for(Int c0=0;c0<3;c0++){
        if(!conj){v[d0*3 + c0] =             ga.g[d0] * buf[ga.ind[d0]*3 + c0] ;}
        if( conj){v[d0*3 + c0] = qlat::qconj(ga.g[d0] * buf[ga.ind[d0]*3 + c0]);}
      }
    }
  });
}

template<typename Ty>
void FieldM_gamma_src(qlat::FieldM<Ty , 12>& vec, ga_M& ga){
  FieldM_gammaT<Ty, 0, false>(vec, ga);
}

template<typename Ty>
void FieldM_gamma_sink(qlat::FieldM<Ty , 12>& vec, ga_M& ga){
  FieldM_gammaT<Ty, 1, false>(vec, ga);
}

template<typename Td, Int dir, bool conj>
void prop4d_src_gammaT(Propagator4dT<Td >& prop, ga_M& ga){
  TIMERA("prop4d_src_gamma");
  // Rowmajor (a,b), b is continues in memory
  qacc_for(isp, Long(prop.geo().local_volume()),{
    qlat::WilsonMatrixT<Td>& v0 =  prop.get_elem_offset(isp);
    qlat::WilsonMatrixT<Td>  v1 = v0;

    for(Int s = 0; s < 4; ++s)
    for(Int c0 = 0;c0< 3 ; c0++)
    for(Int d0 = 0; d0 < 4; ++d0)
    {
      // Source multiply
      if(dir==0)for(Int c1=0;c1<3;c1++){
        if(!conj){v0(s*3 + c0, ga.ind[d0]*3 + c1) = ga.g[d0] * v1(s*3 + c0, d0*3 + c1);}
        if( conj){v0(s*3 + c0, ga.ind[d0]*3 + c1) = qlat::qconj(ga.g[d0] * v1(s*3 + c0, d0*3 + c1));}
      }
      // Sink multiply
      if(dir==1)for(Int c1=0;c1<3;c1++){
        if(!conj){v0(d0*3 + c0, s*3 + c1) = ga.g[d0] * v1(ga.ind[d0]*3 + c0, s*3 + c1);}
        if( conj){v0(d0*3 + c0, s*3 + c1) = qlat::qconj(ga.g[d0] * v1(ga.ind[d0]*3 + c0, s*3 + c1));}
      }
    }
  });
}

template<typename Td>
void prop4d_src_gamma(Propagator4dT<Td >& prop, ga_M& ga, bool conj = false){
  if(!conj){prop4d_src_gammaT<Td, 0, false>(prop, ga);}
  if( conj){prop4d_src_gammaT<Td, 0, true >(prop, ga);}
}
template<typename Td>
void prop4d_sink_gamma(Propagator4dT<Td >& prop, ga_M& ga, bool conj = false){
  if(!conj){prop4d_src_gammaT<Td, 1, false>(prop, ga);}
  if( conj){prop4d_src_gammaT<Td, 1, true >(prop, ga);}
}

template<typename Ty, Int dir, bool conj>
void qprop_src_gamma_T(Ty* res, ga_M& ga, const Long Nsize){
  TIMERA("qprop_src_gamma");
  ////Rowmajor (a,b), b is continues in memory
  ///Ty* res = (Ty*) qlat::get_data(prop1).data()
  //const Long Nsize = prop.geo().local_volume();
  qacc_for(isp, Nsize,{
    Ty src[12*12];
    Ty buf[12*12];
    for(Int i=0;i<12*12;i++){src[i] = res[i*Nsize + isp];}

    for(Int d1 = 0;d1 < 4; d1++)
    for(Int c1 = 0;c1 < 3; c1++)
    for(Int d0 = 0;d0 < 4; d0++)
    for(Int c0 = 0;c0 < 3; c0++)
    {
      /////Source multiply
      if(dir==0)
      {
        buf[(ga.ind[d1]*3 + c1) * 12 + d0*3 + c0] = ga.g[d1] * src[(d1*3 + c1 )*12 + d0*3 + c0];
      }
      /////Sink multiply
      if(dir==1)
      {
        buf[(d1*3 + c1)*12 + d0*3 + c0] = ga.g[d0] * src[(d1*3 + c1)*12 + ga.ind[d0]*3 + c0];
      }
    }

    for(Int i=0;i<12*12;i++){
      if(!conj){res[i*Nsize + isp] = buf[i];}
      if( conj){res[i*Nsize + isp] = qlat::qconj(buf[i]);}
    }
  });
}

template<typename Ty>
void qprop_src_gamma(qpropT& prop, ga_M& ga, bool conj = false){
  Ty* res = (Ty*) qlat::get_data(prop).data();
  const Long Nsize = prop.geo().local_volume();
  if(!conj)qprop_src_gamma_T<Ty, 0, false>(res, ga, Nsize);
  if( conj)qprop_src_gamma_T<Ty, 0, true >(res, ga, Nsize);
}


template<typename Ty>
void qprop_sink_gamma(qpropT& prop, ga_M& ga, bool conj = false){
  Ty* res = (Ty*) qlat::get_data(prop).data();
  const Long Nsize = prop.geo().local_volume();
  if(!conj)qprop_src_gamma_T<Ty, 1, false>(res, ga, Nsize);
  if( conj)qprop_src_gamma_T<Ty, 1, true >(res, ga, Nsize);
}

template<typename Ty>
void Gprop_src_gamma(EigenTy& prop, ga_M& ga, bool conj = false){
  for(unsigned long iv=0;iv<prop.size();iv++)
  {
    Ty* res = (Ty*) qlat::get_data(prop[iv]).data();
    if(!conj)qprop_src_gamma_T<Ty, 0, false>(res, ga, prop[iv].size()/(12*12));
    if( conj)qprop_src_gamma_T<Ty, 0, true >(res, ga, prop[iv].size()/(12*12));
  }
}

template<typename Ty>
void Gprop_sink_gamma(EigenTy& prop, ga_M& ga, bool conj = false){
  for(unsigned long iv=0;iv<prop.size();iv++)
  {
    Ty* res = (Ty*) qlat::get_data(prop[iv]).data();
    if(!conj)qprop_src_gamma_T<Ty, 1, false>(res, ga, prop[iv].size()/(12*12));
    if( conj)qprop_src_gamma_T<Ty, 1, true >(res, ga, prop[iv].size()/(12*12));
  }
}

// can be expanded fields and selected fields
template<class Fieldy>
void fieldG_src_gamma(Fieldy& prop, ga_M& ga, bool conj = false, const bool srcG = true){
  TIMER("fieldG_src_gamma");
  Qassert(prop.initialized);
  Qassert(GetBasicDataType<Fieldy>::get_type_name() != std::string("unknown_type"));
  using D = typename GetBasicDataType<Fieldy>::ElementaryType;
  Qassert(IsBasicTypeReal<D>());

  Qassert(prop.multiplicity % (12 * 12) == 0);
  Qassert(prop.mem_order == QLAT_OUTTER);
  const Int Nvec   = prop.multiplicity / (12 * 12);
  const Long Nd    = prop.field.size() / Nvec;
  const Long Nsize = prop.field.size() / (prop.multiplicity);
  ComplexT<D >* p = (ComplexT<D >*) get_data(prop).data();
  for(Int iv=0;iv<Nvec;iv++){
    ComplexT<D >* res = &p[iv * Nd];
    if(srcG){
      if(!conj)qprop_src_gamma_T<ComplexT<D >, 0, false>(res, ga, Nsize);
      if( conj)qprop_src_gamma_T<ComplexT<D >, 0, true >(res, ga, Nsize);
    }else{
      if(!conj)qprop_src_gamma_T<ComplexT<D >, 1, false>(res, ga, Nsize);
      if( conj)qprop_src_gamma_T<ComplexT<D >, 1, true >(res, ga, Nsize);
    }
  }
}

template<class Fieldy>
void fieldG_sink_gamma(Fieldy& prop, ga_M& ga, bool conj = false){
  fieldG_src_gamma(prop, ga, conj, false);
}

template<class Ty>
void fieldG_src_gammaG(std::vector<FieldG<Ty > > & prop, ga_M& ga, bool conj = false){
  for(unsigned int i=0;i<prop.size();i++){
    fieldG_src_gamma(prop[i], ga, conj, true);
  }
}

template<class Ty>
void fieldG_sink_gammaG(std::vector<FieldG<Ty > >& prop, ga_M& ga, bool conj = false){
  for(unsigned int i=0;i<prop.size();i++){
    fieldG_src_gamma(prop[i], ga, conj, false);
  }
}

template<typename Ty>
void qprop_move_dc_in(Ty* src, const qlat::Geometry& geo, const Int dir = 1)
{
  move_index mv_civ;
  const Long sizeF = geo.local_volume();

  if(dir == 1){mv_civ.move_civ_in( src, src, 1, 12*12, sizeF, 1, true);}
  if(dir == 0){mv_civ.move_civ_out(src, src, 1, sizeF, 12*12, 1, true);}
}

template<typename Ty>
void qprop_move_dc_in(qpropT& src, const Int dir = 1)
{
  Qassert(src.initialized);
  qprop_move_dc_in((Ty*) qlat::get_data(src).data(), src.geo(), dir);
}

template<typename Ty>
void qprop_move_dc_out(Ty* src, const qlat::Geometry& geo)
{
  qprop_move_dc_in(src, geo, 0);
}

template<typename Ty>
void qprop_move_dc_out(qpropT& src)
{
  qprop_move_dc_in(src, 0);
}

template<typename Ty>
void qprop_sub_add(std::vector<qpropT >& res, std::vector< qpropT >& s0, const Ty f0, const Ty f1)
{
  if(s0.size() == 0){res.resize(0); return ;}
  const Int Nvec = s0.size();
  const qlat::Geometry& geo = s0[0].geo();
  const Long Nvol = geo.local_volume();

  init_qpropT(res, Nvec, geo);
  for(Int vi=0;vi<Nvec;vi++)
  {
    Ty* p0 = (Ty* ) qlat::get_data(s0[vi]).data();
    Ty* r0 = (Ty* ) qlat::get_data(res[vi]).data(); 
    for(Int dc=0;dc<12*12;dc++){
      qacc_for(isp, geo.local_volume(),{
        r0[dc*Nvol + isp] = r0[dc*Nvol + isp]*f0 + p0[dc*Nvol + isp] * f1;
      });
    }
  }
}

template<typename Ty>
void qprop_sub_add(std::vector<qpropT >& res, std::vector< qpropT >& s0, std::vector< qpropT >& s1, const Ty f0, const Ty f1)
{
  if(s0.size() == 0){res.resize(0); return ;}
  Qassert(s0.size() == s1.size());
  const Int Nvec = s0.size();
  const qlat::Geometry& geo = s0[0].geo();
  const Long Nvol = geo.local_volume();

  init_qpropT(res, Nvec, geo);
  for(Int vi=0;vi<Nvec;vi++)
  {
    Ty* p0 = (Ty* ) qlat::get_data(s0[vi]).data();
    Ty* p1 = (Ty* ) qlat::get_data(s1[vi]).data(); 
    Ty* r0 = (Ty* ) qlat::get_data(res[vi]).data(); 
    for(Int dc=0;dc<12*12;dc++){
      qacc_for(isp, geo.local_volume(),{
        r0[dc*Nvol + isp] = (p0[dc*Nvol + isp] + p1[dc*Nvol + isp] * f0) * f1;
      });
    }
  }
}

/*
  covariant shifts to 4 directions
*/
template<typename Ty>
void shift_vecs_cov_qpropT(std::vector< std::vector<qpropT > >& res, std::vector< qpropT >& s0, shift_vec& svec,
  std::vector<std::vector<qpropT >>& buf)
{
  if(res.size() != 5){res.resize(5);}
  if(buf.size() != 2){buf.resize(2);}
  qprop_sub_add(res[0], s0, Ty(0.0,0.0), Ty(1.0,0.0) );////equal

  init_qpropT(buf[0], s0.size(), s0[0].geo());
  init_qpropT(buf[1], s0.size(), s0[0].geo());

  for(Int nu = 0; nu < 4 ; nu++)
  {
    shift_vecs_dir_qpropT(s0, buf[0], nu, +1, svec);
    shift_vecs_dir_qpropT(s0, buf[1], nu, -1, svec);
    qprop_sub_add(res[1 + nu], buf[0], buf[1], Ty(-1.0,0.0), Ty(1.0/2.0, 0.0) );////equal
  }
}

template<typename Td>
void prop4d_cps_to_ps(Propagator4dT<Td >& prop, Int dir=0){
  /////sn is -1 for default
  qlat::ComplexT<Td > sn =-1;if(dir == 1){sn= 1;}
  const qlat::ComplexT<Td>sqrt2= qlat::ComplexT<Td>(std::sqrt(2.0), 0.0);

  ////Rowmajor (a,b), b is continues in memory
  qacc_for(isp, prop.geo().local_volume(),{
    qlat::WilsonMatrixT<Td >  v0 = prop.get_elem_offset(isp);
    qlat::WilsonMatrixT<Td >  v1 = prop.get_elem_offset(isp);

    Int dr,d0,d1;
    /////Src rotation
    for (Int s0 = 0; s0 < 4; ++s0)
    for(Int c0 = 0;c0< 3 ; c0++)
    {
      dr=0;d0=1;d1=3;
      for(Int c1=0;c1<3;c1++)v0(s0*3+c0, dr*3+c1) = (-v1(s0*3+c0, d0*3+ c1)*sn + v1(s0*3+c0, d1*3+c1)   )/sqrt2;
      dr=1;d0=0;d1=2;
      for(Int c1=0;c1<3;c1++)v0(s0*3+c0, dr*3+c1) = (+v1(s0*3+c0, d0*3+ c1)*sn - v1(s0*3+c0, d1*3+c1)   )/sqrt2;
      dr=2;d0=1;d1=3;
      for(Int c1=0;c1<3;c1++)v0(s0*3+c0, dr*3+c1) = (-v1(s0*3+c0, d0*3+ c1)    - v1(s0*3+c0, d1*3+c1)*sn)/sqrt2;
      dr=3;d0=0;d1=2;
      for(Int c1=0;c1<3;c1++)v0(s0*3+c0, dr*3+c1) = (+v1(s0*3+c0, d0*3+ c1)    + v1(s0*3+c0, d1*3+c1)*sn)/sqrt2;
    }

    /////Copy previous results
    v1 = v0;
    /////Sink rotation
    for(Int c0 = 0;c0< 3 ; c0++)
    for (Int s0 = 0; s0 < 4; ++s0)
    {
      dr=0;d0=1;d1=3;
      for(Int c1=0;c1<3;c1++)v0(dr*3+c0, s0*3+c1) = (-v1(d0*3+c0, s0*3+c1)*sn + v1(d1*3+c0, s0*3+c1)   )/sqrt2;
      dr=1;d0=0;d1=2;
      for(Int c1=0;c1<3;c1++)v0(dr*3+c0, s0*3+c1) = ( v1(d0*3+c0, s0*3+c1)*sn - v1(d1*3+c0, s0*3+c1)   )/sqrt2;
      dr=2;d0=1;d1=3;
      for(Int c1=0;c1<3;c1++)v0(dr*3+c0, s0*3+c1) = (-v1(d0*3+c0, s0*3+c1)    - v1(d1*3+c0, s0*3+c1)*sn)/sqrt2;
      dr=3;d0=0;d1=2;
      for(Int c1=0;c1<3;c1++)v0(dr*3+c0, s0*3+c1) = ( v1(d0*3+c0, s0*3+c1)    + v1(d1*3+c0, s0*3+c1)*sn)/sqrt2;
    }
    prop.get_elem_offset(isp) = v0;

  });
}

template<typename Td>
void prop4d_ps_to_cps(Propagator4dT<Td >& prop){
  prop4d_cps_to_ps(prop, 1);
}

template<typename Td>
void get_corr_pion(std::vector<qlat::FermionField4dT<Td > > &prop,const Coordinate &x_ini, std::vector<double > &write ){

  const qlat::Geometry& geo = prop[0].geo();

  unsigned long Nvol = geo.local_volume();
  ///int Nt = geo.node_site[3];
  ///Long Nsum = Nvol/Nt;
  Int tini = x_ini[3];

  qlat::vector<qlat::ComplexT<Td> > res;res.resize(Nvol);

  qacc_for(isp, Long(Nvol),{
    qlat::ComplexT<Td> buf(0.0,0.0);

    for(Int dc2=0;dc2<12;dc2++){
      qlat::ComplexT<Td>* a = (qlat::ComplexT<Td>* ) &(prop[dc2].get_elem_offset(isp));
      for(Int dc1=0;dc1<12;dc1++)
      {
        buf+=a[dc1]*qlat::qconj(a[dc1]);
      }
    }
    res[isp] = buf;
    ////src[isp] = buf[isp];
  });

  const Coordinate vg = geo.total_site();
  Int nt = vg[3];
  write.resize(0);
  write.resize(2*nt);qlat::set_zero(write);

  for(unsigned long isp=0;isp<Nvol;isp++){
    Coordinate xl0 = geo.coordinate_from_index(isp);
    Coordinate xg0 = geo.coordinate_g_from_l(xl0);
    Int t = xg0[3];

    Int toff = ((t-tini+nt)%nt);
    write[ toff*2 + 0 ] += res[isp].real();
    write[ toff*2 + 1 ] += res[isp].imag();
  }
  ////May need to be changed for EigenV
  //sum_all_size((double*) &write[0],2*nt);
  sum_all_size((double*) write.data(), 2*nt);
}

template<typename Ty>
void get_src_phase(Ty& phase, const qlat::vector<Int >& nv,
  const Coordinate& pos = Coordinate(), const Coordinate& mom = Coordinate()){
    double p0[3]={2 * QLAT_PI_LOCAL /nv[0],
                  2 * QLAT_PI_LOCAL /nv[1],
                  2 * QLAT_PI_LOCAL /nv[2]};
    double theta= mom[0]*p0[0]*pos[0] + mom[1]*p0[1]*pos[1] + mom[2]*p0[2]*pos[2];
    phase = Ty(std::cos(theta), -1.0*std::sin(theta));
}

template<typename Ty>
void vec_corrE(Ty* srcE, qlat::vector<Ty >& res,qlat::fft_desc_basic &fd,const Int nvec,const Int clear=0,const Coordinate& mom = Coordinate(), const Ty& src_phase = 1.0, const Int t0 = 0){
  TIMER("Reduce vec_corrE");
  Int NTt  = fd.Nv[3];
  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  Ty* src = srcE;
  //int nmass = resE.size()/(NTt*Nxyz);

  ////position p = fd.desc.get_position(0,fd.rank);
  Int t_rank = fd.Pos0[fd.rank][3];
  qlat::vector_gpu<Ty > bufE;
  if(mom != Coordinate())
  {
    bufE.resize(nvec * Nxyz*NTt);
    cpy_data_thread(bufE.data(), srcE, bufE.size(), 1, QTRUE);
    src = bufE.data();

    qlat::vector<double > p0;p0.resize(3);
    for(Int i=0;i<3;i++){p0[i] = 2 * QLAT_PI_LOCAL /fd.nv[i];}

    qlat::vector<Ty > phaseEG;phaseEG.resize(Nxyz);
    //Ty* phaseE = (Ty*) qlat::get_data(phaseEG).data();

    ////qlat::vector_gpu<int > pos_tem;pos_tem.copy_from(fd.Pos0[fd.rank]);int* posP = pos_tem.data();
    /////===may not be consistent for fd definiations under qacc
    /////===slow
    qthread_for(xi, Long(Nxyz),{
      Coordinate pi = fd.coordinate_g_from_index(xi );

      double theta=mom[0]*p0[0]*pi[0]+mom[1]*p0[1]*pi[1]+mom[2]*p0[2]*pi[2];
      phaseEG[xi] = Ty(cos(theta),sin(theta));
    });

    size_t Ns = Nxyz;
    Ns = nvec*NTt*Nxyz;
    qacc_for(i, Long(Ns),{
      LInt mi = i/(NTt*Nxyz);
      LInt ti = (i%(NTt*Nxyz))/Nxyz;
      LInt xi = i%(Nxyz);
      src[(mi*NTt + ti)*Nxyz + xi] = src[(mi*NTt + ti)*Nxyz + xi]*phaseEG[xi];
    });
    /////#endif
  }

  ////TODO
  Int nt = fd.nt;
  //if(clear == 1){res.resize(0);res.resize(nvecs*nt);qlat::set_zero(res);}
  if(clear == 1){if(res.size() != nvec*nt){res.resize(nvec*nt);} qlat::set_zero(res);}
  if(clear == 0){if(res.size() != nvec*nt){qmessage("res size wrong for corr.\n");Qassert(false);}}

  qlat::vector<Ty > tmp;tmp.resize(nvec*NTt);qlat::set_zero(tmp);//tmp.set_zero();
  reduce_vecs(src, tmp.data(), Nxyz, nvec*NTt);

  qlat::vector_gpu<Ty > RES;RES.resize(nvec*nt );RES.set_zero();
  Ty* s1 = RES.data();Ty* s0 = tmp.data();
  Long Ntotal = nvec*NTt;
  qacc_for(mti, Ntotal, {
    Long mi  = mti/NTt;
    Long  ti = mti%NTt;
    s1[mi*nt + (t_rank + ti + nt - t0)%nt ] = s0[mi*NTt + ti] * src_phase;
  });

  //////sum_all_size((Ftype*) (RES.data()), 2*RES.size(), 1);
  sum_all_size(RES.data(), RES.size(), 1);
  cpy_data_thread((Ty*) qlat::get_data(res).data(), RES.data(), RES.size(), 1, QTRUE, 1.0);
}

template<typename Ty>
void vec_corrE(qlat::vector_gpu<Ty >& resE, qlat::vector<Ty >& res,qlat::fft_desc_basic &fd,const Int clear=0,
  const Coordinate& mom = Coordinate(), const Ty& src_phase = 1.0, Int t0=0){
  Int NTt  = fd.Nv[3];
  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  Int nvec = resE.size()/(NTt*Nxyz);
  vec_corrE(resE.data(), res, fd, nvec, clear, mom, src_phase, t0);
}

template<typename Ty>
void vec_corrE(qlat::vector<Ty >& resE, qlat::vector<Ty >& res,qlat::fft_desc_basic &fd,const Int clear=0,
  const Coordinate& mom = Coordinate(), const Ty& src_phase = 1.0, Int t0 = 0){
  Int NTt  = fd.Nv[3];
  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  Int nvec = resE.size()/(NTt*Nxyz);
  ///qlat::vector_gpu<Ty > r0;r0.copy_from(resE);
  Ty* r0 = (Ty*) qlat::get_data(resE).data();
  vec_corrE(r0, res, fd, nvec, clear, mom, src_phase, t0);
}

template<typename Ty>
void shift_result_t(qlat::vector<Ty >& Esrc, Int nt, Int tini){
  if(tini == 0){return ;}
  Long Ntotal = Esrc.size();
  if(Ntotal %(nt) != 0){abort_r("Correlation function size wrong!\n");}
  qlat::vector<Ty > tmp;tmp.resize(Ntotal);
  qacc_for(i, Ntotal, {
    const Int iv = i/nt;
    const Int t  = i%nt;
    tmp[iv*nt + (t - tini + nt)%nt] = Esrc[iv*nt + (t)%nt];
  });
  cpy_data_thread(Esrc.data(), tmp.data(), tmp.size(), 1);
}

//template<typename Ta>
//void ini_propE(EigenMTa &prop,Int nmass, qlat::fft_desc_basic &fd, bool clear = true){
//  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
//  Int NTt  = fd.Nv[3];
//  Int do_resize = 0;
//  if(prop.size() != (LInt) (nmass*12*12*NTt)){do_resize = 1;}
//  for(unsigned int i=0;i<prop.size();i++){if((LInt) prop[i].size() != Nxyz){do_resize=1;}}
//
//  if(do_resize == 1)
//  {
//    for(unsigned int i=0;i<prop.size();i++){prop[i].resize(0);}prop.resize(0);
//    prop.resize(nmass*12*12*NTt);
//    for(unsigned int i=0;i<prop.size();i++){
//      prop[i].resize(Nxyz);
//    }
//  }
//  if(clear){zeroE(prop);}
//}

//template<typename Ty, typename Ta>
//void copy_prop4d_to_propE(EigenMTa &prop, std::vector<Propagator4dT<Ty > > &pV1, qlat::fft_desc_basic &fd){
//  copy_propE(pV1, prop, fd, 0);
//}
//template<typename Ty, typename Ta>
//void copy_propE_to_prop4d(std::vector<Propagator4dT<Ty > > &pV1, EigenMTa &prop, qlat::fft_desc_basic &fd){
//  copy_propE(pV1, prop, fd, 1);
//}
//
//template<typename Ty, typename Ta>
//void copy_prop4d_to_propE(EigenMTa &prop, Propagator4dT<Ty > &pV1, qlat::fft_desc_basic &fd){
//  copy_propE(pV1, prop, fd, 0);
//}
//template<typename Ty, typename Ta>
//void copy_propE_to_prop4d(Propagator4dT<Ty > &pV1, EigenMTa &prop, qlat::fft_desc_basic &fd){
//  copy_propE(pV1, prop, fd, 1);
//}
//

//template<typename Ty>
//void ini_propG(EigenTy& prop, const Long nmass, size_t Nsize, bool clear = true){
//  if(Long(prop.size()) != nmass){prop.resize(nmass);}
//  for(unsigned long i=0;i<prop.size();i++){
//    if(prop[i].size() != Nsize){
//      prop[i].resize(Nsize);
//    }
//    else{
//      if(clear){prop[i].set_zero();}
//    }
//  }
//}

template <typename Ty >
void check_prop_size(EigenTy& prop, fft_desc_basic& fd){
  for(unsigned int i=0;i<prop.size();i++)
  {
    if(prop[0].size() != size_t(fd.Nvol)*12*12)
    {
      qmessage("Size of Prop wrong. \n");
      Qassert(false);
    }
  }
}

template <typename Ty >
void copy_qprop_to_propG(EigenTy& res, std::vector<qlat::FieldM<Ty, 12*12> >& src, const qlat::Geometry& geo, Int GPU = 1, Int dir = 1)
{
  Int nvec = 0;
  if(dir == 1){
    nvec = src.size();
    res.resize(nvec);
  }
  if(dir == 0){
    nvec = res.size();
    src.resize(0);src.resize(nvec);
    for(Int ni=0;ni<nvec;ni++){
      src[ni].init(geo);
      Qassert(res[ni].size() == size_t(12*12*geo.local_volume()));
    }
  }
  if(nvec == 0){return ;}
  
  for(Int ni=0;ni<nvec;ni++)
  {
    if(dir == 1){res[ni].copy_from((Ty*) qlat::get_data(src[ni]).data(), 12*12*geo.local_volume(), GPU);}
    if(dir == 0){res[ni].copy_to((Ty*) qlat::get_data(src[ni]).data(), GPU);}
  }
}

template <typename Ty >
void copy_propG_to_qprop(std::vector<qpropT >& res, EigenTy& src, const qlat::Geometry& geo, Int GPU = 1)
{
  copy_qprop_to_propG(src, res, geo, 0, GPU);
}

template <typename Ty >
void ini_resE(qlat::vector<Ty > &res, Int nmass, qlat::fft_desc_basic &fd){
  Int NTt  = fd.Nv[3];
  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  Int do_resize = 0;
  if((LInt) res.size() != (LInt) nmass*NTt * Nxyz){do_resize=1;}
  if(do_resize == 1)
  {
    res.resize(0);res.resize(nmass*NTt * Nxyz);
  }
  clear_qv(res);
}

inline Coordinate get_src_pos(std::string src_n, Coordinate& off_L, const Geometry& geo){
  std::string noi_name = ssprintf("%s",src_n.c_str()  );

  qlat::FieldM<Complexq, 1> noi;
  noi.init(geo);

  qmessage("Noise %s \n",noi_name.c_str());
  qlat::set_zero(noi);
  load_gwu_noi(noi_name.c_str(), noi);
  Coordinate pos;////qlat::vector<Int > off_L;
  check_noise_pos(noi, pos,off_L);

  return pos;
}

/////momentum related
inline std::vector<double >  hash_mom(const Coordinate& m){
  Qassert(m[3] == 0);
  std::vector<double > q;q.resize(4);
  for(Int i=0;i<4;i++){q[i] = 0;}

  for(Int i=0;i<3;i++)
  {
    q[0] += std::abs(m[i]);
    q[1] += m[i]*m[i];
    q[2] += m[i]*m[i]*m[i]*m[i];
    q[3] += m[i]*m[i]*m[i]*m[i]*m[i]*m[i];
  }
  return q;
}

////0 equal, -1 a<b, +1 a>b
inline Int compare_mom(const Coordinate& a, const Coordinate& b){
  std::vector<double > pa =  hash_mom(a);
  std::vector<double > pb =  hash_mom(b);
  Int equal =  0;
  for(unsigned int i=0;i<pa.size();i++){if(pa[i] != pb[i]){equal = -1;}}
  if(equal == 0){return equal;}

  std::vector<Int > cL = {1, 0, 2, 3};
  for(Int ci=0;ci<4;ci++)
  {
    if(pa[cL[ci]] <  pb[cL[ci]]){equal = -1;break;}
    if(pa[cL[ci]] >  pb[cL[ci]]){equal =  1;break;}
    ////equal then next ci
  }
  return equal;
}

/*
  phases for small number of momenta apply
  exp(sign * ( ( x - offset ) * mom) * 2pi/L)
*/
template<typename Ty>
void get_phases(std::vector<vector_gpu<Ty >>& phases, const std::vector<Coordinate >& momL,
            const Geometry& geo, const int8_t sign = 1, const Coordinate& offset = Coordinate() )
{
  TIMER("get_phases");
  Int Nmom = momL.size();
  phases.resize(Nmom);if(Nmom == 0){return ;}
  qlat::vector<Int > nv, Nv, mv;geo_to_nv(geo, nv, Nv, mv);
  Long vol = Nv[0]*Nv[1]*Nv[2];
  for(Int momi=0;momi<Nmom;momi++){phases[momi].resize(vol);}
  qlat::vector<Ty* > Pres = EigenM_to_pointers(phases);

  qlat::vector<double > p0;p0.resize(3);
  for(Int i=0;i<3;i++){p0[i] = 2 * QLAT_PI_LOCAL /nv[i];}

  ////copy momentum to gpu memery
  qlat::vector<Int > momLV;momLV.resize(momL.size() * 3);
  for(Int momi = 0;momi < Nmom; momi++)
  for(Int i=0;i<3;i++){momLV[momi*3 + i] = momL[momi][i];}
  Int* momLP = momLV.data();
  ////copy momentum to gpu memery

  qacc_for(xi, Long(vol),{
    const Coordinate xl  = geo.coordinate_from_index(xi);
    const Coordinate pos = geo.coordinate_g_from_l(xl);
    for(Int momi = 0;momi < Nmom; momi++)
    {
      double theta = 0.0;
      for(Int i=0;i<3;i++){theta += ( p0[i] * momLP[momi * 3 + i] * (pos[i] - offset[i]) );}
      Pres[momi][xi] = Ty(cos(theta), sign * sin(theta));
    }
  });
}

template<typename Ty>
void get_phases(vector_gpu<Ty >& phases, const Coordinate& mom,
            const Geometry& geo, const int8_t sign = 1, const Coordinate& offset = Coordinate() )
{
  TIMER("get_phases");
  std::vector<Coordinate > momL;momL.resize(1);momL[0] = mom;
  std::vector<vector_gpu<Ty > > phaseL;
  get_phases(phaseL, momL, geo, sign, offset);
  phases.copy_from(phaseL[0]);
}

template <class Ty, Int civ>
void apply_phases(qlat::FieldM<Ty, civ >& src, qlat::FieldM<Ty, civ>* res, vector_gpu<Ty >* phases, const unsigned int Nmom){
  Qassert(src.initialized);
  const Geometry& geo = src.geo();
  //const unsigned int Nmom = phases.size();

  //qlat::vector<Ty* > Pphase = EigenM_to_pointers(phases);
  qlat::vector<Ty* > Pphase;Pphase.resize(Nmom);
  for(unsigned int mi=0;mi<Nmom;mi++){Pphase[mi] = (Ty*) qlat::get_data(phases[mi]).data();}
  ////= EigenM_to_pointers(phases);
  Ty* psrc = (Ty*) qlat::get_data(src).data();

  qlat::vector<Ty* > pres;pres.resize(Nmom);
  /////if(res.size() != Nmom){res.resize(Nmom);}
  for(unsigned int mi=0;mi<Nmom;mi++){
    if(!res[mi].initialized){res[mi].init(geo);}
    pres[mi] = (Ty*) qlat::get_data(res[mi]).data();
  }

  qlat::vector<Int > nv, Nv, mv;geo_to_nv(geo, nv, Nv, mv);
  const Long vol = Nv[0]*Nv[1]*Nv[2];
  const Int   Nt = Nv[3];

  qacc_for(xi, Long(vol),{
    for(unsigned int momi = 0;momi < Nmom; momi++)
    {
      const Ty& ph   = Pphase[momi][xi];
      for(Int ti=0;ti<Nt;ti++)
      for(Int ci=0;ci<civ;ci++){
        pres[momi][(ti*vol+xi)*civ + ci] = psrc[(ti*vol+xi)*civ + ci] * ph;
      }
    }
  });
}

template <class Ty, Int civ>
void apply_phases(qlat::FieldM<Ty, civ >& src, std::vector<qlat::FieldM<Ty, civ>>& res, std::vector<vector_gpu<Ty >>& phases){
  Qassert(src.initialized);
  //const Geometry& geo = src.geo();
  const unsigned int Nmom = phases.size();

  if(res.size() != Nmom){res.resize(Nmom);}
  //for(unsigned int mi=0;mi<Nmom;mi++){
  //  if(!res[mi].initialized){res[mi].init(geo);}
  //  ////pres[mi] = (Ty*) qlat::get_data(res[mi]).data();
  //}
  apply_phases(src, &res[0], &phases[0], Nmom);
}

template <class Ty, Int civ>
void apply_phases(qlat::FieldM<Ty, civ >& src, qlat::FieldM<Ty, civ>& res, vector_gpu<Ty >& phases){
  Qassert(src.initialized);
  //const Geometry& geo = src.geo();
  apply_phases(src, &res, &phases, 1);
}

template <class Ty, class Td>
void apply_phases(Propagator4dT<Td>& src, Propagator4dT<Td>& res, qlat::vector_gpu<Ty >& phases){

  Qassert(src.initialized);
  const Geometry& geo = src.geo();
  if(!res.initialized){res.init(geo);}

  Ty* Pphase = (Ty*) qlat::get_data(phases).data();
  Ty* psrc = (Ty*) qlat::get_data(src).data();
  Ty* pres = (Ty*) qlat::get_data(res).data();

  qlat::vector<Int > nv, Nv, mv;geo_to_nv(geo, nv, Nv, mv);
  const Long vol = Nv[0]*Nv[1]*Nv[2];
  const Int   Nt = Nv[3];
  const Int civ = 12 * 12;

  qacc_for(xi, Long(vol),{
    {
      const Ty& ph   = Pphase[xi];
      for(Int ti=0;ti<Nt;ti++)
      for(Int ci=0;ci<civ;ci++){
        pres[(ti*vol+xi)*civ + ci] = psrc[(ti*vol+xi)*civ + ci] * ph;
      }
    }
  });

}

//template <class Ty, Int civ>
//void apply_phases(qlat::FieldM<Ty, civ >& src, std::vector<qlat::FieldM<Ty, civ>>& res, std::vector<vector_gpu<Ty >>& phases){
//  Qassert(src.initialized);
//  const Geometry& geo = src.geo();
//  const unsigned int Nmom = phases.size();
//
//  qlat::vector<Ty* > Pphase = EigenM_to_pointers(phases);
//  Ty* psrc = (Ty*) qlat::get_data(src).data();
//
//  qlat::vector<Ty* > pres;pres.resize(Nmom);
//  if(res.size() != Nmom){res.resize(Nmom);}
//  for(unsigned int mi=0;mi<Nmom;mi++){
//    if(!res[mi].initialized){res[mi].init(geo);}
//    pres[mi] = (Ty*) qlat::get_data(res[mi]).data();
//  }
//
//  qlat::vector<Int > nv, Nv, mv;geo_to_nv(geo, nv, Nv, mv);
//  const Long vol = Nv[0]*Nv[1]*Nv[2];
//  const Int   Nt = Nv[3];
//
//  qacc_for(xi, Long(vol),{
//    for(Int momi = 0;momi < Nmom; momi++)
//    {
//      const Ty& ph   = Pphase[momi][xi];
//      for(Int ti=0;ti<Nt;ti++)
//      for(Int ci=0;ci<civ;ci++){
//        pres[momi][(ti*vol+xi)*civ + ci] = psrc[(ti*vol+xi)*civ + ci] * ph;
//      }
//    }
//  });
//}

template <class Ty, Int civ>
void multi_factor(qlat::FieldM<Ty, civ >& src, qlat::FieldM<Ty, civ>& res, const Ty& phases){
  Qassert(src.initialized);
  const Geometry& geo = src.geo();
  if(!res.initialized){res.init(geo);}

  const Ty* psrc = (Ty*) qlat::get_data(src).data();
  Ty* pres = (Ty*) qlat::get_data(res).data();

  qacc_for(xi, geo.local_volume(),{
    for(Int ci=0;ci<civ;ci++){
      res[xi*civ + ci] = src[xi*civ + ci] * phases;
    }
  });
}

////phases for momentum data
template<typename Ty>
void get_phases(std::vector<Ty >& phases, Coordinate& pL, const Coordinate& src, const Coordinate& Lat){
  Long vol = pL[0]*pL[1]*pL[2];
  phases.resize(vol);
  for(Int i=0;i<3;i++){Qassert(pL[i] % 2 == 1);}
  #pragma omp parallel for
  for(Long isp =0;isp<vol;isp++){
    Coordinate pos = qlat::coordinate_from_index(isp, pL);

    for(Int i=0;i<3;i++){
      if(pos[i] > pL[i]/2){
        pos[i] = Lat[i] - (pL[i] - pos[i]);
      }
    }
    double v0 = 0.0;
    for(Int i=0;i<3;i++){v0 += (2.0 * QLAT_PI_LOCAL * src[i] * pos[i]/Lat[i]);}

    phases[isp] = Ty(std::cos(v0), -1.0* std::sin(v0));
  }
}

/////V -- 12a x 12b   to   12b x 12a -- V
template<typename Ty>
void copy_qprop_to_propE(std::vector<qlat::vector<Ty > >& Eprop, std::vector<qpropT >& src, Int dir = 1){
  TIMERA("copy_qprop_to_propE");
  const Int nmass = src.size();
  std::vector<Ty* > ps;ps.resize(nmass);
  
  for(Int mi=0;mi<nmass;mi++){
    Qassert(src[mi].initialized);
    ps[mi] = (Ty*) qlat::get_data(src[mi]).data();
  }

  const qlat::Geometry& geo = src[0].geo();
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);
  if(dir == 1){ini_propE(Eprop, nmass, fd);}
  
  ///const Long sizeF = geo.local_volume();
  const Long nvec  = Eprop.size()/nmass;
  const Long sizeF = Eprop[0].size();
  Qassert(nvec * sizeF == 12*12*geo.local_volume());

  ////V x 12 a x 12 b to 12b x 12a x V
  for(Int mi=0;mi<nmass;mi++)
  for(Long i=0;i<nvec;i++)
  {
    if(dir == 1)cpy_data_thread(Eprop[mi*nvec + i].data(), &ps[mi][i*sizeF], sizeF, 1, false);
    if(dir == 0)cpy_data_thread(&ps[mi][i*sizeF], Eprop[mi*nvec + i].data(), sizeF, 1, false);
  }
  qacc_barrier(dummy);
}

template<typename Ty>
void copy_propE_to_qprop(std::vector<qpropT >& src, std::vector<qlat::vector<Ty > >& Eprop){
  copy_qprop_to_propE(Eprop, src, 0);
}

template<typename Ty>
void noise_to_propT(qpropT& prop, qnoiT& noi){
  Qassert(noi.initialized);
  const Geometry& geo = noi.geo();

  if(!prop.initialized){prop.init(geo);}

  Ty* res = (Ty*) qlat::get_data(prop).data();
  Ty* src = (Ty*) qlat::get_data(noi ).data();

  const Long Nvol = geo.local_volume();

  for(Int d0=0;d0<12;d0++){
    cpy_data_thread(&res[(d0*12+d0)*Nvol + 0], src, Nvol, 1, QFALSE);
  }
  qacc_barrier(dummy);

}

template <class Td>
void prop4D_factor(Propagator4dT<Td>& prop, const qlat::ComplexT<Td >& factor)
{
  const size_t Nvol = size_t(prop.geo().local_volume()) * 12 * 12;
  qlat::ComplexT<Td >* src = (qlat::ComplexT<Td >*) qlat::get_data(prop).data();
  qacc_for(isp, Long(Nvol), {
    src[isp] *= factor;
  })
  ////if(qlat::qnorm(factor) >= QLAT_COPY_LIMIT){cpy_data_threadC(src, src, Nvol, 1, true, factor);}
  ////else{;}
}



template <class Ty, Int civ>
void copy_FieldM(std::vector<qlat::FieldM<Ty, civ> >& res, const std::vector<qlat::FieldM<Ty, civ> >& src)
{
  if(src.size() == 0){res.resize(0); return;}
  const Int Nv = src.size();
  for(Int iv=0;iv<Nv;iv++){Qassert(src[iv].initialized);}
  if(res.size() != Nv){res.resize(Nv);}

  const Geometry& geo = src[0].geo();
  const size_t Nvol = size_t(geo.local_volume()) * civ;
  for(Int iv=0;iv<Nv;iv++)
  {
    if(!res[iv].initialized){res[iv].init(geo);}
    Ty* resP = (Ty*) qlat::get_data(res[iv]).data();
    Ty* srcP = (Ty*) qlat::get_data(src[iv]).data();
    cpy_data_thread(resP, srcP, Nvol, 1, false);
  }
  qacc_barrier(dummy);
}


}

#endif
