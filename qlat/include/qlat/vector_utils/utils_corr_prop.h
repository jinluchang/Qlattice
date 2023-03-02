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


template<typename Td, int dir, bool conj>
void prop4d_src_gammaT(Propagator4dT<Td >& prop, ga_M& ga){
  TIMERA("prop4d_src_gamma");
  ////Rowmajor (a,b), b is continues in memory
  qacc_for(isp, long(prop.geo().local_volume()),{
    qlat::WilsonMatrixT<Td>& v0 =  prop.get_elem_offset(isp);
    qlat::WilsonMatrixT<Td>  v1 = v0;

    for(int s = 0; s < 4; ++s)
    for(int c0 = 0;c0< 3 ; c0++)
    for(int d0 = 0; d0 < 4; ++d0)
    {
      /////Source multiply
      if(dir==0)for(int c1=0;c1<3;c1++){
        if(!conj){v0(s*3 + c0, ga.ind[d0]*3 + c1) = ga.g[d0] * v1(s*3 + c0, d0*3 + c1);}
        if( conj){v0(s*3 + c0, ga.ind[d0]*3 + c1) = qlat::qconj(ga.g[d0] * v1(s*3 + c0, d0*3 + c1));}
      }
      /////Sink multiply
      if(dir==1)for(int c1=0;c1<3;c1++){
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

template<typename Ty, int dir, bool conj>
void qprop_src_gamma_T(Ty* res, ga_M& ga, const long Nsize){
  TIMERA("qprop_src_gamma");
  ////Rowmajor (a,b), b is continues in memory
  ///Ty* res = (Ty*) qlat::get_data(prop1).data()
  //const long Nsize = prop.geo().local_volume();
  qacc_for(isp, Nsize,{
    Ty src[12*12];
    Ty buf[12*12];
    for(int i=0;i<12*12;i++){src[i] = res[i*Nsize + isp];}

    for(int d1 = 0;d1 < 4; d1++)
    for(int c1 = 0;c1 < 3; c1++)
    for(int d0 = 0;d0 < 4; d0++)
    for(int c0 = 0;c0 < 3; c0++)
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

    for(int i=0;i<12*12;i++){
      if(!conj){res[i*Nsize + isp] = buf[i];}
      if( conj){res[i*Nsize + isp] = qlat::qconj(buf[i]);}
    }
  });
}

template<typename Ty>
void qprop_src_gamma(qpropT& prop, ga_M& ga, bool conj = false){
  Ty* res = (Ty*) qlat::get_data(prop).data();
  const long Nsize = prop.geo().local_volume();
  if(!conj)qprop_src_gamma_T<Ty, 0, false>(res, ga, Nsize);
  if( conj)qprop_src_gamma_T<Ty, 0, true >(res, ga, Nsize);
}


template<typename Ty>
void qprop_sink_gamma(qpropT& prop, ga_M& ga, bool conj = false){
  Ty* res = (Ty*) qlat::get_data(prop).data();
  const long Nsize = prop.geo().local_volume();
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

template<typename Ty>
void qprop_move_dc_in(Ty* src, const qlat::Geometry &geo, const int dir = 1)
{
  move_index mv_civ;
  const long sizeF = geo.local_volume();

  if(dir == 1){mv_civ.move_civ_in( src, src, 1, 12*12, sizeF, 1, true);}
  if(dir == 0){mv_civ.move_civ_out(src, src, 1, sizeF, 12*12, 1, true);}
}

template<typename Ty>
void qprop_move_dc_in(qpropT& src, const int dir = 1)
{
  qassert(src.initialized);
  qprop_move_dc_in((Ty*) qlat::get_data(src).data(), src.geo(), dir);
}

template<typename Ty>
void qprop_move_dc_out(Ty* src, const qlat::Geometry &geo)
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
  const int Nvec = s0.size();
  const qlat::Geometry &geo = s0[0].geo();
  const long Nvol = geo.local_volume();

  init_qpropT(res, Nvec, geo);
  for(int vi=0;vi<Nvec;vi++)
  {
    Ty* p0 = (Ty* ) qlat::get_data(s0[vi]).data();
    Ty* r0 = (Ty* ) qlat::get_data(res[vi]).data(); 
    for(int dc=0;dc<12*12;dc++){
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
  qassert(s0.size() == s1.size());
  const int Nvec = s0.size();
  const qlat::Geometry &geo = s0[0].geo();
  const long Nvol = geo.local_volume();

  init_qpropT(res, Nvec, geo);
  for(int vi=0;vi<Nvec;vi++)
  {
    Ty* p0 = (Ty* ) qlat::get_data(s0[vi]).data();
    Ty* p1 = (Ty* ) qlat::get_data(s1[vi]).data(); 
    Ty* r0 = (Ty* ) qlat::get_data(res[vi]).data(); 
    for(int dc=0;dc<12*12;dc++){
      qacc_for(isp, geo.local_volume(),{
        r0[dc*Nvol + isp] = (p0[dc*Nvol + isp] + p1[dc*Nvol + isp] * f0) * f1;
      });
    }
  }
}

template<typename Ty>
void shift_vecs_cov_qpropT(std::vector< std::vector<qpropT > >& res, std::vector< qpropT >& s0, shift_vec& svec,
  std::vector<std::vector<qpropT >>& buf)
{
  if(res.size() != 5){res.resize(5);}
  if(buf.size() != 2){buf.resize(2);}
  qprop_sub_add(res[0], s0, Ty(0.0,0.0), Ty(1.0,0.0) );////equal

  init_qpropT(buf[0], s0.size(), s0[0].geo());
  init_qpropT(buf[1], s0.size(), s0[0].geo());

  for(int nu = 0; nu < 4 ; nu++)
  {
    shift_vecs_dir_qpropT(s0, buf[0], nu, +1, svec);
    shift_vecs_dir_qpropT(s0, buf[1], nu, -1, svec);
    qprop_sub_add(res[1 + nu], buf[0], buf[1], Ty(-1.0,0.0), Ty(1.0/2.0, 0.0) );////equal
  }
}


template<typename Td>
void prop4d_cps_to_ps(Propagator4dT<Td >& prop, int dir=0){
  /////sn is -1 for default
  qlat::ComplexT<Td > sn =-1;if(dir == 1){sn= 1;}
  const qlat::ComplexT<Td>sqrt2= qlat::ComplexT<Td>(std::sqrt(2.0), 0.0);

  ////Rowmajor (a,b), b is continues in memory
  qacc_for(isp, prop.geo().local_volume(),{
    qlat::WilsonMatrixT<Td >  v0 = prop.get_elem_offset(isp);
    qlat::WilsonMatrixT<Td >  v1 = prop.get_elem_offset(isp);

    int dr,d0,d1;
    /////Src rotation
    for (int s0 = 0; s0 < 4; ++s0)
    for(int c0 = 0;c0< 3 ; c0++)
    {
      dr=0;d0=1;d1=3;
      for(int c1=0;c1<3;c1++)v0(s0*3+c0, dr*3+c1) = (-v1(s0*3+c0, d0*3+ c1)*sn + v1(s0*3+c0, d1*3+c1)   )/sqrt2;
      dr=1;d0=0;d1=2;
      for(int c1=0;c1<3;c1++)v0(s0*3+c0, dr*3+c1) = (+v1(s0*3+c0, d0*3+ c1)*sn - v1(s0*3+c0, d1*3+c1)   )/sqrt2;
      dr=2;d0=1;d1=3;
      for(int c1=0;c1<3;c1++)v0(s0*3+c0, dr*3+c1) = (-v1(s0*3+c0, d0*3+ c1)    - v1(s0*3+c0, d1*3+c1)*sn)/sqrt2;
      dr=3;d0=0;d1=2;
      for(int c1=0;c1<3;c1++)v0(s0*3+c0, dr*3+c1) = (+v1(s0*3+c0, d0*3+ c1)    + v1(s0*3+c0, d1*3+c1)*sn)/sqrt2;
    }

    /////Copy previous results
    v1 = v0;
    /////Sink rotation
    for(int c0 = 0;c0< 3 ; c0++)
    for (int s0 = 0; s0 < 4; ++s0)
    {
      dr=0;d0=1;d1=3;
      for(int c1=0;c1<3;c1++)v0(dr*3+c0, s0*3+c1) = (-v1(d0*3+c0, s0*3+c1)*sn + v1(d1*3+c0, s0*3+c1)   )/sqrt2;
      dr=1;d0=0;d1=2;
      for(int c1=0;c1<3;c1++)v0(dr*3+c0, s0*3+c1) = ( v1(d0*3+c0, s0*3+c1)*sn - v1(d1*3+c0, s0*3+c1)   )/sqrt2;
      dr=2;d0=1;d1=3;
      for(int c1=0;c1<3;c1++)v0(dr*3+c0, s0*3+c1) = (-v1(d0*3+c0, s0*3+c1)    - v1(d1*3+c0, s0*3+c1)*sn)/sqrt2;
      dr=3;d0=0;d1=2;
      for(int c1=0;c1<3;c1++)v0(dr*3+c0, s0*3+c1) = ( v1(d0*3+c0, s0*3+c1)    + v1(d1*3+c0, s0*3+c1)*sn)/sqrt2;
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

  const qlat::Geometry &geo = prop[0].geo();

  unsigned long Nvol = geo.local_volume();
  ///int Nt = geo.node_site[3];
  ///long Nsum = Nvol/Nt;
  int tini = x_ini[3];

  qlat::vector_acc<qlat::ComplexT<Td> > res;res.resize(Nvol);

  qacc_for(isp, long(Nvol),{
    qlat::ComplexT<Td> buf(0.0,0.0);

    for(int dc2=0;dc2<12;dc2++){
      qlat::ComplexT<Td>* a = (qlat::ComplexT<Td>* ) &(prop[dc2].get_elem_offset(isp));
      for(int dc1=0;dc1<12;dc1++)
      {
        buf+=a[dc1]*qlat::qconj(a[dc1]);
      }
    }
    res[isp] = buf;
    ////src[isp] = buf[isp];
  });

  const Coordinate vg = geo.total_site();
  int nt = vg[3];
  write.resize(0);
  write.resize(2*nt);qlat::set_zero(write);

  for(unsigned long isp=0;isp<Nvol;isp++){
    Coordinate xl0 = geo.coordinate_from_index(isp);
    Coordinate xg0 = geo.coordinate_g_from_l(xl0);
    int t = xg0[3];

    int toff = ((t-tini+nt)%nt);
    write[ toff*2 + 0 ] += res[isp].real();
    write[ toff*2 + 1 ] += res[isp].imag();
  }
  ////May need to be changed for EigenV
  //sum_all_size((double*) &write[0],2*nt);
  sum_all_size((double*) write.data(), 2*nt);
}

template<typename Ty>
void get_src_phase(Ty& phase, const qlat::vector_acc<int >& nv,
  const Coordinate& pos = Coordinate(), const Coordinate& mom = Coordinate()){
    double p0[3]={2*PI/nv[0],
                  2*PI/nv[1],
                  2*PI/nv[2]};
    double theta= mom[0]*p0[0]*pos[0] + mom[1]*p0[1]*pos[1] + mom[2]*p0[2]*pos[2];
    phase = Ty(std::cos(theta), -1.0*std::sin(theta));
}

template<typename Ty>
void vec_corrE(Ty* srcE, qlat::vector_acc<Ty >& res,qlat::fft_desc_basic &fd,const int nvec,const int clear=0,const Coordinate& mom = Coordinate(), const Ty& src_phase = 1.0, const int t0 = 0){
  TIMER("Reduce vec_corrE");
  int NTt  = fd.Nv[3];
  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  Ty* src = srcE;
  //int nmass = resE.size()/(NTt*Nxyz);

  ////position p = fd.desc.get_position(0,fd.rank);
  int t_rank = fd.Pos0[fd.rank][3];
  qlat::vector_gpu<Ty > bufE;
  if(mom != Coordinate())
  {
    bufE.resize(nvec * Nxyz*NTt);
    cpy_data_thread(bufE.data(), srcE, bufE.size(), 1, true);
    src = bufE.data();

    qlat::vector_acc<double > p0;p0.resize(3);
    for(int i=0;i<3;i++){p0[i] = 2*3.1415926535898/fd.nv[i];}

    qlat::vector_acc<Ty > phaseEG;phaseEG.resize(Nxyz);
    //Ty* phaseE = (Ty*) qlat::get_data(phaseEG).data();

    //const qlat::vector_acc<int >& orderN = fd.orderN;
    //const qlat::vector_acc<int >& Nv = fd.Nv;

    ////qlat::vector_gpu<int > pos_tem;pos_tem.copy_from(fd.Pos0[fd.rank]);int* posP = pos_tem.data();
    /////===may not be consistent for fd definiations under qacc
    /////===slow
    qthread_for(xi, long(Nxyz),{
      //int pi[3];
      //pi[orderN[0]] = xi/(Nv[orderN[1]]*Nv[orderN[2]]);
      //pi[orderN[1]] = (xi%(Nv[orderN[1]]*Nv[orderN[2]]))/Nv[orderN[2]];
      //pi[orderN[2]] = xi%Nv[orderN[2]];
      //for(int ptem=0;ptem<3;ptem++){pi[ptem] = pi[ptem] + posP[ptem];}
      Coordinate pi = fd.coordinate_g_from_index(xi );

      double theta=mom[0]*p0[0]*pi[0]+mom[1]*p0[1]*pi[1]+mom[2]*p0[2]*pi[2];
      phaseEG[xi] = Ty(cos(theta),sin(theta));
    });

    size_t Ns = Nxyz;
    Ns = nvec*NTt*Nxyz;
    qacc_for(i, long(Ns),{
      LInt mi = i/(NTt*Nxyz);
      LInt ti = (i%(NTt*Nxyz))/Nxyz;
      LInt xi = i%(Nxyz);
      src[(mi*NTt + ti)*Nxyz + xi] = src[(mi*NTt + ti)*Nxyz + xi]*phaseEG[xi];
    });
    /////#endif
  }

  ////TODO
  int nt = fd.nt;
  //if(clear == 1){res.resize(0);res.resize(nvecs*nt);qlat::set_zero(res);}
  if(clear == 1){if(res.size() != nvec*nt){res.resize(nvec*nt);} qlat::set_zero(res);}
  if(clear == 0){if(res.size() != nvec*nt){print0("res size wrong for corr.");qassert(false);}}

  qlat::vector_acc<Ty > tmp;tmp.resize(nvec*NTt);qlat::set_zero(tmp);//tmp.set_zero();
  reduce_vec(src, tmp.data(), Nxyz, nvec*NTt);

  qlat::vector_gpu<Ty > RES;RES.resize(nvec*nt );RES.set_zero();
  Ty* s1 = RES.data();Ty* s0 = tmp.data();
  long Ntotal = nvec*NTt;
  qacc_for(mti, Ntotal, {
    long mi  = mti/NTt;
    long  ti = mti%NTt;
    s1[mi*nt + (t_rank + ti + nt - t0)%nt ] = s0[mi*NTt + ti] * src_phase;
  });

  //////sum_all_size((Ftype*) (RES.data()), 2*RES.size(), 1);
  sum_all_size(RES.data(), RES.size(), 1);
  cpy_data_thread((Ty*) qlat::get_data(res).data(), RES.data(), RES.size(), 1, true, 1.0);
}

template<typename Ty>
void vec_corrE(qlat::vector_gpu<Ty >& resE, qlat::vector_acc<Ty >& res,qlat::fft_desc_basic &fd,const int clear=0,
  const Coordinate& mom = Coordinate(), const Ty& src_phase = 1.0){
  int NTt  = fd.Nv[3];
  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  int nvec = resE.size()/(NTt*Nxyz);
  vec_corrE(resE.data(), res, fd, nvec, clear, mom, src_phase);
}

template<typename Ty>
void vec_corrE(qlat::vector_acc<Ty >& resE, qlat::vector_acc<Ty >& res,qlat::fft_desc_basic &fd,const int clear=0,
  const Coordinate& mom = Coordinate(), const Ty& src_phase = 1.0){
  int NTt  = fd.Nv[3];
  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  int nvec = resE.size()/(NTt*Nxyz);
  ///qlat::vector_gpu<Ty > r0;r0.copy_from(resE);
  Ty* r0 = (Ty*) qlat::get_data(resE).data();
  vec_corrE(r0, res, fd, nvec, clear, mom, src_phase);
}

void write_pos_to_string(std::string& POS_LIST, Coordinate& pos){
  std::string buf;
  char pnum[500];
  sprintf(pnum," ");
  for(int i=0;i<4;i++){
    sprintf(pnum,"%d ", pos[i]);
    buf += std::string(pnum);
  }
  buf += std::string(" ; ");
  POS_LIST += buf;
}

inline std::vector<Coordinate > string_to_coord(std::string& INFO){
  std::vector<Coordinate > posL ;
  std::vector<std::string > a = stringtolist(INFO);
  qassert(a.size() % 5 == 0);
  int Npos = a.size()/5;
  for(int i=0;i< Npos;i++)
  {
    qassert(a[i*5+4] == std::string(";"));
    Coordinate c;
    for(int j = 0;j<4;j++){c[j] = stringtonum(a[i*5 + j]);}
    posL.push_back(c);
  }
  return posL;
}

template<typename Ty>
void shift_result_t(qlat::vector_acc<Ty >& Esrc, int nt, int tini){
  if(tini == 0){return ;}
  long Ntotal = Esrc.size();
  if(Ntotal %(nt) != 0){abort_r("Correlation function size wrong!\n");}
  qlat::vector_acc<Ty > tmp;tmp.resize(Ntotal);
  qacc_for(i, Ntotal, {
    const int iv = i/nt;
    const int t  = i%nt;
    tmp[iv*nt + (t - tini + nt)%nt] = Esrc[iv*nt + (t)%nt];
  });
  cpy_data_thread(Esrc.data(), tmp.data(), tmp.size(), 1);
}

//template<typename Ta>
//void ini_propE(EigenMTa &prop,int nmass, qlat::fft_desc_basic &fd, bool clear = true){
//  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
//  int NTt  = fd.Nv[3];
//  int do_resize = 0;
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

//////default dir == 0 from prop4d to propE
//template<typename Ty, typename Ta>
//void copy_propE(std::vector<Ty* > &pV1, EigenMTa &prop, qlat::fft_desc_basic &fd, int dir=0){
//  TIMER("Copy prop");
//  qassert(fd.order_ch == 0);
//  const LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
//  const int NTt  = fd.Nv[3];
//  Geometry geo;fd.get_geo(geo);
//  int nmass = 0;
//  if(dir==0){nmass = pV1.size();ini_propE(prop,nmass,fd);}
//  if(dir==1){
//    nmass = prop.size()/(12*12*NTt);
//    qassert(int(pV1.size()) == nmass);
//    //if(pV1.size() != (LInt) nmass){
//    //  pV1.resize(0);
//    //  pV1.resize(nmass);for(int i=0;i<nmass;i++){pV1[i].init(geo);}
//    //}
//  }
//
//  for(int mi = 0;mi < nmass;mi++)
//  {
//    Ty* pv = pV1[mi];
//    /////Propagator4dT<Ty >& pv = pV1[mi];
//    qlat::vector_acc<Ta* > ps = EigenM_to_pointers(prop);
//    qacc_for(isp, long(NTt*Nxyz),{
//      int ti = isp/Nxyz;
//      int xi = isp%Nxyz;
//      /////qlat::WilsonMatrixT<Ty>& v0 =  pv.get_elem_offset(isp);
//      Ty* v0 = &pv[isp * 12 * 12];
//
//      for(int c0 = 0;c0 < 3; c0++)
//      for(int d0 = 0;d0 < 4; d0++)
//      for(int c1 = 0;c1 < 3; c1++)
//      for(int d1 = 0;d1 < 4; d1++)
//      {
//        LInt off = mi*12*12 + (d1*3+c1)*12+d0*3+c0;
//        if(dir==0){ps[off*NTt+ti][xi] = v0[(d0*3 + c0)*12 +  d1*3 + c1];}
//        if(dir==1){v0[(d0*3 + c0)*12 +  d1*3 + c1] = ps[off*NTt+ti][xi];}
//      }
//    });
//  }
//
//}

//template<typename Ty, typename Ta>
//void copy_propE(Propagator4dT<Ty > &pV, EigenMTa &prop, qlat::fft_desc_basic &fd, int dir=0){
//  int nmass = 1;
//  if(dir==0){
//    qassert(pV.initialized);
//    ini_propE(prop,nmass,fd);
//  }
//  if(dir==1){
//    nmass = prop.size()/(12*12*fd.Nv[3]);
//    qassert(nmass == 1);
//    if(!pV.initialized){
//      Geometry geo;fd.get_geo(geo);
//      pV.init(geo);
//    }
//  }
//  std::vector<Ty* > PTy;PTy.resize(1);
//  PTy[0] = (Ty*) (qlat::get_data(pV).data());
//  copy_propE(PTy, prop, fd, dir);
//}

//template<typename Ty, typename Ta>
//void copy_propE(std::vector<Propagator4dT<Ty > > &pV, EigenMTa &prop, qlat::fft_desc_basic &fd, int dir=0){
//  int nmass = 0;
//  if(dir==0){nmass = pV.size();ini_propE(prop,nmass,fd);}
//  if(dir==1){
//    nmass = prop.size()/(12*12*fd.Nv[3]);
//    bool do_init = false;
//    if(pV.size() != (LInt) nmass){do_init = true;}
//    if(pV.size() > 0){if(!pV[0].initialized){do_init = true;}}
//    if(do_init)
//    {
//      pV.resize(0);
//      Geometry geo;fd.get_geo(geo);
//      pV.resize(nmass);for(int i=0;i<nmass;i++){pV[i].init(geo);}
//    }
//  }
//  std::vector<Ty* > PTy;PTy.resize(nmass);
//  for(int im=0;im<nmass;im++){PTy[im] = (Ty*) (qlat::get_data(pV[im]).data());}
//  copy_propE(PTy, prop, fd, dir);
//}
//
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

template<typename Ty>
void ini_propG(EigenTy& prop, const long nmass, size_t Nsize, bool clear = true){
  if(long(prop.size()) != nmass){prop.resize(nmass);}
  for(unsigned long i=0;i<prop.size();i++){
    if(prop[i].size() != Nsize){
      prop[i].resize(Nsize);
    }
    else{
      if(clear){prop[i].set_zero();}
    }
  }
}

template <typename Ty >
void check_prop_size(EigenTy& prop, fft_desc_basic& fd){
  for(unsigned int i=0;i<prop.size();i++)
  {
    if(prop[0].size() != size_t(fd.Nvol)*12*12)
    {
      print0("Size of Prop wrong. \n");
      qassert(false);
    }
  }
}

template <typename Ty >
void copy_qprop_to_propG(EigenTy& res, std::vector<qpropT >& src, const qlat::Geometry &geo, int GPU = 1, int dir = 1)
{
  int nvec = 0;
  if(dir == 1){
    nvec = src.size();
    res.resize(nvec);
  }
  if(dir == 0){
    nvec = res.size();
    src.resize(nvec);
    for(int ni=0;ni<nvec;ni++){
      src[ni].init(geo);
      qassert(res[ni].size() == size_t(12*12*geo.local_volume()));
    }
  }
  if(nvec == 0){return ;}
  
  for(int ni=0;ni<nvec;ni++)
  {
    if(dir == 1){res[ni].copy_from((Complexq*) qlat::get_data(src[ni]).data(), 12*12*geo.local_volume(), GPU);}
    if(dir == 0){res[ni].copy_to((Complexq*) qlat::get_data(src[ni]).data(), GPU);}
  }
}

template <typename Ty >
void copy_propG_to_qprop(std::vector<qpropT >& res, EigenTy& src, const qlat::Geometry &geo, int GPU = 1)
{
  copy_qprop_to_propG(src, res, geo, 0, GPU);
}

template <typename Ty >
void ini_resE(qlat::vector_acc<Ty > &res, int nmass, qlat::fft_desc_basic &fd){
  int NTt  = fd.Nv[3];
  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  int do_resize = 0;
  if((LInt) res.size() != (LInt) nmass*NTt * Nxyz){do_resize=1;}
  if(do_resize == 1)
  {
    res.resize(0);res.resize(nmass*NTt * Nxyz);
  }
  clear_qv(res);
}

inline std::vector<int >  get_map_sec(int dT,int nt){
  std::vector<int > map_sec;map_sec.resize(nt);
  int secN = 2*nt/dT;double lensec = nt/(1.0*secN);
  int tcount = 0;
  int t0 = 0;
  for(int si=0;si<secN;si++)
  {
    for(int t=t0;t < (si+1)*lensec;t++)///boundary with the same sector?
    {
      qassert(t < nt);
      map_sec[t] = si;
      tcount = tcount + 1;
    }
    t0 = tcount;
  }
  return map_sec;

}

template <typename Ty>
void get_num_time(qlat::FieldM<Ty, 1>& noise,int &number_t, int &t_ini){
  qlat::Geometry& geo = noise.geo();
  qlat::vector_acc<int > nv,Nv,mv;
  geo_to_nv(geo, nv, Nv, mv);
  //int nx,ny,nz,nt;
  //nx = nv[0];ny = nv[1];nz = nv[2];nt = nv[3];
  int nt = nv[3];
  LInt Nsite = Nv[0]*Nv[1]*Nv[2]*Nv[3];

  std::vector<double > fullt(nt);for(int ti=0;ti<nt;ti++){fullt[ti]=0.0;}
  for(unsigned int isp=0; isp< Nsite; isp++)
  {
    ////position p = noise.desc->get_position(isp,get_node_rank());
    Coordinate xl0 = geo.coordinate_from_index(isp);
    Coordinate xg0 = geo.coordinate_g_from_l(xl0);
    {
      ///auto tem_source = noise.data[isp];
      auto tem_source =  noise.get_elem_offset(isp);
      if(qnorm(tem_source)>0.01)
      {
        fullt[xg0[3]] = 1.0;
      }
    }
  }
  sum_all_size((double* ) &fullt[0],nt);
  number_t = 0;
  for(int ti=0;ti<nt;ti++){if(fullt[ti]>0.0)number_t += 1;}
  for(int ti=0;ti<nt;ti++){if(fullt[ti]>0.0){t_ini = ti;break;}}
}

inline Coordinate get_src_pos(std::string src_n, Coordinate& off_L, const Geometry &geo){
  char noi_name[500];
  sprintf(noi_name ,"%s",src_n.c_str()  );

  qlat::FieldM<Complexq,1> noi;
  noi.init(geo);

  print0("Noise %s \n",noi_name);
  qlat::set_zero(noi);
  load_gwu_noi(noi_name,noi);
  Coordinate pos;////qlat::vector<int > off_L;
  check_noise_pos(noi, pos,off_L);

  return pos;
}

/////momentum related
inline std::vector<double >  hash_mom(const Coordinate& m){
  qassert(m[3] == 0);
  std::vector<double > q;q.resize(4);
  for(int i=0;i<4;i++){q[i] = 0;}

  for(int i=0;i<3;i++)
  {
    q[0] += std::abs(m[i]);
    q[1] += m[i]*m[i];
    q[2] += m[i]*m[i]*m[i]*m[i];
    q[3] += m[i]*m[i]*m[i]*m[i]*m[i]*m[i];
  }
  return q;
}

////0 equal, -1 a<b, +1 a>b
inline int compare_mom(const Coordinate& a, const Coordinate& b){
  std::vector<double > pa =  hash_mom(a);
  std::vector<double > pb =  hash_mom(b);
  int equal =  0;
  for(unsigned int i=0;i<pa.size();i++){if(pa[i] != pb[i]){equal = -1;}}
  if(equal == 0){return equal;}

  std::vector<int > cL = {1, 0, 2, 3};
  for(int ci=0;ci<4;ci++)
  {
    if(pa[cL[ci]] <  pb[cL[ci]]){equal = -1;break;}
    if(pa[cL[ci]] >  pb[cL[ci]]){equal =  1;break;}
    ////equal then next ci
  }
  return equal;
}

////phases for small number of momenta apply
template<typename Ty>
void get_phases(std::vector<vector_gpu<Ty >>& phases, const std::vector<Coordinate >& momL,
            const Geometry& geo, const char sign = 1, const Coordinate& offset = Coordinate() )
{
  TIMER("get_phases");
  int Nmom = momL.size();
  phases.resize(Nmom);if(Nmom == 0){return ;}
  qlat::vector_acc<int > nv, Nv, mv;geo_to_nv(geo, nv, Nv, mv);
  long vol = Nv[0]*Nv[1]*Nv[2];
  for(int momi=0;momi<Nmom;momi++){phases[momi].resize(vol);}
  qlat::vector_acc<Ty* > Pres = EigenM_to_pointers(phases);

  qlat::vector_acc<double > p0;p0.resize(3);
  for(int i=0;i<3;i++){p0[i] = 2*PI/nv[i];}

  ////copy momentum to gpu memery
  qlat::vector_acc<int > momLV;momLV.resize(momL.size() * 3);
  for(int momi = 0;momi < Nmom; momi++)
  for(int i=0;i<3;i++){momLV[momi*3 + i] = momL[momi][i];}
  int* momLP = momLV.data();
  ////copy momentum to gpu memery

  qacc_for(xi, long(vol),{
    const Coordinate xl  = geo.coordinate_from_index(xi);
    const Coordinate pos = geo.coordinate_g_from_l(xl);
    for(int momi = 0;momi < Nmom; momi++)
    {
      double theta = 0.0;
      for(int i=0;i<3;i++){theta += ( p0[i] * momLP[momi * 3 + i] * (pos[i] - offset[i]) );}
      Pres[momi][xi] = Ty(cos(theta), sign * sin(theta));
    }
  });
}

template <class Ty, int civ>
void apply_phases(qlat::FieldM<Ty, civ >& src, std::vector<qlat::FieldM<Ty, civ>>& res, std::vector<vector_gpu<Ty >>& phases){
  qassert(src.initialized);
  const Geometry& geo = src.geo();
  const unsigned int Nmom = phases.size();

  qlat::vector_acc<Ty* > Pphase = EigenM_to_pointers(phases);
  Ty* psrc = (Ty*) qlat::get_data(src).data();

  qlat::vector_acc<Ty* > pres;pres.resize(Nmom);
  if(res.size() != Nmom){res.resize(Nmom);}
  for(unsigned int mi=0;mi<Nmom;mi++){
    if(!res[mi].initialized){res[mi].init(geo);}
    pres[mi] = (Ty*) qlat::get_data(res[mi]).data();
  }

  qlat::vector_acc<int > nv, Nv, mv;geo_to_nv(geo, nv, Nv, mv);
  const long vol = Nv[0]*Nv[1]*Nv[2];
  const int   Nt = Nv[3];

  qacc_for(xi, long(vol),{
    for(int momi = 0;momi < Nmom; momi++)
    {
      const Ty& ph   = Pphase[momi][xi];
      for(int ti=0;ti<Nt;ti++)
      for(int ci=0;ci<civ;ci++){
        pres[momi][(ti*vol+xi)*civ + ci] = psrc[(ti*vol+xi)*civ + ci] * ph;
      }
    }
  });

}

template <class Ty, int civ>
void multi_factor(qlat::FieldM<Ty, civ >& src, qlat::FieldM<Ty, civ>& res, const Ty& phases){
  qassert(src.initialized);
  const Geometry& geo = src.geo();
  if(!res.initialized){res.init(geo);}

  const Ty* psrc = (Ty*) qlat::get_data(src).data();
  Ty* pres = (Ty*) qlat::get_data(res).data();

  qacc_for(xi, geo.local_volume(),{
    for(int ci=0;ci<civ;ci++){
      res[xi*civ + ci] = src[xi*civ + ci] * phases;
    }
  });
}

////phases for momentum data
template<typename Ty>
void get_phases(std::vector<Ty >& phases, Coordinate& pL, const Coordinate& src, const Coordinate& Lat){
  long vol = pL[0]*pL[1]*pL[2];
  phases.resize(vol);
  for(int i=0;i<3;i++){qassert(pL[i] % 2 == 1);}
  #pragma omp parallel for
  for(long isp =0;isp<vol;isp++){
    Coordinate pos = qlat::coordinate_from_index(isp, pL);

    for(int i=0;i<3;i++){
      if(pos[i] > pL[i]/2){
        pos[i] = Lat[i] - (pL[i] - pos[i]);
      }
    }
    double v0 = 0.0;
    for(int i=0;i<3;i++){v0 += (2.0*PI * src[i] * pos[i]/Lat[i]);}

    phases[isp] = Ty(std::cos(v0), -1.0* std::sin(v0));
  }
}

/////V -- 12a x 12b   to   12b x 12a -- V
template<typename Ty>
void copy_qprop_to_propE(std::vector<qlat::vector_acc<Ty > >& Eprop, std::vector<qpropT >& src, int dir = 1){
  TIMERA("copy_qprop_to_propE");
  const int nmass = src.size();
  std::vector<Ty* > ps;ps.resize(nmass);
  
  for(int mi=0;mi<nmass;mi++){
    qassert(src[mi].initialized);
    ps[mi] = (Ty*) qlat::get_data(src[mi]).data();
  }

  const qlat::Geometry &geo = src[0].geo();
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);
  if(dir == 1){ini_propE(Eprop, nmass, fd);}
  
  ///const long sizeF = geo.local_volume();
  const long nvec  = Eprop.size()/nmass;
  const long sizeF = Eprop[0].size();
  qassert(nvec * sizeF == 12*12*geo.local_volume());

  ////V x 12 a x 12 b to 12b x 12a x V
  for(int mi=0;mi<nmass;mi++)
  for(long i=0;i<nvec;i++)
  {
    if(dir == 1)cpy_data_thread(Eprop[mi*nvec + i].data(), &ps[mi][i*sizeF], sizeF, 1, false);
    if(dir == 0)cpy_data_thread(&ps[mi][i*sizeF], Eprop[mi*nvec + i].data(), sizeF, 1, false);
  }
  qacc_barrier(dummy);
}

template<typename Ty>
void copy_propE_to_qprop(std::vector<qpropT >& src, std::vector<qlat::vector_acc<Ty > >& Eprop){
  copy_qprop_to_propE(Eprop, src, 0);
}

template<typename Ty>
void noise_to_propT(qpropT& prop, qnoiT& noi){
  qassert(noi.initialized);
  const Geometry& geo = noi.geo();

  if(!prop.initialized){prop.init(geo);}

  Ty* res = (Ty*) qlat::get_data(prop).data();
  Ty* src = (Ty*) qlat::get_data(noi ).data();

  const long Nvol = geo.local_volume();

  for(int d0=0;d0<12;d0++){
    cpy_data_thread(&res[(d0*12+d0)*Nvol + 0], src, Nvol, 1, false);
  }
  qacc_barrier(dummy);

}

template <class Td>
void prop4D_factor(Propagator4dT<Td>& prop, const qlat::ComplexT<Td >& factor)
{
  const size_t Nvol = size_t(prop.geo().local_volume()) * 12 * 12;
  qlat::ComplexT<Td >* src = (qlat::ComplexT<Td >*) qlat::get_data(prop).data();
  qacc_for(isp, long(Nvol), {
    src[isp] *= factor;
  })
  ////if(qlat::qnorm(factor) >= QLAT_COPY_LIMIT){cpy_data_threadC(src, src, Nvol, 1, true, factor);}
  ////else{;}
}



template <class Ty, int civ>
void copy_FieldM(std::vector<qlat::FieldM<Ty, civ> >& res, const std::vector<qlat::FieldM<Ty, civ> >& src)
{
  if(src.size() == 0){res.resize(0); return;}
  const int Nv = src.size();
  for(int iv=0;iv<Nv;iv++){qassert(src[iv].initialized);}
  if(res.size() != Nv){res.resize(Nv);}

  const Geometry& geo = src[0].geo();
  const size_t Nvol = size_t(geo.local_volume()) * civ;
  for(int iv=0;iv<Nv;iv++)
  {
    res[iv].init(geo);
    Ty* resP = (Ty*) qlat::get_data(res[iv]).data();
    Ty* srcP = (Ty*) qlat::get_data(src[iv]).data();
    cpy_data_thread(resP, srcP, Nvol, 1, false);
  }
  qacc_barrier(dummy);
}


}

#endif
