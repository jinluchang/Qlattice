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

namespace qlat{

template<typename Ty>
void prop4d_src_gamma(Propagator4dT<Ty >& prop, ga_M& ga,int dir = 0){
  ////Rowmajor (a,b), b is continues in memory
  qacc_for(isp, long(prop.geo().local_volume()),{
    qlat::WilsonMatrixT<Ty>& v0 =  prop.get_elem_offset(isp);
    qlat::WilsonMatrixT<Ty>  v1 = v0;

    for (int s = 0; s < 4; ++s)
    for(int c0 = 0;c0< 3 ; c0++)
    for (int d0 = 0; d0 < 4; ++d0)
    {
      /////Source multiply
      if(dir==0)for(int c1=0;c1<3;c1++)v0(s*3 + c0, ga.ind[d0]*3 + c1) = ga.g[d0] * v1(s*3 + c0, d0*3 + c1);
      /////Sink multiply
      if(dir==1)for(int c1=0;c1<3;c1++)v0(d0*3 + c0, s*3 + c1) = ga.g[d0] * v1(ga.ind[d0]*3 + c0, s*3 + c1);

      ///////Source multiply
      //if(dir==0)for(int c1=0;c1<3;c1++)cp_C(v0(s*3 + c0, ga.ind[d0]*3 + c1), ga.g[d0] * v1(s*3 + c0, d0*3 + c1));
      ///////Sink multiply
      //if(dir==1)for(int c1=0;c1<3;c1++)cp_C(v0(d0*3 + c0, s*3 + c1), ga.g[d0] * v1(ga.ind[d0]*3 + c0, s*3 + c1));
    }
  });
}

template<typename Ty>
void prop4d_sink_gamma(Propagator4dT<Ty >& prop, ga_M& ga){
  prop4d_src_gamma(prop, ga ,1);
}


template<typename Ty>
void prop4d_cps_to_ps(Propagator4dT<Ty >& prop, int dir=0){
  /////sn is -1 for default
  Ty sn =-1;if(dir == 1){sn= 1;}
  const Ty sqrt2= Ty(std::sqrt(2.0), 0.0);

  ////Rowmajor (a,b), b is continues in memory
  qacc_for(isp, prop.geo().local_volume(),{
    qlat::WilsonMatrixT<Ty >  v0 = prop.get_elem_offset(isp);
    qlat::WilsonMatrixT<Ty >  v1 = prop.get_elem_offset(isp);

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

template<typename Ty>
void prop4d_ps_to_cps(Propagator4dT<Ty >& prop){
  prop4d_cps_to_ps(prop, 1);
}

template<typename Ty>
void get_corr_pion(std::vector<qlat::FermionField4dT<Ty > > &prop,const Coordinate &x_ini, std::vector<double > &write ){

  const qlat::Geometry &geo = prop[0].geo();

  unsigned long Nvol = geo.local_volume();
  ///int Nt = geo.node_site[3];
  ///long Nsum = Nvol/Nt;
  int tini = x_ini[3];

  qlat::vector_acc<Ty > res;res.resize(Nvol);

  qacc_for(isp, long(Nvol),{
    Ty buf(0.0,0.0);

    for(int dc2=0;dc2<12;dc2++){
      Ty* a = (Ty* ) &(prop[dc2].get_elem_offset(isp));
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
void vec_corrE(Ty* srcE, qlat::vector_acc<Ty >& res,qlat::fft_desc_basic &fd,const int nvec,const int clear=0,const Coordinate& mom = Coordinate(), const Ty& src_phase = 1.0){
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
    s1[mi*nt + t_rank + ti ] = s0[mi*NTt + ti] * src_phase;
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

template<typename Ta>
void shift_result_t(EigenVTa& Esrc, int nt, int tini){
  if(tini == 0){return ;}
  long Ntotal = Esrc.size();
  if(Ntotal %(nt) != 0){abort_r("Correlation function size wrong!\n");}
  EigenVTa tmp;tmp.resize(Ntotal);
  qacc_for(i, Ntotal, {
    const int iv = i/nt;
    const int t  = i%nt;
    tmp[iv*nt + (t - tini + nt)%nt] = Esrc[iv*nt + (t)%nt];
  });
  cpy_data_thread(Esrc.data(), tmp.data(), tmp.size(), 1);
}

template<typename Ta>
void ini_propE(EigenMTa &prop,int nmass, qlat::fft_desc_basic &fd, bool clear = true){
  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  int NTt  = fd.Nv[3];
  int do_resize = 0;
  if(prop.size() != (LInt) (nmass*12*12*NTt)){do_resize = 1;}
  for(unsigned int i=0;i<prop.size();i++){if((LInt) prop[i].size() != Nxyz){do_resize=1;}}

  if(do_resize == 1)
  {
    for(unsigned int i=0;i<prop.size();i++){prop[i].resize(0);}prop.resize(0);
    prop.resize(nmass*12*12*NTt);
    for(unsigned int i=0;i<prop.size();i++){
      prop[i].resize(Nxyz);
    }
  }
  if(clear){zeroE(prop);}
}

////default dir == 0 from prop4d to propE
template<typename Ty, typename Ta>
void copy_propE(std::vector<Ty* > &pV1, EigenMTa &prop, qlat::fft_desc_basic &fd, int dir=0){
  TIMER("Copy prop");
  qassert(fd.order_ch == 0);
  const LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  const int NTt  = fd.Nv[3];
  Geometry geo;fd.get_geo(geo);
  int nmass = 0;
  if(dir==0){nmass = pV1.size();ini_propE(prop,nmass,fd);}
  if(dir==1){
    nmass = prop.size()/(12*12*NTt);
    qassert(int(pV1.size()) == nmass);
    //if(pV1.size() != (LInt) nmass){
    //  pV1.resize(0);
    //  pV1.resize(nmass);for(int i=0;i<nmass;i++){pV1[i].init(geo);}
    //}
  }

  for(int mi = 0;mi < nmass;mi++)
  {
    Ty* pv = pV1[mi];
    /////Propagator4dT<Ty >& pv = pV1[mi];
    qlat::vector_acc<Ta* > ps = EigenM_to_pointers(prop);
    qacc_for(isp, long(NTt*Nxyz),{
      int ti = isp/Nxyz;
      int xi = isp%Nxyz;
      /////qlat::WilsonMatrixT<Ty>& v0 =  pv.get_elem_offset(isp);
      Ty* v0 = &pv[isp * 12 * 12];

      for(int c0 = 0;c0 < 3; c0++)
      for(int d0 = 0;d0 < 4; d0++)
      for(int c1 = 0;c1 < 3; c1++)
      for(int d1 = 0;d1 < 4; d1++)
      {
        LInt off = mi*12*12 + (d1*3+c1)*12+d0*3+c0;
        if(dir==0){ps[off*NTt+ti][xi] = v0[(d0*3 + c0)*12 +  d1*3 + c1];}
        if(dir==1){v0[(d0*3 + c0)*12 +  d1*3 + c1] = ps[off*NTt+ti][xi];}
      }
    });
  }

}

template<typename Ty, typename Ta>
void copy_propE(Propagator4dT<Ty > &pV, EigenMTa &prop, qlat::fft_desc_basic &fd, int dir=0){
  int nmass = 1;
  if(dir==0){
    qassert(pV.initialized);
    ini_propE(prop,nmass,fd);
  }
  if(dir==1){
    nmass = prop.size()/(12*12*fd.Nv[3]);
    qassert(nmass == 1);
    if(!pV.initialized){
      Geometry geo;fd.get_geo(geo);
      pV.init(geo);
    }
  }
  std::vector<Ty* > PTy;PTy.resize(1);
  PTy[0] = (Ty*) (qlat::get_data(pV).data());
  copy_propE(PTy, prop, fd, dir);
}

template<typename Ty, typename Ta>
void copy_propE(std::vector<Propagator4dT<Ty > > &pV, EigenMTa &prop, qlat::fft_desc_basic &fd, int dir=0){
  int nmass = 0;
  if(dir==0){nmass = pV.size();ini_propE(prop,nmass,fd);}
  if(dir==1){
    nmass = prop.size()/(12*12*fd.Nv[3]);
    bool do_init = false;
    if(pV.size() != (LInt) nmass){do_init = true;}
    if(pV.size() > 0){if(!pV[0].initialized){do_init = true;}}
    if(do_init)
    {
      pV.resize(0);
      Geometry geo;fd.get_geo(geo);
      pV.resize(nmass);for(int i=0;i<nmass;i++){pV[i].init(geo);}
    }
  }
  std::vector<Ty* > PTy;PTy.resize(nmass);
  for(int im=0;im<nmass;im++){PTy[im] = (Ty*) (qlat::get_data(pV[im]).data());}
  copy_propE(PTy, prop, fd, dir);
}

template<typename Ty, typename Ta>
void copy_prop4d_to_propE(EigenMTa &prop, std::vector<Propagator4dT<Ty > > &pV1, qlat::fft_desc_basic &fd){
  copy_propE(pV1, prop, fd, 0);
}
template<typename Ty, typename Ta>
void copy_propE_to_prop4d(std::vector<Propagator4dT<Ty > > &pV1, EigenMTa &prop, qlat::fft_desc_basic &fd){
  copy_propE(pV1, prop, fd, 1);
}

template<typename Ty, typename Ta>
void copy_prop4d_to_propE(EigenMTa &prop, Propagator4dT<Ty > &pV1, qlat::fft_desc_basic &fd){
  copy_propE(pV1, prop, fd, 0);
}
template<typename Ty, typename Ta>
void copy_propE_to_prop4d(Propagator4dT<Ty > &pV1, EigenMTa &prop, qlat::fft_desc_basic &fd){
  copy_propE(pV1, prop, fd, 1);
}

template <typename Ta >
void check_prop_size(EigenMTa &prop){
  int sizep = prop.size();
  if(sizep%(12*12) != 0 or sizep == 0)
  {
    print0("Size of Prop wrong. \n");
    qassert(false);
  }
}

template <typename Ta >
void ini_resE(EigenVTa &res, int nmass, qlat::fft_desc_basic &fd){
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

inline std::vector<int >  get_sec_map(int dT,int nt){
  std::vector<int > map_sec;map_sec.resize(nt);
  int secN = 2*nt/dT;double lensec = nt/(1.0*secN);
  int tcount = 0;
  int t0 = 0;
  for(int si=0;si<secN;si++)
  {
    for(int t=t0;t <= (si+1)*lensec;t++)
    {
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


}

#endif
