#ifndef UTILS_CONSTRUCTION_H
#define UTILS_CONSTRUCTION_H

#pragma once
#include <qlat/qlat.h>
#include "gammas.h"
#include "utils_momentum.h"
#include "fft_desc.h"

namespace qlat{

//template<typename Ty1, typename Ty2>
//inline void cp_C(Ty1 a, Ty2 b)
//{
//  if(typeid(Ty1) == typeid(Ty2)){a = b;}
//  else{
//    a = Ty1(b.real, b.imag);
//  }
//}

void prop4d_src_gamma(Propagator4d& prop, ga_M& ga,int dir = 0)
{
  ////Rowmajor (a,b), b is continues in memory
  qacc_for(isp, long(prop.geo().local_volume()),{ 
    qlat::WilsonMatrix& v0 =  prop.get_elem(isp);
    qlat::WilsonMatrix  v1 = v0;

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

void prop4d_sink_gamma(Propagator4d& prop, ga_M& ga)
{
  prop4d_src_gamma(prop, ga ,1);
}

void prop4d_cps_to_ps(Propagator4d& prop, int dir=0)
{
  /////sn is -1 for default
  double sn =-1;if(dir == 1){sn= 1;}
  const double sqrt2=std::sqrt(2.0);

  ////Rowmajor (a,b), b is continues in memory
  qacc_for(isp, prop.geo().local_volume(),{ 
    qlat::WilsonMatrix  v0 = prop.get_elem(isp);
    qlat::WilsonMatrix  v1 = prop.get_elem(isp);

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
    prop.get_elem(isp) = v0;

  });
}

void prop4d_ps_to_cps(Propagator4d& prop)
{
  prop4d_cps_to_ps(prop, 1);
}

void get_corr_pion(std::vector<qlat::FermionField4dT<qlat::Complex> > &prop,const Coordinate &x_ini, std::vector<double > &write ){

  const qlat::Geometry &geo = prop[0].geo();

  unsigned long Nvol = geo.local_volume();
  ///int Nt = geo.node_site[3];
  ///long Nsum = Nvol/Nt;
  int tini = x_ini[3];

  std::vector<qlat::Complex > res;res.resize(Nvol);

  qacc_for(isp, long(Nvol),{
    qlat::Complex buf(0.0,0.0);

    for(int dc2=0;dc2<12;dc2++){
      qlat::Complex* a = (qlat::Complex* ) &(prop[dc2].get_elem(isp));
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
  sum_all_size((double*) &write[0],2*nt);

}

void vec_corrE(EigenM &resE,Evector &res,qlat::fft_desc_basic &fd,int clear=0,int imom=505050)
{
  int NTt  = fd.Nv[3];
  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  int nmass = resE.size()/(NTt);

  int nt = fd.nt;

  if(clear == 1){res.resize(0);res.resize(nmass*nt);qlat::set_zero(res);}
  if(clear == 0)if(res.size() != nmass*nt){print0("res size wrong for corr.");qassert(false);}

  Evector tmp;tmp.resize(res.size());qlat::set_zero(tmp);

  ////position p = fd.desc.get_position(0,fd.rank);
  int t_rank = fd.Pos0[fd.rank][3];

  if(imom != 505050)
  {
    int mom[3];
    qlat::momentum tem_mom(imom);
    tem_mom.ind_2_mom(mom);
    double p0[3]={2*3.1415926535898/fd.nx,
                  2*3.1415926535898/fd.ny,
                  2*3.1415926535898/fd.nz};

    Evector phaseE;phaseE.resize(Nxyz);

    qacc_for(xi, long(Nxyz),{
      int pi[3];
      pi[fd.orderN[0]] = xi/(fd.Nv[fd.orderN[1]]*fd.Nv[fd.orderN[2]]);
      pi[fd.orderN[1]] = (xi%(fd.Nv[fd.orderN[1]]*fd.Nv[fd.orderN[2]]))/fd.Nv[fd.orderN[2]];
      pi[fd.orderN[2]] = xi%fd.Nv[fd.orderN[2]];
      for(int ptem=0;ptem<3;ptem++){pi[ptem] = pi[ptem] + fd.Pos0[fd.rank][ptem];}

      double theta=mom[0]*p0[0]*pi[0]+mom[1]*p0[1]*pi[1]+mom[2]*p0[2]*pi[2];
      phaseE[xi] = Complexq(cos(theta),sin(theta));
    });
    //////print0("Phase %13.5f \n",phaseE[0].real());

    //#ifndef QLAT_USE_ACC
    //for(int mi=0;mi<nmass;mi++)
    //for(int ti=0;ti<NTt;ti++)
    //{
    //  resE[mi*NTt + ti] = resE[mi*NTt + ti]*phaseE;
    //}
    //#endif

    ////#ifdef QLAT_USE_ACC
    size_t Ns = Nxyz;
    Ns = nmass*NTt*Nxyz;
    qacc_for(i, long(Ns),{
      LInt mi = i/(NTt*Nxyz);
      LInt ti = (i%(NTt*Nxyz))/Nxyz;
      LInt xi = i%(Nxyz);

      resE[mi*NTt + ti][xi] = resE[mi*NTt + ti][xi]*phaseE[xi];
    });
    /////#endif
  }

  for(int mi=0;mi<nmass;mi++)
  for(int ti=0;ti<NTt;ti++)
  {
    //tmp[mi*nt + t_rank + ti ] += resE[mi*NTt + ti].sum();
    tmp[mi*nt + t_rank + ti ] = 0.0;
    for(unsigned int xi=0;xi<Nxyz;xi++){
      tmp[mi*nt + t_rank + ti ] += resE[mi*NTt + ti][xi];
    }
  }

  sum_all_size(reinterpret_cast<Ftype* > (&tmp[0]),2*tmp.size());

  for(int mi=0;mi<nmass;mi++)
  for(int t=0;t<fd.nt;t++)
  {
    res[mi*nt + t ] += tmp[mi*nt + t ];
  }

}

void ini_propE(EigenM &prop,int nmass, qlat::fft_desc_basic &fd)
{
  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  int NTt  = fd.Nv[3];
  prop.resize(0);
  prop.resize(nmass*12*12*NTt);
  for(unsigned int i=0;i<prop.size();i++){
    prop[i].resize(Nxyz);
    //Evector& p = prop[i];
    //qacc_for(isp, long(Nxyz),{ p[isp] = 0.0;});
    //set_zero(prop[i]);
    zeroE(prop[i],0);
  }
}

void copy_propE(std::vector<Propagator4d > &pV1,EigenM &prop, qlat::fft_desc_basic &fd, int dir=0)
{
  TIMER("Copy prop");
  qassert(fd.order_ch == 0);
  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  int NTt  = fd.Nv[3];
  int nmass = 0;
  if(dir==0){nmass = pV1.size();ini_propE(prop,nmass,fd);}
  if(dir==1){
    nmass = prop.size()/(12*12*NTt);
    pV1.resize(0);
    pV1.resize(nmass);for(int i=0;i<nmass;i++){pV1[i].init(*fd.geop);}
  }

  ////////only cpu version
  //#pragma omp parallel for
  //for(long isp=0;isp< long(NTt*Nxyz);isp++)
  //{
  //  int ti = isp/Nxyz;
  //  int xi = isp%Nxyz;

  //  for(int mi = 0;mi < nmass;mi++)
  //  {
  //    qlat::WilsonMatrix& v0 =  pV1[mi].get_elem(isp);

  //    for(int c0 = 0;c0 < 3; c0++)
  //    for(int d0 = 0;d0 < 4; d0++)
  //    for(int c1 = 0;c1 < 3; c1++)
  //    for(int d1 = 0;d1 < 4; d1++)
  //    {
  //      LInt off = mi*12*12 + (d1*3+c1)*12+d0*3+c0;
  //      if(dir==0)prop[off*NTt+ti][xi] = v0(d0*3 + c0, d1*3 + c1);
  //      if(dir==1)v0(d0*3 + c0, d1*3 + c1) = prop[off*NTt+ti][xi];
  //    }
  //  }
  //}

  for(int mi = 0;mi < nmass;mi++)
  {
    Propagator4d& pv = pV1[mi];
    qacc_for(isp, long(NTt*Nxyz),{ 
      int ti = isp/Nxyz;
      int xi = isp%Nxyz;

        qlat::WilsonMatrix& v0 =  pv.get_elem(isp);

        for(int c0 = 0;c0 < 3; c0++)
        for(int d0 = 0;d0 < 4; d0++)
        for(int c1 = 0;c1 < 3; c1++)
        for(int d1 = 0;d1 < 4; d1++)
        {
          LInt off = mi*12*12 + (d1*3+c1)*12+d0*3+c0;
          if(dir==0)prop[off*NTt+ti][xi] = v0(d0*3 + c0, d1*3 + c1);
          if(dir==1)v0(d0*3 + c0, d1*3 + c1) = prop[off*NTt+ti][xi];
        }
    });
  }

}


void check_prop_size(EigenM &prop)
{
  int sizep = prop.size();
  if(sizep%(12*12) != 0 or sizep == 0)
  {
    print0("Size of Prop wrong. \n");
    qassert(false);
    ///shutdown_machine();
    ///abort();
  }
}

void ini_resE(EigenM &res,int nmass, qlat::fft_desc_basic &fd)
{
  ////int Nvol = fd.Nv[fd.orderN[0]]*fd.Nv[fd.orderN[1]]*fd.Nv[fd.orderN[2]];
  int NTt  = fd.Nv[3];
  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  res.resize(0);
  res.resize(nmass*NTt);
  for(unsigned int i=0;i<res.size();i++){
    res[i].resize(Nxyz);
    qacc_for(isp, long(Nxyz),{ res[i][isp] = 0.0;});
  }

}

void clear_vector(EigenM &res,qlat::fft_desc_basic &fd)
{
  ////int NTt  = fd.Nv[3];
  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  for(LInt i=0;i<res.size();i++){
    res[i].resize(Nxyz);
    qacc_for(isp, long(Nxyz),{ res[i][isp] = 0.0;});
    /////res[i].setZero();
  }
}

void clear_qv(qlat::vector<Complexq > &G)
{
  qacc_for(i, long(G.size()),{G[i] = 0.0; });
  /////for(int i = 0;i< G.size(); i++){G[i] = 0.0;}
}


void meson_vectorE(std::vector<Propagator4d > &pV1, std::vector<Propagator4d > &pV2, ga_M &ga1,ga_M &ga2,
        EigenM &res, qlat::fft_desc_basic &fd,int clear=1)
{
  TIMER("Meson_vectorE");
  qassert(fd.order_ch == 0);
  ///////check_prop_size(prop1);check_prop_size(prop2);
  int  NTt  = fd.Nv[3];
  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  int  nmass = pV1.size();
  if(nmass == 0){res.resize(0);return;}

  if(clear == 1){
    ini_resE(res,nmass,fd);
    ////res.resize(nmass*NTt);
    ////clear_vector(res, fd);
  }
  if(res.size()%NTt !=0 or res.size()==0){print0("Size of res wrong. \n");qassert(false);}
  qassert(pV1.size() == pV2.size());

  for(int mi=0;mi<nmass;mi++)
  {
  Propagator4d& pL1 = pV1[mi];
  Propagator4d& pL2 = pV2[mi];

  qacc_for(isp, long(pV1[0].geo().local_volume()),{ 
    int ti = isp/Nxyz;
    int xi = isp%Nxyz;
      //std::complex< Ftype > pres;pres = 0.0;
      Complexq pres;pres = 0.0;
      const qlat::WilsonMatrix& p1 =  pL1.get_elem(isp);
      const qlat::WilsonMatrix& p2 =  pL2.get_elem(isp);
      ///int off1 = mi*12*12 + (d2*3+c2)*12+ga1.ind[d1]*3+c1;
      ///int off2 = mi*12*12 + (ga2.ind[d2]*3+c2)*12+d1*3+c1;
      ///double_complex gi_tem = ga2.g[d2]*ga1.g[d1];
      ///std::complex<Ftype> giE(gi_tem.real,gi_tem.imag);
      ///res[mi*NTt + ti] += (prop1[off1*NTt+ti]*prop2[off2*NTt+ti].conjugate())*(giE);

      for(int d1=0;d1<4;d1++)
      for(int c1=0;c1<3;c1++)
      for(int d2=0;d2<4;d2++)
      {
      Complexq g_tem = ga2.g[d2]*ga1.g[d1];
      for(int c2=0;c2<3;c2++)
      {
        pres += g_tem * 
          p1(ga1.ind[d1]*3+c1,d2*3+c2) * qlat::qconj(p2(d1*3+c1,ga2.ind[d2]*3+c2)) ;

        //pres += 
        //  p1(ga1.ind[d1]*3+c1,d2*3+c2) * qlat::qconj(p2(d1*3+c1,ga2.ind[d2]*3+c2)) 
        //      * g_tem;
      }
      }
      res[mi*NTt + ti][xi%Nxyz] += pres;
  });
  }

}

void meson_vectorE(EigenM &prop1, EigenM &prop2, ga_M &ga1,ga_M &ga2,
        EigenM &res, qlat::fft_desc_basic &fd,int clear=1, int invmode=1)
{
  TIMER("Meson_vectorE");
  check_prop_size(prop1);check_prop_size(prop2);
  ///////check_prop_size(prop1);check_prop_size(prop2);
  int  NTt  = fd.Nv[3];
  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  int  nmass = prop1.size()/(12*12*NTt);
  if(nmass == 0){res.resize(0);return;}

  if(clear == 1){
    ini_resE(res,nmass,fd);
    ////res.resize(nmass*NTt);
    ////clear_vector(res, fd);
  }
  if(res.size()%NTt !=0 or res.size()==0){print0("Size of res wrong. \n");qassert(false);}
  qassert(prop1.size() == prop2.size());

  for(int mi=0;mi<nmass;mi++)
  for(int d2=0;d2<4;d2++)
  for(int c2=0;c2<3;c2++)
  for(int d1=0;d1<4;d1++)
  for(int c1=0;c1<3;c1++)
  {
  int off1 = mi*12*12 + (d2*3+c2)*12+ga1.ind[d1]*3+c1;
  int off2 = mi*12*12 + (ga2.ind[d2]*3+c2)*12+d1*3+c1;

  for(int ti=0;ti<NTt;ti++)
  {
    //double_complex gi_tem = ga2.g[d2]*ga1.g[d1];
    //std::complex<Ftype> giE(gi_tem.real,gi_tem.imag);
    //res[mi*NTt + ti] += (prop1[off1*NTt+ti]*prop2[off2*NTt+ti].conjugate())*(giE);
    EA vp1(&prop1[off1*NTt+ti][0],Nxyz);
    EA vp2(&prop2[off2*NTt+ti][0],Nxyz);
    EA vr0(&res[mi*NTt + ti][0],Nxyz);
    Complexq g_tem = ga2.g[d2]*ga1.g[d1];
    if(invmode == 1)vr0 += vp1*vp2.conjugate() * g_tem;
    if(invmode == 0)vr0 += vp1*vp2             * g_tem;
  }
  }

}


void meson_corrE(std::vector<Propagator4d > &pV1, std::vector<Propagator4d > &pV2,  ga_M &ga1, ga_M &ga2,
  Evector &res, qlat::fft_desc_basic &fd,int clear=1,int imom=505050,int mode_GPU = 0)
{
  ///int NTt  = fd.Nv[3];
  ///LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  int nmass = pV1.size();
  ///int nt = fd.nt;

  EigenM resE;
  ini_resE(resE,nmass,fd);

  if(mode_GPU == 0){
    EigenM prop1;
    EigenM prop2;
    copy_propE(pV1,prop1, fd);
    copy_propE(pV2,prop2, fd);
    meson_vectorE(prop1,prop2,ga1,ga2,resE,fd,1);
  }
  if(mode_GPU == 1){meson_vectorE(pV1,pV2,ga1,ga2,resE,fd,1);}

  vec_corrE(resE,res,fd, clear, imom);
}

void meson_corrE(EigenM &prop1,EigenM &prop2,  ga_M &ga1, ga_M &ga2,
  Evector &res, qlat::fft_desc_basic &fd,int clear=1,int imom=505050)
{
  ///int NTt  = fd.Nv[3];
  ////LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  ///int  nmass = prop1.size()/(12*12*NTt);
  ////int nt = fd.nt;

  EigenM resE;
  ////ini_resE(resE,nmass,fd);

  meson_vectorE(prop1,prop2,ga1,ga2,resE,fd,1);
  vec_corrE(resE,res,fd, clear, imom);
}


///////Proton contractions

void proton_vectorE(EigenM &prop1, EigenM &prop2, EigenM &prop3,
  ga_M &ga2,int ind2, ga_M &ga1,int ind1,
        EigenM &res, fft_desc_basic &fd,int clear=1)
{
  TIMER("Proton_vectorE");
  ////ga2/ind2 for source, gam1/ind1 for sink
  ////"[]|+N" type diagram
  check_prop_size(prop1);check_prop_size(prop2);check_prop_size(prop3);
  int NTt  = fd.Nv[3];
  LInt Nxyz = prop1[0].size();
  int nmass = prop1.size()/(12*12*NTt);
  qassert(prop1.size() == prop2.size());
  qassert(prop1.size() == prop3.size());
  if(clear == 1){ini_resE(res,nmass,fd);}
    
  ////Prop format, src d-4, c-3, sink d-4, c-3, Nt, EigenV<Nxyz>
  ////EA vs0(&outv_sinki0[0],maxN);
  ///int Nt   = prop1.size()/(12*12);
  if(res.size()%NTt !=0 or res.size()==0){print0("Size of res wrong. \n");qassert(false);}

  Ftype epsl[3][3];
  for(int i=0;i<3;i++){epsl[i][i]=0;epsl[i][(i+1)%3]=1;epsl[i][(i+2)%3]=-1;}

  for(int mi=0;mi<nmass;mi++)
  for(int d2=0;d2<4;d2++)
  for(int c21=0;c21<3;c21++)
  for(int ib=1;ib<3;ib++)
  {
    int c22=(c21+ib)%3,c23=(c22+ib)%3;
    for(int d1=0;d1<4;d1++)
    for(int c11=0;c11<3;c11++)
    for(int ia=1;ia<3;ia++)
    {
      int c12=(c11+ia)%3,c13=(c12+ia)%3;

      int m1 = mi*12*12 + (ind2*3+c21)*12+ind1*3+c11;
      int m2 = mi*12*12 + (ga2.ind[d2]*3+c22)*12+d1*3+c12;
      int m3 = mi*12*12 + (d2*3+c23)*12+ga1.ind[d1]*3+c13;

      int n1 = mi*12*12 + (ind2*3+c21)*12+ga1.ind[d1]*3+c11;
      int n2 = mi*12*12 + (ga2.ind[d2]*3+c22)*12+d1*3+c12;
      int n3 = mi*12*12 + (d2*3+c23)*12+ind1*3+c13;

      //double_complex gi_tem = epsl[c11][c12]*epsl[c21][c22]*ga1.g[d1]*ga2.g[d2];
      //std::complex<Ftype> giE(gi_tem.real,gi_tem.imag);
      Complexq giE = epsl[c11][c12]*epsl[c21][c22]*ga1.g[d1]*ga2.g[d2];

      for(int ti=0;ti<NTt;ti++)
      {
        EA vp1(&prop1[m1*NTt+ti][0],Nxyz);
        EA vp2(&prop2[m2*NTt+ti][0],Nxyz);
        EA vp3(&prop3[m3*NTt+ti][0],Nxyz);
        EA vn1(&prop1[n1*NTt+ti][0],Nxyz);
        EA vn2(&prop2[n2*NTt+ti][0],Nxyz);
        EA vn3(&prop3[n3*NTt+ti][0],Nxyz);
        EA vr0(&res[mi*NTt + ti][0],Nxyz);
        vr0 -= (vp1*vp2*vp3 + vn1*vn2*vn3)*giE;
      }
    }
  }

}

std::vector<int >  get_sec_map(int dT,int nt)
{
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

void proton_vectorE(EigenM &prop1, EigenM &prop2, EigenM &prop3,
        EigenM &res, fft_desc_basic &fd, ga_M &ga1,int t0,int dT,int clear=1,int oppo=0)
{
  /////if(clear == 1)clear_vector(res);
  TIMER("Proton_vectorE");
  int NTt  = fd.Nv[3];
  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  int nmass = prop1.size()/(12*12*NTt);
  qassert(prop1.size() == prop2.size());
  qassert(prop1.size() == prop3.size());

  if(clear == 1){ini_resE(res,nmass,fd);}

  int nv = res.size();int Nsize = res[0].size();
  EigenM resE0;resE0.resize(nv);
  EigenM resE1;resE1.resize(nv);
  for(int i=0;i<nv;i++)
  {
    resE0[i].resize(Nsize);
    resE1[i].resize(Nsize);
  }

  //proton_vectorE(prop1,prop2,prop3,ga.ga[1][3],0,ga.ga[1][3],0,resE0,fd,1);
  //proton_vectorE(prop1,prop2,prop3,ga.ga[1][3],1,ga.ga[1][3],1,resE0,fd,0);
  //proton_vectorE(prop1,prop2,prop3,ga.ga[1][3],2,ga.ga[1][3],2,resE1,fd,1);
  //proton_vectorE(prop1,prop2,prop3,ga.ga[1][3],3,ga.ga[1][3],3,resE1,fd,0);

  proton_vectorE(prop1,prop2,prop3,ga1,0,ga1,0,resE0,fd,1);
  proton_vectorE(prop1,prop2,prop3,ga1,1,ga1,1,resE0,fd,0);
  proton_vectorE(prop1,prop2,prop3,ga1,2,ga1,2,resE1,fd,1);
  proton_vectorE(prop1,prop2,prop3,ga1,3,ga1,3,resE1,fd,0);

  std::vector<int > map_sec = get_sec_map(dT,fd.nt);
  ////EA vp1(&prop1[m1*NTt+ti][0],Nxyz);
  ////EA vp2(&prop2[m2*NTt+ti][0],Nxyz);

  int Nt = fd.Nt;

  /////int t0 = 0;
  int t = 0;int nt = fd.nt;
  for(int mi=0;mi<nmass;mi++)
  for(int ti = 0;ti<Nt;ti++)
  {
    t = ti + fd.Pos0[fd.rank][3];
    EA r0(&res[mi*Nt+ti][0],Nxyz);
    EA v0(&resE0[mi*Nt+ti][0],Nxyz);
    EA v1(&resE1[mi*Nt+ti][0],Nxyz);
    if(map_sec[(t-t0+nt)%nt]%2==0)
    {
      if(oppo==0)r0 += v0;
      if(oppo==1)r0 += v1;
      //if(oppo==0)res[mi*Nt+ti] += resE0[mi*Nt+ti];
      //if(oppo==1)res[mi*Nt+ti] += resE1[mi*Nt+ti];
    }
    if(map_sec[(t-t0+nt)%nt]%2==1)
    {
      if(oppo==0)r0 += v1;
      if(oppo==1)r0 += v0;
      //if(oppo==0)res[mi*Nt+ti] += resE1[mi*Nt+ti];
      //if(oppo==1)res[mi*Nt+ti] += resE0[mi*Nt+ti];
    }
  }

}

void proton_corrE(EigenM &prop1, EigenM &prop2, EigenM &prop3,
   ga_M &ga2,int ind2,ga_M &ga1,int ind1,
  Evector &res, fft_desc_basic &fd,int clear=1,int imom=50505)
{
  ///int NTt  = fd.Nv[3];
  ////LInt Nxyz = prop1[0].size();
  ///int nmass = prop1.size()/(12*12*NTt);
  ////int nt = fd.nt;

  EigenM resE;
  ////ini_resE(resE,nmass,fd);

  proton_vectorE(prop1,prop2,prop3,ga2,ind2,ga1,ind1,resE,fd,1);

  vec_corrE(resE,res,fd,clear,imom);
}

void proton_corrE(EigenM &prop1, EigenM &prop2, EigenM &prop3,
 Evector &res, fft_desc_basic &fd, ga_M &ga1,int t0,int dT,int clear=1,int imom=505050)
{
  ///int NTt  = fd.Nv[3];
  ////LInt Nxyz = prop1[0].size();
  ///int nmass = prop1.size()/(12*12*NTt);
  ////int nt = fd.nt;

  EigenM resE;
  ////ini_resE(resE,nmass,fd);
  proton_vectorE(prop1,prop2,prop3,resE,fd, ga1, t0,dT,1);

  vec_corrE(resE,res,fd,clear,imom);
}


/////A source gamma, B sink Gamma, G projections with fermion sign, mL shape of diagram
void baryon_vectorE(EigenM &prop1, EigenM &prop2, EigenM &prop3,
  ga_M &A, ga_M &B, qlat::vector<Complexq > &G, qlat::vector<int > &mL,
        EigenM &res, fft_desc_basic &fd,int clear=1)
{
  TIMER("Proton_vectorE");
  check_prop_size(prop1);check_prop_size(prop2);check_prop_size(prop3);
  int NTt  = fd.Nv[3];
  LInt Nxyz = prop1[0].size();
  int nmass = prop1.size()/(12*12*NTt);
  qassert(prop1.size() == prop2.size());
  qassert(prop1.size() == prop3.size());
  qassert(G.size()  == 16);
  qassert(mL.size() == 3);
  if(clear == 1){ini_resE(res,nmass,fd);}

  if(res.size()%NTt !=0 or res.size()==0){print0("Size of res wrong. \n");qassert(false);}

  Ftype epsl[3][3];
  for(int i=0;i<3;i++){epsl[i][i]=0;epsl[i][(i+1)%3]=1;epsl[i][(i+2)%3]=-1;}

  //mL = {};
  //std::vector<int > mL;mL.resize(3);
  //mL[0] = 0;mL[1] = 1;mL[2] = 2;
  std::vector<int > nmL;nmL.resize(3);
  std::vector<int > bmL;bmL.resize(3);

  for(int massi=0;massi<nmass;massi++)
  {
    for(int a1=0;a1<3;a1++)
    for(int ia=1;ia<3;ia++)
    for(int b1=0;b1<3;b1++)
    for(int ib=1;ib<3;ib++)
    {
      int b2=(b1+ib)%3,b3=(b2+ib)%3;
      int a2=(a1+ia)%3,a3=(a2+ia)%3;
      for(int m2=0;m2<4;m2++)
      for(int m1=0;m1<4;m1++)
      for(int n2=0;n2<4;n2++)
      for(int n1=0;n1<4;n1++)
      {
        Complexq Gv =  G[m1*4+n1];
        double norm = std::sqrt(Gv.real()*Gv.real() + Gv.imag()*Gv.imag());
        if(norm < 1e-20)continue;

        int m3 = A.ind[m2];
        int n3 = B.ind[n2];
        Complexq giE = epsl[a1][a2]*epsl[b1][b2]*A.g[m2]*B.g[n2]*G[m1*4+n1];
        nmL[0] = n1;nmL[1] = n2;nmL[2] = n3;
        bmL[0] = b1;bmL[1] = b2;bmL[2] = b3;
        int nm1 = nmL[mL[0]];
        int nm2 = nmL[mL[1]];
        int nm3 = nmL[mL[2]];

        int bm1 = bmL[mL[0]];
        int bm2 = bmL[mL[1]];
        int bm3 = bmL[mL[2]];

        int o1 = massi*12*12 + (m1*3+a1)*12+(nm1*3+bm1);
        int o2 = massi*12*12 + (m2*3+a2)*12+(nm2*3+bm2);
        int o3 = massi*12*12 + (m3*3+a3)*12+(nm3*3+bm3);

        for(int ti=0;ti<NTt;ti++)
        {
          EA vp1(&prop1[o1*NTt+ti][0],Nxyz);
          EA vp2(&prop2[o2*NTt+ti][0],Nxyz);
          EA vp3(&prop3[o3*NTt+ti][0],Nxyz);
          EA vr0(&res[massi*NTt + ti][0],Nxyz);
          vr0 += vp1*vp2*vp3 * giE;
        }
      }
    }
  }

}


////A source gamma, B sink Gamma, G projections with fermion sign, mL shape of diagram
void baryon_vectorE(EigenM &prop1, EigenM &prop2, EigenM &prop3,
  ga_M &A, ga_M &B, qlat::vector<Complexq > &G, qlat::vector<int > &mL, int insertion,
        EigenM &resP, fft_desc_basic &fd,int clear=1)
{
  TIMER("Proton_vectorE");
  check_prop_size(prop1);check_prop_size(prop2);check_prop_size(prop3);
  int NTt  = fd.Nv[3];
  LInt Nxyz = prop1[0].size();
  int nmass = prop1.size()/(12*12*NTt);
  qassert(prop1.size() == prop2.size());
  qassert(prop1.size() == prop3.size());
  qassert(G.size()  == 16);
  qassert(mL.size() == 3);
  qassert(fd.order_ch == 0);

  if(clear==1){ini_propE(resP,nmass,fd);}
  check_prop_size(resP);

  Ftype epsl[3][3];
  for(int i=0;i<3;i++){epsl[i][i]=0;epsl[i][(i+1)%3]=1;epsl[i][(i+2)%3]=-1;}

  //mL = {};
  //std::vector<int > mL;mL.resize(3);
  //mL[0] = 0;mL[1] = 1;mL[2] = 2;
  std::vector<int > nmL;nmL.resize(3);
  std::vector<int > bmL;bmL.resize(3);

  for(int massi=0;massi<nmass;massi++)
  {
    for(int a1=0;a1<3;a1++)
    for(int ia=1;ia<3;ia++)
    for(int b1=0;b1<3;b1++)
    for(int ib=1;ib<3;ib++)
    {
      int b2=(b1+ib)%3,b3=(b2+ib)%3;
      int a2=(a1+ia)%3,a3=(a2+ia)%3;
      for(int m2=0;m2<4;m2++)
      for(int m1=0;m1<4;m1++)
      for(int n2=0;n2<4;n2++)
      for(int n1=0;n1<4;n1++)
      {
        Complexq Gv =  G[m1*4+n1];
        double norm = std::sqrt(Gv.real()*Gv.real() + Gv.imag()*Gv.imag());
        if(norm < 1e-20)continue;

        int m3 = A.ind[m2];
        int n3 = B.ind[n2];
        Complexq giE = epsl[a1][a2]*epsl[b1][b2]*A.g[m2]*B.g[n2]*G[m1*4+n1];
        nmL[0] = n1;nmL[1] = n2;nmL[2] = n3;
        bmL[0] = b1;bmL[1] = b2;bmL[2] = b3;
        int nm1 = nmL[mL[0]];
        int nm2 = nmL[mL[1]];
        int nm3 = nmL[mL[2]];

        int bm1 = bmL[mL[0]];
        int bm2 = bmL[mL[1]];
        int bm3 = bmL[mL[2]];

        int o1 = massi*12*12 + (m1*3+a1)*12+(nm1*3+bm1);
        int o2 = massi*12*12 + (m2*3+a2)*12+(nm2*3+bm2);
        int o3 = massi*12*12 + (m3*3+a3)*12+(nm3*3+bm3);

        int r0 = 0;
        if(insertion == 0)r0 = massi*12*12 + (m1*3+a1)*12+(nm1*3+bm1);
        if(insertion == 1)r0 = massi*12*12 + (m2*3+a2)*12+(nm2*3+bm2);
        if(insertion == 2)r0 = massi*12*12 + (m3*3+a3)*12+(nm3*3+bm3);

        for(int ti=0;ti<NTt;ti++)
        {
          EA vp1(&prop1[o1*NTt+ti][0],Nxyz);
          EA vp2(&prop2[o2*NTt+ti][0],Nxyz);
          EA vp3(&prop3[o3*NTt+ti][0],Nxyz);
          EA vr0( &resP[r0*NTt+ti][0],Nxyz);
          if(insertion == 0)vr0 += vp2*vp3 * giE;
          if(insertion == 1)vr0 += vp1*vp3 * giE;
          if(insertion == 2)vr0 += vp1*vp2 * giE;
          
        }
      }
    }
  }

}


void baryon_corrE(EigenM &prop1, EigenM &prop2, EigenM &prop3,
   ga_M &ga2,int ind2,ga_M &ga1,int ind1,
  Evector &res, fft_desc_basic &fd,int clear=1,int imom=50505)
{
  int NTt  = fd.Nv[3];
  ////LInt Nxyz = prop1[0].size();
  int nmass = prop1.size()/(12*12*NTt);
  ////int nt = fd.nt;

  EigenM resE;
  ini_resE(resE,nmass,fd);

  qlat::vector<Complexq > G;G.resize(16);
  qlat::vector<int > mL;mL.resize(3);

  clear_qv(G);G[ind2*4 + ind1] = +1.0;
  mL[0] = 0;mL[1] = 1;mL[2] = 2;
  baryon_vectorE(prop1,prop2,prop3, ga2,ga1, G, mL, resE, fd, 1);
  clear_qv(G);G[ind2*4 + ind1] = -1.0;
  mL[0] = 1;mL[1] = 0;mL[2] = 2;
  baryon_vectorE(prop1,prop2,prop3, ga2,ga1, G, mL, resE, fd, 0);

  ////proton_vectorE(prop1,prop2,prop3,ga2,ind2,ga1,ind1,resE,fd,1);

  vec_corrE(resE,res,fd,clear,imom);
}


void Omega_corrE(EigenM &prop1, EigenM &prop2, EigenM &prop3,
   ga_M &ga2,int ind2,ga_M &ga1,int ind1,
  Evector &res, fft_desc_basic &fd,int clear=1,int imom=50505)
{
  int NTt  = fd.Nv[3];
  ///LInt Nxyz = prop1[0].size();
  int nmass = prop1.size()/(12*12*NTt);
  ///int nt = fd.nt;

  EigenM resE;
  ini_resE(resE,nmass,fd);

  qlat::vector<Complexq > G;G.resize(16);
  qlat::vector<int > mL;mL.resize(3);

  std::vector<int > dia;dia.resize(6);
  std::vector<int > sn ;sn.resize(6);
  dia[0] = 9012;sn[0] =  1;
  dia[1] = 9102;sn[1] = -1;
  dia[2] = 9021;sn[2] = -1;
  dia[3] = 9201;sn[3] =  1;
  dia[4] = 9210;sn[4] = -1;
  dia[5] = 9120;sn[5] =  1;

  for(int di=0;di<6;di++)
  {
    clear_qv(G);G[ind2*4 + ind1] = sn[di];
    mL[0] = (dia[di]/100)%10;mL[1] =  (dia[di]%100)/10;mL[2] = dia[di]%10;
    baryon_vectorE(prop1,prop2,prop3, ga2,ga1, G, mL, resE, fd, 0);
  }
     
  ////clear_qv(G);G[ind2*4 + ind1] = +1.0;
  ////mL[0] = 0;mL[1] = 1;mL[2] = 2;
  ////baryon_vectorE(prop1,prop2,prop3, ga2,ga1, G, mL, resE, fd, 0);

  ////clear_qv(G);G[ind2*4 + ind1] = -1.0;
  ////mL[0] = 1;mL[1] = 0;mL[2] = 2;
  ////baryon_vectorE(prop1,prop2,prop3, ga2,ga1, G, mL, resE, fd, 0);

  ////clear_qv(G);G[ind2*4 + ind1] = -1.0;
  ////mL[0] = 0;mL[1] = 2;mL[2] = 1;
  ////baryon_vectorE(prop1,prop2,prop3, ga2,ga1, G, mL, resE, fd, 0);

  ////clear_qv(G);G[ind2*4 + ind1] = +1.0;
  ////mL[0] = 2;mL[1] = 0;mL[2] = 1;
  ////baryon_vectorE(prop1,prop2,prop3, ga2,ga1, G, mL, resE, fd, 0);

  ////clear_qv(G);G[ind2*4 + ind1] = -1.0;
  ////mL[0] = 2;mL[1] = 1;mL[2] = 0;
  ////baryon_vectorE(prop1,prop2,prop3, ga2,ga1, G, mL, resE, fd, 0);

  ////clear_qv(G);G[ind2*4 + ind1] = +1.0;
  ////mL[0] = 1;mL[1] = 2;mL[2] = 0;
  ////baryon_vectorE(prop1,prop2,prop3, ga2,ga1, G, mL, resE, fd, 0);


  ////proton_vectorE(prop1,prop2,prop3,ga2,ind2,ga1,ind1,resE,fd,1);

  vec_corrE(resE,res,fd,clear,imom);
}

}


#endif
