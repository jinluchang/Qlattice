// utils_stagger_contractions.h
// Gen Wang
// Jan. 2022

#ifndef UTILS_STAGGER_CONTRACTIONS_H
#define UTILS_STAGGER_CONTRACTIONS_H

#pragma once

#include "utils_float_type.h"
#include "utils_gammas.h"
#include "utils_fft_desc.h"
#include "utils_reduce_vec.h"
#include "utils_grid_src.h"
#include "utils_construction.h"

namespace qlat{

template<typename Ty , typename Td>
struct stag_inv_buf{

  std::vector<std::vector< colorFT> > bufV;
  qlat::vector_gpu<Ty > tmp;
  std::vector<qlat::vector_gpu<Ty >  > propE;
  std::vector<qlat::vector_gpu<Ty >  > propV;
  std::vector<qlat::vector_gpu<Ty >  > propB;
  std::vector< qlat::vector_gpu<Ty > > propS;
  std::vector<qlat::vector_gpu<Ty >  > prop_shift;
  qlat::vector_gpu<Ty > buf_vec;

  std::vector< GaugeFieldT<Td >   > gfL;

  qlat::FieldM<char, 1> eo;

  //qlat::vector_gpu<Ty > prop_src_gpu;
  int prop_load_src;
  std::vector<qlat::vector_gpu<Ty > > propS_cpu;
  std::vector<qlat::vector_gpu<Ty > > propS_smear;
  std::vector<int > prop_load;
  std::vector<int > prop_load_smear;

  inline void free_buf(){

    for(unsigned int i=0;i<bufV.size();i++)
    {
      bufV[i].resize(0);
    }
    bufV.resize(0);
    tmp.resize(0);
    propE.resize(0);
    propV.resize(0);
    propB.resize(0);
    propS.resize(0);
    prop_shift.resize(0);
    gfL.resize(0);
    buf_vec.resize(0);
    propS_cpu.resize(0);
    prop_load.resize(0);
    ///prop_src_gpu.resize(0);
  }

  inline void init(){
    free_buf();
  }
};

////tr[cf cf^\dagger]
template <class Ty>
void cf_simple_pion(std::vector<colorFT >& cf0, std::vector<colorFT >& cf1, EigenV &corr, qlat::fft_desc_basic &fd,int clear=1, bool print=false, double factor = 1.0)
{
  TIMER("cf_simple_pion");
  qassert(cf0.size() == 3);qassert(cf1.size() == 3);
  int  NTt  = fd.Nv[3];
  ////LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  const Geometry& geo = cf0[0].geo();

  EigenV resV;ini_resE(resV, 1, fd);
  if(resV.size()%NTt !=0 or resV.size()==0){print0("Size of res wrong. \n");qassert(false);}

  const int Dim = 3;
  qlat::vector_acc<Ty* > d0;d0.resize(Dim);
  qlat::vector_acc<Ty* > d1;d1.resize(Dim);
  for(int c=0;c<Dim;c++){d0[c] = (Ty*) qlat::get_data(cf0[c]).data();d1[c] = (Ty*) qlat::get_data(cf1[c]).data();}

  qacc_for(isp, geo.local_volume(), {
    for(int c0=0;c0<Dim;c0++)
    for(int c1=0;c1<Dim;c1++)
    {
      resV[isp] += d0[c0][isp*Dim + c1] * qlat::qconj(d1[c0][isp*Dim + c1]);
    }

  });

  vec_corrE(resV, corr, fd, clear);
  if(print )
  for(int ti=0;ti<fd.nt;ti++)
  {
    auto v0 = corr[ti] * factor;
    print0("ti %5d , cf cf^dagger %+.8e  %+.8e \n", ti, v0.real(), v0.imag());
  }

}

inline void qlat_map_eo_site(qlat::FieldM<char, 1>& eo, const Geometry& geo)
{
  if(eo.initialized)if(eo.geo() == geo){return ;}
  eo.init(geo);
  char* res = (char*) qlat::get_data(eo).data();
  ////only bool is not write thread safe
  qacc_for(isp, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(isp);
    int site_eo = (xl[0] + xl[1] + xl[2] + xl[3]) % 2;
    res[isp] = site_eo;
  });
}

/////src and res can be tthe same pointer
template <class Ty, int civ>
void apply_eo_sign(Ty* sP, Ty* rP, qlat::FieldM<char, 1>& eo, const char dir = 1)
{
  TIMER("apply_eo_sign");
  const Geometry& geo = eo.geo();
  qassert(eo.initialized);
  char* eP = (char*) qlat::get_data(eo).data();
  ///////DATA_TYPE typenum = get_data_type<Ty >();
  qacc_for(isp, geo.local_volume(), {
    qlat::Complex sign = qlat::Complex(-1.0 * dir *(eP[isp]*2 - 1), 0);
    for(int ic=0;ic<civ;ic++){rP[isp*civ+ic] = sign * sP[isp*civ+ic];}
  });
}

/////src and res can be tthe same pointer
template <class Ty, int civ>
void apply_eo_sign(qlat::FieldM<Ty , civ>& src, qlat::FieldM<Ty , civ>& res, qlat::FieldM<char, 1>& eo, const char dir = 1)
{
  TIMER("apply_eo_sign");
  if(!src.initialized or !res.initialized){abort_r("src should be initialized with geo!\n");}
  const Geometry& geo = src.geo();
  if(!eo.initialized){qlat_map_eo_site(eo, geo);}
  Ty*   sP = (Ty*  ) qlat::get_data(src).data();
  Ty*   rP = (Ty*  ) qlat::get_data(res).data();
  apply_eo_sign<Ty, civ>(sP, rP, eo, dir);
}

template <class Ty, int civ>
void apply_eo_sign(std::vector<qlat::FieldM<Ty , civ> >& src, std::vector<qlat::FieldM<Ty , civ> >& res, qlat::FieldM<char, 1>& eo)
{
  qassert(src.size() == res.size());
  for(unsigned int si = 0; si<src.size(); si++){
    apply_eo_sign<Ty, civ>(src[si], res[si], eo);
  }
}



}

#endif

