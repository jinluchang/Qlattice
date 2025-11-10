// utils_corr_seq.h
// Gen Wang
// Oct. 2022

#ifndef UTILS_CORR_SEQ_H
#define UTILS_CORR_SEQ_H

#pragma once

#include "utils_float_type.h"
#include "utils_gammas.h"
#include "utils_fft_desc.h"
#include "utils_sector_funs.h"
#include "utils_reduce_vec.h"
#include "utils_grid_src.h"
#include "utils_corr_prop.h"
#include "utils_eigen_ov.h"
#include "utils_lms_funs.h"
#include "utils_corr_sparse_baryon.h"

namespace qlat{

////simple sequential sources
template <class Ty, Int civ>
void local_sequential_source(qlat::FieldM<Ty, civ>& src, const qlat::vector<Int >& tseq, const Int order_data = 0)
{
  TIMERA("local_sequential_source");
  Qassert(src.initialized);
  const Geometry& geo = src.geo();
  const long V = geo.local_volume();
  const Int Dim = qlat::get_data_size(src) / ( sizeof(Ty) * V );

  qlat::vector<Int > nv, Nv, mv;geo_to_nv(geo, nv, Nv, mv);
  long  Nvol = Nv[0]*Nv[1]*Nv[2];
  Qassert(tseq.size() > 0);
  for(long i=0;i<tseq.size();i++){
    Qassert(tseq[i] >= 0 and tseq[i] < nv[3]);
  }

  Ty* srcP = (Ty*) qlat::get_data(src).data();
  const Int Nt = Nv[3];
  const Int nt = nv[3];
  /////const long  V= Nt * Nvol;

  qacc_for(xi, long(Nvol),{
    for(Int ti = 0;ti < Nt; ti ++)
    {
      const long isp = ti*Nvol + xi;
      const Coordinate xl  = geo.coordinate_from_index(isp);
      const Coordinate pos = geo.coordinate_g_from_l(xl);
      bool find = false;
      for(long tj = 0;tj < tseq.size(); tj++){
        if(pos[3] == (tseq[tj])%nt)
        {
          find = true;break;
        }
      }
      if(find == false and order_data == 0)
      {
        for(Int ic=0;ic<Dim;ic++){srcP[isp*Dim + ic] = 0;}
      }
      if(find == false and order_data == 1)
      {
        for(Int ic=0;ic<Dim;ic++){srcP[ic * Nt * Nvol + isp] = 0;}
      }
    }
  });
}

template <class Ty, Int civ>
void local_sequential_source(qlat::FieldM<Ty, civ>& src, const Int tseq, const Int order_data = 0)
{
  std::vector<Int > tL;tL.resize(1);tL[0] = tseq;
  local_sequential_source(src, tL, order_data);
}

template <class Td>
void local_sequential_source(Propagator4dT<Td >& res, Propagator4dT<Td >& src, const qlat::vector<Int >& tseq, const Int gammai = -1)
{
  TIMERA("local_sequential_source");
  Qassert(src.initialized);
  const Geometry& geo = src.geo();
  const long V = geo.local_volume();
  //const Int Dim = src.multiplicity;
  const Int Dim = qlat::get_data_size(src) / (sizeof(qlat::ComplexT<Td >) * V) ;

  if(!res.initialized){res.init(geo);}
  qlat::vector<Int > nv, Nv, mv;geo_to_nv(geo, nv, Nv, mv);
  long  Nvol = Nv[0]*Nv[1]*Nv[2];
  ///const long Ndata = geo.local_volume() * Dim;
  Qassert(tseq.size() > 0);
  for(long i=0;i<tseq.size();i++){
    Qassert(tseq[i] >= 0 and tseq[i] < nv[3]);
  }
  Qassert(gammai >= -1 and gammai < 16);

  qlat::ComplexT<Td >* srcP = (qlat::ComplexT<Td >*) qlat::get_data(src).data();
  qlat::ComplexT<Td >* resP = (qlat::ComplexT<Td >*) qlat::get_data(res).data();
  ////svecT.shift_vecs_dir(tmp[5], tmp[4],  3, -1  );
  const Int Nt = Nv[3];
  const Int nt = nv[3];
  /////const long  V= Nt * Nvol;

  qacc_for(isp, V, {
    for(Int ic=0;ic<Dim;ic++){resP[isp*Dim + ic] = 0;}
  });
  qacc_for(xi, long(Nvol),{
    for(Int ti = 0;ti < Nt; ti ++)
    {
      const long isp = ti*Nvol + xi;
      const Coordinate xl  = geo.coordinate_from_index(isp);
      const Coordinate pos = geo.coordinate_g_from_l(xl);
      for(long tj = 0;tj < tseq.size(); tj++){
        if(pos[3] == (tseq[tj])%nt)
        {
          for(Int ic=0;ic<Dim;ic++){resP[isp*Dim + ic] = srcP[isp*Dim + ic];}
        }
      }
    }
  });

  if(gammai > -1){
    ga_matrices_cps ga_cps;
    std::vector<ga_M > gL;gL.resize(16);
    {int o=0;
    for(Int i=0;i<6;i++){gL[o] = ga_cps.ga[0][i];o+=1;}
    for(Int i=2;i<6;i++){gL[o] = ga_cps.ga[1][i];o+=1;}
    for(Int i=3;i<6;i++){gL[o] = ga_cps.ga[2][i];o+=1;}
    for(Int i=4;i<6;i++){gL[o] = ga_cps.ga[3][i];o+=1;}
    for(Int i=5;i<6;i++){gL[o] = ga_cps.ga[4][i];o+=1;}}
    ga_M& ga = gL[gammai];
    prop4d_sink_gamma(res, ga );
  }
}

}

#endif
