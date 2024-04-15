// utils_GAUGEFIELD.h
// Gen Wang
// Mar. 2022

#ifndef UTILS_GAUGEFIELD_H
#define UTILS_GAUGEFIELD_H
#pragma once

#include "general_funs.h"
#include <qlat/qcd-topology.h>

//////TODO

namespace qlat
{

/////E is 3x3 array
template <typename Ty >
qacc void su3_one(Ty* E)
{
  for(int i=0;i<9;i++){ E[i] = 0.0; }
  E[0*3 + 0] = 1.0;
  E[1*3 + 1] = 1.0;
  E[2*3 + 2] = 1.0;
}

qacc Long su3_n(const Geometry& geo, const Coordinate& xl)
{
  const int nD = 9;
  return geo.offset_from_coordinate(xl)*nD + 0;
}

qacc Long su3_n(const Geometry& geo, const Coordinate& xl, const int mu)
{
  const int Dim = 4;
  const int nD  = 9;
  qassert(mu >= 0 and mu < Dim);
  return (geo.offset_from_coordinate(xl) * Dim + mu)*nD + 0;
  //if(mu > 0){
  //  return (geo.offset_from_coordinate(xl) * Dim + mu)*nD + 0;
  //}else{
  //  return (geo.offset_from_coordinate(xl)           )*nD + 0;
  //}
}

template <typename Ty, bool clear, bool dag1, bool dag2>
qacc void su3_multi_kerB(Ty* res, Ty* Q1, Ty* Q2, Ty* BUF)
{
  for(int a=0;a<3;a++)
  for(int b=0;b<3;b++)
  {
    Ty& buf = BUF[b*3 + a];buf = 0.0;
    for(int d=0;d<3;d++)
    {
      if(dag1 == false and dag2 == false){buf += Q1[b*3 + d]*Q2[d*3 + a];}
      if(dag1 == false and dag2 == true ){buf += Q1[b*3 + d]*qlat::qconj(Q2[a*3 + d]);}
      if(dag1 == true  and dag2 == false){buf += qlat::qconj(Q1[d*3 + b])*Q2[d*3 + a];}
      if(dag1 == true  and dag2 == true ){buf += qlat::qconj(Q1[d*3 + b])*qlat::qconj(Q2[a*3 + d]);}
    }
  }
  if( clear){for(int a=0;a<9;a++){res[a]   = BUF[a];}}
  if(!clear){for(int a=0;a<9;a++){res[a]  += BUF[a];}}
}

/////3x3 inner 3 is faster and 1,2,3 (row 0) --> 4,5,6 (row 1) --> 7,8,9 (row 2)
template <typename Ty, bool clear, bool dag1, bool dag2>
qacc void su3_multi_ker(Ty* res, Ty* Q1, Ty* Q2)
{
  QLAT_ALIGN(QLAT_ALIGNED_BYTES) Ty BUF[9];
  su3_multi_kerB<Ty, clear, dag1, dag2>(res, Q1, Q2, BUF);
}

template <typename Ty, bool clear>
qacc void su3_multi(Ty* res, Ty* Q1, Ty* Q2)
{
  su3_multi_ker<Ty, true, false, false>(res, Q1, Q2);
}

template <typename Ty>
qacc void su3_dagger(Ty* res, Ty* src)
{
  for(int a=0;a<3;a++)
  for(int b=0;b<3;b++)
  {
    res[a*3 + b] = qlat::qconj(src[b*3 + a]);
  }
}

template <class Ty, Long N>
qacc void normalize_array_c(Ty* src)
{
  RealD norm = 0.0;
  for(Long i=0;i<N;i++){norm = qlat::qconj(src[i]) * src[i];}
  norm = std::sqrt(norm);
  if (not(norm == 1.0))
  for(Long i=0;i<N;i++){src[i] = src[i] / Ty(norm, 0.0);}
}

template <class Ty, Long N>
qacc void orthogonalize_array_c(Ty* p2, Ty* p1)
{
  qlat::ComplexT<RealD > c = 0.0;
  for (Long i = 0; i < N; ++i) {
    c += qconj(p1[i]) * p2[i];
  }
  if (not(c.real() == 0.0)) {
    for (Long i = 0; i < N; ++i) {
      p2[i] -= c * p1[i];
    }
  }
}

template <class Ty>
qacc void cross_product_conj_c(Ty* v3, Ty* v1, Ty* v2)
// v3 = ( v1 x v2 )^*
{
  v3[0] = qconj(v1[1] * v2[2] - v1[2] * v2[1]);
  v3[1] = qconj(v1[2] * v2[0] - v1[0] * v2[2]);
  v3[2] = qconj(v1[0] * v2[1] - v1[1] * v2[0]);
}

template <typename Ty>
qacc void su3_unitarize(Ty* src)
{
  Ty* p1 = src[0 * 3 + 0];// may need to rotate index?
  Ty* p2 = src[1 * 3 + 0];
  Ty* p3 = src[2 * 3 + 0];
  normalize_array_c<Ty, 3>(p1);
  orthogonalize_array_c<Ty, 3>(p2, p1);
  normalize_array_complex<Ty, 3>(p2);
  cross_product_conj_c(p3, p1, p2);
}

template <class Td>
void set_rand_link(GaugeFieldT<Td> &gf, const int seed = -1)
{
  if(seed == -1)
  {
    qacc_for(isp, gf.field.size(), { set_unit(gf.get_elem_offset(isp), 1.0);});
  }else{
    //T* res = (T*) gf.get_elem_offset(0).p;
    const Geometry& geo = gf.geo();
    qlat::ComplexT<Td>* res = (qlat::ComplexT<Td>*) qlat::get_data(gf).data();
    random_Ty(res, geo.local_volume()*geo.multiplicity*sizeof(ColorMatrixT<Td>)/(sizeof(Td)*2), 1, seed);

    //qacc_for(isp, gf.field.size(), { set_unit(gf.get_elem_offset(isp), 1.0);});
    ColorMatrixT<Td> unit;set_unit(unit, 1.0);
    /////TODO This function cannot be done on GPU
    /////Eigen normalize/normalized problem 
    for(Long isp=0;isp<gf.field.size();isp++)
    {
      gf.get_elem_offset(isp) = gf.get_elem_offset(isp) * (1/2.0) + unit;
      unitarize(gf.get_elem_offset(isp));
    }
  }
}

template <class Ta, class Td>
void copy_gf(GaugeFieldT<Ta> &g1, GaugeFieldT<Td> &g0)
{
  TIMER("copy_gf");
  const Geometry geo = g0.geo();
  if(!g1.initialized){g1.init(geo);}

  ///cannot use this due to extended gauge fields
  //const long Ndata = geo.local_volume() * 4 * 9;
  //qlat::ComplexT<Td >* src = (qlat::ComplexT<Td >*) qlat::get_data(g0).data();
  //qlat::ComplexT<Td >* res = (qlat::ComplexT<Td >*) qlat::get_data(g1).data();
  //qlat::cpy_GPU(res, src, Ndata,   1, 1);

  qacc_for(isp, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(isp);
    qlat::ComplexT<Ta >* res = (qlat::ComplexT<Ta >*) g1.get_elem(xl, 0).p;
    qlat::ComplexT<Td >* src = (qlat::ComplexT<Td >*) g0.get_elem(xl, 0).p;
    for(int i=0;i<9*4;i++){res[i] = src[i];}
    //for (int dir = 0; dir < 4; ++dir)
    //{
    //  g1.get_elem(xl, dir) = g0.get_elem(xl, dir);
    //}
  });
}

}

#endif
