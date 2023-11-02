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

/////3x3 inner 3 is faster and 1,2,3 (row 0) --> 4,5,6 (row 1) --> 7,8,9 (row 2)
template <typename Ty, bool clear, bool dag1, bool dag2>
qacc void su3_multi_ker(Ty* res, Ty* Q1, Ty* Q2)
{
  for(int a=0;a<3;a++)
  for(int b=0;b<3;b++)
  {
    Ty buf = 0.0;
    for(int d=0;d<3;d++)
    {
      if(dag1 == false and dag2 == false){buf += Q1[b*3 + d]*Q2[d*3 + a];}
      if(dag1 == false and dag2 == true ){buf += Q1[b*3 + d]*qlat::qconj(Q2[a*3 + d]);}
      if(dag1 == true  and dag2 == false){buf += qlat::qconj(Q1[d*3 + b])*Q2[d*3 + a];}
      if(dag1 == true  and dag2 == true ){buf += qlat::qconj(Q1[d*3 + b])*qlat::qconj(Q2[a*3 + d]);}
    }
    if( clear){res[b*3 + a]  = buf;}
    if(!clear){res[b*3 + a] += buf;}
  }
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

template <class Td>
void copy_gf(GaugeFieldT<Td> &g1, GaugeFieldT<Td> &g0)
{
  const Geometry geo = g0.geo();
  if(!g1.initialized){g1.init(geo);}

  qacc_for(isp, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(isp);
    for (int dir = 0; dir < 4; ++dir)
    {
      g1.get_elem(xl, dir) = g0.get_elem(xl, dir);
    }
  });
}

}

#endif
