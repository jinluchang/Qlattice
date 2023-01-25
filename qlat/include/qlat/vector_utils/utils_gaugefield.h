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
    for(long isp=0;isp<gf.field.size();isp++)
    {
      gf.get_elem_offset(isp) = gf.get_elem_offset(isp) * (1/2.0) + unit;
      unitarize(gf.get_elem_offset(isp));
    }
  }
}


}

#endif
