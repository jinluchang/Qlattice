// utils_GAUGEFIELD.h
// Gen Wang
// Mar. 2022

#ifndef UTILS_GAUGEFIELD_H
#define UTILS_GAUGEFIELD_H
#pragma once

#include "general_funs.h"
#include "utils_copy_data.h"

//////TODO

namespace qlat
{

template <class T>
void set_rand_link(GaugeFieldT<T> &gf, const int seed = -1)
{
  if(seed == -1)
  {
    qacc_for(isp, gf.field.size(), { set_unit(gf.get_elem(isp), 1.0);});
  }else{
    //T* res = (T*) gf.get_elem(0).p;
    const Geometry& geo = gf.geo();
    T* res = (T*) qlat::get_data(gf).data();
    random_Ty(res, geo.local_volume()*geo.multiplicity*sizeof(ColorMatrixT<T>)/sizeof(T), 1, seed);

    //qacc_for(isp, gf.field.size(), { set_unit(gf.get_elem(isp), 1.0);});
    ColorMatrixT<T> unit;set_unit(unit, 1.0);
    /////TODO This function cannot be done on GPU
    /////Eigen normalize/normalized problem 
    for(long isp=0;isp<gf.field.size();isp++)
    {
      gf.get_elem(isp) = gf.get_elem(isp) * (1/2.0) + unit;
      unitarize(gf.get_elem(isp));
    }
  }
}


}

#endif
