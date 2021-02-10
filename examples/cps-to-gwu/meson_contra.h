// io_gwu.h
// Gen Wang
// Feb. 2021

#ifndef meson_contra_h
#define meson_contra_h
#pragma once

#include <qlat/qlat.h>
#include "general_funs.h"

namespace qlat
{

void get_corr_pion(std::vector<qlat::FermionField4dT<qlat::Complex> > &prop,const Coordinate &x_ini, std::vector<double > &write ){

  const qlat::Geometry &geo = prop[0].geo();

  unsigned long Nvol = geo.local_volume();
  int Nt = geo.node_site[3];
  long Nsum = Nvol/Nt;
  int tini = x_ini[3];

  std::vector<qlat::Complex > res;res.resize(Nvol); 

  qacc_for(isp, Nvol,{ 
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
  write.resize(2*nt);set_zero(write);

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

}

#endif
