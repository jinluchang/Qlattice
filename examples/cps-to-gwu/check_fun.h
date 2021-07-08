// io_gwu.h
// Gen Wang
// Jan. 2021

#ifndef Check_FUN_H
#define Check_FUN_H
#pragma once

#include "general_funs.h"

namespace qlat
{

template <class T>
void diff_gauge( GaugeFieldT<T> &g0, GaugeFieldT<T> &g1)
{
  double diff = 0.0;int count_print = 0;
  for (long index = 0; index < g0.geo().local_volume(); ++index) {
    Vector<ColorMatrixT<T> > v0 = g0.get_elems(index);
    Vector<ColorMatrixT<T> > v1 = g1.get_elems(index);
    for(int m = 0; m < 4; ++m){
      const double *p0= (double*) &v0[m](0,0);
      const double *p1= (double*) &v1[m](0,0);
      for(int pi=0;pi<9*2;pi++)
      {
        diff += std::fabs(p0[pi]-p1[pi]);
        if(std::fabs(p0[pi]-p1[pi])>1e-6 and count_print < 100){
          print0("Wrong %.5e %.5e \n",p0[pi],p1[pi]);
          count_print += 1;
        }
      }
    }
  }
  sum_all_size(&diff,0);
  print0("===diff conf %.5e \n",diff);
}

template <class T>
void diff_prop(Propagator4dT<T>& p0, Propagator4dT<T>& p1, double err=1e-15)
{
  int rank = qlat::get_id_node();
  double diffp = 0.0; int countp = 0;
  for (long index = 0; index < p0.geo().local_volume(); ++index) {
    Coordinate xl0 = p0.geo().coordinate_from_index(index);
    Coordinate xg0 = p0.geo().coordinate_g_from_l(xl0);

    qlat::WilsonMatrixT<T>&  s0 =  p0.get_elem(index);
    qlat::WilsonMatrixT<T>&  s1 =  p1.get_elem(index);
    for(int d0=0;d0<12;d0++)
    for(int d1=0;d1<12;d1++)
    {
      T p0 = s0(d0,d1);
      T p1 = s1(d0,d1);

      double diff = 0.0;
      if(fabs(qlat::qnorm(p1)) > 1e-28){diff = std::fabs(qlat::qnorm((p0-p1))/qlat::qnorm(p1));}
      if(fabs(qlat::qnorm(p0)) > 1e-28){diff = std::fabs(qlat::qnorm((p0-p1))/qlat::qnorm(p0));}
      diffp += diff;
      if(diff > err and countp < 64)
      {
        printf("rank %5d, x %3d, y %3d, z %3d, t %3d, d %3d, value %.3e %.3e, %.3e %.3e, %.3e \n",
          rank ,xg0[0],xg0[1],xg0[2],xg0[3],d0*12+d1 ,p0.real(),p0.imag(),p1.real(),p1.imag(), diff);
        countp += 1;
      }
    }

  }
  sum_all_size(&diffp,1);
  MPI_Barrier(get_comm());fflush(stdout);
  print0("==prop diff %.5e \n",diffp/(p0.geo().local_volume()*12*24.0));
  MPI_Barrier(get_comm());fflush(stdout);
}

void print_src(qlat::FieldM<qlat::Complex,1> &noi)
{
  for (long index = 0; index < noi.geo().local_volume(); ++index) {
    Coordinate xl0 = noi.geo().coordinate_from_index(index);
    Coordinate xg0 = noi.geo().coordinate_g_from_l(xl0);

    qlat::Complex tem = noi.get_elems(index)[0]; 
    double sum = 0.0;
    sum += std::fabs(tem.real());
    sum += std::fabs(tem.imag());
    if(sum > 1e-6){
      printf("==src x %3d, y %3d, z %3d, t %3d, value %.5e %.5e \n",xg0[0],xg0[1],xg0[2],xg0[3],tem.real(),tem.imag());
    }

  }

}

}

#endif
