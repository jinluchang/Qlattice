#ifndef utils_low_rho_h
#define utils_low_rho_h
#pragma once


#include <qlat/qlat.h>
#include "io_gwu.h"
namespace qlat{

void get_low_rho(std::vector<qlat::FermionField4dT<qlat::ComplexF> > &eigen,const std::vector<double >  &values,std::vector<Ftype > &Mres,const qlat::Geometry &geo){
  TIMER("Sum vectors");

  qlat::FermionField4dT<qlat::ComplexF> sum_gpu;sum_gpu.init(geo);

  qlat::FermionField4dT<qlat::ComplexF> sum_cpu;sum_cpu.init(geo);

  {TIMER("Sum cpu");
  for(int i=0;i<5;i++)
  for(size_t iv=0;iv<eigen.size();iv++){
  const qlat::FermionField4dT<qlat::ComplexF> &ei = eigen[iv];
  #pragma omp parallel for
  for(int index=0;index<geo.local_volume();index++)
  {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<qlat::WilsonVectorT<qlat::ComplexF > > vs = sum_cpu.get_elems(xl);
    const Vector<qlat::WilsonVectorT<qlat::ComplexF > > v = ei.get_elems_const(xl);
    vs += v;
    vs *= std::cos(i);
  }
  }
  }

  {TIMER("Sum gpu 0 ");
  for(int i=0;i<1;i++)
  for(size_t iv=0;iv<eigen.size();iv++){
  const qlat::FermionField4dT<qlat::ComplexF> &ei = eigen[iv];
  qacc_for(index, geo.local_volume(), {
    //float* p0 = (float*) (qlat::get_data(sum_gpu).data());
    //float* p1 = (float*) (qlat::get_data(ei).data());
    //for(int c=0;c<24;c++)p0[index*24+c]+=p1[index*24 + c];
    //for(int c=0;c<24;c++)p0[index*24+c]*=std::cos(i);

    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<qlat::WilsonVectorT<qlat::ComplexF > > vs = sum_gpu.get_elems(xl);
    const Vector<qlat::WilsonVectorT<qlat::ComplexF > > v = ei.get_elems_const(xl);
    vs += v;
    vs *= std::cos(i);
  });
  }
  }

  sum_gpu.init(geo);
  {TIMER("Sum gpu 1");
  for(int i=0;i<5;i++)
  for(size_t iv=0;iv<eigen.size();iv++){
  const qlat::FermionField4dT<qlat::ComplexF> &ei = eigen[iv];
  qacc_for(index, geo.local_volume(), {
    //float* p0 = (float*) (qlat::get_data(sum_gpu).data());
    //float* p1 = (float*) (qlat::get_data(ei).data());
    //for(int c=0;c<24;c++)p0[index*24+c]+=p1[index*24 + c];
    //for(int c=0;c<24;c++)p0[index*24+c]*=std::cos(i);

    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<qlat::WilsonVectorT<qlat::ComplexF > > vs = sum_gpu.get_elems(xl);
    const Vector<qlat::WilsonVectorT<qlat::ComplexF > > v = ei.get_elems_const(xl);
    vs += v;
    vs *= std::cos(i);
  });
  }
  }


  double sum_diff = 0.0;
  #pragma omp parallel for reduction(+: sum_diff)
  for(int index=0;index<geo.local_volume();index++)
  {
    const Coordinate xl = geo.coordinate_from_index(index);
    float* p0 = (float*) (qlat::get_data(sum_cpu.get_elems(xl)).data());
    float* p1 = (float*) (qlat::get_data(sum_gpu.get_elems(xl)).data());
    for(int i=0;i<24;i++)sum_diff += (p0[i]-p1[i])*(p0[i]-p1[i]);
  };
 
 sum_all_size(&sum_diff,1);
 print0("====sum diff gpu cpu %.3e \n",sum_diff);

}

}

#endif

