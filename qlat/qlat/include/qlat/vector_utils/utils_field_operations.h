// utils_field_operations.h
// Gen Wang
// Feb. 2024

#ifndef UTILS_FIELD_OPERATIONS_H
#define UTILS_FIELD_OPERATIONS_H

#pragma once

#include "utils_float_type.h"

namespace qlat{

template <class T1, class T2, class T3, int civ >
void fields_operations(qlat::FieldM<T1, civ>& pr, qlat::FieldM<T2, civ>& p0, qlat::FieldM<T3, civ>& p1, 
  const T1 f0 = T1(1.0, 0.0), const T1 f1 = T1(1.0, 0.0), const T1 f2 = T1(1.0, 0.0))
{
  Qassert(p0.initialized);
  Qassert(p1.initialized);

  const Geometry& geo = p0.geo();
  if(!pr.initialized){pr.init(geo);}

  const long Ndata = qlat::get_data_size(pr)/ sizeof(T1);
  T1* r0 = (T1*) qlat::get_data(pr).data();
  T2* s0 = (T2*) qlat::get_data(p0).data();
  T3* s1 = (T3*) qlat::get_data(p1).data();
  qacc_for(isp, Ndata, {
    r0[isp] = f0 * r0[isp] + f1 * s0[isp] + f2 * s1[isp];
  });
}

template <class T1, class T2, class T3, int civ >
void fields_operations(std::vector<qlat::FieldM<T1, civ> >& pr, std::vector<qlat::FieldM<T2, civ> >& p0, std::vector<qlat::FieldM<T3, civ> >& p1, const T1 f0 = T1(1.0, 0.0), const T1 f1 = T1(1.0, 0.0), const T1 f2 = T1(1.0, 0.0))
{
  Qassert(p0.size() == p1.size());
  if(p0.size() == 0){return ;}
  const Geometry& geo = p0[0].geo();
  const int Nvec = p0.size();
  if(pr.size() != Nvec){pr.resize(Nvec);}
  for(int vi=0;vi<Nvec;vi++){
    fields_operations(pr[vi], p0[vi], p1[vi], f0, f1, f2);
  }
}

template <class T1, class T2, class T3, int civ >
void fields_additions(std::vector<qlat::FieldM<T1, civ> >& pr, std::vector<qlat::FieldM<T2, civ> >& p0, std::vector<qlat::FieldM<T3, civ> >& p1)
{
  fields_operations(pr, p0, p1, T1( 0.0, 0.0), T2( 1.0, 0.0), T2( 1.0, 0.0));
}

template <class T1, class T2, class T3, int civ >
void fields_subtractions(std::vector<qlat::FieldM<T1, civ> >& pr, std::vector<qlat::FieldM<T2, civ> >& p0, std::vector<qlat::FieldM<T3, civ> >& p1)
{
  fields_operations(pr, p0, p1, T1( 0.0, 0.0), T2( 1.0, 0.0), T2(-1.0, 0.0));
}

template <class T1, class T2, int civ >
void fields_equal(std::vector<qlat::FieldM<T1, civ> >& pr, std::vector<qlat::FieldM<T2, civ> >& ps)
{
  fields_operations(pr, ps, ps, T1( 0.0, 0.0), T2( 1.0, 0.0), T2( 0.0, 0.0));
}

}


#endif

