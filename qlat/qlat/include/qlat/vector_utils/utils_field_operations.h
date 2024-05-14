// utils_field_operations.h
// Gen Wang
// Feb. 2024

#ifndef UTILS_FIELD_OPERATIONS_H
#define UTILS_FIELD_OPERATIONS_H

#pragma once

#include "utils_float_type.h"

namespace qlat{

template <class T1, int civ >
void clear_fields(qlat::FieldM<T1, civ>& pr, int GPU = 1)
{
  TIMER("clear_fields");
  Qassert(pr.initialized);
  Qassert(qlat::get_data_size(pr)%sizeof(double) == 0);
  const Long Ndata = qlat::get_data_size(pr) / sizeof(double);
  double* r0 = (double*) qlat::get_data(pr).data();
  qGPU_for(isp, Ndata, GPU, {
    r0[isp] = 0.0;
  });
}

// can be expanded ones, but only double and float
template <class T1, class T2 >
void copy_fields(T1* pr, const T2* p0, const int civ, const Geometry& geor, const Geometry& geo0)
{
  TIMER("copy_fields");
  //Qassert(geor.multiplicity == geo0.multiplicity and geo0.multiplicity == 1);// to avoid complications for pointer issues
  Qassert(IsBasicDataType<T1>::value and IsBasicDataType<T2>::value);

  using M1 = typename IsBasicDataType<T1>::ElementaryType;
  using M2 = typename IsBasicDataType<T2>::ElementaryType;

  //int mode_copy = 0;
  //int Ndata1    = 1;
  //int Ndata2    = 1;
  int Ndata1 = sizeof(T1) / sizeof(M1);
  int Ndata2 = sizeof(T2) / sizeof(M2);
  Qassert(Ndata1 == Ndata2);
  Ndata1 = Ndata1 * civ;// multiplicity

  qacc_for(isp, geor.local_volume(), {
    const Coordinate xl = geor.coordinate_from_index(isp);
    const Long dr = geor.offset_from_coordinate(xl, 1) * civ + 0;
    const Long d0 = geo0.offset_from_coordinate(xl, 1) * civ + 0;
    void* res = (void*) &pr[dr];//qlat::get_data(pr.get_elems(xl)).data();
    const void* src = (void*) &p0[d0];//qlat::get_data(p0.get_elems(xl)).data();
    copy_double_float((M1*) res, (M2*) src, Ndata1);
  });

  //if(get_data_type_is_double<T1 >()){Ndata1 = sizeof(T1)/sizeof(double);}else{Ndata1 = sizeof(T1)/sizeof(float);}
  //if(get_data_type_is_double<T2 >()){Ndata2 = sizeof(T2)/sizeof(double);}else{Ndata2 = sizeof(T2)/sizeof(float);}

  //if(get_data_type_is_double<T1 >()){Ndata1 = sizeof(T1)/sizeof(double);}else{Ndata1 = sizeof(T1)/sizeof(float);}
  //if(get_data_type_is_double<T2 >()){Ndata2 = sizeof(T2)/sizeof(double);}else{Ndata2 = sizeof(T2)/sizeof(float);}
  //Qassert(Ndata1 == Ndata2);
  //Ndata1 = Ndata1 * civ;// multiplicity

  //if( get_data_type_is_double<T1 >() and  get_data_type_is_double<T2 >()){mode_copy = 0;}
  //if( get_data_type_is_double<T1 >() and !get_data_type_is_double<T2 >()){mode_copy = 1;}
  //if(!get_data_type_is_double<T1 >() and  get_data_type_is_double<T2 >()){mode_copy = 2;}
  //if(!get_data_type_is_double<T1 >() and !get_data_type_is_double<T2 >()){mode_copy = 3;}

  //qacc_for(isp, geor.local_volume(), {
  //  const Coordinate xl = geor.coordinate_from_index(isp);
  //  const Long dr = geor.offset_from_coordinate(xl, 1) * civ + 0;
  //  const Long d0 = geo0.offset_from_coordinate(xl, 1) * civ + 0;
  //  void* res = (void*) &pr[dr];//qlat::get_data(pr.get_elems(xl)).data();
  //  void* src = (void*) &p0[d0];//qlat::get_data(p0.get_elems(xl)).data();
  //  if(mode_copy == 0){copy_double_float((double*) res, (double*) src, Ndata1);}
  //  if(mode_copy == 1){copy_double_float((double*) res, (float* ) src, Ndata1);}
  //  if(mode_copy == 2){copy_double_float((float* ) res, (double*) src, Ndata1);}
  //  if(mode_copy == 3){copy_double_float((float* ) res, (float* ) src, Ndata1);}
  //});
}

template <class T1, class T2, int civ >
void copy_fields(qlat::FieldM<T1, civ>& pr, const qlat::FieldM<T2, civ>& p0)
{
  Qassert(p0.initialized);
  const Geometry& geo = p0.geo();
  if(!pr.initialized){pr.init(geo);}

  //Geometry geor = pr.geo();
  //Geometry geo0 = p0.geo();
  //geor.multiplicity = 1;geo0.multiplicity = 1;
  T1* res = (T1*) qlat::get_data(pr).data();
  const T2* src = (T2*) qlat::get_data(p0).data();
  copy_fields<T1, T2>(res, src, civ, pr.geo(), p0.geo());
}

template <class Ty, int civ>
double fields_quick_checksum(qlat::FieldM<Ty, civ>& fs, const Long block = 128)
{
  TIMER("fields_quick_checksum");
  using D = typename IsBasicDataType<Ty>::ElementaryType;
  const Long Ndata = qlat::get_data_size(fs) / ( 2 * sizeof(D));//complex checksum
  //if(!get_data_type_is_double<Ty >()){Ndata = Ndata * 2;}
  Qassert(Ndata % block == 0);
  qlat::ComplexT<D >  res = 0.0;
  void* src = (void*) qlat::get_data(fs).data();
  res = vec_norm2((qlat::ComplexT<D >*) src, (qlat::ComplexT<D >*)src, Ndata, QMGPU, block);
  //if(get_data_type_is_double<Ty >()){
  //}
  //else{
  //  res = vec_norm2((qlat::ComplexT<float  >*) src, (qlat::ComplexT<float  >*)src, Ndata, QMGPU, block);
  //}
  print0("gauge norm %.15e %.15e \n",  res.real(), res.imag());
  return res.real();
}

template <class T1, class T2, class T3, class Ty, int civ >
void fields_operations(qlat::FieldM<T1, civ>& pr, qlat::FieldM<T2, civ>& p0, qlat::FieldM<T3, civ>& p1, 
  const Ty f0 = Ty(1.0, 0.0), const Ty f1 = Ty(1.0, 0.0), const Ty f2 = Ty(1.0, 0.0))
{
  TIMER("fields_operations");
  Qassert(p0.initialized);
  Qassert(p1.initialized);

  const Geometry& geo = p0.geo();
  if(!pr.initialized){pr.init(geo);}

  //const long Ndata = qlat::get_data_size(pr)/ sizeof(T1);
  //T1* r0 = (T1*) qlat::get_data(pr).data();
  //T2* s0 = (T2*) qlat::get_data(p0).data();
  //T3* s1 = (T3*) qlat::get_data(p1).data();
  //qacc_for(isp, Ndata, {
  //  r0[isp] = f0 * r0[isp] + f1 * s0[isp] + f2 * s1[isp];
  //});

  // geometry can be extended, but not operated
  qacc_for(isp, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(isp);
    T1* r0 = (T1*) pr.get_elems(xl).p;
    T2* s0 = (T2*) p0.get_elems(xl).p;
    T3* s1 = (T2*) p1.get_elems(xl).p;
    for (int m = 0; m < civ; ++m) {
      r0[m] = f0 * r0[m] + f1 * s0[m] + f2 * s1[m];
      //r0[isp] = f0 * r0[isp] + f1 * s0[isp] + f2 * s1[isp];
      //res[m] = src[m];
    }
  });
}

template <class T1, class T2, class T3, class Ty, int civ >
void fields_operations(std::vector<qlat::FieldM<T1, civ> >& pr, std::vector<qlat::FieldM<T2, civ> >& p0, std::vector<qlat::FieldM<T3, civ> >& p1, const Ty f0 = Ty(1.0, 0.0), const Ty f1 = Ty(1.0, 0.0), const Ty f2 = Ty(1.0, 0.0))
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
  fields_operations(pr, p0, p1, qlat::ComplexT<double>( 0.0, 0.0), qlat::ComplexT<double>( 1.0, 0.0), qlat::ComplexT<double>( 1.0, 0.0));
}

template <class T1, class T2, class T3, int civ >
void fields_subtractions(std::vector<qlat::FieldM<T1, civ> >& pr, std::vector<qlat::FieldM<T2, civ> >& p0, std::vector<qlat::FieldM<T3, civ> >& p1)
{
  fields_operations(pr, p0, p1, qlat::ComplexT<double>( 0.0, 0.0), qlat::ComplexT<double>( 1.0, 0.0), qlat::ComplexT<double>(-1.0, 0.0));
}

template <class T1, class T2, int civ >
void fields_equal(std::vector<qlat::FieldM<T1, civ> >& pr, std::vector<qlat::FieldM<T2, civ> >& ps)
{
  fields_operations(pr, ps, ps, qlat::ComplexT<double>( 0.0, 0.0), qlat::ComplexT<double>( 1.0, 0.0), qlat::ComplexT<double>( 0.0, 0.0));
}

}


#endif

