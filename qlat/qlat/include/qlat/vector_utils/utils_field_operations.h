// utils_field_operations.h
// Gen Wang
// Feb. 2024

#ifndef UTILS_FIELD_OPERATIONS_H
#define UTILS_FIELD_OPERATIONS_H

#pragma once

#include "utils_float_type.h"

namespace qlat{

template <class Fieldy >
void clear_fields(Fieldy& pr, int GPU = 1)
{
  TIMER("clear_fields");
  qassert(pr.initialized);
  const Long Nd = GetFieldSize(pr);
  qassert(Nd % sizeof(double) == 0);
  const Long Ndata = Nd / sizeof(double);
  double* r0 = (double*) get_data(pr).data();
  qGPU_for(isp, Ndata, GPU, {
    r0[isp] = 0.0;
  });
}

// can be expanded ones, but only double and float without expanded parts
template <class T1, class T2 >
void copy_fields(T1* pr, const T2* p0, const int civ, const Geometry& geor, const Geometry& geo0)
{
  TIMER("copy_fields");
  qassert(IsBasicTypeReal<T1>() and IsBasicTypeReal<T2>());

  using M1 = typename IsBasicDataType<T1>::ElementaryType;
  using M2 = typename IsBasicDataType<T2>::ElementaryType;

  //int mode_copy = 0;
  //int Ndata1    = 1;
  //int Ndata2    = 1;
  int Ndata1 = sizeof(T1) / sizeof(M1);
  int Ndata2 = sizeof(T2) / sizeof(M2);
  qassert(Ndata1 == Ndata2);
  Ndata1 = Ndata1 * civ;// multiplicity

  qacc_for(isp, geor.local_volume(), {
    const Coordinate xl = geor.coordinate_from_index(isp);
    const Long dr = geor.offset_from_coordinate(xl, 1) * civ + 0;
    const Long d0 = geo0.offset_from_coordinate(xl, 1) * civ + 0;
    void* res = (void*) &pr[dr];//qlat::get_data(pr.get_elems(xl)).data();
    const void* src = (void*) &p0[d0];//qlat::get_data(p0.get_elems(xl)).data();
    copy_double_float((M1*) res, (M2*) src, Ndata1);
  });

}

template <class T1, class T2 >
void copy_fields(qlat::Field<T1>& pr, const qlat::Field<T2>& p0)
{
  qassert(p0.initialized);
  if(!pr.initialized){
    pr.init(p0.geo(), p0.multiplicity);
  }

  T1* res = (T1*) qlat::get_data(pr).data();
  const T2* src = (T2*) qlat::get_data(p0).data();
  qassert(pr.multiplicity == p0.multiplicity);
  const int civ = pr.multiplicity;
  copy_fields<T1, T2>(res, src, civ, pr.geo(), p0.geo());
}

template <class T1, class T2 >
void copy_fieldsG(qlat::FieldG<T1>& pr, const qlat::FieldG<T2>& p0)
{
  copy_fields(pr, p0);
  pr.mem_order = p0.mem_order;
}

template <class T1, class T2 >
void copy_fieldsG(SelectedFieldG<T1>& pr, const SelectedFieldG<T2>& p0)
{
  qassert(p0.initialized);
  if(!pr.initialized){
    pr.init(p0.geo(), p0.n_elems, p0.multiplicity, get_type_mem(p0.field_gpu.GPU));
  }

  T1* res       = (T1*) qlat::get_data(pr).data();
  const T2* src = (T2*) qlat::get_data(p0).data();
  const Long Nd = p0.field.size();
  cpy_GPU(res, src, Nd, pr.field_gpu.GPU, p0.field_gpu.GPU);
  pr.mem_order = p0.mem_order;
}

// norm2 = f^\dagger f 
template <class Fieldy>
double fields_quick_checksum(Fieldy& fs, const Long block = 128, const bool log = true)
{
  TIMER("fields_quick_checksum");
  qassert(fs.initialized);
  //Ty --> double float ...
  qassert(GetBasicDataType<Fieldy>::get_type_name() != std::string("unknown_type"));
  using D = typename GetBasicDataType<Fieldy>::ElementaryType;
  qassert(IsBasicTypeReal<D>());

  const Long Nd = GetFieldSize(fs);
  qassert(Nd % (2 * sizeof(D)) == 0);
  const Long Ndata = Nd / ( 2 * sizeof(D));
  qassert(Ndata % block == 0);
  qlat::ComplexT<D >  res = 0.0;
  void* src = (void*) qlat::get_data(fs).data();
  res = vec_norm2((qlat::ComplexT<D >*) src, (qlat::ComplexT<D >*)src, Ndata, QMGPU, block);
  if(log){
    print0("gauge norm %.15e %.15e \n",  res.real(), res.imag());
  }
  return res.real();
}

template <class T1, class T2, class T3, class Ty >
void fields_operations(qlat::Field<T1>& pr, qlat::Field<T2>& p0, qlat::Field<T3>& p1, 
  const Ty f0 = Ty(1.0, 0.0), const Ty f1 = Ty(1.0, 0.0), const Ty f2 = Ty(1.0, 0.0))
{
  TIMER("fields_operations");
  qassert(p0.initialized);
  qassert(p1.initialized);

  const Geometry& geo = p0.geo();
  const int civ = p0.multiplicity;
  if(!pr.initialized){pr.init(geo, civ);}
  qassert(civ == pr.multiplicity and civ == p1.multiplicity);

  // geometry can be extended, but not operated on local numbers
  qacc_for(isp, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(isp);
    T1* r0 = (T1*) pr.get_elems(xl).p;
    T2* s0 = (T2*) p0.get_elems(xl).p;
    T3* s1 = (T2*) p1.get_elems(xl).p;
    for (int m = 0; m < civ; ++m) {
      r0[m] = f0 * r0[m] + f1 * s0[m] + f2 * s1[m];
    }
  });
}

template <class T1, class T2, class T3, class Ty, int civ >
void fields_operations(std::vector<qlat::FieldM<T1, civ> >& pr, std::vector<qlat::FieldM<T2, civ> >& p0, std::vector<qlat::FieldM<T3, civ> >& p1, const Ty f0 = Ty(1.0, 0.0), const Ty f1 = Ty(1.0, 0.0), const Ty f2 = Ty(1.0, 0.0))
{
  qassert(p0.size() == p1.size());
  if(p0.size() == 0){return ;}
  //const Geometry& geo = p0[0].geo();
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

template <class T1 >
void fields_conj(T1** src, const int nvec, const int civ, const Geometry& geo)
{
  TIMER("fields_conj");
  qassert(IsBasicDataType<T1>::value);

  using M1 = typename IsBasicDataType<T1>::ElementaryType;
  int Ndata = sizeof(T1) / sizeof(M1);

  Ndata = Ndata * civ;// multiplicity

  for(int iv=0;iv<nvec;iv++){
    M1* res = (M1*) src[iv];
    qacc_for(isp, geo.local_volume(), {
      M1* p = &res[isp*Ndata + 0];
      for(int d=0;d<Ndata;d++)
      {
        p[d] = qlat::qconj(p[d]);
      }
    });
  }
}

template <class T1, int civ >
void fields_conj(qlat::FieldM<T1, civ>& pr)
{
  qassert(pr.initialized);
  qlat::vector_acc<T1* > src;src.resize(1);
  src[0] = (T1*) qlat::get_data(pr).data();
  fields_conj(src.data(), 1, civ, pr.geo());
}

template <class T1, int civ >
void fields_conj(std::vector<qlat::FieldM<T1, civ> >& pr)
{
  const int nvec = pr.size();
  if(nvec == 0){return ;}
  qlat::vector_acc<T1* > src;src.resize(nvec);
  for(int si=0;si<nvec;si++)
  {
    qassert(pr[si].initialized);
    src[si] = (T1*) qlat::get_data(pr[si]).data();
  }
  fields_conj(src.data(), nvec, civ, pr[0].geo());
}

template <class M >
void switch_orders(M& f, const QMEM_ORDER order_){
  if(f.mem_order == order_){return ;}
  qassert(order_ == QLAT_DEFAULT or order_== QLAT_OUTTER);
  qassert(f.mem_order == QLAT_DEFAULT or f.mem_order == QLAT_OUTTER);
  move_index mv_idx;
  const Long Nd  = f.get_nsites();
  const Long Dim = f.multiplicity;
  if(f.mem_order == QLAT_DEFAULT and order_ == QLAT_OUTTER){
    mv_idx.move_civ_out(f.field_gpu.p, f.field_gpu.p, 1, Nd, Dim, 1, f.field_gpu.GPU);
    f.mem_order = QLAT_OUTTER;
    return ;
  }

  if(f.mem_order == QLAT_OUTTER and order_ == QLAT_DEFAULT){
    mv_idx.move_civ_in( f.field_gpu.p, f.field_gpu.p, 1, Dim, Nd, 1, f.field_gpu.GPU);
    f.mem_order = QLAT_DEFAULT;
    return ;
  }
}


}


#endif

