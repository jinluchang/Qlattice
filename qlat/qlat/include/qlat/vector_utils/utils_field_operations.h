// utils_field_operations.h
// Gen Wang
// Feb. 2024

#ifndef UTILS_FIELD_OPERATIONS_H
#define UTILS_FIELD_OPERATIONS_H

#pragma once

#include "utils_float_type.h"

namespace qlat{

template <class Fieldy >
void clear_fields(Fieldy& pr, Int GPU = 1)
{
  TIMER("clear_fields");
  Qassert(pr.initialized);
  const Long Nd = GetFieldSize(pr);
  Qassert(Nd % sizeof(RealD) == 0);
  const Long Ndata = Nd / sizeof(RealD);
  double* r0 = (double*) get_data(pr).data();
  qGPU_for(isp, Ndata, GPU, {
    r0[isp] = 0.0;
  });
}

template<class Fieldy>
inline bool is_field_local(Fieldy& field){
  Qassert(field.initialized);
  const Geometry& geo = field.geo();
  Geometry geo_l;Get_geo_local(geo, geo_l);
  if(geo == geo_l){return true;}
  return false;
}

template<class Fieldy>
void init_fields(std::vector<Fieldy >& res, Fieldy& src, const unsigned int size)
{
  Qassert(src.initialized);
  if(res.size() != size){
    res.resize(0);res.resize(size);
  }
  for(unsigned int i=0;i<res.size();i++){
    if(!res[i].initialized){res[i].init(src.geo(), src.multiplicity);}
  }

}

template<class Fieldy>
void init_fieldsG(std::vector<Fieldy >& res, Fieldy& src, const unsigned int size)
{
  init_fields(res, src, size);
  for(unsigned int i=0;i<size;i++){
    res[i].mem_order = src.mem_order;
  }
}

template<class Fieldy>
void init_fieldsG(Fieldy & res, Fieldy& src)
{
  Qassert(src.initialized);
  if(!res.initialized or res.multiplicity == src.multiplicity){
    res.init(src.geo(), src.multiplicity);
  }
  res.mem_order = src.mem_order;
}

// can be expanded ones, but only double and float without expanded parts
template <class T1, class T2 >
void copy_fields(T1* pr, const T2* p0, const Int civ, const Geometry& geor, const Geometry& geo0)
{
  TIMER("copy_fields");
  Qassert(IsBasicTypeReal<T1>() and IsBasicTypeReal<T2>());

  using M1 = typename IsBasicDataType<T1>::ElementaryType;
  using M2 = typename IsBasicDataType<T2>::ElementaryType;

  //int mode_copy = 0;
  //int Ndata1    = 1;
  //int Ndata2    = 1;
  Int Ndata1 = sizeof(T1) / sizeof(M1);
  Int Ndata2 = sizeof(T2) / sizeof(M2);
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
}

template <class T1, class T2 >
void copy_fields(qlat::Field<T1>& pr, const qlat::Field<T2>& p0)
{
  Qassert(p0.initialized);
  if(!pr.initialized){
    pr.init(p0.geo(), p0.multiplicity);
  }

  T1* res = (T1*) qlat::get_data(pr).data();
  const T2* src = (T2*) qlat::get_data(p0).data();
  Qassert(pr.multiplicity == p0.multiplicity);
  const Int civ = pr.multiplicity;
  copy_fields<T1, T2>(res, src, civ, pr.geo(), p0.geo());
}

template <class T1, class T2 >
void copy_fields(std::vector<qlat::Field<T1> >& pr, const std::vector< qlat::Field<T2> >& p0)
{
  if(pr.size() != p0.size()){
    pr.resize(0);
    pr.resize(p0.size());
  }
  for(unsigned int iv=0;iv < p0.size();iv++){
    copy_fields(pr[iv], p0[iv]);
  }
}

template <class T1, class T2 >
void copy_fieldsG(qlat::FieldG<T1>& pr, const qlat::FieldG<T2>& p0)
{
  copy_fields(pr, p0);
  pr.mem_order = p0.mem_order;
}

template <class T1, class T2 >
void copy_fieldsG(std::vector< qlat::FieldG<T1> >& pr, const std::vector< qlat::FieldG<T2> >& p0)
{
  if(pr.size() != p0.size()){
    pr.resize(0);
    pr.resize(p0.size());
  }
  for(unsigned int iv=0;iv < p0.size();iv++){
    if(!pr[iv].initialized){pr[iv].init(p0[iv]);}
    copy_fieldsG(pr[iv], p0[iv]);
  }
}

template <class T1, class T2 >
void copy_fieldsG(SelectedFieldG<T1>& pr, const SelectedFieldG<T2>& p0)
{
  Qassert(p0.initialized);
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
  Qassert(fs.initialized);
  //Ty --> double float ...
  Qassert(GetBasicDataType<Fieldy>::get_type_name() != std::string("unknown_type"));
  using D = typename GetBasicDataType<Fieldy>::ElementaryType;
  Qassert(IsBasicTypeReal<D>());

  const Long Nd = GetFieldSize(fs);
  Qassert(Nd % (2 * sizeof(D)) == 0);
  const Long Ndata = Nd / ( 2 * sizeof(D));
  Qassert(Ndata % block == 0);
  qlat::ComplexT<D >  res = 0.0;
  void* src = (void*) qlat::get_data(fs).data();
  res = vec_norm2((qlat::ComplexT<D >*) src, (qlat::ComplexT<D >*)src, Ndata, QMGPU, block);
  if(log){
    qmessage("gauge norm %.15e %.15e \n",  res.real(), res.imag());
  }
  return res.real();
}

template <class T1, class T2, class T3 >
void fields_operationsGT(std::vector<FieldG<T1> >& pr, std::vector<FieldG<T2> >& p0, std::vector<FieldG<T3> >& p1, 
  const T1 fr = T1(1.0, 0.0), const T1 f0 = T1(1.0, 0.0), const T1 f1 = T1(1.0, 0.0))
{
  TIMER_FLOPS("fields_operationsG");
  Qassert(IsTypeComplex<T1>() and IsTypeComplex<T1>() and IsTypeComplex<T1>());
  const Long nvec = p1.size();
  Qassert(nvec > 0);Qassert(p0.size() == p1.size() and p0.size() == pr.size());
  const Int civ = p0[0].multiplicity;
  const QMEM_ORDER mem_order = p0[0].mem_order;
  for(Long iv=0;iv<nvec;iv++){
    Qassert(p0[iv].initialized and p1[iv].initialized and pr[iv].initialized);
    Qassert(civ == p0[iv].multiplicity and civ == p1[iv].multiplicity and civ == pr[iv].multiplicity);
    if(civ != 1){
      Qassert(mem_order == p0[iv].mem_order and mem_order == p1[iv].mem_order and mem_order == pr[iv].mem_order);
    }
  }
  Qassert(mem_order == QLAT_DEFAULT or mem_order == QLAT_OUTTER);

  const Geometry& geo = p0[0].geo();
  const Long V = geo.local_volume();
  long long Tfloat = 3 * V * civ * sizeof(T1);// only count the memory bandwith
  timer.flops += Tfloat;

  const Geometry& geor = pr[0].geo(); 
  const Geometry& geo0 = p0[0].geo(); 
  const Geometry& geo1 = p1[0].geo(); 
  vector<T1* > Pr;
  vector<T2* > P0;
  vector<T3* > P1;
  Pr.resize(nvec);
  P0.resize(nvec);
  P1.resize(nvec);
  for(Long iv=0;iv<nvec;iv++){
    Pr[iv] = (T1*) get_data(pr[iv]).data();
    P0[iv] = (T2*) get_data(p0[iv]).data();
    P1[iv] = (T3*) get_data(p1[iv]).data();
  }
  const Long Vr = geor.local_volume_expanded();
  const Long V0 = geo0.local_volume_expanded();
  const Long V1 = geo1.local_volume_expanded();

  // geometry can be extended, but only operated on local numbers
  qacc_for(isp, V, {
    const Coordinate xl = geo.coordinate_from_index(isp);
    const Long ir = geor.offset_from_coordinate(xl, 1);
    const Long i0 = geo0.offset_from_coordinate(xl, 1);
    const Long i1 = geo1.offset_from_coordinate(xl, 1);
    if(mem_order == QLAT_DEFAULT)
    {
      for(Int i=0;i<nvec;i++)
      for(Int j=0;j<civ;j++)
      {
        Pr[i][ir*civ + j] = fr * Pr[i][ir*civ + j] + f0 * P0[i][i0*civ + j] + f1 * P1[i][i1*civ + j];
      }
    }
    if(mem_order == QLAT_OUTTER){
      for(Int i=0;i<nvec;i++)
      for(Int j=0;j<civ;j++)
      {
        Pr[i][j*Vr + ir] = fr * Pr[i][j*Vr + ir] + f0 * P0[i][j * V0 + i0] + f1 * P1[i][j * V1 + i1];
      }
    }
  });
}

//template <class Fd1, class Fd2, class Fd3, class Ty >
//void fields_operations(Fd1& pr, Fd2& p0, Fd3& p1, 
//  const Ty f0 = Ty(1.0, 0.0), const Ty f1 = Ty(1.0, 0.0), const Ty f2 = Ty(1.0, 0.0))
//{
//  TIMER("fields_operations");
//  Qassert(p0.initialized);
//  Qassert(p1.initialized);
//
//  const Geometry& geo = p0.geo();
//  const Int civ = p0.multiplicity;
//  if(!pr.initialized){pr.init(geo, civ);}
//  Qassert(civ == pr.multiplicity and civ == p1.multiplicity);
//
//  Qassert(GetBasicDataType<Fd1>::get_type_name() != std::string("unknown_type"));
//  Qassert(GetBasicDataType<Fd2>::get_type_name() != std::string("unknown_type"));
//  Qassert(GetBasicDataType<Fd3>::get_type_name() != std::string("unknown_type"));
//  using D1 = typename GetBasicDataType<Fd1>::ElementaryType; 
//  using D2 = typename GetBasicDataType<Fd2>::ElementaryType; 
//  using D3 = typename GetBasicDataType<Fd3>::ElementaryType; 
//  Qassert(IsBasicTypeReal<D1>() and IsBasicTypeReal<D2>() and IsBasicTypeReal<D3>());
//
//  // geometry can be extended, but only operated on local numbers
//  qacc_for(isp, geo.local_volume(), {
//    const Coordinate xl = geo.coordinate_from_index(isp);
//    ComplexT<D1>* r0 = (ComplexT<D1>*) pr.get_elems(xl).p;
//    ComplexT<D2>* s0 = (ComplexT<D2>*) p0.get_elems(xl).p;
//    ComplexT<D3>* s1 = (ComplexT<D2>*) p1.get_elems(xl).p;
//    for (Int m = 0; m < civ; ++m) {
//      r0[m] = f0 * r0[m] + f1 * s0[m] + f2 * s1[m];
//    }
//  });
//}

template <class T1, class T2, class T3 >
void fields_operations(std::vector<FieldG<T1> >& pr, std::vector<FieldG<T2> >& p0, std::vector<FieldG<T3> >& p1, const T1 f0 = T1(1.0, 0.0), const T1 f1 = T1(1.0, 0.0), const T1 f2 = T1(1.0, 0.0))
{
  fields_operationsGT(pr, p0, p1, f0, f1, f2);
}

template <class T1, class T2, class T3, Int civ >
void fields_operations(std::vector<qlat::FieldM<T1, civ> >& pr, std::vector<qlat::FieldM<T2, civ> >& p0, std::vector<qlat::FieldM<T3, civ> >& p1, const T1 f0 = T1(1.0, 0.0), const T1 f1 = T1(1.0, 0.0), const T1 f2 = T1(1.0, 0.0) )
{
  const Int Nvec = pr.size();
  std::vector<FieldG<T1>> prG;
  std::vector<FieldG<T2>> p1G;
  std::vector<FieldG<T3>> p0G;
  prG.resize(Nvec);
  p1G.resize(Nvec);
  p0G.resize(Nvec);
  for(Int iv=0;iv<Nvec;iv++){
    prG[iv].set_pointer(pr[iv], QLAT_DEFAULT);
    p1G[iv].set_pointer(p1[iv], QLAT_DEFAULT);
    p0G[iv].set_pointer(p0[iv], QLAT_DEFAULT);
  }
  fields_operationsGT(prG, p0G, p1G, f0, f1, f2);
}

template <class T1, class T2, class T3 >
void fields_operations(FieldG<T1>& pr, FieldG<T2>& p0, FieldG<T3>& p1, const T1 f0 = T1(1.0, 0.0), const T1 f1 = T1(1.0, 0.0), const T1 f2 = T1(1.0, 0.0))
{
  std::vector<FieldG<T1>> prG;prG.resize(1);
  std::vector<FieldG<T2>> p1G;p1G.resize(1);
  std::vector<FieldG<T3>> p0G;p0G.resize(1);
  prG[0].set_pointer(pr);
  p1G[0].set_pointer(p1);
  p0G[0].set_pointer(p0);
  fields_operationsGT(prG, p0G, p1G, f0, f1, f2);
}

template <class T1, class T2, class T3, Int civ >
void fields_operations(FieldM<T1, civ>& pr, FieldM<T2, civ>& p0, FieldM<T3, civ>& p1, const T1 f0 = T1(1.0, 0.0), const T1 f1 = T1(1.0, 0.0), const T1 f2 = T1(1.0, 0.0) )
{
  std::vector<FieldG<T1>> prG;prG.resize(1);
  std::vector<FieldG<T2>> p1G;p1G.resize(1);
  std::vector<FieldG<T3>> p0G;p0G.resize(1);
  prG[0].set_pointer(pr, QLAT_DEFAULT);
  p1G[0].set_pointer(p1, QLAT_DEFAULT);
  p0G[0].set_pointer(p0, QLAT_DEFAULT);
  fields_operationsGT(prG, p0G, p1G, f0, f1, f2);
}

//template <class Fd1, class Fd2, class Fd3, class Ty>
//void fields_operations(std::vector<Fd1 >& pr, std::vector<Fd2 >& p0, std::vector<Fd3 >& p1, const Ty f0 = Ty(1.0, 0.0), const Ty f1 = Ty(1.0, 0.0), const Ty f2 = Ty(1.0, 0.0))
//{
//  Qassert(p0.size() == p1.size());
//  if(p0.size() == 0){return ;}
//  //const Geometry& geo = p0[0].geo();
//  const Int Nvec = p0.size();
//  if(Long(pr.size()) != Nvec){pr.resize(Nvec);}
//  for(Int vi=0;vi<Nvec;vi++){
//    fields_operations(pr[vi], p0[vi], p1[vi], f0, f1, f2);
//  }
//}

//template <class T1, class T2, class T3, class Ty>
//void fields_operationsG(std::vector<qlat::FieldG<T1> >& pr, std::vector<qlat::FieldG<T2> >& p0, std::vector<qlat::FieldG<T3> >& p1, const Ty f0 = Ty(1.0, 0.0), const Ty f1 = Ty(1.0, 0.0), const Ty f2 = Ty(1.0, 0.0))
//{
//  Qassert(p0.size() == p1.size());
//  if(p0.size() == 0){return ;}
//  //const Geometry& geo = p0[0].geo();
//  const unsigned int Nvec = p0.size();
//  if(pr.size() != Nvec){pr.resize(Nvec);}
//  for(unsigned int vi=0;vi<Nvec;vi++){
//    fields_operations(pr[vi], p0[vi], p1[vi], f0, f1, f2);
//  }
//}

template <class T1, class T2, class T3, Int civ >
void fields_additions(std::vector<qlat::FieldM<T1, civ> >& pr, std::vector<qlat::FieldM<T2, civ> >& p0, std::vector<qlat::FieldM<T3, civ> >& p1)
{
  fields_operations(pr, p0, p1, T1( 0.0, 0.0), T1( 1.0, 0.0), T1( 1.0, 0.0));
}

template <class T1, class T2, class T3, Int civ >
void fields_subtractions(std::vector<qlat::FieldM<T1, civ> >& pr, std::vector<qlat::FieldM<T2, civ> >& p0, std::vector<qlat::FieldM<T3, civ> >& p1)
{
  fields_operations(pr, p0, p1, T1( 0.0, 0.0), T1( 1.0, 0.0), T1(-1.0, 0.0));
}

template <class T1, class T2, class T3>
void fields_additionsG(std::vector<FieldG<T1> >& pr, std::vector<FieldG<T2> >& p0, std::vector<FieldG<T3> >& p1)
{
  fields_operationsGT(pr, p0, p1, T1(0.0, 0.0), T1( 1.0, 0.0), T1( 1.0, 0.0));
}

template <class T1, class T2, class T3>
void fields_subtractionsG(std::vector<FieldG<T1> >& pr, std::vector<FieldG<T2> >& p0, std::vector<FieldG<T3> >& p1)
{
  fields_operationsGT(pr, p0, p1, T1(0.0, 0.0), T1( 1.0, 0.0), T1(-1.0, 0.0));
}

//template <class Fd1, class Fd2, class Fd3>
//void fields_additionsG(std::vector<Fd1 >& pr, std::vector<Fd2 >& p0, std::vector<Fd3 >& p1)
//{
//  using Ty = ComplexT<double>;
//  fields_operations(pr, p0, p1, Ty(0.0, 0.0), Ty( 1.0, 0.0), Ty( 1.0, 0.0));
//}
//
//template <class Fd1, class Fd2, class Fd3>
//void fields_subtractionsG(std::vector<Fd1 >& pr, std::vector<Fd2 >& p0, std::vector<Fd3 >& p1)
//{
//  using Ty = ComplexT<double>;
//  fields_operations(pr, p0, p1, Ty(0.0, 0.0), Ty( 1.0, 0.0), Ty(-1.0, 0.0));
//}

template <class T1, class T2, Int civ >
void fields_equal(std::vector<qlat::FieldM<T1, civ> >& pr, std::vector<qlat::FieldM<T2, civ> >& ps)
{
  fields_operations(pr, ps, ps, T1( 0.0, 0.0), T1( 1.0, 0.0), T1( 0.0, 0.0));
}

template <class T1 >
void fields_conj(T1** src, const Int nvec, const Int civ, const Geometry& geo)
{
  TIMER("fields_conj");
  Qassert(IsBasicDataType<T1>::value and IsBasicDataType<T1>::is_complex);
  Qassert(IsBasicTypeReal<T1>());// real numbers

  // M1 will be double / float
  using M1 = typename IsBasicDataType<T1>::ElementaryType;
  Int Ndata = sizeof(T1) / sizeof(ComplexT<M1 >);

  Ndata = Ndata * civ;// multiplicity

  for(Int iv=0;iv<nvec;iv++){
    ComplexT<M1 >* res = (ComplexT<M1 >*) src[iv];
    qacc_for(isp, geo.local_volume(), {
      ComplexT<M1 >* p = &res[isp*Ndata + 0];
      for(Int d=0;d<Ndata;d++)
      {
        p[d] = qlat::qconj(p[d]);
      }
    });
  }
}

template <class T1>
void fields_conj(qlat::FieldG<T1>& pr)
{
  Qassert(pr.initialized);
  qlat::vector<T1* > src;src.resize(1);
  src[0] = (T1*) qlat::get_data(pr).data();
  const Int civ = pr.multiplicity;
  fields_conj(src.data(), 1, civ, pr.geo());
}


template <class T1>
void fields_conj(qlat::Field<T1>& pr)
{
  Qassert(pr.initialized);
  qlat::vector<T1* > src;src.resize(1);
  src[0] = (T1*) qlat::get_data(pr).data();
  const Int civ = pr.multiplicity;
  fields_conj(src.data(), 1, civ, pr.geo());
}

template <class T1, Int civ >
void fields_conj(std::vector<qlat::FieldM<T1, civ> >& pr)
{
  const Int nvec = pr.size();
  if(nvec == 0){return ;}
  qlat::vector<T1* > src;src.resize(nvec);
  for(Int si=0;si<nvec;si++)
  {
    Qassert(pr[si].initialized);
    src[si] = (T1*) qlat::get_data(pr[si]).data();
  }
  fields_conj(src.data(), nvec, civ, pr[0].geo());
}

template <class M >
void switch_orders(M& f, const QMEM_ORDER order_){
  if(f.mem_order == order_){return ;}
  Qassert(order_ == QLAT_DEFAULT or order_== QLAT_OUTTER);
  Qassert(f.mem_order == QLAT_DEFAULT or f.mem_order == QLAT_OUTTER);
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

// TODO change it into extended fields
template <class T1, class T2, class Ty >
void fields_operations_localG(std::vector<qlat::FieldG<T1> >& pr, std::vector<qlat::FieldG<T2> >& ps, const Ty f0 = Ty(1.0, 0.0))
{
  TIMER("fields_operations_localG");
  Qassert(pr.size() == ps.size());
  const Long Nsrc = ps.size();
  vector<T1* > rP;rP.resize(Nsrc);
  vector<T2* > sP;sP.resize(Nsrc);
  for(Long i=0;i<Nsrc;i++){
    Qassert(pr[i].initialized and ps[i].initialized);
    Qassert(pr[i].multiplicity == ps[i].multiplicity);
    rP[i] = (T1*) get_data(pr[i]).data();
    sP[i] = (T2*) get_data(ps[i]).data();
  }

  const Geometry& geo = pr[0].geo();
  const Int civ = pr[0].multiplicity;
  // geometry can be extended, but only operated on local numbers
  qacc_for(isp, geo.local_volume(), {
    for(Long si=0;si<Nsrc;si++){
      for (Int m = 0; m < civ; ++m) {
        const Long idx = isp * civ + m;
        rP[si][idx] += f0 * sP[si][idx];
      }
    }
  });
}

}


#endif

