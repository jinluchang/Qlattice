// utils_GAUGE_FIELD.h
// Gen Wang
// Mar. 2022

#ifndef UTILS_GAUGE_FIELD_H
#define UTILS_GAUGE_FIELD_H
#pragma once

#include "general_funs.h"
#include <qlat/qcd-topology.h>

namespace qlat
{

/////E is 3x3 array
template <typename Ty >
qacc void su3_one(Ty* E)
{
  for(Int i=0;i<9;i++){ E[i] = 0.0; }
  E[0*3 + 0] = 1.0;
  E[1*3 + 1] = 1.0;
  E[2*3 + 2] = 1.0;
}

qacc Long su3_6(const Geometry& geo, const Coordinate& xl)
{
  const Int nD = 6;
  return geo.offset_from_coordinate(xl, 1)*nD + 0;
}

qacc Long su3_n(const Geometry& geo, const Coordinate& xl)
{
  const Int nD = 9;
  return geo.offset_from_coordinate(xl, 1)*nD + 0;
}

qacc Long su3_6(const Geometry& geo, const Coordinate& xl, const Int mu)
{
  const Int Dim = 4;
  const Int nD  = 6;
  return (geo.offset_from_coordinate(xl, 1) * Dim + mu)*nD + 0;
  //return geo.offset_from_coordinate(xl, 1)*nD + 0;
}

qacc Long su3_n(const Geometry& geo, const Coordinate& xl, const Int mu)
{
  const Int Dim = 4;
  const Int nD  = 9;
  qassert(mu >= 0 and mu < Dim);
  return (geo.offset_from_coordinate(xl, 1) * Dim + mu)*nD + 0;
  //if(mu > 0){
  //  return (geo.offset_from_coordinate(xl, 1) * Dim + mu)*nD + 0;
  //}else{
  //  return (geo.offset_from_coordinate(xl, 1)           )*nD + 0;
  //}
}

template <typename Ty, bool clear, bool dag1, bool dag2>
qacc void su3_multi_kerB(Ty* res, Ty* Q1, Ty* Q2, Ty* BUF)
{
  for(Int a=0;a<3;a++)
  for(Int b=0;b<3;b++)
  {
    Ty& buf = BUF[b*3 + a];buf = 0.0;
    for(Int d=0;d<3;d++)
    {
      if(dag1 == false and dag2 == false){buf += Q1[b*3 + d]*Q2[d*3 + a];}
      if(dag1 == false and dag2 == true ){buf += Q1[b*3 + d]*qlat::qconj(Q2[a*3 + d]);}
      if(dag1 == true  and dag2 == false){buf += qlat::qconj(Q1[d*3 + b])*Q2[d*3 + a];}
      if(dag1 == true  and dag2 == true ){buf += qlat::qconj(Q1[d*3 + b])*qlat::qconj(Q2[a*3 + d]);}
    }
  }
  if( clear){for(Int a=0;a<9;a++){res[a]   = BUF[a];}}
  if(!clear){for(Int a=0;a<9;a++){res[a]  += BUF[a];}}
}

/////3x3 inner 3 is faster and 1,2,3 (row 0) --> 4,5,6 (row 1) --> 7,8,9 (row 2)
template <typename Ty, bool clear, bool dag1, bool dag2>
qacc void su3_multi_ker(Ty* res, Ty* Q1, Ty* Q2)
{
  Ty BUF[9];
  su3_multi_kerB<Ty, clear, dag1, dag2>(res, Q1, Q2, BUF);
}

template <typename Ty, bool clear>
qacc void su3_multi(Ty* res, Ty* Q1, Ty* Q2)
{
  su3_multi_ker<Ty, true, false, false>(res, Q1, Q2);
}

template <typename Ty>
qacc void su3_dagger(Ty* res, Ty* src)
{
  for(Int a=0;a<3;a++)
  for(Int b=0;b<3;b++)
  {
    res[a*3 + b] = qlat::qconj(src[b*3 + a]);
  }
}

template <class Ty, Long N>
qacc void normalize_array_c(Ty* src)
{
  Ty norm = Ty(0.0, 0.0);
  for(Long i=0;i<N;i++){norm += (qlat::qconj(src[i]) * src[i]).real();}
  norm = qsqrt(norm.real());
  for(Long i=0;i<N;i++){src[i] /= norm;}
}

template <class Ty, Long N>
qacc void orthogonalize_array_c(Ty* p2, Ty* p1)
{
  Ty c = Ty(0.0, 0.0);
  for (Long i = 0; i < N; i++) {
    c += qconj(p1[i]) * p2[i];
  }
  //qmessage("n %.8e %.8e, %.8e %.8e, %.8e %.8e \n", c.real(), c.imag(), p1[0].real(), p1[0].imag(), p2[0].real(), p2[0].imag());
  //if (not(c == 0.0)) 
  {
    for (Long i = 0; i < N; i++) {
      p2[i] -= c * p1[i];
      //p2[i] = p2[i] - c * p1[i];
    }
  }
}

template <class Ty>
qacc void cross_product_conj_c(Ty* v3, Ty* v1, Ty* v2)
// v3 = ( v1 x v2 )^*
{
  v3[0] = qconj(v1[1] * v2[2] - v1[2] * v2[1]);
  v3[1] = qconj(v1[2] * v2[0] - v1[0] * v2[2]);
  v3[2] = qconj(v1[0] * v2[1] - v1[1] * v2[0]);
}

template <typename Ty>
qacc void su3_unitarize_col(Ty* src)
{
  Ty p1[3],p2[3],p3[3];
  for(Int i=0;i<3;i++)
  {
    p1[i] = src[i*3 + 0];
    p2[i] = src[i*3 + 1];
    p3[i] = src[i*3 + 2];
  }

  normalize_array_c<Ty, 3>(p1);
  orthogonalize_array_c<Ty, 3>(p2, p1);
  normalize_array_c<Ty, 3>(p2);
  cross_product_conj_c(p3, p1, p2);
  //normalize_array_c<Ty, 3>(p3);

  for(Int i=0;i<3;i++)
  {
    src[i*3 + 0] = p1[i];
    src[i*3 + 1] = p2[i];
    src[i*3 + 2] = p3[i];
  }
}

template <typename Ty>
qacc void su3_unitarize(Ty* src)
{
  //Ty p1[3],p2[3],p3[3];
  //for(Int i=0;i<3;i++)
  //{
  //  p1[i] = src[i*3 + 0];
  //  p2[i] = src[i*3 + 1];
  //  p3[i] = src[i*3 + 2];
  //}

  Ty* p1 = &src[0 * 3 + 0];// may need to rotate index?
  Ty* p2 = &src[1 * 3 + 0];
  Ty* p3 = &src[2 * 3 + 0];

  normalize_array_c<Ty, 3>(p1);
  orthogonalize_array_c<Ty, 3>(p2, p1);
  normalize_array_c<Ty, 3>(p2);
  cross_product_conj_c(p3, p1, p2);
  //normalize_array_c<Ty, 3>(p3);

  //for(Int i=0;i<3;i++)
  //{
  //  src[i*3 + 0] = p1[i];
  //  src[i*3 + 1] = p2[i];
  //  src[i*3 + 2] = p3[i];
  //}

  //orthogonalize_array_c<Ty, 3>(p2, p1);
  //normalize_array_c<Ty, 3>(p2);
  //cross_product_conj_c(p3, p1, p2);
}

//only touch the diagonal part
template <typename Tb>
qacc void su3_traceless(ComplexT<Tb>* sz)
{
  Int use_type = 0;
  if(is_same<RealDD, Tb>()){use_type = 1;}
  ////==type 0
  if(use_type == 0)
  {
    Tb c_tr = Tb(0.0);
    for(Int a=0;a<3;a++){
      c_tr += sz[a*3 + a].imag();
    }
    Tb one3 = Tb(1.0) / Tb(3.0);
    c_tr = c_tr * one3;

    for(Int a=0;a<3;a++){
      Tb tmp = sz[a*3 + a].imag();
      tmp   -= c_tr;
      sz[a*3 + a] = ComplexT<Tb>(0.0, tmp);
    }
  }

  ////type 1, need to split Y, X when sum is near zero ......
  if(use_type == 1)
  {
    RealDD a[3];
    RealDD b[3];
    RealDD c0 = 0.0;
    RealDD c1 = 0.0;

    for(Int i=0;i<3;i++){
      a[i] = sz[i*3 + i].imag();

      b[i].Y() = a[i].X();
      b[i].X() = 0.0;

      a[i].X() = 0.0;
      c0 += a[i];
      c1 += b[i];
    }

    RealDD one3 = RealDD(1.0) / RealDD(3.0);

    c0  = c0 * one3;
    c1  = c1 * one3;
    c0 += c1;

    //if(double(qfabs(c0)) > 1e-50)
    for(Int a=0;a<3;a++){
      RealDD tmp = sz[a*3 + a].imag();
      tmp   -= c0;
      Tb f   = tmp;
      sz[a*3 + a] = ComplexT<Tb>(0.0, f);
    }
  }

  //RealDD one3;
  //one3.Y() = 0x1.5555555555555p-2;
  //one3.X() = 0x1.5555555555555p-56;
  //c_tr = c_tr / Tb(3.0);

}

template <typename Tb>
qacc void su3_anti_hermition(ComplexT<Tb>* sz)
{
  //for(Int i=0;i<9;i++){BUF[i] = sz[i] ;}
  for(Int a=0;a<3;a++)
  for(Int b=0;b<3;b++)
  {
    if(a > b){
      //Ty tmp = BUF[a*3 + b] - qconj(BUF[b*3 + a]);
      ComplexT<Tb> tmp = sz[a*3 + b] - qconj(sz[b*3 + a]);
      tmp = Tb(0.5) * tmp;
      sz[a*3 + b] =        tmp;
      sz[b*3 + a] = Tb(-1.0) * qconj(tmp);
    }
    // may not need when need traceless
    if(a == b){
      sz[a*3 + b] = ComplexT<Tb>(0.0, sz[a*3 + b].imag());
    }
  }
}

template <typename Tb>
qacc void su3_traceless_anti_hermition(ComplexT<Tb>* sz)
{
  su3_anti_hermition(sz);
  su3_traceless(sz);
}

template <typename Tf>
qacc void su3_reconstruct_row(qlat::ComplexT<Tf >* r)
{
  r[2*3 + 0] = qlat::qconj( r[0*3+1]*r[1*3+2] - r[0*3+2]*r[1*3+1] );
  r[2*3 + 1] = qlat::qconj( r[0*3+2]*r[1*3+0] - r[0*3+0]*r[1*3+2] );
  r[2*3 + 2] = qlat::qconj( r[0*3+0]*r[1*3+1] - r[0*3+1]*r[1*3+0] );
}

template <typename Tf , typename Td>
qacc void su3_reconstruct_row(qlat::ComplexT<Tf >* r, qlat::ComplexT<Td >* s)
{
  for(Int i=0;i<3;i++){r[0*3 + i] = s[0*3+i];}
  for(Int i=0;i<3;i++){r[1*3 + i] = s[1*3+i];}
  su3_reconstruct_row(r);
  //r[2*3 + 0] = qlat::qconj( r[0*3+1]*r[1*3+2] - r[0*3+2]*r[1*3+1] );
  //r[2*3 + 1] = qlat::qconj( r[0*3+2]*r[1*3+0] - r[0*3+0]*r[1*3+2] );
  //r[2*3 + 2] = qlat::qconj( r[0*3+0]*r[1*3+1] - r[0*3+1]*r[1*3+0] );
}

template <typename Tf>
qacc void su3_reconstruct_col(qlat::ComplexT<Tf >* r)
{
  r[0*3 + 2] = qlat::qconj( r[1*3+0]*r[2*3+1] - r[2*3+0]*r[1*3+1] );
  r[1*3 + 2] = qlat::qconj( r[2*3+0]*r[0*3+1] - r[0*3+0]*r[2*3+1] );
  r[2*3 + 2] = qlat::qconj( r[0*3+0]*r[1*3+1] - r[1*3+0]*r[0*3+1] );
}

template <typename Tf , typename Td>
qacc void su3_reconstruct_col(qlat::ComplexT<Tf >* r, qlat::ComplexT<Td >* s)
{
  for(Int i=0;i<3;i++){r[i*3 + 0] = s[0*3+i];}
  for(Int i=0;i<3;i++){r[i*3 + 1] = s[1*3+i];}
  su3_reconstruct_col(r);
}

template <class Td>
void set_rand_link(GaugeFieldT<Td> &gf, const Int seed = -1)
{
  if(seed == -1)
  {
    qacc_for(isp, gf.field.size(), { set_unit(gf.get_elem_offset(isp), 1.0);});
  }else{
    //T* res = (T*) gf.get_elem_offset(0).p;
    const Geometry& geo = gf.geo();
    qlat::ComplexT<Td>* res = (qlat::ComplexT<Td>*) qlat::get_data(gf).data();
    random_Ty(res, geo.local_volume()*gf.multiplicity*sizeof(ColorMatrixT<Td>)/(sizeof(Td)*2), 1, seed);

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

template <class Ta, class Td>
void copy_gf(GaugeFieldT<Ta> &g1, GaugeFieldT<Td> &g0)
{
  TIMER("copy_gf");
  const Geometry& geo = g0.geo();
  if(!g1.initialized){g1.init(geo);}

  ///cannot use this due to extended gauge fields
  //const long Ndata = geo.local_volume() * 4 * 9;
  //qlat::ComplexT<Td >* src = (qlat::ComplexT<Td >*) qlat::get_data(g0).data();
  //qlat::ComplexT<Td >* res = (qlat::ComplexT<Td >*) qlat::get_data(g1).data();
  //qlat::cpy_GPU(res, src, Ndata,   1, 1);

  qacc_for(isp, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(isp);
    qlat::ComplexT<Ta >* res = (qlat::ComplexT<Ta >*) g1.get_elem(xl, 0).p;
    qlat::ComplexT<Td >* src = (qlat::ComplexT<Td >*) g0.get_elem(xl, 0).p;
    for(Int i=0;i<9*4;i++){res[i] = src[i];}
    //for (Int dir = 0; dir < 4; ++dir)
    //{
    //  g1.get_elem(xl, dir) = g0.get_elem(xl, dir);
    //}
  });
}

template <class Td>
void Gauge_antihermition(GaugeFieldT<Td> &gf)
{
  const Geometry& geo = gf.geo();
  const Long V = geo.local_volume();
  qacc_for(isp, V, {
    const Coordinate xl = geo.coordinate_from_index(isp);
    for(Int mu=0;mu<4;mu++)
    {
      qlat::ComplexT<double >* s1  = (qlat::ComplexT<double >*) gf.get_elem(xl, mu).p;
      su3_traceless_anti_hermition(s1);
    }
  });
}

template <class Td>
void Gauge_reconstruct_col(GaugeFieldT<Td> &gf)
{
  const Geometry& geo = gf.geo();
  const Long V = geo.local_volume();
  qacc_for(isp, V, {
    const Coordinate xl = geo.coordinate_from_index(isp);
    for(Int mu=0;mu<4;mu++)
    {
      qlat::ComplexT<double >* s1  = (qlat::ComplexT<double >*) gf.get_elem(xl, mu).p;
      su3_reconstruct_col(s1);
    }
  });
}
}

#endif
