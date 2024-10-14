// utils_gauge_box.h
// Gen Wang
// Apr. 2024

#ifndef UTILS_GAUGE_BOX_H
#define UTILS_GAUGE_BOX_H
#pragma once

#include "general_funs.h"
#include "utils_gaugefield.h"

namespace qlat
{

//u9 usual, u6 reconsturct col 6, uA anti-hermition 6
enum struct Su3{
  u9,
  u6,
  uA
};

//referece from some pointers
template <class Ty, Su3 type>
struct API FieldBoxT {
  // Avoid copy constructor when possible
  // (it is likely not what you think it is)
  // Only used in qacc macros, or if it is already a copy.
  //
  //Su3 type;   // 6, anti-hermition, 9 usual
  int Multiplicity;
  box_acc<Geometry> geo;
  Handle<Ty > v;//pointers to object
  //
  FieldBoxT()
  {
    // TIMER("FieldBoxT::FieldBoxT()")
    qassert(v.p == NULL);
    //type = Su3::u9;
    Multiplicity = 1;
  }

  template <class M, int civ >
  void init(qlat::FieldM<M, civ>& fr)
  {
    // TIMER("FieldBoxT::FieldBoxT(&)")
    //Qassert(get_data_type_is_double<Ty>() == get_data_type_is_double<M>());
    Qassert(Is_data_double<Ty>() == Is_data_double<M>());
    Qassert(fr.multiplicity == civ);
    Qassert(sizeof(M) % sizeof(Ty) == 0);
    const int Nd = sizeof(M) / sizeof(Ty);
    v.p  = (Ty*) qlat::get_data(fr).data();
    geo.set(fr.geo());
    Multiplicity = fr.multiplicity * Nd;

    if(type == Su3::u9){Qassert(Multiplicity % 9 == 0);}
    if(type == Su3::u6){Qassert(Multiplicity % 6 == 0);}
    if(type == Su3::uA){Qassert(Multiplicity % 6 == 0);}
    //geo().multiplicity = 1;
    //print0("multi %5d, geoM %5d \n", Multiplicity, geo().multiplicity);
  }

  template <class M, int civ >
  FieldBoxT(qlat::FieldM<M, civ>& fr)
  {
    init(fr);
  }

  FieldBoxT(Ty* p, const Geometry geo_, char Multiplicity_ = 1)
  {
    // TIMER("FieldBoxT::FieldBoxT(&)")
    init(p, geo_, Multiplicity_);
  }

  void init(Ty* p, const Geometry geo_, char Multiplicity_ = 1)
  {
    Multiplicity = Multiplicity_;
    v.p  = p;
    geo.set(geo_);
    //geo().multiplicity = 1;

    if(type == Su3::u9){Qassert(Multiplicity % 9 == 0);}
    if(type == Su3::u6){Qassert(Multiplicity % 6 == 0);}
    if(type == Su3::uA){Qassert(Multiplicity % 6 == 0);}
  }

  FieldBoxT(const FieldBoxT<Ty, type>& vp)
  {
    // TIMER("FieldBoxT::FieldBoxT(&)")
    //qassert(x.type == type);
    v    = vp.v;
    geo.set(vp.geo());
    Multiplicity = vp.Multiplicity;
  }

  FieldBoxT(FieldBoxT<Ty, type>&& vp) noexcept
  {
    // TIMER("FieldBoxT::FieldBoxT(&&)")
    //type = vp.type;
    //qassert(x.type == type);
    v    = vp.v;
    geo.set(vp.geo());
    Multiplicity = vp.Multiplicity;
  }
  //FieldBoxT(const FieldBoxT<Ty>& vp)
  //{
  //  // TIMER("FieldBoxT::FieldBoxT(std::FieldBoxT&)")
  //  *this = vp;
  //}
  //
  ~FieldBoxT()
  {
    // TIMER("FieldBoxT::~FieldBoxT()")
    clear();
  }
  //
  void init()
  {
    clear();
  }
  //
  void clear()
  {
    v.p = NULL;
  }

  qacc int initialized()
  {
    if(v.p == NULL){return 0;}
    else{return 1;}
  }
  //
  qacc void swap(FieldBoxT<Ty, type>& x)
  {
    //qassert(x.type == type);
    //Su3 tb = x.type;
    //x.type = type;
    //type = tb;

    Handle<Ty > t = v;
    v = x.v;
    x.v = t;

    box_acc<Geometry> g =   geo;
    geo   = x.geo;
    x.geo = g;

    char m = x.Multiplicity;
    Multiplicity = x.Multiplicity;
    x.Multiplicity = m;
  }
  //
  //
  FieldBoxT<Ty, type >& operator=(const FieldBoxT<Ty, type >& vp)
  {
    // TIMER("FieldBoxT::operator=(&)");
    //clear();
    //type = vp.type;
    //qassert(x.type == type);
    v    = vp.v;
    geo.set(vp.geo());
    Multiplicity = vp.Multiplicity;
    return *this;
  }
  FieldBoxT<Ty, type>& operator=(FieldBoxT<Ty, type>&& vp) noexcept
  {
    // TIMER("FieldBoxT::operator=(&&)");
    //type = vp.type;
    //qassert(x.type == type);
    v    = vp.v;
    geo.set(vp.geo());
    Multiplicity = vp.Multiplicity;
    return *this;
  }
  //FieldBoxT<Ty>& operator=(const std::FieldBoxT<Ty>& vp)
  //{
  //  // TIMER("FieldBoxT::operator=(std::FieldBoxT&)");
  //  //clear();
  //  type = vp.type;
  //  v    = vp.v;
  //  geo.set(vp.geo());
  //  Multiplicity = vp.Multiplicity;
  //  //resize(vp.size());
  //  //for (Long i = 0; i < v.n; ++i) {
  //  //  v[i] = vp[i];
  //  //}
  //  return *this;
  //}

  qacc Ty* data(const Coordinate& xl, const Long mu = 0)
  {
    const Long index = geo().offset_from_coordinate(xl, 1);
    if(type == Su3::u9)
    {
      qassert(mu < Multiplicity/9);
      return &v.p[index * Multiplicity + mu*9 + 0];
    }
    if(type == Su3::u6 or type == Su3::uA)
    {
      qassert(mu < Multiplicity/6);
      return &v.p[index * Multiplicity + mu*6 + 0];
    }
    qassert(false);
    return NULL;
  }
  //
  //qacc const M& operator[](const Long i) const { return v[i]; }
  //qacc M& operator[](const Long i) { return v[i]; }
  ////
  //qacc Long size() const { return v.size(); }
  ////
  qacc Ty* data() { return v.p; }
  qacc const Ty* data() const { return v.p; }

  template <class Tb>
  qacc void construct(Tb* r, const Coordinate& xl, const Long mu = 0)
  {
    const Ty* s = data(xl, mu);
    if(type == Su3::u9)
    {
      for(int i=0;i<9;i++){r[i] = s[i];}
      //su3_unitarize(r);
      return ;
    }
    //row col may be exchanged
    if(type == Su3::u6)
    {
      //for(int i=0;i<6;i++){r[i] = s[i];}
      for(int i=0;i<3;i++){r[i*3 + 0] = s[0*3+i];}
      for(int i=0;i<3;i++){r[i*3 + 1] = s[1*3+i];}
      qlat::su3_reconstruct_col(r);
      //su3_unitarize(r);
      return ;
    }
    if(type == Su3::uA)
    {
      r[1*3 + 0] = s[0];
      r[2*3 + 0] = s[1];
      r[2*3 + 1] = s[2];

      r[0*3 + 1] = Tb(-1.0, 0.0) * qconj(s[0]);
      r[0*3 + 2] = Tb(-1.0, 0.0) * qconj(s[1]);
      r[1*3 + 2] = Tb(-1.0, 0.0) * qconj(s[2]);

      r[0*3 + 0] = s[3];
      r[1*3 + 1] = s[4];
      r[2*3 + 2] = s[5];
      return ;
    }
  }

  template <class Tb>
  qacc void copy_from(const Tb* s, const Coordinate& xl, const Long mu = 0)
  {
    Ty* r = data(xl, mu);
    if(type == Su3::u9)
    {
      for(int i=0;i<9;i++){r[i] = s[i];}
      return ;
    }
    if(type == Su3::u6)
    {
      //for(int i=0;i<6;i++){r[i] = s[i];}
      //qlat::su3_reconstruct_col(r);
      for(int i=0;i<3;i++){r[0*3+i] = s[i*3 + 0];}
      for(int i=0;i<3;i++){r[1*3+i] = s[i*3 + 1];}
      return ;
    }
    if(type == Su3::uA)
    {
      QLAT_ALIGN(QLAT_ALIGNED_BYTES) Ty b[9];
      for(int i=0;i<9;i++){b[i] = s[i];}
      //su3_traceless_anti_hermition(b);
      r[0] = b[1*3 + 0];
      r[1] = b[2*3 + 0];
      r[2] = b[2*3 + 1];

      r[3] = b[0*3 + 0];
      r[4] = b[1*3 + 1];
      r[5] = b[2*3 + 2];
      return ;
    }
  }

  template <class Tb>
  qacc void add_from(const Tb* s, const Coordinate& xl, const Long mu = 0)
  {
    Ty* r = data(xl, mu);
    if(type == Su3::u9)
    {
      for(int i=0;i<9;i++){r[i] += s[i];}
      return ;
    }
    if(type == Su3::u6)
    {
      //for(int i=0;i<6;i++){r[i] = s[i];}
      //qlat::su3_reconstruct_col(r);
      for(int i=0;i<3;i++){r[0*3+i] += s[i*3 + 0];}
      for(int i=0;i<3;i++){r[1*3+i] += s[i*3 + 1];}
      return ;
    }
    if(type == Su3::uA)
    {
      QLAT_ALIGN(QLAT_ALIGNED_BYTES) Ty b[9];
      for(int i=0;i<9;i++){b[i] = s[i];}
      //su3_traceless_anti_hermition(b);
      r[0] += b[1*3 + 0];
      r[1] += b[2*3 + 0];
      r[2] += b[2*3 + 1];

      r[3] += b[0*3 + 0];
      r[4] += b[1*3 + 1];
      r[5] += b[2*3 + 2];
      return ;
    }
  }

  void set_zero()
  {
    zero_Ty(data(), geo().local_volume_expanded()*Multiplicity, 1);
  }

};

template <class Ty, Su3 ta, class Tf, Su3 tb>
void copy_ux(FieldBoxT<Ty, ta >& fr, FieldBoxT<Tf, tb >& fs)
{
  //Qassert(fs.geo() == fr.geo());
  Qassert(fs.initialized() and fr.initialized())
  const Geometry& geo = fs.geo();
  const Long V = geo.local_volume();

  int Nr = 4;
  int Nb = 4;
  if(ta == Su3::u9){Nr = fr.Multiplicity / 9;}
  if(ta == Su3::u6){Nr = fr.Multiplicity / 6;}
  if(ta == Su3::uA){Nr = fr.Multiplicity / 6;}

  if(tb == Su3::u9){Nb = fs.Multiplicity / 9;}
  if(tb == Su3::u6){Nb = fs.Multiplicity / 6;}
  if(tb == Su3::uA){Nb = fs.Multiplicity / 6;}
  Qassert(Nr == Nb);

  qacc_for(isp, V, {
    const Coordinate xl = geo.coordinate_from_index(isp);
    QLAT_ALIGN(QLAT_ALIGNED_BYTES) Ty b0[9];
    for(int mu=0;mu<Nr;mu++)
    {
      fs.construct(b0, xl, mu);
      fr.copy_from(b0, xl, mu);
    }
  });
}

template <class Ty, class Tf, Su3 ta>
void copy_fields(FieldBoxT<Ty, ta >& fr, FieldBoxT<Tf, ta >& fs)
{
  if(fr.data() != fs.data())
  {
    copy_fields<Ty, Tf>(fr.data(), fs.data(), fs.Multiplicity, fr.geo(), fs.geo());
  }
}


//template <class Ty, Su3 ta, class Tf, Su3 tb>
//void copy_u9_from_u6(FieldBoxT<Ty, Su3::u9 >& f9, FieldBoxT<Ty, Su3::u6 >& f6)
//{
//  copy_u6_from_u9(f6, f9, 0);
//}

//template <class Ty, class Td>
//void copy_uA_from_u9(FieldBoxT<Ty, Su3::u6 >& fA, FieldBoxT<Ty, Su3::u9 >& f9, int dir = 1)
//{
//  const Geometry& geo = gf.geo();
//  const Long V = geo.local_volume();
//  fu.init(geo);
//  qacc_for(isp, V, {
//    const Coordinate xl = geo.coordinate_from_index(isp);
//    Ty* s2  = (Ty*) fu.get_elems(xl).p;
//    Ty sb[9];
//    Ty BUF[9];
//    for(int mu=0;mu<4;mu++)
//    {
//      qlat::ComplexT<Td >* s1  = (qlat::ComplexT<Td >*) gf.get_elem(xl, mu).p;
//      for(int i=0;i<9;i++){sb[i] = s1[i];}
//      su3_traceless_anti_hermition(sb);
//      s2[mu*6 + 0] = sb[1*3 + 0]; 
//      s2[mu*6 + 1] = sb[2*3 + 0]; 
//      s2[mu*6 + 2] = sb[2*3 + 1]; 
//
//      s2[mu*6 + 3] = sb[0*3 + 0]; 
//      s2[mu*6 + 4] = sb[1*3 + 1]; 
//      s2[mu*6 + 5] = sb[2*3 + 2]; 
//    }
//  });
//}
//
//template <class Ty, class Td>
//void copy_u9_from_uA(GaugeFieldT<Td >& gf, FieldM<Ty, 24>& fu)
//{
//  const Geometry& geo = fu.geo();
//  const Long V = geo.local_volume();
//  gf.init(geo);
//  qacc_for(isp, V, {
//    const Coordinate xl = geo.coordinate_from_index(isp);
//    Ty* sf  = (Ty*) fu.get_elems(xl).p;
//    for(int mu=0;mu<4;mu++)
//    {
//      Ty* s = &sf[mu*6 + 0];
//      qlat::ComplexT<Td >* r  = (qlat::ComplexT<Td >*) gf.get_elem(xl, mu).p;
//      r[1*3 + 0] = s[0];
//      r[2*3 + 0] = s[1];
//      r[2*3 + 1] = s[2];
//
//      r[0*3 + 1] = Ty(-1.0, 0.0) * qconj(s[0]);
//      r[0*3 + 2] = Ty(-1.0, 0.0) * qconj(s[1]);
//      r[1*3 + 2] = Ty(-1.0, 0.0) * qconj(s[2]);
//
//      r[0*3 + 0] = s[3];
//      r[1*3 + 1] = s[4];
//      r[2*3 + 2] = s[5];
//    }
//  });
//}


}

#endif
