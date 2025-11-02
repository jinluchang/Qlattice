// utils_gauge_box.h
// Gen Wang
// Apr. 2024

#ifndef UTILS_GAUGE_BOX_H
#define UTILS_GAUGE_BOX_H
#pragma once

#include "general_funs.h"
#include "utils_gauge_field.h"
#include "utils_field_expand.h"

namespace qlat
{

//u9 usual, u6 reconsturct col 6, uA anti-hermition 6
enum struct Su3 : Int {
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
  Int Multiplicity;
  box<Geometry> geo;
  Handle<Ty > v;//pointers to object
  //
  FieldBoxT()
  {
    // TIMER("FieldBoxT::FieldBoxT()")
    Qassert(v.p == NULL);
    //type = Su3::u9;
    Multiplicity = 9;
  }

  template <class M, Int civ >
  void init(FieldM<M, civ>& fr)
  {
    Qassert(Is_data_double<Ty>() == Is_data_double<M>());
    Qassert(fr.multiplicity == civ);
    Qassert(sizeof(M) % sizeof(Ty) == 0);
    const Int Nd = sizeof(M) / sizeof(Ty);
    v.p  = (Ty*) get_data(fr).data();
    geo.set_view(fr.geo());
    Multiplicity = fr.multiplicity * Nd;

    if(type == Su3::u9){Qassert(Multiplicity % 9 == 0);}
    if(type == Su3::u6){Qassert(Multiplicity % 6 == 0);}
    if(type == Su3::uA){Qassert(Multiplicity % 6 == 0);}
  }
  //
  template <class M, Int civ >
  FieldBoxT(FieldM<M, civ>& fr)
  {
    init(fr);
  }
  //
  FieldBoxT(Ty* p, const Geometry& geo_, int8_t Multiplicity_ = 1)
  {
    init(p, geo_, Multiplicity_);
  }
  // will check whether geo is good to use
  void init(Ty* p, const Geometry& geo_, int8_t Multiplicity_ = 1)
  {
    const MemType gmem = check_mem_type( &geo_);
    const MemType pmem  = check_mem_type(&p);
    if(gmem != MemType::Uvm){
      Qassert(gmem == pmem);
    }
    Multiplicity = Multiplicity_;
    v.p  = p;
    geo.set_view(geo_);
    //
    if(type == Su3::u9){Qassert(Multiplicity % 9 == 0);}
    if(type == Su3::u6){Qassert(Multiplicity % 6 == 0);}
    if(type == Su3::uA){Qassert(Multiplicity % 6 == 0);}
  }

  FieldBoxT(const FieldBoxT<Ty, type>& vp)
  {
    v    = vp.v;
    geo.set_view(vp.geo());
    Multiplicity = vp.Multiplicity;
  }

  FieldBoxT(FieldBoxT<Ty, type>&& vp) noexcept
  {
    v    = vp.v;
    geo.set_view(vp.geo());
    Multiplicity = vp.Multiplicity;
  }
  ~FieldBoxT()
  {
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
  //
  qacc Int initialized()
  {
    if(v.p == NULL){return 0;}
    else{return 1;}
  }
  //
  qacc void swap(FieldBoxT<Ty, type>& x)
  {
    Handle<Ty > t = v;
    v = x.v;
    x.v = t;

    //box<Geometry> g =   geo;
    //geo   = x.geo;
    //x.geo = g;

    box<Geometry> g;

    g.set_view(geo);
    geo.set_view(x.geo);
    x.geo.set_view(g);

    int8_t m = x.Multiplicity;
    Multiplicity = x.Multiplicity;
    x.Multiplicity = m;
  }
  //
  //
  FieldBoxT<Ty, type >& operator=(const FieldBoxT<Ty, type >& vp)
  {
    v    = vp.v;
    geo.set_view(vp.geo());
    Multiplicity = vp.Multiplicity;
    return *this;
  }
  FieldBoxT<Ty, type>& operator=(FieldBoxT<Ty, type>&& vp) noexcept
  {
    v    = vp.v;
    geo.set_view(vp.geo());
    Multiplicity = vp.Multiplicity;
    return *this;
  }
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
  qacc Ty* data() { return v.p; }
  qacc const Ty* data() const { return v.p; }

  qacc Int nc(){
    Int Nc = 0;
    if(type == Su3::u9){Nc = Multiplicity / 9;}
    if(type == Su3::u6){Nc = Multiplicity / 6;}
    if(type == Su3::uA){Nc = Multiplicity / 6;}
    return Nc;
  }

  // print out info
  void show(){
    const Long V    = geo().local_volume();
    const Long V_ex = geo().local_volume_expanded();
    qmessage("FieldM %2d, nc %2d, V %10ld, Vex %10ld, ", Multiplicity, nc(), V, V_ex);
    display_mem_type(v.p);
  }

  template <class Tb>
  qacc void construct(Tb* r, const Coordinate& xl, const Long mu = 0)
  {
    const Ty* s = data(xl, mu);
    if(type == Su3::u9)
    {
      for(Int i=0;i<9;i++){r[i] = s[i];}
      //su3_unitarize(r);
      return ;
    }
    //row col may be exchanged
    if(type == Su3::u6)
    {
      //for(Int i=0;i<6;i++){r[i] = s[i];}
      for(Int i=0;i<3;i++){r[i*3 + 0] = s[0*3+i];}
      for(Int i=0;i<3;i++){r[i*3 + 1] = s[1*3+i];}
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
      for(Int i=0;i<9;i++){r[i] = s[i];}
      return ;
    }
    if(type == Su3::u6)
    {
      //for(Int i=0;i<6;i++){r[i] = s[i];}
      //qlat::su3_reconstruct_col(r);
      for(Int i=0;i<3;i++){r[0*3+i] = s[i*3 + 0];}
      for(Int i=0;i<3;i++){r[1*3+i] = s[i*3 + 1];}
      return ;
    }
    if(type == Su3::uA)
    {
      Ty b[9];
      for(Int i=0;i<9;i++){b[i] = s[i];}
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
      for(Int i=0;i<9;i++){r[i] += s[i];}
      return ;
    }
    if(type == Su3::u6)
    {
      //for(Int i=0;i<6;i++){r[i] = s[i];}
      //qlat::su3_reconstruct_col(r);
      for(Int i=0;i<3;i++){r[0*3+i] += s[i*3 + 0];}
      for(Int i=0;i<3;i++){r[1*3+i] += s[i*3 + 1];}
      return ;
    }
    if(type == Su3::uA)
    {
      Ty b[9];
      for(Int i=0;i<9;i++){b[i] = s[i];}
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
  const Int Nr = fr.nc();
  const Int Ns = fs.nc();
  Qassert(Nr > 0 and Ns > 0 and Nr == Ns);
  //
  //fr.show();
  //fs.show();
  qacc_for(isp, V, {
    Ty b0[9];
    const Coordinate xl = geo.coordinate_from_index(isp);
    for(Int mu=0;mu<Nr;mu++)
    {
      fs.construct(b0, xl, mu);
      fr.copy_from(b0, xl, mu);
    }
  });
}

template <class Ty, Su3 ta, class Tf, Su3 tb>
void copy_uf(FieldBoxT<Ty, ta >& fr, FieldBoxT<Tf, tb >& fs, Int dr, Int ds)
{
  //Qassert(fs.geo() == fr.geo());
  Qassert(fs.initialized() and fr.initialized())
  const Geometry& geo = fs.geo();
  const Long V = geo.local_volume();
  const Int Nr = fr.nc();
  const Int Ns = fs.nc();
  Qassert(dr < Nr and ds < Ns and Nr > 0 and Ns > 0);
  //
  qacc_for(isp, V, {
    Ty b0[9];
    const Coordinate xl = geo.coordinate_from_index(isp);
    fs.construct(b0, xl, ds);
    fr.copy_from(b0, xl, dr);
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

template <class Ty, Su3 type>
void refresh_expanded_GPU(FieldBoxT<Ty, type>& g, const std::string& tag = "", const Int GPU = 1, const QBOOL dummy = QTRUE)
{
  refresh_expanded_GPU(g.data(), g.geo(), g.Multiplicity, tag, GPU, dummy );
}

template <class Ty, Su3 type>
void refresh_expanded_GPU(FieldBoxT<Ty, type>& g, const Int dir, const Int GPU = 1, const QBOOL dummy = QTRUE)
{
  refresh_expanded_GPU(g.data(), g.geo(), g.Multiplicity, dir, GPU, dummy );
}

}

#endif
