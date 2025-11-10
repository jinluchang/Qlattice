// utils_smear_vecs_2link.h
// Gen Wang
// Feb. 2023

#ifndef UTILS_SMEAR_VECS_2link_H
#define UTILS_SMEAR_VECS_2link_H
#pragma once

#include "general_funs.h"
#include <qlat/qcd-utils.h>
#include <qlat/qcd-prop.h>
#include <qlat/qcd-smear.h>
#include "utils_Vec_redistribute.h"
#include "utils_shift_vecs.h"
#include "utils_check_fun.h"
#include "utils_smear_vecs.h"
//#include "utils_construction.h"
 
namespace qlat{

template <class Td>
void smear_propagator_gwu_convension_shift(Propagator4dT<Td>& prop, const GaugeFieldT<Td >& gf, const double width, const Int step)
{
  const Geometry& geo = gf.geo();
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);

  shift_vec svec(fd, true);
  qlat::vector_gpu<qlat::ComplexT<Td > > gauge;
  {
    gauge.resize(2*4* gf.geo().local_volume() * 9);
    extend_links_to_vecs(gauge.data(), gf);
    svec.set_gauge(gauge.data(), 4,12);
  }
  Propagator4dT<Td> buf0;buf0.init(geo);
  Propagator4dT<Td> buf1;buf1.init(geo);
  Propagator4dT<Td> bufc;bufc.init(geo);
  Propagator4dT<Td> bufa;bufa.init(geo);
  qlat::ComplexT<Td > bw   = qlat::ComplexT<Td >( width*width/(4.0*step - 6.0*width*width), 0.0);
  qlat::ComplexT<Td > norm = qlat::ComplexT<Td >( (1 - 3.0*width*width/(2*step)), 0.0);
  //qlat::ComplexT<Td > b0 = qlat::ComplexT<Td >(1.0, 0.0) - norm;
  qlat::ComplexT<Td > b1 = norm * bw;
  ////qmessage("norm %.8e \n", norm.real());

  //qlat::vector_gpu<qlat::ComplexT<Td > > ps;
  //qlat::vector_gpu<qlat::ComplexT<Td > > p0;
  //qlat::vector_gpu<qlat::ComplexT<Td > > p1;
  //qlat::vector_gpu<qlat::ComplexT<Td > > pc;

  prop4D_factor(bufc, qlat::ComplexT<Td >(0.0, 0.0));
  bufc += prop;
  for (Int i = 0; i < step; ++i){
    ///bufa = bufc;
    prop4D_factor(bufa, qlat::ComplexT<Td >(0.0, 0.0));
    for(Int mu=0;mu<3;mu++)
    {
      qlat::ComplexT<Td >* src = (qlat::ComplexT<Td >*) qlat::get_data(bufc).data();
      qlat::ComplexT<Td >* r0  = (qlat::ComplexT<Td >*) qlat::get_data(buf0).data();
      qlat::ComplexT<Td >* r1  = (qlat::ComplexT<Td >*) qlat::get_data(buf1).data();
      ////svec.shift_vecs_dir(bufc, buf0, mu, -1);
      ////svec.shift_vecs_dir(bufc, buf1, mu, +1);
      std::vector<Int > iDir(4);for(Int i=0;i<4;i++){iDir[i] = 0;}
      iDir[mu] = +1;
      svec.shift_vecP(src, r0, iDir, 12*12);
      iDir[mu] = -1;
      svec.shift_vecP(src, r1, iDir, 12*12);

      bufa += buf0;
      bufa += buf1;
    }
    /////for(Int i=0;i<nsites;i++){wmp[i]  = norm*(wmp[i] + bw*buf[i]);}
    prop4D_factor(bufc, norm);
    prop4D_factor(bufa, b1  );
    bufc += bufa;
  }
  //prop4D_factor(prop, qlat::ComplexT<Td >(0.0, 0.0));
  //prop += bufc;
  prop = bufc;
}

template <class Td>
void smear_propagator_gwu_convension_1shift(Propagator4dT<Td>& prop, const GaugeFieldT<Td >& gf, const double width, const Int step)
{
  const Geometry& geo = gf.geo();
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);

  shift_vec svec(fd, true);
  qlat::vector_gpu<qlat::ComplexT<Td > > gauge;
  {
    gauge.resize(2*4* gf.geo().local_volume() * 9);
    extend_links_to_vecs(gauge.data(), gf);
    svec.set_gauge(gauge.data(), 4,12);
  }
  Propagator4dT<Td> buf0;buf0.init(geo);
  Propagator4dT<Td> buf1;buf1.init(geo);
  Propagator4dT<Td> bufc;bufc.init(geo);
  Propagator4dT<Td> bufa;bufa.init(geo);
  qlat::ComplexT<Td > b0 = qlat::ComplexT<Td >( 1.0 - width, 0.0);
  qlat::ComplexT<Td > b1 = qlat::ComplexT<Td >( width/6.0, 0.0);

  prop4D_factor(bufc, qlat::ComplexT<Td >(0.0, 0.0));
  bufc += prop;
  for (Int i = 0; i < step; ++i){
    ///bufa = bufc;
    prop4D_factor(bufa, qlat::ComplexT<Td >(0.0, 0.0));
    for(Int mu=0;mu<3;mu++)
    {
      qlat::ComplexT<Td >* src = (qlat::ComplexT<Td >*) qlat::get_data(bufc).data();
      qlat::ComplexT<Td >* r0  = (qlat::ComplexT<Td >*) qlat::get_data(buf0).data();
      qlat::ComplexT<Td >* r1  = (qlat::ComplexT<Td >*) qlat::get_data(buf1).data();
      std::vector<Int > iDir(4);for(Int i=0;i<4;i++){iDir[i] = 0;}
      iDir[mu] = +1;
      svec.shift_vecP(src, r0, iDir, 12*12);
      iDir[mu] = -1;
      svec.shift_vecP(src, r1, iDir, 12*12);

      bufa += buf0;
      bufa += buf1;
    }
    /////for(Int i=0;i<nsites;i++){wmp[i]  = norm*(wmp[i] + bw*buf[i]);}
    prop4D_factor(bufc, b0);
    prop4D_factor(bufa, b1);
    bufc += bufa;
  }
  ////prop4D_factor(prop, qlat::ComplexT<Td >(0.0, 0.0));
  prop = bufc;
}

template <class Td>
void smear_propagator_gwu_convension_2shift_ori(Propagator4dT<Td>& prop, const GaugeFieldT<Td >& gf, const double width, const Int step)
{
  const Geometry& geo = gf.geo();
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);

  shift_vec svec(fd, true);
  qlat::vector_gpu<qlat::ComplexT<Td > > gauge;
  {
    gauge.resize(2*4* gf.geo().local_volume() * 9);
    extend_links_to_vecs(gauge.data(), gf);
    svec.set_gauge(gauge.data(), 4,12);
  }
  Propagator4dT<Td> buf0;buf0.init(geo);
  Propagator4dT<Td> buf1;buf1.init(geo);
  Propagator4dT<Td> bufc;bufc.init(geo);
  Propagator4dT<Td> bufa;bufa.init(geo);
  qlat::ComplexT<Td > b0 = qlat::ComplexT<Td >( 1.0 - width, 0.0);
  qlat::ComplexT<Td > b1 = qlat::ComplexT<Td >( width/6.0, 0.0);

  prop4D_factor(bufc, qlat::ComplexT<Td >(0.0, 0.0));
  bufc += prop;
  for (Int i = 0; i < step; ++i){
    ///bufa = bufc;
    prop4D_factor(bufa, qlat::ComplexT<Td >(0.0, 0.0));
    for(Int mu=0;mu<3;mu++)
    {
      qlat::ComplexT<Td >* src = (qlat::ComplexT<Td >*) qlat::get_data(bufc).data();
      qlat::ComplexT<Td >* r0  = (qlat::ComplexT<Td >*) qlat::get_data(buf0).data();
      qlat::ComplexT<Td >* r1  = (qlat::ComplexT<Td >*) qlat::get_data(buf1).data();
      std::vector<Int > iDir(4);for(Int i=0;i<4;i++){iDir[i] = 0;}
      iDir[mu] = +2;
      svec.shift_vecP(src, r0, iDir, 12*12);
      iDir[mu] = -2;
      svec.shift_vecP(src, r1, iDir, 12*12);

      bufa += buf0;
      bufa += buf1;
    }
    /////for(Int i=0;i<nsites;i++){wmp[i]  = norm*(wmp[i] + bw*buf[i]);}
    prop4D_factor(bufc, b0);
    prop4D_factor(bufa, b1);
    bufc += bufa;
  }
  ////prop4D_factor(prop, qlat::ComplexT<Td >(0.0, 0.0));
  prop = bufc;
}

template <typename Ta, typename Tb>
void prepare_gauge_buffer(std::vector< GaugeFieldT<Ta>   >& gfL, const GaugeFieldT<Tb >& gf)
{
  TIMERA("prepare_gauge_buffer");
  if(gfL.size() != 8)gfL.resize(8);
  const Geometry& geo = gf.geo();
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);
  Coordinate total_siteh = Coordinate(fd.nx/2, fd.ny/2, fd.nz/2, fd.nt);

  Qassert(fd.Nx % 2 == 0 and fd.Ny % 2 == 0 and fd.Nz % 2 == 0);

  //Geometry geoh;geoh.init(total_siteh);
  const Geometry& geoh = get_geo_cache(total_siteh);

  const Coordinate expan_left( 0, 0, 0, 0);
  const Coordinate expan_right(1, 1, 1, 0);
  const Int dir_limit = 3;
  GaugeField gf1;
  gf1.init(geo_resize(gf.geo(), expan_left, expan_right));
  gf1 = gf;
  refresh_expanded(gf1);

  for(Int z=0;z<2;z++)
  for(Int y=0;y<2;y++)
  for(Int x=0;x<2;x++)
  {
    const long i = (z*2+y)*2+x;
    gfL[i].init();gfL[i].init(geoh);
    GaugeFieldT<Ta >& gfh = gfL[i];
    qlat::set_zero(gfh);

    qacc_for(index, gfh.geo().local_volume(), {
      const Coordinate xlh = gfh.geo().coordinate_from_index(index);
      const Coordinate xl1 = Coordinate(xlh[0]*2 + x, xlh[1]*2 + y, xlh[2]*2 + z, xlh[3]);

      for (Int dir = 0; dir < dir_limit; ++dir) {
        ////auto&
        const Coordinate xld = coordinate_shifts(xl1, dir);
        qlat::ComplexT<Ta>*  res = (qlat::ComplexT<Ta>*) &gfh.get_elem(xlh, dir)(0,0);
        qlat::ComplexT<Tb>*  s0  = (qlat::ComplexT<Tb>*) &gf1.get_elem(xl1, dir)(0,0);
        qlat::ComplexT<Tb>*  s1  = (qlat::ComplexT<Tb>*) &gf1.get_elem(xld, dir)(0,0);
        for(Int a=0;a<3;a++)
        for(Int b=0;b<3;b++)
        {
          res[a*3+b] = 0.0;
          for(Int c=0;c<3;c++)
          {
            res[a*3+b] += s0[a*3+c] * s1[c*3+b];
          }
        }
      }
    });
  }
}

template <class Ty, const Int dir>
void prop_full_copy_to_prop_8(std::vector< qlat::vector_gpu<Ty > >& propL, Ty* prop, const Int Ndata,
  const Geometry& geo, const Geometry& geoh,
  const Int even, const Int idx = 0, Int NdataP = -1)
{
  TIMERA("prop_full_copy_to_prop_8");
  if(propL.size()!= 8 and dir == 1){propL.resize(8);}
  if(dir == 0){Qassert(propL.size() == 8);}
  if(NdataP == -1){NdataP = Ndata;}
  Qassert(NdataP <= Ndata);
  Qassert(idx*NdataP <  Ndata);
  Qassert((idx + 1)*NdataP <= Ndata);

  for(Int z=0;z<2;z++)
  for(Int y=0;y<2;y++)
  for(Int x=0;x<2;x++)
  {
    const long i = (z*2+y)*2+x;
    long writei = i;
    if(even != -1 and even != i){continue;}
    if(even != -1){writei = 0;}

    const size_t Nfull = geoh.local_volume() * Ndata;
    if(dir == 1){if(propL[writei].size() < Nfull){propL[writei].resizeL(Nfull);}}
    if(dir == 0){Qassert(propL[writei].size() >= Nfull);}

    Ty* proph = (Ty*) propL[writei].data();
    qacc_for(index, geoh.local_volume(), {
      const Coordinate xlh = geoh.coordinate_from_index(index);
      long sh = geoh.index_from_coordinate(xlh);

      const Coordinate xl1 = Coordinate(xlh[0]*2 + x, xlh[1]*2 + y, xlh[2]*2 + z, xlh[3]);
      long sp = geo.index_from_coordinate(xl1);
      if(dir == 1)for(Int di=0;di<NdataP;di++){proph[sh*Ndata + idx*NdataP + di] = prop[sp*NdataP + di];}
      if(dir == 0)for(Int di=0;di<NdataP;di++){prop[sp*NdataP + di] = proph[sh*Ndata + idx*NdataP + di];}
    });
  }
}

// wuppertal convension
template <class Ty, class Td, Int c0, Int d0>
void smear_propagator_gwu_convension_2shift_modi(std::vector< qlat::vector_gpu<Ty > >& propL, const Geometry& geo, std::vector< GaugeFieldT<Td> >& gfL, const double width, const Int step, const Int even = -1, const Int force_update = -1, const Int wuppertal_conv = 1, const CoordinateD& mom = CoordinateD())
{
  (void)geo;
  if(width == 0 or step == 0){return ;}
  Qassert(propL.size() == 8);
  Qassert(gfL.size() == 8);
  const Geometry& geoh = gfL[0].geo();
  const Long Ndata = c0 * 3 * d0;

  double factor_sigma = width;//width can be negative for normalization
  if(wuppertal_conv == 1){
    factor_sigma = std::sqrt( 2 * step * width / 3.0);
  }
  for(Int i=0;i<8;i++){
    long writei = i;if(even != -1 and even != i){continue;}if(even != -1){writei = 0;}

    Qassert(Long(propL[writei].size()) >= geoh.local_volume() * Ndata);
    Ty* proph = (Ty*) propL[writei].data();
    smear_propagator_gwu_convension_inner<Ty, c0, d0, Td>(proph, gfL[i], factor_sigma, step, mom, false, 1, i, force_update);
  }
}


template <class Ty, class Td, Int c0, Int d0>
void smear_propagator_gwu_convension_2shift_modi(Ty* prop, const Geometry& geo, std::vector< GaugeFieldT<Td> >& gfL, const double width, const Int step, std::vector< qlat::vector_gpu<Ty > >& propL, const Int even = -1, const Int wuppertal_conv = 1, const CoordinateD& mom = CoordinateD())
{
  Qassert(gfL.size() == 8);
  const Geometry& geoh = gfL[0].geo();
  if(width == 0 or step == 0){return ;}

  const Int Ndata = c0 * 3 * d0;
  prop_full_copy_to_prop_8<Ty, 1>(propL, prop, Ndata, geo, geoh, even);
  smear_propagator_gwu_convension_2shift_modi<Ty, Td, c0, d0>(propL, geo, gfL, width, step, even, -1, wuppertal_conv, mom);
  prop_full_copy_to_prop_8<Ty, 0>(propL, prop, Ndata, geo, geoh, even);
}

// wuppertal convension
template <class Td>
void smear_propagator_wuppertal_convension_2shift(Propagator4dT<Td>& prop, std::vector< GaugeFieldT<Td> >& gfL, const double width, const Int step, std::vector< qlat::vector_gpu<qlat::ComplexT<Td > > >& propL, const Int even = -1, const Int wuppertal_conv = 1, const CoordinateD& mom = CoordinateD())
{
  Qassert(prop.initialized);
  if(width == 0 or step == 0){return ;}
  rotate_prop(prop, 0);
  qlat::ComplexT<Td >* src = (qlat::ComplexT<Td >*) qlat::get_data(prop).data();
  smear_propagator_gwu_convension_2shift_modi<qlat::ComplexT<Td >, Td, 1, 12*4>(src, prop.geo(), gfL, width, step, propL, even, wuppertal_conv, mom);
  rotate_prop(prop, 1);
}

// gwu convension
template <class Td>
void smear_propagator_gwu_convension_2shift(Propagator4dT<Td>& prop, std::vector< GaugeFieldT<Td> >& gfL, const double width, const Int step, std::vector< qlat::vector_gpu<qlat::ComplexT<Td > > >& propL, const Int even = -1, const CoordinateD& mom = CoordinateD())
{
  const Int wuppertal_conv = 0;
  smear_propagator_wuppertal_convension_2shift(prop, gfL, width, step, propL, even, wuppertal_conv, mom);
}

template <class Ty, class Td>
void smear_propagator_gwu_convension_2shift(FieldG<Ty >& prop, std::vector< GaugeFieldT<Td> >& gfL,
         const double width, const Int step, std::vector< qlat::vector_gpu<qlat::ComplexT<Td > > >& propL, const Int even = -1, const CoordinateD& mom = CoordinateD())
{
  if (0 == step) {return;}
  Qassert(prop.initialized and prop.multiplicity == 12 * 12 and prop.mem_order == QLAT_OUTTER);
  const Long Nvol = prop.geo().local_volume();
  Ty* src = (Ty*) qlat::get_data(prop).data();
  move_index mv_civ;int flag = 0;
  flag = 0;mv_civ.dojob(src, src, 1, 12*12, Nvol, flag, 1, false);

  qacc_for(isp, Nvol, {
    Ty buf[12*12];
    for(Int i=0;i<12*12;i++){buf[i] = src[isp*12*12 + i];}
    for(Int d0=0;d0<12*4;d0++)
    for(Int c0=0;c0<   3;c0++)
    {
      src[isp*12*12 + c0*12*4 + d0] = buf[d0*3 + c0];
    }
  });

  smear_propagator_gwu_convension_2shift_modi<Ty, Td, 1, 12*4>(src, prop.geo(), gfL, 
    width, step, propL, even, 0, mom);

  qacc_for(isp, Nvol, {
    Ty buf[12*12];
    for(Int i=0;i<12*12;i++){buf[i] = src[isp*12*12 + i];}
    for(Int d0=0;d0<12*4;d0++)
    for(Int c0=0;c0<   3;c0++)
    {
      src[isp*12*12 + d0*3 + c0] = buf[c0*12*4 + d0];
    }
  });
  flag = 1;mv_civ.dojob(src, src, 1, 12*12, Nvol, flag, 1, false);
}


// gwu convension with all cs outer prop : 4*3 * V * complex
template <class Ty, class Td>
void smear_propagator_gwu_convension_2shift(qpropT& prop, std::vector< GaugeFieldT<Td> >& gfL,
         const double width, const Int step, std::vector< qlat::vector_gpu<qlat::ComplexT<Td > > >& propL, const Int even = -1, const CoordinateD& mom = CoordinateD())
{
  if (0 == step) {return;}
  Qassert(prop.initialized);
  Ty* src = (Ty*) qlat::get_data(prop).data();
  FieldG<Ty> tmp_prop;
  const Long Nd = 12 * 12 * prop.geo().local_volume();
  tmp_prop.set_pointer(src, Nd, prop.geo(), QMGPU, QLAT_OUTTER);
  smear_propagator_gwu_convension_2shift(tmp_prop, gfL, width, step, propL, even, mom);
}

template <typename Ty, typename Td>
void smear_propagator_gwu_convension_2shift_modi(colorFT& prop, std::vector< GaugeFieldT<Td> >& gfL, const double width, const Int step, std::vector< qlat::vector_gpu<Ty > >& propL, const Int even = -1, const Int wuppertal_conv = 1, const CoordinateD& mom = CoordinateD())
{
  Qassert(prop.initialized);
  Ty* src = (Ty*) qlat::get_data(prop).data();
  smear_propagator_gwu_convension_2shift_modi<Ty, Td, 1, 1>(src, prop.geo(), gfL, width, step, propL, even, wuppertal_conv, mom);
}

template <typename Ty, typename Td>
void smear_propagator_gwu_convension_2shift_modi(std::vector<colorFT>& prop, std::vector< GaugeFieldT<Td> >& gfL, const double width, const Int step, std::vector< qlat::vector_gpu<Ty > >& propL, const Int even = -1, const Int wuppertal_conv = 1, const CoordinateD& mom = CoordinateD())
{
  for(unsigned int i=0;i<prop.size();i++){
    smear_propagator_gwu_convension_2shift_modi(prop[i], gfL, width, step, propL, even, wuppertal_conv, mom);
  }
}

}
#endif
