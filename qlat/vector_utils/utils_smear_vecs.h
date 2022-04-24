// utils_smear_vecs.h
// Gen Wang
// Sep. 2021

#ifndef UTILS_SMEAR_VECS_H
#define UTILS_SMEAR_VECS_H
#pragma once

#include "general_funs.h"
#include <qlat/qcd-utils.h>
#include <qlat/qcd-prop.h>
#include <qlat/qcd-smear.h>
#include "utils_shift_vecs.h"
#include "check_fun.h"
 
namespace qlat{

template <class T>
void rotate_prop(Propagator4dT<T>& prop, int dir = 0)
{
  TIMER("Rotate color prop");
  const Geometry& geo = prop.geo();
  long long Nvol =  geo.local_volume();

  qacc_for(isp, long(Nvol), {
    qlat::WilsonMatrixT<T>& v0 =  prop.get_elem(isp);
    qlat::WilsonMatrixT<T>  v1 =  prop.get_elem(isp);

    for(int c0 = 0;c0 < 3; c0++)
    for(int d0 = 0;d0 < 4; d0++)
    for(int c1 = 0;c1 < 3; c1++)
    for(int d1 = 0;d1 < 4; d1++)
    {
      //int a0 = (d0*3+c0)*12+ d1*3+c1;
      //LInt off = ((d0*4+d1)*3 + c0)*3 + c1;
      LInt off = (c0*3+c1)*16+d0*4 + d1;
      int a0 = off/12;
      int a1 = off%12;
      if(dir == 0)v0(a0, a1) = v1(d0*3+c0, d1*3+c1);
      if(dir == 1)v0(d0*3+c0, d1*3+c1) = v1(a0, a1);

      //if(dir==0)pE[off*Nvol + isp] = v0(d0*3 + c0, d1*3 + c1);
      //if(dir==1)v0(d0*3 + c0, d1*3 + c1) = pE[off*Nvol + isp];
    }
  });
}

////TODO need to change other parts for c0, c1
template <class T, class Tg>
void extend_links_to_vecs(qlat::vector_acc<T >& gfE, const GaugeFieldT<Tg >& gf){
  // set_left_expanded_gauge_field(gf1, gf)
  // prop is of qnormal size
  TIMERB("extend_links_to_vecs");
  const Geometry& geo = gf.geo();
  GaugeFieldT<Tg > gf1;
  set_left_expanded_gauge_field(gf1, gf);

  long Nvol = geo.local_volume();
  const int dir_limit = 4;
  gfE.resize(2*dir_limit*Nvol*9);
  qacc_for(index,  geo.local_volume(),{
    for (int dir = -dir_limit; dir < dir_limit; ++dir) {
      const Coordinate xl = geo.coordinate_from_index(index);
      const ColorMatrixT<Tg > link =
          dir >= 0 ? gf1.get_elem(xl, dir)
                   : (ColorMatrixT<Tg >)matrix_adjoint(
                         gf1.get_elem(coordinate_shifts(xl, dir), -dir - 1));
      for(int ci=0; ci<9; ci++){
        ///gfE[(index*dir_limit*2+ (dir+dir_limit))*9 + ci] = link.p[ci];
        gfE[(index*dir_limit*2+ (dir+dir_limit))*9 + (ci%3)*3 + ci/3] = link.p[ci];
      }
    }
  });
}

void get_smear_para(std::string smear_para, double& width, int& step)
{
  if(smear_para == std::string("NONE")){width = 0.0;step = 0;return; }
  std::vector<std::string > temL = stringtolist(smear_para);
  if(temL.size() != 2){abort_r("read smear para wrong! \n");}
  width = stringtodouble(temL[0]);
  step  = stringtonum(temL[1]);
}

template <class T>
void smear_propagator_box1_kernel(Propagator4dT<T>& prop, std::vector<Propagator4dT<T> >& vA,
  const int Bsize, const int dir, qlat::vector_acc<T >& gfE)
{
  TIMER("Shift Kernel");qassert(dir >= 0);
  //std::vector<int > iDir(4);
  //for(int i=0;i<4;i++){iDir[i] = 0;}
  const Geometry& geo = prop.geo();

  const int dir_max = 4;
  std::vector<int > dirL(2);dirL[0] = -dir-1;dirL[1] = dir;
  std::vector<int > dirG(2);dirG[0] = dirL[0] + dir_max;dirG[1] = dirL[1] + dir_max;

  /////Buffer index
  const long Nvol = geo.local_volume();
  std::vector<long > Pdir0;Pdir0.resize(Nvol);
  std::vector<long > Pdir1;Pdir1.resize(Nvol * 2);
  qacc_for(index, geo.local_volume(),{
    const Coordinate xl0 = prop.geo().coordinate_from_index(index);
    Pdir0[index] = prop.geo().offset_from_coordinate(xl0);
    for(int di=0;di<2;di++)
    {
      const Coordinate xl1 = coordinate_shifts(xl0, dirL[di]);
      Pdir1[index*2 + di] = vA[2 + dir*2 + di].geo().offset_from_coordinate(xl1);
    }
  });
  T* prop0P  = (T*) qlat::get_data(prop).data();
  std::vector<T* > vAP;vAP.resize(2);
  std::vector<T* > vBP;vBP.resize(2);

  for(unsigned int i=0;i<2;i++){vAP[i] = (T* ) qlat::get_data(vA[i]).data();}
  for(unsigned int i=0;i<2;i++){vBP[i] = (T* ) qlat::get_data(vA[2 + dir*2 + i]).data();}

  vA[0] = prop;vA[1] = prop;
  for (int i = 0; i < Bsize; ++i) {
    TIMER("Matrix multiply");
    //vp = prop;
    vA[2 + dir*2 + 0] = vA[0];vA[2 + dir*2 + 1] = vA[1];
    refresh_expanded_1(vA[2 + dir*2 + 0]);
    refresh_expanded_1(vA[2 + dir*2 + 1]);

    qacc_for(index,  geo.local_volume(),{

      std::vector<T* > tmp(2);
      for(int di = 0; di < 2; di++){
        tmp[di] = &vAP[di][size_t(Pdir0[index])*12*12];
        Eigen::Matrix<T, 3, 3*16, Eigen::RowMajor>& pAE = *((Eigen::Matrix<T, 3, 3*16, Eigen::RowMajor>*)  tmp[di]);
        const T* wmAP = &vBP[di][size_t(Pdir1[index*2 + di])*12*12];

        const T* lp = &gfE[(index*dir_max*2 + dirG[di])*9];
        Eigen::Matrix<T, 3, 3, Eigen::ColMajor>&     lE = *((Eigen::Matrix<T, 3, 3, Eigen::ColMajor>*) lp);
        Eigen::Matrix<T, 3, 3*16, Eigen::RowMajor>&  pE = *((Eigen::Matrix<T, 3, 3*16, Eigen::RowMajor>*) wmAP);
        ////wmE += (lE * pE);
        pAE = lE*pE;
      }
      T* wmp = &prop0P[size_t(Pdir0[index])*12*12];
      for(int i=0;i<3*3*16;i++){wmp[i] += (tmp[0][i] + tmp[1][i]);}
    });
  }

}


template <class T>
void smear_propagator_box_kernel(Propagator4dT<T>& prop, Propagator4dT<T>& vp, Propagator4dT<T>& vm,
  const int Bsize, const int dir, shift_vec& svec)
{
  TIMER("Shift Kernel");
  std::vector<int > iDir(4);
  for(int i=0;i<4;i++){iDir[i] = 0;}
  const Geometry& geo = prop.geo();
  vp = prop;vm = prop;
  for(int bi=0;bi<Bsize;bi++){
    iDir[dir] =  1;shift_fieldM(svec, vp, vp, iDir);
    iDir[dir] = -1;shift_fieldM(svec, vm, vm, iDir);

    qacc_for(index, geo.local_volume(), {
      const Coordinate xl = geo.coordinate_from_index(index);
      WilsonMatrixT<T>& wm = prop.get_elem(xl);
      wm += vp.get_elem(xl);
      wm += vm.get_elem(xl);
    })
  }
}

////expansion of field as needed, not shift on positions
template <class T, class Tg>
void smear_propagator_box1(Propagator4dT<T>& prop, const GaugeFieldT<Tg >& gf, const int Bsize){
  TIMER_VERBOSE("smear_propagator_box");
  if (Bsize <= 0) {return;}
  const Geometry& geo = prop.geo();
  const int dir_limit = 3;
  std::vector<int > nv,Nv,mv;geo_to_nv(geo, nv, Nv, mv);

  fft_desc_basic fd(geo);
  qlat::vector_acc<T > gfE; 
  extend_links_to_vecs(gfE, gf);

  rotate_prop(prop,0);

  std::vector<Propagator4dT<T> > vL(3);
  for(unsigned int i=0;i<vL.size();i++){vL[i].init(geo);}
  Propagator4dT<T > v1;v1.init(geo);
  Propagator4dT<T > v2;v2.init(geo);

  std::vector<Propagator4dT<T> > vA;vA.resize(2 + 8);
  for(unsigned int i=0;i<2;i++){vA[i].init(geo);}

  for(unsigned int dir=0;dir<4;dir++){
    Coordinate expansion_left( 0, 0, 0, 0);
    Coordinate expansion_right(0, 0, 0, 0);
    Coordinate center(0, 0, 0, 0);
    expansion_left[dir] = 1;expansion_right[dir] = 1;
    const Geometry geoL = geo_resize(prop.geo(), expansion_left, center );
    const Geometry geoR = geo_resize(prop.geo(), center, expansion_right);
    ////TODO Bugs with init for different geos
    vA[2 + dir*2 + 0].init(geoL);
    vA[2 + dir*2 + 1].init(geoR);
  }

  //Coordinate shift;
  for(int dir = 0; dir < dir_limit; ++dir)
  {
    v1 = prop;
    smear_propagator_box1_kernel(v1, vA, Bsize, dir, gfE);
    int d1=(dir+1)%3;int d2=(dir+2)%3;v2 = v1;
    smear_propagator_box1_kernel(v1, vA, Bsize, d1 , gfE);
    smear_propagator_box1_kernel(v2, vA, Bsize, d2 , gfE);
    vL[d1] += v2;
    vL[d2] += v1;
  }
  qlat::set_zero(prop);
  for(int dir = 0; dir < dir_limit; ++dir)
  {
    smear_propagator_box1_kernel(vL[dir], vA, Bsize, dir, gfE);
    prop += vL[dir];
  }
  T* data = (T*) qlat::get_data(prop).data();
  qacc_for(index, geo.local_volume(), {
    for(int c=0;c<12*12;c++){data[index*12*12 + c] *= 1.0/6.0;} 
  });

  rotate_prop(prop,1);
}


template <class T, class Tg>
void smear_propagator_box(Propagator4dT<T>& prop, const GaugeFieldT<Tg >& gf, const int Bsize){
  TIMER_VERBOSE("smear_propagator_box");
  if (Bsize <= 0) {return;}
  /////qassert(!smear_in_time_dir);
  const Geometry& geo = prop.geo();
  //const int dir_limit = smear_in_time_dir ? 4 : 3;
  const int dir_limit = 3;
  std::vector<int > nv,Nv,mv;geo_to_nv(geo, nv, Nv, mv);

  fft_desc_basic fd(geo);
  shift_vec svec(fd, true);
  qlat::vector_acc<T > gfE; 
  extend_links_to_vecs(gfE, gf);
  svec.set_gauge(gfE.data(), 4, 12);

  Propagator4dT<T> vp;vp.init(geo);
  Propagator4dT<T> vm;vm.init(geo);
  std::vector<Propagator4dT<T> > vL(3);
  for(unsigned int i=0;i<vL.size();i++){vL[i].init(geo);}
  Propagator4dT<T> v1;v1.init(geo);
  Propagator4dT<T> v2;v2.init(geo);
  //Coordinate shift;
  for(int dir = 0; dir < dir_limit; ++dir)
  {
    v1 = prop;
    smear_propagator_box_kernel(v1, vp, vm, Bsize, dir, svec);
    int d1=(dir+1)%3;int d2=(dir+2)%3;v2 = v1;
    smear_propagator_box_kernel(v1, vp, vm, Bsize, d1 , svec);
    smear_propagator_box_kernel(v2, vp, vm, Bsize, d2 , svec);
    vL[d1] += v2;
    vL[d2] += v1;
  }
  qlat::set_zero(prop);
  for(int dir = 0; dir < dir_limit; ++dir)
  {
    smear_propagator_box_kernel(vL[dir], vp, vm, Bsize, dir, svec);
    prop += vL[dir];
  }
  T* data = (T*) qlat::get_data(prop).data();
  qacc_for(index, geo.local_volume(), {
    for(int c=0;c<12*12;c++){data[index*12*12 + c] *= 1.0/6.0;} 
  });

}

template <class T >
void smear_propagator1(Propagator4dT<T>& prop, const qlat::vector_acc<T >& gf,
                      const double width, const int step)
{

  const double aw   = 3.0*width*width/(2*step);
  const T bw = width*width/(4.0*step - 6.0*width*width);

  const Geometry& geo = prop.geo();
  const long Nvol = geo.local_volume();
  const Geometry geo1 = geo_resize(geo, Coordinate(1, 1, 1, 0), Coordinate(1, 1, 1, 0));
  const int dir_limit = 3;
  const int dir_max   = 4;
  rotate_prop(prop);

  Propagator4dT<T> prop1;
  prop1.init(geo1);

  /////Buffer index
  std::vector<long > Pdir;Pdir.resize(Nvol);
  std::vector<long > Pdir1;Pdir1.resize(Nvol * dir_limit*2);
  qacc_for(index, geo.local_volume(),{
    const Coordinate xl = geo.coordinate_from_index(index);
    Pdir[index] = prop.geo().offset_from_coordinate(xl);
    for (int dir = -dir_limit; dir < dir_limit; ++dir) {
      const Coordinate xl1 = coordinate_shifts(xl, dir);
      Pdir1[index*dir_limit*2 + (dir + dir_limit)] = prop1.geo().offset_from_coordinate(xl1);
    }
  });
  /////Buffer index

  T* propP  = (T*) qlat::get_data(prop).data();
  T* prop1P = (T*) qlat::get_data(prop1).data();

  for (int i = 0; i < step; ++i) {
    prop1 = prop;
    refresh_expanded_1(prop1);

    {
    TIMER("Matrix multiply");
    qacc_for(index,  geo.local_volume(),{
      T buf[3*3*16];for(int i=0;i<3*3*16;i++){buf[i] = 0;}
      Eigen::Matrix<T, 3, 3*16, Eigen::RowMajor>& wmE = *((Eigen::Matrix<T, 3, 3*16, Eigen::RowMajor>*)  buf);

      for (int dir = -dir_limit; dir < dir_limit; ++dir) {

        const T* wm1p = &prop1P[size_t(Pdir1[index*dir_limit*2 + (dir + dir_limit)])*12*12];
        const T* lp = &gf[(index*dir_max*2 + dir+dir_max)*9];
        Eigen::Matrix<T, 3, 3, Eigen::ColMajor>&     lE = *((Eigen::Matrix<T, 3, 3, Eigen::ColMajor>*) lp);
        Eigen::Matrix<T, 3, 3*16, Eigen::RowMajor>&  pE = *((Eigen::Matrix<T, 3, 3*16, Eigen::RowMajor>*) wm1p);
        wmE += (lE * pE);
      }
      T* wmp = &propP[size_t(Pdir[index])*12*12];
      for(int i=0;i<3*3*16;i++){wmp[i] += bw*buf[i];}
    });
    }
  }

  const T factor = std::pow(1-aw, step); 
  T* res = (T*) qlat::get_data(prop).data();
  qacc_for(index,  geo.local_volume(),{
    for(int i=0;i<12*12;i++){res[index*12*12 + i] *= factor;}
    //const Coordinate xl = geo.coordinate_from_index(index);
    //WilsonMatrixT<T >& wm = prop.get_elem(xl);
    //for(int i=0;i++
    //wm *= factor;
  });

  rotate_prop(prop,1);
}

template <class T, class Tg>
void smear_propagator_gwu_convension(Propagator4dT<T>& prop, const GaugeFieldT<Tg >& gf,
                      const double width, const int step, const int mode = 0)
{
  // gf is left_expanded and refreshed
  // set_left_expanded_gauge_field(gf, gf)
  // prop is of qnormal size
  ////For complex numbers addition and subtraction require two flops, and multiplication and division require six flops
  ///complex multi 6 + plus 2
  #if PRINT_TIMER>4
  TIMER_FLOPS("==smear propagator");
  long long Tfloat = 0;
  ///double mem       = 0.0;

  {
    long long Lat = prop.geo().local_volume();
    int nsrc = 12;
    long long vGb = Lat *nsrc*4;
    int Fcount = 3*(3*6 + 2*2); 
    if(step >= 0){
    int direction   = 6;
    Tfloat = step*direction*vGb*Fcount;
    }else{
    Tfloat = 12*2*int(width)*vGb*Fcount;
    }
    //mem = (Lat*nsrc*12 + Lat*4*9)*8.0;
  }
  timer.flops += Tfloat;
  //print0("Memory size %.3e GB, %.3e Gflop \n", 
  //  mem/(1024.0*1024*1024), Tfloat/(1024.0*1024*1024));
  #else
  TIMER("==smear propagator")
  #endif

  if (0 == step) {
    return;
  }
  if(step >= 0){
    //if(mode == 0){
    //  const double coef =  3.0*width*width/(2*step);
    //  GaugeFieldT<Tg > gf1;
    //  set_left_expanded_gauge_field(gf1, gf);
    //  Propagator4dT<qlat::Complex > propD;propD.init(prop.geo());
    //  propD = prop;
    //  smear_propagator(    propD, gf1, coef , step      );
    //  prop = propD;
    //}
    if(mode == 1){
      qlat::vector_acc<T > gfE; 
      ////TODO Need phase for links
      extend_links_to_vecs(gfE, gf);
      smear_propagator1(prop, gfE, width, step);          
    }
  }

  if(step ==-1){smear_propagator_box(prop, gf,int(width + 0.001));}
  //////slow version, need improvement
  //if(step ==-1){smear_propagator_box1(prop, gf,int(width + 0.001));}
}

//template <class Ty, int civ>
//void smear_fieldM_gwu_convension_simple(std::vector<qlat::FieldM<Ty, civ> >& src, const GaugeFieldT<Ty >& gf1,
//                      const double width, const int step,
//                      const CoordinateD& mom = CoordinateD(),
//                      const bool smear_in_time_dir = false)
//{
//  // gf1 is left_expanded and refreshed
//  // set_left_expanded_gauge_field(gf1, gf)
//  // prop is of qnormal size
//  ////For complex numbers addition and subtraction require two flops, and multiplication and division require six flops
//  ///complex multi 6 + plus 2
//  if(civ % 3 != 0){abort_r("FieldM smearing size wrong! \n");}
//
//  if (0 == step) {
//    return;
//  }
//
//  const double coef =  3.0*width*width/(2*step);
//
//  const Geometry& geo = src[0].geo();
//  const Geometry geo1 =
//      smear_in_time_dir
//          ? geo_resize(geo, 1)
//          : geo_resize(geo, Coordinate(1, 1, 1, 0), Coordinate(1, 1, 1, 0));
//  const int n_avg = smear_in_time_dir ? 8 : 6;
//  const int dir_limit = smear_in_time_dir ? 4 : 3;
//  array<Ty, 8> mom_factors;
//  for (int i = 0; i < 8; ++i) {
//    const int dir = i - 4;
//    const double phase = dir >= 0 ? mom[dir] : -mom[-dir - 1]; 
//    mom_factors[i] = std::polar(coef / n_avg, -phase);
//  }
//
//  qlat::FieldM<Ty, civ> src1;
//  src1.init(geo1);
//  for(LInt iv=0;iv < src.size(); iv++){
//  for (int i = 0; i < step; ++i) {
//    src1 = src[iv];
//    qlat::FieldM<Ty, civ>& tmp = src[iv];
//    src1 = tmp;
//    refresh_expanded_1(src1);
//    qacc_for(index, geo.local_volume(), {
//      const Coordinate xl = geo.coordinate_from_index(index);
//      Ty* wm = (Ty*) qlat::get_data(tmp.get_elems(xl)).data();
//      //Ty* wm = (Ty*) tmp.get_elem(xl).p;
//
//      Eigen::Matrix<Ty, 3, civ/3, Eigen::ColMajor>& wmE = *((Eigen::Matrix<Ty, 3, civ/3, Eigen::ColMajor>*)  wm);
//      wmE *= (1-coef);
//      //WilsonMatrixT<T>& wm = tmp.get_elem(xl);
//      //wm *= 1 - coef;
//      //for(int ci=0;ci<civ;ci++){wm[ci] = (1-coef)*wm[ci];}
//      for (int dir = -dir_limit; dir < dir_limit; ++dir) {
//        const Coordinate xl1 = coordinate_shifts(xl, dir);
//        const ColorMatrixT<Ty> link =
//            dir >= 0
//                ? gf1.get_elem(xl, dir)
//                : (ColorMatrixT<Ty>)matrix_adjoint(gf1.get_elem(xl1, -dir - 1));
//        const Ty* linkP = (Ty*) link.p;
//        //const Ty* linkP = (Ty*) qlat::get_data(link).data();
//        //Ty* src1p = (Ty*) qlat::get_data(src1.get_elem(xl1)).data();
//        Ty* src1p = (Ty*) qlat::get_data(src1.get_elems(xl1)).data();
//
//        Eigen::Matrix<Ty, 3, 3, Eigen::RowMajor>&  lE = *((Eigen::Matrix<Ty, 3, 3, Eigen::RowMajor>*) linkP);
//        Eigen::Matrix<Ty, 3, civ/3, Eigen::ColMajor>&  pE = *((Eigen::Matrix<Ty, 3, civ/3, Eigen::ColMajor>*) src1p);
//        wmE += mom_factors[dir + 4] * (lE * pE);
//
//      }
//    });
//  }
//  }
//
//
//}
//
//
//template <class Ty, int civ>
//void smear_fieldM_gwu_convension(std::vector<qlat::FieldM<Ty, civ> >& src, const GaugeFieldT<Ty >& gf1,
//                      const double width, const int step,
//                      const CoordinateD& mom = CoordinateD(),
//                      const bool smear_in_time_dir = false, int mode = 0)
//{
//
//  TIMER_FLOPS("==smear FieldM");
//  long long Tfloat = 0;
//  double mem       = 0.0;
//  if(civ % 3 != 0){abort_r("FieldM smearing size wrong! \n");}
//  const int bfac = civ/3;
//
//  {long long Lat = src[0].geo().local_volume();
//  int nsrc = src.size();
//  long long vGb = Lat * nsrc * bfac;
//  int Fcount = 3*(3*6 + 2*2); 
//  int direction   = 6;
//  Tfloat = step*direction*vGb*Fcount;
//  mem = (Lat*nsrc*12 + Lat*bfac*9)*8.0;}
//  timer.flops += Tfloat;
//  print0("Memory size %.3e GB, %.3e Gflop \n", 
//    mem/(1024.0*1024*1024), Tfloat/(1024.0*1024*1024));
//
//  if(mode == 0){
//    smear_fieldM_gwu_convension_simple(src, gf1, width, step, mom, smear_in_time_dir);
//    return;}
//
//}





}
#endif
