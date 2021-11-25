// utils_smear_vecs.h
// Gen Wang
// Sep. 2021

#ifndef UTILS_SMEAR_VECS_H
#define UTILS_SMEAR_VECS_H
#pragma once

#include "general_funs.h"
#include <qlat/qcd-utils.h>
#include <qlat/qcd-prop.h>
 
namespace qlat{

void get_smear_para(std::string smear_para, double& width, int& step)
{
  if(smear_para == std::string("NONE")){width = 0.0;step = 0;return; }
  std::vector<std::string > temL = stringtolist(smear_para);
  if(temL.size() != 2){abort_r("read smear para wrong! \n");}
  width = stringtodouble(temL[0]);
  step  = stringtonum(temL[1]);
}

template <class T>
void smear_propagator_gwu_convension(Propagator4dT<T>& prop, const GaugeFieldT<T>& gf1,
                      const double width, const int step,int mode=0)
{
  // gf1 is left_expanded and refreshed
  // set_left_expanded_gauge_field(gf1, gf)
  // prop is of qnormal size
  ////For complex numbers addition and subtraction require two flops, and multiplication and division require six flops
  ///complex multi 6 + plus 2
  TIMER_FLOPS("==smear propagator");
  long long Tfloat = 0;
  double mem       = 0.0;

  {long long Lat = prop.geo().local_volume();
  int nsrc = 12;
  long long vGb = Lat *nsrc*4;
  int Fcount = 3*(3*6 + 2*2); 
  int direction   = 6;
  Tfloat = step*direction*vGb*Fcount;
  mem = (Lat*nsrc*12 + Lat*4*9)*8.0;}
  timer.flops += Tfloat;
  print0("Memory size %.3e GB, %.3e Gflop \n", 
    mem/(1024.0*1024*1024), Tfloat/(1024.0*1024*1024));

  if (0 == step) {
    return;
  }
  const double aw   = 3.0*width*width/(2*step);
  const double coef = aw;
  if(mode == 0)smear_propagator(prop, gf1, coef, step);
  ////if(mode == 1)smear_propagator1(prop, gf1, coef, step);

}

template <class Ty, int civ>
void smear_fieldM_gwu_convension_simple(std::vector<qlat::FieldM<Ty, civ> >& src, const GaugeFieldT<Ty >& gf1,
                      const double width, const int step,
                      const CoordinateD& mom = CoordinateD(),
                      const bool smear_in_time_dir = false)
{
  // gf1 is left_expanded and refreshed
  // set_left_expanded_gauge_field(gf1, gf)
  // prop is of qnormal size
  ////For complex numbers addition and subtraction require two flops, and multiplication and division require six flops
  ///complex multi 6 + plus 2
  if(civ % 3 != 0){abort_r("FieldM smearing size wrong! \n");}

  if (0 == step) {
    return;
  }

  const double aw   = 3.0*width*width/(2*step);
  const double coef = aw;

  const Geometry& geo = src[0].geo();
  const Geometry geo1 =
      smear_in_time_dir
          ? geo_resize(geo, 1)
          : geo_resize(geo, Coordinate(1, 1, 1, 0), Coordinate(1, 1, 1, 0));
  const int n_avg = smear_in_time_dir ? 8 : 6;
  const int dir_limit = smear_in_time_dir ? 4 : 3;
  array<Ty, 8> mom_factors;
  for (int i = 0; i < 8; ++i) {
    const int dir = i - 4;
    const double phase = dir >= 0 ? mom[dir] : -mom[-dir - 1]; 
    mom_factors[i] = std::polar(coef / n_avg, -phase);
  }

  qlat::FieldM<Ty, civ> src1;
  src1.init(geo1);
  for(LInt iv=0;iv < src.size(); iv++){
  for (int i = 0; i < step; ++i) {
    src1 = src[iv];
    qlat::FieldM<Ty, civ>& tmp = src[iv];
    src1 = tmp;
    refresh_expanded_1(src1);
    qacc_for(index, geo.local_volume(), {
      const Coordinate xl = geo.coordinate_from_index(index);
      Ty* wm = (Ty*) qlat::get_data(tmp.get_elems(xl)).data();
      //Ty* wm = (Ty*) tmp.get_elem(xl).p;

      Eigen::Matrix<Ty, 3, civ/3, Eigen::ColMajor>& wmE = *((Eigen::Matrix<Ty, 3, civ/3, Eigen::ColMajor>*)  wm);
      wmE *= (1-coef);
      //WilsonMatrixT<T>& wm = tmp.get_elem(xl);
      //wm *= 1 - coef;
      //for(int ci=0;ci<civ;ci++){wm[ci] = (1-coef)*wm[ci];}
      for (int dir = -dir_limit; dir < dir_limit; ++dir) {
        const Coordinate xl1 = coordinate_shifts(xl, dir);
        const ColorMatrixT<Ty> link =
            dir >= 0
                ? gf1.get_elem(xl, dir)
                : (ColorMatrixT<Ty>)matrix_adjoint(gf1.get_elem(xl1, -dir - 1));
        const Ty* linkP = (Ty*) link.p;
        //const Ty* linkP = (Ty*) qlat::get_data(link).data();
        //Ty* src1p = (Ty*) qlat::get_data(src1.get_elem(xl1)).data();
        Ty* src1p = (Ty*) qlat::get_data(src1.get_elems(xl1)).data();

        Eigen::Matrix<Ty, 3, 3, Eigen::RowMajor>&  lE = *((Eigen::Matrix<Ty, 3, 3, Eigen::RowMajor>*) linkP);
        Eigen::Matrix<Ty, 3, civ/3, Eigen::ColMajor>&  pE = *((Eigen::Matrix<Ty, 3, civ/3, Eigen::ColMajor>*) src1p);
        wmE += mom_factors[dir + 4] * (lE * pE);

      }
    });
  }
  }


}


template <class Ty, int civ>
void smear_fieldM_gwu_convension(std::vector<qlat::FieldM<Ty, civ> >& src, const GaugeFieldT<Ty >& gf1,
                      const double width, const int step,
                      const CoordinateD& mom = CoordinateD(),
                      const bool smear_in_time_dir = false, int mode = 0)
{

  TIMER_FLOPS("==smear FieldM");
  long long Tfloat = 0;
  double mem       = 0.0;
  if(civ % 3 != 0){abort_r("FieldM smearing size wrong! \n");}
  const int bfac = civ/3;

  {long long Lat = src[0].geo().local_volume();
  int nsrc = src.size();
  long long vGb = Lat * nsrc * bfac;
  int Fcount = 3*(3*6 + 2*2); 
  int direction   = 6;
  Tfloat = step*direction*vGb*Fcount;
  mem = (Lat*nsrc*12 + Lat*bfac*9)*8.0;}
  timer.flops += Tfloat;
  print0("Memory size %.3e GB, %.3e Gflop \n", 
    mem/(1024.0*1024*1024), Tfloat/(1024.0*1024*1024));

  if(mode == 0){
    smear_fieldM_gwu_convension_simple(src, gf1, width, step, mom, smear_in_time_dir);
    return;}

}





}
#endif
