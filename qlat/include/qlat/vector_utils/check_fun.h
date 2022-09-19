// Gen Wang
// Jun. 2021

#ifndef Check_FUN_H
#define Check_FUN_H
#pragma once

#include "general_funs.h"
#include "utils_reduce_vec.h"

namespace qlat
{

template <class T>
void diff_gauge( GaugeFieldT<T> &g0, GaugeFieldT<T> &g1)
{
  double diff = 0.0;int count_print = 0;
  for (long index = 0; index < g0.geo().local_volume(); ++index) {
    Vector<ColorMatrixT<T> > v0 = g0.get_elems(index);
    Vector<ColorMatrixT<T> > v1 = g1.get_elems(index);
    for(int m = 0; m < 4; ++m){
      const double *p0= (double*) &v0[m](0,0);
      const double *p1= (double*) &v1[m](0,0);
      for(int pi=0;pi<9*2;pi++)
      {
        diff += std::fabs(p0[pi]-p1[pi]);
        if(std::fabs(p0[pi]-p1[pi])>1e-6 and count_print < 100){
          print0("Wrong %.5e %.5e \n",p0[pi],p1[pi]);
          count_print += 1;
        }
      }
    }
  }
  sum_all_size(&diff,0);
  print0("===diff conf %.5e \n",diff);
}

template <class T>
void diff_prop(Propagator4dT<T>& p0, Propagator4dT<T>& p1, double err=1e-15)
{
  int rank = qlat::get_id_node();
  long MAX_COUNT = 64;
  double diffp = 0.0; long countp = 0;
  for (long index = 0; index < p0.geo().local_volume(); ++index) {
    Coordinate xl0 = p0.geo().coordinate_from_index(index);
    Coordinate xg0 = p0.geo().coordinate_g_from_l(xl0);

    qlat::WilsonMatrixT<T>&  s0 =  p0.get_elem(index);
    qlat::WilsonMatrixT<T>&  s1 =  p1.get_elem(index);
    for(int d0=0;d0<12;d0++)
    for(int d1=0;d1<12;d1++)
    {
      T p0 = s0(d0,d1);
      T p1 = s1(d0,d1);

      double diff = 0.0;
      double n0 = qlat::qnorm(p0);
      double n1 = qlat::qnorm(p1);
      bool checknan = false;
      if(std::isnan(p0.real()) or std::isnan(p0.imag())){checknan = true;}
      if(std::isnan(p1.real()) or std::isnan(p1.imag())){checknan = true;}

      if(n1 > 1e-18){diff = std::fabs(qlat::qnorm(p0-p1)/n1);}
      else{
      if(n0 > 1e-18){diff = std::fabs(qlat::qnorm(p0-p1)/n0);}}
      diffp += diff;
      if((diff > err or checknan) and countp < MAX_COUNT )
      {
        printf("rank %5d, x %5d, y %5d, z %5d, t %5d, d %5d, value %+.3e %+.3e, %+.3e %+.3e, %+.3e \n",
          rank ,xg0[0],xg0[1],xg0[2],xg0[3],d0*12+d1 ,p0.real(),p0.imag(),p1.real(),p1.imag(), diff);
        countp += 1;
      }
    }

  }
  sum_all_size(&diffp,1);
  MPI_Barrier(get_comm());fflush(stdout);
  print0("==prop diff %.5e \n",diffp/(p0.geo().local_volume()*12*24.0));
  MPI_Barrier(get_comm());fflush(stdout);
}

template <class Ty, int civ>
void diff_propT(qlat::FieldM<Ty , civ>& p0, qlat::FieldM<Ty , civ>& p1, double err=1e-15, long MAX_COUNT = 64)
{
  int rank = qlat::get_id_node();
  /////long MAX_COUNT = 64;
  double diffp = 0.0; long countp = 0;
  Ty* s0 = (Ty*) qlat::get_data(p0).data();
  Ty* s1 = (Ty*) qlat::get_data(p1).data();
  for (long index = 0; index < p0.geo().local_volume(); ++index) {
    Coordinate xl0 = p0.geo().coordinate_from_index(index);
    Coordinate xg0 = p0.geo().coordinate_g_from_l(xl0);

    for(int c0=0;c0<civ;c0++)
    {
      Ty& pa = s0[index*civ + c0];
      Ty& pb = s1[index*civ + c0];

      double diff = 0.0;
      double n0 = qlat::qnorm(pa);
      double n1 = qlat::qnorm(pb);
      bool checknan = false;
      if(std::isnan(pa.real()) or std::isnan(pa.imag())){checknan = true;}
      if(std::isnan(pb.real()) or std::isnan(pb.imag())){checknan = true;}

      if(n1 > 1e-18){diff = std::fabs(qlat::qnorm(pa-pb)/n1);}
      else{
      if(n0 > 1e-18){diff = std::fabs(qlat::qnorm(pa-pb)/n0);}}
      diffp += diff;
      if((diff > err or checknan) and countp < MAX_COUNT )
      {
        printf("rank %5d, x %5d, y %5d, z %5d, t %5d, d %5d, value %+.8e %+.8e, %+.3e %+.3e, %+.3e \n",
          rank ,xg0[0],xg0[1],xg0[2],xg0[3], c0 ,pa.real(),pa.imag(),pb.real(),pb.imag(), diff);
        countp += 1;
      }
    }
  }
  sum_all_size(&diffp,1);
  MPI_Barrier(get_comm());fflush(stdout);
  print0("==prop diff %.5e \n",diffp/(p0.geo().local_volume()*qlat::get_num_node()*civ));
  MPI_Barrier(get_comm());fflush(stdout);
}

template <class T, int bfac>
void print_src(qlat::FieldM<T, bfac> &noi)
{
  for (long index = 0; index < noi.geo().local_volume(); ++index) {
    Coordinate xl0 = noi.geo().coordinate_from_index(index);
    Coordinate xg0 = noi.geo().coordinate_g_from_l(xl0);

    for(int bi=0;bi<bfac;bi++){
      qlat::Complex tem = noi.get_elems(index)[bi]; 
      double sum = 0.0;
      sum += std::fabs(tem.real());
      sum += std::fabs(tem.imag());
      if(sum > 1e-6){
        printf("==src x %5d, y %5d, z %5d, t %5d, bi %5d, value %.5e %.5e \n",
          xg0[0],xg0[1],xg0[2],xg0[3], bi,tem.real(),tem.imag());
      }
    }

  }

}

template <class T>
void diff_EigenM(qlat::vector<T >& a, qlat::vector<T >& b, std::string lab)
{
  double diff = 0.0;
  if(a.size() != b.size()){print0("%s size not equal %d, %d\n",
    lab.c_str(), int(a.size()), int(b.size()) );abort_r("");}

  for(long i=0;i < a.size(); i++)
  {
    Complexq tem = a[i] - b[i];
    diff += (tem.real()*tem.real() + tem.imag()*tem.imag());
  }

  sum_all_size(&diff,1);
  print0("%s diff %.5e, count %d, ava %.5e \n",
    lab.c_str(),diff, int(a.size()), diff/a.size());
}

template <class T>
void random_point_src(Propagator4dT<T>& prop, int seed = 0)
{
  set_zero(prop);
  qlat::RngState rs(seed);
  int rank = qlat::get_id_node();
  const qlat::Geometry& geo = prop.geo();
  Coordinate xg;
  for(int i=0;i<4;i++){
    int nl = geo.node_site[i] * geo.geon.size_node[i];
    xg[i] = int(nl*qlat::u_rand_gen(rs));
  }
  qlat::displayln_info(qlat::ssprintf("===Point src x %5d, y %5d, z %5d, t %5d \n",xg[0], xg[1], xg[2], xg[3]));
  const Coordinate xl = geo.coordinate_l_from_g(xg);
  if(geo.is_local(xl)){
    long index = geo.index_from_coordinate(xl);
    qlat::WilsonMatrixT<T>&  s0 =  prop.get_elem(index);
    for(int d0=0;d0<12;d0++)s0(d0,d0) = 1.0;
    printf("===set value, node %d, index %ld \n", rank, index );
  }

}

template <typename Ty0, typename Ty1 , int civ>
void get_mom_apply(qlat::FieldM<Ty0, civ> &src, std::vector<int >& mom, std::vector<Ty1 >& dat, bool dir=true, bool ft4D=false)
{
  Geometry& geo = src.geo();
  std::vector<int > Nv,nv,mv;
  Nv.resize(4);nv.resize(4);mv.resize(4);
  for(int i=0;i<4;i++){Nv[i]=geo.node_site[i];nv[i] = geo.node_site[i] * geo.geon.size_node[i];}
  for(int i=0;i<4;i++){mv[i] = nv[i]/Nv[i];}

  Ty0* P0 = (Ty0*) (qlat::get_data(src).data());int psum = 3;
  if( ft4D){dat.resize(civ);psum = 4;}
  if(!ft4D){dat.resize(nv[3]*civ);psum = 3;}
  for(long di=0;di<long(dat.size());di++){dat[di] = 0.0;}

  double pi = 3.14159265358979323846;
  for(long isp = 0; isp < geo.local_volume(); isp++)
  {
    Coordinate xl = geo.coordinate_from_index(isp);
    Coordinate xg = geo.coordinate_g_from_l(xl);
    double tem = 0.0;
    for(int i=0;i<psum;i++){tem += (2*pi/nv[i]) * xg[i] * mom[i];}
    Ty1 phase = 0.0;
    if(!dir)phase = Ty1(std::cos(tem), std::sin(tem));
    if( dir)phase = Ty1(std::cos(tem), -1.0*std::sin(tem));
    for(long di=0;di<civ;di++)
    {
      if( ft4D)dat[            di] += phase * P0[isp*civ + di];
      if(!ft4D)dat[xg[3]*civ + di] += phase * P0[isp*civ + di];
    }
  }
  sum_all_size(&dat[0], dat.size());

}

template <typename Ty0, typename Ty1 , int civ>
void get_mom_fft(qlat::FieldM<Ty0, civ> &src, std::vector<int >& mom, std::vector<Ty1 >& dat, bool ft4D=false)
{
  Geometry& geo = src.geo();
  std::vector<int > Nv,nv,mv;
  Nv.resize(4);nv.resize(4);mv.resize(4);
  for(int i=0;i<4;i++){Nv[i]=geo.node_site[i];nv[i] = geo.node_site[i] * geo.geon.size_node[i];}
  for(int i=0;i<4;i++){mv[i] = nv[i]/Nv[i];}

  Ty0* P0 = (Ty0*) (qlat::get_data(src).data());
  if( ft4D){dat.resize(civ);}
  if(!ft4D){dat.resize(nv[3]*civ);}
  for(long di=0;di<long(dat.size());di++){dat[di] = 0.0;}

  for(long isp = 0; isp < geo.local_volume(); isp++)
  {
    Coordinate xl = geo.coordinate_from_index(isp);
    Coordinate xg = geo.coordinate_g_from_l(xl);
    if( ft4D)
    if(xg[0] == mom[0] and xg[1] == mom[1] and xg[2] == mom[2] and xg[3] == mom[3])
    {for(int di=0;di<civ;di++){dat[di] += P0[isp*civ + di];}}
    if(!ft4D)
    if(xg[0] == mom[0] and xg[1] == mom[1] and xg[2] == mom[2])
    {for(int di=0;di<civ;di++){dat[xg[3]*civ + di] += P0[isp*civ + di];}}
  }

  sum_all_size(&dat[0], dat.size());

}


template<typename Ty>
void print_sum(const Ty* a, size_t size, std::string stmp=std::string(""), int GPU = 1)
{
  qlat::vector_acc<Ty > sum;sum.resize(1);sum[0] = 0.0;
  Ty* pres = sum.data();
  if(GPU == 1){reduce_vec(a, pres , size, 1);}
  if(GPU == 0){reduce_cpu(a, pres , size, 1);}
  print0("%s, sum %.3e \n", stmp.c_str(), qlat::qnorm(sum[0]) );
}



}

#endif
