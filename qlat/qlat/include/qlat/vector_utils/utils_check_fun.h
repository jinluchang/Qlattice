// Gen Wang
// Jun. 2021

#ifndef UTILS_CHECK_FUN_H
#define UTILS_CHECK_FUN_H
#pragma once

#include "general_funs.h"
#include "utils_reduce_vec.h"
#include "utils_field_operations.h"

namespace qlat
{

template <class Ta, class Tb>
double diff_gauge_extended( GaugeFieldT<Ta> &g0, GaugeFieldT<Tb> &g1)
{
  TIMER("diff_gauge_extended");
  const Geometry& geo = g0.geo();
  const Long V = geo.local_volume_expanded();

  const Ta* p0 = (Ta*) get_data(g0).data();
  const Tb* p1 = (Tb*) get_data(g1).data();

  qlat::vector<Ta > dL;dL.resize(V);
  const Long Nd = 4*9*2;
  qacc_for(index, V, {
    Ta diff = 0.0;
    for(Long m=0;m<Nd;m++){
      diff += qfabs(p0[index*Nd + m] - p1[index*Nd + m]);
    }
    dL[index] = diff;
  });

  double diff = Reduce(dL.data(), dL.size(), true);
  diff = diff/(g0.geo().local_volume()*4*9*2.0);
  qmessage("===diff conf %.5e \n",diff);
  return double(diff);
}


template <class Ta, class Tb>
double diff_gauge( GaugeFieldT<Ta> &g0, GaugeFieldT<Tb> &g1, double err=1e-6)
{
  TIMER("diff_gauge");
  const Geometry& geo = g0.geo();
  double diff = 0.0;int count_print = 0;
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    for(Int m = 0; m < 4; ++m){
      const Ta *p0 = (Ta*)  g0.get_elem(xl, m).p;
      const Tb *p1 = (Tb*)  g1.get_elem(xl, m).p;
      for(Int pi=0;pi<9*2;pi++)
      {
        diff += std::fabs(p0[pi]-p1[pi]);
        if(std::fabs(p0[pi]-p1[pi]) > err and count_print < 300){
          Coordinate xl0 = geo.coordinate_from_index(index);
          Coordinate xg0 = geo.coordinate_g_from_l(xl0);

          qmessage("Wrong %3d %3d %3d %3d, dir %1d, ids %3d, %+.8e %+.8e \n", xg0[0], xg0[1], xg0[2], xg0[3], m, pi,p0[pi],p1[pi]);
          count_print += 1;
        }
      }
    }
  }
  sum_all_size(&diff, 1);
  diff = diff/(g0.geo().local_volume()*4*9*2.0);
  qmessage("==prop diff %.5e \n", diff);
  MPI_Barrier(get_comm());fflush(stdout);
  return diff;
}

template <class Ta, class Tb>
double diff_gauge_GPU( GaugeFieldT<Ta> &g0, GaugeFieldT<Tb> &g1)
{
  TIMER("diff_gauge");
  const Geometry& geo = g0.geo();
  const Long V = geo.local_volume();
  vector<Ta > dL;dL.resize(V);
  qacc_for(index, V, {
    dL[index] = 0.0;
    const Coordinate xl = geo.coordinate_from_index(index);
    for(Int m = 0; m < 4; ++m){
      const Ta *p0 = (Ta*)  g0.get_elem(xl, m).p;
      const Tb *p1 = (Tb*)  g1.get_elem(xl, m).p;
      for(Int pi=0;pi<9*2;pi++)
      {
        dL[index] += qfabs(p0[pi]-p1[pi]);
      }
    }
  });

  double diff = Reduce(dL.data(), dL.size(), true);
  diff = diff / (g0.geo().local_volume()*4*9*2.0);
  qmessage("===diff conf %.5e \n",diff);
  return double(diff);
}

template <class Ta, class Tb>
double diff_fields(Field<Ta>& p0, Field<Tb>& p1)
{
  Qassert(p0.geo() == p1.geo());
  Qassert(p0.multiplicity == p1.multiplicity);
  //int rank = qlat::get_id_node();
  const Int Dim = p0.multiplicity;
  const Long V = p0.geo().local_volume();
  Field<double > fd;fd.init(p0.geo(), 1);
  qacc_for(index, V, {
    Ta* r0 = (Ta*) p0.get_elems(index).p;
    Tb* r1 = (Tb*) p1.get_elems(index).p;
    double d0 = 0.0;
    for(Int i=0;i<Dim;i++){
      d0 += qnorm(r0[i] - r1[i]);
    }
    fd.get_elem(index) = d0;
  });
  double diff = fields_quick_checksum(fd, 8, false);
  return diff;
}

template <class Ta, class Tb>
void diff_prop(Propagator4dT<Ta>& p0, Propagator4dT<Tb>& p1, double err=1e-15)
{
  Int rank = qlat::get_id_node();
  Long MAX_COUNT = 64;
  double diffp = 0.0; Long countp = 0;
  for (Long index = 0; index < p0.geo().local_volume(); ++index) {
    Coordinate xl0 = p0.geo().coordinate_from_index(index);
    Coordinate xg0 = p0.geo().coordinate_g_from_l(xl0);

    qlat::WilsonMatrixT<Ta>&  s0 =  p0.get_elem_offset(index);
    qlat::WilsonMatrixT<Tb>&  s1 =  p1.get_elem_offset(index);
    for(Int d0=0;d0<12;d0++)
    for(Int d1=0;d1<12;d1++)
    {
      qlat::ComplexT<Ta > p0 = s0(d0,d1);
      qlat::ComplexT<Tb > p1 = s1(d0,d1);

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
  qmessage("==prop diff %.5e \n",diffp/(p0.geo().local_volume()*12*24.0));
  MPI_Barrier(get_comm());fflush(stdout);
}

template <class Ty, Int civ>
void print_norm2(qlat::FieldM<Ty , civ>& p0)
{
  const Long V = p0.geo().local_volume() * civ;
  Ty* src = (Ty*) qlat::get_data(p0).data();
  qlat::vector_gpu<Ty > buf;
  buf.copy_from(src, V, 1, 1);
  buf.print_norm2();
}


template <class Ty, Int civ>
void diff_propT(qlat::FieldM<Ty , civ>& p0, qlat::FieldM<Ty , civ>& p1, double err=1e-15, Long MAX_COUNT = 64)
{
  Int rank = qlat::get_id_node();
  /////Long MAX_COUNT = 64;
  double diffp = 0.0; Long countp = 0;
  Ty* s0 = (Ty*) qlat::get_data(p0).data();
  Ty* s1 = (Ty*) qlat::get_data(p1).data();
  for (Long index = 0; index < p0.geo().local_volume(); ++index) {
    Coordinate xl0 = p0.geo().coordinate_from_index(index);
    Coordinate xg0 = p0.geo().coordinate_g_from_l(xl0);

    for(Int c0=0;c0<civ;c0++)
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
  sum_all_size(&diffp, 1, 0);
  MPI_Barrier(get_comm());fflush(stdout);
  qmessage("==prop diff %.5e \n",diffp/(p0.geo().local_volume()*qlat::get_num_node()*civ));
  MPI_Barrier(get_comm());fflush(stdout);
}

template <class T, Int bfac>
void print_src(qlat::FieldM<T, bfac> &noi)
{
  for (Long index = 0; index < noi.geo().local_volume(); ++index) {
    Coordinate xl0 = noi.geo().coordinate_from_index(index);
    Coordinate xg0 = noi.geo().coordinate_g_from_l(xl0);

    for(Int bi=0;bi<bfac;bi++){
      qlat::ComplexD tem = noi.get_elems(index)[bi]; 
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
  if(a.size() != b.size()){qmessage("%s size not equal %d, %d\n",
    lab.c_str(), int(a.size()), int(b.size()) );abort_r("");}

  for(Long i=0;i < a.size(); i++)
  {
    Complexq tem = a[i] - b[i];
    diff += (tem.real()*tem.real() + tem.imag()*tem.imag());
  }

  sum_all_size(&diff,1);
  qmessage("%s diff %.5e, count %d, ava %.5e \n",
    lab.c_str(),diff, int(a.size()), diff/a.size());
}

template <class T>
void random_point_src(Propagator4dT<T>& prop, Int seed = 0)
{
  set_zero(prop);
  qlat::RngState rs(seed);
  Int rank = qlat::get_id_node();
  const qlat::Geometry& geo = prop.geo();
  Coordinate xg;
  for(Int i=0;i<4;i++){
    Int nl = geo.node_site[i] * geo.geon.size_node[i];
    xg[i] = int(nl*qlat::u_rand_gen(rs));
  }
  qlat::displayln_info(qlat::ssprintf("===Point src x %5d, y %5d, z %5d, t %5d \n",xg[0], xg[1], xg[2], xg[3]));
  const Coordinate xl = geo.coordinate_l_from_g(xg);
  if(geo.is_local(xl)){
    Long index = geo.index_from_coordinate(xl);
    qlat::WilsonMatrixT<T>&  s0 =  prop.get_elem_offset(index);
    for(Int d0=0;d0<12;d0++)s0(d0,d0) = 1.0;
    printf("===set value, node %d, index %ld \n", rank, (long)index );
  }

}

template <typename Ty0, typename Ty1 , Int civ>
void get_mom_apply(qlat::FieldM<Ty0, civ> &src, std::vector<Int >& mom, std::vector<Ty1 >& dat, bool dir=true, bool ft4D=false)
{
  Geometry& geo = src.geo();
  std::vector<Int > Nv,nv,mv;
  Nv.resize(4);nv.resize(4);mv.resize(4);
  for(Int i=0;i<4;i++){Nv[i]=geo.node_site[i];nv[i] = geo.node_site[i] * geo.geon.size_node[i];}
  for(Int i=0;i<4;i++){mv[i] = nv[i]/Nv[i];}

  Ty0* P0 = (Ty0*) (qlat::get_data(src).data());int psum = 3;
  if( ft4D){dat.resize(civ);psum = 4;}
  if(!ft4D){dat.resize(nv[3]*civ);psum = 3;}
  for(Long di=0;di<Long(dat.size());di++){dat[di] = 0.0;}

  double pi = 3.14159265358979323846;
  for(Long isp = 0; isp < geo.local_volume(); isp++)
  {
    Coordinate xl = geo.coordinate_from_index(isp);
    Coordinate xg = geo.coordinate_g_from_l(xl);
    double tem = 0.0;
    for(Int i=0;i<psum;i++){tem += (2*pi/nv[i]) * xg[i] * mom[i];}
    Ty1 phase = 0.0;
    if(!dir)phase = Ty1(std::cos(tem), std::sin(tem));
    if( dir)phase = Ty1(std::cos(tem), -1.0*std::sin(tem));
    for(Long di=0;di<civ;di++)
    {
      if( ft4D)dat[            di] += phase * P0[isp*civ + di];
      if(!ft4D)dat[xg[3]*civ + di] += phase * P0[isp*civ + di];
    }
  }
  sum_all_size(&dat[0], dat.size());

}

template <typename Ty0, typename Ty1 , Int civ>
void get_mom_fft(qlat::FieldM<Ty0, civ> &src, std::vector<Int >& mom, std::vector<Ty1 >& dat, bool ft4D=false)
{
  Geometry& geo = src.geo();
  std::vector<Int > Nv,nv,mv;
  Nv.resize(4);nv.resize(4);mv.resize(4);
  for(Int i=0;i<4;i++){Nv[i]=geo.node_site[i];nv[i] = geo.node_site[i] * geo.geon.size_node[i];}
  for(Int i=0;i<4;i++){mv[i] = nv[i]/Nv[i];}

  Ty0* P0 = (Ty0*) (qlat::get_data(src).data());
  if( ft4D){dat.resize(civ);}
  if(!ft4D){dat.resize(nv[3]*civ);}
  for(Long di=0;di<Long(dat.size());di++){dat[di] = 0.0;}

  for(Long isp = 0; isp < geo.local_volume(); isp++)
  {
    Coordinate xl = geo.coordinate_from_index(isp);
    Coordinate xg = geo.coordinate_g_from_l(xl);
    if( ft4D)
    if(xg[0] == mom[0] and xg[1] == mom[1] and xg[2] == mom[2] and xg[3] == mom[3])
    {for(Int di=0;di<civ;di++){dat[di] += P0[isp*civ + di];}}
    if(!ft4D)
    if(xg[0] == mom[0] and xg[1] == mom[1] and xg[2] == mom[2])
    {for(Int di=0;di<civ;di++){dat[xg[3]*civ + di] += P0[isp*civ + di];}}
  }

  sum_all_size(&dat[0], dat.size());

}


template<typename Ty>
void print_sum(const Ty* a, size_t size, std::string stmp=std::string(""), Int GPU = 1)
{
  Ty sum = Reduce(a, size, GPU);
  qmessage("%s, sum %.3e \n", stmp.c_str(), qlat::qnorm(sum) );
}

template <typename Ty>
double check_sum_FieldM(qlat::Field<Ty>& p0)
{
  const Geometry& geo = p0.geo();
  double check_sum = 0.0;
  const Int civ = p0.multiplicity;
  const Ty* src = (Ty*) qlat::get_data(p0).data();
  const Long Nvol = geo.local_volume();
  for(Long index = 0; index < Nvol; ++index){
    Coordinate xl = geo.coordinate_from_index(index);
    Coordinate xg = geo.coordinate_g_from_l(xl);
    double cor = xg[0]*0.5 + xg[1]*1.7 + xg[2]*xg[2]*2.5 + xg[3]*xg[3]*0.2;
    for(Long dc = 0;dc<civ;dc++)
    {
      Ty fac = Ty(std::cos(cor + dc) , std::sin(cor*5.7  + dc*dc) );
      check_sum += qlat::qnorm( src[index*civ + dc] * fac ) ;
    }
  }
  sum_all_size(&check_sum,1);
  MPI_Barrier(get_comm());fflush(stdout);
  return check_sum;
}

template <typename Td>
double check_sum_prop(Propagator4dT<Td >& p0)
{
  ////int rank = qlat::get_id_node();
  const Geometry& geo = p0.geo();
  double check_sum = 0.0;
  const ComplexT<Td>* src = (ComplexT<Td>*) qlat::get_data(p0).data();
  const Long Nvol = geo.local_volume();
  for(Long index = 0; index < Nvol; ++index){
    Coordinate xl = geo.coordinate_from_index(index);
    Coordinate xg = geo.coordinate_g_from_l(xl);
    double cor = xg[0]*0.5 + xg[1]*1.7 + xg[2]*xg[2]*2.5 + xg[3]*xg[3]*0.2;
    for(Long dc = 0;dc<12*12;dc++)
    {
      ComplexT<Td> fac = ComplexT<Td>(std::cos(cor + dc) , std::sin(cor*5.7  + dc*dc) );
      check_sum += qlat::qnorm( src[index*12*12 + dc] * fac ) ;
    }
  }
  sum_all_size(&check_sum,1);
  MPI_Barrier(get_comm());fflush(stdout);
  return check_sum;
}

template <typename Ty>
double check_sum_prop(qpropT& p0)
{
  ////int rank = qlat::get_id_node();
  const Geometry& geo = p0.geo();
  double check_sum = 0.0;
  const Ty* src = (Ty*) qlat::get_data(p0).data();
  const Long Nvol = geo.local_volume();
  for(Long index = 0; index < Nvol; ++index){
    Coordinate xl = geo.coordinate_from_index(index);
    Coordinate xg = geo.coordinate_g_from_l(xl);
    double cor = xg[0]*0.5 + xg[1]*1.7 + xg[2]*xg[2]*2.5 + xg[3]*xg[3]*0.2;
    for(Long dc = 0;dc<12*12;dc++)
    {
      Ty fac = Ty(std::cos(cor + dc) , std::sin(cor*5.7  + dc*dc) );
      check_sum += qlat::qnorm( src[dc*Nvol + index] * fac ) ;
    }
  }
  sum_all_size(&check_sum,1);
  MPI_Barrier(get_comm());fflush(stdout);
  /////qmessage("==prop check_sum %.8e \n", check_sum);
  return check_sum;
}

template <class T, Int civ>
double diff_FieldM(qlat::FieldM<T , civ>& prop0, qlat::FieldM<T , civ>& prop1, double err=1e-15)
{
  Int rank = qlat::get_id_node();
  Long MAX_COUNT = 64;
  double diffp = 0.0; Long countp = 0;
  T* s0 = (T*) qlat::get_data(prop0).data();
  T* s1 = (T*) qlat::get_data(prop1).data();
  for (Long index = 0; index < prop0.geo().local_volume(); ++index) {
    Coordinate xl0 = prop0.geo().coordinate_from_index(index);
    Coordinate xg0 = prop0.geo().coordinate_g_from_l(xl0);
    for(Int d0=0;d0<civ;d0++)
    {
      T p0 = s0[index*civ + d0];
      T p1 = s1[index*civ + d0];

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
          rank ,xg0[0],xg0[1],xg0[2],xg0[3],d0 ,p0.real(),p0.imag(),p1.real(),p1.imag(), diff);
        countp += 1;
      }
    }

  }
  sum_all_size(&diffp,1);
  MPI_Barrier(get_comm());fflush(stdout);
  double diff = diffp/(prop0.geo().local_volume()*civ*2.0);
  qmessage("==prop diff %.5e \n", diff);
  MPI_Barrier(get_comm());fflush(stdout);
  return diff;
}


}

#endif
