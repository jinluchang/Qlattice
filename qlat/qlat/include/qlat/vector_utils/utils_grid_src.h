// utils_grid_src.h
// Gen Wang
// Oct. 2021

#ifndef UTILS_GRID_SRC_H
#define UTILS_GRID_SRC_H

#pragma once

#include "utils_float_type.h"
#include "utils_gammas.h"
#include "utils_momentum.h"
#include "utils_fft_desc.h"
#include "utils_field_gpu.h"

namespace qlat{

template <typename Ty>
void get_num_time(qlat::FieldM<Ty, 1>& noise,Int &number_t, Int &t_ini){
  qlat::Geometry& geo = noise.geo();
  qlat::vector<Int > nv,Nv,mv;
  geo_to_nv(geo, nv, Nv, mv);
  //int nx,ny,nz,nt;
  //nx = nv[0];ny = nv[1];nz = nv[2];nt = nv[3];
  Int nt = nv[3];
  LInt Nsite = Nv[0]*Nv[1]*Nv[2]*Nv[3];

  std::vector<double > fullt(nt);for(Int ti=0;ti<nt;ti++){fullt[ti]=0.0;}
  for(unsigned int isp=0; isp< Nsite; isp++)
  {
    ////position p = noise.desc->get_position(isp,get_node_rank());
    Coordinate xl0 = geo.coordinate_from_index(isp);
    Coordinate xg0 = geo.coordinate_g_from_l(xl0);
    {
      ///auto tem_source = noise.data[isp];
      auto tem_source =  noise.get_elem_offset(isp);
      if(qnorm(tem_source)>0.01)
      {
        fullt[xg0[3]] = 1.0;
      }
    }
  }
  sum_all_size((double* ) &fullt[0],nt);
  number_t = 0;
  for(Int ti=0;ti<nt;ti++){if(fullt[ti]>0.0)number_t += 1;}
  for(Int ti=0;ti<nt;ti++){if(fullt[ti]>0.0){t_ini = ti;break;}}
}

/*
  Check grid src positions and offset
*/
template <typename Ty>
void check_noise_pos(qlat::FieldM<Ty, 1>& noise, Coordinate& pos, Coordinate&off_L,Int printS=0,Int mod=0)
{
  qlat::Geometry& geo = noise.geo();
  qlat::vector<Int > nv,Nv,mv;
  geo_to_nv(geo, nv, Nv, mv);
  Int nx,ny,nz,nt;
  nx = nv[0];ny = nv[1];nz = nv[2];nt = nv[3];
  LInt Nsite = Nv[0]*Nv[1]*Nv[2]*Nv[3];

  std::vector<Int > NL(4);NL[0]=nx;NL[1]=ny;NL[2]=nz;NL[3]=nt;
  Int grid_count = 0;
  std::vector<std::vector<Int > > grid;
  for(Int iL=0;iL<4;iL++){
    grid.push_back(std::vector<Int>(NL[iL]));
    for(LInt giL=0;giL<grid[iL].size();giL++){grid[iL][giL] = 0.0;}
  }
  //grid.push_back(std::vector<double > (nx));
  Int number_t = 1;int t_ini = 0;
  if(mod == 1){get_num_time(noise,number_t,t_ini);}
  for(LInt isp=0; isp< Nsite; isp++)
  {
    Coordinate xl0 = geo.coordinate_from_index(isp);
    Coordinate xg0 = geo.coordinate_g_from_l(xl0);
    ////position p = noise.desc->get_position(isp,get_node_rank());
    //int t = xg0[3];
    //int toff = ((t-tini+nt)%nt);

    {
      auto tem_source =  noise.get_elem_offset(isp);
      ////auto tem_source = noise.data[isp];
      if(qnorm(tem_source)>0.01 and xg0[3] < nt/number_t)
      {
        for(Int i=0;i<4;i++){grid[i][xg0[i]] += 1;}
        ///grid[0][p.x()] += 1.0;
        ///grid[1][p.y()] += 1.0;
        ///grid[2][p.z()] += 1.0;
        ///grid[3][p.t()] += 1.0;
        grid_count = grid_count + 1;
      }
    }
  }
  for(Int iL=0;iL<4;iL++){sum_all_size(&grid[iL][0],NL[iL]);}
  sum_all_size(&grid_count,1);
  ////global_sum_all(&grid_count,1);
  ////off_L.resize(4);
  for(Int oi=0;oi<4;oi++){off_L[oi] = 0;}
  for(Int iL=0;iL<4;iL++)for(Int k=0;k<NL[iL];k++)if(grid[iL][k]>0.0)off_L[iL] += 1;
  //for(Int x=0;x<nx;x++){if(grid[0][x]>0.0)off_L[0] += 1;}
  if(int(grid_count) != off_L[0]*off_L[1]*off_L[2]*off_L[3])
  {
    qmessage("Source Check Failed grid_count %10d, offx %5d, offy %5d, offz %5d, offt %5d!\n",
          int(grid_count),off_L[0],off_L[1],off_L[2],off_L[3]);
    Qassert(false);
    ////shutdown_machine();
    ////abort();
  }

  //int pos = 0;int t_ini = 0;
  ////pos.resize(4);
  for(Int i=0;i<4;i++){
    pos[i] = 0;
    for(Int x=0;x<nv[i];x++){if(grid[i][x] > 0){pos[i] += x;break;}}
  }
  t_ini = pos[3];
  //for(Int x=0;x<nx;x++){if(grid[0][x]>0.0){pos += ((x*100)*100)*1000;break;}}
  //for(Int y=0;y<ny;y++){if(grid[1][y]>0.0){pos += (y*100)*1000;break;}}
  //for(Int z=0;z<nx;z++){if(grid[2][z]>0.0){pos += (z)*1000;break;}}
  //for(Int t=0;t<nt;t++){if(grid[3][t]>0.0){pos += (t);t_ini = t;break;}}

  qmessage("Check T %5d %5d %5d %5d, offx %5d, offy %5d, offz %5d, offt %5d. \n",
    pos[0],pos[1],pos[2],pos[3],
    off_L[0],off_L[1],off_L[2],off_L[3]);

  if(printS == 1)
  {
    for(unsigned int isp=0; isp< Nsite; isp++)
    {
      Coordinate xl0 = geo.coordinate_from_index(isp);
      Coordinate p = geo.coordinate_g_from_l(xl0);
      ////position p = noise.desc->get_position(isp,get_node_rank());
      {
        auto tem_source =  noise.get_elem_offset(isp);
        //auto tem_source = noise.data[isp];
        //if(abs(tem_source)>0.01)
        if(qnorm(tem_source)>0.01)
        {
          printf("Check K %5d %5d %5d %5d node %5d %13.5f %13.5f !\n",p[0],p[1],p[2],p[3],qlat::get_id_node()
            , tem_source.real(), tem_source.imag());
        }
      }
    }
    fflush_MPI();
  }
  if(printS == 2)
  {
    /////int x = pos/10000000;int y = (pos%10000000)/100000;int z = (pos%100000)/1000;int t = pos%1000;
    for(unsigned int isp=0; isp< Nsite; isp++)
    {
      Coordinate xl0 = geo.coordinate_from_index(isp);
      Coordinate p = geo.coordinate_g_from_l(xl0);

      ///position p = noise.desc->get_position(isp,get_node_rank());
      {
        auto tem_source =  noise.get_elem_offset(isp);
        ///auto tem_source = noise.data[isp];
        //if(abs(tem_source)>0.01)
        Int printv = 1;
        for(Int i=0;i<4;i++)
        {
          if(p[i] != pos[i] and p[i] != pos[i] + off_L[i]){printv = 0;}
          /////p[i] = pos[i]
        }
        ///if(p[0] == x or p[0] == x + off_L[0])
        ///if(p[1] == y or p[1] == y + off_L[1])
        ///if(p[2] == z or p[2] == z + off_L[2])
        ///if(p[3] == t or p[3] == t + off_L[3])

        if(printv == 1)
        if(qnorm(tem_source)>0.01)
        {
          printf("Check N %5d %5d %5d %5d node %5d %13.5f %13.5f !\n",p[0],p[1],p[2],p[3],qlat::get_id_node()
            ,tem_source.real(),tem_source.imag());
          ////printf("Check N %5d %5d %5d %5d node %5d %13.5f %13.5f %13.5f !\n",p.x(),p.y(),p.z(),p.t(),get_node_rank(),tem_source.real,tem_source.imag);
        }
      }
    }
    ////Check Source Position
    fflush_MPI();
  }

}

template <typename Ty>
qacc Int check_zn(const Ty& value, const Int Zn = 3){
  Int zk = -1;
  for(Int i=0;i<Zn;i++){
    double eta = 2 * QLAT_PI_LOCAL * (i) / 3;
    const Ty tmp = Ty(std::cos(eta), std::sin(eta));
    if(qnorm(value - tmp) < 1e-8){
      zk = i;
      break;
    }
  }
  return zk;
}

/*
  Get all noise positions
  need_copy data Fieldy on GPU or not, need change into mem_type checks
*/
template <class Fieldy>
void get_noise_pos(Fieldy& noise, std::vector<Coordinate >& grids, std::vector<Int >& Zlist,
 const Int Zn = 3, const Int printS=0, const bool need_copy = false)
{
  Qassert(noise.initialized and noise.multiplicity == 1);
  Qassert(GetBasicDataType<Fieldy>::get_type_name() != std::string("unknown_type"));
  using D = typename GetBasicDataType<Fieldy>::ElementaryType;
  Qassert(IsBasicTypeReal<D>());
  using Ty = ComplexT<D>;

  qlat::Geometry& geo = noise.geo();
  Ty* srcP = (Ty*) get_data(noise).data();
  const Long Nvol = geo.local_volume();
  FieldM<Ty, 1> tmp;
  if(need_copy){
    tmp.init(geo);
    Ty* resP = (Ty*) get_data(tmp).data();
    cpy_GPU(resP, srcP, Nvol);
    srcP = resP;
  }

  qlat::vector<Int > nv,Nv,mv;
  geo_to_nv(geo, nv, Nv, mv);
  //int nx,ny,nz,nt;
  //nx = nv[0];ny = nv[1];nz = nv[2];nt = nv[3];
  LInt Nsite = Nv[0]*Nv[1]*Nv[2]*Nv[3];


  //int t_ini = -1;
  std::vector<Int > pos_list;
  std::vector<Int > pos_Zn;
  for(LInt isp=0; isp< Nsite; isp++)
  {
    Coordinate xl0 = geo.coordinate_from_index(isp);
    Coordinate xg0 = geo.coordinate_g_from_l(xl0);
    {
      auto tem_source =  srcP[isp];
      if(qnorm(tem_source)>0.01)
      {
        for(Int i=0;i<4;i++){ pos_list.push_back(xg0[i]);}
        const Int zk = check_zn(tem_source, Zn);
        pos_Zn.push_back(zk);
        //if(t_ini == -1){
        //  t_ini = xg0[3];
        //}
        //if(t_ini != -1){Qassert(t_ini == xg0[3]);}
      }
    }
  }
  const Int Nmpi  = qlat::get_num_node();
  const Int rank  = qlat::get_id_node();
  //std::vector<Int > tL;tL.resize(Nmpi);
  std::vector<Int > gL;gL.resize(Nmpi);
  for(Int i=0;i<Nmpi;i++){
    //tL[i] = 0;
    gL[i] = 0;
  }
  //tL[rank] = t_ini;
  gL[rank] = pos_list.size() / 4;

  sum_all_size(gL.data(), gL.size());
  //sum_all_size(tL.data(), tL.size());
  Long total_pos = 0;
  Long curr_off   = 0;
  for(Int i=0;i<Nmpi;i++){
    //if(t_ini >= 0 and gL[i] != 0){Qassert(tL[i] == t_ini);}
    total_pos += gL[i];
    if(i < rank){curr_off += gL[i];}
  }
  std::vector<Int > zk_mpi ;zk_mpi.resize(total_pos);
  std::vector<Int > pos_mpi;pos_mpi.resize(total_pos * 4);
  for(Long i=0;i<total_pos*4;i++){pos_mpi[i] = 0;}
  for(Long i=0;i<total_pos;i++){zk_mpi[i] = 0;}

  for(Int i=0;i<gL[rank]*4;i++){
    pos_mpi[curr_off*4 + i] = pos_list[i];
  }
  for(Int i=0;i<gL[rank];i++){
    zk_mpi[curr_off + i] = pos_Zn[i];
  }
 
  sum_all_size(pos_mpi.data(), pos_mpi.size());
  sum_all_size(zk_mpi.data() , zk_mpi.size() );
  grids.resize(total_pos);
  Zlist.resize(total_pos);
  for(Long i=0;i<total_pos;i++)
  {
    Zlist[i] = zk_mpi[i];
    for(Int j=0;j<4;j++){
      grids[i][j] = pos_mpi[i*4 + j];
    }
    //if(t_ini >= 0){Qassert(grids[i][0][3] == t_ini);}
    if(printS >= 1)
    {
      const Coordinate p = grids[i];
      const Int zk = Zlist[i];
      qmessage("Check P %5d %5d %5d %5d, Zn %+3d !\n",p[0],p[1],p[2],p[3], zk);
    }
  }

  if(printS >= 2)
  {
    for(unsigned int isp=0; isp< Nsite; isp++)
    {
      Coordinate xl0 = geo.coordinate_from_index(isp);
      Coordinate p = geo.coordinate_g_from_l(xl0);
      {
        auto tem_source =  noise.get_elem_offset(isp);
        if(qnorm(tem_source)>0.01)
        {
          printf("Check K %5d %5d %5d %5d node %5d %13.5f %13.5f !\n",p[0],p[1],p[2],p[3],qlat::get_id_node()
            , tem_source.real(), tem_source.imag());
        }
      }
    }
    fflush_MPI();
  }
}

template <class Fieldy>
void get_noise_posG(Fieldy& noise, std::vector<Coordinate >& grids, std::vector<Int >& Zlist,
 const Int Zn = 3, const Int printS=0, const bool need_copy = false)
{
  (void)need_copy;
  get_noise_pos(noise, grids, Zlist, Zn, printS, true);
}

template <class Fieldy>
void check_noise_high(Fieldy& noise, std::vector<Int >& sinkt, Long& npoints, const bool need_copy = false)
{
  Qassert(noise.initialized and noise.multiplicity == 1);
  const qlat::Geometry& geo = noise.geo();
  ////const Long Nvol = geo.local_volume();
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);

  Qassert(GetBasicDataType<Fieldy>::get_type_name() != std::string("unknown_type"));
  using D = typename GetBasicDataType<Fieldy>::ElementaryType;
  Qassert(IsBasicTypeReal<D>());
  using Ty = ComplexT<D>;

  Ty* srcP = (Ty*) get_data(noise).data();
  const Long Nvol = geo.local_volume();
  FieldM<Ty, 1> tmp;
  if(need_copy){
    tmp.init(geo);
    Ty* resP = (Ty*) get_data(tmp).data();
    cpy_GPU(resP, srcP, Nvol);
    srcP = resP;
  }

  const Long Nxyz = fd.Nx*fd.Ny*fd.Nz;

  std::vector<double > count;
  count.resize(fd.nt);
  for(Int i=0;i<fd.nt;i++){count[i] = 0;}

  qthread_for(ti, fd.Nt, {
    const Int tg = ti + fd.init;

    for(Long isp=0;isp<Nxyz;isp++)
    {
      const Long off = ti * Nxyz + isp;
      {
        double tem_source =  qlat::qnorm(srcP[off]);
        if(qnorm(tem_source)> 1e-3)
        {
          count[tg] += 1;
        }
      }
    }
  });

  sum_all_size(count.data(), count.size());
  sinkt.resize(0);
  npoints = -1;
  for(Int t=0;t<fd.nt;t++){
    if(count[t] > 0){
      if(npoints == -1){
        npoints = count[t];
      }else{
        Qassert(npoints == count[t]);
      }
      sinkt.push_back(t);
    }
  }

}

/////get positions by spatial setups
inline void grid_list_pos(const Coordinate& off_L, qlat::vector<Long >& Ngrid)
{
  TIMERA("===grid_list_pos===")
  if(off_L.size() != 4){abort_r("dimention of positions wrong!\n ");}
  Int Ntot = off_L[0]*off_L[1]*off_L[2];
  double max_d = 1.0*off_L[0]*off_L[0] + 1.0*off_L[1]*off_L[1] + 1.0*off_L[2]*off_L[2] + 1.0;

  Ngrid.resize(Ntot);
  for(Int j0=0;j0<Ntot;j0++)
  {
    double dis_large = 0;
    for(Int i0=0;i0<Ntot;i0++)
    {
      Int ix0= i0/(off_L[1]*off_L[2]);
      Int iy0= (i0%(off_L[1]*off_L[2]))/off_L[2];
      Int iz0= i0%off_L[2];

      double dis_min = max_d;
      for(Int iv=0;iv<j0;iv++)
      {
        //std::vector<Int > temv = map_index[i0];
        Int isrc = Ngrid[iv];
        Int ix= isrc/(off_L[1]*off_L[2]);
        Int iy= (isrc%(off_L[1]*off_L[2]))/off_L[2];
        Int iz= isrc%off_L[2];
        double dis = 0.0;
        dis += 1.0*(ix0-ix)*(ix0-ix);
        dis += 1.0*(iy0-iy)*(iy0-iy);
        dis += 1.0*(iz0-iz)*(iz0-iz);
        if(dis < dis_min){dis_min = dis;}
        //double dis =
      }
      if(dis_min > dis_large){
        Ngrid[j0] = i0;
        ////Ngrid[j0*2+1] = int(dis_min);
        dis_large = dis_min;
      }
    }
    //count += 1;
  }
}

inline Coordinate get_grid_off(Long j0, const Coordinate& off_L, const Coordinate& pos_ini, const Coordinate& Lat)
{
  if(pos_ini.size() != 4 or off_L.size() != 4){abort_r("dimension of positions wrong!\n ");}
  //////std::vector<Int > pos;pos.resize(4);
  Coordinate pos = pos_ini;
  ////for(Int i=0;i<4;i++){pos[i] = pos_ini[i];}

  Coordinate off_pos = qlat::coordinate_from_index(j0, off_L);
  for(Int i=0;i<4;i++){pos[i] += ((Lat[i]/(off_L[i]))*off_pos[i] ) %(Lat[i]);}

  //int it  = j0/(off_L[0]*off_L[1]*off_L[2]);
  //Long i0 = j0%(off_L[0]*off_L[1]*off_L[2]);
  //int ix= i0/(off_L[1]*off_L[2]);
  //int iy= (i0%(off_L[1]*off_L[2]))/off_L[2];
  //int iz= i0%off_L[2];
  //pos[0] += (Lat[0]/(off_L[0]))*ix;
  //pos[1] += (Lat[1]/(off_L[1]))*iy;
  //pos[2] += (Lat[2]/(off_L[2]))*iz;
  //pos[3] += (Lat[3]/(off_L[3]))*it;
  return pos;
}

/////get positions by spatial and time setups
inline void grid_list_posT(std::vector<PointsSelection >& LMS_points, const Coordinate& off_L, const Coordinate& pos, const Int combineT, const Coordinate& Lat)
{
  TIMERA("===grid_list_posT===")
  qlat::vector<Long > Nfull;
  /////get positions by spatial setups
  grid_list_pos(off_L, Nfull);
  LMS_points.resize(0);
  if(combineT == int(0)){LMS_points.resize(Nfull.size()*off_L[3]);}
  if(combineT == int(1)){LMS_points.resize(Nfull.size()         );}

  Coordinate cur_pos;
  Coordinate cur_off;
  Long ci = 0;
  for(Long gi=0;gi<Nfull.size();gi++){
    cur_pos = get_grid_off(Nfull[gi], off_L, pos, Lat);
    if(combineT == 0){
      for(Int it = 0; it < off_L[3]; it++){
        cur_off = cur_pos;
        cur_off[3] += ((Lat[3]/(off_L[3]))*it)%Lat[3]; ////need be careful not to exceed boundaries
        PointsSelection lms_res;lms_res.resize(1);lms_res[0] = cur_off;
        LMS_points[ci] = lms_res;ci += 1;
      }
    }
    if(combineT == 1){
      PointsSelection lms_res;lms_res.resize(off_L[3]);
      for(Int it = 0; it < off_L[3]; it++){
        cur_off = cur_pos;
        cur_off[3] += ((Lat[3]/(off_L[3]))*it)%Lat[3];
        lms_res[it] = cur_off;
      }
      LMS_points[ci] = lms_res;ci += 1;
    }
  }
  //if(combineT == int(0)){Qassert(Long(LMS_points.size()) == Long(Nfull.size()*off_L[3]));}
  //if(combineT == int(1)){Qassert(Long(LMS_points.size()) == Long(Nfull.size()         ));}
}

template<typename Ty>
void write_grid_point_to_src(FieldG<Ty>& res, const FieldG<Ty>& src, const std::vector<Coordinate >& pos, Int b_size, qlat::fft_desc_basic& fd)
{
  (void)b_size;
  TIMER("write_grid_point_to_src");
  Qassert(src.initialized and src.multiplicity == 1);
  const Geometry& geo = src.geo();
  const Int Ndc = 12 * 12 ;
  if(!res.initialized){res.init(geo, Ndc, QMGPU, QLAT_OUTTER);}
  Qassert(res.mem_order == QLAT_OUTTER and res.multiplicity == Ndc);
  set_zero(res);
  const Ty* srcP = (Ty*) get_data(src).data();
        Ty* resP = (Ty*) get_data(res).data();
  const Long V = geo.local_volume();
  vector<Long > listP;listP.resize(pos.size());
  for(Long gi=0;gi<listP.size();gi++){
    const Coordinate& g = pos[gi];
    if(fd.coordinate_g_is_local(g)){
      listP[gi] = fd.index_l_from_g_coordinate(g);
    }else{
      listP[gi] = -1;
    }
  }
  const Long Np = listP.size();
  qacc_for(li, Np, {
    if(listP[li] >= 0){
      const Long isp = listP[li];
      for(Int dc=0;dc<12;dc++){
        resP[(dc*12 + dc) * V + isp] = srcP[isp];
      }
    }
  });
}

//////assume res have been cleared
//////asuume res on GPU
template<typename Ty>
void write_grid_point_to_src(Ty* res, const Ty* srcP, const PointsSelection& posL, Int b_size, qlat::fft_desc_basic& fd)
{
  TIMERA("===write_grid_point_to_src===")
  /////if(pos.size() != 4){abort_r("dimension of positions wrong!\n ");}
  Int  NTt  = fd.Nv[3];
  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  LInt total = 6*NTt*Nxyz;
  LInt bfac = total/(b_size);
  const Int  Ns = 12;
  if(total % b_size != 0){abort_r("eigen system configurations wrong! \n");}

  Ty phase = 0.0;
  for(Long pi=0;pi<Long(posL.size());pi++){
  const Coordinate& pos = posL[pi];
  if(fd.coordinate_g_is_local(pos)){
    LInt isp = fd.index_l_from_g_coordinate(pos);
    phase = srcP[isp];
    ////printf("src pos %d %d %d %d, real %.3f imag %.3f \n", pos[0], pos[1], pos[2], pos[3], phase.real(), phase.imag());
    qacc_forNB(d0, 12, {
      Int d1 = d0;
      Long mi = d1*NTt*Nxyz + isp;
      Int chi = mi/(total);
      LInt xi = mi%(total);
      Long bi = xi/b_size;
      Long bj = xi%b_size;
      Long off  = (chi*bfac+bi)*Ns*b_size  + d0*b_size + bj;
      res[off] += phase;
    });
    /////phaseN = qlat::qnorm(src[isp]);
  }
  }
  qacc_barrier(dummy);
}



template<typename Ty>
void write_grid_point_to_src(Ty* res, const qnoiT& src, const PointsSelection& posL, Int b_size, qlat::fft_desc_basic& fd)
{
  const Ty* srcP = (Ty*) get_data(src).data();
  write_grid_point_to_src(res, srcP, posL, b_size, fd);
}

template<typename Ty>
void write_grid_point_to_src(Ty* res, const qnoiT& src, const Coordinate& pos, Int b_size, qlat::fft_desc_basic& fd)
{
  const Coordinate total_site_fake;
  PointsSelection posL;posL.init(total_site_fake, 1);posL[0] = pos;
  write_grid_point_to_src(res, src, posL, b_size, fd);
}

void print_psel(PointsSelection& psel)
{
  for(Long i=0;i<psel.size();i++){
    Coordinate& xg0 = psel[i];
    qmessage("x %d, y %d, z %d, t %d\n", xg0[0], xg0[1], xg0[2], xg0[3]);
  }
}

void add_psel(PointsSelection& p0, const PointsSelection& p1)
{
  const Long p0_size = p0.size();
  const Long p1_size = p1.size();
  p0.resize(p0_size + p1_size);
  for (unsigned int i = 0; i < p1.size(); i++) {
    p0[p0_size + i] = p1[i];
  }
}

void vector_to_Coordinate(qlat::vector<Int >& nv, Coordinate& pos, Int dir = 1)
{
  if(dir == 1){Qassert(nv.size() != 4);}
  if(dir == 1){for(Int i=0;i<4;i++){pos[i] = nv[i] ;}}

  if(dir == 0){nv.resize(4);}
  if(dir == 0){for(Int i=0;i<4;i++){nv[i]  = pos[i];}}
}

inline void get_grid_psel(PointsSelection& psel, const Coordinate& nv, const Coordinate& grid, qlat::RngState& rs, Int t0 = -1, const Int even = -1, const Coordinate& ini_ = Coordinate(-1,-1,-1,-1))
{
  /////put seed to all the same as rank 0
  //if(qlat::get_id_node() != 0){seed = 0;}
  //sum_all_size(&seed, 1 );
  //qlat::RngState rs(seed);

  psel.init();

  // Long total = 1;
  for (Int i = 0; i < 4; i++) {
    if (nv[i] < 0 or grid[i] < 0 or nv[i] % grid[i] != 0) {
      qmessage("Grid offset wrong nv[i] %d, grid[i] %d !\n", nv[i], grid[i]);
    }
    // total *= grid[i];
  }

  Coordinate ini;
  if(ini_ != Coordinate(-1,-1,-1,-1)){ini = ini_;}
  else{
    for(Int i=0;i<4;i++){ini[i] = int(qlat::u_rand_gen(rs)*(nv[i]/grid[i]));}
  }
  if(t0 != -1){ini[3] = t0;}
  if(even != -1){
    //// even = (z*2+y)*2+x;
    Int v[3];
    v[2] =  even/4;
    v[1] = (even%4)/2;
    v[0] =  even%2;
    for(Int si=0;si<3;si++)
    {
      Qassert((nv[si]/grid[si]) % 2 == 0);////stagger 8 eo requirements
      if(ini[si]%2 != v[si]){
        ini[si] = (ini[si] + 1 ) % (nv[si]/grid[si]);
      }
    }
  }


  std::vector<Coordinate> psel_xgs;

  for(Int xi=0;xi<grid[0];xi++)
  for(Int yi=0;yi<grid[1];yi++)
  for(Int zi=0;zi<grid[2];zi++)
  for(Int ti=0;ti<grid[3];ti++)
  {
    Coordinate xg;
    Coordinate ci;
    ci[0] = xi;  ci[1] = yi; ci[2] = zi; ci[3] = ti;

    for(Int i=0;i<4;i++){xg[i] = (ini[i] + ci[i]*(nv[i]/grid[i])) % (nv[i]);}

    psel_xgs.push_back(xg);
  }

  const Coordinate total_site_fake;
  psel.init(total_site_fake, psel_xgs);

}

template <typename Ty>
void get_noises_Coordinate(const qlat::FieldM<Ty, 1>& noise, PointsSelection& psel, Int printv = 0)
{
  const qlat::Geometry& geo = noise.geo();
  qlat::vector<Int > nv,Nv,mv;
  geo_to_nv(geo, nv, Nv, mv);
  //int nx,ny,nz,nt;
  LInt Nsite = Nv[0]*Nv[1]*Nv[2]*Nv[3];

  ////PointSelection local_tem;
  std::vector<Int > grid_pos;const Int DIM = 4;

  for(LInt isp=0; isp< Nsite; isp++)
  {
    if(qnorm(noise.get_elem_offset(isp)) > 0.01)
    {
      Coordinate xl0 = geo.coordinate_from_index(isp);
      Coordinate xg0 = geo.coordinate_g_from_l(xl0);
      for(Int i=0;i<DIM;i++){grid_pos.push_back(xg0[i]);}
      if(printv == 2){
        Ty tem = noise.get_elem_offset(isp);
        printf("rank %d, x %d, y %d, z %d, t %d, value %.3e %.3e \n",
            qlat::get_id_node() , xg0[0], xg0[1], xg0[2], xg0[3], tem.real(), tem.imag());
     }
    }
  }

  std::vector<Int > grid_pos_global = sum_local_to_global_vector(grid_pos);
  /////for(unsigned int i=0;i< grid_pos_global.size();i++){qmessage("i %d %d \n",i ,grid_pos_global[i]);}

  Long total = grid_pos_global.size()/DIM;

  psel.resize(total);Coordinate tem;
  for(Long p =0;p < total;p++){
    for(Int i=0;i < DIM; i++ ){tem[i] = grid_pos_global[p*4 + i];}
    psel[p] = tem;
  }

  if(printv >= 1){printf("rank %d, number of non-zeros %ld \n", qlat::get_id_node(), total);}

}

template <class Ty, Int civ>
void get_mix_color_src(qlat::FieldM<Ty , civ>& src, const Coordinate& sp, 
  const std::vector<double >& phases, const FieldSelection& fsel, const Int type_src = 0, Int seed = 0, const Int offT = -1, const Coordinate& offG = Coordinate(1,1,1,1))
{
  TIMERA("get_mix_color_src");
  Qassert(src.initialized);
  const qlat::Geometry& geo = src.geo();
  const Long V_local = geo.local_volume();

  Ty* srcP = (Ty*) qlat::get_data(src).data();
  zero_Ty(srcP, V_local*civ, 0);

  qlat::vector<Ty > color_phases(phases.size());
  Qassert(color_phases.size() >= civ);
  const Int tsrc = sp[3];
  for(unsigned int c=0;c<color_phases.size();c++){
    double r = phases[c];
    color_phases[c] = Ty(std::cos(r), std::sin(r));
  }
  //fft_desc_basic fd(geo);
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);

  if(type_src <= -1) ////point src, with only color zero
  {
    if(fd.coordinate_g_is_local(sp)){
      LInt isp = fd.index_l_from_g_coordinate(sp);
      for(Int c=0;c<civ;c++){
        if(c == -1 * (type_src + 1)){
          srcP[isp*civ + c] = color_phases[c];
        }
      }
    }
  }

  if(type_src == 0) ////point src
  {
    if(fd.coordinate_g_is_local(sp)){
      LInt isp = fd.index_l_from_g_coordinate(sp);
      for(Int c=0;c<civ;c++){
        srcP[isp*civ + c] = color_phases[c];
      }
    }
  }

  if(type_src == 1) ////Wall src
  {
    std::vector<qlat::RngState > rsL;rsL.resize(omp_get_max_threads());
    for(Int is=0;is<omp_get_max_threads();is++)
    {
      rsL[is] = qlat::RngState(seed + qlat::get_id_node()*omp_get_max_threads() + is);
    }
    
    qthread_for(isp, geo.local_volume(), {
      Coordinate xl = geo.coordinate_from_index(isp);
      Coordinate xg = geo.coordinate_g_from_l(xl);
      if(xg[3] == tsrc){
        qlat::RngState& rs = rsL[omp_get_thread_num()];
        for(Int c=0;c<civ;c++){
          double r = 2 * PI * qlat::u_rand_gen(rs);
          srcP[isp*civ + c] = Ty(std::cos(r), std::sin(r));
        }
      }
    });
  }

  if(type_src == 11) ////Wall src tests
  {
    qacc_for(isp, geo.local_volume(), {
      Coordinate xl = geo.coordinate_from_index(isp);
      Coordinate xg = geo.coordinate_g_from_l(xl);
      if(xg[3] == tsrc)
      for(Int c=0;c<civ;c++){
        srcP[isp*civ + c] = color_phases[c];
      }
    });
  }

  if(type_src == 12 or type_src == 13 or type_src == 14) ////Wall src, even or odd
  {
    std::vector<qlat::RngState > rsL;rsL.resize(omp_get_max_threads());
    for(Int is=0;is<omp_get_max_threads();is++)
    {
      rsL[is] = qlat::RngState(seed + qlat::get_id_node()*omp_get_max_threads() + is);
    }

    Int seed_g = seed;if(qlat::get_id_node() != 0){seed_g = 0;}
    sum_all_size(&seed_g, 1);
    qlat::RngState rs = qlat::RngState(seed_g + 127482);///global random
    Int src_eo = int(qlat::u_rand_gen(rs) * 2);
    if(type_src == 13){src_eo = 1;}
    if(type_src == 14){src_eo = 0;}
    
    qthread_for(isp, geo.local_volume(), {
      Coordinate xl = geo.coordinate_from_index(isp);
      Coordinate xg = geo.coordinate_g_from_l(xl);
      Int site_eo = (xl[0] + xl[1] + xl[2] + xl[3]) % 2;
      if(xg[3] == tsrc and site_eo == src_eo){
        qlat::RngState& rs = rsL[omp_get_thread_num()];
        for(Int c=0;c<civ;c++){
          double r = 2 * PI * qlat::u_rand_gen(rs);
          srcP[isp*civ + c] = Ty(std::cos(r), std::sin(r));
        }
      }
    });
  }

  if(type_src == 2) ////sparse src
  {
    std::vector<qlat::RngState > rsL;rsL.resize(omp_get_max_threads());
    for(Int is=0;is<omp_get_max_threads();is++)
    {
      rsL[is] = qlat::RngState(seed + qlat::get_id_node()*omp_get_max_threads() + is);
    }
    
    qthread_for(isp, geo.local_volume(), {
      const Long rank = fsel.f_local_idx.get_elem_offset(isp);
      if(rank >= 0){
        const Coordinate xl  = geo.coordinate_from_index(isp);
        const Coordinate xg  = geo.coordinate_g_from_l(xl);
        if(xg[3] == tsrc){
          qlat::RngState& rs = rsL[omp_get_thread_num()];
          for(Int c=0;c<civ;c++){
            double r = 2 * PI * qlat::u_rand_gen(rs);
            srcP[isp*civ + c] = Ty(std::cos(r), std::sin(r));
          }
        }
      }
    }); 
  }

  if(type_src == 3) ////grid src
  {
    std::vector<qlat::RngState > rsL;rsL.resize(omp_get_max_threads());
    for(Int is=0;is<omp_get_max_threads();is++)
    {
      rsL[is] = qlat::RngState(seed + qlat::get_id_node()*omp_get_max_threads() + is);
    }

    qlat::vector<Int > nv,Nv,mv;
    geo_to_nv(geo, nv, Nv, mv);
    for(Int i=0;i<4;i++){Qassert(nv[i] % offG[i] == 0);}

    qthread_for(isp, geo.local_volume(), {
      const Coordinate xl  = geo.coordinate_from_index(isp);
      const Coordinate xg  = geo.coordinate_g_from_l(xl);
      Int found = 0;
      for(Int i=0;i<4;i++){
        if((xg[i]-sp[i]+nv[i])%offG[i] == 0){
          found += 1;
        }
      };

      if(found == 4)
      {
        qlat::RngState& rs = rsL[omp_get_thread_num()];
        for(Int c=0;c<civ;c++){
          double r = 2 * PI * qlat::u_rand_gen(rs);
          srcP[isp*civ + c] = Ty(std::cos(r), std::sin(r));
        }
      }
    }); 
  }


  if(type_src == 20) ////T grid src, all spatial the same for momenta projections
  {
    Qassert(offT > 0);Qassert(fd.nt % offT == 0);
    Qassert(Long(color_phases.size()) == fd.nt/offT * civ);
    /////qlat::RngState rs = qlat::RngState(seed + type_src*10 + qlat::get_id_node() * 5);
    Coordinate tem = sp;
    for(Int ti = 0; ti < fd.nt/offT; ti ++)
    {
      tem[3] = (sp[3] + ti * offT) % fd.nt;
      if(fd.coordinate_g_is_local(tem)){
        LInt isp = fd.index_l_from_g_coordinate(tem);
        for(Int c=0;c<civ;c++){
          ////double r = 2 * PI * qlat::u_rand_gen(rs);
          ////srcP[isp*civ + c] = Ty(std::cos(r), std::sin(r));;
          srcP[isp*civ + c] = color_phases[ti*civ + c];
        }
      }
    }
  }


}

template <class Ty>
void vec_apply_cut(qlat::vector_gpu<Ty >& res, const Coordinate& sp, const double rmax, const Geometry& geo)
{
  TIMER("vec_apply_cut");
  if(rmax < 0 ){return ;}
  //const Long V_local = geo.local_volume();
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);

  qlat::vector<Int > nv, mv, Nv;
  geo_to_nv(geo, nv, Nv, mv);

  const Int Nxyz = fd.Nx * fd.Ny * fd.Nz;
  const Int Nt   = fd.Nt;
  const double rmax2 = rmax * rmax;
  const Int civ = res.size() / (Nt*Nxyz);

  Ty* srcP = (Ty*) qlat::get_data(res).data();
  qacc_for(xi, Nxyz, {
    Coordinate xl = geo.coordinate_from_index(xi);
    Coordinate xg = geo.coordinate_g_from_l(xl);

    double dis = 0.0;
    for(Int i=0;i<3;i++){
      double tmp = (xg[i] - sp[i] + nv[i])%(nv[i]/2);
      dis += tmp * tmp;
    }

    if(dis > rmax2)
    for(Int ci=0;ci<civ*Nt;ci++)
    {
      srcP[ci*Nxyz + xi] = 0.0;
    }
  });

}

template <class Tr, class Ty, Int civ>
void get_point_color_src(std::vector<qlat::FieldM<Tr , civ> >& srcL, 
  const PointsSelection& grids, const std::vector<Ty >& phases)
{
  TIMER("get_point_color_src");
  Qassert(srcL.size() == civ);
  Qassert(srcL[0].initialized);
  const qlat::Geometry& geo = srcL[0].geo();
  const Long V_local = geo.local_volume();

  std::vector<Tr* > srcP;srcP.resize(srcL.size());
  for(Int ic=0;ic<srcL.size();ic++){
    Qassert(srcL[ic].initialized);
    srcP[ic] = (Tr*) qlat::get_data(srcL[ic]).data();
    zero_Ty(srcP[ic], V_local*civ, 0);
  }

  const fft_desc_basic& fd = get_fft_desc_basic_plan(geo);
  Qassert(grids.size() == phases.size());
  for(unsigned int gi=0;gi<grids.size();gi++){
    const Coordinate& sp = grids[gi];
    if(fd.coordinate_g_is_local(sp)){
      LInt isp = fd.index_l_from_g_coordinate(sp);
      for(Int c=0;c<civ;c++){
        srcP[c][isp*civ + c] = phases[gi];
      }
    }
  }

}

template <class Tr, Int civ>
void get_wall_color_src(std::vector<qlat::FieldM<Tr , civ> >& srcL, const Coordinate& ps, const Int src_eo = -1)
{
  TIMER("get_wall_color_src");
  Qassert(srcL.size() == civ);
  Qassert(srcL[0].initialized);
  const qlat::Geometry& geo = srcL[0].geo();
  const Long V_local = geo.local_volume();
  const Int tini = ps[3];

  std::vector<Tr* > srcP;srcP.resize(srcL.size());
  for(Int ic=0;ic<srcL.size();ic++){
    Qassert(srcL[ic].initialized);
    srcP[ic] = (Tr*) qlat::get_data(srcL[ic]).data();
    zero_Ty(srcP[ic], V_local*civ, 0);
  }

  qthread_for(isp, geo.local_volume(), {
    Coordinate xl   = geo.coordinate_from_index(isp);
    Coordinate p    = geo.coordinate_g_from_l(xl);
    //int site_eo = (xl[0] + xl[1] + xl[2] + xl[3]) % 2;
    Int site_eo = ( (xl[2]%2) * 2 + xl[1]%2 ) * 2 + xl[0]%2;
    if(src_eo == -1){site_eo = src_eo;}
    if(p[3] == tini and site_eo == src_eo)
    {
      for(Int c=0;c<civ;c++){
        srcP[c][isp*civ + c] = Tr(1.0, 0.0);
      }
    }
  });
}

template <class T>
void make_point_prop(Propagator4dT<T>& prop, const Coordinate& sp = Coordinate(0, 0, 0, 0))
{
  const qlat::Geometry& geo = prop.geo();
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);

  const Int civ = 12 * 12;
  qlat::ComplexT<T >* propP = (qlat::ComplexT<T >*) qlat::get_data(prop).data();
  zero_Ty(propP, geo.local_volume()*civ, 0);

  if(fd.coordinate_g_is_local(sp))
  {
    /////Long isp = 0;
    LInt isp = fd.index_l_from_g_coordinate(sp);
    Coordinate xl   = geo.coordinate_from_index(isp);
    Coordinate p    = geo.coordinate_g_from_l(xl);
    if(p[0] == sp[0] and p[1] == sp[1] and p[2] == sp[2] and p[3] == sp[3])
    {
      qlat::WilsonMatrixT<double >& p1 =  prop.get_elem_offset(isp);
      for(Int dc0 =0;dc0<12;dc0++)
      {
        p1(dc0, dc0) = 1.0;
      }
    }
  }
}

template <class Td>
void make_grid_src(Propagator4dT<Td >& src, const Coordinate& sp, const Coordinate& offG = Coordinate(1,1,1,1),  Int seed = 0)
{
  TIMERA("make_grid_src");
  Qassert(src.initialized);
  const qlat::Geometry& geo = src.geo();
  const Long V_local = geo.local_volume();
  const Int civ = 12 * 12;

  qlat::ComplexT<Td >* srcP = (qlat::ComplexT<Td >*) qlat::get_data(src).data();
  zero_Ty(srcP, V_local*civ, 0);

  // fft_desc_basic fd(geo);
  // fft_desc_basic& fd = get_fft_desc_basic_plan(geo);

  std::vector<qlat::RngState > rsL;rsL.resize(omp_get_max_threads());
  for(Int is=0;is<omp_get_max_threads();is++)
  {
    rsL[is] = qlat::RngState(seed + qlat::get_id_node()*omp_get_max_threads() + is);
  }

  qlat::vector<Int > nv,Nv,mv;
  geo_to_nv(geo, nv, Nv, mv);
  for(Int i=0;i<4;i++){Qassert(nv[i] % offG[i] == 0);}
  qmessage("===grid numbers %8d \n", (nv[0]*nv[1]*nv[2]*nv[3])/(offG[0]*offG[1]*offG[2]*offG[3]) );

  qthread_for(isp, geo.local_volume(), {
    const Coordinate xl  = geo.coordinate_from_index(isp);
    const Coordinate xg  = geo.coordinate_g_from_l(xl);
    Int found = 0;
    for(Int i=0;i<4;i++){
      if((xg[i]-sp[i]+nv[i])%offG[i] == 0){
        found += 1;
      }
    };

    if(found == 4)
    {
      qlat::RngState& rs = rsL[omp_get_thread_num()];
      double r = 2 * PI * qlat::u_rand_gen(rs);
      //for(Int dc=0;dc<12;dc++)
      for(Int dc=0;dc<12;dc++)
      {
        srcP[(isp*12+dc)*12 + dc] = qlat::ComplexT<Td >(std::cos(r), std::sin(r));
      }
    }
  }); 

}


//seed to be -1 for all 1 wall src
template <class Td>
void make_volume_src(Propagator4dT<Td >& src, Int seed = 0, Int mix_color = 0, Int mix_spin = 0, Int tini = -1)
{
  TIMERA("make_volume_src");
  Qassert(src.initialized);
  const qlat::Geometry& geo = src.geo();
  const Long V_local = geo.local_volume();
  const Int civ = 12 * 12;

  qlat::ComplexT<Td >* srcP = (qlat::ComplexT<Td >*) qlat::get_data(src).data();
  zero_Ty(srcP, V_local*civ, 0);
  ////fft_desc_basic fd(geo);

  const Int Nthread = omp_get_max_threads();
  std::vector<qlat::RngState > rsL;rsL.resize(Nthread);
  for(Int is=0;is<Nthread;is++)
  {
    rsL[is] = qlat::RngState(seed + qlat::get_id_node()*Nthread + is);
  }

  qthread_for(isp, geo.local_volume(), {
    qlat::RngState& rs = rsL[omp_get_thread_num()];
    Coordinate xl   = geo.coordinate_from_index(isp);
    Coordinate p    = geo.coordinate_g_from_l(xl);
    if(tini == -1 or p[3] == tini)
    {
      if(seed == -1) 
      {   
        for(Int dc=0;dc<12;dc++)
        {   
          srcP[(isp*12+dc)*12 + dc] = qlat::ComplexT<Td >(1.0, 0.0);
        }   
      }else{   

        if(mix_color == 0 and mix_spin == 0)
        {
          ////double r = (2 * PI /3.0 ) * int(3 * qlat::u_rand_gen(rs));
          double r = 2 * PI * qlat::u_rand_gen(rs);
          for(Int dc=0;dc<12;dc++)
          {
            srcP[(isp*12+dc)*12 + dc] = qlat::ComplexT<Td >(std::cos(r), std::sin(r));
          }
        }

        if(mix_color == 1 and mix_spin == 0)
        {
          for(Int c0=0;c0<3;c0++)
          for(Int c1=0;c1<3;c1++)
          {
            ////double r = (2 * PI /3.0 ) * int(3 * qlat::u_rand_gen(rs));
            double r = 2 * PI * qlat::u_rand_gen(rs);
            for(Int d =0;d<4;d++)
            {
              srcP[(isp*12+d*3 + c0)*12 + d*3+c1] = qlat::ComplexT<Td >(std::cos(r), std::sin(r)) / qlat::ComplexT<Td >(3.0, 0.0);
            }
          }
        }

        if(mix_color == 1 and mix_spin == 1)
        {
          for(Int d0 =0;d0<2;d0++)
          for(Int d1 =0;d1<2;d1++)
          for(Int c0=0;c0<3;c0++)
          for(Int c1=0;c1<3;c1++)
          {
            double r = 2 * PI * qlat::u_rand_gen(rs);
            srcP[(isp*12+d0*3 + c0)*12 + d1*3+c1] = qlat::ComplexT<Td >(std::cos(r), std::sin(r))  / qlat::ComplexT<Td >(6.0, 0.0);
          }

          for(Int d0 =2;d0<4;d0++)
          for(Int d1 =2;d1<4;d1++)
          for(Int c0=0;c0<3;c0++)
          for(Int c1=0;c1<3;c1++)
          {
            double r = 2 * PI * qlat::u_rand_gen(rs);
            srcP[(isp*12+d0*3 + c0)*12 + d1*3+c1] = qlat::ComplexT<Td >(std::cos(r), std::sin(r))  / qlat::ComplexT<Td >(6.0, 0.0);
          }

        }
      }
    }
  }); 
}

inline void split_points_simple(std::vector<std::vector<Coordinate > >& res, std::vector<Coordinate >& grids){
  const Long Ng = grids.size();
  res.resize(Ng);
  for(Long gi=0;gi<Ng;gi++){
    res[gi].resize(1);
    res[gi][0] = grids[gi];
  }
}

/*
  sort spatial grid points into sub-grids if necessary
*/
inline void sort_grid_points(std::vector<Coordinate >& res, std::vector<Coordinate >& grids, const Coordinate& Lat, const Int num_points = 0){
  const Long Ng = grids.size();
  const Long max_g = num_points * num_points * num_points;
  res.resize(Ng);

  // do all points if needed points larger than numbers
  if(num_points == 0 or max_g >= Ng){
    for(Long gi=0;gi<Ng;gi++){
      res[gi] = grids[gi];
    }
    return ;
  }

  // Need points to be a grid
  Qassert(Ng % max_g == 0);
  const Long Ngroup = Ng / max_g;
  std::vector<Int > Loff;Loff.resize(3);
  for(Int i=0;i<3;i++){
    Qassert(Lat[i] % num_points == 0);
    Loff[i] = Lat[i] / num_points;
  }

  std::vector<Long > select;select.resize(Ng);
  for(Long gi=0;gi<Ng;gi++){
    select[gi] = Coordinate_to_index(grids[gi], Lat);
  }
  Long count = 0;
  for(Long gi=0;gi<Ngroup;gi++){
    // get the initial points
    Coordinate cur;bool findc = false;
    for(Long gj=0;gj<Ng;gj++){
      if(select[gj] >= 0){
        cur = grids[gj];
        findc = true;
        break;
      }
    }
    Qassert(findc == true);
    for(Int iz=0;iz<num_points;iz++)
    for(Int iy=0;iy<num_points;iy++)
    for(Int ix=0;ix<num_points;ix++)
    {
      Coordinate sp = cur;
      sp[0] = (sp[0] + ix * Loff[0]) % Lat[0];
      sp[1] = (sp[1] + iy * Loff[1]) % Lat[1];
      sp[2] = (sp[2] + iz * Loff[2]) % Lat[2];
      //std::string sm = Coordinate_to_string(sp);
      //qmessage("count %5d %s\n", count, sm.c_str());

      const Long idx = Coordinate_to_index(sp, Lat);
      const Int find = std_find(select, idx);
      Qassert(find >= 0 and find < Ng);
      select[find] = -1;
      //res[count] = sp;
      res[count] = grids[find];
      //for(Int m=0;m<4;m++){res[count][m] = grids[find][m];}
      count += 1;
    }
  }
}

/*
  sort points with grids
    64 points with 8 point in grid 
  num_points : 
    0   simple loop over points
    > 1 groups into subgrid
*/
inline void split_grid_points(std::vector<std::vector<Coordinate > >& res, std::vector<Coordinate >& grids, const Coordinate& Lat, const Int combineT = 0, const Int num_points = 0){
  const Long Ng = grids.size();
  // get all time slice
  std::vector<Int > tsrc;
  std::vector<std::vector<Long > > idx_list;
  std::vector<std::vector<Coordinate > > spx_list;
  for(Long gi=0;gi<Ng;gi++){
    const Int tcur = grids[gi][3];
    if(std_find(tsrc, tcur) == -1){tsrc.push_back(tcur);}
  }
  std::sort(tsrc.begin(), tsrc.end());
  const Long Nt = tsrc.size();
  const Long Ntsrc = combineT == 0 ? 1 : Nt;
  const Long max_g = num_points * num_points * num_points;

  bool do_simple = false;
  if(num_points == 0 and Ntsrc == 1){do_simple = true;}
  if(max_g == Ng and Ntsrc == 1){do_simple = true;}
  if(do_simple){split_points_simple(res, grids);return ;}

  if(combineT != 0){
    Qassert(Ng % Nt == 0);
  }// each time slice much have same number of points
  idx_list.resize(Nt);
  spx_list.resize(Nt);
  for(Long gi=0;gi<Ng;gi++){
    Coordinate sp = grids[gi];
    const Long ti = std_find(tsrc, sp[3]);
    Qassert(ti >= 0 and ti < Nt);
    sp[3] = 0; // all points to zero to compare
    idx_list[ti].push_back(Coordinate_to_index(sp, Lat));
    spx_list[ti].push_back(grids[gi]);
  }

  for(Long ti=0;ti<Nt;ti++){
    std::sort(idx_list[ti].begin(), idx_list[ti].end());
  }
  const Int Nx_grid = idx_list[0].size();

  bool same_points_time = true;
  for(Long ti=1;ti<Nt;ti++){
    Qassert(int(idx_list[ti].size()) == Nx_grid);
    if(same_points_time)
    for(Int i=0;i<Nx_grid;i++){
      if(idx_list[ti][i] != idx_list[0][i]){
        same_points_time = false;
      }
    }
  }

  const Long Nr = Ng / Ntsrc;
  res.resize(Nr);

  std::vector<std::vector<Coordinate > > sort_grids;sort_grids.resize(Nt);
  for(Long ti=0;ti<Nt;ti++){
    sort_grid_points( sort_grids[ti], spx_list[ti], Lat, num_points);
  }

  if(Ntsrc > 1){
    Qassert(Ntsrc == Nt and Nr == Long(sort_grids[0].size()));
    for(Long gi=0;gi<Ng;gi++){
      if(gi >= Nr){continue; }

      res[gi].resize(Nt);
      for(Long ti=0;ti<Nt;ti++){
        res[gi][ti] = sort_grids[ti][gi]; // currently not trying to have max distance
      }

      if(!same_points_time){
        res[gi].resize(Nt);
        for(Long ti=0;ti<Nt;ti++){
          res[gi][ti] = sort_grids[ti][gi]; // currently not trying to have max distance
        }
      }
      if( same_points_time){
        res[gi].resize(Nt);
        Coordinate cur = sort_grids[0][gi];
        for(Long ti=0;ti<Nt;ti++){
          cur[3] = tsrc[ti];
          res[gi][ti] = cur;  // only 0 point will be used for shifts
        }
      }
    }
  }

  if(Ntsrc == 1){
    const Long goffT = Ng / Nt;
    for(Long gi=0;gi<Ng;gi++){
      if(gi > Nr){continue; }
      const Long ti = gi / goffT;
      const Long oi = gi % goffT;
      res[gi].resize(1);
      res[gi][0] = sort_grids[ti][oi]; // currently not trying to have max distance
    }
  }

}

template <typename Td>
void check_prop_noise_positions(qlat::Propagator4dT<Td >& prop, Coordinate& pos, Coordinate&off_L,Int printS=0,Int mod=0)
{
  qlat::Geometry& geo = prop.geo();
  using Ty = ComplexT<Td>;
  qlat::FieldM<Ty, 1> noise;noise.init(geo);
  const Long Nvol = geo.local_volume();
  Ty* pD = (Ty*) get_data(prop ).data();
  Ty* nD = (Ty*) get_data(noise).data();
  qacc_for(isp, Nvol,{
    nD[isp] = pD[isp*12*12 + 0];
  });
  check_noise_pos(noise, pos, off_L, printS, mod);
}


}


#endif

