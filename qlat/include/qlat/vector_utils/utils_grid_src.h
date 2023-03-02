// utils_lms_funs.h
// Gen Wang
// Oct. 2021

#ifndef UTILS_GRID_SRC_H
#define UTILS_GRID_SRC_H

#pragma once

#include "utils_float_type.h"
#include "utils_gammas.h"
#include "utils_momentum.h"
#include "utils_fft_desc.h"

namespace qlat{

template <typename Ty>
void check_noise_pos(qlat::FieldM<Ty, 1>& noise, Coordinate& pos, Coordinate&off_L,int printS=0,int mod=0)
{
  qlat::Geometry& geo = noise.geo();
  qlat::vector_acc<int > nv,Nv,mv;
  geo_to_nv(geo, nv, Nv, mv);
  int nx,ny,nz,nt;
  nx = nv[0];ny = nv[1];nz = nv[2];nt = nv[3];
  LInt Nsite = Nv[0]*Nv[1]*Nv[2]*Nv[3];

  std::vector<int > NL(4);NL[0]=nx;NL[1]=ny;NL[2]=nz;NL[3]=nt;
  int grid_count = 0;
  std::vector<std::vector<int > > grid;
  for(int iL=0;iL<4;iL++){
    grid.push_back(std::vector<int>(NL[iL]));
    for(LInt giL=0;giL<grid[iL].size();giL++){grid[iL][giL] = 0.0;}
  }
  //grid.push_back(std::vector<double > (nx));
  int number_t = 1;int t_ini = 0;
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
        for(int i=0;i<4;i++){grid[i][xg0[i]] += 1;}
        ///grid[0][p.x()] += 1.0;
        ///grid[1][p.y()] += 1.0;
        ///grid[2][p.z()] += 1.0;
        ///grid[3][p.t()] += 1.0;
        grid_count = grid_count + 1;
      }
    }
  }
  for(int iL=0;iL<4;iL++){sum_all_size(&grid[iL][0],NL[iL]);}
  sum_all_size(&grid_count,1);
  ////global_sum_all(&grid_count,1);
  ////off_L.resize(4);
  for(int oi=0;oi<4;oi++){off_L[oi] = 0;}
  for(int iL=0;iL<4;iL++)for(int k=0;k<NL[iL];k++)if(grid[iL][k]>0.0)off_L[iL] += 1;
  //for(int x=0;x<nx;x++){if(grid[0][x]>0.0)off_L[0] += 1;}
  if(int(grid_count) != off_L[0]*off_L[1]*off_L[2]*off_L[3])
  {
    print0("Source Check Failed grid_count %10d, offx %5d, offy %5d, offz %5d, offt %5d!\n",
          int(grid_count),off_L[0],off_L[1],off_L[2],off_L[3]);
    qassert(false);
    ////shutdown_machine();
    ////abort();
  }

  //int pos = 0;int t_ini = 0;
  ////pos.resize(4);
  for(int i=0;i<4;i++){
    pos[i] = 0;
    for(int x=0;x<nv[i];x++){if(grid[i][x] > 0){pos[i] += x;break;}}
  }
  t_ini = pos[3];
  //for(int x=0;x<nx;x++){if(grid[0][x]>0.0){pos += ((x*100)*100)*1000;break;}}
  //for(int y=0;y<ny;y++){if(grid[1][y]>0.0){pos += (y*100)*1000;break;}}
  //for(int z=0;z<nx;z++){if(grid[2][z]>0.0){pos += (z)*1000;break;}}
  //for(int t=0;t<nt;t++){if(grid[3][t]>0.0){pos += (t);t_ini = t;break;}}

  print0("Check T %5d %5d %5d %5d, offx %5d, offy %5d, offz %5d, offt %5d. \n",
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
        int printv = 1;
        for(int i=0;i<4;i++)
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
void check_noise_high(qlat::FieldM<Ty, 1>& noise, std::vector<int >& sinkt, double& factor)
{
  const qlat::Geometry& geo = noise.geo();
  ////const long Nvol = geo.local_volume();
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);

  const long Nxyz = fd.Nx*fd.Ny*fd.Nz;

  std::vector<double > count;
  count.resize(fd.nt);
  for(int i=0;i<fd.nt;i++){count[i] = 0;}

  qthread_for(ti, fd.Nt, {
    const int tg = ti + fd.init;

    for(long isp=0;isp<Nxyz;isp++)
    {
      const long off = ti * Nxyz + isp;
      {
        double tem_source =  qlat::qnorm(noise.get_elem_offset(off));
        if(qnorm(tem_source)> 1e-3)
        {
          count[tg] += 1;
        }
      }
    }
  });

  sum_all_size(count.data(), count.size());
  sinkt.resize(0);
  factor = -1;
  for(int t=0;t<fd.nt;t++){
    if(count[t] > 0){
      if(factor == -1){
        factor = count[t] / fd.vol;
      }else{
        qassert(factor == count[t] / fd.vol);
      }
      sinkt.push_back(t);
    }
  }

}

/////get positions by spatial setups
inline void grid_list_pos(const Coordinate& off_L,qlat::vector_acc<long > &Ngrid)
{
  TIMERA("===grid_list_pos===")
  if(off_L.size() != 4){abort_r("dimention of positions wrong!\n ");}
  int Ntot = off_L[0]*off_L[1]*off_L[2];
  double max_d = 1.0*off_L[0]*off_L[0] + 1.0*off_L[1]*off_L[1] + 1.0*off_L[2]*off_L[2] + 1.0;

  Ngrid.resize(Ntot);
  for(int j0=0;j0<Ntot;j0++)
  {
    double dis_large = 0;
    for(int i0=0;i0<Ntot;i0++)
    {
      int ix0= i0/(off_L[1]*off_L[2]);
      int iy0= (i0%(off_L[1]*off_L[2]))/off_L[2];
      int iz0= i0%off_L[2];

      double dis_min = max_d;
      for(int iv=0;iv<j0;iv++)
      {
        //std::vector<int > temv = map_index[i0];
        int isrc = Ngrid[iv];
        int ix= isrc/(off_L[1]*off_L[2]);
        int iy= (isrc%(off_L[1]*off_L[2]))/off_L[2];
        int iz= isrc%off_L[2];
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

inline Coordinate get_grid_off(long j0, const Coordinate& off_L, const Coordinate& pos_ini, const Coordinate& Lat)
{
  if(pos_ini.size() != 4 or off_L.size() != 4){abort_r("dimension of positions wrong!\n ");}
  //////std::vector<int > pos;pos.resize(4);
  Coordinate pos = pos_ini;
  ////for(int i=0;i<4;i++){pos[i] = pos_ini[i];}

  Coordinate off_pos = qlat::coordinate_from_index(j0, off_L);
  for(int i=0;i<4;i++){pos[i] += ((Lat[i]/(off_L[i]))*off_pos[i] ) %(Lat[i]);}

  //int it  = j0/(off_L[0]*off_L[1]*off_L[2]);
  //long i0 = j0%(off_L[0]*off_L[1]*off_L[2]);
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
inline void grid_list_posT(std::vector<PointSelection > &LMS_points, const Coordinate& off_L, const Coordinate& pos, const int combineT, const Coordinate& Lat)
{
  TIMERA("===grid_list_posT===")
  qlat::vector_acc<long > Nfull;
  /////get positions by spatial setups
  grid_list_pos(off_L, Nfull);

  Coordinate cur_pos;
  Coordinate cur_off;
  for(long gi=0;gi<Nfull.size();gi++){
    cur_pos = get_grid_off(Nfull[gi], off_L, pos, Lat);
    if(combineT == 0){
      for(int it = 0; it < off_L[3]; it++){
        cur_off = cur_pos;
        cur_off[3] += ((Lat[3]/(off_L[3]))*it)%Lat[3]; ////need be careful not to exceed boundaries
        PointSelection lms_res;lms_res.push_back(cur_off);
        LMS_points.push_back(lms_res);
      }
    }
    if(combineT == 1){
      PointSelection lms_res;
      for(int it = 0; it < off_L[3]; it++){
        cur_off = cur_pos;
        cur_off[3] += ((Lat[3]/(off_L[3]))*it)%Lat[3];
        lms_res.push_back(cur_off);
      }
      LMS_points.push_back(lms_res);
    }
  }
  if(combineT == int(0)){qassert(long(LMS_points.size()) == long(Nfull.size()*off_L[3]));}
  if(combineT == int(1)){qassert(long(LMS_points.size()) == long(Nfull.size()         ));}
}


//////assume res have been cleared
//////asuume res on GPU
template<typename Ty>
void write_grid_point_to_src(Ty* res, const qnoiT& src, const PointSelection& posL, int b_size, qlat::fft_desc_basic& fd)
{
  TIMERA("===write_grid_point_to_src===")
  /////if(pos.size() != 4){abort_r("dimension of positions wrong!\n ");}
  int  NTt  = fd.Nv[3];
  LInt Nxyz = fd.Nv[0]*fd.Nv[1]*fd.Nv[2];
  LInt total = 6*NTt*Nxyz;
  LInt bfac = total/(b_size);
  const int  Ns = 12;
  if(total % b_size != 0){abort_r("eigen system configurations wrong! \n");}

  ////Geometry geo;fd.get_geo(geo );
  ////Coordinate xl0 = geo.coordinate_from_index(isp);
  ////Coordinate xg0 = geo.coordinate_g_from_l(xl0);

  Ty phase = 0.0;
  for(long pi=0;pi<long(posL.size());pi++){
  const Coordinate& pos = posL[pi];
  if(fd.coordinate_g_is_local(pos)){
    LInt isp = fd.index_l_from_g_coordinate(pos);
    phase = src.get_elem_offset(isp);
    ////printf("src pos %d %d %d %d, real %.3f imag %.3f \n", pos[0], pos[1], pos[2], pos[3], phase.real(), phase.imag());
    qacc_forNB(d0, 12, {
      int d1 = d0;
      long mi = d1*NTt*Nxyz + isp;
      int chi = mi/(total);
      LInt xi = mi%(total);
      long bi = xi/b_size;
      long bj = xi%b_size;
      long off  = (chi*bfac+bi)*Ns*b_size  + d0*b_size + bj;
      res[off] += phase;
    });
    /////phaseN = qlat::qnorm(src[isp]);
  }
  }
  qacc_barrier(dummy);

  //double flag = qlat::qnorm(phase);
  //sum_all_size(&flag, 1);
  ///////printf("rank %d, flag %.3f \n", fd.rank, flag);
  //if(flag < 1e-15){abort_r("src position phase equal zero!\n ");}

}

template<typename Ty>
void write_grid_point_to_src(Ty* res, const qnoiT& src, const Coordinate& pos, int b_size, qlat::fft_desc_basic& fd)
{
  PointSelection posL;posL.resize(1);posL[0] = pos;
  write_grid_point_to_src(res, src, posL, b_size, fd);
}

void print_psel(PointSelection& psel)
{
  for(unsigned long i=0;i<psel.size();i++){
    Coordinate& xg0 = psel[i];
    print0("x %d, y %d, z %d, t %d\n", xg0[0], xg0[1], xg0[2], xg0[3]);
  }
}

void add_psel(PointSelection& p0, const PointSelection& p1)
{
  for(unsigned int i=0;i<p1.size();i++){p0.push_back(p1[i]);}
}

void vector_to_Coordinate(qlat::vector_acc<int >& nv, Coordinate& pos, int dir = 1)
{
  if(dir == 1){qassert(nv.size() != 4);}
  if(dir == 1){for(int i=0;i<4;i++){pos[i] = nv[i] ;}}

  if(dir == 0){nv.resize(4);}
  if(dir == 0){for(int i=0;i<4;i++){nv[i]  = pos[i];}}
}

inline void get_grid_psel(PointSelection& psel, Coordinate& nv, Coordinate& grid, int t0 = -1, int seed = 123, const int even = -1)
{
  /////put seed to all the same as rank 0
  if(qlat::get_id_node() != 0){seed = 0;}
  sum_all_size(&seed, 1 );

  qlat::RngState rs(seed);
  long total = 1;
  for(int i=0;i<4;i++){
   if(nv[i] < 0 or grid[i] < 0 or nv[i]%grid[i] != 0){print0("Grid offset wrong nv[i] %d, grid[i] %d !\n", nv[i], grid[i]);}
   total *= grid[i];
  }

  Coordinate ini;
  for(int i=0;i<4;i++){ini[i] = int(qlat::u_rand_gen(rs)*(nv[i]/grid[i]));}
  if(t0 != -1){ini[3] = t0;}
  if(even != -1){
    //// even = (z*2+y)*2+x;
    int v[3];
    v[2] =  even/4;
    v[1] = (even%4)/2;
    v[0] =  even%2;
    for(int si=0;si<3;si++)
    {
      qassert((nv[si]/grid[si]) % 2 == 0);////stagger 8 eo requirements
      if(ini[si]%2 != v[si]){
        ini[si] = (ini[si] + 1 ) % (nv[si]/grid[si]);
      }
    }
  }

  psel.resize(0);

  for(int xi=0;xi<grid[0];xi++)
  for(int yi=0;yi<grid[1];yi++)
  for(int zi=0;zi<grid[2];zi++)
  for(int ti=0;ti<grid[3];ti++)
  {
    Coordinate xg;
    Coordinate ci;
    ci[0] = xi;  ci[1] = yi; ci[2] = zi; ci[3] = ti;

    for(int i=0;i<4;i++){xg[i] = ini[i] + ci[i]*(nv[i]/grid[i]);}

    psel.push_back(xg);
  }


}

template <typename Ty>
void get_noises_Coordinate(const qlat::FieldM<Ty, 1>& noise, PointSelection& psel, int printv = 0)
{
  const qlat::Geometry& geo = noise.geo();
  qlat::vector_acc<int > nv,Nv,mv;
  geo_to_nv(geo, nv, Nv, mv);
  //int nx,ny,nz,nt;
  LInt Nsite = Nv[0]*Nv[1]*Nv[2]*Nv[3];

  ////PointSelection local_tem;
  std::vector<int > grid_pos;const int DIM = 4;

  for(LInt isp=0; isp< Nsite; isp++)
  {
    if(qnorm(noise.get_elem_offset(isp)) > 0.01)
    {
      Coordinate xl0 = geo.coordinate_from_index(isp);
      Coordinate xg0 = geo.coordinate_g_from_l(xl0);
      for(int i=0;i<DIM;i++){grid_pos.push_back(xg0[i]);}
      if(printv == 2){
        Ty tem = noise.get_elem_offset(isp);
        printf("rank %d, x %d, y %d, z %d, t %d, value %.3e %.3e \n",
            qlat::get_id_node() , xg0[0], xg0[1], xg0[2], xg0[3], tem.real(), tem.imag());
     }
    }
  }

  std::vector<int > grid_pos_global = sum_local_to_global_vector(grid_pos);
  /////for(unsigned int i=0;i< grid_pos_global.size();i++){print0("i %d %d \n",i ,grid_pos_global[i]);}

  long total = grid_pos_global.size()/DIM;

  psel.resize(total);Coordinate tem;
  for(long p =0;p < total;p++){
    for(int i=0;i < DIM; i++ ){tem[i] = grid_pos_global[p*4 + i];}
    psel[p] = tem;
  }

  if(printv >= 1){printf("rank %d, number of non-zeros %ld \n", qlat::get_id_node(), total);}

}

template <class Ty, int civ>
void get_mix_color_src(qlat::FieldM<Ty , civ>& src, const Coordinate& sp, 
  const std::vector<double >& phases, const FieldSelection& fsel, const int type_src = 0, int seed = 0, const int offT = -1)
{
  TIMER("get_mix_color_src");
  qassert(src.initialized);
  const qlat::Geometry& geo = src.geo();
  const long V_local = geo.local_volume();

  Ty* srcP = (Ty*) qlat::get_data(src).data();
  zero_Ty(srcP, V_local*civ, 0);

  qlat::vector_acc<Ty > color_phases(phases.size());
  const int tsrc = sp[3];
  for(unsigned int c=0;c<color_phases.size();c++){
    double r = phases[c];
    color_phases[c] = Ty(std::cos(r), std::sin(r));
  }
  fft_desc_basic fd(geo);

  if(type_src == 0) ////point src
  {
    if(fd.coordinate_g_is_local(sp)){
      LInt isp = fd.index_l_from_g_coordinate(sp);
      for(int c=0;c<civ;c++){
        srcP[isp*civ + c] = color_phases[c];
      }
    }
  }

  if(type_src == 1) ////Wall src
  {
    std::vector<qlat::RngState > rsL;rsL.resize(omp_get_max_threads());
    for(int is=0;is<omp_get_max_threads();is++)
    {
      rsL[is] = qlat::RngState(seed + qlat::get_id_node()*omp_get_max_threads() + is);
    }
    
    qthread_for(isp, geo.local_volume(), {
      Coordinate xl = geo.coordinate_from_index(isp);
      Coordinate xg = geo.coordinate_g_from_l(xl);
      if(xg[3] == tsrc){
        qlat::RngState& rs = rsL[omp_get_thread_num()];
        for(int c=0;c<civ;c++){
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
      for(int c=0;c<civ;c++){
        srcP[isp*civ + c] = color_phases[c];
      }
    });
  }

  if(type_src == 2) ////sparse src
  {
    std::vector<qlat::RngState > rsL;rsL.resize(omp_get_max_threads());
    for(int is=0;is<omp_get_max_threads();is++)
    {
      rsL[is] = qlat::RngState(seed + qlat::get_id_node()*omp_get_max_threads() + is);
    }
    
    qthread_for(isp, geo.local_volume(), {
      const long rank = fsel.f_local_idx.get_elem_offset(isp);
      if(rank >= 0){
        const Coordinate xl  = geo.coordinate_from_index(isp);
        const Coordinate xg  = geo.coordinate_g_from_l(xl);
        if(xg[3] == tsrc){
          qlat::RngState& rs = rsL[omp_get_thread_num()];
          for(int c=0;c<civ;c++){
            double r = 2 * PI * qlat::u_rand_gen(rs);
            srcP[isp*civ + c] = Ty(std::cos(r), std::sin(r));
          }
        }
      }
    }); 
  }

  if(type_src == 20) ////T grid src, all spatial the same for momenta projections
  {
    qassert(offT > 0);qassert(fd.nt % offT == 0);
    qassert(long(color_phases.size()) == fd.nt/offT * civ);
    /////qlat::RngState rs = qlat::RngState(seed + type_src*10 + qlat::get_id_node() * 5);
    Coordinate tem = sp;
    for(int ti = 0; ti < fd.nt/offT; ti ++)
    {
      tem[3] = (sp[3] + ti * offT) % fd.nt;
      if(fd.coordinate_g_is_local(tem)){
        LInt isp = fd.index_l_from_g_coordinate(tem);
        for(int c=0;c<civ;c++){
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
  const long V_local = geo.local_volume();
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);

  qlat::vector_acc<int > nv, mv, Nv;
  geo_to_nv(geo, nv, Nv, mv);

  const int Nxyz = fd.Nx * fd.Ny * fd.Nz;
  const int Nt   = fd.Nt;
  const double rmax2 = rmax * rmax;
  const int civ = res.size() / (Nt*Nxyz);

  Ty* srcP = (Ty*) qlat::get_data(res).data();
  qacc_for(xi, Nxyz, {
    Coordinate xl = geo.coordinate_from_index(xi);
    Coordinate xg = geo.coordinate_g_from_l(xl);

    double dis = 0.0;
    for(int i=0;i<3;i++){
      double tmp = (xg[i] - sp[i] + nv[i])%(nv[i]/2);
      dis += tmp * tmp;
    }

    if(dis > rmax2)
    for(int ci=0;ci<civ*Nt;ci++)
    {
      srcP[ci*Nxyz + xi] = 0.0;
    }
  });

}

template <class Ty, int civ>
void get_point_color_src(std::vector<qlat::FieldM<Ty , civ> >& srcL, 
  const PointSelection& grids, const std::vector<Ty >& phases)
{
  TIMER("get_point_color_src");
  qassert(srcL.size() == civ);
  qassert(srcL[0].initialized);
  const qlat::Geometry& geo = srcL[0].geo();
  const long V_local = geo.local_volume();

  std::vector<Ty* > srcP;srcP.resize(srcL.size());
  for(int ic=0;ic<srcL.size();ic++){
    qassert(srcL[ic].initialized);
    srcP[ic] = (Ty*) qlat::get_data(srcL[ic]).data();
    zero_Ty(srcP[ic], V_local*civ, 0);
  }

  const fft_desc_basic& fd = get_fft_desc_basic_plan(geo);
  qassert(grids.size() == phases.size());
  for(unsigned int gi=0;gi<grids.size();gi++){
    const Coordinate& sp = grids[gi];
    if(fd.coordinate_g_is_local(sp)){
      LInt isp = fd.index_l_from_g_coordinate(sp);
      for(int c=0;c<civ;c++){
        srcP[c][isp*civ + c] = phases[gi];
      }
    }
  }

}



}


#endif

