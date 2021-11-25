// utils_lms_funs.h
// Gen Wang
// Oct. 2021

#ifndef UTILS_GRID_SRC_H
#define UTILS_GRID_SRC_H

#pragma once

#include "float_type.h"
#include "gammas.h"
#include "utils_momentum.h"

namespace qlat{

template <typename Ty>
void check_noise_pos(qlat::FieldM<Ty, 1>& noise, Coordinate& pos,qlat::vector_acc<int > &off_L,int printS=0,int mod=0)
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
      auto tem_source =  noise.get_elem(isp);
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
  off_L.resize(4);
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
        auto tem_source =  noise.get_elem(isp);
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
        auto tem_source =  noise.get_elem(isp);
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

inline void grid_list_pos(qlat::vector_acc<int >& off_L,qlat::vector_acc<int > &Ngrid)
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

inline Coordinate get_grid_off(int i0, qlat::vector_acc<int >& off_L, Coordinate& pos_ini, qlat::fft_desc_basic& fd)
{
  if(pos_ini.size() != 4 or off_L.size() != 4){abort_r("dimension of positions wrong!\n ");}
  //////std::vector<int > pos;pos.resize(4);
  Coordinate pos = pos_ini;
  ////for(int i=0;i<4;i++){pos[i] = pos_ini[i];}

  int ix= i0/(off_L[1]*off_L[2]);
  int iy= (i0%(off_L[1]*off_L[2]))/off_L[2];
  int iz= i0%off_L[2];

  pos[0] += (fd.nx/(off_L[0]))*ix;
  pos[1] += (fd.ny/(off_L[1]))*iy;
  pos[2] += (fd.nz/(off_L[2]))*iz;
  return pos;
}

//////assume res have been cleared
template<typename Ty>
void write_grid_point_to_src(Ty* res, const qnoiT& src, const Coordinate& pos, int b_size, qlat::fft_desc_basic& fd)
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
  if(fd.coordinate_g_is_local(pos)){
    LInt isp = fd.index_l_from_g_coordinate(pos);
    phase = src.get_elem(isp);
    ////print0("position %d %d %d %d !\n", cur_pos[0], cur_pos[1], cur_pos[2], cur_pos[3]);
    ////printf("src pos %d %d %d %d, real %.3f imag %.3f \n", pos[0], pos[1], pos[2], pos[3], phase.real(), phase.imag());
    for(int d0=0;d0<12;d0++){
      int d1 = d0;

      long mi = d1*NTt*Nxyz + isp;
      int chi = mi/(total);
      LInt xi = mi%(total);
      long bi = xi/b_size;
      long bj = xi%b_size;

      long off  = (chi*bfac+bi)*Ns*b_size  + d0*b_size + bj;
      res[off] += phase;
    }
    /////phaseN = qlat::qnorm(src[isp]);
  }

  //double flag = qlat::qnorm(phase);
  //sum_all_size(&flag, 1);
  ///////printf("rank %d, flag %.3f \n", fd.rank, flag);
  //if(flag < 1e-15){abort_r("src position phase equal zero!\n ");}

}



}


#endif

