#pragma once
#ifndef fft_desc_h
#define fft_desc_h

#include <stdlib.h> 
#include <time.h>
#include <stdio.h>
#include <algorithm>

#include <complex>
#include <cmath>

#include "general_funs.h"

namespace qlat{

class fft_desc_basic
{
public:
  LInt noden;int rank;int Nmpi;
  int nx,ny,nz,nt;
  LInt vol,Nvol;
  int inix,iniy,iniz,init;
  int Nx,Ny,Nz,Nt;
  int mx,my,mz,mt;

  //std::vector< std::vector<int> > ranged;
  qlat::vector<int> Nv,nv,mv;
  std::vector<qlat::vector<int> > Pos0;
  std::vector<qlat::vector<int>  > mi_list;

  int order_ch;
  int f0,f1,f2,f3;
  qlat::vector<int > orderN;
  std::vector<qlat::vector<LInt > > maporder_Nitoi;
  std::vector<qlat::vector<LInt > > maporder_itoNi;
  std::vector<qlat::vector<LInt > > maporder_Nitoipos;

  const qlat::Geometry* geop;

  fft_desc_basic(const qlat::Geometry& geo,int order_ch_or=0, int check_mem=1)
  {
  
    inix = 1;iniy = 1;iniz = 1;init = 1;
    Nx = 1;  Ny = 1;  Nz = 1;  Nt = 1;
    mx = 1;  my = 1;  mz = 1;  mt = 1;
    ////Order change of memory within a node
    order_ch = order_ch_or;
  
    geop = &geo;
    ///noden = desc->sites_on_node;
    ///rank = get_node_rank();
    //vol = nx*ny*nz;Nmpi = vol*nt/noden;
    Nmpi  = qlat::get_num_node();
    rank  = qlat::get_id_node();

    Nv.resize(4);nv.resize(4);mv.resize(4);
    for(int i=0;i<4;i++){Nv[i]=geo.node_site[i];nv[i] = geo.node_site[i] * geo.geon.size_node[i];}
    nx = nv[0];ny = nv[1];nz = nv[2];nt = nv[3];
    Nx = Nv[0];Ny = Nv[1];Nz = Nv[2];Nt = Nv[3];
    ///Nv[0] = Nx;Nv[1] = Ny;Nv[2] = Nz;Nv[3] = Nt;
    ///nv[0] = nx;nv[1] = ny;nv[2] = nz;nv[3] = nt;
    ///nx=desc->nx;ny=desc->ny;nz=desc->nz;nt=desc->nt;

    mx = nx/Nx;my = ny/Ny;mz = nz/Nz;mt = nt/Nt;
    mv[0] = mx;mv[1] = my;mv[2] = nz;mv[3] = mt;
    vol   =  nx*ny*nz;////Only the spatial volume
    Nvol  =  Nx*Ny*Nz*Nt;
    noden =  Nx*Ny*Nz*Nt;
    print0("=Need=====Size of dim,   x %5d  y %5d  z %5d  t %5d nvol %ld \n",nx,ny,nz,nt,vol*nt);
    print0("=Need=====Size of dim,  Nx %5d Ny %5d Nz %5d Nt %5d Nvol %ld \n",Nx,Ny,Nz,Nt,Nvol);

    qlat::Coordinate ts = geo.coordinate_from_index(0);
    qlat::Coordinate gs = geo.coordinate_g_from_l(ts);
    inix = gs[0];iniy = gs[1];iniz = gs[2];init = gs[3];

    ////Do not assume 0 is the initial positions
    std::vector<int> Pos0_tem;Pos0_tem.resize(Nmpi*4);
    for(unsigned int i=0;i<Pos0_tem.size();i++){Pos0_tem[i] = 0;}
    for(unsigned int i=0;i<4;i++)Pos0_tem[rank*4 + i] = gs[i];
    sum_all_size(&Pos0_tem[0], Nmpi*4);
    Pos0.resize(Nmpi);
    for(int ri=0;ri<Nmpi;ri++)
    {
      Pos0[ri].resize(4);
      for(int i=0;i<4;i++){Pos0[ri][i] = Pos0_tem[ri*4+i];}
    }
  
    /////Default memory order, nt, nz, ny, nx
    f0 = nt;f1 = nz;f2 = ny;f3 = nx;
    orderN.resize(3);orderN[0] = 2;orderN[1] = 1;orderN[2] = 0;

    /////Define fatest order to inner loop
    if(order_ch == 1){
      int maxD = 0;int maxN = Nv[maxD];
      for(int i=0;i<3;i++){if(Nv[i]>maxN){maxD = i;maxN = Nv[i];}}
      orderN[2] = maxD;
      int ci = 0;
      for(int i=0;i<3;i++){if(i != maxD){orderN[ci] = i;ci++;}}
      f1 = nv[orderN[0]];f2 = nv[orderN[1]]; f3 = nv[orderN[2]];
    }
    print0("=Order====Size of dim,  M0 %5d M1 %5d M2 %5d, d0 %2d d1 %2d d2 %2d. \n",
            Nv[orderN[0]],Nv[orderN[1]],Nv[orderN[2]],orderN[0],orderN[1],orderN[2]);


    /////Check assumption, Initial point is the smallest point
    if(check_mem == 1){
      std::vector< std::vector<int> > ranged;
      ranged.resize(4);
      for(size_t isp= 0; isp < size_t(noden); isp++)
      {
        //////position p = desc->get_position(i,rank);
        qlat::Coordinate t0 = geo.coordinate_from_index(isp);
        qlat::Coordinate g0 = geo.coordinate_g_from_l(t0);
        for(int i=0;i<4;i++){ranged[i].push_back(g0[i]);}
      }
  
      for(int i=0;i<4;i++){
        std::sort(ranged[i].begin(), ranged[i].end());
        ranged[i].erase( unique( ranged[i].begin(), ranged[i].end() ), ranged[i].end() );
      }

      int flagC = 0;
      for(unsigned int i=0;i<4;i++){unsigned int tem = Nv[i];if(tem != ranged[i].size()){flagC = 1;}}
      for(unsigned int i=0;i<4;i++){if(gs[i] != ranged[i][0]){flagC = 1;}}
      sum_all_size(&flagC,1);
      if(flagC>0){
        for(int i=0;i<4;i++)print0("%5d %5d, ", Nv[i], int(ranged[i].size()));
        print0("\n");
        for(int i=0;i<4;i++)print0("%5d %5d, ", gs[i], int(ranged[i][0]));
        print0("\n");
        abort_r("Layout not continuous in x! \n");
      }

      /////Check continues between nodes
      std::vector<std::vector<LInt > > map_Nitoi;
      map_Nitoi.resize(Nmpi);
      for(int ri=0;ri<Nmpi;ri++){
        map_Nitoi[ri].resize(noden);
      }

      for(size_t isp=0;isp<size_t(noden);isp++){
        qlat::Coordinate t0 = geo.coordinate_from_index(isp);
        qlat::Coordinate g0 = geo.coordinate_g_from_l(t0);
        LInt offv = ((g0[3]*nz+g0[2])*ny+g0[1])*nx+g0[0];
        map_Nitoi[rank][isp] = offv;
      }

      for(int ri=0;ri<Nmpi;ri++)MPI_Bcast(&map_Nitoi[ri][0], noden*sizeof(LInt), MPI_CHAR, ri, get_comm());

      int flag = 0;
      for(int i=0;i<Nmpi;i++)
      for(size_t isp=0;isp<size_t(noden/Nx);isp++){
        LInt inix = map_Nitoi[i][isp*Nx+ 0];
        for(int xi=0;xi<Nx;xi++)if(map_Nitoi[i][isp*Nx+xi] != inix + xi){flag=1;}
      }
      sum_all_size(&flag,1);
      if(flag>0){abort_r("Layout not continuous in x! \n");}
    }

    ////////position p;
  
    maporder_itoNi.resize(Nmpi);
    maporder_Nitoi.resize(Nmpi);
    maporder_Nitoipos.resize(Nmpi);
    std::vector<int > cur0(4),curi(4);
    for(int ri=0;ri<Nmpi;ri++)
    {
      maporder_itoNi[ri].resize(noden);
      maporder_Nitoi[ri].resize(noden);
      maporder_Nitoipos[ri].resize(noden);
      for(size_t i=0;i < size_t(noden);i++)
      {
        maporder_itoNi[ri][i] = 0;
        maporder_Nitoi[ri][i] = 0;
        maporder_Nitoipos[ri][i] = 0;
      }
    }

    cur0[0] = Pos0[rank][0];
    cur0[1] = Pos0[rank][1];
    cur0[2] = Pos0[rank][2];
    cur0[3] = Pos0[rank][3];
  
    for(size_t i= 0; i < size_t(noden); i++)
    {
      qlat::Coordinate t0 = geo.coordinate_from_index(i);
      qlat::Coordinate p = geo.coordinate_g_from_l(t0);
      ///////Take care of this when porting code
      curi[0] = p[0] - cur0[0];
      curi[1] = p[1] - cur0[1];
      curi[2] = p[2] - cur0[2];

      curi[3] = p[3] - cur0[3];
      LInt newi  =  ((curi[3]*Nv[orderN[0]]+curi[orderN[0]])*Nv[orderN[1]]+curi[orderN[1]])*Nv[orderN[2]]+curi[orderN[2]];
      maporder_itoNi[rank][i] = newi;
      maporder_Nitoi[rank][newi] = i;
  
      LInt Newvi = ((curi[3]*nv[orderN[0]]+p[orderN[0]])*nv[orderN[1]]+p[orderN[1]])*nv[orderN[2]]+p[orderN[2]];
      maporder_Nitoipos[rank][newi] = Newvi;
    }

    for(int ri=0;ri<Nmpi;ri++)
    {
      sum_all_size(&maporder_itoNi[ri][0], noden);
      sum_all_size(&maporder_Nitoi[ri][0], noden);
      sum_all_size(&maporder_Nitoipos[ri][0], noden);
    }
  
    std::vector<int> mi_list_tem;
    mi_list_tem.resize(mt*mz*my*mx);
    for(unsigned int i=0;i<mi_list_tem.size();i++){mi_list_tem[i] = 0;}
    qlat::Coordinate p_tem;
    for(int tmi=0;tmi<mt;tmi++)
    {
      for(int zi=0;zi<mz;zi++)
      for(int yi=0;yi<my;yi++)
      for(int xi=0;xi<mx;xi++)
      {
        LInt off = (zi*my+yi)*mx+xi;
        p_tem[0] = xi*Nx;
        p_tem[1] = yi*Ny;
        p_tem[2] = zi*Nz;
        p_tem[3] = tmi*Nt;
        int flagE = 1;
        for(int i=0;i<4;i++){if(p_tem[i] != gs[i]){flagE = 0;}}
        if(flagE == 1){mi_list_tem[tmi*mz*my*mx + off] = rank;}
      }
    }
    sum_all_size(&mi_list_tem[0], mi_list_tem.size());
    std::vector<int> cp_l;
    for(unsigned int i=0;i<cp_l.size();i++){cp_l[i] = mi_list_tem[i];}
    std::sort(cp_l.begin(), cp_l.end());
    for(unsigned int i=0;i<cp_l.size();i++){
      unsigned int v = cp_l[i];
      //if(v < 0 or v > Nmpi){abort_r("Node map Duplicated! \n");}
      if(v != i ){abort_r("Node map Duplicated! \n");}
    }

    mi_list.resize(mt);
    for(int tmi=0;tmi<mt;tmi++)
    {
      mi_list[tmi].resize(mz*my*mx);
      for(int off=0;off<mi_list[tmi].size();off++){
        mi_list[tmi][off] = mi_list_tem[tmi*mz*my*mx + off];
      }
    }
  
  }

  long get_mi_curr(int dir=3);


  ~fft_desc_basic()
  {
    Pos0.resize(0);
    mi_list.resize(0);
    orderN.resize(0);
    maporder_Nitoi.resize(0);
    maporder_itoNi.resize(0);
    maporder_Nitoipos.resize(0);
  }
};

/////Get the index for current node in 3D or 4D
/////fft redistribution needed
long fft_desc_basic::get_mi_curr(int dir)
{
  int tmi = Pos0[rank][3]/Nt;
  int zmi = Pos0[rank][2]/Nz;
  int ymi = Pos0[rank][1]/Ny;
  int xmi = Pos0[rank][0]/Nx;

  long re = (zmi*my + ymi)*mx + xmi;

  if(dir == 4){re += tmi*mz*my*mx;}

  return re;
}

}

#endif
