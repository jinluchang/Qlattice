#pragma once
#ifndef utils_fft_desc_h
#define utils_fft_desc_h

#include <stdlib.h> 
#include <time.h>
#include <stdio.h>
#include <algorithm>

#include <complex>
#include <cmath>

#include "utils_mpi.h"

namespace qlat{

struct fft_desc_basic
{
  LInt noden;int rank;int Nmpi;
  Int nx,ny,nz,nt;
  LInt vol,Nvol;
  Int inix,iniy,iniz,init;
  Int Nx,Ny,Nz,Nt;
  Int mx,my,mz,mt;

  qlat::vector<Int> iniv,Nv,nv,mv;     ///4 

  std::vector<std::vector<Int> > Pos0;     ////size Nmpi * 4
  std::vector<std::vector<Int>  > mi_list; ////mt , mx*my*mz
  box<Geometry > geoB;

  Int variable_set;

  Int order_ch;
  std::vector<Int > orderN;   ////4

  fft_desc_basic(Int order_ch_or=0)
  {
    TIMERA("Create fft_desc_basic");
  
    variable_set = -1;
    inix =-1;iniy =-1;iniz =-1;init =-1;
    Nx =-1;  Ny =-1;  Nz =-1;  Nt =-1;
    mx =-1;  my =-1;  mz =-1;  mt =-1;
  
    Nmpi  = qlat::get_num_node();
    rank  = qlat::get_id_node();
    order_ch = order_ch_or;

    ////Need set
    nv.resize(4);Nv.resize(4);iniv.resize(4);
  }

  fft_desc_basic(const qlat::Geometry& geo,Int order_ch_or=0)
  {
    TIMER("Create fft_desc_basic");
  
    variable_set = -1;
    inix =-1;iniy =-1;iniz =-1;init =-1;
    Nx =-1;  Ny =-1;  Nz =-1;  Nt =-1;
    mx =-1;  my =-1;  mz =-1;  mt =-1;
    ////Order change of memory within a node
  
    Nmpi  = qlat::get_num_node();
    rank  = qlat::get_id_node();

    Qassert(geo.node_site != qlat::Coordinate());

    ////set nv, Nv, mv, iniv
    geo_to_nv(geo, nv, Nv, mv);

    qlat::Coordinate ts = geo.coordinate_from_index(0);
    qlat::Coordinate gs = geo.coordinate_g_from_l(ts);
    iniv.resize(4);for(unsigned int i=0;i<4;i++){iniv[i] = gs[i];}

    geoB.set_view(get_geo_cache(geo));

    order_ch = order_ch_or;
    set_variable();
    ////check_mem();
  }

  inline void set_variable();
  inline void print_info();
  inline Long get_mi_curr(Int dir=3);
  inline Long get_tmi_curr();
  inline Long get_xyzmi_curr();
  inline void check_mem();
  inline Coordinate coordinate_l_from_index(LInt isp) const;
  inline Coordinate coordinate_g_from_g_index(LInt isp) const ;
  inline Coordinate coordinate_g_from_index(LInt isp, Int rank_set = -1) const;
  inline LInt index_g_from_local(LInt isp, Int rank_set = -1) const;
  //////LInt index_g_from_g_coordinate(std::vector<Int > pos);
  inline LInt index_g_from_g_coordinate(Int t, Int z, Int y, Int x) const;
  inline LInt index_g_from_g_coordinate(const Coordinate& g0) const ;
  inline bool coordinate_g_is_local(const Coordinate& g0) const ;

  inline LInt index_l_from_g_coordinate(const Coordinate& pos, Int rank_set = -1) const;

  inline const Geometry& geo(){
    return geoB();
    //const Coordinate total_site = Coordinate(nx, ny, nz, nt);
    //return qlat::get_geo_cache(total_site);
  }

  inline size_t get_prop_size(){
    size_t Nsize = 1;
    Nsize = noden;
    Nsize = 12*12*Nsize;
    return Nsize;
  }

  ~fft_desc_basic()
  {
    Pos0.resize(0);
    mi_list.resize(0);
    orderN.resize(0);
    Nv.resize(0);
    nv.resize(0);
    mv.resize(0);
  }
};

inline void copy_fft_desc(fft_desc_basic& fd1, const fft_desc_basic& fd0)
{
  fd1.nv = fd0.nv;fd1.Nv = fd0.Nv;fd1.iniv = fd0.iniv;
  fd1.set_variable();
}

/////Get the index for current node in 3D or 4D
/////fft redistribution needed
inline Long fft_desc_basic::get_mi_curr(Int dir)
{
  Int tmi = Pos0[rank][3]/Nt;
  Int zmi = Pos0[rank][2]/Nz;
  Int ymi = Pos0[rank][1]/Ny;
  Int xmi = Pos0[rank][0]/Nx;

  Long re = (zmi*my + ymi)*mx + xmi;
  if(dir == 4){
    if(Nv[0] == nv[0]){re += tmi*mz*my*mx;}
    if(Nv[0] != nv[0]){re = re*(mt) + tmi;}
  }

  return re;
}

inline Long fft_desc_basic::get_tmi_curr()
{
  Int tmi = Pos0[rank][3]/Nt;
  return tmi;
}

inline Long fft_desc_basic::get_xyzmi_curr()
{
  return get_mi_curr(3);
}


inline Coordinate fft_desc_basic::coordinate_l_from_index(LInt isp) const
{
  if(variable_set == -1){abort_r("fft_desc not set! \n");}

  ////std::vector<Int > g0;g0.resize(4);for(Int i=0;i<4;i++){g0[i] = 0;}
  Coordinate g0;for(Int i=0;i<4;i++){g0[i] = 0;}
  g0[3] += isp/(Nv[0]*Nv[1]*Nv[2]);Long tem = isp%(Nv[0]*Nv[1]*Nv[2]);
  g0[orderN[0]] += tem/(Nv[orderN[1]]*Nv[orderN[2]]);tem = tem%(Nv[orderN[1]]*Nv[orderN[2]]);
  g0[orderN[1]] += tem/(Nv[orderN[2]]);
  g0[orderN[2]] += tem%(Nv[orderN[2]]);

  return g0;
}

inline Coordinate fft_desc_basic::coordinate_g_from_g_index(LInt isp) const
{
  if(variable_set == -1){abort_r("fft_desc not set! \n");}

  ////std::vector<Int > g0;g0.resize(4);
  Coordinate g0;for(Int i=0;i<4;i++){g0[i] = 0;}
  g0[3] += isp/(nv[0]*nv[1]*nv[2]);Long tem = isp%(nv[0]*nv[1]*nv[2]);
  g0[orderN[0]] += tem/(nv[orderN[1]]*nv[orderN[2]]);tem = tem%(nv[orderN[1]]*nv[orderN[2]]);
  g0[orderN[1]] += tem/(nv[orderN[2]]);
  g0[orderN[2]] += tem%(nv[orderN[2]]);
  return g0;
}


inline LInt fft_desc_basic::index_g_from_g_coordinate(const Coordinate& g0)  const
{
  if(variable_set == -1){abort_r("fft_desc not set! \n");}
  /////if(g0.size() != 4){abort_r("input pos wrong! \n");}

  LInt index = ((g0[3]*Nv[orderN[0]] + g0[orderN[0]])*Nv[orderN[1]] + g0[orderN[1]])*Nv[orderN[2]] + g0[orderN[2]];
  return index;
}

inline bool fft_desc_basic::coordinate_g_is_local(const Coordinate& g0) const
{
  /////if(g0.size() != 4){abort_r("input pos wrong! \n");}
  bool is_local = true;
  for(Int i=0;i<4;i++)
  {
    if(g0[i] < Pos0[rank][i] or g0[i] >= (Pos0[rank][i]+Nv[i])){is_local = false;}
  }
  return is_local;
}


inline LInt fft_desc_basic::index_g_from_g_coordinate(Int t, Int z, Int y, Int x) const
{
  if(variable_set == -1){abort_r("fft_desc not set! \n");}
  /////std::vector<Int > g0;g0.resize(4);
  Coordinate g0;
  g0[0] = x;g0[1] = y;g0[2] = z;g0[3] = t;
  return index_g_from_g_coordinate(g0);
}

inline Coordinate fft_desc_basic::coordinate_g_from_index(LInt isp, Int rank_set) const
{
  if(variable_set == -1){abort_r("fft_desc not set! \n");}

  Int rank_cur = rank;if(rank_set != -1){rank_cur = rank_set;}
  ////std::vector<Int > g0;g0.resize(4);
  Coordinate g0;/////for(Int i=0;i<4;i++){g0[i] = 0;}
  for(Int i=0;i<4;i++){g0[i] = Pos0[rank_cur][i];}

  g0[3] += isp/(Nv[0]*Nv[1]*Nv[2]);Long tem = isp%(Nv[0]*Nv[1]*Nv[2]);
  g0[orderN[0]] += tem/(Nv[orderN[1]]*Nv[orderN[2]]);tem = tem%(Nv[orderN[1]]*Nv[orderN[2]]);
  g0[orderN[1]] += tem/(Nv[orderN[2]]);
  g0[orderN[2]] += tem%(Nv[orderN[2]]);

  return g0;
}

inline LInt fft_desc_basic::index_l_from_g_coordinate(const Coordinate& pos, Int rank_set) const
{
  if(variable_set == -1){abort_r("fft_desc not set! \n");}

  Int rank_cur = rank;if(rank_set != -1){rank_cur = rank_set;}
  std::vector<Int > p;p.resize(4);
  for(Int i=0;i<4;i++){
    p[i] = pos[i] - Pos0[rank_cur][i];
    if(p[i] < 0 or p[i] >= Nv[i]){abort_r("pos not local!\n");}
  }

  LInt index = ((p[3]*Nv[orderN[0]]+p[orderN[0]])*Nv[orderN[1]]+p[orderN[1]])*Nv[orderN[2]]+p[orderN[2]];

  return index;
}

inline LInt fft_desc_basic::index_g_from_local(LInt isp, Int rank_set) const
{
  if(variable_set == -1){abort_r("fft_desc not set! \n");}
  const Coordinate p = this->coordinate_g_from_index(isp, rank_set);
  LInt newi = ((p[3]*nv[orderN[0]]+p[orderN[0]])*nv[orderN[1]]+p[orderN[1]])*nv[orderN[2]]+p[orderN[2]];
  return newi;
}

inline void fft_desc_basic::print_info(){
  if(variable_set == -1){abort_r("fft_desc not set! \n");}
  qmessage("=Need=====Size of dim,   x %10d  y %10d  z %10d  t %10d nvol %ld \n",nx,ny,nz,nt,vol*nt);
  qmessage("=Need=====Size of dim,  Nx %10d Ny %10d Nz %10d Nt %10d Nvol %ld \n",Nx,Ny,Nz,Nt,Nvol);
  qmessage("=Order====Size of dim,  M0 %10d M1 %10d M2 %10d, d0 %2d d1 %2d d2 %2d. \n",
          Nv[orderN[0]],Nv[orderN[1]],Nv[orderN[2]],orderN[0],orderN[1],orderN[2]);

  fflush_MPI();
}

/////mi_list, same t0 rank list
inline void fft_desc_basic::set_variable()
{
  TIMERA("fft_desc_basic::set_variable");
  for(Int i=0;i<4;i++){Qassert(iniv[i] >= 0);Qassert(nv[i] > 0);Qassert(Nv[i] > 0);}
  mv.resize(4);for(Int i=0;i<4;i++){mv[i] = nv[i]/Nv[i];Qassert(mv[i] > 0);}

  nx = nv[0];ny = nv[1];nz = nv[2];nt = nv[3];
  Nx = Nv[0];Ny = Nv[1];Nz = Nv[2];Nt = Nv[3];
  mx = mv[0];my = mv[1];mz = mv[2];mt = mv[3];

  vol   =  nx*ny*nz;////Only the spatial volume
  Nvol  =  Nx*Ny*Nz*Nt;
  noden =  Nx*Ny*Nz*Nt;

  inix = iniv[0];iniy = iniv[1];iniz = iniv[2];init = iniv[3];

  ////Do not assume 0 is the initial positions
  std::vector<Int> Pos0_tem;Pos0_tem.resize(Nmpi*4);
  for(unsigned int i=0;i<Pos0_tem.size();i++){Pos0_tem[i] = 0;}
  for(unsigned int i=0;i<4;i++)Pos0_tem[rank*4 + i] = iniv[i];
  sum_all_size(Pos0_tem.data(), Nmpi*4);
  Pos0.resize(Nmpi);
  for(Int ri=0;ri<Nmpi;ri++)
  {
    Pos0[ri].resize(4);
    for(Int i=0;i<4;i++){Pos0[ri][i] = Pos0_tem[ri*4+i];}
  }

  orderN.resize(3);orderN[0] = 2;orderN[1] = 1;orderN[2] = 0;

  if(order_ch == 1){
    Int maxD = 0;int maxN = Nv[maxD];
    for(Int i=0;i<3;i++){if(Nv[i]>maxN){maxD = i;maxN = Nv[i];}}
    orderN[2] = maxD;
    Int ci = 0;
    for(Int i=0;i<3;i++){if(i != maxD){orderN[ci] = i;ci++;}}
    ////f1 = nv[orderN[0]];f2 = nv[orderN[1]]; f3 = nv[orderN[2]];
  }

  ////////position p;

  ////depend on rank, geo.coordinate_g_from_l or gs
  std::vector<Int> mi_list_tem;
  mi_list_tem.resize(mt*mz*my*mx);
  for(unsigned int i=0;i<mi_list_tem.size();i++){mi_list_tem[i] = 0;}
  std::vector<Int > p_tem(4);
  /////qlat::Coordinate p_tem;
  for(Int tmi=0;tmi<mt;tmi++)
  {
    for(Int zi=0;zi<mz;zi++)
    for(Int yi=0;yi<my;yi++)
    for(Int xi=0;xi<mx;xi++)
    {
      LInt off = (zi*my+yi)*mx+xi;
      p_tem[0] = xi*Nx;
      p_tem[1] = yi*Ny;
      p_tem[2] = zi*Nz;
      p_tem[3] = tmi*Nt;
      Int flagE = 1;
      for(Int i=0;i<4;i++){if(p_tem[i] != Pos0[rank][i]){flagE = 0;}}
      if(flagE == 1){mi_list_tem[tmi*mz*my*mx + off] = rank;}
    }
  }
  sum_all_size(&mi_list_tem[0], mi_list_tem.size());
  std::vector<Int> cp_l(mi_list_tem.size());
  for(unsigned int i=0;i<cp_l.size();i++){cp_l[i] = mi_list_tem[i];}
  std::sort(cp_l.begin(), cp_l.end());
  for(unsigned int i=0;i<cp_l.size();i++){
    unsigned int v = cp_l[i];
    //if(v < 0 or v > Nmpi){abort_r("Node map Duplicated! \n");}
    if(v != i ){abort_r("Node map Duplicated! \n");}
  }

  mi_list.resize(mt);
  for(Int tmi=0;tmi<mt;tmi++)
  {
    mi_list[tmi].resize(mz*my*mx);
    for(LInt off=0;off<mi_list[tmi].size();off++){
      mi_list[tmi][off] = mi_list_tem[tmi*mz*my*mx + off];
    }
  }
  variable_set = 1;
}

inline void fft_desc_basic::check_mem()
{
  if(variable_set == -1){abort_r("fft_desc not set! \n");}

  /////Check assumption, Initial point is the smallest point
  std::vector< std::vector<Int> > ranged;
  ranged.resize(4);
  for(size_t isp= 0; isp < size_t(noden); isp++)
  {
    Coordinate g0 = coordinate_g_from_index(isp);
    for(Int i=0;i<4;i++){ranged[i].push_back(g0[i]);}
  }

  for(Int i=0;i<4;i++){
    std::sort(ranged[i].begin(), ranged[i].end());
    ranged[i].erase( unique( ranged[i].begin(), ranged[i].end() ), ranged[i].end() );
  }

  Int flagC = 0;
  for(unsigned int i=0;i<4;i++){unsigned int tem = Nv[i];if(tem != ranged[i].size()){flagC = 1;}}
  for(unsigned int i=0;i<4;i++){if(Pos0[rank][i] != ranged[i][0]){flagC = 1;}}
  sum_all_size(&flagC,1);
  if(flagC>0){
    for(Int i=0;i<4;i++)qmessage("%5d %5d, ", Nv[i], int(ranged[i].size()));
    qmessage("\n");
    for(Int i=0;i<4;i++)qmessage("%5d %5d, ", Pos0[rank][i], int(ranged[i][0]));
    qmessage("\n");
    abort_r("Layout not continuous in x, A! \n");
  }

  /////Check continues between nodes
  std::vector<std::vector<LInt > > map_Nitoi;
  map_Nitoi.resize(Nmpi);
  for(Int ri=0;ri<Nmpi;ri++){
    map_Nitoi[ri].resize(noden);
  }

  for(size_t isp=0;isp<size_t(noden);isp++){
    Coordinate g0 = coordinate_g_from_index(isp);
    LInt offv = ((g0[3]*nz+g0[2])*ny+g0[1])*nx+g0[0];
    map_Nitoi[rank][isp] = offv;
  }

  for(Int ri=0;ri<Nmpi;ri++)MPI_Bcast(&map_Nitoi[ri][0], noden*sizeof(LInt), MPI_CHAR, ri, get_comm());

  Int flag = 0;
  for(Int i=0;i<Nmpi;i++)
  for(size_t isp=0;isp<size_t(noden/Nx);isp++){
    LInt inix = map_Nitoi[i][isp*Nx+ 0];
    for(Int xi=0;xi<Nx;xi++)if(map_Nitoi[i][isp*Nx+xi] != inix + xi){flag=1;}
  }
  sum_all_size(&flag,1);
  if(flag>0){abort_r("Layout not continuous in x, B! \n");}

}

///mode 0, mxyz ->  Nt ; mode 1, Nt -> mxyz
inline void desc_xyz_in_one(fft_desc_basic& fd, const Geometry& geo, Int mode = 1){
  qlat::vector<Int> Nv,nv,mv;
  geo_to_nv(geo, nv, Nv, mv);
  /////dim of lat
  fd.nv[0] = nv[0];
  fd.nv[1] = nv[1];
  fd.nv[2] = nv[2];
  fd.nv[3] = nv[3]*mv[2]*mv[1]*mv[0];
  /////dim of lat on local
  fd.Nv[0] = nv[0];
  fd.Nv[1] = nv[1];
  fd.Nv[2] = nv[2];
  fd.Nv[3] = Nv[3];

  qlat::Coordinate ts = geo.coordinate_from_index(0);
  qlat::Coordinate gs = geo.coordinate_g_from_l(ts);

  ///////dim of MPI lat
  //fd.mv[0] = 1;
  //fd.mv[1] = 1;
  //fd.mv[2] = mv[2]*mv[1]*mv[0];
  //fd.mv[3] = mv[3];
  fd.iniv.resize(4);
  //int t0 = fd.iniv[3];
  for(unsigned int i=0;i<4;i++){fd.iniv[i] = 0;} 
  Coordinate tnode = coor_node_from_id_node(qlat::get_id_node());
  Int tA = (tnode[2]*mv[1] + tnode[1])*mv[0] + tnode[0];
  if(mode == 0){fd.iniv[3] = tA*nv[3] + gs[3]; }
  if(mode == 1){fd.iniv[3] = gs[3] * (mv[2]*mv[1]*mv[0]) + tA * Nv[3]; }
  fd.set_variable();
}

/////fft desc buffers related
struct FFTdescKey {
  Coordinate total_site;
  Int order_ch_or;
  FFTdescKey(const Geometry& geo, Int order_ch_or_ = 0)
  {
    total_site = geo.total_site();
    order_ch_or = order_ch_or_;
  }

};
inline bool operator<(const FFTdescKey& x, const FFTdescKey& y)
{
  if(x.total_site < y.total_site ){  return true;}
  if(y.total_site < x.total_site ){  return false;}
  if(x.order_ch_or < y.order_ch_or ){  return true;}
  if(y.order_ch_or < x.order_ch_or ){  return false;}
  return false;
}

inline Cache<FFTdescKey, fft_desc_basic >& get_fft_desc_basic_cache()
{
  static Cache<FFTdescKey, fft_desc_basic > cache("FFTdescKey", 16);
  return cache;
}

inline fft_desc_basic& get_fft_desc_basic_plan(const FFTdescKey& fkey)
{
  if (!get_fft_desc_basic_cache().has(fkey)) {
    //Geometry geo;geo.init(fkey.total_site);
    const Geometry& geo = get_geo_cache(fkey.total_site);
    get_fft_desc_basic_cache()[fkey] = fft_desc_basic(geo, fkey.order_ch_or);
  }
  fft_desc_basic& buf = get_fft_desc_basic_cache()[fkey];
  return buf;
}

inline fft_desc_basic& get_fft_desc_basic_plan(const Geometry& geo)
{
  FFTdescKey fkey(geo);
  return get_fft_desc_basic_plan(fkey);
}




}

#endif
