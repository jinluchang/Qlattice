// utils_sector_funs.h
// Gen Wang
// Jun. 2025

#ifndef UTILS_SECTOR_FUNS_H
#define UTILS_SECTOR_FUNS_H

#pragma once

#include "utils_float_type.h"
#include "utils_gammas.h"

namespace qlat{

/*
  Create sections with tini and dT
  Buffer sends for each sections
*/
inline std::vector<Int >  get_map_sec(Int dT,Int nt){
  std::vector<Int > map_sec;map_sec.resize(nt);
  Qassert(nt % dT == 0);
  Int secN = 2 * nt / dT;double lensec = nt/(1.0*secN);
  Int tcount = 0;
  Int t0 = 0;
  for(Int si=0;si<secN;si++)
  {
    for(Int t=t0;t < (si+1)*lensec;t++)///boundary with the same sector?
    {
      Qassert(t < nt);
      map_sec[t] = si;
      tcount = tcount + 1;
    }
    t0 = tcount;
  }
  return map_sec;
}

inline void get_map_sec(vector<Int >& map_sec, Int tini, Int dT,Int nt, bool message = false){
  std::vector<Int > map_sec_0 = get_map_sec(dT, nt);
  map_sec.resize(map_sec_0.size());
  std::string m0 = "t ";
  std::string m1 = "s ";
  //shift section number to 0 for first ordered time slice
  for(unsigned int t=0;t<map_sec.size();t++)
  {
    map_sec[t] = map_sec_0[ ( t - tini + nt)%nt];
    m0 += ssprintf("%3d ", t);
    m1 += ssprintf("%3d ", map_sec[t]);
  }
  if(message){qmessage("%s\n%s\n", m0.c_str(), m1.c_str());}
}

inline void get_src_times(vector<Int>& src_t, vector<Int>& src_t_order, vector<Int >& map_sec, const Int tini, const Int dT){
  Qassert(map_sec.size() != 0);
  const Int nt = map_sec.size();
  if(src_t.size() != nt){src_t.resize(nt);}
  const Long Nsrc_perT   = nt / dT;
  const Int Nsec = 2 * Nsrc_perT;

  //const Int tl = (tk + (Nsrc_perT-1) * dT) % nt;// last time
  const Int tk = (tini + dT)%dT;// first time 
  const Int sec0 = map_sec[tk];
  src_t_order.resize(Nsrc_perT);
  for(Int src_num=0;src_num < Nsrc_perT; src_num++){
    const Int src_time   = (src_num * dT + tk)%nt;// src time slice
    src_t_order[src_num] = src_time;
  }
  for(Int t=0;t<nt;t++){
    const Int seci = map_sec[t];// sector number shift to zero
    const Int secj = ( seci - sec0 + Nsec ) % Nsec;// sector number shift to zero
    const Int src_num = (( secj + 1 ) / 2 ) % Nsrc_perT;
    const Int src_time   = (src_num * dT + tk)%nt;// src time slice

    const Int dis0 = (src_time - t + nt ) % nt;
    const Int dis1 = dis0 >  dT/2 ? nt - dis0 : dis0;
    //qmessage("t %5d, src_time %5d, dis1 %5d, dT %d \n", t, src_time, dis1, dT);
    Qassert(dis1 <= dT/2);
    src_t[t] = src_time;
  }
}

struct sec_list{
  box<Geometry> geoB;

  Int tini;
  Int dT;
  Int nt;

  Int Nt;
  Int Nsec;
  bool antiP;
  std::vector<Int > Loop_sec;// defines the calculation order of sections start from 0, gap of Nsec / 2
  /*
    get sector number from time
    map_sec  [nt]  : global variable to check t is within which sections
    has_sec [Nsec] : local variable to check whether a sector is within nodes
    anti_sign : anti periodic signs for baryon
    src_t     : the source position for each time slice
    src_t_r_order : get src number of each time slice
    src_t_order : order the sources
    map to sections
  */
  vector<Int > map_sec;
  vector<Int > has_sec;
  vector<signed char > anti_sign;
  vector<Int > src_t;
  vector<Int > src_t_order;
  vector<Int > src_t_r_order;

  ////communicater with the same init
  ////for spatial global sum
  MPI_Comm xyz_comm;
  ////for time direction global sum
  std::vector< MPI_Comm > t_comm;
  ////rank within each mt
  std::vector<Int > map_mpi_xyz;
  std::vector<std::vector<Int > > map_mpi_t;

  sec_list(){
    tini = -1;
    dT   =  0;
    nt   =  0;
  }

  inline std::vector<Int > get_sinkt(const Int sep)
  {
    Qassert(nt != 0);
    std::vector<Int > sinkt;
    // src sink seperation must be less than half T
    Qassert(sep < (dT/2));
    for(Int si=0;si<Nsec;si++)
    {
      Int tsink = 0;
      if(si%2 == 0){tsink = (tini + dT*(si/2) + sep)%nt;}
      if(si%2 == 1){tsink = (tini + dT*(si/2 + 1) - sep)%nt;}
      sinkt.push_back(tsink);
    }
    std::sort(sinkt.begin(), sinkt.end());// sort the sink time slices
    return sinkt;
  }

  inline void setup(){
    Qassert(nt != 0);
    get_src_times(src_t, src_t_order, map_sec, tini, dT);
    src_t_r_order.resize(nt);
    for(Int t=0;t<nt;t++){
      const Int src_time = src_t[t];
      for(unsigned int si=0;si<src_t_order.size();si++)
      {
        if(src_t_order[si] == src_time){
          src_t_r_order[t] = si;
          break;
        }
      }
    }
    if(anti_sign.size() != nt){anti_sign.resize(nt);}
    const Long Nsrc_perT   = nt / dT;
    const Int tk = (tini + dT)%dT;// first time 
    const Int tl = (tk + (Nsrc_perT-1) * dT) % nt;// last time
    for(Int t=0;t<nt;t++){
      anti_sign[t] = 1;
      const Int src_time = src_t[t];
      if(antiP){
        // dT/2 not included in forward section
        // if src is first time and time cross boundary
        if(src_time == tk and t > tk and t >= tk + int(dT/2)){anti_sign[t] = -1;}

        // if src is last time and time cross boundary
        if(src_time == tl and t < tl and t <  tl - int(dT/2)){anti_sign[t] = -1;}
      }
    }
  }

  inline void init(const Geometry& geo_, const Int tini_, const Int dT_, const bool antiP_ = false, const bool message = true){
    TIMERA("Initialize sec_list ");
    geoB.set(geo_);
    fft_desc_basic& fd = get_fft_desc_basic_plan(geoB());
    nt = fd.nt;
    dT = dT_;
    antiP = antiP_;

    Qassert(nt >= dT);
    Qassert(nt % dT == 0);
    Qassert(dT % 2  == 0);

    Nt = fd.Nt;
    ////total seci numbers
    Nsec = 2 * nt / dT;
    Loop_sec.resize(Nsec);
    for(Int si=0;si<Nsec/2;si++)
    {
      Loop_sec[si*2 + 0] = si;
      Loop_sec[si*2 + 1] = si + Nsec/2;
    }

    update(tini_, message);
  }

  inline void update(Int tini_, const bool message = true){
    TIMERA("update sec_list ");
    if(tini == tini_ % dT){return ;}
    tini = tini_ % dT;

    fft_desc_basic& fd = get_fft_desc_basic_plan(geoB());
    ////from t (without tini) to seci
    get_map_sec(map_sec, tini, dT, nt, message);
    //std::vector<Int > map_sec_0 = get_map_sec(dT, nt);
    //map_sec.resize(map_sec_0.size());
    //for(unsigned int t=0;t<map_sec.size();t++)
    //{
    //  map_sec[t] = map_sec_0[ ( t - tini + nt)%nt];
    //  if(message)qmessage("t %5d, sec %5d \n", t, map_sec[t]);
    //}
    has_sec.resize(Nsec); //// local variable to check whether a sector is within nodes
    for(Int si=0;si<Nsec;si++)
    {
      has_sec[si] = 0;
      for(Int t0=0;t0<fd.Nt;t0++)
      {
        if( map_sec[t0 + fd.init] == si){has_sec[si] = 1;}
      }
    }

    Int Nmpi = fd.Nmpi;
    map_mpi_xyz.resize(Nmpi);
    for(unsigned long mapi=0;mapi<map_mpi_xyz.size();mapi++){
      map_mpi_xyz[mapi] = 0;
    }

    Int color_xyz = fd.init/fd.Nt;
    MPI_Comm_split(get_comm() ,color_xyz, fd.rank, &xyz_comm);
    {
      Int int_tem = -1;
      MPI_Comm_rank(xyz_comm, &int_tem);
      map_mpi_xyz[fd.rank] = int_tem;
    }
    sum_all_size((int*) (&map_mpi_xyz[0]),Nmpi);

    //  split bcase for each sectors
    t_comm.resize(Nsec);
    for(Int s0=0;s0<Nsec;s0++){
      Int si = Loop_sec[s0];
      Int color_t = fd.get_xyzmi_curr();
      color_t = color_t * 2 + has_sec[si];
      //printf("seci %5d, rank %5d, color %5d \n", si, fd.rank, color_t);
      MPI_Comm_split(get_comm() ,color_t, fd.rank, &t_comm[si]);
      Qassert(&t_comm[si] != NULL);
    }

    map_mpi_t.resize(Nsec);
    {
      Int int_tem = -1;
      for(Int s0=0;s0<Nsec;s0++){
        Int si = Loop_sec[s0];
        map_mpi_t[si].resize(Nmpi);
        for(unsigned int mi=0;mi<map_mpi_t[si].size();mi++)
        {
          map_mpi_t[si][mi] = 0;
        }
        MPI_Comm_rank(t_comm[si], &int_tem);
        map_mpi_t[si][fd.rank] = int_tem;
      }
    }
    for(Int s0=0;s0<Nsec;s0++){
      Int si = Loop_sec[s0];
      sum_all_size((int*) (&map_mpi_t[si][0]),Nmpi);
    }
    setup();
  }

  /*
    gsum xyz
    then bcast with from the host_ti
    not sure whether it's necessary since the data is usually small ?
  */
  template <typename Ty>
  void bcast_sink_vecs(std::vector<qlat::vector_gpu<Ty > >& data, const std::vector< Int >& host_ti)
  {
    TIMERA("bcast_sink_vecs");
    Qassert(nt != 0 and data.size() == host_ti.size());
    fft_desc_basic& fd = get_fft_desc_basic_plan(geoB());
  
    /*
      spatial sum only on the host ranks
      within xyz_comm
    */
    for(size_t i=0;i<data.size();i++){
      Qassert(host_ti[i] < fd.nt);
      if(host_ti[i]/fd.Nt == fd.init/fd.Nt) /// host t could be not the initial time of the node
      {
        sum_all_size( data[i].data(), data[i].size(), data[i].GPU, &xyz_comm);
      }
    }
    //  MPI_Barrier(get_comm());
  
    /*
      bcast to other mt nodes from the host ranks
      assuming host_ti is ordered to avoid sector conflicts 
    */
    for(size_t i=0;i<data.size();i++){
      const Int si   = map_sec[host_ti[i]];   //  current sectors of the data
      const Int tmi0 = host_ti[i]/fd.Nt;
      const Int g_rank = fd.mi_list[tmi0][fd.get_xyzmi_curr()];  /// each xyz do the bcast to each other
      const Int rank = map_mpi_t[si][g_rank];   ////host rank for bcast data

      if(i != data.size() - 1){
        if(host_ti[i] >= host_ti[i + 1]){
          print_info();
          for(size_t i=0;i<data.size();i++)
          {
            printf("rank %5d, datai %5d, host_ti %5d, sec %5d has sec %5d \n", fd.rank, int(i), host_ti[i], si, has_sec[si]);
          }
          Qassert(false);
        }
        Qassert( host_ti[i] < host_ti[i + 1] );
      }


      ////printf("rank %5d, local %5d, seci %5d \n", fd.rank, rank, si);
      /* 
        bcast from the data sector
        bcast only to the needed sectors
      */
      if(has_sec[si] == 1){  
        bcast_all_size(data[i].data(), data[i].size(), rank, data[i].GPU, &t_comm[si]);
      }
    }
    MPI_Barrier(get_comm());
  }

  inline void print_info(){
    Qassert(nt != 0);
    for(Int t=0;t<nt;t++){
      qmessage("t %5d, seci %5d, src %5d, anti %5d \n", t, map_sec[t], src_t[t], anti_sign[t]);
    }
  }
};


}

#endif
