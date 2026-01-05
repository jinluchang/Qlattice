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
        bcast_all_size(data[i].data(), data[i].size(), rank, &t_comm[si]);
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

/*
  get src phase off set for different sections
    [Nsrc, Nmom]
*/
template <typename Ty>
inline void get_src_phases(std::vector<vector<Ty > >& src_phases, const std::vector<std::vector<Coordinate > >& pos_new, const std::vector<Coordinate >& sink_momL, const int sign_sink, const int NsecT, const fft_desc_basic& fd){
  TIMER("src phase calculate");
  src_phases.resize(0);
  const int Npsrc = pos_new.size();
  const int Nmom  = sink_momL.size();
  src_phases.resize(Npsrc * sink_momL.size());
  qlat::vector<double > p0;p0.resize(3);
  for(int i=0;i<3;i++){p0[i] = 2 * QLAT_PI_LOCAL / fd.nv[i];}

  for(int srci=0;srci<Npsrc;srci++)
  {
    for(int momi=0;momi<Nmom;momi++){
      ////const int& NsecT = sec.Nsec / 2;// number of src points
      const int  Npsec = pos_new[srci].size();
      Qassert(Npsec == 1 or Npsec == NsecT);
      src_phases[srci*Nmom + momi].resize(NsecT);
      for(int si=0;si < NsecT;si++){
        //pos_new[srci].size(); // src points from each sections
        const Coordinate& sp       = Npsec == 1 ? pos_new[srci][0] : pos_new[srci][si];
        double theta = 0.0;
        // TODO check extra minus sign here
        for(int i=0;i<3;i++){theta += ( p0[i] * sink_momL[momi][i] * ( -1 * sp[i]) );} 
        src_phases[srci*Nmom + momi][si] = Ty(cos(theta), sign_sink * sin(theta));
      }
    }
  }
}

/* 
  shift result vectors based on sections
    [Nsrc_ops (), Nops (Nmom - Nlh - 32)]
  posL : positions of sources, with each section one source 
  Will also apply anti-periodic signs
  with_phases : simple shift or with momL phases and anti-periodic signs
*/
template <typename Ty>
void shift_results_with_phases(std::vector<FieldG<Ty>>& curr_res, const Long Nops, 
  const std::vector<std::vector<Coordinate > >& posL, const std::vector<Coordinate >& sink_momL, 
  sec_list& sec, std::vector<std::vector<FieldG<Ty>>>& curr_shift_buf, std::vector<bool >& baryon_list, const int with_phases = 1){
  TIMER("shift_results_with_phases");
  //const Long Nsrc_ops = pos_new.size();
  //Qassert(Long(curr_res.size()) == Nsrc_ops * Nops);
  Qassert(curr_res.size() % Nops == 0);
  const Long Nsrc_tot = curr_res.size() / Nops;
  const Long Nmom = sink_momL.size();
  Qassert(Nsrc_tot % posL.size() == 0);
  Qassert(Nops % Nmom == 0);

  // Nlh , 32
  const int Nrest = Nops / Nmom;
  // possible mesons / baryons
  const Long Nsrc = posL.size();
  const int Nsrc_ops = Nsrc_tot / Nsrc;
  Qassert(Nsrc_ops == Long(baryon_list.size()));

  std::vector<std::vector<Coordinate > > pos_new;
  pos_new.resize(Nsrc_tot);
  for(Long i=0;i<Nsrc_tot;i++){
    const std::vector<Coordinate >& s = posL[i / Nsrc_ops];
    const int Np = s.size();
    bool same_spatial = true;
    Coordinate s0 = s[0];s0[3] = 0;
    for(int j=1;j<Np;j++){
      Coordinate s1 = s[j];s1[3] = 0;
      if(s0 != s1){same_spatial = false;break;}
    }

    /*
      Inputs of posL should be consistent with section source positions
    */
    for(int j=0;j<Np;j++){
      Qassert(sec.src_t_order[j] == s[j][3]);
    }

    if(!same_spatial){
      pos_new[i].resize(Np);
      for(int j=0;j<Np;j++){
        pos_new[i][j] = s[j];
      }
    }
    // same spatial, then only one shift
    if( same_spatial){
      pos_new[i].resize(1);
      pos_new[i][0] = s[0];
    }

  }

  Qassert(curr_res[0].initialized);
  const Geometry& geo = curr_res[0].geo();
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);

  const int sign_sink = -1;// DEBUG
  const int Nsec = sec.Nsec;Qassert(Nsec % 2 == 0);
  //const vector<int>& map_sec = sec.map_sec;
  if(with_phases == 1){
  std::vector<vector<Ty > > src_phases;src_phases.resize(0);
  get_src_phases(src_phases, pos_new, sink_momL, sign_sink, Nsec / 2, fd);

  // apply source phases and anti-periodic signs
  {
    TIMER("Apply source phases");
    const Long Nvol = geo.local_volume();
    //const int Nsrc_perT = Nsec / 2;
    const int init = fd.init;

    vector<int>& src_t_r_order = sec.src_t_r_order;
    vector<signed char>& anti_sign = sec.anti_sign;
    //vector<int>& src_t = sec.src_t;

    vector<Ty* > resP;resP.resize(Nrest);
    for(Long srci=0;srci<Nsrc_tot;srci++)
    {
      const bool is_baryon = baryon_list[srci % Nsrc_ops];
      const bool antiP = sec.antiP;
      for(Long momi=0;momi<Nmom;momi++)
      {
        vector<Ty >& sphase = src_phases[srci*Nmom + momi];
        //const Ty p0 = phaseP[isp % Nxyz] * sphase[srci];
        for(Long ri=0;ri<Nrest;ri++){
          resP[ri] = (Ty*) get_data(curr_res[srci*Nops + momi * Nrest + ri]).data(); 
        }

        qacc_for(isp, Nvol, {
          const Coordinate xl = geo.coordinate_from_index(isp);
          const int ti   = xl[3] + init;
          //const int seci = map_sec[ti];
          // get the source number
          //const int src_num = ((seci + 1) / 2) % (Nsrc_perT);
          const int src_num = src_t_r_order[ti];
          Ty ph = sphase[src_num];
          // anti-periodic signs
          if(is_baryon and antiP){
            ph = ph * Ty(anti_sign[ti], 0.0);
          }
          for(Long ri=0;ri<Nrest;ri++){
            resP[ri][isp] = resP[ri][isp] * ph;
          }
        });
      }
    }
  }
  }

  bool need_buf = false;
  for(Long si=0;si<Nsrc_tot;si++){
    if(pos_new[si].size() != 1){need_buf = true;}
  }

  std::vector<FieldG<Ty> > shift_res;shift_res.resize(Nops);
  if(need_buf){
    curr_shift_buf.resize(2);
    for(int i=0;i<2;i++){
      if(Long(curr_shift_buf[i].size()) != Nops){
        curr_shift_buf[i].resize(0);
        curr_shift_buf[i].resize(Nops);
      }
      for(Long j=0;j<Nops;j++){
        curr_shift_buf[i][j].init_size(curr_res[0]);
      }
    }
  }

  const int max_shift_grid = 128;

  for(Long si=0;si<Nsrc_tot;si++){
    // TODO check change it for different section T
    if(pos_new[si].size() == 1){
      for(Long j=0;j<Nops;j++){
        shift_res[j].set_pointer( curr_res[si*Nops + j] );
      }
      shift_fields_grid(shift_res, shift_res, pos_new[si][0], -1, max_shift_grid);
    }
    else
    {
      Qassert(false);
      //for(Long j=0;j<Nops;j++){
      //  shift_res[j].set_pointer( curr_res[si*Nops + j] );
      //}
      //copy_fieldsG(curr_shift_buf[0], shift_res);
      //const int Nsrc_perT = Nsec / 2;
      //Qassert(Nsrc_perT == int(pos_new[si].size()));
      //const int init = fd.init;
      //const Long Nvol = geo.local_volume();

      //vector<Ty*> resP = FieldG_to_pointers(shift_res);
      //vector<Ty*> srcP = FieldG_to_pointers(curr_shift_buf[1]);
      //vector<int>& src_t_order = sec.src_t_order;
      //vector<int>& src_t = sec.src_t;

      //const int stime = pos_new[si][0][3];//shift only first time slice to zero
      //// copy only the time slice start with srci 
      //for(int srci=0;srci<Nsrc_perT;srci++){
      //  Coordinate sp_curr = pos_new[si][srci];sp_curr[3] = stime;
      //  shift_fields_grid(curr_shift_buf[0], curr_shift_buf[1], sp_curr, -1, max_shift_grid);
      //  //const int nt = fd.nt;
      //  qacc_for(isp, Nvol, {
      //    const Coordinate xl = geo.coordinate_from_index(isp);
      //    const int ti      = xl[3] + init;
      //    //const int seci    = map_sec[(ti + stime)%nt];// time already shifted
      //    //const int src_num = ((seci + 1) / 2) % (Nsrc_perT);
      //    //if(src_num == srci)
      //    if(src_t[ti] == src_t_order[srci])
      //    {
      //      for(Long op=0;op<Nops;op++){
      //        resP[op][isp] = srcP[op][isp];
      //      }
      //    }
      //  });
      //}
    }
  }
}

// 2pt shifts
template <typename Ty>
void shift_results_with_phases_2pt(qlat::vector_gpu<Ty >& res, const Geometry& geo, 
  std::vector<Coordinate >& pos_src, sec_list& sec, std::vector<std::vector<FieldG<Ty>>>& shift_buf){
  TIMER("shift_results_with_phases_2pt");
  const Long Vol = geo.local_volume();
  Qassert(Long(res.size()) == 32 * Vol);
  std::vector<bool > baryon_list;baryon_list.resize(32);
  for(int i= 0;i<16;i++){baryon_list[i] = false;}
  for(int i=16;i<32;i++){baryon_list[i] = true ;}

  std::vector<std::vector<Coordinate > > posL;posL.resize(1);
  for(int i=0;i<1;i++){
    posL[i].resize(pos_src.size());
    for(unsigned int j=0;j<pos_src.size();j++){
      posL[i][j] = pos_src[j];
    }
  }
  std::vector<Coordinate > sink_momL;sink_momL.resize(1);
  sink_momL[0] = Coordinate(0, 0, 0, 0);

  std::vector<FieldG<Ty>> curr;curr.resize(32);
  for(int i=0;i<32;i++){
    Ty* srcp = &res[i * Vol];
    curr[i].set_pointer(srcp, Vol, geo, QMGPU, QLAT_OUTTER);
  }

  // apply baryon antiperiodic signs
  const bool antiP = sec.antiP;
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);
  if(antiP){
    TIMER("Apply 2pt antiP");
    const int init = fd.init;

    vector<signed char>& anti_sign = sec.anti_sign;
    for(Long srci=0;srci<32;srci++)
    {
      const bool is_baryon = baryon_list[srci];
      Ty* resP = &res[srci * Vol];

      qacc_for(isp, Vol, {
        const Coordinate xl = geo.coordinate_from_index(isp);
        const int ti   = xl[3] + init;
        if(is_baryon){
          resP[isp] = resP[isp] * Ty(anti_sign[ti], 0.0);
        }
      });
    }
  }

  const int Nops = 32;
  std::vector<bool > baryon_list_;baryon_list_.resize(1);
  baryon_list_[0] = false;
  shift_results_with_phases(curr, Nops, posL, sink_momL, sec, shift_buf, baryon_list_, 0);
}

template <typename Ty>
void shift_results_with_phases_2pt(qlat::vector_gpu<Ty >& res, const Geometry& geo, 
  qlat::PointsSelection& pL, sec_list& sec, std::vector<std::vector<FieldG<Ty>>>& shift_buf){
  std::vector<Coordinate > pos_src;
  pos_src.resize(pL.size());
  for(unsigned int i=0;i<pos_src.size();i++){pos_src[i] = pL[i];}
  shift_results_with_phases_2pt(res, geo, pos_src, sec, shift_buf);
}

}

#endif
