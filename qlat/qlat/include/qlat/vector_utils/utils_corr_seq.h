// utils_corr_seq.h
// Gen Wang
// Oct. 2022

#ifndef UTILS_CORR_SEQ_H
#define UTILS_CORR_SEQ_H

#pragma once

#include "utils_float_type.h"
#include "utils_gammas.h"
#include "utils_fft_desc.h"
#include "utils_reduce_vec.h"
#include "utils_grid_src.h"
#include "utils_corr_prop.h"

////#ifdef QLAT_USE_ACC
////#define REDUCEFAC  1
////#else
////#define REDUCEFAC 16
////#endif

namespace qlat{

struct sec_list{
  Geometry geo;

  int tini;
  int dT;
  int nt;

  int Nt;
  int Nsec;
  std::vector<int > Loop_sec;
  ////get sector number from time
  std::vector<int > map_sec;
  std::vector<int > has_sec;

  ////communicater with the same init
  ////for spatial global sum
  MPI_Comm xyz_comm;
  ////for time direction global sum
  std::vector< MPI_Comm > t_comm;
  ////rank within each mt
  std::vector<int > map_mpi_xyz;
  std::vector<std::vector<int > > map_mpi_t;

  sec_list(){
    tini = -1;
    dT   =  0;
    nt   =  0;
  }

  inline std::vector<int > get_sinkt(const int sep)
  {
    Qassert(nt != 0);
    std::vector<int > sinkt;
    Qassert(sep < (dT/2));
    for(int si=0;si<Nsec;si++)
    {
      int tsink = 0;
      if(si%2 == 0){tsink = (tini + dT*(si/2) + sep)%nt;}
      if(si%2 == 1){tsink = (tini + dT*(si/2 + 1) - sep)%nt;}
      sinkt.push_back(tsink);
    }
    return sinkt;
  }

  inline void init(Geometry& geo_, int tini_, int dT_){
    TIMERA("Initialize sec_list ");
    geo = geo_;
    fft_desc_basic& fd = get_fft_desc_basic_plan(geo);
    nt = fd.nt;
    dT = dT_;

    Qassert(nt >= dT);
    Qassert(nt % dT == 0);
    Qassert(dT % 2  == 0);

    Nt = fd.Nt;
    ////total seci numbers
    Nsec = 2 * nt / dT;
    Loop_sec.resize(Nsec);
    for(int si=0;si<Nsec/2;si++)
    {
      Loop_sec[si*2 + 0] = si;
      Loop_sec[si*2 + 1] = si + Nsec/2;
    }

    update(tini_);
  }

  inline void update(int tini_){
    TIMERA("update sec_list ");
    if(tini == tini_){return ;}
    tini = tini_;

    fft_desc_basic& fd = get_fft_desc_basic_plan(geo);
    ////from t (without tini) to seci
    std::vector<int > map_sec_0 = get_map_sec(dT, nt);
    map_sec.resize(map_sec_0.size());
    for(unsigned int t=0;t<map_sec.size();t++)
    {
      map_sec[t] = map_sec_0[ ( t - tini + nt)%nt];
      print0("t %5d, sec %5d \n", t, map_sec[t]);
    }
    has_sec.resize(Nsec); //// local variable to check whether a sector is within nodes
    for(int si=0;si<Nsec;si++)
    {
      has_sec[si] = 0;
      for(int t0=0;t0<fd.Nt;t0++)
      {
        if( map_sec[t0 + fd.init] == si){has_sec[si] = 1;}
      }
    }

    int Nmpi = fd.Nmpi;
    map_mpi_xyz.resize(Nmpi);
    for(unsigned long mapi=0;mapi<map_mpi_xyz.size();mapi++){
      map_mpi_xyz[mapi] = 0;
    }

    int color_xyz = fd.init/fd.Nt;
    MPI_Comm_split(get_comm() ,color_xyz, fd.rank, &xyz_comm);
    {
      int int_tem = -1;
      MPI_Comm_rank(xyz_comm, &int_tem);
      map_mpi_xyz[fd.rank] = int_tem;
    }
    sum_all_size((int*) (&map_mpi_xyz[0]),Nmpi);

    ////split bcase for each sectors
    t_comm.resize(Nsec);
    for(int s0=0;s0<Nsec;s0++){
      int si = Loop_sec[s0];
      int color_t = fd.get_xyzmi_curr();
      color_t = color_t * 2 + has_sec[si];
      //printf("seci %5d, rank %5d, color %5d \n", si, fd.rank, color_t);
      MPI_Comm_split(get_comm() ,color_t, fd.rank, &t_comm[si]);
      Qassert(&t_comm[si] != NULL);
    }

    map_mpi_t.resize(Nsec);
    {
      int int_tem = -1;
      for(int s0=0;s0<Nsec;s0++){
        int si = Loop_sec[s0];
        map_mpi_t[si].resize(Nmpi);
        for(unsigned int mi=0;mi<map_mpi_t[si].size();mi++)
        {
          map_mpi_t[si][mi] = 0;
        }
        MPI_Comm_rank(t_comm[si], &int_tem);
        map_mpi_t[si][fd.rank] = int_tem;
      }
    }
    for(int s0=0;s0<Nsec;s0++){
      int si = Loop_sec[s0];
      sum_all_size((int*) (&map_mpi_t[si][0]),Nmpi);
    }
  }

  ////gsum xyz and bcast with from the host_ti
  template <typename Ty>
  void bcast_sink_vecs(std::vector<qlat::vector_gpu<Ty > >& data, const std::vector< int >& host_ti)
  {
    TIMERA("bcast_sink_vecs");
    Qassert(nt != 0);
    fft_desc_basic& fd = get_fft_desc_basic_plan(geo);
  
    ////spatial sum only on the host ranks
    for(size_t i=0;i<data.size();i++){
      Qassert(host_ti[i] < fd.nt);
      if(host_ti[i]/fd.Nt == fd.init/fd.Nt) /// host t could be not the initial time of the node
      {
        sum_all_size( data[i].data(), data[i].size(), data[i].GPU, &xyz_comm);
      }
    }
    /////MPI_Barrier(get_comm());
  
    ////bcast to other mt nodes from the host ranks
    ////assuming host_ti is ordered to avoid sector conflicts 
    for(size_t i=0;i<data.size();i++){
      int si   = map_sec[host_ti[i]];   ////current sectors of the data
      int tmi0 = host_ti[i]/fd.Nt;
      int g_rank = fd.mi_list[tmi0][fd.get_xyzmi_curr()];  /// each xyz do the bcast to each other
      const int rank = map_mpi_t[si][g_rank];   ////host rank for bcast data
      ////printf("rank %5d, local %5d, seci %5d \n", fd.rank, rank, si);
      if(has_sec[si] == 1){   ////bcast only to the needed sectors
        bcast_all_size(data[i].data(), data[i].size(), rank, data[i].GPU, &t_comm[si]);
      }
    }
    MPI_Barrier(get_comm());
  }

};

template <typename Ty>
void seq_high_single(qpropT& src, qpropT& sinkH, qpropT& noiseH, const std::vector<int >& sinkt, qpropT& res,
  std::vector<Coordinate >& sink_mom, 
  sec_list& sec, int clear=1)
{
  const Geometry& geo = src.geo();
  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);
  const long Nvol = geo.local_volume();
  const long Nxyz = fd.Nx * fd.Ny * fd.Nz;

  Qassert(src.initialized);
  Qassert(sinkH.initialized);
  Qassert(noiseH.initialized);
  if(!res.initialized){res.init(geo);}
  if(clear == 1){qlat::set_zero(res);}

  std::vector<vector_gpu<Ty > > phases;
  get_phases(phases, sink_mom, geo);

  Ty* srcP = (Ty*) qlat::get_data(src).data();
  Ty* sinP = (Ty*) qlat::get_data(sinkH).data();
  Ty* noiP = (Ty*) qlat::get_data(noiseH).data();
  Ty* resP = (Ty*) qlat::get_data(res).data();

  //std::vector<vector_gpu<Ty > > buf_vec;
  //buf_vec.resize(Nsrc * Nsink * Nmom);
  //for(unsigned long i=0;i<buf_vec.size();i++){buf_vec.resize(Nsum * 12 * 12);}
  //qlat::vector_acc<Ty* > buf_vecP = EigenM_to_pointers(buf_vec);

  qlat::vector_gpu<Ty > buf_vec;
  buf_vec.resize(12*12*Nxyz);
  Ty* buf_vecP = buf_vec.data();

  std::vector< qlat::vector_gpu<Ty > > sum_d;
  sum_d.resize(sinkt.size());
  for(unsigned int i=0;i<sum_d.size();i++){sum_d[i].resize(12*12);sum_d[i].set_zero();}

  for(unsigned int hosti=0;hosti<sinkt.size();hosti++)
  {
    const int t0 = sinkt[hosti];
    Qassert( t0 >= 0 and t0 < fd.nt );
    if(t0 < fd.init or t0 >= fd.init + fd.Nt){continue ;}
    const int ti = t0 - fd.init;
    ////buf_vec.set_zero();
    qacc_for(isp, Nxyz, {
      const long off = ti*Nxyz + isp;
      Ty* s0 = &srcP[(0*12 + 0)*Nvol + off];
      Ty* s1 = &noiP[(0*12 + 0)*Nvol + off];
      //Ty* r0 = &buf_vecP[(0*12 + 0)*Nvol + off];
      Ty* r0 = &buf_vecP[(0*12 + 0)*Nxyz + isp];

      for(int dc1=0;dc1<12;dc1++)
      for(int dc2=0;dc2<12;dc2++)
      {
        Ty buf = 0;
        for(int dcs=0;dcs<12;dcs++)
        {
          buf += s0[(dc1*12 + dcs)*Nvol] * s1[(dcs*12 + dc2) *Nvol];
          //buf += s0[(dc1*12 + dcs)*Nvol] * qlat::qconj(s1[(dcs*12 + dc2) *Nvol]);
        }
        r0[(dc1*12 + dc2)*Nxyz] = buf;
      }
    });
    reduce_vecs(buf_vec.data(), sum_d[hosti].data(), buf_vec.size()/(12*12), 12*12);
  }

  sec.bcast_sink_vecs(sum_d, sinkt);
  ////for(int i=0;i<sum_d.size();i++){sum_all_size(sum_d[i].data(), sum_d[i].size());}

  for(int ti=0;ti<fd.Nt;ti++){
    const int ta = ti + fd.init;
    //printf("time %d \n", ta);
    const int seci = sec.map_sec[ta];
    int mapi = -1;
    for(unsigned int si=0;si<sinkt.size();si++)
    {
      ////printf("si %d %d %d \n", si, sinkt[si], int(sec.map_sec.size()));
      if(sec.map_sec[sinkt[si]] == seci)
      {
        if(mapi == -1){mapi = si;}
        else{
          /////print0("Dumplicate %d %d %d \n", mapi, sinkt[si], seci);
          Qassert(false);
        }
      }
    }
    Qassert(mapi >= 0);

    Ty* alpha = sum_d[mapi].data();

    qacc_for(isp, Nxyz,{
      const long off = ti*Nxyz + isp;
      Ty* r0 = &resP[(0*12 + 0)*Nvol + off];
      Ty* s0 = &sinP[(0*12 + 0)*Nvol + off];

      for(int dc1=0;dc1<12;dc1++)
      for(int dc2=0;dc2<12;dc2++)
      {
        Ty buf = 0;
        for(int dcs=0;dcs<12;dcs++)
        {
          buf += alpha[dc1*12 + dcs] * s0[(dcs*12 + dc2)*Nvol] ;
          //buf += alpha[dc1*12 + dcs] * qlat::qconj(s0[(dc2*12 + dcs)*Nvol] ) ;
          //buf += alpha[dc1*12 + dcs] * s0[(dcs*12 + dc2)*Nvol] ;
          //buf += alpha[dcs*12 + dc1] * s0[(dcs*12 + dc2)*Nvol] ;
        }
        r0[(dc1*12 + dc2)*Nvol] = buf;
      }
    });
  }

}

//////sinkH reversed already, with conjugate
//////noiseH smeared already
//////(src * noise) * mom X sinkH
//template <typename Ty, int Nmom>
//void baryon_vectorE(std::vector<qpropT >& src, std::vector<qpropT >& sinkH,
//  std::vector<qpropT >& noiseH, std::vector<int >& sinkT,
//  std::vector<int > sum_list, const int sum_max,
//  std::vector<Coordinate >& sink_mom, 
//  std::vector<qpropT >& res, std::vector<qpropT >& buf,
//  ///std::vector<qpropT >& resP, ga_M &A, ga_M &B, qlat::vector_acc<Ty > &G, qlat::vector_acc<int > &mL, int insertion,int clear=1
//  const int tini, sec_list& sec, int clear=1)
//{
//  fft_desc_basic& fd = get_fft_desc_basic_plan(geo);
//  sec.update(tini);
//  Qassert(sinkT.size()  == sinkH.size());
//  Qassert(noiseH.size() == sinkH.size());
//  const int Nsrc  = src.size();
//  const int Nsink = sum_max;
//  ///const int Nmom  = sink_mom.size()
//  Qassert(Nmom == sink_mom.size()); ////boost buffers
//  const long Nvol = geo.local_volume();
//  const long Nxyz = fd.Nx * fd.Ny * fd.Nz;
//
//  Qassert(Nmom != 0 and Nsrc != 0 and Nsink != 0);
//  const int Nres  = Nsrc * Nsink * Nmom;
//
//  Qassert(src[0].initialized);
//  const Geometry& geo = src[0].geo();
//
//  init_qpropT(res, Nres, geo);
//  if(clear == 1){clear_qpropT(res);}
//
//
//  std::vector<vector_gpu<Ty > > phases;
//  get_phases(phases, sink_mom, geo);
//
//  std::vector<qpropT > srcM;
//  init_qpropT(srcM, Nsrc, geo);
//
//  qlat::vector_acc<Ty* > srcP   = EigenM_to_pointers(src);
//  qlat::vector_acc<Ty* > srcMP  = EigenM_to_pointers(srcM);
//  qlat::vector_acc<Ty* > sinkHP = EigenM_to_pointers(sinkH);
//  qlat::vector_acc<Ty* > noiseHP = EigenM_to_pointers(noiseH);
//  qlat::vector_acc<Ty* > phasesP = EigenM_to_pointers(phases);
//
//  std::vector<vector_gpu<Ty > > buf_vec;
//  buf_vec.resize(Nsrc * Nsink * Nmom);
//
//  Qassert(Nxyz % REDUCEFAC == 0);
//  const long Nsum = Nvol / REDUCEFAC;
//  for(unsigned long i=0;i<buf_vec.size();i++){buf_vec.resize(Nsum * 12 * 12);}
//  qlat::vector_acc<Ty* > buf_vecP = EigenM_to_pointers(buf_vec);
//
//  for(int momi=0;momi<Nmom;momi++)
//  {
//    qacc_for(isp, Nsum,{
//      const long off  = isp*REDUCEFAC + 0;
//      const long xoff = off%Nxyz;
//
//      Ty* m0 = &phasesP[momi][xoff];
//
//      for(int iv=0;iv<Nsrc;iv++)
//      for(int dc=0;dc<12*12;dc++)
//      for(int i=0;i<REDUCEFAC;i++){
//        srcMP[iv][dc*Nvol + off + i] = srcP[iv][dc*Nvol + off + i] * m0[i];
//      }
//    })
//
//    qacc_for(isp, Nsum,{
//      const long off  = isp*REDUCEFAC + 0;
//      const long xoff  = off%Nxyz;
//
//      //Ty* m0[Nmom];
//      //for(int mi=0;mi<Nmom;mi++){m0[mi] = &phasesP[momi][xoff];}
//
//      //for(int i=0;i<12*12;i++){buf[i] = 0;}
//
//      for(int srci=0;srci<Nsrc;srci++)
//      for(int sinki=0;sinki<Nsink;sinki++)
//      {
//        Ty* s0 = &srcMP[srci][(0*12 + 0)*Nvol + off];
//        Ty* s1 = &noiseHP[sinki][(0*12 + 0) *Nvol + off ];
//        Ty* r0 = &buf_vecP[srci*Nsrc + Nsinki][(0*12 + 0)* Nsum + isp];
//
//        for(int dc1=0;dc1<12;dc1++)
//        for(int dc2=0;dc2<12;dc2++)
//        {
//          Ty buf = 0;
//          for(int dcs=0;dcs<12;dcs++)
//          {
//            for(int i=0;i<REDUCEFAC;i++){
//              buf += s0[(dc1*12 + dcs)*Nvol + i] * s1[(dcs*12 + dc2) *Nvol + i ];
//            }
//          }
//          r0[(dc1*12 + dc2)*Nsum + isp] = buf;
//        }
//      }
//    });
//    qlat::vector_acc<Ty* > buf_vecP = EigenM_to_pointers(buf_vec);
//  }
//
//  qlat::vector_acc<Ty* > buf_vecP = EigenM_to_pointers(buf_vec);
//  qlat::vector_gpu<Ty* >
//
//}

////simple sequential sources
template <class Ty, int civ>
void local_sequential_source(qlat::FieldM<Ty, civ>& src, const qlat::vector_acc<int >& tseq)
{
  TIMERA("local_sequential_source");
  Qassert(src.initialized);
  const Geometry& geo = src.geo();
  const long V = geo.local_volume();
  const long Ndata = qlat::get_data_size(src) / sizeof(Ty);
  const int Dim = Ndata / V;

  qlat::vector_acc<int > nv, Nv, mv;geo_to_nv(geo, nv, Nv, mv);
  long  Nvol = Nv[0]*Nv[1]*Nv[2];
  Qassert(tseq.size() > 0);
  for(long i=0;i<tseq.size();i++){
    Qassert(tseq[i] >= 0 and tseq[i] < nv[3]);
  }

  Ty* srcP = (Ty*) qlat::get_data(src).data();
  const int Nt = Nv[3];
  const int nt = nv[3];
  /////const long  V= Nt * Nvol;

  qacc_for(xi, long(Nvol),{
    for(int ti = 0;ti < Nt; ti ++)
    {
      const long isp = ti*Nvol + xi;
      const Coordinate xl  = geo.coordinate_from_index(isp);
      const Coordinate pos = geo.coordinate_g_from_l(xl);
      bool find = false;
      for(long tj = 0;tj < tseq.size(); tj++){
        if(pos[3] == (tseq[tj])%nt)
        {
          find = true;break;
        }
      }
      if(find == false)
      {
        for(int ic=0;ic<Dim;ic++){srcP[isp*Dim + ic] = 0;}
      }
    }
  });
}

template <class Ty, int civ>
void local_sequential_source(qlat::FieldM<Ty, civ>& src, const int tseq)
{
  std::vector<int > tL;tL.resize(1);tL[0] = tseq;
  local_sequential_source(src, tL);
}

template <class Td>
void local_sequential_source(Propagator4dT<Td >& res, Propagator4dT<Td >& src, const qlat::vector_acc<int >& tseq, const int gammai = -1)
{
  TIMERA("local_sequential_source");
  Qassert(src.initialized);
  const Geometry& geo = src.geo();
  const long V = geo.local_volume();
  const long Ndata = qlat::get_data_size(src) / sizeof(qlat::ComplexT<Td>);
  const int Dim = Ndata / V;

  if(!res.initialized){res.init(geo);}
  qlat::vector_acc<int > nv, Nv, mv;geo_to_nv(geo, nv, Nv, mv);
  long  Nvol = Nv[0]*Nv[1]*Nv[2];
  ///const long Ndata = geo.local_volume() * Dim;
  Qassert(tseq.size() > 0);
  for(long i=0;i<tseq.size();i++){
    Qassert(tseq[i] >= 0 and tseq[i] < nv[3]);
  }
  Qassert(gammai >= -1 and gammai < 16);

  qlat::ComplexT<Td >* srcP = (qlat::ComplexT<Td >*) qlat::get_data(src).data();
  qlat::ComplexT<Td >* resP = (qlat::ComplexT<Td >*) qlat::get_data(res).data();
  ////svecT.shift_vecs_dir(tmp[5], tmp[4],  3, -1  );
  const int Nt = Nv[3];
  const int nt = nv[3];
  /////const long  V= Nt * Nvol;

  qacc_for(isp, V, {
    for(int ic=0;ic<Dim;ic++){resP[isp*Dim + ic] = 0;}
  });
  qacc_for(xi, long(Nvol),{
    for(int ti = 0;ti < Nt; ti ++)
    {
      const long isp = ti*Nvol + xi;
      const Coordinate xl  = geo.coordinate_from_index(isp);
      const Coordinate pos = geo.coordinate_g_from_l(xl);
      for(long tj = 0;tj < tseq.size(); tj++){
        if(pos[3] == (tseq[tj])%nt)
        {
          for(int ic=0;ic<Dim;ic++){resP[isp*Dim + ic] = srcP[isp*Dim + ic];}
        }
      }
    }
  });

  if(gammai > -1){
    ga_matrices_cps ga_cps;
    std::vector<ga_M > gL;gL.resize(16);
    {int o=0;
    for(int i=0;i<6;i++){gL[o] = ga_cps.ga[0][i];o+=1;}
    for(int i=2;i<6;i++){gL[o] = ga_cps.ga[1][i];o+=1;}
    for(int i=3;i<6;i++){gL[o] = ga_cps.ga[2][i];o+=1;}
    for(int i=4;i<6;i++){gL[o] = ga_cps.ga[3][i];o+=1;}
    for(int i=5;i<6;i++){gL[o] = ga_cps.ga[4][i];o+=1;}}
    ga_M& ga = gL[gammai];
    prop4d_sink_gamma(res, ga );
  }
}


}

////#undef REDUCEFAC

#endif
