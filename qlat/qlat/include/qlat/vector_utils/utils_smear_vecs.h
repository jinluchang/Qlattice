// utils_smear_vecs.h
// Gen Wang
// Sep. 2021

#ifndef UTILS_SMEAR_VECS_H
#define UTILS_SMEAR_VECS_H
#pragma once

#include "general_funs.h"
#include <qlat/qcd-utils.h>
#include <qlat/qcd-prop.h>
#include <qlat/qcd-smear.h>
#include "utils_Vec_redistribute.h"
#include "utils_shift_vecs.h"
#include "utils_eo_copies.h"
#include "utils_check_fun.h"
#include "utils_field_operations.h"
#include "utils_field_expand.h"
 
namespace qlat{

///////////psrc order, bfac, c, d0, t,z,y,x
#ifdef QLAT_USE_ACC
template <class T, Int bfac, Int d0, Int dirL>
__global__ void gauss_smear_global4(T* pres, const T* psrc, const T* gf, const T bw, const T norm,
  const Long Nvol, const Long* map_bufD, const Long* map_final)
{

  const Int offB = dirL == 3? 1 : 3;  // for bank conflicts ?
  //__shared__ T ls[(dirL*2)*3*3 + offB];
  //__shared__ T ps[(dirL*2)*bfac*3*d0 + offB];
  __shared__ T ls[3][(2*dirL)*3 + offB];
  __shared__ T ps[d0*bfac][2*dirL*3 + offB];

  ///tid should larger than 9
  unsigned int   tid   =  threadIdx.x;
  unsigned long  count =  blockIdx.y*gridDim.x + blockIdx.x;
  unsigned int   ns = blockDim.x;
  const Int dir_max = 4;
  const Int nsites = bfac * 3 * d0;

  if(count < Nvol){

  ///////(2*dirL) --> c0 , c1
  ///////TODO need to check dirL is not 3
  unsigned int off = tid;
  {
    const T* gf_t = &gf[count*(2*dir_max)*9];
    //#pragma unroll
    while(off < (2*dirL)*9){
      //int dir = off/9;int c = off%9;
      //ls[(c/3)*(2*dirL)*3 + dir*3 + c%3] = gf[count*(2*dirL)*9 + off];off += ns;
      //ls[off] = gf_t[off];
      const Int dir = off/9;
      const Int c0 = (off/3)%3;
      const Int c1 =  off%3;
      const Int dirM = dir - dirL;
      //ls[c1*(2*dirL)*3 + dir*3 + c0 ] = gf_t[(dirM + dir_max)*9 + c1*3 + c0];
      ls[c1][dir*3 + c0 ] = gf_t[(dirM + dir_max)*9 + c1*3 + c0];
      off += ns;
    }
  }

  const Long res_off = map_final[count];
  const T* wo = &psrc[res_off*nsites];
  T* wm       = &pres[res_off*nsites];

  ///////(2*dir) --> bi, c1 , d0
  for (Int dir = -dirL; dir < dirL; ++dir){
    ////Long src_off = map_bufD[(dir+4)*Nvol + count];
    Long src_off = map_bufD[count* dirL*2 + (dir + dirL)];
    const T* src_t = &psrc[src_off*nsites];
    //T* res_t     = &ps[(dir+dirL)*bfac*3*d0];
    ////T* res_t     = &ps[dir+dirL];
    off = tid;
    //#pragma unroll
    while(off < nsites){
      //res_t[off] = src_t[off]; off+=ns;
      const unsigned int bi = off/(3*d0);
      const unsigned int c  = (off%(3*d0))/d0;
      const unsigned int di = off%d0;
      //ps[(bi*d0 + di)*(2*dirL*3) + (dir+dirL)*3 + c] = src_t[off];
      ps[(bi*d0 + di)][(dir+dirL)*3 + c] = src_t[off];
      off+=ns;
    }
  }
  __syncthreads();

  T tmp = 0.0;
  off = tid;
  while(off < nsites){
    //unsigned int c0 = off/(bfac*d0);
    //unsigned int bi = (off%(bfac*d0))/d0;
    //unsigned int di = off%d0;
    unsigned int bi = off/(3*d0);
    unsigned int c0 = (off%(3*d0))/d0;
    unsigned int di = off%d0;

    //T* a = &ls[c0*(2*dirL*3)];
    //T* b = &ps[(bi*d0 + di)*(2*dirL*3)];

    T* a = ls[c0];
    T* b = ps[(bi*d0 + di)];

    tmp = 0.0;

    //#pragma unroll
    for(Int dir=0;dir<(2*dirL)*3; dir++)
    {
      //tmp += b[dir*bfac*d0] * a[dir];
      tmp += b[dir] * a[dir];
    }

    //wm[(bi*3 + c0)*d0 + di] = norm*(wo[(bi*3 + c0)*d0 + di ] + bw*tmp);
    wm[off] = norm*(wo[off] + bw*tmp);
    off += ns;
  }

  }

}
#endif

inline void get_mapvq_each(const std::vector<CommPackInfo> &pack_infos, std::vector<qlat::vector<Long > >& mapvq, Int dir=0)
{
  std::vector<std::vector<Long > > mapv(2);//mapv.resize(4*pack_infos.size());
  ////std::vector<Long > factorL;
  for (Long i = 0; i < (Long)pack_infos.size(); ++i) {
    const CommPackInfo& cpi = pack_infos[i];
    //Long bufi = cpi.buffer_idx + 0;
    //Long fi   = cpi.offset     + 0;
    for(Long j=0; j< cpi.size ; j++)
    {
      ////factorL.push_back(cpi.size);
      Long bufi = cpi.buffer_idx + j;
      Long fi   = cpi.offset     + j;
      //if(dir == 1){
      //  qmessage("transfer ori %8ld, buffer %8ld \n", fi, bufi);
      //}
      ////qmessage("size %d, group %d, bufi %d, fi %d \n", int(i), int(cpi.size), int(bufi), int(fi));
      if(dir == 0){mapv[0].push_back(bufi);mapv[1].push_back(  fi);}
      if(dir == 1){mapv[0].push_back(  fi);mapv[1].push_back(bufi);}
    }
  }
  //std::vector<std::vector<Long > > mapv_copy(3);
  //sort_vectors_by_axis(mapv, mapv_copy, 1);
  mapvq.resize(2);
  for(Int j=0;j<2;j++)
  {
    mapvq[j].resize( mapv[j].size() );
    for(unsigned long i=0;i<mapv[j].size();i++)
    {
      mapvq[j][i] = mapv[j][i];
    }
  }

  //mapvq.resize(2);for(unsigned int i=0;i<2;i++){mapvq[i].copy_from(mapv[i], 1);}
}

////  === output
////  qlat::vector<Long > local_map_typeA0;  //// i to count, Nvol_ext, original positions to count new positions
////  qlat::vector<Long > local_map_typeA1;  //// reverse
////  qlat::vector<Long > map_index_typeA0;  //// index to index_typeA
////  qlat::vector<Long > map_index_typeA1;  //// from count to writi (new indexings)
////  qlat::vector<Long > map_index_typeAL;  //// count to index, loop mappings
////  qlat::vector<Long > map_bufD_typeA  ;  //// Nvol * dirL * 2, count*dirL*2 + (dir+dirL), direction indexings
inline void get_maps_hoppings(const Geometry& geo, const Geometry& geo_ext, const Int dirL,
  qlat::vector<Long >& local_map_typeA0,
  qlat::vector<Long >& local_map_typeA1,
  qlat::vector<Long >& map_index_typeA0,
  qlat::vector<Long >& map_index_typeA1,
  qlat::vector<Long >& map_index_typeAL,
  qlat::vector<Long >& map_bufD_typeA  ,
  std::vector< qlat::vector<Long > >& copy_extra_index,
  std::vector<Long >& pos_typeA   )
{
  TIMER("get_maps_hoppings");

  (void)local_map_typeA1;

  Qassert(dirL == 3 or dirL == 4);

  qlat::vector<Int> Nv,nv,mv;
  geo_to_nv(geo, nv, Nv, mv);
  const Long Nvol     = geo.local_volume();
  const Long Nvol_ext = geo_ext.local_volume_expanded();
  pos_typeA.resize(2);

  //std::vector<qlat::vector<Long > > map_bufV;
  qlat::vector<Long > map_bufD;

  /////const Int dir_max = 4;
  qlat::vector<Long > map_index_typeO_0;map_index_typeO_0.resize(Nvol_ext);
  qlat::vector<Long > map_index_typeO_1;map_index_typeO_1.resize(Nvol    );
  for(Long i=0;i<Nvol_ext;i++){map_index_typeO_0[i] = -1;}
  for(Long i=0;i<Nvol    ;i++){map_index_typeO_1[i] = -1;}

  //for(unsigned int i=0;i<map_bufV.size();i++){map_bufV[i].resize(Nvol       );}

  #pragma omp parallel for
  for(Long index=0;index<Nvol;index++){
    const Coordinate xl = geo.coordinate_from_index(index);
    Long index_ext = geo_ext.offset_from_coordinate(xl, 1);
    map_index_typeO_0[index_ext] = index    ;
    map_index_typeO_1[index    ] = index_ext;
    //if(index == 0){
    //  qmessage("==rank %3d, x y z t, %3d %3d %3d %3d, ext %ld \n", qlat::get_id_node(), xl[0], xl[1], xl[2], xl[3], index_ext);
    //}
    //qmessage("index ext %8ld index %8ld\n", index_ext, index);
    //map_bufV[0][index]  = geo_ext.offset_from_coordinate(xl, 1);
    //map_bufV[1][index]  = index;
    //////qmessage(" mappings %ld %ld \n",map_bufV[0][index], map_bufV[1][index]);
  }

  map_bufD.resize(Nvol*dirL*2);

  std::vector<Long > need_pos;need_pos.resize(Nvol_ext); ////for checkings only
  for(Long i=0;i<Long(need_pos.size());i++){need_pos[i] = 0;}

  //std::vector<Int > xlE;xlE.resize(8);
  for(Int dir=-dirL;dir<dirL; dir++){
  #pragma omp parallel for
  for(Long index=0;index<Nvol;index++)
  {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xl1 = coordinate_shifts(xl, dir);
    Long index_ext = geo_ext.offset_from_coordinate(xl1, 1);
    map_bufD[index*dirL*2 + (dir+dirL)] = index_ext;
    need_pos[ map_bufD[index*dirL*2 + (dir+dirL)] ] += 1;
    //if(map_index_typeO_0[index_ext] == -1)
    //{
    //  ////const Coordinate xl = geo_ext.coordinate_from_index(i);///this function is wrong!
    //  if(xl[3] == 0){
    //    xlE[xl1[2]] += 1;
    //    //qmessage("rank %3d, x y z t, %3d %3d %3d %3d, ext %3d %3d %3d %3d , index %ld ext %ld\n", qlat::get_id_node(), 
    //    //  xl[0], xl[1], xl[2], xl[3],
    //    //  xl1[0], xl1[1], xl1[2], xl1[3],
    //    //  index, index_ext
    //    //);
    //  }
    //}

  }}
  //for(Int i=0;i<8;i++)
  //{
  //  qmessage("i %3d, count %5d \n", i, xlE[i]);
  //}

  //Qassert(false)

  //Long Nneed = 0;
  //for(Long i=0;i<Long(need_pos.size());i++){
  //  if(need_pos[i] > 0){Nneed += 1;}
  //}
  //qmessage("rank %3d, Nneed %ld, total %ld \n", qlat::get_id_node(), Nneed - Nvol, Nvol_ext - Nvol );

  QLAT_PUSH_DIAGNOSTIC_DISABLE_DANGLING_REF;
  const CommPlan& plan = get_comm_plan(set_marks_field_1, "", geo_ext, 1);
  QLAT_DIAGNOSTIC_POP;

  std::vector<qlat::vector<Long > > mapvq_send;
  std::vector<qlat::vector<Long > > mapvq_recv;
  get_mapvq_each(plan.send_pack_infos, mapvq_send, 0);
  get_mapvq_each(plan.recv_pack_infos, mapvq_recv, 1);

  qlat::vector<Long > send_buffer_index;
  qlat::vector<Long > recv_buffer_index;

  ////Qassert(false);
  ////qmessage("rank %5d, send %ld, recv %ld \n", qlat::get_id_node(), Long(plan.total_send_size), Long(plan.total_recv_size));

  //{
  //  std::vector<Long > buf_s;
  //  std::vector<Long > buf_r;
  //  buf_s.resize(plan.total_send_size);
  //  buf_r.resize(plan.total_recv_size);
  //  for(unsigned long i=0;i<buf_s.size();i++){buf_s[i] = -1;}
  //  for(unsigned long i=0;i<buf_r.size();i++){buf_r[i] = -1;}

  //  Long count_send = 0;
  //  Long count_recv = 0;
  //  {
  //    ////TIMER("refresh_expanded-comm-init");
  //    for (size_t i = 0; i < plan.recv_msg_infos.size(); ++i) {
  //      const CommMsgInfo& cmi = plan.recv_msg_infos[i];
  //      count_recv += cmi.size;
  //      for(Long i=0;i<cmi.size;i++){buf_r[cmi.buffer_idx + i] = 1;}
  //    }
  //    for (size_t i = 0; i < plan.send_msg_infos.size(); ++i) {
  //      const CommMsgInfo& cmi = plan.send_msg_infos[i];
  //      count_send += cmi.size;
  //      for(Long i=0;i<cmi.size;i++){buf_s[cmi.buffer_idx + i] = 1;}
  //    }
  //  }
  //  //////qmessage("Actual send %ld, recv %ld \n", count_send, count_recv);
  //}


  std::vector<std::vector<Long> > check_send;check_send.resize(Nvol_ext);
  std::vector<std::vector<Long> > check_recv;check_recv.resize(Nvol_ext);
  //for(Long i=0;i<Nvol_ext;i++){
  //  check_send[i] = -1;check_recv[i] = -1;
  //}

  ///#pragma omp parallel for
  ////same data may be send to different nodex
  Long count_s = 0;
  Long count_r = 0;
  ////#pragma omp parallel for
  for(Long i=0;i<Long(mapvq_send[0].size());i++)
  {
    Long ri = mapvq_send[0][i];
    Long si = mapvq_send[1][i];
    ////send_buffer_index[ri] = si;
    //if(need_pos[si])
    //Qassert(check_send[si] == -1)
    check_send[si].push_back(ri);
    count_s += 1;
  }
  ////#pragma omp parallel for
  for(Long i=0;i<Long(mapvq_recv[0].size());i++)
  {
    Long ri = mapvq_recv[1][i];
    Long si = mapvq_recv[0][i];
    //recv_buffer_index[ri] = si;
    ////qmessage("%ld %ld, max %ld %ld \n", ri, si, Long(recv_buffer_index.size()), Nvol_ext);
    Qassert(check_send[si].size() == 0);
    check_recv[si].push_back(ri);
    count_r += 1;
  }
  //qmessage("plan write send %ld, recv %ld \n", count_s, count_r);
  ////Qassert(false)
  /////QLAT_VEC_CKPOINT

  qlat::FieldM<int8_t, 1> eo;
  qlat::qlat_map_eo_site(eo, geo);
  int8_t* eo_int8_t = (int8_t*) qlat::get_data(eo).data();

  std::vector<Long > local_map_typeA0_e;local_map_typeA0_e.resize(Nvol_ext); ////local_map0[count] == original positions
  std::vector<Long > local_map_typeA0_o;local_map_typeA0_o.resize(Nvol_ext); ////local_map0[count] == original positions
  //std::vector<Long > local_map_typeA0;
  //std::vector<Long > local_map_typeA1;
  //std::vector<Long > local_map_typeA2;local_map_typeA1.resize(Nvol_ext);
  #pragma omp parallel for
  for(Long i=0;i<Nvol_ext;i++){
    local_map_typeA0_e[i] = -1;
    local_map_typeA0_o[i] = -1;
  }

  //std::sort(check_send.begin(), check_send.end());
  //std::sort(check_recv.begin(), check_recv.end());

  //for(Long i=0;i<Nvol_ext;i++)
  //{
  //  if(check_send[i] != -1){
  //    qmessage("index ext %8ld, send %8ld \n", i, check_send[i]);
  //  }
  //}

  //for(Long i=0;i<Nvol_ext;i++)
  //{
  //  if(check_recv[i] != -1){
  //    qmessage("index ext %8ld, recv %8ld \n", i, check_recv[i]);
  //  }
  //}

  //QLAT_VEC_CKPOINT

  Long count_e = 0;
  Long count_o = 0;
  //std::vector<Int > xlC;xlC.resize(8);
  for(Long i=0;i<Nvol_ext;i++)
  {
    if(check_send[i].size() == 0 and check_recv[i].size() == 0)
    {
      Long index = map_index_typeO_0[i];
      //if(index == -1)
      //if(index == -1 and need_pos[i] == -1)
      //if(index == -1 and need_pos[i] > 0)
      //{
      //  xlC[0] += 1;
      //  //const Coordinate xl = geo_ext.coordinate_from_index(i);///this function is wrong!
      //  //if(xl[3] == 0){
      //  //  xlC[xl[2]] += 1;
      //  //  qmessage("rank %3d, x y z t, %3d %3d %3d %3d \n", qlat::get_id_node(), xl[0] - 1, xl[1] - 1, xl[2], xl[3]);
      //  //}
      //}
      Qassert(not (index == -1 and need_pos[i] > 0));
      if(index == -1){continue ;} ////corners 
      //#if GET_MAPS_DEBUG==1
      //////qmessage("index %ld %ld \n", index, i);
      //Qassert(index != -1);
      //#endif
      if(eo_int8_t[index] == 0){
        local_map_typeA0_e[i] = count_e;
        count_e += 1;
      }

      if(eo_int8_t[index] == 1){
        local_map_typeA0_o[i] = count_o;
        count_o += 1;
      }
      ////local_map1[i    ] = count;
      //count += 1;
    }
  }
  //for(Int i=0;i<8;i++)
  //{
  //  qmessage("i %3d, count %5d \n", i, xlC[i]);
  //}

  Long Nvol_sum = count_e + count_o + count_s + count_r;
  if(Nvol_sum < Nvol_ext){
    ////qmessage("==Check %ld %ld \n", Nvol_ext, Nvol_sum);
    Nvol_sum = Nvol_ext;
  }
  local_map_typeA0.resize(Nvol_sum); ////local_map0[count] == original positions
  //local_map_typeA1.resize(Nvol_sum);
  #pragma omp parallel for
  for(Long i=0;i<Nvol_sum;i++){
    local_map_typeA0[i]   = -1;
    //local_map_typeA1[i]   = -1;
  }

  //qmessage("Vol %ld, extra %ld, full %ld \n", Nvol, Nvol_ext, Nvol_sum);

  const Int use_eo = 0;

  ////write index with eo splitted
  Long count = 0;
  for(Long i=0;i<Nvol_ext;i++)
  {
    if(local_map_typeA0_e[i] != -1)
    {
      ////local_map_typeA0_e[i]
      Qassert(local_map_typeA0[i] == -1);
      //Qassert(local_map_typeA1[count] == -1);
      local_map_typeA0[i    ] = count;
      //local_map_typeA1[count] = i    ;
      count += 1;
    }

    if(use_eo == 0)
    if(local_map_typeA0_o[i] != -1)
    {
      Qassert(local_map_typeA0[i] == -1);
      //Qassert(local_map_typeA1[count] == -1);
      ////local_map_typeA0_e[i]
      local_map_typeA0[i    ] = count;
      //local_map_typeA1[count] = i    ;
      count += 1;
    }
  }

  if(use_eo == 1)
  for(Long i=0;i<Nvol_ext;i++)
  {
    if(local_map_typeA0_o[i] != -1)
    {
      Qassert(local_map_typeA0[i] == -1);
      //Qassert(local_map_typeA1[count] == -1);
      ////local_map_typeA0_e[i]
      local_map_typeA0[i    ] = count;
      //local_map_typeA1[count] = i    ;
      count += 1;
    }
  }
  //qmessage("total local %8ld, even %8ld, odd %8ld \n", count, count_e, count_o);

  Long pos_send = 0;
  Long pos_recv = 0;
  pos_send = count;

  Long max_pos = pos_send;
  std::vector<std::vector<Long > > copy_extra;////send inner boundaries copy
  //std::vector<Long > copy_extra0;////send inner boundaries copy
  //std::vector<std::vector<Long > > copy_extra1;////send inner boundaries copy

  copy_extra.resize(Nvol_ext);
  /////Lont count_copy = 0;
  for(Long i=0;i<Nvol_ext;i++)
  {
    /////if(check_send[i].size() != 0)
    for(unsigned long j=0;j < check_send[i].size();j++)
    {
      Long send = pos_send + check_send[i][j];
      copy_extra[i].push_back(send);
      if(send + 1 >= max_pos){max_pos = send + 1;}
      if(local_map_typeA0[i] == -1){
        local_map_typeA0[i] = send;
      }
      //else{
      //  Qassert(map_index_typeO_0[i] == -1);////must be boundaries
      //}
      /////local_map_typeA0[i    ] = send;
      //Qassert(local_map_typeA0[i] == -1);
      Qassert(send < Nvol_sum);
      //Qassert(local_map_typeA1[send] == -1);
      //local_map_typeA1[send ] = i;
      count += 1;
    }
  }
  ////Qassert(false);

  std::vector<Long > copy_extraL;////send inner boundaries copy
  std::vector<std::vector<Long > > copy_extra_pos;////send inner boundaries copy
  Long count_sum = 0;
  for(Long i=0;i<Nvol_ext;i++)
  { 
    if(copy_extra[i].size() > 1){
      copy_extraL.push_back(copy_extra[i][0]);
      std::vector<Long > pos;
      for(unsigned int j=1;j<copy_extra[i].size();j++)
      {
        pos.push_back(copy_extra[i][j]);
        count_sum += 1;
      }
      copy_extra_pos.push_back(pos);
    }
  }
  std::vector<Long > sortL = get_sort_index(copy_extraL.data(), copy_extraL.size());
  ////std::vector< qlat::vector<Long > > copy_extra_index;
  copy_extra_index.resize(3);

  qlat::vector<Long >& copy_extra_index0 = copy_extra_index[0];
  qlat::vector<Long >& copy_extra_index1 = copy_extra_index[1];
  qlat::vector<Long >& copy_extra_posL   = copy_extra_index[2];
  copy_extra_index0.resize(sortL.size());
  copy_extra_index1.resize(sortL.size() + 1);
  copy_extra_posL.resize(count_sum);

  //for(unsigned long isp=0;isp<copy_extra_index0.size();isp++){copy_extra_index[isp] = copy_extra_index0[isp];}
  //for(unsigned long isp=0;isp<copy_extra_posL0.size() ;isp++){copy_extra_posL[ isp] = copy_extra_posL0[ isp];}

  Long counti = 0;
  for(unsigned long si=0;si<sortL.size();si++)
  {
    Long li = sortL[si];
    Long send = copy_extraL[li];
    std::vector<Long >& pos = copy_extra_pos[li];
    copy_extra_index0[si ] = send;
    copy_extra_index1[si ] = counti;
    //qmessage("copy %ld, c0 %ld ", send, counti);
    for(unsigned long j=0;j<pos.size();j++)
    {
      copy_extra_posL[counti] = pos[j];
      counti += 1;
    }
    //qmessage(" c1 %ld \n", counti);
  }
  copy_extra_index1[sortL.size() ] = counti;
  /////Qassert(false);

  //Long counti = 0;
  //std::vector<Long > copy_extra_index0;
  //std::vector<Long > copy_extra_posL0;
  //Long counti = 0;
  //for(Long i=0;i<Nvol_ext;i++)
  //{
  //  if(copy_extra[i].size() > 1){
  //    copy_extra_index0.push_back(copy_extra[i][0]);///index of the initial value
  //    ////copy_extra_index0.push_back(copy_extra[i].size() - 1);
  //    copy_extra_index0.push_back(counti);////index of initial positions
  //    qmessage("copy %ld, c0 %ld ", copy_extra[i][0], counti);
  //    for(unsigned int j=1;j<copy_extra[i].size();j++)
  //    {
  //      copy_extra_posL0.push_back(copy_extra[i][j]);
  //      counti += 1;
  //    }
  //    copy_extra_index0.push_back(counti);////index of final positions
  //    qmessage(" c1 %ld \n", counti);
  //  }
  //}

  Qassert(max_pos >= count);////some may be empty
  ////qmessage("actual send %ld, max send %ld \n", count - pos_send, max_pos - pos_send);
  pos_recv = max_pos;
  //Long count0 = count;

  for(Long i=0;i<Nvol_ext;i++)
  {
    ////if(check_recv[i] != -1)
    for(unsigned long j=0;j < check_recv[i].size();j++)
    {
      Long recv = pos_recv + check_recv[i][j];
      if(recv + 1 >= max_pos){max_pos = recv + 1;}
      /////local_map_typeA0[i    ] = recv;
      Qassert(local_map_typeA0[i] == -1);
      //Qassert(local_map_typeA1[recv] == -1);
      //if(local_map_typeA0[i] == -1){
      //  local_map_typeA0[i] = recv;
      //}else{
      //  Qassert(map_index_typeO_0[i] == -1);////must be boundaries
      //}
      local_map_typeA0[i] = recv;
      //local_map_typeA1[recv ] = i    ;
      count += 1;
    }
  }

  Qassert(max_pos >= count);////some may be empty
  /////qmessage("actual recv %ld, max recv %ld \n", count - pos_recv, max_pos - pos_recv);

  pos_typeA[0] = pos_send;
  pos_typeA[1] = pos_recv;

  ////Long count_end = max_pos;
  //for(Long i=pos_send;i<pos_recv;i++){ 
  //  if(local_map_typeA1[i] != -1){
  //  Qassert(local_map_typeA0[local_map_typeA1[i]] < pos_recv );}
  //}
  //for(Long i=pos_recv;i<count_end;i++){
  //  if(local_map_typeA1[i] != -1){
  //  Qassert(local_map_typeA0[local_map_typeA1[i]] < count_end);}
  //}
  //for(Long i=0;i<count_end;i++){
  //  if(local_map_typeA0[i] == -1 or local_map_typeA1[i] == -1)
  //  {
  //    qmessage("layout wrong %ld !\n", i);
  //  }
  //}
  //qmessage("send_pos %ld %ld recv_pos %ld %ld, final %ld %ld %ld, vol %ld \n", 
  //  pos_send, pos_recv - pos_send, pos_recv, count_end - pos_recv, count0, count, count_end, Nvol_ext);
  ////Qassert(false);

  ////map to original vol index
  //qlat::vector<Long > map_index_typeA0;
  //qlat::vector<Long > map_index_typeA1;
  //qlat::vector<Long > map_index_typeAL;
  map_index_typeA0.resize(Nvol);
  map_index_typeA1.resize(Nvol_sum);
  map_index_typeAL.resize(Nvol_sum);
  #pragma omp parallel for
  for(Long i=0;i<Nvol;i++){    map_index_typeA0[i] = -1;}
  #pragma omp parallel for
  for(Long i=0;i<Nvol_sum;i++){map_index_typeA1[i] = -1;}
  #pragma omp parallel for
  for(Long i=0;i<Nvol_sum;i++){    map_index_typeAL[i] = -1;}

  #pragma omp parallel for
  for(Long index = 0;index  < Nvol; index++){
    Long index_ext   = map_index_typeO_1[index];
    Long index_typeA = local_map_typeA0[index_ext]; ////memory positions in typeA format
    Qassert(index_typeA != -1);

    Qassert(map_index_typeA0[index      ] == -1);
    Qassert(map_index_typeA1[index_typeA] == -1);
    map_index_typeA0[index      ] = index_typeA;
    map_index_typeA1[index_typeA] = index;
  }

  count = 0;
  for(Long index_typeA = 0;index_typeA  < Nvol_ext; index_typeA++){
    Long index = map_index_typeA1[index_typeA];
    if(index != -1){
      map_index_typeAL[count] = index;
      count += 1;
      Qassert(count <= Nvol);
    }
  }
  //qmessage("indexed count %8ld, vol %8ld \n", count, Nvol);

  map_index_typeA1.resize(Nvol);////reuse map_index_typeA1

  //////map_bufD[index*dirL*2 + (dir+dirL)] --> needed Nvol_ext positions
  //qlat::vector<Long > map_bufD_typeA;
  map_bufD_typeA.resize(Nvol * dirL * 2 );
  #pragma omp parallel for
  for(Long i=0;i<Long(map_bufD_typeA.size());i++){map_bufD_typeA[i] = -1;}
  ////count = index_typeAL

  #pragma omp parallel for
  for(Long count = 0;count < Nvol;count++){
    Long index = map_index_typeAL[count];
    map_index_typeA1[count] = map_index_typeA0[index];
    ////qmessage("=count %8ld, iwrite %8ld \n", count, map_index_typeA0[index]);
    ////Long index_typeA = map_index_typeA1[index];
    for(Int dir=-dirL;dir<dirL; dir++)
    {
      Long index_ext = map_bufD[index*dirL*2 + (dir+dirL)]; ////memory postions in index_ext
      Long index_typeA = local_map_typeA0[index_ext];
      Qassert(index_typeA != -1);
      map_bufD_typeA[count*dirL*2 + (dir+dirL)] = index_typeA;
    }
  }

  #pragma omp parallel for
  for(Long i=0;i<Nvol * dirL * 2;i++){Qassert(map_bufD_typeA[i] != -1);}
  //count = 0;
  //for(Long i=0;i<Nvol_ext;i++){
  //  if(map_index_typeA1[i]!=-1){count += 1;}
  //}
  //Qassert(count == Nvol);

  count = 0;
  for(Long i=0;i<Nvol;i++){
    if(map_index_typeA0[i]!=-1){count += 1;}
  }
  Qassert(count == Nvol);

}

////Ty must be complex
template <typename Ty >
struct smear_fun{

  Geometry geo;
  Geometry geo_ext;
  Int dirL;
  bool smear_in_time_dir;

  ////shift_vec *svec;
  Vec_redistribute *vec_rot;
  move_index mv_idx;

  /////for shifters and geo
  fft_desc_basic fd;
  fft_desc_basic fd_new;

  //std::vector<qlat::vector_gpu<Long > > mapvq_send;
  //std::vector<qlat::vector_gpu<Long > > mapvq_recv;
  ////int CONTINUS;
  qlat::vector<MPI_Request> send_reqs;
  qlat::vector<MPI_Request> recv_reqs;
  ///unsigned int bfac;
  ////unsigned int Tsize;
  Long Nvol;
  Long Nvol_ext;

  //std::vector<qlat::vector<Long > > map_bufV;
  //qlat::vector<Long > map_bufD;
  qlat::vector<Int> Nv,nv,mv;

  qlat::vector<Long > map_index_typeI ;
  qlat::vector<Long > map_index_typeA0;
  qlat::vector<Long > map_index_typeA1;
  qlat::vector<Long > map_index_typeAL;
  qlat::vector<Long > map_bufD_typeA  ;

  Int use_gauge_mapping;////hack for box smearings, gaussian is 1

  std::vector< qlat::vector<Long > > copy_extra_index;
  std::vector<Long> pos_typeA;

  unsigned int NVmpi;
  Int groupP;

  Int redistribution_copy;

  qlat::vector_gpu<Ty > send_buffer;
  qlat::vector_gpu<Ty > recv_buffer;
  //qlat::vector_gpu<Ty > gauge_buf;
  //std::string prop_name;
  std::string prop_buf0_name;
  std::string prop_buf1_name;
  std::string gauge_buf_name_base;
  std::string gauge_buf_name;
  Int gauge_dup;
  //qlat::vector_gpu<Ty > prop;
  //qlat::vector_gpu<Ty > prop_buf;
  std::vector<qlat::vector_gpu<Ty > > vL;

  bool  gauge_setup_flag;
  void* gauge_check;
  crc32_t gauge_checksum;
  /////buffers

  bool  mem_setup_flag;
  qlat::vector<qlat::ComplexD > mom_factors;

  //template <class Td>
  //void setup(const GaugeFieldT<Td >& gf, const CoordinateD& mom, const bool smear_in_time_dir){
  //  qlat::vector<qlat::ComplexD > momF(8);
  //  for (Int i = 0; i < 8; ++i) {
  //    const Int dir = i - 4;
  //    const double phase = dir >= 0 ? mom[dir] : -mom[-dir - 1];
  //    momF[i] = std::polar(1.0, -phase);
  //  }
  //  gauge_setup(gf, momF);
  //}

  smear_fun(const Geometry& geo_, const bool smear_in_time_dir_){
    Nvol  = 0;
    redistribution_copy = 0;
    NVmpi = 0;
    groupP = 0;

    prop_buf0_name = std::string("Smear_prop_buf0");
    prop_buf1_name = std::string("Smear_prop_buf1");
    gauge_buf_name_base = std::string("Gauge_buf_name");
    gauge_dup = -1;

    gauge_check = NULL;
    gauge_checksum = 0;
    gauge_setup_flag = false;

    mem_setup_flag = false;
    dirL = 3;

    fft_desc_basic fd_(geo_);
    copy_fft_desc(fd, fd_);
    copy_fft_desc(fd_new, fd_);
    NVmpi = fd.mz*fd.my*fd.mx;

    smear_in_time_dir = smear_in_time_dir_;
    init_mem(geo_, smear_in_time_dir);

    ////Box smearing buffer size
    vL.resize(8);
    vec_rot = NULL;
    use_gauge_mapping = -1;
    /////CONTINUS = 1;
  }

  inline void check_setup(){
    if(Nv.size() == 0){Qassert(false);}
    if(Nvol == 0){Qassert(false);}
    if(Nvol_ext == 0){Qassert(false);}
    const Int GPU = 1;
    qlat::vector_gpu<Ty >& gauge_buf = get_vector_gpu_plan<Ty >(0, gauge_buf_name, GPU);
    if(gauge_buf.size() == 0){Qassert(false);}
    if(gauge_check == NULL){Qassert(false);}
    ////if(gauge_checksum == 0){Qassert(false);}
    ////if(map_bufV.size() == 0){Qassert(false);}
    ////if(map_bufD.size() == 0){Qassert(false);}
    if(vL.size() != 8){Qassert(false);}
    if(mem_setup_flag == false){Qassert(false);}
    ////if(gauge.size() != Nvol*4*2*9){Qassert(false);}
  };

  inline void get_mapvq(const std::vector<CommPackInfo> &pack_infos, std::vector<qlat::vector_gpu<Long > >& mapvq, Int dir=0)
  {
    std::vector<std::vector<Long > > mapv(3);//mapv.resize(4*pack_infos.size());
    ///std::vector<Long > factorL;
    for (Long i = 0; i < (Long)pack_infos.size(); ++i) {
      const CommPackInfo& cpi = pack_infos[i];
      Long bufi = cpi.buffer_idx + 0;
      Long fi   = cpi.offset     + 0;

      if(dir == 0){mapv[0].push_back(bufi);mapv[1].push_back(  fi);mapv[2].push_back(  cpi.size);}
      if(dir == 1){mapv[0].push_back(  fi);mapv[1].push_back(bufi);mapv[2].push_back(  cpi.size);}

      //for(Long j=0; j< cpi.size ; j++)
      //{
      //  ////factorL.push_back(cpi.size);
      //  Long bufi = cpi.buffer_idx + j;
      //  Long fi   = cpi.offset     + j;
      //  qmessage("size %d, group %d, bufi %d, fi %d \n", int(i), int(cpi.size), int(bufi), int(fi));
      //  if(dir == 0){mapv[0].push_back(bufi);mapv[1].push_back(  fi);}
      //  if(dir == 1){mapv[0].push_back(  fi);mapv[1].push_back(bufi);}
      //}
    }

    //for(Long j=0;j<mapv[0].size();j++){
    //  qmessage("%8ld %8ld %8ld \n", mapv[0][j], mapv[1][j], mapv[2][j]);
    //}

    //std::vector<Long > index = get_sort_index(&mapv[1][0], mapv[1].size());

    std::vector<std::vector<Long > > mapv_copy(3);
    sort_vectors_by_axis(mapv, mapv_copy, 1);

    //qmessage("\n\n\n");
    //for(Long j=0;j<mapv[0].size();j++){
    //  qmessage("%8ld %8ld %8ld \n", mapv_copy[0][j], mapv_copy[1][j], mapv_copy[2][j]);
    //}
 
    ////how to check continus
    //Long maxF = *max_element(std::begin(factorL), std::end(factorL)); 
    //int C0 = 1;
    //for(Int fi=maxF;fi>0;fi--){
    //  Int wrong = 0;
    //  for(unsigned int f=0;f<factorL.size();f++){
    //    if(factorL[f] % fi != 0){wrong = 1;} 
    //  }
    //  if(wrong == 0){C0 = fi;break;}
    //}
    //qmessage("continus %d \n", C0 );
    //abort_r();
    //return C0;

    mapvq.resize(3);for(unsigned int i=0;i<3;i++){mapvq[i].copy_from(mapv_copy[i], 1);}
  }

  /////index for MPI buffers
  inline void init_mem(const Geometry& geo_, bool smear_in_time_dir_ = false)
  {
    if(geo == geo_ and smear_in_time_dir == smear_in_time_dir_ and mem_setup_flag == true){return ;}
    geo = geo_;
    ////move to default geo
    geo.resize(Coordinate(0, 0, 0, 0), Coordinate(0, 0, 0, 0));
    //geo.multiplicity = 1;
    geo.eo=0;
    Geometry geo1 = geo;
    if(!smear_in_time_dir){geo1.resize(Coordinate(1, 1, 1, 0), Coordinate(1, 1, 1, 0));}
    if( smear_in_time_dir){geo1.resize(Coordinate(1, 1, 1, 1), Coordinate(1, 1, 1, 1));}
    geo_ext = geo1;
    QLAT_PUSH_DIAGNOSTIC_DISABLE_DANGLING_REF;
    const CommPlan& plan = get_comm_plan(set_marks_field_1, "", geo_ext, 1);
    QLAT_DIAGNOSTIC_POP;
    //bfac  = bfac_a;
    ////Tsize = Tsize_a;

    geo_to_nv(geo, nv, Nv, mv);
    Nvol     = geo.local_volume();
    Nvol_ext = geo_ext.local_volume_expanded();
    //Nv.resize(4);nv.resize(4);mv.resize(4);
    //for(Int i=0;i<4;i++){Nv[i]=geo.node_site[i];nv[i] = geo.node_site[i] * geo.geon.size_node[i];}
    //for(Int i=0;i<4;i++){mv[i] = nv[i]/Nv[i];}
    if(!smear_in_time_dir){dirL = 3;}
    if( smear_in_time_dir){dirL = 4;}

    //get_mapvq(plan.send_pack_infos, mapvq_send, 0);
    //get_mapvq(plan.recv_pack_infos, mapvq_recv, 1);

    qlat::vector<Long > local_map_typeA0;
    qlat::vector<Long > local_map_typeA1;
    get_maps_hoppings(geo, geo_ext, dirL,
      local_map_typeA0, local_map_typeA1, map_index_typeA0, map_index_typeA1, map_index_typeAL, map_bufD_typeA, copy_extra_index, pos_typeA);

    map_index_typeI.resize(Nvol);
    for(Long i=0;i<Long(Nvol);i++){map_index_typeI[i] = i;}

    //if(C0 != C1){abort_r("QLAT SMEARING CONTINUES FAILED!");}
    //CONTINUS = C0;

    send_reqs.resize(plan.send_msg_infos.size());
    recv_reqs.resize(plan.recv_msg_infos.size());

    /////const Int dir_max = 4;
    //map_bufV.resize(2);
    //for(unsigned int i=0;i<map_bufV.size();i++){map_bufV[i].resize(Nvol       );}
    //#pragma omp parallel for
    //for(Long index=0;index<Nvol;index++){
    //  const Coordinate xl = geo.coordinate_from_index(index);
    //  map_bufV[0][index]  = geo_ext.offset_from_coordinate(xl, 1);
    //  map_bufV[1][index]  = index;
    //  //////qmessage(" mappings %ld %ld \n",map_bufV[0][index], map_bufV[1][index]);
    //}

    //////abort_r();

    //map_bufD.resize(Nvol*dirL*2);
    //for(Int dir=-dirL;dir<dirL; dir++)
    //#pragma omp parallel for
    //for(Long index=0;index<Nvol;index++)
    //{
    //  const Coordinate xl = geo.coordinate_from_index(index);
    //  const Coordinate xl1 = coordinate_shifts(xl, dir);
    //  map_bufD[index*dirL*2 + (dir+dirL)] = geo_ext.offset_from_coordinate(xl1, 1);
    //}

    mem_setup_flag = true;
  }

  /////define the redistributed smearing kernels
  inline void init_distribute()
  {
    TIMERA("Gauge init_distribute");
    Qassert(dirL==3);
    Qassert(gauge_setup_flag);
    if(redistribution_copy == 1){return ;}

    ///////Redistribute the data
    ////if(int(NVmpi) < bfac){return ;}

    check_setup();
    //groupP = (bfac+NVmpi-1)/NVmpi;
    //qmessage("====Vec setup, NVmpi %d, groupP %d \n", NVmpi, groupP);
    const Int GPU = 1;
    qlat::vector_gpu<Ty >& gauge_buf = get_vector_gpu_plan<Ty >(0, gauge_buf_name, GPU);

    vec_rot = new Vec_redistribute(fd);
    /////update gf to each MPI core
    qlat::vector_gpu<Ty > gfT;gfT.resize(NVmpi*gauge_buf.size());
    qlat::vector_gpu<Ty > gfT_buf;gfT_buf.resize(NVmpi*gauge_buf.size());
    const Int dir_max = 4;
    const size_t Ndata = gauge_buf.size();

    for(Long vi=0;vi<NVmpi;vi++){cpy_data_thread(&(gfT.data()[vi*Ndata]), gauge_buf.data(), Ndata);}

    vec_rot->reorder(gfT.data(), gfT_buf.data(), 1, (dir_max*2)*9 ,   0);
    gauge_buf.copy_from(gfT);

    Nvol     = Nv[3] * nv[2]*nv[1]*nv[0];
    Nvol_ext = Nv[3] * nv[2]*nv[1]*nv[0];

    ///===new varaibles
    //local_map_typeA0.resize(Nvol);
    //local_map_typeA1.resize(Nvol);
    map_index_typeA0.resize(Nvol);
    map_index_typeA1.resize(Nvol);
    map_index_typeAL.resize(Nvol);
    for(Long i=0;i<Long(Nvol);i++)
    {
      //local_map_typeA0[i] = i;
      //local_map_typeA1[i] = i;
      map_index_typeA0[i] = i;
      map_index_typeA1[i] = i;
      map_index_typeAL[i] = i;
    }

    copy_extra_index.resize(3);
    map_bufD_typeA.resize(Nvol*dirL*2);
    ///===new varaibles

    std::vector<Int > Nn;Nn.resize(4);
    for(Int i=0;i<3;i++){Nn[i] = nv[i];}
    Nn[3] = Nv[3];

    //map_bufV.resize(2);
    //for(unsigned int i=0;i<map_bufV.size();i++){map_bufV[i].resize(Nvol       );}
    //#pragma omp parallel for
    //for(Long index=0;index<Nvol;index++){
    //  map_bufV[0][index] = index;
    //  map_bufV[1][index] = index;
    //}

    ////map_bufD.resize(Nvol*dirL*2);
    for(Int dir=-dirL;dir<dirL; dir++)
    #pragma omp parallel for
    for(Long index=0;index<Nvol;index++)
    {
      std::vector<Int > xl;xl.resize(4);
      xl[3] =  index/(Nn[2]*Nn[1]*Nn[0]);
      xl[2] = (index%(Nn[2]*Nn[1]*Nn[0]))/(Nn[1]*Nn[0]);
      xl[1] = (index/(Nn[0]))%Nn[1];
      xl[0] =  index%(Nn[0]);

      Int di = 0;
      if(dir >= 0){di=dir   ;xl[di] = (xl[di]+Nn[di]+1)%Nn[di];}
      if(dir <  0){di=-dir-1;xl[di] = (xl[di]+Nn[di]-1)%Nn[di];}

      //int ti =  index/(nv[2]*nv[1]*nv[0]);
      //int zi = (index%(nv[2]*nv[1]*nv[0]))/(nv[1]*nv[0]);
      //int yi = (index/(nv[0]))%nv[1];
      //int xi =  index%(nv[0]);
      //const Coordinate xl(xi,yi,zi,ti);
      //const Coordinate xl1 = coordinate_shifts(xl, dir);

      ///map_bufD[(dir+4)*Nvol + index] = ((xl[3]*Nn[2]+xl[2])*Nn[1] + xl[1])*Nn[0] + xl[0];
      //map_bufD[index*dirL*2 + (dir+dirL)] = geo_ext.offset_from_coordinate(xl1, 1);
      const Long pos = ((xl[3]*Nn[2]+xl[2])*Nn[1] + xl[1])*Nn[0] + xl[0]; 
      //map_bufD[index*dirL*2 + (dir+dirL)] = pos;
      map_bufD_typeA[index*dirL*2 + (dir+dirL)] = pos;
    }

    //geo_ext.resize(Coordinate(0, 0, 0, 0), Coordinate(0, 0, 0, 0));

    ///////fd update for box smearing
    desc_xyz_in_one(fd_new, geo);

    redistribution_copy = 1;
  }

  //  on CPU or Nvidia GPU, could avoid copy data to buffer
  //  on GPU AMD hip, somehow necessary
  template<int bfac>
  void refresh_expanded_GPU(Ty* f, Int GPU = 1)
  {
    (void)GPU;
    check_setup();
    if(redistribution_copy == 1){return ;}
    QLAT_PUSH_DIAGNOSTIC_DISABLE_DANGLING_REF;
    const CommPlan& plan = get_comm_plan(set_marks_field_1, "", geo_ext, 1);
    QLAT_DIAGNOSTIC_POP;
    const Long total_bytes =
        (plan.total_recv_size + plan.total_send_size) * bfac * sizeof(Ty);
    if (0 == total_bytes) {return;}
    TIMER_FLOPS("refresh_expanded");
    timer.flops += total_bytes / 2;

    #ifdef QLAT_USE_ACC
    const Long Nvol = map_index_typeAL.size();
    send_buffer.resizeL((pos_typeA[1] - pos_typeA[0])*bfac);
    recv_buffer.resizeL((Nvol - pos_typeA[1])*bfac);
    #endif

    ////copy extra send buf
    //copy_extra_index
    qlat::vector<Long >& copy_extra_index0 = copy_extra_index[0];
    qlat::vector<Long >& copy_extra_index1 = copy_extra_index[1];
    qlat::vector<Long >& copy_extra_posL   = copy_extra_index[2];
    //copy_extra_index0.resize(sortL.size());
    //copy_extra_index1.resize(sortL.size());
    //copy_extra_posL.resize(count_sum);

    Ty* s_buf = &f[pos_typeA[0] * bfac];
    Ty* r_buf = &f[pos_typeA[1] * bfac];

    Ty* S_buf = s_buf;
    Ty* R_buf = r_buf;

    const Long shift = pos_typeA[0];
    #ifdef QLAT_USE_ACC
    Ty* s_copy = send_buffer.data();
    cpy_GPU(send_buffer.data(), s_buf, send_buffer.size());

    S_buf = send_buffer.data();
    R_buf = recv_buffer.data();
    #else
    Ty* s_copy = &f[shift*bfac];
    #endif

    const Long Ncopy = copy_extra_index0.size();
    if(Ncopy != 0)
    {TIMERA("smear Copy extra");
    qacc_for(isp, Ncopy, {
      const Long send   = copy_extra_index0[isp] - shift;
      const Long count0 = copy_extra_index1[isp];
      const Long count1 = copy_extra_index1[isp + 1];
      const Long loop = count1 - count0;
      const Long* cP  = &copy_extra_posL[count0];
      Ty buf[bfac];
      for(Long bi=0;bi<bfac;bi++){buf[bi] = s_copy[send*bfac + bi];}

      //////Ty buf = f[send*bfac + bi];
      for(Long j=0;j<loop;j++)
      for(Long bi=0;bi<bfac;bi++)
      {
        ////Qassert(cP[j] >= pos_typeA[0] and cP[j] < pos_typeA[1]);
        s_copy[(cP[j]-shift)*bfac + bi] = buf[bi];
      }
    });
    }

    //cpy_data_from_group(send_buffer.data(), f , mapvq_send[0].data(), mapvq_send[1].data(), mapvq_send[2].data(), mapvq_send[0].size(), bfac, GPU);
    {
      TIMER_FLOPS("refresh_expanded-comm");
      timer.flops +=
          (plan.total_recv_size + plan.total_send_size) * bfac * sizeof(Ty) / 2;
      {
        ////TIMER("refresh_expanded-comm-init");
        const Int mpi_tag = QLAT_VECTOR_UTILS_MPI_TAG;
        for (size_t i = 0; i < plan.recv_msg_infos.size(); ++i) {
          const CommMsgInfo& cmi = plan.recv_msg_infos[i];
          Long count = cmi.size * bfac * sizeof(Ty) / sizeof(RealD);
          //qmessage("==recv comm %ld\n", count);
          //MPI_Irecv(&recv_buffer[cmi.buffer_idx*bfac], count, MPI_DOUBLE,
          //          cmi.id_node, mpi_tag, get_comm(), &recv_reqs[i]);

          MPI_Irecv(&R_buf[cmi.buffer_idx*bfac], count, MPI_DOUBLE,
                    cmi.id_node, mpi_tag, get_comm(), &recv_reqs[i]);
        }
        for (size_t i = 0; i < plan.send_msg_infos.size(); ++i) {
          const CommMsgInfo& cmi = plan.send_msg_infos[i];
          Long count = cmi.size * bfac * sizeof(Ty) / sizeof(RealD);
          //qmessage("==send comm %ld\n", count);
          //MPI_Isend(&send_buffer[cmi.buffer_idx*bfac], count, MPI_DOUBLE,
          //          cmi.id_node, mpi_tag, get_comm(), &send_reqs[i]);

          MPI_Isend(&S_buf[cmi.buffer_idx*bfac], count, MPI_DOUBLE,
                    cmi.id_node, mpi_tag, get_comm(), &send_reqs[i]);
        }
      }
      MPI_Waitall(recv_reqs.size(), recv_reqs.data(), MPI_STATUS_IGNORE);
      MPI_Waitall(send_reqs.size(), send_reqs.data(), MPI_STATUS_IGNORE);
    }

    #ifdef QLAT_USE_ACC
    cpy_GPU(r_buf, recv_buffer.data(), recv_buffer.size());
    #endif

    ////sync_node();
    //cpy_data_from_group(f, recv_buffer.data(), mapvq_recv[0].data(), mapvq_recv[1].data(), mapvq_recv[2].data(), mapvq_recv[0].size(), bfac, GPU);
  }

  inline void gauge_buf_name_append(Int buf)
  {
    if(gauge_dup == buf or buf == -1){return ;}
    //char tmp[100];
    std::string tmp = ssprintf("_%d", buf);
    gauge_buf_name_base += std::string( tmp );
    gauge_dup = buf;
    //qlat::vector_gpu<Ty >& gauge_buf = get_vector_gpu_plan<Ty >(0, gauge_buf_name, GPU);
  }

  //void gauge_setup(qlat::vector_gpu<T >& gfE, const GaugeFieldT<Tg >& gf,
  //           const qlat::vector<qlat::ComplexD >& momF = qlat::vector<qlat::ComplexD >(), bool force_update = false);
  template <class Td>
  void gauge_setup(const GaugeFieldT<Td >& gf, const CoordinateD& mom_, const Int force_update = 0, const Int use_gauge_mapping_ = 1)
  {
    Qassert(mem_setup_flag == true);
    qlat::vector<qlat::ComplexD > momF(8);
    for (Int i = 0; i < 8; ++i) {
      const Int dir = i - 4;
      const double phase = dir >= 0 ? mom_[dir] : -mom_[-dir - 1];
      momF[i] = std::polar(1.0, -phase);
    }

    const Int GPU = 1;
    qlat::vector_gpu<Ty >& gauge_buf = get_vector_gpu_plan<Ty >(0, gauge_buf_name, GPU);

    bool update = false;
    if(gauge_setup_flag == false){ update = true;}
    if(force_update == 1){ update = true;}
    if(gauge_buf.size() == 0){  update = true;}
    if(gauge_check != (void*) qlat::get_data(gf).data()){update = true;}
    if(mom_factors.size() != 8){update = true;}
    else{for(Int i=0;i<momF.size();i++){if(momF[i]!=mom_factors[i]){update = true;}}}
    if(use_gauge_mapping != use_gauge_mapping_){update = true;}

    ////for_update == -1; do not check gauge sum

    //const Long gf_data_size = get_expanded_data_size(gf) / sizeof(Td);
    const Long gf_data_size = GetFieldSize(gf) / sizeof(Td);
    if(force_update == 0){
      crc32_t tmp_gauge_checksum = quick_checksum((Td*) qlat::get_data(gf).data(), gf_data_size );
      if(gauge_checksum != tmp_gauge_checksum ){update = true;}
    }

    if(update){
      TIMERA("gauge setup");
      crc32_t tmp_gauge_checksum = quick_checksum((Td*) qlat::get_data(gf).data(), gf_data_size );
      ////qmessage("==Gauge checksum %X", &tmp_gauge_checksum);

      ///clear previous cache
      {
        qlat::vector_gpu<Ty >& gauge_buf = get_vector_gpu_plan<Ty >(0, gauge_buf_name, GPU);
        gauge_buf.resize(0);
      }

      ////rename gauge_buf_name
      //std::sscanf(gauge_buf_name.c_str(), "%s_%X", gauge_buf_name_base.c_str(), &tmp_gauge_checksum);
      std::string tem = ssprintf("%X", tmp_gauge_checksum);
      gauge_buf_name = gauge_buf_name_base + std::string("_") + std::string(tem);
      qlat::vector_gpu<Ty >& gauge_buf = get_vector_gpu_plan<Ty >(0, gauge_buf_name, GPU);

      Qassert(geo.total_site() == gf.geo().total_site());
      gauge_buf.resizeL(2*4* gf.geo().local_volume() * 9);
      mom_factors = momF;
      if(use_gauge_mapping_ == 1){
        extend_links_to_vecs(gauge_buf.data(), gf, momF, map_index_typeAL);
      }else{
        extend_links_to_vecs(gauge_buf.data(), gf, momF);
      }

      use_gauge_mapping = use_gauge_mapping_;
      gauge_setup_flag = true;
      gauge_check = (void*) qlat::get_data(gf).data();
      gauge_checksum = tmp_gauge_checksum;
      ///Need to redistribute when copied
      redistribution_copy = 0 ;
    }
    ////qmessage("==Gauge buf name %s \n", gauge_buf_name.c_str());
  }

  inline void clear_mem(){
    #ifdef __CLEAR_SMEAR_MEM__
    const Int GPU = 1;
    safe_free_vector_gpu_plan<Ty >(prop_buf0_name, GPU);
    safe_free_vector_gpu_plan<Ty >(prop_buf1_name, GPU);
    for(Int i=0;i<vL.size();i++){vL[i].resize(0);}
    mv_idx.free_mem();
    #endif
  }

  inline void clear_buf_mem(){
    const Int GPU = 1;
    ////clear gauges and need to reload
    qlat::vector_gpu<Ty >& gauge_buf = get_vector_gpu_plan<Ty >(0, gauge_buf_name, GPU);
    gauge_buf.resize(0);

    qlat::vector_gpu<Ty >& prop_buf0 = get_vector_gpu_plan<Ty >(0, prop_buf0_name, GPU);
    qlat::vector_gpu<Ty >& prop_buf1 = get_vector_gpu_plan<Ty >(0, prop_buf1_name, GPU);
    prop_buf0.resize(0);
    prop_buf1.resize(0);
  }

  inline void clear_all_mem(){
    clear_mem();
    clear_buf_mem();
    if(vec_rot != NULL){delete vec_rot; vec_rot = NULL;}

  }

  ~smear_fun(){
    clear_all_mem();
  }

};

/////smear plan buffer related
struct SmearPlanKey {
  Coordinate total_site;
  bool smear_in_time_dir;
  Int bfac;
  Int civ;
  Int dup;
  DATA_TYPE prec;
};

inline bool operator<(const SmearPlanKey& x, const SmearPlanKey& y)
{
  if(x.total_site < y.total_site ){  return true;}
  if(y.total_site < x.total_site ){  return false;}

  if(x.smear_in_time_dir < y.smear_in_time_dir ){  return true;}
  if(y.smear_in_time_dir < x.smear_in_time_dir ){  return false;}

  if(x.bfac < y.bfac ){return true;}
  if(y.bfac < x.bfac ){return false;}

  if(x.civ  < y.civ  ){return true;}
  if(y.civ  < x.civ  ){return false;}

  if(x.prec < y.prec ){return true;}
  if(y.prec < x.prec ){return false;}

  if(x.dup < y.dup ){return true;}
  if(y.dup < x.dup ){return false;}

  return false;
}

struct smear_fun_copy{
  DATA_TYPE prec;
  void* smfP;
  bool is_copy;  // do not free memory if is_copy=true

  smear_fun_copy(){smfP = NULL;is_copy = false;prec = ComplexD_TYPE;}
  smear_fun_copy(const smear_fun_copy& smf)
  {
    #ifndef QLAT_USE_ACC
    Qassert(false);
    #endif
    is_copy = true;
    smfP = smf.smfP;
    prec = smf.prec;
  }
  smear_fun_copy(smear_fun_copy&& smf) noexcept
  {
    is_copy = smf.is_copy;
    smfP = smf.smfP;
    prec = smf.prec;
    smf.is_copy = true;
  }

  inline void swap(smear_fun_copy& x)
  {
    Qassert(not is_copy);
    Qassert(not x.is_copy);
    void* tmp = smfP;
    smfP   = x.smfP;
    x.smfP = tmp;
    DATA_TYPE tmp_prec = prec;
    prec   = x.prec;
    x.prec = tmp_prec;
  }
  //
  const smear_fun_copy& operator=(const smear_fun_copy& smf)
  {
    Qassert(not is_copy);
    delete_pointer();
    prec = smf.prec;

    //smear_fun<T > smf;
    //smf.setup(gf, mom, smear_in_time_dir);
    if(smf.smfP == NULL){abort_r("smear fun point NULL!\n");}
    ////qmessage("Smear fun copy \n");
    //smf.is_copy = true;
    //smfP = smf.smfP;

    if(prec == ComplexD_TYPE ){
      const smear_fun<ComplexD  >* a =   ((const smear_fun<ComplexD  >*) smf.smfP);
      smfP = (void*) (new smear_fun<ComplexD  >(a->geo, a->smear_in_time_dir));
    }
    else if(prec == ComplexF_TYPE ){
      const smear_fun<ComplexF >* a =   ((const smear_fun<ComplexF >*) smf.smfP);
      smfP = (void*) (new smear_fun<ComplexF >(a->geo, a->smear_in_time_dir));
    }
    else{qmessage("Only ComplexD and ComplexF supported for smearing! \n");Qassert(false);}

    /////swap(*this, smf);
    return *this;
  }

  inline void delete_pointer()
  {
    if(smfP != NULL){
    if(prec == ComplexD_TYPE  ){delete ((smear_fun<ComplexD  >*)smfP);}
    if(prec == ComplexF_TYPE ){delete ((smear_fun<ComplexF >*)smfP);}
    smfP = NULL;
    }
  }

  void clear(){
    if(is_copy == false){delete_pointer();}
  }
  ~smear_fun_copy(){clear();}

};


inline smear_fun_copy make_smear_plan(const Geometry& geo, const bool smear_in_time_dir, const DATA_TYPE prec)
{
  TIMER_VERBOSE("make_smear_plan");
  smear_fun_copy st;st.prec = prec;
  if(prec == ComplexD_TYPE ){     st.smfP = (void*) (new smear_fun<ComplexD >(geo, smear_in_time_dir));}
  else if(prec == ComplexF_TYPE){st.smfP = (void*) (new smear_fun<ComplexF>(geo, smear_in_time_dir));}
  else{qmessage("Only ComplexD and ComplexF supported for smearing! \n");Qassert(false);}
  return st;
}

inline smear_fun_copy make_smear_plan(const SmearPlanKey& skey)
{
  Geometry geo;geo.init(skey.total_site);
  return make_smear_plan(geo, skey.smear_in_time_dir, skey.prec);
}

inline Cache<SmearPlanKey, smear_fun_copy >& get_smear_plan_cache()
{
  static Cache<SmearPlanKey, smear_fun_copy > cache("SmearPlanKeyCache", 64);
  return cache;
}

inline void clear_smear_plan_cache()
{
  get_smear_plan_cache().clear();
}

inline const smear_fun_copy& get_smear_plan(const SmearPlanKey& skey)
{
  if (!get_smear_plan_cache().has(skey)) {
    get_smear_plan_cache()[skey] = make_smear_plan(skey);
  }
  ////smear_fun_copy& = get_smear_plan_cache()[skey];
  ////qmessage("setup %5d %5d \n", skey.civ, skey.bfac);
  return get_smear_plan_cache()[skey];
}

template<typename Ty, Int bfac, Int civ>
inline SmearPlanKey get_smear_plan_key(const Geometry& geo, const bool smear_in_time_dir, Int dup = -1)
{
  SmearPlanKey skey;
  //skey.geo  = geo.total_site();
  skey.total_site  = geo.total_site();
  skey.bfac = bfac;
  skey.civ  = civ;
  skey.smear_in_time_dir = smear_in_time_dir;
  skey.prec = get_data_type<Ty >();
  skey.dup  = dup;
  return skey;
}
/////smear plan buffer related

////TODO need to change other parts for c0, c1
template <class T, class Td>
void extend_links_to_vecs(T* resE, const GaugeFieldT<Td >& gf, const qlat::vector<qlat::ComplexD >& mom_factors=qlat::vector<qlat::ComplexD >(), const qlat::vector<Long >& index_map = qlat::vector<Long >(0)){
  TIMERB("extend_links_to_vecs");
  const Geometry& geo = gf.geo();
  GaugeFieldT<Td > gf1;
  Coordinate expan_left(1, 1, 1, 1);
  Coordinate expan_zero(0, 0, 0, 0);
  const Geometry& geo_ex = geo_resize(gf.geo(), expan_left, expan_zero);
  gf1.init(geo_ex);
  copy_fields(gf1, gf);
  {
  TIMER("extend_links_to_vecs expand");
  refresh_expanded_GPU(gf1);
  }

  //set_left_expanded_gauge_field(gf1, gf);

  const Int dir_limit = 4;
  ////set up mom factors
  qlat::ComplexD* momF = NULL;
  qlat::vector<qlat::ComplexD > buf_mom;buf_mom.resize(8);
  for(Int i=0;i<buf_mom.size();i++){buf_mom[i] = 1;}
  if(mom_factors.size() == 0){momF = (qlat::ComplexD*) qlat::get_data(buf_mom).data();}
  if(mom_factors.size() != 0){
    Qassert(mom_factors.size()==8);
    momF = (qlat::ComplexD*) qlat::get_data(mom_factors).data();
  }

  const Long* index_mapT = (Long*) index_map.data();
  qlat::vector<Long > index_mapR;
  const Long Nvol = geo.local_volume();
  if(index_map.size() == 0){
    index_mapR.resize(Nvol);
    qacc_for(i, Nvol, {index_mapR[i] = i;});
    //for(Long i=0;i<Nvol;i++){index_mapR[i] = i;}
    index_mapT = index_mapR.data();
  }
  ////set up mom factors
  //mom_factors()[dir + 4];
  qacc_for(count,  geo.local_volume(),{
    Long index = index_mapT[count];
    const Coordinate xl = geo.coordinate_from_index(index);
    for (Int dir = -dir_limit; dir < dir_limit; ++dir) {
      const ColorMatrixT<Td > link =
          dir >= 0 ? gf1.get_elem(xl, dir)
                   : (ColorMatrixT<Td >)matrix_adjoint(
                         gf1.get_elem(coordinate_shifts(xl, dir), -dir - 1));
      ////convention used in gauss_smear_kernel CPU and gauss_smear_global4
      ////also in multiply_gauge of utils_shift_vecs.h
      for(Int ci=0; ci<9; ci++){
        //resE[(index*dir_limit*2+ (dir+dir_limit))*9 + (ci%3)*3 + ci/3] = momF[dir + 4] * link.p[ci];
        resE[(count*dir_limit*2+ (dir+dir_limit))*9 + ci] = momF[dir + 4] * link.p[ci];
      }
    }
  })
}

template <class T, class Td>
void extend_links_to_vecs(qlat::vector_gpu<T >& gfE, const GaugeFieldT<Td >& gf, const qlat::vector<qlat::ComplexD >& mom_factors=qlat::vector<qlat::ComplexD >()){
  gfE.resize(2*4* gf.geo().local_volume() * 9);
  extend_links_to_vecs(gfE.data(), gf, mom_factors);
}

/////T is double or float
template <class T>
void rotate_Vec_prop(Propagator4dT<T>& prop, qlat::vector_gpu<ComplexT<T> > &propT, unsigned int NVmpi, unsigned int groupP, Int dir = 0)
{
  TIMER("Rotate Vec prop");
  ///unsigned int NVmpi = fd.mz*fd.my*fd.mx;
  ///int groupP = (12+NVmpi-1)/NVmpi;
  Qassert(groupP <= 12);Qassert(groupP * NVmpi >= 12);
  Qassert(groupP >   0);

  const Geometry& geo = prop.geo();
  const Long Nvol =  geo.local_volume();
  if(dir == 0)propT.resize(NVmpi*Nvol*groupP*12);

  qacc_for(index, Long(Nvol), {
    qlat::WilsonMatrixT<T>& v0 =  prop.get_elem_offset(index);

    for(Int c1 = 0;c1 < 3; c1++)
    for(Int d1 = 0;d1 < 4; d1++)
    for(Int c0 = 0;c0 < 3; c0++)
    for(Int d0 = 0;d0 < 4; d0++)
    {
      Int off1 = c1*4+d1;
      Int n0 = off1/groupP;
      Int n1 = off1%groupP;

      Int off0 = c0*4+d0;
      //LInt off = (c0*3+c1)*16+d0*4 + d1;
      //Long offP = n0*Nvol*groupP*12 + index*groupP*12 + n1*12 + off0;
      Long offP = ((n0*Nvol + index)*groupP + n1)*12 + off0;
      if(dir == 0){propT[offP] = v0(d0*3+c0, d1*3+c1);}
      if(dir == 1){v0(d0*3+c0, d1*3+c1) = propT[offP];}
    }
  });


}

inline void get_smear_para(std::string smear_para, double& width, Int& step)
{
  if(smear_para == std::string("NONE")){width = 0.0;step = 0;return; }
  std::vector<std::string > temL = stringtolist(smear_para);
  if(temL.size() != 2){abort_r("read smear para wrong! \n");}
  width = stringtodouble(temL[0]);
  step  = stringtonum(temL[1]);
}

template <class T, Int bfac, Int civ>
void smear_propagator_box_kernel(qlat::vector_gpu<T >& prop, qlat::vector_gpu<T >& vp, qlat::vector_gpu<T >& vm,
  const Int Bsize, const Int dir, shift_vec& svec)
{
  //return ;
  TIMERC("Shift Kernel");
  std::vector<Int > iDir(4);
  for(Int i=0;i<4;i++){iDir[i] = 0;}

  vp.copy_from(prop);vm.copy_from(prop);
  for(Int bi=0;bi<Bsize;bi++){
    iDir[dir] =  1;svec.shift_vecP(vp.p, vp.p, iDir , bfac*3*civ);
    iDir[dir] = -1;svec.shift_vecP(vm.p, vm.p, iDir , bfac*3*civ);
    prop += vp;
    prop += vm;
  }
}

template <class T, Int bfac, Int civ>
void smear_propagator_box(T* src, const Int Bsize, smear_fun<T >& smf){
  if (Bsize <= 0) {return;}
  /////Qassert(!smear_in_time_dir);
  smf.check_setup();
  Long Nvol = smf.Nvol;
  //const Geometry& geo = smf.geo;
  //const Int dir_limit = smear_in_time_dir ? 4 : 3;
  const Int nsites = bfac*3*civ;
  const Int dir_limit = 3;
  /////std::vector<Int > nv,Nv,mv;geo_to_nv(geo, nv, Nv, mv);

  const Int GPU = 1;
  qlat::vector_gpu<T >& gauge_buf = get_vector_gpu_plan<T >(0, smf.gauge_buf_name, GPU);

  //qlat::vector_gpu<T > gfE; 
  //extend_links_to_vecs(gfE, gf);
  ////fft_desc_basic fd(geo);
  shift_vec svec(smf.fd_new, true);
  {
    TIMERC("Box smear setup Gauge");
    svec.set_gauge(gauge_buf.data(), bfac, civ);
  }
  //rotate_prop(prop,0);

  size_t Ndata = size_t(Nvol) * nsites;
  ////T* pdata = (T*) qlat::get_data(prop).data();

  for(unsigned int i=0;i<smf.vL.size();i++){smf.vL[i].resize(Ndata);smf.vL[i].set_zero();}
  qlat::vector_gpu<T >& vprop = smf.vL[0];
  qlat::vector_gpu<T >& v1 = smf.vL[1];
  qlat::vector_gpu<T >& v2 = smf.vL[2];
  qlat::vector_gpu<T >& vp = smf.vL[3];
  qlat::vector_gpu<T >& vm = smf.vL[4];
  std::vector<qlat::vector_gpu<T >*  > vL(3);
  for(unsigned int i=0;i<vL.size();i++){vL[i] = &smf.vL[i + 5];}

  //qlat::vector_gpu<T > vprop; vprop.resize(Ndata);
  //qlat::vector_gpu<T > v1; v1.resize(Ndata);
  //qlat::vector_gpu<T > v2; v2.resize(Ndata);
  //qlat::vector_gpu<T > vp; vp.resize(Ndata);
  //qlat::vector_gpu<T > vm; vm.resize(Ndata);
  //std::vector<qlat::vector_gpu<T >  > vL(3);
  //for(unsigned int i=0;i<vL.size();i++){vL[i].resize(Ndata);vL[i].set_zero();}

  cpy_data_thread(vprop.p, src, Ndata);

  for(Int dir = 0; dir < dir_limit; ++dir)
  {
    v1.copy_from(vprop);
    smear_propagator_box_kernel<T,bfac,civ>(v1, vp, vm, Bsize, dir, svec);
    Int d1=(dir+1)%3;int d2=(dir+2)%3;v2.copy_from(v1);
    smear_propagator_box_kernel<T,bfac,civ>(v1, vp, vm, Bsize, d1 , svec);
    smear_propagator_box_kernel<T,bfac,civ>(v2, vp, vm, Bsize, d2 , svec);
    (*vL[d1]) += v2;
    (*vL[d2]) += v1;
  }
  vprop.set_zero();
  for(Int dir = 0; dir < dir_limit; ++dir)
  {
    smear_propagator_box_kernel<T,bfac,civ>(*vL[dir], vp, vm, Bsize, dir, svec);
    vprop += (*vL[dir]);
  }

  T* data = vprop.p;
  qacc_for(index, Nvol, {
    for(Int c=0;c<nsites;c++){src[index*nsites + c] = data[index*nsites + c] * (1.0/6.0);} 
  });

  smf.clear_mem();
}

template <class T, Int bfac, Int civ>
void gauss_smear_kernel(T* src, const double width, const Int step, const T norm, smear_fun<T >& smf)
{
  double bw_fac = width*width/(4.0*step - 6.0*width*width);
  ////match convention to qlat
  if(smf.dirL == 4 ){bw_fac = bw_fac * 3.0/4.0;}
  const T bw = bw_fac;

  const Long Nvol = smf.Nvol;
  const Int nsites = bfac*3*civ;

  const Int GPU = 1;
  qlat::vector_gpu<T >& gauge_buf = get_vector_gpu_plan<T >(0, smf.gauge_buf_name, GPU);
  qlat::vector_gpu<T >& prop_buf0 = get_vector_gpu_plan<T >(0, smf.prop_buf0_name, GPU);
  qlat::vector_gpu<T >& prop_buf1 = get_vector_gpu_plan<T >(0, smf.prop_buf1_name, GPU);

  smf.check_setup();

  const T* gf = (T*) gauge_buf.data();
  const Long Nvol_sum = smf.map_index_typeAL.size();////only Nvol have data, just save Nvol_sum
  prop_buf0.resizeL(Nvol_sum * bfac * 3* civ);
  prop_buf1.resizeL(Nvol_sum * bfac * 3* civ);
  ////cpy_data_thread(prop, src, Ndata);
  ////copy src to buf0

  ////if distributed, then no need to copy and src already in prop_buf0, not expansion also

  //prop_buf1.copy_from(prop_buf0); /// copy src

  ////int dirL = 3;
  #ifdef QLAT_USE_ACC
  size_t Ndata = smf.Nvol * bfac * 3* civ;
  Int nt  = 3*3*9;
  /////if(3*bfac*civ <= 36){ nt =         36;}
  if(3*bfac*civ <= nt){ nt =         3*bfac*civ;}
  if(bfac*civ <=  6){ nt = 3*bfac*civ;}
  //std::string val = get_env(std::string("q_smear_gpu_thread"));
  //if(val != ""){nt = stringtonum(val);}

  //int nt = 32;
  dim3 dimBlock(nt, 1, 1);
  const Long sn = Nvol;
  dim3 dimGrid( sn, 1, 1);

  #endif
  //const Long* map_final = &smf.map_bufV[0][0];

  //const Long* Pdir1 = (Long*) qlat::get_data(smf.map_bufD).data();
  const Long* Pdir1 = (Long*) qlat::get_data(smf.map_bufD_typeA).data();
  const Long* map_final = &smf.map_index_typeA1[0];

  T* prop_src = (T*) prop_buf0.data();
  T* prop_tmp = (T*) prop_buf0.data();
  T* prop_res = (T*) prop_buf1.data();

  if(smf.redistribution_copy == 0)
  {
  cpy_data_from_index(prop_src, src, &smf.map_index_typeA0[0], &smf.map_index_typeI[0], Nvol, nsites, GPU);
  //cpy_data_from_index(prop_src, src, &smf.map_bufV[0][0], &smf.map_bufV[1][0], Nvol, nsites, GPU);
  }

  //const Long Ndata = smf.Nvol * nsites;
  ////cpy_data_from_index(prop_res, src, &smf.map_bufV[1][0], &smf.map_bufV[1][0], Nvol, nsites, GPU);
  //size_t Ndata = smf.Nvol * bfac * 3* civ;
  //cpy_data_thread(prop_res, src, Ndata);
  //cpy_data_thread(prop_res, prop_src, Ndata);

  for (Int i = 0; i < step; ++i) {
    //{TIMERC("Copy prop");
    //cpy_data_from_index(prop_src, prop_res, &smf.map_bufV[0][0], &smf.map_bufV[0][0], Nvol, nsites, GPU);}

    {TIMERC("Communication");
    smf.template refresh_expanded_GPU<nsites>(prop_src, GPU);}

    {TIMERC("Matrix multiply");
    #ifdef QLAT_USE_ACC
    if(smf.dirL==3){gauss_smear_global4<T, bfac, civ, 3><<< dimGrid, dimBlock >>>(prop_res, prop_src, gf, bw, norm, Nvol, Pdir1, map_final);}
    if(smf.dirL==4){gauss_smear_global4<T, bfac, civ, 4><<< dimGrid, dimBlock >>>(prop_res, prop_src, gf, bw, norm, Nvol, Pdir1, map_final);}
    qacc_barrier(dummy);
    #else

    const Int dir_max  = 4;
    const Int dir_limit = smf.dirL;
    ////const Long* map_count = smf.map_index_typeAL[0];
    qthread_for(count,  Nvol,{
      //const Long index  = smf.map_index_typeAL[count];
      //const Long iwrite = smf.map_index_typeA0[index];
      const Long iwrite = map_final[count];
      QLAT_ALIGN(QLAT_ALIGNED_BYTES) T  buf[nsites];

      for(Int i=0;i<nsites;i++){ 
        buf[i] = 0; 
      }
      for (Int dir = -dir_limit; dir < dir_limit; ++dir) {
        ////const Long index_dir = smf.local_map_typeA0[Pdir1[index*dir_limit*2 + (dir + dir_limit)]];
        const Long index_dir = Pdir1[count*dir_limit*2 + (dir + dir_limit)];
        const T* wm1p = &prop_src[size_t(index_dir)*nsites];
        const T* lp = &gf[(count*dir_max*2 + dir + dir_max)*9];
        //const T* l0 = &gf[(count*dir_max*2 + dir + dir_max)*9];
        //for(Int c0=0;c0<3;c0++)
        //for(Int c1=0;c1<3;c1++)
        //{
        //  lbuf[c0*3+c1] = l0[c0*3+c1];
        //}

        for(Int bi=0;bi<bfac;bi++){
          Eigen::Matrix<T, 3, 3, Eigen::ColMajor>&     lE = *((Eigen::Matrix<T, 3, 3, Eigen::ColMajor>*) lp);
          //Eigen::Matrix<T, 3, 3, Eigen::RowMajor>&     lE = *((Eigen::Matrix<T, 3, 3, Eigen::RowMajor>*) lp);
          Eigen::Matrix<T, civ, 3, Eigen::ColMajor>& wmE = *((Eigen::Matrix<T, civ, 3, Eigen::ColMajor>*)  &buf[bi*3*civ]);

          //Eigen::Matrix<T, civ, 3, Eigen::ColMajor>&  pE = *((Eigen::Matrix<T, civ, 3, Eigen::ColMajor>*) &wm1p[bi*3*civ]);
          //wmE += ( pE * lE);

          if(civ != 1){
            Eigen::Matrix<T, civ, 3, Eigen::ColMajor>&  pE = *((Eigen::Matrix<T, civ, 3, Eigen::ColMajor>*) &wm1p[bi*3*civ]);
            wmE += ( pE * lE);
          }
          if(civ == 1){
            Eigen::Matrix<T, civ, 3, Eigen::RowMajor>&  pE = *((Eigen::Matrix<T, civ, 3, Eigen::RowMajor>*) &wm1p[bi*3*civ]);
            wmE += ( pE * lE);
          }////avoid mix of Col and Row when civ == 1
        }
      }

      const T* wo = &prop_src[iwrite*nsites];
            T* wp = &prop_res[iwrite*nsites];
      //const T* wo = &prop_src[map_final[index]*nsites];
      //T* wp = &prop_res[map_final[index]*nsites];

      for(Int ic=0;ic<nsites;ic++){wp[ic]  = norm*(wo[ic] + bw*buf[ic]);}
    });

    #endif
    }

    //cpy_data_thread(prop_src, prop_res, Ndata);

    //////switch buffers

    prop_tmp = prop_src;
    prop_src = prop_res;
    prop_res = prop_tmp;

    /////cpy_data_thread(prop_res, src, Ndata);
  }

  if(smf.redistribution_copy == 0){
  cpy_data_from_index(src, prop_src, &smf.map_index_typeI[0], &smf.map_index_typeA0[0], Nvol, nsites, GPU);
  //cpy_data_from_index(src, prop_src, &smf.map_bufV[1][0], &smf.map_bufV[0][0], Nvol, nsites, GPU);
  }

  smf.clear_mem();
}

template <class T, Int bfac, Int civ>
void smear_kernel_sort(T* src, const double width_in, const Int step, smear_fun<T >& smf)
{
  ////width < 0, then additional normalization
  if(step >= 0){
    double sm_factor = 1.0;
    double width = width_in;
    if(step > 0 and width < 0){
      width = -1.0 * width;
      sm_factor = std::pow( std::pow( 1.0 * width, 3), 1.0/step );
    }
    const double norm   = (1 - 3.0*width*width/(2*step)) * sm_factor;
    gauss_smear_kernel<T, bfac, civ>(src, width, step, norm, smf);
  }

  if(step ==-1){
    if(smf.smear_in_time_dir == true){abort_r("Box smearing not supported for t smearings \n");}
    smear_propagator_box<T, bfac, civ>(src,int(width_in + 0.001), smf);
  }
}

template <class T>
void smear_kernel(T* src, const double width, const Int step, smear_fun<T >& smf, Int bfac, Int civ)
{
  bool cfind = false;
  {
  TIMERC("==smear_kernel");
  #define smear_macros(ba,da) if(bfac == ba and civ ==  da and cfind == false){cfind = true; \
    smear_kernel_sort<T, ba, da >(src, width, step, smf);}

  /////max size of a prop
  ///////macros for color 3 * dirac 4 in inner prop
  smear_macros(   1,  4);
  smear_macros(   2,  4);
  smear_macros(   3,  4);
  smear_macros(   4,  4);
  smear_macros(   5,  4);
  smear_macros(   6,  4);
  smear_macros(   7,  4);
  smear_macros(   8,  4);
  smear_macros(   9,  4);
  smear_macros(  10,  4);
  smear_macros(  11,  4);
  smear_macros(  12,  4);

  smear_macros(   9,  3);
  smear_macros(   9,  2);

  ///////macros for inner color 3 and all outter prop
  smear_macros(   1,  1);
  smear_macros(   2,  1);
  smear_macros(   3,  1);
  smear_macros(   4,  1);
  smear_macros(   5,  1);
  smear_macros(   6,  1);
  smear_macros(   7,  1);
  smear_macros(   8,  1);
  smear_macros(   9,  1);
  smear_macros(  10,  1);
  smear_macros(  11,  1);
  smear_macros(  12,  1);

  smear_macros(   6,  3);
  smear_macros(  24,  1);
  smear_macros(  48,  1);

  #ifndef QLAT_USE_ACC
  //smear_macros( 384 , 1);
  smear_macros(  12,  12);
  smear_macros(   1, 384);
  smear_macros(   1,  96);
  #endif
  smear_macros(   1, 48);
  smear_macros(   1, 27);
  smear_macros(   1, 24);
  smear_macros(   1, 18);
  smear_macros(   1, 16);
  smear_macros(   1, 12);
  smear_macros(   1,  9);
  smear_macros(   1,  8);
  smear_macros(   1,  6);
  smear_macros(   1,  3);
  smear_macros(   1,  2);

  smear_macros(   4, 12);
  smear_macros(   4,  3);
  smear_macros(   4,  6);
  smear_macros(   4,  8);
  smear_macros(   4,  2);
  smear_macros(   3, 16);

  qacc_barrier(dummy);
  #undef smear_macros
  }
  if(cfind == false){qmessage("Case bfac %5d, civ %5d \n", bfac, civ); abort_r("smear kernel not found!\n");}
  ////Qassert(cfind);
  //qmessage("==check Case bfac %5d, civ %5d \n", bfac, civ);
}

////4*3_0   4*3_1 --> 3_0 --> 4*4*3_1
template <class T>
void rotate_prop(Propagator4dT<T>& prop, Int dir = 0)
{
  TIMERB("rotate_prop");
  const Geometry& geo = prop.geo();
  const Long Nvol =  geo.local_volume();

  ComplexT<T>* src =  (ComplexT<T>*) qlat::get_data(prop).data();
  qacc_for(index, Long(Nvol), {
    ComplexT<T> buf[12*12];
    ComplexT<T>* res = &src[index*12*12];
    for(Int i=0;i<12*12;i++){buf[i] = res[i];}

    for(Int d1 = 0;d1 < 4; d1++)
    for(Int c1 = 0;c1 < 3; c1++)
    for(Int d0 = 0;d0 < 4; d0++)
    for(Int c0 = 0;c0 < 3; c0++)
    {
      Int off0 = c1*4*4*3 + ((d1*4+d0)*3+c0);
      Int off1 = (d1*3+c1)*12 + d0*3 + c0;
      if(dir == 0){res[off0] = buf[off1];}
      if(dir == 1){res[off1] = buf[off0];}
    }
  });
}

/*
  Td is double or float
  gassian normalization -1, meson w^{-3 * 2}, baryon w^{-3 * 3}
*/
template <class Ty, Int c0,Int d0, class Td>
void smear_propagator_gwu_convension_inner(
    Ty* prop, const GaugeFieldT<Td >& gf,
    const double width,
    const Int step,
    const CoordinateD& mom = CoordinateD(),
    const bool smear_in_time_dir = false,
    const Int mode = 1,
    const Int dup = -1,
    const Int force_update = 0
    )
{
  TIMER_FLOPS("smear_propagator_gwu_convension_inner");
  Long Tfloat = 0;
  ///double mem       = 0.0;
  Qassert(gf.initialized);
  const Geometry& geo = gf.geo();
  Qassert(geo.eo == 0);
  #ifdef QLAT_USE_ACC
  if(c0*d0 > 48){abort_r("d0 should be smaller than 48 for gpu mem. \n");}
  #endif

  Int use_gauge_mapping = 1;
  if(step < 0){ use_gauge_mapping = 0;}////hack for box smearings

  {
    const Long Lat = geo.local_volume();
    Int nsrc = c0*d0;
    const Long vGb = Lat *nsrc;
    const Int n_avg = smear_in_time_dir ? 8 : 6;
    Int Fcount = 3*(3*6 + 2*2); 
    if(step >= 0){
    Tfloat = step*n_avg*vGb*Fcount;
    }else{
    //Tfloat = c0*d0*2*int(qlat::qnorm(width))*vGb*Fcount;
    /////qmessage("%3d dir %d vol %d Flops %d \n", int(width), 6, int(vGb), int(Fcount));
    Tfloat = 4 * 6*int(width) * vGb * Fcount;
    }
  }
  timer.flops += Tfloat;

  if (0 == step) {return;}
  /////smear_fun<T > smf(gf.geo(), smear_in_time_dir);

  Ty* src = prop;
  //Ty* src = (Ty*) qlat::get_data(prop).data();
  /////rotate_prop(prop, 0);

  SmearPlanKey skey = get_smear_plan_key<Ty, c0, d0>(geo, smear_in_time_dir, dup);
  const smear_fun_copy& smf_copy = get_smear_plan(skey);
  smear_fun<Ty >& smf = *((smear_fun<Ty >*) (smf_copy.smfP));
  ////hack to reuse links
  smf.gauge_buf_name_append(skey.dup);

  ////fft_desc_basic fd(geo);
  ////Vec_redistribute vec_large(fd);

  Long Nvol_pre = geo.local_volume();
  bool reorder = false;
  Int groupP = (d0+smf.NVmpi-1)/smf.NVmpi;

  ////check on parameters for smearings
  bool ini_check = false;
  ////why do I need d0%smf.NVmpi == 0 previously
  //if(smf.NVmpi <= d0 or smf.NVmpi != 1 and (d0%smf.NVmpi == 0))

  ////do this check only on GPU, CPU set it by hand
  #ifdef QLAT_USE_ACC
  if(smf.NVmpi <= d0 or smf.NVmpi != 1)
  {
    //ini_check = true;
    ini_check = false;////disable rotation by default
  }
  #endif
  std::string val = get_env(std::string("q_smear_do_rotate"));
  if(val != ""){ini_check = stringtonum(val);}
  if(ini_check){
    if(smf.NVmpi == 1){ini_check = false;};
    Qassert(groupP * smf.NVmpi >= d0);
  }

  if(ini_check and smear_in_time_dir == false and mode == 1){reorder = true;}
  Int Nprop = smf.NVmpi * c0*3 * groupP;

  const Int GPU = 1;
  qlat::vector_gpu<Ty >& prop_buf0 = get_vector_gpu_plan<Ty >(0, smf.prop_buf0_name, GPU);
  qlat::vector_gpu<Ty >& prop_buf1 = get_vector_gpu_plan<Ty >(0, smf.prop_buf1_name, GPU);

  if(reorder){use_gauge_mapping = 0;}////hack for reordering
  smf.gauge_setup(gf, mom, force_update, use_gauge_mapping );

  if(reorder ){
    /////qmessage("====Vec setup, NVmpi %d, groupP %d \n", smf.NVmpi, groupP);
    smf.init_distribute();
    //smf.prop.resize(    Nvol_pre * Nprop );
    //smf.prop_buf.resize(Nvol_pre * Nprop );
    size_t psize = Nvol_pre * Nprop;
    prop_buf0.resizeL(psize);
    prop_buf1.resizeL(psize);
    Qassert(!smear_in_time_dir);

    Ty* res = prop_buf0.data();prop_buf0.set_zero();

    ////need to check this rotation
    cpy_GPU2D(res, src, Long(d0), Long(Nvol_pre*c0*3), Long(smf.NVmpi*groupP), Long(d0), QMGPU, QMGPU);
    //void copy_buffers_vecs(Ty *res, Ty *src,Long N0, Long N1, Long Ncopy, size_t Vol, Int GPU = 1)
    //copy_buffers_vecs(res, src, smf.NVmpi*groupP, d0, d0, Nvol_pre*c0*3, 1);
    smf.mv_idx.dojob(res, res, 1, smf.NVmpi, Nvol_pre*c0*3, 1,   groupP, true);
    {TIMERC("Vec prop");smf.vec_rot->reorder(res, prop_buf1.data(), 1, c0*3*groupP ,   0);}
    src = prop_buf0.data();
  }

  ////qmessage("===Case %d %d \n", c0, grouP);
  if( reorder)smear_kernel(src, width, step, smf,  c0, groupP);
  if(!reorder)smear_kernel(src, width, step, smf,  c0, d0);

  if(reorder ){
    Ty* res = prop_buf0.data();
    {TIMERC("Vec prop");smf.vec_rot->reorder(res, prop_buf1.data(), 1, c0*3*groupP , 100);}
    smf.mv_idx.dojob(res, res, 1, smf.NVmpi, Nvol_pre*c0*3, 0,  groupP , true);
    cpy_GPU2D(prop, prop_buf0.data(), Long(d0), Long(Nvol_pre*c0*3), Long(d0), Long(smf.NVmpi*groupP), QMGPU, QMGPU);
    ////cpy_data_thread(prop, smf.prop.data(), smf.prop.size(), 1);
    //copy_buffers_vecs(prop, prop_buf0.data(), d0, smf.NVmpi*groupP, d0, Nvol_pre*c0*3, 1);
  }
  /////rotate_prop(prop, 1);

  //#if PRINT_TIMER>3
  //qmessage("====Vec setup, c0 %d, d0 %d, NVmpi %d, groupP %d , reorder %d \n", c0 , d0, smf.NVmpi, groupP, int(reorder));
  //#endif

}

//Calculate the radiuse of a source
template <typename Td>
double source_radius(Propagator4dT<Td >& prop, const Int tsrc = 0)
{
  double radius = 0; 
  double rho = 0; 
  double rho_x2 = 0; 
  const Geometry& geo = prop.geo();
  const Long Nvol = geo.local_volume();
  std::vector<Int > nv, Nv, mv;
  geo_to_nv(geo, nv, Nv, mv);

  //position p;
  //int tem_x =  (*source.vec[0]).desc->nx;
  //qmessage("tem_x---------%3d\n",tem_x);
  for(Long isp = 0; isp < Nvol; isp ++)
  {
    double rx2 = 0; 
    Coordinate xl   = geo.coordinate_from_index(isp);
    Coordinate p    = geo.coordinate_g_from_l(xl);
    if(p[3] == tsrc)
    {    
      const qlat::WilsonMatrixT<Td >& p1 =  prop.get_elem_offset(isp);

      double x = 0.0; 
      double y = 0.0; 
      double z = 0.0; 
      if(p[0] >= int(nv[0]/2)){x = -nv[0] + p[0];}else{x = p[0];}
      if(p[1] >= int(nv[1]/2)){y = -nv[1] + p[1];}else{y = p[1];}
      if(p[2] >= int(nv[2]/2)){z = -nv[2] + p[2];}else{z = p[2];}

      rx2 = x*x + y*y + z*z; 
      qlat::ComplexT<Td > tem = 0.0;
      for(Int dc0 =0;dc0<12;dc0++)
      for(Int dc1 =0;dc1<12;dc1++)
      {    
        tem += qlat::qnorm( p1(dc0, dc1) );
      }    
      double tem_rho = tem.real();
      rho_x2 += rx2*(tem_rho);
      rho += (tem_rho);
    }
  }

  sum_all_size(&rho, 1);
  sum_all_size(&rho_x2, 1);
  //global_sum_all(&rho,1);
  //global_sum_all(&rho_x2,1);
  ////radius = sqrt(rho_x2/rho)*sqrt(2.0/3.0);
  radius = sqrt(rho_x2/rho);
  //Which is 1/\sigma in exp(-x^2/\sigma^2)
  return radius;
}

template <class Ty, class Td>
void smear_propagator_gwu_convension(FieldG<Ty >& prop, const GaugeFieldT<Td >& gf,
                      const double width, const Int step, const CoordinateD& mom = CoordinateD(), const bool smear_in_time_dir = false, const Int mode = 1, const Int dup = -1, const Int force_update = 0)
{
  if (0 == step) {return;}
  Qassert(prop.initialized and prop.multiplicity == 12 * 12 and prop.mem_order == QLAT_OUTTER);
  const Long Nvol = prop.geo().local_volume();
  Ty* src = (Ty*) qlat::get_data(prop).data();
  move_index mv_civ;int flag = 0;
  flag = 0;mv_civ.dojob(src, src, 1, 12*12, Nvol, flag, 1, true);

  qacc_for(isp, Nvol, {
    Ty buf[12*12];
    for(Int i=0;i<12*12;i++){buf[i] = src[isp*12*12 + i];}
    for(Int d0=0;d0<12*4;d0++)
    for(Int c0=0;c0<   3;c0++)
    {
      src[isp*12*12 + c0*12*4 + d0] = buf[d0*3 + c0];
    }
  });

  smear_propagator_gwu_convension_inner<Ty, 1, 12*4, Td>(src, gf, width, step, mom, smear_in_time_dir, mode, dup, force_update);

  qacc_for(isp, Nvol, {
    Ty buf[12*12];
    for(Int i=0;i<12*12;i++){buf[i] = src[isp*12*12 + i];}
    for(Int d0=0;d0<12*4;d0++)
    for(Int c0=0;c0<   3;c0++)
    {
      src[isp*12*12 + d0*3 + c0] = buf[c0*12*4 + d0];
    }
  });
  flag = 1;mv_civ.dojob(src, src, 1, 12*12, Nvol, flag, 1, true);
}

template <class Ty, class Td>
void smear_propagator_gwu_convension(qpropT& prop, const GaugeFieldT<Td >& gf,
                      const double width, const Int step, const CoordinateD& mom = CoordinateD(), const bool smear_in_time_dir = false, const Int mode = 1, const Int dup = -1, const Int force_update = 0)
{
  if (0 == step) {return;}
  Qassert(prop.initialized);
  Ty* src = (Ty*) qlat::get_data(prop).data();
  FieldG<Ty> tmp_prop;
  const Long Nd = 12 * 12 * prop.geo().local_volume();
  tmp_prop.set_pointer(src, Nd, prop.geo(), QMGPU, QLAT_OUTTER);
  smear_propagator_gwu_convension(tmp_prop, gf, width, step, mom, smear_in_time_dir, mode, dup, force_update);
}

template <class Ty, class Td>
void smear_propagator_gwu_convension(qlat::FieldM<ComplexT<Ty> , 12>& prop, const GaugeFieldT<Td >& gf,
                      const double width, const Int step, const CoordinateD& mom = CoordinateD(), const bool smear_in_time_dir = false, const Int mode = 1, const Int dup = -1, const Int force_update = 0)
{
  if (0 == step) {return;}
  ComplexT<Ty>* src = (ComplexT<Ty>*) qlat::get_data(prop).data();
  qacc_for(isp, prop.geo().local_volume(), {
    ComplexT<Ty> buf[12];
    for(Int i=0;i<12;i++){buf[i] = src[isp*12 + i];}
    for(Int d0=0;d0<4;d0++)
    for(Int c0=0;c0<   3;c0++)
    {
      src[isp*12 + c0*4 + d0] = buf[d0*3 + c0];
    }
  });

  smear_propagator_gwu_convension_inner<ComplexT<Ty>, 1, 4, Td>(src, gf, width, step, mom, smear_in_time_dir, mode, dup, force_update);
  qacc_for(isp, prop.geo().local_volume(), {
    ComplexT<Ty> buf[12];
    for(Int i=0;i<12;i++){buf[i] = src[isp*12 + i];}
    for(Int d0=0;d0<4;d0++)
    for(Int c0=0;c0<   3;c0++)
    {
      src[isp*12 + d0*3 + c0] = buf[c0*4 + d0];
    }
  });
}


template <class Ty, class Td>
void smear_propagator_gwu_convension(Propagator4dT<Ty>& prop, const GaugeFieldT<Td >& gf,
                      const double width, const Int step, const CoordinateD& mom = CoordinateD(), const bool smear_in_time_dir = false, const Int mode = 1, const Int dup = -1, const Int force_update = 0)
{
  if (0 == step) {return;}
  rotate_prop(prop, 0);
  ComplexT<Ty>* src = (ComplexT<Ty>*) qlat::get_data(prop).data();
  smear_propagator_gwu_convension_inner<ComplexT<Ty>, 1, 12*4, Td>(src, gf, width, step, mom, smear_in_time_dir, mode, dup, force_update);
  rotate_prop(prop, 1);
}

template <class Ty, class Td>
void smear_propagator_wuppertal_convension(Propagator4dT<Ty>& prop, const GaugeFieldT<Td >& gf,
                      const double width, const Int step, const CoordinateD& mom = CoordinateD(), const bool smear_in_time_dir = false, const Int mode = 1, const Int dup = -1, const Int force_update = 0)
{
  if (0 == step) {return;}
  rotate_prop(prop, 0);
  ComplexT<Ty>* src = (ComplexT<Ty>*) qlat::get_data(prop).data();
  double factor_sigma = std::sqrt( 2 * step * width / 3.0);
  smear_propagator_gwu_convension_inner<ComplexT<Ty>, 1, 12*4, Td>(src, gf, factor_sigma, step, mom, smear_in_time_dir, mode, dup, force_update);
  rotate_prop(prop, 1);
}

template <class T, class Td>
void prop_smear_qlat_convension(Propagator4dT<T>& prop,
                                const GaugeFieldT<Td>& gf, const double coef,
                                const Int step,
                                const CoordinateD& mom = CoordinateD(),
                                const bool smear_in_time_dir = false,
                                const Int mode = 1, const Int dup = -1,
                                const Int force_update = 0)
{
  if (coef <= 0) {
    return;
  }
  double width = 0.0;
  if (step < 0) {
    width = coef;
  }  ////box smearings
  if (step >= 0) {
    width = std::sqrt(coef * 2 * step / (3.0));
  }  ////gauss smearings
  smear_propagator_gwu_convension(prop, gf, width, step, mom, smear_in_time_dir,
                                  mode, dup, force_update);
}

}
#endif


