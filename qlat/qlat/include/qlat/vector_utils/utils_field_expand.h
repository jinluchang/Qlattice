// utils_field_expand.h
// Gen Wang
// Apr. 2024

#ifndef UTILS_FIELD_EXPAND_H
#define UTILS_FIELD_EXPAND_H
#pragma once

#include <string.h>
#include <sys/resource.h>
#include <mpi.h>
#include <time.h>
#include <typeinfo>

#include "utils_float_type.h"
#include <type_traits>

#include <iterator>
#include "utils_read_txt.h"
#include "utils_vector_GPU.h"
#include "utils_mpi.h"

namespace qlat
{

void setup_expand(const Geometry& geo, const Int multiplicity, qlat::vector_acc<Long>& pack_send,
                  qlat::vector_acc<Long>& pack_recv,
                  const SetMarksField& set_marks_field = set_marks_field_all,
                  const std::string& tag = std::string(""));

void set_marks_field_dir(CommMarks& marks, const Geometry& geo,
                         const Int multiplicity, const std::string& tag);

struct expand_index_buf {
  //Geometry geo; //make a copy of geo if needed
  qlat::vector_acc<Long >  pack_send;
  qlat::vector_acc<Long >  pack_recv;
  Int Multiplicity;//record numbers
  //const long Nindex;
  expand_index_buf()
  {
    pack_send.resize(0);
    pack_recv.resize(0);
  }

  expand_index_buf(const Geometry& geo_, const Int multiplicity, const std::string& tag)
  {
    //geo = geo_;
    //const Int multiplicity = 1;//always 1 for the buffers reuse
    Multiplicity = multiplicity;
    if(tag == std::string("")){
      setup_expand(geo_, multiplicity, pack_send, pack_recv, set_marks_field_all, tag);
    }else{
      setup_expand(geo_, multiplicity, pack_send, pack_recv, set_marks_field_dir, tag);
    }
  }

  ~expand_index_buf()
  {
    pack_send.resize(0);
    pack_recv.resize(0);
  }
};

/////buffers for expand index
struct expand_index_Key {
  Geometry geo;
  Int multiplicity;
  std::string tag;
  //Coordinate total_site;
  //Coordinate expansion_left ;
  //Coordinate expansion_right;
  expand_index_Key(const Geometry& geo_, const Int multiplicity_, const std::string& tag_)
  {
    geo = geo_;
    multiplicity = multiplicity_;
    tag = tag_;
  }
};

bool compare_geo(const Geometry& g0, const Geometry& g1);

bool Compare_geo(const Geometry& g0, const Geometry& g1);

bool compare_less(const Geometry& g0, const Geometry& g1);

bool operator<(const expand_index_Key& x, const expand_index_Key& y);

inline Cache<expand_index_Key, expand_index_buf >& get_expand_index_buf_cache()
{
  static Cache<expand_index_Key, expand_index_buf > cache("expand_index_Key", 64);
  return cache;
}

inline expand_index_buf& get_expand_index_buf_plan(const expand_index_Key& ekey)
{
  if (!get_expand_index_buf_cache().has(ekey)) {
    get_expand_index_buf_cache()[ekey] = expand_index_buf(ekey.geo, ekey.multiplicity, ekey.tag);
  }
  expand_index_buf& buf = get_expand_index_buf_cache()[ekey];
  return buf;
}

inline expand_index_buf& get_expand_index_buf_plan(const Geometry& geo, const Int multiplicity, const std::string& tag)
{
  expand_index_Key ekey(geo, multiplicity, tag);
  return get_expand_index_buf_plan(ekey);
}


template <class M>
void refresh_expanded_GPUT(M* res, const Geometry& geo, const int MULTI, 
  const SetMarksField& set_marks_field = set_marks_field_all, const std::string& tag = std::string(""), int GPU = 1)
{
  Qassert(sizeof(M) % sizeof(double) == 0);
  //const int MULTI = geo.multiplicity;
  const Long mpi_size = MULTI * sizeof(M)/sizeof(double);
  const Int multiplicity = 1;//always 1 for the buffers reuse

  // Geometry geo1 = geo;
  //Qassert(geo1 == geo);
  // geo1.multiplicity = 1;// no multiplicity to save buffers

  const CommPlan& plan = get_comm_plan(set_marks_field, tag, geo, multiplicity);
  const Long total_bytes =
      (plan.total_recv_size + plan.total_send_size) * sizeof(M);
  if (0 == total_bytes) {
    return;
  }
  TIMER_FLOPS("refresh_expanded_GPU");
  timer.flops += total_bytes * MULTI / 2;

  std::vector<MPI_Request> reqs_send;
  std::vector<MPI_Request> reqs_recv;

  qlat::vector_gpu<char >& sbuf = qlat::get_vector_gpu_plan<char >(0, std::string("general_buf0"), GPU);
  qlat::vector_gpu<char >& rbuf = qlat::get_vector_gpu_plan<char >(0, std::string("general_buf1"), GPU);
  expand_index_buf& ebuf = get_expand_index_buf_plan(geo, multiplicity, tag);
  Qassert(ebuf.pack_send.size() == 2*plan.total_send_size and ebuf.pack_recv.size() == 2*plan.total_recv_size);

  const Long Nsend = plan.total_send_size ;
  const Long Nrecv = plan.total_recv_size ;

  sbuf.resizeL(Nsend * MULTI * sizeof(M) / sizeof(char));
  rbuf.resizeL(Nrecv * MULTI * sizeof(M) / sizeof(char));

  M* sP = (M*) &sbuf[0];
  M* rP = (M*) &rbuf[0];

  ////setup reciev
  const int mpi_tag = 10;
  for (size_t i = 0; i < plan.recv_msg_infos.size(); ++i) {
    const CommMsgInfo& cmi = plan.recv_msg_infos[i]; 
    mpi_irecv(&rP[cmi.buffer_idx * MULTI], cmi.size * mpi_size, MPI_DOUBLE,
              cmi.id_node, mpi_tag, get_comm(), reqs_recv);
  }

  //qlat::vector_acc<long > pack_infos;
  Long* pack_send = (Long*) &ebuf.pack_send[0];
  Long* pack_recv = (Long*) &ebuf.pack_recv[0];

  qGPU_for(isp, Nsend, GPU, {
    Long ri = pack_send[isp * 2 + 0] * MULTI;//offset with multivec
    Long si = pack_send[isp * 2 + 1] * MULTI;
    for(int n=0;n<MULTI;n++){sP[ri + n] = res[si + n];}
  });

  { 
    //TIMER("refresh_expanded-comm-init");
    for (size_t i = 0; i < plan.send_msg_infos.size(); ++i) {
      const CommMsgInfo& cmi = plan.send_msg_infos[i];
      mpi_isend(&sP[cmi.buffer_idx * MULTI], cmi.size * mpi_size, MPI_DOUBLE,
                cmi.id_node, mpi_tag, get_comm(), reqs_send);
    }
  }

  mpi_waitall(reqs_recv);////receive done and write
  qGPU_for(isp, Nrecv, GPU, {
    const Long ri = pack_recv[isp * 2 + 0] * MULTI;
    const Long si = pack_recv[isp * 2 + 1] * MULTI;
    for(int n=0;n<MULTI;n++){res[ri + n] = rP[si + n];}
  });

  mpi_waitall(reqs_send);

  //safe_free_vector_gpu_plan<char >(std::string("general_buf0"), GPU);
  //safe_free_vector_gpu_plan<char >(std::string("general_buf1"), GPU);
}

template <class M>
void refresh_expanded_GPU(M* res, const Geometry& geo, const int MULTI, const std::string& tag, int GPU = 1)
{
  if(tag == ""){
    refresh_expanded_GPUT(res, geo, MULTI, set_marks_field_all, std::string(""), GPU);
  }else{
    refresh_expanded_GPUT(res, geo, MULTI, set_marks_field_dir, tag, GPU);
  }
}

template <class M>
void refresh_expanded_GPU(M* res, const Geometry& geo, const int MULTI, int dir = -1000, int GPU = 1)
{
  Qassert(dir == -1000 or (dir >= -3-1 and dir <= 3));
  if(dir == -1000){
    refresh_expanded_GPUT(res, geo, MULTI, set_marks_field_all, std::string(""), GPU);
  }
  if(dir >= -3-1 and dir <= 3)
  {
    std::string tagL[8] = {"dirmt", "dirmz", "dirmy", "dirmx", "dirx", "diry", "dirz", "dirt"};
    refresh_expanded_GPUT(res, geo, MULTI, set_marks_field_dir, tagL[dir + 3 + 1], GPU);
  }

  //left or right expand
  if(dir == 500 or dir == 501)
  {
    std::string tagL[2] = {"dirL", "dirR"};
    refresh_expanded_GPUT(res, geo, MULTI, set_marks_field_dir, tagL[dir % 500], GPU);
  }
}

template <class M>
void refresh_expanded_GPU(Field<M>& f, int dir = -1000, int GPU = 1)
{
  M* res = (M*) qlat::get_data(f).data();
  refresh_expanded_GPU(res, f.geo(), f.multiplicity, dir, GPU);
}


}

#endif
