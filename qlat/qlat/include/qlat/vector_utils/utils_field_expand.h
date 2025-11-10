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

#define MAX_EXPAND_BUF 64

namespace qlat
{

void setup_expand(const Geometry& geo, const Int multiplicity, qlat::vector<Long>& pack_send,
                  qlat::vector<Long>& pack_recv,
                  const SetMarksField& set_marks_field = set_marks_field_all,
                  const std::string& tag = std::string(""));
void set_marks_field_dir(CommMarks& marks, const Geometry& geo,
                         const Int multiplicity, const std::string& tag);

struct expand_index_buf {
  qlat::vector<Long >  pack_send;
  qlat::vector<Long >  pack_recv;
  Int Multiplicity;//record numbers
  //const long Nindex;
  expand_index_buf()
  {
    pack_send.resize(0);
    pack_recv.resize(0);
  }
  //
  expand_index_buf(const Geometry& geo_, const Int multiplicity, const std::string& tag)
  {
    //const Int multiplicity = 1;//always 1 for the buffers reuse
    Multiplicity = multiplicity;
    if(tag == std::string("")){
      setup_expand(geo_, multiplicity, pack_send, pack_recv, set_marks_field_all, tag);
    }else{
      setup_expand(geo_, multiplicity, pack_send, pack_recv, set_marks_field_dir, tag);
    }
  }
  //
  ~expand_index_buf()
  {
    pack_send.resize(0);
    pack_recv.resize(0);
  }
};

// buffers for expand index
struct expand_index_Key {
  Geometry geo;
  Int multiplicity;
  std::string tag;
  expand_index_Key(const Geometry& geo_, const Int multiplicity_, const std::string& tag_)
  {
    geo = geo_;
    multiplicity = multiplicity_;
    tag = tag_;
  }
};

bool operator<(const expand_index_Key& x, const expand_index_Key& y);

inline Cache<expand_index_Key, expand_index_buf >& get_expand_index_buf_cache()
{
  static Cache<expand_index_Key, expand_index_buf > cache("expand_index_Key", MAX_EXPAND_BUF + 1);
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

inline void clear_expand_plan_cache()
{
  get_expand_index_buf_cache().clear();
}

// MPIs for expansion
inline std::vector<MPI_Request>& get_expand_buffer_reqs_send()
{
  static std::vector<MPI_Request> reqs_send;
  return reqs_send;
}

inline std::vector<MPI_Request>& get_expand_buffer_reqs_recv()
{
  static std::vector<MPI_Request> reqs_recv;
  return reqs_recv;
}

/*
  buffer for  field expansions
  index  : store the inx for send, recv copies, will be reused for other buffers
  bs, br : buffers for MPIs
  pres   : result memory locations
  dsize  : the length of each elements copied
*/
struct expand_field_buffer {
  Geometry geo;
  std::string tag;
  //expand_index_buf* index;
  vector<int8_t> bs;
  vector<int8_t> br;
  int8_t* pres;
  size_t dsize;// in terms of size of char
  bool initialized;
  //
  // clear the plan
  inline void init(){
    //index = NULL;
    pres  = NULL;
    dsize = 0;
    tag = std::string("");
    initialized = false;
  }
  //
  expand_field_buffer()
  {
    init();
    bs.resize(0);
    br.resize(0);
  }
  //
  inline Long get_flops(){
    if(!initialized){return 0;}
    expand_index_buf& ebuf = get_expand_index_buf_plan(geo, 1, tag);
    //const expand_index_buf& ebuf = *index;
    const Long Nsend = ebuf.pack_send.size() / 2;
    const Long Nrecv = ebuf.pack_recv.size() / 2;
    return (Nsend + Nrecv) * dsize / 2;
  }
  //
  template<typename T>
  void init(const Geometry& geo_, const std::string& tag_, T* res, const Int MULTI, const Int GPU){
    dsize = MULTI * sizeof(T);
    Qassert(dsize % sizeof(int8_t) == 0);
    //
    geo = geo_;
    tag = tag_;
    const expand_index_buf& ebuf = get_expand_index_buf_plan(geo, 1, tag);
    const Long Nsend = ebuf.pack_send.size() / 2;
    const Long Nrecv = ebuf.pack_recv.size() / 2;
    const Long Ns = Nsend * dsize / sizeof(int8_t);
    const Long Nr = Nrecv * dsize / sizeof(int8_t);
    //
    Qassert(GPU == 1 or GPU == 0);
    //
    if(GPU == 1){
      if(bs.mem_type != MemType::Acc or bs.size() < Ns or br.size() < Nr){
        bs.set_mem_type(MemType::Acc);
        br.set_mem_type(MemType::Acc);
        bs.resize(Ns);br.resize(Nr);
      }
    }
    //
    if(GPU == 0){
      if(bs.mem_type != MemType::Cpu or bs.size() < Ns or br.size() < Nr){
        bs.set_mem_type(MemType::Cpu);
        br.set_mem_type(MemType::Cpu);
        bs.resize(Ns);br.resize(Nr);
      }
    }
    //
    Qassert(bs.mem_type == br.mem_type);
    //
    pres = (int8_t*) res;
    //index = &ebuf;
    initialized = true;
  }
  //
  inline Long buf_size(){
    return (bs.size() + br.size()) * sizeof(char);
  }
  //
  inline void buf_clean(){
    bs.resize(0);
    br.resize(0);
  }
  //
  // type 0 : first copy, 1 : second copy
  template<typename T>
  void excute_copyT(Int type, Int MULTI )
  {
    TIMER_FLOPS("expand_field_buffer excute_copy");
    Qassert(initialized);
    Qassert(bs.size() > 0 and br.size() > 0  and bs.mem_type == br.mem_type);
    Int GPU = 1;
    if(bs.mem_type == MemType::Cpu or bs.mem_type == MemType::Comm){GPU = 0;}
    //
    T* sP  = (T*) bs.data();
    const T* rP  = (T*) br.data();
    T* res = (T*) pres;
    //
    const expand_index_buf& ebuf = get_expand_index_buf_plan(geo, 1, tag);
    //const expand_index_buf& ebuf = *index;
    const Long* pack_send = (Long*) &ebuf.pack_send[0];
    const Long* pack_recv = (Long*) &ebuf.pack_recv[0];
    //
    const Long Nsend = ebuf.pack_send.size() / 2;
    const Long Nrecv = ebuf.pack_recv.size() / 2;
    //
    if(type == 0){
      qGPU_for(isp, Nsend, GPU, {
        const Long ri = pack_send[isp * 2 + 0] * MULTI;//offset with MULTI
        const Long si = pack_send[isp * 2 + 1] * MULTI;
        for(Int n=0;n<MULTI;n++){sP[ri + n] = res[si + n];}
      });
      timer.flops += Nsend * MULTI * sizeof(T);
    }
    //
    if(type == 1){
      qGPU_for(isp, Nrecv, GPU, {
        const Long ri = pack_recv[isp * 2 + 0] * MULTI;
        const Long si = pack_recv[isp * 2 + 1] * MULTI;
        for(Int n=0;n<MULTI;n++){res[ri + n] = rP[si + n];}
      });
      init();
      timer.flops += Nrecv * MULTI * sizeof(T);
    }
  }
  //
  inline void excute_copy(Int type)
  {
    if(!initialized){return ;}
    bool find = false;
    if(!find and dsize % sizeof(int64_t) == 0){excute_copyT<int64_t>(type,  dsize / sizeof(int64_t));find = true;}
    if(!find and dsize % sizeof(int32_t) == 0){excute_copyT<int32_t>(type,  dsize / sizeof(int32_t));find = true;}
    if(!find and dsize % sizeof(int16_t) == 0){excute_copyT<int16_t>(type,  dsize / sizeof(int16_t));find = true;}
    if(!find and dsize % sizeof(int8_t ) == 0){excute_copyT<int8_t >(type,  dsize / sizeof(int8_t ));find = true;}
    if(!find and dsize % sizeof(float  ) == 0){excute_copyT<float  >(type,  dsize / sizeof(float  ));find = true;}
    if(!find and dsize % sizeof(char   ) == 0){excute_copyT<char   >(type,  dsize / sizeof(char   ));find = true;}
    Qassert(find == true);
  }
};

inline std::vector<expand_field_buffer>& get_expand_buffer_full_list()
{
  static std::vector<expand_field_buffer > bufs;
  if(bufs.size() != MAX_EXPAND_BUF){bufs.resize(MAX_EXPAND_BUF);};
  return bufs;
}

inline void clean_expand_buffer_vecs(){
  //
  static const Long max_bytes  = Long( get_env_double_default("q_expand_buffer_size", 0.5) * 1024 * 1024 * 1024 );
  std::vector<expand_field_buffer>& bufs = get_expand_buffer_full_list();
  Long total_size = 0;
  bool clean = false;
  //
  for(unsigned int i=0;i<bufs.size();i++){
    if(!bufs[i].initialized){continue;}
    total_size += bufs[i].buf_size();
    if(total_size > max_bytes){clean = true;}
    if(clean){
      TIMER("clean_expand_buffer_vecs");
      bufs[i].buf_clean();
    }
  }
}

inline void expand_buffer_wait_mpi(){
  TIMER_FLOPS("expand_buffer_wait_mpi");
  std::vector<expand_field_buffer>& bufs = get_expand_buffer_full_list();
  std::vector<MPI_Request>& reqs_recv = get_expand_buffer_reqs_recv();
  std::vector<MPI_Request>& reqs_send = get_expand_buffer_reqs_send();
  static Int do_comm_step = qlat::get_env_long_default(std::string("qlat_field_expand_buffer_AMD"), 0);
  Qassert(do_comm_step == 1 or do_comm_step == 0);
  if(do_comm_step == 1){
    mpi_waitall(reqs_recv);//receive done and write
    for(Int i=0;i<MAX_EXPAND_BUF;i++){
      timer.flops += bufs[i].get_flops();
      bufs[i].excute_copy(1);
    }
    mpi_waitall(reqs_send);//send    done and write
  }
  // copy after communications to avoid bad news on AMD GPUs...
  if(do_comm_step == 0){
    mpi_waitall(reqs_recv);//receive done and write
    mpi_waitall(reqs_send);//send    done and write
    for(Int i=0;i<MAX_EXPAND_BUF;i++){
      timer.flops += bufs[i].get_flops();
      bufs[i].excute_copy(1);
    }
  }
  reqs_send.resize(0);
  reqs_recv.resize(0);
  clean_expand_buffer_vecs();
}

template <class T>
inline expand_field_buffer& expand_buffer(const Geometry& geo, 
  const std::string& tag, T* res, const Int MULTI,  const Int GPU, Int& idx){
  std::vector<expand_field_buffer>& bufs = get_expand_buffer_full_list();
  // check if res is already in buffer or not
  static Int do_comm_default = qlat::get_env_long_default(std::string("qlat_field_expand_buffer_asyn"), 1);
  Qassert(do_comm_default == 0 or do_comm_default == 1);
  Int do_comm = 1 - do_comm_default;
  for(Int i=0;i<MAX_EXPAND_BUF;i++){
    if(bufs[i].initialized){
      if(&bufs[i].pres[0] == (int8_t*) &res[0]){
        do_comm = 1;
      }
    }
  }
  // if not found then do all mpi and return 0
  if(do_comm == 0){
    for(Int i=0;i<MAX_EXPAND_BUF;i++){
      if(!bufs[i].initialized){
        bufs[i].init(geo, tag, res, MULTI, GPU);
        idx = i;
        return bufs[i];
      }
    }
    do_comm = 1;
  }
  if(do_comm == 1)
  {
    expand_buffer_wait_mpi();
    idx = 0;
    bufs[idx].init(geo, tag, res, MULTI, GPU);
    return bufs[idx];
  }
  // qmessage("Compile with larger MAX_EXPAND_BUF!\n");
  // won't reach here to avoid warning
  Qassert(do_comm == 0);
  return bufs[0];
}

template <class M>
void refresh_expanded_GPUT_v0(M* res, const Geometry& geo, const Int MULTI, 
  const SetMarksField& set_marks_field = set_marks_field_all, const std::string& tag = std::string(""), Int GPU = 1, const QBOOL dummy = QTRUE){
  Qassert(sizeof(M) % sizeof(RealD) == 0);
  (void) dummy;
  //const Int MULTI = geo.multiplicity;
  const Long mpi_size = MULTI * sizeof(M)/sizeof(RealD);
  const Int multiplicity = 1;//always 1 for the buffers reuse
  //
  const CommPlan& plan = get_comm_plan(set_marks_field, tag, geo, multiplicity);
  const Long total_bytes =
      (plan.total_recv_size + plan.total_send_size) * sizeof(M);
  if (0 == total_bytes) {
    return;
  }
  TIMER_FLOPS("refresh_expanded_GPU");
  timer.flops += total_bytes * MULTI / 2;
  //
  std::vector<MPI_Request> reqs_send;
  std::vector<MPI_Request> reqs_recv;
  //
  qlat::vector_gpu<int8_t >& sbuf = qlat::get_vector_gpu_plan<int8_t >(0, std::string("general_buf0"), GPU);
  qlat::vector_gpu<int8_t >& rbuf = qlat::get_vector_gpu_plan<int8_t >(0, std::string("general_buf1"), GPU);
  expand_index_buf& ebuf = get_expand_index_buf_plan(geo, multiplicity, tag);
  Qassert(ebuf.pack_send.size() == 2*plan.total_send_size and ebuf.pack_recv.size() == 2*plan.total_recv_size);
  //
  const Long Nsend = plan.total_send_size ;
  const Long Nrecv = plan.total_recv_size ;
  //
  sbuf.resizeL(Nsend * MULTI * sizeof(M) / sizeof(int8_t));
  rbuf.resizeL(Nrecv * MULTI * sizeof(M) / sizeof(int8_t));
  //
  M* sP = (M*) &sbuf[0];
  M* rP = (M*) &rbuf[0];
  //
  ////setup reciev
  const Int mpi_tag = QLAT_VECTOR_UTILS_MPI_TAG;
  for (size_t i = 0; i < plan.recv_msg_infos.size(); ++i) {
    const CommMsgInfo& cmi = plan.recv_msg_infos[i]; 
    mpi_irecv(&rP[cmi.buffer_idx * MULTI], cmi.size * mpi_size, MPI_DOUBLE,
              cmi.id_node, mpi_tag, get_comm(), reqs_recv);
  }
  //
  //qlat::vector<long > pack_infos;
  Long* pack_send = (Long*) &ebuf.pack_send[0];
  Long* pack_recv = (Long*) &ebuf.pack_recv[0];
  //
  qGPU_for(isp, Nsend, GPU, {
    Long ri = pack_send[isp * 2 + 0] * MULTI;//offset with multivec
    Long si = pack_send[isp * 2 + 1] * MULTI;
    for(Int n=0;n<MULTI;n++){sP[ri + n] = res[si + n];}
  });
  //
  { 
    //TIMER("refresh_expanded-comm-init");
    for (size_t i = 0; i < plan.send_msg_infos.size(); ++i) {
      const CommMsgInfo& cmi = plan.send_msg_infos[i];
      mpi_isend(&sP[cmi.buffer_idx * MULTI], cmi.size * mpi_size, MPI_DOUBLE,
                cmi.id_node, mpi_tag, get_comm(), reqs_send);
    }
  }
  //
  mpi_waitall(reqs_recv);////receive done and write
  qGPU_for(isp, Nrecv, GPU, {
    const Long ri = pack_recv[isp * 2 + 0] * MULTI;
    const Long si = pack_recv[isp * 2 + 1] * MULTI;
    for(Int n=0;n<MULTI;n++){res[ri + n] = rP[si + n];}
  });
  //
  mpi_waitall(reqs_send);
}

template <class M>
void refresh_expanded_GPUT(M* res, const Geometry& geo, const Int MULTI, 
  const SetMarksField& set_marks_field = set_marks_field_all, const std::string& tag = std::string(""), Int GPU = 1, const QBOOL dummy = QTRUE){
  Qassert(sizeof(M) % sizeof(RealD) == 0);
  //const Int MULTI = geo.multiplicity;
  const Long mpi_size = MULTI * sizeof(M)/sizeof(RealD);
  const Int multiplicity = 1;//always 1 for the buffers reuse
  //
  const CommPlan& plan = get_comm_plan(set_marks_field, tag, geo, multiplicity);
  const Long total_bytes =
      (plan.total_recv_size + plan.total_send_size) * sizeof(M);
  if (0 == total_bytes) {
    if(dummy == QTRUE){expand_buffer_wait_mpi();}
    return;
  }
  //
  TIMER_FLOPS("refresh_expanded_GPU");
  timer.flops += total_bytes * MULTI / 2;
  //
  expand_index_buf& ebuf = get_expand_index_buf_plan(geo, multiplicity, tag);
  Qassert(ebuf.pack_send.size() == 2*plan.total_send_size and ebuf.pack_recv.size() == 2*plan.total_recv_size);
  //
  std::vector<MPI_Request>& reqs_send = get_expand_buffer_reqs_send();
  std::vector<MPI_Request>& reqs_recv = get_expand_buffer_reqs_recv();
  //
  Int idx_buf = -1;
  expand_field_buffer& efb = expand_buffer(geo, tag, res, MULTI, GPU, idx_buf);
  Qassert(idx_buf >= 0);
  //
  M* sP = (M*) efb.bs.data();
  M* rP = (M*) efb.br.data();
  //
  ////setup reciev
  const Int mpi_tag = QLAT_VECTOR_UTILS_MPI_TAG + idx_buf;
  for (size_t i = 0; i < plan.recv_msg_infos.size(); ++i) {
    const CommMsgInfo& cmi = plan.recv_msg_infos[i]; 
    mpi_irecv(&rP[cmi.buffer_idx * MULTI], cmi.size * mpi_size, MPI_DOUBLE,
              cmi.id_node, mpi_tag, get_comm(), reqs_recv);
  }
  //
  efb.excute_copy(0);
  //
  {
    for (size_t i = 0; i < plan.send_msg_infos.size(); ++i) {
      const CommMsgInfo& cmi = plan.send_msg_infos[i];
      mpi_isend(&sP[cmi.buffer_idx * MULTI], cmi.size * mpi_size, MPI_DOUBLE,
                cmi.id_node, mpi_tag, get_comm(), reqs_send);
    }
  }
  //
  if(dummy == QTRUE){expand_buffer_wait_mpi();}
}

template <class M>
void refresh_expanded_GPU(M* res, const Geometry& geo, const Int MULTI, const std::string& tag, const Int GPU = 1, const QBOOL dummy = QTRUE){
  std::vector<std::string > tagL = stringtolist(tag);
  //
  if(tagL.size() == 0 or tagL[0] == std::string(""))
  {
    // will do corners
    refresh_expanded_GPUT(res, geo, MULTI, set_marks_field_all, std::string(""), GPU, dummy);
  }else{
    for(unsigned int itag=0;itag<tagL.size();itag++){
      const std::string& t= tagL[itag];
      refresh_expanded_GPUT(res, geo, MULTI, set_marks_field_dir, t, GPU, dummy);
    }
  }
}

inline std::vector<std::string > expand_tags(){
  std::vector<std::string > tagL = {"dirmt", "dirmz", "dirmy", "dirmx", "dirx", "diry", "dirz", "dirt"};
  return tagL;
}

template <class M>
void refresh_expanded_GPU(M* res, const Geometry& geo, const Int MULTI, Int dir = -1000, const Int GPU = 1, const QBOOL dummy = QTRUE){
  Qassert(dir == -1000 or dir == -100 or (dir >= -3-1 and dir <= 3));
  // with corner
  if(dir == -1000){
    refresh_expanded_GPUT(res, geo, MULTI, set_marks_field_all, std::string(""), GPU, dummy);
  }
  // without corner
  if(dir == -100){
    refresh_expanded_GPUT(res, geo, MULTI, set_marks_field_dir, std::string(""), GPU, dummy);
  }
  //
  if(dir >= -3-1 and dir <= 3)
  {
    //std::string tagL[8] = {"dirmt", "dirmz", "dirmy", "dirmx", "dirx", "diry", "dirz", "dirt"};
    std::vector<std::string > tagL = expand_tags();
    refresh_expanded_GPUT(res, geo, MULTI, set_marks_field_dir, tagL[dir + 3 + 1], GPU, dummy);
  }
  //
  //left or right expand
  if(dir == 500 or dir == 501)
  {
    std::string tagL[2] = {"dirL", "dirR"};
    refresh_expanded_GPUT(res, geo, MULTI, set_marks_field_dir, tagL[dir % 500], GPU, dummy);
  }
}

template <class M>
void refresh_expanded_GPU(Field<M>& f, Int dir = -1000, const Int GPU = 1, const QBOOL dummy = QTRUE){
  M* res = (M*) qlat::get_data(f).data();
  refresh_expanded_GPU(res, f.geo(), f.multiplicity, dir, GPU, dummy);
}


}

#endif
