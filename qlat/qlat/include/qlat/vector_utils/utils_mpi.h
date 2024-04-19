// utils_mpi.h
// Gen Wang
// Jan. 2021

#ifndef UTILS_MPI_H
#define UTILS_MPI_H
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

namespace qlat
{

template<typename Iy>
void reduce_MPI_type(Iy num, MPI_Datatype& curr, unsigned int& size)
{
  if(num <= 0){curr = MPI_BYTE; size = 1;return;}
  //if(num%(sizeof(ComplexT<double>)) == 0){curr = MPI::DOUBLE_COMPLEX ; size=sizeof( ComplexT<double> );return;}
  //if(num%(sizeof(ComplexT<float >)) == 0){curr = MPI::COMPLEX        ; size=sizeof( ComplexT<float > );return;}

  if(num%(sizeof(std::int64_t )) == 0){curr = MPI_INT64_T ; size=sizeof(std::int64_t );return;}
  if(num%(sizeof(std::int32_t )) == 0){curr = MPI_INT32_T ; size=sizeof(std::int32_t );return;}
  if(num%(sizeof(std::int16_t )) == 0){curr = MPI_INT16_T ; size=sizeof(std::int16_t );return;}
  if(num%(sizeof(std::int8_t  )) == 0){curr = MPI_INT8_T  ; size=sizeof(std::int8_t  );return;}
}

// template<class M>
// struct get_MPI_Type{ static MPI_Datatype c=MPI_LOGICAL; static constexpr int size = 0;};
//
// template<>
// struct get_MPI_Type<char    >{      static MPI_Datatype c=MPI_CHAR; static constexpr int size = 1;};
// //template<>
// //struct get_MPI_Type<unsigned char>{ static MPI_Datatype c=MPI_UNSIGNED_CHAR; static constexpr int size = 1;};
// template<>
// struct get_MPI_Type<int8_t  >{ static MPI_Datatype c=MPI_INT8_T  ; static constexpr int size = 1;};
// template<>
// struct get_MPI_Type<int16_t >{ static MPI_Datatype c=MPI_INT16_T ; static constexpr int size = 1;};
// template<>
// struct get_MPI_Type<int32_t >{ static MPI_Datatype c=MPI_INT32_T ; static constexpr int size = 1;};
// template<>
// struct get_MPI_Type<int64_t >{ static MPI_Datatype c=MPI_INT64_T ; static constexpr int size = 1;};
// template<>
// struct get_MPI_Type<uint8_t >{ static MPI_Datatype c=MPI_UINT8_T ; static constexpr int size = 1;};
// template<>
// struct get_MPI_Type<uint16_t>{ static MPI_Datatype c=MPI_UINT16_T; static constexpr int size = 1;};
// template<>
// struct get_MPI_Type<uint32_t>{ static MPI_Datatype c=MPI_UINT32_T; static constexpr int size = 1;};
// template<>
// struct get_MPI_Type<uint64_t>{ static MPI_Datatype c=MPI_UINT64_T; static constexpr int size = 1;};
// template<>
// struct get_MPI_Type<RealF   >{ static MPI_Datatype c=MPI_FLOAT   ; static constexpr int size = 1;};
// template<>
// struct get_MPI_Type<RealD   >{ static MPI_Datatype c=MPI_DOUBLE  ; static constexpr int size = 1;};
// template<>
// struct get_MPI_Type<RealDD >{ static MPI_Datatype c=MPI_DOUBLE  ; static constexpr int size = 2;};

template <class M>
int set_mpi_type(MPI_Datatype& mpi_type)
// set mpi_type and return size
{
  mpi_type = MPI_LOGICAL;
  if (get_id_node() == 0) {
    printf("Type not found !!!! \n");
  }
  const int set_mpi_type_size = 0;
  qassert(set_mpi_type_size != 0);
  return 0;
}

template <>
inline int set_mpi_type<RealDD>(MPI_Datatype& mpi_type)
{
  mpi_type = MPI_DOUBLE;
  return 2;
}

template <>
inline int set_mpi_type<RealD>(MPI_Datatype& mpi_type)
{
  mpi_type = MPI_DOUBLE;
  return 1;
}

template <>
inline int set_mpi_type<RealF>(MPI_Datatype& mpi_type)
{
  mpi_type = MPI_FLOAT ;
  return 1;
}

template <>
inline int set_mpi_type<uint64_t>(MPI_Datatype& mpi_type)
{
  mpi_type = MPI_UINT64_T;
  return 1;
}

template <>
inline int set_mpi_type<uint32_t>(MPI_Datatype& mpi_type)
{
  mpi_type = MPI_UINT32_T;
  return 1;
}

template <>
inline int set_mpi_type<uint16_t>(MPI_Datatype& mpi_type)
{
  mpi_type = MPI_UINT16_T;
  return 1;
}

template <>
inline int set_mpi_type<uint8_t>(MPI_Datatype& mpi_type)
{
  mpi_type = MPI_UINT8_T;
  return 1;
}

template <>
inline int set_mpi_type<int64_t>(MPI_Datatype& mpi_type)
{
  mpi_type = MPI_INT64_T;
  return 1;
}

template <>
inline int set_mpi_type<int32_t>(MPI_Datatype& mpi_type)
{
  mpi_type = MPI_INT32_T;
  return 1;
}

template <>
inline int set_mpi_type<int16_t>(MPI_Datatype& mpi_type)
{
  mpi_type = MPI_INT16_T;
  return 1;
}

template <>
inline int set_mpi_type<int8_t>(MPI_Datatype& mpi_type)
{
  mpi_type = MPI_INT8_T;
  return 1;
}

template <>
inline int set_mpi_type<char>(MPI_Datatype& mpi_type)
{
  mpi_type = MPI_CHAR;
  return 1;
}

template<class M>
unsigned int get_mpi_type(MPI_Datatype& curr)
{
  using D = typename IsBasicDataType<M>::ElementaryType;
  const int Nsize = set_mpi_type<D>(curr);
  Qassert(sizeof(D) % Nsize == 0);
  const int size = sizeof(D) / Nsize;
  //if(size == 0)if(get_id_node()== 0){printf("Type not found !!!! \n");}
  //qassert(size != 0);
  return size;
}

// template<class M>
// unsigned int get_MPI_type(MPI_Datatype& curr)
// {
//   using D = typename IsBasicDataType<M>::ElementaryType;
//   curr = get_MPI_Type<D>::c;
//   Qassert(sizeof(D) % get_MPI_Type<D>::size == 0);
//   const int size = sizeof(D) / get_MPI_Type<D>::size;
//   //int size = get_MPI_Type<D>::size * sizeof(D);
//   //int size = get_MPI_Type<D>::size * sizeof(M) / sizeof(D);
//   if(size == 0)if(get_id_node()== 0){printf("Type not found !!!! \n");}
//   qassert(size != 0);
//   //size = size * sizeof(D);
//   return size;
//   
//   //curr = MPI_BYTE;unsigned int size = 1;
//   //DATA_TYPE typenum = get_data_type<M >();
//   //if(typenum == INVALID_TYPE){
//   //  if(get_id_node()== 0){printf("Type not found !!!! \n");}Qassert(false); return 0;
//   //}
// 
//   //int dtype = typenum % MAXTYPE;
//   //if(dtype <= FLOATIND + 3){
// 
//   //  size = typenum/MAXTYPE;
// 
//   //  if(dtype == 0){curr =  MPI_CHAR                 ; return size ;}
//   //  if(dtype == 1){curr =  MPI_UNSIGNED_CHAR        ; return size ;}
//   //  if(dtype == 2){curr =  MPI_SHORT                ; return size ;}
//   //  if(dtype == 3){curr =  MPI_UNSIGNED_SHORT       ; return size ;}
//   //  if(dtype == 4){curr =  MPI_INT                  ; return size ;}
//   //  if(dtype == 5){curr =  MPI_UNSIGNED             ; return size ;}
//   //  if(dtype == 6){curr =  MPI_LONG                 ; return size ;}
//   //  if(dtype == 7){curr =  MPI_UNSIGNED_LONG        ; return size ;}
//   //  if(dtype == 8){curr =  MPI_LONG_LONG            ; return size ;}
//   //  if(dtype == 9){curr =  MPI_UNSIGNED_LONG_LONG   ; return size ;}
//   //  if(dtype ==10){curr =  MPI_INT8_T               ; return size ;}
//   //  if(dtype ==11){curr =  MPI_UINT8_T              ; return size ;}
//   //  if(dtype ==12){curr =  MPI_INT16_T              ; return size ;}
//   //  if(dtype ==13){curr =  MPI_UINT16_T             ; return size ;}
//   //  if(dtype ==14){curr =  MPI_INT32_T              ; return size ;}
//   //  if(dtype ==15){curr =  MPI_UINT32_T             ; return size ;}
//   //  if(dtype ==16){curr =  MPI_INT64_T              ; return size ;}
//   //  if(dtype ==17){curr =  MPI_UINT64_T             ; return size ;}
// 
//   //  if(dtype ==FLOATIND+0){curr =  MPI_DOUBLE               ; return size ;}
//   //  if(dtype ==FLOATIND+1){curr =  MPI_FLOAT                ; return size ;}
//   //  if(dtype ==FLOATIND+2){curr =  MPI_C_DOUBLE_COMPLEX     ; return size ;}
//   //  if(dtype ==FLOATIND+3){curr =  MPI_C_FLOAT_COMPLEX      ; return size ;}
//   //}
//   //else{
//   //  if( get_data_type_is_double<M >()){curr = MPI_C_DOUBLE_COMPLEX; size = ComplexD_TYPE/MAXTYPE ;return size ;}
//   //  if(!get_data_type_is_double<M >()){curr = MPI_C_FLOAT_COMPLEX ; size = ComplexF_TYPE/MAXTYPE;return size ;}
//   //}
// 
//   //if(get_id_node()== 0){printf("Type not found !!!! \n");}Qassert(false);
//   //return 0;
// }

template<typename Ty>
void bcast_all_size(Ty *src, Long size, int root, int GPU=0, MPI_Comm* commp=NULL)
{
  TIMER("bcast_all_size");
  if(size == 0){return ;}
  (void) GPU;

  MPI_Datatype curr = MPI_DOUBLE;unsigned int M_size = sizeof(double);
  M_size = get_mpi_type<Ty >(curr);

  Qassert(sizeof(Ty)%M_size == 0);int M_fac = sizeof(Ty)/M_size;
  ////printf("size %5d %5d, type %d \n", int(size), int(fac), int(sizeof(Ty)));

  if(commp == NULL){
    MPI_Bcast(src, size * M_fac, curr, root, get_comm());
  }
  else{
    MPI_Bcast(src, size * M_fac, curr, root, *commp);
  }
}


template<typename Ty>
void sum_all_size(Ty *src,Ty *sav,Long size, int GPU=0, const MPI_Comm* commp=NULL)
{
  //TIMER("global sum sum_all_size");
  if(size == 0){return ;}
  if(qlat::get_num_node() == 1){
    if(src == sav){return;}
    if(src != sav){
      cpy_data_thread(sav, src, size, GPU);return;}
  }

  const int iomp = omp_get_thread_num(); ////each thread will have it's own buf
  Ty* buf_res;int GPU_set = GPU;
  #ifndef QLAT_USE_ACC
  GPU_set = 0;
  #endif
  VectorGPUKey gkey(size_t(size)*sizeof(Ty), ssprintf("sum_all_size_buf_%d", iomp), GPU_set); ////read buffers for global sum
  if(src == sav){
    const vector_gpu<char >& tmp = get_vector_gpu_plan<char >(gkey);
    buf_res = (Ty*) tmp.p;
  }else{buf_res = sav;}////small modify for pointers

  MPI_Datatype curr = MPI_DOUBLE;unsigned int M_size = sizeof(double);
  M_size = get_mpi_type<Ty >(curr);

  Qassert(sizeof(Ty)%M_size == 0);
  const int M_fac = sizeof(Ty)/M_size;
  //print0("mpi size %5d, M_fac %5d, Ty %5d \n", int(size), M_fac, int(sizeof(Ty)) );

  Ty* tem_src = NULL; Ty* tem_res = NULL;
  std::vector<Ty > tem_sHIP,tem_rHIP;
  bool do_copy = false;

  #ifdef __NO_GPU_DIRECT__
  #ifdef QLAT_USE_ACC
  if(GPU == 1){do_copy = true;}
  #endif
  #endif

  if(do_copy == false){tem_src = src;tem_res = buf_res;}
  if(do_copy == true ){
    tem_sHIP.resize(size);tem_rHIP.resize(size);

    cpy_data_thread(&tem_sHIP[0], src, size, 3);
    tem_src = &tem_sHIP[0];tem_res = &tem_rHIP[0];
  }
  
  if(commp == NULL){MPI_Allreduce(tem_src,tem_res, size * M_fac, curr, MPI_SUM, get_comm());}
  else{MPI_Allreduce(tem_src,tem_res, size * M_fac, curr, MPI_SUM, *commp);}

  if(do_copy == true){
    cpy_data_thread(buf_res, &tem_rHIP[0], size, 2);
  }


  if(src == sav)
  {
    cpy_data_thread(sav, buf_res, size, GPU);
  }
  if(src != sav){safe_free_vector_gpu_plan<char>(gkey);}
}

template<typename Ty>
void sum_all_size(Ty *src,Long size, int GPU=0, const MPI_Comm* commp=NULL)
{
  sum_all_size(src,src,size, GPU, commp);
}


inline void abort_sum(double flag, std::string stmp=std::string(""))
{
  sum_all_size(&flag,1);
  if(flag > 0)
  {
    abort_r(stmp);
  }
}

inline void fflush_MPI(){
  MPI_Barrier(get_comm());
  fflush(stdout);
}

//////"INT_MAX"
//////offset by number of char
template<typename Iy0, typename Iy1>
void MPI_Alltoallv_Send_Recv(char* src, Iy0* send, Iy1* spls, char* res, Iy0* recv, Iy1* rpls, const MPI_Comm& comm)
{
  int num_node;MPI_Comm_size(comm, &num_node);
  int id_node;MPI_Comm_rank(comm, &id_node);
  std::vector<MPI_Request> send_reqs(num_node);
  int mpi_tag = id_node;
  int c1 = 0;

  /////===get proper M_size
  MPI_Datatype curr = MPI_BYTE;unsigned int M_size = 1;unsigned int M_tem = 1;
  for(int n = 0; n < num_node; n++){
    if(send[n]!= 0){
      reduce_MPI_type(send[n], curr, M_tem);
      if(M_size == 1){M_size = M_tem;}
      else{if(M_tem != M_size){curr = MPI_BYTE;M_size = 1;break;}}
    }
  }
  /////

  for(int n = 0; n < num_node; n++){
    if(send[n]!=0){MPI_Isend(&src[spls[n]], int(send[n]/M_size), curr, n, mpi_tag + n, comm, &send_reqs[c1]);c1 += 1;}
  }

  for(int n = 0; n < num_node; n++){
    if(recv[n]!=0){MPI_Recv( &res[rpls[n]], int(recv[n]/M_size), curr, n, mpi_tag + n, comm, MPI_STATUS_IGNORE);}
  }
  if(c1!=0){MPI_Waitall(c1, send_reqs.data(), MPI_STATUS_IGNORE);}
}

template<typename Ty>
void MPI_Alltoallv_mode(Ty* src0, int* send, int* spls, Ty* res0, int* recv, int* rpls, const MPI_Comm& comm, int mode=0, int GPU = 0)
{
  (void)GPU;
  Ty* src = NULL;Ty* res = NULL;

  std::vector<Ty > tem_src,tem_res;
  bool do_copy = false;
  #ifdef __NO_GPU_DIRECT__
  #ifdef QLAT_USE_ACC
  if(GPU == 1){do_copy = true;}
  #endif
  #endif

  if(do_copy == false){src = src0; res = res0;}

  ////resize buffers
  Long max_src = 0;
  Long max_res = 0;
  if(do_copy == true){
    int num_node;MPI_Comm_size(comm, &num_node);
    for(int n = 0; n < num_node; n++){
      Long cur_size = spls[n]/sizeof(Ty) + send[n]/sizeof(Ty);
      if(cur_size > max_src){max_src = cur_size;}
      cur_size = rpls[n]/sizeof(Ty) + recv[n]/sizeof(Ty);
      if(cur_size > max_res){max_res = cur_size;}
    }

    tem_src.resize(max_src);tem_res.resize(max_res);
    cpy_data_thread(&tem_src[0], src0, max_src, 3);
    cpy_data_thread(&tem_res[0], res0, max_res, 3);
    src = &tem_src[0]; res = &tem_res[0];
  }

  if(mode == 0){
    MPI_Alltoallv(src, send, spls, MPI_CHAR,
                  res, recv, rpls, MPI_CHAR, comm);
  }
  if(mode == 1){
    MPI_Alltoallv_Send_Recv((char*) src, send, spls, (char*) res, recv, rpls, comm);

    //int num_node;MPI_Comm_size(comm, &num_node);
    //int id_node;MPI_Comm_rank(comm, &id_node);
    //std::vector<MPI_Request> send_reqs(num_node);
    //int mpi_tag = id_node;
    //int c1 = 0;
    //for(int n = 0; n < num_node; n++){
    //  if(send[n]!=0){MPI_Isend(&src[spls[n]/sizeof(Ty)], send[n], MPI_CHAR, n, mpi_tag + n, comm, &send_reqs[c1]);c1 += 1;}
    //}

    //for(int n = 0; n < num_node; n++){
    //  if(recv[n]!=0){MPI_Recv( &res[rpls[n]/sizeof(Ty)], recv[n], MPI_CHAR, n, mpi_tag + n, comm, MPI_STATUS_IGNORE);}
    //}
    //if(c1!=0){MPI_Waitall(c1, send_reqs.data(), MPI_STATUS_IGNORE);}
  }

  if(do_copy == true){cpy_data_thread(res0, &tem_res[0], max_res, 2);}
}


////Need add explanations
////sum all src or bcast each node data to others
template<typename Ty>
void Redistribute_all_Nt(Ty *src,Long size,const qlat::Geometry &geo, int GPU=0)
{
  if(qlat::get_num_node() == 1){return;}
  int Nt = geo.node_site[3];
  int Nmpi  = qlat::get_num_node();

  const Coordinate vg = geo.total_site();
  const int nt = vg[3];

  int mt = nt/Nt;
  if(mt != Nmpi){print0("Not supported !");Qassert(false);return;}

  /////int rank  = qlat::get_id_node();
  Long size_c = sizeof(Ty)*size/mt;

  std::vector<int > send,recv,spls,rpls;
  send.resize(Nmpi);
  recv.resize(Nmpi);
  spls.resize(Nmpi);
  rpls.resize(Nmpi);

  for(int ri=0;ri<Nmpi;ri++)
  {
    send[ri] = size_c;
    spls[ri] = size_c*ri;

    recv[ri] = size_c;
    rpls[ri] = size_c*ri;
  }

  //Ty* buf;
  //if(GPU == 0){buf = (Ty *)aligned_alloc_no_acc(size*sizeof(Ty));}
  qlat::vector_gpu<Ty > buf; buf.resize(size, GPU);

  {
  ////TIMER("MPI call CPU");
  MPI_Alltoallv(src   ,(int*) &send[0],(int*) &spls[0], MPI_CHAR,
            buf.data(),(int*) &recv[0],(int*) &rpls[0], MPI_CHAR, get_comm());
  }

  cpy_data_thread(src, buf.data(), size, GPU);
}

////sum all src or bcast each node data to others
template<typename Ty>
void Bcast_all_Nt(Ty *src,Long size,const qlat::Geometry &geo)
{
  if(qlat::get_num_node() == 1){return;}
  int Nt = geo.node_site[3];
  int Nmpi  = qlat::get_num_node();

  const Coordinate vg = geo.total_site();
  const int nt = vg[3];

  /////if Nmpi is not equal to mt
  if(nt/Nt != Nmpi){
    sum_all_size(src,size);
    return;
  }

  int mt = nt/Nt;
  int rank  = qlat::get_id_node();
  Long size_c = sizeof(Ty)*size/mt;

  unsigned short t0 = 0;
  {
    Coordinate xl = geo.coordinate_from_index(0);
    xl[3] = 0;
    Coordinate xg = geo.coordinate_g_from_l(xl);
    t0 = xg[3];
  }

  std::vector<int > send,recv,spls,rpls;
  send.resize(Nmpi);
  recv.resize(Nmpi);
  spls.resize(Nmpi);
  rpls.resize(Nmpi);

  std::fill(send.begin(), send.end(), 0);
  std::fill(recv.begin(), recv.end(), 0);
  std::fill(spls.begin(), spls.end(), 0);
  std::fill(rpls.begin(), rpls.end(), 0);

  for(int ti=0;ti<mt;ti++){
    int tini = ti*Nt;
    if(t0 == tini){
      for(int ri=0;ri<Nmpi;ri++)if(ri != rank)
      {
        send[ri] = size_c;
        spls[ri] = size_c*ti;
      }
    }

    if(t0 != tini){
      int ri_recv = ti;
      recv[ri_recv] = size_c;
      rpls[ri_recv] = size_c*ti;
    }
  }

  MPI_Alltoallv(src,(int*) &send[0],(int*) &spls[0], MPI_CHAR,
                src,(int*) &recv[0],(int*) &rpls[0], MPI_CHAR, get_comm());
}

///num to be zero for nodes
template<typename Ty>
void sum_value_mpi(Ty& num)
{
  //int Nmpi  = qlat::get_num_node();
  //int rank  = qlat::get_id_node();
  Ty buf = num;
  Long nvalue = 0;
  if(std::fabs(num) > 1e-30){nvalue = 1;}
  sum_all_size(&buf, 1);
  sum_all_size(&nvalue, 1);
  if(nvalue != 0){buf = buf/nvalue;}
  num = buf;
}

inline int get_mpi_id_node_close()
{
  int globalRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &globalRank);
  //Qassert(globalRank == get_id_node());
  // node local comm
  MPI_Comm nodeComm;
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, globalRank,
                      MPI_INFO_NULL, &nodeComm);

  // id within the node
  int localRank;
  MPI_Comm_rank(nodeComm, &localRank);
  //if (0 == get_id_node()) {
  //  Qassert(localRank == 0);
  //}
  //return 0;
  // number of process in this node
  int localSize;
  MPI_Comm_size(nodeComm, &localSize);
  // comm across node (each node select one process with the same local rank)
  MPI_Comm masterComm;
  MPI_Comm_split(MPI_COMM_WORLD, localRank, globalRank, &masterComm);
  // id across node
  int masterRank;
  MPI_Comm_rank(masterComm, &masterRank);
  // size of each master comm
  int masterSize;
  MPI_Comm_size(masterComm, &masterSize);
  // calculate number of node
  Long num_of_node = masterSize;
  MPI_Bcast(&num_of_node, 1, MPI_LONG, 0, nodeComm);
  // calculate id of node (master rank of the 0 local rank process)
  Long id_of_node = masterRank;
  MPI_Bcast(&id_of_node, 1, MPI_LONG, 0, nodeComm);
  Qassert(id_of_node < num_of_node);
  // calculate number of processes for each node
  std::vector<Long> n0(num_of_node, 0);
  std::vector<Long> n1(num_of_node, 0);
  n0[id_of_node] = 1;
  /////glb_sum(get_data(num_process_for_each_node));
  /////MPI_Allreduce(get_data(num_process_for_each_node), get_data(num_process_for_each_node), )

  //if(commp == NULL){MPI_Allreduce(tem_src,tem_res, size * fac, curr, MPI_SUM, get_comm());}
  MPI_Allreduce(&n0[0], &n1[0], num_of_node, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  ////printf("id_of_node %5d, total %5d \n", int(id_of_node), int(n1[id_of_node]));

  int id_node_local = localRank;
  for (Long i = 0; i < id_of_node; ++i) {
    id_node_local += n1[i];
  }

  //// calculate the number of master comm (the maximum in num_process_for_each_node)
  //Long num_of_master_comm = 0;
  //for (Long i = 0; i < (Long)num_process_for_each_node.size(); ++i) {
  //  if (num_process_for_each_node[i] > num_of_master_comm) {
  //    num_of_master_comm = num_process_for_each_node[i];
  //  }
  //}
  //// calculate the id of the master comm (same as local rank)
  //Long id_of_master_comm = localRank;
  //Qassert(id_of_master_comm < num_of_master_comm);
  //// calculate number of processes for each masterComm
  //std::vector<Long> num_process_for_each_master_comm(num_of_master_comm, 0);
  //num_process_for_each_master_comm[id_of_master_comm] = 1;
  //glb_sum(get_data(num_process_for_each_master_comm));
  //Qassert(num_process_for_each_master_comm[id_of_master_comm] == masterSize);
  //// calculate id_node_in_shuffle
  // calculate the list of id_node for each id_node_in_shuffle
  //std::vector<Long> list_long(get_num_node(), 0);
  //list_long[id_node_in_shuffle] = get_id_node();
  //glb_sum(get_data(list_long));
  //std::vector<int> list(get_num_node(), 0);
  //for (Long i = 0; i < get_num_node(); ++i) {
  //  list[i] = list_long[i];
  //}
  //// checking
  //Qassert(list[0] == 0);
  //for (Long i = 0; i < get_num_node(); ++i) {
  //  Qassert(0 <= list[i]);
  //  Qassert(list[i] < get_num_node());
  //  for (Long j = 0; j < i; ++j) {
  //    Qassert(list[i] != list[j]);
  //  }
  //}
  //return list;
  return id_node_local;
}


inline void geo_to_nv(const qlat::Geometry& geo, std::vector<int >& nv, std::vector<int > &Nv, std::vector<int > &mv)
{
  Nv.resize(4);nv.resize(4);mv.resize(4);
  for(int i=0;i<4;i++){Nv[i]=geo.node_site[i];nv[i] = geo.node_site[i] * geo.geon.size_node[i];}
  for(int i=0;i<4;i++){mv[i] = nv[i]/Nv[i];}
}
inline void geo_to_nv(const qlat::Geometry& geo, qlat::vector_acc<int >& nv, qlat::vector_acc<int > &Nv, qlat::vector_acc<int > &mv){
  Nv.resize(4);nv.resize(4);mv.resize(4);
  for(int i=0;i<4;i++){Nv[i]=geo.node_site[i];nv[i] = geo.node_site[i] * geo.geon.size_node[i];}
  for(int i=0;i<4;i++){mv[i] = nv[i]/Nv[i];}
}

//inline void setup_expand(const Geometry& geo, qlat::vector_acc<Long >& pack_send, qlat::vector_acc<Long >& pack_recv)
//{
//  const CommPlan& plan = get_comm_plan(set_marks_field_all, "", geo);
//  const Long Nsend = plan.total_send_size;
//  const Long Nrecv = plan.total_recv_size;
//  /////printf("setup %8d %8d \n", int(Nsend * 2), int(Nrecv * 2));
//
//  pack_send.resize( Nsend * 2 );
//  pack_recv.resize( Nrecv * 2 );
//
//  {
//  TIMER("refresh setup index");
//  Long cur = 0;
//  for (Long i = 0; i < (Long)plan.send_pack_infos.size(); ++i){
//    const CommPackInfo& cpi = plan.send_pack_infos[i];
//    #pragma omp parallel for
//    for(int off=0;off<cpi.size;off++){
//      pack_send[(cur+off)*2 + 0] = cpi.buffer_idx + off;
//      pack_send[(cur+off)*2 + 1] = cpi.offset + off;
//    }
//    cur += cpi.size;
//  }
//
//       cur = 0;
//  for (Long i = 0; i < (Long)plan.recv_pack_infos.size(); ++i){
//    const CommPackInfo& cpi = plan.recv_pack_infos[i];
//    #pragma omp parallel for
//    for(int off=0;off<cpi.size;off++){
//      pack_recv[(cur + off)*2 + 0] = cpi.offset + off;
//      pack_recv[(cur + off)*2 + 1] = cpi.buffer_idx + off;
//    }
//    cur += cpi.size;
//  }
//  }
//}
//
//struct expand_index_buf {
//  //Geometry geo; //make a copy of geo if needed
//  qlat::vector_acc<Long >  pack_send;
//  qlat::vector_acc<Long >  pack_recv;
//  //const long Nindex;
//  expand_index_buf()
//  {
//    pack_send.resize(0);
//    pack_recv.resize(0);
//  }
//
//  expand_index_buf(const Geometry& geo_)
//  {
//    //geo = geo_;
//    setup_expand(geo_, pack_send, pack_recv);
//  }
//
//  ~expand_index_buf()
//  {
//    pack_send.resize(0);
//    pack_recv.resize(0);
//  }
//};
//
///////buffers for expand index
//struct expand_index_Key {
//  Geometry geo;
//  //Coordinate total_site;
//  //Coordinate expansion_left ;
//  //Coordinate expansion_right;
//  expand_index_Key(const Geometry& geo_)
//  {
//    geo = geo_;
//    //total_site      = geo.total_site();
//    //expansion_left  = geo.expansion_left;
//    //expansion_right = geo.expansion_right;
//  }
//};
//
//inline bool compare_geo(const Geometry& g0, const Geometry& g1, const int with_multi = 1)
//{
//  int equal = 1;
//  if(g0.initialized           != g1.initialized ){ return 0; }
//  if(g0.eo                    != g1.eo ){ return 0; }
//  if(g0.is_only_local         != g1.is_only_local    ){ return 0; }
//
//  if(g0.geon                  != g1.geon ){ return 0; }
//  if(g0.node_site             != g1.node_site    ){ return 0; }
//  if(g0.node_site_expanded    != g1.node_site_expanded    ){ return 0; }
//
//  if(g0.expansion_left        != g1.expansion_left  ){ return 0; }
//  if(g0.expansion_right       != g1.expansion_right ){ return 0; }
//
//  //if(g0.total_site()    != g1.total_site()    ){ return 0; }
//
//  if(with_multi){
//    if(g0.multiplicity  != g1.multiplicity    ){ return 0; }
//  }
//  return equal;
//}
//
//inline bool Compare_geo(const Geometry& g0, const Geometry& g1)
//{
//  return compare_geo(g0, g1, 0);
//}
//
//inline bool compare_less(const Geometry& g0, const Geometry& g1, const int with_multi = 1)
//{
//  if(g0.total_site()    < g1.total_site()    ){  return true;}
//  if(g1.total_site()    < g0.total_site()    ){  return false;}
//
//  if(g0.expansion_left  < g1.expansion_left  ){  return true;}
//  if(g1.expansion_left  < g0.expansion_left  ){  return false;}
//
//  if(g0.expansion_right < g1.expansion_right ){  return true;}
//  if(g1.expansion_right < g0.expansion_right ){  return false;}
//
//  if(with_multi){
//    if(g0.multiplicity    < g1.multiplicity    ){  return true;}
//    if(g1.multiplicity    < g0.multiplicity    ){  return false;}
//  }
//
//  return false;
//
//}
//
//inline bool operator<(const expand_index_Key& x, const expand_index_Key& y)
//{
//  return compare_less(x.geo, y.geo, 1);
//  //if(x.geo.total_site < y.geo.total_site ){  return true;}
//  //if(y.geo.total_site < x.geo.total_site ){  return false;}
//
//  //if(x.geo.expansion_left < y.geo.expansion_left ){  return true;}
//  //if(y.geo.expansion_left < x.geo.expansion_left ){  return false;}
//
//  //if(x.geo.expansion_right < y.geo.expansion_right ){  return true;}
//  //if(y.geo.expansion_right < x.geo.expansion_right ){  return false;}
//
//  //if(x.geo.multiplicity < y.geo.multiplicity ){  return true;}
//  //if(y.geo.multiplicity < x.geo.multiplicity ){  return false;}
//
//
//  //return false;
//}
//
//inline Cache<expand_index_Key, expand_index_buf >& get_expand_index_buf_cache()
//{
//  static Cache<expand_index_Key, expand_index_buf > cache("expand_index_Key", 64);
//  return cache;
//}
//
//inline expand_index_buf& get_expand_index_buf_plan(const expand_index_Key& ekey)
//{
//  if (!get_expand_index_buf_cache().has(ekey)) {
//    //Geometry geo(ekey.total_site, 1);
//    //Geometry geo_ext = geo_resize(geo, ekey.expansion_left, ekey.expansion_right);
//    get_expand_index_buf_cache()[ekey] = expand_index_buf(ekey.geo);
//  }
//  expand_index_buf& buf = get_expand_index_buf_cache()[ekey];
//  return buf;
//}
//
//inline expand_index_buf& get_expand_index_buf_plan(const Geometry& geo)
//{
//  expand_index_Key ekey(geo);
//  return get_expand_index_buf_plan(ekey);
//}
//
//template <class M>
//void refresh_expanded_GPU(Field<M>& f, int GPU = 1)
//{
//  const CommPlan& plan = get_comm_plan(set_marks_field_all, "", f.geo());
//  const Long total_bytes =
//      (plan.total_recv_size + plan.total_send_size) * sizeof(M);
//  if (0 == total_bytes) {
//    return;
//  }
//  TIMER_FLOPS("refresh_expanded_GPU");
//  timer.flops += total_bytes / 2;
//
//  std::vector<MPI_Request> reqs_send;
//  std::vector<MPI_Request> reqs_recv;
//
//  qlat::vector_gpu<char >& sbuf = qlat::get_vector_gpu_plan<char >(0, std::string("general_buf0"), GPU);
//  qlat::vector_gpu<char >& rbuf = qlat::get_vector_gpu_plan<char >(0, std::string("general_buf1"), GPU);
//  //qlat::vector_gpu<char >& pack_buf = qlat::get_vector_gpu_plan<char >(0, std::string("general_buf2"), -1);
//  expand_index_buf& ebuf = get_expand_index_buf_plan(f.geo());
//  //expand_index_buf ebuf(f.geo());
//  //printf("send %8d %8d \n", int(ebuf.pack_send.size()), int( 2*plan.total_send_size));
//  //printf("recv %8d %8d \n", int(ebuf.pack_recv.size()), int( 2*plan.total_recv_size));
//  Qassert(ebuf.pack_send.size() == 2*plan.total_send_size and ebuf.pack_recv.size() == 2*plan.total_recv_size);
//  //expand_index_buf ebuf(f.geo());
//
//  const Long Nsend = plan.total_send_size;
//  const Long Nrecv = plan.total_recv_size;
//
//  sbuf.resizeL(Nsend * sizeof(M) / sizeof(char));
//  rbuf.resizeL(Nrecv * sizeof(M) / sizeof(char));
//  //pack_buf.resizeL( 2 * plan.total_send_size * sizeof(Long) / sizeof(char));
//  //pack_buf.resizeL( 2 * plan.total_recv_size * sizeof(Long) / sizeof(char));
//
//  M* sP = (M*) &sbuf[0];
//  M* rP = (M*) &rbuf[0];
//
//  Qassert(sizeof(M) % sizeof(double) == 0);
//  ////setup reciev
//  const int mpi_tag = 10;
//  for (size_t i = 0; i < plan.recv_msg_infos.size(); ++i) {
//    const CommMsgInfo& cmi = plan.recv_msg_infos[i]; 
//    mpi_irecv(&rP[cmi.buffer_idx], cmi.size * sizeof(M)/sizeof(double), MPI_DOUBLE,
//              cmi.id_node, mpi_tag, get_comm(), reqs_recv);
//  }
//
//  //qlat::vector_acc<long > pack_infos;
//  Long* pack_send = (Long*) &ebuf.pack_send[0];
//  Long* pack_recv = (Long*) &ebuf.pack_recv[0];
//
//  //{
//  //TIMER("refresh setup index");
//  //////pack_infos.resize(2 * plan.total_send_size);
//  //Long cur = 0;
//  //for (Long i = 0; i < (Long)plan.send_pack_infos.size(); ++i){
//  //  const CommPackInfo& cpi = plan.send_pack_infos[i];
//  //  #pragma omp parallel for
//  //  for(int off=0;off<cpi.size;off++){
//  //    pack_infos[(cur+off)*2 + 0] = cpi.buffer_idx + off;
//  //    pack_infos[(cur+off)*2 + 1] = cpi.offset + off;
//  //    //send_pack_infos[i*3 + 2] = cur;
//  //  }
//  //  cur += cpi.size;
//  //}
//  //}
//
//  qGPU_for(isp, Nsend, GPU, {
//    Long ri = pack_send[isp* 2 + 0];
//    Long si = pack_send[isp* 2 + 1];
//    sP[ri] = f.get_elem_offset(si);
//  });
//
//  { 
//    //TIMER("refresh_expanded-comm-init");
//    for (size_t i = 0; i < plan.send_msg_infos.size(); ++i) {
//      const CommMsgInfo& cmi = plan.send_msg_infos[i];
//      mpi_isend(&sP[cmi.buffer_idx], cmi.size * sizeof(M)/sizeof(double), MPI_DOUBLE,
//                cmi.id_node, mpi_tag, get_comm(), reqs_send);
//    }
//  }
//
//  //{
//  //TIMER("refresh setup index");
//  /////const Long Nsize = 0;
//  /////pack_infos.resize(2 * plan.total_recv_size);
//  //Long cur = 0;
//  //for (Long i = 0; i < (Long)plan.recv_pack_infos.size(); ++i){
//  //  const CommPackInfo& cpi = plan.recv_pack_infos[i];
//  //  #pragma omp parallel for
//  //  for(int off=0;off<cpi.size;off++){
//  //    pack_infos[(cur + off)*2 + 0] = cpi.offset + off;
//  //    pack_infos[(cur + off)*2 + 1] = cpi.buffer_idx + off;
//  //  }
//  //  cur += cpi.size;
//  //}
//  //}
//
//  mpi_waitall(reqs_recv);////receive done and write
//  qGPU_for(isp, Nrecv, GPU, {
//    const Long ri = pack_recv[isp* 2 + 0];
//    const Long si = pack_recv[isp* 2 + 1];
//    f.get_elem_offset(ri) = rP[si];
//  });
//
//  mpi_waitall(reqs_send);
//  //safe_free_vector_gpu_plan<char >(std::string("general_buf0"), GPU);
//  //safe_free_vector_gpu_plan<char >(std::string("general_buf1"), GPU);
//}


}

#endif
