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
#include "utils_field_gpu.h"

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

template<typename D>
bool IsBasicTypeReal(){
  // RealD, RealF, RealDD
  using M1 = typename IsBasicDataType<D>::ElementaryType;
  std::string type = IsBasicDataType<M1>::get_type_name();
  if(type == std::string("RealF") or type == std::string("RealD") or type == std::string("RealDD")){
    return true;
  }else{
    return false;
  }
}

template<typename Ty>
bool IsTypeComplex(){
  std::string type = IsBasicDataType<Ty>::get_type_name();
  if(type == std::string("ComplexD") or type == std::string("ComplexF") or type == std::string("ComplexDD")){
    return true;
  }else{
    return false;
  }
}

template <class M>
struct GetBasicDataType {
  static const std::string get_type_name() { return "unknown_type"; }
  using ElementaryType = M;
};

template <typename N >
struct GetBasicDataType<SelectedField<N > > {
  static const std::string get_type_name() {
    return IsBasicDataType<N>::get_type_name();
  }
  using ElementaryType = typename IsBasicDataType<N>::ElementaryType;
};

template <typename N >
struct GetBasicDataType<SelectedFieldG<N > > {
  static const std::string get_type_name() {
    return IsBasicDataType<N>::get_type_name();
  }
  using ElementaryType = typename IsBasicDataType<N>::ElementaryType;
};

template <typename N >
struct GetBasicDataType<Field<N > > {
  static const std::string get_type_name() {
    return IsBasicDataType<N>::get_type_name();
  }
  using ElementaryType = typename IsBasicDataType<N>::ElementaryType;
};

template <typename N >
struct GetBasicDataType<FieldG<N > > {
  static const std::string get_type_name() {
    return IsBasicDataType<N>::get_type_name();
  }
  using ElementaryType = typename IsBasicDataType<N>::ElementaryType;
};

template <typename N, Int civ >
struct GetBasicDataType<FieldM<N, civ > > {
  static const std::string get_type_name() {
    return IsBasicDataType<N>::get_type_name();
  }
  using ElementaryType = typename IsBasicDataType<N>::ElementaryType;
};

template <typename D>
struct GetBasicDataType<GaugeFieldT<D> > {
  static const std::string get_type_name() {
    return IsBasicDataType<ColorMatrixT<D>>::get_type_name();
  }
  using ElementaryType = typename IsBasicDataType<ColorMatrixT<D>>::ElementaryType;
};

template <typename D>
struct GetBasicDataType<Propagator4dT<D> > {
  static const std::string get_type_name() {
    return IsBasicDataType<WilsonMatrixT<D>>::get_type_name();
  }
  using ElementaryType = typename IsBasicDataType<WilsonMatrixT<D>>::ElementaryType;
};

template <class M>
qacc Long GetFieldSize(const Field<M>& f)
{
  return f.field.size() * sizeof(M);
}

template <class M>
qacc Long GetFieldSize(const SelectedField<M>& f)
{
  return f.field.size() * sizeof(M);
}

// template<class M>
// struct get_MPI_Type{ static MPI_Datatype c=MPI_LOGICAL; static constexpr Int size = 0;};
//
// template<>
// struct get_MPI_Type<char    >{      static MPI_Datatype c=MPI_CHAR; static constexpr Int size = 1;};
// //template<>
// //struct get_MPI_Type<unsigned char>{ static MPI_Datatype c=MPI_UNSIGNED_CHAR; static constexpr Int size = 1;};
// template<>
// struct get_MPI_Type<int8_t  >{ static MPI_Datatype c=MPI_INT8_T  ; static constexpr Int size = 1;};
// template<>
// struct get_MPI_Type<int16_t >{ static MPI_Datatype c=MPI_INT16_T ; static constexpr Int size = 1;};
// template<>
// struct get_MPI_Type<int32_t >{ static MPI_Datatype c=MPI_INT32_T ; static constexpr Int size = 1;};
// template<>
// struct get_MPI_Type<int64_t >{ static MPI_Datatype c=MPI_INT64_T ; static constexpr Int size = 1;};
// template<>
// struct get_MPI_Type<uint8_t >{ static MPI_Datatype c=MPI_UINT8_T ; static constexpr Int size = 1;};
// template<>
// struct get_MPI_Type<uint16_t>{ static MPI_Datatype c=MPI_UINT16_T; static constexpr Int size = 1;};
// template<>
// struct get_MPI_Type<uint32_t>{ static MPI_Datatype c=MPI_UINT32_T; static constexpr Int size = 1;};
// template<>
// struct get_MPI_Type<uint64_t>{ static MPI_Datatype c=MPI_UINT64_T; static constexpr Int size = 1;};
// template<>
// struct get_MPI_Type<RealF   >{ static MPI_Datatype c=MPI_FLOAT   ; static constexpr Int size = 1;};
// template<>
// struct get_MPI_Type<RealD   >{ static MPI_Datatype c=MPI_DOUBLE  ; static constexpr Int size = 1;};
// template<>
// struct get_MPI_Type<RealDD >{ static MPI_Datatype c=MPI_DOUBLE  ; static constexpr Int size = 2;};

template <class M>
Int set_mpi_type(MPI_Datatype& mpi_type)
// set mpi_type and return size
{
  mpi_type = MPI_LOGICAL;
  if (get_id_node() == 0) {
    printf("Type not found !!!! \n");
  }
  const Int set_mpi_type_size = 0;
  Qassert(set_mpi_type_size != 0);
  return 0;
}

template <>
inline Int set_mpi_type<RealDD>(MPI_Datatype& mpi_type)
{
  mpi_type = MPI_DOUBLE;
  return 2;
}

template <>
inline Int set_mpi_type<RealD>(MPI_Datatype& mpi_type)
{
  mpi_type = MPI_DOUBLE;
  return 1;
}

template <>
inline Int set_mpi_type<RealF>(MPI_Datatype& mpi_type)
{
  mpi_type = MPI_FLOAT ;
  return 1;
}

template <>
inline Int set_mpi_type<uint64_t>(MPI_Datatype& mpi_type)
{
  mpi_type = MPI_UINT64_T;
  return 1;
}

template <>
inline Int set_mpi_type<uint32_t>(MPI_Datatype& mpi_type)
{
  mpi_type = MPI_UINT32_T;
  return 1;
}

template <>
inline Int set_mpi_type<uint16_t>(MPI_Datatype& mpi_type)
{
  mpi_type = MPI_UINT16_T;
  return 1;
}

template <>
inline Int set_mpi_type<uint8_t>(MPI_Datatype& mpi_type)
{
  mpi_type = MPI_UINT8_T;
  return 1;
}

template <>
inline Int set_mpi_type<int64_t>(MPI_Datatype& mpi_type)
{
  mpi_type = MPI_INT64_T;
  return 1;
}

template <>
inline Int set_mpi_type<int32_t>(MPI_Datatype& mpi_type)
{
  mpi_type = MPI_INT32_T;
  return 1;
}

template <>
inline Int set_mpi_type<int16_t>(MPI_Datatype& mpi_type)
{
  mpi_type = MPI_INT16_T;
  return 1;
}

template <>
inline Int set_mpi_type<int8_t>(MPI_Datatype& mpi_type)
{
  mpi_type = MPI_INT8_T;
  return 1;
}

template <>
inline Int set_mpi_type<char>(MPI_Datatype& mpi_type)
{
  mpi_type = MPI_CHAR;
  return 1;
}

template<class M>
unsigned int get_mpi_type(MPI_Datatype& curr)
{
  using D = typename IsBasicDataType<M>::ElementaryType;
  const Int Nsize = set_mpi_type<D>(curr);
  Qassert(sizeof(D) % Nsize == 0);
  const Int size = sizeof(D) / Nsize;
  return size;
}

template<typename Ty>
void bcast_all_size(Ty *src, Long size, Int root, Int GPU=0, MPI_Comm* commp=NULL)
{
  TIMER("bcast_all_size");
  if(size == 0){return ;}
  (void) GPU;

  MPI_Datatype curr = MPI_DOUBLE;unsigned int M_size = sizeof(RealD);
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
void sum_all_size(Ty *src,Ty *sav,Long size, Int GPU=0, const MPI_Comm* commp=NULL)
{
  TIMERB("global sum sum_all_size");
  if(size == 0){return ;}
  if(qlat::get_num_node() == 1){
    if(src == sav){return;}
    if(src != sav){
      cpy_data_thread(sav, src, size, GPU);return;}
  }

  //omp_get_thread_num
  Qassert(omp_get_num_threads() == 1);// need to check whether it works...
  const Int iomp = omp_get_thread_num(); ////each thread will have it's own buf
  Ty* buf_res;int GPU_set = GPU;
  #ifndef QLAT_USE_ACC
  GPU_set = 0;
  #endif
  VectorGPUKey gkey(size_t(size)*sizeof(Ty)/sizeof(int8_t), ssprintf("sum_all_size_buf_%d", iomp), GPU_set); ////read buffers for global sum
  if(src == sav){
    const vector_gpu<int8_t >& tmp = get_vector_gpu_plan<int8_t >(gkey);
    buf_res = (Ty*) tmp.data();
  }else{buf_res = sav;}////small modify for pointers

  MPI_Datatype curr = MPI_DOUBLE;unsigned int M_size = sizeof(RealD);
  M_size = get_mpi_type<Ty >(curr);

  Qassert(sizeof(Ty)%M_size == 0);
  const Int M_fac = sizeof(Ty)/M_size;
  //qmessage("mpi size %5d, M_fac %5d, Ty %5d, int8_t %5d \n", int(size), M_fac, int(sizeof(Ty)), int(sizeof(int8_t)) );

  //bool do_copy = false;
  const bool do_copy = true   ; // always copy for global sum, AMD machine has some issue
  //H100 need to avoid direct gpu collective communications 
  //collective communications directly between GPUs need NCCL for this if necessary...
  //bool copy_sum_gpu  = GPU_set;
  bool copy_sum_gpu  = false;

  #ifdef __NO_GPU_DIRECT__
  #ifdef QLAT_USE_ACC
  if(GPU == 1){copy_sum_gpu = false;}
  #endif
  #endif

  //copy_sum_gpu = false;

  VectorGPUKey gkey0(0, ssprintf("sum_all_size_buf0_%d", iomp), copy_sum_gpu);
  VectorGPUKey gkey1(0, ssprintf("sum_all_size_buf1_%d", iomp), copy_sum_gpu);
  qlat::vector_gpu<int8_t >& tem_sHIP = get_vector_gpu_plan<int8_t >(gkey0);
  qlat::vector_gpu<int8_t >& tem_rHIP = get_vector_gpu_plan<int8_t >(gkey1);

  Ty* tem_src = NULL; Ty* tem_res = NULL;
  //#ifdef QLAT_USE_ACC
  //#ifdef __NVCC__
  //std::vector<Ty > tem_sHIP,tem_rHIP;
  //if(do_copy == true){tem_sHIP.resize(size);tem_rHIP.resize(size);}
  //#else
  //qlat::vector_gpu<Ty > tem_sHIP,tem_rHIP;
  //do_copy = true; // AMD copies always copy to global sum
  //if(do_copy == true){tem_sHIP.resizeL(size);tem_rHIP.resizeL(size);}
  //#endif
  //#else
  //std::vector<Ty > tem_sHIP,tem_rHIP;
  //if(do_copy == true){tem_sHIP.resize(size);tem_rHIP.resize(size);}
  //#endif

  if(do_copy == false){tem_src = src;tem_res = buf_res;}
  if(do_copy == true ){
    tem_sHIP.resizeL(size * sizeof(Ty)/sizeof(int8_t));
    tem_rHIP.resizeL(size * sizeof(Ty)/sizeof(int8_t));
    cpy_GPU((Ty*) tem_sHIP.data(), src, size, copy_sum_gpu, GPU_set);
    //cpy_data_thread(tem_sHIP.data(), src, size, 3);
    tem_src = (Ty*) tem_sHIP.data();tem_res = (Ty*) tem_rHIP.data();
  }
  
  if(commp == NULL){MPI_Allreduce(tem_src,tem_res, size * M_fac, curr, MPI_SUM, get_comm());}
  else{MPI_Allreduce(tem_src,tem_res, size * M_fac, curr, MPI_SUM, *commp);}

  //MPI_Request request;
  //if(commp == NULL){MPI_Iallreduce(tem_src,tem_res, size * M_fac, curr, MPI_SUM, get_comm(), &request);}
  //else{MPI_Iallreduce(tem_src,tem_res, size * M_fac, curr, MPI_SUM, *commp, &request);}
  //MPI_Wait(&request, MPI_STATUS_IGNORE);

  if(do_copy == true){
    cpy_GPU(sav, (Ty*) tem_rHIP.data(), size, GPU_set, copy_sum_gpu);
    //cpy_data_thread(buf_res, tem_rHIP.data(), size, 2);
  }
  //qmessage("===check iomp %d, GPU_set %d, do_copy %d, cpu %d \n", iomp, GPU_set, int(do_copy), int(copy_sum_gpu));
  if(src == sav and do_copy == false)
  {
    cpy_GPU(sav, buf_res, size, GPU, GPU_set);
    //cpy_data_thread(sav, buf_res, size, GPU);
  }

  if(src != sav){
    safe_free_vector_gpu_plan<int8_t>(gkey);
    safe_free_vector_gpu_plan<int8_t>(gkey0);
    safe_free_vector_gpu_plan<int8_t>(gkey1);
  }
}

template<typename Ty>
void sum_all_size(Ty *src,Long size, Int GPU=0, const MPI_Comm* commp=NULL)
{
  sum_all_size(src,src,size, GPU, commp);
}


inline void abort_sum(RealD flag, std::string stmp=std::string(""))
{
  sum_all_size(&flag,1);
  if(flag > 0)
  {
    abort_r(stmp);
  }
}

inline void fflush_MPI(){
  qacc_barrier(dummy);
  MPI_Barrier(get_comm());
  fflush(stdout);
}

//////"INT_MAX"
//////offset by number of int8_t
template<typename Iy0, typename Iy1>
void MPI_Alltoallv_Send_Recv(char* src, Iy0* send, Iy1* spls, char* res, Iy0* recv, Iy1* rpls, const MPI_Comm& comm)
{
  Int num_node;MPI_Comm_size(comm, &num_node);
  Int id_node;MPI_Comm_rank(comm, &id_node);
  std::vector<MPI_Request> send_reqs(num_node);
  std::vector<MPI_Request> recv_reqs(num_node);
  //int mpi_tag = id_node;
  const Int mpi_tag = QLAT_VECTOR_UTILS_MPI_TAG;

  /////===get proper M_size
  MPI_Datatype curr = MPI_BYTE;unsigned int M_size = 1;unsigned int M_tem = 1;
  for(Int n = 0; n < num_node; n++){
    if(send[n]!= 0){
      reduce_MPI_type(send[n], curr, M_tem);
      if(M_size == 1){M_size = M_tem;}
      else{if(M_tem != M_size){curr = MPI_BYTE;M_size = 1;break;}}
    }
  }
  /////

  Int c1 = 0;
  Int c2 = 0;
  for(Int n = 0; n < num_node; n++){
    if(send[n]!=0){MPI_Isend(&src[spls[n]], int(send[n]/M_size), curr, n, mpi_tag, comm, &send_reqs[c1]);c1 += 1;}
  }

  for(Int n = 0; n < num_node; n++){
    if(recv[n]!=0){MPI_Irecv( &res[rpls[n]], int(recv[n]/M_size), curr, n, mpi_tag, comm, &recv_reqs[c2]);c2 += 1;}
  }
  if(c2!=0){MPI_Waitall(c2, recv_reqs.data(), MPI_STATUS_IGNORE);}
  if(c1!=0){MPI_Waitall(c1, send_reqs.data(), MPI_STATUS_IGNORE);}
}

template<typename Ty>
void MPI_Alltoallv_mode(Ty* src0, Int* send, Int* spls, Ty* res0, Int* recv, Int* rpls, const MPI_Comm& comm, Int mode=0, Int GPU = 0)
{
  (void)GPU;
  Ty* src = NULL;Ty* res = NULL;

  std::vector<Ty > tem_src,tem_res;
  //bool do_copy = false;
  //collactive behavior need to be done on CPU ...
  bool do_copy = true;
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
    Int num_node;MPI_Comm_size(comm, &num_node);
    for(Int n = 0; n < num_node; n++){
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
  }

  if(do_copy == true){cpy_data_thread(res0, &tem_res[0], max_res, 2);}
}

////Need add explanations
////sum all src or bcast each node data to others
template<typename Ty>
void Redistribute_all_Nt(Ty *src,Long size,const qlat::Geometry& geo, Int GPU=0)
{
  if(qlat::get_num_node() == 1){return;}
  Int Nt = geo.node_site[3];
  Int Nmpi  = qlat::get_num_node();

  const Coordinate vg = geo.total_site();
  const Int nt = vg[3];

  Int mt = nt/Nt;
  if(mt != Nmpi){qmessage("Not supported !");Qassert(false);return;}

  /////int rank  = qlat::get_id_node();
  Long size_c = sizeof(Ty)*size/mt;

  std::vector<Int > send,recv,spls,rpls;
  send.resize(Nmpi);
  recv.resize(Nmpi);
  spls.resize(Nmpi);
  rpls.resize(Nmpi);

  for(Int ri=0;ri<Nmpi;ri++)
  {
    send[ri] = size_c;
    spls[ri] = size_c*ri;

    recv[ri] = size_c;
    rpls[ri] = size_c*ri;
  }

  //Ty* buf;
  //if(GPU == 0){buf = (Ty *)aligned_alloc_no_acc(size*sizeof(Ty));}
  vector_gpu<Ty > buf; buf.resize(size, GPU);

  {
  ////TIMER("MPI call CPU");
  MPI_Alltoallv(src   ,(int*) &send[0],(int*) &spls[0], MPI_CHAR,
            buf.data(),(int*) &recv[0],(int*) &rpls[0], MPI_CHAR, get_comm());
  }

  cpy_data_thread(src, buf.data(), size, GPU);
}

////sum all src or bcast each node data to others
template<typename Ty>
void Bcast_all_Nt(Ty *src,Long size,const Geometry& geo)
{
  if(qlat::get_num_node() == 1){return;}
  Int Nt = geo.node_site[3];
  Int Nmpi  = qlat::get_num_node();

  const Coordinate vg = geo.total_site();
  const Int nt = vg[3];

  /////if Nmpi is not equal to mt
  if(nt/Nt != Nmpi){
    sum_all_size(src,size);
    return;
  }

  Int mt = nt/Nt;
  Int rank  = qlat::get_id_node();
  Long size_c = sizeof(Ty)*size/mt;

  unsigned short t0 = 0;
  {
    Coordinate xl = geo.coordinate_from_index(0);
    xl[3] = 0;
    Coordinate xg = geo.coordinate_g_from_l(xl);
    t0 = xg[3];
  }

  std::vector<Int > send,recv,spls,rpls;
  send.resize(Nmpi);
  recv.resize(Nmpi);
  spls.resize(Nmpi);
  rpls.resize(Nmpi);

  std::fill(send.begin(), send.end(), 0);
  std::fill(recv.begin(), recv.end(), 0);
  std::fill(spls.begin(), spls.end(), 0);
  std::fill(rpls.begin(), rpls.end(), 0);

  for(Int ti=0;ti<mt;ti++){
    Int tini = ti*Nt;
    if(t0 == tini){
      for(Int ri=0;ri<Nmpi;ri++)if(ri != rank)
      {
        send[ri] = size_c;
        spls[ri] = size_c*ti;
      }
    }

    if(t0 != tini){
      Int ri_recv = ti;
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

inline Int get_mpi_id_node_close()
{
  Int globalRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &globalRank);
  //Qassert(globalRank == get_id_node());
  // node local comm
  MPI_Comm nodeComm;
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, globalRank,
                      MPI_INFO_NULL, &nodeComm);

  // id within the node
  Int localRank;
  MPI_Comm_rank(nodeComm, &localRank);
  //if (0 == get_id_node()) {
  //  Qassert(localRank == 0);
  //}
  //return 0;
  // number of process in this node
  Int localSize;
  MPI_Comm_size(nodeComm, &localSize);
  // comm across node (each node select one process with the same local rank)
  MPI_Comm masterComm;
  MPI_Comm_split(MPI_COMM_WORLD, localRank, globalRank, &masterComm);
  // id across node
  Int masterRank;
  MPI_Comm_rank(masterComm, &masterRank);
  // size of each master comm
  Int masterSize;
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

  Int id_node_local = localRank;
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
  //std::vector<Int> list(get_num_node(), 0);
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


}

#endif
