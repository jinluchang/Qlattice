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
#include<type_traits>

#include <iterator>
#include "utils_read_txt.h"
#include "utils_vector_GPU.h"

namespace qlat
{

template<typename Iy>
void reduce_MPI_type(Iy num, MPI_Datatype& curr, unsigned int& size)
{
  if(num <= 0){curr = MPI_BYTE; size = 1;return;}
  //if(num%(sizeof(std::complex<double>)) == 0){curr = MPI::DOUBLE_COMPLEX ; size=sizeof( std::complex<double> );return;}
  //if(num%(sizeof(std::complex<float >)) == 0){curr = MPI::COMPLEX        ; size=sizeof( std::complex<float > );return;}

  if(num%(sizeof(std::int64_t )) == 0){curr = MPI_INT64_T ; size=sizeof(std::int64_t );return;}
  if(num%(sizeof(std::int32_t )) == 0){curr = MPI_INT32_T ; size=sizeof(std::int32_t );return;}
  if(num%(sizeof(std::int16_t )) == 0){curr = MPI_INT16_T ; size=sizeof(std::int16_t );return;}
  if(num%(sizeof(std::int8_t  )) == 0){curr = MPI_INT8_T  ; size=sizeof(std::int8_t  );return;}
}

template<class M>
unsigned int get_MPI_type(MPI_Datatype& curr)
{
  curr = MPI_BYTE;unsigned int size = 1;
  DATA_TYPE typenum = get_data_type<M >();
  if(typenum == INVALID_TYPE){
    if(get_id_node()== 0){printf("Type not found !!!! \n");}qassert(false); return 0;
  }

  int dtype = typenum % MAXTYPE;
  if(dtype <= FLOATIND + 3){

    size = typenum/MAXTYPE;

    if(dtype == 0){curr =  MPI_CHAR                 ; return size ;}
    if(dtype == 1){curr =  MPI_UNSIGNED_CHAR        ; return size ;}
    if(dtype == 2){curr =  MPI_SHORT                ; return size ;}
    if(dtype == 3){curr =  MPI_UNSIGNED_SHORT       ; return size ;}
    if(dtype == 4){curr =  MPI_INT                  ; return size ;}
    if(dtype == 5){curr =  MPI_UNSIGNED             ; return size ;}
    if(dtype == 6){curr =  MPI_LONG                 ; return size ;}
    if(dtype == 7){curr =  MPI_UNSIGNED_LONG        ; return size ;}
    if(dtype == 8){curr =  MPI_LONG_LONG            ; return size ;}
    if(dtype == 9){curr =  MPI_UNSIGNED_LONG_LONG   ; return size ;}
    if(dtype ==10){curr =  MPI_INT8_T               ; return size ;}
    if(dtype ==11){curr =  MPI_UINT8_T              ; return size ;}
    if(dtype ==12){curr =  MPI_INT16_T              ; return size ;}
    if(dtype ==13){curr =  MPI_UINT16_T             ; return size ;}
    if(dtype ==14){curr =  MPI_INT32_T              ; return size ;}
    if(dtype ==15){curr =  MPI_UINT32_T             ; return size ;}
    if(dtype ==16){curr =  MPI_INT64_T              ; return size ;}
    if(dtype ==17){curr =  MPI_UINT64_T             ; return size ;}
    
    if(dtype ==FLOATIND+0){curr =  MPI_DOUBLE               ; return size ;}
    if(dtype ==FLOATIND+1){curr =  MPI_FLOAT                ; return size ;}
    if(dtype ==FLOATIND+2){curr =  MPI_C_DOUBLE_COMPLEX     ; return size ;}
    if(dtype ==FLOATIND+3){curr =  MPI_C_FLOAT_COMPLEX      ; return size ;}
  }
  else{
    if( get_data_type_is_double<M >()){curr = MPI_C_DOUBLE_COMPLEX; size = Complex_TYPE/MAXTYPE ;return size ;}
    if(!get_data_type_is_double<M >()){curr = MPI_C_FLOAT_COMPLEX ; size = ComplexF_TYPE/MAXTYPE;return size ;}
  }

  if(get_id_node()== 0){printf("Type not found !!!! \n");}qassert(false);
  return 0;

}

template<typename Ty>
void bcast_all_size(Ty *src, long size, int root, int GPU=0, MPI_Comm* commp=NULL)
{
  TIMER("bcast_all_size");
  if(size == 0){return ;}
  (void) GPU;

  MPI_Datatype curr = MPI_DOUBLE;unsigned int M_size = sizeof(double);
  M_size = get_MPI_type<Ty >(curr);

  qassert(sizeof(Ty)%M_size == 0);int fac = sizeof(Ty)/M_size;
  ////printf("size %5d %5d, type %d \n", int(size), int(fac), int(sizeof(Ty)));

  if(commp == NULL){
    MPI_Bcast(src, size * fac, curr, root, get_comm());
  }
  else{
    MPI_Bcast(src, size * fac, curr, root, *commp);
  }
}


template<typename Ty>
void sum_all_size(Ty *src,Ty *sav,long size, int GPU=0, MPI_Comm* commp=NULL)
{
  TIMER("global sum sum_all_size");
  if(size == 0){return ;}
  if(qlat::get_num_node() == 1){
    if(src == sav){return;}
    if(src != sav){
      cpy_data_thread(sav, src, size, GPU, true);return;}
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
  M_size = get_MPI_type<Ty >(curr);

  qassert(sizeof(Ty)%M_size == 0);int fac = sizeof(Ty)/M_size;


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

    cpy_data_thread(&tem_sHIP[0], src, size, 3, true);
    tem_src = &tem_sHIP[0];tem_res = &tem_rHIP[0];
  }
  
  if(commp == NULL){MPI_Allreduce(tem_src,tem_res, size * fac, curr, MPI_SUM, get_comm());}
  else{MPI_Allreduce(tem_src,tem_res, size * fac, curr, MPI_SUM, *commp);}

  if(do_copy == true){
    cpy_data_thread(buf_res, &tem_rHIP[0], size, 2, true);
  }


  if(src == sav)
  {
    cpy_data_thread(sav, buf_res, size, GPU, true);
  }
  if(src != sav){safe_free_vector_gpu_plan<char>(gkey);}
}

template<typename Ty>
void sum_all_size(Ty *src,long size, int GPU=0, MPI_Comm* commp=NULL)
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
void MPI_Alltoallv_Send_Recv(char* src, Iy0* send, Iy1* spls, char* res, Iy0* recv, Iy1* rpls, MPI_Comm& comm)
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
void MPI_Alltoallv_mode(Ty* src0, int* send, int* spls, Ty* res0, int* recv, int* rpls, MPI_Comm& comm, int mode=0, int GPU = 0)
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
  long max_src = 0;
  long max_res = 0;
  if(do_copy == true){
    int num_node;MPI_Comm_size(comm, &num_node);
    for(int n = 0; n < num_node; n++){
      long cur_size = spls[n]/sizeof(Ty) + send[n]/sizeof(Ty);
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
void Redistribute_all_Nt(Ty *src,long size,const qlat::Geometry &geo, int GPU=0)
{
  if(qlat::get_num_node() == 1){return;}
  int Nt = geo.node_site[3];
  int Nmpi  = qlat::get_num_node();

  const Coordinate vg = geo.total_site();
  const int nt = vg[3];

  int mt = nt/Nt;
  if(mt != Nmpi){print0("Not supported !");qassert(false);return;}

  /////int rank  = qlat::get_id_node();
  long size_c = sizeof(Ty)*size/mt;

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
  //if(GPU == 1){gpuMalloc(buf, size, Ty);}
  qlat::vector_gpu<Ty > buf; buf.resize(size, GPU);

  {
  ////TIMER("MPI call CPU");
  MPI_Alltoallv(src   ,(int*) &send[0],(int*) &spls[0], MPI_CHAR,
            buf.data(),(int*) &recv[0],(int*) &rpls[0], MPI_CHAR, get_comm());
  }

  cpy_data_thread(src, buf.data(), size, GPU, true);
}

////sum all src or bcast each node data to others
template<typename Ty>
void Bcast_all_Nt(Ty *src,long size,const qlat::Geometry &geo)
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
  long size_c = sizeof(Ty)*size/mt;

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
  long nvalue = 0;
  if(std::fabs(num) > 1e-30){nvalue = 1;}
  sum_all_size(&buf, 1);
  sum_all_size(&nvalue, 1);
  if(nvalue != 0){buf = buf/nvalue;}
  num = buf;
}




}

#endif
