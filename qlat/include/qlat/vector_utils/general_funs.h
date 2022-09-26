// general_funs.h
// Gen Wang
// Jan. 2021

#ifndef GENERAL_FUNS_H
#define GENERAL_FUNS_H
#pragma once


#include <string.h>
#include <sys/resource.h>
#include <mpi.h>
#include <time.h>
#include <typeinfo>

#include "utils_float_type.h"
#include<type_traits>

#include <iterator>
#include <sys/sysinfo.h>
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
void sum_all_size(Ty *src,Ty *sav,long size, int GPU=0, MPI_Comm* commp=NULL)
{
  TIMER("global sum sum_all_size");
  if(size == 0){return ;}
  qlat::vector_gpu<Ty > res;//// buf.resize(size, GPU);
  ///Ty *res;/////qlat::vector<Ty >buf;
  if(src == sav){
    res.resize(size, GPU);
    //if(GPU == 0){res = (Ty *)aligned_alloc_no_acc(size*sizeof(Ty));}
    //else{gpuMalloc(res, size, Ty);}
  }else{res.p = sav;}////small modify for pointers
  if(qlat::get_num_node() == 1){
    if(src == sav){return;}
    if(src != sav){
      cpy_data_thread(sav, src, size, GPU, true);
      ////#ifdef QLAT_USE_ACC
      ////if(GPU==0)memcpy(sav,src,size*sizeof(Ty));
      ////if(GPU==1){cudaMemcpy(sav, src, size*sizeof(Ty), cudaMemcpyDeviceToDevice);}
      ////#else
      ////memcpy(sav, src, size*sizeof(Ty));
      ////#endif
    return;}
  }

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

  if(do_copy == false){tem_src = src;tem_res = res.data();}
  if(do_copy == true ){
    tem_sHIP.resize(size);tem_rHIP.resize(size);

    cpy_data_thread(&tem_sHIP[0], src, size, 3, true);
    tem_src = &tem_sHIP[0];tem_res = &tem_rHIP[0];
  }
  
  if(commp == NULL){MPI_Allreduce(tem_src,tem_res, size * fac, curr, MPI_SUM, get_comm());}
  else{MPI_Allreduce(tem_src,tem_res, size * fac, curr, MPI_SUM, *commp);}

  if(do_copy == true){
    cpy_data_thread(res.data(), &tem_rHIP[0], size, 2, true);
  }


  if(src == sav)
  {
    cpy_data_thread(sav, res.data(), size, GPU, true);
    //#ifdef QLAT_USE_ACC
    //if(GPU==0){memcpy(sav,res.data(),size*sizeof(Ty));}  ////free(res);
    //if(GPU==1){cudaMemcpy(sav, res.data(), size*sizeof(Ty), cudaMemcpyDeviceToDevice);} ///gpuFree(res);res = NULL;
    //#else
    //memcpy(sav,res.data(),size*sizeof(Ty));////free(res);res = NULL;
    //#endif
  }
  if(src != sav){res.p = NULL;}
}

template<typename Ty>
void sum_all_size(Ty *src,long size, int GPU=0, MPI_Comm* commp=NULL)
{
  sum_all_size(src,src,size, GPU, commp);
}

template<typename Ty>
void Bcast_all_Nt(Ty *src,long size,const qlat::Geometry &geo)
{
  if(qlat::get_num_node() == 1){return;}
  int Nt = geo.node_site[3];
  int Nmpi  = qlat::get_num_node();

  const Coordinate vg = geo.total_site();
  const int nt = vg[3];

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

////Need add explanations
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

  #ifdef QLAT_USE_ACC
  if(GPU == 0){
    #pragma omp parallel for
    for(long isp=0;isp<size;isp++){src[isp] = buf[isp];}
    ///delete [] buf;
    //free(buf);
  }
  if(GPU == 1){
    /////qacc_for(isp, size,{ src[isp] = buf[isp];});
    cudaMemcpy(src, buf.data(), size*sizeof(Ty), cudaMemcpyDeviceToDevice);
    ////gpuFree(buf);
  }
  #else
  #pragma omp parallel for
  for(long isp=0;isp<size;isp++){src[isp] = buf[isp];}
  //delete [] buf;
  //free(buf);
  #endif

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

////Only cpu verstion
////flag = 1 --> biva * sizeF * civ * size_inner --> biva * civ * sizeF * size_inner
inline void reorder_civ(char* src,char* res,int biva,int civ,size_t sizeF,int flag,int size_inner)
{
  //TIMER("reorder_civ vectors char");
  if(biva == 0 or civ == 0 or sizeF == 0 or size_inner == 0){return ;}

  std::vector<char* > psrc;psrc.resize(civ);
  std::vector<char > tmp;tmp.resize(biva*sizeF*civ*size_inner);
  if(size_inner <= 1){abort_r("size_innter too small ! \n");return;}

  if(flag == 1){memcpy((char*)&tmp[0],(char*)&src[0],sizeof(char)*biva*sizeF*civ*size_inner);}
 
  for(size_t bi=0;bi<size_t(biva);bi++)
  {
    for(int ci=0;ci<civ;ci++)
    {
      if(flag==0)psrc[ci] = &src[(bi*sizeF*civ + ci*sizeF + 0)*size_inner];
      if(flag==1)psrc[ci] = &res[(bi*sizeF*civ + ci*sizeF + 0)*size_inner];
    }

    #pragma omp parallel for
    for(LInt si=0;si<sizeF;si++)
    for(int ci=0;ci<civ;ci++)
    {
      if(flag==0){
        memcpy(&tmp[(bi*sizeF*civ + si*civ + ci)*size_inner],&psrc[ci][si*size_inner],sizeof(char)*size_inner);
      }
      if(flag==1){
        memcpy(&psrc[ci][si*size_inner],&tmp[(bi*sizeF*civ + si*civ + ci)*size_inner],sizeof(char)*size_inner);
      }
    }
  }
 
  if(flag == 0){memcpy((char*)&res[0],(char*)&tmp[0],biva*sizeF*civ*size_inner);}
}

///flag = 1 --> biva * sizeF * civ * size_inner --> biva * civ * sizeF * size_inner
#ifdef QLAT_USE_ACC
template <typename Ty, bool flag, int Threads, int Biva>
__global__ void move_index_global(Ty* src, Ty* res, long sizeF, int civ, int inner)
{
  __shared__ Ty buf[Threads*Biva];

  int    tid = threadIdx.x;
  long s0    = blockIdx.x*blockDim.x;

  int Total = Threads*civ*inner;
  if(s0 + Threads > sizeF){Total = (sizeF - s0) * civ*inner;}

  int nB    = (Total + Threads-1)/Threads;
  int nC    = (Total + Biva*Threads-1)/(Biva*Threads);

  int ci, si, i0;
  long z0 = 0;long off = 0;long off1 = 0;
  for(int ni=0;ni < nC; ni++)
  {
    if(z0 >= Total){break;}
    if(flag){
    off = z0 + tid;
    for(int xi=0;xi<Biva;xi++)
    {
      if(off < Total){buf[xi*Threads + tid] = src[s0*civ*inner + off];off += Threads;}
    }
    __syncthreads();
    }

    off = tid;
    for(int xi=0;xi<nB;xi++)
    {
      ci = off/(Threads*inner);
      si = (off/inner)%Threads;
      i0 = off%inner;

      off1 = (si*civ + ci)*inner + i0 - z0;
      if(off1 >= 0)
      if((off1 < Threads*Biva) and (off1 < (Total - z0)) )
      {
        if( flag){res[(ci*sizeF+s0+si)*inner + i0] = buf[off1];}
        if(!flag){buf[off1] = src[(ci*sizeF+s0+si)*inner + i0];}
      }
      off += Threads;
    }
    __syncthreads();

    if(!flag){
    off = z0 + tid;
    for(int xi=0;xi<Biva;xi++)
    {
      if(off < Total){res[s0*civ*inner + off] = buf[xi*Threads + tid];off += Threads;}
    }
    __syncthreads();
    }

    z0 += Threads*Biva;
  }

}
#endif

////TODO change into Ty*
struct move_index
{
  //bool GPU;

  qlat::vector_gpu<char > buf;
  ////size_t buf_size;
  //qlat::vector<char* > pciv;

  //move_index(bool GPU_set=false){
  //  #ifndef QLAT_USE_ACC
  //  GPU = false;
  //  #else
  //  GPU = GPU_set;
  //  #endif
  //  buf = NULL;
  //  buf_size = 0;
  //}

  //void set_mem(int civ, size_t Bsize)
  //{
  //  TIMERA("move_index set_mem");
  //  if(buf_size != Bsize){
  //    free_mem();
  //    if(GPU){gpuMalloc(buf, Bsize, char);}
  //    else{buf = (void *)aligned_alloc_no_acc(Bsize);}
  //    buf_size = Bsize;
  //  }
  //  //////psrc.resize(civ);
  //}

  /////order follow src memory order
  template <typename Ty >
  void move_civ_out(Ty* src,Ty* res,int biva, long sizeF,int civ, int size_inner, bool GPU = false)
  {
    dojob(src, res, biva, civ, sizeF, 1, size_inner, GPU);
  }

  /////order follow src memory order
  template <typename Ty >
  void move_civ_in(Ty* src,Ty* res,int biva, int civ, long sizeF, int size_inner, bool GPU = false)
  {
    dojob(src, res, biva, civ, sizeF, 0, size_inner, GPU);
  }

  ////flag = 1 --> biva * sizeF * civ * size_inner --> biva * civ * sizeF * size_inner
  template <typename Ty >
  void dojob(Ty* src,Ty* res,int biva,int civ,long sizeF,int flag, int size_inner, bool GPU = false)
  {
  if(biva == 0 or civ == 0 or sizeF == 0 or size_inner == 0){return ;}
  /////size_t sizeF = sizeF0;

  ////size_t bufN = biva*civ*size_inner*sizeof(Ty)*sizeF;
  size_t Off = civ*sizeF*size_inner;
  #if PRINT_TIMER>5
  TIMER_FLOPS("reorder index");
  timer.flops += biva*Off*sizeof(Ty);
  #endif

  ////TIMERB("reorder index");
  if(size_inner < 1){qlat::displayln_info(qlat::ssprintf("size_inner too small %d !\n", size_inner));
    MPI_Barrier(get_comm());fflush(stdout);qassert(false);
  }

  if(src == res){buf.resize(Off*sizeof(Ty), GPU);}
  //pciv.resize(civ);
  Ty* s0;Ty *s1;
  //#ifdef QLAT_USE_ACC
  //if(GPU)
  if(src == res)if((Off*sizeof(Ty)) % sizeof(qlat::ComplexF) != 0){
    qlat::displayln_info(qlat::ssprintf("size not divided by 16, too small. \n"));qassert(false);}
  ///#endif

  for(int bi=0;bi<biva;bi++){
    s0 = &src[bi*Off];
    if(src == res){s1 = (Ty*)buf.data();}else{s1 = (Ty*) &res[bi*Off];}
    #ifdef QLAT_USE_ACC
    if(GPU){

      {
      const int Threads = 32;const int Biva =  (16*16+sizeof(Ty)-1)/sizeof(Ty);
      long Nb = (sizeF + Threads -1)/Threads;
      dim3 dimBlock(    Threads,    1, 1);
      dim3 dimGrid(     Nb,    1, 1);
      if(flag==0)move_index_global<Ty, false , Threads, Biva><<< dimGrid, dimBlock >>>(s0, s1, sizeF, civ, size_inner);
      if(flag==1)move_index_global<Ty, true  , Threads, Biva><<< dimGrid, dimBlock >>>(s0, s1, sizeF, civ, size_inner);
      qacc_barrier(dummy);
      }

      if(src == res){
      long Nvol = long(Off*sizeof(Ty)/sizeof(qlat::ComplexF));
      cpy_data_thread((qlat::ComplexF*) &res[bi*Off], (qlat::ComplexF*) s1, Nvol, 1);
      }

    continue ;}
    #endif

    #pragma omp parallel for
    for(long   si=0;si<sizeF;si++)
    for(int    ci=0;ci<civ;ci++)
    {
      Ty* p0=NULL;Ty* p1=NULL;
      if(flag == 1){
        p0 = (Ty*) &s0[(si*civ   + ci)*size_inner];
        p1 = (Ty*) &s1[(ci*sizeF + si)*size_inner];
      }
      if(flag == 0){
        p0 = (Ty*) &s0[(ci*sizeF + si)*size_inner];
        p1 = (Ty*) &s1[(si*civ   + ci)*size_inner];
      }
      memcpy(p1, p0, sizeof(Ty)*size_inner);
    }

    if(src == res){
      long Nvol = long(Off*sizeof(Ty)/sizeof(qlat::ComplexF));
      cpy_data_thread((qlat::ComplexF*) &res[bi*Off], (qlat::ComplexF*) s1, Nvol, 0);
    }

  }

  }

  void free_mem(){
    buf.resize(0);
  }

  ~move_index(){
    free_mem();
  }

};


inline void set_GPU(){
  #ifdef QLAT_USE_ACC
  int num_node;MPI_Comm_size(get_comm(), &num_node);
  int id_node;MPI_Comm_rank(get_comm(), &id_node);

  int num_gpus = 0;
  cudaGetDeviceCount(&num_gpus);
  ////cudaDeviceReset();
  cudaSetDevice(id_node % num_gpus);
  int gpu_id = -1; 
  cudaGetDevice(&gpu_id);
  //printf("CPU node %d (of %d) uses CUDA device %d\n", id_node, num_node, gpu_id);
  fflush(stdout);
  MPI_Barrier(get_comm());
  #endif

}

inline void set_GPU_threads(int mode=0){
  //////Set up gpu map to cpu
  (void)mode;
  #ifdef QLAT_USE_ACC
  int num_gpus = 0;
  cudaGetDeviceCount(&num_gpus);
  cudaDeviceReset();
  if(mode == 0){
  #pragma omp parallel
  {
    unsigned int cpu_thread_id = omp_get_thread_num();
    unsigned int num_cpu_threads = omp_get_num_threads();
    cudaSetDevice(cpu_thread_id % num_gpus);
    int gpu_id = -1; 
    cudaGetDevice(&gpu_id);
    printf("CPU thread %d (of %d) uses CUDA device %d\n", cpu_thread_id, num_cpu_threads, gpu_id);
  }}

  if(mode == 1){
  #pragma omp parallel
  {
    unsigned int cpu_thread_id = omp_get_thread_num();
    unsigned int num_cpu_threads = omp_get_num_threads();
    if(cpu_thread_id%num_gpus != 0){print0("Wrong mapping of omp !\n");qassert(false);}
    int Nthreads = cpu_thread_id/num_gpus;

    cudaSetDevice(cpu_thread_id / Nthreads);
    int gpu_id = -1; 
    printf("CPU thread %d (of %d) uses CUDA device %d\n", cpu_thread_id, num_cpu_threads, gpu_id);
    cudaGetDevice(&gpu_id);
  }}

  #endif

}

inline int init_mpi_thread(int* argc, char **argv[], int mode = 3)
{
  int provided;
  if(mode == 0)MPI_Init_thread(argc, argv, MPI_THREAD_SINGLE, &provided);
  if(mode == 1)MPI_Init_thread(argc, argv, MPI_THREAD_FUNNELED, &provided);
  if(mode == 2)MPI_Init_thread(argc, argv, MPI_THREAD_SERIALIZED, &provided);
  if(mode == 3)MPI_Init_thread(argc, argv, MPI_THREAD_MULTIPLE, &provided);

  int num_node;
  MPI_Comm_size(MPI_COMM_WORLD, &num_node);
  int id_node;
  MPI_Comm_rank(MPI_COMM_WORLD, &id_node);
  if (0 == id_node) {
    displayln("qlat::begin(): " +
              ssprintf("MPI Initialized. num_node = %d", num_node));
  }

  return num_node;
}

inline void begin_thread(
    int* argc, char** argv[],
    const std::vector<Coordinate>& size_node_list = std::vector<Coordinate>())
// begin Qlat and initialize a new comm
{
  const int num_node = init_mpi_thread(argc, argv);
  Coordinate size_node;
  for (int i = 0; i < (int)size_node_list.size(); ++i) {
    size_node = size_node_list[i];
    if (num_node == product(size_node)) {
      break;
    }   
  }
  if (num_node != product(size_node)) {
    size_node = plan_size_node(num_node);
  }
  begin_comm(get_comm(), size_node);
}

inline void print_mem_info(std::string stmp = "")
{
  print0("%s, ",stmp.c_str());
  #ifdef QLAT_USE_ACC
  double freeD = 0;double totalD=0;
  size_t freeM = 0;size_t totalM = 0;
  cudaMemGetInfo(&freeM,&totalM);
  freeD = freeM*pow(0.5,30);totalD = totalM*pow(0.5,30);
  #endif
  struct sysinfo s_info;
  sysinfo(&s_info);
  #ifdef QLAT_USE_ACC
  print0("===CPU free %.3e GB, total %.3e GB; GPU free %.3e GB, total %.3e GB. \n"
          , s_info.freeram*pow(0.5,30),s_info.totalram*pow(0.5,30),freeD,totalD);

  #else
  print0("===CPU free %.3e GB, total %.3e GB. \n"
          , s_info.freeram*pow(0.5,30),s_info.totalram*pow(0.5,30));
  #endif
}


inline int read_vector(const char *filename, std::vector<double > &dat)
{
  int prods = 0; 
  unsigned long Vsize = 0; 
  ////{synchronize();fflush(stdout);}

  if(qlat::get_id_node() == 0)
  {
    Vsize = get_file_size_o(filename);
    if(Vsize == 0){prods = 0;}else{prods=Vsize;}
    Vsize = Vsize/8;
  }
  sum_all_size(&prods,1);
  if(prods==0)return prods;

  sum_all_size((int*)&Vsize,1);
  dat.resize(Vsize);

  if(qlat::get_id_node() == 0)
  {
    FILE* filer = fopen(filename, "rb");
    unsigned long count = 1024*1024;
    unsigned long sizec = 0; 
    unsigned long offr  = 0; 
    for(unsigned int iv=0;iv<Vsize;iv++)
    {    
      if(offr >= Vsize*8)break;
      char* buf = (char *)&dat[offr/8];
      if((offr + count) <= (Vsize*8)){sizec = count;}
      else{sizec = Vsize*8 - offr;}

      fread(buf, 1, sizec, filer);
      offr = offr + sizec;
    }    

    fclose(filer);
  }

  sum_all_size(&dat[0],Vsize);
  return prods;

}

template<typename Yl>
void p_vector(const Yl teml)
{
  std::cout << teml << " ";
}

template<typename Ty>
void p_vector(const std::vector<Ty> teml)
{
  for(int i=0;i< teml.size();i++)
  {
    p_vector(teml[i]);
  }
  std::cout << std::endl;
}

template<typename Ty>
void p_vector(const qlat::vector<Ty> teml)
{
  for(int i=0;i< teml.size();i++)
  {
    p_vector(teml[i]);
  }
  std::cout << std::endl;
}

template<typename Ty>
inline void random_Ty(Ty* a, long N0,int GPU=0, int seed = 0, const int mode = 0)
{
  if(N0 == 0)return;
  qlat::RngState rs(qlat::get_id_node() + 1 + seed);

  double ini = qlat::u_rand_gen(rs);
  #ifdef QLAT_USE_ACC
  if(GPU == 1){
    size_t bfac = size_t(std::sqrt(N0));
    qacc_for(isp, long(N0/bfac + 1),{
     for(size_t i=0;i<size_t(bfac);i++){
      size_t off = isp*bfac + i;
      if(off < size_t(N0)){
        if(mode==0){a[off] = Ty(std::cos((ini+isp)*0.5) , (5.0/(isp+1))*ini*0.1);}
        if(mode==1){a[off] = Ty(std::cos((ini+isp)*0.5) , sin(5.0*(isp+1))*ini*0.1);}
      }
    }
    });
    return ;
  }
  #endif

  #pragma omp parallel for
  for(size_t isp=0;isp< size_t(N0);isp++)
  {
     if(mode==0){a[isp] = Ty(std::cos((ini+isp)*0.5) , (5.0/(isp+1))*ini*0.1);}
     if(mode==1){a[isp] = Ty(std::cos((ini+isp)*0.5) , sin(5.0*(isp+1))*ini*0.1);}
  }

}

template<typename Ty>
inline void random_EigenM(qlat::vector_acc<Ty >& a,int GPU=0, int seed = 0)
{
  Ty* buf = a.data();
  random_Ty(buf, a.size(), GPU, seed);
}

template<typename Ty>
inline void random_EigenM(std::vector<qlat::vector_acc<Ty > >& a, int GPU=0, int seed = 0)
{
  int N0 = a.size();if(N0 == 0)return ;
  for(size_t i=0;i < size_t(N0);i++){random_EigenM(a[i], GPU,  seed + i);}
}

template<typename Ty>
inline void zeroE(qlat::vector_acc<Ty >& a,int GPU=0, bool dummy=true)
{
  zero_Ty(a.data(), a.size(), GPU, dummy);
}

template<typename Ty>
inline void zeroE(std::vector<qlat::vector_acc<Ty > >& a,int GPU=0, bool dummy=true)
{
  for(LInt iv=0;iv<a.size();iv++){zeroE(a[iv], GPU, false);}
  if (dummy) {
    qacc_barrier(dummy);
  }
}

template<typename Ty, int civ>
inline void random_FieldM(qlat::FieldM<Ty , civ>& a,int GPU=0, int seed = 0)
{
  qassert(a.initialized);
  const Geometry& geo = a.geo();
  Ty* buf = (Ty*) qlat::get_data(a).data();
  random_Ty(buf, geo.local_volume() * civ, GPU, seed);
}

template<typename Ty, int civ>
inline Ty norm_FieldM(qlat::FieldM<Ty , civ>& a)
{
  qassert(a.initialized);
  const Geometry& geo = a.geo();
  Ty* buf = (Ty*) qlat::get_data(a).data();
  const long  V = geo.local_volume() ;
  qlat::vector_gpu<Ty > tmp;tmp.resize(V*civ);
  Ty* srcP = (Ty* ) qlat::get_data(a).data();
  Ty* resP = tmp.data();
  qacc_for(isp,V*civ,{ resP[isp] = srcP[isp];});
  Ty  norm = tmp.norm();
  return norm;
  //print0("norm %.3e %.3e \n", norm.real(), norm.imag());
}

template <class T>
void random_prop(Propagator4dT<T >& prop, int seed = -1)
{
  qassert(prop.initialized);
  ////print0("print time %.3f\n", tm.tv_sec);
  int rand_seed = qlat::get_id_node() + 1;
  if(seed == -1){timeval tm;gettimeofday(&tm, NULL);rand_seed += int(tm.tv_sec);}else{rand_seed += seed;}

  qlat::RngState rs(rand_seed);
  double ini = qlat::u_rand_gen(rs);

  /////int dir_limit = 4;
  const Geometry& geo = prop.geo();

  qacc_for(isp,  geo.local_volume(),{
    qlat::WilsonMatrixT<T >& v0 =  prop.get_elem(isp);
    //for(int ci=0;ci<12*12;ci++){
    //  //v0.p[ci] = T(std::cos((ini+isp + ci*2)*0.5 + (isp+ ci%4)/5) , ((5.0+ci)/(isp+1))*ini*0.1 + (isp*2 + ci%3)/5); 
    //  v0.p[ci] = T(std::cos((ini+isp + ci*2)*0.5 + ci) , ((5.0+ci)/(isp+1))*ini*0.1 + 0.2); 
    //}
    for(int ci=0;ci<12*12;ci++){
      v0.p[ci] = (ci/(12*12.0))*T(std::cos((ini+isp + ci*2)*0.5 + ci) , (ci+(5.0+ci)/(isp+1))*ini*0.1 + 0.2); 
    }
  }); 
}

template <class T>
void random_link(GaugeFieldT<T> &gf, const int seed = -1)
{
  if(seed == -1)
  {
    qacc_for(isp, gf.field.size(), { set_unit(gf.get_elem(isp), 1.0);});
  }else{
    //T* res = (T*) gf.get_elem(0).p;
    const Geometry& geo = gf.geo();
    T* res = (T*) qlat::get_data(gf).data();
    random_Ty(res, geo.local_volume()*geo.multiplicity*sizeof(ColorMatrixT<T>)/sizeof(T), 1, seed);

    //qacc_for(isp, gf.field.size(), { set_unit(gf.get_elem(isp), 1.0);});
    ColorMatrixT<T> unit;set_unit(unit, 1.0);
    /////TODO This function cannot be done on GPU
    /////Eigen normalize/normalized problem
    for(long isp=0;isp<gf.field.size();isp++)
    {
      gf.get_elem(isp) = gf.get_elem(isp) * (1/2.0) + unit;
      unitarize(gf.get_elem(isp));
    }
  }
}


inline size_t get_threads(size_t thread, size_t Total, int def=1)
{
  for(size_t temb=thread;temb<Total;temb++)
  {
    if(Total%temb == 0)
    {
      return temb;
    }
  }
  return def;
}

inline std::vector<unsigned int > get_num_power(const size_t x,const std::vector<unsigned int >& a)
{
  std::vector<unsigned int > re; re.resize(a.size());
  for(unsigned int ai=0;ai<a.size();ai++)
  {
    unsigned int    fac   = a[ai];
    unsigned int    count = 0;
    size_t use   = x;
    for(unsigned int j = 0;j < 300; j++)
    {
      if(use % fac == 0){
        count += 1;
        use = use/fac;}
      else{
        re[ai] = count;
        break;
      }
    }
  }
  size_t test = 1;
  for(unsigned int ai=0;ai<a.size();ai++){
    test = test * std::pow(a[ai], re[ai]);
  }
  assert(test == x);
  return re;
}

//////Even node in xyzT directions
inline Coordinate spread_even(const int n, const Coordinate& Lat, const std::vector<unsigned int >& a)
{
  std::vector<unsigned int > Mpow = get_num_power(n, a);
  std::vector<std::vector<unsigned int > > Lpow;
  Lpow.resize(4);for(int i=0;i<4;i++){Lpow[i] =  get_num_power(Lat[i], a);}
  Coordinate re;
  for(int i=0;i<4;i++)re[i] = 1;

  std::vector<int > nL(4);

  for(LInt i=0;i<Mpow.size();i++){
    int fac = a[i];
    int num = Mpow[i];
    if(Mpow[i] != 0){
      int suma = Lpow[0][i] + Lpow[1][i] + Lpow[2][i] + Lpow[3][i] ;
      assert(int(Mpow[i]) <= suma);

      for(int j=0;j<4;j++){nL[j] = Lpow[j][i];}
      //std::vector<int>::iterator result = std::max_element(nL.begin(), nL.end());
      //int pos_max = std::distance(nL.begin(), result);
      for(unsigned int ni = 0; ni< Mpow[i] + 1 ;ni++)
      {
        for(int j=3;j>=0;j--)
        {
          if(nL[j] > 0){num -= 1; nL[j] -= 1; re[j] *= fac;}
          if(num == 0){break;}
        }
        if(num == 0){break;}
      }
    }
  }

  assert(re[0]*re[1]*re[2]*re[3] == n);
  return re;
}


//////Most power in T direction
inline Coordinate spread_powT(const int n, const Coordinate& Lat, const std::vector<unsigned int >& a)
{
  std::vector<unsigned int > Mpow = get_num_power(n, a);
  std::vector<std::vector<unsigned int > > Lpow;
  Lpow.resize(4);for(int i=0;i<4;i++){Lpow[i] =  get_num_power(Lat[i], a);}
  Coordinate re;
  for(int i=0;i<4;i++)re[i] = 1;

  for(LInt i=0;i<Mpow.size();i++){
    int num = Mpow[i];
    if(num != 0){
      int suma = Lpow[0][i] + Lpow[1][i] + Lpow[2][i] + Lpow[3][i] ;
      assert(num <= suma);
      unsigned int tem = num;
      for(unsigned int ni=0;ni<4;ni++){
        if(tem >= Lpow[4-ni-1][i]){
          re[4-ni-1] *= int(std::pow(a[i], Lpow[4-ni-1][i]));
          tem = tem - Lpow[4-ni-1][i];
        }
        else{
          re[4-ni-1] *= (unsigned int)(std::pow(a[i], tem));
          tem = 0;
        }
        ////if(tem >= Lpow[4-ni-1][i]){}
        if(tem == 0)break;
      }
    }
  }

  assert(re[0]*re[1]*re[2]*re[3] == n);
  return re;

}

inline Coordinate guess_nodeL(int n, const Coordinate& Lat, const int mode = 0)
{

  std::vector<unsigned int > a;a.resize(8);
  a[0] = 2;a[1] = 3;a[2] = 5;a[3] = 7;a[4] =11;a[5] =13;a[6] =17;a[7] =19;
  Coordinate re;
  if(mode == 0){re = spread_powT(n, Lat,a);}
  if(mode == 1){re = spread_even(n, Lat,a);}
  return re;
}

inline void add_nodeL(std::vector<Coordinate>& size_node_list)
{
  size_node_list.push_back(Coordinate(1, 1, 1,  1));
  size_node_list.push_back(Coordinate(1, 2, 1,  1));
  //size_node_list.push_back(Coordinate(1, 1, 1,  2));
  //size_node_list.push_back(Coordinate(1, 2, 2,  1));
  size_node_list.push_back(Coordinate(1, 2, 1,  2));
  size_node_list.push_back(Coordinate(2, 2, 2,  1));
  size_node_list.push_back(Coordinate(2, 2, 2,  2));
  size_node_list.push_back(Coordinate(1, 1, 3,  1));
  ////size_node_list.push_back(Coordinate(1, 1, 2,  2));
  size_node_list.push_back(Coordinate(1, 1, 1,  4));
  size_node_list.push_back(Coordinate(1, 1, 3,  2));
  size_node_list.push_back(Coordinate(1, 2, 3,  1));
  size_node_list.push_back(Coordinate(1, 2, 4,  1));
  size_node_list.push_back(Coordinate(1, 1, 1, 12));
  size_node_list.push_back(Coordinate(1, 1, 1, 16));
  size_node_list.push_back(Coordinate(1, 1, 6,  4));
  size_node_list.push_back(Coordinate(1, 1, 1, 32));
  size_node_list.push_back(Coordinate(2, 4, 4, 2 ));

  size_node_list.push_back(Coordinate(4, 4, 8, 16));
  size_node_list.push_back(Coordinate(4, 8, 8, 16));
  size_node_list.push_back(Coordinate(4, 4, 4,  6));
  size_node_list.push_back(Coordinate(4, 4, 4, 12));
  size_node_list.push_back(Coordinate(4, 4, 8, 12));
  size_node_list.push_back(Coordinate(4, 8, 8, 16));
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


////mode_dis % 1 == 0, t in diff node, mode_dis % 2 == 1, t in single node
////mode_dis < 2, T, mode_dis >= 2 even
inline void begin_Lat(int* argc, char** argv[], inputpara& in, int read_Lat = 0){
  if(read_Lat >= 0)
  {
    int n_node = init_mpi(argc, argv);
    in.load_para(*argc, *argv);
    Coordinate Lat(in.nx, in.ny, in.nz, in.nt);
    Coordinate spreadT;
    if(in.layout != std::string("NONE")){
      spreadT = string_to_Coordinate(in.layout);
      if(spreadT[0] * spreadT[1] * spreadT[2] * spreadT[3] != n_node)
      {printf("Wrong input layout\n");abort_r();}
    }
    else{
      if(in.mode_dis >= 0 and in.mode_dis < 2){spreadT = guess_nodeL(n_node, Lat, 0);}
      if(in.mode_dis >= 2 and in.mode_dis < 4){spreadT = guess_nodeL(n_node, Lat, 1);}
    }
    ///3D begin
    ////begin_comm(MPI_COMM_WORLD , spreadT);

    ///4D begin
    int id_node, n;
    MPI_Comm_size(MPI_COMM_WORLD, &n);
    MPI_Comm_rank(MPI_COMM_WORLD, &id_node);
    int t =  id_node/(spreadT[0]*spreadT[1]*spreadT[2]);
    int z = (id_node/(spreadT[0]*spreadT[1]))%(spreadT[2]);
    int y = (id_node/(spreadT[0]))%(spreadT[1]);
    int x = (id_node%(spreadT[0]));
    ///int new_id = ((z*spreadT[1] + y)*spreadT[0] + x)*spreadT[3] + t;
    int new_id = ((x*spreadT[1] + y)*spreadT[2] + z)*spreadT[3] + t;
    if(in.mode_dis % 2 == 0)begin(id_node, spreadT);
    if(in.mode_dis % 2 == 1)begin(new_id, spreadT);
  }

  if(read_Lat == -1)
  {
    std::vector<Coordinate> size_node_list;
    add_nodeL(size_node_list);
    begin(argc, argv, size_node_list);
    in.load_para(*argc, *argv);
  }

  set_GPU();

  omp_set_num_threads(omp_get_max_threads());
  print0("===nthreads %8d %8d, max %8d \n",qlat::qacc_num_threads(),omp_get_num_threads(),omp_get_max_threads());

  fflush_MPI();
  print_mem_info();

}

inline std::vector<long > job_create(long total, long each)
{
  if(total < 1 or each < 1){
    print0("Give me valid job types total %ld, each %ld \n", total, each);
    abort_r();}
  /////std::vector<long > a = job_create(total, each);
  std::vector<long > a;a.resize(0);
  long jobN  = (total + each - 1)/each;
  int i0 = 0; int dj = each;
  for(int ji = 0; ji < jobN ; ji++)
  {
    if(i0 >= total){break;}
    if(i0 + dj > total){dj = total - i0;}
    a.push_back(i0);
    a.push_back(dj);
    i0 += dj;
  }

  return a;
}

template<typename Ty>
inline void allocate_buf(std::vector<qlat::vector_gpu<Ty > > & buf, size_t n0, size_t n1)
{
  TIMERA("CUDA Buf mem allocation");
  buf.resize(n0);
  for(LInt i=0;i<buf.size();i++){
    buf[i].resize(n1);
  }
}

template<typename Ty>
inline void allocate_buf(std::vector<qlat::vector_acc<Ty > > & buf, size_t n0, size_t n1)
{
  TIMERA("CUDA Buf mem allocation");
  buf.resize(0);
  buf.resize(n0);
  for(LInt i=0;i<buf.size();i++){
    buf[i].resize(0);
    buf[i].resize(n1);
  }
}

template<typename Ty>
inline void allocate_buf(std::vector<Ty* >& buf, size_t n0, size_t n1)
{
  TIMERA("CUDA Buf mem allocation");
  buf.resize(n0);
  for(LInt i=0;i<buf.size();i++){
    gpuMalloc(buf[i], n1, Ty);
  }
}

template<typename Ty>
inline Ty inv_self(const Ty& lam, double m, double rho,int one_minus_halfD=1)
{
  std::complex<double > tem(lam.real(),lam.imag());
  std::complex<double > v0 = (one_minus_halfD>0)?(1.0-tem/2.0)/(rho*tem+m*(1.0-tem/2.0)):1.0/(rho*tem+m*(1.0-tem/2.0));
  Ty res(v0.real(),v0.imag());
  return res;
}

template<typename Ty>
qlat::vector_acc<Ty* > EigenM_to_pointers(std::vector<qlat::vector_gpu<Ty > >& src)
{
  qlat::vector_acc< Ty* >  res;
  res.resize(src.size());
  for(LInt iv=0;iv<src.size();iv++)
  {
    res[iv] = src[iv].data();
  }
  return res;
}

template<typename Ty>
qlat::vector_acc<Ty* > EigenM_to_pointers(std::vector<qlat::vector_acc<Ty > >& src)
{
  qlat::vector_acc<Ty* >  res;
  res.resize(src.size());
  for(LInt iv=0;iv<src.size();iv++)
  {
    res[iv] = src[iv].data();
  }
  return res;
}

/////Ty should have size(), resize(), and data()
template<class Ty>
Ty sum_local_to_global_vector(Ty src, MPI_Comm* commp=NULL)
{
  ////int Nt = geo.node_site[3];
  int Nmpi  = qlat::get_num_node();
  int rank  = qlat::get_id_node();
  if(commp != NULL){MPI_Comm_size(*commp, &Nmpi);MPI_Comm_rank(*commp, &rank);}

  std::vector<long > size_global;size_global.resize(Nmpi);for(unsigned long i=0;i<size_global.size();i++){size_global[i]=0;}

  size_global[rank] = src.size();

  sum_all_size(size_global.data(), size_global.size(), 0, commp);

  long total = 0;long current = 0;
  for(int i=0;i<Nmpi;i++){total += size_global[i];if(i < rank){current += size_global[i];}}

  /////for(unsigned int i=0;i<size_global.size();i++){printf("rank %d, size %ld \n", rank, size_global[i]);}
  /////printf("rank %d, Total %ld, current %ld \n", rank, total, current);

  Ty res;res.resize(total);
  for(unsigned long pos=current;pos<src.size();pos++){res[pos] = src[pos - current];}

  sum_all_size(res.data(), res.size(), 0, commp);

  return res;

}

inline std::vector<int > num_to_site(const long num, const std::vector<int > key_T)
{
  qassert(key_T.size() > 0);
  int dim = key_T.size();
  std::vector<int > site;site.resize(dim);
  for(int iv=0;iv<dim;iv++){site[iv] = 0;}

  long tem_i = num;
  for(int Ni=0; Ni < dim; Ni++)
  {
    long N_T = 1;
    for(int numi = Ni+1;numi < dim;numi++)
    {
      N_T = N_T*key_T[numi];
    }
    site[Ni] = tem_i/N_T;
    tem_i = tem_i%N_T;
  }
  return site;
}

inline std::vector<long > random_list(const long n, const long m, const int seed)
{
  std::vector<long > r;
  std::vector<long > b;
  qlat::RngState rs(13021 + seed);
  for(long ri = 0; ri < n; ri++){r.push_back(ri);}
  if(m >= n or m <  0){return r;}
  if(m == 0){r.resize(0);return r;}

  for(long ri=0;ri<m;ri++){
    long u = int(qlat::u_rand_gen(rs) * r.size());
    b.push_back(r[u]);
    r.erase(r.begin()+ u);
  }
  return b;

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
