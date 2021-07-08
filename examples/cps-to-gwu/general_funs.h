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

#include "float_type.h"
#include<type_traits>

#include <iterator>
#include <sys/sysinfo.h>

////#include <cuda_runtime.h>
//
////#include <Eigen/Dense>
////#include "fftw3.h"
////#include "fftw3-mpi.h"

namespace qlat
{

#ifndef QLAT_USE_ACC
#ifdef __GNUC__
///////Multiply of different types of complex

//inline std::complex<double> operator*(const std::complex<float> &a, const double &b) {
//    return std::complex<double>(a.real()*b, a.imag()*b);
//}
//
//inline std::complex<double> operator/(const std::complex<float> &a, const double &b) {
//    return std::complex<double>(a.real()/b, a.imag()/b);
//}

inline std::complex<double> operator*(const std::complex<float> &a, const std::complex<double > &b) {
    return std::complex<double>(a.real() * b.real() - a.imag()*b.imag(), a.imag()*b.real() + a.real()*b.imag());
}

inline std::complex<double> operator*(const std::complex<double > &b, const std::complex<float> &a) {
    return std::complex<double>(a.real() * b.real() - a.imag()*b.imag(), a.imag()*b.real() + a.real()*b.imag());
}
#endif
#endif


inline unsigned int get_node_rank_funs0()
{
  int rank;
  //MPI_Comm_rank(get_comm(), &rank);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
}

#define print0 if(qlat::get_id_node() == 0) printf

inline unsigned long get_file_size_o(const char *filename)
{
  std::ifstream File(filename);
  if(!File.is_open()){print0("file is not exist\n");return 0;}
  unsigned long Begin = File.tellg();
  File.seekg(0, std::ios_base::end);
  unsigned long End = File.tellg();
  File.close();
  return End-Begin;
}

inline size_t get_file_size_MPI(const char *filename)
{
  size_t sizen = 0;
  if(qlat::get_id_node()==0){
    std::ifstream File(filename);
    if(!File.is_open()){print0("file is not exist\n");sizen = 0;}
    else{
      unsigned long Begin = File.tellg();
      File.seekg(0, std::ios_base::end);
      unsigned long End = File.tellg();
      File.close();
      sizen = End-Begin;
    }
  }
  MPI_Bcast(&sizen, sizeof(size_t), MPI_CHAR, 0, get_comm());
  return sizen;
}

template<typename Ty>
void sum_all_size_threads(Ty *src,Ty *sav,long size)
{
  ////if(qlat::get_num_node()==1)return;
  Ty *res;
  if(src == sav){
    res=new Ty[size];
  }else{res = sav;}
  if(qlat::get_num_node() == 1){
    if(src == sav){return;}
    if(src != sav){memcpy(sav,src,size*sizeof(Ty));return;}
  }

  int provided;int status = MPI_Query_thread(&provided);
  if(provided == MPI_THREAD_SINGLE)
  {
  if(std::is_same<Ty, unsigned long>::value)MPI_Allreduce(src,res, size, MPI_UNSIGNED_LONG, MPI_SUM, get_comm());
  if(std::is_same<Ty, int>::value)MPI_Allreduce(src,res, size, MPI_INT   , MPI_SUM, get_comm());
  if(std::is_same<Ty, float>::value)MPI_Allreduce(src,res, size, MPI_FLOAT , MPI_SUM, get_comm());
  if(std::is_same<Ty, double>::value)MPI_Allreduce(src,res, size, MPI_DOUBLE, MPI_SUM, get_comm());
  }

  if(provided == MPI_THREAD_MULTIPLE)
  {
    int Nv = omp_get_num_threads();
    int ny = (size+Nv-1)/Nv;
    #pragma omp parallel
    for(int iv=0;iv<Nv;iv++){
      int sizeI = ny;if(iv==Nv-1){sizeI = size - iv*ny;}
      int off0  = iv*ny;
      Ty *src0 = &src[off0];
      Ty *res0 = &res[off0];
    if(std::is_same<Ty, unsigned long>::value)MPI_Allreduce(src0,res0, size, MPI_UNSIGNED_LONG, MPI_SUM, get_comm());
    if(std::is_same<Ty, int>::value)MPI_Allreduce(src0,res0, size, MPI_INT   , MPI_SUM, get_comm());
    if(std::is_same<Ty, float>::value)MPI_Allreduce(src0,res0, size, MPI_FLOAT , MPI_SUM, get_comm());
    if(std::is_same<Ty, double>::value)MPI_Allreduce(src0,res0, size, MPI_DOUBLE, MPI_SUM, get_comm());
    }
  }

  if(src == sav){
    memcpy(sav,res,size*sizeof(Ty));
    delete []res;
  }
}

template<typename Ty>
void get_MPI_type(Ty& a, MPI_Datatype& curr, unsigned int& size, int mode = 0)
{
  curr = MPI_BYTE; size = 1;
  bool is_same_v = false;

  is_same_v = std::is_same< Ty , int>::value;
  if(is_same_v){curr = MPI_INT;size=sizeof(int);return ;}

  is_same_v = std::is_same< Ty , unsigned int>::value;
  if(is_same_v){curr = MPI_UNSIGNED;size=sizeof(unsigned int);return ;}

  is_same_v = std::is_same< Ty , long>::value;
  if(is_same_v){curr = MPI_LONG;size=sizeof(long);return ;}

  is_same_v = std::is_same< Ty , unsigned long>::value;
  if(is_same_v){curr = MPI_UNSIGNED_LONG;size=sizeof(unsigned long);return ;}

  is_same_v = std::is_same< Ty , float>::value;
  if(is_same_v){curr = MPI_FLOAT;size=sizeof(float);return ;}

  is_same_v = std::is_same< Ty , double>::value;
  if(is_same_v){curr = MPI_DOUBLE;size=sizeof(double);return ;}

  is_same_v = std::is_same< Ty , std::int8_t>::value;
  if(is_same_v){curr = MPI_INT8_T;size=sizeof(std::int8_t);return ;}

  is_same_v = std::is_same< Ty , std::int16_t>::value;
  if(is_same_v){curr = MPI_INT16_T;size=sizeof(std::int16_t);return ;}

  is_same_v = std::is_same< Ty , std::int32_t>::value;
  if(is_same_v){curr = MPI_INT32_T;size=sizeof(std::int32_t);return ;}

  is_same_v = std::is_same< Ty , std::int64_t>::value;
  if(is_same_v){curr = MPI_INT64_T;size=sizeof(std::int64_t);return ;}


  if(mode == 1){print0("Type not found !!!! \n");qassert(false);}

  /////Complex types sum and equal
  if(sizeof(Ty)== 8){curr = MPI_FLOAT  ;size = sizeof(float) ;return ;}
  if(sizeof(Ty)==16){curr = MPI_DOUBLE ;size = sizeof(double) ;return ;}

}


template<typename Ty>
void sum_all_size(Ty *src,Ty *sav,long size)
{
  Ty *res;
  if(src == sav){res = new Ty[size];}else{res = sav;}
  if(qlat::get_num_node() == 1){
    if(src == sav){return;}
    if(src != sav){memcpy(sav,src,size*sizeof(Ty));return;}
  }

  //if(std::is_same<Ty, unsigned long>::value)MPI_Allreduce(src,res, size, MPI_UNSIGNED_LONG, MPI_SUM, get_comm());
  //if(std::is_same<Ty, int>::value)MPI_Allreduce(src,res, size, MPI_INT   , MPI_SUM, get_comm());
  //if(std::is_same<Ty, float>::value)MPI_Allreduce(src,res, size, MPI_FLOAT , MPI_SUM, get_comm());
  //if(std::is_same<Ty, double>::value)MPI_Allreduce(src,res, size, MPI_DOUBLE, MPI_SUM, get_comm());

  MPI_Datatype curr = MPI_DOUBLE;unsigned int M_size = sizeof(double);Ty atem=0;
  get_MPI_type(atem, curr, M_size, 1);

  MPI_Allreduce(src,res, size, curr, MPI_SUM, get_comm());

  if(src == sav)
  {
    memcpy(sav,res,size*sizeof(Ty));
    delete []res;
  }
}

template<typename Ty>
void sum_all_size(Ty *src,long size)
{
  sum_all_size(src,src,size);
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

template<typename Ty>
void Redistribute_all_Nt(Ty *src,long size,const qlat::Geometry &geo)
{
  if(qlat::get_num_node() == 1){return;}
  int Nt = geo.node_site[3];
  int Nmpi  = qlat::get_num_node();

  const Coordinate vg = geo.total_site();
  const int nt = vg[3];

  int mt = nt/Nt;
  if(mt != Nmpi){print0("Note supported !");qassert(false);return;}

  int rank  = qlat::get_id_node();
  long size_c = sizeof(Ty)*size/mt;

  //unsigned short t0 = 0;
  //{
  //  Coordinate xl = geo.coordinate_from_index(0);
  //  xl[3] = 0;
  //  Coordinate xg = geo.coordinate_g_from_l(xl);
  //  t0 = xg[3];
  //}

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

  std::vector<Ty >buf;buf.resize(size);
  {
  ////TIMER("MPI call CPU");
  MPI_Alltoallv(src,(int*) &send[0],(int*) &spls[0], MPI_CHAR,
            &buf[0],(int*) &recv[0],(int*) &rpls[0], MPI_CHAR, get_comm());
  }
  #pragma omp parallel for
  for(long isp=0;isp<size;isp++){src[isp] = buf[isp];}

  /////memcpy(src,&buf[0],size*sizeof(Ty));

  //qlat::vector<Ty >buf;buf.resize(size);
  //{
  //  Ty *res;res = &buf[0];
  //  TIMER("MPI call CPU");
  //  int* sendp = (int*) &send[0];
  //  int* splsp = (int*) &spls[0];
  //  int* recvp = (int*) &recv[0];
  //  int* rplsp = (int*) &rpls[0];
  //  #pragma acc host_data use_device (src, res,sendp,splsp,recvp,rplsp)
  //  MPI_Alltoallv(src,sendp,recvp, MPI_CHAR,
  //                res,recvp,rplsp, MPI_CHAR, get_comm());
  //}
  //qacc_for(isp, size,{ src[isp] = buf[isp];});

  ////cudaMemcpy(src,&buf[0], size*sizeof(Ty), cudaMemcpyDeviceToDevice);

  ////memcpy(src,&buf[0],size*sizeof(Ty));

  //Ty *buf0;cudaMalloc((void **)&buf0, size*sizeof(Ty));
  ////cudaMemcpy(buf0,src, size*sizeof(Ty), cudaMemcpyHostToDevice);
  //{TIMER("MPI call GPU");
  //MPI_Alltoallv(src,(int*) &send[0],(int*) &spls[0], MPI_CHAR,
  //          &buf0[0],(int*) &recv[0],(int*) &rpls[0], MPI_CHAR, get_comm());}
  ////cudaMemcpy(src, buf1, size*sizeof(Ty), cudaMemcpyDeviceToHost);
  //cudaMemcpy(src, buf0, size*sizeof(Ty), cudaMemcpyDeviceToDevice);
  //cudaFree(buf0);



  //Ty *buf0;cudaMalloc((void **)&buf0, size*sizeof(Ty));
  ////cudaMemcpy(buf0,src, size*sizeof(Ty), cudaMemcpyHostToDevice);
  //{TIMER("MPI call GPU");
  //MPI_Alltoallv(src,(int*) &send[0],(int*) &spls[0], MPI_CHAR,
  //          &buf0[0],(int*) &recv[0],(int*) &rpls[0], MPI_CHAR, get_comm());}
  ////cudaMemcpy(src, buf1, size*sizeof(Ty), cudaMemcpyDeviceToHost);
  //cudaMemcpy(src, buf0, size*sizeof(Ty), cudaMemcpyDeviceToDevice);
  //cudaFree(buf0);


  //Ty *buf0;cudaMalloc((void **)&buf0, size*sizeof(Ty));
  //Ty *buf1;cudaMalloc((void **)&buf1, size*sizeof(Ty));
  ////cudaMemcpy(buf0,src, size*sizeof(Ty), cudaMemcpyHostToDevice);
  //cudaMemcpy(buf0,src, size*sizeof(Ty), cudaMemcpyDeviceToDevice);
  //{
  //TIMER("MPI call GPU");
  //int* sendp = (int*) &send[0];
  //int* splsp = (int*) &spls[0];
  //int* recvp = (int*) &recv[0];
  //int* rplsp = (int*) &rpls[0];

  //#pragma acc host_data use_device (src, res,sendp,splsp,recvp,rplsp)
  //MPI_Alltoallv(&buf0[0],(int*) &sendp[0],(int*) &splsp[0], MPI_CHAR,
  //          &buf1[0],(int*) &recvp[0],(int*) &rplsp[0], MPI_CHAR, get_comm());
  //}
  //cudaMemcpy(src, buf1, size*sizeof(Ty), cudaMemcpyDeviceToDevice);
  ////cudaMemcpy(src, buf1, size*sizeof(Ty), cudaMemcpyHostToHost);
  //cudaFree(buf0);
  //cudaFree(buf1);

  //cudaError_t status = cudaMallocHost((void**)&h_aPinned, bytes);
  //if (status != cudaSuccess)printf("===Error allocating pinned host memory\n");
  ////cudaMemcpy(val_host, val_device, sizeof(float), cudaMemcpyDeviceToHost);
  ////cudaMemcpy(val_host, val_device, sizeof(float), cudaMemcpyHostToDevice);

}

inline void abort_r(std::string stmp)
{
  print0("%s\n",stmp.c_str());
  MPI_Barrier(get_comm());
  fflush(stdout);
  ////MPI_Finalize();
  qlat::end();
  abort();
}

inline void abort_sum(double flag)
{
  sum_all_size(&flag,1);
  if(flag > 0)
  {
    abort_r("");
    //MPI_Barrier(get_comm());
    //fflush(stdout);
    //MPI_Finalize();
    //abort();
  }
}

inline void fflush_MPI(){
  MPI_Barrier(get_comm());
  fflush(stdout);
}

inline void reorder_civ(char* src,char* res,int biva,int civ,int sizeF,int flag,int size_inner)
{
  //TIMER("reorder_civ vectors char");
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
    for(int si=0;si<sizeF;si++)
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

////////Assuming M must be field
//template<typename M,int DIMN>
//void touchv(qlat::FieldM<M, DIMN> &a0){
//  TIMER("Touch GPU mem");
//  const qlat::Geometry &geo = a0.geo();
//  qacc_for(index, geo.local_volume(), {
//    for(int di=0;di<DIMN;di++)
//    {
//      a0.get_elems(index)[di] = a0.get_elems(index)[di]*2;
//      a0.get_elems(index)[di] = a0.get_elems(index)[di]/2;
//    }
//  });
//}
//template<class M,int DIMN >
//void touchM(qlat::FieldM<M, DIMN> &a0){
//  TIMER("Touch GPU mem");
//  const qlat::Geometry &geo = a0.geo();
//  int an = a0.get_elems(0)[0].em().cols();
//  int bn = a0.get_elems(0)[0].em().rows();
//  qacc_for(index, geo.local_volume(), {
//    for(int di=0;di<DIMN;di++)
//    for(int ai=0;ai<an;ai++)
//    for(int bi=0;bi<bn;bi++)
//    {
//      a0.get_elems(index)[di].em().col(ai)[bi] = a0.get_elems(index)[di].em().col(ai)[bi]*2;
//      a0.get_elems(index)[di].em().col(ai)[bi] = a0.get_elems(index)[di].em().col(ai)[bi]/2;
//    }
//  });
//}


inline void set_GPU(){
  //////Set up gpu map to cpu
  #ifdef QLAT_USE_ACC
  //int rank, local_rank, local_size;
  //MPI_Comm local_comm;
  //MPI_Comm_rank(get_comm(), &rank);
  ////MPI_Comm_split_type(get_comm(), MPI_COMM_TYPE_SHARED, rank,  MPI_INFO_NULL, &local_comm);
  ////MPI_Comm_size(local_comm, &local_size);
  ////MPI_Comm_rank(local_comm, &local_rank);
  ////cudaSetDevice(local_rank%local_size);
  //cudaSetDevice(local_rank%local_size);

  int num_node;MPI_Comm_size(get_comm(), &num_node);
  int id_node;MPI_Comm_rank(get_comm(), &id_node);

  int num_gpus = 0;
  cudaGetDeviceCount(&num_gpus);
  ////cudaDeviceReset();
  cudaSetDevice(id_node % num_gpus);
  int gpu_id = -1; 
  cudaGetDevice(&gpu_id);
  printf("CPU node %d (of %d) uses CUDA device %d\n", id_node, num_node, gpu_id);
  fflush(stdout);
  MPI_Barrier(get_comm());
  #endif

}

inline void set_GPU_threads(int mode=0){
  //////Set up gpu map to cpu
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

inline int init_mpi_thread(int* argc, char **argv[])
{
  int provided;
  MPI_Init_thread(argc, argv, MPI_THREAD_MULTIPLE, &provided);

  int num_node;
  MPI_Comm_size(get_comm(), &num_node);
  int id_node;
  MPI_Comm_rank(get_comm(), &id_node);
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
  ////double length = (geo.local_volume()*pow(0.5,30))*12*sizeof(Complexq);
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


inline std::vector<std::string > stringtolist(std::string &tem_string)
{
  std::istringstream iss(tem_string);
  std::vector<std::string> results((std::istream_iterator<std::string>(iss)),std::istream_iterator<std::string>());
  return results;
}

inline double stringtodouble(std::string &tem_string)
{
  //double use = atof(tem_string.c_str());
  double use = 0.0;
  if(tem_string!="_NONE_")
  {
    use = atof(tem_string.c_str());
  }
  return use;
}

inline int stringtonum(std::string &tem_string)
{
  //int t_Total = 0;
  //if(tem_string!="_NONE_")
  //{
  //  int tem_length = strlen(tem_string.c_str()); 
  //  for(int i=0;i<tem_length;i++){t_Total = t_Total+(tem_string.c_str()[i]-'0')*std::pow(10,tem_length-i-1);};
  //}
  //return t_Total;
  double use = 0.0;
  if(tem_string!="_NONE_")
  {
    use = atof(tem_string.c_str());
  }
  return int(use);

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
    //std::vector<double > tem;
    //tem.resize(3);
    unsigned long count = 1024*1024;
    unsigned long sizec = 0; 
    unsigned long offr  = 0; 

    //char* buf = (char *)&dat[0];
    //fread(buf, 1, Vsize*8, filer);

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


inline void read_input(const char *filename,std::vector<std::vector<std::string > > &read_f)
{
  read_f.resize(0);
  if(get_file_size_o(filename) == 0){return;}
  FILE* filer = fopen(filename, "r");
  char sTemp[300],tem[300];
  ///std::string s0(sTemp);

  int count_line = 0;
  while(!feof(filer)){
    fgets(tem, 300, filer);
    if(std::string(tem).size() >= 2){
      sprintf(sTemp,"%s",tem);
      std::string s0(sTemp);
      std::vector<std::string > resv = stringtolist(s0);
      read_f.push_back(resv);
      count_line += 1;
    }
  }
  fclose(filer);
}

////Bcast conf_l from zero rank
inline void bcast_vstring(std::vector<std::string> &conf_l, const int Host_rank = 0){

  int rank = get_node_rank_funs0();
  ////Bcast strings
  size_t sizen = 0;if(rank == Host_rank)sizen = conf_l.size();
  MPI_Bcast(&sizen, sizeof(size_t), MPI_CHAR, Host_rank, get_comm());
  if(rank != Host_rank)conf_l.resize(sizen);
  for(unsigned int is=0;is<conf_l.size();is++){
    if(rank == Host_rank)sizen = conf_l[is].size();
    MPI_Bcast(&sizen, sizeof(size_t), MPI_CHAR, Host_rank, get_comm());

    if(rank != Host_rank)conf_l[is].resize(sizen);
    MPI_Bcast(&conf_l[is][0], sizen, MPI_CHAR, Host_rank, get_comm());
  }
  ////Bcast strings

}

struct inputpara{
  std::vector<std::vector<std::string > > read_f;
  int bSize;
  int bSum;
  int cutN;
  std::string lat;
  int icfg;
  int save_prop;

  int bini, ncut0, ncut1;

  int nx;
  int ny;
  int nz;
  int nt;
  int nvec;
  int ionum;
  int nmass;
  std::string Ename;
  std::string Pname;
  std::string paraI;

  int debuga;

  ~inputpara(){
    for(unsigned int is=0;is<read_f.size();is++){
     for(unsigned int ic=0;ic<read_f[is].size();ic++){read_f[is][ic].resize(0);}
     read_f[is].resize(0);
    }
    read_f.resize(0);
  }

  int find_para(const std::string &str2, int &res){
    for(unsigned int is=0;is<read_f.size();is++){
      ////std::string str2("bSize");
      std::size_t found = read_f[is][0].find(str2);
      if(found != std::string::npos and read_f[is].size() >= 2){
        res = stringtonum(read_f[is][1]);
        if(get_node_rank_funs0() == 0)printf("  %20s %10d \n",str2.c_str(), res);
        return 1;
      }
    }
    return 0;
  }

  int find_para(const std::string &str2, std::string &res){
    for(unsigned int is=0;is<read_f.size();is++){
      ////std::string str2("bSize");
      std::size_t found = read_f[is][0].find(str2);
      if(found != std::string::npos and read_f[is].size() >= 2){
        if(read_f[is].size()==2){res = read_f[is][1];}else{
          res = "";
          for(unsigned int temi=1;temi<read_f[is].size();temi++){
            res += read_f[is][temi] + " ";}
        }
        ////res = read_f[is][1];
        if(get_node_rank_funs0() == 0)printf("  %20s %s \n",str2.c_str(), res.c_str());
        return 1;
      }
    }
    return 0;
  }

  void load_para(const char *filename){
    int rank = get_node_rank_funs0();
    if(rank == 0)read_input(filename, read_f);
    ////===Bcast read_f;
    size_t sizen = 0;if(rank == 0)sizen = read_f.size();
    MPI_Bcast(&sizen, sizeof(size_t), MPI_CHAR, 0, MPI_COMM_WORLD);
    if(rank != 0)read_f.resize(sizen);

    for(unsigned int is=0;is<read_f.size();is++)
    {
      if(rank == 0)sizen = read_f[is].size();
      MPI_Bcast(&sizen, sizeof(size_t), MPI_CHAR, 0, MPI_COMM_WORLD);
      if(rank != 0)read_f[is].resize(sizen);
      for(unsigned int ic=0;ic<read_f[is].size();ic++)
      {
        if(rank == 0)sizen = read_f[is][ic].size();
        MPI_Bcast(&sizen, sizeof(size_t), MPI_CHAR, 0, MPI_COMM_WORLD);
        if(rank != 0)read_f[is][ic].resize(sizen);
        MPI_Bcast(&read_f[is][ic][0], sizen, MPI_CHAR, 0, MPI_COMM_WORLD);
      }
    }
    ////===Bcast read_f;

    if(get_node_rank_funs0() == 0)printf("========Start input \n");
    if(find_para(std::string("bSize"),bSize)==0)bSize = 32;
    if(find_para(std::string("bSum"),bSum)==0)bSum  = 512;
    if(find_para(std::string("cutN"),cutN)==0)cutN  = 8;
    if(find_para(std::string("lat"),lat)==0)lat  = std::string("24D");
    if(find_para(std::string("icfg"),icfg)==0)icfg  = 0;

    if(find_para(std::string("bini"),bini)==0)bini = -1;
    if(find_para(std::string("ncut0"),ncut0)==0)ncut0 = 30;
    if(find_para(std::string("ncut1"),ncut1)==0)ncut1 = 30;

    if(find_para(std::string("nx"),nx)==0)nx  = 0;
    if(find_para(std::string("ny"),ny)==0)ny  = 0;
    if(find_para(std::string("nz"),nz)==0)nz  = 0;
    if(find_para(std::string("nt"),nt)==0)nt  = 0;
    if(find_para(std::string("nvec"),nvec)==0)nvec  = 0;
    if(find_para(std::string("ionum"),ionum)==0)ionum  = 0;
    if(find_para(std::string("Ename"),Ename)==0)Ename  = std::string("NONE");
    if(find_para(std::string("Pname"),Pname)==0)Pname  = std::string("NONE");
    if(find_para(std::string("paraI"),paraI)==0)paraI  = std::string("NONE");

    if(find_para(std::string("nmass"),nmass)==0)nmass  = 0;
    if(find_para(std::string("debuga"),debuga)==0)debuga  = 0;
    if(find_para(std::string("save_prop"),save_prop)==0)save_prop  = 0;
    if(get_node_rank_funs0() == 0)printf("========End   input \n");
  }

  void load_para(int argc, char* argv[]){
    for(int i=1;i<argc;i++){
      std::string str=argv[i];
      std::size_t found = str.find(std::string("txt"));
      if(found != std::string::npos)
      {
        load_para(str.c_str());
        return;
      }
    }
  
    //////if(get_file_size_MPI(std::string("input.txt").c_str()) == 0)return;
    load_para("input.txt");return;

  }

};


template<typename Yl>
void p_vector(const Yl teml)
{
  std::cout << teml << " ";
};

template<typename Ty>
void p_vector(const std::vector<Ty> teml)
{
  for(int i=0;i< teml.size();i++)
  {
    p_vector(teml[i]);
  }
  std::cout << std::endl;
  //std::vector<std::vector<std::vector<int> > > c;
  //std::cout << typeid(c).name() << std::endl;
  //if(namev.compare("St6vectorIiSaIiEE") == 0)
};


inline void ran_EigenM(EigenM& a, int cpu=1)
{
  int N0 = a.size();if(N0 == 0)return ;
  //int N1 = a[0].size();
  timeval tm;gettimeofday(&tm, NULL);
  ////print0("print time %.3f\n", tm.tv_sec);
  qlat::RngState rs(qlat::get_id_node() + 1 +int(tm.tv_sec));
  double ini = qlat::u_rand_gen(rs);

  if(cpu == 0){
  for(long i=0;i<N0;i++){
  EigenV& tem = a[i];
  qacc_for(j, tem.size() ,{
    tem[j] = Complexq(std::cos((ini+i)*0.5) , (9.0/(j+1))*ini*0.1);
  });
  }
  //qacc_for(isp, long(N0),{
  //  long i = isp;
  //  for(long j=0;j<a[i].size();j++)
  //  {   a[i][j] = Complexq(std::cos((ini+i)*0.5) , (9.0/(j+1))*ini*0.1);}
  //});
  //qacc_for(isp, long(N0),{
  //  long i = isp;
  //  for(long j=0;j<a[i].size();j++)
  //  {   a[i][j] = Complexq(std::cos((ini+i)*0.5) , (9.0/(j+1))*ini*0.1);}
  //});
  }
  else{
  #pragma omp parallel for
  for(size_t i=0;i < size_t(N0);i++)
  {
    for(size_t j = 0; j< size_t(a[i].size()); j++)
    {
        a[i][j] = Complexq(std::cos((ini+i)*0.5) , (9.0/(j+1))*ini*0.1);
    }
  }}
}

inline void ran_EigenM(Evector& a,int cpu=1)
{
  if(a.size() == 0)return;
  timeval tm;gettimeofday(&tm, NULL);
  ////print0("print time %.3f\n", tm.tv_sec);
  qlat::RngState rs(qlat::get_id_node() + 1 +int(tm.tv_sec));

  double ini = qlat::u_rand_gen(rs);

  if(cpu == 0){
  size_t bfac = size_t(std::sqrt(a.size()));
  qacc_for(isp, long(a.size()/bfac + 1),{
   for(size_t i=0;i<size_t(bfac);i++){  
    size_t off = isp*bfac + i;
    if(off < size_t(a.size())){
      a[isp] = Complexq(std::cos((ini+isp)*0.5) , (5.0/(isp+1))*ini*0.1);
    }
  }
  });}else{
  #pragma omp parallel for
  for(size_t isp=0;isp< size_t(a.size());isp++)
  {
          a[isp] = Complexq(std::cos((ini+isp)*0.5) , (5.0/(isp+1))*ini*0.1);
  }}
}

inline void zeroE(EigenM& a,int cpu=1, bool dummy=true)
{
  int N0 = a.size();if(N0 == 0)return ;
  int N1 = a[0].size();

  #ifdef QLAT_USE_ACC
  if(cpu == 0){
    for(int i=0;i<N0;i++){cudaMemsetAsync(&a[i][0], 0, a[i].size()*sizeof(Complexq));}
    if(dummy)qacc_barrier(dummy);
    return ;
  }
  #endif

  #pragma omp parallel for
  for(long i=0;i < a.size();i++)
  for(long j=0;j < a[i].size();j++)
  {
    a[i][j] = 0.0;
  }

  //qacc_for(isp, long(N0*N1),{
  //  long i = isp/N1;long j = isp%N1;a[i][j] = 0.0;
  //});

}

inline void zeroE(Evector& a,int cpu=1, bool dummy=true)
{
  #ifdef QLAT_USE_ACC
  if(cpu == 0){
    cudaMemsetAsync(&a[0], 0, a.size()*sizeof(Complexq));
    if(dummy)qacc_barrier(dummy);
    return ;
  }
  #endif

  #pragma omp parallel for
  for(long isp=0;isp<a.size();isp++){  a[isp] = 0.0;}

  //qacc_for(isp, long(a.size()),{ 
  //  a[isp] = 0.0;
  //});
}

template<typename Ty>
inline void copy_data(Ty* res, Ty* src, size_t size, int cpu=1, bool dummy=true)
{
  #ifdef QLAT_USE_ACC
  if(cpu == 0){
    cudaMemcpyAsync(res, src, size*sizeof(Ty), cudaMemcpyDeviceToDevice);
    if(dummy)qacc_barrier(dummy);
    return ;
  }
  #endif

  #pragma omp parallel for
  for(size_t isp=0;isp<size;isp++){ res[isp] = src[isp];}
}

template <class T>
void random_link(GaugeFieldT<T >& g)
{
  timeval tm;gettimeofday(&tm, NULL);
  ////print0("print time %.3f\n", tm.tv_sec);
  qlat::RngState rs(qlat::get_id_node() + 1 +int(tm.tv_sec));
  double ini = qlat::u_rand_gen(rs);

  int dir_limit = 4;
  const Geometry& geo = g.geo();

  qacc_for(isp,  geo.local_volume(),{
    for (int dir = 0; dir < dir_limit; ++dir) {
      const Coordinate xl = geo.coordinate_from_index(isp);
      ColorMatrixT<T >& link = g.get_elem(xl, dir);
      for(int ci=0; ci<9; ci++){
        link.p[ci] = T(std::cos((ini+isp)*0.5) , ((5.0+ci)/(isp+1))*ini*0.1);;
      } 
    }   
  }); 

}

template <class T>
void random_prop(Propagator4dT<T >& prop)
{
  timeval tm;gettimeofday(&tm, NULL);
  ////print0("print time %.3f\n", tm.tv_sec);
  qlat::RngState rs(qlat::get_id_node() + 1 +int(tm.tv_sec));
  double ini = qlat::u_rand_gen(rs);

  int dir_limit = 4;
  const Geometry& geo = prop.geo();

  qacc_for(isp,  geo.local_volume(),{
    qlat::WilsonMatrixT<T >& v0 =  prop.get_elem(isp);
    for(int ci=0;ci<12*12;ci++){
      v0.p[ci] = T(std::cos((ini+isp)*0.5) , ((5.0+ci)/(isp+1))*ini*0.1); 
    }
  }); 
}



//qacc void zeroE(std::vector<Complexq >& a)
//{
//  for(unsigned int isp=0;isp<a.size();isp++)
//  {
//    a[isp] = 0.0;
//  }
//}

inline void diff_EigenM(EigenV& a, EigenV& b, std::string lab)
{
  double diff = 0.0;
  if(a.size() != b.size()){print0("%s size not equal %d, %d\n",
    lab.c_str(), int(a.size()), int(b.size()) );abort_r("");}

  for(long i=0;i < a.size(); i++)
  { 
    Complexq tem = a[i] - b[i];
    diff += (tem.real()*tem.real() + tem.imag()*tem.imag());
  } 

  sum_all_size(&diff,1);
  print0("%s diff %.5e, count %d, ava %.5e \n",
    lab.c_str(),diff, int(a.size()), diff/a.size());
}

size_t get_threads(size_t thread, size_t Total, int def=1)
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

std::vector<unsigned int > get_num_power(const size_t x,const std::vector<unsigned int >& a)
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
  qassert(test == x);
  return re;
}

//////Most power in T direction
Coordinate spread_powT(const int n, const Coordinate& Lat, const std::vector<unsigned int >& a)
{
  std::vector<unsigned int > Mpow = get_num_power(n, a);
  std::vector<std::vector<unsigned int > > Lpow;
  Lpow.resize(4);for(int i=0;i<4;i++){Lpow[i] =  get_num_power(Lat[i], a);}
  Coordinate re;
  for(int i=0;i<4;i++)re[i] = 1;

  for(int i=0;i<Mpow.size();i++){
    int num = Mpow[i];
    if(num != 0){
      int suma = Lpow[0][i] + Lpow[1][i] + Lpow[2][i] + Lpow[3][i] ;
      qassert(num < suma);
      int tem = num;
      for(int ni=0;ni<4;ni++){
        if(tem >= Lpow[4-ni-1][i]){
          re[4-ni-1] *= int(std::pow(a[i], Lpow[4-ni-1][i]));
          tem = tem - Lpow[4-ni-1][i];
        }
        else{
          re[4-ni-1] *= int(std::pow(a[i], tem));
          tem = 0;
        }
        ////if(tem >= Lpow[4-ni-1][i]){
        if(tem == 0)break;
      }
    }
  }
  return re;

}

Coordinate guess_nodeL(int n, const Coordinate& Lat, const int mode = 0)
{
  //MPI_Comm comm;comm = MPI_COMM_WORLD;
  //int rank, n;
  //MPI_Comm_rank(comm, &rank);
  //MPI_Comm_size(comm, &n);
  //int nt = Lat[3];

  std::vector<unsigned int > a;a.resize(5);
  a[0] = 2;a[1] = 3;a[2] = 5;a[3] = 7;a[4] =11;
  Coordinate re = spread_powT(n, Lat,a);
  return re;
  ///////If mode 0, the T first
  //size_t get_threads(size_t thread, size_t Total, int def=1)
  //if(rank < nt){int gt = get_threads(rank, nt, 1);}

  //Coordinate(1, 1, 1,  1)
  //size_node_list.push_back()

}

#ifdef QLAT_USE_ACC
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      qlat::displayln(qlat::ssprintf("cuda error: %s %s %d\n", cudaGetErrorString(code), file, line ));
      qassert(false);
   }
}
#endif

}

#endif
