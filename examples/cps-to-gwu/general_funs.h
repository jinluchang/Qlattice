// general_funs.h
// Gen Wang
// Jan. 2021

#ifndef general_funs_h
#define general_funs_h
#pragma once

#include <string.h>
#include <sys/resource.h>
#include <mpi.h>
#include <time.h>
#include <typeinfo>
#include <qlat/qlat.h>

#include <iterator>
#include <sys/sysinfo.h>

//#include <cuda_runtime.h>

//#include <Eigen/Dense>
//#include "fftw3.h"
//#include "fftw3-mpi.h"


#define Enablefloat 1

#define large_vuse Elarge_vector
#if Enablefloat==0
#define Complexq qlat::Complex
#define Ftype double
#define Ctype double
#define CMPI                  MPI_DOUBLE
#define FFTtype               fftw
#define FFTtype_malloc        fftw_malloc
#define FFTtype_complex       fftw_complex
#define FFTtype_plan          fftw_plan
#define FFTtype_plan_dft_3d   fftw_plan_dft_3d
#define FFTtype_plan_many_dft fftw_plan_many_dft
#define FFTtype_execute       fftw_execute
#define FFTtype_free          fftw_free
#define FFTtype_destroy_plan  fftw_destroy_plan
#define qcd_vec               qcd::vector
#define EigenM Eigen::Matrix< std::complex< double >, Eigen::Dynamic, Eigen::Dynamic ,Eigen::RowMajor>
#define EigenMD Eigen::Matrix< std::complex< double >, Eigen::Dynamic, Eigen::Dynamic ,Eigen::RowMajor>
#define EigenV Eigen::VectorXcd
#define Evector Eigen::Array<std::complex<double >,Eigen::Dynamic,1>
#define EM Eigen::Map<Eigen::Array<std::complex<double >,Eigen::Dynamic,1> >
#define FFT_init_threads fftw_init_threads
#define FFT_plan_with_nthreads fftw_plan_with_nthreads

#define FFT_mpi_init              fftw_mpi_init
#define FFT_mpi_local_size_many   fftw_mpi_local_size_many
#define FFT_alloc_complex         fftw_alloc_complex
#define FFT_mpi_plan_many_dft     fftw_mpi_plan_many_dft
#endif

#if Enablefloat==1
#define Complexq qlat::ComplexF
#define Ftype float
#define Ctype float
#define CMPI                  MPI_FLOAT
#define FFTtype               fftwf
#define FFTtype_malloc        fftwf_malloc
#define FFTtype_complex       fftwf_complex
#define FFTtype_plan          fftwf_plan
#define FFTtype_plan_dft_3d   fftwf_plan_dft_3d
#define FFTtype_plan_many_dft fftwf_plan_many_dft
#define FFTtype_execute       fftwf_execute
#define FFTtype_free          fftwf_free
#define FFTtype_destroy_plan  fftwf_destroy_plan
#define qcd_vec               qcd::fvector
#define EigenM Eigen::Matrix< std::complex< float>, Eigen::Dynamic, Eigen::Dynamic ,Eigen::RowMajor>
#define EigenMD Eigen::Matrix< std::complex< float>, Eigen::Dynamic, Eigen::Dynamic ,Eigen::RowMajor>
#define EigenV Eigen::VectorXcf
#define Evector Eigen::Array<std::complex<float >,Eigen::Dynamic,1>
//#define EM Eigen::Map<Eigen::Array<std::complex<float >,Eigen::Dynamic,1> >
#define EM Eigen::Map<Eigen::Array<std::complex<float >,Eigen::Dynamic,1> >
#define FFT_init_threads fftwf_init_threads
#define FFT_plan_with_nthreads fftwf_plan_with_nthreads

#define FFT_mpi_init              fftwf_mpi_init
#define FFT_mpi_local_size_many   fftwf_mpi_local_size_many
#define FFT_alloc_complex         fftwf_alloc_complex
#define FFT_mpi_plan_many_dft     fftwf_mpi_plan_many_dft
#endif

#define ComtypeT std::complex

#define LInt unsigned long int


namespace qlat
{

inline unsigned int get_node_rank_funs()
{
  int rank;
  MPI_Comm_rank(get_comm(), &rank);
  return rank;
}

#define print0 if(get_node_rank_funs() == 0) printf

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
  if(get_node_rank_funs()==0){
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
void sum_all_size(Ty *src,Ty *sav,long size)
{
  Ty *res;
  if(src == sav){res = new Ty[size];}else{res = sav;}
  if(qlat::get_num_node() == 1){
    if(src == sav){return;}
    if(src != sav){memcpy(sav,src,size*sizeof(Ty));return;}
  }

  if(std::is_same<Ty, unsigned long>::value)MPI_Allreduce(src,res, size, MPI_UNSIGNED_LONG, MPI_SUM, get_comm());
  if(std::is_same<Ty, int>::value)MPI_Allreduce(src,res, size, MPI_INT   , MPI_SUM, get_comm());
  if(std::is_same<Ty, float>::value)MPI_Allreduce(src,res, size, MPI_FLOAT , MPI_SUM, get_comm());
  if(std::is_same<Ty, double>::value)MPI_Allreduce(src,res, size, MPI_DOUBLE, MPI_SUM, get_comm());

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

inline void print_mem_info()
{
  ////double length = (geo.local_volume()*pow(0.5,30))*12*sizeof(Complexq);
  size_t freeM = 0;size_t totalM = 0;
  double freeD = 0;double totalD=0;
  #ifdef QLAT_USE_ACC
  cudaMemGetInfo(&freeM,&totalM);
  freeD = freeM*pow(0.5,30);totalD = totalM*pow(0.5,30);
  #endif
  struct sysinfo s_info;
  sysinfo(&s_info);
  print0("===CPU free %.3e GB, total %.3e GB; GPU free %.3e GB, total %.3e GB. \n"
          ,s_info.totalram*pow(0.5,30), s_info.freeram*pow(0.5,30),freeD,totalD);
}


std::vector<std::string > stringtolist(std::string &tem_string)
{
  std::istringstream iss(tem_string);
  std::vector<std::string> results((std::istream_iterator<std::string>(iss)),std::istream_iterator<std::string>());
  return results;
}

double stringtodouble(std::string &tem_string)
{
  //double use = atof(tem_string.c_str());
  double use = 0.0;
  if(tem_string!="_NONE_")
  {
    use = atof(tem_string.c_str());
  }
  return use;
}

int stringtonum(std::string &tem_string)
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

int read_vector(const char *filename, std::vector<double > &dat)
{
  int prods = 0; 
  unsigned long Vsize = 0; 
  ////{synchronize();fflush(stdout);}

  if(get_node_rank_funs() == 0)
  {
    Vsize = get_file_size_o(filename);
    if(Vsize == 0){prods = 0;}else{prods=Vsize;}
    Vsize = Vsize/8;
  }
  sum_all_size(&prods,1);
  if(prods==0)return prods;

  sum_all_size((int*)&Vsize,1);
  dat.resize(Vsize);

  if(get_node_rank_funs() == 0)
  {
    FILE* filer = fopen(filename, "rb");
    //std::vector<double > tem;
    //tem.resize(3);
    unsigned long count = 1024*1024;
    unsigned long sizec = 0; 
    unsigned long offr  = 0; 

    //char* buf = (char *)&dat[0];
    //fread(buf, 1, Vsize*8, filer);

    for(int iv=0;iv<Vsize;iv++)
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
void bcast_vstring(std::vector<std::string> &conf_l, const int Host_rank = 0){

  int rank = get_node_rank_funs();
  ////Bcast strings
  size_t sizen = 0;if(rank == Host_rank)sizen = conf_l.size();
  MPI_Bcast(&sizen, sizeof(size_t), MPI_CHAR, Host_rank, get_comm());
  if(rank != Host_rank)conf_l.resize(sizen);
  for(int is=0;is<conf_l.size();is){
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
  int conf;
  int save_prop;
  ~inputpara(){
    for(int is=0;is<read_f.size();is++){
     for(int ic=0;ic<read_f[is].size();ic++){read_f[is][ic].resize(0);}
     read_f[is].resize(0);
    }
    read_f.resize(0);
  }

  int find_para(const std::string &str2, int &res){
    for(int is=0;is<read_f.size();is++){
      ////std::string str2("bSize");
      std::size_t found = read_f[is][0].find(str2);
      if(found != std::string::npos and read_f[is].size() >= 2){
        res = stringtonum(read_f[is][1]);
        print0("  %10s %10d \n",str2.c_str(), res);
        return 1;
      }
    }
    return 0;
  }

  int find_para(const std::string &str2, std::string &res){
    for(int is=0;is<read_f.size();is++){
      ////std::string str2("bSize");
      std::size_t found = read_f[is][0].find(str2);
      if(found != std::string::npos and read_f[is].size() >= 2){
        res = read_f[is][1];
        print0("  %10s %10s \n",str2.c_str(), res.c_str());
        return 1;
      }
    }
    return 0;
  }

  void load_para(const char *filename){
    int rank = get_node_rank_funs();
    if(rank == 0)read_input(filename, read_f);
    ////===Bcast read_f;
    size_t sizen = 0;if(rank == 0)sizen = read_f.size();
    MPI_Bcast(&sizen, sizeof(size_t), MPI_CHAR, 0, get_comm());
    if(rank != 0)read_f.resize(sizen);

    for(int is=0;is<read_f.size();is++)
    {
      if(rank == 0)sizen = read_f[is].size();
      MPI_Bcast(&sizen, sizeof(size_t), MPI_CHAR, 0, get_comm());
      if(rank != 0)read_f[is].resize(sizen);
      for(int ic=0;ic<read_f[is].size();ic++)
      {
        if(rank == 0)sizen = read_f[is][ic].size();
        MPI_Bcast(&sizen, sizeof(size_t), MPI_CHAR, 0, get_comm());
        if(rank != 0)read_f[is][ic].resize(sizen);
        MPI_Bcast(&read_f[is][ic][0], sizen, MPI_CHAR, 0, get_comm());
      }
    }
    ////===Bcast read_f;

    print0("========Start input \n");
    if(find_para(std::string("bSize"),bSize)==0)bSize = 32;
    if(find_para(std::string("bSum"),bSum)==0)bSum  = 512;
    if(find_para(std::string("cutN"),cutN)==0)cutN  = 8;
    if(find_para(std::string("lat"),lat)==0)lat  = std::string("24D");
    if(find_para(std::string("conf"),conf)==0)conf  = 0;
    if(find_para(std::string("save_prop"),save_prop)==0)save_prop  = 0;
    print0("========End   input \n");
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
  
    if(get_file_size_MPI(std::string("input.txt").c_str()) == 0)return;
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


////inutpara input_global;

}

#endif
