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

//#include <Eigen/Dense>
//#include "fftw3.h"
//#include "fftw3-mpi.h"

#define Enablefloat 1

#define large_vuse Elarge_vector
#if Enablefloat==0
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
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
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
  MPI_Bcast(&sizen, sizeof(size_t), MPI_CHAR, 0, MPI_COMM_WORLD);
  return sizen;
}

//inline void sum_all_size(unsigned long int *value,int size)
//{
//  unsigned long int *work=new unsigned long int[size];
//  MPI_Allreduce(value, work, size, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
//  memcpy(value,work,size*sizeof(unsigned long int));
//  delete []work;
//}

template<typename Ty>
void sum_all_size(Ty *value,int size)
{
  Ty *work=new Ty[size];

  if(std::is_same<Ty, unsigned long>::value)MPI_Allreduce(value, work, size, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  if(std::is_same<Ty, int>::value)MPI_Allreduce(value, work, size, MPI_INT   , MPI_SUM, MPI_COMM_WORLD);
  if(std::is_same<Ty, float>::value)MPI_Allreduce(value, work, size, MPI_FLOAT , MPI_SUM, MPI_COMM_WORLD);
  if(std::is_same<Ty, double>::value)MPI_Allreduce(value, work, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  memcpy(value,work,size*sizeof(Ty));

  delete []work;
}

inline void abort_r(std::string stmp)
{
  print0("%s\n",stmp.c_str());
  MPI_Barrier(MPI_COMM_WORLD);
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
    //MPI_Barrier(MPI_COMM_WORLD);
    //fflush(stdout);
    //MPI_Finalize();
    //abort();
  }
}

inline void fflush_MPI(){
  MPI_Barrier(MPI_COMM_WORLD);
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

}

#endif
