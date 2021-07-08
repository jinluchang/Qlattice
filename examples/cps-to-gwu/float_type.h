// float_type.h
// Gen Wang
// Jun. 2021

#ifndef FLOAT_TYPE_H
#define FLOAT_TYPE_H
#pragma once

#include <string.h>
#include <sys/resource.h>
#include <mpi.h>
#include <time.h>
#include <typeinfo>
#include <qlat/qcd.h>

#include <iterator>
#include <sys/sysinfo.h>


////#include <cuda_runtime.h>
//
////#include <Eigen/Dense>
////#include "fftw3.h"
////#include "fftw3-mpi.h"

#define LInt unsigned long

#define Enablefloat 0

#define large_vuse Elarge_vector
#if Enablefloat==0
#define Complexq qlat::Complex
#define Ftype double
//#define CMPI                  MPI_DOUBLE
//#define FFTtype               fftw
//#define FFTtype_malloc        fftw_malloc
//#define FFTtype_complex       fftw_complex
//#define FFTtype_plan          fftw_plan
//#define FFTtype_plan_dft_3d   fftw_plan_dft_3d
//#define FFTtype_plan_many_dft fftw_plan_many_dft
//#define FFTtype_execute       fftw_execute
//#define FFTtype_free          fftw_free
//#define FFTtype_destroy_plan  fftw_destroy_plan
//
//#define FFT_init_threads fftw_init_threads
//#define FFT_plan_with_nthreads fftw_plan_with_nthreads
//
//#define FFT_mpi_init              fftw_mpi_init
//#define FFT_mpi_local_size_many   fftw_mpi_local_size_many
//#define FFT_alloc_complex         fftw_alloc_complex
//#define FFT_mpi_plan_many_dft     fftw_mpi_plan_many_dft

#endif

#if Enablefloat==1
#define Complexq qlat::ComplexF
#define Ftype float
//#define CMPI                  MPI_FLOAT
//#define FFTtype               fftwf
//#define FFTtype_malloc        fftwf_malloc
//#define FFTtype_complex       fftwf_complex
//#define FFTtype_plan          fftwf_plan
//#define FFTtype_plan_dft_3d   fftwf_plan_dft_3d
//#define FFTtype_plan_many_dft fftwf_plan_many_dft
//#define FFTtype_execute       fftwf_execute
//#define FFTtype_free          fftwf_free
//#define FFTtype_destroy_plan  fftwf_destroy_plan
//
//#define FFT_init_threads fftwf_init_threads
//#define FFT_plan_with_nthreads fftwf_plan_with_nthreads
//
//#define FFT_mpi_init              fftwf_mpi_init
//#define FFT_mpi_local_size_many   fftwf_mpi_local_size_many
//#define FFT_alloc_complex         fftwf_alloc_complex
//#define FFT_mpi_plan_many_dft     fftwf_mpi_plan_many_dft

#endif

////#define Evector qlat::vector<std::complex< Ftype > >
////#define EA Eigen::Map<Eigen::Array<std::complex<Ftype >,Eigen::Dynamic,1 > >
////#define EM Eigen::Map< Eigen::Matrix< std::complex<Ftype >, Eigen::Dynamic, Eigen::Dynamic ,Eigen::RowMajor> > 

#define Evector qlat::vector<Complexq >
#define EigenV  qlat::vector<Complexq >
#define EA Eigen::Map<Eigen::Array<Complexq ,Eigen::Dynamic,1 > >
#define EM Eigen::Map< Eigen::Matrix<Complexq , Eigen::Dynamic, Eigen::Dynamic ,Eigen::RowMajor> > 
#define EMC Eigen::Map< Eigen::Matrix<Complexq , Eigen::Dynamic, Eigen::Dynamic ,Eigen::ColMajor> > 
////#define EM Eigen::Map< Eigen::Matrix<std::complex<Ftype > , Eigen::Dynamic, Eigen::Dynamic ,Eigen::RowMajor> > 
////#define EMC Eigen::Map< Eigen::Matrix<std::complex<Ftype > , Eigen::Dynamic, Eigen::Dynamic ,Eigen::ColMajor> > 

#define EigenM qlat::vector<Evector >
////#define EigenM std::vector<Evector >

#define EG  Eigen::Matrix<Complexq , Eigen::Dynamic, Eigen::Dynamic ,Eigen::RowMajor>  
#define EGC Eigen::Matrix<Complexq , Eigen::Dynamic, Eigen::Dynamic ,Eigen::ColMajor> 




#endif
