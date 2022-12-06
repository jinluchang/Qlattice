// utils_construction.h
// Gen Wang
// Jul. 2021

#ifndef UTILS_CONSTRUCTION_H
#define UTILS_CONSTRUCTION_H

#pragma once

#include "utils_float_type.h"
#include "utils_gammas.h"
#include "utils_fft_desc.h"
#include "utils_reduce_vec.h"
#include "utils_grid_src.h"
#include "utils_io_vec.h"


#ifdef QLAT_USE_ACC
#define USEKERNEL 1
#define USEGLOBAL 1
#define USEQACC   1
#else
#define USEKERNEL 0
#define USEGLOBAL 0
#define USEQACC   0
#endif

#define EigenTy std::vector<qlat::vector_gpu<Ty > >

///#define EigenMTa std::vector<qlat::vector_acc<Ta > >
//#define EigenVTa qlat::vector_acc<Ta >
#define EAy   Eigen::Map<Eigen::Array<Ty ,Eigen::Dynamic,1 > >
//#define EAa   Eigen::Map<Eigen::Array<Ta ,Eigen::Dynamic,1 > >

#include "utils_corr_prop.h"
#include "utils_corr_meson.h"
#include "utils_corr_baryon.h"
#include "utils_corr_seq.h"

#undef  EigenTy
///#undef  EigenMTa
//#undef  EigenVTa
#undef  EAy
//#undef  EAa


#endif

