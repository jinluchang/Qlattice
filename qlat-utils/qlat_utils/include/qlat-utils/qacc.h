#pragma once

// From https://github.com/paboyle/Grid/blob/develop/Grid/threads/Accelerator.h
// Orignal author: Peter Boyle <paboyle@ph.ed.ac.uk>

#ifndef QLAT_NO_ACC

// QLAT_USE_ACC : whether to use GPU or not
#ifdef __NVCC__
#define QLAT_USE_ACC
#endif

#ifdef __HIPCC__
#define QLAT_USE_ACC
#endif

// if QLAT_IN_ACC : code compile on device,  without formated message on error, on host with message
#ifdef __CUDA_ARCH__
#define QLAT_IN_ACC
#endif

#ifdef __HIP_DEVICE_COMPILE__
#define QLAT_IN_ACC
#endif

#endif

#include <qlat-utils/env.h>
#include <qlat-utils/qacc-translator.h>

#include <string>

namespace qlat
{  //

#ifdef QLAT_USE_ACC

#define qacc_no_inline __host__ __device__

#define qacc __host__ __device__ inline

#else

#define qacc_no_inline

#define qacc inline

#endif

enum struct MemType : Int {
  Cpu,      // CPU main memory
  Acc,      // Accelerator
  Uvm,      // Uniform virtual memory
  Comm,     // For communication on CPU
  CommAcc,  // For communication on ACC
  SIZE,
};

std::string show(const MemType mem_type);

MemType read_mem_type(const std::string& mem_type_str);

API inline MemType& get_default_mem_type()
{
  static MemType mem_type =
      read_mem_type(get_env_default("q_default_mem_type", "uvm"));
  return mem_type;
}

}  // namespace qlat
