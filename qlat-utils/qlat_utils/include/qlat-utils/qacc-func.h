#pragma once

// From https://github.com/paboyle/Grid/blob/develop/Grid/threads/Accelerator.h
// Orignal author: Peter Boyle <paboyle@ph.ed.ac.uk>

#include <qlat-utils/qacc.h>
#include <qlat-utils/timer.h>

namespace qlat
{  //

API inline int& qacc_num_threads()
// qlat parameter
{
  static int nt = get_env_long_default("q_acc_num_threads", 32);
  return nt;
}

#ifdef QLAT_USE_ACC

inline void gpuErr(qacc_Error err, const char *file, int line)
{
  if (qacc_Success != err) {
    qlat::displayln(
        qlat::ssprintf("qacc_barrier: ACC error %s from '%s' Line %d.",
                       qacc_GetErrorString(err), file, line));
    qassert(false);
  }
}

#ifdef __NVCC__

#define gpuErrCheck(ans) { gpuErr((ans), __FILE__, __LINE__); }

#else

#define gpuErrCheck(ans) { (ans);gpuErr(qacc_GetLastError(), __FILE__, __LINE__); }

#endif

//inline void qacc_DeviceSynchronize()
//{
//  cudaDeviceSynchronize();
//  qacc_Error err = qacc_GetLastError();
//  if (qacc_Success != err) {
//    qlat::displayln(
//        qlat::ssprintf("qacc_barrier: ACC error %s from '%s' Line %d.",
//                       qacc_GetErrorString(err), __FILE__, __LINE__));
//    qassert(false);
//  }
//}

//inline void qacc_DeviceSynchronize()
//{
//  qacc_Error err = hipDeviceSynchronize();
//  if (qacc_Success != err) {
//    qlat::displayln(
//        qlat::ssprintf("qacc_barrier: ACC error %s from '%s' Line %d.",
//                       qacc_GetErrorString(err), __FILE__, __LINE__));
//    qassert(false);
//  }
//}

inline void qacc_DeviceSynchronize()
{
#ifdef __NVCC__
  gpuErrCheck(cudaDeviceSynchronize());
#else
  gpuErrCheck(hipDeviceSynchronize());
#endif
}

#endif

API inline MemType check_mem_type(void* ptr)
{
  MemType mem_type;
#ifdef QLAT_USE_ACC
  bool find = false;
  qacc_PointerAttributes attr;
  gpuErrCheck(qacc_PointerGetAttributes(&attr, ptr));
  if (attr.type == qacc_MemoryTypeHost) {
    mem_type = MemType::Cpu;find = true;
  }
  if (attr.type == qacc_MemoryTypeDevice) {
    mem_type = MemType::Acc;find = true;
  }
  if (attr.type == qacc_MemoryTypeManaged) {
    mem_type = MemType::Uvm;find = true;
  }
  if(!find){assert(false);}
#else
  (void)ptr;
  mem_type = MemType::Cpu;
#endif
  return mem_type;
}


#define qfor(iter, num, ...)                                   \
  {                                                            \
    const qlat::Long QACC_NUM_VALUE = (num);                   \
    for (qlat::Long iter = 0; iter < QACC_NUM_VALUE; ++iter) { \
      {__VA_ARGS__};                                           \
    }                                                          \
  }

#define q_do_pragma(x) _Pragma(#x)

#define qthread_for(iter, num, ...)                          \
  {                                                          \
    const qlat::Long QACC_NUM_VALUE = (num);                 \
    q_do_pragma(omp parallel for schedule(static))           \
    for (qlat::Long iter = 0; iter < QACC_NUM_VALUE; ++iter) \
    {                                                        \
      {__VA_ARGS__};                                         \
    }                                                        \
  }

#ifdef QLAT_USE_ACC

#ifndef __CUDACC_EXTENDED_LAMBDA__
#error "please compile with --expt-extended-lambda"
#endif

#define qacc_continue return

#define qacc_forNB(iter, num, ...)                                             \
  {                                                                            \
    const qlat::Long QACC_NUM_VALUE = (num);                                   \
    if (QACC_NUM_VALUE != 0) {                                                 \
      auto QACC_FOR_LOOP_LAMBDA =                                              \
          [=] __host__ __device__(qlat::Long iter) mutable { {__VA_ARGS__}; }; \
      const int QACC_NUM_THREADS = qlat::qacc_num_threads();                   \
      dim3 CUDA_THREADS_DIM3(QACC_NUM_THREADS, 1, 1);                          \
      dim3 CUDA_BLOCKS_DIM3(                                                   \
          (QACC_NUM_VALUE + QACC_NUM_THREADS - 1) / QACC_NUM_THREADS, 1, 1);   \
      qacc_Error CUDA_LAST_ERROR_OBJECT = qacc_GetLastError();                 \
      if (qacc_Success != CUDA_LAST_ERROR_OBJECT) {                            \
        qerr(qlat::ssprintf("qacc_for: ACC error %s from '%s' Line %d.",       \
                            qacc_GetErrorString(CUDA_LAST_ERROR_OBJECT),       \
                            __FILE__, __LINE__));                              \
      }                                                                        \
      qlambda_apply<<<CUDA_BLOCKS_DIM3, CUDA_THREADS_DIM3>>>(                  \
          QACC_NUM_VALUE, QACC_FOR_LOOP_LAMBDA);                               \
    }                                                                          \
  }

template <typename Lambda>
__global__ void qlambda_apply(Long num, Lambda lam)
{
  Long x = threadIdx.x + blockDim.x * blockIdx.x;
  if (x < num) {
    lam(x);
  }
}

#define qacc_barrier(dummy)   \
  {                           \
    qacc_DeviceSynchronize(); \
  }

#define qacc_for(iter, num, ...) \
  qacc_forNB(iter, num, {__VA_ARGS__}) qacc_barrier(dummy)

#else

#define qacc_continue continue

#define qacc_forNB(iter, num, ...) qthread_for(iter, num, {__VA_ARGS__})

#define qacc_barrier(dummy)

#define qacc_for(iter, num, ...) \
  qacc_forNB(iter, num, {__VA_ARGS__}) qacc_barrier(dummy);

#endif

#define qmem_for(iter, num, mem_type, ...)                              \
  {                                                                     \
    qlat::MemType QACC_EFF_MEM_TYPE = qlat::get_eff_mem_type(mem_type); \
    if (QACC_EFF_MEM_TYPE == qlat::MemType::Cpu) {                      \
      qthread_for(iter, (num), {__VA_ARGS__});                          \
    } else {                                                            \
      qacc_for(iter, (num), {__VA_ARGS__});                             \
    }                                                                   \
  }

}  // namespace qlat
