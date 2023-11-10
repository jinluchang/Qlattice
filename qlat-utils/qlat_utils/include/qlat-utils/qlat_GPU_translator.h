//// https://github.com/ROCm-Developer-Tools/HIPIFY/blob/amd-staging/docs/tables/CUDA_Runtime_API_functions_supported_by_HIP.md
//// https://github.com/ROCm-Developer-Tools/HIPIFY/blob/amd-staging/docs/tables/CUFFT_API_supported_by_HIP.md

#ifndef QLAT_GPU_TRANSLATOR_H
#define QLAT_GPU_TRANSLATOR_H
#pragma once

#ifdef QLAT_USE_ACC

#ifdef __NVCC__

#define qlat_GPU_DeviceProp                 cudaDeviceProp                
#define qlat_GPU_DeviceReset                cudaDeviceReset               
#define qlat_GPU_DeviceSetCacheConfig       cudaDeviceSetCacheConfig      
#define qlat_GPU_DeviceSynchronize          cudaDeviceSynchronize         
#define qlat_GPU_ErrorCudartUnloading       cudaErrorCudartUnloading      
#define qlat_GPU_Error                      cudaError                     
#define qlat_GPU_Free                       cudaFree                      
#define qlat_GPU_FuncCachePreferL1          cudaFuncCachePreferL1         
#define qlat_GPU_FuncCachePreferNone        cudaFuncCachePreferNone       
#define qlat_GPU_GetDevice                  cudaGetDevice                 
#define qlat_GPU_GetDeviceCount             cudaGetDeviceCount            
#define qlat_GPU_GetDeviceFlags             cudaGetDeviceFlags            
#define qlat_GPU_GetDeviceProperties        cudaGetDeviceProperties       
#define qlat_GPU_GetErrorString             cudaGetErrorString            
#define qlat_GPU_GetLastError               cudaGetLastError              
#define qlat_GPU_Malloc                     cudaMalloc                    
#define qlat_GPU_Malloc3D                   cudaMalloc3D                  
#define qlat_GPU_MallocManaged              cudaMallocManaged             
#define qlat_GPU_MemAdvise                  cudaMemAdvise                 
#define qlat_GPU_MemAdviseSetReadMostly     cudaMemAdviseSetReadMostly    
#define qlat_GPU_MemAdviseUnsetReadMostly   cudaMemAdviseUnsetReadMostly  
#define qlat_GPU_MemGetInfo                 cudaMemGetInfo                
#define qlat_GPU_MemPrefetchAsync           cudaMemPrefetchAsync          
#define qlat_GPU_Memcpy2DAsync              cudaMemcpy2DAsync             
#define qlat_GPU_MemcpyAsync                cudaMemcpyAsync               
#define qlat_GPU_MemcpyDeviceToDevice       cudaMemcpyDeviceToDevice      
#define qlat_GPU_MemcpyDeviceToHost         cudaMemcpyDeviceToHost        
#define qlat_GPU_MemcpyHostToDevice         cudaMemcpyHostToDevice        
#define qlat_GPU_MemcpyHostToHost           cudaMemcpyHostToHost          
#define qlat_GPU_MemcpyKind                 cudaMemcpyKind                
#define qlat_GPU_MemcpyToSymbol             cudaMemcpyToSymbol            
#define qlat_GPU_MemsetAsync                cudaMemsetAsync
#define qlat_GPU_SetDevice                  cudaSetDevice                 
#define qlat_GPU_StreamCreate               cudaStreamCreate              
#define qlat_GPU_Stream_t                   cudaStream_t                  
#define qlat_GPU_Success                    cudaSuccess   
#define qlat_GPU_StreamDestroy              cudaStreamDestroy

#else

#include <hip/hip_runtime.h>
#include <hip/hip_runtime_api.h>

#define qlat_GPU_DeviceProp                 hipDeviceProp_t
#define qlat_GPU_DeviceReset                hipDeviceReset
#define qlat_GPU_DeviceSetCacheConfig       hipDeviceSetCacheConfig
#define qlat_GPU_DeviceSynchronize          hipDeviceSynchronize
#define qlat_GPU_ErrorCudartUnloading       hipErrorDeinitialized
#define qlat_GPU_Error                      hipError_t
#define qlat_GPU_Free                       hipFree  
#define qlat_GPU_FuncCachePreferL1          hipFuncCachePreferL1
#define qlat_GPU_FuncCachePreferNone        hipFuncCachePreferNone
#define qlat_GPU_GetDevice                  hipGetDevice
#define qlat_GPU_GetDeviceCount             hipGetDeviceCount
#define qlat_GPU_GetDeviceFlags             hipGetDeviceFlags
#define qlat_GPU_GetDeviceProperties        hipGetDeviceProperties
#define qlat_GPU_GetErrorString             hipGetErrorString
#define qlat_GPU_GetLastError               hipGetLastError
#define qlat_GPU_Malloc                     hipMalloc
#define qlat_GPU_Malloc3D                   hipMalloc3D
#define qlat_GPU_MallocManaged              hipMallocManaged
#define qlat_GPU_MemAdvise                  hipMemAdvise
#define qlat_GPU_MemAdviseSetReadMostly     hipMemAdviseSetReadMostly
#define qlat_GPU_MemAdviseUnsetReadMostly   hipMemAdviseUnsetReadMostly
#define qlat_GPU_MemGetInfo                 hipMemGetInfo
#define qlat_GPU_MemPrefetchAsync           hipMemPrefetchAsync
#define qlat_GPU_Memcpy2DAsync              hipMemcpy2DAsync
#define qlat_GPU_MemcpyAsync                hipMemcpyAsync
#define qlat_GPU_MemcpyDeviceToDevice       hipMemcpyDeviceToDevice
#define qlat_GPU_MemcpyDeviceToHost         hipMemcpyDeviceToHost
#define qlat_GPU_MemcpyHostToDevice         hipMemcpyHostToDevice
#define qlat_GPU_MemcpyHostToHost           hipMemcpyHostToHost
#define qlat_GPU_MemcpyKind                 hipMemcpyKind
#define qlat_GPU_MemcpyToSymbol             hipMemcpyToSymbol
#define qlat_GPU_MemsetAsync                hipMemsetAsync
#define qlat_GPU_SetDevice                  hipSetDevice
#define qlat_GPU_StreamCreate               hipStreamCreate
#define qlat_GPU_Stream_t                   hipStream_t
#define qlat_GPU_Success                    hipSuccess
#define qlat_GPU_StreamDestroy              hipStreamDestroy

#endif

#endif

#endif
