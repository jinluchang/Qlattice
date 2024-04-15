//// https://github.com/ROCm-Developer-Tools/HIPIFY/blob/amd-staging/docs/tables/CUDA_Runtime_API_functions_supported_by_HIP.md
//// https://github.com/ROCm-Developer-Tools/HIPIFY/blob/amd-staging/docs/tables/CUFFT_API_supported_by_HIP.md

#ifndef QLAT_GPU_TRANSLATOR_H
#define QLAT_GPU_TRANSLATOR_H
#pragma once

#ifdef QLAT_USE_ACC

#ifdef __NVCC__

#define qacc_DeviceProp                 cudaDeviceProp                
#define qacc_DeviceReset                cudaDeviceReset               
#define qacc_DeviceSetCacheConfig       cudaDeviceSetCacheConfig      
#define qacc_DeviceSynchronize          cudaDeviceSynchronize         
#define qacc_ErrorCudartUnloading       cudaErrorCudartUnloading      
#define qacc_Error                      cudaError                     
#define qacc_Free                       cudaFree                      
#define qacc_FuncCachePreferL1          cudaFuncCachePreferL1         
#define qacc_FuncCachePreferNone        cudaFuncCachePreferNone       
#define qacc_GetDevice                  cudaGetDevice                 
#define qacc_GetDeviceCount             cudaGetDeviceCount            
#define qacc_GetDeviceFlags             cudaGetDeviceFlags            
#define qacc_GetDeviceProperties        cudaGetDeviceProperties       
#define qacc_GetErrorString             cudaGetErrorString            
#define qacc_GetLastError               cudaGetLastError              
#define qacc_Malloc                     cudaMalloc                    
#define qacc_Malloc3D                   cudaMalloc3D                  
#define qacc_MallocManaged              cudaMallocManaged             
#define qacc_MemAdvise                  cudaMemAdvise                 
#define qacc_MemAdviseSetReadMostly     cudaMemAdviseSetReadMostly    
#define qacc_MemAdviseUnsetReadMostly   cudaMemAdviseUnsetReadMostly  
#define qacc_MemGetInfo                 cudaMemGetInfo                
#define qacc_MemPrefetchAsync           cudaMemPrefetchAsync          
#define qacc_Memcpy2DAsync              cudaMemcpy2DAsync             
#define qacc_MemcpyAsync                cudaMemcpyAsync               
#define qacc_MemcpyDeviceToDevice       cudaMemcpyDeviceToDevice      
#define qacc_MemcpyDeviceToHost         cudaMemcpyDeviceToHost        
#define qacc_MemcpyHostToDevice         cudaMemcpyHostToDevice        
#define qacc_MemcpyHostToHost           cudaMemcpyHostToHost          
#define qacc_MemcpyKind                 cudaMemcpyKind                
#define qacc_MemcpyToSymbol             cudaMemcpyToSymbol            
#define qacc_MemsetAsync                cudaMemsetAsync
#define qacc_SetDevice                  cudaSetDevice                 
#define qacc_StreamCreate               cudaStreamCreate              
#define qacc_Stream_t                   cudaStream_t                  
#define qacc_Success                    cudaSuccess   
#define qacc_StreamDestroy              cudaStreamDestroy

#else

#include <hip/hip_runtime.h>
#include <hip/hip_runtime_api.h>

#define qacc_DeviceProp                 hipDeviceProp_t
#define qacc_DeviceReset                hipDeviceReset
#define qacc_DeviceSetCacheConfig       hipDeviceSetCacheConfig
#define qacc_DeviceSynchronize          hipDeviceSynchronize
#define qacc_ErrorCudartUnloading       hipErrorDeinitialized
#define qacc_Error                      hipError_t
#define qacc_Free                       hipFree  
#define qacc_FuncCachePreferL1          hipFuncCachePreferL1
#define qacc_FuncCachePreferNone        hipFuncCachePreferNone
#define qacc_GetDevice                  hipGetDevice
#define qacc_GetDeviceCount             hipGetDeviceCount
#define qacc_GetDeviceFlags             hipGetDeviceFlags
#define qacc_GetDeviceProperties        hipGetDeviceProperties
#define qacc_GetErrorString             hipGetErrorString
#define qacc_GetLastError               hipGetLastError
#define qacc_Malloc                     hipMalloc
#define qacc_Malloc3D                   hipMalloc3D
#define qacc_MallocManaged              hipMallocManaged
#define qacc_MemAdvise                  hipMemAdvise
#define qacc_MemAdviseSetReadMostly     hipMemAdviseSetReadMostly
#define qacc_MemAdviseUnsetReadMostly   hipMemAdviseUnsetReadMostly
#define qacc_MemGetInfo                 hipMemGetInfo
#define qacc_MemPrefetchAsync           hipMemPrefetchAsync
#define qacc_Memcpy2DAsync              hipMemcpy2DAsync
#define qacc_MemcpyAsync                hipMemcpyAsync
#define qacc_MemcpyDeviceToDevice       hipMemcpyDeviceToDevice
#define qacc_MemcpyDeviceToHost         hipMemcpyDeviceToHost
#define qacc_MemcpyHostToDevice         hipMemcpyHostToDevice
#define qacc_MemcpyHostToHost           hipMemcpyHostToHost
#define qacc_MemcpyKind                 hipMemcpyKind
#define qacc_MemcpyToSymbol             hipMemcpyToSymbol
#define qacc_MemsetAsync                hipMemsetAsync
#define qacc_SetDevice                  hipSetDevice
#define qacc_StreamCreate               hipStreamCreate
#define qacc_Stream_t                   hipStream_t
#define qacc_Success                    hipSuccess
#define qacc_StreamDestroy              hipStreamDestroy

#endif

#endif

#endif
