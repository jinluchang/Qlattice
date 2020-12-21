#pragma once

// From https://github.com/paboyle/Grid/blob/develop/Grid/threads/Accelerator.h
// Orignal author: Peter Boyle <paboyle@ph.ed.ac.uk>

#ifdef __NVCC__
#define QLAT_USE_ACC
#endif

namespace qlat
{  //

#ifdef QLAT_USE_ACC

#define qacc_no_inline __host__ __device__

#define qacc __host__ __device__ inline

#else

#define qacc_no_inline

#define qacc inline

#endif

}  // namespace qlat
