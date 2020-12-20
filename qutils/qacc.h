#pragma once

namespace qlat
{  //

#ifdef QLAT_USE_GPU

#define qacc_no_inline __host__ __device__

#define qacc __host__ __device__ inline

#else

#define qacc_no_inline

#define qacc inline

#endif

}  // namespace qlat
