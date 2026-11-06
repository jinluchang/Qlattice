// Precision-preserving arithmetic for the opt-in local CG and simple stencil.
#ifndef UTILS_CG_PRECISION_H
#define UTILS_CG_PRECISION_H
#pragma once

#include "general_funs.h"
#include <algorithm>
#include <cmath>
#include <type_traits>
#include <vector>

namespace qlat {

// RealDD stores the leading component in Y() and the trailing component in
// X().  The installed device operator+= does not always normalize the trailing
// component correctly.  Use TwoSum explicitly in this independent solver;
// existing solver arithmetic is deliberately untouched.
qacc void cg_two_sum(double a, double b, double& sum, double& error) {
  sum = a + b;
  const double delta = sum - a;
  error = (a - (sum - delta)) + (b - delta);
}

template<class Real> qacc Real cg_real_add(const Real& a, const Real& b) {
  return a + b;
}
template<> qacc RealDD cg_real_add<RealDD>(const RealDD& a, const RealDD& b) {
  double high, low, tail, tail_error;
  cg_two_sum(a.Y(), b.Y(), high, low);
  cg_two_sum(a.X(), b.X(), tail, tail_error);
  low += tail;
  cg_two_sum(high, low, tail, low);
  low += tail_error;
  RealDD value;
  cg_two_sum(tail, low, value.Y(), value.X());
  return value;
}
template<class Real> qacc Real cg_real_neg(const Real& a) { return -a; }
template<> qacc RealDD cg_real_neg<RealDD>(const RealDD& a) {
  RealDD value;
  value.Y() = -a.Y();
  value.X() = -a.X();
  return value;
}
template<class Real> qacc Real cg_real_sub(const Real& a, const Real& b) {
  return cg_real_add(a, cg_real_neg(b));
}
template<class Real> qacc Real cg_real_mul(const Real& a, const Real& b) {
  return a * b;
}
template<class Real> qacc Real cg_real_div(const Real& a, const Real& b) {
  return a / b;
}
template<class Real> qacc bool cg_real_zero(const Real& a) {
  return a == Real(0.0);
}
template<> qacc bool cg_real_zero<RealDD>(const RealDD& a) {
  return a.Y() == 0.0 && a.X() == 0.0;
}
template<class Real> qacc bool cg_real_finite(const Real& a) {
  return std::isfinite(a);
}
template<> qacc bool cg_real_finite<RealDD>(const RealDD& a) {
  return std::isfinite(a.Y()) && std::isfinite(a.X());
}
template<class Real> qacc bool cg_real_positive(const Real& a) {
  return Real(0.0) < a;
}
template<class Real> qacc bool cg_real_leq(const Real& a, const Real& b) {
  return !(b < a);
}

// Do not use qnorm or ComplexT<RealDD>'s fused helpers here: qnorm returns
// double and some installed complex helpers use the unnormalized DD addition.
template<class Real> qacc ComplexT<Real> cg_complex_add(
    const ComplexT<Real>& a, const ComplexT<Real>& b) {
  return ComplexT<Real>(cg_real_add(a.real(), b.real()),
                        cg_real_add(a.imag(), b.imag()));
}
template<class Real> qacc ComplexT<Real> cg_complex_sub(
    const ComplexT<Real>& a, const ComplexT<Real>& b) {
  return ComplexT<Real>(cg_real_sub(a.real(), b.real()),
                        cg_real_sub(a.imag(), b.imag()));
}
template<class Real> qacc ComplexT<Real> cg_complex_scale(
    const Real& a, const ComplexT<Real>& b) {
  return ComplexT<Real>(cg_real_mul(a, b.real()), cg_real_mul(a, b.imag()));
}
template<class Real> qacc ComplexT<Real> cg_complex_mul(
    const ComplexT<Real>& a, const ComplexT<Real>& b) {
  return ComplexT<Real>(
    cg_real_sub(cg_real_mul(a.real(), b.real()), cg_real_mul(a.imag(), b.imag())),
    cg_real_add(cg_real_mul(a.real(), b.imag()), cg_real_mul(a.imag(), b.real())));
}
template<class Real> qacc ComplexT<Real> cg_complex_conj_mul(
    const ComplexT<Real>& a, const ComplexT<Real>& b) {
  return ComplexT<Real>(
    cg_real_add(cg_real_mul(a.real(), b.real()), cg_real_mul(a.imag(), b.imag())),
    cg_real_sub(cg_real_mul(a.real(), b.imag()), cg_real_mul(a.imag(), b.real())));
}
template<class Real> qacc bool cg_complex_finite(const ComplexT<Real>& a) {
  return cg_real_finite(a.real()) && cg_real_finite(a.imag());
}
template<class Real> qacc Real cg_real_dot_term(
    const ComplexT<Real>& a, const ComplexT<Real>& b) {
  return cg_real_add(cg_real_mul(a.real(), b.real()),
                     cg_real_mul(a.imag(), b.imag()));
}

#ifdef QLAT_USE_ACC
template<class Real> __global__ void cg_typed_dot_partial(
    const ComplexT<Real>* a, const ComplexT<Real>* b, Real* partial, Long length) {
  __shared__ Real sums[256];
  Real sum = Real(0.0);
  for (Long i = Long(blockIdx.x) * 256 + threadIdx.x;
       i < length; i += Long(gridDim.x) * 256)
    sum = cg_real_add(sum, cg_real_dot_term(a[i], b[i]));
  sums[threadIdx.x] = sum;
  __syncthreads();
  for (Int step = 128; step; step >>= 1) {
    if (threadIdx.x < step)
      sums[threadIdx.x] = cg_real_add(sums[threadIdx.x], sums[threadIdx.x + step]);
    __syncthreads();
  }
  if (threadIdx.x == 0) partial[blockIdx.x] = sums[0];
}
template<class Real> __global__ void cg_typed_sum_partial(
    const Real* partial, Real* result, Int length) {
  __shared__ Real sums[256];
  Real sum = Real(0.0);
  for (Int i = threadIdx.x; i < length; i += 256)
    sum = cg_real_add(sum, partial[i]);
  sums[threadIdx.x] = sum;
  __syncthreads();
  for (Int step = 128; step; step >>= 1) {
    if (threadIdx.x < step)
      sums[threadIdx.x] = cg_real_add(sums[threadIdx.x], sums[threadIdx.x + step]);
    __syncthreads();
  }
  if (threadIdx.x == 0) result[0] = sums[0];
}
#endif

// A local, single-participating-GPU reduction. Both GPU trees and the returned
// host coefficient retain Real (including RealDD low components and float
// accumulation). The simple staggered backend performs its rank guard before
// using this helper. No MPI reduction or conversion to double occurs here.
template<class Real> struct CGTypedReduction {
  static_assert(std::is_same_v<Real, float> || std::is_same_v<Real, double> ||
                std::is_same_v<Real, RealDD>);
  vector_gpu<Real> partial, result;
  Int blocks = 0;

  void prepare(Long length) {
    Qassert(length >= 0);
    blocks = Int(std::min<Long>(256, std::max<Long>(1, (length + 255) / 256)));
    partial.resize(blocks, 1);
    result.resize(1, 1);
  }

  Real dot(const ComplexT<Real>* a, const ComplexT<Real>* b, Long length) {
    Qassert(length >= 0);
    if (length == 0) return Real(0.0);
    const Int needed = Int(std::min<Long>(256, (length + 255) / 256));
    if (blocks != needed || partial.size() != size_t(needed) || result.size() != 1)
      prepare(length);
#ifdef QLAT_USE_ACC
    cg_typed_dot_partial<Real><<<blocks, 256>>>(a, b, partial.data(), length);
    Qassert(qacc_GetLastError() == qacc_Success);
    cg_typed_sum_partial<Real><<<1, 256>>>(partial.data(), result.data(), blocks);
    Qassert(qacc_GetLastError() == qacc_Success);
    Real host;
    cpy_GPU(&host, result.data(), 1, 0, 1);
    return host;
#else
    Real sum = Real(0.0);
    for (Long i = 0; i < length; ++i)
      sum = cg_real_add(sum, cg_real_dot_term(a[i], b[i]));
    return sum;
#endif
  }
};

#ifdef QLAT_USE_ACC
// Separate kernels keep the established single-RHS reduction unchanged.
template<class Real> __global__ void cg_typed_batch_dot_partial(
    const ComplexT<Real>* a, const ComplexT<Real>* b, Real* partial, Long length) {
  __shared__ Real sums[256];
  const Long base=Long(blockIdx.y)*length;
  Real sum=Real(0.0);
  for(Long i=Long(blockIdx.x)*256+threadIdx.x;
      i<length;i+=Long(gridDim.x)*256)
    sum=cg_real_add(sum,cg_real_dot_term(a[base+i],b[base+i]));
  sums[threadIdx.x]=sum;
  __syncthreads();
  for(Int step=128;step;step>>=1) {
    if(threadIdx.x<step)
      sums[threadIdx.x]=cg_real_add(sums[threadIdx.x],sums[threadIdx.x+step]);
    __syncthreads();
  }
  if(threadIdx.x==0)
    partial[Long(blockIdx.y)*gridDim.x+blockIdx.x]=sums[0];
}
template<class Real> __global__ void cg_typed_batch_sum_partial(
    const Real* partial, Real* result, Int length) {
  __shared__ Real sums[256];
  const Long base=Long(blockIdx.y)*length;
  Real sum=Real(0.0);
  for(Int i=threadIdx.x;i<length;i+=256)
    sum=cg_real_add(sum,partial[base+i]);
  sums[threadIdx.x]=sum;
  __syncthreads();
  for(Int step=128;step;step>>=1) {
    if(threadIdx.x<step)
      sums[threadIdx.x]=cg_real_add(sums[threadIdx.x],sums[threadIdx.x+step]);
    __syncthreads();
  }
  if(threadIdx.x==0) result[blockIdx.y]=sums[0];
}
#endif

// One typed result per contiguous RHS. Both GPU stages use grid.y for RHSs;
// the arithmetic and the host transfer retain RealDD trailing components.
template<class Real> struct CGTypedBatchReduction {
  static_assert(std::is_same_v<Real,float> || std::is_same_v<Real,double> ||
                std::is_same_v<Real,RealDD>);
  vector_gpu<Real> partial,result;
  Int blocks=0,batch_rhs=0;

  void prepare(Long length,Int rhs) {
    Qassert(length>=0 && rhs>0);
    blocks=Int(std::min<Long>(256,std::max<Long>(1,(length+255)/256)));
    batch_rhs=rhs;
    partial.resize(Long(blocks)*rhs,1);
    result.resize(rhs,1);
  }

  std::vector<Real> dot(const ComplexT<Real>* a,const ComplexT<Real>* b,
                        Long length,Int rhs) {
    Qassert(length>=0 && rhs>0);
    std::vector<Real> host(rhs,Real(0.0));
    if(length==0) return host;
    const Int needed=Int(std::min<Long>(256,(length+255)/256));
    if(blocks!=needed || batch_rhs!=rhs ||
       partial.size()!=size_t(Long(needed)*rhs) || result.size()!=size_t(rhs))
      prepare(length,rhs);
#ifdef QLAT_USE_ACC
    cg_typed_batch_dot_partial<Real><<<dim3(blocks,rhs),256>>>(
      a,b,partial.data(),length);
    Qassert(qacc_GetLastError()==qacc_Success);
    cg_typed_batch_sum_partial<Real><<<dim3(1,rhs),256>>>(
      partial.data(),result.data(),blocks);
    Qassert(qacc_GetLastError()==qacc_Success);
    cpy_GPU(host.data(),result.data(),rhs,0,1);
#else
    for(Int h=0;h<rhs;++h) {
      const Long base=Long(h)*length;
      Real sum=Real(0.0);
      for(Long i=0;i<length;++i)
        sum=cg_real_add(sum,cg_real_dot_term(a[base+i],b[base+i]));
      host[h]=sum;
    }
#endif
    return host;
  }
};

}  // namespace qlat
#endif
