// Optional reduction arithmetic for the custom staggered multi-shift solver.
#ifndef UTILS_MCG_REDUCTION_H
#define UTILS_MCG_REDUCTION_H
#pragma once

// Custom-backend environment controls, both read/validated once per process:
//   qlat_mcg_reduction_precision=0: double; =1: RealDD products and sums.
//   qlat_mcg_color_sum=0: component sums; =1: scalar color-first site sums.
// Both accept only 0/1. Defaults retain the ordinary double reduction path;
// optional buffers are allocated lazily. Enabling these changes arithmetic
// and cost: RealDD expands operations/MPI work; color-first adds site scratch
// traffic and a synchronized pass. Results enter the existing recurrence as
// double. No setting promotes vectors or the QUDA operator. Generic solver
// callers use the default path unless their matrix opts into these controls.
// See utils_staggered_multishift_mcg.h for the full flag/benchmark reference.

#include "general_funs.h"
#include "utils_mcg_debug.h"
#include <algorithm>
#include <stdexcept>
#include <type_traits>

namespace qlat {
// The custom MCG reduction controls are read once per process.  This keeps the
// reduction path stable for the lifetime of a solver and avoids changing the
// arithmetic halfway through an HMC trajectory.  Precision 0 is the historic
// double path and 1 uses a scalar double-double (RealDD) accumulator.
struct MSCGReductionConfig {
  Int precision = 0;
  bool color_sum = false;
};

inline const MSCGReductionConfig& get_mscg_reduction_config() {
  static const MSCGReductionConfig config = [] {
    MSCGReductionConfig value;
    const Long precision = qlat::get_env_long_default(std::string("qlat_mcg_reduction_precision"), 0);
    const Long color = qlat::get_env_long_default(std::string("qlat_mcg_color_sum"), 0);
    if (precision < 0 || precision > 1)
      throw std::invalid_argument("qlat_mcg_reduction_precision must be 0 (double) or 1 (RealDD)");
    if (color < 0 || color > 1)
      throw std::invalid_argument("qlat_mcg_color_sum must be 0 or 1");
    value.precision = Int(precision);
    value.color_sum = color != 0;
    return value;
  }();
  return config;
}

// Accumulator arithmetic; preserve product and summation precision
qacc double mscg_acc_add(double a, double b) {
  return a + b;
}

// Keep both components normalized at each reduction node, on host and device.
// Use an explicit TwoSum rather than relying on host/device-specific RealDD
// operator+= implementations in the external Qlattice installation.
qacc void mscg_two_sum(double a, double b, double& sum, double& error) {
  sum = a + b;
  const double delta = sum - a;
  error = (a - (sum - delta)) + (b - delta);
}

qacc RealDD mscg_acc_add(const RealDD& a, const RealDD& b) {
  double high, low, tail, tail_error;
  mscg_two_sum(a.Y(), b.Y(), high, low);
  mscg_two_sum(a.X(), b.X(), tail, tail_error);
  low += tail;
  mscg_two_sum(high, low, tail, low);
  low += tail_error;
  RealDD result;
  mscg_two_sum(tail, low, result.Y(), result.X());
  return result;
}

template <class Acc, class T> qacc Acc mscg_dot_term(const T& a, const T& b) {
  const Acc ar = Acc(a.real());
  const Acc ai = Acc(a.imag());
  const Acc br = Acc(b.real());
  const Acc bi = Acc(b.imag());
  return mscg_acc_add(ar * br, ai * bi);
}

template <class Acc, class T> qacc Acc mscg_norm_term(const T& a) {
  const Acc ar = Acc(a.real());
  const Acc ai = Acc(a.imag());
  return mscg_acc_add(ar * ar, ai * ai);
}

// QUDA's axpyCGNorm uses the difference of the two rounded residuals,
// not -alpha*Ap. Subtract in vector precision before promoting the products.
template <class Acc, class T, class Real>
qacc void mscg_cg_residual_term(T& residual, const T& ap, Real alpha, Acc& norm, Acc& projection) {
  const T old = residual;
  const T next = old - alpha * ap;
  const T delta = next - old;
  residual = next;
  norm = mscg_norm_term<Acc>(next);
  projection = mscg_dot_term<Acc>(next, delta);
}

// Each entry is Re(conj(a)*b). Promote components before multiplication, so
// RealDD retains product roundoff as well as summation roundoff. ComplexT<RealDD>
// inputs are supported without narrowing their components before the product.
#ifdef QLAT_USE_ACC
template <class Acc, class T>
__global__ void mscg_precise_dot_partial(const T* a, const T* b, Acc* partial, Long length) {
  __shared__ Acc sums[256];
  const Long base = Long(blockIdx.y) * length;
  Acc sum = Acc(0.0);
  for (Long i = Long(blockIdx.x) * 256 + threadIdx.x; i < length; i += Long(gridDim.x) * 256)
    sum = mscg_acc_add(sum, mscg_dot_term<Acc>(a[base + i], b[base + i]));
  sums[threadIdx.x] = sum;
  __syncthreads();
  for (int s = 128; s; s >>= 1) {
    if (threadIdx.x < s)
      sums[threadIdx.x] = mscg_acc_add(sums[threadIdx.x], sums[threadIdx.x + s]);
    __syncthreads();
  }
  if (threadIdx.x == 0)
    partial[blockIdx.y * gridDim.x + blockIdx.x] = sums[0];
}

template <class Acc> __global__ void mscg_sum_partial(const Acc* input, Acc* partial, Long length) {
  __shared__ Acc sums[256];
  const Long base = Long(blockIdx.y) * length;
  Acc sum = Acc(0.0);
  for (Long i = Long(blockIdx.x) * 256 + threadIdx.x; i < length; i += Long(gridDim.x) * 256)
    sum = mscg_acc_add(sum, input[base + i]);
  sums[threadIdx.x] = sum;
  __syncthreads();
  for (int s = 128; s; s >>= 1) {
    if (threadIdx.x < s)
      sums[threadIdx.x] = mscg_acc_add(sums[threadIdx.x], sums[threadIdx.x + s]);
    __syncthreads();
  }
  if (threadIdx.x == 0)
    partial[blockIdx.y * gridDim.x + blockIdx.x] = sums[0];
}

template <class Acc, class T, class Real>
__global__ void mscg_precise_residual_partial(T* residual, const T* ap, const double* alpha, Acc* partial,
                                              Long length) {
  __shared__ Acc sums[256];
  const Int rhs = blockIdx.y;
  const Long base = Long(rhs) * length;
  const Real a = Real(alpha[rhs]);
  Acc sum = Acc(0.0);
  for (Long i = Long(blockIdx.x) * 256 + threadIdx.x; i < length; i += Long(gridDim.x) * 256) {
    const Long index = base + i;
    residual[index] -= a * ap[index];
    sum = mscg_acc_add(sum, mscg_norm_term<Acc>(residual[index]));
  }
  sums[threadIdx.x] = sum;
  __syncthreads();
  for (int s = 128; s; s >>= 1) {
    if (threadIdx.x < s)
      sums[threadIdx.x] = mscg_acc_add(sums[threadIdx.x], sums[threadIdx.x + s]);
    __syncthreads();
  }
  if (threadIdx.x == 0)
    partial[rhs * gridDim.x + blockIdx.x] = sums[0];
}

template <class Acc, class T, class Real>
__global__ void mscg_cg_residual_partial(T* residual, const T* ap, const double* alpha, Acc* partial,
                                         Long length, Int rhs) {
  __shared__ Acc norms[256], projections[256];
  const Int h = blockIdx.y;
  const Long base = Long(h) * length;
  const Real a = Real(alpha[h]);
  Acc norm = Acc(0.0), projection = Acc(0.0);
  for (Long i = Long(blockIdx.x) * 256 + threadIdx.x; i < length; i += Long(gridDim.x) * 256) {
    Acc n, p;
    mscg_cg_residual_term<Acc>(residual[base + i], ap[base + i], a, n, p);
    norm = mscg_acc_add(norm, n);
    projection = mscg_acc_add(projection, p);
  }
  norms[threadIdx.x] = norm;
  projections[threadIdx.x] = projection;
  __syncthreads();
  for (int s = 128; s; s >>= 1) {
    if (threadIdx.x < s) {
      norms[threadIdx.x] = mscg_acc_add(norms[threadIdx.x], norms[threadIdx.x + s]);
      projections[threadIdx.x] = mscg_acc_add(projections[threadIdx.x], projections[threadIdx.x + s]);
    }
    __syncthreads();
  }
  if (threadIdx.x == 0) {
    partial[Long(h) * gridDim.x + blockIdx.x] = norms[0];
    partial[Long(rhs + h) * gridDim.x + blockIdx.x] = projections[0];
  }
}
#endif

inline void mscg_sum_accumulator(RealDD* values, Int count) {
  sum_all_mpi(values, count);
}
inline void mscg_sum_accumulator(double* values, Int count) {
  sum_all_size(values, count, 0);
}

// A reusable reduction with explicitly selected accumulation precision. It is
// separate from MSCGReduction so the old default kernels and other callers keep
// their arithmetic and launch order. Only the final globally summed value is
// narrowed to double for the existing host recurrence.
// Lazy reusable reduction storage and host-facing double results
template <class Acc> struct MSCGTypedReduction {
  static_assert(std::is_same<Acc, double>::value || std::is_same<Acc, RealDD>::value, "Unsupported solver precision");
  MCGDebugBuffer<Acc> partial, result, host, site_values;
  MCGDebugBuffer<double> output;
  Int blocks = 0, mrh = 0;
  Long allocations = 0;

  template <class V> void resize_one(V& v, Long n, Int memory, const char* name) {
    if (v.size() != size_t(n) || v.GPU != memory) {
      v.resize(n, memory);
      v.debug_name(name);
      ++allocations;
    }
  }

  void prepare(Long length, Int rhs, bool color_sum) {
    MCGDebugPhase phase("reduction_prepare");
    Qassert(length > 0 && rhs > 0 && (!color_sum || length % 3 == 0));
    const Long count = color_sum ? length / 3 : length;
    blocks = std::min<Long>(256, (count + 255) / 256);
    mrh = rhs;
    resize_one(partial, Long(blocks) * rhs, 1, "reduction.partial");
    resize_one(result, rhs, 1, "reduction.result");
    resize_one(host, rhs, 0, "reduction.host");
    resize_one(output, rhs, 0, "reduction.output");
    if (color_sum)
      resize_one(site_values, count * rhs, 1, "reduction.site_values");
  }

  const double* finish() {
    MCGDebugPhase phase("reduction_finish");
#ifdef QLAT_USE_ACC
    // Reuse the same tree for stage two, keeping Acc in shared memory too.
    mscg_sum_partial<<<dim3(1, mrh), 256>>>(partial.data(), result.data(), blocks);
    Qassert(qacc_GetLastError() == qacc_Success);
    cpy_GPU(host.data(), result.data(), mrh, 0, 1);
#endif
    {
      MCGDebugPhase mpi_phase("host_mpi_reduction");
      mscg_sum_accumulator(host.data(), mrh);
    }
    for (Int h = 0; h < mrh; ++h)
      output[h] = double(host[h]);
    return output.data();
  }

  void reduce_sites(Long sites) {
#ifdef QLAT_USE_ACC
    mscg_sum_partial<<<dim3(blocks, mrh), 256>>>(site_values.data(), partial.data(), sites);
    Qassert(qacc_GetLastError() == qacc_Success);
#else
    for (Int h = 0; h < mrh; ++h) {
      Acc sum = Acc(0.0);
      for (Long site = 0; site < sites; ++site)
        sum = mscg_acc_add(sum, site_values[h * sites + site]);
      host[h] = sum;
    }
#endif
  }

  template <class T, class Layout>
  const double* dot(const T* a, const T* b, Long length, Int rhs, bool color_sum, const Layout& layout,
                    bool public_vectors = false) {
    prepare(length, rhs, color_sum);
    if (color_sum) {
      const Long sites = length / 3;
      Acc* values = site_values.data();
      // Sum scalar color contributions, never the complex colors themselves.
      qacc_for(i, sites * rhs, {
        const Long site = i % sites, base = (i / sites) * length;
        Acc sum = Acc(0.0);
        for (Int color = 0; color < 3; ++color) {
          const Long index = base + (public_vectors ? layout.template public_index<T>(site, color, sites)
                                                    : layout.storage_index(site, color, sites));
          sum = mscg_acc_add(sum, mscg_dot_term<Acc>(a[index], b[index]));
        }
        values[i] = sum;
      });
      reduce_sites(sites);
    } else {
#ifdef QLAT_USE_ACC
      mscg_precise_dot_partial<Acc><<<dim3(blocks, rhs), 256>>>(a, b, partial.data(), length);
      Qassert(qacc_GetLastError() == qacc_Success);
#else
      for (Int h = 0; h < rhs; ++h) {
        Acc sum = Acc(0.0);
        for (Long i = 0; i < length; ++i)
          sum = mscg_acc_add(sum, mscg_dot_term<Acc>(a[h * length + i], b[h * length + i]));
        host[h] = sum;
      }
#endif
    }
    return finish();
  }

  template <class T, class Real, class Layout>
  const double* residual_norm(T* residual, const T* ap, const double* alpha, Long length, Int rhs,
                              bool color_sum, const Layout& layout) {
    prepare(length, rhs, color_sum);
    if (color_sum) {
      const Long sites = length / 3;
      Acc* values = site_values.data();
      // Residual vectors retain their configured precision and update formula.
      // This pass fuses the update and the three-color norm at each site.
      qacc_for(i, sites * rhs, {
        const Long site = i % sites, h = i / sites, base = h * length;
        const Real a = Real(alpha[h]);
        Acc sum = Acc(0.0);
        for (Int color = 0; color < 3; ++color) {
          const Long index = base + layout.storage_index(site, color, sites);
          residual[index] -= a * ap[index];
          sum = mscg_acc_add(sum, mscg_norm_term<Acc>(residual[index]));
        }
        values[i] = sum;
      });
      reduce_sites(sites);
    } else {
#ifdef QLAT_USE_ACC
      mscg_precise_residual_partial<Acc, T, Real>
          <<<dim3(blocks, rhs), 256>>>(residual, ap, alpha, partial.data(), length);
      Qassert(qacc_GetLastError() == qacc_Success);
#else
      for (Int h = 0; h < rhs; ++h) {
        const Real a = Real(alpha[h]);
        Acc sum = Acc(0.0);
        for (Long i = 0; i < length; ++i) {
          const Long index = h * length + i;
          residual[index] -= a * ap[index];
          sum = mscg_acc_add(sum, mscg_norm_term<Acc>(residual[index]));
        }
        host[h] = sum;
      }
#endif
    }
    return finish();
  }
};

// Two globally reduced channels: [norm per RHS, signed projection per RHS].
// The base reducer retains Acc through both GPU trees and the MPI reduction;
// only its final host-facing output is narrowed to double. No old-vector copy
// or temporary full-volume vector is needed, including with color grouping.
template <class Acc> struct MSCGPairedReduction : MSCGTypedReduction<Acc> {
  template <class T, class Real, class Layout>
  const double* residual_cg_norm(T* residual, const T* ap, const double* alpha, Long length, Int rhs,
                                 bool color_sum, const Layout& layout) {
    this->prepare(length, 2 * rhs, color_sum);
    if (color_sum) {
      const Long sites = length / 3;
      Acc* values = this->site_values.data();
      qacc_for(i, sites * rhs, {
        const Long site = i % sites, h = i / sites, base = h * length;
        const Real a = Real(alpha[h]);
        Acc norm = Acc(0.0), projection = Acc(0.0);
        for (Int color = 0; color < 3; ++color) {
          const Long index = base + layout.storage_index(site, color, sites);
          Acc n, p;
          mscg_cg_residual_term<Acc>(residual[index], ap[index], a, n, p);
          norm = mscg_acc_add(norm, n);
          projection = mscg_acc_add(projection, p);
        }
        values[i] = norm;
        values[Long(rhs) * sites + i] = projection;
      });
      this->reduce_sites(sites);
    } else {
#ifdef QLAT_USE_ACC
      mscg_cg_residual_partial<Acc, T, Real>
          <<<dim3(this->blocks, rhs), 256>>>(residual, ap, alpha, this->partial.data(), length, rhs);
      Qassert(qacc_GetLastError() == qacc_Success);
#else
      for (Int h = 0; h < rhs; ++h) {
        const Real a = Real(alpha[h]);
        Acc norm = Acc(0.0), projection = Acc(0.0);
        for (Long i = 0; i < length; ++i) {
          const Long index = Long(h) * length + i;
          Acc n, p;
          mscg_cg_residual_term<Acc>(residual[index], ap[index], a, n, p);
          norm = mscg_acc_add(norm, n);
          projection = mscg_acc_add(projection, p);
        }
        this->host[h] = norm;
        this->host[rhs + h] = projection;
      }
#endif
    }
    return this->finish();
  }
};
} // namespace qlat
#endif
