// Multi-shift CG: selectable public-vector and Krylov precision.
#ifndef UTILS_MULTI_SHIFT_CG_H
#define UTILS_MULTI_SHIFT_CG_H
#pragma once

// Tuning and diagnostics are configured by MSCGOptions below. Defaults keep
// fused residual/norm updates, block 128, shift grouping 6, temporal window
// max(1, shifts - 1), scoped qacc threads 128/64, and profiling/stopping observers disabled.
// True-residual evaluation is optional and defaults to off. Larger temporal
// windows allocate history lazily and retain it for reuse. Set temporal_batch
// positive to override the automatic window; finished lanes leave replay after
// their final pending update. D/F public
// calls default to float solution updates backed by caller output storage;
// public doubles are exported mass-by-mass after the solve. History copies
// and replay launches stop when no nonseed updates remain. Small scheduling
// buffers are cached with the vectors; fully active windows retain the dense
// kernel. Host-only execution resolves automatic to 1. See
// utils_staggered_multishift_mcg.h for environment controls,
// benchmark input keys, precision-specific evidence and observation costs.

#include "general_funs.h"
#include <qlat-utils/qacc.h>
#include "utils_geo_cache.h"
#include "utils_mcg_reduction.h"
#include <memory>
#include <limits>
#include <type_traits>
#include <chrono>
#include <map>
#include <optional>

namespace qlat {
// Host-only observers. Iterations are one-based; -1 means no stopping poll.
// Host observations and caller options
enum class MSCGStopReason : char {
  target = 'T',
  epsilon = 'E',
  zero_source = 'Z',
  exact_zero = 'X',
  initial_target = 'I',
  fixed_limit = 'F',
  iteration_limit = 'L'
};

inline const char* mscg_stop_reason_name(MSCGStopReason reason) {
  switch (reason) {
  case MSCGStopReason::target:
    return "target";
  case MSCGStopReason::epsilon:
    return "epsilon";
  case MSCGStopReason::zero_source:
    return "zero-source";
  case MSCGStopReason::exact_zero:
    return "exact-zero";
  case MSCGStopReason::initial_target:
    return "initial-target";
  case MSCGStopReason::fixed_limit:
    return "fixed-limit";
  case MSCGStopReason::iteration_limit:
    return "iteration-limit";
  }
  return "unknown";
}

struct MSCGStoppingChecks {
  Int previous_iteration = -1, iteration = -1;
  double previous_ratio = std::numeric_limits<double>::quiet_NaN();
  double ratio = std::numeric_limits<double>::quiet_NaN();

  void record(Int step, double residual2, double threshold2) {
    previous_iteration = iteration;
    previous_ratio = ratio;
    iteration = step;
    // Preserve a zero/underflowed threshold as NaN or Inf in the diagnostic.
    ratio = threshold2 > 0  ? residual2 / threshold2
            : residual2 > 0 ? std::numeric_limits<double>::infinity()
                            : std::numeric_limits<double>::quiet_NaN();
  }
};

struct MSCGStoppingDiagnostics {
  Int maxiter = 0;
  std::vector<double> source_norm2, seed_threshold2, shift_threshold2;
  std::vector<MSCGStoppingChecks> seed_checks, shift_checks;
  std::vector<MSCGStopReason> seed_reason, shift_reason;
  std::vector<char> shift_retired;
  // Shift-major/RHS-major. The seed entry counts active recurrence updates;
  // higher entries count nonzero coefficients in the vector-update precision.
  std::vector<Int> effective_updates;
};

struct MSCGStats {
  Int iterations = 0, refine_iterations = 0, matrix_calls = 0;
  Long reductions = 0, seed_updates = 0, shift_updates = 0, residual_checks = 0;
  Long workspace_allocations = 0, matrix_allocations = 0;
  std::vector<Int> shift_iterations;
  std::vector<double> true_residual, iterative_residual;
  std::map<std::string, double> host_seconds, device_seconds;
  // Local conventional arithmetic counts (not hardware instruction counts).
  // The staggered adapter supplies operator FLOPs and synchronized wall time.
  double modeled_vector_flops = 0, modeled_operator_flops = 0, solve_seconds = 0;
  bool borrowed_solutions = false;
  Int accumulation_bytes = 0;
  double retirement_epsilon = 0;
  MSCGStoppingDiagnostics stopping;
};
// Inherit is resolved by the staggered entry points. Generic matrix callers
// treat it as requested-target retirement and the ordinary norm-ratio beta.
enum class MSCGShiftStopPolicy { inherit, requested, cpu_epsilon, precision_epsilon };
enum class MSCGBetaPolicy { inherit, norm_ratio, quda };
enum class MSCGAccumulationPolicy { inherit, single, double_precision };

struct MSCGOptions {
  bool fixed_iterations = false, profile = false;
  // Unset inherits the staggered process policy, or false for generic matrices.
  // Explicit true/false overrides the process policy.
  std::optional<bool> compute_true_res;
  double refine_tolerance = 0;
  Int refine_maxiter = -1;
  MSCGStats* stats = nullptr;
  // Diagnostic/tuning controls; precision still comes exclusively from qinv.
  bool generic_vector_layout = false;
  bool fuse_residual_norm = true;
  Int vector_block = 128, shift_batch = 6;
  Int temporal_batch = 0; // 0 selects max(1, shifts - 1).
  Int qacc_threads_single = 128, qacc_threads_double = 64;
  bool stopping_diagnostics = false;
  MSCGShiftStopPolicy shift_stop = MSCGShiftStopPolicy::inherit;
  MSCGBetaPolicy beta = MSCGBetaPolicy::inherit;
  MSCGAccumulationPolicy accumulation = MSCGAccumulationPolicy::inherit;
  bool print_mem = false; // Enabled by qlat_mcg_print_mem for custom staggered solves.
};

// Scoped launch settings and reduction primitives
class ScopedQaccNumThreads {
  Int& setting;
  Int previous;

public:
  explicit ScopedQaccNumThreads(Int threads) : setting(qacc_num_threads()), previous(setting) {
    Qassert(threads > 0 && threads <= 1024);
    setting = threads;
  }

  ScopedQaccNumThreads(const ScopedQaccNumThreads&) = delete;
  ScopedQaccNumThreads& operator=(const ScopedQaccNumThreads&) = delete;

  ~ScopedQaccNumThreads() {
    setting = previous;
  }
};

// Local, reusable two-stage reductions. Float components are promoted BEFORE
// products; neither stage uses float sums or QUDA BLAS. One grid-y row per RHS.
#ifdef QLAT_USE_ACC
template <class T> __global__ void mscg_dot_partial(const T* a, const T* b, double* partial, Long length) {
  __shared__ double sums[256];
  const Long base = Long(blockIdx.y) * length;
  double sum = 0;
  for (Long i = Long(blockIdx.x) * 256 + threadIdx.x; i < length; i += Long(gridDim.x) * 256) {
    const double ar = double(a[base + i].real()), ai = double(a[base + i].imag());
    const double br = double(b[base + i].real()), bi = double(b[base + i].imag());
    sum += ar * br + ai * bi;
  }
  sums[threadIdx.x] = sum;
  __syncthreads();
  for (int s = 128; s; s >>= 1) {
    if (threadIdx.x < s)
      sums[threadIdx.x] += sums[threadIdx.x + s];
    __syncthreads();
  }
  if (threadIdx.x == 0)
    partial[blockIdx.y * gridDim.x + blockIdx.x] = sums[0];
}

template <class T, class Real>
__global__ void mscg_residual_norm_partial(T* residual, const T* ap, const double* alpha, double* partial,
                                           Long length) {
  __shared__ double sums[256];
  const Int rhs = blockIdx.y;
  const Long base = Long(rhs) * length;
  const Real a = Real(alpha[rhs]);
  double sum = 0;
  for (Long i = Long(blockIdx.x) * 256 + threadIdx.x; i < length; i += Long(gridDim.x) * 256) {
    const Long index = base + i;
    residual[index] -= a * ap[index];
    const double re = double(residual[index].real()), im = double(residual[index].imag());
    sum += re * re + im * im;
  }
  sums[threadIdx.x] = sum;
  __syncthreads();
  for (int s = 128; s; s >>= 1) {
    if (threadIdx.x < s)
      sums[threadIdx.x] += sums[threadIdx.x + s];
    __syncthreads();
  }
  if (threadIdx.x == 0)
    partial[rhs * gridDim.x + blockIdx.x] = sums[0];
}

static __global__ void mscg_dot_finish(const double* partial, double* result, Int blocks) {
  __shared__ double sums[256];
  double sum = 0;
  for (Int i = threadIdx.x; i < blocks; i += 256)
    sum += partial[blockIdx.x * blocks + i];
  sums[threadIdx.x] = sum;
  __syncthreads();
  for (int s = 128; s; s >>= 1) {
    if (threadIdx.x < s)
      sums[threadIdx.x] += sums[threadIdx.x + s];
    __syncthreads();
  }
  if (threadIdx.x == 0)
    result[blockIdx.x] = sums[0];
}
#endif
struct MSCGReduction {
  MCGDebugBuffer<double> partial, result, host;
  Int blocks = 0, mrh = 0;
  Long allocations = 0;

  void resize(Long length, Int rhs) {
    const Int nb = std::min<Long>(256, (length + 255) / 256);
    if (nb == blocks && rhs == mrh)
      return;
    blocks = nb;
    mrh = rhs;
    partial.resize(Long(nb) * rhs);
    result.resize(rhs);
    host.resize(rhs, 0);
    partial.debug_name("reduction.partial");
    result.debug_name("reduction.result");
    host.debug_name("reduction.host");
    allocations += 3;
  }

  const double* finish(Int rhs) {
    MCGDebugPhase phase("reduction_finish");
#ifdef QLAT_USE_ACC
    mscg_dot_finish<<<rhs, 256>>>(partial.data(), result.data(), blocks);
    Qassert(qacc_GetLastError() == qacc_Success);
    cpy_GPU(host.data(), result.data(), rhs, 0, 1);
#endif
    {
      MCGDebugPhase mpi_phase("host_mpi_reduction");
      sum_all_size(host.data(), rhs, 0);
    }
    return host.data();
  }

  template <class T> const double* dot(const T* a, const T* b, Long length, Int rhs) {
    resize(length, rhs);
#ifdef QLAT_USE_ACC
    mscg_dot_partial<<<dim3(blocks, rhs), 256>>>(a, b, partial.data(), length);
    Qassert(qacc_GetLastError() == qacc_Success);
#else
    for (Int h = 0; h < rhs; ++h) {
      double sum = 0;
      for (Long i = 0; i < length; ++i) {
        const T aa = a[h * length + i], bb = b[h * length + i];
        sum += double(aa.real()) * double(bb.real()) + double(aa.imag()) * double(bb.imag());
      }
      host[h] = sum;
    }
#endif
    return finish(rhs);
  }

  template <class T, class Real>
  const double* residual_norm(T* residual, const T* ap, const double* alpha, Long length, Int rhs) {
    resize(length, rhs);
#ifdef QLAT_USE_ACC
    mscg_residual_norm_partial<T, Real>
        <<<dim3(blocks, rhs), 256>>>(residual, ap, alpha, partial.data(), length);
    Qassert(qacc_GetLastError() == qacc_Success);
#else
    for (Int h = 0; h < rhs; ++h) {
      const Real a = Real(alpha[h]);
      double sum = 0;
      for (Long i = 0; i < length; ++i) {
        residual[h * length + i] -= a * ap[h * length + i];
        const double re = double(residual[h * length + i].real()),
                     im = double(residual[h * length + i].imag());
        sum += re * re + im * im;
      }
      host[h] = sum;
    }
#endif
    return finish(rhs);
  }
};

// Profiling is opt-in. Events are retained across calls; no event creation in
// warmed profiled solves. Headline timings must run with profiling disabled.
// Optional phase profiling
struct MSCGProfile {
  struct Record {
    std::string name;
#ifdef QLAT_USE_ACC
    qacc_Event_t begin = nullptr, end = nullptr;
#endif
  };

  std::vector<Record> records;
  size_t used = 0;
  bool enabled = false;
  std::map<std::string, double> host, device;

  ~MSCGProfile() {
#ifdef QLAT_USE_ACC
    for (auto& r : records) {
      qacc_Errchk(qacc_EventDestroy(r.begin));
      qacc_Errchk(qacc_EventDestroy(r.end));
    }
#endif
  }

  void reset(bool on) {
    enabled = on;
    used = 0;
    host.clear();
    device.clear();
  }

  template <class F> void run(const char* name, F&& f) {
    MCGDebugPhase debug_phase(name);
    if (!enabled) {
      f();
      return;
    }
    if (used == records.size()) {
      records.emplace_back();
#ifdef QLAT_USE_ACC
      Qassert(qacc_EventCreate(&records.back().begin) == qacc_Success);
      Qassert(qacc_EventCreate(&records.back().end) == qacc_Success);
#endif
    }
    const size_t index = used++;
    records[index].name = name;
#ifdef QLAT_USE_ACC
    qacc_Errchk(qacc_EventRecord(records[index].begin));
#endif
    auto start = std::chrono::steady_clock::now();
    f();
    host[name] += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
#ifdef QLAT_USE_ACC
    qacc_Errchk(qacc_EventRecord(records[index].end));
#endif
  }

  void finish() {
#ifdef QLAT_USE_ACC
    if (enabled && used) {
      qacc_Errchk(qacc_EventSynchronize(records[used - 1].end));
      for (size_t i = 0; i < used; ++i) {
        float ms = 0;
        qacc_Errchk(qacc_EventElapsedTime(&ms, records[i].begin, records[i].end));
        device[records[i].name] += ms * 1e-3;
      }
    }
#endif
  }
};

// One group of shifted vectors with unapplied updates. Masks cover at most
// shift_batch lanes, not the total shift count. steps excludes unused tail slots.
struct MSCGShiftWork {
  Int rhs = 0, first = 0, mask = 0, steps = 0;
};

// Layout policies map a logical site/color separately for solution and Krylov
// storage. Non-QUDA matrix callers retain an identical contiguous layout.
// Vector layout policies and CUDA update kernels
struct MSCGFlatLayout {
  qacc Long double_index(Long site, Int color, Long) const {
    return 3 * site + color;
  }

  qacc Long storage_index(Long site, Int color, Long) const {
    return 3 * site + color;
  }

  template <class U> qacc Long public_index(Long site, Int color, Long sites) const {
    return double_index(site, color, sites);
  }
};

template <int Order> qacc Long mscg_native_index(Long site, Int color, Long sites) {
  static_assert(Order >= 2 && Order % 2 == 0, "Native complex packing must be positive and even");
  constexpr Int width = Order / 2, full = 3 / width, remainder = 3 - full * width;
  const Int group = color / width, lane = color % width;
  const Int pack = (remainder == 0 || group < full) ? width : remainder;
  return Long(group) * sites * width + site * pack + lane;
}

// Staggered search vectors retain Krylov arithmetic independently of the
// solution accumulator. Generic layouts retain their established arithmetic.
template <class T, class U, class Layout, class Scalar>
qacc T mscg_search_update(T residual, T old, Scalar z, Scalar beta) {
  if constexpr (requires { Layout::krylov_arithmetic; }) {
    using Real = decltype(T{}.real());
    return Real(z) * residual + Real(beta) * old;
  } else {
    return T(z * U(residual) + beta * U(old));
  }
}

template <bool Borrowed, class U>
qacc U* mscg_solution_row(U* const* rows, Int shift, Int rhs, Int mrh, Long length) {
  if constexpr (Borrowed)
    return rows[Long(shift) * mrh + rhs];
  else
    return rows[shift] + Long(rhs) * length;
}
#ifdef QLAT_USE_ACC
template <class T, class U, class Layout>
__global__ void mscg_initialize_vectors(const U* b, T* r, T* p, T* const* searches, Int shifts, Long sites,
                                        Layout layout) {
  const Long site = Long(blockIdx.x) * blockDim.x + threadIdx.x;
  if (site >= sites)
    return;
  const Long base = Long(blockIdx.y) * sites * 3;
#pragma unroll
  for (Int color = 0; color < 3; ++color) {
    const Long it = base + layout.storage_index(site, color, sites);
    const T value = T(b[base + layout.template public_index<U>(site, color, sites)]);
    r[it] = value;
    p[it] = value;
    for (Int j = 1; j < shifts; ++j)
      searches[j][it] = value;
  }
}

template <class T, class U, bool Borrowed = false, class Layout>
__global__ void mscg_seed_vectors(U* x, T* p, const T* r, const double* c, Int mrh, Long sites,
                                  Layout layout, U* const* rows = nullptr) {
  const Long site = Long(blockIdx.x) * blockDim.x + threadIdx.x;
  if (site >= sites)
    return;
  const Int h = blockIdx.y;
  const Long base = Long(h) * sites * 3;
  using Scalar = decltype(U{}.real());
  const Scalar alpha = Scalar(c[h]);
  const Scalar beta = Scalar(c[mrh + h]);
  if constexpr (Borrowed)
    x = rows[h];
#pragma unroll
  for (Int color = 0; color < 3; ++color) {
    const Long it = base + layout.storage_index(site, color, sites);
    const Long id = (Borrowed ? 0 : base) + layout.template public_index<U>(site, color, sites);
    const U old = U(p[it]);
    x[id] += U(alpha * old);
    p[it] = mscg_search_update<T, U, Layout>(r[it], p[it], Scalar(1), beta);
  }
}

template <int Batch, class T, class U, bool Borrowed = false, class Layout>
__global__ void mscg_shift_vectors(U* const* solutions, T* const* searches, const T* r, const double* c,
                                   Int shifts, Int mrh, Long sites, Layout layout) {
  using Scalar = decltype(U{}.real());
  const Long site = Long(blockIdx.x) * blockDim.x + threadIdx.x;
  if (site >= sites)
    return;
  const Int h = blockIdx.y, first = 1 + blockIdx.z * Batch;
  const Long base = Long(h) * sites * 3, S = Long(shifts) * mrh;
#pragma unroll
  for (Int color = 0; color < 3; ++color) {
    const Long it = base + layout.storage_index(site, color, sites);
    const Long id = layout.template public_index<U>(site, color, sites);
    const U residual = U(r[it]);
#pragma unroll
    for (Int lane = 0; lane < Batch; ++lane) {
      const Int j = first + lane;
      if (j < shifts) {
        const Long v = 2 * mrh + Long(j) * mrh + h;
        const Scalar residual_coefficient = Scalar(c[v]);
        if (residual_coefficient != Scalar(0)) {
          const Scalar solution_coefficient = Scalar(c[v + S]);
          const Scalar search_coefficient = Scalar(c[v + 2 * S]);
          const U old = U(searches[j][it]);
          mscg_solution_row<Borrowed>(solutions, j, h, mrh, 3 * sites)[id] += U(solution_coefficient * old);
          searches[j][it] = mscg_search_update<T, U, Layout>(T(residual), searches[j][it],
                                                           residual_coefficient, search_coefficient);
        }
      }
    }
  }
}

#endif  // QLAT_USE_ACC

// The temporal coefficient history must use the arithmetic precision of the
// solution vectors.  In the all-single path, keeping this history in double
// silently promotes every shifted update to FP64. Mixed D/F solves use float
// history for float accumulation and double history for double accumulation.
template <class U> struct mscg_history_scalar;

template <> struct mscg_history_scalar<ComplexF> {
  using type = float;
};

template <> struct mscg_history_scalar<ComplexD> {
  using type = double;
};
template <class U> using mscg_history_scalar_t = typename mscg_history_scalar<U>::type;

#ifdef QLAT_USE_ACC
template <class U>
__global__ void mscg_convert_coeff_history(const double* src, mscg_history_scalar_t<U>* dst, Long n) {
  for (Long i = Long(blockIdx.x) * blockDim.x + threadIdx.x; i < n; i += Long(gridDim.x) * blockDim.x)
    dst[i] = mscg_history_scalar_t<U>(src[i]);
}
// Window 0 consumes the runtime step count; 4/8/16 retain bounded specializations.
template <int Window, int Batch, class T, class U, bool Borrowed = false, class Layout>
__global__ void
mscg_shift_vectors_temporal(U* const* solutions, T* const* searches, const T* residual_history,
                            const mscg_history_scalar_t<U>* coefficient_history, Int steps, Int shifts,
                            Int mrh, Long sites, Long coeff_stride, Layout layout) {
  const Long site = Long(blockIdx.x) * blockDim.x + threadIdx.x;
  if (site >= sites)
    return;
  const Int h = blockIdx.y, first = 1 + blockIdx.z * Batch;
  const Long L = 3 * sites, base = Long(h) * L, S = Long(shifts) * mrh;
  for (Int color = 0; color < 3; ++color) {
    const Long it = base + layout.storage_index(site, color, sites);
    const Long id = layout.template public_index<U>(site, color, sites);
    // Keep each shifted solution/search value live across the temporal window.
    // This changes the dominant mass-vector traffic from one read/write pair
    // per iteration to one load and one store per window.  The lane loop is
    // compile-time unrolled so these small arrays become registers rather than
    // dynamically indexed local memory.
    U x_reg[Batch];
    T p_reg[Batch];
#pragma unroll
    for (Int lane = 0; lane < Batch; ++lane) {
      const Int j = first + lane;
      if (j < shifts) {
        x_reg[lane] = mscg_solution_row<Borrowed>(solutions, j, h, mrh, L)[id];
        p_reg[lane] = searches[j][it];
      } else {
        x_reg[lane] = U(0.0);
        p_reg[lane] = T(0.0);
      }
    }
    for (Int t = 0; t < steps && (Window == 0 || t < Window); ++t) {
      const U residual = U(residual_history[Long(t) * L * mrh + it]);
      const mscg_history_scalar_t<U>* c = coefficient_history + Long(t) * coeff_stride;
#pragma unroll
      for (Int lane = 0; lane < Batch; ++lane) {
        const Int j = first + lane;
        if (j < shifts) {
          const Long v = 2 * mrh + Long(j) * mrh + h;
          if (c[v] != 0) {
            const U old = U(p_reg[lane]);
            x_reg[lane] += U(c[v + S] * old);
            p_reg[lane] = mscg_search_update<T, U, Layout>(T(residual), p_reg[lane], c[v], c[v + 2 * S]);
          }
        }
      }
    }
#pragma unroll
    for (Int lane = 0; lane < Batch; ++lane) {
      const Int j = first + lane;
      if (j < shifts) {
        mscg_solution_row<Borrowed>(solutions, j, h, mrh, L)[id] = x_reg[lane];
        searches[j][it] = p_reg[lane];
      }
    }
  }
}

// Sparse replay is used only after a lane/RHS becomes inactive. Test the
// coefficients before reading history: unused RHS slots are deliberately not
// copied. In particular, a retired lane never loads/stores its vectors again
// once its last pending update has been applied.
template <bool Temporal, int Batch, class T, class U, bool Borrowed = false, class Layout>
__global__ void mscg_shift_vectors_pending(
    U* const* solutions, T* const* searches, const T* residuals,
    const mscg_history_scalar_t<U>* history, const double* coefficients,
    const MSCGShiftWork* work, Int shifts, Int mrh, Long sites, Long stride, Layout layout) {
  const Long site = Long(blockIdx.x) * blockDim.x + threadIdx.x;
  if (site >= sites)
    return;
  const auto entry = work[blockIdx.y];
  const Int h = entry.rhs, first = entry.first;
  const Long L = 3 * sites, base = Long(h) * L, S = Long(shifts) * mrh;
  using Scalar = mscg_history_scalar_t<U>;
  for (Int color = 0; color < 3; ++color) {
    const Long it = base + layout.storage_index(site, color, sites);
    const Long id = layout.template public_index<U>(site, color, sites);
    U x_reg[Batch];
    T p_reg[Batch];
#pragma unroll
    for (Int lane = 0; lane < Batch; ++lane) {
      x_reg[lane] = U(0.0);
      p_reg[lane] = T(0.0);
      if (entry.mask & (1 << lane)) {
        x_reg[lane] = mscg_solution_row<Borrowed>(solutions, first + lane, h, mrh, L)[id];
        p_reg[lane] = searches[first + lane][it];
      }
    }
    for (Int t = 0; t < entry.steps; ++t) {
      Scalar a[Batch], b[Batch], beta[Batch];
      bool needed = false;
#pragma unroll
      for (Int lane = 0; lane < Batch; ++lane) {
        a[lane] = 0;
        if (entry.mask & (1 << lane)) {
          const Long v = 2 * mrh + Long(first + lane) * mrh + h;
          if constexpr (Temporal) {
            const auto* c = history + Long(t) * stride;
            a[lane] = c[v];
            b[lane] = c[v + S];
            beta[lane] = c[v + 2 * S];
          } else {
            a[lane] = Scalar(coefficients[v]);
            b[lane] = Scalar(coefficients[v + S]);
            beta[lane] = Scalar(coefficients[v + 2 * S]);
          }
          needed = needed || a[lane] != Scalar(0);
        }
      }
      if (!needed)
        continue;
      const U residual = U(residuals[(Temporal ? Long(t) * L * mrh : 0) + it]);
#pragma unroll
      for (Int lane = 0; lane < Batch; ++lane) {
        if (a[lane] != Scalar(0)) {
          const U old = U(p_reg[lane]);
          x_reg[lane] += U(b[lane] * old);
          p_reg[lane] = mscg_search_update<T, U, Layout>(T(residual), p_reg[lane], a[lane], beta[lane]);
        }
      }
    }
#pragma unroll
    for (Int lane = 0; lane < Batch; ++lane) {
      if (entry.mask & (1 << lane)) {
        mscg_solution_row<Borrowed>(solutions, first + lane, h, mrh, L)[id] = x_reg[lane];
        searches[first + lane][it] = p_reg[lane];
      }
    }
  }
}

#endif

// Cached solver workspace
struct MSCGKey {
  Geometry geo;
  Int shifts, mrh;
  bool borrowed = false;
};

inline bool operator<(const MSCGKey& a, const MSCGKey& b) {
  if (Compare_geo_less(a.geo, b.geo))
    return true;
  if (Compare_geo_less(b.geo, a.geo))
    return false;
  if (a.shifts != b.shifts)
    return a.shifts < b.shifts;
  if (a.mrh != b.mrh)
    return a.mrh < b.mrh;
  return a.borrowed < b.borrowed;
}

template <class T, class U = ComplexD> struct mscg_workspace {
  MCGDebugBuffer<U> source;
  MCGDebugBuffer<T> r;
  std::vector<MCGDebugBuffer<U>> xs;
  std::vector<MCGDebugBuffer<T>> ps;
  std::vector<U*> solution_rows;
  MCGDebugBuffer<double> coefficients;
  MCGDebugBuffer<U*> solution_ptrs;
  MCGDebugBuffer<T*> search_ptrs;
  MCGDebugBuffer<T> residual_history;
  MCGDebugBuffer<mscg_history_scalar_t<U>> coefficient_history;
  MCGDebugBuffer<MSCGShiftWork> pending_work;
  std::vector<MSCGShiftWork> window_work, replay_work;
  std::vector<char> history_rhs;
  MCGDebugBuffer<ComplexD> refine_r, refine_p, refine_ap;
  MCGDebugBuffer<ComplexD> public_source;
  MSCGReduction reduction;
  MSCGTypedReduction<double> color_reduction;
  MSCGTypedReduction<RealDD> dd_reduction;
  MSCGPairedReduction<double> cg_reduction;
  MSCGPairedReduction<RealDD> cg_dd_reduction;
  MSCGProfile profile;
  Long resize_changes = 0;
  bool borrowed_solutions = false;

  void debug_contract(Long n, Int rhs, Int shifts, Int window) {
    if (!mcg_debug_active()) return;
    mcg_debug_require(source.size()==size_t(n) && r.size()==size_t(n),"workspace source/residual sizes");
    mcg_debug_require(xs.size()==size_t(shifts) && ps.size()==size_t(shifts),"workspace shift vectors");
    mcg_debug_require(solution_ptrs.size()==size_t(borrowed_solutions?Long(shifts)*rhs:shifts) &&
      search_ptrs.size()==size_t(shifts),"workspace pointer-table lengths");
    mcg_debug_require(coefficients.size()>=size_t(2*rhs+3*Long(shifts)*rhs),"coefficient capacity");
    mcg_debug_require(pending_work.size()>=size_t(Long(shifts-1)*rhs),"replay schedule capacity");
    if (window>1 && shifts>1) {
      mcg_debug_require(residual_history.size()>=size_t(Long(window)*n),"residual history capacity");
      mcg_debug_require(coefficient_history.size()>=size_t(Long(window)*(2*rhs+3*Long(shifts)*rhs)),
        "coefficient history capacity");
    }
    source.debug_name("workspace.source");r.debug_name("workspace.residual");
    coefficients.debug_name("workspace.coefficients");solution_ptrs.debug_name("workspace.solution_ptrs");
    search_ptrs.debug_name("workspace.search_ptrs");pending_work.debug_name("workspace.pending_work");
    residual_history.debug_name("workspace.residual_history");coefficient_history.debug_name("workspace.coefficient_history");
    for(Int j=0;j<shifts;++j) {
      mcg_debug_require((borrowed_solutions || xs[j].size()==size_t(n)) &&
        (!j || ps[j].size()==size_t(n)),"solution/search vector capacities");
      xs[j].debug_name("workspace.solution."+std::to_string(j));
      ps[j].debug_name("workspace.search."+std::to_string(j));
    }
  }

  void resize(Int shifts, Long n, Int mrh, Int temporal_batch = 1, bool borrowed = false) {
    borrowed_solutions = borrowed;
    auto resize_one = [&](auto& v, Long size) {
      if (v.size() != size) {
        v.resize(size);
        ++resize_changes;
      }
    };
    resize_one(source, n);
    resize_one(r, n);
    const Long before_vectors = resize_changes;
    const Long table_size = borrowed ? Long(shifts) * mrh : shifts;
    const bool new_tables = solution_ptrs.size() != size_t(table_size);
    xs.resize(shifts);
    ps.resize(shifts);
    for (Int j = 0; j < shifts; ++j) {
      if (!borrowed)
        resize_one(xs[j], n);
      if (j)
        resize_one(ps[j], n);
    }
    resize_one(solution_ptrs, table_size);
    resize_one(search_ptrs, shifts);
    if (new_tables || resize_changes != before_vectors) {
      std::vector<U*> xp(table_size, nullptr);
      std::vector<T*> pp(shifts, nullptr);
      for (Int j = 0; j < shifts; ++j) {
        if (!borrowed)
          xp[j] = xs[j].data();
        if (j)
          pp[j] = ps[j].data();
      }
#ifdef QLAT_USE_ACC
      Qassert(qacc_Memcpy(solution_ptrs.data(), xp.data(), table_size * sizeof(U*), qacc_MemcpyHostToDevice) ==
              qacc_Success);
      Qassert(qacc_Memcpy(search_ptrs.data(), pp.data(), shifts * sizeof(T*), qacc_MemcpyHostToDevice) ==
              qacc_Success);
#else
      std::copy(xp.begin(), xp.end(), solution_ptrs.data());
      std::copy(pp.begin(), pp.end(), search_ptrs.data());
#endif
    }
    resize_one(coefficients, 2 * mrh + 3 * shifts * mrh);
    // Maximum grouping is one shift per entry. Retain these small scheduling
    // buffers when a caller changes the grouping or temporal window.
    const Long work_size = Long(shifts - 1) * mrh;
    if (pending_work.size() < work_size)
      resize_one(pending_work, work_size);
    window_work.reserve(work_size);
    replay_work.reserve(work_size);
    history_rhs.resize(mrh);
    // Batch 1 and single-shift solves never use temporal histories. Grow only
    // when needed, so switching windows cannot shrink and reallocate a cache
    // entry on every call. History indexing uses the current solve's stride.
    if (shifts > 1 && temporal_batch > 1) {
      const Long residual_size = Long(temporal_batch) * n;
      const Long coefficient_size = Long(temporal_batch) * (2 * mrh + 3 * shifts * mrh);
      if (residual_history.size() < residual_size)
        resize_one(residual_history, residual_size);
      if (coefficient_history.size() < coefficient_size)
        resize_one(coefficient_history, coefficient_size);
    }
    const Long before = reduction.allocations;
    reduction.resize(n / mrh, mrh);
    resize_changes += reduction.allocations - before;
    debug_contract(n,mrh,shifts,temporal_batch);
  }

  void prepare_refinement(Long n) {
    for (auto* v : {&refine_r, &refine_p, &refine_ap})
      if (v->size() != n) {
        v->resize(n);
        ++resize_changes;
      }
  }
};

// Call-specific output pointers must never survive in a cached workspace.
// Normal release stays asynchronous; finish() completes it before returning.
template <class T, class U> class MSCGSolutionRows {
  mscg_workspace<T, U>& w;
  bool active = true;

public:
  template <class V>
  MSCGSolutionRows(mscg_workspace<T, U>& workspace, const std::vector<V*>& rows) : w(workspace) {
    MCGDebugPhase phase("borrowed_pointer_table_upload");
    w.solution_rows.resize(rows.size());
    for (size_t i = 0; i < rows.size(); ++i)
      w.solution_rows[i] = reinterpret_cast<U*>(rows[i]);
#ifdef QLAT_USE_ACC
    Qassert(qacc_MemcpyAsync(w.solution_ptrs.data(), w.solution_rows.data(), rows.size() * sizeof(U*),
                            qacc_MemcpyHostToDevice) == qacc_Success);
#else
    std::copy(w.solution_rows.begin(), w.solution_rows.end(), w.solution_ptrs.data());
#endif
  }

  void release() {
    if (!active) return;
    MCGDebugPhase phase("borrowed_pointer_table_release");
#ifdef QLAT_USE_ACC
    Qassert(qacc_MemsetAsync(w.solution_ptrs.data(), 0, w.solution_rows.size() * sizeof(U*)) == qacc_Success);
#else
    std::fill(w.solution_ptrs.data(), w.solution_ptrs.data() + w.solution_rows.size(), nullptr);
#endif
    w.solution_rows.clear();
    active = false;
  }

  ~MSCGSolutionRows() {
    if (active) {
      release();
      // Exception cleanup must also complete use of caller-owned buffers.
      qacc_barrier(dummy);
    }
  }
};

template <class T, class U = ComplexD>
inline Cache<MSCGKey, std::shared_ptr<mscg_workspace<T, U>>>& get_mscg_cache() {
  static Cache<MSCGKey, std::shared_ptr<mscg_workspace<T, U>>> cache("MSCGWorkspace", 8);
  return cache;
}

// Resolve before allocating histories. A single shift has no replay work.
inline Int resolve_mscg_temporal_batch(Int requested, size_t shifts) {
  if (requested < 0)
    throw std::invalid_argument("temporal_batch must be 0 (automatic) or positive");
  if (shifts == 0 || shifts > size_t(std::numeric_limits<Int>::max()))
    throw std::invalid_argument("shift count must be positive and fit Int");
#ifndef QLAT_USE_ACC
  if (requested > 1)
    throw std::invalid_argument("host-only MCG supports temporal_batch 0 or 1");
  return 1;
#else
  if (shifts == 1)
    return 1;
  return requested > 0 ? requested : Int(shifts - 1);
#endif
}

template <class T, class U = ComplexD>
inline std::shared_ptr<mscg_workspace<T, U>> get_mscg_workspace(const Geometry& geo, Int shifts, Int mrh,
                                                                Long n, Int temporal_batch = 1,
                                                                bool borrowed = false) {
  MSCGKey key{geo, shifts, mrh, borrowed};
  MCGDebugPhase phase("workspace_acquire");
  auto& cache = get_mscg_cache<T, U>();
  if (!cache.has(key))
    cache[key] = std::make_shared<mscg_workspace<T, U>>();
  auto w = cache[key];
  w->resize(shifts, n, mrh, temporal_batch, borrowed);
  return w;
}

// multi/multi_precise are massless; shifts are absolute 4*m*m. An optional
// multi_seed folds the lowest shift into the operator. Independent scalar
// recurrences per RHS. Cache entries are not safe for concurrent solves.
// Multishift recurrence and optional true-residual/refinement pass
template <class M, class Tv> struct multi_shift_cg {
  using T = typename std::remove_pointer<decltype(std::declval<M>().buf[0].data())>::type;
  using Real = typename std::conditional<sizeof(T) == sizeof(ComplexD), double, float>::type;
  M* mat;
  Int iterations = 0, refine_iterations = 0, multi_count = 0;
  Long reduction_count = 0, seed_update_count = 0, shifted_update_count = 0;
  Long residual_check_count = 0, result_copy_count = 0, workspace_resize_changes = 0;
  double modeled_vector_flops = 0;
  double retirement_epsilon = 0;
  std::vector<double> true_residual, iterative_residual;
  std::vector<Int> shift_iterations;
  MSCGStoppingDiagnostics stopping;

  explicit multi_shift_cg(M& m) : mat(&m) {}

  Int internal_qacc_threads(const MSCGOptions& options) const {
    return sizeof(T) == sizeof(ComplexD) ? options.qacc_threads_double : options.qacc_threads_single;
  }

  // Raw outputs are the solution storage, including independently allocated RHS rows.
  template <class U>
  void call_raw(const std::vector<U*>& res, const std::vector<U*>& src, const std::vector<double>& shifts,
                double tolerance, Int maxiter, const MSCGOptions& options = {}) {
    Qassert(!src.empty() && res.size() == src.size() * shifts.size());
    MSCGOptions opt = options;
    opt.temporal_batch = resolve_mscg_temporal_batch(options.temporal_batch, shifts.size());
    mat->set_mrh(src.size());
    const Long N = mat->get_vsize(), L = N / src.size();
    auto& cache = get_mscg_cache<T, U>();
    MSCGKey key{mat->qinvP->geo(), Int(shifts.size()), Int(src.size()), true};
    const Long before = cache.has(key) ? cache[key]->resize_changes : 0;
    const auto acquire_start =
        options.profile ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
    auto ws = get_mscg_workspace<T, U>(mat->qinvP->geo(), shifts.size(), src.size(), N, opt.temporal_batch, true);
    const double acquire_seconds =
        options.profile
            ? std::chrono::duration<double>(std::chrono::steady_clock::now() - acquire_start).count()
            : 0;
    auto& w = *ws;
    w.profile.reset(opt.profile);
    w.profile.run("input_conversion", [&] {
      ScopedQaccNumThreads threads(internal_qacc_threads(options));
      for (Int h = 0; h < Int(src.size()); ++h) {
        U* s = src[h];
        U* b = w.source.data() + h * L;
        qGPU_for(i, L, 1, { b[i] = s[i]; });
      }
    });
    MSCGSolutionRows<T, U> rows(w, res);
    solve<true>(w, shifts, tolerance, maxiter, src.size(), opt);
    rows.release();
    workspace_resize_changes = w.resize_changes - before;
    finish(w, opt);
    if (opt.profile && opt.stats)
      opt.stats->host_seconds["workspace_acquire"] = acquire_seconds;
  }

  // Float solutions borrow the first sizeof(ComplexF)*L bytes of each double
  // result. RHS pointers need not be contiguous. No result-sized allocation is
  // made; the dead residual vector stages one mass before widening its output.
  void call_raw_mixed(const std::vector<ComplexD*>& res, const std::vector<ComplexD*>& src,
                      const std::vector<double>& shifts, double tolerance, Int maxiter,
                      const MSCGOptions& options) {
    static_assert(std::is_same_v<T, ComplexF>);
    Qassert(!src.empty() && res.size() == src.size() * shifts.size());
    MSCGOptions opt = options;
    opt.temporal_batch = resolve_mscg_temporal_batch(opt.temporal_batch, shifts.size());
    mat->set_mrh(src.size());
    const Int rhs = Int(src.size());
    const Long N = mat->get_vsize(), L = N / rhs, sites = L / 3;
    auto& cache = get_mscg_cache<T, ComplexF>();
    MSCGKey key{mat->qinvP->geo(), Int(shifts.size()), rhs, true};
    const Long before = cache.has(key) ? cache[key]->resize_changes : 0;
    auto ws = get_mscg_workspace<T, ComplexF>(mat->qinvP->geo(), shifts.size(), rhs, N,
                                             opt.temporal_batch, true);
    auto& w = *ws;
    w.profile.reset(opt.profile);
    const auto layout = mat->vector_layout(opt.generic_vector_layout);
    const bool check = opt.compute_true_res.value_or(false) || opt.refine_tolerance > 0;
    if (check && w.public_source.size() != size_t(N)) {
      w.public_source.resize(N);
      ++w.resize_changes;
    }
    w.profile.run("input_conversion", [&] {
      ScopedQaccNumThreads threads(opt.qacc_threads_single);
      for (Int h = 0; h < rhs; ++h) {
        const ComplexD* input = src[h];
        ComplexF* b = w.source.data() + Long(h) * L;
        qGPU_for(i, L, 1, {
          const Long site = i / 3;
          const Int color = Int(i % 3);
          b[layout.storage_index(site, color, sites)] =
              ComplexF(input[layout.template public_index<ComplexD>(site, color, sites)]);
        });
        if (check) {
          ComplexD* saved = w.public_source.data() + Long(h) * L;
          qGPU_for(i, L, 1, {
            const Long site = i / 3;
            const Int color = Int(i % 3);
            const Long id = layout.template public_index<ComplexD>(site, color, sites);
            saved[id] = input[id];
          });
        }
      }
    });
    MSCGSolutionRows<T, ComplexF> rows(w, res);
    solve<true, true>(w, shifts, tolerance, maxiter, rhs, opt);
    w.profile.run("output_conversion", [&] {
      ScopedQaccNumThreads threads(opt.qacc_threads_single);
      for (size_t j = 0; j < shifts.size(); ++j) {
        for (Int h = 0; h < rhs; ++h) {
          const ComplexF* packed = w.solution_rows[j * rhs + h];
          ComplexF* scratch = w.r.data() + Long(h) * L;
          qGPU_for(i, L, 1, { scratch[i] = packed[i]; });
        }
        for (Int h = 0; h < rhs; ++h) {
          const ComplexF* scratch = w.r.data() + Long(h) * L;
          ComplexD* output = res[j * rhs + h];
          qGPU_for(i, L, 1, {
            const Long site = i / 3;
            const Int color = Int(i % 3);
            output[layout.template public_index<ComplexD>(site, color, sites)] =
                ComplexD(scratch[layout.storage_index(site, color, sites)]);
          });
          ++result_copy_count;
        }
      }
    });
    if (check)
      check_exported_double(w, res, shifts, tolerance, maxiter, rhs, opt);
    rows.release();
    workspace_resize_changes = w.resize_changes - before;
    finish(w, opt);
  }

  // Public-double observation/internal correction after a borrowed float solve.
  // The original source is saved only when this optional work is requested.
  void check_exported_double(mscg_workspace<T, ComplexF>& w, const std::vector<ComplexD*>& res,
                             const std::vector<double>& shifts, double tolerance, Int maxiter,
                             Int rhs, const MSCGOptions& opt) {
    const Long N = mat->get_vsize(), L = N / rhs;
    w.prepare_refinement(N);
    auto* b = w.public_source.data();
    auto* r = w.refine_r.data();
    auto* p = w.refine_p.data();
    auto* ap = w.refine_ap.data();
    const auto layout = mat->vector_layout(opt.generic_vector_layout);
    const auto config = get_mscg_reduction_config();
    auto dot = [&](const ComplexD* a, const ComplexD* c) {
      const double* values = config.precision ? w.dd_reduction.dot(a, c, L, rhs, config.color_sum, layout, true)
                             : config.color_sum ? w.color_reduction.dot(a, c, L, rhs, true, layout, true)
                                                : w.reduction.dot(a, c, L, rhs);
      ++reduction_count;
      modeled_vector_flops += 4.0 * N;
      return std::vector<double>(values, values + rhs);
    };
    auto precise = [&] {
      mat->multi_precise(ap, p);
      ++multi_count;
    };
    const auto b2 = dot(b, b);
    const Int limit = opt.refine_maxiter >= 0 ? opt.refine_maxiter : maxiter;
    ScopedQaccNumThreads threads(opt.qacc_threads_double);
    w.profile.run("public_checks_refinement", [&] {
      for (size_t j = 0; j < shifts.size(); ++j) {
        const double sigma = shifts[j];
        Int used = 0;
        while (true) {
          for (Int h = 0; h < rhs; ++h) {
            const ComplexD* output = res[j * rhs + h];
            ComplexD* ph = p + Long(h) * L;
            qGPU_for(i, L, 1, { ph[i] = output[i]; });
          }
          precise();
          qGPU_for(i, N, 1, { r[i] = b[i] - ap[i] - sigma * p[i]; p[i] = r[i]; });
          modeled_vector_flops += 6.0 * N;
          auto rr = dot(r, r);
          std::vector<char> done(rhs);
          bool all = true;
          const double target = opt.refine_tolerance > 0 ? std::min(tolerance, opt.refine_tolerance) : tolerance;
          for (Int h = 0; h < rhs; ++h) {
            const double residual = b2[h] == 0 ? std::sqrt(rr[h]) : std::sqrt(rr[h] / b2[h]);
            true_residual[j * rhs + h] = residual;
            ++residual_check_count;
            done[h] = b2[h] == 0 || residual <= target;
            all = all && done[h];
            if (done[h]) {
              ComplexD* ph = p + Long(h) * L;
              qGPU_for(i, L, 1, { ph[i] = 0.0; });
            }
          }
          if (opt.refine_tolerance <= 0 || all || used >= limit)
            break;
          for (Int pass = 0; pass < 32 && used < limit; ++pass) {
            precise();
            qGPU_for(i, N, 1, { ap[i] += sigma * p[i]; });
            auto pap = dot(p, ap);
            for (Int h = 0; h < rhs; ++h)
              if (!done[h]) {
                Qassert(pap[h] > 0 && std::isfinite(pap[h]));
                const double alpha = rr[h] / pap[h];
                ComplexD* xh = res[j * rhs + h];
                ComplexD *rh = r + Long(h) * L, *ph = p + Long(h) * L, *ah = ap + Long(h) * L;
                qGPU_for(i, L, 1, { xh[i] += alpha * ph[i]; rh[i] -= alpha * ah[i]; });
                modeled_vector_flops += 8.0 * L;
              }
            auto rn = dot(r, r);
            bool finished = true;
            for (Int h = 0; h < rhs; ++h)
              if (!done[h]) {
                const double beta = rn[h] / rr[h];
                done[h] = rn[h] <= target * target * b2[h];
                finished = finished && done[h];
                const double keep = done[h] ? 0.0 : 1.0;
                ComplexD *rh = r + Long(h) * L, *ph = p + Long(h) * L;
                qGPU_for(i, L, 1, { ph[i] = keep * (rh[i] + beta * ph[i]); });
                modeled_vector_flops += 6.0 * L;
              }
            rr = rn;
            ++used;
            ++refine_iterations;
            if (finished)
              break;
          }
        }
      }
    });
    if (opt.print_mem)
      print_mem_info(refine_iterations > 0 ? "MCG refinement done" : "MCG residual checks done");
  }

  // Retain the vector_cs entry point; its source/results now follow double I/O.
  void call(std::vector<Tv>& x, Tv& src, const std::vector<double>& shifts, double tolerance, Int maxiter,
            Int ir = 0, Int is = 0, Int mrh = 1, double refine_tolerance = 0, Int refine_maxiter = -1,
            const std::vector<double>& tolerances = {}, const MSCGOptions& options = {}) {
    static_assert(std::is_same<decltype(std::declval<Tv>().get_pointer_b(0, 0)), ComplexD*>::value,
                  "MCG public input/output must use double complex vectors");
    Qassert(tolerances.empty() || tolerances.size() == shifts.size());
    MSCGOptions opt = options;
    opt.temporal_batch = resolve_mscg_temporal_batch(options.temporal_batch, shifts.size());
    mat->set_mrh(mrh);
    const Long N = mat->get_vsize();
    auto& cache = get_mscg_cache<T>();
    MSCGKey key{mat->qinvP->geo(), Int(shifts.size()), mrh};
    const Long before = cache.has(key) ? cache[key]->resize_changes : 0;
    const auto acquire_start =
        options.profile ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
    auto ws = get_mscg_workspace<T>(mat->qinvP->geo(), shifts.size(), mrh, N, opt.temporal_batch);
    auto& w = *ws;
    const double acquire_seconds =
        options.profile
            ? std::chrono::duration<double>(std::chrono::steady_clock::now() - acquire_start).count()
            : 0;
    if (refine_tolerance > 0)
      opt.refine_tolerance = refine_tolerance;
    if (refine_maxiter >= 0)
      opt.refine_maxiter = refine_maxiter;
    Qassert(src.nsum == N / mrh && is >= 0 && is + mrh <= src.nvec);
    w.profile.reset(options.profile);
    {
      ScopedQaccNumThreads threads(opt.qacc_threads_double);
      src.copy_to_group(w.source.data(), is, mrh);
    }
    solve(w, shifts, tolerance, maxiter, mrh, opt, tolerances);
    x.resize(shifts.size());
    for (size_t j = 0; j < shifts.size(); ++j) {
      if (!x[j].initialized)
        x[j].resize(ir + mrh, src);
      {
        ScopedQaccNumThreads threads(opt.qacc_threads_double);
        x[j].copy_from_group(w.xs[j].data(), ir, mrh);
      }
      ++result_copy_count;
    }
    workspace_resize_changes = w.resize_changes - before;
    finish(w, opt);
    if (opt.profile && opt.stats)
      opt.stats->host_seconds["workspace_acquire"] = acquire_seconds;
  }

  template <class U> void finish(mscg_workspace<T, U>& w, const MSCGOptions& options) {
    MCGDebugPhase phase("solver_finish");
    qacc_barrier(dummy);
    w.profile.finish();
    if (options.stats) {
      auto& s = *options.stats;
      s.iterations = iterations;
      s.refine_iterations = refine_iterations;
      s.matrix_calls = multi_count;
      s.reductions = reduction_count;
      s.seed_updates = seed_update_count;
      s.shift_updates = shifted_update_count;
      s.residual_checks = residual_check_count;
      s.workspace_allocations = workspace_resize_changes;
      s.shift_iterations = shift_iterations;
      s.true_residual = true_residual;
      s.iterative_residual = iterative_residual;
      s.stopping = stopping;
      s.host_seconds = w.profile.host;
      s.device_seconds = w.profile.device;
      s.modeled_vector_flops = modeled_vector_flops;
      s.borrowed_solutions = w.borrowed_solutions;
      s.accumulation_bytes = sizeof(U) / 2;
      s.retirement_epsilon = retirement_epsilon;
      s.modeled_operator_flops = 0;
      s.solve_seconds = 0;
    }
  }

  template <bool Borrowed = false, bool PackedFloat = false, class U>
  void solve(mscg_workspace<T, U>& w, const std::vector<double>& shifts, double tolerance, Int maxiter,
             Int mrh, const MSCGOptions& opt, const std::vector<double>& targets = {}) {
    Qassert(mrh > 0 && !shifts.empty() && tolerance > 0 && maxiter >= 0);
    for (size_t j = 0; j < shifts.size(); ++j) {
      Qassert(std::isfinite(shifts[j]) && shifts[j] >= 0);
      Qassert(j == 0 || shifts[j] >= shifts[j - 1]);
    }
    iterations = refine_iterations = multi_count = 0;
    reduction_count = seed_update_count = shifted_update_count = 0;
    residual_check_count = result_copy_count = 0;
    modeled_vector_flops = 0;
    const Long N = mat->get_vsize(), L = N / mrh;
    const size_t J = shifts.size(), S = J * mrh;
    w.debug_contract(N,mrh,J,opt.temporal_batch);
    Qassert(L % 3 == 0 && (opt.vector_block == 64 || opt.vector_block == 128 || opt.vector_block == 256));
    Qassert(opt.shift_batch == 1 || opt.shift_batch == 2 || opt.shift_batch == 4 || opt.shift_batch == 6);
    Qassert(opt.temporal_batch > 0);
    Qassert(opt.qacc_threads_single > 0 && opt.qacc_threads_single <= 1024 && opt.qacc_threads_double > 0 &&
            opt.qacc_threads_double <= 1024);
    const Long sites = L / 3;
    Qassert(opt.shift_stop == MSCGShiftStopPolicy::inherit ||
            opt.shift_stop == MSCGShiftStopPolicy::requested ||
            opt.shift_stop == MSCGShiftStopPolicy::cpu_epsilon ||
            opt.shift_stop == MSCGShiftStopPolicy::precision_epsilon);
    Qassert(opt.beta == MSCGBetaPolicy::inherit || opt.beta == MSCGBetaPolicy::norm_ratio ||
            opt.beta == MSCGBetaPolicy::quda);
    const bool epsilon_retirement = opt.shift_stop == MSCGShiftStopPolicy::cpu_epsilon ||
                                    opt.shift_stop == MSCGShiftStopPolicy::precision_epsilon;
    retirement_epsilon = !epsilon_retirement ? 0 :
        opt.shift_stop == MSCGShiftStopPolicy::precision_epsilon ? std::numeric_limits<Real>::epsilon()
                                                               : std::numeric_limits<double>::epsilon();
    const bool quda_beta = opt.beta == MSCGBetaPolicy::quda;
    const auto layout = [&]() {
      if constexpr (requires { mat->vector_layout(opt.generic_vector_layout); })
        return mat->vector_layout(opt.generic_vector_layout);
      else
        return MSCGFlatLayout{};
    }();
    // Only the staggered adapter opts into the process-wide controls. Generic
    // matrix callers and independent benchmark/reference reductions retain the
    // historic arithmetic unless they explicitly call MSCGTypedReduction.
    const MSCGReductionConfig reduction_config = [&] {
      if constexpr (requires { M::staggered_reductions; })
        return M::staggered_reductions ? get_mscg_reduction_config() : MSCGReductionConfig{};
      else
        return MSCGReductionConfig{};
    }();
    if (reduction_config.precision || reduction_config.color_sum)
      w.profile.run("reduction_prepare", [&] {
        if (reduction_config.precision == 1) {
          const Long before = w.dd_reduction.allocations;
          w.dd_reduction.prepare(L, mrh, reduction_config.color_sum);
          w.resize_changes += w.dd_reduction.allocations - before;
        } else if (reduction_config.color_sum) {
          const Long before = w.color_reduction.allocations;
          w.color_reduction.prepare(L, mrh, true);
          w.resize_changes += w.color_reduction.allocations - before;
        }
      });
    if (quda_beta)
      w.profile.run("cg_reduction_prepare", [&] {
        auto prepare = [&](auto& reduction) {
          const Long before = reduction.allocations;
          reduction.prepare(L, 2 * mrh, reduction_config.color_sum);
          w.resize_changes += reduction.allocations - before;
        };
        if (reduction_config.precision == 1)
          prepare(w.cg_dd_reduction);
        else
          prepare(w.cg_reduction);
      });
    // Collective even on cache reuse: allocation counts can differ by rank.
    if (opt.print_mem)
      print_mem_info("MCG workspace ready");
    vector<T*> mp;
    mat->set_buf_pointers(mp, mrh);
    U* b = w.source.data();
    T* r = w.r.data();
    T* p = mp[0];
    T* ap = mp[1];
    auto dot = [&](auto* a, auto* c, bool public_vectors = false) {
      const double* out = nullptr;
      w.profile.run("reductions", [&] {
        if (reduction_config.precision == 1)
          out = w.dd_reduction.dot(a, c, L, mrh, reduction_config.color_sum, layout, public_vectors);
        else if (reduction_config.color_sum)
          out = w.color_reduction.dot(a, c, L, mrh, true, layout, public_vectors);
        else
          out = w.reduction.dot(a, c, L, mrh);
      });
      modeled_vector_flops += 4.0 * N;
      ++reduction_count;
      return std::vector<double>(out, out + mrh);
    };
    std::vector<double> b2 = dot(b, b, true), rr = b2, aold(mrh, 1), bold(mrh, 0);
    std::vector<double> z(S, 1), zo(S, 1), stop(J, tolerance), coeff(2 * mrh + 3 * S);
    if (!targets.empty())
      stop = targets;
    for (double target : stop)
      Qassert(std::isfinite(target) && target > 0);
    std::vector<double> retirement_threshold2(J * mrh);
    for (size_t j = 0; j < J; ++j)
      for (Int h = 0; h < mrh; ++h) {
        const double target = j && epsilon_retirement ? retirement_epsilon : stop[j];
        retirement_threshold2[j * mrh + h] = target * target * b2[h];
      }
    std::vector<char> done(mrh, 0), active(S, 1);
    true_residual.assign(S, std::numeric_limits<double>::quiet_NaN());
    iterative_residual.assign(S, 0);
    shift_iterations.assign(S, 0);
    for (size_t j = 0; j < J; ++j)
      for (Int h = 0; h < mrh; ++h)
        iterative_residual[j * mrh + h] = b2[h] == 0 ? 0.0 : 1.0;
    w.profile.run("workspace_prepare", [&] {
#ifdef QLAT_USE_ACC
      mscg_initialize_vectors<T, U>
          <<<dim3((sites + opt.vector_block - 1) / opt.vector_block, mrh), opt.vector_block>>>(
              b, r, p, w.search_ptrs.data(), J, sites, layout);
      Qassert(qacc_GetLastError() == qacc_Success);
      qacc_barrier(dummy);
#else
      ScopedQaccNumThreads threads(internal_qacc_threads(opt));
      qGPU_for(i,N/3,0,{
        const Long site=i%sites,base=(i/sites)*L;
        for(Int color=0;color<3;++color) {
          const Long it=base+layout.storage_index(site,color,sites);
          r[it]=p[it]=T(b[base+layout.template public_index<U>(site,color,sites)]);
        }
      });
      for(size_t j=1;j<J;++j) { T* pj=w.ps[j].data();qGPU_for(i,N,0,{ pj[i]=r[i]; }); }
#endif
      for (size_t j = 0; j < J; ++j) {
        if constexpr (Borrowed) {
          for (Int h = 0; h < mrh; ++h) {
            U* xx = w.solution_rows[j * mrh + h];
            qGPU_for(i, L, 1, { xx[i] = U(0); });
          }
        } else {
          U* xx = w.xs[j].data();
          ScopedQaccNumThreads threads(internal_qacc_threads(opt));
          qGPU_for(i, N, 1, { xx[i] = 0.0; });
        }
      }
    });
    for (Int h = 0; h < mrh; ++h)
      done[h] = b2[h] == 0 || (!opt.fixed_iterations && b2[h] < stop[0] * stop[0] * b2[h]);
    stopping = {};
    if (opt.stopping_diagnostics) {
      stopping.maxiter = maxiter;
      stopping.source_norm2 = b2;
      stopping.seed_threshold2.resize(mrh);
      stopping.shift_threshold2.resize(S);
      stopping.seed_checks.resize(mrh);
      stopping.shift_checks.resize(S);
      const auto limit = opt.fixed_iterations ? MSCGStopReason::fixed_limit : MSCGStopReason::iteration_limit;
      stopping.seed_reason.assign(mrh, limit);
      stopping.shift_reason.assign(S, limit);
      stopping.shift_retired.assign(S, 0);
      stopping.effective_updates.assign(S, 0);
      for (Int h = 0; h < mrh; ++h) {
        stopping.seed_threshold2[h] = stop[0] * stop[0] * b2[h];
        for (size_t j = 0; j < J; ++j)
          stopping.shift_threshold2[j * mrh + h] = retirement_threshold2[j * mrh + h];
        if (done[h])
          stopping.seed_reason[h] = b2[h] == 0 ? MSCGStopReason::zero_source : MSCGStopReason::initial_target;
      }
    }
    Int pending_temporal = 0;
    const Int shift_groups = (J - 1 + opt.shift_batch - 1) / opt.shift_batch;
    w.window_work.resize(Long(shift_groups) * mrh);
    for (Int g = 0; g < shift_groups; ++g)
      for (Int h = 0; h < mrh; ++h)
        w.window_work[Long(g) * mrh + h] = {h, 1 + g * opt.shift_batch, 0, 0};
    bool dense_window = true;
    double* c = w.coefficients.data();
    // Launch one shifted-vector replay for the pending temporal history.
    auto launch_shift_window = [&](Int steps) {
      MCGDebugPhase phase("shift_replay");
      if (mcg_debug_active()) {
        mcg_debug_replay_steps(steps,opt.temporal_batch);
        for (const auto& e:w.window_work) if (e.mask) {
          mcg_debug_require(e.rhs>=0 && e.rhs<mrh && e.first>=1 && size_t(e.first)<J,
            "replay bounds: RHS/first shift");
          mcg_debug_require(e.steps>0 && e.steps<=steps,"replay bounds: entry history steps");
          const Int lanes=std::min<Int>(opt.shift_batch,J-e.first);
          mcg_debug_require((e.mask & ~((1<<lanes)-1))==0,"replay bounds: shift mask");
        }
      }
#ifdef QLAT_USE_ACC
      if (J <= 1)
        return;
      w.replay_work.clear();
      for (const auto& entry : w.window_work)
        if (entry.mask)
          w.replay_work.push_back(entry);
      mcg_debug_require(w.replay_work.size()<=w.pending_work.size(),"replay bounds: schedule capacity");
      if (w.replay_work.empty())
        return;
      if (!dense_window)
        w.profile.run("shift_schedule_copy", [&] {
          Qassert(qacc_MemcpyAsync(w.pending_work.data(), w.replay_work.data(),
                                  w.replay_work.size() * sizeof(MSCGShiftWork),
                                  qacc_MemcpyHostToDevice) == qacc_Success);
        });
      const dim3 grid((sites + opt.vector_block - 1) / opt.vector_block,
                      dense_window ? mrh : w.replay_work.size(), dense_window ? shift_groups : 1);
      auto launch = [&]<int Batch>() {
        if (!dense_window) {
          auto sparse = [&]<bool Temporal>() {
            mscg_shift_vectors_pending<Temporal, Batch, T, U, Borrowed><<<grid, opt.vector_block>>>(
                w.solution_ptrs.data(), w.search_ptrs.data(),
                Temporal ? w.residual_history.data() : r, w.coefficient_history.data(), c,
                w.pending_work.data(), J, mrh, sites, 2 * mrh + 3 * S, layout);
          };
          if (opt.temporal_batch == 1)
            sparse.template operator()<false>();
          else
            sparse.template operator()<true>();
        } else if (opt.temporal_batch == 1)
          mscg_shift_vectors<Batch, T, U, Borrowed><<<grid, opt.vector_block>>>(
              w.solution_ptrs.data(), w.search_ptrs.data(), r, c, J, mrh, sites, layout);
        else {
          const Long stride = 2 * mrh + 3 * S;
          auto replay = [&]<int Window>() {
            mscg_shift_vectors_temporal<Window, Batch, T, U, Borrowed><<<grid, opt.vector_block>>>(
                w.solution_ptrs.data(), w.search_ptrs.data(), w.residual_history.data(),
                w.coefficient_history.data(), steps, J, mrh, sites, stride, layout);
          };
          if (opt.temporal_batch == 4)
            replay.template operator()<4>();
          else if (opt.temporal_batch == 8)
            replay.template operator()<8>();
          else if (opt.temporal_batch == 16)
            replay.template operator()<16>();
          else
            replay.template operator()<0>(); // Runtime window; same per-step arithmetic.
        }
      };
      w.profile.run("shift_kernel", [&] {
        switch (opt.shift_batch) {
        case 1:
          launch.template operator()<1>();
          break;
        case 2:
          launch.template operator()<2>();
          break;
        case 4:
          launch.template operator()<4>();
          break;
        case 6:
          launch.template operator()<6>();
          break;
        }
        Qassert(qacc_GetLastError() == qacc_Success);
      });
      // Keep the kernel launch and its completion point visible separately
      // in the optional profile while retaining the barrier required by the
      // host-side recurrence and the following seed update.
      w.profile.run("shift_barrier", [&] { qacc_barrier(dummy); });
#else
      Qassert(opt.temporal_batch == 1);
      if (J <= 1)
        return;
      using Scalar = decltype(U{}.real());
      ScopedQaccNumThreads threads(internal_qacc_threads(opt));
      for (size_t j = 1; j < J; ++j) {
        T* pj = w.ps[j].data();
        const Long off = 2 * mrh + j * mrh;
        qGPU_for(i, N / 3, 0, {
          const Long site = i % sites, h = i / sites, base = h * L, v = off + h;
          U* xx = Borrowed ? w.solution_rows[j * mrh + h] : w.xs[j].data() + base;
          const Scalar residual_coefficient = Scalar(c[v]);
          if (residual_coefficient != Scalar(0))
            for (Int color = 0; color < 3; ++color) {
              const Long it = base + layout.storage_index(site, color, sites),
                         id = layout.template public_index<U>(site, color, sites);
              const Scalar solution_coefficient = Scalar(c[v + S]);
              const Scalar search_coefficient = Scalar(c[v + 2 * S]);
              const U old = U(pj[it]);
              const U residual = U(r[it]);
              xx[id] += U(solution_coefficient * old);
              pj[it] = mscg_search_update<T, U, decltype(layout)>(T(residual), pj[it],
                                                               residual_coefficient, search_coefficient);
            }
        });
      }
#endif
      for (auto& entry : w.window_work)
        entry.mask = entry.steps = 0;
      dense_window = true;
    };
    for (Int k = 0; k < maxiter; ++k) {
      bool all = true;
      for (auto d : done)
        if (!d)
          all = false;
      if (all)
        break;
      if (mcg_debug_active()) mcg_debug_state().iteration=k+1;
      // Use the original solver RHS count, even after some streams finish.
      // Residual norms are computed every iteration; replay size is independent.
      const bool check_stopping = (k + 1) % mrh == 0;
      w.profile.run("dslash", [&] {
        if constexpr (requires { mat->multi_seed(ap, p, shifts[0]); })
          mat->multi_seed(ap, p, shifts[0]);
        else
          mat->multi(ap, p);
      });
      ++multi_count;
      if constexpr (!requires { mat->multi_seed(ap, p, shifts[0]); }) {
        const Real sigma = Real(shifts[0]);
        w.profile.run("seed_shift", [&] {
          ScopedQaccNumThreads threads(internal_qacc_threads(opt));
          qGPU_for(i, N, 1, { ap[i] += sigma * p[i]; });
        });
      }
      // Count the Xpay once, including when it is fused into multi_seed.
      modeled_vector_flops += 4.0 * N;
      auto pap = dot(p, ap);
      Int effective_nonseed = 0;
      bool future_nonseed = false;
      std::fill(w.history_rhs.begin(), w.history_rhs.end(), 0);
      w.profile.run("host_recurrence", [&] {
        for (Int h = 0; h < mrh; ++h) {
          if (done[h])
            coeff[h] = 0;
          else {
            Qassert(std::isfinite(pap[h]) && pap[h] > 0);
            coeff[h] = rr[h] / pap[h];
          }
        }
      });
      w.profile.run("coefficient_copy", [&] {
        // In the temporal path, the residual kernel and following reduction
        // are ordered on the default stream.  The reduction's host readback
        // completes before coeff is updated again, so an extra device-wide
        // barrier after this upload is unnecessary.
        cpy_GPU(c, coeff.data(), mrh, 1, 0, opt.temporal_batch == 1 ? QTRUE : QFALSE);
      });
      std::vector<double> rn, beta_numerator;
      if (quda_beta) {
        // The modified beta needs both rounded residuals, so use the paired
        // update/reduction even when ordinary norm fusion was disabled.
        w.profile.run("residual_update_cg_norm", [&] {
          const double* values = reduction_config.precision == 1
                                     ? w.cg_dd_reduction.template residual_cg_norm<T, Real>(
                                           r, ap, c, L, mrh, reduction_config.color_sum, layout)
                                     : w.cg_reduction.template residual_cg_norm<T, Real>(
                                           r, ap, c, L, mrh, reduction_config.color_sum, layout);
          rn.assign(values, values + mrh);
          beta_numerator.assign(values + mrh, values + 2 * mrh);
        });
        ++reduction_count;
        modeled_vector_flops += 14.0 * N; // residual axpy, norm, difference and real dot
      } else if (opt.fuse_residual_norm) {
        w.profile.run("residual_update_norm", [&] {
          const double* values = nullptr;
          if (reduction_config.precision == 1)
            values = w.dd_reduction.template residual_norm<T, Real>(r, ap, c, L, mrh,
                                                                    reduction_config.color_sum, layout);
          else if (reduction_config.color_sum)
            values = w.color_reduction.template residual_norm<T, Real>(r, ap, c, L, mrh, true, layout);
          else
            values = w.reduction.template residual_norm<T, Real>(r, ap, c, L, mrh);
          rn.assign(values, values + mrh);
        });
        ++reduction_count;
        modeled_vector_flops += 8.0 * N; // residual axpy + norm
      } else {
        w.profile.run("residual_update", [&] {
          ScopedQaccNumThreads threads(internal_qacc_threads(opt));
          qGPU_for(i, N, 1, { r[i] -= Real(c[i / L]) * ap[i]; });
        });
        modeled_vector_flops += 4.0 * N;
        rn = dot(r, r);
      }
      w.profile.run("host_recurrence", [&] {
        for (Int h = 0; h < mrh; ++h) {
          if (done[h]) {
            coeff[mrh + h] = 0;
            continue;
          }
          Qassert(std::isfinite(rn[h]) && rn[h] >= 0);
          Int remaining = 0;
          const double a = coeff[h], beta = (quda_beta ? beta_numerator[h] : rn[h]) / rr[h];
          Qassert(std::isfinite(beta));
          coeff[mrh + h] = beta;
          for (size_t j = 0; j < J; ++j) {
            const size_t v = j * mrh + h;
            if (!active[v]) {
              coeff[2 * mrh + v] = 0;
              continue;
            }
            double zn = 1;
            if (j) {
              const double den =
                  a * bold[h] * (zo[v] - z[v]) + zo[v] * aold[h] * (1 + (shifts[j] - shifts[0]) * a);
              Qassert(std::isfinite(den) && den != 0);
              zn = z[v] * zo[v] * aold[h] / den;
            }
            const double ratio = zn / z[v];
            coeff[2 * mrh + v] = zn;
            coeff[2 * mrh + S + v] = a * ratio;
            coeff[2 * mrh + 2 * S + v] = beta * ratio * ratio;
            zo[v] = z[v];
            z[v] = zn;
            iterative_residual[v] = std::sqrt(rn[h] * zn * zn / b2[h]);
            // Match the actual shifted-kernel branch, including FP32
            // coefficient underflow; the seed kernel runs for every RHS.
            using UpdateReal = decltype(U{}.real());
            if (j && UpdateReal(zn) != UpdateReal(0)) {
              modeled_vector_flops += 10.0 * L;
              auto& entry = w.window_work[((j - 1) / opt.shift_batch) * mrh + h];
              entry.mask |= 1 << ((j - 1) % opt.shift_batch);
              entry.steps = pending_temporal + 1;
              w.history_rhs[h] = 1;
              ++effective_nonseed;
            }
            shift_iterations[v] = k + 1;
            ++shifted_update_count;
            // Keep this boundary iteration's coefficients for its vector update
            // (or temporal replay); inactive shifts are masked next iteration.
            if (!opt.fixed_iterations && check_stopping && j && rn[h] * zn * zn < retirement_threshold2[v])
              active[v] = 0;
            if (j && active[v])
              ++remaining;
            if (opt.stopping_diagnostics) {
              if (!j || UpdateReal(zn) != UpdateReal(0))
                ++stopping.effective_updates[v];
              if (!opt.fixed_iterations && check_stopping && j) {
                stopping.shift_checks[v].record(k + 1, rn[h] * zn * zn, stopping.shift_threshold2[v]);
                if (!active[v]) {
                  stopping.shift_reason[v] =
                      epsilon_retirement ? MSCGStopReason::epsilon : MSCGStopReason::target;
                  stopping.shift_retired[v] = 1;
                }
              }
            }
          }
          rr[h] = rn[h];
          aold[h] = a;
          bold[h] = beta;
          ++seed_update_count;
          if (rn[h] == 0 || (!opt.fixed_iterations && check_stopping && rn[h] < stop[0] * stop[0] * b2[h]))
            done[h] = 1;
          if (!done[h] && remaining)
            future_nonseed = true;
          if (opt.stopping_diagnostics) {
            if (!opt.fixed_iterations && check_stopping)
              stopping.seed_checks[h].record(k + 1, rn[h], stopping.seed_threshold2[h]);
            if (done[h])
              stopping.seed_reason[h] = rn[h] == 0 ? MSCGStopReason::exact_zero : MSCGStopReason::target;
          }
        }
      });
      w.profile.run("coefficient_copy", [&] {
        cpy_GPU(c, coeff.data(), coeff.size(), 1, 0, opt.temporal_batch == 1 ? QTRUE : QFALSE);
      });
      w.profile.run("shift_updates", [&] {
        if (J <= 1)
          return;
        if (!effective_nonseed) {
          if (pending_temporal) {
            launch_shift_window(pending_temporal);
            pending_temporal = 0;
          }
          return;
        }
        dense_window = dense_window && effective_nonseed == Long(J - 1) * mrh;
        if (opt.temporal_batch == 1) {
          launch_shift_window(1);
          return;
        }
        mcg_debug_require(pending_temporal>=0 && pending_temporal<opt.temporal_batch,
          "history bounds: temporal slot");
        const Long coeff_size = 2 * mrh + 3 * S;
        // These copies and the following temporal kernel all use the default
        // stream.  Keep their ordering through stream semantics, but avoid a
        // device-wide synchronize after every history copy; the kernel launch
        // at the end of a window already provides the required completion
        // point before the next host-side recurrence consumes the results.
        w.profile.run("shift_history_copy", [&] {
          // Copy contiguous live RHS ranges. Slots for other RHS are never
          // read by sparse replay, including when a RHS finishes mid-window.
          for (Int h = 0; h < mrh;) {
            if (!w.history_rhs[h]) {
              ++h;
              continue;
            }
            const Int first = h++;
            while (h < mrh && w.history_rhs[h])
              ++h;
            cpy_GPU(w.residual_history.data() + Long(pending_temporal) * N + Long(first) * L,
                    r + Long(first) * L, Long(h - first) * L, 1, 1, QFALSE);
          }
          const Long offset = Long(pending_temporal) * coeff_size;
          if constexpr (std::is_same_v<mscg_history_scalar_t<U>, double>) {
            cpy_GPU(w.coefficient_history.data() + offset, c, coeff_size, 1, 1, QFALSE);
          } else {
            // The recurrence keeps its host-facing coefficients in double.
            // Convert only the temporal history for the all-FP32 path so that
            // the shifted update performs its arithmetic in single precision.
#ifdef QLAT_USE_ACC
            mscg_convert_coeff_history<U><<<1, 128>>>(c, w.coefficient_history.data() + offset, coeff_size);
            Qassert(qacc_GetLastError() == qacc_Success);
#else
            for (Long i = 0; i < coeff_size; ++i)
              w.coefficient_history[offset + i] = mscg_history_scalar_t<U>(c[i]);
#endif
          }
        });
        ++pending_temporal;
        if (pending_temporal == opt.temporal_batch || k + 1 == maxiter || !future_nonseed) {
          launch_shift_window(pending_temporal);
          pending_temporal = 0;
        }
      });
      U* x0 = Borrowed ? nullptr : w.xs[0].data();
      w.profile.run("seed_update", [&] {
#ifdef QLAT_USE_ACC
        mscg_seed_vectors<T, U, Borrowed>
            <<<dim3((sites + opt.vector_block - 1) / opt.vector_block, mrh), opt.vector_block>>>(
                x0, p, r, c, mrh, sites, layout, w.solution_ptrs.data());
        Qassert(qacc_GetLastError() == qacc_Success);
        qacc_barrier(dummy);
#else
        ScopedQaccNumThreads threads(internal_qacc_threads(opt));
        qGPU_for(i,N/3,0,{
          const Long site=i%sites,h=i/sites,base=h*L;
          U* xx = Borrowed ? w.solution_rows[h] : x0 + base;
          for(Int color=0;color<3;++color) {
            const Long it=base+layout.storage_index(site,color,sites),id=layout.template public_index<U>(site,color,sites);
            const U old=U(p[it]);xx[id]+=U(c[h]*old);
            p[it]=mscg_search_update<T,U,decltype(layout)>(r[it],p[it],1.0,c[mrh+h]);
          }
        });
#endif
      });
      modeled_vector_flops += 8.0 * N; // seed solution/search update, including masked RHS
      // Mask finished RHS streams without stopping the remaining MRH sources.
      for (Int h = 0; h < mrh; ++h)
        if (done[h]) {
          T* ph = p + h * L;
          {
            ScopedQaccNumThreads threads(internal_qacc_threads(opt));
            qGPU_for(i, L, 1, { ph[i] = 0.0; });
          }
          coeff[h] = 0;
          for (size_t j = 0; j < J; ++j)
            coeff[2 * mrh + j * mrh + h] = 0;
        }
      iterations = k + 1;
    }
    if (pending_temporal > 0) {
      w.profile.run("shift_updates", [&] { launch_shift_window(pending_temporal); });
      pending_temporal = 0;
    }
    // Freeze base-MCG decisions before optional true-residual/internal-CG work
    // reuses done/rr. Retired shifts retain their own target stopping checks.
    if (opt.stopping_diagnostics)
      for (size_t j = 0; j < J; ++j)
        for (Int h = 0; h < mrh; ++h) {
          const size_t v = j * mrh + h;
          if (active[v])
            stopping.shift_reason[v] = stopping.seed_reason[h];
        }
    // Only packed float-in-double outputs require export before observation.
    // Direct same-type output rows can be checked/refined here.
    if constexpr (PackedFloat)
      return;
    if constexpr (!std::is_same<U, ComplexD>::value) {
      Qassert(opt.refine_tolerance <= 0);
      if (opt.compute_true_res.value_or(false)) {
        for (size_t j = 0; j < J; ++j) {
          for (Int h = 0; h < mrh; ++h) {
            const U* xx = Borrowed ? w.solution_rows[j * mrh + h] : w.xs[j].data() + Long(h) * L;
            T* ph = p + Long(h) * L;
            qGPU_for(i, L, 1, { ph[i] = xx[i]; });
          }
          mat->multi(ap, p);
          ++multi_count;
          const float sigma = float(shifts[j]);
          qGPU_for(i, N, 1, { r[i] = b[i] - ap[i] - sigma * p[i]; });
          modeled_vector_flops += 6.0 * N;
          auto rn = dot(r, r);
          for (Int h = 0; h < mrh; ++h) {
            true_residual[j * mrh + h] = b2[h] == 0 ? std::sqrt(rn[h]) : std::sqrt(rn[h] / b2[h]);
            ++residual_check_count;
          }
        }
        if (opt.print_mem)
          print_mem_info(refine_iterations > 0 ? "MCG refinement done" : "MCG residual checks done");
      }
      return;
    } else {
      if (!opt.compute_true_res.value_or(false) && opt.refine_tolerance <= 0)
        return;
      w.prepare_refinement(N);
      ComplexD *dr = w.refine_r.data(), *dp = w.refine_p.data(), *da = w.refine_ap.data();
      auto precise = [&](ComplexD* out, ComplexD* in) {
        w.profile.run("precise_dslash", [&] {
          if constexpr (requires { mat->multi_precise(out, in); })
            mat->multi_precise(out, in);
          else if constexpr (sizeof(T) == sizeof(ComplexD)) {
            qGPU_for(i, N, 1, { p[i] = in[i]; });
            mat->multi(ap, p);
            qGPU_for(i, N, 1, { out[i] = ap[i]; });
          } else
            Qassert(false); // Float matrices must supply a precise operator.
        });
        ++multi_count;
      };
      const Int limit = opt.refine_maxiter >= 0 ? opt.refine_maxiter : maxiter;
      ScopedQaccNumThreads refinement_threads(opt.qacc_threads_double);
      w.profile.run("residual_checks_refinement", [&] {
        for (size_t j = 0; j < J; ++j) {
          const double sigma = shifts[j];
          Int used = 0;
          while (true) {
            for (Int h = 0; h < mrh; ++h) {
              const ComplexD* xx = Borrowed ? w.solution_rows[j * mrh + h] : w.xs[j].data() + Long(h) * L;
              ComplexD* ph = dp + Long(h) * L;
              qGPU_for(i, L, 1, { ph[i] = xx[i]; });
            }
            precise(da, dp);
            qGPU_for(i, N, 1, {
              dr[i] = b[i] - da[i] - sigma * dp[i];
              dp[i] = dr[i];
            });
            modeled_vector_flops += 6.0 * N;
            rr = dot(dr, dr, true);
            bool all = true;
            for (Int h = 0; h < mrh; ++h) {
              const double target =
                  opt.refine_tolerance > 0 ? std::min(opt.refine_tolerance, stop[j]) : stop[j];
              true_residual[j * mrh + h] = b2[h] == 0 ? std::sqrt(rr[h]) : std::sqrt(rr[h] / b2[h]);
              ++residual_check_count;
              done[h] = b2[h] == 0 || true_residual[j * mrh + h] <= target;
              if (!done[h])
                all = false;
              else {
                ComplexD* ph = dp + h * L;
                qGPU_for(i, L, 1, { ph[i] = 0.0; });
              }
            }
            if (opt.refine_tolerance <= 0 || all || used >= limit)
              break;
            for (Int pass = 0; pass < 32 && used < limit; ++pass) {
              precise(da, dp);
              qGPU_for(i, N, 1, { da[i] += sigma * dp[i]; });
              modeled_vector_flops += 4.0 * N;
              auto pap = dot(dp, da, true);
              for (Int h = 0; h < mrh; ++h)
                if (!done[h]) {
                  Qassert(pap[h] > 0 && std::isfinite(pap[h]));
                  const double a = rr[h] / pap[h];
                  ComplexD* xh = Borrowed ? w.solution_rows[j * mrh + h] : w.xs[j].data() + Long(h) * L;
                  ComplexD *rh = dr + h * L, *ph = dp + h * L, *ah = da + h * L;
                  qGPU_for(i, L, 1, {
                    xh[i] += a * ph[i];
                    rh[i] -= a * ah[i];
                  });
                  modeled_vector_flops += 8.0 * L;
                }
              auto rn = dot(dr, dr, true);
              bool finished = true;
              for (Int h = 0; h < mrh; ++h)
                if (!done[h]) {
                  const double beta = rn[h] / rr[h], target = std::min(opt.refine_tolerance, stop[j]);
                  done[h] = rn[h] <= target * target * b2[h];
                  if (!done[h])
                    finished = false;
                  const double keep = done[h] ? 0 : 1;
                  ComplexD *rh = dr + h * L, *ph = dp + h * L;
                  qGPU_for(i, L, 1, { ph[i] = keep * (rh[i] + beta * ph[i]); });
                  modeled_vector_flops += 6.0 * L;
                }
              rr = rn;
              ++used;
              ++refine_iterations;
              if (finished)
                break;
            }
          }
        }
      });
      if (opt.print_mem)
        print_mem_info(refine_iterations > 0 ? "MCG refinement done" : "MCG residual checks done");
    }
  }
};
} // namespace qlat
#endif
