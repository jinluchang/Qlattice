#ifndef UTILS_STAGGERED_MULTISHIFT_MCG_H
#define UTILS_STAGGERED_MULTISHIFT_MCG_H
#pragma once

/*
 * Staggered multi-shift CG: controls, diagnostics and benchmark reference
 * ===================================================================
 * Public offsets are 4*m*m; output index = shift*number_of_rhs + rhs.
 * Precision follows qinv: D/D = double I/O and Krylov, D/F = double I/O
 * with float Krylov, F/F = float I/O and Krylov. Recurrence scalars remain
 * double for the custom backend. No control below silently changes qinv.
 *
 * Environment control                 Default  Meaning / affected backends
 * ----------------------------------  -------  ---------------------------
 * qlat_quda_mcg_type                   0        Backend: 0 custom MCG,
 *   1 QUDA multishift, 2 independent QUDA CG, 3 Richardson, 4 heavy-ball,
 *   5 standalone typed CG, 6 standalone typed multishift CG. Types 5/6
 *   require double public vectors and one participating GPU (DD/D/F).
 * qlat_multi_refine                    0        Custom: 0 no correction;
 *   1/2 apply QUDA CG corrections, including required residual checks.
 *   QUDA backend 1 retains its existing meanings. Types 2-6 bypass it.
 * qlat_mcg_compute_true_res            0        Custom check-only pass, 0/1.
 *   No redundant pre-correction pass when custom refinement is enabled.
 * qlat_multi_refine_iter               1000     QUDA correction limit >= 0;
 *   zero skips correction and creates no correction context.
 * qlat_mcg_validate_layout             0        0/1 automatic GPU layout
 *   probes for types 0 and 2-4. Cheap contract checks are always active.
 *   Explicit validate_staggered_mcg_ordering*() still requests the probe.
 * qlat_mcg_signature                   0        Any nonzero value enables
 *   gauge/source/output fingerprints for types 0 and 2-6. Expensive: host
 *   buffers, device-to-host copies, CPU calculations, MPI sums and logging.
 * qlat_quda_verbos                     0        -1 suppresses MCG logging;
 *   >= 0 and -2 print summaries; >= 2 also records stopping diagnostics;
 *   >= 3 retains exhaustive legacy records. Output/host work has a cost.
 * qlat_mcg_reduction_precision         0        0 double, 1 RealDD products
 *   and accumulation for type 0; the recurrence receives double scalars.
 * qlat_mcg_color_sum                   0        0/1 scalar color-first sums
 *   for type 0; adds per-site scratch traffic and a synchronized pass.
 * qlat_mcg_shift_stop                  1        0 requested target,
 *   1 precision epsilon retirement of higher shifts (types 0/6): FLT_EPSILON
 *   for float Krylov vectors, DBL_EPSILON for double. Not true accuracy.
 *   Fixed-step mode disables target retirement.
 * qlat_mcg_mixed_accumulation          0        Mixed D/F mode 0/1 only:
 *   0 backend default (custom single, QUDA double); 1 force single solution
 *   accumulation; 2 force double. Dslash/Krylov/reductions remain unchanged.
 * qlat_mcg_beta                        0        0 norm ratio, 1 historical
 *   QUDA residual-difference beta (types 0/6); changes recurrence arithmetic.
 * qlat_mcg_single_rhs                  0        0 batched, 1 independent
 *   RHS (types 0/3/4/6). Types 2/5 are always independent. Types 0/6 poll
 *   every original batch-width iterations, even after some RHS finish.
 * qlat_mcg_maxiter                     0        0 preserves caller limit;
 *   positive caps types 0 and 2-6. Type 1 uses caller niter directly.
 * qlat_mcg_max_temporal_batch          0        Type 0 only: positive caps
 *   automatic AND explicit temporal windows in double-history steps; double
 *   the cap when both histories are single. Unset/0 leaves windows unchanged.
 * qlat_mcg_print_mem                   0        0/1 print_mem_info at custom
 *   entry, matrix/workspace setup, optional refinement, and solve exit.
 * qlat_mcg_wait_quda                   0        0 original async behavior;
 *   1 checked QUDA completion at custom handoffs and completion of correction
 *   wrapper copies. Independent of QUDA_ENABLE_P2P and QUDA_ENABLE_GDR.
 * qlat_quda_mcg_type_rho               0.1      Finite nonnegative damping
 *   for types 3/4 only; no automatic spectral stabilization.
 * qlat_mcg_type5_precision             D        DD/D/F (case insensitive)
 *   selects complete standalone arithmetic for types 5/6, not qinv sloppy.
 * qlat_mcg_type5_err                   unset    Finite positive override
 *   for all type-5/6 targets. Squaring tiny values can underflow to zero;
 *   e.g. 1e-200 does not establish an attainable 1e-200 true residual.
 * qlat_quda_stag_mrh                   9        Native QUDA source-group
 *   limit, distinct from custom MCG batching and temporal replay.
 *
 * Lifetimes: backend, verbosity and native source-group limit are read at
 * inverter construction. Refinement mode is read per public call and its
 * limit per correction. Other controls are cached once per process on
 * first relevant use; type-5 controls are first read on type-5/6 use.
 * Restart the process to change cached settings. Existing parsing remains
 * unchanged; the new validate-layout switch strictly accepts 0 or 1.
 * Explicit MSCGOptions beta/retirement policies override the environment.
 *
 * C++ tuning/observation controls (MSCGOptions in utils_multi_shift_cg.h):
 *   fuse_residual_norm=true; vector_block=128 (64/128/256);
 *   shift_batch=6 (1/2/4/6); temporal_batch=0 means max(1, shifts-1);
 *     positive values explicitly override the automatic window;
 *   qacc_threads_single=128, qacc_threads_double=64 (1..1024, scoped);
 *   generic_vector_layout=false (reference accessor); profile=false;
 *   stopping_diagnostics=false (also enabled by qinv verbosity >= 2);
 *   accumulation=inherit/single/double_precision; public D/F mode 0 defaults
 *     to single solution updates and borrows caller output storage;
 *   fixed_iterations=false; compute_true_res unset inherits the custom
 *   process flag (default false); explicit true/false overrides it.
 *   refine_tolerance=0; refine_maxiter=-1 (caller limit); stats=nullptr.
 *   Direct _mcg refinement uses refine_tolerance, independently of the public
 *   dispatcher's qlat_multi_refine. Refinement always performs needed checks.
 * No MCG profiling events or detailed phase timing maps are populated when
 * profile=false. Stopping observers add host work only. Disabled signatures
 * and probes add no diagnostic GPU/MPI work. Ordinary solver reductions,
 * transfers, completion synchronization and summary accounting still run.
 *
 * Default temporal window is exactly max(1, supplied shifts-1), excluding
 * the seed. Single-shift calls use 1 and allocate no temporal history.
 * Positive C++ overrides are accepted without rounding; an explicitly set
 * qlat_mcg_max_temporal_batch caps them along with automatic windows, with
 * twice the cap for F/F and D/F single accumulation. D/D and D/F double
 * accumulation use the exported cap. Explicit step requests are not doubled. 4/8/16
 * retain specialized replay kernels, other windows use runtime replay.
 * Host-only automatic execution uses 1. Negative windows are rejected.
 * History memory grows with window*RHS*vector length and retains the largest
 * capacity until cache cleanup. No silent memory-based fallback is applied.
 * Completed shift/RHS lanes leave replay after their last pending update;
 * empty groups are not launched. History is saved only for RHS with effective
 * updates and stops entirely after the final nonseed flush, including for
 * window 1. Recurrence arithmetic and retirement decisions are unchanged.
 * Ordinary mode-2 corrections retain QUDA CG, source/solution buffers and
 * reduction storage per inverter/geometry/public precision. Sources are
 * refreshed and normalized once per call/RHS, then reused across masses.
 * Cached Dirac objects rebind current resident precise links on each call;
 * no precision/gauge reload or inverter-field rebuild is needed. Capacity
 * grows to the largest RHS count and is released by inverter/cache cleanup.
 * Split-grid, deflation, resident/chronological features, nonnative/padded
 * fields and nonstandard normalization use the legacy correction interface,
 * which retains its setup costs. QUDA still uses internal temporary/pool
 * operations: zero MCG-owned buffer growth is not zero total memory activity.
 * Keep fusion, block 128, shift batch 6, qacc threads 128/64, double
 * reductions and batched RHS. Use an explicit smaller window when needed.
 * Select with two synchronized, unprofiled sweeps at equal accuracy/work:
 * require >= 5% gain in both sweeps, prefer the smallest window within 2%
 * of the fastest, then validate convergence and available memory separately.
 * Audit on RTX 5060 Ti, 16^3*32, 7 shifts, fixed 64, RHS 1/4 selected:
 * D/D window 8, D/F window 4, F/F window 8/16 (19-40% faster than window 1).
 * F/F 24^3*48, RHS 8, 17 shifts selected 16 (921.4 -> 302.4 ms).
 * These equal-work harness tests include a duplicate shift; they do not
 * establish convergence gains or replace workload-specific validation.
 * Corrected RTX 5060 Ti / 24^3*48 / 17-mass F/F fixed-64 tests measured
 * batch 1 -> 16: 155.5 -> 69.4 ms (RHS 1), 987.8 -> 328.0 ms (RHS 8).
 * Earlier FP32 temporal tables used double history and are superseded.
 * Isolated D/D or D/F replay gains are not full-solver convergence gains.
 * See /Shared/mcg_mrh/{kernel_fusion_trial,full_cg_retest,
 * full_cg_fp32_temporal_nosync2_trial}/REPORT.md. Compare synchronized full
 * calls separately from instrumented phases and use the same accuracy work.
 * Unique-mass convergence/refinement baseline for the former automatic default:
 * multishift_benchmark/nonseed_default_20261028/REPORT.md. Mode 2 meets
 * D/D 1e-10, D/F 1e-8 and F/F 1e-5 on the tested workloads. Refinement is now
 * off by default; independently validate uncorrected solves when true accuracy
 * is required.
 * Completed-shift/correction-reuse audit (same GPU, large F/F/RHS 8/17 shifts):
 * mode 0: 3756.2 -> 2929.8 ms; mode 2: 4119.4 -> 3274.8 ms, with the same
 * corrected true residual 9.45e-6 (target 1e-5). Refined QUDA: 8001.3 ms.
 * Replay launches drop 52 -> 30, with no history copies/replays after the
 * final required nonseed flush. Warmed correction traces have no CUDA buffer
 * allocation/free calls; QUDA internal pool/host temporary activity remains.
 * Small D/D cases are largely unchanged; refined one-RHS D/F still favors
 * QUDA. See multishift_benchmark/retired_reuse_20261028/REPORT.md for timings,
 * exact accuracy checks, cache lifetime, sanitizer and disabled-test evidence.
 * Precision-default audit (2026-10-29):
 *   custom D/F converts the source once, accumulates in float in the first
 *   half (bytes) of each caller double output, then stages/widens one mass
 *   at a time through the dead residual buffer. Output rows may be separate
 *   allocations. No additional mass-sized solution buffers are allocated.
 *   QUDA D/F retains double solution accumulation by default. The optional
 *   single override uses QUDA's partial accumulator/reliable-update semantics.
 *   congrad.c uses DBL_EPSILON literally; this adapter follows its squared
 *   norm test/poll cadence with epsilon selected by Krylov precision instead.
 *   The seed keeps the requested target and ends the whole solve; nonseed
 *   shifts retire at epsilon or finish with the seed, whichever occurs first.
 *   RTX 5060 Ti, F/F 24^3*48, RHS 8, 17 unique masses: approximately 3.12 s
 *   custom / 7.74 s QUDA uncorrected (both miss the true 1e-5 target);
 *   corrected 3.43 s / 8.01 s, both meet target (2.34x). D/F at true 1e-8
 *   still favors QUDA after corrections. Full data and reproduction:
 *   multishift_benchmark/precision_defaults_20261029/REPORT.md.
 *
 * Benchmark input keys below belong to control_staggered_multishift_cg.cpp,
 * not the production environment/API. Defaults refer to the current driver.
 *   precision=64 selects double/float Krylov for double I/O; fp32_vectors=0
 *     selects the normal driver, =1 selects all-single fixed-work solves.
 *   mrh=1; masses="0.01 0.02 0.03 0.04"; niter=64; cg_err=1e-30;
 *     fixed_iterations=1; repeats=5; warmup=2. Normal driver additionally
 *     uses refine_err=0, refine_iter=10000 and zero_source=0 (1 mixed,
 *     2 all zero). Use fixed_iterations=0 for convergence measurements.
 *   reduction_check=1, layout_check=1: normal-driver checks outside timing;
 *     gauge_check=0: cached/fresh gauge test; generic_layout=0 selects the
 *     C++ generic_vector_layout option. Tuning keys vector_block=128,
 *     shift_batch=6, fuse_residual_norm=1 and qacc_threads_single/double
 *     =128/64 map to the C++ options. profile=1 adds a SEPARATE profiling
 *     pass; headline repetitions always have C++ profile=false.
 *   temporal_batch=1 in the existing F/F driver overrides the automatic
 *     production default; that driver still accepts only 1/4/8/16. The
 *     normal D/D and D/F driver ignores this key: use MSCGOptions in a
 *     caller/harness to test automatic or arbitrary positive windows.
 *     fp32_memory_light=0 optionally
 *     separates F/F custom/QUDA workspace residency when set to 1.
 *   reduction_only=0, backend_check=0, simple_check=0, type6_check=0:
 *     nonzero selects a dedicated suite instead of the solver benchmark.
 *   seed_simple_trial=0: standalone seed comparison; seed_repeats=21,
 *     seed_warmup=7, seed_batch=16. Retired seed_update_trial,
 *     seed_schedule_trial and seed_streaming_trial are rejected when set.
 *   shift_temporal_trial=0: standalone replay/matrix experiment;
 *     temporal_checks=1, temporal_only_checks=0, temporal_warmup=5,
 *     temporal_repeats=21, temporal_profile=0. Here temporal_batch=0 means
 *     all variants; 4/8/16 selects a replay/matrix pair. This is a separate
 *     flat-vector experiment, not a production CG tuning selector.
 */

#include "utils_quda_inverter.h"
#include "utils_vector_cs.h"
#include "utils_multi_shift_cg.h"
#include "utils_mcg_signature.h"
#include "utils_staggered_simple.h"
#include "utils_multi_shift_cg_typed.h"
#include <cctype>
#include <chrono>
#include <cstdlib>
#include <limits>
#include <map>
#include <memory>
#include <color_spinor_field_order.h>
#include <sstream>
#include <stdexcept>

namespace qlat {

// Process settings and policy resolution
inline bool staggered_mcg_print_mem_enabled() {
  static const bool enabled = [] {
    const Long value = qlat::get_env_long_default("qlat_mcg_print_mem", 0);
    if (value < 0 || value > 1)
      throw std::invalid_argument("qlat_mcg_print_mem must be 0 or 1");
    return value != 0;
  }();
  return enabled;
}

// Optional check-only work. Refinement always retains its required checks.
inline bool staggered_mcg_compute_true_res_enabled() {
  static const bool enabled = [] {
    const Long value = qlat::get_env_long_default("qlat_mcg_compute_true_res", 0);
    if (value < 0 || value > 1)
      throw std::invalid_argument("qlat_mcg_compute_true_res must be 0 or 1");
    return value != 0;
  }();
  return enabled;
}

// Automatic GPU probes are diagnostic work, independent of the inexpensive
// layout contract checks required by every solve. Explicit validate_* calls
// below still request a once-per-public-type probe, regardless of this flag.
inline bool staggered_mcg_layout_probe_enabled() {
  static const bool enabled = [] {
    const Long value = qlat::get_env_long_default("qlat_mcg_validate_layout", 0);
    if (value < 0 || value > 1)
      throw std::invalid_argument("qlat_mcg_validate_layout must be 0 or 1");
    return value != 0;
  }();
  return enabled;
}

struct StaggeredMCGAlgorithmConfig {
  MSCGShiftStopPolicy shift_stop = MSCGShiftStopPolicy::precision_epsilon;
  MSCGBetaPolicy beta = MSCGBetaPolicy::norm_ratio;
  bool single_rhs = false;
  Int max_temporal_batch = 0;
};

inline const StaggeredMCGAlgorithmConfig& get_staggered_mcg_algorithm_config() {
  static const StaggeredMCGAlgorithmConfig config = [] {
    const Long stop = qlat::get_env_long_default("qlat_mcg_shift_stop", 1);
    const Long beta = qlat::get_env_long_default("qlat_mcg_beta", 0);
    const Long single_rhs = qlat::get_env_long_default("qlat_mcg_single_rhs", 0);
    const Long max_temporal_batch = qlat::get_env_long_default("qlat_mcg_max_temporal_batch", 0);
    if (stop < 0 || stop > 1)
      throw std::invalid_argument("qlat_mcg_shift_stop must be 0 (requested) or 1 (precision epsilon)");
    if (beta < 0 || beta > 1)
      throw std::invalid_argument("qlat_mcg_beta must be 0 (norm ratio) or 1 (QUDA modified beta)");
    if (single_rhs < 0 || single_rhs > 1)
      throw std::invalid_argument("qlat_mcg_single_rhs must be 0 (batched RHS) or 1 (single RHS)");
    if (max_temporal_batch < 0 || max_temporal_batch > std::numeric_limits<Int>::max())
      throw std::invalid_argument("qlat_mcg_max_temporal_batch must be 0 (no cap) or a positive value fitting Int");
    StaggeredMCGAlgorithmConfig value;
    value.shift_stop = stop ? MSCGShiftStopPolicy::precision_epsilon : MSCGShiftStopPolicy::requested;
    value.beta = beta ? MSCGBetaPolicy::quda : MSCGBetaPolicy::norm_ratio;
    value.single_rhs = single_rhs != 0;
    value.max_temporal_batch = Int(max_temporal_batch);
    return value;
  }();
  return config;
}

// Custom handoff policy is independent of QUDA's P2P/GDR configuration.
// Unset/zero retains the original asynchronous path.
inline bool staggered_mcg_wait_for_quda_enabled() {
  static const bool enabled = [] {
    const Long value = qlat::get_env_long_default("qlat_mcg_wait_quda", 0);
    if (value < 0 || value > 1)
      throw std::invalid_argument("qlat_mcg_wait_quda must be 0 or 1");
    return value != 0;
  }();
  return enabled;
}

inline void staggered_mcg_wait_for_quda() {
#ifdef QLAT_USE_ACC
  if (staggered_mcg_wait_for_quda_enabled())
    qacc_Errchk(qacc_DeviceSynchronize());
#endif
}

inline bool staggered_mcg_single_accumulation(const MSCGOptions& options) {
  return options.accumulation == MSCGAccumulationPolicy::single ||
         options.accumulation == MSCGAccumulationPolicy::inherit;
}

// Resolve accumulation options before determining the history precision.
// D/F single accumulation borrows float outputs even though public I/O is double.
template <class Ty>
inline Int staggered_mcg_effective_temporal_cap(QudaPrecision krylov, const MSCGOptions& options) {
  const Int cap = get_staggered_mcg_algorithm_config().max_temporal_batch;
  const bool single_histories = krylov == QUDA_SINGLE_PRECISION &&
      (std::is_same_v<Ty, ComplexF> || staggered_mcg_single_accumulation(options));
  if (single_histories && cap > 0)
    return cap > std::numeric_limits<Int>::max() / 2 ? std::numeric_limits<Int>::max() : 2 * cap;
  return cap;
}

template <class Ty>
inline Int resolve_staggered_mcg_temporal_batch(Int requested, size_t shifts, QudaPrecision krylov,
                                                const MSCGOptions& options) {
  // Validate first, retaining host-only and single-shift behavior.
  const Int window = resolve_mscg_temporal_batch(requested, shifts);
  const Int cap = staggered_mcg_effective_temporal_cap<Ty>(krylov, options);
  return cap > 0 ? std::min(window, cap) : window;
}

// Explicit caller policies override process settings; generic solver callers
// never read these environment variables. Resolve before solving AND logging.
inline MSCGOptions resolve_staggered_mcg_options(MSCGOptions options) {
  const auto& config = get_staggered_mcg_algorithm_config();
  options.print_mem = staggered_mcg_print_mem_enabled();
  const bool compute_true_res = staggered_mcg_compute_true_res_enabled();
  if (!options.compute_true_res.has_value())
    options.compute_true_res = compute_true_res;
  if (options.shift_stop == MSCGShiftStopPolicy::inherit)
    options.shift_stop = config.shift_stop;
  if (options.beta == MSCGBetaPolicy::inherit)
    options.beta = config.beta;
  const auto accumulation = get_mcg_mixed_accumulation();
  if (options.accumulation == MSCGAccumulationPolicy::inherit && accumulation)
    options.accumulation = accumulation == 1 ? MSCGAccumulationPolicy::single
                                              : MSCGAccumulationPolicy::double_precision;
  return options;
}

// Keep native double and Krylov indices independent. The generic path uses
// QUDA's accessor, including native packing introduced by later QUDA versions.
// Native layout contracts and optional startup probes
template <class T> struct staggered_mcg_layout {
  static constexpr bool krylov_arithmetic = true;
  using Real = std::conditional_t<sizeof(T) == sizeof(ComplexD), double, float>;
  static constexpr Int double_order = quda::colorspinor::get_vector_order<double>(6);
  static constexpr Int storage_order = quda::colorspinor::get_vector_order<Real>(6);
  quda::colorspinor::AccessorCB<double, 1, 3, 1, QUDA_NATIVE_FIELD_ORDER> d;
  quda::colorspinor::AccessorCB<Real, 1, 3, 1, QUDA_NATIVE_FIELD_ORDER> t;
  bool generic = false;

  qacc Long double_index(Long site, Int color, Long sites) const {
    if constexpr (double_order == 2 || double_order == 4)
      if (!generic)
        return mscg_native_index<double_order>(site, color, sites);
    return d.index(0, site, 0, color, 0, sites);
  }

  qacc Long storage_index(Long site, Int color, Long sites) const {
    if constexpr (storage_order == 2 || storage_order == 4)
      if (!generic)
        return mscg_native_index<storage_order>(site, color, sites);
    return t.index(0, site, 0, color, 0, sites);
  }

  template <class U> qacc Long public_index(Long site, Int color, Long sites) const {
    if constexpr (sizeof(U) == sizeof(T))
      return storage_index(site, color, sites);
    else
      return double_index(site, color, sites);
  }
};

inline void check_staggered_mcg_param(const quda::ColorSpinorParam& p, QudaPrecision precision, Long sites) {
  auto fail = [&](const char* reason) {
    std::ostringstream out;
    out << "MCG ordering validation: " << reason << " (field_order=" << int(p.fieldOrder)
        << ", precision=" << int(p.Precision()) << ", configured D/F=" << QUDA_ORDER_DOUBLE << "/"
        << QUDA_ORDER_SINGLE << ")";
    throw std::runtime_error(out.str());
  };
  if (p.fieldOrder != QUDA_NATIVE_FIELD_ORDER)
    fail("unsupported/unknown device field order");
  if (p.Precision() != precision || p.location != QUDA_CUDA_FIELD_LOCATION)
    fail("unexpected precision/location");
  if (p.nSpin != 1 || p.nColor != 3 || p.nVec != 1 || p.nVec_actual != 1)
    fail("expected one staggered vector");
  if (p.siteSubset != QUDA_PARITY_SITE_SUBSET || p.siteOrder != QUDA_EVEN_ODD_SITE_ORDER)
    fail("unexpected parity/site ordering");
  if (p.pad != 0 || p.is_composite || p.is_component)
    fail("padded/composite fields are unsupported");
  if (p.nDim != 4)
    fail("expected four-dimensional staggered geometry");
  Long volume = 1;
  for (Int j = 0; j < p.nDim; ++j) {
    if (p.x[j] <= 0)
      fail("invalid field extent");
    volume *= p.x[j];
  }
  if (volume != sites)
    fail("field geometry does not match qinv parity volume");
}

// Expensive double-input native-layout self-test, including mixed D/F.
inline void probe_staggered_mcg_ordering(quda_inverter& qinv, Int mrh) {
  const Long sites = qinv.geo().local_volume() / 2, L = sites * 3, N = L * mrh;
  const auto start = std::chrono::steady_clock::now();
  MCGDebugBuffer<ComplexD> nativeD(N);
  MCGDebugBuffer<ComplexF> nativeF(N), converted(N);
  std::vector<ComplexD> host(L), back(L), rawD(N);
  std::vector<ComplexF> rawF(N), rawConverted(N);
  staggered_mcg_layout<ComplexF> layout;
  auto value = [&](Int h, Long site, Int color) {
    const Long index = Long(h) * L + site * 3 + color;
    // Pair of exactly representable floats, unique throughout the buffer.
    return ComplexD(double(index % 1048576 + 1), -double(index / 1048576 + 1));
  };
  auto mismatch = [&](const char* stage, Int h, Long site, Int color, ComplexD expected, ComplexD actual) {
    std::ostringstream out;
    out << "MCG ordering validation failed at " << stage << ", rhs=" << h << ", site=" << site
        << ", color=" << color << ", expected=" << expected << ", actual=" << actual
        << ", configured D/F=" << QUDA_ORDER_DOUBLE << '/' << QUDA_ORDER_SINGLE
        << ", effective D/F=" << layout.double_order << '/' << layout.storage_order;
    throw std::runtime_error(out.str());
  };
  for (Int h = 0; h < mrh; ++h) {
    for (Long site = 0; site < sites; ++site)
      for (Int color = 0; color < 3; ++color)
        host[site * 3 + color] = value(h, site, color);
    auto hp = qinv.cs_gpuD;
    hp.location = QUDA_CPU_FIELD_LOCATION;
    hp.mem_type = QUDA_MEMORY_HOST;
    hp.fieldOrder = QUDA_SPACE_SPIN_COLOR_FIELD_ORDER;
    hp.create = QUDA_REFERENCE_FIELD_CREATE;
    hp.v = host.data();
    quda::ColorSpinorField hf(hp);
    hp.v = back.data();
    quda::ColorSpinorField output(hp);
    auto dp = qinv.cs_gpuD, fp = qinv.cs_gpuF;
    dp.create = fp.create = QUDA_REFERENCE_FIELD_CREATE;
    dp.v = nativeD.data() + h * L;
    fp.v = nativeF.data() + h * L;
    quda::ColorSpinorField df(dp), ff(fp);
    // copy() uses Bytes(). Reject alignment padding before a tightly packed
    // RHS buffer can be overrun; reference construction itself never writes.
    if (df.Bytes() != size_t(L) * sizeof(ComplexD) || ff.Bytes() != size_t(L) * sizeof(ComplexF) ||
        df.Stride() != size_t(sites) || ff.Stride() != size_t(sites))
      throw std::runtime_error(
          "MCG ordering validation: native field bytes/stride exceed tightly packed RHS buffers");
    auto check_raw = [&](const char* stage, bool single) {
      if (single)
        cpy_GPU(rawF.data() + h * L, nativeF.data() + h * L, L, 0, 1);
      else
        cpy_GPU(rawD.data() + h * L, nativeD.data() + h * L, L, 0, 1);
      for (Long site = 0; site < sites; ++site)
        for (Int color = 0; color < 3; ++color) {
          const ComplexD expected = value(h, site, color);
          const ComplexD actual = single ? ComplexD(rawF[h * L + layout.storage_index(site, color, sites)])
                                         : rawD[h * L + layout.double_index(site, color, sites)];
          if (actual != expected)
            mismatch(stage, h, site, color, expected, actual);
        }
    };
    df.copy(hf);
    ff.copy(hf);
    check_raw("host-to-native-double", false);
    check_raw("host-to-native-float", true);
    for (const auto* field : {&df, &ff}) {
      output.copy(*field);
      for (Long site = 0; site < sites; ++site)
        for (Int color = 0; color < 3; ++color)
          if (back[site * 3 + color] != value(h, site, color))
            mismatch("native-to-host", h, site, color, value(h, site, color), back[site * 3 + color]);
    }
    ff.copy(df);
    check_raw("QUDA-double-to-float", true);
    df.copy(ff);
    check_raw("QUDA-float-to-double", false);
  }
  // Compare specialized and generic initialization against QUDA's raw result.
  for (bool generic : {false, true}) {
    layout.generic = generic;
#ifdef QLAT_USE_ACC
    mscg_initialize_vectors<<<dim3((sites + 127) / 128, mrh), 128>>>(
        nativeD.data(), converted.data(), nativeF.data(), static_cast<ComplexF* const*>(nullptr), 1, sites,
        layout);
    Qassert(qacc_GetLastError() == qacc_Success);
    qacc_barrier(dummy);
#else
    for (Int h = 0; h < mrh; ++h)
      for (Long site = 0; site < sites; ++site)
        for (Int color = 0; color < 3; ++color)
          converted[h * L + layout.storage_index(site, color, sites)] =
              ComplexF(nativeD[h * L + layout.double_index(site, color, sites)]);
#endif
    cpy_GPU(rawConverted.data(), converted.data(), N, 0, 1);
    for (Int h = 0; h < mrh; ++h)
      for (Long site = 0; site < sites; ++site)
        for (Int color = 0; color < 3; ++color) {
          const ComplexD actual = ComplexD(rawConverted[h * L + layout.storage_index(site, color, sites)]);
          if (actual != value(h, site, color))
            mismatch(generic ? "MCG-generic-conversion" : "MCG-specialized-conversion", h, site, color,
                     value(h, site, color), actual);
        }
  }
  if (qinv.quda_verbos >= 2)
    qmessage("BENCH_LAYOUT,%d,%d,%d,%d,%d,%.12e,1\n", QUDA_ORDER_DOUBLE, QUDA_ORDER_SINGLE,
             layout.double_order, layout.storage_order, mrh,
             std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count());
}

// All-single input also probes its native public layout and both initialization
// paths, without constructing a double-I/O field or changing qinv precision.
inline void probe_staggered_mcg_ordering_single(quda_inverter& qinv, Int mrh) {
  const Long sites = qinv.geo().local_volume() / 2, L = sites * 3, N = L * mrh;
  const auto start = std::chrono::steady_clock::now();
  MCGDebugBuffer<ComplexF> native(N), converted(N);
  std::vector<ComplexF> host(L), back(L), raw(N);
  staggered_mcg_layout<ComplexF> layout;
  auto value = [&](Int h, Long site, Int color) {
    const Long index = Long(h) * L + site * 3 + color;
    return ComplexF(float(index % 1048576 + 1), -float(index / 1048576 + 1));
  };
  auto check_raw = [&] {
    for (Int h = 0; h < mrh; ++h)
      for (Long site = 0; site < sites; ++site)
        for (Int color = 0; color < 3; ++color)
          if (raw[h * L + layout.storage_index(site, color, sites)] != value(h, site, color))
            throw std::runtime_error("MCG ordering validation: single native layout mismatch");
  };
  for (Int h = 0; h < mrh; ++h) {
    for (Long site = 0; site < sites; ++site)
      for (Int color = 0; color < 3; ++color)
        host[site * 3 + color] = value(h, site, color);
    auto hp = qinv.cs_gpuF;
    hp.location = QUDA_CPU_FIELD_LOCATION;
    hp.mem_type = QUDA_MEMORY_HOST;
    hp.fieldOrder = QUDA_SPACE_SPIN_COLOR_FIELD_ORDER;
    hp.create = QUDA_REFERENCE_FIELD_CREATE;
    hp.v = host.data();
    quda::ColorSpinorField hf(hp);
    hp.v = back.data();
    quda::ColorSpinorField output(hp);
    auto fp = qinv.cs_gpuF;
    fp.create = QUDA_REFERENCE_FIELD_CREATE;
    fp.v = native.data() + h * L;
    quda::ColorSpinorField ff(fp);
    if (ff.Bytes() != size_t(L) * sizeof(ComplexF) || ff.Stride() != size_t(sites))
      throw std::runtime_error(
          "MCG ordering validation: native field bytes/stride exceed tightly packed RHS buffers");
    ff.copy(hf);
    output.copy(ff);
    for (Long i = 0; i < L; ++i)
      if (back[i] != host[i])
        throw std::runtime_error("MCG ordering validation: single native round-trip mismatch");
  }
  cpy_GPU(raw.data(), native.data(), N, 0, 1);
  check_raw();
  for (bool generic : {false, true}) {
    layout.generic = generic;
#ifdef QLAT_USE_ACC
    mscg_initialize_vectors<<<dim3((sites + 127) / 128, mrh), 128>>>(
        native.data(), converted.data(), native.data(), static_cast<ComplexF* const*>(nullptr), 1, sites,
        layout);
    Qassert(qacc_GetLastError() == qacc_Success);
    qacc_barrier(dummy);
#else
    for (Int h = 0; h < mrh; ++h)
      for (Long site = 0; site < sites; ++site)
        for (Int color = 0; color < 3; ++color)
          converted[h * L + layout.storage_index(site, color, sites)] =
              native[h * L + layout.template public_index<ComplexF>(site, color, sites)];
#endif
    cpy_GPU(raw.data(), converted.data(), N, 0, 1);
    check_raw();
  }
  if (qinv.quda_verbos >= 2)
    qmessage("BENCH_LAYOUT,%d,%d,%d,%d,%d,%.12e,1\n", QUDA_ORDER_DOUBLE, QUDA_ORDER_SINGLE,
             layout.double_order, layout.storage_order, mrh,
             std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count());
}

template <class U> inline bool staggered_mcg_ordering_self_test_once(quda_inverter& qinv, Int mrh) {
  bool probed = false;
  static const bool ordering_validated = [&] {
    if constexpr (std::is_same_v<U, ComplexD>)
      probe_staggered_mcg_ordering(qinv, mrh);
    else
      probe_staggered_mcg_ordering_single(qinv, mrh);
    probed = true;
    return true;
  }();
  (void)ordering_validated;
  return probed;
}

// Required per-call contract checks; these do not allocate or launch GPU work.
inline void check_staggered_mcg_ordering(quda_inverter& qinv, Int mrh = 1) {
  const Long sites = qinv.geo().local_volume() / 2;
  if (mrh <= 0 || qinv.inv_param.dirac_order != QUDA_INTERNAL_DIRAC_ORDER ||
      qinv.inv_param.cpu_prec != QUDA_DOUBLE_PRECISION || qinv.inv_param.cuda_prec != QUDA_DOUBLE_PRECISION ||
      (qinv.inv_param.cuda_prec_sloppy != QUDA_SINGLE_PRECISION &&
       qinv.inv_param.cuda_prec_sloppy != QUDA_DOUBLE_PRECISION))
    throw std::runtime_error(
        "MCG ordering validation: expected double native I/O and float/double sloppy precision");
  check_staggered_mcg_param(qinv.cs_gpuD, QUDA_DOUBLE_PRECISION, sites);
  check_staggered_mcg_param(qinv.cs_gpuF, QUDA_SINGLE_PRECISION, sites);
}

inline void check_staggered_mcg_ordering_single(quda_inverter& qinv, Int mrh = 1) {
  const Long sites = qinv.geo().local_volume() / 2;
  if (mrh <= 0 || qinv.inv_param.dirac_order != QUDA_INTERNAL_DIRAC_ORDER ||
      qinv.inv_param.cpu_prec != QUDA_SINGLE_PRECISION || qinv.inv_param.cuda_prec != QUDA_SINGLE_PRECISION ||
      qinv.inv_param.cuda_prec_sloppy != QUDA_SINGLE_PRECISION)
    throw std::runtime_error(
        "MCG ordering validation: expected single native I/O and single sloppy precision");
  check_staggered_mcg_param(qinv.cs_gpuF, QUDA_SINGLE_PRECISION, sites);
}

// Explicit validation requests the expensive probe once per public type and
// process, including across inverters, RHS changes and cache clears. Return
// true only for the call that actually performed it.
inline bool validate_staggered_mcg_ordering(quda_inverter& qinv, Int mrh = 1) {
  check_staggered_mcg_ordering(qinv, mrh);
  return staggered_mcg_ordering_self_test_once<ComplexD>(qinv, mrh);
}

inline bool validate_staggered_mcg_ordering_single(quda_inverter& qinv, Int mrh = 1) {
  check_staggered_mcg_ordering_single(qinv, mrh);
  return staggered_mcg_ordering_self_test_once<ComplexF>(qinv, mrh);
}

template <class U> inline bool check_staggered_mcg_ordering_for_solve(quda_inverter& qinv, Int mrh) {
  if constexpr (std::is_same_v<U, ComplexD>)
    check_staggered_mcg_ordering(qinv, mrh);
  else
    check_staggered_mcg_ordering_single(qinv, mrh);
  if (!staggered_mcg_layout_probe_enabled())
    return false;
  return staggered_mcg_ordering_self_test_once<U>(qinv, mrh);
}

// Local matrix follows utils_stagger_eigen.h's MRH Dslash/DslashXpay sequence.
// quda_inverter manages resident links and Dirac operators. This adapter owns
// only spinor views/workspace; refresh() prepares the inverter before a solve.
// QUDA matrix adapter, reusable fields and cache lifetime
template <class T> class staggered_mcg_matrix {
  std::unique_ptr<staggered_mcg_matrix<ComplexD>> precise_matrix;
  std::vector<std::vector<quda::ColorSpinorField>> fields;
  Int mrh = 0;
  Long length = 0;
  bool precise;
  // A ComplexF Krylov matrix can be used either by the all-single public
  // path (F/F) or by the mixed double-I/O path (D/F).  Keep that distinction
  // explicit because the ColorSpinor layout contract is different: D/F still
  // requires QUDA's double native field and single sloppy field.
  bool native_double;

public:
  static constexpr bool staggered_reductions = true;
  quda_inverter* qinvP;
  std::vector<MCGDebugBuffer<T>> buf;
  Long allocations = 0;

  staggered_mcg_matrix(quda_inverter* q, Int rhs, bool precise_ = false, bool native_double_ = false)
      : precise(precise_), native_double(native_double_), qinvP(q) {
    if constexpr (sizeof(T) == sizeof(ComplexD)) {
      check_staggered_mcg_ordering_for_solve<ComplexD>(*qinvP, rhs);
    } else if (native_double) {
      check_staggered_mcg_ordering_for_solve<ComplexD>(*qinvP, rhs);
    } else {
      check_staggered_mcg_ordering_for_solve<ComplexF>(*qinvP, rhs);
    }
    set_mrh(rhs);
    refresh();
  }

  ~staggered_mcg_matrix() {
    // Reference fields must be released before their backing buffers.
    fields.clear();
  }

  void refresh(bool prepare = true) {
    MCGDebugPhase phase("adapter_refresh");
    const auto& param = qinvP->inv_param;
    const QudaPrecision expected =
        sizeof(T) == sizeof(ComplexD) ? QUDA_DOUBLE_PRECISION : QUDA_SINGLE_PRECISION;
    Qassert((precise ? param.cuda_prec : param.cuda_prec_sloppy) == expected);
    if (prepare) {
      MCGDebugPhase gauge_phase("resident_gauge_prepare");
      qinvP->prepare_staggered_dslash();
    }
    if (precise_matrix)
      precise_matrix->refresh(false);
  }

  staggered_mcg_layout<T> vector_layout(bool generic = false) const {
    staggered_mcg_layout<T> result;
    result.generic = generic;
    return result;
  }

  Int get_mrh() const {
    return mrh;
  }

  Long get_vsize() const {
    return length * mrh;
  }

  void set_mrh(Int rhs) {
    MCGDebugPhase phase("matrix_buffers");
    Qassert(rhs > 0);
    const Long l = qinvP->geo().local_volume() * 3 / 2;
    if (mrh == rhs && length == l)
      return;
    mrh = rhs;
    length = l;
    fields.clear();
    buf.resize(3);
    fields.resize(3);
    quda::ColorSpinorParam param = sizeof(T) == sizeof(ComplexD) ? qinvP->cs_gpuD : qinvP->cs_gpuF;
    param.create = QUDA_REFERENCE_FIELD_CREATE;
    for (Int j = 0; j < 3; ++j) {
      buf[j].resize(length * mrh);
      ++allocations;
      buf[j].debug_name("matrix."+std::to_string(j));
      mcg_debug_device_pointer(buf[j].data(),"matrix buffer");
      fields[j].resize(mrh);
      for (Int h = 0; h < mrh; ++h) {
        param.v = buf[j].data() + h * length;
        fields[j][h] = quda::ColorSpinorField(param);
        mcg_debug_require(fields[j][h].data()==param.v &&
          size_t(h+1)*size_t(length)<=buf[j].size(),"QUDA field RHS offset/range");
        if (mcg_debug_active()) {
          auto& d=mcg_debug_state();d.open();
          std::fprintf(d.events,"MCG_FIELD solve=%ld buffer=%d rhs=%d address=%p offset=%zu bytes=%zu stride=%zu\n",
            d.solve,j,h,param.v,size_t(h)*size_t(length)*sizeof(T),fields[j][h].Bytes(),fields[j][h].Stride());
        }
        if (fields[j][h].Bytes() != size_t(length) * sizeof(T) || fields[j][h].Stride() != size_t(length / 3))
          throw std::runtime_error(
              "MCG ordering validation: native field bytes/stride exceed tightly packed RHS buffers");
      }
    }
    if (precise_matrix)
      precise_matrix->set_mrh(rhs);
  }

  void set_buf_pointers(vector<T*>& pointers, Int rhs) {
    set_mrh(rhs);
    pointers.resize(3);
    for (Int j = 0; j < 3; ++j)
      pointers[j] = buf[j].data();
  }

  void multi(T* out, T* in) {
    Qassert(out == buf[1].data() && in == buf[0].data());
    auto& dirac = qinvP->staggered_dslash(precise);
    if (mcg_debug_active()) ++mcg_debug_state().operators;
    {
      MCGDebugPhase phase("quda_dslash_odd");
      dirac.Dslash(fields[2], fields[0], QUDA_ODD_PARITY);
    }
    {
      MCGDebugPhase phase("quda_dslash_even");
      dirac.DslashXpay(fields[1], fields[2], QUDA_EVEN_PARITY, fields[0], 0.0);
    }
    // Optional explicit completion before returning to custom kernels/copies.
    {
      MCGDebugPhase phase("quda_completion_handoff");
      staggered_mcg_wait_for_quda();
    }
  }

  void multi_seed(T* out, T* in, double sigma0) {
    Qassert(out == buf[1].data() && in == buf[0].data());
    Qassert(std::isfinite(sigma0) && sigma0 >= 0);
    // Match QUDA's staggered multi-shift driver, including the sqrt/multiply
    // rounding of 4*mass*mass, without changing the inverter-owned mass.
    const double mass = std::sqrt(sigma0 / 4.0);
    auto& dirac = qinvP->staggered_dslash(precise);
    if (mcg_debug_active()) ++mcg_debug_state().operators;
    {
      MCGDebugPhase phase("quda_dslash_odd");
      dirac.Dslash(fields[2], fields[0], QUDA_ODD_PARITY);
    }
    {
      MCGDebugPhase phase("quda_dslash_even");
      dirac.DslashXpay(fields[1], fields[2], QUDA_EVEN_PARITY, fields[0], 4 * mass * mass);
    }
    {
      MCGDebugPhase phase("quda_completion_handoff");
      staggered_mcg_wait_for_quda();
    }
  }

  void multi_precise(ComplexD* out, ComplexD* in) {
    MCGDebugPhase phase("precise_operator_conversion");
    if (!precise_matrix)
      precise_matrix = std::make_unique<staggered_mcg_matrix<ComplexD>>(qinvP, mrh, true);
    auto& m = *precise_matrix;
    const Long n = get_vsize();
    ComplexD* p = m.buf[0].data();
    ComplexD* ap = m.buf[1].data();
    qGPU_for(i, n, 1, { p[i] = in[i]; });
    m.multi(ap, p);
    qGPU_for(i, n, 1, { out[i] = ap[i]; });
  }

  Long allocation_count() const {
    return allocations + (precise_matrix ? precise_matrix->allocation_count() : 0);
  }
};

template <class T>
inline std::map<std::pair<quda_inverter*, Int>, std::shared_ptr<staggered_mcg_matrix<T>>>&
staggered_mcg_matrix_cache() {
  static std::map<std::pair<quda_inverter*, Int>, std::shared_ptr<staggered_mcg_matrix<T>>> c;
  return c;
}

inline void clear_staggered_mcg_correction_cache(quda_inverter* qinv = nullptr);

// Clear before inverter destruction/replacement. Gauge updates need no clearing:
// each adapter entry prepares the inverter's operators against current links.
inline void clear_staggered_multishift_mcg_cache() {
  MCGDebugMaintenance activity("cache_cleanup");
  MCGDebugPhase phase("cache_cleanup");
  clear_staggered_mcg_correction_cache();
  staggered_mcg_matrix_cache<ComplexD>().clear();
  staggered_mcg_matrix_cache<ComplexF>().clear();
  get_mscg_cache<ComplexD>().clear();
  get_mscg_cache<ComplexF>().clear();
  get_mscg_cache<ComplexF, ComplexF>().clear();
}

// Release only the cached QUDA matrix wrappers owned by one inverter.  The
// wrappers contain reference ColorSpinorFields; leaving them
// in the process-wide cache until endQuda makes QUDA report their host-side
// object allocations as leaks. The workspace cache contains
// only qlat buffers and is intentionally retained for reuse by other
// inverters.
inline void clear_staggered_multishift_mcg_cache(quda_inverter* qinv) {
  MCGDebugMaintenance activity("cache_cleanup");
  MCGDebugPhase phase("cache_cleanup");
  clear_staggered_mcg_correction_cache(qinv);
  auto erase = [qinv](auto& cache) {
    for (auto it = cache.begin(); it != cache.end();) {
      if (it->first.first == qinv)
        it = cache.erase(it);
      else
        ++it;
    }
  };
  erase(staggered_mcg_matrix_cache<ComplexD>());
  erase(staggered_mcg_matrix_cache<ComplexF>());
}

inline void install_staggered_multishift_mcg_cache_cleanup(quda_inverter& qinv) {
  if (!qinv.staggered_mcg_cache_cleanup)
    qinv.staggered_mcg_cache_cleanup = [&qinv] { clear_staggered_multishift_mcg_cache(&qinv); };
}

template <class T, class U>
inline void staggered_multishift_mcg_dispatch(quda_inverter& qinv, const std::vector<U*>& res,
                                              const std::vector<U*>& src, const std::vector<double>& masses2,
                                              double err, Int niter, const MSCGOptions& options) {
  install_staggered_multishift_mcg_cache_cleanup(qinv);
  const auto check_start =
      options.profile ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
  const bool probed = check_staggered_mcg_ordering_for_solve<U>(qinv, src.size());
  const double check_time =
      options.profile ? std::chrono::duration<double>(std::chrono::steady_clock::now() - check_start).count()
                      : 0;
  const auto start =
      options.profile ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
  auto& cache = staggered_mcg_matrix_cache<T>();
  const auto key = std::make_pair(&qinv, Int(src.size()));
  auto it = cache.find(key);
  Long before = it == cache.end() ? 0 : it->second->allocation_count();
  const bool create = it == cache.end();
  {
    MCGDebugPhase phase(create ? "adapter_create" : "adapter_reuse");
    if (create) {
      // U selects the public vector precision.  D/F therefore has T=ComplexF
      // but still uses double native I/O; pass that fact to the matrix adapter
      // so its layout validation checks the correct QUDA fields.
      constexpr bool native_double = std::is_same_v<U, ComplexD>;
      it =
          cache.emplace(key, std::make_shared<staggered_mcg_matrix<T>>(&qinv, src.size(), false, native_double))
              .first;
    } else {
      it->second->set_mrh(src.size());
      it->second->refresh();
    }
  }
  const double matrix_time =
      options.profile ? std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() : 0;
  if (options.print_mem)
    print_mem_info("MCG matrix ready");
  multi_shift_cg<staggered_mcg_matrix<T>, vector_cs<U>> solver(*it->second);
  if constexpr (std::is_same_v<U, ComplexD> && std::is_same_v<T, ComplexF>) {
    if (staggered_mcg_single_accumulation(options))
      solver.call_raw_mixed(res, src, masses2, err, niter, options);
    else
      solver.call_raw(res, src, masses2, err, niter, options);
  } else
    solver.call_raw(res, src, masses2, err, niter, options);
  if (options.stats) {
    if (options.profile) {
      options.stats->host_seconds[probed ? "ordering_validation" : "ordering_param_check"] = check_time;
      options.stats->host_seconds[create ? "matrix_create" : "matrix_refresh"] = matrix_time;
    }
    options.stats->matrix_allocations = it->second->allocation_count() - before;
    // Two parity Dslashes: 570 (staggered) or 1146 (Asqtad) per full-volume
    // site/RHS, excluding reconstruction and communication arithmetic.
    const double site_flops = qinv.inv_param.dslash_type == QUDA_ASQTAD_DSLASH ? 1146.0 : 570.0;
    options.stats->modeled_operator_flops =
        site_flops * qinv.geo().local_volume() * src.size() * options.stats->matrix_calls;
  }
}

// Combine independent single-RHS solves without changing the public shift-major
// indexing. Internal iterations remain the maximum; headline work is summed.
inline void merge_staggered_mcg_rhs_stats(MSCGStats& total, const MSCGStats& one, size_t rhs, size_t h) {
  total.iterations = std::max(total.iterations, one.iterations);
  total.refine_iterations += one.refine_iterations;
  total.matrix_calls += one.matrix_calls;
  total.reductions += one.reductions;
  total.seed_updates += one.seed_updates;
  total.shift_updates += one.shift_updates;
  total.residual_checks += one.residual_checks;
  total.workspace_allocations += one.workspace_allocations;
  total.matrix_allocations += one.matrix_allocations;
  total.modeled_vector_flops += one.modeled_vector_flops;
  total.modeled_operator_flops += one.modeled_operator_flops;
  total.solve_seconds += one.solve_seconds;
  total.borrowed_solutions = one.borrowed_solutions;
  total.accumulation_bytes = one.accumulation_bytes;
  total.retirement_epsilon = one.retirement_epsilon;
  for (const auto& [name, seconds] : one.host_seconds)
    total.host_seconds[name] += seconds;
  for (const auto& [name, seconds] : one.device_seconds)
    total.device_seconds[name] += seconds;
  auto scatter = [rhs, h](auto& out, const auto& in) {
    out.resize(in.size() * rhs);
    for (size_t j = 0; j < in.size(); ++j)
      out[j * rhs + h] = in[j];
  };
  scatter(total.shift_iterations, one.shift_iterations);
  scatter(total.iterative_residual, one.iterative_residual);
  scatter(total.true_residual, one.true_residual);
  auto& d = total.stopping;
  const auto& s = one.stopping;
  d.maxiter = s.maxiter;
  scatter(d.source_norm2, s.source_norm2);
  scatter(d.seed_threshold2, s.seed_threshold2);
  scatter(d.shift_threshold2, s.shift_threshold2);
  scatter(d.seed_checks, s.seed_checks);
  scatter(d.shift_checks, s.shift_checks);
  scatter(d.seed_reason, s.seed_reason);
  scatter(d.shift_reason, s.shift_reason);
  scatter(d.shift_retired, s.shift_retired);
  scatter(d.effective_updates, s.effective_updates);
}

// Preserve only inputs overlapped by outputs, for work spanning multiple raw
// solves (split RHS or external correction). Ordinary disjoint calls allocate nothing.
template <class Ty>
inline void preserve_staggered_mcg_aliased_sources(const std::vector<Ty*>& res, std::vector<Ty*>& src,
                                                  Long length,
                                                  std::vector<MCGDebugBuffer<Ty>>& saved) {
  const uintptr_t bytes = size_t(length) * sizeof(Ty);
  saved.resize(src.size());
  for (size_t h = 0; h < src.size(); ++h) {
    const uintptr_t a = reinterpret_cast<uintptr_t>(src[h]);
    bool overlap = false;
    for (const auto* row : res) {
      const uintptr_t b = reinterpret_cast<uintptr_t>(row);
      overlap = overlap || (a <= b ? b - a < bytes : a - b < bytes);
    }
    if (overlap) {
      saved[h].resize(length);
      const Ty* input = src[h];
      Ty* copy = saved[h].data();
      qGPU_for(i, length, 1, { copy[i] = input[i]; });
      src[h] = copy;
    }
  }
}

// Implementation shared by both public entries; neither nests a public logger.
// Optional experiment: split public sources into independent single-RHS solves.
template <class Ty>
void staggered_multishift_even_mcg_impl(quda_inverter& qinv, std::vector<Ty*> res, std::vector<Ty*> src,
                                        const std::vector<double>& masses2, double err, Int niter,
                                        const MSCGOptions& options = {}) {
  static_assert(std::is_same<Ty, ComplexD>::value || std::is_same<Ty, ComplexF>::value,
                "MCG public input/output must be ComplexD or ComplexF");
  MCGDebugSolve debug(src.size(),masses2.size(),sizeof(Ty),int(qinv.inv_param.cuda_prec_sloppy)*2,
    0);
  mcg_debug_public_buffers(res,src,masses2.size(),qinv.geo().local_volume()*3/2);
  if (mcg_debug_active())
    mcg_debug_state().window=resolve_mscg_temporal_batch(options.temporal_batch,masses2.size());
  Qassert(!src.empty() && !masses2.empty() && res.size() == src.size() * masses2.size());
  Qassert(qinv.inv_param.use_sloppy_partial_accumulator == 0 && qinv.inv_param.dagger == QUDA_DAG_NO);
  if (get_staggered_mcg_algorithm_config().single_rhs && src.size() > 1) {
    std::vector<MCGDebugBuffer<Ty>> saved_sources;
    preserve_staggered_mcg_aliased_sources(res, src, qinv.geo().local_volume()*3/2, saved_sources);
    if (options.stats)
      *options.stats = {};
    std::vector<Ty*> single_res(masses2.size());
    for (size_t h = 0; h < src.size(); ++h) {
      MCGDebugSubcall debug_subcall(1,h);
      for (size_t j = 0; j < masses2.size(); ++j)
        single_res[j] = res[j * src.size() + h];
      MSCGStats one;
      MSCGOptions single = options;
      single.stats = options.stats ? &one : nullptr;
      staggered_multishift_even_mcg_impl(qinv, single_res, std::vector<Ty*>{src[h]}, masses2, err, niter,
                                         single);
      if (options.stats)
        merge_staggered_mcg_rhs_stats(*options.stats, one, src.size(), h);
    }
    return;
  }
  if constexpr (std::is_same<Ty, ComplexD>::value) {
    Qassert(qinv.inv_param.cpu_prec == QUDA_DOUBLE_PRECISION &&
            qinv.inv_param.cuda_prec == QUDA_DOUBLE_PRECISION);
    Qassert(qinv.inv_param.cuda_prec_sloppy == QUDA_DOUBLE_PRECISION ||
            qinv.inv_param.cuda_prec_sloppy == QUDA_SINGLE_PRECISION);
    if (qinv.inv_param.cuda_prec_sloppy == QUDA_DOUBLE_PRECISION)
      staggered_multishift_mcg_dispatch<ComplexD, ComplexD>(qinv, res, src, masses2, err, niter, options);
    else
      staggered_multishift_mcg_dispatch<ComplexF, ComplexD>(qinv, res, src, masses2, err, niter, options);
  } else {
    Qassert(qinv.inv_param.cpu_prec == QUDA_SINGLE_PRECISION &&
            qinv.inv_param.cuda_prec == QUDA_SINGLE_PRECISION &&
            qinv.inv_param.cuda_prec_sloppy == QUDA_SINGLE_PRECISION);
    staggered_multishift_mcg_dispatch<ComplexF, ComplexF>(qinv, res, src, masses2, err, niter, options);
  }
}

// Host diagnostics and public-vector fingerprints
inline const char* staggered_mcg_precision_name(QudaPrecision p) {
  return p == QUDA_DOUBLE_PRECISION ? "double" : p == QUDA_SINGLE_PRECISION ? "single" : "other";
}

// Shared across public-vector types; IDs count diagnostic custom completions.
inline Long next_staggered_mcg_diagnostic_id() {
  static Long id = 0;
  return ++id;
}

inline void emit_staggered_mcg_signature_gauge(quda_inverter& qinv, Long call) {
  const auto& links = qinv.qinv_twist_link == QTRUE ? qinv.quda_gf_default_twist : qinv.quda_gf_default;
  const auto& geo = qinv.geo();
  const Coordinate dims = geo.total_site();
  const MCGSignature sig = mcg_signature_gauge(geo, links.data(), links.size());
  qmessage("  MCG_GAUGE_SIG call=%ld source=%s elements=%ld dims=%dx%dx%dx%d l2=%.17e coord_l2=%.17e "
           "probe0=%.17e probe1=%.17e probe2=%.17e\n",
           long(call), qinv.qinv_twist_link == QTRUE ? "qinv.quda_gf_default_twist" : "qinv.quda_gf_default",
           long(links.size()), dims[0], dims[1], dims[2], dims[3], sig.l2, sig.coord_l2, sig.probe0,
           sig.probe1, sig.probe2);
}

template <class Ty>
inline void emit_staggered_mcg_signature_sources(quda_inverter& qinv, const std::vector<Ty*>& src,
                                                 Long call) {
  staggered_mcg_layout<Ty> layout;
  const Long sites = qinv.geo().local_volume() / 2;
  for (Int h = 0; h < Int(src.size()); ++h)
    mcg_signature_print_vector(call, "src", h, -1,
                               mcg_signature_public_vector(qinv.geo(), layout, src[h], sites));
}

template <class Ty>
inline void emit_staggered_mcg_signature_outputs(quda_inverter& qinv, const std::vector<Ty*>& res,
                                                 const std::vector<double>& offsets, Long call,
                                                 const char* stage) {
  staggered_mcg_layout<Ty> layout;
  const Long sites = qinv.geo().local_volume() / 2;
  const Int rhs = res.size() / offsets.size();
  for (size_t j = 0; j < offsets.size(); ++j)
    for (Int h = 0; h < rhs; ++h)
      mcg_signature_print_vector(call, stage, h, Long(j),
                                 mcg_signature_public_vector(qinv.geo(), layout, res[j * rhs + h], sites));
}

template <class Ty>
void log_staggered_multishift_mcg(quda_inverter& qinv, const std::vector<double>& offsets, Int rhs,
                                  double target, const MSCGOptions& options, const MSCGStats& stats,
                                  double total_seconds, bool corrected = false, Int correction_iter = 0,
                                  double correction_seconds = 0, double correction_gflops = 0,
                                  Long signature_call = -1, Int requested_batch = -1) {
  MCGDebugPhase phase("solve_logging_and_collectives");
  if (qinv.quda_verbos < 0 && qinv.quda_verbos != -2)
    return;
  double seconds = stats.solve_seconds;
  glb_max(seconds);
  glb_max(total_seconds);
  double cost = (stats.modeled_operator_flops + stats.modeled_vector_flops) / 1e9;
  sum_all_size(&cost, 1, 0);
  const char* io = std::is_same_v<Ty, ComplexD> ? "double" : "single";
  const char* accumulation = stats.accumulation_bytes == 4 ? "single" : "double";
  bool finite = true, met = true, zero = stats.iterations == 0;
  for (size_t i = 0; i < stats.iterative_residual.size(); ++i) {
    const double r = stats.iterative_residual[i];
    // The seed controls global termination; higher shifts either retire at
    // their own threshold or stop with the seed. Neither is an iteration cap.
    finite = finite && std::isfinite(r);
    if (i < size_t(rhs))
      met = met && r <= target;
    zero = zero && r == 0;
  }
  const bool checked = options.compute_true_res.value_or(false) || options.refine_tolerance > 0;
  if (checked) {
    Qassert(stats.true_residual.size() == offsets.size() * rhs);
    for (double r : stats.true_residual)
      finite = finite && std::isfinite(r);
  }
  const char* status = !finite                    ? "nonfinite-residual"
                       : zero                     ? "zero-source"
                       : options.fixed_iterations ? "fixed-iterations"
                       : stats.iterations == 0    ? "no-iterations"
                       : met                      ? "iterative-seed-target-reached"
                                                  : "iteration-limit";
  // Aggregate endpoint estimates over RHS without hiding NaN/Inf values.
  auto worst = [&](size_t j) {
    double value = 0;
    for (Int h = 0; h < rhs; ++h) {
      const double r = stats.iterative_residual[j * rhs + h];
      if (!std::isfinite(r))
        return r;
      value = std::max(value, r);
    }
    return value;
  };
  qmessage("  Donec : %ld iter + %d cor : m%d / %.6f secs; ", long(stats.seed_updates),
           stats.refine_iterations + correction_iter, rhs, total_seconds);
  if (seconds > 0 && std::isfinite(seconds))
    qmessage("%.3f GFLOPs/GPU", cost / (seconds * get_num_node()));
  else
    qmessage("n/a GFLOPs/GPU");
  qmessage(", Cost %.6f GFLOP; massl %.6g; res lh(%ld)%s %.1e / %.1e", cost, std::sqrt(offsets.front() / 4.0),
           long(offsets.size()), corrected ? " pre-cor" : "", worst(0), worst(offsets.size() - 1));
  if (!finite)
    qmessage("; nonfinite-residual");
  qmessage("\n");
  if (corrected)
    qmessage("  MCG rate/cost: modeled, MCG wall %.6f secs; headline wall includes QUDA-CG corrections\n",
             seconds);
  if (qinv.quda_verbos < 2)
    return;
  const Long call = signature_call >= 0 ? signature_call : next_staggered_mcg_diagnostic_id();
  const auto& d = stats.stopping;
  Qassert(d.seed_checks.size() == size_t(rhs) && d.shift_checks.size() == offsets.size() * rhs);
  const bool single_rhs = get_staggered_mcg_algorithm_config().single_rhs;
  const Int rhs_batch = single_rhs ? 1 : rhs;
  qmessage("  MCG_CONFIG call=%ld rhs=%d shifts=%ld target=%.17e limit=%d mode=%s batch=%d "
           "requested_batch=%d max_temporal_batch=%d effective_max_temporal_batch=%d wait_quda=%d single_rhs=%d "
           "rhs_batch=%d check_interval=%d shift_stop=%s beta=%s seed_operator=quda-folded io=%s krylov=%s "
           "update=%s accum=%s borrowed=%d scalar=double true=%s%s offsets=[",
           long(call), rhs, long(offsets.size()), target, d.maxiter,
           options.fixed_iterations ? "fixed" : "target", options.temporal_batch,
           requested_batch < 0 ? options.temporal_batch : requested_batch,
           get_staggered_mcg_algorithm_config().max_temporal_batch,
           staggered_mcg_effective_temporal_cap<Ty>(qinv.inv_param.cuda_prec_sloppy, options),
           int(staggered_mcg_wait_for_quda_enabled()), int(single_rhs), rhs_batch,
           rhs_batch,
           options.fixed_iterations ? "disabled"
           : options.shift_stop == MSCGShiftStopPolicy::precision_epsilon ? "precision-epsilon"
           : options.shift_stop == MSCGShiftStopPolicy::cpu_epsilon ? "double-epsilon"
                                                                    : "requested",
           options.beta == MSCGBetaPolicy::quda ? "quda" : "norm-ratio", io,
           staggered_mcg_precision_name(qinv.inv_param.cuda_prec_sloppy), accumulation,
           accumulation, int(stats.borrowed_solutions),
           checked ? "checked" : "not-evaluated",
           corrected ? " residuals=pre-QUDA-CG" : "");
  for (size_t j = 0; j < offsets.size(); ++j)
    qmessage("%s%.17e", j ? "," : "", offsets[j]);
  qmessage("]\n");
  const auto& reduction_config = get_mscg_reduction_config();
  if (reduction_config.precision || reduction_config.color_sum)
    qmessage("  MCG_REDUCTION call=%ld accum=%s color_sum=%d recurrence=double; cost=conventional FLOPs\n",
             long(call), reduction_config.precision ? "RealDD" : "double", int(reduction_config.color_sum));
  if (call == 1)
    qmessage("  MCG_CODES T=target E=epsilon Z=zero-source X=exact-zero I=initial-target F=fixed-limit "
             "L=iteration-limit; headline=sum seed iterations across RHS; polls=iteration:residual2/threshold2, pass<1; "
             "shift/RHS indices zero-based\n");
  for (Int h = 0; h < rhs; ++h) {
    qmessage("  MCG_MASS call=%ld rhs=%d entries=[", long(call), h);
    for (size_t j = 0; j < offsets.size(); ++j) {
      const size_t v = j * rhs + h;
      if (j)
        qmessage(",");
      qmessage("%ld:%.12e:%d:%c:%.12e:", long(j), std::sqrt(offsets[j] / 4.0), stats.shift_iterations[v],
               char(d.shift_reason[v]), stats.iterative_residual[v]);
      if (!checked)
        qmessage("n/a:n/a");
      else {
        const double true_res = stats.true_residual[v];
        if (std::isfinite(true_res))
          qmessage("%.12e:%d", true_res, int(true_res <= target));
        else
          qmessage("n/a:0");
      }
    }
    qmessage("]%s\n", corrected ? " residuals=pre-QUDA-CG" : "");
  }
  auto print_checks = [](const MSCGStoppingChecks& checks) {
    if (checks.previous_iteration < 0)
      qmessage("n/a");
    else
      qmessage("%d:%.17e", checks.previous_iteration, checks.previous_ratio);
    qmessage("->");
    if (checks.iteration < 0)
      qmessage("n/a");
    else
      qmessage("%d:%.17e", checks.iteration, checks.ratio);
  };
  for (Int h = 0; h < rhs; ++h) {
    qmessage("  MCG_SEED call=%ld rhs=%d iter=%d stop=%s b2=%.17e threshold2=%.17e res=%.12e polls=",
             long(call), h, stats.shift_iterations[h], mscg_stop_reason_name(d.seed_reason[h]),
             d.source_norm2[h], d.seed_threshold2[h], stats.iterative_residual[h]);
    print_checks(d.seed_checks[h]);
    if (checked)
      qmessage(" true%s=%.12e", stats.refine_iterations ? "(after-internal-CG)" : "", stats.true_residual[h]);
    else
      qmessage(" true=not-evaluated");
    qmessage("\n  MCG_ITERS call=%ld rhs=%d counts=[", long(call), h);
    for (size_t j = 0; j < offsets.size(); ++j) {
      const size_t v = j * rhs + h;
      qmessage("%s%d%c", j ? "," : "", stats.shift_iterations[v], char(d.shift_reason[v]));
    }
    qmessage("]");
    bool differing_updates = false, retired = false;
    for (size_t j = 1; j < offsets.size(); ++j) {
      const size_t v = j * rhs + h;
      if (d.effective_updates[v] != stats.shift_iterations[v]) {
        qmessage("%s%ld:%d(zero-coeff=%d)", differing_updates ? "," : " effective=[", long(j),
                 d.effective_updates[v], stats.shift_iterations[v] - d.effective_updates[v]);
        differing_updates = true;
      }
      retired = retired || d.shift_retired[v];
    }
    if (differing_updates)
      qmessage("]");
    qmessage("\n");
    if (retired) {
      qmessage("  MCG_TARGET call=%ld rhs=%d checks=[", long(call), h);
      bool first = true;
      for (size_t j = 1; j < offsets.size(); ++j) {
        const size_t v = j * rhs + h;
        if (!d.shift_retired[v])
          continue;
        qmessage("%s%ld(threshold2=%.17e)@", first ? "" : ",", long(j), d.shift_threshold2[v]);
        print_checks(d.shift_checks[v]);
        first = false;
      }
      qmessage("]\n");
    }
    if (checked) {
      double worst_true = 0;
      size_t worst_shift = 0;
      bool nonfinite = false, target_met = true;
      qmessage("  MCG_TRUE call=%ld rhs=%d res=[", long(call), h);
      for (size_t j = 0; j < offsets.size(); ++j) {
        const double r = stats.true_residual[j * rhs + h];
        qmessage("%s%.17e", j ? "," : "", r);
        target_met = target_met && std::isfinite(r) && r <= target;
        if (!std::isfinite(r)) {
          if (!nonfinite) {
            worst_true = r;
            worst_shift = j;
            nonfinite = true;
          }
        } else if (!nonfinite && r > worst_true) {
          worst_true = r;
          worst_shift = j;
        }
      }
      qmessage("] worst_shift=%ld worst=%.17e target_met=%d%s\n", long(worst_shift), worst_true,
               int(target_met), corrected ? " residuals=pre-QUDA-CG" : "");
    }
  }
  if (qinv.quda_verbos >= 3) {
    qmessage("Custom multishift: shifts=%ld, RHS=%d, min mass=%.12e, max mass=%.12e, target=%.6e, mode=%s, "
             "temporal batch=%d\n",
             long(offsets.size()), rhs, std::sqrt(offsets.front() / 4.0), std::sqrt(offsets.back() / 4.0),
             target, options.fixed_iterations ? "fixed-iterations" : "residual-target",
             options.temporal_batch);
    if (reduction_config.precision || reduction_config.color_sum)
      qmessage("  precision: source/solution=%s, MCG operator/Krylov=%s, solution updates=%s, reductions=%s, "
               "recurrence=double, true-residual operator=%s\n",
               io, staggered_mcg_precision_name(qinv.inv_param.cuda_prec_sloppy), accumulation,
               reduction_config.precision ? "RealDD" : "double", checked ? io : "not evaluated");
    else
      qmessage("  precision: source/solution=%s, MCG operator/Krylov=%s, solution updates=%s, "
               "reductions/recurrence=double, true-residual operator=%s\n",
               io, staggered_mcg_precision_name(qinv.inv_param.cuda_prec_sloppy), accumulation,
               checked ? io : "not evaluated");
    for (Int h = 0; h < rhs; ++h)
      for (Int end = 0; end < (offsets.size() == 1 ? 1 : 2); ++end) {
        const size_t j = end ? offsets.size() - 1 : 0, v = j * rhs + h;
        qmessage("  MCG%s rhs=%d, %s mass=%.12e, shift=%ld, %d iterations, relative residual: iterated=%.12e",
                 corrected ? " before correction" : "", h,
                 offsets.size() == 1 ? "lowest/highest"
                 : end               ? "highest"
                                     : "lowest",
                 std::sqrt(offsets[j] / 4.0), long(j), stats.shift_iterations[v],
                 stats.iterative_residual[v]);
        if (checked)
          qmessage(", true%s=%.12e", stats.refine_iterations ? " (after internal CG)" : "",
                   stats.true_residual[v]);
        else
          qmessage(", true=not evaluated");
        qmessage("\n");
      }
    // Retain the exhaustive legacy records at level 3+.
    for (size_t j = 0; j < offsets.size(); ++j)
      for (Int h = 0; h < rhs; ++h) {
        const size_t v = j * rhs + h;
        qmessage("  MCG_SHIFT_ITER,%ld,%d,%d,%.12e\n", long(j), h, stats.shift_iterations[v],
                 stats.iterative_residual[v]);
      }
  }
  qmessage("  MCG status=%s; %d internal-CG iter, %d QUDA-CG correction iter; MCG wall %.6f secs; modeled "
           "Cost %.6f GFLOP\n",
           status, stats.refine_iterations, correction_iter, seconds, cost);
  if (corrected)
    qmessage("  QUDA-CG correction: precision=%s, %d aggregate shift/RHS iterations, native %.6f secs, "
             "native cost %.6f GFLOP; total wall includes correction setup and precision reloads\n",
             io, correction_iter, correction_seconds, correction_gflops);
}

// Same six inputs as get_staggered_multishift_even. Optional controls do not
// override precision: qinv selects D/D, D/F, or F/F.
template <class Ty>
void get_staggered_multishift_even_mcg(quda_inverter& qinv, std::vector<Ty*> res, std::vector<Ty*> src,
                                       const std::vector<double>& offsets, double err, Int niter,
                                       const MSCGOptions& options = {}) {
  if (staggered_mcg_print_mem_enabled())
    print_mem_info("MCG begin");
  MCGDebugSolve debug(src.size(),offsets.size(),sizeof(Ty),int(qinv.inv_param.cuda_prec_sloppy)*2,
    0);
  mcg_debug_public_buffers(res,src,offsets.size(),qinv.geo().local_volume()*3/2);
  MSCGStats local_stats;
  MSCGOptions opt = resolve_staggered_mcg_options(options);
  opt.temporal_batch = resolve_staggered_mcg_temporal_batch<Ty>(
      opt.temporal_batch, offsets.size(), qinv.inv_param.cuda_prec_sloppy, opt);
  if (mcg_debug_active()) mcg_debug_state().window=opt.temporal_batch;
  if (!opt.stats)
    opt.stats = &local_stats;
  opt.stopping_diagnostics = opt.stopping_diagnostics || qinv.quda_verbos >= 2;
  Long signature_call = -1;
  if (mcg_signature_enabled()) {
    mcg_signature_config_once();
    signature_call = next_staggered_mcg_diagnostic_id();
    MCGDebugPhase phase("signature_inputs");
    emit_staggered_mcg_signature_gauge(qinv, signature_call);
    emit_staggered_mcg_signature_sources(qinv, src, signature_call);
  }
  qacc_barrier(dummy);
  const auto start = std::chrono::steady_clock::now();
  staggered_multishift_even_mcg_impl(qinv, res, src, offsets, err, niter, opt);
  qacc_barrier(dummy);
  opt.stats->solve_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  if (signature_call >= 0)
    emit_staggered_mcg_signature_outputs(qinv, res, offsets, signature_call, "mcg");
  log_staggered_multishift_mcg<Ty>(qinv, offsets, src.size(), err, opt, *opt.stats, opt.stats->solve_seconds,
                                   false, 0, 0, 0, signature_call, options.temporal_batch);
  if (opt.print_mem)
    print_mem_info("MCG end");
}

// Compatibility corrections can change QUDA's resident sloppy gauges. Reload on entry
// and on exit, and discard cached Dirac objects before their gauges disappear.
// The caller's complete configuration is restored; only solve counters escape.
// QUDA corrections and scoped inverter-state restoration
class StaggeredCorrectionPrecisionGuard {
  quda_inverter& q;
  QudaInvertParam inv;
  QudaGaugeParam gauge;
  Int prec;
  double mass;

  void reload() {
    clear_staggered_mcg_correction_cache();
    staggered_mcg_matrix_cache<ComplexD>().clear();
    staggered_mcg_matrix_cache<ComplexF>().clear();
    auto& links = q.qinv_twist_link == QTRUE ? q.quda_gf_default_twist : q.quda_gf_default;
    Qassert(links.size() != 0);
    q.setup_link_p(links.data(), 0, 0, q.qinv_twist_link);
    q.alloc_csfield_gpu();
    q.setup_mat_mass(mass, true);
  }

public:
  StaggeredCorrectionPrecisionGuard(quda_inverter& q_, QudaPrecision p)
      : q(q_), inv(q.inv_param), gauge(q.gauge_param), prec(q.prec_type_check), mass(q.mass_mat) {
    q.inv_param.cpu_prec = p;
    q.inv_param.cuda_prec = q.inv_param.cuda_prec_sloppy = p;
    q.inv_param.cuda_prec_refinement_sloppy = p;
    q.inv_param.cuda_prec_precondition = q.inv_param.cuda_prec_eigensolver = p;
    q.gauge_param.cuda_prec = q.gauge_param.cuda_prec_sloppy = p;
    q.gauge_param.cuda_prec_refinement_sloppy = p;
    q.gauge_param.cuda_prec_precondition = q.gauge_param.cuda_prec_eigensolver = p;
    reload();
  }

  ~StaggeredCorrectionPrecisionGuard() {
    q.inv_param = inv;
    q.gauge_param = gauge;
    reload();
    q.inv_param = inv;
    q.gauge_param = gauge;
    q.prec_type_check = prec;
  }
};

// Level 2 supplies bounded wrapper diagnostics, without native CG progress.
// Scope only the native call so setup/reload and caller verbosity are retained.
class StaggeredCorrectionVerbosityGuard {
  quda_inverter& q;
  QudaVerbosity global, inverter;
  bool enabled;

public:
  explicit StaggeredCorrectionVerbosityGuard(quda_inverter& q_)
      : q(q_), global(::getVerbosity()), inverter(q.inv_param.verbosity), enabled(q.quda_verbos == 2) {
    if (enabled) {
      ::setVerbosity(QUDA_SUMMARIZE);
      q.inv_param.verbosity = QUDA_SUMMARIZE;
    }
  }

  StaggeredCorrectionVerbosityGuard(const StaggeredCorrectionVerbosityGuard&) = delete;
  StaggeredCorrectionVerbosityGuard& operator=(const StaggeredCorrectionVerbosityGuard&) = delete;

  ~StaggeredCorrectionVerbosityGuard() {
    if (enabled) {
      q.inv_param.verbosity = inverter;
      ::setVerbosity(global);
    }
  }
};

template <class Ty>
void staggered_multishift_cg_correction_legacy(quda_inverter& qinv, const std::vector<Ty*>& res,
                                        const std::vector<Ty*>& src, const std::vector<double>& offsets,
                                        double err) {
  const Int limit = qlat::get_env_long_default("qlat_multi_refine_iter", 1000);
  if (limit < 0)
    throw std::invalid_argument("qlat_multi_refine_iter must be nonnegative");
  double seconds = 0, gflops = 0;
  Int iterations = 0;
  const Long L = qinv.geo().local_volume() * 3 / 2;
  // QUDA CG rejects zero RHS. The MCG result for these RHS stays zero.
  MSCGReduction reduction;
  std::vector<Int> active;
  for (Int h = 0; h < Int(src.size()); ++h) {
    const double norm = *reduction.dot(src[h], src[h], L, 1);
    if (!std::isfinite(norm))
      throw std::invalid_argument("nonfinite MCG source");
    if (norm > 0)
      active.push_back(h);
  }
  if (!active.empty() && limit > 0) {
    constexpr QudaPrecision precision =
        std::is_same_v<Ty, ComplexD> ? QUDA_DOUBLE_PRECISION : QUDA_SINGLE_PRECISION;
    StaggeredCorrectionPrecisionGuard guard(qinv, precision);
    auto& p = qinv.inv_param;
    p.inv_type = QUDA_CG_INVERTER;
    p.solution_type = QUDA_MATPC_SOLUTION;
    p.solve_type = QUDA_DIRECT_PC_SOLVE;
    p.matpc_type = QUDA_MATPC_EVEN_EVEN;
    p.input_location = p.output_location = QUDA_CUDA_FIELD_LOCATION;
    p.num_offset = 0;
    p.compute_true_res = 1;
    p.use_init_guess = QUDA_USE_INIT_GUESS_YES;
    p.tol = err;
    p.maxiter = limit;
    p.residual_type = QUDA_L2_RELATIVE_RESIDUAL;
    p.reliable_delta = 0.1;
    p.use_sloppy_partial_accumulator = 0;
    qlat::vector<void*> sp(1), rp(1);
    for (size_t j = 0; j < offsets.size(); ++j) {
      p.mass = std::sqrt(offsets[j] / 4.0);
      for (Int h : active) {
        sp[0] = src[h];
        rp[0] = res[j * src.size() + h];
        {
          StaggeredCorrectionVerbosityGuard verbosity(qinv);
          qinv.callMultiSrcQuda(rp, sp);
          staggered_mcg_wait_for_quda();
        }
        seconds += p.secs;
        gflops += p.gflops;
        iterations += p.iter;
        if (qinv.quda_verbos >= 2)
          qmessage("MCG_QUDA_CORRECTION,shift=%ld,rhs=%d,rhs_batch=1,precision=%d/%d/%d/"
                   "%d,iter=%d,secs=%.12e,true=%.12e,target=%.12e\n",
                   long(j), h, int(p.cpu_prec), int(p.cuda_prec), int(p.cuda_prec_sloppy),
                   int(p.cuda_prec_refinement_sloppy), p.iter, p.secs, p.true_res[0], err);
      }
    }
  }
  qinv.inv_param.secs = seconds;
  qinv.inv_param.gflops = gflops;
  qinv.inv_param.iter = iterations;
}

// QUDA's public invert interface constructs matrices, fields and a new solver
// per call. Keep the same QUDA CG algorithm, but retain its setup for ordinary
// staggered corrections. All matrices use the resident precise gauge; no global
// precision change or resident-link replacement is necessary.
template <class Ty> class StaggeredMCGCorrection {
  struct Source {
    MCGDebugBuffer<Ty> data;
    quda::ColorSpinorField field;
  };
  Geometry geo;
  quda::ColorSpinorParam field_param;
  QudaDslashType dslash_type;
  QudaDagType dagger;
  std::vector<std::unique_ptr<Source>> sources;
  MCGDebugBuffer<Ty> solution;
  quda::ColorSpinorField x;
  // Destruction order matters: CG holds references to parameters and matrix,
  // and DiracM holds a reference to Dirac.
  std::unique_ptr<quda::Dirac> dirac;
  std::unique_ptr<quda::DiracM> matrix;
  std::unique_ptr<quda::SolverParam> parameters;
  std::unique_ptr<quda::CG> cg;
  quda::TimeProfile profile{"MCG correction", false};
  std::vector<Int> active;
  std::vector<double> scale, source_norm;

  static quda::ColorSpinorParam native_param(const quda_inverter& q) {
    if constexpr (std::is_same_v<Ty, ComplexD>)
      return q.cs_gpuD;
    else
      return q.cs_gpuF;
  }

public:
  MSCGReduction reduction;
  Long allocations = 0, solver_creations = 0, dirac_creations = 0;

  explicit StaggeredMCGCorrection(const quda_inverter& q)
      : geo(q.geo()), field_param(native_param(q)), dslash_type(q.inv_param.dslash_type),
        dagger(q.inv_param.dagger) {
    field_param.create = QUDA_REFERENCE_FIELD_CREATE;
  }

  bool compatible(const quda_inverter& q) const {
    const auto p = native_param(q);
    if (Compare_geo_less(geo, q.geo()) || Compare_geo_less(q.geo(), geo) ||
        dslash_type != q.inv_param.dslash_type || dagger != q.inv_param.dagger ||
        p.Precision() != field_param.Precision() || p.fieldOrder != field_param.fieldOrder ||
        p.nDim != field_param.nDim || p.pad != field_param.pad ||
        p.siteSubset != field_param.siteSubset || p.siteOrder != field_param.siteOrder ||
        p.nSpin != field_param.nSpin || p.nColor != field_param.nColor ||
        p.nVec != field_param.nVec || p.nVec_actual != field_param.nVec_actual)
      return false;
    for (Int d = 0; d < p.nDim; ++d)
      if (p.x[d] != field_param.x[d])
        return false;
    return true;
  }

  void apply(quda_inverter& q, const std::vector<Ty*>& res, const std::vector<Ty*>& src,
             const std::vector<double>& offsets, double err, Int limit) {
    const Long L = q.geo().local_volume() * 3 / 2;
    const Long before = allocations, reduction_before = reduction.allocations;
    active.clear();
    active.reserve(src.size());
    scale.resize(src.size());
    source_norm.resize(src.size());
    // Preserve rejection of nonfinite sources and avoid QUDA's zero-source
    // rejection. This reduction's storage survives subsequent corrections.
    for (Int h = 0; h < Int(src.size()); ++h) {
      const double norm = *reduction.dot(src[h], src[h], L, 1);
      if (!std::isfinite(norm))
        throw std::invalid_argument("nonfinite MCG source");
      if (norm > 0)
        active.push_back(h);
    }
    allocations += reduction.allocations - reduction_before;
    if (active.empty()) {
      q.inv_param.secs = q.inv_param.gflops = 0;
      q.inv_param.iter = 0;
      return;
    }
    QudaInvertParam p = q.inv_param;
    constexpr QudaPrecision precision =
        std::is_same_v<Ty, ComplexD> ? QUDA_DOUBLE_PRECISION : QUDA_SINGLE_PRECISION;
    p.cpu_prec = p.cuda_prec = p.cuda_prec_sloppy = precision;
    p.cuda_prec_refinement_sloppy = p.cuda_prec_precondition = p.cuda_prec_eigensolver = precision;
    p.inv_type = QUDA_CG_INVERTER;
    p.solution_type = QUDA_MATPC_SOLUTION;
    p.solve_type = QUDA_DIRECT_PC_SOLVE;
    p.matpc_type = QUDA_MATPC_EVEN_EVEN;
    p.input_location = p.output_location = QUDA_CUDA_FIELD_LOCATION;
    p.num_src = 1;
    p.num_offset = 0;
    p.compute_true_res = 1;
    p.use_init_guess = QUDA_USE_INIT_GUESS_YES;
    p.tol = err;
    p.maxiter = limit;
    p.residual_type = QUDA_L2_RELATIVE_RESIDUAL;
    p.reliable_delta = 0.1;
    p.use_sloppy_partial_accumulator = 0;
    p.iter = 0;
    p.secs = p.gflops = 0;
    // Always rebind, even if a replacement gauge reuses the same address.
    // Mass is the only changing scalar in this even-even staggered operator.
    quda::DiracParam dp;
    quda::setDiracParam(dp, &p, true);
    if (!dirac) {
      dirac.reset(quda::Dirac::create(dp));
      matrix = std::make_unique<quda::DiracM>(*dirac);
      ++dirac_creations;
    } else
      dirac->updateFields(dp.gauge, dp.fatGauge, dp.longGauge, dp.clover);
    if (!parameters)
      parameters = std::make_unique<quda::SolverParam>(p);
    else
      *parameters = quda::SolverParam(p);
    if (!cg) {
      cg = std::make_unique<quda::CG>(*matrix, *matrix, *matrix, *matrix, *parameters);
      ++solver_creations;
    }
    auto reference = [&](Ty* pointer) {
      auto fp = field_param;
      fp.v = pointer;
      quda::ColorSpinorField field(fp);
      if (field.Bytes() != size_t(L) * sizeof(Ty) || field.Stride() != size_t(L / 3))
        throw std::runtime_error("MCG correction: padded native field cannot reference a packed vector");
      return field;
    };
    if (solution.size() != L) {
      solution.resize(L);
      solution.debug_name("correction.solution");
      ++allocations;
      x = reference(solution.data());
    }
    if (sources.size() < src.size())
      sources.resize(src.size());
    // Normalize the source once per RHS, with the same QUDA BLAS operation as
    // its public interface. Copies of already normalized data can then serve
    // every mass without additional source conversions or allocations.
    for (Int h : active) {
      if (!sources[h]) {
        sources[h] = std::make_unique<Source>();
        sources[h]->data.resize(L);
        sources[h]->data.debug_name("correction.source."+std::to_string(h));
        sources[h]->field = reference(sources[h]->data.data());
        ++allocations;
      }
      auto& b = sources[h]->field;
      cpy_GPU(sources[h]->data.data(), src[h], L, 1, 1,
              staggered_mcg_wait_for_quda_enabled() ? QTRUE : QFALSE);
      scale[h] = 1;
      if (p.solver_normalization == QUDA_SOURCE_NORMALIZATION) {
        source_norm[h] = std::sqrt(double(quda::blas::norm2(b)));
        scale[h] = 1.0 / source_norm[h];
        quda::blas::ax(scale[h], b);
        staggered_mcg_wait_for_quda();
      }
    }
    struct Verbosity {
      explicit Verbosity(QudaVerbosity value) { ::pushVerbosity(value); }
      ~Verbosity() { ::popVerbosity(); }
    } verbosity(q.quda_verbos == 2 ? QUDA_SUMMARIZE : p.verbosity);
    double seconds = 0, gflops = 0;
    Int iterations = 0;
    for (size_t j = 0; j < offsets.size(); ++j) {
      p.mass = std::sqrt(offsets[j] / 4.0);
      dirac->setMass(p.mass);
      for (Int h : active) {
        p.iter = parameters->iter = 0;
        p.secs = p.gflops = 0;
        {
          quda::pushProfile timer(profile, &p);
          cpy_GPU(solution.data(), res[j * src.size() + h], L, 1, 1,
                  staggered_mcg_wait_for_quda_enabled() ? QTRUE : QFALSE);
          if (p.solver_normalization == QUDA_SOURCE_NORMALIZATION)
            quda::blas::ax(scale[h], x);
          quda::cvector_ref<quda::ColorSpinorField> out(x);
          quda::cvector_ref<const quda::ColorSpinorField> in(sources[h]->field);
          { MCGDebugPhase phase("quda_correction_cg");(*cg)(out, in); }
          staggered_mcg_wait_for_quda();
          parameters->updateInvertParam(p);
          if (p.solver_normalization == QUDA_SOURCE_NORMALIZATION) {
            quda::blas::ax(source_norm[h], x);
            staggered_mcg_wait_for_quda();
          }
          cpy_GPU(res[j * src.size() + h], solution.data(), L, 1, 1,
                  staggered_mcg_wait_for_quda_enabled() ? QTRUE : QFALSE);
        }
        seconds += p.secs;
        gflops += p.gflops;
        iterations += p.iter;
        if (q.quda_verbos >= 2)
          qmessage("MCG_QUDA_CORRECTION,shift=%ld,rhs=%d,rhs_batch=1,precision=%d/%d/%d/"
                   "%d,iter=%d,secs=%.12e,true=%.12e,target=%.12e,cached=1\n",
                   long(j), h, int(p.cpu_prec), int(p.cuda_prec), int(p.cuda_prec_sloppy),
                   int(p.cuda_prec_refinement_sloppy), p.iter, p.secs, p.true_res[0], err);
      }
    }
    q.inv_param.secs = seconds;
    q.inv_param.gflops = gflops;
    q.inv_param.iter = iterations;
    if (q.quda_verbos >= 2)
      qmessage("MCG_QUDA_CORRECTION_CACHE,allocations=%ld,growth=%ld,solvers=%ld,dirac=%ld,"
               "source_capacity=%ld,gauge_reloads=0\n",
               long(allocations), long(allocations - before), long(solver_creations),
               long(dirac_creations), long(sources.size()));
  }
};

template <class Ty>
inline std::map<quda_inverter*, std::unique_ptr<StaggeredMCGCorrection<Ty>>>&
staggered_mcg_correction_cache() {
  static std::map<quda_inverter*, std::unique_ptr<StaggeredMCGCorrection<Ty>>> cache;
  return cache;
}

inline void clear_staggered_mcg_correction_cache(quda_inverter* qinv) {
  auto clear = [qinv](auto& cache) {
    if (qinv)
      cache.erase(qinv);
    else
      cache.clear();
  };
  clear(staggered_mcg_correction_cache<ComplexD>());
  clear(staggered_mcg_correction_cache<ComplexF>());
}

// The general interface still handles split-grid redistribution and optional
// resident/deflation/chronological features. Never silently bypass them.
inline bool staggered_mcg_correction_reusable(const QudaInvertParam& p) {
  return p.split_grid[0] == 1 && p.split_grid[1] == 1 && p.split_grid[2] == 1 && p.split_grid[3] == 1 &&
         !p.eig_param && !p.preconditioner && !p.chrono_use_resident && !p.chrono_make_resident &&
         !p.use_resident_solution && !p.make_resident_solution &&
         p.dirac_order == QUDA_INTERNAL_DIRAC_ORDER &&
         p.mass_normalization == QUDA_MASS_NORMALIZATION &&
         (p.solver_normalization == QUDA_SOURCE_NORMALIZATION ||
          p.solver_normalization == QUDA_DEFAULT_NORMALIZATION) &&
         p.distance_pc_alpha0 == 0;
}

template <class Ty> inline bool staggered_mcg_correction_native(const quda_inverter& q) {
  constexpr auto precision = std::is_same_v<Ty, ComplexD> ? QUDA_DOUBLE_PRECISION : QUDA_SINGLE_PRECISION;
  const auto& p = std::is_same_v<Ty, ComplexD> ? q.cs_gpuD : q.cs_gpuF;
  const size_t bytes = size_t(q.geo().local_volume() * 3 / 2) * sizeof(Ty);
  return q.inv_param.cuda_prec == precision && p.Precision() == precision &&
         p.location == QUDA_CUDA_FIELD_LOCATION && p.fieldOrder == QUDA_NATIVE_FIELD_ORDER &&
         p.siteSubset == QUDA_PARITY_SITE_SUBSET && p.siteOrder == QUDA_EVEN_ODD_SITE_ORDER &&
         p.nDim == 4 && p.nSpin == 1 && p.nColor == 3 && p.nVec == 1 && p.nVec_actual == 1 &&
         p.pad == 0 && !p.is_composite && !p.is_component && ALIGNMENT_ADJUST(bytes) == bytes;
}

template <class Ty>
void staggered_multishift_cg_correction(quda_inverter& qinv, const std::vector<Ty*>& res,
                                        const std::vector<Ty*>& src, const std::vector<double>& offsets,
                                        double err) {
  MCGDebugPhase phase("correction");
  const Long requested = qlat::get_env_long_default("qlat_multi_refine_iter", 1000);
  if (requested < 0 || requested > std::numeric_limits<Int>::max())
    throw std::invalid_argument("qlat_multi_refine_iter must be nonnegative and fit Int");
  if (requested == 0) {
    qinv.inv_param.secs = qinv.inv_param.gflops = 0;
    qinv.inv_param.iter = 0;
    return;
  }
  if (!staggered_mcg_correction_reusable(qinv.inv_param) || !staggered_mcg_correction_native<Ty>(qinv)) {
    staggered_multishift_cg_correction_legacy(qinv, res, src, offsets, err);
    return;
  }
  install_staggered_multishift_mcg_cache_cleanup(qinv);
  auto& context = staggered_mcg_correction_cache<Ty>()[&qinv];
  if (!context || !context->compatible(qinv))
    context = std::make_unique<StaggeredMCGCorrection<Ty>>(qinv);
  context->apply(qinv, res, src, offsets, err, Int(requested));
}

// Read once across input precisions, as with the custom reduction controls.
// Unset or zero leaves the caller's iteration limit unchanged.
// Experimental backend configuration and execution
inline Int get_staggered_mcg_maxiter() {
  static const Int maxiter = [] {
    const Long value = qlat::get_env_long_default(std::string("qlat_mcg_maxiter"), 0);
    if (value < 0 || value > std::numeric_limits<Int>::max())
      throw std::invalid_argument(
          "qlat_mcg_maxiter must be 0 (caller limit) or a positive value fitting Int");
    return Int(value);
  }();
  return maxiter;
}

// Fixed-function diagnostic backends use this process-static damping value.
// It intentionally has no effect on the custom or legacy QUDA backends.
inline double get_staggered_mcg_type_rho() {
  static const double rho = [] {
    const double value = qlat::get_env_double_default("qlat_quda_mcg_type_rho", 0.1);
    if (!std::isfinite(value) || value < 0)
      throw std::invalid_argument("qlat_quda_mcg_type_rho must be finite and nonnegative");
    return value;
  }();
  return rho;
}

// A test-only tolerance override shared by DD/D/F and read on first type-5/6
// use. Unset preserves the caller's tolerance; zero is not a disable sentinel.
struct StaggeredType5ErrorConfig {
  bool overridden = false;
  double value = std::numeric_limits<double>::quiet_NaN();
};

inline const StaggeredType5ErrorConfig& get_staggered_type5_error_config() {
  static const StaggeredType5ErrorConfig config = [] {
    StaggeredType5ErrorConfig result;
    const char* setting = std::getenv("qlat_mcg_type5_err");
    if (!setting)
      return result;
    char* end = nullptr;
    const double value = std::strtod(setting, &end);
    const bool parsed = end != setting;
    while (*end && std::isspace(static_cast<unsigned char>(*end)))
      ++end;
    if (!parsed || *end || !std::isfinite(value) || value <= 0)
      throw std::invalid_argument(
          "qlat_mcg_type5_err must be a finite positive number, or unset to use the input err");
    result.overridden = true;
    result.value = value;
    return result;
  }();
  return config;
}

inline double get_staggered_type5_error(double requested) {
  const auto& config = get_staggered_type5_error_config();
  const double effective = config.overridden ? config.value : requested;
  if (!std::isfinite(effective) || effective <= 0)
    throw std::invalid_argument("Standalone type 5/6 requires a finite positive effective tolerance");
  return effective;
}

struct StaggeredMCGTypeStats {
  Int type = 0, steps = 0;
  double rho = 0, seconds = 0, gflops = 0;
  QudaPrecision io_precision = QUDA_DOUBLE_PRECISION;
  QudaPrecision krylov_precision = QUDA_DOUBLE_PRECISION;
  const char* internal_precision = nullptr;
  double requested_err = std::numeric_limits<double>::quiet_NaN();
  double effective_err = std::numeric_limits<double>::quiet_NaN();
  bool err_overridden = false;
  bool shift_epsilon = false;
  Int rhs_batch = 1, check_interval = 1;
  Long mass_rhs_iterations = 0, matrix_calls = 0;
  double retirement_target = std::numeric_limits<double>::quiet_NaN();
  std::vector<Int> last_poll;
  std::vector<char> retired;
  Long total_iterations = 0;
  std::vector<double> iterative_residual;
  std::vector<Int> iterations;
  std::vector<char> stop;
  // Post-function residuals (type 5 uses DD on exported double fields), never recurrence estimates.
  std::vector<double> residual;
};

inline Int get_staggered_mcg_experiment_iterations(Int niter) {
  const Int cap = get_staggered_mcg_maxiter();
  return cap > 0 ? std::min(niter, cap) : niter;
}

template <class Ty>
inline void validate_staggered_mcg_type_input(quda_inverter& qinv, const std::vector<Ty*>& res,
                                              const std::vector<Ty*>& src, const std::vector<double>& offsets,
                                              double err, Int niter, Int backend = -1) {
  static_assert(std::is_same_v<Ty, ComplexD> || std::is_same_v<Ty, ComplexF>);
  if (backend < 0)
    backend = qinv.mcg_type;
  const double target = backend >= 5 ? get_staggered_type5_error(err) : err;
  if (src.empty() || offsets.empty() || src.size() > size_t(std::numeric_limits<Int>::max()) ||
      offsets.size() > std::numeric_limits<size_t>::max() / src.size() ||
      res.size() != src.size() * offsets.size() || niter < 0 || !std::isfinite(target) || target <= 0)
    throw std::invalid_argument("MCG replacement: invalid vector counts, target or iteration limit");
  for (double offset : offsets)
    if (!std::isfinite(offset) || offset < 0)
      throw std::invalid_argument("MCG replacement: offsets must be finite and nonnegative");
  if (backend == 6 && (offsets.size() > size_t(std::numeric_limits<Int>::max()) ||
                       !std::is_sorted(offsets.begin(), offsets.end())))
    throw std::invalid_argument("Standalone type 6 requires nondecreasing offsets fitting Int");
  for (auto p : src)
    if (!p)
      throw std::invalid_argument("MCG replacement: null source");
  for (auto p : res)
    if (!p)
      throw std::invalid_argument("MCG replacement: null output");
  Qassert(qinv.inv_param.dagger == QUDA_DAG_NO);
  if (backend >= 5) {
    if constexpr (!std::is_same_v<Ty, ComplexD>)
      throw std::invalid_argument("Standalone type 5/6 requires double public input/output");
    if (qinv.inv_param.cpu_prec != QUDA_DOUBLE_PRECISION ||
        qinv.inv_param.cuda_prec != QUDA_DOUBLE_PRECISION ||
        qinv.inv_param.dirac_order != QUDA_INTERNAL_DIRAC_ORDER)
      throw std::invalid_argument("Standalone type 5/6 requires native double public fields");
    check_staggered_mcg_param(qinv.cs_gpuD, QUDA_DOUBLE_PRECISION, qinv.geo().local_volume() / 2);
  } else {
    check_staggered_mcg_ordering_for_solve<Ty>(qinv, src.size());
  }
}

// Optional observation after all output vectors have been produced. Use the
// public/main precision and resident precise gauge even for a D/F recurrence.
// Neither this operator work nor its reductions enter time or cost accounting.
template <class Ty>
inline void measure_staggered_mcg_type_residuals(quda_inverter& qinv, const std::vector<Ty*>& res,
                                                 const std::vector<Ty*>& src,
                                                 const std::vector<double>& offsets,
                                                 StaggeredMCGTypeStats& stats) {
  staggered_mcg_matrix<Ty> matrix(&qinv, 1, true, std::is_same_v<Ty, ComplexD>);
  const Long L = matrix.get_vsize();
  Ty* x = matrix.buf[0].data();
  Ty* ap = matrix.buf[1].data();
  MSCGReduction reduction;
  for (size_t h = 0; h < src.size(); ++h) {
    Ty* b = src[h];
    const double b2 = reduction.dot(b, b, L, 1)[0];
    for (size_t j = 0; j < offsets.size(); ++j) {
      const size_t v = j * src.size() + h;
      Ty* out = res[v];
      qGPU_for(i, L, 1, { x[i] = out[i]; });
      matrix.multi_seed(ap, x, offsets[j]);
      qGPU_for(i, L, 1, { ap[i] = b[i] - ap[i]; });
      const double r2 = reduction.dot(ap, ap, L, 1)[0];
      // Preserve NaN/Inf rather than clamping them into a passing residual.
      stats.residual[v] = std::sqrt(b2 > 0 ? r2 / b2 : r2);
    }
  }
}

template <class Ty>
inline StaggeredMCGTypeStats
run_staggered_mcg_type2(quda_inverter& qinv, const std::vector<Ty*>& res, const std::vector<Ty*>& src,
                        const std::vector<double>& offsets, double err, Int niter) {
  quda_inverter::ResidentGaugeMutation resident_mutation;
  Qassert(!src.empty() && !offsets.empty() && res.size() == src.size() * offsets.size());
  StaggeredMCGTypeStats stats;
  stats.type = 2;
  stats.steps = get_staggered_mcg_experiment_iterations(niter);
  stats.rho = std::numeric_limits<double>::quiet_NaN();
  const Int rhs = src.size();
  const size_t total = offsets.size() * size_t(rhs);
  stats.iterations.assign(total, 0);
  stats.stop.assign(total, 'Q');
  stats.residual.assign(total, std::numeric_limits<double>::quiet_NaN());
  QudaInvertParam param = qinv.inv_param;
  param.inv_type = QUDA_CG_INVERTER;
  param.solution_type = QUDA_MATPC_SOLUTION;
  param.solve_type = QUDA_DIRECT_PC_SOLVE;
  param.matpc_type = QUDA_MATPC_EVEN_EVEN;
  param.input_location = QUDA_CUDA_FIELD_LOCATION;
  param.output_location = QUDA_CUDA_FIELD_LOCATION;
  param.use_init_guess = QUDA_USE_INIT_GUESS_NO;
  param.num_offset = 0;
  param.num_src = 1;
  param.preserve_source = QUDA_PRESERVE_SOURCE_YES;
  param.inv_type_precondition = QUDA_INVALID_INVERTER;
  param.deflate = QUDA_BOOLEAN_FALSE;
  param.make_resident_solution = param.use_resident_solution = 0;
  param.chrono_make_resident = param.chrono_use_resident = 0;
  param.dagger = QUDA_DAG_NO;
  param.tol = err;
  param.maxiter = stats.steps;
  param.residual_type = QUDA_L2_RELATIVE_RESIDUAL;
  param.reliable_delta = 0;
  param.use_alternative_reliable = 0;
  param.use_sloppy_partial_accumulator = 0;
  param.solution_accumulator_pipeline = 1;
  param.pipeline = 0;
  // Diagnostic residuals are evaluated separately, outside QUDA's timing.
  param.compute_true_res = 0;
  stats.io_precision = param.cuda_prec;
  stats.krylov_precision = param.cuda_prec_sloppy;
  const Long L = qinv.geo().local_volume() * 3 / 2;
  MSCGReduction reduction;
  std::vector<bool> zero(rhs);
  for (Int h = 0; h < rhs; ++h)
    zero[h] = reduction.dot(src[h], src[h], L, 1)[0] == 0;
  double seconds = 0, gflops = 0;
  Long iterations = 0;
  for (size_t j = 0; j < offsets.size(); ++j) {
    Qassert(std::isfinite(offsets[j]) && offsets[j] >= 0);
    param.mass = std::sqrt(offsets[j] / 4.0);
    for (Int h = 0; h < rhs; ++h) {
      const size_t v = j * size_t(rhs) + h;
      if (zero[h] || stats.steps == 0) {
        Ty* out = res[v];
        qGPU_for(i, L, 1, { out[i] = Ty(0); });
        stats.stop[v] = zero[h] ? 'Z' : 'L';
        continue;
      }
      param.secs = param.gflops = 0;
      param.iter = 0;
      qacc_barrier(dummy);
      const auto start = std::chrono::steady_clock::now();
      invertQuda(res[j * size_t(rhs) + h], src[h], &param);
      qacc_barrier(dummy);
      seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
      // The public single-CG API exposes the iteration count, but neither the
      // iterative residual nor an unambiguous stopping reason. Q records that
      // limitation instead of inferring convergence from a later true residual.
      stats.iterations[v] = param.iter;
      iterations += param.iter;
      gflops += param.gflops;
    }
  }
  stats.total_iterations = iterations;
  stats.seconds = seconds;
  stats.gflops = gflops;
  return stats;
}

// Read once: public fields remain double; this selects the complete local
// operator/CG arithmetic, independently of QUDA's sloppy-precision setting.
enum class StaggeredType5Precision { DD, D, F };

inline StaggeredType5Precision get_staggered_type5_precision() {
  static const auto precision = [] {
    const char* setting = std::getenv("qlat_mcg_type5_precision");
    std::string value = setting ? setting : "D";
    for (char& c : value)
      c = char(std::toupper(static_cast<unsigned char>(c)));
    if (value == "DD")
      return StaggeredType5Precision::DD;
    if (value == "D")
      return StaggeredType5Precision::D;
    if (value == "F")
      return StaggeredType5Precision::F;
    throw std::invalid_argument("qlat_mcg_type5_precision must be DD, D, or F");
  }();
  return precision;
}

// Each batch shares a lowest-mass Krylov stream across all shifts. Public
// inputs/outputs remain double; all recurrence state stays in the chosen Real.
template <class Real>
inline StaggeredMCGTypeStats
run_staggered_mcg_type6_impl(quda_inverter& qinv, const std::vector<ComplexD*>& res,
                             const std::vector<ComplexD*>& src, const std::vector<double>& offsets,
                             double err, Int niter) {
  using C = ComplexT<Real>;
  validate_staggered_mcg_type_input(qinv, res, src, offsets, err, niter, 6);
  const double effective = get_staggered_type5_error(err);
  const auto& config = get_staggered_mcg_algorithm_config();
  const bool overridden = get_staggered_type5_error_config().overridden;
  const bool epsilon = !overridden &&
                       (config.shift_stop == MSCGShiftStopPolicy::cpu_epsilon ||
                        config.shift_stop == MSCGShiftStopPolicy::precision_epsilon);
  const Int rhs = Int(src.size()), J = Int(offsets.size()), batch = config.single_rhs ? 1 : rhs;
  StaggeredSimpleOperator<Real> matrix(qinv);
  const Long L = matrix.get_vsize(), N = L * rhs;
  std::vector<Real> absolute(J), delta(J);
  for (Int j = 0; j < J; ++j) {
    matrix.set_shift(offsets[j]);
    absolute[j] = matrix.shift;
  }
  for (Int j = 0; j < J; ++j)
    delta[j] = cg_real_sub(absolute[j], absolute[0]);
  matrix.set_shift(offsets[0]);
  matrix.prepare_batch(batch);
  vector_gpu<C> b(N);
  vector_gpu<ComplexD> mask(N);
  staggered_mcg_layout<ComplexD> layout;
  for (Int h = 0; h < rhs; ++h) {
    matrix.import_public(b.data() + Long(h) * L, src[h], layout);
    const ComplexD* input = src[h];
    ComplexD* row = mask.data() + Long(h) * L;
    qacc_for(i, L, { row[i] = ComplexD(input[i].real() != 0.0 || input[i].imag() != 0.0 ? 1.0 : 0.0, 0.0); });
  }
  CGTypedBatchReduction<double> public_reduction;
  public_reduction.prepare(L, rhs);
  const auto nonzero = public_reduction.dot(mask.data(), mask.data(), L, rhs);
  TypedMultiShiftCG<Real> solver;
  solver.prepare(L, batch, J);
  StaggeredMCGTypeStats stats;
  stats.type = 6;
  stats.steps = get_staggered_mcg_experiment_iterations(niter);
  stats.rhs_batch = batch;
  stats.check_interval = batch;
  stats.rho = std::numeric_limits<double>::quiet_NaN();
  stats.requested_err = err;
  stats.effective_err = effective;
  stats.err_overridden = overridden;
  stats.shift_epsilon = epsilon;
  stats.retirement_target =
      epsilon ? (config.shift_stop == MSCGShiftStopPolicy::precision_epsilon
                     ? double(std::numeric_limits<Real>::epsilon())
                     : std::numeric_limits<double>::epsilon())
              : effective;
  stats.io_precision = QUDA_DOUBLE_PRECISION;
  stats.krylov_precision = std::is_same_v<Real, float> ? QUDA_SINGLE_PRECISION : QUDA_DOUBLE_PRECISION;
  stats.internal_precision = std::is_same_v<Real, RealDD> ? "DD" : std::is_same_v<Real, double> ? "D" : "F";
  const size_t S = size_t(J) * rhs;
  stats.iterations.resize(S);
  stats.stop.resize(S);
  stats.iterative_residual.resize(S);
  stats.last_poll.resize(S);
  stats.retired.resize(S);
  stats.residual.assign(S, std::numeric_limits<double>::quiet_NaN());
  for (Int first = 0; first < rhs; first += batch) {
    std::vector<char> original(batch);
    for (Int h = 0; h < batch; ++h)
      original[h] = nonzero[first + h] > 0;
    qacc_barrier(dummy);
    const auto start = std::chrono::steady_clock::now();
    const auto result =
        solver.solve(matrix, b.data() + Long(first) * L, delta, effective, stats.retirement_target,
                     stats.steps, batch, config.beta == MSCGBetaPolicy::quda, original, epsilon);
    qacc_barrier(dummy);
    stats.seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    stats.total_iterations += result.iterations;
    stats.matrix_calls += result.matrix_calls;
    stats.gflops += (result.vector_flops +
                     double(result.matrix_calls) * batch * (570.0 * qinv.geo().local_volume() + 4.0 * L)) *
                    1e-9;
    for (Int j = 0; j < J; ++j)
      for (Int h = 0; h < batch; ++h) {
        const size_t from = size_t(j) * batch + h, to = size_t(j) * rhs + first + h;
        stats.iterations[to] = result.shift_iterations[from];
        stats.stop[to] = result.shift_reason[from];
        stats.iterative_residual[to] = result.iterative_residual[from];
        stats.last_poll[to] = result.last_poll[from];
        stats.retired[to] = result.retired[from];
        stats.mass_rhs_iterations += result.shift_iterations[from];
        matrix.export_public(res[to], solver.solution_data() + from * L, layout);
      }
  }
  return stats;
}

template <class Ty>
inline StaggeredMCGTypeStats
run_staggered_mcg_type6(quda_inverter& qinv, const std::vector<Ty*>& res, const std::vector<Ty*>& src,
                        const std::vector<double>& offsets, double err, Int niter) {
  if constexpr (!std::is_same_v<Ty, ComplexD>)
    throw std::invalid_argument("Standalone type 6 requires double public input/output");
  else {
    switch (get_staggered_type5_precision()) {
    case StaggeredType5Precision::DD:
      return run_staggered_mcg_type6_impl<RealDD>(qinv, res, src, offsets, err, niter);
    case StaggeredType5Precision::D:
      return run_staggered_mcg_type6_impl<double>(qinv, res, src, offsets, err, niter);
    case StaggeredType5Precision::F:
      return run_staggered_mcg_type6_impl<float>(qinv, res, src, offsets, err, niter);
    }
    throw std::logic_error("invalid standalone precision");
  }
}

template <class Real>
inline StaggeredMCGTypeStats
run_staggered_mcg_type5_impl(quda_inverter& qinv, const std::vector<ComplexD*>& res,
                             const std::vector<ComplexD*>& src, const std::vector<double>& offsets,
                             double err, Int niter) {
  using C = ComplexT<Real>;
  const double effective_err = get_staggered_type5_error(err);
  StaggeredSimpleOperator<Real> matrix(qinv);
  const Long L = matrix.get_vsize();
  vector_gpu<C> b(L), x(L), r(L), p(L), ap(L);
  vector_gpu<ComplexD> public_nonzero_mask(L);
  CGTypedReduction<Real> reduction;
  reduction.prepare(L);
  CGTypedReduction<double> public_reduction;
  public_reduction.prepare(L);
  cg_iter<StaggeredSimpleOperator<Real>> solver(&matrix, &matrix);
  staggered_mcg_layout<ComplexD> layout;
  StaggeredMCGTypeStats stats;
  stats.type = 5;
  stats.requested_err = err;
  stats.effective_err = effective_err;
  stats.err_overridden = get_staggered_type5_error_config().overridden;
  stats.steps = get_staggered_mcg_experiment_iterations(niter);
  stats.rho = std::numeric_limits<double>::quiet_NaN();
  stats.io_precision = QUDA_DOUBLE_PRECISION;
  stats.krylov_precision = std::is_same_v<Real, float> ? QUDA_SINGLE_PRECISION : QUDA_DOUBLE_PRECISION;
  stats.internal_precision = std::is_same_v<Real, RealDD> ? "DD" : std::is_same_v<Real, double> ? "D" : "F";
  const size_t total = offsets.size() * src.size();
  stats.iterations.assign(total, 0);
  stats.stop.assign(total, 'L');
  stats.iterative_residual.assign(total, std::numeric_limits<double>::quiet_NaN());
  stats.residual.assign(total, std::numeric_limits<double>::quiet_NaN());
  for (size_t h = 0; h < src.size(); ++h) {
    // A boolean mask avoids underflow when testing a tiny public source before
    // float conversion. Norm-zero after conversion must not become stop=Z.
    const ComplexD* input = src[h];
    ComplexD* mask = public_nonzero_mask.data();
    qacc_for(i, L,
             { mask[i] = ComplexD(input[i].real() != 0.0 || input[i].imag() != 0.0 ? 1.0 : 0.0, 0.0); });
    const bool original_nonzero = public_reduction.dot(mask, mask, L) > 0.0;
    matrix.import_public(b.data(), src[h], layout);
    for (size_t j = 0; j < offsets.size(); ++j) {
      matrix.set_shift(offsets[j]);
      qacc_barrier(dummy);
      const auto start = std::chrono::steady_clock::now();
      const auto one =
          solver.template solve_typed<Real>(x.data(), b.data(), r.data(), p.data(), ap.data(), L, stats.steps,
                                            effective_err, original_nonzero, &reduction);
      qacc_barrier(dummy);
      stats.seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
      const size_t v = j * src.size() + h;
      stats.iterations[v] = one.iterations;
      stats.total_iterations += one.iterations;
      stats.stop[v] = one.reason == TypedCGStop::zero_source           ? 'Z'
                      : one.reason == TypedCGStop::numerical_breakdown ? 'B'
                      : one.reason == TypedCGStop::converged ? (cg_real_zero(one.residual_norm2) ? 'X' : 'T')
                                                             : 'L';
      stats.iterative_residual[v] = one.relative_residual();
      // Conventional real arithmetic, independent of DD instruction expansion.
      stats.gflops +=
          (double(one.matrix_calls) * (570.0 * qinv.geo().local_volume() + 4.0 * L) + one.vector_flops) *
          1e-9;
      matrix.export_public(res[v], x.data(), layout);
    }
  }
  return stats;
}

template <class Ty>
inline StaggeredMCGTypeStats
run_staggered_mcg_type5(quda_inverter& qinv, const std::vector<Ty*>& res, const std::vector<Ty*>& src,
                        const std::vector<double>& offsets, double err, Int niter) {
  if constexpr (!std::is_same_v<Ty, ComplexD>) {
    throw std::invalid_argument("Standalone type 5 requires double public input/output");
  } else {
    switch (get_staggered_type5_precision()) {
    case StaggeredType5Precision::DD:
      return run_staggered_mcg_type5_impl<RealDD>(qinv, res, src, offsets, err, niter);
    case StaggeredType5Precision::D:
      return run_staggered_mcg_type5_impl<double>(qinv, res, src, offsets, err, niter);
    case StaggeredType5Precision::F:
      return run_staggered_mcg_type5_impl<float>(qinv, res, src, offsets, err, niter);
    }
    throw std::logic_error("invalid standalone precision");
  }
}

// Observe the EXPORTED double solution using the full saved double links with
// DD arithmetic. This independently measures both recurrence and export error.
inline void measure_staggered_mcg_type5_residuals(quda_inverter& qinv, const std::vector<ComplexD*>& res,
                                                  const std::vector<ComplexD*>& src,
                                                  const std::vector<double>& offsets,
                                                  StaggeredMCGTypeStats& stats) {
  using C = ComplexT<RealDD>;
  StaggeredSimpleOperator<RealDD> matrix(qinv);
  const Long L = matrix.get_vsize();
  vector_gpu<C> b(L), x(L), ap(L);
  CGTypedReduction<RealDD> reduction;
  reduction.prepare(L);
  staggered_mcg_layout<ComplexD> layout;
  for (size_t h = 0; h < src.size(); ++h) {
    matrix.import_public(b.data(), src[h], layout);
    const RealDD b2 = reduction.dot(b.data(), b.data(), L);
    for (size_t j = 0; j < offsets.size(); ++j) {
      const size_t v = j * src.size() + h;
      matrix.set_shift(offsets[j]);
      matrix.import_public(x.data(), res[v], layout);
      matrix.multi(ap.data(), x.data());
      C* residual = ap.data();
      const C* source = b.data();
      qacc_for(i, L, { residual[i] = cg_complex_sub(source[i], residual[i]); });
      const RealDD r2 = reduction.dot(ap.data(), ap.data(), L);
      stats.residual[v] = std::sqrt(double(cg_real_positive(b2) ? cg_real_div(r2, b2) : r2));
    }
  }
}

template <class T, class Ty>
inline StaggeredMCGTypeStats
run_staggered_mcg_type_fixed_impl(quda_inverter& qinv, const std::vector<Ty*>& res,
                                  const std::vector<Ty*>& src, const std::vector<double>& offsets, Int niter,
                                  Int type, bool native_double) {
  Qassert(!src.empty() && !offsets.empty() && res.size() == src.size() * offsets.size());
  const Int rhs = src.size(), steps = get_staggered_mcg_experiment_iterations(niter);
  const double rho = get_staggered_mcg_type_rho();
  staggered_mcg_matrix<T> matrix(&qinv, rhs, false, native_double);
  vector<T*> pointers;
  matrix.set_buf_pointers(pointers, rhs);
  T* x = pointers[0];
  T* ap = pointers[1];
  const Long N = matrix.get_vsize(), L = N / rhs, sites = L / 3;
  vector_gpu<T> source(N);
  vector_gpu<T> momentum;
  if (type == 4)
    momentum.resize(N);
  T* b = source.data();
  T* pvec = type == 4 ? momentum.data() : nullptr;
  staggered_mcg_layout<T> layout = matrix.vector_layout();
  for (Int h = 0; h < rhs; ++h) {
    Ty* in = src[h];
    const Long base = Long(h) * L;
    qGPU_for(i, sites * 3, 1, {
      const Long site = i / 3;
      const Int color = Int(i % 3);
      const Long it = base + layout.storage_index(site, color, sites);
      const Long id = layout.template public_index<Ty>(site, color, sites);
      b[it] = T(in[id]);
    });
  }
  StaggeredMCGTypeStats stats;
  stats.type = type;
  stats.steps = steps;
  stats.rho = rho;
  stats.io_precision = qinv.inv_param.cuda_prec;
  stats.krylov_precision = qinv.inv_param.cuda_prec_sloppy;
  const size_t total = offsets.size() * size_t(rhs);
  stats.iterations.assign(total, steps);
  stats.stop.assign(total, 'F');
  stats.residual.assign(total, std::numeric_limits<double>::quiet_NaN());
  using Real = typename staggered_mcg_layout<T>::Real;
  const Real trho = Real(rho), tbeta = Real(0.5);
  for (size_t j = 0; j < offsets.size(); ++j) {
    Qassert(std::isfinite(offsets[j]) && offsets[j] >= 0);
    qacc_barrier(dummy);
    const auto start = std::chrono::steady_clock::now();
    // Every mass defines an independent function of b, with fresh state.
    qGPU_for(i, N, 1, {
      x[i] = T(0);
      if (type == 4)
        pvec[i] = T(0);
    });
    for (Int k = 0; k < steps; ++k) {
      matrix.multi_seed(ap, x, offsets[j]);
      qGPU_for(i, N, 1, {
        const T r = b[i] - ap[i];
        if (type == 3) {
          x[i] += trho * r;
        } else {
          pvec[i] = r + tbeta * pvec[i];
          x[i] += trho * pvec[i];
        }
      });
    }
    for (Int h = 0; h < rhs; ++h) {
      Ty* out = res[j * size_t(rhs) + h];
      const Long base = Long(h) * L;
      qGPU_for(i, sites * 3, 1, {
        const Long site = i / 3;
        const Int color = Int(i % 3);
        const Long it = base + layout.storage_index(site, color, sites);
        const Long id = layout.template public_index<Ty>(site, color, sites);
        out[id] = Ty(x[it]);
      });
    }
    qacc_barrier(dummy);
    stats.seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  }
  stats.total_iterations = Long(steps) * offsets.size() * rhs;
  // Conventional operator model (same convention as type 0), plus the massive
  // xpay and real vector updates. No setup, diagnostics or reductions counted.
  const double site_flops = qinv.inv_param.dslash_type == QUDA_ASQTAD_DSLASH ? 1146.0 : 570.0;
  stats.gflops = double(steps) * offsets.size() *
                 (site_flops * qinv.geo().local_volume() * rhs + (type == 3 ? 10.0 : 14.0) * N) / 1e9;
  return stats;
}

template <class Ty>
inline StaggeredMCGTypeStats
run_staggered_mcg_type_fixed(quda_inverter& qinv, const std::vector<Ty*>& res, const std::vector<Ty*>& src,
                             const std::vector<double>& offsets, Int niter, Int type) {
  Qassert(type == 3 || type == 4);
  if (get_staggered_mcg_algorithm_config().single_rhs && src.size() > 1) {
    StaggeredMCGTypeStats total;
    total.type = type;
    total.steps = get_staggered_mcg_experiment_iterations(niter);
    total.rho = get_staggered_mcg_type_rho();
    total.iterations.resize(res.size());
    total.stop.resize(res.size());
    total.residual.resize(res.size());
    std::vector<Ty*> one_res(offsets.size());
    for (size_t h = 0; h < src.size(); ++h) {
      for (size_t j = 0; j < offsets.size(); ++j)
        one_res[j] = res[j * src.size() + h];
      const auto one =
          run_staggered_mcg_type_fixed(qinv, one_res, std::vector<Ty*>{src[h]}, offsets, niter, type);
      total.seconds += one.seconds;
      total.gflops += one.gflops;
      total.total_iterations += one.total_iterations;
      for (size_t j = 0; j < offsets.size(); ++j) {
        const size_t v = j * src.size() + h;
        total.iterations[v] = one.iterations[j];
        total.stop[v] = one.stop[j];
        total.residual[v] = one.residual[j];
      }
    }
    return total;
  }
  if constexpr (std::is_same_v<Ty, ComplexD>) {
    Qassert(qinv.inv_param.cpu_prec == QUDA_DOUBLE_PRECISION &&
            qinv.inv_param.cuda_prec == QUDA_DOUBLE_PRECISION);
    Qassert(qinv.inv_param.cuda_prec_sloppy == QUDA_DOUBLE_PRECISION ||
            qinv.inv_param.cuda_prec_sloppy == QUDA_SINGLE_PRECISION);
    if (qinv.inv_param.cuda_prec_sloppy == QUDA_DOUBLE_PRECISION)
      return run_staggered_mcg_type_fixed_impl<ComplexD, Ty>(qinv, res, src, offsets, niter, type, true);
    return run_staggered_mcg_type_fixed_impl<ComplexF, Ty>(qinv, res, src, offsets, niter, type, true);
  } else {
    Qassert(qinv.inv_param.cpu_prec == QUDA_SINGLE_PRECISION &&
            qinv.inv_param.cuda_prec == QUDA_SINGLE_PRECISION &&
            qinv.inv_param.cuda_prec_sloppy == QUDA_SINGLE_PRECISION);
    return run_staggered_mcg_type_fixed_impl<ComplexF, Ty>(qinv, res, src, offsets, niter, type, false);
  }
}

template <class Ty>
inline void log_staggered_mcg_type(quda_inverter& qinv, const StaggeredMCGTypeStats& stats,
                                   const std::vector<double>& offsets, Int rhs, double target, Long call) {
  if (stats.type < 2 || stats.type > 6 || (qinv.quda_verbos < 0 && qinv.quda_verbos != -2))
    return;
  const char* method = stats.type == 2   ? "independent-cg"
                       : stats.type == 3 ? "richardson"
                       : stats.type == 4 ? "heavy-ball"
                       : stats.type == 5 ? "standalone-cg"
                                         : "standalone-multishift-cg";
  double seconds = stats.seconds;
  glb_max(seconds);
  double cost = stats.gflops;
  if (stats.type == 3 || stats.type == 4)
    sum_all_size(&cost, 1, 0);
  if (stats.type == 6)
    qmessage("  Donec : %ld batch iter + 0 cor / %.6f secs; Cost %.6f GFLOP; type=6 method=%s "
             "mass_rhs_iterations=%ld matrix_calls=%ld\n",
             long(stats.total_iterations), seconds, cost, method, long(stats.mass_rhs_iterations),
             long(stats.matrix_calls));
  else
    qmessage("  Donec : %ld total mass/RHS iter + 0 cor / %.6f secs; Cost %.6f GFLOP; type=%d method=%s\n",
             long(stats.total_iterations), seconds, cost, stats.type, method);
  if (qinv.quda_verbos < 2)
    return;
  if (stats.type == 6) {
    const auto& config = get_staggered_mcg_algorithm_config();
    qmessage("  MCG_TYPE_CONFIG call=%ld type=6 method=standalone-multishift-cg operator=simple-staggered-ee "
             "steps=%d rhs=%d shifts=%ld rhs_batch=%d check_gap=%d io=double krylov=%s gauge=%s reduction=%s "
             "coefficients=%s requested_err=%.17e effective_err=%.17e err_source=%s shift_stop=%s "
             "retirement_target=%.17e beta=%s batch_iterations=%ld mass_rhs_iterations=%ld matrix_calls=%ld "
             "refine=bypassed reliable_updates=0 true=post-export-standalone-DD T=target E=epsilon "
             "I=initial-target X=exact-zero Z=zero-source L=iteration-limit B=numerical-breakdown\n",
             long(call), stats.steps, rhs, long(offsets.size()), stats.rhs_batch, stats.check_interval,
             stats.internal_precision, stats.internal_precision, stats.internal_precision,
             stats.internal_precision, stats.requested_err, stats.effective_err,
             stats.err_overridden ? "env" : "input",
             stats.err_overridden  ? "override"
             : stats.shift_epsilon ? "cpu-epsilon"
                                   : "requested",
             stats.retirement_target, config.beta == MSCGBetaPolicy::quda ? "quda" : "norm-ratio",
             long(stats.total_iterations), long(stats.mass_rhs_iterations), long(stats.matrix_calls));
    for (Int h = 0; h < rhs; ++h)
      for (size_t j = 0; j < offsets.size(); ++j) {
        const size_t v = j * size_t(rhs) + h;
        const double r = stats.residual[v];
        const bool met = stats.stop[v] != 'B' && std::isfinite(r) && r <= stats.effective_err;
        qmessage(
            "  MCG_TYPE_RESULT call=%ld type=6 rhs=%d shift=%ld mass=%.12e iter=%d stop=%c iter_res=%.12e "
            "true_res=%.12e target_met=%d retirement_target=%.17e last_poll=%d retired=%d\n",
            long(call), h, long(j), std::sqrt(offsets[j] / 4.0), stats.iterations[v], stats.stop[v],
            stats.iterative_residual[v], r, int(met), j ? stats.retirement_target : stats.effective_err,
            stats.last_poll[v], int(stats.retired[v]));
      }
    return;
  }
  if (stats.type == 5) {
    qmessage("  MCG_TYPE_CONFIG call=%ld type=5 method=standalone-cg operator=simple-staggered-ee steps=%d "
             "rhs=%d shifts=%ld rhs_batch=1 io=double krylov=%s gauge=%s reduction=%s coefficients=%s "
             "requested_err=%.17e effective_err=%.17e err_source=%s check_gap=1 refine=bypassed "
             "reliable_updates=0 true=post-export-standalone-DD T=target X=exact-zero Z=zero-source "
             "L=iteration-limit B=numerical-breakdown\n",
             long(call), stats.steps, rhs, long(offsets.size()), stats.internal_precision,
             stats.internal_precision, stats.internal_precision, stats.internal_precision,
             stats.requested_err, stats.effective_err, stats.err_overridden ? "env" : "input");
    for (Int h = 0; h < rhs; ++h)
      for (size_t j = 0; j < offsets.size(); ++j) {
        const size_t v = j * size_t(rhs) + h;
        const double r = stats.residual[v];
        const bool met = stats.stop[v] != 'B' && std::isfinite(r) && r <= stats.effective_err;
        qmessage("  MCG_TYPE_RESULT call=%ld type=5 rhs=%d shift=%ld mass=%.12e iter=%d stop=%c "
                 "iter_res=%.12e true_res=%.12e target_met=%d\n",
                 long(call), h, long(j), std::sqrt(offsets[j] / 4.0), stats.iterations[v], stats.stop[v],
                 stats.iterative_residual[v], r, int(met));
      }
    return;
  }
  qmessage("  MCG_TYPE_CONFIG call=%ld type=%d method=%s rho=%.17e steps=%d rhs=%d shifts=%ld rhs_batch=%d "
           "io=%s krylov=%s refine=bypassed reliable_updates=0 true=post-function-public-precision "
           "Q=quda-stop-unavailable F=fixed-limit Z=zero-source L=iteration-limit\n",
           long(call), stats.type, method, stats.rho, stats.steps, rhs, long(offsets.size()),
           stats.type == 2 || get_staggered_mcg_algorithm_config().single_rhs ? 1 : rhs,
           staggered_mcg_precision_name(stats.io_precision),
           staggered_mcg_precision_name(stats.krylov_precision));
  for (Int h = 0; h < rhs; ++h)
    for (size_t j = 0; j < offsets.size(); ++j) {
      const size_t v = j * size_t(rhs) + h;
      const double r = stats.residual[v];
      const bool met = std::isfinite(r) && r <= target;
      qmessage("  MCG_TYPE_RESULT call=%ld type=%d rhs=%d shift=%ld mass=%.12e iter=%d stop=%c iter_res=n/a "
               "true_res=%.12e target_met=%d\n",
               long(call), stats.type, h, long(j), std::sqrt(offsets[j] / 4.0), stats.iterations[v],
               stats.stop[v], r, int(met));
    }
}

// Public dispatcher. Offset convention is sigma=4*m*m and results are
// shift-major then RHS-major, exactly as in the original QUDA entry point.
// Public backend dispatch
template <class Ty>
void get_staggered_multishift_even(quda_inverter& qinv, std::vector<Ty*> res, std::vector<Ty*> src,
                                   const std::vector<double>& offsets, double err, Int niter) {
  if (qinv.mcg_type < 0 || qinv.mcg_type > 6)
    throw std::invalid_argument(
        "qlat_quda_mcg_type must be 0 (custom), 1 (QUDA), 2 (independent CG), 3 (Richardson), 4 "
        "(heavy-ball), 5 (standalone CG), or 6 (standalone multishift CG)");
  if (qinv.mcg_type == 1) {
    // Direct legacy QUDA dispatch: keep pointer order and arguments unchanged.
    get_staggered_multishift_even_quda(qinv, res, src, offsets, err, niter);
    return;
  }
  if (qinv.mcg_type >= 2) {
    if (qinv.mcg_type >= 5)
      check_staggered_simple_single_gpu_once();
    validate_staggered_mcg_type_input(qinv, res, src, offsets, err, niter);
    const Int experiment_niter = get_staggered_mcg_experiment_iterations(niter);
    Long signature_call = -1;
    if (mcg_signature_enabled()) {
      mcg_signature_config_once();
      signature_call = next_staggered_mcg_diagnostic_id();
      emit_staggered_mcg_signature_gauge(qinv, signature_call);
      emit_staggered_mcg_signature_sources(qinv, src, signature_call);
    }
    StaggeredMCGTypeStats stats =
        qinv.mcg_type == 6   ? run_staggered_mcg_type6(qinv, res, src, offsets, err, experiment_niter)
        : qinv.mcg_type == 5 ? run_staggered_mcg_type5(qinv, res, src, offsets, err, experiment_niter)
        : qinv.mcg_type == 2
            ? run_staggered_mcg_type2(qinv, res, src, offsets, err, experiment_niter)
            : run_staggered_mcg_type_fixed(qinv, res, src, offsets, experiment_niter, qinv.mcg_type);
    qacc_barrier(dummy);
    qinv.inv_param.secs = stats.seconds;
    qinv.inv_param.gflops = stats.gflops;
    // The public counter is Int; the diagnostic retains the full work count.
    qinv.inv_param.iter = Int(std::min<Long>(stats.total_iterations, std::numeric_limits<Int>::max()));
    if (qinv.quda_verbos >= 2) {
      if (qinv.mcg_type >= 5) {
        if constexpr (std::is_same_v<Ty, ComplexD>)
          measure_staggered_mcg_type5_residuals(qinv, res, src, offsets, stats);
      } else
        measure_staggered_mcg_type_residuals(qinv, res, src, offsets, stats);
    }
    if (signature_call >= 0)
      emit_staggered_mcg_signature_outputs(qinv, res, offsets, signature_call,
                                           qinv.mcg_type == 2   ? "type2"
                                           : qinv.mcg_type == 3 ? "type3"
                                           : qinv.mcg_type == 4 ? "type4"
                                           : qinv.mcg_type == 5 ? "type5"
                                                                : "type6");
    const Long call = signature_call >= 0 ? signature_call : next_staggered_mcg_diagnostic_id();
    log_staggered_mcg_type<Ty>(qinv, stats, offsets, src.size(), err, call);
    return;
  }
  if (staggered_mcg_print_mem_enabled())
    print_mem_info("MCG begin");
  MCGDebugSolve debug(src.size(),offsets.size(),sizeof(Ty),int(qinv.inv_param.cuda_prec_sloppy)*2,
    0);
  mcg_debug_public_buffers(res,src,offsets.size(),qinv.geo().local_volume()*3/2);
  // Apply an optional environment cap while preserving early stopping.
  const Int experiment_maxiter = get_staggered_mcg_maxiter();
  if (experiment_maxiter > 0)
    niter = std::min(niter, experiment_maxiter);
  const Int refine = qlat::get_env_long_default("qlat_multi_refine", 0);
  if (refine < 0 || refine > 2)
    throw std::invalid_argument("qlat_multi_refine must be 0, 1 or 2");
  MSCGOptions options = resolve_staggered_mcg_options({});
  options.temporal_batch = resolve_staggered_mcg_temporal_batch<Ty>(
      options.temporal_batch, offsets.size(), qinv.inv_param.cuda_prec_sloppy, options);
  if (mcg_debug_active()) mcg_debug_state().window=options.temporal_batch;
  MSCGStats stats;
  options.stats = &stats;
  options.stopping_diagnostics = qinv.quda_verbos >= 2;
  // Corrections perform their own required residual checks. Avoid a redundant
  // custom pre-correction observation even when check-only work is enabled.
  if (refine != 0)
    options.compute_true_res = false;
  options.refine_tolerance = 0;
  std::vector<MCGDebugBuffer<Ty>> saved_sources;
  if (refine != 0)
    preserve_staggered_mcg_aliased_sources(res, src, qinv.geo().local_volume()*3/2, saved_sources);
  Long signature_call = -1;
  if (mcg_signature_enabled()) {
    mcg_signature_config_once();
    signature_call = next_staggered_mcg_diagnostic_id();
    MCGDebugPhase phase("signature_inputs");
    emit_staggered_mcg_signature_gauge(qinv, signature_call);
    emit_staggered_mcg_signature_sources(qinv, src, signature_call);
  }
  qacc_barrier(dummy);
  const auto start = std::chrono::steady_clock::now();
  staggered_multishift_even_mcg_impl(qinv, res, src, offsets, err, niter, options);
  qacc_barrier(dummy);
  stats.solve_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  if (signature_call >= 0)
    emit_staggered_mcg_signature_outputs(qinv, res, offsets, signature_call, "mcg");
  const auto correction_start = std::chrono::steady_clock::now();
  qinv.inv_param.secs = 0;
  qinv.inv_param.gflops = 0;
  qinv.inv_param.iter = 0;
  if (refine != 0)
    staggered_multishift_cg_correction(qinv, res, src, offsets, err);
  const Int correction_iter = qinv.inv_param.iter;
  const double correction_seconds = qinv.inv_param.secs, correction_gflops = qinv.inv_param.gflops;
  qacc_barrier(dummy);
  const double correction_wall =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - correction_start).count();
  qinv.inv_param.iter += stats.iterations;
  if (signature_call >= 0) {
    if (refine != 0)
      emit_staggered_mcg_signature_outputs(qinv, res, offsets, signature_call, "final");
  }
  // Include all correction setup, checks and solve overhead in custom wall time.
  // gflops counts only the QUDA corrections, not the custom MCG work.
  qinv.inv_param.secs = stats.solve_seconds + correction_wall;
  log_staggered_multishift_mcg<Ty>(qinv, offsets, src.size(), err, options, stats, qinv.inv_param.secs,
                                   refine != 0, correction_iter, correction_seconds, correction_gflops,
                                   signature_call, 0);
  if (options.print_mem)
    print_mem_info("MCG end");
}
} // namespace qlat
#endif
