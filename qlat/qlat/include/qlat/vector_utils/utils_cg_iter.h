// utils_cg_iter.h
// Gen Wang
// Jun. 2023

#ifndef UTILS_CG_ITER_H
#define UTILS_CG_ITER_H

#pragma once
#include "utils_float_type.h"
#include "utils_cg_precision.h"
#include <limits>
#include <stdexcept>

namespace qlat
{

/*
  M is the method which provides Ty bufs, and multi on Ty; have nothing to do
  with vector_cs Tv is the src and result vectors, could have different
  precision as Ty
*/

template <class M, typename Tv>
inline double check_residue(Tv& res, Tv& src, M& mat, const Int ir = 0,
                            const Int is = 0)
{
  const long Ndata = mat.get_vsize();
  Qassert(long(src.nsum) == Ndata);
  Int GPU = 1;
  auto* ab = mat.mat_src[0];
  auto* ax = mat.mat_src[1];
  auto* ak = mat.mat_src[2];
  //
  src.copy_to(ab, is, true);
  res.copy_to(ax, ir, true);
  mat.multi(ak, ax);
  qGPU_for(isp, Ndata, GPU, { ax[isp] = ak[isp] - ab[isp]; });
  auto residu = vec_norm2(ax, ax, Ndata);
  qmessage("QLAT::mat residue %.8e \n", residu.real());
  return residu.real();
}

// Opt-in precision-preserving CG status.  These types are independent of the
// historic iter()/call() entry points below, whose arithmetic remains intact.
enum class TypedCGStop {
  converged,
  zero_source,
  iteration_limit,
  numerical_breakdown
};

inline const char* typed_cg_stop_name(TypedCGStop reason)
{
  switch (reason) {
  case TypedCGStop::converged: return "tolerance";
  case TypedCGStop::zero_source: return "zero_source";
  case TypedCGStop::iteration_limit: return "iteration_limit";
  case TypedCGStop::numerical_breakdown: return "numerical_breakdown";
  }
  return "numerical_breakdown";
}

template<class Real> struct TypedCGResult {
  Int iterations = 0;
  Int matrix_calls = 0;
  double vector_flops = 0;
  TypedCGStop reason = TypedCGStop::iteration_limit;
  Real source_norm2 = Real(0.0);
  Real residual_norm2 = Real(0.0);

  bool converged() const
  {
    return reason == TypedCGStop::converged || reason == TypedCGStop::zero_source;
  }
  double relative_residual() const
  {
    if (reason == TypedCGStop::zero_source) return 0.0;
    if (!cg_real_finite(source_norm2) || !cg_real_positive(source_norm2) ||
        !cg_real_finite(residual_norm2) || residual_norm2 < Real(0.0))
      return std::numeric_limits<double>::infinity();
    return std::sqrt(double(cg_real_div(residual_norm2, source_norm2)));
  }
};

template <class M>
struct cg_iter {
  // double resid;
  // Int niter;
  Int debug_mat;
  Int multi_count;
  double resid;
  Int niter;
  M* matP;
  M* matS;
  double it_res;
  Int verbose;
  // qlat::vector<qlat::Complex > RkL;

  cg_iter(M* matP_, M* matS_, const Int debug_mat_ = 0)
  {
    // qmessage("setup chebyshev \n");
    matP = matP_;
    matS = matS_;
    // buf.resize(4);
    multi_count = 0;
    //
    /// RkL.resize(2);
    it_res = 0.0;
    resid = 0.0;
    niter = 0;
    verbose = qlat::get_env_long_default(std::string("qlat_cg_verbose"), 1);
    ;
    debug_mat = debug_mat_;
  }

  // Classic real-coefficient CG from a zero guess on raw complex pointers.
  // All vectors, products, tree reductions, alpha/beta coefficients, and the
  // squared-residual comparison remain in Real.  Only diagnostics returned by
  // relative_residual() narrow to double.  Caller-owned buffers and an optional
  // preallocated reduction keep setup outside the caller's solve timing.
  //
  // original_source_nonzero must describe the public source before conversion:
  // a nonzero double source which becomes zero (or has a zero norm) in float is
  // numerical breakdown, not successful zero-source convergence.
  template<class Real>
  inline TypedCGResult<Real> solve_typed(
      ComplexT<Real>* x, const ComplexT<Real>* source,
      ComplexT<Real>* r, ComplexT<Real>* p, ComplexT<Real>* ap,
      Long length, Int max_iterations, double relative_tolerance,
      bool original_source_nonzero = false,
      CGTypedReduction<Real>* prepared_reduction = nullptr)
  {
    if (length < 0 || max_iterations < 0 ||
        !std::isfinite(relative_tolerance) || relative_tolerance < 0.0)
      throw std::invalid_argument("typed CG requires nonnegative length, cap, and finite tolerance");
    Qassert(matP != nullptr);
    Qassert(x != source && r != source && p != source && ap != source);
    Qassert(x != r && x != p && x != ap && r != p && r != ap && p != ap);
    TypedCGResult<Real> stats;
    multi_count = 0;
    niter = max_iterations;
    resid = relative_tolerance;
    it_res = 0.0;
    CGTypedReduction<Real> local_reduction;
    CGTypedReduction<Real>& reduction = prepared_reduction ?
      *prepared_reduction : local_reduction;
    if (!prepared_reduction) reduction.prepare(length);
    qacc_for(i, length, {
      x[i] = ComplexT<Real>(Real(0.0), Real(0.0));
      r[i] = source[i];
      p[i] = source[i];
    });
    Real rho = reduction.dot(r, r, length);
    stats.vector_flops += 4.0 * length;
    stats.source_norm2 = rho;
    stats.residual_norm2 = rho;
    const Real tolerance = Real(relative_tolerance);
    const Real threshold = cg_real_mul(cg_real_mul(tolerance, tolerance), rho);
    stats.vector_flops += 2.0;
    if (!cg_real_finite(rho) || rho < Real(0.0) ||
        !cg_real_finite(tolerance) || !cg_real_finite(threshold)) {
      stats.reason = TypedCGStop::numerical_breakdown;
    } else if (cg_real_zero(rho)) {
      stats.reason = original_source_nonzero ? TypedCGStop::numerical_breakdown :
                                              TypedCGStop::zero_source;
    } else if (cg_real_leq(rho, threshold)) {
      stats.reason = TypedCGStop::converged;
    } else {
      for (Int iteration = 0; iteration < max_iterations; ++iteration) {
        matP->multi(ap, p);
        ++stats.matrix_calls;
        const Real curvature = reduction.dot(p, ap, length);
        stats.vector_flops += 4.0 * length;
        if (!cg_real_finite(curvature) || !cg_real_positive(curvature)) {
          stats.reason = TypedCGStop::numerical_breakdown;
          break;
        }
        const Real alpha = cg_real_div(rho, curvature);
        stats.vector_flops += 1.0;
        if (!cg_real_finite(alpha)) {
          stats.reason = TypedCGStop::numerical_breakdown;
          break;
        }
        qacc_for(i, length, {
          const ComplexT<Real> next_x = cg_complex_add(x[i], cg_complex_scale(alpha, p[i]));
          x[i] = next_x;
          // An overflowing solution must not appear converged merely because
          // the recurrence residual cancels. Propagate it into the norm pass.
          r[i] = cg_complex_finite(next_x) ?
            cg_complex_sub(r[i], cg_complex_scale(alpha, ap[i])) :
            ComplexT<Real>(Real(NAN), Real(0.0));
        });
        ++stats.iterations;
        stats.vector_flops += 8.0 * length;
        const Real next_rho = reduction.dot(r, r, length);
        stats.vector_flops += 4.0 * length;
        stats.residual_norm2 = next_rho;
        if (!cg_real_finite(next_rho) || next_rho < Real(0.0)) {
          stats.reason = TypedCGStop::numerical_breakdown;
          break;
        }
        // Exact zero is a valid solution even when a tiny tolerance squared
        // underflows in the selected working precision.
        if (cg_real_zero(next_rho) || cg_real_leq(next_rho, threshold)) {
          stats.reason = TypedCGStop::converged;
          break;
        }
        const Real beta = cg_real_div(next_rho, rho);
        stats.vector_flops += 1.0;
        if (!cg_real_finite(beta)) {
          stats.reason = TypedCGStop::numerical_breakdown;
          break;
        }
        qacc_for(i, length, {
          p[i] = cg_complex_add(r[i], cg_complex_scale(beta, p[i]));
        });
        stats.vector_flops += 4.0 * length;
        rho = next_rho;
      }
    }
    multi_count = stats.matrix_calls;
    it_res = stats.relative_residual();
    return stats;
  }

  // template<typename Tv, typename Ty>
  // inline void exit(Tv& dst, Ty* x, const Int ir = 0, qlat::vector<Ty >& RkL){
  //   ////multi_count
  //   dst.copy_from(x, ir, true);
  //   qmessage("QLAT::CG converged with resid %.8e, request %.8e, iter %5d, max
  //   %5d \n",
  //     RkL[1].real(), resid, multi_count, niter);
  // }

  //   may need to re-start as below GPT
  //   algorithms/inverter/defect_correcting.py #     lhs^{(0)} = inner_mat^{-1}
  //   rhs #     lhs^{(1)} = lhs^{(0)} - inner_mat^{-1} (outer_mat lhs^{(0)} -
  //   rhs) #     lhs^{(2)} = lhs^{(1)} - inner_mat^{-1} (outer_mat lhs^{(1)} -
  //   rhs)

  ////x is the results with preconditioning, p contains the initial src
  template <typename Ty>
  inline void iter(Ty* r, Ty* pA, Ty* p, Ty* x, const Int niter_local,
                   double resid_local)
  {
    M& mat = *matP;
    const long Ndata = mat.get_vsize();
    //
    qlat::vector<Ty> RkL;
    RkL.resize(2);
    Ty& Rk = RkL[0];
    Ty& Rkp1 = RkL[1];
    it_res = 0.0;
    //
    const Int GPU = 1;
    VectorGPUKey gkey(0, ssprintf("vec_norm2_buf"), GPU);
    vector_gpu<int8_t>& Buf = get_vector_gpu_plan<int8_t>(gkey);
    Buf.resizeL(size_t(Ndata) * sizeof(Ty));
    Ty* buf = (Ty*)Buf.data();
    //
    mat.multi(pA, x);
    multi_count += 1;
    qGPU_for(isp, Ndata, GPU, {
      r[isp] = p[isp] - pA[isp];
      p[isp] = r[isp];
      buf[isp] = qconj(r[isp]) * r[isp];
    });
    //
    Rkp1 = Reduce(buf, Ndata);
    Rk = Rkp1;
    if (std::sqrt(Rk.real()) < resid_local) {
      it_res = std::sqrt(Rkp1.real());
      return;
    }
    //
    for (Int ni = 0; ni < niter_local; ni++) {
      mat.multi(pA, p);
      multi_count += 1;
      Ty pAp = vec_norm2(p, pA, Ndata);
      Ty alphak = Rk / pAp;
      qGPU_for(isp, Ndata, GPU, {
        r[isp] -= alphak * pA[isp];
        buf[isp] = qconj(r[isp]) * r[isp];
        x[isp] += alphak * p[isp];
      });
      Rkp1 = Reduce(buf, Ndata);
      if (std::sqrt(Rkp1.real()) < resid_local) {
        it_res = std::sqrt(Rkp1.real());
        return;
      }
      //
      Ty betak = Rkp1 / Rk;
      qGPU_for(isp, Ndata, GPU, {
        // r[isp]  -= alphak * pA[isp];
        // buf[isp] = qconj(r[isp]) * r[isp];
        // x[isp]  += alphak * p[isp];
        p[isp] = r[isp] + betak * p[isp];
      });
      Rk = Rkp1;
    }
    it_res = std::sqrt(Rkp1.real());
    return;
  }

  template <typename Tv>
  inline void call(Tv& dst, Tv& src, double resid_, const Int niter_,
                   const Int ir = 0, const Int is = 0)
  {
    if (!dst.initialized) {
      Qassert(ir == 0);
      dst.resize(1, src);
    }
    M& mat = *matP;
    resid = resid_;
    niter = niter_;
    Qassert(mat.get_mrh() == 1);
    multi_count = 0;
    using Td = typename GetBasicDataType<Tv>::ElementaryType;
    vector<ComplexT<Td>*> mat_src;
    mat.set_buf_pointers(mat_src, 1);
    // auto& mat_src = mat.mat_src;
    Qassert(mat_src.size() >= 4);
    Qassert(mat.get_vsize() == long(src.nsum));
    /////const long Ndata = src.nsum;
    //
    auto* x = mat_src[0];
    auto* p = mat_src[1];
    auto* pA = mat_src[2];
    auto* r = mat_src[3];  //// 4 will be used for multi buffers
    //
    dst.copy_to(x, ir, true);
    src.copy_to(p, is, true);
    //
    ////x is the results with preconditioning, p contains the initial src
    Int niter_local = niter;
    double resid_local = resid;
    iter(r, pA, p, x, niter_local, resid_local);
    // exit(dst, x, ir, RkL);
    dst.copy_from(x, ir, true);
    if (verbose >= 1)
      qmessage(
          "QLAT::CG converged with resid %.8e, request %.8e, iter %5d, max %5d "
          "\n",
          it_res, resid, multi_count, niter);
    //
    // VectorGPUKey gkey(0, ssprintf("vec_norm2_buf"), GPU);
    // vector_gpu<int8_t >& Buf = get_vector_gpu_plan<int8_t >(gkey);
    // Buf.resizeL(size_t(Ndata)* sizeof(Ty));
    // Ty* buf = (Ty*) Buf.data();
    //
    ///////mat.multi(dst, src, ir, is); ////M mult(res, src)
    // mat.multi(pA , x);multi_count += 1;
    //
    // qGPU_for(isp, Ndata, {
    //   r[isp]   = p[isp] - pA[isp];
    //   p[isp]   = r[isp];
    //   buf[isp] = qconj(r[isp]) * r[isp];
    // });
    //
    // Rk = Reduce(buf, Ndata );
    // if(Rk < resid){exit(dst, x, ir);return ;}
    //
    // Rkp1 = Rk;
    //
    // for(Int ni=0;ni<niter;ni++)
    //{
    //   mat.multi(pA , p);multi_count += 1;
    //   Ty pAp = vec_norm2(p, pA, Ndata);
    //   Ty alphak = Rk / pAp;
    //   qGPU_for(isp, Ndata, {
    //     r[isp]  -= alphak * pA[isp];
    //     buf[isp] = qconj(r[isp]) * r[isp];
    //     x[isp]  += alphak * p[isp];
    //   });
    //   Rkp1 = Reduce(buf, Ndata );
    //   if(Rkp1 < resid){exit(dst, x, ir);return ;}
    //
    //   Ty betak = Rkp1 / Rk;
    //   qGPU_for(isp, Ndata, {
    //     //r[isp]  -= alphak * pA[isp];
    //     //buf[isp] = qconj(r[isp]) * r[isp];
    //     //x[isp]  += alphak * p[isp];
    //     p[k] = r[isp] + betak * p[k];
    //   });
    //   Rk = Rkp1;
    // }
    //
    // exit(dst, x, ir);return ;
  }
};

}  // namespace qlat

#endif
