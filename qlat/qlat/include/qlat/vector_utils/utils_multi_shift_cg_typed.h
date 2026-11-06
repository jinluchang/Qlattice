#ifndef UTILS_MULTI_SHIFT_CG_TYPED_H
#define UTILS_MULTI_SHIFT_CG_TYPED_H
#pragma once

#include "utils_cg_precision.h"
#include <limits>
#include <stdexcept>
#include <vector>

namespace qlat {

// Local multishift CG. Unlike the historical custom solver, every recurrence
// scalar and reduction stays in Real, including DD low components. There is
// no residual replacement, restart, correction, or implicit precision floor.
template<class Real> struct TypedMultiShiftCGResult {
  Int iterations=0, matrix_calls=0;
  double vector_flops=0;
  std::vector<Int> shift_iterations,last_poll;
  std::vector<char> shift_reason,retired;
  std::vector<double> iterative_residual;
  std::vector<Real> source_norm2,residual_norm2,threshold2;
};

template<class Real> struct TypedMultiShiftCG {
  using C=ComplexT<Real>;
  Long length=0;
  Int rhs=0,shifts=0;
  vector_gpu<C> xs,ps,r,next_r,ap,difference,trial_x,trial_p,invalid;
  vector_gpu<Real> coefficients;
  vector_gpu<Int> enabled;
  CGTypedBatchReduction<Real> reduction,validity_reduction;

  void prepare(Long n,Int nrhs,Int nshift) {
    if(n<=0 || nrhs<=0 || nshift<=0 ||
       size_t(nrhs)>size_t(std::numeric_limits<Int>::max())/size_t(nshift) ||
       size_t(n)>size_t(std::numeric_limits<Long>::max())/size_t(nrhs)/size_t(nshift))
      throw std::invalid_argument("typed multishift CG requires positive fitting workspace dimensions");
    length=n;rhs=nrhs;shifts=nshift;
    const Long N=n*nrhs,S=Long(nrhs)*nshift;
    for(auto* v:{&xs,&ps,&trial_x,&trial_p,&invalid}) v->resize(n*S);
    for(auto* v:{&r,&next_r,&ap,&difference}) v->resize(N);
    coefficients.resize(3*S+nrhs);enabled.resize(S+nrhs);
    reduction.prepare(n,nrhs);validity_reduction.prepare(n,Int(S));
  }
  C* solution_data() {return xs.data();}

  template<class Matrix>
  TypedMultiShiftCGResult<Real> solve(Matrix& matrix,const C* source,
      const std::vector<Real>& delta,double target,double retirement_target,
      Int cap,Int check_gap,bool quda_beta,const std::vector<char>& original_nonzero,
      bool epsilon_retirement=false) {
    if(!source || length<=0 || rhs<=0 || shifts<=0 || delta.size()!=size_t(shifts) ||
       original_nonzero.size()!=size_t(rhs) || cap<0 || check_gap<=0 ||
       !std::isfinite(target) || target<=0 || !std::isfinite(retirement_target) || retirement_target<=0)
      throw std::invalid_argument("invalid typed multishift CG inputs");
    for(Int j=0;j<shifts;++j)
      if(!cg_real_finite(delta[j]) || delta[j]<Real(0) || (j==0 && !cg_real_zero(delta[j])) ||
         (j && delta[j]<delta[j-1]))
        throw std::invalid_argument("typed multishift CG requires ordered nonnegative relative shifts starting at zero");
    const Int R=rhs;
    const Long L=length,N=L*R,S=Long(R)*shifts;
    C* x=xs.data();C* p=ps.data();C* residual=r.data();C* nr=next_r.data();C* applied=ap.data();
    C* diff=difference.data();C* tx=trial_x.data();C* tp=trial_p.data();C* bad=invalid.data();
    Real* c=coefficients.data();Int* flags=enabled.data();
    qacc_for(i,N,{residual[i]=source[i];});
    qacc_for(i,L*S,{x[i]=C(Real(0),Real(0));p[i]=source[i%N];});
    auto rho=reduction.dot(residual,residual,L,rhs);
    const auto b2=rho;
    const Real tolerance=Real(target),retirement=Real(retirement_target);
    const Real t2=cg_real_mul(tolerance,tolerance),rt2=cg_real_mul(retirement,retirement);
    TypedMultiShiftCGResult<Real> out;
    out.source_norm2=b2;out.residual_norm2.resize(S);out.threshold2.resize(S);
    out.shift_iterations.assign(S,0);out.last_poll.assign(S,-1);
    out.shift_reason.assign(S,'L');out.retired.assign(S,0);out.iterative_residual.assign(S,1);
    out.vector_flops=4.0*N+2;
    std::vector<char> alive(rhs,1),active(S,1);
    std::vector<Real> z(S,Real(1)),previous_z(S,Real(1)),previous_alpha(rhs,Real(1)),previous_beta(rhs,Real(0));
    std::vector<Real> host_coeff(3*S+rhs,Real(0));
    std::vector<Int> host_flags(S+rhs,0);
    auto finish_seed=[&](Int h,char reason,Int stopping_iteration=-1) {
      alive[h]=0;
      for(Int j=0;j<shifts;++j) {
        const Long v=Long(j)*rhs+h;
        if(active[v]) {
          active[v]=0;out.shift_reason[v]=reason;
          if(stopping_iteration>=0) out.shift_iterations[v]=stopping_iteration;
        }
      }
    };
    auto relative=[&](const Real& r2,const Real& norm) {
      if(!cg_real_finite(r2) || r2<Real(0) || !cg_real_positive(norm))
        return std::numeric_limits<double>::infinity();
      return std::sqrt(double(cg_real_div(r2,norm)));
    };
    for(Int h=0;h<rhs;++h) {
      for(Int j=0;j<shifts;++j) {
        const Long v=Long(j)*rhs+h;
        out.threshold2[v]=cg_real_mul(j?rt2:t2,b2[h]);out.residual_norm2[v]=b2[h];
      }
      if(!cg_real_finite(b2[h]) || b2[h]<Real(0) || !cg_real_finite(t2) || !cg_real_finite(rt2) ||
         !cg_real_finite(out.threshold2[h]) || !cg_real_finite(out.threshold2[S-rhs+h])) {
        finish_seed(h,'B');
      } else if(cg_real_zero(b2[h])) {
        finish_seed(h,original_nonzero[h]?'B':'Z');
        for(Int j=0;j<shifts;++j) out.iterative_residual[Long(j)*rhs+h]=original_nonzero[h]?
          std::numeric_limits<double>::infinity():0.0;
      } else if(b2[h]<out.threshold2[h]) finish_seed(h,'I');
      host_flags[S+h]=alive[h];
    }
    cpy_GPU(flags,host_flags.data(),S+rhs,1,0);
    qacc_for(i,N,{if(!flags[S+i/L]) p[i]=C(Real(0),Real(0));});
    for(Int k=0;k<cap;++k) {
      bool any=false;for(char a:alive) any=any || a;if(!any) break;
      matrix.multi_batch(applied,p,rhs);++out.matrix_calls;out.iterations=k+1;
      auto curvature=reduction.dot(p,applied,L,rhs);out.vector_flops+=4.0*N;
      for(Int h=0;h<rhs;++h) {
        Real alpha=Real(0);
        if(alive[h]) {
          if(!cg_real_finite(curvature[h]) || !cg_real_positive(curvature[h])) finish_seed(h,'B',k+1);
          else {alpha=cg_real_div(rho[h],curvature[h]);if(!cg_real_finite(alpha)) finish_seed(h,'B',k+1);}
          out.vector_flops+=1;
        }
        host_coeff[3*S+h]=alpha;host_flags[S+h]=alive[h];
      }
      cpy_GPU(c,host_coeff.data(),3*S+rhs,1,0);cpy_GPU(flags,host_flags.data(),S+rhs,1,0);
      qacc_for(i,N,{
        if(flags[S+i/L]) {
          nr[i]=cg_complex_sub(residual[i],cg_complex_scale(c[3*S+i/L],applied[i]));
          if(quda_beta) diff[i]=cg_complex_sub(nr[i],residual[i]);
        } else {nr[i]=C(Real(0),Real(0));diff[i]=nr[i];}
      });
      out.vector_flops+=4.0*N;
      auto next_rho=reduction.dot(nr,nr,L,rhs);out.vector_flops+=4.0*N;
      std::vector<Real> beta_numerator;
      if(quda_beta) {beta_numerator=reduction.dot(nr,diff,L,rhs);out.vector_flops+=6.0*N;}
      const bool poll=(k+1)%check_gap==0;
      std::vector<Real> beta(rhs,Real(0)),new_z=z;
      for(Int h=0;h<rhs;++h) if(alive[h]) {
        if(!cg_real_finite(next_rho[h]) || next_rho[h]<Real(0)) {finish_seed(h,'B',k+1);continue;}
        beta[h]=cg_real_div(quda_beta?beta_numerator[h]:next_rho[h],rho[h]);out.vector_flops+=1;
        if(!cg_real_finite(beta[h])) {finish_seed(h,'B',k+1);continue;}
        const Real alpha=host_coeff[3*S+h];
        for(Int j=0;j<shifts;++j) {
          const Long v=Long(j)*rhs+h;if(!active[v]) continue;
          Real zn=Real(1),ratio=Real(1);
          if(j) {
            const Real den=cg_real_add(
              cg_real_mul(cg_real_mul(alpha,previous_beta[h]),cg_real_sub(previous_z[v],z[v])),
              cg_real_mul(cg_real_mul(previous_z[v],previous_alpha[h]),
                cg_real_add(Real(1),cg_real_mul(delta[j],alpha))));
            if(!cg_real_finite(den) || cg_real_zero(den) || cg_real_zero(z[v])) {
              active[v]=0;out.shift_reason[v]='B';out.shift_iterations[v]=k+1;continue;
            }
            zn=cg_real_div(cg_real_mul(cg_real_mul(z[v],previous_z[v]),previous_alpha[h]),den);
            if(!cg_real_finite(zn) || (cg_real_zero(zn) && !cg_real_zero(next_rho[h]))) {
              active[v]=0;out.shift_reason[v]='B';out.shift_iterations[v]=k+1;continue;
            }
            ratio=cg_real_div(zn,z[v]);out.vector_flops+=14;
          }
          const Real as=cg_real_mul(alpha,ratio),bs=cg_real_mul(beta[h],cg_real_mul(ratio,ratio));
          if(!cg_real_finite(ratio) || !cg_real_finite(as) || !cg_real_finite(bs)) {
            active[v]=0;out.shift_reason[v]='B';out.shift_iterations[v]=k+1;continue;
          }
          new_z[v]=zn;host_coeff[v]=zn;host_coeff[S+v]=as;host_coeff[2*S+v]=bs;
          out.vector_flops+=4;
        }
      }
      for(Long v=0;v<S;++v) host_flags[v]=active[v];
      cpy_GPU(c,host_coeff.data(),3*S+rhs,1,0);cpy_GPU(flags,host_flags.data(),S+rhs,1,0);
      // Validate a whole mass/RHS candidate before committing it, so a failed
      // shift cannot corrupt its saved solution or the other shifted streams.
      qacc_for(i,L*S,{
        const Long v=i/L,h=v%R;
        if(flags[v]) {
          tx[i]=cg_complex_add(x[i],cg_complex_scale(c[S+v],p[i]));
          tp[i]=cg_complex_add(cg_complex_scale(c[v],nr[h*L+i%L]),cg_complex_scale(c[2*S+v],p[i]));
          bad[i]=C(Real(cg_complex_finite(tx[i]) && cg_complex_finite(tp[i])?0:1),Real(0));
        } else bad[i]=C(Real(0),Real(0));
      });
      const auto invalid_count=validity_reduction.dot(bad,bad,L,Int(S));out.vector_flops+=4.0*L*S;
      for(Int h=0;h<rhs;++h) if(alive[h] && !cg_real_zero(invalid_count[h])) finish_seed(h,'B',k+1);
      for(Long v=0;v<S;++v) if(active[v]) {
        if(!cg_real_zero(invalid_count[v])) {active[v]=0;out.shift_reason[v]='B';out.shift_iterations[v]=k+1;}
        else {
          out.shift_iterations[v]=k+1;
          out.residual_norm2[v]=cg_real_mul(cg_real_mul(next_rho[v%rhs],new_z[v]),new_z[v]);
          out.iterative_residual[v]=relative(out.residual_norm2[v],b2[v%rhs]);
          out.vector_flops+=10.0*L;
        }
      }
      for(Long v=0;v<S;++v) host_flags[v]=active[v];
      for(Int h=0;h<rhs;++h) host_flags[S+h]=alive[h];
      cpy_GPU(flags,host_flags.data(),S+rhs,1,0);
      qacc_for(i,L*S,{if(flags[i/L]) {x[i]=tx[i];p[i]=tp[i];}});
      qacc_for(i,N,{if(flags[S+i/L]) residual[i]=nr[i];});
      // Poll after the boundary update. Strict '<' intentionally preserves a
      // squared threshold of zero rather than silently replacing it by epsilon.
      for(Int h=0;h<rhs;++h) if(alive[h]) {
        for(Int j=0;j<shifts;++j) {
          const Long v=Long(j)*rhs+h;if(!active[v]) continue;
          previous_z[v]=z[v];z[v]=new_z[v];
          if(poll) {
            out.last_poll[v]=k+1;
            if(j && out.residual_norm2[v]<out.threshold2[v]) {
              active[v]=0;out.retired[v]=1;out.shift_reason[v]=epsilon_retirement?'E':'T';
            }
          }
        }
        rho[h]=next_rho[h];previous_alpha[h]=host_coeff[3*S+h];previous_beta[h]=beta[h];
        if(cg_real_zero(rho[h])) finish_seed(h,'X');
        else if(poll && rho[h]<out.threshold2[h]) finish_seed(h,'T');
      }
      for(Int h=0;h<rhs;++h) host_flags[S+h]=alive[h];
      cpy_GPU(flags,host_flags.data(),S+rhs,1,0);
      qacc_for(i,N,{if(!flags[S+i/L]) p[i]=C(Real(0),Real(0));});
    }
    return out;
  }
};

} // namespace qlat
#endif
