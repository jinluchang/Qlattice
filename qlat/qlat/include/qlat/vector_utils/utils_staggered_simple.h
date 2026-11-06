#ifndef UTILS_STAGGERED_SIMPLE_H
#define UTILS_STAGGERED_SIMPLE_H
#pragma once

#include "utils_quda_inverter.h"
#include "utils_cg_precision.h"
#include <comm_quda.h>
#include <stdexcept>

namespace qlat {

// This operator deliberately has no halo exchange. Check participating ranks,
// rather than qacc_GetDeviceCount(): unused visible GPUs are harmless.
inline void check_staggered_simple_single_gpu_once() {
  static const bool checked=[] {
    if(get_num_node()!=1 || quda::comm_size()!=1)
      abort_r("Standalone staggered Dslash requires one participating Qlat rank and one QUDA rank (single GPU); halo exchange is not implemented.");
    return true;
  }();
  (void)checked;
}

// One checkerboard site's full lexicographic index. With even X[0], adjacent
// lexicographic sites have opposite parity and site=full/2 is QUDA's CB index.
struct StaggeredSimpleGeometry {
  Int x[4]={0,0,0,0};
  Long sites=0;
  qacc Long full_index(Long site,Int parity) const {
    Long full=2*site,rem=full;Int sum=0;
    for(Int mu=0;mu<4;++mu) { sum+=Int(rem%x[mu]);rem/=x[mu]; }
    return full+((sum&1)!=parity);
  }
  qacc Long neighbor(Long full,Int mu,Int direction) const {
    Long stride=1;
    for(Int nu=0;nu<mu;++nu) stride*=x[nu];
    const Int coordinate=Int((full/stride)%x[mu]);
    if(direction>0) return full+stride*(coordinate==x[mu]-1?1-x[mu]:1);
    return full+stride*(coordinate==0?x[mu]-1:-1);
  }
};

// An independent, full-storage nearest-neighbor staggered operator:
// D v(x) = sum_mu [U_mu(x)v(x+mu)-U_mu(x-mu)^dagger v(x-mu)],
// A_m = 4*m*m - D_eo D_oe. Saved qinv links use MILC checkerboard order
// [parity][site][direction][row][column] and ALREADY include staggered phases,
// boundary factors and setup-time link scaling. Coordinate wrapping therefore
// adds no phase. Internal vectors use [color][checkerboard site]; public native
// QUDA vectors must enter/leave through import_public/export_public.
template<class Real> struct StaggeredSimpleOperator {
  static_assert(std::is_same<Real,float>::value || std::is_same<Real,double>::value || std::is_same<Real,RealDD>::value, "Unsupported solver precision");
  using Complex=ComplexT<Real>;
  using real_type=Real;
  StaggeredSimpleGeometry geometry;
  vector_gpu<Complex> links,odd,even;
  vector_gpu<Complex> batch_odd,batch_even;
  Int batch_rhs=0;
  vector_gpu<Complex> cg_buffers;
  Real shift=Real(0.0);
  double effective_shift=0;

  explicit StaggeredSimpleOperator(quda_inverter& qinv) {
    check_staggered_simple_single_gpu_once();
    if(qinv.inv_param.dslash_type!=QUDA_STAGGERED_DSLASH)
      throw std::invalid_argument("Standalone staggered operator supports ordinary nearest-neighbor staggered links only (ASQTAD is unsupported)");
    if(qinv.gauge_param.gauge_order!=QUDA_MILC_GAUGE_ORDER ||
       qinv.gauge_param.cpu_prec!=QUDA_DOUBLE_PRECISION)
      throw std::invalid_argument("Standalone staggered operator requires saved double MILC gauge links");
    if(qinv.inv_param.dagger!=QUDA_DAG_NO)
      throw std::invalid_argument("Standalone staggered operator requires dagger=no");
    const auto& input=qinv.qinv_twist_link==QTRUE?qinv.quda_gf_default_twist:qinv.quda_gf_default;
    initialize(qinv.geo(),input.data(),Long(input.size()));
    for(Int mu=0;mu<4;++mu) if(qinv.X[mu]!=geometry.x[mu])
      throw std::invalid_argument("Standalone gauge dimensions do not match qinv");
  }
  // Raw saved-gauge interface also supports standalone reference applications.
  StaggeredSimpleOperator(const Geometry& geo,const ComplexD* source,Long elements) {
    check_staggered_simple_single_gpu_once();initialize(geo,source,elements);
  }
  void initialize(const Geometry& geo,const ComplexD* source,Long elements) {
    const Coordinate dims=geo.node_site,total=geo.total_site();Long volume=1;
    for(Int mu=0;mu<4;++mu) {
      if(dims[mu]<=0 || dims[mu]%2 || dims[mu]!=total[mu])
        throw std::invalid_argument("Standalone staggered operator requires four even unpartitioned extents");
      geometry.x[mu]=dims[mu];volume*=dims[mu];
    }
    if(volume!=geo.local_volume() || !source || elements!=volume*36)
      throw std::invalid_argument("Standalone staggered operator requires a complete four-dimensional saved gauge");
    geometry.sites=volume/2;
    links.resize(elements);odd.resize(get_vsize());even.resize(get_vsize());
    import_gauge(source,elements);
  }
  void import_gauge(const ComplexD* source,Long count) {
    if(count!=Long(links.size()))
      throw std::invalid_argument("Standalone staggered gauge import has an incorrect element count");
    Complex* dest=links.data();
    qacc_for(i,count,{ dest[i]=Complex(Real(source[i].real()),Real(source[i].imag())); });
  }
  Long get_vsize() const { return geometry.sites*3; }
  Int get_mrh() const { return 1; }
  void set_shift(double offset) {
    if(!std::isfinite(offset) || offset<0)
      throw std::invalid_argument("Standalone staggered operator offset must be finite and nonnegative");
    const double mass=std::sqrt(offset/4.0);
    effective_shift=4.0*mass*mass;
    shift=Real(effective_shift);
  }
  void set_buf_pointers(vector<Complex*>& buffers,Int mrh=1) {
    if(mrh!=1) throw std::invalid_argument("Standalone staggered CG solves one RHS at a time");
    cg_buffers.resize(4*get_vsize());buffers.resize(4);
    for(Int i=0;i<4;++i) buffers[i]=cg_buffers.data()+i*get_vsize();
  }
  template<class Layout> void import_public(Complex* out,const ComplexD* in,const Layout& layout) const {
    const Long sites=geometry.sites;
    qacc_for(i,3*sites,{
      const Int color=Int(i/sites);const Long site=i%sites;
      const ComplexD value=in[layout.template public_index<ComplexD>(site,color,sites)];
      out[i]=Complex(Real(value.real()),Real(value.imag()));
    });
  }
  template<class Layout> void export_public(ComplexD* out,const Complex* in,const Layout& layout) const {
    const Long sites=geometry.sites;
    qacc_for(i,3*sites,{
      const Int color=Int(i/sites);const Long site=i%sites;
      out[layout.template public_index<ComplexD>(site,color,sites)]=ComplexD(double(in[i].real()),double(in[i].imag()));
    });
  }
  // output_parity identifies destination sites, source has opposite parity.
  // dslash itself requires disjoint buffers; multi() supports arbitrary aliases.
  void dslash(Complex* out,const Complex* in,Int output_parity) const {
    if(output_parity!=0 && output_parity!=1)
      throw std::invalid_argument("Standalone Dslash output parity must be even (0) or odd (1)");
    if(out==in) throw std::invalid_argument("Standalone parity Dslash requires disjoint input/output buffers");
    const auto g=geometry;const Long sites=g.sites;const Complex* gauge=links.data();
    qacc_for(i,3*sites,{
      const Int row=Int(i/sites);const Long site=i%sites;
      const Long full=g.full_index(site,output_parity);
      Complex sum(Real(0.0),Real(0.0));
      for(Int mu=0;mu<4;++mu) {
        const Long forward=g.neighbor(full,mu,1),backward=g.neighbor(full,mu,-1);
        const Long here=Long(output_parity)*sites+site;
        const Long back=Long(1-output_parity)*sites+backward/2;
        for(Int col=0;col<3;++col) {
          const Complex u=gauge[(here*4+mu)*9+row*3+col];
          const Complex v=gauge[(back*4+mu)*9+col*3+row];
          const Complex ud(v.real(),cg_real_neg(v.imag()));
          sum=cg_complex_add(sum,cg_complex_mul(u,in[Long(col)*sites+forward/2]));
          sum=cg_complex_sub(sum,cg_complex_mul(ud,in[Long(col)*sites+backward/2]));
        }
      }
      out[i]=sum;
    });
  }
  void multi(Complex* out,const Complex* in) {
    dslash(odd.data(),in,1);dslash(even.data(),odd.data(),0);
    const Complex* dd=even.data();const Real mass_term=shift;const Long count=get_vsize();
    qacc_for(i,count,{ out[i]=cg_complex_sub(cg_complex_scale(mass_term,in[i]),dd[i]); });
  }
  // Batch storage is independent of the scalar path. Allocate once before
  // entering a solve; all batch applications reuse these buffers.
  void prepare_batch(Int rhs) {
    if(rhs<=0)
      throw std::invalid_argument("Standalone staggered batch requires a positive RHS count");
    const Long count=Long(rhs)*get_vsize();
    batch_odd.resize(count);batch_even.resize(count);batch_rhs=rhs;
  }
  // Vectors use [RHS][color][checkerboard site], while every RHS shares links.
  void dslash_batch(Complex* out,const Complex* in,Int output_parity,Int rhs) const {
    if(output_parity!=0 && output_parity!=1)
      throw std::invalid_argument("Standalone batched Dslash output parity must be even (0) or odd (1)");
    if(rhs<=0)
      throw std::invalid_argument("Standalone staggered batch requires a positive RHS count");
    if(out==in)
      throw std::invalid_argument("Standalone batched parity Dslash requires disjoint input/output buffers");
    const auto g=geometry;const Long sites=g.sites,length=3*sites;
    const Complex* gauge=links.data();
    qacc_for(i,Long(rhs)*length,{
      const Long base=(i/length)*length,index=i%length;
      const Int row=Int(index/sites);const Long site=index%sites;
      const Long full=g.full_index(site,output_parity);
      Complex sum(Real(0.0),Real(0.0));
      for(Int mu=0;mu<4;++mu) {
        const Long forward=g.neighbor(full,mu,1),backward=g.neighbor(full,mu,-1);
        const Long here=Long(output_parity)*sites+site;
        const Long back=Long(1-output_parity)*sites+backward/2;
        for(Int col=0;col<3;++col) {
          const Complex u=gauge[(here*4+mu)*9+row*3+col];
          const Complex v=gauge[(back*4+mu)*9+col*3+row];
          const Complex ud(v.real(),cg_real_neg(v.imag()));
          sum=cg_complex_add(sum,cg_complex_mul(u,in[base+Long(col)*sites+forward/2]));
          sum=cg_complex_sub(sum,cg_complex_mul(ud,in[base+Long(col)*sites+backward/2]));
        }
      }
      out[i]=sum;
    });
  }
  void multi_batch(Complex* out,const Complex* in,Int rhs) {
    const Long count=Long(rhs)*get_vsize();
    if(rhs<=0 || rhs>batch_rhs || batch_odd.size()<size_t(count) ||
       batch_even.size()<size_t(count))
      throw std::invalid_argument("Standalone staggered batch must be prepared before application");
    dslash_batch(batch_odd.data(),in,1,rhs);
    dslash_batch(batch_even.data(),batch_odd.data(),0,rhs);
    Complex* result=batch_odd.data();const Complex* dd=batch_even.data();
    const Real mass_term=shift;
    qacc_for(i,count,{ result[i]=cg_complex_sub(cg_complex_scale(mass_term,in[i]),dd[i]); });
    // Finish in scratch before copying, including for overlapping public views.
    qacc_for(i,count,{ out[i]=result[i]; });
  }
};

} // namespace qlat
#endif
