// utils_cg_iter.h
// Gen Wang
// Jun. 2023

#ifndef UTILS_CG_ITER_H
#define UTILS_CG_ITER_H

#pragma once
#include "utils_float_type.h"

namespace qlat{

/*
  M is the method which provides Ty bufs, and multi on Ty; have nothing to do with vector_cs
  Tv is the src and result vectors, could have different precision as Ty
*/

template<class M, typename Tv>
inline double check_residue(Tv& res, Tv& src, M& mat, const Int ir = 0, const Int is = 0)
{
  const long Ndata = mat.get_vsize();
  Qassert(long(src.nsum) == Ndata);
  Int  GPU = 1;
  auto* ab = mat.mat_src[0];
  auto* ax = mat.mat_src[1];
  auto* ak = mat.mat_src[2];

  src.copy_to(ab, is, true);
  res.copy_to(ax, ir, true);
  mat.multi(ak, ax);
  qGPU_for(isp, Ndata, GPU, {ax[isp] = ak[isp] - ab[isp];});
  auto residu = vec_norm2(ax, ax, Ndata);
  qmessage("QLAT::mat residue %.8e \n", residu.real() );
  return residu.real();
}

template<class M>
struct cg_iter{
  //double resid;
  //Int niter;
  Int debug_mat;
  Int multi_count;
  double resid;
  Int niter;
  M* matP;
  M* matS;
  double it_res;
  Int verbose;
  //qlat::vector<qlat::Complex > RkL;
 
  cg_iter(M* matP_, M* matS_, const Int debug_mat_= 0)
  {
    //qmessage("setup chebyshev \n");
    matP = matP_;
    matS = matS_;
    //buf.resize(4);
    multi_count = 0;
    //
    ///RkL.resize(2);
    it_res = 0.0;
    resid = 0.0;
    niter = 0;
    verbose = qlat::get_env_long_default(std::string("qlat_cg_verbose"), 1);;
    debug_mat = debug_mat_;
  }

  //template<typename Tv, typename Ty>
  //inline void exit(Tv& dst, Ty* x, const Int ir = 0, qlat::vector<Ty >& RkL){
  //  ////multi_count
  //  dst.copy_from(x, ir, true);
  //  qmessage("QLAT::CG converged with resid %.8e, request %.8e, iter %5d, max %5d \n", 
  //    RkL[1].real(), resid, multi_count, niter);
  //}

  //   may need to re-start as below GPT algorithms/inverter/defect_correcting.py
  //   #     lhs^{(0)} = inner_mat^{-1} rhs
  //   #     lhs^{(1)} = lhs^{(0)} - inner_mat^{-1} (outer_mat lhs^{(0)} - rhs)
  //   #     lhs^{(2)} = lhs^{(1)} - inner_mat^{-1} (outer_mat lhs^{(1)} - rhs)
 
  ////x is the results with preconditioning, p contains the initial src
  template<typename Ty>
  inline void iter( Ty* r, Ty* pA, Ty* p, Ty* x , const Int niter_local, double resid_local){
    M& mat  = *matP;
    const long Ndata = mat.get_vsize();
    //
    qlat::vector<Ty > RkL;RkL.resize(2);
    Ty& Rk = RkL[0]; Ty& Rkp1 = RkL[1];
    it_res = 0.0;
    //
    const Int GPU = 1;
    VectorGPUKey gkey(0, ssprintf("vec_norm2_buf"), GPU);
    vector_gpu<int8_t >& Buf = get_vector_gpu_plan<int8_t >(gkey);
    Buf.resizeL(size_t(Ndata)* sizeof(Ty));
    Ty* buf = (Ty*) Buf.data();
    //
    mat.multi(pA , x);multi_count += 1;
    qGPU_for(isp, Ndata, GPU, {
      r[isp]   = p[isp] - pA[isp];
      p[isp]   = r[isp];
      buf[isp] = qconj(r[isp]) * r[isp];
    });
    //
    Rkp1 = Reduce(buf, Ndata );
    Rk   = Rkp1;
    if(std::sqrt(Rk.real()) < resid_local){it_res = std::sqrt(Rkp1.real());return ;}
    //
    for(Int ni=0;ni<niter_local;ni++)
    {
      mat.multi(pA , p);multi_count += 1;
      Ty pAp = vec_norm2(p, pA, Ndata);
      Ty alphak = Rk / pAp;
      qGPU_for(isp, Ndata, GPU, {
        r[isp]  -= alphak * pA[isp];
        buf[isp] = qconj(r[isp]) * r[isp];
        x[isp]  += alphak * p[isp];
      });
      Rkp1 = Reduce(buf, Ndata );
      if(std::sqrt(Rkp1.real()) < resid_local){it_res = std::sqrt(Rkp1.real());return ;}

      Ty betak = Rkp1 / Rk;
      qGPU_for(isp, Ndata, GPU, {
        //r[isp]  -= alphak * pA[isp];
        //buf[isp] = qconj(r[isp]) * r[isp];
        //x[isp]  += alphak * p[isp];
        p[isp] = r[isp] + betak * p[isp];
      }); 
      Rk = Rkp1;
    }
    it_res = std::sqrt(Rkp1.real());return ;
  }

  template<typename Tv>
  inline void call(Tv& dst, Tv& src, double resid_, const Int niter_, const Int ir = 0, const Int is = 0 ){
    if(!dst.initialized){Qassert(ir == 0);dst.resize(1, src);}
    M& mat  = *matP;resid = resid_; niter = niter_;
    Qassert(mat.get_mrh() == 1);
    multi_count = 0;
    using Td = typename GetBasicDataType<Tv>::ElementaryType;
    vector<ComplexT<Td>* > mat_src;
    mat.set_buf_pointers(mat_src, 1);
    //auto& mat_src = mat.mat_src;
    Qassert(mat_src.size()  >= 4);
    Qassert(mat.get_vsize() == long(src.nsum));
    /////const long Ndata = src.nsum;
    //
    auto* x  = mat_src[0];
    auto* p  = mat_src[1];
    auto* pA = mat_src[2];
    auto* r  = mat_src[3]; //// 4 will be used for multi buffers
    //
    dst.copy_to(x, ir, true);
    src.copy_to(p, is, true);
    //
    ////x is the results with preconditioning, p contains the initial src
    Int niter_local = niter; double resid_local = resid;
    iter( r, pA, p, x , niter_local, resid_local);
    //exit(dst, x, ir, RkL);
    dst.copy_from(x, ir, true);
    if(verbose >= 1)
    qmessage("QLAT::CG converged with resid %.8e, request %.8e, iter %5d, max %5d \n", 
      it_res, resid, multi_count, niter);
    //
    //VectorGPUKey gkey(0, ssprintf("vec_norm2_buf"), GPU);
    //vector_gpu<int8_t >& Buf = get_vector_gpu_plan<int8_t >(gkey);
    //Buf.resizeL(size_t(Ndata)* sizeof(Ty));
    //Ty* buf = (Ty*) Buf.data();
    //
    ///////mat.multi(dst, src, ir, is); ////M mult(res, src)
    //mat.multi(pA , x);multi_count += 1;
    //
    //qGPU_for(isp, Ndata, {
    //  r[isp]   = p[isp] - pA[isp];
    //  p[isp]   = r[isp];
    //  buf[isp] = qconj(r[isp]) * r[isp];
    //});
    //
    //Rk = Reduce(buf, Ndata );
    //if(Rk < resid){exit(dst, x, ir);return ;}
    //
    //Rkp1 = Rk;
    //
    //for(Int ni=0;ni<niter;ni++)
    //{
    //  mat.multi(pA , p);multi_count += 1;
    //  Ty pAp = vec_norm2(p, pA, Ndata);
    //  Ty alphak = Rk / pAp;
    //  qGPU_for(isp, Ndata, {
    //    r[isp]  -= alphak * pA[isp];
    //    buf[isp] = qconj(r[isp]) * r[isp];
    //    x[isp]  += alphak * p[isp];
    //  });
    //  Rkp1 = Reduce(buf, Ndata );
    //  if(Rkp1 < resid){exit(dst, x, ir);return ;}
    //
    //  Ty betak = Rkp1 / Rk;
    //  qGPU_for(isp, Ndata, {
    //    //r[isp]  -= alphak * pA[isp];
    //    //buf[isp] = qconj(r[isp]) * r[isp];
    //    //x[isp]  += alphak * p[isp];
    //    p[k] = r[isp] + betak * p[k];
    //  }); 
    //  Rk = Rkp1;
    //}
    //
    //exit(dst, x, ir);return ;
  }
};

}

#endif

