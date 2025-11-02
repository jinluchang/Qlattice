// utils_vector_GPU.h
// Gen Wang
// Sep. 2021

#ifndef UTILS_VECTOR_GPU_H
#define UTILS_VECTOR_GPU_H

#pragma once
#include <qlat/qcd.h>
#include "utils_float_type.h"
#include "utils_read_txt.h"
#include "utils_COPY_data.h"

////needed for norm calculation
#include "utils_reduce_vec.h"

namespace qlat{

template <typename Ty >
struct vector_gpu{
  ///Vector<Ty > v;
  Ty*    p;
  size_t n;
  Int GPU;///1 for GPU, 0 for CPU, -1 for unified
  bool is_copy;

  vector_gpu()
  {
    p = NULL; n = 0;GPU = 1;is_copy = false;
  }

  vector_gpu(const size_t n_set, const Int GPU_set = 1)
  {
    p = NULL; n = 0;GPU = 1;is_copy = false;
    resize(n_set, GPU_set);
  }

  inline void resize(const size_t n_set, const Int GPU_set)
  {
    Qassert(not is_copy);
    Qassert(GPU_set == -1 or GPU_set == 0 or GPU_set == 1);
    if(n_set == 0){clean_mem(); return ;}
    if((n != n_set) or (GPU != GPU_set))
    {
      /////qmessage("Reisze n %d, n_set %d !\n", int(n), int(n_set));
      clean_mem();
      n = n_set;
      GPU    = GPU_set;
      gpuMalloc(p, size_t(n), Ty, GPU);
      set_zero();
    }
    ////qacc_barrier(dummy);
  }

  vector_gpu(const vector_gpu<Ty>& vp)
  {
    #ifndef QLAT_USE_ACC
        Qassert(false);
    #endif
    is_copy = true;
    p = vp.p;
    n = vp.n;
    GPU = vp.GPU;
  }

  vector_gpu(vector_gpu<Ty>&& vp) noexcept
  {
    is_copy = vp.is_copy;
    p = vp.p;
    n = vp.n;
    GPU = vp.GPU;
    vp.is_copy = true;
  }

  inline void resize(const size_t n_set)
  {
    Int GPU_tem = GPU;
    resize(n_set, GPU_tem);
  }

  inline void resizeL(const size_t n_set, const Int GPU_ = -2)
  {
    Int GPU_cur = GPU;if(GPU_ != -2){GPU_cur = GPU_;}
    if(GPU_cur != GPU){resize(n_set, GPU_cur); return ;}
    if(n < n_set){resize(n_set, GPU_cur);}
  }

  qacc size_t size() const{return n;}
  qacc Ty* data(){return p;}
  qacc const Ty* data() const { return p; }
  qacc const Ty& operator[](const size_t i) const { return p[i]; }
  qacc Ty& operator[](const size_t i) { return p[i]; }

  inline void set_zero(QBOOL dummy=QTRUE)
  {
    zero_Ty(p, n, GPU, dummy);
  }

  inline void set_zero_pt(QBOOL dummy=QTRUE)
  {
    (void)dummy;
    Ty* pr = p;
    Int GPU_ = GPU;
    Long n_ = n;
    qGPU_for(isp, n_, GPU_, {pr[isp] = NULL;});
  }

  inline void clean_mem(){
    free_buf(p, GPU);
    p = NULL;n = 0;
  }
  void clear(){
    Qassert(not is_copy);
    clean_mem();
  }

  void clear_copy(){
    if(is_copy){
      p = NULL;
      n = 0   ;
      GPU = 0;
      is_copy = false;
    }
    else{
      clear();
    }
  }

  ~vector_gpu(){
    if (not is_copy){
      clean_mem();
    }
  }

  //template <class T >
  //vector_gpu<Ty>& operator=(const vector_gpu<T >& vp)
  //{
  //  qmessage("NO SUPPORT yet!\n");
  //  Qassert(false);
  //  ////bool tem_GPU = vp.GPU;
  //  resize(vp.size(), vp.GPU);

  //  Int mode_cpu = 0;
  //  if(vp.GPU == false and GPU == false){mode_cpu =  0;}
  //  if(vp.GPU == true  and GPU == true ){mode_cpu =  1;}

  //  cpy_data_thread(p, vp.p, n, mode_cpu, true);

  //  return *this;
  //}

  /////template <class T >
  const vector_gpu<Ty>& operator=(const vector_gpu<Ty >& vp)
  {
    Qassert(not is_copy);
    resize(vp.size(), vp.GPU);
    //int mode_cpu = 0;
    //if(vp.GPU == false and GPU == false){mode_cpu =  0;}
    //if(vp.GPU == true  and GPU == true ){mode_cpu =  1;}
    //cpy_data_thread(p, vp.p, n, mode_cpu, true);
    cpy_GPU(p, vp.p, n, GPU, vp.GPU, QTRUE);
    return *this;
  }

  ////
  template <class T >
  void copy_from(const T* src, size_t Ndata, Int GPU_set = -2, Int GPU_ori = 0)
  {
    //Qassert(GPU_ori == 0 or GPU_ori == 1);
    Qassert(GPU_ori == 0 or GPU_ori == 1 or GPU_ori == -1);
    Int GPU_src = GPU_ori;
    //if(GPU_ori == 0 ){GPU_src = false ;}
    //if(GPU_ori == 1 ){GPU_src = true  ;}

    Int tem_GPU = GPU;
    if(GPU_set != -2){tem_GPU = GPU_set;}
    resize(Ndata, tem_GPU);
  
    //int mode_cpu = 0;
    //if(GPU_src == false and GPU == false){mode_cpu =  0;} // host to host
    //if(GPU_src == true  and GPU == true ){mode_cpu =  1;} // device to device
    //if(GPU_src == false and GPU == true ){mode_cpu =  2;} // host to device
    //if(GPU_src == true  and GPU == false){mode_cpu =  3;} // device to host
    //cpy_data_thread(p, src, n, mode_cpu, true);

    cpy_GPU(p, src, n, GPU, GPU_src, QTRUE);

  }

  template <class T >
  void copy_from(const vector_gpu<T >& vp, Int GPU_set = -2)
  {
    copy_from(vp.p, vp.n, GPU_set, int(vp.GPU));
  }

  template <class T >
  void copy_from(const std::vector<T >& vp, Int GPU_set = -2)
  {
    Int GPU_mem = 0;
    if(GPU_set == -2){GPU_mem = GPU;}else{GPU_mem = GPU_set;}
    copy_from(&vp[0], vp.size(), GPU_mem, 0);
  }


  template <class T >
  void copy_from(qlat::vector<T >& vp, Int GPU_set = -2, Int GPU_ori = 0)
  {
    Int GPU_mem = 0;
    if(GPU_set == -2){GPU_mem = GPU;}else{GPU_mem = GPU_set;}
    T* src = (T*) qlat::get_data(vp).data();
    copy_from(src, vp.size(), GPU_mem, GPU_ori);
  }

  template <class T >
  void copy_to(T* res, Int GPU_ori = -2)
  {
    ////int mode_cpu = 0;
    Int GPU_set = GPU_ori;if(GPU_ori == -2){GPU_set  =  1;}
    //if(GPU == false and GPU_set == 0){mode_cpu =  0;} // host to host
    //if(GPU == true  and GPU_set == 1){mode_cpu =  1;} // device to device
    //if(GPU == true  and GPU_set == 0){mode_cpu =  3;} // device to host
    //if(GPU == false and GPU_set == 1){mode_cpu =  2;} // host to device
    //cpy_data_thread(res, p, n, mode_cpu, true);
    cpy_GPU(res, p, n, GPU_set, GPU, QTRUE);
  }

  template <class T >
  void copy_to(vector_gpu<T >& vp, Int GPU_set = -2)
  {
    Int tem_GPU =  GPU;
    if(GPU_set != -2){tem_GPU = GPU_set;}
    //if(GPU_set == -2){tem_GPU = GPU;}
    //else{tem_GPU = GPU_set;}
    //if(GPU_set == -1 ){tem_GPU = true ;}
    //if(GPU_set ==  0 ){tem_GPU = false ;}
    //if(GPU_set ==  1 ){tem_GPU = true  ;}
    vp.resize(size(), tem_GPU);
    copy_to(vp.p, GPU_set);
  }

  inline Ty norm2()
  {
    qlat::vector_gpu<Ty > copy;copy.resize(n, GPU);
    Ty* res = copy.data();Ty* src = p;
    //if(GPU){
    //  qacc_for(isp, Long(n),    {res[isp] = qlat::qconj(src[isp]) * src[isp];});
    //}
    //else{
    //  qthread_for(isp, Long(n), {res[isp] = qlat::qconj(src[isp]) * src[isp];});
    //}
    qGPU_for(isp, Long(n), GPU, { res[isp] = qlat::qconj(src[isp]) * src[isp]; });
    qlat::vector<Ty > tmp;tmp.resize(2);tmp[0] = 0;
    reduce_vecs(res, tmp.data(), n, 1, GPU);
    glb_sum(tmp[0]);
    //MPI_Datatype curr = MPI_DOUBLE;unsigned int M_size = sizeof(double);
    //M_size = get_MPI_type<Ty >(curr);
    //MPI_Allreduce(&tmp[0], &tmp[1], sizeof(Ty)/M_size, curr, MPI_SUM, get_comm());
    return tmp[0];
  }
  
  inline void print_norm2(std::string prec = std::string("%.8e"))
  {
    Ty normC = norm2();
    std::string ktem = ssprintf( prec.c_str(), normC.real());
    qmessage("==norm %s \n", ktem.c_str());
  }

  template <class T >
  void swap(std::vector<T >& vp)
  {
    Qassert(not is_copy);
    Qassert(not vp.is_copy);
    Ty*  p_tmp   = vp.p;
    size_t n_tmp = vp.n;
    Int GPU_tmp = vp.GPU;

    ////copy to vp
    vp.p = p;
    vp.n = n;
    vp.GPU = GPU;

    ////copy to self
    p = p_tmp;
    n = n_tmp;
    GPU = GPU_tmp;
  }

  template <class T >
  void copy_to(std::vector<T >& vp)
  {
    vp.resize(n);
    //int mode_cpu = 0;
    //if(GPU == false){mode_cpu =  0;}
    //if(GPU == true ){mode_cpu =  3;} // device to host
    //cpy_data_thread(&vp[0], p, n, mode_cpu, true);
    cpy_GPU(&vp[0], p, n, 0, GPU, QTRUE);
  }

  template <class T >
  void copy_to(qlat::vector<T >& vp, Int GPU_ori = 0)
  {
    vp.resize(size());
    T* res = (T*) qlat::get_data(vp).data();
    copy_to(res, GPU_ori);
  }

  template <class T >
  const vector_gpu<Ty>& operator+=(const vector_gpu<T >& vp)
  {
    Qassert(GPU == vp.GPU and n == vp.n);
    //int mode_cpu = 0;
    //if(vp.GPU == false and GPU == false){mode_cpu =  0;}
    //if(vp.GPU == true  and GPU == true ){mode_cpu =  1;}
    //cpy_data_thread(p, vp.p, n, mode_cpu, true, 1.0);
    //cpy_GPU(p, vp.p,n, GPU, vp.GPU, true, 1.0);
    const Long N = n;
    const Int GPU_ = GPU;
    T*  r = p;
    Ty* s = vp.p;
    qGPU_for(isp, N, GPU_, { r[isp] += s[isp]; } );
    return *this;
  }

  template <class T >
  const vector_gpu<Ty>& operator-=(const vector_gpu<T >& vp)
  {
    Qassert(GPU == vp.GPU and n == vp.n);
    //int mode_cpu = 0;
    //if(vp.GPU == false and GPU == false){mode_cpu =  0;}
    //if(vp.GPU == true  and GPU == true ){mode_cpu =  1;}
    //cpy_data_thread(p, vp.p, n, mode_cpu, true, -1.0);
    //cpy_GPU(p, vp.p,n, GPU, vp.GPU, true, -1.0);
    const Long N = n;
    const Int GPU_ = GPU;
    T*  r = p;
    Ty* s = vp.p;
    qGPU_for(isp, N, GPU_, { r[isp] -= s[isp]; } );
    return *this;
  }

  template <class T >
  const vector_gpu<Ty>& operator*=(const T& f)
  {
    Ty* data = p;
    const Int GPU_ = GPU;
    qGPU_for(isp, Long(n), GPU_, {
      data[isp] = f * data[isp];
    });
    return *this;
  }

};

////y = alpha x + beta y
template <typename Ty, typename T >
void qblas_xAXPY(vector_gpu<Ty>& y, vector_gpu<Ty>& x, T& alpha, T& beta = 0.0)
{
  Qassert(x.GPU == y.GPU and x.n == y.n);
  Ty* xp = x.p;
  Ty* yp = y.p;
  const Long N = x.n;
  const Int GPU = x.GPU;
  if(beta != 0.0){
    qGPU_for(isp, N, GPU, {
      yp[isp]  = alpha * xp[isp] + beta * yp[isp];
    });
  }else{
    qGPU_for(isp, N, GPU, {
      yp[isp]  = alpha * xp[isp];
    });
  }
}

////dot = x^T y
template <typename Ty >
Ty qblas_dot(vector_gpu<Ty>& y, vector_gpu<Ty>& x, bool conj = true)
{
  Qassert(x.GPU == y.GPU and x.n == y.n);
  Ty* xp = x.p;
  Ty* yp = y.p;
  const Long N = x.n;
  const Int GPU = x.GPU;
  vector_gpu<Ty > buf;buf.resize(N, GPU);
  if(conj){
    qGPU_for(isp, N, GPU, {
      buf[isp] = qconj(xp[isp]) * yp[isp];
    });
  }else{
    qGPU_for(isp, N, GPU, {
      buf[isp] = xp[isp] * yp[isp];
    });
  }
  return buf.norm2();
}

template <typename Ty >
qacc Vector<Ty> get_data(vector_gpu<Ty>& vec)
{
  return Vector<Ty>(vec.data(), vec.size());
}

template <typename Ty >
qacc void set_zero(vector_gpu<Ty>& vec)
{
  vec.set_zero();
}

/////vector buffers related
struct VectorGPUKey {
  std::string tag;
  size_t size;
  Int GPU;
  VectorGPUKey()
  {
    size = 0; GPU = false;tag = std::string("");
  }
  VectorGPUKey(size_t size_, const std::string tag_, Int GPU_)
  {
    size = size_; GPU = GPU_;tag = tag_;
  }

};

inline bool operator<(const VectorGPUKey& x, const VectorGPUKey& y)
{
  if(x.GPU < y.GPU ){  return true;}
  if(y.GPU < x.GPU ){  return false;}

  if(x.tag.compare(y.tag) < 0){return true;}
  if(y.tag.compare(x.tag) < 0){return false;}

  return false;
}

template <typename Ty >
inline Cache<VectorGPUKey, vector_gpu<Ty > >& get_vector_gpu_cache()
{
  static Cache<VectorGPUKey, vector_gpu<Ty > > cache("VectorGPUKey", 128);
  return cache;
}

template <typename Ty >
inline vector_gpu<Ty >& get_vector_gpu_plan(const VectorGPUKey& gkey)
{
  if (!get_vector_gpu_cache<Ty>().has(gkey)) {
    get_vector_gpu_cache<Ty>()[gkey] = vector_gpu<Ty >(); 
  }
  vector_gpu<Ty >& buf = get_vector_gpu_cache<Ty>()[gkey];

  buf.GPU = gkey.GPU;
  if(buf.size() < gkey.size){
    buf.resize(gkey.size);
  }
  return buf;
}

template <typename Ty >
inline void safe_free_vector_gpu_plan(const VectorGPUKey& gkey, const bool zero = false)
{
  if ( get_vector_gpu_cache<Ty>().has(gkey)) {
    vector_gpu<Ty >& buf = get_vector_gpu_cache<Ty>()[gkey];
    size_t MAX = MAX_VECTOR_GPU_BUF;
    std::string val = get_env(std::string("q_mem_vec_gpu"));
    if(val != ""){MAX = stringtonum(val);}
    MAX = MAX * 1024 * 1024;  ////to MB

    if(!zero)
    if(buf.size() * sizeof(Ty) > MAX){
      buf.resize(MAX);
    }
    if(zero){buf.resize(0);}
  }
}

template <typename Ty >
inline vector_gpu<Ty >& get_vector_gpu_plan(size_t vol, const std::string& info, const Int GPU)
{
  VectorGPUKey gkey(vol, info, GPU);
  vector_gpu<Ty >& buf = get_vector_gpu_plan<Ty >(gkey);
  return buf;
}

template <typename Ty >
inline void safe_free_vector_gpu_plan(const std::string& info, const Int GPU, const bool zero = false)
{
  VectorGPUKey gkey(0, info, GPU);
  safe_free_vector_gpu_plan<Ty >(gkey, zero);
}

inline void clear_vector_gpu_cache()
{
  get_vector_gpu_cache<int8_t >().clear();
  get_vector_gpu_cache<float >().clear();
  get_vector_gpu_cache<double >().clear();
  get_vector_gpu_cache<qlat::ComplexT<float> >().clear();
  get_vector_gpu_cache<qlat::ComplexT<double> >().clear();
 
  //Cache<VectorGPUKey, vector_gpu<int8_t > >& c0 = get_vector_gpu_cache<int8_t >();
  //Cache<VectorGPUKey, vector_gpu<float > >& c1 = get_vector_gpu_cache<float >();
  //Cache<VectorGPUKey, vector_gpu<double > >& c2 = get_vector_gpu_cache<double >();
  //Cache<VectorGPUKey, vector_gpu<qlat::ComplexT<float > > >& c3 = get_vector_gpu_cache<qlat::ComplexT<float> >();
  //Cache<VectorGPUKey, vector_gpu<qlat::ComplexT<double> > >& c4 = get_vector_gpu_cache<qlat::ComplexT<double> >();
  //c0.clear();
  //c1.clear();
  //c2.clear();
  //c3.clear();
  //c4.clear();
}


}


#endif
