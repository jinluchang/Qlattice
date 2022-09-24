// utils_vector_GPU.h
// Gen Wang
// Sep. 2021

#ifndef UTILS_VECTOR_GPU_H
#define UTILS_VECTOR_GPU_H

#pragma once
#include <qlat/qcd.h>
#include "utils_float_type.h"
#include "utils_copy_data.h"

////needed for norm calculation
#include "utils_reduce_vec.h"

namespace qlat{

template <typename Ty >
struct vector_gpu{
  ///Vector<Ty > v;
  Ty*    p;
  size_t n;
  bool GPU;

  vector_gpu<Ty >()
  {
    p = NULL; n = 0;GPU = true;
  }

  vector_gpu<Ty >(const size_t n_set, const bool GPU_set = true)
  {
    p = NULL; n = 0;GPU = true;
    resize(n_set, GPU_set);
  }

  void resize(const size_t n_set, const bool GPU_set)
  {
    if(n_set == 0){clean_mem(); return ;}
    if((n != n_set) or (GPU != GPU_set))
    {
      /////print0("Reisze n %d, n_set %d !\n", int(n), int(n_set));
      clean_mem();
      n = n_set;
      GPU    = GPU_set;
      if(GPU){
        #ifdef QLAT_USE_ACC
        gpuErrchk(cudaMalloc(&p, n*sizeof(Ty)));
        #else
        p = (Ty*) aligned_alloc_no_acc( n*sizeof(Ty));
        #endif
      }
      else{p = (Ty*) aligned_alloc_no_acc(n*sizeof(Ty));}

      set_zero();
    }
    ////qacc_barrier(dummy);
  }

  void resize(const size_t n_set)
  {
    bool GPU_tem = GPU;
    resize(n_set, GPU_tem);
  }

  inline size_t size() const{return n;}
  inline Ty* data(){return p;}
  inline const Ty* data() const { return p; }

  inline const Ty& operator[](const size_t i) const { return p[i]; }
  inline Ty& operator[](const size_t i) { return p[i]; }

  void set_zero(bool dummy = true)
  {
    zero_Ty(p, n, GPU, dummy);
  }

  void clean_mem(){
    free_buf(p, GPU);
    p = NULL;n = 0;
  }
  void clear(){clean_mem();}

  ~vector_gpu(){
    clean_mem();
  }

  template <class T >
  vector_gpu<Ty>& operator=(const vector_gpu<T >& vp)
  {
    print0("NO SUPPORT yet!\n");
    qassert(false);
    ////bool tem_GPU = vp.GPU;
    resize(vp.size(), vp.GPU);

    int mode_cpu = 0;
    if(vp.GPU == false and GPU == false){mode_cpu =  0;}
    if(vp.GPU == true  and GPU == true ){mode_cpu =  1;}

    cpy_data_thread(p, vp.p, n, mode_cpu, true);

    return *this;
  }


  ////
  template <class T >
  void copy_from(const T* src, size_t Ndata, int GPU_set = -1, int GPU_ori = 0)
  {
    bool tem_GPU = GPU;
    bool GPU_src = true;
    qassert(GPU_ori == 0 or GPU_ori == 1);
    if(GPU_ori == 0 ){GPU_src = false ;}
    if(GPU_ori == 1 ){GPU_src = true  ;}

    if(GPU_set == -1){tem_GPU = GPU_src;}
    if(GPU_set == 0 ){tem_GPU = false ;}
    if(GPU_set == 1 ){tem_GPU = true  ;}
    resize(Ndata, tem_GPU);
  
    int mode_cpu = 0;
    if(GPU_src == false and GPU == false){mode_cpu =  0;} // host to host
    if(GPU_src == true  and GPU == true ){mode_cpu =  1;} // device to device
    if(GPU_src == false and GPU == true ){mode_cpu =  2;} // host to device
    if(GPU_src == true  and GPU == false){mode_cpu =  3;} // device to host
  
    cpy_data_thread(p, src, n, mode_cpu, true);

    //bool tem_GPU =  GPU;
    //if(GPU_set == -1){tem_GPU = true  ;}
    //if(GPU_set == 0 ){tem_GPU = false ;}
    //if(GPU_set == 1 ){tem_GPU = true  ;}
    //resize(Ndata, tem_GPU);
  
    //int mode_cpu = 0;
    //if(GPU == false){mode_cpu =  3;} // device to host
    //if(GPU == true ){mode_cpu =  1;} // device to device
    //cpy_data_thread(p, &vp[0], n, mode_cpu, true);
  }

  template <class T >
  void copy_from(const vector_gpu<T >& vp, int GPU_set = -1)
  {
    copy_from(vp.p, vp.n, GPU_set, int(vp.GPU));
    //bool tem_GPU =  GPU;
    //if(GPU_set == -1){tem_GPU = vp.GPU;}
    //if(GPU_set == 0 ){tem_GPU = false ;}
    //if(GPU_set == 1 ){tem_GPU = true  ;}
    //resize(vp.size(), tem_GPU);
  
    //int mode_cpu = 0;
    //if(vp.GPU == false and GPU == false){mode_cpu =  0;}
    //if(vp.GPU == true  and GPU == true ){mode_cpu =  1;}
    //if(vp.GPU == false and GPU == true ){mode_cpu =  2;} // host to device
    //if(vp.GPU == true  and GPU == false){mode_cpu =  3;} // device to host
  
    //cpy_data_thread(p, vp.p, n, mode_cpu, true);
  }


  template <class T >
  void copy_from(const std::vector<T >& vp, int GPU_set = -1)
  {
    int GPU_mem = 0;
    if(GPU_set == -1){GPU_mem = 1;}else{GPU_mem = GPU_set;}
    copy_from(&vp[0], vp.size(), GPU_mem, 0);
    //bool tem_GPU =  GPU;
    ////////if(GPU_set == -1){tem_GPU = true  ;}
    //if(GPU_set == 0 ){tem_GPU = false ;}
    //if(GPU_set == 1 ){tem_GPU = true  ;}
    //resize(vp.size(), tem_GPU);
  
    //int mode_cpu = 0;
    //if(GPU == false){mode_cpu =  0;}
    //if(GPU == true ){mode_cpu =  2;} // host to device
    //cpy_data_thread(p, &vp[0], n, mode_cpu, true);
  }


  template <class T >
  void copy_from(qlat::vector_acc<T >& vp, int GPU_set = -1, int GPU_ori = 0)
  {
    int GPU_mem = 0;
    if(GPU_set == -1){GPU_mem = 1;}else{GPU_mem = GPU_set;}
    T* src = (T*) qlat::get_data(vp).data();
    copy_from(src, vp.size(), GPU_mem, GPU_ori);
    //bool tem_GPU =  true;
    //if(GPU_set == -1){tem_GPU = GPU  ;}
    //if(GPU_set == 0 ){tem_GPU = false ;}
    //if(GPU_set == 1 ){tem_GPU = true  ;}
    //resize(vp.size(), tem_GPU);
  
    //int mode_cpu = 0;
    //if(GPU_ori == 0 and GPU == false){mode_cpu =  0;}
    //if(GPU_ori == 1 and GPU == true ){mode_cpu =  1;}
    //if(GPU_ori == 0 and GPU == true ){mode_cpu =  2;} // host to device
    //if(GPU_ori == 1 and GPU == false){mode_cpu =  3;} // device to host
    //T* src = (T*) qlat::get_data(vp).data();
    //cpy_data_thread(p, src, n, mode_cpu, true);
  }

  //////
  //template <class T >
  //void copy_from_GPU(const T* vp, size_t Ndata, int GPU_set = -1)
  //{
  //  bool tem_GPU =  GPU;
  //  if(GPU_set == -1){tem_GPU = true  ;}
  //  if(GPU_set == 0 ){tem_GPU = false ;}
  //  if(GPU_set == 1 ){tem_GPU = true  ;}
  //  resize(Ndata, tem_GPU);
  //
  //  int mode_cpu = 0;
  //  if(GPU == false){mode_cpu =  3;} // device to host
  //  if(GPU == true ){mode_cpu =  1;} // device to device
  //  cpy_data_thread(p, &vp[0], n, mode_cpu, true);
  //}

  //template <class T >
  //void copy_from_CPU(const T* vp, size_t Ndata, int GPU_set = -1)
  //{
  //  bool tem_GPU =  GPU;
  //  if(GPU_set == -1){tem_GPU = true  ;}
  //  if(GPU_set == 0 ){tem_GPU = false ;}
  //  if(GPU_set == 1 ){tem_GPU = true  ;}
  //  resize(Ndata, tem_GPU);
  //
  //  int mode_cpu = 0;
  //  if(GPU == false){mode_cpu =  0;} // host to host
  //  if(GPU == true ){mode_cpu =  2;} // host to device
  //  cpy_data_thread(p, &vp[0], n, mode_cpu, true);
  //}


  template <class T >
  void copy_to(T* res, int GPU_ori = -1)
  {
    int mode_cpu = 0;
    int GPU_set = 0;if(GPU_ori == -1){GPU_set  =  1;}
    if(GPU == false and GPU_set == 0){mode_cpu =  0;} // host to host
    if(GPU == true  and GPU_set == 1){mode_cpu =  1;} // device to device
    if(GPU == true  and GPU_set == 0){mode_cpu =  3;} // device to host
    if(GPU == false and GPU_set == 1){mode_cpu =  2;} // host to device
    cpy_data_thread(res, p, n, mode_cpu, true);
  }


  template <class T >
  void copy_to(vector_gpu<T >& vp, int GPU_set = -1)
  {
    bool tem_GPU =  GPU;
    if(GPU_set == -1){tem_GPU = GPU;}
    if(GPU_set == 0 ){tem_GPU = false ;}
    if(GPU_set == 1 ){tem_GPU = true  ;}
    vp.resize(size(), tem_GPU);
    copy_to(vp.p, GPU_set);
  
    //int mode_cpu = 0;
    //if(GPU == false and vp.GPU == false){mode_cpu =  0;}
    //if(GPU == true  and vp.GPU == true ){mode_cpu =  1;}
    //if(GPU == true  and vp.GPU == false){mode_cpu =  3;} // device to host
    //if(GPU == false and vp.GPU == true ){mode_cpu =  2;} // host to device
    //cpy_data_thread(vp.p, p, n, mode_cpu, true);
  }

  inline Ty norm()
  {
    qlat::vector_acc<Ty > tmp;tmp.resize(1);tmp[0] = 0;
    qlat::vector_gpu<Ty > copy;copy.resize(n, GPU);
    Ty* res = copy.data();Ty* src = p;
    if(GPU){
      qacc_for(isp, long(n),    {res[isp] = qlat::qconj(src[isp]) * src[isp];});
    }
    else{
      qthread_for(isp, long(n), {res[isp] = qlat::qconj(src[isp]) * src[isp];});
    }
    if(GPU == true ){reduce_vec(res, tmp.data(), n, 1);}
    if(GPU == false){reduce_cpu(res, tmp.data(), n, 1);}
    glb_sum(tmp[0]);
    return tmp[0];
  }

  template <class T >
  void copy_to(std::vector<T >& vp)
  {
    vp.resize(n);
    int mode_cpu = 0;
    if(GPU == false){mode_cpu =  0;}
    if(GPU == true ){mode_cpu =  3;} // device to host
    cpy_data_thread(&vp[0], p, n, mode_cpu, true);
  }

  template <class T >
  void copy_to(qlat::vector_acc<T >& vp, int GPU_ori = 0)
  {
    vp.resize(size());
    T* res = (T*) qlat::get_data(vp).data();
    copy_to(res, GPU_ori);
  }

  template <class T >
  const vector_gpu<Ty>& operator+=(const vector_gpu<T >& vp)
  {
    qassert(GPU == vp.GPU and n == vp.n);
    int mode_cpu = 0;
    if(vp.GPU == false and GPU == false){mode_cpu =  0;}
    if(vp.GPU == true  and GPU == true ){mode_cpu =  1;}
    cpy_data_thread(p, vp.p, n, mode_cpu, true, 1.0);
    return *this;
  }

  template <class T >
  const vector_gpu<Ty>& operator-=(const vector_gpu<T >& vp)
  {
    qassert(GPU == vp.GPU and n == vp.n);
    int mode_cpu = 0;
    if(vp.GPU == false and GPU == false){mode_cpu =  0;}
    if(vp.GPU == true  and GPU == true ){mode_cpu =  1;}
    cpy_data_thread(p, vp.p, n, mode_cpu, true, -1.0);
    return *this;
  }

};

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

}


#endif
