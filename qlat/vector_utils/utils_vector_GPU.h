// utils_vector_GPU.h
// Gen Wang
// Sep. 2021

#ifndef UTILS_VECTOR_GPU_H
#define UTILS_VECTOR_GPU_H

#pragma once
#include <qlat/qcd.h>
#include "float_type.h"
#include "utils_copy_data.h"

namespace qlat{

template <typename Ty >
struct vector_gpu
{
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
      else{p = (Ty*) alloc_mem_alloc_no_acc(n*sizeof(Ty));}
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
  const vector_gpu<Ty>& operator=(const vector_gpu<T >& vp)
  {
    ////bool tem_GPU = vp.GPU;
    resize(vp.size(), vp.GPU);

    int mode_cpu = 0;
    if(vp.GPU == false and GPU == false){mode_cpu =  0;}
    if(vp.GPU == true  and GPU == true ){mode_cpu =  1;}

    cpy_data_thread(p, vp.p, n, mode_cpu, true);

    return *this;
  }

  template <class T >
  void copy_from(vector_gpu<T >& vp, int GPU_set = -1)
  {
    bool tem_GPU =  GPU;
    if(GPU_set == -1){tem_GPU = vp.GPU;}
    if(GPU_set == 0 ){tem_GPU = false ;}
    if(GPU_set == 1 ){tem_GPU = true  ;}
    resize(vp.size(), tem_GPU);
  
    int mode_cpu = 0;
    if(vp.GPU == false and GPU == false){mode_cpu =  0;}
    if(vp.GPU == true  and GPU == true ){mode_cpu =  1;}
    if(vp.GPU == false and GPU == true ){mode_cpu =  2;} // host to device
    if(vp.GPU == true  and GPU == false){mode_cpu =  3;} // device to host
  
    cpy_data_thread(p, vp.p, n, mode_cpu, true);
  }

  template <class T >
  void copy_to(vector_gpu<T >& vp, int GPU_set = -1)
  {
    bool tem_GPU =  GPU;
    if(GPU_set == -1){tem_GPU = GPU;}
    if(GPU_set == 0 ){tem_GPU = false ;}
    if(GPU_set == 1 ){tem_GPU = true  ;}
    vp.resize(size(), tem_GPU);
  
    int mode_cpu = 0;
    if(GPU == false and vp.GPU == false){mode_cpu =  0;}
    if(GPU == true  and vp.GPU == true ){mode_cpu =  1;}
    if(GPU == true  and vp.GPU == false){mode_cpu =  3;} // device to host
    if(GPU == false and vp.GPU == true ){mode_cpu =  2;} // host to device
  
    cpy_data_thread(vp.p, p, n, mode_cpu, true);
  }


  template <class T >
  void copy_from(const std::vector<T >& vp, int GPU_set = -1)
  {
    bool tem_GPU =  GPU;
    //////if(GPU_set == -1){tem_GPU = true  ;}
    if(GPU_set == 0 ){tem_GPU = false ;}
    if(GPU_set == 1 ){tem_GPU = true  ;}
    resize(vp.size(), tem_GPU);
  
    int mode_cpu = 0;
    if(GPU == false){mode_cpu =  0;}
    if(GPU == true ){mode_cpu =  2;} // host to device
    cpy_data_thread(p, &vp[0], n, mode_cpu, true);
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
  void copy_from(qlat::vector_acc<T >& vp, int GPU_set = -1, int GPU_ori = 0)
  {
    bool tem_GPU =  true;
    if(GPU_set == -1){tem_GPU = GPU  ;}
    if(GPU_set == 0 ){tem_GPU = false ;}
    if(GPU_set == 1 ){tem_GPU = true  ;}
    resize(vp.size(), tem_GPU);
  
    int mode_cpu = 0;
    if(GPU_ori == 0 and GPU == false){mode_cpu =  0;}
    if(GPU_ori == 1 and GPU == true ){mode_cpu =  1;}
    if(GPU_ori == 0 and GPU == true ){mode_cpu =  2;} // host to device
    if(GPU_ori == 1 and GPU == false){mode_cpu =  3;} // device to host
    T* src = (T*) qlat::get_data(vp).data();
    cpy_data_thread(p, src, n, mode_cpu, true);
  }

  template <class T >
  void copy_to(qlat::vector_acc<T >& vp, int GPU_ori = 0)
  {
    vp.resize(size());
  
    int mode_cpu = 0;
    if(GPU == false and GPU_ori == 0){mode_cpu =  0;}
    if(GPU == true  and GPU_ori == 1){mode_cpu =  1;}
    if(GPU == true  and GPU_ori == 0){mode_cpu =  3;} // device to host
    if(GPU == false and GPU_ori == 1){mode_cpu =  2;} // host to device
    T* res = (T*) qlat::get_data(vp).data();
  
    cpy_data_thread(res, p, n, mode_cpu, true);
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


}


#endif
