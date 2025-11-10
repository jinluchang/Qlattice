// utils_field_gpu.h
// Gen Wang
// May. 2025

#ifndef UTILS_FIELD_GPU_H
#define UTILS_FIELD_GPU_H

#pragma once

#include "utils_vector_GPU.h"

/*
  Field and SelectedField generations
    1. memeory could be pure GPU/CPU/Mananged
    2. Constract from a memeory pointer
    3. Fix `civ` for contraction vectorizations
    4. Dynamical group of vectors Int boxes
    5. Different memeory order support
    6. Unified interface for Field and SelectedField
*/

namespace qlat
{

/* 
  View the field to field_gpu
  field is always is_copy
*/
template <class M>
struct API FieldG : Field<M> {
  // vector<M > field; need to be copy
  vector_gpu<M> field_gpu;
  QMEM_ORDER mem_order;

  void set_ghost_field_zero(){
    // resize to 0 if not a copy
    if(!Field<M>::field.is_copy){
      Field<M>::field.resize(0);
    }
    Field<M>::field.is_copy = true;
    Field<M>::field.v.p = NULL;
    Field<M>::field.v.n = 0;
    Field<M>::initialized = false;
  }

  void set_ghost_field_pointer(){
    set_ghost_field_zero();
    Field<M>::field.is_copy = true;
    Field<M>::field.v.p = field_gpu.p;
    Field<M>::field.v.n = field_gpu.n;
    Field<M>::initialized = true;
    const Long Ve = Field<M>::geo().local_volume_expanded();
    Qassert(field_gpu.n % Ve == 0);
    Field<M>::multiplicity = field_gpu.n / Ve;
  }

  // set the fields from qacc fields and qacc into is_copy
  void set_ghost_field(QMEM GPU, bool copy_data = false){
    Qassert(Field<M>::field.size() != 0);
    const Long n = Field<M>::field.size();
    field_gpu.resize(n, GPU);
    if(copy_data){
      cpy_GPU(field_gpu.p, Field<M>::field.v.p, n, GPU, 1);
    }
    set_ghost_field_pointer();
  }

  Int get_memtype(){
    return field_gpu.GPU;
  }

  Long get_nsites(){
    return Field<M>::geo().local_volume_expanded();
  }

  Long get_nelems(){
    return Field<M>::geo().local_volume_expanded() * Field<M>::multiplicity;
  }

  // resize all to zero
  // reset field to default not copy one
  void clear_copy()
  {
    if(Field<M>::initialized){
      field_gpu.clear_copy();
      Field<M>::field.is_copy = false;
      Field<M>::field.v.p     = NULL;
      Field<M>::field.v.n     = 0;
      Field<M>::initialized   = false;
    }
  }

  void init() { 
    clear_copy();
    Field<M>::init(); 
    set_ghost_field_zero();
    mem_order = QLAT_DEFAULT;
  }

  // QLAT_OUTTER
  void init(const Geometry& geo_, const Int multiplicity_, QMEM GPU = QMGPU, QMEM_ORDER mem_order_ = QLAT_DEFAULT)
  {
    clear_copy();
    Field<M>::init(geo_, multiplicity_);
    set_ghost_field(GPU);
    mem_order = mem_order_;
  }

  // initialize only not initilized
  void init_zero(const Geometry& geo_, const Int multiplicity_, QMEM GPU = QMGPU, QMEM_ORDER mem_order_ = QLAT_DEFAULT)
  {
    bool need_init = false;
    if(Field<M>::initialized){
      Qassert(field_gpu.GPU == GPU);
      if(Field<M>::geo() != geo_ or Field<M>::multiplicity != multiplicity_ or field_gpu.GPU != GPU or mem_order_ != mem_order){
        need_init = true;
      }else{
        // check parameters
        Field<M>::init_zero(geo_, multiplicity_);
      }
    }else{
      need_init = true;
    }
    if(need_init){
      init(geo_, multiplicity_, QMGPU, mem_order_);
    }
  }

  // only structures
  void init_size(const FieldG<M>& f)
  {
    if(!Field<M>::initialized or Field<M>::geo() != f.geo() or Field<M>::multiplicity != f.multiplicity or field_gpu.GPU != f.field_gpu.GPU or mem_order != f.mem_order)
    {
      QMEM Gtype = get_type_mem(f.field_gpu.GPU);
      init(f.geo(), f.multiplicity, Gtype, f.mem_order);
    }
  }

  // copy the content
  void init(const Field<M>& f, QMEM GPU = QMGPU, QMEM_ORDER mem_order_ = QLAT_DEFAULT)
  {
    clear_copy();
    Field<M>::init(f);
    set_ghost_field(GPU, true);
    mem_order = mem_order_;
  }

  // construct from vector_gpu
  void set_pointer(M* srcp, const Long Nd, const Geometry& geo_, const QMEM GPU = QMGPU,
    const QMEM_ORDER order_ = QLAT_DEFAULT)
  {
    TIMER("FieldG set_pointer");
    Qassert(Nd % geo_.local_volume_expanded() == 0);
    // clear current field
    clear_copy();
    //Field<M>::geo.set(geo_);
    Field<M>::geo.set_view(geo_);
    field_gpu.p = srcp;
    field_gpu.n = Nd;
    field_gpu.GPU = GPU;
    field_gpu.is_copy = true;
    mem_order = order_;
    set_ghost_field_pointer();
  }

  void set_pointer(FieldG<M>& src, const Long size = 0, const Long offset = 0)
  {
    TIMER("FieldG set_pointer");
    const Geometry& geo = src.geo();
    Long Nd = src.field_gpu.n;
    if(size > 0){Nd = size;}
    Qassert(offset + Nd <= Long(src.field_gpu.n));
    clear_copy();
    //Field<M>::geo.set(geo);
    Field<M>::geo.set_view(geo);
    field_gpu.p = &src.field_gpu.p[offset];
    field_gpu.n = Nd;
    field_gpu.GPU = src.field_gpu.GPU;
    field_gpu.is_copy = true;
    mem_order = src.mem_order;
    set_ghost_field_pointer();
  }

  void set_pointer(Field<M>& src, QMEM_ORDER mem_order_, const Long size = 0, const Long offset = 0)
  {
    TIMER("FieldG set_pointer");
    const Geometry& geo = src.geo();
    M* p = (M*) get_data(src).data();
    Long Nd = src.field.size();
    if(size > 0){Nd = size;}
    Qassert(offset + Nd <= src.field.size());
    clear_copy();
    //Field<M>::geo.set(geo);
    Field<M>::geo.set_view(geo);
    field_gpu.p = &p[offset];
    field_gpu.n = Nd;
    field_gpu.GPU = QMGPU;
    field_gpu.is_copy = true;
    mem_order = mem_order_;
    set_ghost_field_pointer();
  }

  // construct from vector_gpu
  void set_pointer(const vector_gpu<M>& src, const Geometry& geo_, 
    const QMEM_ORDER order_ = QLAT_DEFAULT)
  {
    set_pointer(src.p, src.n, src.GPU, geo_, order_);
  }

  void set_zero()
  {
    field_gpu.set_zero();
  }

  //
  FieldG() { init(); }
  FieldG(const FieldG<M>&) = default;
  FieldG(FieldG<M>&&) noexcept = default;
  FieldG<M>& operator=(FieldG<M>&&) noexcept = default;
  //
};

/* 
  View the field to field_gpu
  field is always is_copy
*/
template <class M>
struct API SelectedFieldG : SelectedField<M> {
  // vector<M > field; need to be copy
  vector_gpu<M> field_gpu;
  QMEM_ORDER mem_order;

  void set_ghost_field_zero(){
    // resize to 0 if not a copy
    if(!SelectedField<M>::field.is_copy){
      SelectedField<M>::field.resize(0);
    }
    SelectedField<M>::field.is_copy = true;
    SelectedField<M>::field.v.p = NULL;
    SelectedField<M>::field.v.n = 0;
    SelectedField<M>::initialized  = false;
  }

  void set_ghost_field_pointer(){
    set_ghost_field_zero();
    SelectedField<M>::field.is_copy = true;
    SelectedField<M>::field.v.p = field_gpu.p;
    SelectedField<M>::field.v.n = field_gpu.n;
    SelectedField<M>::initialized = true;
    const Long Ve = SelectedField<M>::n_elems;
    Qassert(field_gpu.n % Ve == 0);
    SelectedField<M>::multiplicity = field_gpu.n / Ve;
  }

  // set the fields from qacc fields and qacc into is_copy
  void set_ghost_field(QMEM GPU, bool copy_data = false){
    Qassert(SelectedField<M>::field.size() != 0);
    const Long n = SelectedField<M>::field.size();
    field_gpu.resize(n, GPU);
    if(copy_data){
      cpy_GPU(field_gpu.p, SelectedField<M>::field.v.p, n, GPU, 1);
    }
    set_ghost_field_pointer();
  }

  Long get_nsites(){
    return SelectedField<M>::n_elems;
  }

  Long get_nelems(){
    return SelectedField<M>::n_elems * SelectedField<M>::multiplicity;
  }

  // resize all to zero
  // reset field to default not copy one
  void clear_copy()
  {
    if(SelectedField<M>::initialized){
      field_gpu.clear_copy();
      SelectedField<M>::field.is_copy = false;
      SelectedField<M>::field.v.p     = NULL;
      SelectedField<M>::field.v.n     = 0;
      SelectedField<M>::initialized   = false;
      SelectedField<M>::n_elems       = 0;
      SelectedField<M>::multiplicity  = 0;
    }
  }

  void init() {
    clear_copy();
    SelectedField<M>::init(); 
    set_ghost_field_zero();
    mem_order = QLAT_DEFAULT;
  }

  void init(const Geometry& geo_, const Long n_elems_, const Int multiplicity_, QMEM GPU = QMGPU, QMEM_ORDER mem_order_ = QLAT_DEFAULT)
  {
    clear_copy();
    SelectedField<M>::init(geo_, n_elems_, multiplicity_);
    set_ghost_field(GPU);
    mem_order = mem_order_;
  }

  void init(const FieldSelection& fsel, const Int multiplicity_, QMEM GPU = QMGPU, QMEM_ORDER mem_order_ = QLAT_DEFAULT)
  {
    clear_copy();
    SelectedField<M>::init(fsel, multiplicity_);
    set_ghost_field(GPU);
    mem_order = mem_order_;
  }

  // initialize only not initilized, could change parameters compare to initial ones
  void init_zero(const Geometry& geo_, const Long n_elems_, const Int multiplicity_, QMEM GPU = QMGPU, QMEM_ORDER mem_order_ = QLAT_DEFAULT)
  {
    bool need_init = false;
    if(SelectedField<M>::initialized){
      Qassert(field_gpu.GPU == GPU);
      if(SelectedField<M>::geo() != geo_ or SelectedField<M>::n_elems != n_elems_ or SelectedField<M>::multiplicity != multiplicity_ or field_gpu.GPU != GPU or mem_order_ != mem_order){
        need_init = true;
      }else{ 
        // check parameters
        SelectedField<M>::init_zero(geo_, n_elems_, multiplicity_, mem_order_);
      }
    }else{
      need_init = true;
    }
    if(need_init){
      init(geo_, n_elems_, multiplicity_, QMGPU, mem_order_);
    }
  }

  void init_zero(const FieldSelection& fsel, const Int multiplicity_, QMEM GPU = QMGPU)
  {
    if(SelectedField<M>::initialized){
      SelectedField<M>::init_zero(fsel, multiplicity_);
      Qassert(field_gpu.GPU == GPU);
    }else{
      init(fsel, multiplicity_, QMGPU);
    }
  }

  // construct from vector_gpu
  void set_pointer(M* srcp, const Long n_elems_, const Int multiplicity_ , const Geometry& geo_, const QMEM GPU = QMGPU,
    const QMEM_ORDER order_ = QLAT_DEFAULT)
  {
    TIMER("SelectedFieldG set_pointer");
    // clear current field
    const Long Nd = n_elems_ * multiplicity_;
    clear_copy();
    SelectedField<M>::geo.set(geo_);
    SelectedField<M>::n_elems = n_elems_;
    field_gpu.p = srcp;
    field_gpu.n = Nd;
    field_gpu.GPU = GPU;
    field_gpu.is_copy = true;
    mem_order = order_;
    set_ghost_field_pointer();
  }

  // construct from vector_gpu
  void set_pointer(const vector_gpu<M>& src, const Long n_elems_, const Geometry& geo_, 
    const QMEM_ORDER order_ = QLAT_DEFAULT)
  {
    Qassert(src.n % n_elems_ == 0);
    const Long multiplicity_ = src.n / n_elems_;
    set_pointer(src.p, n_elems_, multiplicity_, src.GPU, geo_, order_);
  }

  void set_zero()
  {
    field_gpu.set_zero();
  }

  //
  SelectedFieldG() { init(); }
  SelectedFieldG(const SelectedFieldG<M>&) = default;
  SelectedFieldG(SelectedFieldG<M>&&) noexcept = default;
  SelectedFieldG<M>& operator=(SelectedFieldG<M>&&) noexcept = default;
  //
};


template <typename M >
inline void set_zero(FieldG<M>& f)
{
  f.set_zero();
}

template <typename M >
inline void set_zero(SelectedFieldG<M>& f)
{
  f.set_zero();
}

template <typename M >
inline QMEM_ORDER get_mem_order(FieldG<M>& f)
{
  return f.mem_order;
}

template <typename M >
inline QMEM_ORDER get_mem_order(SelectedFieldG<M>& f)
{
  return f.mem_order;
}

template <typename M >
inline QMEM_ORDER get_mem_order(Field<M>& f)
{
  Qassert(f.initialized);
  return QLAT_DEFAULT;
}

template <class M>
inline void set_field(Field<M >& res, M* src, const Long n, const Geometry& geo,
  const MemType& type   = MemType::Acc, 
  const MemOrder& order =MemOrder::TZYXM)
{
  if(res.initialized){
    res.init();
  }
  res.initialized = true;
  //res.geo.set_view(geo);

  res.geo.set(geo);
  const Long Nd = geo.local_volume_expanded();
  Qassert(n % Nd == 0);
  const Long inner = n / Nd;

  res.multiplicity = inner;
  res.mem_order    = order;
  
  //Vector<M> vec;
  //vec.p = src;
  //vec.n =   n;

  //vector<M> field;
  //field.set_view(vec);
  //field.mem_type = type;

  res.field.is_copy  = true;
  res.field.mem_type = type;
  res.field.v.p = src;
  res.field.v.n = n;

  //res.field.set_view(field);
}

}

#endif
