#pragma once

#include <qlat-utils/config.h>
#include <qlat-utils/handle.h>
#include <qlat-utils/qacc-func.h>
#include <qlat-utils/utils.h>
#include <unistd.h>

#include <cstdlib>
#include <cstring>
#include <unordered_map>

#ifdef QLAT_USE_MALLOC_STATS
#include <malloc.h>
#endif

namespace qlat
{  //

void clear_mem_cache();

void clear_all_caches();

// --------------------

API inline Long& get_alignment(const MemType mem_type)
// qlat parameter
//
// Should NOT change in the middle of the run.
//
// mem_type identified as same type by `is_same_mem_type` should have the same
// chunk_size.
{
#ifdef QLAT_USE_ACC
  static Long v = get_env_long_default("q_alignment", 256);
  static Long v_cpu = get_env_long_default("q_alignment_cpu", v);
  static Long v_acc = get_env_long_default("q_alignment_acc", v);
  static Long v_uvm = get_env_long_default("q_alignment_uvm", v);
  static Long v_comm = get_env_long_default("q_alignment_comm", v);
  static Long v_comm_acc = get_env_long_default("q_alignment_comm_acc", v);
  if (mem_type == MemType::Cpu) {
    return v_cpu;
  } else if (mem_type == MemType::Acc) {
    return v_acc;
  } else if (mem_type == MemType::Uvm) {
    return v_uvm;
  } else if (mem_type == MemType::Comm) {
    return v_comm;
  } else if (mem_type == MemType::CommAcc) {
    return v_comm_acc;
  } else {
    qassert(false);
    return v;
  }
#else
  static Long v = get_env_long_default("q_alignment", 256);
  static Long v_comm = get_env_long_default("q_alignment_comm", v);
  if (mem_type == MemType::Cpu) {
    return v;
  } else if (mem_type == MemType::Acc) {
    return v;
  } else if (mem_type == MemType::Uvm) {
    return v;
  } else if (mem_type == MemType::Comm) {
    return v_comm;
  } else if (mem_type == MemType::CommAcc) {
    return v_comm;
  } else {
    qassert(false);
    return v;
  }
#endif
}

API inline Long& get_mem_chunk_size(const MemType mem_type)
// qlat parameter
//
// Should NOT change in the middle of the run.
//
// mem_type identified as same type by `is_same_mem_type` should have the same
// chunk_size.
{
#ifdef QLAT_USE_ACC
  static Long v = get_env_long_default("q_mem_chunk_size", 512);
  static Long v_cpu = get_env_long_default("q_mem_chunk_size_cpu", v);
  static Long v_acc = get_env_long_default("q_mem_chunk_size_acc", 1);
  static Long v_uvm = get_env_long_default("q_mem_chunk_size_uvm", v);
  static Long v_comm = get_env_long_default("q_mem_chunk_size_comm", 1);
  static Long v_comm_acc = get_env_long_default("q_mem_chunk_size_comm_acc", 1);
  if (mem_type == MemType::Cpu) {
    return v_cpu;
  } else if (mem_type == MemType::Acc) {
    return v_acc;
  } else if (mem_type == MemType::Uvm) {
    return v_uvm;
  } else if (mem_type == MemType::Comm) {
    return v_comm;
  } else if (mem_type == MemType::CommAcc) {
    return v_comm_acc;
  } else {
    qassert(false);
    return v;
  }
#else
  static Long v = get_env_long_default("q_mem_chunk_size", 512);
  static Long v_comm = get_env_long_default("q_mem_chunk_size_comm", 1);
  if (mem_type == MemType::Cpu) {
    return v;
  } else if (mem_type == MemType::Acc) {
    return v;
  } else if (mem_type == MemType::Uvm) {
    return v;
  } else if (mem_type == MemType::Comm) {
    return v_comm;
  } else if (mem_type == MemType::CommAcc) {
    return v_comm;
  } else {
    qassert(false);
    return v;
  }
#endif
}

API inline bool& get_mem_type_comm_use_acc()
// qlat parameter
//
// Should NOT change in the middle of the run.
{
  static bool v = get_env_default("q_mem_type_comm_use_acc", "true") == "true";
  return v;
}

inline bool is_same_mem_type(const MemType t1, const MemType t2)
{
  if (t1 == t2) {
    return true;
  } else {
#ifdef QLAT_USE_ACC
    if (get_mem_type_comm_use_acc()) {
      return false;
    } else if (t1 == MemType::Comm) {
      return t2 == MemType::CommAcc;
    } else if (t1 == MemType::CommAcc) {
      return t2 == MemType::Comm;
    } else {
      return false;
    }
#else
    if (t1 == MemType::Comm or t1 == MemType::CommAcc) {
      return (t2 == MemType::Comm or t2 == MemType::CommAcc);
    } else if (t2 == MemType::Comm or t2 == MemType::CommAcc) {
      return (t1 == MemType::Comm or t1 == MemType::CommAcc);
    } else {
      return true;
    }
#endif
  }
}

inline MemType get_eff_mem_type(const MemType mem_type)
// Return only Cpu, Acc, Uvm
{
  if (mem_type == MemType::Comm) {
    return MemType::Cpu;
  } else if (mem_type == MemType::CommAcc) {
    if (get_mem_type_comm_use_acc()) {
      return MemType::Acc;
    } else {
      return MemType::Cpu;
    }
  }
  return mem_type;
}

// assume always have Acc memeory to compare even on CPU
inline MemType get_comm_mem_type(const MemType mem_type)
{
  if(!get_mem_type_comm_use_acc()){return MemType::Cpu;}
  else{
    MemType b = get_eff_mem_type(mem_type);
    if(b == MemType::Cpu){return MemType::Cpu;}
    else{return MemType::Acc;}
  }
}

inline bool is_same_comm_mem_type(const MemType t1, const MemType t2)
{
  if(get_comm_mem_type(t1) == get_comm_mem_type(t2)){
    return true;
  }
  return false;
}

inline Long get_chunked_mem_size(const Long chunk_size, const Long min_size)
{
  const Long n_elem = 1 + (min_size - 1) / chunk_size;
  const Long size = n_elem * chunk_size;
  return size;
}

API inline Long& get_mem_cache_max_size(const MemType mem_type)
// qlat parameter
//
// environment variable in unit in MB
{
  static Long v = get_env_long_default("q_mem_cache_max_size", 512);
  static Long v_cpu =
      get_env_long_default("q_mem_cache_max_size_cpu", v) * 1024L * 1024L;
  static Long v_acc =
      get_env_long_default("q_mem_cache_max_size_acc", v) * 1024L * 1024L;
  static Long v_uvm =
      get_env_long_default("q_mem_cache_max_size_uvm", v) * 1024L * 1024L;
  static Long v_comm =
      get_env_long_default("q_mem_cache_max_size_comm", v) * 1024L * 1024L;
  static Long v_comm_acc =
      get_env_long_default("q_mem_cache_max_size_comm_acc", v) * 1024L * 1024L;
  if (mem_type == MemType::Cpu) {
    return v_cpu;
  } else if (mem_type == MemType::Acc) {
    return v_acc;
  } else if (mem_type == MemType::Uvm) {
    return v_uvm;
  } else if (mem_type == MemType::Comm) {
    return v_comm;
  } else if (mem_type == MemType::CommAcc) {
    return v_comm_acc;
  } else {
    qassert(false);
    return v;
  }
}

API inline Long& get_alloc_mem_max_size()
// qlat parameter
// unit in MB
{
  static Long max_size =
      get_env_long_default("q_alloc_mem_max_size", 256L * 1024L) * 1024L *
      1024L;
  return max_size;
}

struct MemoryStats {
  Long alloc[static_cast<Int>(MemType::SIZE)];
  Long cache[static_cast<Int>(MemType::SIZE)];
  //
  MemoryStats() { init(); }
  //
  void init()
  {
    for (Int i = 0; i < static_cast<Int>(MemType::SIZE); ++i) {
      alloc[i] = 0;
      cache[i] = 0;
    }
  }
  //
  Long total()
  // total memory include the size of the cache
  {
    Long total = 0;
    for (Int i = 0; i < static_cast<Int>(MemType::SIZE); ++i) {
      total += alloc[i] + cache[i];
    }
    return total;
  }
};

API inline MemoryStats& get_mem_stats()
{
  static MemoryStats ms;
  return ms;
}

void* alloc_mem_alloc(const Long size, const MemType mem_type);

void free_mem_free(void* ptr, const Long size, const MemType mem_type);

struct API MemCache {
  MemType mem_type;
  Long mem_cache_max_size;
  Long mem_cache_size;
  std::unordered_multimap<Long, void*> db;
  //
  MemCache(const MemType mem_type_ = get_default_mem_type(),
           const Long mem_cache_max_size_ = -1)
  {
    init(mem_type_, mem_cache_max_size_);
  }
  ~MemCache() { gc(); }
  //
  void init(const MemType mem_type_ = get_default_mem_type(),
            const Long mem_cache_max_size_ = -1);
  void add(void* ptr, const Long size);
  void* del(const Long size);
  void gc();
};

std::vector<MemCache> mk_mem_cache_vec();

API inline MemCache& get_mem_cache(
    const MemType mem_type = get_default_mem_type())
{
  static std::vector<MemCache> cache_vec = mk_mem_cache_vec();
  return cache_vec[static_cast<Int>(mem_type)];
}

void* alloc_mem(const Long min_size, const MemType mem_type);

void free_mem(void* ptr, const Long min_size, const MemType mem_type);

void copy_mem(void* dst, const MemType mem_type_dst, const void* src,
              const MemType mem_type_src, const Long size);

void set_mem(void* ptr, const Int v, const Long size, const MemType mem_type);

// --------------------

inline void displayln_malloc_stats()
{
#ifdef QLAT_USE_MALLOC_STATS
  malloc_stats();
#endif
}

// --------------------

template <class M>
struct API vector {
  // Avoid copy constructor when possible
  // (it is likely not what you think it is)
  // Only used in qacc macros, or if it is already a copy.
  //
  Vector<M> v;
  mutable MemType mem_type;  // if place data on accelerator memory
  bool is_copy;              // do not free memory if is_copy=true
  //
  vector(const MemType mem_type_ = get_default_mem_type())
  {
    qassert(v.p == NULL);
    mem_type = mem_type_;
    is_copy = false;
  }
  vector(const Long size, const MemType mem_type_ = get_default_mem_type())
  {
    // TIMER("vector::vector(size)")
    qassert(v.p == NULL);
    mem_type = mem_type_;
    is_copy = false;
    resize(size);
  }
  vector(const std::vector<M>& vp,
         const MemType mem_type_ = get_default_mem_type())
  {
    // TIMER("vector::vector(std::vector&)")
    qassert(v.p == NULL);
    mem_type = mem_type_;
    is_copy = false;
    *this = vp;
  }
  vector(const vector<M>& vp)
  {
    // TIMER("vector::vector(&)")
#ifndef QLAT_USE_ACC
    qassert(vp.is_copy);
#endif
    is_copy = true;
    const MemType eff_mem_type = get_eff_mem_type(vp.mem_type);
    qassert(eff_mem_type == MemType::Uvm or eff_mem_type == MemType::Acc);
    mem_type = vp.mem_type;
    v = vp.v;
  }
  vector(vector<M>&& vp) noexcept
  {
    // TIMER("vector::vector(&&)")
    is_copy = vp.is_copy;
    mem_type = vp.mem_type;
    v = vp.v;
    vp.is_copy = true;
  }
  //
  ~vector()
  {
    // TIMER("vector::~vector()")
    if (not is_copy) {
      clear();
    }
  }
  //
  void init()
  // does not change mem_type
  {
    if (not is_copy) {
      clear();
    } else {
      v = Vector<M>();
    }
    is_copy = false;
  }
  //
  void clear()
  // does not change mem_type
  {
    qassert(not is_copy);
    if (v.p != NULL) {
      free_mem(v.p, v.n * sizeof(M), mem_type);
    }
    v = Vector<M>();
    qassert(v.p == NULL);
  }
  //
  void set_mem_type(const MemType mem_type_) const
  {
    qassert(not is_copy);
    if (NULL == v.p) {
      mem_type = mem_type_;
      return;
    }
    if (is_same_mem_type(mem_type, mem_type_)) {
      mem_type = mem_type_;
      return;
    }
    {
      TIMER("vector::set_mem_type");
      vector<M> vec(mem_type_);
      vec = *this;
      std::swap(mem_type, vec.mem_type);
      std::swap(v.p, vec.v.p);
      qassert(mem_type == mem_type_);
    }
  }
  MemType get_mem_type() const
  {
    return mem_type;
  }
  //
  void set_view(const Vector<M>& vec) { set_view(vec, mem_type); }
  void set_view(const Vector<M>& vec, const MemType mem_type_)
  // does not change mem_type
  {
    init();
    set_mem_type(mem_type_);
    is_copy = true;
    v = vec;
  }
  void set_view(const vector<M>& vec)
  {
    init();
    is_copy = true;
    mem_type = vec.mem_type;
    v = vec.v;
  }
  //
  template <class N>
  void set_view_cast(const Vector<N>& vec)
  // does not change mem_type
  {
    init();
    is_copy = true;
    v.set_cast(vec.v);
  }
  template <class N>
  void set_view_cast(const vector<N>& vec)
  {
    init();
    is_copy = true;
    mem_type = vec.mem_type;
    v.set_cast(vec.v);
  }
  //
  void resize(const Long size) { resize(size, mem_type); }
  void resize(const Long size, const MemType mem_type_)
  {
    qassert(not is_copy);
    qassert(0 <= size);
    if (v.p == NULL) {
      v.p = (M*)alloc_mem(size * sizeof(M), mem_type_);
      v.n = size;
      mem_type = mem_type_;
      return;
    }
    if (size == v.n) {
      if (mem_type != mem_type_) {
        set_mem_type(mem_type_);
      }
      return;
    }
    if (size == 0) {
      clear();
      mem_type = mem_type_;
      return;
    }
    {
      TIMER("vector::resize");
      vector<M> vp(size, mem_type_);
      qswap(*this, vp);
      const Long n_min = std::min(size, vp.v.n);
      copy_mem(v.p, mem_type, vp.v.p, vp.mem_type, n_min * sizeof(M));
    }
  }
  //
  void resize_zero(const Long size){ resize_zero(size, mem_type); }
  void resize_zero(const Long size, const MemType mem_type_)
  {
    resize(0);
    if(size == 0){return ;}
    resize(size, mem_type_);
    set_zero(*this);
  }
  //
  vector<M>& operator=(const vector<M>& vp)
  // does not change mem_type
  {
    // TIMER("vector::operator=(&)");
    qassert(not is_copy);
    clear();
    resize(vp.size());
    copy_mem(v.p, mem_type, vp.v.p, vp.mem_type, v.n * sizeof(M));
    return *this;
  }
  vector<M>& operator=(vector<M>&& vp) noexcept
  {
    // TIMER("vector::operator=(&&)");
    is_copy = vp.is_copy;
    mem_type = vp.mem_type;
    v = vp.v;
    vp.is_copy = true;
    return *this;
  }
  vector<M>& operator=(const std::vector<M>& vp)
  {
    // TIMER("vector::operator=(std::vector&)");
    qassert(not is_copy);
    clear();
    resize(vp.size());
    copy_mem(v.p, mem_type, vp.data(), MemType::Cpu, v.n * sizeof(M));
    return *this;
  }
  //
  qacc const M& operator[](const Long i) const { return v[i]; }
  qacc M& operator[](const Long i) { return v[i]; }
  //
  qacc Long size() const { return v.size(); }
  //
  qacc M* data() { return v.data(); }
  qacc const M* data() const { return v.data(); }
};

template <class M>
void clear(vector<M>& v)
{
  v.clear();
}

template <class M>
qacc void qswap(vector<M>& v1, vector<M>& v2)
{
  qassert(not v1.is_copy);
  qassert(not v2.is_copy);
  std::swap(v1.mem_type, v2.mem_type);
  std::swap(v1.v, v2.v);
}

template <class M, class N>
qacc void qswap_cast(vector<M>& v1, vector<N>& v2)
{
  qassert(not v1.is_copy);
  qassert(not v2.is_copy);
  const Long data_size1 = v2.v.data_size();
  const Long data_size2 = v1.v.data_size();
  const Long size1 = data_size1 / sizeof(M);
  const Long size2 = data_size2 / sizeof(N);
  qassert(size1 * (Long)sizeof(M) == data_size1);
  qassert(size2 * (Long)sizeof(N) == data_size2);
  void* p_tmp = (void*)v1.v.p;
  v1.v.p = (M*)v2.v.p;
  v2.v.p = (N*)p_tmp;
  v1.v.n = size1;
  v2.v.n = size2;
  std::swap(v1.mem_type, v2.mem_type);
}

// --------------------

template <class M>
qacc bool operator==(const vector<M>& v1, const vector<M>& v2)
{
  if (v1.size() != v2.size()) {
    return false;
  }
  const Int cmp = std::memcmp(v1.data(), v2.data(), v1.size() * sizeof(M));
  return cmp == 0;
}

template <class M>
qacc bool operator!=(const vector<M>& v1, const vector<M>& v2)
{
  return not(v1 == v2);
}

// --------------------

template <class M>
struct IsDataVectorType<vector<M>> {
  using DataType = M;
  using BasicDataType = typename IsDataValueType<DataType>::BasicDataType;
  using ElementaryType = typename IsDataValueType<DataType>::ElementaryType;
  static constexpr bool value = is_data_value_type<DataType>();
};

template <class M>
qacc Vector<M> get_data(const vector<M>& v)
{
  return v.v;
}

template <class M>
qacc Vector<Char> get_data_char(const vector<M>& v)
{
  Vector<Char> vc((Char*)v.v.p, v.v.data_size());
  return vc;
}

template <class M>
qacc void set_zero(vector<M>& xx)
{
  Vector<Char> vec = get_data_char(xx);
#ifdef QLAT_IN_ACC
  memset(vec.data(), 0, vec.size());
#else
  set_mem(vec.data(), 0, vec.size(), xx.mem_type);
#endif
}

// --------------------

template <class M>
struct API box {
  //
  // like a one element vector
  //
  // Avoid copy constructor when possible
  // (it is likely not what you think it is)
  // Only used in qacc macros, or if it is already a copy.
  //
  bool is_copy;              // do not free memory if is_copy=true
  mutable MemType mem_type;  // if place data on accelerator memory
  Handle<M> v;
  //
  box()
  {
    // TIMER("box::box()");
    qassert(v.p == NULL);
    is_copy = false;
    mem_type = get_default_mem_type();
  }
  box(const box<M>& vp)
  {
    // TIMER("box::box(&)");
#ifndef QLAT_USE_ACC
    qassert(vp.is_copy);
#endif
    is_copy = true;
    const MemType eff_mem_type = get_eff_mem_type(vp.mem_type);
    qassert(eff_mem_type == MemType::Uvm or eff_mem_type == MemType::Acc);
    mem_type = vp.mem_type;
    v = vp.v;
  }
  box(box<M>&& vp) noexcept
  {
    // TIMER("box::box(&&)");
    is_copy = vp.is_copy;
    mem_type = vp.mem_type;
    v = vp.v;
    vp.is_copy = true;
  }
  box(const M& x, MemType mem_type_ = get_default_mem_type())
  {
    // TIMER("box::box(x)");
    qassert(v.p == NULL);
    is_copy = false;
    mem_type = mem_type_;
    set(x);
  }
  //
  ~box()
  {
    // TIMER("box::~box()");
    if (not is_copy) {
      clear();
    }
  }
  //
  void init()
  // does not change mem_type
  {
    if (not is_copy) {
      clear();
    } else {
      v = Handle<M>();
    }
    is_copy = false;
  }
  //
  void clear()
  {
    qassert(not is_copy);
    if (v.p != NULL) {
      free_mem(v.p, sizeof(M), mem_type);
    }
    v = Handle<M>();
    qassert(v.p == NULL);
  }
  //
  void set_mem_type(const MemType mem_type_) const
  {
    qassert(not is_copy);
    if (NULL == v.p) {
      mem_type = mem_type_;
      return;
    }
    if (is_same_mem_type(mem_type, mem_type_)) {
      mem_type = mem_type_;
      return;
    }
    {
      TIMER("box::set_mem_type");
      box<M> b;
      b.mem_type = mem_type_;
      b = *this;
      std::swap(mem_type, b.mem_type);
      std::swap(v.p, b.v.p);
      qassert(mem_type == mem_type_);
    }
  }
  //
  void set_view(const M& x)
  // does not change mem_type
  {
    qassert(mem_type != MemType::Acc);
    init();
    is_copy = true;
    // Be cautious about the const property
    // 改不改靠自觉
    v.p = const_cast<M*>(&x);
  }
  void set_view(const Handle<M> h)
  // does not change mem_type
  {
    init();
    is_copy = true;
    v = h;
  }
  void set_view(const box<M>& b)
  {
    init();
    is_copy = true;
    mem_type = b.mem_type;
    v = b.v;
  }
  //
  void alloc()
  {
    qassert(not is_copy);
    if (v.p == NULL) {
      v.p = (M*)alloc_mem(sizeof(M), mem_type);
    }
  }
  //
  void set(const M& x)
  {
    alloc();
    copy_mem(v.p, mem_type, &x, MemType::Cpu, sizeof(M));
  }
  //
  qacc M get() const
  {
    qassert(not v.null());
#ifdef QLAT_USE_ACC
#ifdef QLAT_IN_ACC
    assert(mem_type == MemType::Acc or mem_type == MemType::Uvm);
    return v();
#else
    M x;
    copy_mem(&x, MemType::Cpu, v.p, mem_type, sizeof(M));
    return x;
#endif
#else
    return v();
#endif
  }
  //
  box<M>& operator=(const box<M>& vp)
  {
    // TIMER("box::operator=(&)");
    qassert(not is_copy);
    alloc();
    copy_mem(v.p, mem_type, vp.v.p, vp.mem_type, sizeof(M));
    return *this;
  }
  box<M>& operator=(box<M>&& vp) noexcept
  {
    // TIMER("box::operator=(&&)");
    is_copy = vp.is_copy;
    mem_type = vp.mem_type;
    v = vp.v;
    vp.is_copy = true;
    return *this;
  }
  //
  qacc const M& operator()() const { return v(); }
  qacc M& operator()() { return v(); }
  //
  qacc bool null() const { return v.null(); }
};

template <class M>
void clear(box<M>& v)
{
  v.clear();
}

template <class M>
qacc void qswap(box<M>& v1, box<M>& v2)
{
  qassert(not v1.is_copy);
  qassert(not v2.is_copy);
  std::swap(v1.mem_type, v2.mem_type);
  std::swap(v1.v, v2.v);
}

// --------------------

template <class M>
struct IsDataVectorType<box<M>> {
  using DataType = M;
  using BasicDataType = typename IsDataValueType<DataType>::BasicDataType;
  using ElementaryType = typename IsDataValueType<DataType>::ElementaryType;
  static constexpr bool value = is_data_value_type<DataType>();
};

template <class M>
qacc Vector<M> get_data(const box<M>& v)
{
  if (not v.null()) {
    return get_data(v());
  } else {
    return Vector<M>();
  }
}

}  // namespace qlat
