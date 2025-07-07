#pragma once

#include <qlat-utils/config.h>
#include <qlat-utils/handle.h>
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

enum struct MemType : Int {
  Cpu,   // CPU main memory
  Acc,   // Accelerator
  Comm,  // For communication on CPU
  Uvm,   // Uniform virtual memory
  SIZE,
};

std::string show(const MemType mem_type);

MemType read_mem_type(const std::string& mem_type_str);

API inline Long& get_alignment()
// qlat parameter
//
// Should NOT change in the middle of the run.
{
  static Long alignment = 256;
  return alignment;
}

inline Long get_aligned_mem_size(const Long alignment, const Long min_size)
{
  const Long n_elem = 1 + (min_size - 1) / alignment;
  const Long size = n_elem * alignment;
  return size;
}

API inline Long& get_mem_cache_max_size(const MemType mem_type = MemType::Cpu)
// qlat parameter
// unit in MB
{
  static Long max_size = get_env_long_default("q_mem_cache_max_size", 512);
  static Long max_size_cpu =
      get_env_long_default("q_mem_cache_cpu_max_size", max_size) * 1024L *
      1024L;
  static Long max_size_acc =
      get_env_long_default("q_mem_cache_acc_max_size", max_size) * 1024L *
      1024L;
  static Long max_size_comm =
      get_env_long_default("q_mem_cache_comm_max_size", max_size) * 1024L *
      1024L;
  static Long max_size_uvm =
      get_env_long_default("q_mem_cache_uvm_max_size", max_size) * 1024L *
      1024L;
  if (mem_type == MemType::Cpu) {
    return max_size_cpu;
  } else if (mem_type == MemType::Acc) {
    return max_size_acc;
  } else if (mem_type == MemType::Comm) {
    return max_size_comm;
  } else if (mem_type == MemType::Uvm) {
    return max_size_uvm;
  } else {
    qassert(false);
    return max_size;
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
  MemCache(const MemType mem_type_ = MemType::Cpu,
           const Long mem_cache_max_size_ = -1)
  {
    init(mem_type_, mem_cache_max_size_);
  }
  ~MemCache() { gc(); }
  //
  void init(const MemType mem_type_ = MemType::Cpu,
            const Long mem_cache_max_size_ = -1);
  void add(void* ptr, const Long size);
  void* del(const Long size);
  void gc();
};

std::vector<MemCache> mk_mem_cache_vec();

API inline MemCache& get_mem_cache(const MemType mem_type = MemType::Cpu)
{
  static std::vector<MemCache> cache_vec = mk_mem_cache_vec();
  return cache_vec[static_cast<Int>(mem_type)];
}

void* alloc_mem(const Long min_size, const MemType mem_type);

void free_mem(void* ptr, const Long min_size, const MemType mem_type);

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
  MemType mem_type;  // if place data on accelerator memory
  bool is_copy;  // do not free memory if is_copy=true
  //
  vector()
  {
    // TIMER("vector::vector()")
    mem_type = MemType::Cpu;
    qassert(v.p == NULL);
    is_copy = false;
  }
  vector(const vector<M>& vp)
  {
    // TIMER("vector::vector(&)")
#ifndef QLAT_USE_ACC
    qassert(vp.is_copy);
#endif
    is_copy = true;
    qassert(vp.mem_type == MemType::Uvm);
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
  vector(const Long size)
  {
    // TIMER("vector::vector(size)")
    qassert(v.p == NULL);
    is_copy = false;
    mem_type = MemType::Cpu;
    resize(size);
  }
  vector(const Long size, const M& x)
  {
    // TIMER("vector::vector(size,x)")
    qassert(v.p == NULL);
    is_copy = false;
    mem_type = MemType::Cpu;
    resize(size, x);
  }
  vector(const std::vector<M>& vp)
  {
    // TIMER("vector::vector(std::vector&)")
    is_copy = false;
    mem_type = MemType::Cpu;
    *this = vp;
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
    }
    is_copy = false;
  }
  //
  void clear()
  {
    qassert(not is_copy);
    if (v.p != NULL) {
      free_mem(v.p, v.n * sizeof(M), mem_type);
    }
    v = Vector<M>();
    qassert(v.p == NULL);
  }
  //
  void set_mem_type(const MemType mem_type_)
  {
    qassert(not is_copy);
    if (mem_type == mem_type_) {
      return;
    }
    if (NULL == v.p) {
      mem_type = mem_type_;
      return;
    }
    vector<M> vec;
    vec.set_mem_type(mem_type_);
    vec = *this;
    swap(vec);
    qassert(mem_type == mem_type_);
  }
  //
  qacc void swap(vector<M>& x)
  {
    qassert(not is_copy);
    qassert(not x.is_copy);
    const MemType tb = x.mem_type;
    x.mem_type = mem_type;
    mem_type = tb;
    Vector<M> t = v;
    v = x.v;
    x.v = t;
  }
  //
  void set_view(const Vector<M>& vec)
  // does not change mem_type
  {
    init();
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
  void resize(const Long size)
  {
    qassert(not is_copy);
    qassert(0 <= size);
    if (v.p == NULL) {
      v.p = (M*)alloc_mem(size * sizeof(M), mem_type);
      v.n = size;
    } else if (v.n != size) {
      vector<M> vp;
      vp.set_mem_type(mem_type);
      vp.v = v;
      v.p = (M*)alloc_mem(size * sizeof(M), mem_type);
      v.n = size;
      if (size <= vp.v.n) {
        std::memcpy((void*)v.p, (void*)vp.v.p, size * sizeof(M));
      } else {
        std::memcpy((void*)v.p, (void*)vp.v.p, vp.v.n * sizeof(M));
      }
    }
  }
  void resize(const Long size, const M& x)
  {
    qassert(not is_copy);
    qassert(0 <= size);
    if (v.p == NULL) {
      v.p = (M*)alloc_mem(size * sizeof(M), mem_type);
      v.n = size;
      for (Long i = 0; i < v.n; ++i) {
        v[i] = x;
      }
    } else if (v.n != size) {
      vector<M> vp;
      vp.set_mem_type(mem_type);
      vp.v = v;
      v.p = (M*)alloc_mem(size * sizeof(M), mem_type);
      v.n = size;
      if (size <= vp.v.n) {
        std::memcpy(v.p, vp.v.p, v.n * sizeof(M));
      } else {
        std::memcpy(v.p, vp.v.p, vp.v.n * sizeof(M));
        for (Long i = vp.v.n; i < v.n; ++i) {
          v[i] = x;
        }
      }
    }
  }
  //
  vector<M>& operator=(const vector<M>& vp)
  {
    // TIMER("vector::operator=(&)");
    qassert(not is_copy);
    clear();
    resize(vp.size());
    for (Long i = 0; i < v.n; ++i) {
      v[i] = vp[i];
    }
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
    for (Long i = 0; i < v.n; ++i) {
      v[i] = vp[i];
    }
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
struct API vector_acc : vector<M> {
  // Avoid copy constructor when possible
  // (it is likely not what you think it is)
  // Only used in qacc macros.
  //
  using vector<M>::v;
  using vector<M>::is_copy;
  using vector<M>::mem_type;
  using vector<M>::resize;
  //
  vector_acc()
  {
    // TIMER("vector_acc::vector_acc()");
    qassert(v.p == NULL);
    is_copy = false;
    mem_type = MemType::Uvm;
  }
  vector_acc(const vector_acc<M>& vp)
  {
    // TIMER("vector_acc::vector_acc(&)");
#ifndef QLAT_USE_ACC
    qassert(vp.is_copy);
#endif
    is_copy = true;
    qassert(vp.mem_type == MemType::Uvm);
    mem_type = vp.mem_type;
    v = vp.v;
  }
  vector_acc(vector_acc<M>&& vp) noexcept
  {
    // TIMER("vector_acc::vector_acc(&&)")
    // qassert(vp.mem_type == MemType::Uvm);
    is_copy = vp.is_copy;
    mem_type = vp.mem_type;
    v = vp.v;
    vp.is_copy = true;
  }
  vector_acc(const Long size)
  {
    // TIMER("vector_acc::vector_acc(size)");
    qassert(v.p == NULL);
    is_copy = false;
    mem_type = MemType::Uvm;
    resize(size);
  }
  vector_acc(const Long size, const M& x)
  {
    // TIMER("vector_acc::vector_acc(size,x)");
    qassert(v.p == NULL);
    is_copy = false;
    mem_type = MemType::Uvm;
    resize(size, x);
  }
  vector_acc(const std::vector<M>& vp)
  {
    // TIMER("vector_acc::vector_acc(std::vector&)");
    is_copy = false;
    mem_type = MemType::Uvm;
    *this = vp;
  }
  //
  vector_acc<M>& operator=(const vector_acc<M>& vp)
  {
    // TIMER("vector_acc::operator=(&)");
    vector<M>::operator=(vp);
    return *this;
  }
  vector_acc<M>& operator=(const vector<M>& vp)
  {
    // TIMER("vector_acc::operator=(&)");
    vector<M>::operator=(vp);
    return *this;
  }
  vector_acc<M>& operator=(vector_acc<M>&& vp) noexcept
  {
    // TIMER("vector_acc::operator=(&&)");
    vector<M>::operator=(std::move(vp));
    return *this;
  }
  vector_acc<M>& operator=(vector<M>&& vp) noexcept
  {
    // TIMER("vector_acc::operator=(&&)");
    vector<M>::operator=(std::move(vp));
    return *this;
  }
  vector_acc<M>& operator=(const std::vector<M>& vp)
  {
    // TIMER("vector_acc::operator=(std::vector&)");
    vector<M>::operator=(vp);
    return *this;
  }
};

template <class M>
void clear(vector<M>& v)
{
  v.clear();
}

template <class M>
qacc void qswap(vector<M>& v1, vector<M>& v2)
{
  v1.swap(v2);
}

// --------------------

template <class M>
qacc bool operator==(const vector<M>& v1, const vector<M>& v2)
{
  if (v1.size() != v2.size()) {
    return false;
  }
  const int cmp = std::memcmp(v1.data(), v2.data(), v1.size() * sizeof(M));
  return cmp == 0;
}

template <class M>
qacc bool operator==(const vector_acc<M>& v1, const vector_acc<M>& v2)
{
  if (v1.size() != v2.size()) {
    return false;
  }
  const int cmp = std::memcmp(v1.data(), v2.data(), v1.size() * sizeof(M));
  return cmp == 0;
}

template <class M>
qacc bool operator!=(const vector<M>& v1, const vector<M>& v2)
{
  return not(v1 == v2);
}

template <class M>
qacc bool operator!=(const vector_acc<M>& v1, const vector_acc<M>& v2)
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
struct IsDataVectorType<vector_acc<M>> {
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
qacc Vector<M> get_data(const vector_acc<M>& v)
{
  return get_data((const vector<M>&)v);
}

// --------------------

template <class M>
struct API box {
  //
  // like a one element vector
  //
  // Avoid copy constructor when possible
  // (it is likely not be what you think it is)
  // Only used in qacc macros, or if it is already a copy.
  //
  bool is_copy;  // do not free memory if is_copy=true
  MemType mem_type;  // if place data on accelerator memory
  Handle<M> v;
  //
  box()
  {
    // TIMER("box::box()");
    qassert(v.p == NULL);
    is_copy = false;
    mem_type = MemType::Cpu;
  }
  box(const box<M>& vp)
  {
    // TIMER("box::box(&)");
#ifndef QLAT_USE_ACC
    qassert(vp.is_copy);
#endif
    is_copy = true;
    qassert(vp.mem_type == MemType::Uvm);
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
  box(const M& x)
  {
    // TIMER("box::box(x)");
    qassert(v.p == NULL);
    is_copy = false;
    mem_type = MemType::Uvm;
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
  void set_mem_type(const MemType mem_type_)
  {
    qassert(not is_copy);
    if (mem_type == mem_type_) {
      return;
    }
    if (NULL == v.p) {
      mem_type = mem_type_;
      return;
    }
    box<M> b;
    b.set_mem_type(mem_type_);
    b = *this;
    swap(b);
    qassert(mem_type == mem_type_);
  }
  //
  qacc void swap(box<M>& x)
  {
    qassert(not is_copy);
    qassert(not x.is_copy);
    const MemType tb = x.mem_type;
    x.mem_type = mem_type;
    mem_type = tb;
    Handle<M> t = v;
    v = x.v;
    x.v = t;
  }
  //
  void set_view(const M& x)
  // does not change mem_type
  {
    init();
    is_copy = true;
    v.p = &x;
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
  void set(const M& x)
  {
    qassert(not is_copy);
    if (v.p == NULL) {
      v.p = (M*)alloc_mem(sizeof(M), mem_type);
    }
    v() = x;
  }
  //
  box<M>& operator=(const box<M>& vp)
  {
    // TIMER("box::operator=(&)");
    qassert(not is_copy);
    set(vp());
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
struct API box_acc : box<M> {
  //
  // like a one element vector_acc
  //
  // Avoid copy constructor when possible
  // (it is likely not be what you think it is)
  // Only used in qacc macros, or if it is already a copy.
  //
  using box<M>::v;
  using box<M>::is_copy;
  using box<M>::mem_type;
  using box<M>::set;
  //
  box_acc()
  {
    // TIMER("box_acc::box_acc()");
    qassert(v.p == NULL);
    is_copy = false;
    mem_type = MemType::Uvm;
  }
  box_acc(const box_acc<M>& vp)
  {
    // TIMER("box::box(&)");
#ifndef QLAT_USE_ACC
    qassert(vp.is_copy);
#endif
    is_copy = true;
    qassert(vp.mem_type == MemType::Uvm);
    mem_type = vp.mem_type;
    v = vp.v;
  }
  box_acc(box_acc<M>&& vp) noexcept
  {
    // TIMER("box_acc::box_acc(&&)");
    // qassert(vp.mem_type == MemType::Uvm);
    is_copy = vp.is_copy;
    mem_type = vp.mem_type;
    v = vp.v;
    vp.is_copy = true;
  }
  box_acc(const M& x)
  {
    // TIMER("box_acc::box_acc(x)");
    qassert(v.p == NULL);
    is_copy = false;
    mem_type = MemType::Uvm;
    set(x);
  }
  //
  box_acc<M>& operator=(const box_acc<M>& vp)
  {
    // TIMER("box_acc::operator=(&)");
    box<M>::operator=(vp);
    return *this;
  }
  box_acc<M>& operator=(const box<M>& vp)
  {
    // TIMER("box_acc::operator=(&)");
    box<M>::operator=(vp);
    return *this;
  }
  box_acc<M>& operator=(box_acc<M>&& vp) noexcept
  {
    // TIMER("box_acc::operator=(&&)");
    box<M>::operator=(std::move(vp));
    return *this;
  }
  box_acc<M>& operator=(box<M>&& vp) noexcept
  {
    // TIMER("box_acc::operator=(&&)");
    box<M>::operator=(std::move(vp));
    return *this;
  }
};

template <class M>
void clear(box<M>& v)
{
  v.clear();
}

template <class M>
qacc void qswap(box<M>& v1, box<M>& v2)
{
  v1.swap(v2);
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
struct IsDataVectorType<box_acc<M>> {
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

template <class M>
qacc Vector<M> get_data(const box_acc<M>& v)
{
  return get_data((const box<M>&)v);
}

}  // namespace qlat
