// Opt-in diagnostics for the custom MCG adapter. No QUDA policy changes.
#ifndef UTILS_MCG_DEBUG_H
#define UTILS_MCG_DEBUG_H
#pragma once
#include "general_funs.h"
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <exception>
#include <limits>
#include <utility>
#include <algorithm>

namespace qlat {
inline int mcg_debug_uncaught_exceptions() noexcept {
#if defined(__cpp_lib_uncaught_exceptions) && __cpp_lib_uncaught_exceptions >= 201411L
  return std::uncaught_exceptions();
#else
  // C++14 exposes only an unwinding flag. This detects a new exception when
  // a phase/solve starts outside unwinding, as normal solver calls do. It
  // cannot distinguish nested exceptions in a caller's unwinding destructor.
  return std::uncaught_exception() ? 1 : 0;
#endif
}

inline int mcg_debug_level() {
  static const int value = [] {
    const char* s = std::getenv("qlat_mcg_debug");
    if (!s) return 0;
    if (std::strcmp(s,"0") == 0) return 0;
    if (std::strcmp(s,"1") == 0) return 1;
    if (std::strcmp(s,"2") == 0) return 2;
    std::fprintf(stderr,"qlat_mcg_debug must be 0, 1, or 2\n");
    int initialized=0; MPI_Initialized(&initialized);
    if (initialized) MPI_Abort(MPI_COMM_WORLD,97);
    std::abort();
  }();
  return value;
}
struct MCGDebugState {
  bool active=false;
  long solve=0, next_solve=0, subcall=-1, rhs=0, masses=0, io=0, krylov=0;
  long window=0, iteration=0, operators=0;
  int rank=0;
  FILE* events=nullptr;
  char* current=nullptr;
  char record[2048]{};
  void open() {
    if (events) return;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    static const std::string prefix=[] {
      const char* s=std::getenv("qlat_mcg_debug_prefix");
      return std::string(s && *s?s:"mcg-debug");
    }();
    const std::string base=prefix+".rank"+std::to_string(rank)+".pid"+std::to_string(getpid());
    events=std::fopen((base+".log").c_str(),"w");
    int fd=::open((base+".current").c_str(),O_RDWR|O_CREAT|O_TRUNC,0600);
    if (!events || fd<0 || ftruncate(fd,2048)!=0) {
      std::fprintf(stderr,"MCG_DEBUG cannot create %s: %s\n",base.c_str(),std::strerror(errno));
      MPI_Abort(MPI_COMM_WORLD,97); std::abort();
    }
    void* mapping=mmap(nullptr,2048,PROT_READ|PROT_WRITE,MAP_SHARED,fd,0);
    close(fd);
    if (mapping==MAP_FAILED) {
      std::fprintf(stderr,"MCG_DEBUG mmap failed: %s\n",std::strerror(errno));
      MPI_Abort(MPI_COMM_WORLD,97); std::abort();
    }
    current=static_cast<char*>(mapping);
    std::setvbuf(events,nullptr,_IOLBF,0);
  }
  bool sampled() const { return iteration<=2 || iteration%1024==0; }
  void mark(const char* phase,const char* state,bool force=false) {
    if (!active) return;
    open();
    std::snprintf(record,sizeof(record),
      "MCG_DEBUG rank=%d pid=%ld solve=%ld subcall=%ld rhs=%ld masses=%ld io=%ld krylov=%ld "
      "window=%ld iter=%ld operators=%ld phase=%s state=%s\n",
      rank,long(getpid()),solve,subcall,rhs,masses,io,krylov,window,iteration,operators,phase,state);
    // MAP_SHARED retains the latest phase even if MPI or the allocator aborts.
    // A reader during execution may see a partial update; read again when stopped.
    std::memset(current,0,2048);
    std::memcpy(current,record,std::strlen(record));
    if (force || sampled()) std::fputs(record,events);
  }
};
inline MCGDebugState& mcg_debug_state() {
  // Intentionally process-lived: caches can be destroyed after MPI_Finalize.
  static auto* state=new MCGDebugState;
  return *state;
}
inline bool mcg_debug_active() { return mcg_debug_level() && mcg_debug_state().active; }
[[noreturn]] inline void mcg_debug_fail(const std::string& message) {
  auto& s=mcg_debug_state();
  std::fprintf(stderr,"MCG_DEBUG_ERROR rank=%d %s\n%s",s.rank,message.c_str(),s.record);
  if (s.events) {
    std::fprintf(s.events,"MCG_DEBUG_ERROR %s\n%s",message.c_str(),s.record);
    std::fflush(s.events);
  }
  if (s.current) msync(s.current,2048,MS_SYNC);
  int initialized=0,finalized=0;MPI_Initialized(&initialized);MPI_Finalized(&finalized);
  if (initialized && !finalized) MPI_Abort(MPI_COMM_WORLD,97);
  std::abort();
}
inline void mcg_debug_require(bool ok,const char* message) {
  if (mcg_debug_active() && !ok) mcg_debug_fail(message);
}
inline void mcg_debug_complete() {
#ifdef QLAT_USE_ACC
  // No guard copies or diagnostic kernels may be submitted before this returns.
  const auto error=qacc_DeviceSynchronize();
  if (error!=qacc_Success) mcg_debug_fail("GPU completion failed, code="+std::to_string(int(error)));
#endif
}
inline void mcg_debug_replay_steps(long steps,long window) {
  mcg_debug_require(steps>0 && steps<=window,"replay bounds: steps outside temporal window");
}
struct MCGDebugMaintenance {
  bool activate=false;
  explicit MCGDebugMaintenance(const char* phase) {
    if (!mcg_debug_level()) return;
    auto& s=mcg_debug_state();activate=!s.active;s.active=true;s.mark(phase,"enter",true);
  }
  ~MCGDebugMaintenance() { if(activate) mcg_debug_state().active=false; }
};
struct MCGDebugGuard {
  unsigned char *allocation,*payload;
  size_t bytes;
  int memory;
  std::string name;
};
inline std::vector<MCGDebugGuard>& mcg_debug_guards() {
  static auto* guards=new std::vector<MCGDebugGuard>;
  return *guards;
}
inline void mcg_debug_copy(void* dst,const void* src,size_t bytes,bool to_device) {
#ifdef QLAT_USE_ACC
  const auto error=qacc_Memcpy(dst,src,bytes,to_device?qacc_MemcpyHostToDevice:qacc_MemcpyDeviceToHost);
  if (error!=qacc_Success) mcg_debug_fail("guard copy failed, code="+std::to_string(int(error)));
#else
  (void)to_device;std::memcpy(dst,src,bytes);
#endif
}
inline void mcg_debug_check_guards() {
  if (!mcg_debug_active() || mcg_debug_level()<2) return;
  mcg_debug_complete();
  unsigned char bytes[256];
  for (const auto& g:mcg_debug_guards()) {
    for (int side=0;side<2;++side) {
      const void* p=side?g.payload+g.bytes:g.allocation;
      if (g.memory) mcg_debug_copy(bytes,p,256,false);
      else std::memcpy(bytes,p,256);
      for (size_t i=0;i<256;++i) if (bytes[i]!=0xa7)
        mcg_debug_fail("canary damaged: buffer="+g.name+" address="+
          std::to_string(reinterpret_cast<uintptr_t>(g.payload))+" bytes="+std::to_string(g.bytes)+
          " side="+(side?"trailing":"leading")+" guard_byte="+std::to_string(i));
    }
  }
}
struct MCGDebugPhase {
  const char* name;
  bool enabled;
  int exceptions;
  explicit MCGDebugPhase(const char* n):name(n),enabled(mcg_debug_active()),exceptions(enabled?mcg_debug_uncaught_exceptions():0) {
    if (enabled) mcg_debug_state().mark(name,"enter");
  }
  ~MCGDebugPhase() {
    if (!enabled) return;
    if (mcg_debug_uncaught_exceptions()>exceptions) mcg_debug_fail(std::string("exception in ")+name);
    mcg_debug_check_guards();
    mcg_debug_state().mark(name,"complete");
  }
};
struct MCGDebugSolve {
  bool enabled=false,owner=false;
  int exceptions=0;
  explicit MCGDebugSolve(long rhs,long masses,long io,long krylov,long window) {
    if (!mcg_debug_level()) return;
    enabled=true;auto& s=mcg_debug_state();owner=!s.active;
    if (!owner) return;
    s.active=true;s.solve=++s.next_solve;s.subcall=-1;s.rhs=rhs;s.masses=masses;s.io=io;
    s.krylov=krylov;s.window=window;s.iteration=0;s.operators=0;
    exceptions=mcg_debug_uncaught_exceptions();s.mark("public_solve","enter",true);
  }
  ~MCGDebugSolve() {
    if (!enabled || !owner) return;
    if (mcg_debug_uncaught_exceptions()>exceptions) mcg_debug_fail("exception escaping custom MCG solve");
    mcg_debug_check_guards();
    auto& s=mcg_debug_state();s.mark("public_solve","complete",true);s.active=false;
  }
};
struct MCGDebugSubcall {
  long old_rhs=0,old_subcall=0;
  explicit MCGDebugSubcall(long rhs,long subcall) {
    if (!mcg_debug_active()) return;
    auto& s=mcg_debug_state();old_rhs=s.rhs;old_subcall=s.subcall;s.rhs=rhs;s.subcall=subcall;
  }
  ~MCGDebugSubcall() {
    if (!mcg_debug_active()) return;
    auto& s=mcg_debug_state();s.rhs=old_rhs;s.subcall=old_subcall;
  }
};
inline void mcg_debug_device_pointer(const void* p,const char* name) {
  if (!mcg_debug_active()) return;
  if (!p) mcg_debug_fail(std::string("null pointer: ")+name);
#ifdef QLAT_USE_ACC
  qacc_PointerAttributes attr{};int device=-1;
  if (qacc_GetDevice(&device)!=qacc_Success || qacc_PointerGetAttributes(&attr,p)!=qacc_Success)
    mcg_debug_fail(std::string("cannot query pointer: ")+name);
  if ((attr.type!=qacc_MemoryTypeDevice && attr.type!=qacc_MemoryTypeManaged) || attr.device!=device)
    mcg_debug_fail(std::string("pointer is not on the current device: ")+name);
#endif
}
template<class T> inline void mcg_debug_public_buffers(const std::vector<T*>& out,
    const std::vector<T*>& in,size_t masses,size_t length) {
  if (!mcg_debug_active()) return;
  mcg_debug_require(!in.empty() && masses>0 && masses<=SIZE_MAX/in.size() && out.size()==masses*in.size(),
                    "public pointer-table sizes");
  mcg_debug_require(in.size()<=size_t(std::numeric_limits<Int>::max()) &&
    masses<=size_t(std::numeric_limits<Int>::max()),"RHS/mass count exceeds Int");
  mcg_debug_require(length>0 && length<=SIZE_MAX/sizeof(T) &&
    length<=size_t(std::numeric_limits<Long>::max())/in.size(),"public vector byte count overflow");
  const size_t bytes=length*sizeof(T);
  for (auto* p:in) mcg_debug_device_pointer(p,"source");
  for (size_t i=0;i<out.size();++i) {
    mcg_debug_device_pointer(out[i],"output");
    const uintptr_t a=reinterpret_cast<uintptr_t>(out[i]);
    mcg_debug_require(a<=UINTPTR_MAX-bytes,"output address range overflow");
    for (size_t j=0;j<i;++j) {
      const uintptr_t b=reinterpret_cast<uintptr_t>(out[j]);
      mcg_debug_require(a+bytes<=b || b+bytes<=a,"overlapping output ranges");
    }
  }
  // Source/output aliasing is supported: input is saved before outputs change.
}

// Derived only to retain compatibility with existing vector_gpu references.
// Levels 0/1 delegate to the original allocator, with identical logical sizes.
// Level 2 offsets the payload by 256 bytes; every QUDA alignment is preserved.
template<class T> struct MCGDebugBuffer:vector_gpu<T> {
  using Base=vector_gpu<T>;
  unsigned char* guarded=nullptr;
  std::string label="custom";
  MCGDebugBuffer()=default;
  explicit MCGDebugBuffer(size_t n,int memory=1) { resize(n,memory); }
  MCGDebugBuffer(const MCGDebugBuffer&)=delete;
  MCGDebugBuffer& operator=(const MCGDebugBuffer&)=delete;
  MCGDebugBuffer(MCGDebugBuffer&& other) noexcept:Base(std::move(other)),
    guarded(std::exchange(other.guarded,nullptr)),label(std::move(other.label)) {
    other.p=nullptr;other.n=0;other.is_copy=false;
  }
  ~MCGDebugBuffer() { clear(); }
  void clear() {
    if (guarded) {
      if (mcg_debug_active()) mcg_debug_check_guards();
      auto& guards=mcg_debug_guards();
      guards.erase(std::remove_if(guards.begin(),guards.end(),[&](const auto& g){return g.allocation==guarded;}),guards.end());
      free_buf(guarded,this->GPU);guarded=nullptr;this->p=nullptr;this->n=0;this->is_copy=false;
    } else if (!this->is_copy) Base::clean_mem();
  }
  void resize(size_t n,int memory) {
    if (n==this->n && memory==this->GPU) return;
    clear();
    if (!n) return;
    if (mcg_debug_level()<2) { Base::resize(n,memory);return; }
    // Include buffers warmed by direct helper calls before the public solve.
    MCGDebugMaintenance activity("guarded_allocation");
    mcg_debug_require(n<=(SIZE_MAX-767)/sizeof(T),"guarded allocation byte count overflow");
    mcg_debug_complete();
    const size_t bytes=n*sizeof(T),total=(bytes+512+255)/256*256;
    gpuMalloc(guarded,total,unsigned char,memory);
    this->p=reinterpret_cast<T*>(guarded+256);this->n=n;this->GPU=memory;this->is_copy=true;
    unsigned char canary[256];std::memset(canary,0xa7,256);
    if (memory) {
      mcg_debug_copy(guarded,canary,256,true);mcg_debug_copy(guarded+256+bytes,canary,256,true);
    } else {
      std::memcpy(guarded,canary,256);std::memcpy(guarded+256+bytes,canary,256);
    }
    Base::set_zero();
    mcg_debug_guards().push_back({guarded,guarded+256,bytes,memory,label});
  }
  void resize(size_t n) { resize(n,this->GPU); }
  void debug_name(const std::string& name) {
    if (!mcg_debug_active()) return;
    label=name;
    for (auto& g:mcg_debug_guards()) if (g.allocation==guarded) g.name=name;
    if (mcg_debug_active() && this->size()) {
      auto& s=mcg_debug_state();s.open();
      std::fprintf(s.events,"MCG_BUFFER solve=%ld name=%s address=%p elements=%zu bytes=%zu memory=%d guarded=%d\n",
        s.solve,name.c_str(),static_cast<void*>(this->p),this->n,this->n*sizeof(T),this->GPU,int(guarded!=nullptr));
    }
  }
};
} // namespace qlat
#endif
