// vim: set ts=2 sw=2 expandtab:

#pragma once

#include <map>
#include <set>
#include <utility>

#include <qlat-utils/timer.h>
#include <qlat-utils/vector.h>

namespace qlat
{  //

struct CacheBase {
  std::string name;
  long limit;
  long buffer_size;  // number of empty slots created by gc()
  //
  virtual ~CacheBase(){};
  //
  virtual long size() = 0;
  virtual void clear() = 0;
  //
  void init(const long limit_, const long buffer_size_ = 0)
  // If buffer_size_ == 0, then buffer_size = limit - limit / 2, which is also
  // the maximally allowed value
  {
    if (limit_ <= 0) {
      qassert(limit_ == 0);
      limit = 0;
      buffer_size = 1;
    } else {
      qassert(limit_ >= 1);
      limit = limit_;
      buffer_size = limit - limit / 2;
      if (buffer_size_ > 0) {
        qassert(buffer_size_ <= buffer_size);
        buffer_size = buffer_size_;
      }
      qassert(buffer_size >= 1);
    }
  }
  //
  virtual void resize(const long limit_, const long buffer_size_ = 0)
  {
    init(limit_, buffer_size_);
  }
  //
  virtual std::string show_info()
  {
    return ssprintf("Cache: '%s': %ld / %ld [%ld]", name.c_str(), size(), limit,
                    buffer_size);
  }
};

API inline std::set<CacheBase*>& get_all_caches()
{
  static std::set<CacheBase*> all_caches;
  return all_caches;
}

inline std::vector<std::string> get_all_caches_info()
{
  TIMER("get_all_caches_info");
  std::vector<std::string> infos;
  for (auto it = get_all_caches().cbegin(); it != get_all_caches().cend();
       ++it) {
    infos.push_back((*it)->show_info());
  }
  return infos;
}

inline void clear_all_caches()
// (1) Clear all cache using CacheBase (or Cache which internally uses
// CacheBase).
//
// (2) Clear memory cache (for vector and vector_acc memory freed but not
// cached).
{
  TIMER_VERBOSE("clear_all_caches");
  for (auto it = get_all_caches().cbegin(); it != get_all_caches().cend();
       ++it) {
    (*it)->clear();
  }
  clear_mem_cache();
}

template <class K, class M>
struct Cache : CacheBase {
  std::map<K, std::pair<long, M> > m;
  long idx;
  //
  Cache(const std::string& name_ = "Cache", const long limit_ = 16,
        const long buffer_size_ = 0)
  {
    init(name_, limit_, buffer_size_);
  }
  //
  ~Cache() noexcept
  {
    assert(qlat::has(get_all_caches(), (CacheBase*)this));
    get_all_caches().erase(this);
  }
  //
  void init(const std::string& name_, const long limit_,
            const long buffer_size_ = 0)
  {
    name = name_;
    CacheBase::init(limit_, buffer_size_);
    clear();
    if (not qlat::has(get_all_caches(), (CacheBase*)this)) {
      get_all_caches().insert(this);
    }
  }
  //
  bool has(const K& key) const { return qlat::has(m, key); }
  //
  long erase(const K& key) { return m.erase(key); }
  //
  M& operator[](const K& key)
  {
    typename std::map<K, std::pair<long, M> >::iterator it = m.find(key);
    if (it != m.end()) {
      if ((it->second).first != idx - 1) {
        (it->second).first = idx;
        idx += 1;
      }
      return (it->second).second;
    } else {
      gc();
      displayln_info(0, show_info() + " to add");
      std::pair<long, M>& v = m[key];
      v.first = idx;
      idx += 1;
      return v.second;
    }
  }
  //
  void gc()
  {
    if (limit > 0 and (long) m.size() >= limit) {
      TIMER_VERBOSE("Cache::gc");
      displayln_info(0, show_info() + " before gc");
      std::vector<long> idxes;
      for (typename std::map<K, std::pair<long, M> >::iterator it = m.begin();
           it != m.end(); ++it) {
        const long i = (it->second).first;
        idxes.push_back(i);
      }
      std::sort(idxes.begin(), idxes.end());
      qassert((long)m.size() == limit);
      qassert((long)m.size() > buffer_size);
      const long threshhold = idxes[buffer_size];
      std::vector<K> to_free;
      for (typename std::map<K, std::pair<long, M> >::iterator it = m.begin();
           it != m.end(); ++it) {
        const K& k = it->first;
        const long i = (it->second).first;
        if (i <= threshhold) {
          to_free.push_back(k);
        }
      }
      for (size_t i = 0; i < to_free.size(); ++i) {
        m.erase(to_free[i]);
      }
      displayln_info(0, show_info() + " after gc");
    }
  }
  //
  long size() { return m.size(); }
  //
  void clear()
  {
    TIMER_VERBOSE("Cache::clear");
    idx = 0;
    if (m.size() > 0) {
      displayln_info(0, show_info() + " clear");
      m.clear();
    }
  }
};

}  // namespace qlat
