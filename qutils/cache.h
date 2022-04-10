// vim: set ts=2 sw=2 expandtab:

#pragma once

#include <map>
#include <set>
#include <utility>

#include <qutils/timer.h>

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
  virtual std::string show_info()
  {
    return ssprintf("Cache: '%s': %ld / %ld [%ld]", name.c_str(), size(), limit,
                    buffer_size);
  }
};

inline std::set<CacheBase*>& get_all_caches()
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
{
  TIMER_VERBOSE("clear_all_caches");
  for (auto it = get_all_caches().cbegin(); it != get_all_caches().cend();
       ++it) {
    (*it)->clear();
  }
}

template <class K, class M>
struct Cache : CacheBase {
  std::map<K, std::pair<long, M> > m;
  long idx;
  //
  Cache(const std::string& name_ = "Cache", const long limit_ = 16,
        const long buffer_size_ = 8)
  {
    init(name_, limit_, buffer_size_);
  }
  //
  ~Cache()
  {
    qassert(qlat::has(get_all_caches(), (CacheBase*)this));
    get_all_caches().erase(this);
  }
  //
  void init(const std::string& name_, const long limit_,
            const long buffer_size_)
  {
    name = name_;
    if (limit_ <= 0) {
      limit = 0;
      buffer_size = 1;
    } else {
      limit = std::max(limit_, (long)1);
      buffer_size = std::min(buffer_size_, limit - limit / 2);
      qassert(buffer_size >= 1);
    }
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
      displayln_info(
          ssprintf("%s: to add %d / %d.", name.c_str(), m.size() + 1, limit));
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
      displayln_info(0, show_info());
      m.clear();
    }
  }
};

}  // namespace qlat
