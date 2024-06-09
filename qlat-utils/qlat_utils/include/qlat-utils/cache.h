// vim: set ts=2 sw=2 expandtab:

#pragma once

#include <map>
#include <set>
#include <utility>

#include <qlat-utils/timer.h>
#include <qlat-utils/vector.h>

namespace qlat
{  //

struct API CacheBase {
  std::string name;
  Long limit;
  Long buffer_size;  // number of empty slots created by gc()
  //
  virtual ~CacheBase(){};
  //
  virtual Long size() = 0;
  virtual void clear() = 0;
  //
  void init(const Long limit_, const Long buffer_size_ = 0)
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
  virtual void resize(const Long limit_, const Long buffer_size_ = 0)
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

template <class K, class M>
struct API Cache : CacheBase {
  std::map<K, std::pair<Long, M> > m;
  Long idx;
  //
  Cache(const std::string& name_ = "Cache", const Long limit_ = 16,
        const Long buffer_size_ = 0)
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
  void init(const std::string& name_, const Long limit_,
            const Long buffer_size_ = 0)
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
  Long erase(const K& key) { return m.erase(key); }
  //
  M& operator[](const K& key)
  {
    typename std::map<K, std::pair<Long, M> >::iterator it = m.find(key);
    if (it != m.end()) {
      if ((it->second).first != idx - 1) {
        (it->second).first = idx;
        idx += 1;
      }
      return (it->second).second;
    } else {
      gc();
      displayln_info(0, show_info() + " to add");
      std::pair<Long, M>& v = m[key];
      v.first = idx;
      idx += 1;
      return v.second;
    }
  }
  //
  void gc()
  {
    if (limit > 0 and (Long) m.size() >= limit) {
      TIMER_VERBOSE("Cache::gc");
      displayln_info(0, show_info() + " before gc");
      std::vector<Long> idxes;
      for (typename std::map<K, std::pair<Long, M> >::iterator it = m.begin();
           it != m.end(); ++it) {
        const Long i = (it->second).first;
        idxes.push_back(i);
      }
      std::sort(idxes.begin(), idxes.end());
      qassert((Long)m.size() == limit);
      qassert((Long)m.size() > buffer_size);
      const Long threshhold = idxes[buffer_size];
      std::vector<K> to_free;
      for (typename std::map<K, std::pair<Long, M> >::iterator it = m.begin();
           it != m.end(); ++it) {
        const K& k = it->first;
        const Long i = (it->second).first;
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
  Long size() { return m.size(); }
  //
  void clear()
  {
    idx = 0;
    if (m.size() > 0) {
      TIMER_VERBOSE("Cache::clear");
      displayln_info(0, show_info() + " clear");
      m.clear();
    }
  }
};

}  // namespace qlat
