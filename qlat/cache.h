// vim: set ts=2 sw=2 expandtab:

#pragma once

#include <qlat/config.h>

#include <map>
#include <utility>

QLAT_START_NAMESPACE

template <class K, class M>
struct Cache {
  std::string name;
  std::map<K, std::pair<long, M> > m;
  long idx;
  long limit;
  //
  Cache(const std::string& name_ = "Cache")
  {
    name = name_;
    idx = 0;
    limit = 16;
  }
  Cache(const std::string& name_, const long limit_)
  {
    name = name_;
    idx = 0;
    limit = limit_;
  }
  //
  bool has(const K& key) const
  {
    typename std::map<K, std::pair<long, M> >::const_iterator it = m.find(key);
    return it != m.end();
  }
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
      displayln_info(ssprintf("%s::%s: to add %d / %d.", cname().c_str(),
                              name.c_str(), m.size() + 1, limit));
      std::pair<long, M>& v = m[key];
      v.first = idx;
      idx += 1;
      return v.second;
    }
  }
  //
  void gc()
  {
    if (m.size() >= limit) {
      TIMER_VERBOSE("Cache::gc");
      displayln_info(ssprintf("%s::%s: before gc: %d / %d.", cname().c_str(),
                              name.c_str(), m.size(), limit));
      std::vector<int> idxes;
      for (typename std::map<K, std::pair<long, M> >::iterator it = m.begin();
           it != m.end(); ++it) {
        const K& k = it->first;
        const long i = (it->second).first;
        idxes.push_back(i);
      }
      std::sort(idxes.begin(), idxes.end());
      const long threshhold = idxes[m.size() / 2];
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
      displayln_info(ssprintf("%s::%s:  after gc: %d / %d.", cname().c_str(),
                              name.c_str(), m.size(), limit));
    }
  }
  //
  void clear() { m.clear(); }
};

QLAT_END_NAMESPACE
