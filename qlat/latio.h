#pragma once

#include <qlat/utils.h>
#include <qutils/lat-io.h>
#include <qutils/qutils.h>

namespace qutils
{  //

inline void clear(latio::LatData& ld)
{
  clear(ld.info);
  clear(ld.res);
}

}  // namespace qutils

QLAT_START_NAMESPACE

inline LatDim lat_dim_re_im()
{
  LatDim dim;
  dim.name = "re-im";
  dim.size = 2;
  dim.indices.resize(2);
  dim.indices[0] = "re";
  dim.indices[1] = "im";
  return dim;
}

inline LatDim lat_dim_number(const std::string& name, const long start,
                             const long end, const long inc = 1)
{
  using namespace qshow;
  LatDim dim;
  dim.name = name;
  for (long i = start; i <= end; i += inc) {
    dim.size += 1;
    dim.indices.push_back(show(i));
  }
  return dim;
}

inline long lat_dim_idx(const LatDim& dim, const std::string& idx)
{
  assert(dim.indices.size() <= dim.size);
  for (long i = 0; i < dim.indices.size(); ++i) {
    if (idx == dim.indices[i]) {
      return i;
    }
  }
  const long i = -read_long(idx) - 1;
  assert(dim.indices.size() <= i and i < dim.size);
  return i;
}

inline long lat_dim_idx(const LatDim& dim, const long& idx)
{
  assert(dim.indices.size() <= dim.size);
  assert(0 <= idx and idx < dim.size);
  return idx;
}

template <class VecS>
inline long lat_data_offset(const LatInfo& info, const VecS& idx)
// will return offset at the level the idx specify
// VecS can be std::vector<std::string> or std::vector<long>
// or can be std::array of certain length
{
  assert(idx.size() <= info.size());
  long ret = 0;
  for (int i = 0; i < idx.size(); ++i) {
    const long k = lat_dim_idx(info[i], idx[i]);
    ret = ret * info[i].size + k;
  }
  return ret;
}

inline bool is_lat_info_complex(const LatInfo& info)
{
  assert(info.size() >= 1);
  const LatDim& dim = info.back();
  if (dim.name != "re-im" or dim.size != 2) {
    return false;
  } else if (dim.indices.size() != 2) {
    return false;
  } else if (dim.indices[0] != "re" or dim.indices[1] != "im") {
    return false;
  } else {
    return true;
  }
}

inline bool operator==(const LatDim& d1, const LatDim& d2)
{
  return d1.name == d2.name and d1.size == d2.size and d1.indices == d2.indices;
}

inline bool operator!=(const LatDim& d1, const LatDim& d2)
{
  return not(d1 == d2);
}

inline bool is_matching(const LatData& ld1, const LatData& ld2)
{
  return ld1.res.size() == ld2.res.size() and ld1.info == ld2.info;
}

inline void set_zero(LatData& ld)
{
  std::memset(ld.res.data(), 0, ld.res.size() * sizeof(double));
}

inline const LatData& operator+=(LatData& ld, const LatData& ld1)
{
  assert(is_matching(ld, ld1));
  for (long i = 0; i < ld.res.size(); ++i) {
    ld.res[i] += ld1.res[i];
  }
  return ld;
}

inline const LatData& operator-=(LatData& ld, const LatData& ld1)
{
  assert(is_matching(ld, ld1));
  for (long i = 0; i < ld.res.size(); ++i) {
    ld.res[i] -= ld1.res[i];
  }
  return ld;
}

inline const LatData& operator*=(LatData& ld, const double factor)
{
  for (long i = 0; i < ld.res.size(); ++i) {
    ld.res[i] *= factor;
  }
  return ld;
}

inline std::string idx_name(const LatDim& dim, const long idx)
{
  if (idx < dim.indices.size()) {
    return dim.indices[idx];
  } else {
    return show(-idx-1);
  }
}

inline void print(const LatData& ld)
{
  const LatInfo& info = ld.info;
  display(ssprintf("%s", show(info).c_str()));
  std::vector<long> idx(info.size(), 0);
  for (long k = 0; k < lat_data_size(info); ++k) {
    for (int a = 0; a < info.size(); ++a) {
      display(ssprintf("%s[%8s] ", info[a].name.c_str(), idx_name(info[a], idx[a]).c_str()));
    }
    display(ssprintf("%24.17E\n", ld.res[lat_data_offset(info, idx)]));
    idx[info.size() - 1] += 1;
    for (int a = info.size() - 1; a > 0; --a) {
      if (idx[a] == info[a].size) {
        idx[a] = 0;
        idx[a-1] += 1;
      }
    }
  }
}

// qlat specific below

template <class VecS>
inline Vector<double> lat_data_get(LatData& ld, const VecS& idx)
{
  const long offset = lat_data_offset(ld.info, idx);
  const long size = lat_data_size(ld.info, idx.size());
  qassert(offset * size + size <= ld.res.size());
  Vector<double> ret(&ld.res[offset * size], size);
  return ret;
}

template <class VecS>
inline Vector<double> lat_data_get_const(const LatData& ld, const VecS& idx)
// Be cautious about the const property
// 改不改靠自觉
{
  const long offset = lat_data_offset(ld.info, idx);
  const long size = lat_data_size(ld.info, idx.size());
  qassert(offset * size + size <= ld.res.size());
  Vector<double> ret(&ld.res[offset * size], size);
  return ret;
}

template <class VecS>
inline Vector<Complex> lat_data_complex_get(LatData& ld, const VecS& idx)
{
  qassert(is_lat_info_complex(ld.info));
  qassert(idx.size() < ld.info.size());
  const long offset = lat_data_offset(ld.info, idx);
  const long size = lat_data_size(ld.info, idx.size());
  qassert(size % 2 == 0);
  qassert(offset * size + size <= ld.res.size());
  Vector<Complex> ret((Complex*)&ld.res[offset * size], size / 2);
  return ret;
}

template <class VecS>
inline Vector<Complex> lat_data_complex_get_const(const LatData& ld,
                                                  const VecS& idx)
// Be cautious about the const property
// 改不改靠自觉
{
  qassert(is_lat_info_complex(ld.info));
  qassert(idx.size() < ld.info.size());
  const long offset = lat_data_offset(ld.info, idx);
  const long size = lat_data_size(ld.info, idx.size());
  qassert(size % 2 == 0);
  qassert(offset * size + size <= ld.res.size());
  Vector<Complex> ret((Complex*)&ld.res[offset * size], size / 2);
  return ret;
}

inline int lat_data_glb_sum(LatData& ld)
{
  return glb_sum_double_vec(get_data(ld.res));
}

QLAT_END_NAMESPACE
