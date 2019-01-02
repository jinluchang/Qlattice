#pragma once

#include <qutils/lat-io.h>
#include <qutils/qutils.h>
#include <qlat/utils.h>

namespace qutils {

inline void clear(latio::LatData& ld)
{
  clear(ld.info);
  clear(ld.res);
}

}

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

inline LatDim lat_dim_number(const std::string& name, const long start, const long end, const long inc = 1)
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
  qassert(dim.indices.size() <= dim.size);
  for (long i = 0; i < dim.indices.size(); ++i) {
    if (idx == dim.indices[i]) {
      return i;
    }
  }
  const long i = -read_long(idx)-1;
  qassert(dim.indices.size() <= i and i < dim.size);
  return i;
}

inline long lat_dim_idx(const LatDim& dim, const long& idx)
{
  qassert(dim.indices.size() <= dim.size);
  qassert(0 <= idx and idx < dim.size);
  return idx;
}

template <class M>
inline void set_size(std::vector<M>& vec, const long size)
{
  vec.resize(size);
}

template <class M, unsigned long N>
inline void set_size(std::array<M,N>& vec, const long size)
{
  qassert(N == size);
}

template <class Vec, class VecS>
inline Vec lat_data_idx(const LatInfo& info, const VecS& idx)
  // Vec can be std::vector<long>
  // VecS can be std::vector<std::string> or std::vector<long>
  // or both can be std::array of same length
{
  qassert(idx.size() <= info.size());
  Vec ret;
  set_size(ret, idx.size());
  for (long i = 0; i < ret.size(); ++i) {
    ret[i] = lat_dim_idx(info[i], idx[i]);
  }
  return ret;
}

template <class VecS>
inline long lat_data_offset(const LatInfo& info, const VecS& idx)
  // will return offset at the level the idx specify
  // VecS can be std::vector<std::string> or std::vector<long>
  // or can be std::array of certain length
{
  qassert(idx.size() <= info.size());
  long ret = 0;
  for (int i = 0; i < idx.size(); ++i) {
    const long k = lat_dim_idx(info[i], idx[i]);
    ret = ret * info[i].size + k;
  }
  return ret;
}

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

inline bool is_lat_info_complex(const LatInfo& info)
{
  qassert(info.size() >= 1);
  const LatDim& dim = info.back();
  if (dim.name != "re-im" or dim.size != 2) {
    return false;
  } else if (dim.indices.size() !=2 ) {
    return false;
  } else if (dim.indices[0] != "re" or dim.indices[1] != "im") {
    return false;
  } else {
    return true;
  }
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
  Vector<Complex> ret((Complex*)&ld.res[offset * size], size/2);
  return ret;
}

template <class VecS>
inline Vector<Complex> lat_data_complex_get_const(const LatData& ld, const VecS& idx)
  // Be cautious about the const property
  // 改不改靠自觉
{
  qassert(is_lat_info_complex(ld.info));
  qassert(idx.size() < ld.info.size());
  const long offset = lat_data_offset(ld.info, idx);
  const long size = lat_data_size(ld.info, idx.size());
  qassert(size % 2 == 0);
  qassert(offset * size + size <= ld.res.size());
  Vector<Complex> ret((Complex*)&ld.res[offset * size], size/2);
  return ret;
}

QLAT_END_NAMESPACE
