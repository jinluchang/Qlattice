#pragma once

#include <qlat/utils.h>
#include <qutils/lat-io.h>
#include <qutils/qutils.h>

QLAT_START_NAMESPACE

inline bool operator==(const LatDim& d1, const LatDim& d2)
{
  return d1.name == d2.name and d1.size == d2.size and d1.indices == d2.indices;
}

inline bool operator!=(const LatDim& d1, const LatDim& d2)
{
  return not(d1 == d2);
}

inline void set_zero(LatData& ld)
{
  std::memset(ld.res.data(), 0, ld.res.size() * sizeof(double));
}

inline bool is_matching(const LatData& ld1, const LatData& ld2)
{
  return ld1.res.size() == ld2.res.size() and ld1.info == ld2.info;
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

inline LatData operator+(const LatData& ld, const LatData& ld1)
{
  LatData ret = ld;
  ret += ld1;
  return ret;
}

inline LatData operator-(const LatData& ld, const LatData& ld1)
{
  LatData ret = ld;
  ret -= ld1;
  return ret;
}

inline LatData operator*(const LatData& ld, const double factor)
{
  LatData ret = ld;
  for (long i = 0; i < ret.res.size(); ++i) {
    ret.res[i] *= factor;
  }
  return ret;
}

inline LatData operator*(const double factor, const LatData& ld)
{
  return ld * factor;
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
