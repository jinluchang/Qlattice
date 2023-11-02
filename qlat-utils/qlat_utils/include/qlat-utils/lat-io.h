#pragma once

#include <qlat-utils/timer.h>
#include <qlat-utils/crc32.h>
#include <qlat-utils/qutils-io.h>
#include <qlat-utils/qutils.h>
#include <qlat-utils/qutils-vec.h>
#include <qlat-utils/show.h>
#include <qlat-utils/qar.h>
#include <qlat-utils/qar-cache.h>
#include <stdint.h>
#include <zlib.h>

#include <cassert>
#include <string>
#include <vector>

namespace qlat
{  //

const std::string lat_data_header = "#!/usr/bin/env lat-io-glimpse\n";

struct API LatDim {
  std::string name;
  Long size;                         // size of this dimension
  std::vector<std::string> indices;  // indices names
                                     // (default: "-1", "-2", "-3", ...)
                                     // If indices.size() < size then example
                                     // will be indices[0], indices[1], "-3",
                                     // "-4", ... Won't check duplication (when
                                     // assess, search from left to right)
  //
  LatDim() { size = 0; }
};

typedef std::vector<LatDim> LatInfo;

bool is_lat_info_complex(const LatInfo& info);

LatDim lat_dim_re_im();

LatDim lat_dim_number(const std::string& name, const Long start, const Long end,
                      const Long inc = 1);

std::string show(const LatDim& dim);

std::string show(const LatInfo& info);

LatDim read_lat_dim(const std::string& str);

LatInfo read_lat_info(const std::string& str);

struct API LatData {
  LatInfo info;
  std::vector<double> res;
  //
  LatData() { init(); };
  //
  void init()
  {
    clear(info);
    clear(res);
    res.resize(1);
    res[0] = 0.0;
  }
  //
  void load(QFile& qfile);
  void load(const std::string& fn)
  {
    QFile qfile = qfopen(fn, "r");
    load(qfile);
    qfclose(qfile);
  }
  //
  void save(QFile& qfile) const;
  void save(const std::string& fn) const
  {
    QFile qfile = qfopen(fn + ".partial", "w");
    save(qfile);
    qfclose(qfile);
    qrename(fn + ".partial", fn);
  };
  //
  bool is_complex() const;
  int ndim() const;
  double* data() { return res.data(); }
  const double* data() const { return res.data(); }
};

inline bool LatData::is_complex() const { return is_lat_info_complex(info); }

inline int LatData::ndim() const
{
  if (is_lat_info_complex(info)) {
    return info.size() - 1;
  } else {
    return info.size();
  }
}

inline bool is_initialized(const LatData& ld) { return ld.res.size() > 0; }

inline bool is_zero(const LatData& ld)
{
  if (ld.info.size() > 0) {
    return false;
  }
  if (ld.res.size() > 1) {
    return false;
  }
  qassert(ld.res.size() == 1);
  return ld.res[0] == 0.0;
}

inline Long lat_info_size(const LatInfo& info, const int level = 0)
{
  Long total = 1;
  for (int i = level; i < (int)info.size(); ++i) {
    total *= info[i].size;
  }
  return total;
}

inline Long lat_data_size(const LatData& ld, const int level = 0)
{
  return lat_info_size(ld.info, level);
}

inline void lat_data_alloc(LatData& ld)
{
  ld.res.resize(lat_data_size(ld));
}

inline Vector<double> get_data(const LatData& ld) { return get_data(ld.res); }

inline void clear(LatData& ld) { ld.init(); }

std::string show_double(const LatData& ld);

std::string show_complex(const LatData& ld);

std::string show(const LatData& ld);

void print(const LatData& ld);

const LatData& operator*=(LatData& ld, const double factor);

const LatData& operator*=(LatData& ld, const ComplexD& factor);

LatData operator*(const LatData& ld, const double factor);

LatData operator*(const double factor, const LatData& ld);

LatData operator*(const LatData& ld, const ComplexD& factor);

LatData operator*(const ComplexD& factor, const LatData& ld);

const LatData& operator+=(LatData& ld, const LatData& ld1);

const LatData& operator-=(LatData& ld, const LatData& ld1);

LatData operator+(const LatData& ld, const LatData& ld1);

LatData operator-(const LatData& ld, const LatData& ld1);

LatData operator-(const LatData& ld);

template <class VecS>
LatDim lat_dim_string(const std::string& name, const VecS& indices)
{
  LatDim dim;
  dim.name = name;
  dim.size = indices.size();
  for (int i = 0; i < dim.size; ++i) {
    dim.indices.push_back(indices[i]);
  }
  return dim;
}

inline Long lat_dim_idx(const LatDim& dim, const std::string& idx)
{
  if ((Long)dim.indices.size() == 0) {
    Long i = read_long(idx);
    if (i >= 0) {
      qassert(i < dim.size);
      return i;
    } else {
      i = -i - 1;
      qassert(i < dim.size);
      return i;
    }
  } else if ((Long)dim.indices.size() == dim.size) {
    for (Long i = 0; i < (Long)dim.indices.size(); ++i) {
      if (idx == dim.indices[i]) {
        return i;
      }
    }
    std::vector<Long> possible_i_vec;
    for (Long i = 0; i < (Long)dim.indices.size(); ++i) {
      if (dim.indices[i].find(idx) != std::string::npos) {
        possible_i_vec.push_back(i);
      }
    }
    qassert(possible_i_vec.size() == 1);
    return possible_i_vec[0];
  } else {
    qassert((Long)dim.indices.size() <= dim.size);
    for (Long i = 0; i < (Long)dim.indices.size(); ++i) {
      if (idx == dim.indices[i]) {
        return i;
      }
    }
    const Long i = -read_long(idx) - 1;
    qassert((Long)dim.indices.size() <= i and i < dim.size);
    return i;
  }
}

inline Long lat_dim_idx(const LatDim& dim, const Long& idx)
{
  qassert((Long)dim.indices.size() <= dim.size);
  qassert(0 <= idx and idx < dim.size);
  return idx;
}

template <class VecS>
Long lat_info_offset(const LatInfo& info, const VecS& idx)
// will return offset at the level the idx specify
// VecS can be std::vector<std::string> or std::vector<Long>
// or can be array of certain length
{
  qassert((Long)idx.size() <= (Long)info.size());
  Long ret = 0;
  for (int i = 0; i < (int)idx.size(); ++i) {
    const Long k = lat_dim_idx(info[i], idx[i]);
    ret = ret * info[i].size + k;
  }
  return ret;
}

template <class VecS>
Long lat_data_offset(const LatData& ld, const VecS& idx)
// will return offset at the level the idx specify
// VecS can be std::vector<std::string> or std::vector<Long>
// or can be array of certain length
{
  return lat_info_offset(ld.info, idx);
}

inline std::string idx_name(const LatDim& dim, const Long idx)
{
  if (0 == (Long)dim.indices.size()) {
    return show(idx);
  } else if (idx < (Long)dim.indices.size()) {
    return dim.indices[idx];
  } else {
    return show(-idx - 1);
  }
}

inline bool operator==(const LatDim& d1, const LatDim& d2)
{
  // return d1.name == d2.name and d1.size == d2.size and d1.indices ==
  // d2.indices;
  return d1.name == d2.name and d1.size == d2.size;
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

template <class VecS>
Vector<double> lat_data_get(LatData& ld, const VecS& idx)
{
  const Long offset = lat_data_offset(ld, idx);
  const Long size = lat_data_size(ld, idx.size());
  qassert(offset * size + size <= (Long)ld.res.size());
  Vector<double> ret(&ld.res[offset * size], size);
  return ret;
}

template <class VecS>
Vector<double> lat_data_get_const(const LatData& ld, const VecS& idx)
// Be cautious about the const property
// 改不改靠自觉
{
  const Long offset = lat_data_offset(ld, idx);
  const Long size = lat_data_size(ld, idx.size());
  qassert(offset * size + size <= (Long)ld.res.size());
  Vector<double> ret(&ld.res[offset * size], size);
  return ret;
}

template <class VecS>
Vector<ComplexD> lat_data_complex_get(LatData& ld, const VecS& idx)
{
  qassert(is_lat_info_complex(ld.info));
  qassert((Long)idx.size() < (Long)ld.info.size());
  const Long offset = lat_data_offset(ld, idx);
  const Long size = lat_data_size(ld, idx.size());
  qassert(size % 2 == 0);
  qassert(offset * size + size <= (Long)ld.res.size());
  Vector<ComplexD> ret((ComplexD*)&ld.res[offset * size], size / 2);
  return ret;
}

template <class VecS>
Vector<ComplexD> lat_data_complex_get_const(const LatData& ld, const VecS& idx)
// Be cautious about the const property
// 改不改靠自觉
{
  qassert(is_lat_info_complex(ld.info));
  qassert((Long)idx.size() < (Long)ld.info.size());
  const Long offset = lat_data_offset(ld, idx);
  const Long size = lat_data_size(ld, idx.size());
  qassert(size % 2 == 0);
  qassert(offset * size + size <= (Long)ld.res.size());
  Vector<ComplexD> ret((ComplexD*)&ld.res[offset * size], size / 2);
  return ret;
}

template <class VecS>
Vector<ComplexD> lat_data_cget(LatData& ld, const VecS& idx)
{
  return lat_data_complex_get(ld, idx);
}

template <class VecS>
Vector<ComplexD> lat_data_cget_const(const LatData& ld, const VecS& idx)
// Be cautious about the const property
// 改不改靠自觉
{
  return lat_data_complex_get_const(ld, idx);
}

inline Vector<double> lat_data_get(LatData& ld)
{
  array<int, 0> idx;
  return lat_data_get(ld, idx);
}

inline Vector<double> lat_data_get_const(const LatData& ld)
// Be cautious about the const property
// 改不改靠自觉
{
  array<int, 0> idx;
  return lat_data_get_const(ld, idx);
}

inline Vector<ComplexD> lat_data_cget(LatData& ld)
{
  array<int, 0> idx;
  return lat_data_cget(ld, idx);
}

inline Vector<ComplexD> lat_data_cget_const(const LatData& ld)
// Be cautious about the const property
// 改不改靠自觉
{
  array<int, 0> idx;
  return lat_data_cget_const(ld, idx);
}

inline double qnorm(const LatData& ld) { return qnorm(ld.res); }

// -------------------

inline void lat_data_save_info(const std::string& path, const LatData& ld)
{
  TIMER("lat_data_save_info");
  if (get_id_node() == 0) {
    ld.save(path);
  }
}

}  // namespace qlat
