#pragma once

#include <qlat-utils/crc32.h>
#include <qlat-utils/mpi-auto.h>
#include <qlat-utils/qar.h>
#include <qlat-utils/show.h>
#include <qlat-utils/timer.h>
#include <qlat-utils/utils-io.h>
#include <qlat-utils/utils-vec.h>
#include <qlat-utils/utils.h>
#include <stdint.h>
#include <zlib.h>

#include <cassert>
#include <string>
#include <vector>

namespace qlat
{  //

const std::string lat_data_header = "#!/usr/bin/env lat-io-glimpse\n";
// Recommended file extension ".lat"

const std::string lat_data_int_header = "#!/usr/bin/env lat-io-int-glimpse\n";
// Recommended file extension ".lati"

const std::string lat_data_long_header = "#!/usr/bin/env lat-io-long-glimpse\n";
// Recommended file extension ".latl"

const std::string lat_data_real_f_header = "#!/usr/bin/env lat-io-real-f-glimpse\n";
// Recommended file extension ".latf"

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

inline Long lat_info_size(const LatInfo& info, const Int level = 0)
{
  Long total = 1;
  for (Int i = level; i < (Int)info.size(); ++i) {
    total *= info[i].size;
  }
  return total;
}

template <class VecS>
LatDim lat_dim_string(const std::string& name, const VecS& indices)
{
  LatDim dim;
  dim.name = name;
  dim.size = indices.size();
  for (Int i = 0; i < dim.size; ++i) {
    dim.indices.push_back(indices[i]);
  }
  return dim;
}

Long lat_dim_idx(const LatDim& dim, const std::string& idx);

inline Long lat_dim_idx(const LatDim& dim, const Long& idx)
{
  Qassert((Long)dim.indices.size() <= dim.size);
  Qassert(0 <= idx and idx < dim.size);
  return idx;
}

template <class VecS>
Long lat_info_offset(const LatInfo& info, const VecS& idx)
// will return offset at the level the idx specify
// VecS can be std::vector<std::string> or std::vector<Long>
// or can be array of certain length
{
  Qassert((Long)idx.size() <= (Long)info.size());
  Long ret = 0;
  for (Int i = 0; i < (Int)idx.size(); ++i) {
    const Long k = lat_dim_idx(info[i], idx[i]);
    ret = ret * info[i].size + k;
  }
  return ret;
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

// ------------------------------------------

template <class T>
struct API LatDataT {
  LatInfo info;
  std::vector<T> res;
  //
  LatDataT() { init(); };
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
  void load_str(std::string& content)
  // Allow to destroy `content` to be more efficient.
  {
    QFile qfile =
        qfopen(QFileType::String, "/ load LatData /", QFileMode::Read, content);
    load(qfile);
    qfclose(qfile);
  }
  std::string save_str() const
  {
    QFile qfile =
        qfopen(QFileType::String, "/ save LatData /", QFileMode::Write);
    save(qfile);
    const std::string ret = qfile.content();
    qfclose(qfile);
    return ret;
  }
  //
  bool is_complex() const;
  Int ndim() const;
  T* data() { return res.data(); }
  const T* data() const { return res.data(); }
};

template <class T>
bool LatDataT<T>::is_complex() const
{
  return is_lat_info_complex(info);
}

template <class T>
Int LatDataT<T>::ndim() const
{
  if (is_lat_info_complex(info)) {
    return info.size() - 1;
  } else {
    return info.size();
  }
}

template <class T>
struct IsDataVectorType<LatDataT<T>> {
  using DataType = T;
  using BasicDataType = typename IsDataValueType<DataType>::BasicDataType;
  using ElementaryType = typename IsDataValueType<DataType>::ElementaryType;
  static constexpr bool value = is_data_value_type<DataType>();
};

// ------------------------------------------

template <class T>
bool is_initialized(const LatDataT<T>& ld)
{
  return ld.res.size() > 0;
}

template <class T>
void clear(LatDataT<T>& ld)
{
  ld.init();
}

template <class T>
Long lat_data_size(const LatDataT<T>& ld, const Int level = 0)
{
  return lat_info_size(ld.info, level);
}

template <class T>
void lat_data_alloc(LatDataT<T>& ld)
{
  ld.res.resize(lat_data_size(ld));
}

template <class T>
bool is_zero(const LatDataT<T>& ld)
{
  if (ld.info.size() > 0) {
    return false;
  }
  if (ld.res.size() > 1) {
    return false;
  }
  Qassert(ld.res.size() == 1);
  return ld.res[0] == 0;
}

template <class T>
Vector<T> get_data(const LatDataT<T>& ld)
{
  return get_data(ld.res);
}

template <class T, class VecS>
Long lat_data_offset(const LatDataT<T>& ld, const VecS& idx)
// will return offset at the level the idx specify
// VecS can be std::vector<std::string> or std::vector<Long>
// or can be array of certain length
{
  return lat_info_offset(ld.info, idx);
}

template <class T>
void set_zero(LatDataT<T>& ld)
{
  memset(ld.res.data(), 0, ld.res.size() * sizeof(T));
}

template <class T>
bool is_matching(const LatDataT<T>& ld1, const LatDataT<T>& ld2)
{
  return ld1.res.size() == ld2.res.size() and ld1.info == ld2.info;
}

template <class T, class VecS>
Vector<T> lat_data_get(LatDataT<T>& ld, const VecS& idx)
{
  const Long offset = lat_data_offset(ld, idx);
  const Long size = lat_data_size(ld, idx.size());
  Qassert(offset * size + size <= (Long)ld.res.size());
  Vector<T> ret(&ld.res[offset * size], size);
  return ret;
}

template <class T, class VecS>
Vector<T> lat_data_get_const(const LatDataT<T>& ld, const VecS& idx)
// Be cautious about the const property
// 改不改靠自觉
{
  const Long offset = lat_data_offset(ld, idx);
  const Long size = lat_data_size(ld, idx.size());
  Qassert(offset * size + size <= (Long)ld.res.size());
  Vector<T> ret(&ld.res[offset * size], size);
  return ret;
}

template <class T>
Vector<T> lat_data_get(LatDataT<T>& ld)
{
  array<Int, 0> idx;
  return lat_data_get(ld, idx);
}

template <class T>
Vector<T> lat_data_get_const(const LatDataT<T>& ld)
// Be cautious about the const property
// 改不改靠自觉
{
  array<Int, 0> idx;
  return lat_data_get_const(ld, idx);
}

template <class T>
Int glb_sum(LatDataT<T>& ld)
{
  return glb_sum_vec(get_data(ld.res));
}

template <class T>
Int bcast(LatDataT<T>& ld, const Int root)
{
  TIMER("bcast(ld,root)");
  if (1 == get_num_node()) {
    return 0;
  }
  Int ret = 0;
  std::string info_str;
  if (get_id_node() == root) {
    info_str = show(ld.info);
  }
  ret = bcast_val(info_str, root);
  if (ret != 0) {
    return ret;
  }
  if (get_id_node() == root) {
    Qassert((Long)ld.res.size() == lat_data_size(ld));
  } else {
    ld.info = read_lat_info(info_str);
    lat_data_alloc(ld);
  }
  ret = bcast_val(ld.res, root);
  return ret;
}

template <class T>
void lat_data_load_sync_node(LatDataT<T>& ld, const std::string& path)
{
  TIMER("lat_data_load_sync_node(ld)");
  ld.init();
  if (get_id_node() == 0) {
    ld.load(path);
  }
  Int ret = bcast(ld, 0);
  Qassert(ret == 0);
}

template <class T>
void lat_data_save_info(const std::string& path, const LatDataT<T>& ld)
{
  TIMER("lat_data_save_info");
  if (get_id_node() == 0) {
    ld.save(path);
  }
}

// ------------------------------------------

struct API LatDataInt : LatDataT<Int> {
};

template <>
void LatDataT<Int>::load(QFile& qfile);

template <>
void LatDataT<Int>::save(QFile& qfile) const;

inline LatDataInt lat_data_int_load_sync_node(const std::string& path)
{
  LatDataInt ld;
  lat_data_load_sync_node(ld, path);
  return ld;
}

template <>
struct IsDataVectorType<LatDataInt> {
  using DataType = Int;
  using BasicDataType = typename IsDataValueType<DataType>::BasicDataType;
  using ElementaryType = typename IsDataValueType<DataType>::ElementaryType;
  static constexpr bool value = is_data_value_type<DataType>();
};

// ------------------------------------------

struct API LatDataLong : LatDataT<Long> {
};

template <>
void LatDataT<Long>::load(QFile& qfile);

template <>
void LatDataT<Long>::save(QFile& qfile) const;

inline LatDataLong lat_data_long_load_sync_node(const std::string& path)
{
  LatDataLong ld;
  lat_data_load_sync_node(ld, path);
  return ld;
}

template <>
struct IsDataVectorType<LatDataLong> {
  using DataType = Long;
  using BasicDataType = typename IsDataValueType<DataType>::BasicDataType;
  using ElementaryType = typename IsDataValueType<DataType>::ElementaryType;
  static constexpr bool value = is_data_value_type<DataType>();
};

// ------------------------------------------

struct API LatData;
struct API LatDataRealF;

// ------------------------------------------

struct API LatData : LatDataT<RealD> {
  LatData() = default;
  LatData(const LatData&) = default;
  LatData(LatData&&) noexcept = default;
  LatData& operator=(const LatData&) = default;
  LatData& operator=(LatData&&) noexcept = default;
  LatData(const LatDataRealF& ld);
  LatData& operator=(const LatDataRealF& ld);
};

template <>
void LatDataT<RealD>::load(QFile& qfile);

template <>
void LatDataT<RealD>::save(QFile& qfile) const;

inline LatData lat_data_load_sync_node(const std::string& path)
{
  LatData ld;
  lat_data_load_sync_node(ld, path);
  return ld;
}

template <>
struct IsDataVectorType<LatData> {
  using DataType = RealD;
  using BasicDataType = typename IsDataValueType<DataType>::BasicDataType;
  using ElementaryType = typename IsDataValueType<DataType>::ElementaryType;
  static constexpr bool value = is_data_value_type<DataType>();
};

// ------------------------------------------

std::string show(const LatData& ld);

std::string show_real(const LatData& ld);

std::string show_complex(const LatData& ld);

void print(const LatData& ld);

const LatData& operator*=(LatData& ld, const RealD factor);

const LatData& operator*=(LatData& ld, const ComplexD& factor);

LatData operator*(const LatData& ld, const RealD factor);

LatData operator*(const RealD factor, const LatData& ld);

LatData operator*(const LatData& ld, const ComplexD& factor);

LatData operator*(const ComplexD& factor, const LatData& ld);

const LatData& operator+=(LatData& ld, const LatData& ld1);

const LatData& operator-=(LatData& ld, const LatData& ld1);

LatData operator+(const LatData& ld, const LatData& ld1);

LatData operator-(const LatData& ld, const LatData& ld1);

LatData operator-(const LatData& ld);

template <class VecS>
Vector<ComplexD> lat_data_cget(LatData& ld, const VecS& idx)
{
  Qassert(is_lat_info_complex(ld.info));
  Qassert((Long)idx.size() < (Long)ld.info.size());
  const Long offset = lat_data_offset(ld, idx);
  const Long size = lat_data_size(ld, idx.size());
  Qassert(size % 2 == 0);
  Qassert(offset * size + size <= (Long)ld.res.size());
  Vector<ComplexD> ret((ComplexD*)&ld.res[offset * size], size / 2);
  return ret;
}

template <class VecS>
Vector<ComplexD> lat_data_cget_const(const LatData& ld, const VecS& idx)
// Be cautious about the const property
// 改不改靠自觉
{
  Qassert(is_lat_info_complex(ld.info));
  Qassert((Long)idx.size() < (Long)ld.info.size());
  const Long offset = lat_data_offset(ld, idx);
  const Long size = lat_data_size(ld, idx.size());
  Qassert(size % 2 == 0);
  Qassert(offset * size + size <= (Long)ld.res.size());
  Vector<ComplexD> ret((ComplexD*)&ld.res[offset * size], size / 2);
  return ret;
}

inline Vector<ComplexD> lat_data_cget(LatData& ld)
{
  array<Int, 0> idx;
  return lat_data_cget(ld, idx);
}

inline Vector<ComplexD> lat_data_cget_const(const LatData& ld)
// Be cautious about the const property
// 改不改靠自觉
{
  array<Int, 0> idx;
  return lat_data_cget_const(ld, idx);
}

inline RealD qnorm(const LatData& ld) { return qnorm(ld.res); }

// ------------------------------------------

struct API LatDataRealF : LatDataT<RealF> {
  LatDataRealF() = default;
  LatDataRealF(const LatDataRealF&) = default;
  LatDataRealF(LatDataRealF&&) noexcept = default;
  LatDataRealF& operator=(const LatDataRealF&) = default;
  LatDataRealF& operator=(LatDataRealF&&) noexcept = default;
  LatDataRealF(const LatData& ld);
  LatDataRealF& operator=(const LatData& ld);
};

template <>
void LatDataT<RealF>::load(QFile& qfile);

template <>
void LatDataT<RealF>::save(QFile& qfile) const;

inline LatDataRealF lat_data_real_f_load_sync_node(const std::string& path)
{
  LatDataRealF ld;
  lat_data_load_sync_node(ld, path);
  return ld;
}

template <>
struct IsDataVectorType<LatDataRealF> {
  using DataType = RealF;
  using BasicDataType = typename IsDataValueType<DataType>::BasicDataType;
  using ElementaryType = typename IsDataValueType<DataType>::ElementaryType;
  static constexpr bool value = is_data_value_type<DataType>();
};

}  // namespace qlat
