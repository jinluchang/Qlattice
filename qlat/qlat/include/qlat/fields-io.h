// vim: set ts=2 sw=2 expandtab:

#pragma once

// Copyright (c) 2016 Luchang Jin
// All rights reserved.

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// File format should be compatible with Christoph Lehner's file format.

// Add fields index, which stores fn_list, and offsets_map for each Fields file.
// index.qar
// Names will be "0", "1", "2", ...
// Content will be
/*
LENGTH_OF_FN
[newline character]
FN
[newline character]
DENSE_OR_SPARSE
[newline character]
OFFSET_START_0 OFFSET_START_1 OFFSET_START_2 ... OFFSET_START_N
[newline character]
OFFSET_STOP_0 OFFSET_STOP_1 OFFSET_STOP_2 ... OFFSET_STOP_N
[newline character]
*/
// (DENSE_OR_SPARSE will be "dense" or "sparse", not including quote.)

#include <errno.h>
#include <qlat-utils/qar.h>
#include <qlat/field-io.h>
#include <qlat/selected-field.h>
#include <qlat/selected-points.h>

namespace qlat
{  //

// ---------------------------------------------------

struct API BitSet {
  std::vector<unsigned char> bytes;
  size_t N;   // number of uncompressed elements
  size_t cN;  // number of compressed elements
  //
  BitSet()
  {
    N = 0;
    cN = 0;
  }
  BitSet(const size_t N_)
  {
    N = N_;
    cN = 0;
    clear(bytes);
    bytes.resize(1 + (N - 1) / 8, 0);
  }
  //
  void set(const void* pbytes, const size_t nbytes);
  //
  void set(const size_t idx, const bool v)
  {
    qassert(idx < N);
    if (v) {
      if (!get(idx)) {
        cN++;
      }
      bytes[idx / 8] |= 1 << (idx % 8);
    } else {
      if (get(idx)) {
        qassert(cN > 0);
        cN--;
      }
      bytes[idx / 8] &= ~(1 << (idx % 8));
    }
  }
  //
  bool get(const size_t idx) const
  {
    qassert(idx < N);
    return (bytes[idx / 8] & (1 << (idx % 8))) != 0;
  }
  //
  void set_f_rank(FieldRank& f_rank, const int64_t rank = 0);
  //
  bool check_f_rank(const FieldRank& f_rank);
  //
  void compress(const void* src, void* dst, size_t block_size) const;
  //
  void decompress(const void* src, void* dst, size_t block_size) const;
  //
  template <class M>
  std::vector<char> compress(const Vector<M>& src) const;
  //
  template <class M>
  std::vector<char> compress_selected(const Vector<M>& src) const;
};

template <class M>
std::vector<char> BitSet::compress(const Vector<M>& src) const
{
  qassert(src.size() % N == 0);
  size_t sz_block = sizeof(M) * src.size() / N;
  size_t sz_compressed = sz_block * cN;
  std::vector<char> dst(sz_compressed + bytes.size(), 0);
  memcpy(&dst[0], &bytes[0], bytes.size());
  compress((void*)src.data(), (void*)&dst[bytes.size()], sz_block);
  return dst;
}

template <class M>
std::vector<char> BitSet::compress_selected(const Vector<M>& src) const
{
  size_t sz_compressed = src.size() * sizeof(M);
  std::vector<char> dst(sz_compressed + bytes.size(), 0);
  memcpy(&dst[0], &bytes[0], bytes.size());
  if (src.size() == 0) {
    qassert(cN == 0);
  } else {
    qassert(src.size() % cN == 0);
    memcpy((void*)&dst[bytes.size()], src.data(), sz_compressed);
  }
  return dst;
}

std::vector<char> bitset_decompress(const std::vector<char>& data,
                                    const Long local_volume);

BitSet mk_bitset_from_field_rank(const FieldRank& f_rank);

// ---------------------------------------------------

struct API FieldsSegmentInfo {
  Long offset_start;
  Long offset_stop;
  bool is_sparse_field;
  //
  void init()
  {
    offset_start = 0;
    offset_stop = 0;
    is_sparse_field = false;
  }
  void init(const Long offset_start_, const Long offset_stop_,
            const bool is_sparse_field_)
  {
    offset_start = offset_start_;
    offset_stop = offset_stop_;
    is_sparse_field = is_sparse_field_;
  }
  //
  FieldsSegmentInfo() { init(); }
  FieldsSegmentInfo(const Long offset_start_, const Long offset_stop_,
                    const bool is_sparse_field_)
  {
    init(offset_start_, offset_stop_, is_sparse_field_);
  }
};

struct API FieldsWriter {
  //
  // should only use ShuffledFieldsWriter
  //
  std::string path;
  GeometryNode geon;
  QFile qfile;
  bool is_little_endian;  // should be true
  //
  std::vector<std::string> fn_list;
  std::map<std::string, FieldsSegmentInfo> offsets_map;
  Long max_offset;
  //
  FieldsWriter() { init(); }
  //
  void init();
  void init(const std::string& path_, const GeometryNode& geon_,
            const bool is_append = false);
  //
  void close() { qfclose(qfile); }
};

struct API FieldsReader {
  //
  // should only use ShuffledFieldsReader
  //
  std::string path;
  GeometryNode geon;
  QFile qfile;
  bool is_little_endian;  // should be true
  //
  bool is_read_through;
  std::vector<std::string> fn_list;
  std::map<std::string, FieldsSegmentInfo> offsets_map;
  Long max_offset;
  //
  Long file_size;
  //
  FieldsReader() { init(); }
  //
  void init();
  void init(const std::string& path_, const GeometryNode& geon_);
  //
  void close() { qfclose(qfile); }
};

void fields_writer_dirs_geon_info(const GeometryNode& geon,
                                  const std::string& path,
                                  const mode_t mode = default_dir_mode());

Coordinate shuffled_fields_reader_size_node_info(const std::string& path);

void mkfile(FieldsReader& fr);

std::string get_file_path(FieldsReader& fr);

template <class M>
void convert_endian(Vector<M> data, const bool is_little_endian)
{
  if (is_little_endian) {
    to_from_little_endian(data);
  } else {
    to_from_big_endian(data);
  }
}

Long qfwrite_convert_endian(void* ptr, const size_t size, const size_t nmemb,
                            QFile& qfile, const bool is_little_endian);

Long write(FieldsWriter& fw, const std::string& fn, const Geometry& geo,
           const Vector<char> data, const bool is_sparse_field = false);

Long qfread_convert_endian(void* ptr, const size_t size, const size_t nmemb,
                           QFile& qfile, const bool is_little_endian);

bool read_tag(FieldsReader& fr, std::string& fn, Coordinate& total_site,
              crc32_t& crc, int64_t& data_len, bool& is_sparse_field);

Long read_data(FieldsReader& fr, std::vector<char>& data,
               const int64_t data_len, const crc32_t crc);

Long read_next(FieldsReader& fr, std::string& fn, Coordinate& total_site,
               std::vector<char>& data, bool& is_sparse_field);

Long read_skip_next(FieldsReader& fr, std::string& fn);

void read_through(FieldsReader& fr);

bool does_file_exist(const FieldsReader& fr, const std::string& fn);

bool does_file_exist(const FieldsWriter& fr, const std::string& fn);

Long read(FieldsReader& fr, const std::string& fn, Coordinate& total_site,
          std::vector<char>& data, bool& is_sparse_field);

Long read_skip(FieldsReader& fr, const std::string& fn);

Long check_file(FieldsReader& fr, const std::string& fn, const bool is_check_data);

Int flush(FieldsWriter& fw);

// -----------------

template <class M>
Long write(FieldsWriter& fw, const std::string& fn, const Field<M>& field)
// field already have endianess converted correctly
{
  TIMER_FLOPS("write(fw,fn,field)");
  const Geometry& geo = field.geo();
  const Vector<M> v = get_data(field);
  const Vector<char> data((const char*)v.data(), v.data_size());
  const Long total_bytes = write(fw, fn, geo, data, false);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
Long write(FieldsWriter& fw, const std::string& fn, const SelectedField<M>& sf,
           const BitSet& bs)
// field already have endianness converted correctly
{
  TIMER_FLOPS("write(fw,fn,sf,bs)");
  const Geometry& geo = sf.geo();
  const std::vector<char> data = bs.compress_selected(get_data(sf));
  const Long total_bytes = write(fw, fn, geo, get_data(data), true);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
void set_field_from_data(Field<M>& field, const GeometryNode& geon,
                         const Coordinate& total_site,
                         const std::vector<char>& data,
                         const bool is_sparse_field)
{
  TIMER("set_field_from_data(field,geon,total_site,data)");
  const Coordinate node_site = total_site / geon.size_node;
  const Long local_volume = product(node_site);
  ConstHandle<std::vector<char>> hdata(data);
  std::vector<char> dc_data;
  if (is_sparse_field) {
    dc_data = bitset_decompress(data, local_volume);
    hdata.init(dc_data);
  }
  if (hdata().size() == 0) {
    field.init();
    return;
  }
  const Long local_data_size = hdata().size();
  const Long site_data_size = local_data_size / local_volume;
  qassert(site_data_size % sizeof(M) == 0);
  const Int multiplicity = site_data_size / sizeof(M);
  Geometry geo;
  geo.init(geon, node_site);
  field.init();
  field.init(geo, multiplicity);
  Vector<M> fv = get_data(field);
  qassert(fv.data_size() == (Long)hdata().size());
  memcpy(fv.data(), hdata().data(), fv.data_size());
}

template <class M>
void set_field_from_data(SelectedField<M>& sf, FieldRank& f_rank,
                         const std::vector<char>& data)
// obtain f_rank from data
{
  TIMER("set_field_from_data(sf,f_rank,data)");
  const Geometry& geo = f_rank.geo();
  const Coordinate& node_site = geo.node_site;
  const Long local_volume = product(node_site);
  const size_t N = local_volume;
  const size_t nbytes = 1 + (N - 1) / 8;
  BitSet bs(N);
  bs.set(&data[0], nbytes);
  bs.set_f_rank(f_rank);
  const Long n_elems = bs.cN;
  if (n_elems == 0) {
    sf.init();
    return;
  }
  const size_t sz_compressed = data.size() - nbytes;
  qassert(sz_compressed % n_elems == 0);
  const size_t sz_block = sz_compressed / n_elems;
  qassert(sz_block % sizeof(M) == 0);
  const Int multiplicity = sz_block / sizeof(M);
  const Vector<M> fdata((const M*)&data[nbytes], n_elems * multiplicity);
  sf.init(geo, n_elems, multiplicity);
  assign(get_data(sf), fdata);
}

template <class M>
void set_field_from_data(SelectedField<M>& sf, const std::vector<char>& data,
                         const FieldSelection& fsel)
// fsel must be a subset of the actual data
{
  TIMER("set_field_from_data(sf,data,fsel)");
  const Geometry& geo = fsel.f_rank.geo();
  FieldSelection fsel_data;
  fsel_data.f_rank.init(geo);
  SelectedField<M> sf_data;
  set_field_from_data(sf_data, fsel_data.f_rank, data);
  update_field_selection(fsel_data);
  qassert(is_containing(fsel_data, fsel));
  if (fsel.n_elems == 0) {
    sf.init();
  } else {
    set_selected_field(sf, sf_data, fsel, fsel_data);
  }
}

template <class M>
Long read(FieldsReader& fr, const std::string& fn, Field<M>& field)
// field endianess not converted at all
{
  TIMER_FLOPS("read(fr,fn,field)");
  Coordinate total_site;
  std::vector<char> data;
  bool is_sparse_field = false;
  const Long total_bytes = read(fr, fn, total_site, data, is_sparse_field);
  if (0 == total_bytes) {
    return 0;
  }
  set_field_from_data(field, fr.geon, total_site, data, is_sparse_field);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
Long read(FieldsReader& fr, const std::string& fn, SelectedField<M>& sf,
          FieldRank& f_rank)
// field endianess not converted at all
// f_rank does not need to be initialized
{
  TIMER_FLOPS("read(fr,fn,sf,f_rank)");
  Coordinate total_site;
  std::vector<char> data;
  bool is_sparse_field = false;
  const Long total_bytes = read(fr, fn, total_site, data, is_sparse_field);
  if (0 == total_bytes) {
    return 0;
  }
  if (not is_sparse_field) {
    return 0;
  }
  Geometry geo;
  geo.init(fr.geon, total_site / fr.geon.size_node);
  f_rank.init(geo, 1);
  qassert(f_rank.geo().is_only_local);
  set_field_from_data(sf, f_rank, data);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
Long read(FieldsReader& fr, const std::string& fn, const FieldSelection& fsel,
          SelectedField<M>& sf)
// field endianess not converted at all
// fsel must match the actual data
// (code will verify & will fail if not match)
{
  TIMER_FLOPS("read(fr,fn,fsel,sf)");
  Coordinate total_site;
  std::vector<char> data;
  bool is_sparse_field = false;
  const Long total_bytes = read(fr, fn, total_site, data, is_sparse_field);
  if (0 == total_bytes) {
    return 0;
  }
  if (not is_sparse_field) {
    return 0;
  }
  qassert(total_site == fsel.f_rank.geo().total_site());
  set_field_from_data(sf, data, fsel);
  timer.flops += total_bytes;
  return total_bytes;
}

// -----------------

struct API ShuffledBitSet {
  FieldSelection fsel;
  ShufflePlan sp;
  std::vector<FieldSelection> fsels;
  std::vector<BitSet> vbs;
};

ShuffledBitSet mk_shuffled_bitset(const FieldRank& f_rank,
                                  const Coordinate& new_size_node);

ShuffledBitSet mk_shuffled_bitset(const FieldSelection& fsel,
                                  const Coordinate& new_size_node);

ShuffledBitSet mk_shuffled_bitset(const Coordinate& total_site,
                                  const PointsSelection& psel,
                                  const Coordinate& new_size_node);

ShuffledBitSet mk_shuffled_bitset(const FieldRank& f_rank,
                                  const PointsSelection& psel,
                                  const Coordinate& new_size_node);

// -----------------

struct API ShuffledFieldsWriter {
  std::string path;
  Coordinate new_size_node;
  std::vector<FieldsWriter> fws;
  QarFile qar_index;
  Long qar_index_idx;
  //
  ShuffledFieldsWriter() { init(); }
  ShuffledFieldsWriter(const std::string& path_,
                       const Coordinate& new_size_node_,
                       const bool is_append = false,
                       const bool is_removing_old = false)
  // interface function
  {
    init(path_, new_size_node_, is_append, is_removing_old);
  }
  //
  ~ShuffledFieldsWriter() { close(); }
  //
  void init();
  void init(const std::string& path_, const Coordinate& new_size_node_,
            const bool is_append = false, const bool is_removing_old = false);
  //
  void close();
};

struct API ShuffledFieldsReader {
  std::string path;
  Coordinate new_size_node;
  std::vector<FieldsReader> frs;
  //
  ShuffledFieldsReader()
  // interface function
  {
    init();
  }
  ShuffledFieldsReader(const std::string& path_,
                       const Coordinate& new_size_node_ = Coordinate())
  // interface function
  {
    init(path_, new_size_node_);
  }
  //
  ~ShuffledFieldsReader() { close(); }
  //
  void init();
  void init(const std::string& path_,
            const Coordinate& new_size_node_ = Coordinate());
  //
  void close();
};

// --------------------------------------------

std::string show(const ShuffledFieldsWriter& sfw);

using ShuffledFieldsWriterMap = std::map<Long, Handle<ShuffledFieldsWriter>>;

API inline ShuffledFieldsWriterMap& get_all_shuffled_fields_writer()
{
  static ShuffledFieldsWriterMap sfwm;
  return sfwm;
}

std::vector<std::string> show_all_shuffled_fields_writer();

void add_shuffled_fields_writer(ShuffledFieldsWriter& sfw);

void remove_shuffled_fields_writer(ShuffledFieldsWriter& sfw);

void close_all_shuffled_fields_writer();

// --------------------------------------------

std::string show(const ShuffledFieldsReader& sfr);

using ShuffledFieldsReaderMap = std::map<Long, Handle<ShuffledFieldsReader>>;

API inline ShuffledFieldsReaderMap& get_all_shuffled_fields_reader()
{
  static ShuffledFieldsReaderMap sfrm;
  return sfrm;
}

std::vector<std::string> show_all_shuffled_fields_reader();

void add_shuffled_fields_reader(ShuffledFieldsReader& sfr);

void remove_shuffled_fields_reader(ShuffledFieldsReader& sfr);

void close_all_shuffled_fields_reader();

// --------------------------------------------

typedef Cache<std::string, ShuffledFieldsReader> ShuffledFieldsReaderCache;

typedef Cache<std::string, ShuffledFieldsWriter> ShuffledFieldsWriterCache;

API inline ShuffledFieldsReaderCache& get_shuffled_fields_reader_cache()
{
  static ShuffledFieldsReaderCache cache("ShuffledFieldsReaderCache", 4, 1);
  return cache;
}

API inline ShuffledFieldsWriterCache& get_shuffled_fields_writer_cache()
{
  static ShuffledFieldsWriterCache cache("ShuffledFieldsWriterCache", 4, 1);
  return cache;
}

ShuffledFieldsReader& get_shuffled_fields_reader(
    const std::string& path, const Coordinate& new_size_node = Coordinate());

ShuffledFieldsWriter& get_shuffled_fields_writer(
    const std::string& path,
    const Coordinate& new_size_node = Coordinate(2, 2, 2, 4),
    const bool is_append = false, const bool is_removing_old = false);

inline void close_shuffled_fields_reader(const std::string& path)
{
  TIMER("close_all_shuffled_fields_reader");
  get_shuffled_fields_reader_cache().erase(path);
}

inline void close_shuffled_fields_writer(const std::string& path)
{
  TIMER("close_all_shuffled_fields_writer");
  get_shuffled_fields_writer_cache().erase(path);
}

inline void clear_shuffled_fields_reader_cache()
{
  TIMER_VERBOSE("clear_shuffled_fields_reader_cache");
  get_shuffled_fields_reader_cache().clear();
}

inline void clear_shuffled_fields_writer_cache()
{
  TIMER_VERBOSE("clear_shuffled_fields_writer_cache");
  get_shuffled_fields_writer_cache().clear();
}

// --------------------------------------------

Long flush(ShuffledFieldsWriter& sfw);

void read_through_sync_node(ShuffledFieldsReader& sfr);

bool does_file_exist_sync_node(const ShuffledFieldsReader& sfr,
                               const std::string& fn);

bool does_file_exist_sync_node(const ShuffledFieldsWriter& sfw,
                               const std::string& fn);

bool is_sparse_field_sync_node(const ShuffledFieldsReader& sfr, const std::string& fn);

bool check_file_sync_node(ShuffledFieldsReader& sfr, const std::string& fn,
                          const bool is_check_data,
                          std::vector<Long>& final_offsets);

std::vector<std::string> list_fields(const ShuffledFieldsReader& sfr, bool is_skipping_check = false);

std::vector<std::string> list_fields(const ShuffledFieldsWriter& sfw, bool is_skipping_check = false);

Int truncate_fields_sync_node(const std::string& path,
                              const std::vector<std::string>& fns_keep,
                              const Coordinate& new_size_node = Coordinate());

void properly_truncate_fields_sync_node(
    std::vector<std::string>& fn_list,
    std::vector<std::vector<FieldsSegmentInfo>>& offsets_list, const std::string& path,
    const bool is_check_all, const bool is_only_check,
    const Coordinate& new_size_node);

std::vector<std::string> properly_truncate_fields_sync_node(
    const std::string& path, const bool is_check_all = false,
    const bool is_only_check = false,
    const Coordinate& new_size_node = Coordinate());

bool has_duplicates(const ShuffledFieldsReader& sfr);

// -----------------

std::string show_field_index(const std::string& fn,
                             const std::vector<FieldsSegmentInfo>& offsets);

void parse_field_index(std::string& fn, std::vector<FieldsSegmentInfo>& all_offsets,
                       const std::string& field_index_content);

std::vector<FieldsSegmentInfo> collect_fields_offsets(const ShuffledFieldsWriter& sfw,
                                         const std::string& fn);

std::vector<FieldsSegmentInfo> collect_fields_offsets(const ShuffledFieldsReader& sfr,
                                         const std::string& fn);

bool populate_fields_offsets(ShuffledFieldsReader& sfr, const std::string& fn,
                             const std::vector<FieldsSegmentInfo>& offsets);

void load_all_fields_index(std::vector<std::string>& fn_list,
                           std::vector<std::vector<FieldsSegmentInfo>>& all_offsets_list,
                           QarFile& qar_index);

void save_fields_index(QarFile& qar_index, const Long idx,
                       const std::string& fn, const std::vector<FieldsSegmentInfo>& offsets);

void save_fields_index(ShuffledFieldsWriter& sfw, const std::string& fn);

void fields_build_index(const ShuffledFieldsReader& sfr);

void fields_build_index(const std::string& path);

// -----------------

template <class M>
Long write(ShuffledFieldsWriter& sfw, const std::string& fn,
           const Field<M>& field)
// interface function
{
  TIMER_VERBOSE_FLOPS("write(sfw,fn,field)");
  displayln_info(0, fname + ssprintf(": writing field with fn='%s' from '%s'.",
                                     fn.c_str(), sfw.path.c_str()));
  std::vector<Field<M>> fs;
  shuffle_field(fs, field, sfw.new_size_node);
  qassert(fs.size() == sfw.fws.size());
  Long total_bytes = 0;
  const Int num_node = product(sfw.new_size_node);
  const Int n_cycle = std::max(1, num_node / dist_write_par_limit());
  for (Int cycle = 0; cycle < n_cycle; cycle++) {
    for (Int i = 0; i < (int)fs.size(); ++i) {
      const Int id_node = sfw.fws[i].geon.id_node;
      qassert(0 <= id_node && id_node < num_node);
      if (id_node % n_cycle == cycle) {
        total_bytes += write(sfw.fws[i], fn, fs[i]);
      }
    }
    SYNC_NODE();
  }
  glb_sum(total_bytes);
  save_fields_index(sfw, fn);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
Long write(ShuffledFieldsWriter& sfw, const std::string& fn,
           const ShuffledBitSet& sbs, const SelectedField<M>& sf)
// interface function
// sbs must match the actual data
{
  TIMER_VERBOSE_FLOPS("write(sfw,fn,sbs,sf)");
  displayln_info(
      0, fname + ssprintf(": writing sparse field with fn='%s' from '%s'.",
                          fn.c_str(), sfw.path.c_str()));
  std::vector<SelectedField<M>> sfs;
  shuffle_field(sfs, sf, sbs.sp);
  qassert(sfs.size() == sfw.fws.size());
  qassert(sbs.vbs.size() == sfw.fws.size());
  Long total_bytes = 0;
  const Int num_node = product(sfw.new_size_node);
  const Int n_cycle = std::max(1, num_node / dist_write_par_limit());
  for (Int cycle = 0; cycle < n_cycle; cycle++) {
    for (Int i = 0; i < (int)sfs.size(); ++i) {
      const Int id_node = sfw.fws[i].geon.id_node;
      qassert(0 <= id_node && id_node < num_node);
      if (id_node % n_cycle == cycle) {
        total_bytes += write(sfw.fws[i], fn, sfs[i], sbs.vbs[i]);
      }
    }
    SYNC_NODE();
  }
  glb_sum(total_bytes);
  save_fields_index(sfw, fn);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
void set_field_info_from_fields(Coordinate& total_site, Int& multiplicity,
                                std::vector<Field<M>>& fs,
                                const ShuffledFieldsReader& sfr)
{
  TIMER("set_field_info_from_fields(total_site,multiplicity,fs,sfr)");
  total_site = Coordinate();
  multiplicity = 0;
  std::vector<Long> available_nodes(product(sfr.new_size_node), 0);
  for (Int i = 0; i < (int)fs.size(); ++i) {
    const Int id_node = sfr.frs[i].geon.id_node;
    qassert(0 <= id_node and id_node < (int)available_nodes.size());
    if (fs[i].initialized) {
      available_nodes[id_node] = get_id_node() + 1;
    }
  }
  glb_sum(available_nodes);
  Int id_node_first_available = 0;
  Int id_node_bcast_from = 0;
  for (Int i = 0; i < (int)available_nodes.size(); ++i) {
    if (available_nodes[i] > 0) {
      id_node_first_available = i;
      id_node_bcast_from = available_nodes[i] - 1;
      break;
    }
  }
  for (Int i = 0; i < (int)fs.size(); ++i) {
    const Int id_node = sfr.frs[i].geon.id_node;
    if (id_node == id_node_first_available) {
      total_site = fs[i].geo().total_site();
      multiplicity = fs[i].multiplicity;
      qassert(get_id_node() == id_node_bcast_from);
    }
  }
  bcast(get_data_one_elem(total_site), id_node_bcast_from);
  bcast(get_data_one_elem(multiplicity), id_node_bcast_from);
  for (Int i = 0; i < (int)fs.size(); ++i) {
    if (not fs[i].initialized) {
      const GeometryNode& geon = sfr.frs[i].geon;
      const Coordinate node_site = total_site / geon.size_node;
      Geometry geo;
      geo.init(geon, node_site);
      fs[i].init(geo, multiplicity);
      set_zero(fs[i]);
    }
  }
}

template <class M>
void set_field_info_from_fields(Coordinate& total_site, Int& multiplicity,
                                std::vector<SelectedField<M>>& sfs,
                                const ShuffledFieldsReader& sfr)
{
  TIMER("set_field_info_from_fields(total_site,multiplicity,sfs,sfr)");
  total_site = Coordinate();
  multiplicity = 0;
  std::vector<Long> available_nodes(product(sfr.new_size_node), 0);
  for (Int i = 0; i < (int)sfs.size(); ++i) {
    const Int id_node = sfr.frs[i].geon.id_node;
    qassert(0 <= id_node and id_node < (int)available_nodes.size());
    if (is_initialized(sfs[i])) {
      available_nodes[id_node] = get_id_node() + 1;
    }
  }
  glb_sum(available_nodes);
  Int id_node_first_available = 0;
  Int id_node_bcast_from = 0;
  for (Int i = 0; i < (int)available_nodes.size(); ++i) {
    if (available_nodes[i] > 0) {
      id_node_first_available = i;
      id_node_bcast_from = available_nodes[i] - 1;
      break;
    }
  }
  for (Int i = 0; i < (int)sfs.size(); ++i) {
    const Int id_node = sfr.frs[i].geon.id_node;
    if (id_node == id_node_first_available) {
      qassert(is_initialized(sfs[i]));
      total_site = sfs[i].geo().total_site();
      multiplicity = sfs[i].multiplicity;
      qassert(get_id_node() == id_node_bcast_from);
    }
  }
  bcast(get_data_one_elem(total_site), id_node_bcast_from);
  bcast(get_data_one_elem(multiplicity), id_node_bcast_from);
  for (Int i = 0; i < (int)sfs.size(); ++i) {
    if (not sfs[i].initialized) {
      const GeometryNode& geon = sfr.frs[i].geon;
      const Coordinate node_site = total_site / geon.size_node;
      Geometry geo;
      geo.init(geon, node_site);
      sfs[i].init(geo, 0, multiplicity);
      set_zero(sfs[i]);
    }
  }
}

template <class M>
Long read(ShuffledFieldsReader& sfr, const std::string& fn, Field<M>& field)
// interface function
{
  TIMER_VERBOSE_FLOPS("read(sfr,fn,field)");
  Long total_bytes = 0;
  displayln_info(0, fname + ssprintf(": reading field with fn='%s' from '%s'.",
                                     fn.c_str(), sfr.path.c_str()));
  std::vector<Field<M>> fs(sfr.frs.size());
  Long zero_size_count = 0;
  const Int num_node = product(sfr.new_size_node);
  const Int n_cycle = std::max(1, num_node / dist_read_par_limit());
  for (Int cycle = 0; cycle < n_cycle; cycle++) {
    for (Int i = 0; i < (int)fs.size(); ++i) {
      const Int id_node = sfr.frs[i].geon.id_node;
      qassert(0 <= id_node && id_node < num_node);
      if (id_node % n_cycle == cycle) {
        const Long bytes = read(sfr.frs[i], fn, fs[i]);
        if (0 == bytes) {
          zero_size_count += 1;
          qassert(0 == total_bytes);
        } else {
          total_bytes += bytes;
        }
      }
    }
    SYNC_NODE();
  }
  glb_sum(total_bytes);
  if (0 != zero_size_count) {
    qassert(0 == total_bytes);
  }
  if (0 == total_bytes) {
    qwarn(fname + ssprintf(": total_bytes=%ld (fn='%s', sfr.path='%s')",
                           total_bytes, fn.c_str(), sfr.path.c_str()));
    return 0;
  }
  Coordinate total_site;
  Int multiplicity = 0;
  set_field_info_from_fields(total_site, multiplicity, fs, sfr);
  Geometry geo;
  geo.init(total_site);
  field.init(geo, multiplicity);
  shuffle_field_back(field, fs, sfr.new_size_node);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
Long read(ShuffledFieldsReader& sfr, const std::string& fn,
          SelectedField<M>& sf, FieldSelection& fsel)
// interface function
{
  TIMER_VERBOSE_FLOPS("read(sfr,fn,sf,fsel)")
  Long total_bytes = 0;
  displayln_info(0, fname + ssprintf(": reading field with fn='%s' from '%s'.",
                                     fn.c_str(), sfr.path.c_str()));
  std::vector<SelectedField<M>> sfs(sfr.frs.size());
  std::vector<Field<int64_t>> f_rank_s(sfr.frs.size());
  Long zero_size_count = 0;
  const Int num_node = product(sfr.new_size_node);
  const Int n_cycle = std::max(1, num_node / dist_read_par_limit());
  for (Int cycle = 0; cycle < n_cycle; cycle++) {
    for (Int i = 0; i < (int)sfs.size(); ++i) {
      const Int id_node = sfr.frs[i].geon.id_node;
      qassert(0 <= id_node && id_node < num_node);
      if (id_node % n_cycle == cycle) {
        FieldRank& f_rank = static_cast<FieldRank&>(f_rank_s[i]);
        const Long bytes = read(sfr.frs[i], fn, sfs[i], f_rank);
        if (0 == bytes) {
          zero_size_count += 1;
          qassert(0 == total_bytes);
        } else {
          total_bytes += bytes;
        }
      }
    }
    SYNC_NODE();
  }
  glb_sum(total_bytes);
  if (0 != zero_size_count) {
    qassert(0 == total_bytes);
  }
  if (0 == total_bytes) {
    qwarn(fname + ssprintf(": total_bytes=%ld (fn='%s', sfr.path='%s')",
                           total_bytes, fn.c_str(), sfr.path.c_str()));
    return 0;
  }
  Coordinate total_site;
  Int multiplicity = 0;
  set_field_info_from_fields(total_site, multiplicity, sfs, sfr);
  const Geometry geo(total_site);
  fsel.f_rank.init(geo, 1);
  shuffle_field_back(fsel.f_rank, f_rank_s, sfr.new_size_node);
  update_field_selection(fsel);
  sf.init(fsel, multiplicity);
  shuffle_field_back(sf, sfs, sfr.new_size_node, fsel);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
Long read(ShuffledFieldsReader& sfr, const std::string& fn,
          const ShuffledBitSet& sbs, SelectedField<M>& sf)
// interface function
// sbs must match (can be subset of) the actual data
// (code will verify & will fail if not match)
{
  TIMER_VERBOSE_FLOPS("read(sfr,fn,sbs,sf)");
  sf.init();
  Long total_bytes = 0;
  displayln_info(
      0, fname + ssprintf(": reading sparse field with fn='%s' from '%s'.",
                          fn.c_str(), sfr.path.c_str()));
  std::vector<SelectedField<M>> sfs(sfr.frs.size());
  Long zero_size_count = 0;
  const Int num_node = product(sfr.new_size_node);
  const Int n_cycle = std::max(1, num_node / dist_read_par_limit());
  for (Int cycle = 0; cycle < n_cycle; cycle++) {
    for (Int i = 0; i < (int)sfs.size(); ++i) {
      const Int id_node = sfr.frs[i].geon.id_node;
      qassert(0 <= id_node && id_node < num_node);
      if (id_node % n_cycle == cycle) {
        const Long bytes = read(sfr.frs[i], fn, sbs.fsels[i], sfs[i]);
        if (0 == bytes) {
          zero_size_count += 1;
          qassert(0 == total_bytes);
        } else {
          total_bytes += bytes;
        }
      }
    }
    SYNC_NODE();
  }
  glb_sum(total_bytes);
  if (0 != zero_size_count) {
    qassert(0 == total_bytes);
  }
  if (0 == total_bytes) {
    qwarn(fname + ssprintf(": total_bytes=%ld (fn='%s', sfr.path='%s')",
                           total_bytes, fn.c_str(), sfr.path.c_str()));
    return 0;
  }
  Coordinate total_site;
  Int multiplicity = 0;
  set_field_info_from_fields(total_site, multiplicity, sfs, sfr);
  qassert(total_site != Coordinate());
  qassert(multiplicity > 0);
  sf.init(sbs.fsel, multiplicity);
  shuffle_field_back(sf, sfs, sbs.sp);
  qassert(is_consistent(sf, sbs.fsel));
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
Long write_float_from_double(ShuffledFieldsWriter& sfw, const std::string& fn,
                             const Field<M>& field)
// interface function
{
  TIMER_VERBOSE_FLOPS("write_float_from_double(sfw,fn,field)");
  Field<float> ff;
  convert_field_float_from_double(ff, field);
  to_from_little_endian(get_data(ff));
  const Long total_bytes = write(sfw, fn, ff);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
Long write_float_from_double(ShuffledFieldsWriter& sfw, const std::string& fn,
                             const ShuffledBitSet& sbs,
                             const SelectedField<M>& sf)
// interface function
// sbs must match the actual data
{
  TIMER_VERBOSE_FLOPS("write_float_from_double(sfw,fn,sbs,sf)");
  SelectedField<float> sff;
  convert_field_float_from_double(sff, sf);
  to_from_little_endian(get_data(sff));
  const Long total_bytes = write(sfw, fn, sbs, sff);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
Long read_double_from_float(ShuffledFieldsReader& sfr, const std::string& fn,
                            Field<M>& field)
// interface function
{
  TIMER_VERBOSE_FLOPS("read_double_from_float(sfr,fn,field)");
  Field<float> ff;
  const Long total_bytes = read(sfr, fn, ff);
  if (total_bytes == 0) {
    return 0;
  } else {
    to_from_little_endian(get_data(ff));
    convert_field_double_from_float(field, ff);
    timer.flops += total_bytes;
    return total_bytes;
  }
}

template <class M>
Long read_double_from_float(ShuffledFieldsReader& sfr, const std::string& fn,
                            SelectedField<M>& sf, FieldSelection& fsel)
// interface function
{
  TIMER_VERBOSE_FLOPS("read_double_from_float(sfr,fn,sf,fsel)");
  sf.init();
  SelectedField<float> sff;
  const Long total_bytes = read(sfr, fn, sff, fsel);
  if (total_bytes == 0) {
    return 0;
  } else {
    to_from_little_endian(get_data(sff));
    convert_field_double_from_float(sf, sff);
    timer.flops += total_bytes;
    return total_bytes;
  }
}

template <class M>
Long read_double_from_float(ShuffledFieldsReader& sfr, const std::string& fn,
                            const ShuffledBitSet& sbs, SelectedField<M>& sf)
// interface function
// sbs must match the actual data
{
  TIMER_VERBOSE_FLOPS("read_double_from_float(sfr,fn,sbs,sf)");
  sf.init();
  SelectedField<float> sff;
  const Long total_bytes = read(sfr, fn, sbs, sff);
  if (total_bytes == 0) {
    return 0;
  } else {
    to_from_little_endian(get_data(sff));
    convert_field_double_from_float(sf, sff);
    timer.flops += total_bytes;
    return total_bytes;
  }
}

// -----------------

template <class M>
Long read_field(Field<M>& field, const std::string& path, const std::string& fn)
// interface function
{
  TIMER_VERBOSE("read_field(field,path,fn)");
  ShuffledFieldsReader& sfr = get_shuffled_fields_reader(path);
  return read(sfr, fn, field);
}

template <class M>
Long read(SelectedField<M>& sf, FieldSelection& fsel, const std::string& path,
          const std::string& fn)
// interface function
{
  TIMER_VERBOSE("read_field(sf,fsel,path,fn)");
  ShuffledFieldsReader& sfr = get_shuffled_fields_reader(path);
  return read(sfr, fn, sf, fsel);
}

template <class M>
Long read_field(SelectedField<M>& sf, const std::string& path,
                const std::string& fn, const ShuffledBitSet& sbs)
// interface function
// sbs must match (can be subset of) the actual data
// (code will verify & will fail if not match)
{
  TIMER_VERBOSE("read_field(sf,path,fn,sbs)");
  ShuffledFieldsReader& sfr = get_shuffled_fields_reader(path);
  return read(sfr, fn, sbs, sf);
}

template <class M>
Long read_field_double_from_float(Field<M>& field, const std::string& path,
                                  const std::string& fn)
// interface function
{
  TIMER_VERBOSE("read_field_double_from_float(field,path,fn)");
  ShuffledFieldsReader& sfr = get_shuffled_fields_reader(path);
  return read_double_from_float(sfr, fn, field);
}

template <class M>
Long read_field_double_from_float(SelectedField<M>& sf, const std::string& path,
                                  const std::string& fn,
                                  const ShuffledBitSet& sbs)
// interface function
// sbs must match the actual data
{
  TIMER_VERBOSE("read_field_double_from_float(sf,path,fn,sbs)");
  ShuffledFieldsReader& sfr = get_shuffled_fields_reader(path);
  return read_double_from_float(sfr, fn, sbs, sf);
}

template <class M>
Long write_field(const Field<M>& field, const std::string& path,
                 const std::string& fn)
// interface function
{
  TIMER_VERBOSE_FLOPS("write_field(field,path,fn)");
  ShuffledFieldsWriter& sfw = get_shuffled_fields_writer(path);
  const Long total_bytes = write(sfw, fn, field);
  flush(sfw);
  return total_bytes;
}

template <class M>
Long write_field(const SelectedField<M>& sf, const ShuffledBitSet& sbs,
                 const std::string& path, const std::string& fn)
// interface function
{
  TIMER_VERBOSE_FLOPS("write_field(sf,sbs,path,fn,is_append)");
  ShuffledFieldsWriter& sfw = get_shuffled_fields_writer(path);
  const Long total_bytes = write(sfw, fn, sbs, sf);
  flush(sfw);
  return total_bytes;
}

template <class M>
Long write_field_float_from_double(const Field<M>& field,
                                   const std::string& path,
                                   const std::string& fn)
// interface function
{
  TIMER_VERBOSE_FLOPS("write_field(field,path,fn)");
  ShuffledFieldsWriter& sfw = get_shuffled_fields_writer(path);
  const Long total_bytes = write_float_from_double(sfw, fn, field);
  flush(sfw);
  return total_bytes;
}

template <class M>
Long write_field_float_from_double(const SelectedField<M>& sf,
                                   const ShuffledBitSet& sbs,
                                   const std::string& path,
                                   const std::string& fn)
// interface function
{
  TIMER_VERBOSE_FLOPS("write_field(sf,sbs,path,fn,is_append)");
  ShuffledFieldsWriter& sfw = get_shuffled_fields_writer(path);
  const Long total_bytes = write_float_from_double(sfw, fn, sbs, sf);
  flush(sfw);
  return total_bytes;
}

// -----------------

std::vector<std::string> list_fields(const std::string& path);

bool does_file_exist_sync_node(const std::string& path, const std::string& fn);

// --------------------

#ifdef QLAT_INSTANTIATE_FIELDS_IO
#define QLAT_EXTERN
#else
#define QLAT_EXTERN extern
#endif

#define QLAT_EXTERN_TEMPLATE(TYPENAME)                                         \
                                                                               \
  QLAT_EXTERN template std::vector<char> BitSet::compress(                     \
      const Vector<TYPENAME>& src) const;                                      \
                                                                               \
  QLAT_EXTERN template std::vector<char> BitSet::compress_selected(            \
      const Vector<TYPENAME>& src) const;                                      \
                                                                               \
  QLAT_EXTERN template Long write(FieldsWriter& fw, const std::string& fn,     \
                                  const Field<TYPENAME>& field);               \
                                                                               \
  QLAT_EXTERN template Long write(FieldsWriter& fw, const std::string& fn,     \
                                  const SelectedField<TYPENAME>& sf,           \
                                  const BitSet& bs);                           \
                                                                               \
  QLAT_EXTERN template void set_field_from_data(                               \
      Field<TYPENAME>& field, const GeometryNode& geon,                        \
      const Coordinate& total_site, const std::vector<char>& data,             \
      const bool is_sparse_field);                                             \
                                                                               \
  QLAT_EXTERN template void set_field_from_data(                               \
      SelectedField<TYPENAME>& sf, FieldRank& f_rank,                          \
      const std::vector<char>& data);                                          \
                                                                               \
  QLAT_EXTERN template void set_field_from_data(SelectedField<TYPENAME>& sf,   \
                                                const std::vector<char>& data, \
                                                const FieldSelection& fsel);   \
                                                                               \
  QLAT_EXTERN template Long read(FieldsReader& fr, const std::string& fn,      \
                                 Field<TYPENAME>& field);                      \
                                                                               \
  QLAT_EXTERN template Long read(FieldsReader& fr, const std::string& fn,      \
                                 SelectedField<TYPENAME>& sf,                  \
                                 FieldRank& f_rank);                           \
                                                                               \
  QLAT_EXTERN template Long read(FieldsReader& fr, const std::string& fn,      \
                                 const FieldSelection& fsel,                   \
                                 SelectedField<TYPENAME>& sf);                 \
                                                                               \
  QLAT_EXTERN template Long write(ShuffledFieldsWriter& sfw,                   \
                                  const std::string& fn,                       \
                                  const Field<TYPENAME>& field);               \
                                                                               \
  QLAT_EXTERN template Long write(                                             \
      ShuffledFieldsWriter& sfw, const std::string& fn,                        \
      const ShuffledBitSet& sbs, const SelectedField<TYPENAME>& sf);           \
                                                                               \
  QLAT_EXTERN template void set_field_info_from_fields(                        \
      Coordinate& total_site, Int& multiplicity,                               \
      std::vector<Field<TYPENAME>>& fs, const ShuffledFieldsReader& sfr);      \
                                                                               \
  QLAT_EXTERN template void set_field_info_from_fields(                        \
      Coordinate& total_site, Int& multiplicity,                               \
      std::vector<SelectedField<TYPENAME>>& sfs,                               \
      const ShuffledFieldsReader& sfr);                                        \
                                                                               \
  QLAT_EXTERN template Long read(ShuffledFieldsReader& sfr,                    \
                                 const std::string& fn,                        \
                                 Field<TYPENAME>& field);                      \
                                                                               \
  QLAT_EXTERN template Long read(                                              \
      ShuffledFieldsReader& sfr, const std::string& fn,                        \
      SelectedField<TYPENAME>& sf, FieldSelection& fsel);                      \
                                                                               \
  QLAT_EXTERN template Long read(                                              \
      ShuffledFieldsReader& sfr, const std::string& fn,                        \
      const ShuffledBitSet& sbs, SelectedField<TYPENAME>& sf);                 \
                                                                               \
  QLAT_EXTERN template Long write_float_from_double(                           \
      ShuffledFieldsWriter& sfw, const std::string& fn,                        \
      const Field<TYPENAME>& field);                                           \
                                                                               \
  QLAT_EXTERN template Long write_float_from_double(                           \
      ShuffledFieldsWriter& sfw, const std::string& fn,                        \
      const ShuffledBitSet& sbs, const SelectedField<TYPENAME>& sf);           \
                                                                               \
  QLAT_EXTERN template Long read_double_from_float(ShuffledFieldsReader& sfr,  \
                                                   const std::string& fn,      \
                                                   Field<TYPENAME>& field);    \
                                                                               \
  QLAT_EXTERN template Long read_double_from_float(                            \
      ShuffledFieldsReader& sfr, const std::string& fn,                        \
      SelectedField<TYPENAME>& sf, FieldSelection& fsel);                      \
                                                                               \
  QLAT_EXTERN template Long read_double_from_float(                            \
      ShuffledFieldsReader& sfr, const std::string& fn,                        \
      const ShuffledBitSet& sbs, SelectedField<TYPENAME>& sf);                 \
                                                                               \
  QLAT_EXTERN template Long read_field(                                        \
      Field<TYPENAME>& field, const std::string& path, const std::string& fn); \
                                                                               \
  QLAT_EXTERN template Long read_field(                                        \
      SelectedField<TYPENAME>& sf, const std::string& path,                    \
      const std::string& fn, const ShuffledBitSet& sbs);                       \
                                                                               \
  QLAT_EXTERN template Long read_field_double_from_float(                      \
      Field<TYPENAME>& field, const std::string& path, const std::string& fn); \
                                                                               \
  QLAT_EXTERN template Long read_field_double_from_float(                      \
      SelectedField<TYPENAME>& sf, const std::string& path,                    \
      const std::string& fn, const ShuffledBitSet& sbs);                       \
                                                                               \
  QLAT_EXTERN template Long write_field(const Field<TYPENAME>& field,          \
                                        const std::string& path,               \
                                        const std::string& fn);                \
                                                                               \
  QLAT_EXTERN template Long write_field(                                       \
      const SelectedField<TYPENAME>& sf, const ShuffledBitSet& sbs,            \
      const std::string& path, const std::string& fn);                         \
                                                                               \
  QLAT_EXTERN template Long write_field_float_from_double(                     \
      const Field<TYPENAME>& field, const std::string& path,                   \
      const std::string& fn);                                                  \
                                                                               \
  QLAT_EXTERN template Long write_field_float_from_double(                     \
      const SelectedField<TYPENAME>& sf, const ShuffledBitSet& sbs,            \
      const std::string& path, const std::string& fn)

QLAT_CALL_WITH_TYPES(QLAT_EXTERN_TEMPLATE);
#undef QLAT_EXTERN_TEMPLATE

#undef QLAT_EXTERN

}  // namespace qlat
