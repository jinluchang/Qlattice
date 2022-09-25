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

#include <errno.h>
#include <qlat/selected-field.h>
#include <qlat/selected-points.h>
#include <qlat-utils/qar-cache.h>

namespace qlat
{  //

struct BitSet {
  std::vector<unsigned char> bytes;
  size_t N; // number of uncompressed elements
  size_t cN; // number of compressed elements
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
  void set(const void* pbytes, const size_t nbytes)
  {
    qassert(nbytes == bytes.size());
    memcpy(&bytes[0], pbytes, nbytes);
    cN = 0;
    for (size_t i = 0; i < N; i++)
      if (get(i)) cN++;
  }
  //
  void set(const size_t idx, const bool v)
  {
    qassert(idx < N);
    if (v) {
      if (!get(idx)) cN++;
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
  void set_f_rank(FieldM<int64_t, 1>& f_rank, const int64_t rank = 0)
  {
    TIMER("BitSet::set_f_rank")
    const Geometry& geo = f_rank.geo();
    qassert(geo.local_volume() == (long)N);
    qassert(geo.is_only_local());
    qassert(geo.multiplicity == 1);
    for (size_t i = 0; i < N; i++) {
      if (get(i)) {
        f_rank.get_elem(i) = rank;
      } else {
        f_rank.get_elem(i) = -1;
      }
    }
  }
  bool check_f_rank(const FieldM<int64_t, 1>& f_rank)
  {
    TIMER("BitSet::check_f_rank")
    const Geometry& geo = f_rank.geo();
    qassert(geo.is_only_local());
    qassert(geo.multiplicity == 1);
    if (not(geo.local_volume() == (long)N)) {
      return false;
    }
    for (size_t i = 0; i < N; i++) {
      if (get(i)) {
        if (not(f_rank.get_elem(i) >= 0)) {
          return false;
        }
      } else {
        if (not(f_rank.get_elem(i) == -1)) {
          return false;
        }
      }
    }
    return true;
  }
  //
  void compress(const void* src, void* dst, size_t block_size) const
  {
    unsigned char* psrc = (unsigned char*)src;
    unsigned char* pdst = (unsigned char*)dst;
    size_t j = 0;
    for (size_t i = 0; i < N; i++) {
      if (get(i)) {
        memcpy(&pdst[j * block_size], &psrc[i * block_size], block_size);
        j += 1;
      }
    }
  }
  //
  void decompress(const void* src, void* dst, size_t block_size) const
  {
    unsigned char* psrc = (unsigned char*)src;
    unsigned char* pdst = (unsigned char*)dst;
    size_t j = 0;
    for (size_t i = 0; i < N; i++) {
      if (get(i)) {
        memcpy(&pdst[i * block_size], &psrc[j * block_size], block_size);
        j += 1;
      } else {
        memset(&pdst[i * block_size], 0, block_size);
      }
    }
  }
  //
  template <class M>
  std::vector<char> compress(const Vector<M>& src) const
  {
    qassert(src.size() % N == 0);
    size_t sz_block = sizeof(M) * src.size() / N;
    size_t sz_compressed = sz_block * cN;
    std::vector<char> dst(sz_compressed + bytes.size(), 0);
    memcpy(&dst[0], &bytes[0], bytes.size());
    compress((void*)src.data(), (void*)&dst[bytes.size()], sz_block);
    return dst;
  }
  //
  template <class M>
  std::vector<char> compress_selected(const Vector<M>& src) const
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
};

inline std::vector<char> bitset_decompress(const std::vector<char>& data,
                                           const long local_volume)
{
  TIMER("bitset_decompress");
  const size_t N = local_volume;
  const size_t nbytes = 1 + (N - 1) / 8;
  BitSet bs(N);
  bs.set(&data[0], nbytes);
  const size_t sz_compressed = data.size() - nbytes;
  size_t sz_block = 0;
  if (bs.cN == 0) {
    qassert(sz_compressed == 0);
    std::vector<char> ret;
    return ret;
  } else {
    qassert(sz_compressed % bs.cN == 0);
    sz_block = sz_compressed / bs.cN;
  }
  std::vector<char> ret(sz_block * local_volume, 0);
  bs.decompress(&data[nbytes], &ret[0], sz_block);
  return ret;
}

inline BitSet mk_bitset_from_field_rank(const FieldM<int64_t, 1>& f_rank,
                                        const int64_t n_per_tslice = -1)
{
  TIMER("mk_bitset_from_field_rank");
  const Geometry& geo = f_rank.geo();
  BitSet bs(geo.local_volume());
  qassert(geo.is_only_local());
  for (long index = 0; index < geo.local_volume(); ++index) {
    const int64_t rank = f_rank.get_elem(index);
    if (0 <= rank and (rank < n_per_tslice or n_per_tslice == -1)) {
      bs.set(index, true);
    } else {
      bs.set(index, false);
    }
  }
  return bs;
}

inline void fields_writer_dirs_geon_info(const GeometryNode& geon,
                                         const std::string& path,
                                         const mode_t mode = default_dir_mode())
{
  TIMER("fields_writer_dirs_geon_info");
  dist_mkdir(path, geon.num_node, mode);
  const std::string fn = path + "/geon-info.txt";
  QFile qfile = qfopen(fn, "w");
  qwrite_data(ssprintf("geon.num_node = %d\n", geon.num_node), qfile);
  qwrite_data(ssprintf("geon.size_node[0] = %d\n", geon.size_node[0]), qfile);
  qwrite_data(ssprintf("geon.size_node[1] = %d\n", geon.size_node[1]), qfile);
  qwrite_data(ssprintf("geon.size_node[2] = %d\n", geon.size_node[2]), qfile);
  qwrite_data(ssprintf("geon.size_node[3] = %d\n", geon.size_node[3]), qfile);
}

inline Coordinate shuffled_fields_reader_size_node_info(const std::string& path)
{
  TIMER("shuffled_fields_reader_size_node_info");
  Coordinate size_node;
  Coordinate node_site;
  if (get_id_node() == 0) {
    const std::string fn = path + "/geon-info.txt";
    const std::vector<std::string> lines = qgetlines(fn);
    int num_node;
    reads(num_node, info_get_prop(lines, "geon.num_node = "));
    for (int i = 0; i < 4; ++i) {
      reads(size_node[i],
            info_get_prop(lines, ssprintf("geon.size_node[%d] = ", i)));
    }
    qassert(num_node == product(size_node));
  }
  bcast(get_data(size_node));
  return size_node;
}

struct FieldsWriter {
  //
  // should only use ShuffledFieldsWriter
  //
  std::string path;
  GeometryNode geon;
  QFile qfile;
  bool is_little_endian;  // should be true
  //
  FieldsWriter()
  {
    init();
  }
  //
  void init()
  {
    path = "";
    geon.init();
    qfile.init();
    is_little_endian = true;
  }
  void init(const std::string& path_, const GeometryNode& geon_,
            const bool is_append = false)
  {
    path = path_;
    geon = geon_;
    qfile.init();
    if (geon.id_node == 0) {
      if (does_file_exist(path + ".partial")) {
        if (is_append) {
          qwarn(ssprintf("FieldsWriter: Cannot append '%s.partial' exists.", path.c_str()));
          qassert(false);
        }
        qremove_all(path + ".partial");
      }
      displayln(0, "FieldsWriter: open '" + path + "'.");
      if (is_append and does_file_exist(path)) {
        qassert(does_file_exist(path + "/geon-info.txt"));
      } else {
        fields_writer_dirs_geon_info(geon, path);
      }
    }
    qfile = qfopen(dist_file_name(path, geon.id_node, geon.num_node), is_append ? "a" : "w");
    qassert(not qfile.null());
  }
  //
  void close() { qfile.close(); }
};

struct FieldsReader {
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
  std::map<std::string, long> offsets_map;
  long max_offset;
  //
  FieldsReader()
  {
    init();
  }
  //
  void init()
  {
    path = "";
    geon.init();
    qfile.init();
    is_little_endian = true;
    is_read_through = false;
    fn_list.clear();
    offsets_map.clear();
    max_offset = 0;
  }
  void init(const std::string& path_, const GeometryNode& geon_)
  {
    path = path_;
    geon = geon_;
    qfile.init();
    if (geon.id_node == 0) {
      displayln(0, "FieldsReader: open '" + path + "'.");
    }
    qfile = qfopen(dist_file_name(path, geon.id_node, geon.num_node), "r");
    if (qfile.null()) {
      is_read_through = true;
    } else {
      is_read_through = false;
    }
    fn_list.clear();
    offsets_map.clear();
    max_offset = 0;
  }
  //
  void close() { qfile.close(); }
};

inline void mkfile(FieldsReader& fr)
// create the file (open with appending)
// does not open the file
{
  if (fr.qfile.null() and fr.path != "") {
    fr.qfile =
        qfopen(dist_file_name(fr.path, fr.geon.id_node, fr.geon.num_node), "a");
    fr.qfile.close();
  }
}

inline std::string get_file_path(FieldsReader& fr)
{
  if (fr.path != "") {
    return dist_file_name(fr.path, fr.geon.id_node, fr.geon.num_node);
  }
  return "";
}

inline long get_file_size(FieldsReader& fr)
// the file must be opened for reading
// will restore the position.
{
  if (not fr.qfile.null()) {
    const long pos = qftell(fr.qfile);
    qfseek(fr.qfile, 0L, SEEK_END);
    const long sz = qftell(fr.qfile);
    qfseek(fr.qfile, pos, SEEK_SET);
    return sz;
  } else {
    return -1;
  }
}

template <class M>
void convert_endian_32(Vector<M> data, const bool is_little_endian)
{
  if (is_little_endian) {
    to_from_little_endian_32(data);
  } else {
    to_from_big_endian_32(data);
  }
}

template <class M>
void convert_endian_64(Vector<M> data, const bool is_little_endian)
{
  if (is_little_endian) {
    to_from_little_endian_64(data);
  } else {
    to_from_big_endian_64(data);
  }
}

inline void qfwrite_convert_endian(void* ptr, const size_t size,
                                  const size_t nmemb, QFile& qfile,
                                  const bool is_little_endian)
{
  if (size == 4) {
    convert_endian_32(Vector<int32_t>((int32_t*)ptr, nmemb), is_little_endian);
  } else if (size == 8) {
    convert_endian_64(Vector<int64_t>((int64_t*)ptr, nmemb), is_little_endian);
  } else {
    qassert(false);
  }
  qfwrite(ptr, size, nmemb, qfile);
  if (size == 4) {
    convert_endian_32(Vector<int32_t>((int32_t*)ptr, nmemb), is_little_endian);
  } else if (size == 8) {
    convert_endian_64(Vector<int64_t>((int64_t*)ptr, nmemb), is_little_endian);
  } else {
    qassert(false);
  }
}

inline long write(FieldsWriter& fw, const std::string& fn, const Geometry& geo,
                  const Vector<char> data, const bool is_sparse_field = false)
{
  TIMER("write(fw,fn,geo,data)");
  // first write tag
  int32_t tag_len = fn.size() + 1;  // fn is the name of the field (say prop1)
  qfwrite_convert_endian(&tag_len, 4, 1, fw.qfile, fw.is_little_endian);
  qfwrite(fn.c_str(), tag_len, 1, fw.qfile);
  //
  // then write crc
  crc32_t crc = crc32_par(data);
  qfwrite_convert_endian(&crc, 4, 1, fw.qfile, fw.is_little_endian);
  //
  // then write geometry info
  int32_t nd = 4;  // <- number of dimensions of field, typically 4
  qfwrite_convert_endian(&nd, 4, 1, fw.qfile, fw.is_little_endian);
  //
  std::vector<int32_t> gd(4, 0);
  std::vector<int32_t> num_procs(4, 0);
  //
  const Coordinate total_site = geo.total_site();
  const Coordinate size_node = geo.geon.size_node;
  for (int mu = 0; mu < 4; ++mu) {
    gd[mu] = total_site[mu];
    num_procs[mu] = size_node[mu];
  }
  //
  if (is_sparse_field) {
    gd[0] = -gd[0];
    // note that gd[0] is negative which I use as a flag to tell my reader that
    // data is compressed.  For positive gd[0] my reader assumes a dense field
    // and skips the decompress.  in this way the format is backwards compatible
    // with my old one
  }
  //
  qfwrite_convert_endian(&gd[0], 4, nd, fw.qfile, fw.is_little_endian);
  qfwrite_convert_endian(&num_procs[0], 4, nd, fw.qfile, fw.is_little_endian);
  //
  // then data size
  int64_t data_len = data.size();
  qfwrite_convert_endian(&data_len, 8, 1, fw.qfile, fw.is_little_endian);
  //
  // then write data
  qfwrite(&data[0], data_len, 1, fw.qfile);
  //
  return data_len;
}

inline long qfread_convert_endian(void* ptr, const size_t size,
                                 const size_t nmemb, QFile& qfile,
                                 const bool is_little_endian)
{
  if (qfile.null()) {
    return 0;
  }
  const long total_nmemb = qfread(ptr, size, nmemb, qfile);
  if (size == 4) {
    convert_endian_32(Vector<int32_t>((int32_t*)ptr, nmemb), is_little_endian);
  } else if (size == 8) {
    convert_endian_64(Vector<int64_t>((int64_t*)ptr, nmemb), is_little_endian);
  } else {
    qassert(false);
  }
  return total_nmemb;
}

inline bool read_tag(FieldsReader& fr, std::string& fn, Coordinate& total_site,
                     crc32_t& crc, int64_t& data_len, bool& is_sparse_field)
{
  TIMER("read_tag(fr,fn,total_site,crc,data_len,is_sparse_field)");
  fn = "";
  total_site = Coordinate();
  crc = 0;
  data_len = 0;
  is_sparse_field = false;
  //
  if (fr.qfile.null()) {
    qwarn(ssprintf("read_tag: fr.qfile.null()==true fn='%s'", get_file_path(fr).c_str()));
    return false;
  }
  //
  const long offset_initial = qftell(fr.qfile);
  //
  // first read tag
  int32_t tag_len = 0;
  if (1 != qfread_convert_endian(&tag_len, 4, 1, fr.qfile, fr.is_little_endian)) {
    fr.is_read_through = true;
    return false;
  }
  if (not (tag_len > 0)) {
    qwarn(ssprintf("read_tag: tag_len <= 0 fn='%s'", get_file_path(fr).c_str()));
    fr.is_read_through = true;
    return false;
  }
  std::vector<char> fnv(tag_len);
  if (1 != qfread(fnv.data(), tag_len, 1, fr.qfile)) {
    qwarn(ssprintf("read_tag: fn='%s'", get_file_path(fr).c_str()));
    fr.is_read_through = true;
    return false;
  }
  fn = std::string(fnv.data(), tag_len - 1);
  //
  if (has(fr.offsets_map, fn)) {
    qassert(fr.offsets_map[fn] == offset_initial);
  } else {
    fr.fn_list.push_back(fn);
    fr.offsets_map[fn] = offset_initial;
  }
  //
  // then read crc
  if (1 != qfread_convert_endian(&crc, 4, 1, fr.qfile, fr.is_little_endian)) {
    qwarn(ssprintf("read_tag: fn='%s'", get_file_path(fr).c_str()));
    fr.is_read_through = true;
    return false;
  }
  //
  // then read geometry info
  int32_t nd = 0;
  if (1 != qfread_convert_endian(&nd, 4, 1, fr.qfile, fr.is_little_endian)) {
    qwarn(ssprintf("read_tag: fn='%s'", get_file_path(fr).c_str()));
    fr.is_read_through = true;
    return false;
  }
  if (not (4 == nd)) {
    qwarn(ssprintf("read_tag: fn='%s'", get_file_path(fr).c_str()));
    fr.is_read_through = true;
    return false;
  }
  //
  std::vector<int32_t> gd(4, 0);
  std::vector<int32_t> num_procs(4, 0);
  if (4 != qfread_convert_endian(&gd[0], 4, 4, fr.qfile, fr.is_little_endian)) {
    qwarn(ssprintf("read_tag: fn='%s'", get_file_path(fr).c_str()));
    fr.is_read_through = true;
    return false;
  }
  if (4 !=
      qfread_convert_endian(&num_procs[0], 4, 4, fr.qfile, fr.is_little_endian)) {
    qwarn(ssprintf("read_tag: fn='%s'", get_file_path(fr).c_str()));
    fr.is_read_through = true;
    return false;
  }
  //
  if (gd[0] < 0) {
    is_sparse_field = true;
    gd[0] = -gd[0];
  }
  //
  const Coordinate size_node = fr.geon.size_node;
  for (int mu = 0; mu < 4; ++mu) {
    total_site[mu] = gd[mu];
    if (not(size_node[mu] == (int)num_procs[mu])) {
      qwarn(ssprintf("read_tag: fn='%s'", get_file_path(fr).c_str()));
      fr.is_read_through = true;
      return false;
    }
  }
  //
  // then read data size
  if (1 != qfread_convert_endian(&data_len, 8, 1, fr.qfile, fr.is_little_endian)) {
    qwarn(ssprintf("read_tag: fn='%s'", get_file_path(fr).c_str()));
    fr.is_read_through = true;
    return false;
  }
  if (not (data_len > 0)) {
    qwarn(ssprintf("read_tag: fn='%s'", get_file_path(fr).c_str()));
    fr.is_read_through = true;
    return false;
  }
  //
  const long final_offset = qftell(fr.qfile) + data_len;
  if (final_offset > fr.max_offset) {
    fr.max_offset = final_offset;
  }
  //
  if (fr.geon.id_node == 0) {
    displayln(
        0, fname + ssprintf(": '%s' from '%s'.", fn.c_str(), fr.path.c_str()));
  }
  return true;
}

inline long read_data(FieldsReader& fr, std::vector<char>& data,
                      const int64_t data_len, const crc32_t crc)
// return data_len (if not successful then return 0)
{
  TIMER_FLOPS("read_data(fr,fn,geo,data)");
  clear(data);
  data.resize(data_len, 0);
  if (fr.qfile.null()) {
    qwarn(ssprintf("read_data: file does not exist fn='%s'", get_file_path(fr).c_str()));
    return 0;
  }
  const long read_data_all = qfread(&data[0], data_len, 1, fr.qfile);
  if (not (1 == read_data_all)) {
    qwarn(ssprintf("read_data: data not complete fn='%s'", get_file_path(fr).c_str()));
    fr.is_read_through = true;
    return 0;
  }
  crc32_t crc_read = crc32_par(get_data(data));
  if (not(crc_read == crc)) {
    qwarn(ssprintf("read_data: crc does not match fn='%s'", get_file_path(fr).c_str()));
    return 0;
  }
  timer.flops += data_len;
  return data_len;
}

inline long read_next(FieldsReader& fr, std::string& fn, Coordinate& total_site,
                      std::vector<char>& data, bool& is_sparse_field)
{
  TIMER_FLOPS("read_next(fr,fn,geo,data)");
  crc32_t crc = 0;
  int64_t data_len = 0;
  const bool is_ok =
      read_tag(fr, fn, total_site, crc, data_len, is_sparse_field);
  const long total_bytes = is_ok ? read_data(fr, data, data_len, crc) : 0;
  timer.flops += total_bytes;
  return total_bytes;
}

inline void read_through(FieldsReader& fr)
{
  Coordinate total_site;
  bool is_sparse_field;
  while (true) {
    std::string fn_read = "";
    crc32_t crc = 0;
    int64_t data_len = 0;
    const bool is_ok =
        read_tag(fr, fn_read, total_site, crc, data_len, is_sparse_field);
    if (is_ok) {
      qfseek(fr.qfile, data_len, SEEK_CUR);
    } else {
      errno = 0;
      return;
    }
  }
}

inline bool does_file_exist(FieldsReader& fr, const std::string& fn)
{
  TIMER("does_file_exist(fr,fn,site)");
  if (fr.offsets_map.count(fn) == 1) {
    return true;
  } else if (fr.is_read_through) {
    return false;
  } else if (fr.qfile.null()) {
    return false;
  } else {
    const int ret = qfseek(fr.qfile, fr.max_offset, SEEK_SET);
    if (ret != 0) {
      return false;
    }
  }
  Coordinate total_site;
  bool is_sparse_field;
  while (true) {
    std::string fn_read = "";
    crc32_t crc = 0;
    int64_t data_len = 0;
    const bool is_ok =
        read_tag(fr, fn_read, total_site, crc, data_len, is_sparse_field);
    if (is_ok) {
      if (fn == fn_read) {
        return true;
      } else {
        qfseek(fr.qfile, data_len, SEEK_CUR);
      }
    } else {
      return false;
    }
  }
}

inline long read(FieldsReader& fr, const std::string& fn,
                 Coordinate& total_site, std::vector<char>& data,
                 bool& is_sparse_field)
{
  TIMER_FLOPS("read(fr,fn,site,data)");
  if (not does_file_exist(fr, fn)) {
    return 0;
  }
  qassert(fr.offsets_map.count(fn) == 1);
  qfseek(fr.qfile, fr.offsets_map[fn], SEEK_SET);
  std::string fn_r;
  const long total_bytes =
      read_next(fr, fn_r, total_site, data, is_sparse_field);
  qassert(fn == fn_r);
  return total_bytes;
}

inline long check_file(FieldsReader& fr, const std::string& fn)
// return final offset of the data
// if check_file fail, return 0
{
  TIMER_FLOPS("check_file(fr,fn)");
  Coordinate total_site;
  std::vector<char> data;
  bool is_sparse_field;
  const long total_bytes = read(fr, fn, total_site, data, is_sparse_field);
  if (total_bytes > 0) {
    return qftell(fr.qfile);
  } else {
    return 0;
  }
}

template <class M>
long write(FieldsWriter& fw, const std::string& fn, const Field<M>& field)
// field already have endianess converted correctly
{
  TIMER_FLOPS("write(fw,fn,field)");
  const Geometry& geo = field.geo();
  const Vector<M> v = get_data(field);
  const Vector<char> data((const char*)v.data(), v.data_size());
  const long total_bytes = write(fw, fn, geo, data, false);
  timer.flops += total_bytes;
  return total_bytes;
}

/*
template <class M>
long write(FieldsWriter& fw, const std::string& fn, const Field<M>& field,
           const BitSet& bs)
// field already have endianess converted correctly
{
  TIMER_FLOPS("write(fw,fn,field,bs)");
  const Geometry& geo = field.geo();
  const std::vector<char> data = bs.compress(get_data(field));
  const long total_bytes = write(fw, fn, geo, get_data(data), true);
  timer.flops += total_bytes;
  return total_bytes;
}
*/

template <class M>
long write(FieldsWriter& fw, const std::string& fn, const SelectedField<M>& sf,
           const BitSet& bs)
// field already have endianness converted correctly
{
  TIMER_FLOPS("write(fw,fn,sf,bs)");
  const Geometry& geo = sf.geo();
  const std::vector<char> data = bs.compress_selected(get_data(sf));
  const long total_bytes = write(fw, fn, geo, get_data(data), true);
  timer.flops += total_bytes;
  return total_bytes;
}

inline int flush(FieldsWriter& fw)
{
  TIMER("flush(fw)");
  return qfflush(fw.qfile);
}

template <class M>
void set_field_from_data(Field<M>& field, const GeometryNode& geon,
                         const Coordinate& total_site,
                         const std::vector<char>& data,
                         const bool is_sparse_field)
{
  TIMER("set_field_from_data");
  const Coordinate node_site = total_site / geon.size_node;
  const long local_volume = product(node_site);
  ConstHandle<std::vector<char> > hdata(data);
  std::vector<char> dc_data;
  if (is_sparse_field) {
    dc_data = bitset_decompress(data, local_volume);
    hdata.init(dc_data);
  }
  if (hdata().size() == 0) {
    field.init();
    return;
  }
  const long local_data_size = hdata().size();
  const long site_data_size = local_data_size / local_volume;
  qassert(site_data_size % sizeof(M) == 0);
  const int multiplicity = site_data_size / sizeof(M);
  Geometry geo;
  geo.init(geon, node_site, multiplicity);
  field.init();
  field.init(geo);
  Vector<M> fv = get_data(field);
  qassert(fv.data_size() == (long)hdata().size());
  memcpy(fv.data(), hdata().data(), fv.data_size());
}

template <class M>
void set_field_from_data(SelectedField<M>& sf, FieldM<int64_t, 1>& f_rank,
                         const std::vector<char>& data)
{
  TIMER("set_field_from_data");
  const Geometry& geo = f_rank.geo();
  const Coordinate& node_site = geo.node_site;
  const long local_volume = product(node_site);
  const size_t N = local_volume;
  const size_t nbytes = 1 + (N - 1) / 8;
  BitSet bs(N);
  bs.set(&data[0], nbytes);
  bs.set_f_rank(f_rank);
  const long n_elems = bs.cN;
  if (n_elems == 0) {
    sf.init();
    return;
  }
  const size_t sz_compressed = data.size() - nbytes;
  qassert(sz_compressed % n_elems == 0);
  const size_t sz_block = sz_compressed / n_elems;
  qassert(sz_block % sizeof(M) == 0);
  const int multiplicity = sz_block / sizeof(M);
  const Vector<M> fdata((const M*)&data[nbytes], n_elems * multiplicity);
  sf.init(geo, n_elems, multiplicity);
  assign(get_data(sf), fdata);
}

template <class M>
void set_field_from_data(SelectedField<M>& sf, const std::vector<char>& data,
                         const FieldSelection& fsel)
// fsel must match the actual data
{
  TIMER("set_field_from_data");
  const Geometry& geo = fsel.f_rank.geo();
  const Coordinate& node_site = geo.node_site;
  const long local_volume = product(node_site);
  const size_t N = local_volume;
  const size_t nbytes = 1 + (N - 1) / 8;
  BitSet bs(N);
  bs.set(&data[0], nbytes);
  qassert(bs.check_f_rank(fsel.f_rank));
  const long n_elems = fsel.n_elems;
  qassert(n_elems == (long)bs.cN);
  if (n_elems == 0) {
    sf.init();
    return;
  }
  const size_t sz_compressed = data.size() - nbytes;
  qassert(sz_compressed % n_elems == 0);
  const size_t sz_block = sz_compressed / n_elems;
  qassert(sz_block % sizeof(M) == 0);
  const int multiplicity = sz_block / sizeof(M);
  const Vector<M> fdata((const M*)&data[nbytes], n_elems * multiplicity);
  sf.init(fsel, multiplicity);
  assign(get_data(sf), fdata);
}

template <class M>
long read(FieldsReader& fr, const std::string& fn, Field<M>& field)
// field endianess not converted at all
{
  TIMER_FLOPS("read(fr,fn,field)");
  Coordinate total_site;
  std::vector<char> data;
  bool is_sparse_field = false;
  const long total_bytes = read(fr, fn, total_site, data, is_sparse_field);
  if (0 == total_bytes) {
    return 0;
  }
  set_field_from_data(field, fr.geon, total_site, data, is_sparse_field);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
long read(FieldsReader& fr, const std::string& fn, SelectedField<M>& sf,
          FieldM<int64_t, 1>& f_rank)
// field endianess not converted at all
// f_rank does not need to be initialized
{
  TIMER_FLOPS("read(fr,fn,sf,f_rank)");
  Coordinate total_site;
  std::vector<char> data;
  bool is_sparse_field = false;
  const long total_bytes = read(fr, fn, total_site, data, is_sparse_field);
  if (0 == total_bytes) {
    return 0;
  }
  qassert(is_sparse_field);
  Geometry geo;
  geo.init(fr.geon, total_site / fr.geon.size_node, 1);
  f_rank.init(geo);
  qassert(f_rank.geo().is_only_local());
  set_field_from_data(sf, f_rank, data);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
long read(FieldsReader& fr, const std::string& fn, const FieldSelection& fsel,
          SelectedField<M>& sf)
// field endianess not converted at all
// fsel must match the actual data
// (code will verify & will fail if not match)
{
  TIMER_FLOPS("read(fr,fn,fsel,sf)");
  Coordinate total_site;
  std::vector<char> data;
  bool is_sparse_field = false;
  const long total_bytes = read(fr, fn, total_site, data, is_sparse_field);
  if (0 == total_bytes) {
    return 0;
  }
  qassert(is_sparse_field);
  qassert(total_site == fsel.f_rank.geo().total_site());
  set_field_from_data(sf, data, fsel);
  timer.flops += total_bytes;
  return total_bytes;
}

struct ShuffledBitSet {
  FieldSelection fsel;
  ShufflePlan sp;
  std::vector<FieldSelection> fsels;
  std::vector<BitSet> vbs;
};

inline ShuffledBitSet mk_shuffled_bitset(const FieldM<int64_t, 1>& f_rank,
                                         const int64_t n_per_tslice,
                                         const Coordinate& new_size_node)
// do not enforce n_per_tslice
{
  TIMER("mk_shuffled_bitset(f_rank,n_per_tslice,new_size_node)");
  std::vector<Field<int64_t> > fs_rank;
  shuffle_field(fs_rank, f_rank, new_size_node);
  ShuffledBitSet sbs;
  set_field_selection(sbs.fsel, f_rank, n_per_tslice);
  sbs.sp = make_shuffle_plan(sbs.fsels, sbs.fsel, new_size_node);
  sbs.vbs.resize(fs_rank.size());
  for (int i = 0; i < (int)fs_rank.size(); ++i) {
    FieldM<int64_t, 1> fs_rank_i;
    fs_rank_i.init(fs_rank[i]);
    sbs.vbs[i] = mk_bitset_from_field_rank(fs_rank_i);
  }
  return sbs;
}

inline ShuffledBitSet mk_shuffled_bitset(const FieldSelection& fsel,
                                         const Coordinate& new_size_node)
// interface function
// do not enforce fsel.n_per_tslice
{
  TIMER_VERBOSE("mk_shuffled_bitset(fsel,new_size_node)");
  return mk_shuffled_bitset(fsel.f_rank, fsel.n_per_tslice, new_size_node);
}

inline ShuffledBitSet mk_shuffled_bitset(const Coordinate& total_site,
                                         const std::vector<Coordinate>& xgs,
                                         const Coordinate& new_size_node)
{
  TIMER("mk_shuffled_bitset");
  FieldM<int64_t, 1> f_rank;
  mk_field_selection(f_rank, total_site, xgs);
  return mk_shuffled_bitset(f_rank, 0, new_size_node);
}

inline ShuffledBitSet mk_shuffled_bitset(const FieldM<int64_t, 1>& f_rank,
                                         const std::vector<Coordinate>& xgs,
                                         const Coordinate& new_size_node)
{
  TIMER_VERBOSE("mk_shuffled_bitset");
  const Geometry& geo = f_rank.geo();
  FieldM<int64_t, 1> f_rank_combined;
  f_rank_combined = f_rank;
  const Coordinate total_site = geo.total_site();
  const long spatial_vol = total_site[0] * total_site[1] * total_site[2];
#pragma omp parallel for
  for (long i = 0; i < (long)xgs.size(); ++i) {
    const Coordinate xl = geo.coordinate_l_from_g(xgs[i]);
    if (geo.is_local(xl)) {
      f_rank_combined.get_elem(xl) = spatial_vol + i;
    }
  }
  return mk_shuffled_bitset(f_rank_combined, 0, new_size_node);
}

struct ShuffledFieldsWriter;

typedef std::map<long, Handle<ShuffledFieldsWriter> > ShuffledFieldsWriterMap;

API inline ShuffledFieldsWriterMap& get_all_shuffled_fields_writer()
{
  static ShuffledFieldsWriterMap sfwm;
  return sfwm;
}

inline void add_shuffled_fields_writer(ShuffledFieldsWriter& sfw)
{
  ShuffledFieldsWriterMap& sfwm = get_all_shuffled_fields_writer();
  const long key = (long)&sfw;
  qassert(not has(sfwm, key));
  sfwm[key] = Handle<ShuffledFieldsWriter>(sfw);
}

inline void remove_shuffled_fields_writer(ShuffledFieldsWriter& sfw)
{
  ShuffledFieldsWriterMap& sfwm = get_all_shuffled_fields_writer();
  const long key = (long)&sfw;
  if (has(sfwm, key)) {
    sfwm.erase(key);
  }
}

struct ShuffledFieldsWriter {
  std::string path;
  Coordinate new_size_node;
  std::vector<FieldsWriter> fws;
  //
  ShuffledFieldsWriter() { init(); }
  ShuffledFieldsWriter(const std::string& path_,
                       const Coordinate& new_size_node_,
                       const bool is_append = false)
  // interface function
  {
    init(path_, new_size_node_, is_append);
  }
  //
  ~ShuffledFieldsWriter() { close(); }
  //
  void init()
  // interface function
  {
    close();
    path = "";
    new_size_node = Coordinate();
  }
  void init(const std::string& path_, const Coordinate& new_size_node_,
            const bool is_append = false)
  // interface function
  {
    TIMER_VERBOSE("ShuffledFieldsWriter::init")
    init();
    path = path_;
    if (is_append and does_file_exist_sync_node(path + ".partial")) {
      qwarn(ssprintf("ShuffledFieldsWriter: Cannot append '%s.partial' exists.",
                     path.c_str()));
      qassert(false);
    }
    if (is_append and does_file_exist_sync_node(path)) {
      qassert(does_file_exist_sync_node(path + "/geon-info.txt"));
      new_size_node = shuffled_fields_reader_size_node_info(path);
      if (new_size_node_ != Coordinate() and new_size_node_ != new_size_node) {
        qwarn(ssprintf(
            "ShuffledFieldsWriter::init(p,sn,app): WARNING: new_size_node do "
            "not match. file=%s argument=%s . Will use the new_size_node from "
            "the existing file.",
            show(new_size_node).c_str(), show(new_size_node_).c_str()));
      }
    } else {
      new_size_node = new_size_node_;
    }
    std::vector<GeometryNode> geons = make_dist_io_geons(new_size_node);
    fws.resize(geons.size());
    for (int i = 0; i < (int)geons.size(); ++i) {
        fws[i].init(path, geons[i], is_append);
    }
    sync_node();
    add_shuffled_fields_writer(*this);
  }
  //
  void close()
  // interface function
  {
    TIMER_VERBOSE("ShuffledFieldsWriter::close")
    remove_shuffled_fields_writer(*this);
    clear(fws);
    sync_node();
  }
};

inline void close_all_all_shuffled_fields_writer()
// Force close all the ShuffledFieldsWriter.
// Only call this when quitting the program (e.g. in qquit(msg)).
{
  TIMER_VERBOSE("close_all_all_shuffled_fields_writer");
  ShuffledFieldsWriterMap& sfwm = get_all_shuffled_fields_writer();
  std::vector<Handle<ShuffledFieldsWriter> > sfwv;
  for (auto it = sfwm.begin(); it != sfwm.end(); ++it) {
    sfwv.push_back(it->second);
  }
  for (long i = 0; i < (long)sfwv.size(); ++i) {
    sfwv[i]().close();
  }
  qassert(sfwm.size() == 0);
}

struct ShuffledFieldsReader {
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
  void init()
  // interface function
  {
    path = "";
    new_size_node = Coordinate();
    clear(frs);
  }
  void init(const std::string& path_,
            const Coordinate& new_size_node_ = Coordinate())
  // interface function
  {
    init();
    path = path_;
    if (does_file_exist_qar_sync_node(path + "/geon-info.txt")) {
      new_size_node = shuffled_fields_reader_size_node_info(path);
    } else {
      qassert(new_size_node_ != Coordinate());
      new_size_node = new_size_node_;
    }
    std::vector<GeometryNode> geons = make_dist_io_geons(new_size_node);
    frs.resize(geons.size());
    for (int i = 0; i < (int)geons.size(); ++i) {
      frs[i].init(path, geons[i]);
    }
  }
};

template <class M>
long write(ShuffledFieldsWriter& sfw, const std::string& fn,
           const Field<M>& field)
// interface function
{
  TIMER_VERBOSE_FLOPS("write(sfw,fn,field)");
  displayln_info(0, fname + ssprintf(": writing field with fn='%s' from '%s'.",
                                     fn.c_str(), sfw.path.c_str()));
  std::vector<Field<M> > fs;
  shuffle_field(fs, field, sfw.new_size_node);
  qassert(fs.size() == sfw.fws.size());
  long total_bytes = 0;
  for (int i = 0; i < (int)fs.size(); ++i) {
    total_bytes += write(sfw.fws[i], fn, fs[i]);
  }
  glb_sum(total_bytes);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
long write(ShuffledFieldsWriter& sfw, const std::string& fn,
           const SelectedField<M>& sf, const ShuffledBitSet& sbs)
// interface function
// sbs must match the actual data
{
  TIMER_VERBOSE_FLOPS("write(sfw,fn,sf,sbs)");
  displayln_info(
      0, fname + ssprintf(": writing sparse field with fn='%s' from '%s'.",
                          fn.c_str(), sfw.path.c_str()));
  std::vector<SelectedField<M> > sfs;
  shuffle_field(sfs, sf, sbs.sp);
  qassert(sfs.size() == sfw.fws.size());
  qassert(sbs.vbs.size() == sfw.fws.size());
  long total_bytes = 0;
  for (int i = 0; i < (int)sfs.size(); ++i) {
    total_bytes += write(sfw.fws[i], fn, sfs[i], sbs.vbs[i]);
  }
  glb_sum(total_bytes);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
long write(ShuffledFieldsWriter& sfw, const std::string& fn,
           const Field<M>& field, const ShuffledBitSet& sbs)
// interface function
{
  TIMER_VERBOSE_FLOPS("write(sfw,fn,field,sbs)");
  SelectedField<M> sf;
  set_selected_field(sf, field, sbs.fsel);
  return write(sfw, fn, sf, sbs);
}

inline long flush(ShuffledFieldsWriter& sfw)
// interface function
{
  TIMER_VERBOSE("flush(sfw)");
  long ret = 0;
  for (int i = 0; i < (int)sfw.fws.size(); ++i) {
    ret += flush(sfw.fws[i]);
  }
  glb_sum(ret);
  return ret;
}

inline void read_through_sync_node(ShuffledFieldsReader& sfr)
{
  TIMER_VERBOSE("read_through_sync_node(sfr)");
  for (int i = 0; i < (int)sfr.frs.size(); ++i) {
    read_through(sfr.frs[i]);
  }
  sync_node();
}

inline bool does_file_exist_sync_node(ShuffledFieldsReader& sfr,
                                      const std::string& fn)
// interface function
{
  TIMER_VERBOSE("does_file_exist_sync_node(sfr,fn)");
  long total_counts = 0;
  displayln_info(0, fname + ssprintf(": check fn='%s' from '%s'.", fn.c_str(),
                                     sfr.path.c_str()));
  for (int i = 0; i < (int)sfr.frs.size(); ++i) {
    if (does_file_exist(sfr.frs[i], fn)) {
      total_counts += 1;
    }
  }
  glb_sum(total_counts);
  if (total_counts == 0) {
    return false;
  } else {
    qassert(total_counts == product(sfr.new_size_node));
    return true;
  }
}

inline bool check_file_sync_node(ShuffledFieldsReader& sfr,
                                 const std::string& fn,
                                 std::vector<long>& final_offsets)
// interface function
// set final_offsets to be the files position after loading the data ``fn'' (zero if failed for that file)
// return if data is loaded successfully
{
  TIMER_VERBOSE("check_file_sync_node(sfr,fn)");
  displayln_info(0, fname + ssprintf(": reading field with fn='%s' from '%s'.",
                                     fn.c_str(), sfr.path.c_str()));
  clear(final_offsets);
  final_offsets.resize(sfr.frs.size(), 0);
  long total_failed_counts = 0;
  for (int i = 0; i < (int)sfr.frs.size(); ++i) {
    final_offsets[i] = check_file(sfr.frs[i], fn);
    if (final_offsets[i] == 0) {
      total_failed_counts += 1;
    }
  }
  glb_sum(total_failed_counts);
  bool ret = total_failed_counts == 0;
  displayln_info(0, fname + ssprintf(": check=%s field fn='%s' from '%s'.",
                                     ret ? "true" : "false", fn.c_str(),
                                     sfr.path.c_str()));
  return ret;
}

inline std::vector<std::string> list_fields(ShuffledFieldsReader& sfr)
// interface function
{
  TIMER_VERBOSE("list_fields");
  read_through_sync_node(sfr);
  std::vector<std::string> ret;
  if (0 == get_id_node()) {
    qassert(sfr.frs.size() > 0);
    FieldsReader& fr = sfr.frs[0];
    qassert(fr.is_read_through);
    ret = fr.fn_list;
  }
  bcast(ret);
  return ret;
}

inline int truncate_fields_sync_node(
    const std::string& path, const std::vector<std::string>& fns_keep,
    const Coordinate& new_size_node = Coordinate())
{
  TIMER_VERBOSE("truncate_fields_sync_node");
  ShuffledFieldsReader sfr;
  sfr.init(path, new_size_node);
  const std::vector<std::string> fns = list_fields(sfr);
  if (fns.size() < fns_keep.size()) {
    qwarn(fname + ssprintf(": fns.size()=%ld fns_keep.size()=%ld", fns.size(),
                           fns_keep.size()));
    return 1;
  }
  for (long i = 0; i < (long)fns_keep.size(); ++i) {
    if (fns[i] != fns_keep[i]) {
      qwarn(fname + ssprintf(": fns[i]='%s' fns_keep[i]='%s'", fns[i].c_str(),
                             fns_keep[i].c_str()));
      return 2;
    }
  }
  std::vector<long> final_offsets(sfr.frs.size(), 0);
  if (fns_keep.size() >= 1) {
    const std::string& fn_last = fns_keep.back();
    const bool is_fn_last_valid =
        check_file_sync_node(sfr, fn_last, final_offsets);
    if (not is_fn_last_valid) {
      qwarn(fname + ssprintf(": fn_last='%s' check failed", fn_last.c_str()));
      return 2;
    }
  }
  for (int i = 0; i < (int)sfr.frs.size(); ++i) {
    FieldsReader& fr = sfr.frs[i];
    const std::string path_file = get_file_path(fr);
    const long file_size = get_file_size(fr);
    fr.close();
    const long final_offset = final_offsets[i];
    if (file_size != final_offset) {
      displayln_info(
          0,
          fname +
              ssprintf(
                  ": Truncate '%s': final_offset=%ld, original file_size=%ld.",
                  path_file.c_str(), final_offset, file_size));
      if (file_size < 0) {
        mkfile(fr);
      }
      const bool b = qtruncate(path_file, final_offset);
      qassert(b);
    }
  }
  return 0;
}

inline std::vector<std::string> properly_truncate_fields_sync_node(
    const std::string& path, const bool is_check_all = false,
    const bool is_only_check = false,
    const Coordinate& new_size_node = Coordinate())
// interface function
// return available fns
{
  TIMER_VERBOSE("properly_truncate_fields_sync_node");
  std::vector<std::string> fns;
  if (not does_file_exist_qar_sync_node(path + "/geon-info.txt")) {
    displayln_info(0, fname + ssprintf(": '%s' does not exist.", path.c_str()));
    return fns;
  }
  ShuffledFieldsReader sfr;
  sfr.init(path, new_size_node);
  fns = list_fields(sfr);
  std::vector<long> last_final_offsets(sfr.frs.size(), 0);
  long last_idx = -1;
  if (is_check_all) {
    for (long i = 0; i < (long)fns.size(); ++i) {
      const std::string& fn = fns[i];
      std::vector<long> final_offsets;
      const bool b = check_file_sync_node(sfr, fn, final_offsets);
      if (b) {
        last_final_offsets = final_offsets;
        last_idx = i;
      } else {
        break;
      }
    }
  } else {
    for (long i = (long)fns.size() - 1; i >= 0; i -= 1) {
      const std::string& fn = fns[i];
      std::vector<long> final_offsets;
      const bool b = check_file_sync_node(sfr, fn, final_offsets);
      if (b) {
        last_final_offsets = final_offsets;
        last_idx = i;
        break;
      }
    }
  }
  for (int i = 0; i < (int)sfr.frs.size(); ++i) {
    FieldsReader& fr = sfr.frs[i];
    const std::string path_file = get_file_path(fr);
    const long file_size = get_file_size(fr);
    fr.close();
    const long final_offset = last_final_offsets[i];
    if (file_size != final_offset) {
      if (is_only_check) {
        qwarn(fname +
              ssprintf(
                  ": Error in '%s': final_offset=%ld, original file_size=%ld.",
                  path_file.c_str(), final_offset, file_size));
      } else {
        qwarn(fname +
              ssprintf(
                  ": Truncate '%s': final_offset=%ld, original file_size=%ld.",
                  path_file.c_str(), final_offset, file_size));
        if (file_size < 0) {
          mkfile(fr);
        }
        const bool b = qtruncate(path_file, final_offset);
        qassert(b);
      }
    }
  }
  errno = 0;
  fns.resize(last_idx + 1);
  for (long i = 0; i < (long)fns.size(); ++i) {
    const std::string& fn = fns[i];
    displayln_info(0, fname + ssprintf(": i=%5ld fn='%s'", i, fn.c_str()));
  }
  displayln_info(
      0, fname + ssprintf(": fns.size()=%5ld '%s'", fns.size(), path.c_str()));
  sync_node();
  return fns;
}

template <class M>
void set_field_info_from_fields(Coordinate& total_site, int& multiplicity,
                                std::vector<Field<M> >& fs,
                                const ShuffledFieldsReader& sfr)
{
  TIMER_VERBOSE("set_field_info_from_fields");
  total_site = Coordinate();
  multiplicity = 0;
  std::vector<long> available_nodes(product(sfr.new_size_node), 0);
  for (int i = 0; i < (int)fs.size(); ++i) {
    const int id_node = sfr.frs[i].geon.id_node;
    qassert(0 <= id_node and id_node < (int)available_nodes.size());
    if (fs[i].initialized) {
      available_nodes[id_node] = get_id_node() + 1;
    }
  }
  glb_sum(available_nodes);
  int id_node_first_available = 0;
  int id_node_bcast_from = 0;
  for (int i = 0; i < (int)available_nodes.size(); ++i) {
    if (available_nodes[i] > 0) {
      id_node_first_available = i;
      id_node_bcast_from = available_nodes[i] - 1;
      break;
    }
  }
  for (int i = 0; i < (int)fs.size(); ++i) {
    const int id_node = sfr.frs[i].geon.id_node;
    if (id_node == id_node_first_available) {
      total_site = fs[i].geo().total_site();
      multiplicity = fs[i].geo().multiplicity;
      qassert(get_id_node() == id_node_bcast_from);
    }
  }
  bcast(get_data_one_elem(total_site), id_node_bcast_from);
  bcast(get_data_one_elem(multiplicity), id_node_bcast_from);
  for (int i = 0; i < (int)fs.size(); ++i) {
    if (not fs[i].initialized) {
      const GeometryNode& geon = sfr.frs[i].geon;
      const Coordinate node_site = total_site / geon.size_node;
      Geometry geo;
      geo.init(geon, node_site, multiplicity);
      fs[i].init(geo);
      set_zero(fs[i]);
    }
  }
}

template <class M>
void set_field_info_from_fields(Coordinate& total_site, int& multiplicity,
                                std::vector<SelectedField<M> >& sfs,
                                const ShuffledFieldsReader& sfr)
{
  TIMER_VERBOSE("set_field_info_from_fields");
  total_site = Coordinate();
  multiplicity = 0;
  std::vector<long> available_nodes(product(sfr.new_size_node), 0);
  for (int i = 0; i < (int)sfs.size(); ++i) {
    const int id_node = sfr.frs[i].geon.id_node;
    qassert(0 <= id_node and id_node < (int)available_nodes.size());
    if (is_initialized(sfs[i])) {
      available_nodes[id_node] = get_id_node() + 1;
    }
  }
  glb_sum(available_nodes);
  int id_node_first_available = 0;
  int id_node_bcast_from = 0;
  for (int i = 0; i < (int)available_nodes.size(); ++i) {
    if (available_nodes[i] > 0) {
      id_node_first_available = i;
      id_node_bcast_from = available_nodes[i] - 1;
      break;
    }
  }
  for (int i = 0; i < (int)sfs.size(); ++i) {
    const int id_node = sfr.frs[i].geon.id_node;
    if (id_node == id_node_first_available) {
      qassert(is_initialized(sfs[i]));
      total_site = sfs[i].geo().total_site();
      multiplicity = sfs[i].geo().multiplicity;
      qassert(get_id_node() == id_node_bcast_from);
    }
  }
  bcast(get_data_one_elem(total_site), id_node_bcast_from);
  bcast(get_data_one_elem(multiplicity), id_node_bcast_from);
  for (int i = 0; i < (int)sfs.size(); ++i) {
    if (not sfs[i].initialized) {
      const GeometryNode& geon = sfr.frs[i].geon;
      const Coordinate node_site = total_site / geon.size_node;
      Geometry geo;
      geo.init(geon, node_site, multiplicity);
      sfs[i].init(geo, 0, multiplicity);
      set_zero(sfs[i]);
    }
  }
}

template <class M>
long read(ShuffledFieldsReader& sfr, const std::string& fn, Field<M>& field)
// interface function
{
  TIMER_VERBOSE_FLOPS("read(sfr,fn,field)");
  long total_bytes = 0;
  displayln_info(0, fname + ssprintf(": reading field with fn='%s' from '%s'.",
                                     fn.c_str(), sfr.path.c_str()));
  std::vector<Field<M> > fs(sfr.frs.size());
  long zero_size_count = 0;
  for (int i = 0; i < (int)fs.size(); ++i) {
    const long bytes = read(sfr.frs[i], fn, fs[i]);
    if (0 == bytes) {
      zero_size_count += 1;
      qassert(0 == total_bytes);
    } else {
      total_bytes += bytes;
    }
  }
  glb_sum(total_bytes);
  if (0 != zero_size_count) {
    qassert(0 == total_bytes);
  }
  if (0 == total_bytes) {
    return 0;
  }
  Coordinate total_site;
  int multiplicity = 0;
  set_field_info_from_fields(total_site, multiplicity, fs, sfr);
  Geometry geo;
  geo.init(total_site, multiplicity);
  field.init(geo);
  shuffle_field_back(field, fs, sfr.new_size_node);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
long read(ShuffledFieldsReader& sfr, const std::string& fn,
          SelectedField<M>& sf, FieldSelection& fsel)
// interface function
{
  TIMER_VERBOSE_FLOPS("read(sfr,fn,sf,fsel)")
  long total_bytes = 0;
  displayln_info(0, fname + ssprintf(": reading field with fn='%s' from '%s'.",
                                     fn.c_str(), sfr.path.c_str()));
  std::vector<SelectedField<M> > sfs(sfr.frs.size());
  std::vector<Field<int64_t> > f_rank_s(sfr.frs.size());
  long zero_size_count = 0;
  for (int i = 0; i < (int)sfs.size(); ++i) {
    FieldM<int64_t, 1>& f_rank = static_cast<FieldM<int64_t, 1>&>(f_rank_s[i]);
    const long bytes = read(sfr.frs[i], fn, sfs[i], f_rank);
    if (0 == bytes) {
      zero_size_count += 1;
      qassert(0 == total_bytes);
    } else {
      total_bytes += bytes;
    }
  }
  glb_sum(total_bytes);
  if (0 != zero_size_count) {
    qassert(0 == total_bytes);
  }
  if (0 == total_bytes) {
    return 0;
  }
  Coordinate total_site;
  int multiplicity = 0;
  set_field_info_from_fields(total_site, multiplicity, sfs, sfr);
  const Geometry geo(total_site, 1);
  fsel.f_rank.init(geo);
  shuffle_field_back(fsel.f_rank, f_rank_s, sfr.new_size_node);
  update_field_selection(fsel);
  sf.init(fsel, multiplicity);
  shuffle_field_back(sf, sfs, sfr.new_size_node, fsel);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
long read(ShuffledFieldsReader& sfr, const std::string& fn,
          const ShuffledBitSet& sbs, SelectedField<M>& sf)
// interface function
// sbs must match the actual data
// (code will verify & will fail if not match)
{
  TIMER_VERBOSE_FLOPS("read(sfr,fn,sbs,sf)");
  sf.init();
  long total_bytes = 0;
  displayln_info(
      0, fname + ssprintf(": reading sparse field with fn='%s' from '%s'.",
                          fn.c_str(), sfr.path.c_str()));
  std::vector<SelectedField<M> > sfs(sfr.frs.size());
  long zero_size_count = 0;
  for (int i = 0; i < (int)sfs.size(); ++i) {
    const long bytes = read(sfr.frs[i], fn, sbs.fsels[i], sfs[i]);
    if (0 == bytes) {
      zero_size_count += 1;
      qassert(0 == total_bytes);
    } else {
      total_bytes += bytes;
    }
  }
  glb_sum(total_bytes);
  if (0 != zero_size_count) {
    qassert(0 == total_bytes);
  }
  if (0 == total_bytes) {
    return 0;
  }
  Coordinate total_site;
  int multiplicity = 0;
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
long write_float_from_double(ShuffledFieldsWriter& sfw, const std::string& fn,
                             const Field<M>& field)
// interface function
{
  TIMER_VERBOSE_FLOPS("write_float_from_double(sfw,fn,field)");
  Field<float> ff;
  convert_field_float_from_double(ff, field);
  to_from_little_endian_32(get_data(ff));
  const long total_bytes = write(sfw, fn, ff);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
long write_float_from_double(ShuffledFieldsWriter& sfw, const std::string& fn,
                             const SelectedField<M>& sf, const ShuffledBitSet& sbs)
// interface function
// sbs must match the actual data
{
  TIMER_VERBOSE_FLOPS("write_float_from_double(sfw,fn,sf,sbs)");
  SelectedField<float> sff;
  convert_field_float_from_double(sff, sf);
  to_from_little_endian_32(get_data(sff));
  const long total_bytes = write(sfw, fn, sff, sbs);
  timer.flops += total_bytes;
  return total_bytes;
}

template <class M>
long write_float_from_double(ShuffledFieldsWriter& sfw, const std::string& fn,
                             const Field<M>& field, const ShuffledBitSet& sbs)
// interface function
{
  TIMER_VERBOSE_FLOPS("write_float_from_double(sfw,fn,field,sbs)");
  SelectedField<M> sf;
  set_selected_field(sf, field, sbs.fsel);
  return write_float_from_double(sfw, fn, sf, sbs);
}

template <class M>
long read_double_from_float(ShuffledFieldsReader& sfr, const std::string& fn,
                            Field<M>& field)
// interface function
{
  TIMER_VERBOSE_FLOPS("read_double_from_float(sfr,fn,field)");
  Field<float> ff;
  const long total_bytes = read(sfr, fn, ff);
  if (total_bytes == 0) {
    return 0;
  } else {
    to_from_little_endian_32(get_data(ff));
    convert_field_double_from_float(field, ff);
    timer.flops += total_bytes;
    return total_bytes;
  }
}

template <class M>
long read_double_from_float(ShuffledFieldsReader& sfr, const std::string& fn,
                            SelectedField<M>& sf, FieldSelection& fsel)
// interface function
{
  TIMER_VERBOSE_FLOPS("read_double_from_float(sfr,fn,sf,fsel)");
  sf.init();
  SelectedField<float> sff;
  const long total_bytes = read(sfr, fn, sff, fsel);
  if (total_bytes == 0) {
    return 0;
  } else {
    to_from_little_endian_32(get_data(sff));
    convert_field_double_from_float(sf, sff);
    timer.flops += total_bytes;
    return total_bytes;
  }
}

template <class M>
long read_double_from_float(ShuffledFieldsReader& sfr, const std::string& fn,
                            const ShuffledBitSet& sbs, SelectedField<M>& sf)
// interface function
// sbs must match the actual data
{
  TIMER_VERBOSE_FLOPS("read_double_from_float(sfr,fn,sbs,sf)");
  sf.init();
  SelectedField<float> sff;
  const long total_bytes = read(sfr, fn, sbs, sff);
  if (total_bytes == 0) {
    return 0;
  } else {
    to_from_little_endian_32(get_data(sff));
    convert_field_double_from_float(sf, sff);
    timer.flops += total_bytes;
    return total_bytes;
  }
}

typedef Cache<std::string, ShuffledFieldsReader> ShuffledFieldsReaderCache;

API inline ShuffledFieldsReaderCache& get_shuffled_fields_reader_cache()
{
  static ShuffledFieldsReaderCache cache("ShuffledFieldsReaderCache", 4, 1);
  return cache;
}

inline ShuffledFieldsReader& get_shuffled_fields_reader(
    const std::string& path, const Coordinate& new_size_node = Coordinate())
{
  TIMER("get_shuffled_fields_reader");
  ShuffledFieldsReader& sfr = get_shuffled_fields_reader_cache()[path];
  if (sfr.path == "") {
    sfr.init(path, new_size_node);
  }
  return sfr;
}

template <class M>
long read_field(Field<M>& field, const std::string& path, const std::string& fn)
// interface function
{
  TIMER_VERBOSE("read_field(field,path,fn)");
  ShuffledFieldsReader& sfr = get_shuffled_fields_reader(path);
  return read(sfr, fn, field);
}

template <class M>
long read_field(SelectedField<M>& sf, const std::string& path, const std::string& fn, const ShuffledBitSet& sbs)
// interface function
// sbs must match the actual data
{
  TIMER_VERBOSE("read_field(sf,path,fn,sbs)");
  ShuffledFieldsReader& sfr = get_shuffled_fields_reader(path);
  return read(sfr, fn, sbs, sf);
}

template <class M>
long read_field_double_from_float(Field<M>& field, const std::string& path,
                                  const std::string& fn)
// interface function
{
  TIMER_VERBOSE("read_field_double_from_float(field,path,fn)");
  ShuffledFieldsReader& sfr = get_shuffled_fields_reader(path);
  return read_double_from_float(sfr, fn, field);
}

template <class M>
long read_field_double_from_float(SelectedField<M>& sf,
                                  const std::string& path,
                                  const std::string& fn,
                                  const ShuffledBitSet& sbs)
// interface function
// sbs must match the actual data
{
  TIMER_VERBOSE("read_field_double_from_float(sf,path,fn,sbs)");
  ShuffledFieldsReader& sfr = get_shuffled_fields_reader(path);
  return read_double_from_float(sfr, fn, sbs, sf);
}

inline std::vector<std::string> list_fields(const std::string& path)
{
  ShuffledFieldsReader& sfr = get_shuffled_fields_reader(path);
  return list_fields(sfr);
}

inline bool does_file_exist_sync_node(const std::string& path,
                                      const std::string& fn)
{
  TIMER_VERBOSE("does_file_exist_sync_node(path,fn)");
  ShuffledFieldsReader& sfr = get_shuffled_fields_reader(path);
  return does_file_exist_sync_node(sfr, fn);
}

// old code

/*
template <class M>
long read_next(FieldsReader& fr, std::string& fn, Field<M>& field)
// field endianess not converted at all
{
  TIMER_FLOPS("read_next(fr,fn,field)");
  Coordinate total_site;
  std::vector<char> data;
  bool is_sparse_field = false;
  const long total_bytes = read_next(fr, fn, total_site, data, is_sparse_field);
  if (0 == total_bytes) {
    return 0;
  }
  set_field_from_data(field, fr.geon, total_site, data, is_sparse_field);
  timer.flops += total_bytes;
  return total_bytes;
}
*/

/*
template <class M>
long read_next(ShuffledFieldsReader& sfr, std::string& fn, Field<M>& field)
// interface function
{
  TIMER_VERBOSE_FLOPS("read_next(sfr,fn,field)");
  fn = "";
  long total_bytes = 0;
  std::vector<Field<M> > fs(sfr.frs.size());
  for (int i = 0; i < (int)fs.size(); ++i) {
    std::string fni;
    const long bytes = read_next(sfr.frs[i], fni, fs[i]);
    if (0 == bytes) {
      qassert(0 == total_bytes);
    } else {
      total_bytes += bytes;
    }
    const int id_node = sfr.frs[i].geon.id_node;
    if (id_node == 0) {
      fn = fni;
      qassert(get_id_node() == 0);
    }
  }
  bcast(fn);
  glb_sum(total_bytes);
  if (0 == total_bytes) {
    fn = "";
    return 0;
  }
  displayln_info(0, fname +
                 ssprintf(": read the next field with fn='%s'.", fn.c_str()));
  Coordinate total_site;
  int multiplicity = 0;
  set_field_info_from_fields(total_site, multiplicity, fs, sfr);
  Geometry geo;
  geo.init(total_site, multiplicity);
  field.init(geo);
  shuffle_field_back(field, fs, sfr.new_size_node);
  timer.flops += total_bytes;
  return total_bytes;
}
*/

/*
template <class M>
long read_next_double_from_float(ShuffledFieldsReader& sfr, std::string& fn,
                                 Field<M>& field)
// interface function
{
  TIMER_VERBOSE_FLOPS("read_next_double_from_float(sfr,fn,field)");
  Field<float> ff;
  const long total_bytes = read_next(sfr, fn, ff);
  if (total_bytes == 0) {
    return 0;
  } else {
    to_from_little_endian_32(get_data(ff));
    convert_field_double_from_float(field, ff);
    timer.flops += total_bytes;
    return total_bytes;
  }
}
*/

}  // namespace qlat
