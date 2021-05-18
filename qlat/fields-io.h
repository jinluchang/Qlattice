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

#include <qlat/selected-field.h>
#include <qlat/selected-points.h>

namespace qlat
{  //

struct BitSet {
  std::vector<unsigned char> bytes;
  size_t N, cN;
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
  FILE* fp = qopen(fn, "w");
  displayln(ssprintf("geon.num_node = %d", geon.num_node), fp);
  displayln(ssprintf("geon.size_node[0] = %d", geon.size_node[0]), fp);
  displayln(ssprintf("geon.size_node[1] = %d", geon.size_node[1]), fp);
  displayln(ssprintf("geon.size_node[2] = %d", geon.size_node[2]), fp);
  displayln(ssprintf("geon.size_node[3] = %d", geon.size_node[3]), fp);
  qclose(fp);
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
  FILE* fp;
  bool is_little_endian;  // should be true
  //
  FieldsWriter()
  {
    fp = NULL;
    init();
  }
  //
  ~FieldsWriter() { close(); }
  //
  void init()
  {
    close();
    path = "";
    geon.init();
    qassert(fp == NULL);
    is_little_endian = true;
  }
  void init(const std::string& path_, const GeometryNode& geon_,
            const bool is_append = false)
  {
    close();
    path = path_;
    geon = geon_;
    qassert(fp == NULL);
    if (geon.id_node == 0) {
      if (does_file_exist(path + ".partial")) {
        if (is_append) {
          displayln(ssprintf("FieldsWriter: Cannot append '%s.partial' exists.", path.c_str()));
          qassert(false);
        }
        qremove_all(path + ".partial");
      }
      displayln("FieldsWriter: open '" + path + ".partial'.");
      if (is_append and does_file_exist(path)) {
        qassert(does_file_exist(path + "/geon-info.txt"));
        qrename(path, path + ".partial");
      } else {
        fields_writer_dirs_geon_info(geon, path + ".partial");
      }
    }
    fp = dist_open(path + ".partial", geon.id_node, geon.num_node,
                   is_append ? "a" : "w");
    qassert(NULL != fp);
  }
  //
  void close()
  {
    const bool is_need_close = fp != NULL;
    if (is_need_close and geon.id_node == 0) {
      displayln("FieldsWriter: close '" + path + ".partial" + "'.");
    }
    dist_close(fp);
    if (is_need_close and geon.id_node == 0) {
      qrename(path + ".partial", path);
    }
    qassert(fp == NULL);
  }
};

struct FieldsReader {
  //
  // should only use ShuffledFieldsReader
  //
  std::string path;
  GeometryNode geon;
  FILE* fp;
  bool is_little_endian;  // should be true
  //
  bool is_read_through;
  std::vector<std::string> fn_list;
  std::map<std::string, long> offsets_map;
  long max_offset;
  //
  FieldsReader()
  {
    fp = NULL;
    init();
  }
  //
  ~FieldsReader() { close(); }
  //
  void init()
  {
    close();
    path = "";
    geon.init();
    qassert(fp == NULL);
    is_little_endian = true;
    is_read_through = false;
    fn_list.clear();
    offsets_map.clear();
    max_offset = 0;
  }
  void init(const std::string& path_, const GeometryNode& geon_)
  {
    close();
    path = path_;
    geon = geon_;
    qassert(fp == NULL);
    if (geon.id_node == 0) {
      displayln("FieldsReader: open '" + path + "'.");
    }
    fp = dist_open(path, geon.id_node, geon.num_node, "r");
    qassert(NULL != fp);
    is_read_through = false;
    fn_list.clear();
    offsets_map.clear();
    max_offset = 0;
  }
  //
  std::string path_file()
  {
    if (path != "") {
      return dist_file_name(path, geon.id_node, geon.num_node);
    }
    return "";
  }
  //
  void close()
  {
    if (fp != NULL and geon.id_node == 0) {
      displayln("FieldsReader: close '" + path + "'.");
    }
    dist_close(fp);
    qassert(fp == NULL);
  }
};

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

inline void fwrite_convert_endian(void* ptr, const size_t size,
                                  const size_t nmemb, FILE* fp,
                                  const bool is_little_endian)
{
  if (size == 4) {
    convert_endian_32(Vector<int32_t>((int32_t*)ptr, nmemb), is_little_endian);
  } else if (size == 8) {
    convert_endian_64(Vector<int64_t>((int64_t*)ptr, nmemb), is_little_endian);
  } else {
    qassert(false);
  }
  fwrite(ptr, size, nmemb, fp);
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
  fwrite_convert_endian(&tag_len, 4, 1, fw.fp, fw.is_little_endian);
  fwrite(fn.c_str(), tag_len, 1, fw.fp);
  //
  // then write crc
  crc32_t crc = crc32_par(data);
  fwrite_convert_endian(&crc, 4, 1, fw.fp, fw.is_little_endian);
  //
  // then write geometry info
  int32_t nd = 4;  // <- number of dimensions of field, typically 4
  fwrite_convert_endian(&nd, 4, 1, fw.fp, fw.is_little_endian);
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
  fwrite_convert_endian(&gd[0], 4, nd, fw.fp, fw.is_little_endian);
  fwrite_convert_endian(&num_procs[0], 4, nd, fw.fp, fw.is_little_endian);
  //
  // then data size
  int64_t data_len = data.size();
  fwrite_convert_endian(&data_len, 8, 1, fw.fp, fw.is_little_endian);
  //
  // then write data
  fwrite(&data[0], data_len, 1, fw.fp);
  //
  return data_len;
}

inline long fread_convert_endian(void* ptr, const size_t size,
                                 const size_t nmemb, FILE* fp,
                                 const bool is_little_endian)
{
  const long total_nmemb = fread(ptr, size, nmemb, fp);
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
  TIMER("read_tag(fr,fn,geo)");
  fn = "";
  total_site = Coordinate();
  crc = 0;
  data_len = 0;
  is_sparse_field = false;
  //
  const long offset_initial = ftell(fr.fp);
  //
  // first read tag
  int32_t tag_len = 0;
  if (1 != fread_convert_endian(&tag_len, 4, 1, fr.fp, fr.is_little_endian)) {
    fr.is_read_through = true;
    return false;
  }
  if (not (tag_len > 0)) {
    qwarn("read_tag: tag_len <= 0");
    fr.is_read_through = true;
    return false;
  }
  std::vector<char> fnv(tag_len);
  if (1 != fread(fnv.data(), tag_len, 1, fr.fp)) {
    qwarn("read_tag:");
    fr.is_read_through = true;
    return false;
  }
  fn = std::string(fnv.data(), tag_len - 1);
  //
  fr.fn_list.push_back(fn);
  fr.offsets_map[fn] = offset_initial;
  //
  // then read crc
  if (1 != fread_convert_endian(&crc, 4, 1, fr.fp, fr.is_little_endian)) {
    qwarn("read_tag:");
    fr.is_read_through = true;
    return false;
  }
  //
  // then read geometry info
  int32_t nd = 0;
  if (1 != fread_convert_endian(&nd, 4, 1, fr.fp, fr.is_little_endian)) {
    qwarn("read_tag:");
    fr.is_read_through = true;
    return false;
  }
  if (not (4 == nd)) {
    qwarn("read_tag:");
    fr.is_read_through = true;
    return false;
  }
  //
  std::vector<int32_t> gd(4, 0);
  std::vector<int32_t> num_procs(4, 0);
  if (4 != fread_convert_endian(&gd[0], 4, 4, fr.fp, fr.is_little_endian)) {
    qwarn("read_tag:");
    fr.is_read_through = true;
    return false;
  }
  if (4 !=
      fread_convert_endian(&num_procs[0], 4, 4, fr.fp, fr.is_little_endian)) {
    qwarn("read_tag:");
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
      qwarn("read_tag:");
      fr.is_read_through = true;
      return false;
    }
  }
  //
  // then read data size
  if (1 != fread_convert_endian(&data_len, 8, 1, fr.fp, fr.is_little_endian)) {
    qwarn("read_tag:");
    fr.is_read_through = true;
    return false;
  }
  if (not (data_len > 0)) {
    qwarn("read_tag:");
    fr.is_read_through = true;
    return false;
  }
  //
  const long final_offset = ftell(fr.fp) + data_len;
  if (final_offset > fr.max_offset) {
    fr.max_offset = final_offset;
  }
  //
  if (fr.geon.id_node == 0) {
    displayln(fname +
              ssprintf(": '%s' from '%s'.", fn.c_str(), fr.path.c_str()));
  }
  return true;
}

inline long read_data(FieldsReader& fr, std::vector<char>& data,
                      const int64_t data_len, const crc32_t crc)
{
  TIMER_FLOPS("read_data(fr,fn,geo,data)");
  clear(data);
  data.resize(data_len, 0);
  const long read_data_all = fread(&data[0], data_len, 1, fr.fp);
  if (not (1 == read_data_all)) {
    qwarn("read_data: data not complete");
    fr.is_read_through = true;
    return 0;
  }
  crc32_t crc_read = crc32_par(get_data(data));
  if (not(crc_read == crc)) {
    qwarn("read_data: crc does not match");
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
      fseek(fr.fp, data_len, SEEK_CUR);
    } else {
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
  } else {
    const int ret = fseek(fr.fp, fr.max_offset, SEEK_SET);
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
        fseek(fr.fp, data_len, SEEK_CUR);
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
  fseek(fr.fp, fr.offsets_map[fn], SEEK_SET);
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
  if (not does_file_exist(fr, fn)) {
    return false;
  }
  qassert(fr.offsets_map.count(fn) == 1);
  fseek(fr.fp, fr.offsets_map[fn], SEEK_SET);
  std::string fn_r;
  Coordinate total_site;
  std::vector<char> data;
  bool is_sparse_field;
  const long total_bytes =
      read_next(fr, fn_r, total_site, data, is_sparse_field);
  qassert(fn == fn_r);
  if (total_bytes > 0) {
    return ftell(fr.fp);
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
// field already have endianess converted correctly
{
  TIMER_FLOPS("write(fw,fn,sf,bs)");
  const Geometry& geo = sf.geo();
  const std::vector<char> data = bs.compress_selected(get_data(sf));
  const long total_bytes = write(fw, fn, geo, get_data(data), true);
  timer.flops += total_bytes;
  return total_bytes;
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
void set_field_from_data_fsel(SelectedField<M>& sf,
                              const std::vector<char>& data,
                              const FieldSelection& fsel)
{
  TIMER("set_field_from_data_fsel");
  const Geometry& geo = fsel.f_rank.geo();
  const long n_elems = fsel.n_elems;
  if (n_elems == 0) {
    sf.init();
    return;
  }
  const Coordinate& node_site = geo.node_site;
  const long local_volume = product(node_site);
  const size_t N = local_volume;
  const size_t nbytes = 1 + (N - 1) / 8;
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
long read(FieldsReader& fr, const std::string& fn, const FieldSelection& fsel,
          SelectedField<M>& sf)
// field endianess not converted at all
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
  set_field_from_data_fsel(sf, data, fsel);
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

inline ShuffledFieldsWriterMap& get_all_shuffled_fields_writer()
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
    init();
    path = path_;
    if (is_append and does_file_exist_sync_node(path + ".partial")) {
      displayln_info(
          ssprintf("ShuffledFieldsWriter: Cannot append '%s.partial' exists.",
                   path.c_str()));
      qassert(false);
    }
    if (is_append and does_file_exist_sync_node(path)) {
      qassert(does_file_exist_sync_node(path + "/geon-info.txt"));
      new_size_node = shuffled_fields_reader_size_node_info(path);
      if (new_size_node_ != Coordinate() and new_size_node_ != new_size_node) {
        displayln_info(ssprintf(
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
      if (geons[i].id_node == 0) {
        fws[i].init(path, geons[i], is_append);
      }
    }
    sync_node();
    for (int i = 0; i < (int)geons.size(); ++i) {
      if (geons[i].id_node != 0) {
        fws[i].init(path, geons[i], is_append);
      }
    }
    add_shuffled_fields_writer(*this);
  }
  //
  void close()
  // interface function
  {
    remove_shuffled_fields_writer(*this);
    std::vector<GeometryNode> geons = make_dist_io_geons(new_size_node);
    for (int i = 0; i < (int)fws.size(); ++i) {
      if (geons[i].id_node != 0) {
        fws[i].close();
      }
    }
    sync_node();
    for (int i = 0; i < (int)fws.size(); ++i) {
      if (geons[i].id_node == 0) {
        fws[i].close();
      }
    }
    clear(fws);
  }
};

inline void close_all_all_shuffled_fields_writer()
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
    if (does_file_exist_sync_node(path + "/geon-info.txt")) {
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
  displayln_info(fname + ssprintf(": writing field with fn='%s' from '%s'.",
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
{
  TIMER_VERBOSE_FLOPS("write(sfw,fn,sf,sbs)");
  displayln_info(fname +
                 ssprintf(": writing sparse field with fn='%s' from '%s'.",
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
  displayln_info(fname + ssprintf(": check fn='%s' from '%s'.", fn.c_str(),
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

inline std::vector<long> check_file_sync_node(ShuffledFieldsReader& sfr,
                                              const std::string& fn)
// interface function
// return final offset if all succeed
// return empty vector if any fail
{
  TIMER_VERBOSE("check_file_sync_node(sfr,fn)");
  std::vector<long> ret(sfr.frs.size(), 0);
  displayln_info(fname + ssprintf(": reading field with fn='%s' from '%s'.",
                                  fn.c_str(), sfr.path.c_str()));
  long total_failed_counts = 0;
  for (int i = 0; i < (int)sfr.frs.size(); ++i) {
    ret[i] = check_file(sfr.frs[i], fn);
    if (ret[i] == 0) {
      total_failed_counts += 1;
    }
  }
  glb_sum(total_failed_counts);
  if (total_failed_counts > 0) {
    clear(ret);
  }
  displayln_info(fname + ssprintf(": check=%s field fn='%s' from '%s'.",
                                  ret.size() == 0 ? "false" : "true",
                                  fn.c_str(), sfr.path.c_str()));
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

inline std::vector<std::string> properly_truncate_fields_sync_node(
    const std::string& path, const bool is_check_all = false,
    const Coordinate& new_size_node = Coordinate())
// interface function
// return available fns
{
  TIMER_VERBOSE("properly_truncate_fields_sync_node");
  std::vector<std::string> fns;
  if (not does_file_exist_sync_node(path)) {
    displayln_info(
        fname + ssprintf(": '%s' does not exist.", path.c_str()));
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
      const std::vector<long> final_offsets = check_file_sync_node(sfr, fn);
      if (0 != final_offsets.size()) {
        last_final_offsets = final_offsets;
        last_idx = i;
      } else {
        break;
      }
    }
  } else {
    for (long i = (long)fns.size() - 1; i >= 0; i -= 1) {
      const std::string& fn = fns[i];
      const std::vector<long> final_offsets = check_file_sync_node(sfr, fn);
      if (0 != final_offsets.size()) {
        last_final_offsets = final_offsets;
        last_idx = i;
        break;
      }
    }
  }
  for (int i = 0; i < (int)sfr.frs.size(); ++i) {
    FieldsReader& fr = sfr.frs[i];
    const std::string path_file = fr.path_file();
    const long final_offset = last_final_offsets[i];
    fr.close();
    displayln_info(fname + ssprintf(": '%s' final_offset=%ld",
                                    path_file.c_str(), final_offset));
    const bool b = qtruncate(path_file, final_offset);
    qassert(b);
  }
  fns.resize(last_idx + 1);
  for (long i = 0; i < (long)fns.size(); ++i) {
    const std::string& fn = fns[i];
    displayln_info(fname + ssprintf(": i=%5ld fn='%s'", i, fn.c_str()));
  }
  displayln_info(fname +
                 ssprintf(": fns.size()=%5ld '%s'", fns.size(), path.c_str()));
  sync_node();
  return fns;
}

template <class M>
long read(ShuffledFieldsReader& sfr, const std::string& fn, Field<M>& field)
// interface function
{
  TIMER_VERBOSE_FLOPS("read(sfr,fn,field)");
  long total_bytes = 0;
  displayln_info(fname + ssprintf(": reading field with fn='%s' from '%s'.",
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
long read(ShuffledFieldsReader& sfr, const std::string& fn, const ShuffledBitSet& sbs, SelectedField<M>& sf)
// interface function
{
  TIMER_VERBOSE_FLOPS("read(sfr,fn,sbs,sf)");
  sf.init();
  long total_bytes = 0;
  displayln_info(fname +
                 ssprintf(": reading sparse field with fn='%s' from '%s'.",
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
                            const ShuffledBitSet& sbs, SelectedField<M>& sf)
// interface function
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

inline ShuffledFieldsReaderCache& get_shuffled_fields_reader_cache()
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
  displayln_info(fname +
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
