// vim: set ts=2 sw=2 expandtab:

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
  qassert(sz_compressed % bs.cN == 0);
  const size_t sz_block = sz_compressed / bs.cN;
  std::vector<char> ret(sz_block * local_volume, 0);
  bs.decompress(&data[nbytes], &ret[0], sz_block);
  return ret;
}

inline BitSet mk_bitset_from_field_rank(const FieldM<int64_t, 1>& f_rank,
                                        const int64_t n_per_tslice)
{
  TIMER("mk_bitset_from_field_rank");
  const Geometry& geo = f_rank.geo;
  BitSet bs(geo.local_volume());
  qassert(geo.is_only_local());
  for (long index = 0; index < geo.local_volume(); ++index) {
    const int64_t rank = f_rank.get_elem(index);
    if (0 <= rank and rank < n_per_tslice) {
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
  void init(const std::string& path_, const GeometryNode& geon_)
  {
    close();
    path = path_;
    geon = geon_;
    qassert(fp == NULL);
    if (geon.id_node == 0) {
      if (does_file_exist(path + ".partial")) {
        qremove_all(path + ".partial");
      }
      displayln("FieldsWriter: open '" + path + ".partial'.");
      fields_writer_dirs_geon_info(geon, path + ".partial");
    }
    fp = dist_open(path + ".partial", geon.id_node, geon.num_node, "w");
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

inline void write(FieldsWriter& fw, const std::string& fn, const Geometry& geo,
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
  // first read tag
  int32_t tag_len = 0;
  if (1 != fread_convert_endian(&tag_len, 4, 1, fr.fp, fr.is_little_endian)) {
    return false;
  }
  std::vector<char> fnv(tag_len);
  if (1 != fread(fnv.data(), tag_len, 1, fr.fp)) {
    qassert(false);
  }
  fn = std::string(fnv.data(), tag_len - 1);
  //
  // then read crc
  if (1 != fread_convert_endian(&crc, 4, 1, fr.fp, fr.is_little_endian)) {
    qassert(false);
  }
  //
  // then read geometry info
  int32_t nd = 0;
  if (1 != fread_convert_endian(&nd, 4, 1, fr.fp, fr.is_little_endian)) {
    qassert(false);
  }
  qassert(4 == nd);
  //
  std::vector<int32_t> gd(4, 0);
  std::vector<int32_t> num_procs(4, 0);
  if (4 != fread_convert_endian(&gd[0], 4, 4, fr.fp, fr.is_little_endian)) {
    qassert(false);
  }
  if (4 !=
      fread_convert_endian(&num_procs[0], 4, 4, fr.fp, fr.is_little_endian)) {
    qassert(false);
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
    qassert(size_node[mu] == (int)num_procs[mu]);
  }
  //
  // then read data size
  if (1 != fread_convert_endian(&data_len, 8, 1, fr.fp, fr.is_little_endian)) {
    qassert(false);
  }
  return true;
}

inline void read_data(FieldsReader& fr, std::vector<char>& data,
                      const int64_t data_len, const crc32_t crc)
{
  TIMER("read_data(fr,fn,geo,data)");
  clear(data);
  data.resize(data_len, 0);
  const long read_data_all = fread(&data[0], data_len, 1, fr.fp);
  qassert(1 == read_data_all);
  crc32_t crc_read = crc32_par(get_data(data));
  qassert(crc_read == crc);
}

inline bool read_next(FieldsReader& fr, std::string& fn, Coordinate& total_site,
                      std::vector<char>& data, bool& is_sparse_field)
{
  TIMER("read_next(fr,fn,geo,data)");
  crc32_t crc = 0;
  int64_t data_len = 0;
  const bool is_ok =
      read_tag(fr, fn, total_site, crc, data_len, is_sparse_field);
  if (is_ok) {
    read_data(fr, data, data_len, crc);
  }
  return is_ok;
}

inline bool read(FieldsReader& fr, const std::string& fn,
                 Coordinate& total_site, std::vector<char>& data,
                 bool& is_sparse_field)
{
  TIMER("read(fr,fn,geo,data)");
  bool has_restart = false;
  while (true) {
    std::string fn_read = "";
    crc32_t crc = 0;
    int64_t data_len = 0;
    const bool is_ok =
        read_tag(fr, fn_read, total_site, crc, data_len, is_sparse_field);
    if (is_ok) {
      if (fn == fn_read) {
        read_data(fr, data, data_len, crc);
        return true;
      } else {
        fseek(fr.fp, data_len, SEEK_CUR);
      }
    } else {
      if (has_restart) {
        return false;
      } else {
        rewind(fr.fp);
        has_restart = true;
      }
    }
  }
}

template <class M>
void write(FieldsWriter& fw, const std::string& fn, const Field<M>& field)
// field already have endianess converted correctly
{
  TIMER("write(fw,fn,field)");
  const Geometry& geo = field.geo;
  const Vector<M> v = get_data(field);
  const Vector<char> data((const char*)v.data(), v.data_size());
  write(fw, fn, geo, data, false);
}

template <class M>
void write(FieldsWriter& fw, const std::string& fn, const Field<M>& field,
           const BitSet& bs)
// field already have endianess converted correctly
{
  TIMER("write(fw,fn,field,bs)");
  const Geometry& geo = field.geo;
  const std::vector<char> data = bs.compress(get_data(field));
  write(fw, fn, geo, get_data(data), true);
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
bool read_next(FieldsReader& fr, std::string& fn, Field<M>& field)
// field endianess not converted at all
{
  TIMER("read_next(fw,fn,field)");
  Coordinate total_site;
  std::vector<char> data;
  bool is_sparse_field = false;
  const bool is_ok = read_next(fr, fn, total_site, data, is_sparse_field);
  if (not is_ok) {
    return false;
  }
  set_field_from_data(field, fr.geon, total_site, data, is_sparse_field);
  return true;
}

template <class M>
bool read(FieldsReader& fr, const std::string& fn, Field<M>& field)
// field endianess not converted at all
{
  TIMER("read(fw,fn,field)");
  Coordinate total_site;
  std::vector<char> data;
  bool is_sparse_field = false;
  const bool is_ok = read(fr, fn, total_site, data, is_sparse_field);
  if (not is_ok) {
    return false;
  }
  set_field_from_data(field, fr.geon, total_site, data, is_sparse_field);
  return true;
}

typedef std::vector<BitSet> ShuffledBitSet;

inline ShuffledBitSet mk_shuffled_bitset(const FieldM<int64_t, 1>& f_rank,
                                         const int64_t n_per_tslice,
                                         const Coordinate& new_size_node)
{
  TIMER("mk_shuffled_bitset(f_rank,n_per_tslice,new_size_node)");
  std::vector<Field<int64_t> > fs_rank;
  shuffle_field(fs_rank, f_rank, new_size_node);
  ShuffledBitSet sbs(fs_rank.size());
  for (int i = 0; i < (int)fs_rank.size(); ++i) {
    FieldM<int64_t, 1> fs_rank_i;
    fs_rank_i.init(fs_rank[i]);
    sbs[i] = mk_bitset_from_field_rank(fs_rank_i, n_per_tslice);
  }
  return sbs;
}

inline ShuffledBitSet mk_shuffled_bitset(const FieldSelection& fsel,
                                         const Coordinate& new_size_node)
{
  TIMER_VERBOSE("mk_shuffled_bitset(fsel,new_size_node)");
  return mk_shuffled_bitset(fsel.f_rank, fsel.n_per_tslice, new_size_node);
}

struct ShuffledFieldsWriter {
  std::string path;
  Coordinate new_size_node;
  std::vector<FieldsWriter> fws;
  //
  ShuffledFieldsWriter() { init(); }
  ShuffledFieldsWriter(const std::string& path_,
                       const Coordinate& new_size_node_)
  // interface function
  {
    init(path_, new_size_node_);
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
  void init(const std::string& path_, const Coordinate& new_size_node_)
  // interface function
  {
    init();
    path = path_;
    new_size_node = new_size_node_;
    std::vector<GeometryNode> geons = make_dist_io_geons(new_size_node);
    fws.resize(geons.size());
    for (int i = 0; i < (int)geons.size(); ++i) {
      if (geons[i].id_node == 0) {
        fws[i].init(path, geons[i]);
      }
    }
    sync_node();
    for (int i = 0; i < (int)geons.size(); ++i) {
      if (geons[i].id_node != 0) {
        fws[i].init(path, geons[i]);
      }
    }
  }
  //
  void close()
  // interface function
  {
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
  ShuffledFieldsReader(const std::string& path_)
  // interface function
  {
    init(path_);
  }
  //
  void init()
  // interface function
  {
    path = "";
    new_size_node = Coordinate();
    clear(frs);
  }
  void init(const std::string& path_)
  // interface function
  {
    init();
    path = path_;
    new_size_node = shuffled_fields_reader_size_node_info(path);
    std::vector<GeometryNode> geons = make_dist_io_geons(new_size_node);
    frs.resize(geons.size());
    for (int i = 0; i < (int)geons.size(); ++i) {
      frs[i].init(path, geons[i]);
    }
  }
};

template <class M>
void write(ShuffledFieldsWriter& sfw, const std::string& fn,
           const Field<M>& field)
// interface function
{
  TIMER_VERBOSE("write(sfw,fn,field)");
  displayln_info(fname +
                 ssprintf(": writting field with fn='%s'.", fn.c_str()));
  std::vector<Field<M> > fs;
  shuffle_field(fs, field, sfw.new_size_node);
  qassert(fs.size() == sfw.fws.size());
  for (int i = 0; i < (int)fs.size(); ++i) {
    write(sfw.fws[i], fn, fs[i]);
  }
}

template <class M>
void write(ShuffledFieldsWriter& sfw, const std::string& fn,
           const Field<M>& field, const ShuffledBitSet& sbs)
// interface function
{
  TIMER_VERBOSE("write(sfw,fn,field)");
  displayln_info(fname +
                 ssprintf(": writting sparse field with fn='%s'.", fn.c_str()));
  std::vector<Field<M> > fs;
  shuffle_field(fs, field, sfw.new_size_node);
  qassert(fs.size() == sfw.fws.size());
  qassert(sbs.size() == sfw.fws.size());
  for (int i = 0; i < (int)fs.size(); ++i) {
    write(sfw.fws[i], fn, fs[i], sbs[i]);
  }
}

template <class M>
bool read_next(ShuffledFieldsReader& sfr, std::string& fn, Field<M>& field)
// interface function
{
  TIMER_VERBOSE("read_next(sfr,fn,field)");
  fn = "";
  Coordinate total_site;
  int multiplicity = 0;
  bool is_ok = true;
  std::vector<Field<M> > fs(sfr.frs.size());
  for (int i = 0; i < (int)fs.size(); ++i) {
    std::string fni;
    is_ok = is_ok and read_next(sfr.frs[i], fni, fs[i]);
    if (sfr.frs[i].geon.id_node == 0) {
      fn == fni;
      total_site = fs[i].geo.total_site();
      multiplicity = fs[i].geo.multiplicity;
      qassert(get_id_node() == 0);
    }
  }
  bcast(get_data_one_elem(is_ok));
  if (not is_ok) {
    fn = "";
    return false;
  }
  displayln_info(fname +
                 ssprintf(": read the next field with fn='%s'.", fn.c_str()));
  bcast(fn);
  bcast(get_data_one_elem(total_site));
  bcast(get_data_one_elem(multiplicity));
  Geometry geo;
  geo.init(total_site, multiplicity);
  field.init(geo);
  shuffle_field_back(field, fs, sfr.new_size_node);
  return true;
}

template <class M>
bool read(ShuffledFieldsReader& sfr, const std::string& fn, Field<M>& field)
// interface function
{
  TIMER_VERBOSE("read(sfr,fn,field)");
  Coordinate total_site;
  int multiplicity = 0;
  bool is_ok = true;
  displayln_info(fname + ssprintf(": reading field with fn='%s'.", fn.c_str()));
  std::vector<Field<M> > fs(sfr.frs.size());
  for (int i = 0; i < (int)fs.size(); ++i) {
    is_ok = is_ok and read(sfr.frs[i], fn, fs[i]);
    if (sfr.frs[i].geon.id_node == 0) {
      total_site = fs[i].geo.total_site();
      multiplicity = fs[i].geo.multiplicity;
      qassert(get_id_node() == 0);
    }
  }
  bcast(get_data_one_elem(is_ok));
  if (not is_ok) {
    return false;
  }
  bcast(get_data_one_elem(total_site));
  bcast(get_data_one_elem(multiplicity));
  Geometry geo;
  geo.init(total_site, multiplicity);
  field.init(geo);
  shuffle_field_back(field, fs, sfr.new_size_node);
  return true;
}

}  // namespace qlat
