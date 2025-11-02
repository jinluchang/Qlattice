#define QLAT_INSTANTIATE_FIELDS_IO

#include <qlat/fields-io.h>

namespace qlat
{  //

void BitSet::set(const void* pbytes, const size_t nbytes)
{
  Qassert(nbytes == bytes.size());
  memcpy(&bytes[0], pbytes, nbytes);
  cN = 0;
  for (size_t i = 0; i < N; i++) {
    if (get(i)) {
      cN++;
    }
  }
}

void BitSet::set_f_rank(FieldRank& f_rank, const int64_t rank)
{
  TIMER("BitSet::set_f_rank")
  const Geometry& geo = f_rank.geo();
  const Int multiplicity = f_rank.multiplicity;
  Qassert(geo.local_volume() == (Long)N);
  Qassert(geo.is_only_local);
  Qassert(multiplicity == 1);
  for (size_t i = 0; i < N; i++) {
    if (get(i)) {
      f_rank.get_elem(i) = rank;
    } else {
      f_rank.get_elem(i) = -1;
    }
  }
}

bool BitSet::check_f_rank(const FieldRank& f_rank)
{
  TIMER("BitSet::check_f_rank")
  const Geometry& geo = f_rank.geo();
  const Int multiplicity = f_rank.multiplicity;
  Qassert(geo.is_only_local);
  Qassert(multiplicity == 1);
  if (not(geo.local_volume() == (Long)N)) {
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

void BitSet::compress(const void* src, void* dst, size_t block_size) const
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

void BitSet::decompress(const void* src, void* dst, size_t block_size) const
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

std::vector<char> bitset_decompress(const std::vector<char>& data,
                                    const Long local_volume)
{
  TIMER("bitset_decompress");
  const size_t N = local_volume;
  const size_t nbytes = 1 + (N - 1) / 8;
  BitSet bs(N);
  bs.set(&data[0], nbytes);
  const size_t sz_compressed = data.size() - nbytes;
  size_t sz_block = 0;
  if (bs.cN == 0) {
    Qassert(sz_compressed == 0);
    std::vector<char> ret;
    return ret;
  } else {
    Qassert(sz_compressed % bs.cN == 0);
    sz_block = sz_compressed / bs.cN;
  }
  std::vector<char> ret(sz_block * local_volume, 0);
  bs.decompress(&data[nbytes], &ret[0], sz_block);
  return ret;
}

BitSet mk_bitset_from_field_rank(const FieldRank& f_rank)
{
  TIMER("mk_bitset_from_field_rank");
  const Geometry& geo = f_rank.geo();
  BitSet bs(geo.local_volume());
  Qassert(geo.is_only_local);
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const int64_t rank = f_rank.get_elem(index);
    if (0 <= rank) {
      bs.set(index, true);
    } else {
      bs.set(index, false);
    }
  }
  return bs;
}

void fields_writer_dirs_geon_info(const GeometryNode& geon,
                                  const std::string& path, const mode_t mode)
{
  TIMER("fields_writer_dirs_geon_info");
  dist_mkdir(path, geon.num_node, mode);
  const std::string fn = path + "/geon-info.txt";
  std::string contents;
  contents += ssprintf("geon.num_node = %d\n", geon.num_node);
  contents += ssprintf("geon.size_node[0] = %d\n", geon.size_node[0]);
  contents += ssprintf("geon.size_node[1] = %d\n", geon.size_node[1]);
  contents += ssprintf("geon.size_node[2] = %d\n", geon.size_node[2]);
  contents += ssprintf("geon.size_node[3] = %d\n", geon.size_node[3]);
  qtouch(fn, contents);
}

Coordinate shuffled_fields_reader_size_node_info(const std::string& path)
{
  TIMER("shuffled_fields_reader_size_node_info");
  Coordinate size_node;
  Coordinate node_site;
  if (get_id_node() == 0) {
    const std::string fn = path + "/geon-info.txt";
    const std::vector<std::string> lines = qgetlines(fn);
    Int num_node;
    reads(num_node, info_get_prop(lines, "geon.num_node = "));
    for (Int i = 0; i < 4; ++i) {
      reads(size_node[i],
            info_get_prop(lines, ssprintf("geon.size_node[%d] = ", i)));
    }
    Qassert(num_node == product(size_node));
  }
  bcast(get_data(size_node));
  return size_node;
}

void FieldsWriter::init()
{
  path = "";
  geon.init();
  qfile.init();
  is_little_endian = true;
  fn_list.clear();
  offsets_map.clear();
  max_offset = 0;
}

void FieldsWriter::init(const std::string& path_, const GeometryNode& geon_,
                        const bool is_append)
{
  path = path_;
  geon = geon_;
  qfile.init();
  if (geon.id_node == 0) {
    displayln(0, "FieldsWriter: open '" + path + "'.");
    if (is_append and does_file_exist(path)) {
      if (not does_file_exist(path + "/geon-info.txt")) {
        qwarn(
            "FieldsWriter::init: " +
            ssprintf("geon-info.txt does not exist! Will create again. path='%s'", path.c_str()));
        fields_writer_dirs_geon_info(geon, path);
      }
    } else {
      fields_writer_dirs_geon_info(geon, path);
    }
  }
  qfile = qfopen(dist_file_name(path, geon.id_node, geon.num_node),
                 is_append ? "a" : "w");
  Qassert(not qfile.null());
}

void FieldsReader::init()
{
  path = "";
  geon.init();
  qfile.init();
  is_little_endian = true;
  is_read_through = false;
  fn_list.clear();
  offsets_map.clear();
  max_offset = 0;
  file_size = 0;
}

void FieldsReader::init(const std::string& path_, const GeometryNode& geon_)
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
    file_size = 0;
  } else {
    is_read_through = false;
    file_size = qfile.size();
  }
  fn_list.clear();
  offsets_map.clear();
  max_offset = 0;
}

void mkfile(FieldsReader& fr)
// create the file (open with appending and then close)
// does not open the file
{
  if (fr.qfile.null() and fr.path != "") {
    fr.qfile =
        qfopen(dist_file_name(fr.path, fr.geon.id_node, fr.geon.num_node), "a");
    qfclose(fr.qfile);
  }
}

std::string get_file_path(FieldsReader& fr)
{
  if (fr.path != "") {
    return dist_file_name(fr.path, fr.geon.id_node, fr.geon.num_node);
  }
  return "";
}

Long qfwrite_convert_endian(void* ptr, const size_t size, const size_t nmemb,
                            QFile& qfile, const bool is_little_endian)
{
  if (size == 4) {
    convert_endian(Vector<int32_t>((int32_t*)ptr, nmemb), is_little_endian);
  } else if (size == 8) {
    convert_endian(Vector<int64_t>((int64_t*)ptr, nmemb), is_little_endian);
  } else {
    Qassert(false);
  }
  const Long ret = qfwrite(ptr, size, nmemb, qfile);
  if (size == 4) {
    convert_endian(Vector<int32_t>((int32_t*)ptr, nmemb), is_little_endian);
  } else if (size == 8) {
    convert_endian(Vector<int64_t>((int64_t*)ptr, nmemb), is_little_endian);
  } else {
    Qassert(false);
  }
  return ret;
}

Long write(FieldsWriter& fw, const std::string& fn, const Geometry& geo,
           const Vector<char> data, const bool is_sparse_field)
{
  TIMER("write(fw,fn,geo,data)");
  Long n_elem_write;
  Long offset;
  // get initial offset
  const Long offset_start = fw.qfile.tell();
  offset = offset_start;
  // first write tag
  int32_t tag_len = fn.size() + 1;  // fn is the name of the field (say prop1)
  n_elem_write =
      qfwrite_convert_endian(&tag_len, 4, 1, fw.qfile, fw.is_little_endian);
  Qassert(n_elem_write == 1);
  offset += 4;
  n_elem_write = qfwrite(fn.c_str(), tag_len, 1, fw.qfile);
  Qassert(n_elem_write == 1);
  offset += tag_len;
  //
  // then write crc
  crc32_t crc = crc32_par(data);
  n_elem_write =
      qfwrite_convert_endian(&crc, 4, 1, fw.qfile, fw.is_little_endian);
  Qassert(n_elem_write == 1);
  offset += 4;
  //
  // then write geometry info
  int32_t nd = 4;  // <- number of dimensions of field, typically 4
  n_elem_write =
      qfwrite_convert_endian(&nd, 4, 1, fw.qfile, fw.is_little_endian);
  Qassert(n_elem_write == 1);
  offset += 4;
  //
  std::vector<int32_t> gd(4, 0);
  std::vector<int32_t> num_procs(4, 0);
  //
  const Coordinate total_site = geo.total_site();
  const Coordinate size_node = geo.geon.size_node;
  for (Int mu = 0; mu < 4; ++mu) {
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
  n_elem_write =
      qfwrite_convert_endian(&gd[0], 4, nd, fw.qfile, fw.is_little_endian);
  Qassert(n_elem_write == nd);
  offset += 4 * nd;
  n_elem_write = qfwrite_convert_endian(&num_procs[0], 4, nd, fw.qfile,
                                        fw.is_little_endian);
  Qassert(n_elem_write == nd);
  offset += 4 * nd;
  //
  // then data size
  int64_t data_len = data.size();
  n_elem_write =
      qfwrite_convert_endian(&data_len, 8, 1, fw.qfile, fw.is_little_endian);
  Qassert(n_elem_write == 1);
  offset += 8;
  //
  // then write data
  n_elem_write = qfwrite(&data[0], data_len, 1, fw.qfile);
  Qassert(n_elem_write == 1);
  offset += data_len;
  //
  const Long offset_stop = fw.qfile.tell();
  Qassert(offset_stop == offset);
  //
  // register file
  fw.fn_list.push_back(fn);
  fw.offsets_map[fn] =
      FieldsSegmentInfo(offset_start, offset_stop, is_sparse_field);
  //
  return data_len;
}

Long qfread_convert_endian(void* ptr, const size_t size, const size_t nmemb,
                           QFile& qfile, const bool is_little_endian)
{
  if (qfile.null()) {
    return 0;
  }
  const Long total_nmemb = qfread(ptr, size, nmemb, qfile);
  if (size == 4) {
    convert_endian(Vector<int32_t>((int32_t*)ptr, nmemb), is_little_endian);
  } else if (size == 8) {
    convert_endian(Vector<int64_t>((int64_t*)ptr, nmemb), is_little_endian);
  } else {
    Qassert(false);
  }
  return total_nmemb;
}

bool read_tag(FieldsReader& fr, std::string& fn, Coordinate& total_site,
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
    qwarn(ssprintf("read_tag: fr.qfile.null()==true fn='%s'",
                   get_file_path(fr).c_str()));
    return false;
  }
  //
  const Long offset_initial = qftell(fr.qfile);
  //
  // first read tag
  int32_t tag_len = 0;
  if (1 !=
      qfread_convert_endian(&tag_len, 4, 1, fr.qfile, fr.is_little_endian)) {
    fr.is_read_through = true;
    return false;
  }
  if (not(tag_len > 0)) {
    qwarn(
        ssprintf("read_tag: tag_len <= 0 fn='%s'", get_file_path(fr).c_str()));
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
  if (not(4 == nd)) {
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
  if (4 != qfread_convert_endian(&num_procs[0], 4, 4, fr.qfile,
                                 fr.is_little_endian)) {
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
  for (Int mu = 0; mu < 4; ++mu) {
    total_site[mu] = gd[mu];
    if (not(size_node[mu] == (int)num_procs[mu])) {
      qwarn(ssprintf("read_tag: fn='%s'", get_file_path(fr).c_str()));
      fr.is_read_through = true;
      return false;
    }
  }
  //
  // then read data size
  if (1 !=
      qfread_convert_endian(&data_len, 8, 1, fr.qfile, fr.is_little_endian)) {
    qwarn(ssprintf("read_tag: fn='%s'", get_file_path(fr).c_str()));
    fr.is_read_through = true;
    return false;
  }
  if (not(data_len > 0)) {
    qwarn(ssprintf("read_tag: fn='%s'", get_file_path(fr).c_str()));
    fr.is_read_through = true;
    return false;
  }
  //
  const Long offset_final = qftell(fr.qfile) + data_len;
  //
  const FieldsSegmentInfo fsinfo = FieldsSegmentInfo(offset_initial, offset_final, is_sparse_field);
  //
  if (has(fr.offsets_map, fn)) {
    if (fr.offsets_map[fn].offset_start != offset_initial) {
      qwarn(
          ssprintf("read_tag: fn='%s' of '%s' appeared twice! recorded=%ld ; "
                   "newly read=%ld ; id_node=%d.",
                   fn.c_str(), fr.path.c_str(), fr.offsets_map[fn].offset_start,
                   offset_initial, fr.geon.id_node));
      if (offset_initial > fr.offsets_map[fn].offset_start) {
        fr.fn_list.push_back(fn);
        fr.offsets_map[fn] = fsinfo;
      }
    }
  } else {
    fr.fn_list.push_back(fn);
    fr.offsets_map[fn] = fsinfo;
  }
  //
  if (offset_final > fr.max_offset) {
    fr.max_offset = offset_final;
  }
  //
  if (fr.geon.id_node == 0) {
    displayln_info(
        0, fname + ssprintf(": '%s' from '%s'.", fn.c_str(), fr.path.c_str()));
  }
  return true;
}

Long read_data(FieldsReader& fr, std::vector<char>& data,
               const int64_t data_len, const crc32_t crc)
// return data_len (if not successful then return 0)
{
  TIMER_FLOPS("read_data(fr,fn,geo,data)");
  clear(data);
  data.resize(data_len, 0);
  if (fr.qfile.null()) {
    qwarn(ssprintf("read_data: file does not exist fn='%s'",
                   get_file_path(fr).c_str()));
    return 0;
  }
  const Long read_data_all = qfread(&data[0], data_len, 1, fr.qfile);
  if (not(1 == read_data_all)) {
    qwarn(ssprintf("read_data: data not complete fn='%s'",
                   get_file_path(fr).c_str()));
    fr.is_read_through = true;
    return 0;
  }
  crc32_t crc_read = crc32_par(get_data(data));
  if (not(crc_read == crc)) {
    qwarn(ssprintf("read_data: crc does not match fn='%s'",
                   get_file_path(fr).c_str()));
    return 0;
  }
  timer.flops += data_len;
  return data_len;
}

Long read_next(FieldsReader& fr, std::string& fn, Coordinate& total_site,
               std::vector<char>& data, bool& is_sparse_field)
{
  TIMER_FLOPS("read_next(fr,fn,geo,data)");
  crc32_t crc = 0;
  int64_t data_len = 0;
  const bool is_ok =
      read_tag(fr, fn, total_site, crc, data_len, is_sparse_field);
  const Long total_bytes = is_ok ? read_data(fr, data, data_len, crc) : 0;
  timer.flops += total_bytes;
  return total_bytes;
}

Long read_skip_next(FieldsReader& fr, std::string& fn)
// return offset of the end of this data segment
// return -1 if failed.
{
  Coordinate total_site;
  bool is_sparse_field;
  fn = "";
  crc32_t crc = 0;
  int64_t data_len = 0;
  const bool is_ok =
      read_tag(fr, fn, total_site, crc, data_len, is_sparse_field);
  if (not is_ok) {
    errno = 0;
    fn = "";
    return -1;
  }
  const Int ret = qfseek(fr.qfile, data_len, SEEK_CUR);
  if (ret != 0) {
    errno = 0;
    return -1;
  }
  const Long offset_final = fr.qfile.tell();
  if (offset_final > fr.file_size) {
    return -1;
  }
  return offset_final;
}

void read_through(FieldsReader& fr)
{
  TIMER("read_through(fr)");
  Qassert(not fr.qfile.null());
  Qassert(fr.max_offset <= fr.file_size);
  const Int code = qfseek(fr.qfile, fr.max_offset, SEEK_SET);
  Qassert(code == 0);
  while (true) {
    std::string fn_read = "";
    const Long offset_final = read_skip_next(fr, fn_read);
    if (offset_final < 0) {
      Qassert(offset_final == -1);
      fr.is_read_through = true;
      return;
    }
  }
}

bool does_file_exist(const FieldsReader& fr, const std::string& fn)
{
  TIMER("does_file_exist(fr,fn,site)");
  Qassert(fr.is_read_through);
  return has(fr.offsets_map, fn);
}

bool does_file_exist(const FieldsWriter& fw, const std::string& fn)
{
  TIMER("does_file_exist(fw,fn,site)");
  return has(fw.offsets_map, fn);
}

Long read(FieldsReader& fr, const std::string& fn, Coordinate& total_site,
          std::vector<char>& data, bool& is_sparse_field)
{
  TIMER_FLOPS("read(fr,fn,site,data)");
  if (not does_file_exist(fr, fn)) {
    return 0;
  }
  Qassert(has(fr.offsets_map, fn));
  qfseek(fr.qfile, fr.offsets_map[fn].offset_start, SEEK_SET);
  std::string fn_r;
  const Long total_bytes =
      read_next(fr, fn_r, total_site, data, is_sparse_field);
  Qassert(fn == fn_r);
  return total_bytes;
}

Long read_skip(FieldsReader& fr, const std::string& fn)
// return offset of the end of this data segment
// return -1 if failed.
{
  TIMER_FLOPS("read_skip(fr,fn)");
  if (not does_file_exist(fr, fn)) {
    return -1;
  }
  Qassert(has(fr.offsets_map, fn));
  qfseek(fr.qfile, fr.offsets_map[fn].offset_start, SEEK_SET);
  std::string fn_r;
  const Long offset_final = read_skip_next(fr, fn_r);
  if (fn != fn_r) {
    return -1;
  }
  return offset_final;
}

Long check_file(FieldsReader& fr, const std::string& fn, const bool is_check_data)
// return final offset of the data
// if check_file fail, return -1
{
  TIMER_FLOPS("check_file(fr,fn)");
  if (not is_check_data) {
    return read_skip(fr, fn);
  }
  Coordinate total_site;
  std::vector<char> data;
  bool is_sparse_field;
  const Long total_bytes = read(fr, fn, total_site, data, is_sparse_field);
  if (total_bytes > 0) {
    return qftell(fr.qfile);
  } else {
    return -1;
  }
}

Int flush(FieldsWriter& fw)
{
  TIMER("flush(fw)");
  return qfflush(fw.qfile);
}

// ------------------------

ShuffledBitSet mk_shuffled_bitset(const FieldRank& f_rank,
                                  const Coordinate& new_size_node)
{
  TIMER("mk_shuffled_bitset(f_rank,new_size_node)");
  std::vector<Field<int64_t> > fs_rank;
  shuffle_field(fs_rank, f_rank, new_size_node);
  ShuffledBitSet sbs;
  set_field_selection(sbs.fsel, f_rank);
  sbs.sp = make_shuffle_plan(sbs.fsels, sbs.fsel, new_size_node);
  sbs.vbs.resize(fs_rank.size());
  for (Int i = 0; i < (int)fs_rank.size(); ++i) {
    FieldRank fs_rank_i;
    fs_rank_i.init(fs_rank[i]);
    sbs.vbs[i] = mk_bitset_from_field_rank(fs_rank_i);
  }
  return sbs;
}

ShuffledBitSet mk_shuffled_bitset(const FieldSelection& fsel,
                                  const Coordinate& new_size_node)
// interface function
{
  TIMER_VERBOSE("mk_shuffled_bitset(fsel,new_size_node)");
  return mk_shuffled_bitset(fsel.f_rank, new_size_node);
}

ShuffledBitSet mk_shuffled_bitset(const Coordinate& total_site,
                                  const PointsSelection& psel,
                                  const Coordinate& new_size_node)
{
  TIMER("mk_shuffled_bitset(total_site,psel,new_size_node)");
  FieldRank f_rank;
  mk_field_selection(f_rank, total_site, psel);
  return mk_shuffled_bitset(f_rank, new_size_node);
}

ShuffledBitSet mk_shuffled_bitset(const FieldRank& f_rank,
                                  const PointsSelection& psel,
                                  const Coordinate& new_size_node)
{
  TIMER_VERBOSE("mk_shuffled_bitset(f_rank,psel,new_size_node)");
  FieldRank f_rank_combined;
  f_rank_combined = f_rank;
  add_field_selection(f_rank_combined, psel);
  return mk_shuffled_bitset(f_rank_combined, new_size_node);
}

// ------------------------

void ShuffledFieldsWriter::init()
// interface function
{
  close();
  path = "";
  new_size_node = Coordinate();
}

void ShuffledFieldsWriter::init(const std::string& path_,
                                const Coordinate& new_size_node_,
                                const bool is_append,
                                const bool is_removing_old)
// interface function
{
  TIMER_VERBOSE("ShuffledFieldsWriter::init(p,sn,app)")
  remove_entry_directory_cache(path_);
  init();
  path = path_;
  new_size_node = new_size_node_;
  std::vector<std::string> fn_list;
  std::vector<std::vector<FieldsSegmentInfo>> offsets_list;
  if (is_append) {
    if (does_file_exist_sync_node(path + "/geon-info.txt")) {
      new_size_node = shuffled_fields_reader_size_node_info(path);
      if (new_size_node_ != Coordinate() and new_size_node_ != new_size_node) {
        if (get_id_node() == 0) {
          qwarn(
              fname +
              ssprintf(
                  ": WARNING: new_size_node do not match. file=%s argument=%s "
                  ". Will use the new_size_node from the existing file.",
                  show(new_size_node).c_str(), show(new_size_node_).c_str()));
        }
      }
    } else {
      if (does_file_exist_sync_node(path)) {
        if (get_id_node() == 0) {
          qwarn(fname +
                ssprintf(": WARNING: path='%s' exists but "
                         "'geon-info.txt' is not present in the folder.",
                         path.c_str()));
        }
      }
    }
    properly_truncate_fields_sync_node(fn_list, offsets_list, path, false,
                                       false, new_size_node);
    if (get_id_node() == 0) {
      qar_index.init(path + "/index.qar", QFileMode::Append);
    }
  } else {
    if (does_file_exist_sync_node(path + "/geon-info.txt")) {
      if (is_removing_old) {
        if (get_id_node() == 0) {
          qwarn(fname + ssprintf(": 'geon-info.txt' exist. Removing '%s'.",
                                 path.c_str()));
        }
        qremove_all_sync_node(path);
      } else {
        qerr(fname +
             ssprintf(": cannot open for write '%s/geon-info.txt' exist",
                      path.c_str()));
      }
    }
    if (does_file_exist_sync_node(path)) {
      if (is_removing_old) {
        if (get_id_node() == 0) {
          qwarn(fname +
                ssprintf(": directory exist. Removing '%s'.", path.c_str()));
        }
        qremove_all_sync_node(path);
      } else {
        qerr(fname +
             ssprintf(": cannot open for write '%s' exist", path.c_str()));
      }
    }
    if (get_id_node() == 0) {
      qar_index.init(path + "/index.qar", QFileMode::Write);
    }
  }
  if (get_id_node() == 0) {
    qar_index.flush();
  }
  qar_index_idx = list(qar_index).size();
  std::vector<GeometryNode> geons = make_dist_io_geons(new_size_node);
  fws.resize(geons.size());
  for (Int i = 0; i < (int)geons.size(); ++i) {
    const GeometryNode& geon = geons[i];
    if (geon.id_node == 0) {
      fws[i].init(path, geon, is_append);
    }
  }
  SYNC_NODE();
  for (Int i = 0; i < (int)geons.size(); ++i) {
    const GeometryNode& geon = geons[i];
    if (geon.id_node != 0) {
      fws[i].init(path, geon, is_append);
    }
  }
  for (Int i = 0; i < (int)fws.size(); ++i) {
    FieldsWriter& fw = fws[i];
    fw.fn_list = fn_list;
    for (Long j = 0; j < (Long)fn_list.size(); ++j) {
      const std::string& fn = fn_list[j];
      fw.offsets_map[fn] = offsets_list[j][i];
    }
  }
  add_shuffled_fields_writer(*this);
}

void ShuffledFieldsWriter::close()
// interface function
{
  remove_shuffled_fields_writer(*this);
  if (fws.size() > 0) {
    TIMER_VERBOSE("ShuffledFieldsWriter::close");
    for (Long i = 0; i < (Long)fws.size(); ++i) {
      fws[i].close();
    }
    clear(fws);
  }
  qar_index.close();
  qar_index_idx = 0;
}

void ShuffledFieldsReader::init()
// interface function
{
  close();
  path = "";
  new_size_node = Coordinate();
}

void ShuffledFieldsReader::init(const std::string& path_,
                                const Coordinate& new_size_node_)
// interface function
{
  TIMER_VERBOSE("ShuffledFieldsReader::init(path,new_size_node)")
  init();
  path = path_;
  if (does_file_exist_qar_sync_node(path + "/geon-info.txt")) {
    new_size_node = shuffled_fields_reader_size_node_info(path);
  } else {
    Qassert(new_size_node_ != Coordinate());
    new_size_node = new_size_node_;
  }
  std::vector<GeometryNode> geons = make_dist_io_geons(new_size_node);
  frs.resize(geons.size());
  for (Int i = 0; i < (int)geons.size(); ++i) {
    frs[i].init(path, geons[i]);
  }
  if (does_file_exist_qar_sync_node(path + "/index.qar")) {
    QarFile qar_index;
    if (get_id_node() == 0) {
      qar_index.init(path + "/index.qar", QFileMode::Read);
    }
    std::vector<std::string> fn_list;
    std::vector<std::vector<FieldsSegmentInfo>> all_offsets_list;
    load_all_fields_index(fn_list, all_offsets_list, qar_index);
    qar_index.close();
    for (Long i = 0; i < (Long)fn_list.size(); ++i) {
      const std::string& fn = fn_list[i];
      const std::vector<FieldsSegmentInfo>& all_offsets = all_offsets_list[i];
      populate_fields_offsets(*this, fn, all_offsets);
    }
  }
  read_through_sync_node(*this);
  add_shuffled_fields_reader(*this);
}

void ShuffledFieldsReader::close()
// interface function
{
  remove_shuffled_fields_reader(*this);
  if (frs.size() > 0) {
    TIMER_VERBOSE("ShuffledFieldsReader::close")
    for (Long i = 0; i < (Long)frs.size(); ++i) {
      frs[i].close();
    }
    clear(frs);
  }
}

// ------------------------

std::string show(const ShuffledFieldsWriter& sfw)
{
  return ssprintf("ShuffledFieldsWriter(path='%s',new_size_node=%s)",
                  sfw.path.c_str(), show(sfw.new_size_node).c_str());
}

void add_shuffled_fields_writer(ShuffledFieldsWriter& sfw)
{
  ShuffledFieldsWriterMap& sfwm = get_all_shuffled_fields_writer();
  const Long key = (Long)&sfw;
  Qassert(not has(sfwm, key));
  sfwm[key] = Handle<ShuffledFieldsWriter>(sfw);
}

void remove_shuffled_fields_writer(ShuffledFieldsWriter& sfw)
{
  ShuffledFieldsWriterMap& sfwm = get_all_shuffled_fields_writer();
  const Long key = (Long)&sfw;
  if (has(sfwm, key)) {
    sfwm.erase(key);
  }
}

void close_all_shuffled_fields_writer()
// Force close all the ShuffledFieldsWriter.
// Only call this when quitting the program (e.g. in qquit(msg)).
{
  TIMER_VERBOSE("close_all_shuffled_fields_writer");
  ShuffledFieldsWriterMap& sfwm = get_all_shuffled_fields_writer();
  std::vector<Handle<ShuffledFieldsWriter>> sfwv;
  for (auto it = sfwm.cbegin(); it != sfwm.cend(); ++it) {
    sfwv.push_back(it->second);
  }
  for (Long i = 0; i < (Long)sfwv.size(); ++i) {
    sfwv[i]().close();
  }
  Qassert(sfwm.size() == 0);
  SYNC_NODE();
}

std::vector<std::string> show_all_shuffled_fields_writer()
{
  std::vector<std::string> ret;
  const ShuffledFieldsWriterMap& sfwm = get_all_shuffled_fields_writer();
  for (auto it = sfwm.cbegin(); it != sfwm.cend(); ++it) {
    const ShuffledFieldsWriter& sfw = (it->second)();
    ret.push_back(show(sfw));
  }
  return ret;
}

// ------------------------

std::string show(const ShuffledFieldsReader& sfr)
{
  return ssprintf("ShuffledFieldsReader(path='%s',new_size_node=%s)",
                  sfr.path.c_str(), show(sfr.new_size_node).c_str());
}

void add_shuffled_fields_reader(ShuffledFieldsReader& sfr)
{
  ShuffledFieldsReaderMap& sfrm = get_all_shuffled_fields_reader();
  const Long key = (Long)&sfr;
  Qassert(not has(sfrm, key));
  sfrm[key] = Handle<ShuffledFieldsReader>(sfr);
}

void remove_shuffled_fields_reader(ShuffledFieldsReader& sfr)
{
  ShuffledFieldsReaderMap& sfrm = get_all_shuffled_fields_reader();
  const Long key = (Long)&sfr;
  if (has(sfrm, key)) {
    sfrm.erase(key);
  }
}

void close_all_shuffled_fields_reader()
// Force close all the ShuffledFieldsReader.
// Only call this when quitting the program (e.g. in qquit(msg)).
{
  TIMER_VERBOSE("close_all_shuffled_fields_reader");
  ShuffledFieldsReaderMap& sfrm = get_all_shuffled_fields_reader();
  std::vector<Handle<ShuffledFieldsReader>> sfrv;
  for (auto it = sfrm.cbegin(); it != sfrm.cend(); ++it) {
    sfrv.push_back(it->second);
  }
  for (Long i = 0; i < (Long)sfrv.size(); ++i) {
    sfrv[i]().close();
  }
  Qassert(sfrm.size() == 0);
  SYNC_NODE();
}

std::vector<std::string> show_all_shuffled_fields_reader()
{
  std::vector<std::string> ret;
  const ShuffledFieldsReaderMap& sfrm = get_all_shuffled_fields_reader();
  for (auto it = sfrm.cbegin(); it != sfrm.cend(); ++it) {
    const ShuffledFieldsReader& sfr = (it->second)();
    ret.push_back(show(sfr));
  }
  return ret;
}

// ------------------------

Long flush(ShuffledFieldsWriter& sfw)
// interface function
{
  TIMER_VERBOSE("flush(sfw)");
  Long ret = 0;
  for (Int i = 0; i < (int)sfw.fws.size(); ++i) {
    ret += flush(sfw.fws[i]);
  }
  if (get_id_node() == 0) {
    sfw.qar_index.flush();
  }
  glb_sum(ret);
  return ret;
}

void read_through_sync_node(ShuffledFieldsReader& sfr)
{
  TIMER_VERBOSE("read_through_sync_node(sfr)");
  for (Int i = 0; i < (int)sfr.frs.size(); ++i) {
    read_through(sfr.frs[i]);
  }
  SYNC_NODE();
}

bool does_file_exist_sync_node(const ShuffledFieldsReader& sfr, const std::string& fn)
// interface function
{
  TIMER("does_file_exist_sync_node(sfr,fn)");
  Long total_counts = 0;
  displayln_info(1, fname + ssprintf(": check fn='%s' from '%s'.", fn.c_str(),
                                     sfr.path.c_str()));
  for (Int i = 0; i < (int)sfr.frs.size(); ++i) {
    if (does_file_exist(sfr.frs[i], fn)) {
      total_counts += 1;
    }
  }
  glb_sum(total_counts);
  if (total_counts == 0) {
    return false;
  } else {
    Qassert(total_counts == product(sfr.new_size_node));
    return true;
  }
}

bool does_file_exist_sync_node(const ShuffledFieldsWriter& sfw,
                               const std::string& fn)
// interface function
{
  TIMER("does_file_exist_sync_node(sfw,fn)");
  Long total_counts = 0;
  displayln_info(1, fname + ssprintf(": check fn='%s' from '%s'.", fn.c_str(),
                                     sfw.path.c_str()));
  for (Int i = 0; i < (int)sfw.fws.size(); ++i) {
    if (does_file_exist(sfw.fws[i], fn)) {
      total_counts += 1;
    }
  }
  glb_sum(total_counts);
  if (total_counts == 0) {
    return false;
  } else {
    Qassert(total_counts == product(sfw.new_size_node));
    return true;
  }
}

bool is_sparse_field_sync_node(const ShuffledFieldsReader& sfr,
                               const std::string& fn)
// interface function
{
  TIMER("is_sparse_field_sync_node(sfr,fn)");
  Long total_counts = 0;
  displayln_info(1, fname + ssprintf(": check fn='%s' from '%s'.", fn.c_str(),
                                     sfr.path.c_str()));
  for (Int i = 0; i < (int)sfr.frs.size(); ++i) {
    const FieldsReader& fr = sfr.frs[i];
    if (has(fr.offsets_map, fn)) {
      const FieldsSegmentInfo& fsinfo = fr.offsets_map.at(fn);
      if (fsinfo.is_sparse_field) {
        total_counts += 1;
      }
    }
  }
  glb_sum(total_counts);
  if (total_counts == 0) {
    return false;
  } else {
    Qassert(total_counts == product(sfr.new_size_node));
    return true;
  }
}

bool check_file_sync_node(ShuffledFieldsReader& sfr, const std::string& fn,
                          const bool is_check_data,
                          std::vector<Long>& final_offsets)
// interface function
// set final_offsets to be the files position after loading the data ``fn''
// (zero if failed for that file) return if data is loaded successfully
{
  TIMER_VERBOSE("check_file_sync_node(sfr,fn,is_check_data,final_offsets)");
  displayln_info(0, fname + ssprintf(": reading field with fn='%s' from '%s'.",
                                     fn.c_str(), sfr.path.c_str()));
  clear(final_offsets);
  final_offsets.resize(sfr.frs.size(), 0);
  Long total_failed_counts = 0;
  for (Int i = 0; i < (int)sfr.frs.size(); ++i) {
    final_offsets[i] = check_file(sfr.frs[i], fn, is_check_data);
    if (final_offsets[i] == -1) {
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

std::vector<std::string> list_fields(const ShuffledFieldsReader& sfr, bool is_skipping_check)
// interface function
{
  TIMER_VERBOSE("list_fields(sfr)");
  std::vector<std::string> ret;
  if (0 == get_id_node()) {
    Qassert(sfr.frs.size() > 0);
    const FieldsReader& fr = sfr.frs[0];
    Qassert(fr.is_read_through);
    ret = fr.fn_list;
  }
  bcast(ret);
  if (not is_skipping_check) {
    for (Int i = 0; i < (int)sfr.frs.size(); ++i) {
      const FieldsReader& fr = sfr.frs[i];
      Qassert(fr.is_read_through);
      Qassert(fr.fn_list.size() == ret.size());
      for (Long j = 0; j < (Long)ret.size(); ++j) {
        Qassert(ret[j] == fr.fn_list[j]);
      }
    }
  }
  return ret;
}

std::vector<std::string> list_fields(const ShuffledFieldsWriter& sfw, bool is_skipping_check)
// interface function
{
  TIMER_VERBOSE("list_fields(sfw)");
  std::vector<std::string> ret;
  if (0 == get_id_node()) {
    Qassert(sfw.fws.size() > 0);
    const FieldsWriter& fw = sfw.fws[0];
    ret = fw.fn_list;
  }
  bcast(ret);
  if (not is_skipping_check) {
    for (Int i = 0; i < (int)sfw.fws.size(); ++i) {
      const FieldsWriter& fw = sfw.fws[i];
      Qassert(fw.fn_list.size() == ret.size());
      for (Long j = 0; j < (Long)ret.size(); ++j) {
        Qassert(ret[j] == fw.fn_list[j]);
      }
    }
  }
  return ret;
}

// ------------------------

Int truncate_fields_sync_node(const std::string& path,
                              const std::vector<std::string>& fns_keep,
                              const Coordinate& new_size_node)
{
  TIMER_VERBOSE("truncate_fields_sync_node");
  ShuffledFieldsReader sfr;
  sfr.init(path, new_size_node);
  const std::vector<std::string> fns = list_fields(sfr, true);
  if (fns.size() < fns_keep.size()) {
    qwarn(fname + ssprintf(": fns.size()=%ld fns_keep.size()=%ld", fns.size(),
                           fns_keep.size()));
    return 1;
  }
  for (Long i = 0; i < (Long)fns_keep.size(); ++i) {
    if (fns[i] != fns_keep[i]) {
      qwarn(fname + ssprintf(": fns[i]='%s' fns_keep[i]='%s'", fns[i].c_str(),
                             fns_keep[i].c_str()));
      return 2;
    }
  }
  std::vector<Long> final_offsets(sfr.frs.size(), 0);
  if (fns_keep.size() >= 1) {
    const std::string& fn_last = fns_keep.back();
    const bool is_fn_last_valid =
        check_file_sync_node(sfr, fn_last, true, final_offsets);
    if (not is_fn_last_valid) {
      qwarn(fname + ssprintf(": fn_last='%s' check failed", fn_last.c_str()));
      return 2;
    }
  }
  for (Int i = 0; i < (int)sfr.frs.size(); ++i) {
    FieldsReader& fr = sfr.frs[i];
    const std::string path_file = get_file_path(fr);
    const Long file_size = fr.file_size;
    fr.close();
    const Long final_offset = final_offsets[i];
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
      const Int b = qtruncate(path_file, final_offset);
      Qassert(b == 0);
      fr.fn_list = fns_keep;
      fr.max_offset = final_offset;
    }
  }
  fields_build_index(sfr);
  return 0;
}

void properly_truncate_fields_sync_node(
    std::vector<std::string>& fn_list,
    std::vector<std::vector<FieldsSegmentInfo>>& offsets_list, const std::string& path,
    const bool is_check_all, const bool is_only_check,
    const Coordinate& new_size_node)
// interface function
// offsets_list.size() == fn_list.size()
// offsets_list[0].size() == sfr.frs.size()
{
  TIMER_VERBOSE("properly_truncate_fields_sync_node");
  fn_list.clear();
  offsets_list.clear();
  if (not does_file_exist_qar_sync_node(path + "/geon-info.txt")) {
    displayln_info(0, fname + ssprintf(": '%s' does not exist.", path.c_str()));
    return;
  }
  ShuffledFieldsReader sfr;
  sfr.init(path, new_size_node);
  fn_list = list_fields(sfr, true);
  std::vector<Long> last_final_offsets(sfr.frs.size(), 0);
  Long last_idx = -1;
  if (is_check_all) {
    for (Long i = 0; i < (Long)fn_list.size(); ++i) {
      const std::string& fn = fn_list[i];
      std::vector<Long> final_offsets;
      const bool b = check_file_sync_node(sfr, fn, true, final_offsets);
      if (b) {
        last_final_offsets = final_offsets;
        last_idx = i;
      } else {
        break;
      }
    }
  } else {
    for (Long i = (Long)fn_list.size() - 1; i >= 0; i -= 1) {
      const std::string& fn = fn_list[i];
      std::vector<Long> final_offsets;
      const bool b = check_file_sync_node(sfr, fn, false, final_offsets);
      if (b) {
        last_final_offsets = final_offsets;
        last_idx = i;
        break;
      }
    }
  }
  fn_list.resize(last_idx + 1);
  offsets_list.resize(fn_list.size());
  for (Int j = 0; j < (int)sfr.frs.size(); ++j) {
    FieldsReader& fr = sfr.frs[j];
    fr.fn_list.resize(last_idx + 1);
    Qassert(fr.fn_list == fn_list);
    for (Long i = 0; i < (Long)fn_list.size(); ++i) {
      const std::string& fn = fn_list[i];
      const FieldsSegmentInfo offset = fr.offsets_map[fn];
      offsets_list[i].push_back(offset);
    }
  }
  if (not is_only_check) {
    fields_build_index(sfr);
  }
  for (Int i = 0; i < (int)sfr.frs.size(); ++i) {
    FieldsReader& fr = sfr.frs[i];
    const std::string path_file = get_file_path(fr);
    const Long file_size = fr.file_size;
    fr.close();
    const Long final_offset = last_final_offsets[i];
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
        const Int b = qtruncate(path_file, final_offset);
        Qassert(b == 0);
      }
    }
  }
  sfr.close();
  errno = 0;
  for (Long i = 0; i < (Long)fn_list.size(); ++i) {
    const std::string& fn = fn_list[i];
    displayln_info(0, fname + ssprintf(": i=%5ld fn='%s'", i, fn.c_str()));
  }
  displayln_info(0, fname + ssprintf(": fn_list.size()=%5ld '%s'",
                                     fn_list.size(), path.c_str()));
  SYNC_NODE();
}

std::vector<std::string> properly_truncate_fields_sync_node(
    const std::string& path, const bool is_check_all, const bool is_only_check,
    const Coordinate& new_size_node)
// interface function
// return available fns
{
  std::vector<std::string> fn_list;
  std::vector<std::vector<FieldsSegmentInfo>> offsets_list;
  properly_truncate_fields_sync_node(fn_list, offsets_list, path, is_check_all,
                                     is_only_check, new_size_node);
  return fn_list;
}

bool has_duplicates(const ShuffledFieldsReader& sfr)
{
  TIMER("has_duplicates(sfr)");
  const std::vector<std::string> fn_list = list_fields(sfr);
  Long count = 0;
  for (Int j = 0; j < (int)sfr.frs.size(); ++j) {
    const FieldsReader& fr = sfr.frs[j];
    Qassert(fr.is_read_through);
    Qassert(fr.fn_list.size() == fn_list.size());
    Qassert(fr.offsets_map.size() <= fn_list.size());
    count += fn_list.size() - fr.offsets_map.size();
  }
  glb_sum(count);
  return count > 0;
}

// ------------------------

std::string show_field_index(const std::string& fn, const std::vector<FieldsSegmentInfo>& all_offsets)
// all_offsets.size() == total number of part-files.
{
  QFile qfile = qfopen(QFileType::String, "show_field_index", QFileMode::Write);
  qfile.write_data(ssprintf("%ld\n", (long)fn.size()));
  qfile.write_data(fn);
  qfile.write_data("\n");
  const bool is_sparse_field = all_offsets[0].is_sparse_field;
  for (Long i = 0; i < (Long)all_offsets.size(); ++i) {
    Qassert(is_sparse_field == all_offsets[i].is_sparse_field);
  }
  if (is_sparse_field) {
    qfile.write_data("sparse");
  } else {
    qfile.write_data("dense");
  }
  qfile.write_data("\n");
  for (Long i = 0; i < (Long)all_offsets.size(); ++i) {
    qfile.write_data(ssprintf("%ld", (long)all_offsets[i].offset_start));
    if (i != (Long)all_offsets.size() - 1) {
      qfile.write_data(" ");
    }
  }
  qfile.write_data("\n");
  for (Long i = 0; i < (Long)all_offsets.size(); ++i) {
    qfile.write_data(ssprintf("%ld", (long)all_offsets[i].offset_stop));
    if (i != (Long)all_offsets.size() - 1) {
      qfile.write_data(" ");
    }
  }
  qfile.write_data("\n");
  const std::string ret = qfile.content();
  qfile.close();
  return ret;
}

void parse_field_index(std::string& fn, std::vector<FieldsSegmentInfo>& all_offsets, const std::string& field_index_content)
{
  fn.clear();
  all_offsets.clear();
  Long cur = 0;
  std::string fn_len_str;
  bool b;
  b = parse_line(fn_len_str, cur, field_index_content);
  Qassert(b);
  const Long fn_len = read_long(fn_len_str);
  b = parse_len(fn, cur, field_index_content, fn_len);
  Qassert(b);
  b = parse_literal(cur, field_index_content, '\n');
  Qassert(b);
  std::string is_sparse_field_str;
  b = parse_line(is_sparse_field_str, cur, field_index_content);
  Qassert(b);
  bool is_sparse_field;
  if (is_sparse_field_str == "sparse\n") {
    is_sparse_field = true;
  } else if (is_sparse_field_str == "dense\n") {
    is_sparse_field = false;
  } else {
    qerr(ssprintf("parse_field_index: is_sparse_field_str='%s'.",
                  is_sparse_field_str.c_str()));
  }
  std::string offsets_start_str;
  b = parse_line(offsets_start_str, cur, field_index_content);
  Qassert(b);
  std::vector<Long> all_offsets_start = read_longs(offsets_start_str);
  Qassert(all_offsets_start.size() > 0);
  std::string offsets_stop_str;
  b = parse_line(offsets_stop_str, cur, field_index_content);
  Qassert(b);
  std::vector<Long> all_offsets_stop = read_longs(offsets_stop_str);
  Qassert(all_offsets_stop.size() > 0);
  Qassert(all_offsets_stop.size() == all_offsets_start.size());
  b = parse_end(cur, field_index_content);
  Qassert(b);
  all_offsets.resize(all_offsets_start.size());
  for (Long i = 0; i < (Long)all_offsets_start.size(); ++i) {
    all_offsets[i] = FieldsSegmentInfo(all_offsets_start[i],
                                       all_offsets_stop[i], is_sparse_field);
  }
}

std::vector<FieldsSegmentInfo> collect_fields_offsets(const ShuffledFieldsWriter& sfw, const std::string& fn)
{
  TIMER("collect_fields_offsets(sfw,fn)");
  Long num_files = sfw.fws.size();
  glb_sum(num_files);
  std::vector<FieldsSegmentInfo> all_offsets(num_files);
  Vector<char> all_offsets_view((char*)all_offsets.data(),
                                all_offsets.size() * sizeof(FieldsSegmentInfo));
  set_zero(all_offsets_view);
  for (Int i = 0; i < (int)sfw.fws.size(); ++i) {
    const FieldsWriter& fw = sfw.fws[i];
    Qassert(fw.geon.num_node == num_files);
    Qassert(fw.geon.id_node >= 0);
    Qassert(fw.geon.id_node < num_files);
    if (has(fw.offsets_map, fn)) {
      all_offsets[fw.geon.id_node] = fw.offsets_map.at(fn);
    }
  }
  glb_sum(all_offsets_view);
  return all_offsets;
}

std::vector<FieldsSegmentInfo> collect_fields_offsets(const ShuffledFieldsReader& sfr, const std::string& fn)
{
  TIMER("collect_fields_offsets(sfr,fn)");
  Long num_files = sfr.frs.size();
  glb_sum(num_files);
  std::vector<FieldsSegmentInfo> all_offsets(num_files);
  Vector<char> all_offsets_view((char*)all_offsets.data(),
                                all_offsets.size() * sizeof(FieldsSegmentInfo));
  set_zero(all_offsets_view);
  for (Int i = 0; i < (int)sfr.frs.size(); ++i) {
    const FieldsReader& fr = sfr.frs[i];
    Qassert(fr.geon.num_node == num_files);
    Qassert(fr.geon.id_node >= 0);
    Qassert(fr.geon.id_node < num_files);
    if (has(fr.offsets_map, fn)) {
      all_offsets[fr.geon.id_node] = fr.offsets_map.at(fn);
    }
  }
  glb_sum(all_offsets_view);
  return all_offsets;
}

bool populate_fields_offsets(ShuffledFieldsReader& sfr, const std::string& fn,
                             const std::vector<FieldsSegmentInfo>& all_offsets)
// should populate fields offsets sequentially in correct order
{
  TIMER("populate_fields_offsets(sfr,fn)");
  const Long num_files = all_offsets.size();
  Long count = 0;
  for (Int i = 0; i < (int)sfr.frs.size(); ++i) {
    FieldsReader& fr = sfr.frs[i];
    Qassert(fr.geon.num_node == num_files);
    Qassert(fr.geon.id_node >= 0);
    Qassert(fr.geon.id_node < num_files);
    if (has(fr.offsets_map, fn)) {
      qwarn(fname + ssprintf(": duplicate fn='%s' id_node='%d'", fn.c_str(),
                             fr.geon.id_node));
    }
    const Long offset_start = all_offsets[fr.geon.id_node].offset_start;
    const Long offset_stop = all_offsets[fr.geon.id_node].offset_stop;
    Qassert(offset_stop > offset_start);
    if (fr.file_size < offset_stop) {
      qwarn(fname + ssprintf(": file_size not large enough. fn='%s' "
                             "id_node='%d' file_size='%ld' offset_stop='%ld'",
                             fn.c_str(), fr.geon.id_node, (long)fr.file_size,
                             (long)offset_stop));
    } else if (fr.max_offset != offset_start) {
      qwarn(fname + ssprintf(": fr.max_offset does not match. fn='%s' "
                             "id_node='%d' fr.max_offset='%ld' offset='%ld'",
                             fn.c_str(), fr.geon.id_node, (long)fr.max_offset,
                             (long)offset_start));
    } else {
      count += 1;
    }
  }
  glb_sum(count);
  if (count != num_files) {
    if (get_id_node() == 0) {
      qwarn(fname + ssprintf(": '%s' ; count=%ld ; num_files=%ld.", fn.c_str(),
                             (long)count, (long)num_files));
    }
    return false;
  }
  for (Int i = 0; i < (int)sfr.frs.size(); ++i) {
    FieldsReader& fr = sfr.frs[i];
    fr.fn_list.push_back(fn);
    fr.offsets_map[fn] = all_offsets[fr.geon.id_node];
    fr.max_offset = all_offsets[fr.geon.id_node].offset_stop;
  }
  return true;
}

void load_all_fields_index(std::vector<std::string>& fn_list,
                           std::vector<std::vector<FieldsSegmentInfo>>& all_offsets_list,
                           QarFile& qar_index)
{
  TIMER("load_all_fields_index(fn_list,all_offsets_list,qar_index)");
  fn_list.clear();
  all_offsets_list.clear();
  Long num_files = 0;
  if (get_id_node() == 0) {
    const std::vector<std::string> qar_idx_list = list(qar_index);
    for (Long i = 0; i < (Long)qar_idx_list.size(); ++i) {
      const std::string& qar_idx = qar_idx_list[i];
      const std::string field_index_content = read_data(qar_index, qar_idx);
      std::string fn;
      std::vector<FieldsSegmentInfo> all_offsets;
      parse_field_index(fn, all_offsets, field_index_content);
      Qassert((Long)all_offsets.size() > 0);
      if (num_files == 0) {
        num_files = all_offsets.size();
      } else {
        Qassert(num_files == (Long)all_offsets.size());
      }
      fn_list.push_back(fn);
      all_offsets_list.push_back(all_offsets);
    }
    Qassert(all_offsets_list.size() == fn_list.size());
  }
  bcast(fn_list);
  bcast(num_files);
  all_offsets_list.resize(fn_list.size());
  for (Long i = 0; i < (Long)all_offsets_list.size(); ++i) {
    std::vector<FieldsSegmentInfo>& all_offsets = all_offsets_list[i];
    all_offsets.resize(num_files);
    Vector<char> v((char*)all_offsets.data(),
                   all_offsets.size() * sizeof(FieldsSegmentInfo));
    bcast(v);
  }
}

void save_fields_index(QarFile& qar_index, const Long idx,
                       const std::string& fn, const std::vector<FieldsSegmentInfo>& all_offsets)
{
  TIMER("save_fields_index(qar_index,idx,fn,all_offsets)");
  if (get_id_node() == 0) {
    const std::string fn_index = show_field_index(fn, all_offsets);
    write_from_data(qar_index, show(idx), "", fn_index);
    qar_index.flush();
  }
}

void save_fields_index(ShuffledFieldsWriter& sfw, const std::string& fn)
{
  TIMER("save_fields_index(sfw,fn)");
  const std::vector<FieldsSegmentInfo> all_offsets = collect_fields_offsets(sfw, fn);
  save_fields_index(sfw.qar_index, sfw.qar_index_idx, fn, all_offsets);
  sfw.qar_index_idx += 1;
}

void fields_build_index(const ShuffledFieldsReader& sfr)
{
  TIMER_VERBOSE("fields_build_index(sfr)");
  const std::vector<std::string> fn_list = list_fields(sfr);
  QarFile qar_index;
  if (get_id_node() == 0) {
    qar_index.init(sfr.path + "/index.qar", QFileMode::Write);
  }
  for (Long i = 0; i < (Long)fn_list.size(); ++i) {
    const std::string& fn = fn_list[i];
    const std::vector<FieldsSegmentInfo> all_offsets = collect_fields_offsets(sfr, fn);
    save_fields_index(qar_index, i, fn, all_offsets);
  }
  qar_index.close();
}

void fields_build_index(const std::string& path)
{
  TIMER_VERBOSE("fields_build_index(path)");
  ShuffledFieldsReader sfr(path);
  fields_build_index(sfr);
  sfr.close();
}

// ------------------------

ShuffledFieldsReader& get_shuffled_fields_reader(
    const std::string& path, const Coordinate& new_size_node)
{
  TIMER("get_shuffled_fields_reader");
  ShuffledFieldsReader& sfr = get_shuffled_fields_reader_cache()[path];
  if (get_shuffled_fields_writer_cache().has(path)) {
    // warn if there is a writer in cache
    if (get_id_node() == 0) {
      qwarn(fname + ssprintf(": path='%s' is in writer cache.", path.c_str()));
    }
  }
  if (sfr.path == "") {
    sfr.init(path, new_size_node);
  }
  return sfr;
}

ShuffledFieldsWriter& get_shuffled_fields_writer(
    const std::string& path, const Coordinate& new_size_node,
    const bool is_append, const bool is_removing_old)
{
  TIMER("get_shuffled_fields_writer");
  ShuffledFieldsWriter& sfw = get_shuffled_fields_writer_cache()[path];
  if (get_shuffled_fields_reader_cache().has(path)) {
    if (get_id_node() == 0) {
      qwarn(
          fname +
          ssprintf(": path='%s' is in reader cache. Remove this cached reader.",
                   path.c_str()));
    }
    close_shuffled_fields_reader(path);
  }
  if (sfw.path == "") {
    sfw.init(path, new_size_node, is_append, is_removing_old);
  }
  return sfw;
}

std::vector<std::string> list_fields(const std::string& path)
{
  ShuffledFieldsReader& sfr = get_shuffled_fields_reader(path);
  return list_fields(sfr);
}

bool does_file_exist_sync_node(const std::string& path, const std::string& fn)
{
  TIMER_VERBOSE("does_file_exist_sync_node(path,fn)");
  ShuffledFieldsReader& sfr = get_shuffled_fields_reader(path);
  return does_file_exist_sync_node(sfr, fn);
}

}  // namespace qlat
