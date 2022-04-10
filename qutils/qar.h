#pragma once

#include <stdint.h>
#include <qutils/qutils-vec.h>
#include <qutils/qutils-io.h>

#include <cassert>
#include <string>
#include <vector>
#include <utility>

namespace qlat
{  //

struct QFile {
  // Interface to FILE* which allow a view of a portion of the file specified by
  // offset_start and offset_end.
  // The view can be nested.
  //
  std::string path;
  std::string mode;  // can be "r", "a", "w"
  FILE* fp;
  //
  QFile* parent;  // If parent == NULL, then this QFile own the fp pointer and
                  // will be responsible for close it.
  long number_of_child;  // Can only close the FILE when number_of_child == 0.
  //
  bool is_eof;  // the eof state of QFile.
                // NOTE: may not match with eof state of fp.
  long pos;  // position of the QFile. (correspond to position of fp should be
             // pos + offset_start).
             // NOTE: actual fp position may be adjust elsewhere and does not
             // match this pos.
  //
  long offset_start;  // start offset of fp for QFile
  long offset_end;    // end offset of fp for QFile (-1 if not limit, useful
                      // when writing)
  //
  QFile()
  {
    fp = NULL;
    parent = NULL;
    number_of_child = 0;
    init();
  }
  QFile(const std::string& path_, const std::string& mode_)
  {
    fp = NULL;
    parent = NULL;
    number_of_child = 0;
    init(path_, mode_);
  }
  QFile(QFile& qfile, const long q_offset_start, const long q_offset_end)
  {
    fp = NULL;
    parent = NULL;
    number_of_child = 0;
    init(qfile, q_offset_start, q_offset_end);
  }
  //
  QFile(const QFile&) = delete;
  //
  QFile(QFile&& qfile) noexcept
  {
    fp = NULL;
    parent = NULL;
    number_of_child = 0;
    init();
    swap(qfile);
  }
  //
  ~QFile() { close(); }
  //
  void init()
  {
    close();
    path = "";
    mode = "";
    is_eof = false;
    pos = 0;
    offset_start = 0;
    offset_end = -1;
  }
  void init(const std::string& path_, const std::string& mode_)
  {
    close();
    path = path_;
    mode = mode_;
    if (verbose_level() > 0) {
      displayln(
          ssprintf("QFile: open '%s' with '%s'.", path.c_str(), mode.c_str()));
    }
    fp = qopen(path, mode);
    if (fp == NULL) {
      qwarn(
          ssprintf("QFile: open '%s' with '%s'.", path.c_str(), mode.c_str()));
    }
    is_eof = false;
    pos = ftell(fp);
    offset_start = 0;
    offset_end = -1;
  }
  void init(QFile& qfile, const long q_offset_start, const long q_offset_end)
  // Become a child of qfile.
  // NOTE: q_offset_start and q_offset_end are relative offset for qfile not the
  // absolute offset for qfile.fp .
  // q_offset_end == -1 means no additional limit
  // NOTE: Initial position set to be 0. Does not perform fseek to appropriate
  // position.
  {
    close();
    if (qfile.null()) {
      return;
    }
    qfile.number_of_child += 1;
    path = qfile.path;
    mode = qfile.mode;
    fp = qfile.fp;
    parent = &qfile;
    qassert(q_offset_start >= 0);
    is_eof = false;
    pos = 0;
    offset_start = qfile.offset_start + q_offset_start;
    if (q_offset_end == -1) {
      offset_end = qfile.offset_end;
    } else {
      qassert(q_offset_end >= q_offset_start);
      offset_end = qfile.offset_start + q_offset_end;
      if (qfile.offset_end != -1) {
        qassert(offset_end <= qfile.offset_end);
      }
    }
    if (mode == "r" and offset_end != -1) {
      const int code = fseek(fp, offset_end, SEEK_SET);
      if (code != 0) {
        qwarn(ssprintf("QFile: '%s' with '%s' offset=%ld,%ld failed.",
                       path.c_str(), mode.c_str(), offset_start, offset_end));
        close();
      }
    }
  }
  //
  void close()
  {
    // to close the file, it cannot have any child
    qassert(number_of_child == 0);
    if (NULL == parent) {
      if (fp != NULL) {
        if (verbose_level() > 0) {
          displayln(ssprintf("QFile: close '%s' with '%s'.", path.c_str(),
                             mode.c_str()));
        }
        qclose(fp);
      }
    } else {
      fp = NULL;
      (*parent).number_of_child -= 1;
      parent = NULL;
    }
    qassert(fp == NULL);
    qassert(parent == NULL);
  }
  //
  void swap(QFile& qfile)
  {
    // cannot swap if has child
    qassert(number_of_child == 0);
    qassert(qfile.number_of_child == 0);
    std::swap(path, qfile.path);
    std::swap(mode, qfile.mode);
    std::swap(fp, qfile.fp);
    std::swap(parent, qfile.parent);
    std::swap(number_of_child, qfile.number_of_child);
    std::swap(is_eof, qfile.is_eof);
    std::swap(pos, qfile.pos);
    std::swap(offset_start, qfile.offset_start);
    std::swap(offset_end, qfile.offset_end);
  }
  //
  bool null() { return fp == NULL; }
};

inline void qswap(QFile& qfile1, QFile& qfile2)
// interface function
{
  qfile1.swap(qfile2);
}

inline bool qfeof(const QFile& qfile)
// interface function
{
  return qfile.is_eof;
}

inline long qftell(const QFile& qfile)
// interface function
{
  return qfile.pos;
}

inline int qfseek(QFile& qfile, const long q_offset, const int whence)
// interface function
// Always call fseek and adjust qfile.is_eof and qfile.pos
// qfile.pos will be set to the actual QFile position after qfseek.
// return 0 if successful
{
  qfile.is_eof = false;
  int ret = 0;
  if (SEEK_SET == whence) {
    const long offset = qfile.offset_start + q_offset;
    ret = fseek(qfile.fp, offset, SEEK_SET);
  } else if (SEEK_CUR == whence) {
    ret = fseek(qfile.fp, qfile.offset_start + qfile.pos + q_offset, SEEK_SET);
  } else if (SEEK_END == whence) {
    if (qfile.offset_end == -1) {
      ret = fseek(qfile.fp, q_offset, SEEK_END);
    } else {
      const long offset = qfile.offset_end + q_offset;
      ret = fseek(qfile.fp, offset, SEEK_SET);
    }
  } else {
    qassert(false);
  }
  qfile.pos = ftell(qfile.fp) - qfile.offset_start;
  qassert(qfile.pos >= 0);
  if (qfile.offset_end != -1) {
    qassert(qfile.offset_start + qfile.pos <= qfile.offset_end);
  }
  return ret;
}

inline long qfread(void* ptr, const long size, const long nmemb, QFile& qfile)
// interface function
// Only read portion of data if not enough content in qfile.
{
  if (0 == size or 0 == nmemb) {
    return 0;
  }
  qassert(size > 0);
  qassert(nmemb > 0);
  const int code = qfseek(qfile, qfile.pos, SEEK_SET);
  qassert(code == 0);
  long actual_nmemb = 0;
  if (qfile.offset_end != -1) {
    const long remaining_size =
        qfile.offset_end - qfile.offset_start - qfile.pos;
    qassert(remaining_size >= 0);
    const long target_nmemb = std::min(remaining_size / size, nmemb);
    actual_nmemb = std::fread(ptr, size, target_nmemb, qfile.fp);
    qassert(actual_nmemb == target_nmemb);
    qfile.pos += target_nmemb * size;
    qassert(qfile.pos == ftell(qfile.fp) - qfile.offset_start);
    if (target_nmemb < nmemb) {
      qfile.is_eof = true;
    } else {
      qassert(target_nmemb == nmemb);
      qfile.is_eof = false;
    }
  } else {
    actual_nmemb = std::fread(ptr, size, nmemb, qfile.fp);
    qfile.pos = ftell(qfile.fp) - qfile.offset_start;
    qfile.is_eof = feof(qfile.fp) != 0;
  }
  return actual_nmemb;
}

inline long qfwrite(const void* ptr, const long size, const long nmemb,
                    QFile& qfile)
// interface function
// Crash if no enough space
{
  if (0 == size or 0 == nmemb) {
    return 0;
  }
  qassert(size > 0);
  qassert(nmemb > 0);
  const int code = qfseek(qfile, qfile.pos, SEEK_SET);
  qassert(code == 0);
  if (qfile.offset_end != -1) {
    const long remaining_size =
        qfile.offset_end - qfile.offset_start - qfile.pos;
    qassert(remaining_size >= size * nmemb);
  }
  const long actual_nmemb = std::fwrite(ptr, size, nmemb, qfile.fp);
  qassert(actual_nmemb == nmemb);
  qfile.pos = ftell(qfile.fp) - qfile.offset_start;
  qassert(qfile.pos >= 0);
  if (qfile.offset_end != -1) {
    qassert(qfile.offset_start + qfile.pos <= qfile.offset_end);
  }
  return actual_nmemb;
}

inline std::string qgetline(QFile& qfile)
// interface function
// read an entire line including the final '\n' char.
{
  const int code = qfseek(qfile, qfile.pos, SEEK_SET);
  qassert(code == 0);
  char* lineptr = NULL;
  size_t n = 0;
  const long size = getline(&lineptr, &n, qfile.fp);
  qfile.is_eof = feof(qfile.fp) != 0;
  if (size > 0) {
    std::string ret;
    const long pos = ftell(qfile.fp) - qfile.offset_start;
    qassert(pos >= 0);
    if (qfile.offset_end != -1 and
        qfile.offset_start + pos > qfile.offset_end) {
      qfseek(qfile, 0, SEEK_END);
      qfile.is_eof = true;
      const long size_truncate =
          size - (qfile.offset_start + pos - qfile.offset_end);
      qassert(size_truncate >= 0);
      ret = std::string(lineptr, size_truncate);
    } else {
      qfile.pos = pos;
      ret = std::string(lineptr, size);
    }
    std::free(lineptr);
    return ret;
  } else {
    std::free(lineptr);
    return std::string();
  }
}

inline std::vector<std::string> qgetlines(QFile& qfile)
// interface function
{
  std::vector<std::string> ret;
  while (not qfeof(qfile)) {
    ret.push_back(qgetline(qfile));
  }
  return ret;
}

template <class M>
long qwrite_data(const Vector<M>& v, QFile& qfile)
// interface function
{
  TIMER_FLOPS("qwrite_data(v,qfile)");
  const long data_size = sizeof(M) * qfwrite((void*)v.p, sizeof(M), v.n, qfile);
  timer.flops += data_size;
  return data_size;
}

inline long qwrite_data(const std::string& line, QFile& qfile)
// interface function
{
  return qwrite_data(get_data(line), qfile);
}

inline long qwrite_data(const std::vector<std::string>& lines, QFile& qfile)
// interface function
{
  long total_bytes = 0;
  for (long i = 0; i < (long)lines.size(); ++i) {
    total_bytes += qwrite_data(lines[i], qfile);
  }
  return total_bytes;
}

template <class M>
long qread_data(const Vector<M>& v, QFile& qfile)
// interface function
{
  TIMER_FLOPS("qread_data(v,qfile)");
  const long data_size = sizeof(M) * qfread((void*)v.p, sizeof(M), v.n, qfile);
  timer.flops += data_size;
  return data_size;
}

template <class M>
long qread_data_all(std::vector<M>& v, QFile& qfile)
// interface function
// Read all the remaining data.
// (Remaining size must be multiple of sizeof(M) otherwise will fail.)
// return total bytes read.
{
  TIMER_FLOPS("qread_data_all(v,qfile)");
  const long pos_initial = qftell(qfile);
  qfseek(qfile, 0, SEEK_END);
  const long pos_final = qftell(qfile);
  const long data_size = pos_final - pos_initial;
  const long n = data_size / sizeof(M);
  qassert(data_size == sizeof(M) * n);
  qfseek(qfile, pos_initial, SEEK_SET);
  v.resize(n);
  const long data_size_read = qread_data(get_data(v), qfile);
  qassert(data_size_read == data_size);
  timer.flops += data_size;
  return data_size;
}

inline long& write_from_qfile_chunk_size()
// qlat parameter
// size in bytes
{
  static long size = 512 * 1024;
  return size;
}

inline long write_from_qfile(QFile& qfile_out, QFile& qfile_in)
{
  TIMER_FLOPS("write_from_qfile(qfile_out,qfile_in)");
  const long chunk_size = write_from_qfile_chunk_size();
  std::vector<char> buf(chunk_size);
  long total_bytes = 0;
  while (not qfeof(qfile_in)) {
    const long size = qread_data(get_data(buf), qfile_in);
    qassert(size <= chunk_size);
    const long size_out = qwrite_data(get_data(get_data(buf), size), qfile_out);
    qassert(size_out == size);
    total_bytes += size;
  }
  timer.flops += total_bytes;
  return total_bytes;
}

const std::string qar_header = "#!/usr/bin/env qar-glimpse\n\n";

struct QarSegmentInfo {
  long offset;
  long offset_fn;
  long offset_info;
  long offset_data;
  long offset_end;
  long fn_len;
  long info_len;
  long data_len;
  //
  QarSegmentInfo() { init(); }
  //
  void init() { set_zero(get_data_one_elem(*this)); }
};

struct QarFile {
  QFile qfile;
  //
  bool is_read_through;
  std::vector<std::string> fn_list;
  std::map<std::string, QarSegmentInfo> qsinfo_map;
  long max_offset;  // when read, maximum offset reached so far
  //
  long current_write_segment_offset;  // when write, the offset of the beginning
                                      // of the current working segment.
  long current_write_segment_offset_data;  // when write, the offset of the
                                           // beginning of the data of the
                                           // current working segment.
  long current_write_segment_offset_header_len;
  long current_write_segment_offset_fn_len;
  long current_write_segment_offset_info_len;
  long current_write_segment_offset_data_len;
  //
  QarFile() { init(); }
  QarFile(const std::string& path, const std::string& mode)
  {
    init(path, mode);
  }
  //
  void init()
  {
    qfile.init();
    is_read_through = false;
    fn_list.clear();
    qsinfo_map.clear();
    max_offset = 0;
    current_write_segment_offset = -1;
    current_write_segment_offset_data = -1;
    current_write_segment_offset_header_len = -1;
    current_write_segment_offset_fn_len = -1;
    current_write_segment_offset_info_len = -1;
    current_write_segment_offset_data_len = -1;
  }
  void init(const std::string& path, const std::string& mode)
  {
    init();
    qfile.init(path, mode);
    if (qfile.null()) {
      return;
    }
    if (mode == "w") {
      qfwrite(qar_header.data(), qar_header.size(), 1, qfile);
    } else if (mode == "r") {
      std::vector<char> check_line(qar_header.size(), 0);
      const long qfread_check_len =
          qfread(check_line.data(), qar_header.size(), 1, qfile);
      if (not(qfread_check_len == 1 and
              std::string(check_line.data(), check_line.size()) ==
                  qar_header)) {
        qfile.close();
      };
      max_offset = qftell(qfile);
      QarSegmentInfo qsinfo;
      qsinfo.offset_end = qftell(qfile);
      qsinfo_map[""] = qsinfo;
    }
  }
  //
  void close() { init(); }
  //
  bool null() { return qfile.null(); }
};

inline bool read_qar_segment_info(QarFile& qar, QarSegmentInfo& qsinfo)
// Initial pos: beginning of the segment, just before FILE-HEADER.
// Final pos: at the end of the segment, at the beginning of the next segment.
// Return true if read successfully (also qfseek to the beginning of the next
// segment).
{
  qassert(qar.qfile.mode == "r");
  set_zero(get_data_one_elem(qsinfo));
  if (qar.qfile.null()) {
    qwarn(ssprintf("read_tag: fn='%s' pos=%ld.", qar.qfile.path.c_str(),
                   qftell(qar.qfile)));
    return false;
  }
  qsinfo.offset = qftell(qar.qfile);
  const std::string header_prefix = "QAR-FILE ";
  const std::string header = qgetline(qar.qfile);
  if (header.size() == 0) {
    qar.is_read_through = true;
    return false;
  }
  if (header.size() <= header_prefix.size()) {
    qwarn(ssprintf("read_tag: fn='%s' pos=%ld.", qar.qfile.path.c_str(),
                   qftell(qar.qfile)));
    qar.is_read_through = true;
    return false;
  }
  if (header.substr(0, header_prefix.size()) != header_prefix) {
    qwarn(ssprintf("read_tag: fn='%s' pos=%ld.", qar.qfile.path.c_str(),
                   qftell(qar.qfile)));
    qar.is_read_through = true;
    return false;
  }
  const std::vector<long> len_vec = read_longs(header.substr(header_prefix.size()));
  if (len_vec.size() != 3) {
    qwarn(ssprintf("read_tag: fn='%s' pos=%ld.", qar.qfile.path.c_str(),
                   qftell(qar.qfile)));
    qar.is_read_through = true;
    return false;
  }
  qsinfo.fn_len = len_vec[0];
  qsinfo.info_len = len_vec[1];
  qsinfo.data_len = len_vec[2];
  qsinfo.offset_fn = qftell(qar.qfile);
  qsinfo.offset_info = qsinfo.offset_fn + qsinfo.fn_len + 1;
  qsinfo.offset_data = qsinfo.offset_info + qsinfo.info_len + 1;
  qsinfo.offset_end = qsinfo.offset_data + qsinfo.data_len + 2;
  const int code = qfseek(qar.qfile, qsinfo.offset_end, SEEK_SET);
  if (code != 0) {
    qwarn(ssprintf("read_tag: fn='%s' pos=%ld offset_end=%ld.",
                   qar.qfile.path.c_str(), qftell(qar.qfile),
                   qsinfo.offset_end));
    qar.is_read_through = true;
    return false;
  }
  if (qar.max_offset < qsinfo.offset_end) {
    qar.max_offset = qsinfo.offset_end;
  }
  return true;
}

inline void read_fn(QarFile& qar, std::string& fn, const QarSegmentInfo& qsinfo)
// interface function
{
  qassert(qar.qfile.mode == "r");
  qassert(not qar.qfile.null());
  std::vector<char> data(qsinfo.fn_len);
  const int code = qfseek(qar.qfile, qsinfo.offset_fn, SEEK_SET);
  qassert(code == 0);
  if (1 != qfread(data.data(), qsinfo.fn_len, 1, qar.qfile)) {
    qassert(false);
  }
  fn = std::string(data.data(), qsinfo.fn_len);
}

inline void read_info(QarFile& qar, std::string& info,
                      const QarSegmentInfo& qsinfo)
// interface function
{
  qassert(qar.qfile.mode == "r");
  qassert(not qar.qfile.null());
  std::vector<char> data(qsinfo.info_len);
  const int code = qfseek(qar.qfile, qsinfo.offset_info, SEEK_SET);
  qassert(code == 0);
  if (1 != qfread(data.data(), qsinfo.info_len, 1, qar.qfile)) {
    qassert(false);
  }
  info = std::string(data.data(), qsinfo.info_len);
}

inline void get_qfile_of_data(QarFile& qar, QFile& qfile,
                              const QarSegmentInfo& qsinfo)
// interface function
// set qfile to be a qfile containing the data specified by qsinfo.
// qfile initial pos is zero
{
  qassert(qar.qfile.mode == "r");
  qassert(not qar.qfile.null());
  qfile.init(qar.qfile, qsinfo.offset_data,
             qsinfo.offset_data + qsinfo.data_len);
  qassert(not qfile.null());
}

inline bool read_next(QarFile& qar, std::string& fn, QFile& qfile)
// interface function
// Initial pos of qar should be at the beginning of a segment.
{
  qassert(qar.qfile.mode == "r");
  fn = std::string();
  qfile.init();
  QarSegmentInfo qsinfo;
  if (not read_qar_segment_info(qar, qsinfo)) {
    return false;
  }
  read_fn(qar, fn, qsinfo);
  if (not has(qar.qsinfo_map, fn)) {
    qar.fn_list.push_back(fn);
    qar.qsinfo_map[fn] = qsinfo;
  } else {
    qassert(qar.qsinfo_map[fn].offset == qsinfo.offset);
  }
  get_qfile_of_data(qar, qfile, qsinfo);
  const int code = qfseek(qar.qfile, qsinfo.offset_end, SEEK_SET);
  qassert(code == 0);
  return true;
}

inline void read_through(QarFile& qar)
{
  qassert(qar.qfile.mode == "r");
  std::string fn;
  QFile qfile;
  while (true) {
    const bool b = read_next(qar, fn, qfile);
    if (not b) {
      break;
    }
  }
}

inline bool read(QarFile& qar, const std::string& fn, QFile& qfile_in)
// interface function
{
  qassert(qar.qfile.mode == "r");
  qassert(fn != "");
  if (has(qar.qsinfo_map, fn)) {
    const QarSegmentInfo& qsinfo = qar.qsinfo_map[fn];
    get_qfile_of_data(qar, qfile_in, qsinfo);
    return true;
  }
  std::string fn_read;
  while (true) {
    const bool b = read_next(qar, fn_read, qfile_in);
    if (not b) {
      return false;
    }
    if (fn == fn_read) {
      return true;
    }
  }
}

inline std::vector<std::string> list(QarFile& qar)
// interface function
{
  if (qar.null()) {
    return std::vector<std::string>();
  }
  qassert(qar.qfile.mode == "r");
  read_through(qar);
  return qar.fn_list;
}

inline std::vector<std::string> list_qar(const std::string& path)
{
  TIMER_VERBOSE("list_qar");
  QarFile qar(path, "r");
  return list(qar);
}

inline int truncate_qar_file(const std::string& path,
                             const std::vector<std::string>& fns_keep)
// interface function
// return nonzero if failed.
// return 0 if truncated successfully.
// if fns_keep is empty, the resulting qar file should have and only have qar_header.
{
  TIMER_VERBOSE("truncate_qar_file");
  QarFile qar(path, "r");
  if (qar.null()) {
    if (fns_keep.size() == 0) {
      qar.init(path, "w");
      return 0;
    } else {
      qwarn(fname + ssprintf(": fns_keep.size()=%ld", fns_keep.size()));
      return 1;
    }
  }
  const std::vector<std::string> fns = list(qar);
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
  if (fns_keep.size() > 0) {
  }
  std::string fn_last = "";
  if (fns_keep.size() > 0) {
    fn_last = fns_keep.back();
  }
  const long offset_final = qar.qsinfo_map[fn_last].offset_end;
  qar.close();
  const bool b = qtruncate(path, offset_final);
  if (not b) {
    qwarn(fname +
          ssprintf(": fns.size()=%ld fns_keep.size()=%ld offset_final=%ld",
                   fns.size(), fns_keep.size(), offset_final));
    return 3;
  }
  return 0;
}

inline std::vector<std::string> properly_truncate_qar_file(const std::string& path)
// interface function
// The resulting qar file should at least have qar_header.
// Should call this function before append.
{
  std::vector<std::string> fns_keep;
  QarFile qar(path, "r");
  if (qar.null()) {
    qar.init(path, "w");
    qar.close();
    return fns_keep;
  }
  fns_keep = list(qar);
  std::string fn_last = "";
  if (fns_keep.size() > 0) {
    fn_last = fns_keep.back();
  }
  const long offset_final = qar.qsinfo_map[fn_last].offset_end;
  qar.close();
  const bool b = qtruncate(path, offset_final);
  qassert(b);
  return fns_keep;
}

inline void write_start(QarFile& qar, const std::string& fn,
                        const std::string& info, QFile& qfile_out,
                        const long data_len = -1, const long header_len = 60)
// interface function
// Initial pos should be the end of the qar
// Set the qfile_out to be a writable QFile to qar.
// When the final size of qfile_out is unknow (data_len == -1), header_len is
// reserved for header.
// Should call write_end(qar) after writing to qfile_out is finished.
{
  qassert(qar.current_write_segment_offset = -1);
  qar.current_write_segment_offset = qftell(qar.qfile);
  qfseek(qar.qfile, 0, SEEK_END);
  qassert(qftell(qar.qfile) == qar.current_write_segment_offset);
  const std::string header_prefix = "QAR-FILE ";
  std::string header;
  header = ssprintf("%ld %ld %ld", fn.size(), info.size(), data_len);
  if (data_len < 0) {
    qassert(data_len == -1);
    qassert(header_len - (long)header_prefix.size() >= (long)header.size());
    const std::string header_pad(
        header_len - header_prefix.size() - header.size(), ' ');
    header = header_pad + header;
  }
  header = header_prefix + header;
  std::string meta;
  meta += header;
  meta += "\n";
  meta += fn;
  meta += "\n";
  meta += info;
  meta += "\n";
  qwrite_data(meta, qar.qfile);
  qar.current_write_segment_offset_data = qftell(qar.qfile);
  qar.current_write_segment_offset_header_len = header.size();
  qar.current_write_segment_offset_fn_len = fn.size();
  qar.current_write_segment_offset_info_len = info.size();
  qar.current_write_segment_offset_data_len = data_len;
  const long offset_start = qar.current_write_segment_offset_data;
  const long offset_end = data_len == -1 ? -1 : offset_start + data_len;
  qfile_out.init(qar.qfile, offset_start, offset_end);
}

inline void write_end(QarFile& qar)
// interface function
// Call after finish writing to a qfile set by write_start (qfile should be
// closed already).
// Use the end of file as the end of the data.
// Will check / add data_len information in header.
// Finally, will write "\n\n" after the end of file.
{
  qfseek(qar.qfile, 0, SEEK_END);
  const long offset_end = qftell(qar.qfile);
  qassert(qar.current_write_segment_offset >= 0);
  qassert(qar.current_write_segment_offset_data >=
          qar.current_write_segment_offset);
  qfseek(qar.qfile, qar.current_write_segment_offset, SEEK_SET);
  long data_len = qar.current_write_segment_offset_data_len;
  if (data_len >= 0) {
    qassert(qar.current_write_segment_offset_data + data_len == offset_end);
  } else {
    qassert(data_len == -1);
    const long header_len = qar.current_write_segment_offset_header_len;
    const long fn_len = qar.current_write_segment_offset_fn_len;
    const long info_len = qar.current_write_segment_offset_info_len;
    data_len = offset_end - qar.current_write_segment_offset_data;
    qassert(data_len >= 0);
    const std::string header_prefix = "QAR-FILE ";
    std::string header = ssprintf("%ld %ld %ld", fn_len, info_len, data_len);
    qassert(header_len - (long)header_prefix.size() >= (long)header.size());
    const std::string header_pad(
        header_len - header_prefix.size() - header.size(), ' ');
    header = header_prefix + header_pad + header;
    qassert(header.size() == header_len);
    qfseek(qar.qfile, qar.current_write_segment_offset, SEEK_SET);
    qwrite_data(header, qar.qfile);
  }
  qfseek(qar.qfile, 0, SEEK_END);
  qwrite_data("\n\n", qar.qfile);
  qar.current_write_segment_offset = -1;
}

inline long write_from_qfile(QarFile& qar, const std::string& fn,
                             const std::string& info, QFile& qfile_in)
// interface function
// Write content (start from the current position) of qfile_in to qar.
// NOTE: different from write_start and write_end functions!
{
  TIMER_FLOPS("write_from_qfile");
  const long offset_start = qftell(qfile_in);
  qfseek(qfile_in, 0, SEEK_END);
  const long offset_end = qftell(qfile_in);
  qfseek(qfile_in, offset_start, SEEK_SET);
  const long data_len = offset_end - offset_start;
  qassert(data_len >= 0);
  QFile qfile_out;
  write_start(qar, fn, info, qfile_out, data_len);
  const long total_bytes = write_from_qfile(qfile_out, qfile_in);
  write_end(qar);
  timer.flops += total_bytes;
  return total_bytes;
}

inline int qar_create(const std::string& path_qar,
                      const std::string& path_folder_,
                      const bool is_remove_folder_after = false)
// interface function
// return 0 if successful.
{
  TIMER_VERBOSE_FLOPS("qar_create");
  const std::string path_folder = remove_trailing_slashes(path_folder_);
  if (not is_directory(path_folder)) {
    qwarn(fname + ssprintf(": '%s' '%s' no folder.", path_qar.c_str(),
                           path_folder.c_str()));
    return 1;
  }
  if (does_file_exist(path_qar)) {
    qwarn(fname + ssprintf(": '%s' '%s' qar already exist.", path_qar.c_str(),
                           path_folder.c_str()));
    return 2;
  }
  if (does_file_exist(path_qar + ".acc")) {
    qwarn(fname + ssprintf(": '%s' '%s' qar.acc already exist.",
                           path_qar.c_str(), path_folder.c_str()));
    return 3;
  }
  const std::string path_prefix = path_folder + "/";
  const long path_prefix_len = path_prefix.size();  // including the final "/"
  const std::vector<std::string> contents = qls_all(path_folder);
  std::vector<std::string> reg_files;
  for (long i = 0; i < (long)contents.size(); ++i) {
    const std::string path = contents[i];
    qassert(path.substr(0, path_prefix_len) == path_prefix);
    if (not is_directory(path)) {
      if (not is_regular_file(path)) {
        qwarn(fname + ssprintf(": '%s' '%s' '%s' not regular file.",
                               path_qar.c_str(), path_folder.c_str(),
                               path.c_str()));
        return 4;
      }
      reg_files.push_back(path);
    }
  }
  QarFile qar(path_qar + ".acc", "w");
  for (long i = 0; i < (long)reg_files.size(); ++i) {
    const std::string path = reg_files[i];
    const std::string fn = path.substr(path_prefix_len);
    if (verbose_level() > 0) {
      displayln(fname + ssprintf(": '%s' '%s'/'%s' %ld/%ld.", path_qar.c_str(),
                                 path_prefix.c_str(), fn.c_str(), i + 1,
                                 reg_files.size()));
    }
    QFile qfile_in(path, "r");
    timer.flops += write_from_qfile(qar, fn, "", qfile_in);
  }
  qar.close();
  qrename(path_qar + ".acc", path_qar);
  if (is_remove_folder_after) {
    qremove_all(path_folder);
  }
  return 0;
}

inline int qar_extract(const std::string& path_qar,
                       const std::string& path_folder_,
                       const bool is_remove_qar_after = false)
// interface function
// return 0 if successful.
{
  TIMER_VERBOSE_FLOPS("qar_extract");
  const std::string path_folder = remove_trailing_slashes(path_folder_);
  if (not does_file_exist(path_qar)) {
    qwarn(fname + ssprintf(": '%s' '%s' qar does not exist.", path_qar.c_str(),
                           path_folder.c_str()));
    return 1;
  }
  if (does_file_exist(path_folder)) {
    qwarn(fname + ssprintf(": '%s' '%s' folder exist.", path_qar.c_str(),
                           path_folder.c_str()));
    return 2;
  }
  if (does_file_exist(path_folder + ".acc")) {
    qwarn(fname + ssprintf(": '%s' '%s' qar.acc already exist.",
                           path_qar.c_str(), path_folder.c_str()));
    return 3;
  }
  QarFile qar(path_qar, "r");
  const std::vector<std::string> contents = list(qar);
  std::set<std::string> dirs;
  qmkdir_p(path_folder);
  dirs.insert(".");
  for (long i = 0; i < (long)contents.size(); ++i) {
    const std::string& fn = contents[i];
    const std::string dn = dirname(fn);
    if (not has(dirs, dn)) {
      const int code = qmkdir_p(path_folder + "/" + dn);
      qassert(code == 0);
      dirs.insert(dn);
    }
    QFile qfile_in;
    const bool b = read(qar, fn, qfile_in);
    qassert(b);
    QFile qfile_out(path_folder + "/" + fn, "w");
    qassert(not qfile_out.null());
    timer.flops += write_from_qfile(qfile_out, qfile_in);
  }
  qar.close();
  qrename(path_qar + ".acc", path_qar);
  if (is_remove_qar_after) {
    qremove(path_qar);
  }
  return 0;
}

}  // namespace qlat
