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
  long pos;  // position of the Qfile. (correspond to position of fp should be
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
    displayln(
        ssprintf("QFile: open '%s' with '%s'.", path.c_str(), mode.c_str()));
    fp = qopen(path, mode);
    qassert(NULL != fp);
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
  // NOTE: Initial position set to be 0. Does not perform fseek to appropriate position.
  {
    close();
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
  }
  //
  void close()
  {
    // to close the file, it cannot have any child
    qassert(number_of_child == 0);
    if (NULL == parent) {
      displayln(
          ssprintf("QFile: close '%s' with '%s'.", path.c_str(), mode.c_str()));
      qclose(fp);
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
};

inline void qswap(QFile& qfile1, QFile& qfile2) { qfile1.swap(qfile2); }

inline bool qfeof(const QFile& qfile) { return qfile.is_eof; }

inline long qftell(const QFile& qfile) { return qfile.pos; }

inline int qfseek(QFile& qfile, const long q_offset, const int whence)
// Always call fseek and adjust qfile.is_eof and qfile.pos
// qfile.pos will be set to the actual QFile position after qfseek.
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
    const long remaining_size = qfile.offset_end - qfile.offset_start - qfile.pos;
    qassert(remaining_size >= 0);
    const long target_nmemb = std::min(remaining_size / size, nmemb);
    if (target_nmemb == 0) {
      return 0;
    }
    long actual_nmemb = std::fread(ptr, size, target_nmemb, qfile.fp);
    if (target_nmemb < nmemb) {
      qfseek(qfile, 0, SEEK_END);
      qfile.is_eof = true;
    }
    qassert(actual_nmemb == target_nmemb);
  } else {
    long actual_nmemb = std::fread(ptr, size, nmemb, qfile.fp);
    qfile.pos = ftell(qfile.fp) - qfile.offset_start;
    qfile.is_eof = feof(qfile.fp) != 0;
  }
  return actual_nmemb;
}

inline long qfwrite(const void* ptr, const long size, const long nmemb,
                    QFile& qfile)
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
    const long remaining_size = qfile.offset_end - qfile.offset_start - qfile.pos;
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
      const long size_truncate = size - (qfile.offset_start + pos - qfile.offset_end);
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
{
  std::vector<std::string> ret;
  while (not qfeof(qfile)) {
    ret.push_back(qgetline(qfile));
  }
  return ret;
}

const std::string qar_header = "#!/usr/bin/env qar-glimpse\n\n";

struct QarFile {
  QFile qfile;
  //
  bool is_read_through;
  std::vector<std::string> fn_list;
  std::map<std::string, long> offsets_map;
  long max_offset; // maximum offset reached so far
  //
  QarWriter()
  {
    init();
  }
  //
  void init()
  {
    qfile.init();
    is_read_through = false;
    fn_list.clear();
    offsets_map.clear();
    max_offset = 0;
  }
  void init(const std::string& path, const std::string& mode)
  {
    qfile.init(path, mode);
    if (mode == "w") {
      qfwrite(qar_header.data(), qar_header.size(), 1, qfile);
    } else if (mode == "r") {
      std::vector<char> check_line(qar_header.size(), 0);
      const long qfread_check_len =
          qfread(check_line.data(), qar_header.size(), 1, qfile);
      qassert(qfread_check_len == 1);
      qassert(std::string(check_line.data(), check_line.size()) == qar_header);
    }
  }
};

}  // namespace qlat
