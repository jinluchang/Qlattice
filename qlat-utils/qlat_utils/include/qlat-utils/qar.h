#pragma once

#include <stdint.h>
#include <qlat-utils/env.h>
#include <qlat-utils/qutils-vec.h>
#include <qlat-utils/qutils-io.h>

#include <cassert>
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <memory>

namespace qlat
{  //

struct QFileInternal;

typedef std::map<long, std::weak_ptr<QFileInternal> > QFileMap;

API inline QFileMap& get_all_qfile()
// get_all_qfile()[key] -> std::weak_ptr<QFileInternal>
// key = (long)&qfile_internal
{
  static QFileMap all_qfile;
  return all_qfile;
}

struct API QFile {
  // Interface to FILE* which allow a view of a portion of the file specified by
  // offset_start and offset_end.
  // The view can be nested.
  //
  std::shared_ptr<QFileInternal> p;
  //
  QFile() { init(); }
  QFile(const std::weak_ptr<QFileInternal>& wp) { init(wp); }
  QFile(const std::string& path, const std::string& mode) { init(path, mode); }
  QFile(const QFile& qfile, const long q_offset_start, const long q_offset_end)
  {
    init(qfile, q_offset_start, q_offset_end);
  }
  //
  void init() { p = nullptr; }
  void init(const std::weak_ptr<QFileInternal>& wp)
  {
    p = std::shared_ptr<QFileInternal>(wp);
  }
  void init(const std::string& path, const std::string& mode);
  void init(const QFile& qfile, const long q_offset_start,
            const long q_offset_end);
  //
  void close();
  //
  bool null() const { return p == nullptr; }
  //
  const std::string& path() const;
  //
  const std::string& mode() const;
  //
  FILE* get_fp() const;
};

inline void add_qfile(const QFile& qfile)
{
  QFileMap& qfile_map = get_all_qfile();
  const long key = (long)qfile.p.get();
  qassert(not has(qfile_map, key));
  qfile_map[key] = qfile.p;
}

inline void remove_qfile(const QFileInternal& qfile_internal)
{
  QFileMap& qfile_map = get_all_qfile();
  const long key = (long)&qfile_internal;
  qassert(has(qfile_map, key));
  qfile_map.erase(key);
}

inline QFile get_qfile(const QFileInternal& qfile_internal)
{
  QFileMap& qfile_map = get_all_qfile();
  const long key = (long)&qfile_internal;
  qassert(has(qfile_map, key));
  return QFile(qfile_map[key]);
}

struct QFileInternal {
  // Interface to FILE* which allow a view of a portion of the file specified by
  // offset_start and offset_end.
  // The view can be nested.
  //
  std::string path;
  std::string mode;  // can be "r", "a", "w"
  FILE* fp;
  //
  QFile parent;  // If parent.null(), then this QFileInternal own the fp pointer
                 // and will be responsible for close it.
  long number_of_child;  // Can close the FILE only when number_of_child == 0.
  //
  bool is_eof;  // the eof state of QFileInternal.
                // NOTE: may not match with eof state of fp.
  long pos;     // position of the QFileInternal. (correspond to position of fp
                // should be pos + offset_start).
  // NOTE: Actual fp position may be adjust elsewhere and does not
  // match this pos. When performing operations, always fseek fp to
  // location indicated by pos first.
  //
  long offset_start;  // start offset of fp for QFileInternal
  long offset_end;    // end offset of fp for QFileInternal (-1 if not limit,
                      // useful when writing)
  //
  QFileInternal()
  {
    fp = NULL;
    number_of_child = 0;
    init();
  }
  QFileInternal(const std::string& path_, const std::string& mode_)
  {
    fp = NULL;
    number_of_child = 0;
    init(path_, mode_);
  }
  QFileInternal(const QFile& qfile, const long q_offset_start,
                const long q_offset_end)
  {
    fp = NULL;
    number_of_child = 0;
    init(qfile, q_offset_start, q_offset_end);
  }
  //
  QFileInternal(const QFileInternal&) = delete;
  //
  QFileInternal(QFileInternal&& qfile) noexcept
  {
    fp = NULL;
    number_of_child = 0;
    init();
    swap(qfile);
  }
  //
  ~QFileInternal()
  {
    close();
    remove_qfile(*this);
  }
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
    displayln_info(
        1, ssprintf("QFile: open '%s' with '%s'.", path.c_str(), mode.c_str()));
    if (mode == "r" and (not is_regular_file(path))) {
      qwarn(ssprintf("QFile: open '%s' with '%s' not regular file.",
                     path.c_str(), mode.c_str()));
    }
    fp = qopen(path, mode);
    if (fp == NULL) {
      qwarn(ssprintf("QFile: open '%s' with '%s' failed.", path.c_str(),
                     mode.c_str()));
    } else {
      pos = ftell(fp);
    }
    is_eof = false;
    offset_start = 0;
    offset_end = -1;
  }
  void init(const QFile& qfile, const long q_offset_start,
            const long q_offset_end)
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
    qfile.p->number_of_child += 1;
    path = qfile.p->path;
    mode = qfile.p->mode;
    fp = qfile.p->fp;
    parent = qfile;
    qassert(q_offset_start >= 0);
    is_eof = false;
    pos = 0;
    offset_start = qfile.p->offset_start + q_offset_start;
    if (q_offset_end == -1) {
      offset_end = qfile.p->offset_end;
    } else {
      qassert(q_offset_end >= q_offset_start);
      offset_end = qfile.p->offset_start + q_offset_end;
      if (qfile.p->offset_end != -1) {
        qassert(offset_end <= qfile.p->offset_end);
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
    if (parent.null()) {
      if (fp != NULL) {
        displayln_info(1, ssprintf("QFile: close '%s' with '%s'.", path.c_str(),
                                   mode.c_str()));
        qfclose(fp);
      }
    } else {
      fp = NULL;
      parent.p->number_of_child -= 1;
      parent = QFile();
    }
    qassert(fp == NULL);
    qassert(parent.null());
  }
  //
  void swap(QFileInternal& qfile)
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
  bool null() const { return fp == NULL; }
};

inline void QFile::init(const std::string& path, const std::string& mode)
{
  if (p == nullptr) {
    p = std::shared_ptr<QFileInternal>(new QFileInternal());
    add_qfile(*this);
  }
  p->init(path, mode);
}

inline void QFile::init(const QFile& qfile, const long q_offset_start,
                        const long q_offset_end)
{
  if (p == nullptr) {
    p = std::shared_ptr<QFileInternal>(new QFileInternal());
    add_qfile(*this);
  }
  p->init(qfile, q_offset_start, q_offset_end);
}

inline void QFile::close()
{
  if (p != nullptr) {
    p->close();
    p = nullptr;
  }
}

inline const std::string& QFile::path() const { return p->path; }

inline const std::string& QFile::mode() const { return p->mode; }

inline FILE* QFile::get_fp() const { return p->fp; }

inline void qfclose(QFile& qfile)
// interface function
{
  qfile.close();
}

inline void qswap(QFile& qfile1, QFile& qfile2)
// interface function
{
  std::swap(qfile1, qfile2);
}

inline bool qfeof(const QFile& qfile)
// interface function
{
  qassert(not qfile.null());
  return qfile.p->is_eof;
}

inline long qftell(const QFile& qfile)
// interface function
{
  qassert(not qfile.null());
  return qfile.p->pos;
}

inline int qfflush(const QFile& qfile)
// interface function
{
  qassert(not qfile.null());
  return fflush(qfile.get_fp());
}

int qfseek(const QFile& qfile, const long q_offset, const int whence);

long qfread(void* ptr, const long size, const long nmemb, const QFile& qfile);

long qfwrite(const void* ptr, const long size, const long nmemb,
             const QFile& qfile);

int qvfscanf(const QFile& qfile, const char* fmt, va_list args);

int qfscanf(const QFile& qfile, const char* fmt, ...);

std::string qgetline(const QFile& qfile);

inline std::vector<std::string> qgetlines(const QFile& qfile)
// interface function
{
  qassert(not qfile.null());
  std::vector<std::string> ret;
  while (not qfeof(qfile)) {
    ret.push_back(qgetline(qfile));
  }
  return ret;
}

inline long qfile_remaining_size(const QFile& qfile)
// interface function
// qfile should have definite size.
// return the remaining size of qfile (start from the current position).
{
  TIMER_FLOPS("qfile_remaining_size");
  qassert(not qfile.null());
  const long offset_start = qftell(qfile);
  qfseek(qfile, 0, SEEK_END);
  const long offset_end = qftell(qfile);
  qfseek(qfile, offset_start, SEEK_SET);
  const long data_len = offset_end - offset_start;
  qassert(data_len >= 0);
  return data_len;
}

template <class M>
long qwrite_data(const Vector<M>& v, const QFile& qfile)
// interface function
{
  TIMER_FLOPS("qwrite_data(v,qfile)");
  qassert(not qfile.null());
  const long data_size = sizeof(M) * qfwrite((void*)v.p, sizeof(M), v.n, qfile);
  timer.flops += data_size;
  return data_size;
}

inline long qwrite_data(const std::string& line, const QFile& qfile)
// interface function
{
  qassert(not qfile.null());
  return qwrite_data(get_data(line), qfile);
}

inline long qwrite_data(const std::vector<std::string>& lines,
                        const QFile& qfile)
// interface function
{
  qassert(not qfile.null());
  long total_bytes = 0;
  for (long i = 0; i < (long)lines.size(); ++i) {
    total_bytes += qwrite_data(lines[i], qfile);
  }
  return total_bytes;
}

template <class M>
long qread_data(const Vector<M>& v, const QFile& qfile)
// interface function
{
  TIMER_FLOPS("qread_data(v,qfile)");
  qassert(not qfile.null());
  const long data_size = sizeof(M) * qfread((void*)v.p, sizeof(M), v.n, qfile);
  timer.flops += data_size;
  return data_size;
}

template <class M>
long qread_data_all(std::vector<M>& v, const QFile& qfile)
// interface function
// Read all the remaining data.
// (Remaining size must be multiple of sizeof(M) otherwise will fail.)
// return total bytes read.
{
  TIMER_FLOPS("qread_data_all(v,qfile)");
  qassert(not qfile.null());
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

inline long qvfprintf(const QFile& qfile, const char* fmt, va_list args)
{
  const std::string str = vssprintf(fmt, args);
  return qwrite_data(str, qfile);
}

inline long qfprintf(const QFile& qfile, const char* fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  return qvfprintf(qfile, fmt, args);
}

long write_from_qfile(const QFile& qfile_out, const QFile& qfile_in);

// -------------------

const std::string qar_header = "#!/usr/bin/env qar-glimpse\n\n";

const std::string qar_idx_header = "#!/usr/bin/env qar-idx-glimpse\n\n";

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

struct QarFileVolInternal;

struct API QarFileVol {
  std::shared_ptr<QarFileVolInternal> p;
  //
  QarFileVol() { init(); }
  QarFileVol(const std::string& path, const std::string& mode)
  {
    init(path, mode);
  }
  QarFileVol(const QFile& qfile) { init(qfile); }
  //
  void init() { p = nullptr; }
  void init(const std::string& path, const std::string& mode)
  {
    init(QFile(path, mode));
  }
  void init(const QFile& qfile);
  //
  void close();  // may not close the underling file, only release the pointer
                 // from this QarFileVol
  //
  bool null() const { return p == nullptr; }
  //
  const std::string& path() const;
  //
  const std::string& mode() const;
  //
  const QFile& qfile() const;
};

struct QarFileVolInternal {
  QFile qfile;
  //
  bool is_read_through;
  std::vector<std::string> fn_list;
  std::map<std::string, QarSegmentInfo> qsinfo_map;
  std::set<std::string> directories;
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
  QarFileVolInternal() { init(); }
  QarFileVolInternal(const QFile& qfile) { init(qfile); }
  //
  void init()
  {
    qfile.init();
    is_read_through = false;
    fn_list.clear();
    qsinfo_map.clear();
    directories.clear();
    max_offset = 0;
    current_write_segment_offset = -1;
    current_write_segment_offset_data = -1;
    current_write_segment_offset_header_len = -1;
    current_write_segment_offset_fn_len = -1;
    current_write_segment_offset_info_len = -1;
    current_write_segment_offset_data_len = -1;
  }
  void init(const QFile& qfile_)
  {
    init();
    qfile = qfile_;
    if (qfile.null()) {
      return;
    }
    if (mode() == "w") {
      qfwrite(qar_header.data(), qar_header.size(), 1, qfile);
    } else if (mode() == "r") {
      std::vector<char> check_line(qar_header.size(), 0);
      const long qfread_check_len =
          qfread(check_line.data(), qar_header.size(), 1, qfile);
      if (not(qfread_check_len == 1 and
              std::string(check_line.data(), check_line.size()) ==
                  qar_header)) {
        qfile.close();
      };
      max_offset = qftell(qfile);
      directories.insert("");
    }
  }
  //
  void close() { init(); }
  //
  bool null() const { return qfile.null(); }
  //
  const std::string& path() const { return qfile.path(); }
  //
  const std::string& mode() const { return qfile.mode(); }
};

inline void QarFileVol::init(const QFile& qfile)
{
  if (p == nullptr) {
    p = std::shared_ptr<QarFileVolInternal>(new QarFileVolInternal());
  }
  p->init(qfile);
}

inline void QarFileVol::close()
{
  if (p != nullptr) {
    p->close();
    p = nullptr;
  }
}

inline const std::string& QarFileVol::path() const { return p->path(); }

inline const std::string& QarFileVol::mode() const { return p->mode(); }

inline const QFile& QarFileVol::qfile() const { return p->qfile; }

void register_file(const QarFileVol& qar, const std::string& fn,
                   const QarSegmentInfo& qsinfo);

bool read_qar_segment_info(QarFileVolInternal& qar, QarSegmentInfo& qsinfo);

std::string read_fn(const QarFileVol& qar, const QarSegmentInfo& qsinfo);

void read_info(const QarFileVol& qar, std::string& info,
               const QarSegmentInfo& qsinfo);

QFile get_qfile_of_data(const QarFileVol& qar, const QarSegmentInfo& qsinfo);

QFile read_next(const QarFileVol& qar, std::string& fn);

void read_through(const QarFileVol& qar);

QFile read(const QarFileVol& qar, const std::string& fn);

bool has_regular_file(const QarFileVol& qar, const std::string& fn);

bool has(const QarFileVol& qar, const std::string& fn);

std::vector<std::string> list(const QarFileVol& qar);

void write_start(const QarFileVol& qar, const std::string& fn,
                 const std::string& info, QFile& qfile_out,
                 const long data_len = -1, const long header_len = 60);

void write_end(const QarFileVol& qar);

long write_from_qfile(const QarFileVol& qar, const std::string& fn,
                      const std::string& info, const QFile& qfile_in);

int truncate_qar_file(const std::string& path,
                      const std::vector<std::string>& fns_keep);

std::vector<std::string> properly_truncate_qar_file(const std::string& path);

}  // namespace qlat
