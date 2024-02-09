#pragma once

#include <qlat-utils/cache.h>
#include <qlat-utils/crc32.h>
#include <qlat-utils/env.h>
#include <qlat-utils/types.h>
#include <qlat-utils/utils-io.h>
#include <qlat-utils/utils-vec.h>
#include <stdint.h>

#include <cassert>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace qlat
{  //

struct QFileInternal;

typedef std::map<Long, std::weak_ptr<QFileInternal>> QFileMap;

API inline QFileMap& get_all_qfile()
// get_all_qfile()[key] -> std::weak_ptr<QFileInternal>
// key = (Long)&qfile_internal
{
  static QFileMap all_qfile;
  return all_qfile;
}

// ---------------------

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
  QFile(const QFile& qfile, const Long q_offset_start, const Long q_offset_end)
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
  void init(const QFile& qfile, const Long q_offset_start,
            const Long q_offset_end);
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

// ---------------------

inline void add_qfile(const QFile& qfile)
{
  QFileMap& qfile_map = get_all_qfile();
  const Long key = (Long)qfile.p.get();
  qassert(not has(qfile_map, key));
  qfile_map[key] = qfile.p;
}

inline void remove_qfile(const QFileInternal& qfile_internal)
{
  QFileMap& qfile_map = get_all_qfile();
  const Long key = (Long)&qfile_internal;
  qassert(has(qfile_map, key));
  qfile_map.erase(key);
}

inline QFile get_qfile(const QFileInternal& qfile_internal)
{
  QFileMap& qfile_map = get_all_qfile();
  const Long key = (Long)&qfile_internal;
  qassert(has(qfile_map, key));
  return QFile(qfile_map[key]);
}

// ---------------------

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
  Long number_of_child;  // Can close the FILE only when number_of_child == 0.
  //
  bool is_eof;  // the eof state of QFileInternal.
                // NOTE: may not match with eof state of fp.
  Long pos;     // position of the QFileInternal. (correspond to position of fp
                // should be pos + offset_start).
  // NOTE: Actual fp position may be adjust elsewhere and does not
  // match this pos. When performing operations, always fseek fp to
  // location indicated by pos first.
  //
  Long offset_start;  // start offset of fp for QFileInternal
  Long offset_end;    // end offset of fp for QFileInternal (-1 if not limit,
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
  QFileInternal(const QFile& qfile, const Long q_offset_start,
                const Long q_offset_end)
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
  void init(const std::string& path_, const std::string& mode_);
  void init(const QFile& qfile, const Long q_offset_start,
            const Long q_offset_end);
  // Become a child of qfile.
  // NOTE: q_offset_start and q_offset_end are relative offset for qfile not the
  // absolute offset for qfile.fp .
  // q_offset_end == -1 means no additional limit
  // NOTE: Initial position set to be 0. Does not perform fseek to appropriate
  // position.
  //
  void close();
  //
  void swap(QFileInternal& qfile);
  //
  bool null() const { return fp == NULL; }
};

inline const std::string& QFile::path() const { return p->path; }

inline const std::string& QFile::mode() const { return p->mode; }

inline FILE* QFile::get_fp() const { return p->fp; }

// ---------------------

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

inline Long qftell(const QFile& qfile)
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

int qfseek(const QFile& qfile, const Long q_offset, const int whence);

Long qfread(void* ptr, const Long size, const Long nmemb, const QFile& qfile);

Long qfwrite(const void* ptr, const Long size, const Long nmemb,
             const QFile& qfile);

int qvfscanf(const QFile& qfile, const char* fmt, va_list args);

int qfscanf(const QFile& qfile, const char* fmt, ...);

std::string qgetline(const QFile& qfile);

std::vector<std::string> qgetlines(const QFile& qfile);

Long qfile_size(const QFile& qfile);

Long qfile_remaining_size(const QFile& qfile);

template <class M>
Long qwrite_data(const Vector<M>& v, const QFile& qfile)
// interface function
{
  TIMER_FLOPS("qwrite_data(v,qfile)");
  qassert(not qfile.null());
  const Long data_size = sizeof(M) * qfwrite((void*)v.p, sizeof(M), v.n, qfile);
  timer.flops += data_size;
  return data_size;
}

Long qwrite_data(const std::string& line, const QFile& qfile);

Long qwrite_data(const std::vector<std::string>& lines, const QFile& qfile);

template <class M>
Long qread_data(const Vector<M>& v, const QFile& qfile)
// interface function
{
  TIMER_FLOPS("qread_data(v,qfile)");
  qassert(not qfile.null());
  const Long data_size = sizeof(M) * qfread((void*)v.p, sizeof(M), v.n, qfile);
  timer.flops += data_size;
  return data_size;
}

template <class M>
Long qread_data_all(std::vector<M>& v, const QFile& qfile)
// interface function
// Read all the remaining data.
// (Remaining size must be multiple of sizeof(M) otherwise will fail.)
// return total bytes read.
{
  TIMER_FLOPS("qread_data_all(v,qfile)");
  qassert(not qfile.null());
  const Long pos_initial = qftell(qfile);
  qfseek(qfile, 0, SEEK_END);
  const Long pos_final = qftell(qfile);
  const Long data_size = pos_final - pos_initial;
  const Long n = data_size / sizeof(M);
  qassert(data_size == sizeof(M) * n);
  qfseek(qfile, pos_initial, SEEK_SET);
  v.resize(n);
  const Long data_size_read = qread_data(get_data(v), qfile);
  qassert(data_size_read == data_size);
  timer.flops += data_size;
  return data_size;
}

inline Long qvfprintf(const QFile& qfile, const char* fmt, va_list args)
{
  const std::string str = vssprintf(fmt, args);
  return qwrite_data(str, qfile);
}

inline Long qfprintf(const QFile& qfile, const char* fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  return qvfprintf(qfile, fmt, args);
}

Long write_from_qfile(const QFile& qfile_out, const QFile& qfile_in);

// -------------------

const std::string qar_header = "#!/usr/bin/env qar-glimpse\n\n";

const std::string qar_idx_header = "#!/usr/bin/env qar-idx-glimpse\n\n";

struct QarSegmentInfo {
  Long offset;
  Long offset_fn;
  Long offset_info;
  Long offset_data;
  Long offset_end;
  Long fn_len;
  Long info_len;
  Long data_len;
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
  Long max_offset;  // when read, maximum offset reached so far
  //
  Long current_write_segment_offset;  // when write, the offset of the beginning
                                      // of the current working segment.
  Long current_write_segment_offset_data;  // when write, the offset of the
                                           // beginning of the data of the
                                           // current working segment.
  Long current_write_segment_offset_header_len;
  Long current_write_segment_offset_fn_len;
  Long current_write_segment_offset_info_len;
  Long current_write_segment_offset_data_len;
  //
  QarFileVolInternal() { init(); }
  QarFileVolInternal(const QFile& qfile) { init(qfile); }
  //
  void init();
  void init(const QFile& qfile_);
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
                 const Long data_len = -1, const Long header_len = 60);

void write_end(const QarFileVol& qar);

Long write_from_qfile(const QarFileVol& qar, const std::string& fn,
                      const std::string& info, const QFile& qfile_in);

int truncate_qar_file(const std::string& path,
                      const std::vector<std::string>& fns_keep);

std::vector<std::string> properly_truncate_qar_file(const std::string& path);

// -------------------

bool does_regular_file_exist_qar(const std::string& path);

bool does_file_exist_qar(const std::string& path);

QFile qfopen(const std::string& path, const std::string& mode);

std::string qcat(const std::string& path);

int qar_build_index(const std::string& path_qar);

int qar_create(const std::string& path_qar, const std::string& path_folder_,
               const bool is_remove_folder_after = false);

int qar_extract(const std::string& path_qar, const std::string& path_folder_,
                const bool is_remove_qar_after = false);

int qcopy_file(const std::string& path_src, const std::string& path_dst);

std::string mk_key_from_qar_path(const std::string& path);

std::vector<std::string> list_qar(const std::string& path);

DataTable qload_datatable_serial(QFile& qfile);

DataTable qload_datatable_par(QFile& qfile);

DataTable qload_datatable_serial(const std::string& path);

DataTable qload_datatable_par(const std::string& path);

DataTable qload_datatable(const std::string& path, const bool is_par = false);

int qtouch(const std::string& path);

int qtouch(const std::string& path, const std::string& content);

int qtouch(const std::string& path, const std::vector<std::string>& content);

int qappend(const std::string& path, const std::string& content);

crc32_t compute_crc32(QFile& qfile);

crc32_t compute_crc32(const std::string& path);

std::vector<std::pair<std::string, crc32_t>> check_all_files_crc32(
    const std::string& path);

void check_all_files_crc32_info(const std::string& path);

// -------------------

struct API QarFile;

int save_qar_index(const QarFile& qar, const std::string& fn);

int save_qar_index_info(const QarFile& qar, const std::string& fn);

int parse_qar_index(const QarFile& qar, const std::string& qar_index_content);

int load_qar_index(const QarFile& qar, const std::string& fn);

// -------------------

inline std::string qar_file_multi_vol_suffix(const Long i)
{
  if (i == 0) {
    return "";
  } else {
    return ssprintf(".v%ld", i);
  }
  qassert(false);
  return "";
}

struct API QarFile : std::vector<QarFileVol> {
  // Only for reading
  QarFile() { init(); }
  QarFile(const std::string& path_qar, const std::string& mode)
  {
    init(path_qar, mode);
  }
  //
  void init()
  {
    std::vector<QarFileVol>& v = *this;
    qlat::clear(v);
  }
  void init(const std::string& path_qar, const std::string& mode);
  //
  void close() { init(); }
  //
  bool null() const { return size() == 0; }
};

std::vector<std::string> list(const QarFile& qar);

bool has_regular_file(const QarFile& qar, const std::string& fn);

bool has(const QarFile& qar, const std::string& fn);

QFile read(const QarFile& qar, const std::string& fn);

// -------------------

API inline Cache<std::string, QarFile>& get_qar_read_cache()
// key should be the path prefix of the contents of the qar file.
// Note: key should end with '/'.
{
  static Cache<std::string, QarFile> cache("QarReadCache", 64, 1);
  return cache;
}

std::string mk_new_qar_read_cache_key(const QarFile& qar,
                                      const std::string& key,
                                      const std::string& path);

std::string mk_new_qar_read_cache_key(const std::string& path);

std::string get_qar_read_cache_key(const std::string& path);

// -------------------

std::vector<std::string> qgetlines(const std::string& fn);

// -------------------

std::string show_file_crc32(const std::pair<std::string, crc32_t>& fcrc);

std::string show_files_crc32(
    const std::vector<std::pair<std::string, crc32_t>>& fcrcs);

std::pair<std::string, crc32_t> check_file_crc32(const std::string& fn);

// -------------------

std::string qcat_info(const std::string& path);

int qtouch_info(const std::string& path);

int qtouch_info(const std::string& path, const std::string& content);

int qtouch_info(const std::string& path,
                const std::vector<std::string>& content);

int qappend_info(const std::string& path, const std::string& content);

int qar_build_index_info(const std::string& path_qar);

int qar_create_info(const std::string& path_qar,
                    const std::string& path_folder_,
                    const bool is_remove_folder_after = false);

int qar_extract_info(const std::string& path_qar,
                     const std::string& path_folder_,
                     const bool is_remove_qar_after = false);

int qcopy_file_info(const std::string& path_src, const std::string& path_dst);

// -------------------

bool does_regular_file_exist_qar_sync_node(const std::string& fn);

bool does_file_exist_qar_sync_node(const std::string& fn);

int qar_create_sync_node(const std::string& path_qar,
                         const std::string& path_folder_,
                         const bool is_remove_folder_after = false);

int qar_extract_sync_node(const std::string& path_qar,
                          const std::string& path_folder_,
                          const bool is_remove_qar_after = false);

int qcopy_file_sync_node(const std::string& path_src,
                         const std::string& path_dst);

}  // namespace qlat
