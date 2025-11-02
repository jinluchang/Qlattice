#pragma once

#include <qlat-utils/cache.h>
#include <qlat-utils/crc32.h>
#include <qlat-utils/env.h>
#include <qlat-utils/mpi-auto.h>
#include <qlat-utils/timer.h>
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

enum struct QFileMode : Int {
  Read,
  Write,
  Append,
};

enum struct QFileType : Int {
  CFile,
  String,
};

// ---------------------

std::string show(const QFileMode mode);

QFileMode read_qfile_mode(const std::string& mode);

std::string show(const QFileType ftype);

QFileType read_qfile_type(const std::string& ftype);

// ---------------------

struct QFileBase {
  virtual void init() = 0;
  virtual void close() = 0;
  //
  virtual QFileType ftype() const = 0;
  virtual const std::string& path() const = 0;
  virtual QFileMode mode() const = 0;
  virtual bool null() const = 0;
  virtual Long size() const = 0;
  virtual bool eof() const = 0;
  virtual Long tell() const = 0;
  virtual Int flush() const = 0;
  virtual Int seek(const Long offset, const Int whence) = 0;
  virtual Long read(void* ptr, const Long size, const Long nmemb) = 0;
  virtual Long write(const void* ptr, const Long size, const Long nmemb) = 0;
  virtual const std::string& content() = 0;
  //
  virtual ~QFileBase();
  //
  virtual Long remaining_size();
  //
  virtual Long read_data(Vector<Char> v);
  virtual std::string read_all();
  virtual std::string cat();
  virtual std::string getline();
  virtual std::vector<std::string> getlines();
  //
  virtual Long write_data(Vector<Char> v);
  virtual Long write_data(const std::string& v);
  virtual Long write_data(const std::vector<std::string>& v);
  virtual Int append(const std::string& content);
  virtual Int append(const std::vector<std::string>& content);
  //
  virtual Long vprintf(const char* fmt, va_list args);
  virtual Long printf(const char* fmt, ...);
};

// ---------------------

Long write_from_qfile(QFileBase& qfile_out, QFileBase& qfile_in);

// ---------------------

struct QFileObjCFile : QFileBase {
  // can not copy
  //
  std::string path_v;
  QFileMode mode_v;
  FILE* fp;
  Long file_size;
  //
  QFileObjCFile(const std::string& path_, const QFileMode mode_);
  void init(const std::string& path_, const QFileMode mode_);
  //
  QFileObjCFile();
  ~QFileObjCFile();
  QFileObjCFile(const QFileObjCFile&) = delete;
  QFileObjCFile& operator=(const QFileObjCFile&) = delete;
  //
  void init();
  void close();
  QFileType ftype() const;
  const std::string& path() const;
  QFileMode mode() const;
  bool null() const;
  Long size() const;
  bool eof() const;
  Long tell() const;
  Int flush() const;
  Int seek(const Long offset, const Int whence);
  Long read(void* ptr, const Long size, const Long nmemb);
  Long write(const void* ptr, const Long size, const Long nmemb);
  const std::string& content();
};

// ---------------------

struct QFileObjString : QFileBase {
  // can not copy
  //
  std::string path_v;
  QFileMode mode_v;
  bool is_null;
  std::string content_v;
  Long pos;
  bool is_eof;
  //
  QFileObjString(const std::string& path_, const QFileMode mode_);
  QFileObjString(const std::string& path_, const QFileMode mode_,
                 std::string& content_);
  void init(const std::string& path_, const QFileMode mode_);
  void init(const std::string& path_, const QFileMode mode_,
            std::string& content_);
  //
  QFileObjString();
  ~QFileObjString();
  QFileObjString(const QFileObjString&) = delete;
  QFileObjString& operator=(const QFileObjString&) = delete;
  //
  void init();
  void close();
  QFileType ftype() const;
  const std::string& path() const;
  QFileMode mode() const;
  bool null() const;
  Long size() const;
  bool eof() const;
  Long tell() const;
  Int flush() const;
  Int seek(const Long offset, const Int whence);
  Long read(void* ptr, const Long size, const Long nmemb);
  Long write(const void* ptr, const Long size, const Long nmemb);
  const std::string& content();
};

// ---------------------

struct QFileObj : QFileBase {
  // Interface to a `QFileBase` object which allow a view of a portion of the
  // file specified by offset_start and offset_end. The view can be nested.
  //
  std::shared_ptr<QFileBase> fp;
  //
  std::shared_ptr<QFileObj>
      parent;  // If parent.null(), then the `fp` pointer own the `QFileBase`
               // object and is responsible for closing it.
  Long number_of_child;  // Can close `fp` only when number_of_child == 0.
  //
  Long pos;     // position of this `QFileObj`. (correspond to position of `fp`
                // should be pos + offset_start).
  bool is_eof;  // the eof state of QFileObj.
                // NOTE: may not match with eof state of fp.
  // NOTE: Actual fp position may be adjust elsewhere and does not
  // match this pos. When performing operations, always fseek fp to
  // location indicated by pos first.
  //
  Long offset_start;  // start offset of `fp` for `QFileObj`
  Long offset_end;    // end offset of `fp` for `QFileObj` (-1 if not limit,
                      // useful when writing)
  //
  Long file_size;
  //
  QFileObj(const QFileType ftype_, const std::string& path_,
           const QFileMode mode_);
  QFileObj(const QFileType ftype_, const std::string& path_,
           const QFileMode mode_, std::string& content_);
  QFileObj(const std::shared_ptr<QFileObj>& qfile, const Long q_offset_start,
           const Long q_offset_end);
  void init(const QFileType ftype_, const std::string& path_,
            const QFileMode mode_);
  void init(const QFileType ftype_, const std::string& path_,
            const QFileMode mode_, std::string& content_);
  void init(const std::shared_ptr<QFileObj>& qfile, const Long q_offset_start,
            const Long q_offset_end);
  //
  QFileObj(QFileObj&& qfile) noexcept;
  //
  QFileObj();
  ~QFileObj();
  QFileObj(const QFileObj&) = delete;
  QFileObj& operator=(const QFileObj&) = delete;
  //
  void init();
  void close();
  QFileType ftype() const;
  const std::string& path() const;
  QFileMode mode() const;
  bool null() const;
  Long size() const;
  bool eof() const;
  Long tell() const;
  Int flush() const;
  Int seek(const Long q_offset, const Int whence);
  Long read(void* ptr, const Long size, const Long nmemb);
  Long write(const void* ptr, const Long size, const Long nmemb);
  const std::string& content();
};

using QFileMap = std::map<Long, std::weak_ptr<QFileObj>>;

API inline QFileMap& get_all_qfile()
// get_all_qfile()[key] -> std::weak_ptr<QFileObj>
// key = (Long)&qfile_internal
{
  static QFileMap all_qfile;
  return all_qfile;
}

// ---------------------

struct API QFile : QFileBase {
  // Smart pointer to `QFileObj` which allow a view of a portion of the file
  // specified by offset_start and offset_end. The view can be nested.
  //
  std::shared_ptr<QFileObj> p;
  //
  QFile(const std::weak_ptr<QFileObj>& wp);
  QFile(const std::string& path, const QFileMode mode);
  QFile(const QFileType ftype, const std::string& path, const QFileMode mode);
  QFile(const QFileType ftype, const std::string& path, const QFileMode mode,
        std::string& content);
  QFile(const QFile& qfile, const Long q_offset_start, const Long q_offset_end);
  //
  void init(const std::weak_ptr<QFileObj>& wp);
  void init(const std::string& path, const QFileMode mode);
  void init(const QFileType ftype, const std::string& path,
            const QFileMode mode);
  void init(const QFileType ftype, const std::string& path,
            const QFileMode mode, std::string& content);
  void init(const QFile& qfile, const Long q_offset_start,
            const Long q_offset_end);
  //
  QFile();
  //
  void init();
  void close();
  QFileType ftype() const;
  const std::string& path() const;
  QFileMode mode() const;
  bool null() const;
  Long size() const;
  bool eof() const;
  Long tell() const;
  Int flush() const;
  Int seek(const Long offset, const Int whence);
  Long read(void* ptr, const Long size, const Long nmemb);
  Long write(const void* ptr, const Long size, const Long nmemb);
  const std::string& content();
};

// ---------------------

Int clean_up_qfile_map();

std::vector<std::string> show_all_qfile();

std::string show(const QFileObj& qfile);

void qswap(QFileObj& qfile1, QFileObj& qfile2);

// ---------------------

std::string show(const QFile& qfile);

void qswap(QFile& qfile1, QFile& qfile2);

QFile qfopen(const std::string& path, const QFileMode mode);

QFile qfopen(const QFileType ftype, const std::string& path,
             const QFileMode mode);

QFile qfopen(const QFileType ftype, const std::string& path,
             const QFileMode mode, std::string& content);

QFile qfopen(const std::string& path, const std::string& mode);

void qfclose(QFile& qfile);

bool qfeof(const QFile& qfile);

Long qftell(const QFile& qfile);

Int qfflush(const QFile& qfile);

Int qfseek(QFile& qfile, const Long offset, const Int whence);

Int qfseek_set(QFile& qfile, const Long offset);

Int qfseek_end(QFile& qfile, const Long offset);

Int qfseek_cur(QFile& qfile, const Long offset);

Long qfile_size(QFile& qfile);

Long qfile_remaining_size(QFile& qfile);

// ---------------------

Long qfread(void* ptr, const Long size, const Long nmemb, QFile& qfile);

Long qfwrite(const void* ptr, const Long size, const Long nmemb,
             QFile& qfile);

Long qvfprintf(QFile& qfile, const char* fmt, va_list args);

Long qfprintf(QFile& qfile, const char* fmt, ...);

std::string qgetline(QFile& qfile);

std::vector<std::string> qgetlines(QFile& qfile);

template <class M>
Long qwrite_data(const Vector<M>& v, QFile& qfile)
{
  return qwrite_data(get_data_char(v), qfile);
}

Long qwrite_data(const Vector<Char>& v, QFile& qfile);

Long qwrite_data(const std::string& v, QFile& qfile);

Long qwrite_data(const std::vector<std::string>& v, QFile& qfile);

template <class M>
Long qread_data(const Vector<M>& v, QFile& qfile)
// interface function
{
  return qread_data(get_data_char(v), qfile);
}

Long qread_data(const Vector<Char>& v, QFile& qfile);

template <class M>
Long qread_data_all(std::vector<M>& v, QFile& qfile)
// interface function
// Read all the remaining data.
// (Remaining size must be multiple of sizeof(M) otherwise will fail.)
// return total bytes read.
{
  TIMER_FLOPS("qread_data_all(v,qfile)");
  qassert(not qfile.null());
  const Long data_size = qfile_remaining_size(qfile);
  qassert(data_size >= 0);
  const Long n = data_size / sizeof(M);
  qassert(data_size == sizeof(M) * n);
  v.resize(n);
  const Long data_size_read = qread_data(get_data_char(v), qfile);
  qassert(data_size_read == data_size);
  timer.flops += data_size;
  return data_size;
}

std::string qcat(QFile& qfile);

Int qappend(QFile& qfile, const std::string& content);

Int qappend(QFile& qfile, const std::vector<std::string>& content);

// -------------------

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
  //
  void update_offset();
  bool check_offset() const;
};

struct QarFileIndex {
  std::vector<Long> vol_idx_vec;
  std::vector<std::string> fn_vec;
  std::vector<QarSegmentInfo> qsinfo_vec;
  //
  QarFileIndex() { init(); }
  //
  void init()
  {
    vol_idx_vec.clear();
    fn_vec.clear();
    qsinfo_vec.clear();
  }
  //
  bool check() const
  {
    return (fn_vec.size() == vol_idx_vec.size()) and
           (qsinfo_vec.size() == vol_idx_vec.size());
  }
  //
  Long size() const
  {
    qassert(check());
    return fn_vec.size();
  }
};

struct QarFileVolObj {
  QFile qfile;
  //
  bool is_read_through;  // false if in write/append mode
  //
  std::vector<std::string> fn_list;  // update through register_file
  std::map<std::string, QarSegmentInfo>
      qsinfo_map;                     // update through register_file
  std::set<std::string> directories;  // update through register_file
  Long max_offset;  // when read, maximum offset reached so far, update through
                    // register_file
  //
  std::string current_write_segment_fn;  // when write the fn of the current
                                         // working segment
  QarSegmentInfo current_write_segment_info;  // when write, the info of the
                                              // current working segment
  //
  QarFileVolObj() { init(); }
  QarFileVolObj(const std::string& path, const QFileMode mode)
  {
    init(path, mode);
  }
  QarFileVolObj(const QFile& qfile) { init(qfile); }
  ~QarFileVolObj() { close(); };
  //
  void init();
  void init(const std::string& path, const QFileMode mode,
            const Long vol_idx = 0,
            const QarFileIndex& qar_index = QarFileIndex());
  void init(const QFile& qfile_);
  //
  void close();
  //
  bool null() const { return qfile.null(); }
  //
  Int flush() const { return qfile.flush(); }
  //
  const std::string& path() const { return qfile.path(); }
  //
  QFileMode mode() const { return qfile.mode(); }
  //
  Long size() const { return fn_list.size(); }
};

struct API QarFileVol {
  std::shared_ptr<QarFileVolObj> p;
  //
  QarFileVol() { init(); }
  QarFileVol(const std::string& path, const QFileMode mode)
  {
    init(path, mode);
  }
  QarFileVol(const QFile& qfile) { init(qfile); }
  //
  void init() { p = nullptr; }
  void init(const std::string& path, const QFileMode mode,
            const Long vol_idx = 0,
            const QarFileIndex& qar_index = QarFileIndex());
  void init(const QFile& qfile);
  //
  void close();
  bool null() const;
  Int flush() const;
  //
  const std::string& path() const;
  //
  QFileMode mode() const;
  //
  Long size() const;
  //
  QFile& qfile() const;
};

// -------------------

std::vector<std::string> list(const QarFileVol& qar);

bool has_regular_file(const QarFileVol& qar, const std::string& fn);

bool has(const QarFileVol& qar, const std::string& fn);

QFile read(const QarFileVol& qar, const std::string& fn);

std::string read_data(const QarFileVol& qar, const std::string& fn);

std::string read_info(const QarFileVol& qar, const std::string& fn);

void write_start(const QarFileVol& qar, const std::string& fn,
                 const std::string& info, QFile& qfile_out,
                 const Long data_len = -1, const Long header_len = 60);

void write_end(const QarFileVol& qar);

Long write_from_qfile(const QarFileVol& qar, const std::string& fn,
                      const std::string& info, const QFile& qfile_in);

Long write_from_data(const QarFileVol& qar, const std::string& fn,
                     const std::string& info, const Vector<Char> data);

// -------------------

Int truncate_qar_vol_file(const std::string& path,
                          const std::vector<std::string>& fns_keep);

void properly_truncate_qar_vol_file(
    std::vector<std::string>& fn_list,
    std::map<std::string, QarSegmentInfo>& qsinfo_map,
    std::set<std::string>& directories, Long& max_offset,
    const std::string& path, const Long vol_idx, const QarFileIndex& qar_index,
    const bool is_only_check = false);

std::vector<std::string> properly_truncate_qar_vol_file(
    const std::string& path, const bool is_only_check = false);

// -------------------

struct API QarFile : std::vector<QarFileVol> {
  std::string path;
  QFileMode mode;
  //
  Long qar_index_size_saved;
  //
  QarFile()
  {
    qar_index_size_saved = 0;
    init();
  }
  QarFile(const std::string& path_qar, const QFileMode mode)
  {
    qar_index_size_saved = 0;
    init(path_qar, mode);
  }
  //
  ~QarFile() { init(); };
  //
  void init();
  void init(const std::string& path_qar, const QFileMode mode);
  //
  void close();
  //
  bool null() const { return size() == 0; }
  //
  Int flush() const;
  //
  Long index_size() const;
  //
  void save_index(const Long max_diff = 0);
};

std::string show(const QarFile& qar);

using QarFileMap = std::map<Long, Handle<QarFile>>;

API inline QarFileMap& get_all_qar_file_map()
// Note: key be the location of the QarFile converted to Long
{
  static QarFileMap qar_map;
  return qar_map;
}

std::vector<std::string> show_all_qar_file();

void add_qar_file(QarFile& qar);

void remove_qar_file(QarFile& qar);

void close_all_qar_file();

// -------------------

API inline Cache<std::string, QarFile>& get_qar_read_cache()
// key should be the path prefix of the contents of the qar file.
// Note: key should end with '/'.
{
  static Cache<std::string, QarFile> cache("QarReadCache", 64, 1);
  return cache;
}

void update_qar_cache_due_to_change_of_qar_file(const std::string& path);

void update_qar_cache_due_to_change_of_directory(const std::string& path);

// -------------------

std::vector<std::string> list(const QarFile& qar);

bool has_regular_file(const QarFile& qar, const std::string& fn);

bool has(const QarFile& qar, const std::string& fn);

QFile read(const QarFile& qar, const std::string& fn);

std::string read_data(const QarFile& qar, const std::string& fn);

std::string read_info(const QarFile& qar, const std::string& fn);

bool verify_index(const QarFile& qar);

Long write_from_qfile(QarFile& qar, const std::string& fn,
                      const std::string& info, QFile& qfile_in);

Long write_from_data(QarFile& qar, const std::string& fn,
                     const std::string& info, const Vector<Char> data);

Long write_from_data(QarFile& qar, const std::string& fn,
                     const std::string& info, const std::string& data);

Long write_from_data(QarFile& qar, const std::string& fn,
                     const std::string& info,
                     const std::vector<std::string>& data);

// -------------------

std::vector<std::string> properly_truncate_qar_file(
    const std::string& path, const bool is_only_check = false);

// -------------------

std::string show_qar_index(const QarFile& qar);

Int save_qar_index(const QarFile& qar, const std::string& fn);

Int parse_qar_index(QarFileIndex& qar_index,
                    const std::string& qar_index_content);

void install_qar_index(const QarFileVol& qar, const Long vol_idx,
                       const QarFileIndex& qar_index);

Int read_qar_index(const QarFile& qar, const std::string& qar_index_content);

// -------------------

bool does_regular_file_exist_qar(const std::string& path);

bool does_file_exist_qar(const std::string& path);

Int qar_build_index(const std::string& path_qar);

Int qar_create(const std::string& path_qar, const std::string& path_folder_,
               const bool is_remove_folder_after = false);

Int qar_extract(const std::string& path_qar, const std::string& path_folder_,
                const bool is_remove_qar_after = false);

Int qcopy_file(const std::string& path_src, const std::string& path_dst);

std::vector<std::string> list_qar(const std::string& path);

std::string qcat(const std::string& path);

std::vector<std::string> qgetlines(const std::string& fn);

Int qtouch(const std::string& path);

Int qtouch(const std::string& path, const std::string& content);

Int qtouch(const std::string& path, const std::vector<std::string>& content);

Int qappend(const std::string& path, const std::string& content);

Int qappend(const std::string& path, const std::vector<std::string>& content);

// -------------------

DataTable qload_datatable_serial(QFile& qfile);

DataTable qload_datatable_par(QFile& qfile);

DataTable qload_datatable_serial(const std::string& path);

DataTable qload_datatable_par(const std::string& path);

DataTable qload_datatable(const std::string& path, const bool is_par = false);

// -------------------

crc32_t compute_crc32(QFile& qfile);

crc32_t compute_crc32(const std::string& path);

std::pair<std::string, crc32_t> check_file_crc32(const std::string& fn);

std::vector<std::pair<std::string, crc32_t>> check_all_files_crc32(
    const std::string& path);

std::string show_file_crc32(const std::pair<std::string, crc32_t>& fcrc);

std::string show_files_crc32(
    const std::vector<std::pair<std::string, crc32_t>>& fcrcs);

void check_all_files_crc32_info(const std::string& path);

// -------------------

Int qar_build_index_info(const std::string& path_qar);

Int qar_create_info(const std::string& path_qar,
                    const std::string& path_folder_,
                    const bool is_remove_folder_after = false);

Int qar_extract_info(const std::string& path_qar,
                     const std::string& path_folder_,
                     const bool is_remove_qar_after = false);

Int qcopy_file_info(const std::string& path_src, const std::string& path_dst);

std::string qcat_info(const std::string& path);

Int qtouch_info(const std::string& path);

Int qtouch_info(const std::string& path, const std::string& content);

Int qtouch_info(const std::string& path,
                const std::vector<std::string>& content);

Int qappend_info(const std::string& path, const std::string& content);

Int qappend_info(const std::string& path,
                 const std::vector<std::string>& content);

// -------------------

bool does_regular_file_exist_qar_sync_node(const std::string& fn);

bool does_file_exist_qar_sync_node(const std::string& fn);

Int qar_create_sync_node(const std::string& path_qar,
                         const std::string& path_folder_,
                         const bool is_remove_folder_after = false);

Int qar_extract_sync_node(const std::string& path_qar,
                          const std::string& path_folder_,
                          const bool is_remove_qar_after = false);

Int qcopy_file_sync_node(const std::string& path_src,
                         const std::string& path_dst);

std::string qcat_sync_node(const std::string& path);

DataTable qload_datatable_sync_node(const std::string& path,
                                    const bool is_par = false);

}  // namespace qlat
