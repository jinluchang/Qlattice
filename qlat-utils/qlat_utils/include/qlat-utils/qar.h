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

struct QFileObj;

struct API QFile {
  // Interface to FILE* which allow a view of a portion of the file specified by
  // offset_start and offset_end.
  // The view can be nested.
  //
  std::shared_ptr<QFileObj> p;
  //
  QFile() { init(); }
  QFile(const std::weak_ptr<QFileObj>& wp) { init(wp); }
  QFile(const std::string& path, const std::string& mode) { init(path, mode); }
  QFile(const QFile& qfile, const Long q_offset_start, const Long q_offset_end)
  {
    init(qfile, q_offset_start, q_offset_end);
  }
  //
  void init();
  void init(const std::weak_ptr<QFileObj>& wp);
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

struct QFileObj {
  // Interface to FILE* which allow a view of a portion of the file specified by
  // offset_start and offset_end.
  // The view can be nested.
  //
  std::string path;
  std::string mode;  // can be "r", "a", "w"
  FILE* fp;
  //
  QFile parent;  // If parent.null(), then this QFileObj own the fp pointer
                 // and will be responsible for close it.
  Long number_of_child;  // Can close the FILE only when number_of_child == 0.
  //
  bool is_eof;  // the eof state of QFileObj.
                // NOTE: may not match with eof state of fp.
  Long pos;     // position of the QFileObj. (correspond to position of fp
                // should be pos + offset_start).
  // NOTE: Actual fp position may be adjust elsewhere and does not
  // match this pos. When performing operations, always fseek fp to
  // location indicated by pos first.
  //
  Long offset_start;  // start offset of fp for QFileObj
  Long offset_end;    // end offset of fp for QFileObj (-1 if not limit,
                      // useful when writing)
  //
  QFileObj();
  QFileObj(const std::string& path_, const std::string& mode_);
  QFileObj(const QFile& qfile, const Long q_offset_start,
           const Long q_offset_end);
  //
  QFileObj(const QFileObj&) = delete;
  //
  QFileObj(QFileObj&& qfile) noexcept;
  //
  ~QFileObj();
  //
  void init();
  void init(const std::string& path_, const std::string& mode_);
  void init(const QFile& qfile, const Long q_offset_start,
            const Long q_offset_end);
  //
  void close();
  //
  void swap(QFileObj& qfile);
  //
  bool null() const { return fp == NULL; }
};

using QFileMap = std::map<Long, std::weak_ptr<QFileObj>>;

API inline QFileMap& get_all_qfile()
// get_all_qfile()[key] -> std::weak_ptr<QFileObj>
// key = (Long)&qfile_internal
{
  static QFileMap all_qfile;
  return all_qfile;
}

std::vector<std::string> show_all_qfile();

// ---------------------

std::string show(const QFileObj& qfile);

void qswap(QFileObj& qfile1, QFileObj& qfile2);

// ---------------------

std::string show(const QFile& qfile);

void qswap(QFile& qfile1, QFile& qfile2);

QFile qfopen(const std::string& path, const std::string& mode);

void qfclose(QFile& qfile);

bool qfeof(const QFile& qfile);

Long qftell(const QFile& qfile);

int qfflush(const QFile& qfile);

int qfseek(const QFile& qfile, const Long q_offset, const int whence);

Long qfile_size(const QFile& qfile);

Long qfile_remaining_size(const QFile& qfile);

// ---------------------

Long qfread(void* ptr, const Long size, const Long nmemb, const QFile& qfile);

Long qfwrite(const void* ptr, const Long size, const Long nmemb,
             const QFile& qfile);

int qvfscanf(const QFile& qfile, const char* fmt, va_list args);

int qfscanf(const QFile& qfile, const char* fmt, ...);

Long qvfprintf(const QFile& qfile, const char* fmt, va_list args);

Long qfprintf(const QFile& qfile, const char* fmt, ...);

std::string qgetline(const QFile& qfile);

std::vector<std::string> qgetlines(const QFile& qfile);

template <class M>
Long qwrite_data(const Vector<M>& v, const QFile& qfile)
{
  return qwrite_data(get_data_char(v), qfile);
}

Long qwrite_data(const Vector<char>& v, const QFile& qfile);

Long qwrite_data(const std::string& line, const QFile& qfile);

Long qwrite_data(const std::vector<std::string>& lines, const QFile& qfile);

template <class M>
Long qread_data(const Vector<M>& v, const QFile& qfile)
// interface function
{
  return qread_data(get_data_char(v), qfile);
}

Long qread_data(const Vector<char>& v, const QFile& qfile);

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
  const Long data_size_read = qread_data(get_data_char(v), qfile);
  qassert(data_size_read == data_size);
  timer.flops += data_size;
  return data_size;
}

Long write_from_qfile(const QFile& qfile_out, const QFile& qfile_in);

// -------------------

struct QarFileVolObj;

struct API QarFileVol {
  std::shared_ptr<QarFileVolObj> p;
  //
  QarFileVol() { init(); }
  QarFileVol(const std::string& path, const std::string& mode)
  {
    init(path, mode);
  }
  QarFileVol(const QFile& qfile) { init(qfile); }
  //
  void init() { p = nullptr; }
  void init(const std::string& path, const std::string& mode);
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
  bool check_offset();
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
  QarFileVolObj(const std::string& path, const std::string& mode)
  {
    init(path, mode);
  }
  QarFileVolObj(const QFile& qfile) { init(qfile); }
  //
  void init();
  void init(const std::string& path, const std::string& mode);
  void init(const QFile& qfile_);
  //
  void close();
  //
  bool null() const { return qfile.null(); }
  //
  const std::string& path() const { return qfile.path(); }
  //
  const std::string& mode() const { return qfile.mode(); }
};

// -------------------

QFile read_next(const QarFileVol& qar, std::string& fn);

void read_through(const QarFileVol& qar);

QFile read(const QarFileVol& qar, const std::string& fn);

void read_info(const QarFileVol& qar, std::string& info,
               const QarSegmentInfo& qsinfo);

bool has_regular_file(const QarFileVol& qar, const std::string& fn);

bool has(const QarFileVol& qar, const std::string& fn);

std::vector<std::string> list(const QarFileVol& qar);

void write_start(const QarFileVol& qar, const std::string& fn,
                 const std::string& info, QFile& qfile_out,
                 const Long data_len = -1, const Long header_len = 60);

void write_end(const QarFileVol& qar);

Long write_from_qfile(const QarFileVol& qar, const std::string& fn,
                      const std::string& info, const QFile& qfile_in);

Long write_from_data(const QarFileVol& qar, const std::string& fn,
                     const std::string& info, const Vector<char> data);

// -------------------

int truncate_qar_vol_file(const std::string& path,
                          const std::vector<std::string>& fns_keep);

void properly_truncate_qar_vol_file(
    std::vector<std::string>& fn_list,
    std::map<std::string, QarSegmentInfo>& qsinfo_map,
    std::set<std::string>& directories, Long& max_offset,
    const std::string& path, const bool is_check_all = false,
    const bool is_only_check = false);

std::vector<std::string> properly_truncate_qar_vol_file(
    const std::string& path, const bool is_check_all = false,
    const bool is_only_check = false);

// -------------------

struct API QarFile;

std::vector<std::string> show_qar_index(const QarFile& qar,
                                        const std::string& fn);

int save_qar_index(const QarFile& qar, const std::string& fn);

int save_qar_index_info(const QarFile& qar, const std::string& fn);

int parse_qar_index(std::vector<Long>& vol_idx_vec,
                    std::vector<std::string>& fn_vec,
                    std::vector<QarSegmentInfo>& qsinfo_vec,
                    const std::string& qar_index_content);

int parse_qar_index(const QarFile& qar, const std::string& qar_index_content);

int load_qar_index(const QarFile& qar, const std::string& fn);

// -------------------

struct API QarFile : std::vector<QarFileVol> {
  // Only for reading
  QarFile() { init(); }
  QarFile(const std::string& path_qar, const std::string& mode)
  {
    init(path_qar, mode);
  }
  //
  void init();
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

// -------------------

// -------------------

bool does_regular_file_exist_qar(const std::string& path);

bool does_file_exist_qar(const std::string& path);

int qar_build_index(const std::string& path_qar);

int qar_create(const std::string& path_qar, const std::string& path_folder_,
               const bool is_remove_folder_after = false);

int qar_extract(const std::string& path_qar, const std::string& path_folder_,
                const bool is_remove_qar_after = false);

int qcopy_file(const std::string& path_src, const std::string& path_dst);

std::vector<std::string> list_qar(const std::string& path);

std::string qcat(const std::string& path);

std::vector<std::string> qgetlines(const std::string& fn);

int qtouch(const std::string& path);

int qtouch(const std::string& path, const std::string& content);

int qtouch(const std::string& path, const std::vector<std::string>& content);

int qappend(const std::string& path, const std::string& content);

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

int qar_build_index_info(const std::string& path_qar);

int qar_create_info(const std::string& path_qar,
                    const std::string& path_folder_,
                    const bool is_remove_folder_after = false);

int qar_extract_info(const std::string& path_qar,
                     const std::string& path_folder_,
                     const bool is_remove_qar_after = false);

int qcopy_file_info(const std::string& path_src, const std::string& path_dst);

std::string qcat_info(const std::string& path);

int qtouch_info(const std::string& path);

int qtouch_info(const std::string& path, const std::string& content);

int qtouch_info(const std::string& path,
                const std::vector<std::string>& content);

int qappend_info(const std::string& path, const std::string& content);

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
