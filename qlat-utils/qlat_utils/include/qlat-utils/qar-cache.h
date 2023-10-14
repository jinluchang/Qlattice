#pragma once

#include <qlat-utils/env.h>
#include <qlat-utils/cache.h>
#include <qlat-utils/crc32.h>
#include <qlat-utils/qar.h>
#include <qlat-utils/types.h>

namespace qlat
{  //

bool does_regular_file_exist_qar(const std::string& path);

bool does_file_exist_qar(const std::string& path);

QFile qfopen(const std::string& path, const std::string& mode);

std::string qcat(const std::string& path);

void qar_build_index(const std::string& path_qar);

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

std::vector<std::pair<std::string, crc32_t> > check_all_files_crc32(
    const std::string& path);

void check_all_files_crc32_info(const std::string& path);

// -------------------

struct API QarFile;

void save_qar_index(const QarFile& qar, const std::string& fn);

void save_qar_index_info(const QarFile& qar, const std::string& fn);

void parse_qar_index(const QarFile& qar, const std::string& qar_index_content);

void load_qar_index(const QarFile& qar, const std::string& fn);

// -------------------

inline std::string qar_file_multi_vol_suffix(const long i)
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
  void init(const std::string& path_qar, const std::string& mode)
  {
    init();
    if (mode == "r") {
      // maximally 1024 * 1024 * 1024 volumes
      for (long iv = 0; iv < 1024 * 1024 * 1024; ++iv) {
        const std::string path_qar_v = path_qar + qar_file_multi_vol_suffix(iv);
        if (not does_regular_file_exist_qar(path_qar_v)) {
          break;
        }
        push_back(qfopen(path_qar_v, mode));
        if (back().null()) {
          pop_back();
          break;
        }
      }
      load_qar_index(*this, path_qar + ".idx");
    } else {
      qassert(false);
    }
  }
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

inline std::vector<std::string> qgetlines(const std::string& fn)
{
  QFile qfile = qfopen(fn, "r");
  qassert(not qfile.null());
  std::vector<std::string> lines = qgetlines(qfile);
  qfclose(qfile);
  return lines;
}

// -------------------

inline std::string qcat_info(const std::string& path)
{
  TIMER("qcat_info");
  if (0 == get_id_node()) {
    return qcat(path);
  } else {
    return std::string();
  }
}

inline std::string show_file_crc32(const std::pair<std::string, crc32_t>& fcrc)
{
  return ssprintf("%08X  fn='%s'", fcrc.second, fcrc.first.c_str());
}

inline std::string show_files_crc32(
    const std::vector<std::pair<std::string, crc32_t> >& fcrcs)
{
  std::ostringstream out;
  for (long i = 0; i < (long)fcrcs.size(); ++i) {
    out << ssprintf("%5ld ", i) << show_file_crc32(fcrcs[i]) << std::endl;
  }
  return out.str();
}

inline std::pair<std::string, crc32_t> check_file_crc32(const std::string& fn)
{
  TIMER_VERBOSE("check_file_crc32");
  std::pair<std::string, crc32_t> p;
  p.first = fn;
  p.second = compute_crc32(fn);
  displayln_info(show_file_crc32(p));
  return p;
}

inline int qtouch_info(const std::string& path)
{
  TIMER("qtouch_info");
  if (0 == get_id_node()) {
    return qtouch(path);
  } else {
    return 0;
  }
}

inline int qtouch_info(const std::string& path, const std::string& content)
{
  TIMER("qtouch_info");
  if (0 == get_id_node()) {
    return qtouch(path, content);
  } else {
    return 0;
  }
}

inline int qappend_info(const std::string& path, const std::string& content)
{
  TIMER("qappend_info");
  if (0 == get_id_node()) {
    return qappend(path, content);
  } else {
    return 0;
  }
}

}  // namespace qlat
