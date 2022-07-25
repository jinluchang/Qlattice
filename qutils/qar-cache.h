#pragma once

#include <qutils/qar.h>
#include <qutils/cache.h>

namespace qlat
{  //

inline Cache<std::string, QarFileMultiVol>& get_qar_read_cache()
// key should be the path prefix of the contents of the qar file.
// Note: key should end with '/'.
{
  static Cache<std::string, QarFileMultiVol> cache("QarReadCache", 16, 1);
  return cache;
}

inline std::string get_qar_read_cache_key(const std::string& path)
// return key of get_qar_read_cache() that may contain path
// return empty string if no cached key is found.
// Note: key should end with '/'.
// Steps:
// (1) Search in Cache. Return if found matching key.
// (2) If not found, check if path exists. If exists, return empty key.
// (3) If does not exist, try to find qar file yet to be in cache.
{
  TIMER("get_qar_read_cache_key");
  Cache<std::string, QarFileMultiVol>& cache = get_qar_read_cache();
  for (auto it = cache.m.cbegin(); it != cache.m.cend(); ++it) {
    const std::string& key = it->first;
    if (key == path.substr(0, key.size())) {
      return key;
    }
  }
  if (does_file_exist(path)) {
    return std::string();
  }
  std::string path_dir = remove_trailing_slashes(path);
  while (true) {
    if (path_dir == "/" or path_dir == ".") {
      return std::string();
    }
    if (does_file_exist(path_dir + ".qar")) {
      const std::string key = path_dir + "/";
      qassert(not cache.has(key));
      QarFileMultiVol& qar = cache[key];
      qar.init(path_dir + ".qar", "r");
      qassert(not qar.null());
      return key;
    }
    path_dir = dirname(path_dir);
  }
  qassert(false);
  return std::string();
}

inline bool does_file_exist_qar(const std::string& path)
// interface function
// Note: should only check file, not directory.
{
  TIMER("does_file_exist_qar");
  const std::string key = get_qar_read_cache_key(path);
  if (key == "") {
    return does_file_exist(path);
  }
  qassert(key == path.substr(0, key.size()));
  const std::string fn = path.substr(key.size());
  QarFileMultiVol& qar = get_qar_read_cache()[key];
  return has(qar, fn);
}

inline bool does_file_or_directory_exist_qar(const std::string& path)
// interface function
{
  TIMER("does_file_or_directory_exist_qar");
  if (does_file_exist_qar(path)) {
    return true;
  }
  const std::string path_test = remove_trailing_slashes(path) + "/TEST-FILE";
  const std::string key = get_qar_read_cache_key(path_test);
  qassert(key == path_test.substr(0, key.size()));
  const std::string fn =
      (remove_trailing_slashes(path) + "/").substr(key.size());
  if (fn == "") {
    return true;
  }
  QarFileMultiVol& qar = get_qar_read_cache()[key];
  const std::vector<std::string> fn_list = list(qar);
  for (long i = 0; i < (long)fn_list.size(); ++i) {
    const std::string& fni = fn_list[i];
    if (fni.compare(0, fn.size(), fn) == 0) {
      return true;
    }
  }
  return false;
}

inline void qopen(QFile& qfile, const std::string& path,
                  const std::string& mode)
// interface function
// qfile.null() == true if qopen failed.
// Will open files in qar for read
// Will create directories needed for write / append
{
  TIMER("qopen(path,mode,qfile)");
  if (mode == "r") {
    const std::string key = get_qar_read_cache_key(path);
    if (key == "") {
      qfile.init(path, mode);
      return;
    } else {
      qassert(key == path.substr(0, key.size()));
      const std::string fn = path.substr(key.size());
      QarFileMultiVol& qar = get_qar_read_cache()[key];
      read(qar, fn, qfile);
    }
  } else if (mode == "w" or mode == "a") {
    const std::string path_dir = dirname(path);
    qmkdir_p(path_dir);
    qfile.init(path, mode);
  } else {
    qassert(false);
  }
}

// -------------------

inline crc32_t compute_crc32(QFile& qfile)
// interface function
// compute_crc32 for all the remaining data.
{
  TIMER_VERBOSE_FLOPS("compute_crc32");
  qassert(not qfile.null());
  const size_t chunk_size = 16 * 1024 * 1024;
  std::vector<char> data(chunk_size);
  crc32_t crc = 0;
  while (true) {
    const long size = qread_data(get_data(data), qfile);
    timer.flops += size;
    if (size == 0) {
      break;
    }
    crc = crc32_par(crc, Vector<char>(data.data(), size));
  }
  return crc;
}

inline crc32_t compute_crc32(const std::string& path)
{
  QFile qfile;
  qopen(qfile, path, "r");
  return compute_crc32(qfile);
}

inline std::vector<std::string> qgetlines(const std::string& fn)
{
  QFile qfile;
  qopen(qfile, fn, "r");
  qassert(not qfile.null());
  std::vector<std::string> lines = qgetlines(qfile);
  return lines;
}

inline std::string qcat(const std::string& path)
{
  TIMER("qcat");
  QFile qfile;
  qopen(qfile, path, "r");
  if (qfile.null()) {
    return "";
  }
  qfseek(qfile, 0, SEEK_END);
  const long length = qftell(qfile);
  qfseek(qfile, 0, SEEK_SET);
  std::string ret(length, 0);
  const long length_actual = qfread(&ret[0], 1, length, qfile);
  qassert(length == length_actual);
  return ret;
}

inline int qtouch(const std::string& path)
// return 0 if success
{
  TIMER("qtouch");
  QFile qfile;
  qopen(qfile, path, "w");
  if (qfile.null()) {
    return 1;
  }
  return 0;
}

inline int qtouch(const std::string& path, const std::string& content)
{
  TIMER("qtouch");
  QFile qfile;
  qopen(qfile, path + ".partial", "w");
  if (qfile.null()) {
    return 1;
  }
  const long total_bytes = qwrite_data(content, qfile);
  qassert(total_bytes == content.size());
  qfile.close();
  return qrename(path + ".partial", path);
}

inline int qappend(const std::string& path, const std::string& content)
{
  TIMER("qappend");
  QFile qfile;
  qopen(qfile, path, "a");
  if (qfile.null()) {
    return 1;
  }
  const long total_bytes = qwrite_data(content, qfile);
  qassert(total_bytes == content.size());
  qfile.close();
  return 0;
}

}  // namespace qlat
