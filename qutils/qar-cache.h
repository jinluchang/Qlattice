#pragma once

#include <qutils/qar.h>
#include <qutils/cache.h>

namespace qlat
{  //

inline Cache<std::string, QarFile>& get_qar_read_cache()
// key should be the path prefix of the contents of the qar file.
// Note: key should end with '/'.
{
  static Cache<std::string, QarFile> cache("QarReadCache", 8, 1);
  return cache;
}

inline std::string get_qar_read_cache_key(const std::string& path)
// return key of get_qar_read_cache() that may contain path
// return empty string if no cached key is found.
// Note: key should end with '/'.
{
  TIMER("get_qar_read_cache_key");
  Cache<std::string, QarFile>& cache = get_qar_read_cache();
  for (auto it = cache.m.cbegin(); it != cache.m.cend(); ++it) {
    const std::string& key = it->first;
    if (key == path.substr(0, key.size())) {
      return key;
    }
  }
  std::string path_dir = dirname(path);
  while (true) {
    if (path_dir == "/" or path_dir == ".") {
      return std::string();
    }
    if (does_file_exist(path_dir + ".qar")) {
      const std::string key = path_dir + "/";
      qassert(not cache.has(key));
      QarFile& qar = cache[key];
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
  QarFile& qar = get_qar_read_cache()[key];
  return has(qar, fn);
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
      QarFile& qar = get_qar_read_cache()[key];
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

}  // namespace qlat
