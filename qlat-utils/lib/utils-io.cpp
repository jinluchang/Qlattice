#include <qlat-utils/qutils-io.h>

namespace qlat
{  //

std::vector<std::string> qls_aux(const std::string& path,
                                 const bool is_sort = true)
{
  std::vector<std::string> contents;
  DIR* dir = opendir(path.c_str());
  if (dir == NULL) {
    return contents;
  }
  struct dirent* d;
  while ((d = readdir(dir)) != NULL) {
    if (!strcmp(d->d_name, ".") || !strcmp(d->d_name, "..")) {
      continue;
    }
    contents.push_back(path + "/" + d->d_name);
  }
  closedir(dir);
  if (is_sort) {
    std::sort(contents.begin(), contents.end());
  }
  return contents;
}

std::vector<std::string> qls_all_aux(const std::string& path,
                                     const bool is_folder_before_files = false,
                                     const bool is_sort = true)
// list all files and folder in path (not including it self)
{
  std::vector<std::string> all_contents;
  if (not is_directory(path)) {
    return all_contents;
  }
  const std::vector<std::string> contents = qls_aux(path, is_sort);
  for (long i = 0; i < (long)contents.size(); ++i) {
    const std::string& path_i = contents[i];
    if (not is_directory(path_i)) {
      all_contents.push_back(path_i);
    } else {
      if (is_folder_before_files) {
        all_contents.push_back(path_i);
        vector_append(all_contents,
                      qls_all_aux(path_i, is_folder_before_files, is_sort));
      } else {
        // default behavior
        vector_append(all_contents,
                      qls_all_aux(path_i, is_folder_before_files, is_sort));
        all_contents.push_back(path_i);
      }
    }
  }
  return all_contents;
}

int qremove_all_aux(const std::string& path)
{
  if (not is_directory(path)) {
    return qremove(path);
  } else {
    int ret = 0;
    const std::vector<std::string> paths = qls_aux(path);
    for (long i = 0; i < (long)paths.size(); ++i) {
      ret += qremove_all_aux(paths[i]);
    }
    return ret + qremove(path);
  }
}

// --------------------------------

void flush() { fflush(get_output_file()); }

std::string dirname(const std::string& fn)
// try to follow libgen.h version see man 3 dirname
{
  long cur = fn.size() - 1;
  // remove trailing '/'
  while (cur > 0 and fn[cur] == '/') {
    cur -= 1;
  }
  if (cur < 0) {
    return ".";
  } else if (cur == 0) {
    if (fn[cur] == '/') {
      return "/";
    } else {
      return ".";
    }
  } else {
    // remove last component
    while (cur >= 0 and fn[cur] != '/') {
      cur -= 1;
    }
    if (cur < 0) {
      return ".";
    } else {
      // remove trailing '/'
      while (cur > 0 and fn[cur] == '/') {
        cur -= 1;
      }
      return std::string(fn, 0, cur + 1);
    }
  }
  qassert(false);
  return std::string();
}

std::string basename(const std::string& fn)
// try to follow libgen.h version see man 3 basename
{
  long cur = fn.size() - 1;
  // remove trailing '/'
  while (cur > 0 and fn[cur] == '/') {
    cur -= 1;
  }
  if (cur < 0) {
    return "";
  } else if (cur == 0) {
    if (fn[cur] == '/') {
      return "/";
    } else {
      return std::string(fn, 0, cur + 1);
    }
  } else {
    const long pos_stop = cur + 1;
    // skip last component
    while (cur >= 0 and fn[cur] != '/') {
      cur -= 1;
    }
    return std::string(fn, cur + 1, pos_stop);
  }
  qassert(false);
  return std::string();
}

std::string remove_trailing_slashes(const std::string& fn)
// remove trailing slashes (but won't remove the first slash)
// e.g.
// remove_trailing_slashes("/home/") = "/home"
// remove_trailing_slashes("//") = "/"
{
  long cur = fn.size() - 1;
  while (cur > 0 and fn[cur] == '/') {
    cur -= 1;
  }
  return std::string(fn, 0, cur + 1);
}

int qrename(const std::string& old_path, const std::string& new_path)
{
  TIMER("qrename");
  displayln(0,
            ssprintf("qrename: '%s' '%s'", old_path.c_str(), new_path.c_str()));
  return rename(old_path.c_str(), new_path.c_str());
}

int qrename_info(const std::string& old_path, const std::string& new_path)
{
  TIMER("qrename_info");
  if (0 == get_id_node()) {
    return qrename(old_path, new_path);
  } else {
    return 0;
  }
}

std::vector<std::string> qls(const std::string& path, const bool is_sort)
{
  return qls_aux(remove_trailing_slashes(path), is_sort);
}

std::vector<std::string> qls_all(const std::string& path,
                                 const bool is_folder_before_files,
                                 const bool is_sort)
// list files before its folder
{
  return qls_all_aux(remove_trailing_slashes(path), is_folder_before_files,
                     is_sort);
}

bool does_file_exist(const std::string& fn)
{
  TIMER("does_file_exist")
  struct stat sb;
  return 0 == stat(fn.c_str(), &sb);
}

bool is_directory(const std::string& fn)
{
  struct stat sb;
  if (0 != stat(fn.c_str(), &sb)) {
    return false;
  }
  return S_ISDIR(sb.st_mode);
}

bool is_regular_file(const std::string& fn)
{
  struct stat sb;
  if (0 != stat(fn.c_str(), &sb)) {
    return false;
  }
  return S_ISREG(sb.st_mode);
}

bool is_directory_cache(const std::string& dir_)
// Check the existence of directory from the root directory of dir
// Cache all intermediate results
{
  TIMER("is_directory_cache");
  const std::string dir = remove_trailing_slashes(dir_) + "/";
  Cache<std::string, bool>& cache = get_is_directory_cache();
  for (long cur = 0; cur < (long)dir.size(); ++cur) {
    if (dir[cur] != '/') {
      continue;
    }
    const std::string subdir = dir.substr(0, cur + 1);
    bool b;
    if (cache.has(subdir)) {
      b = cache[subdir];
    } else {
      b = is_directory(subdir);
      cache[subdir] = b;
    }
    if (not b) {
      return false;
    }
  }
  return true;
}

bool is_regular_file_cache(const std::string& fn)
{
  TIMER("is_regular_file_cache");
  const std::string dir = dirname(fn);
  if (not is_directory_cache(dir)) {
    return false;
  }
  return is_regular_file(fn);
}

bool does_file_exist_cache(const std::string& fn)
{
  TIMER("does_file_exist_cache");
  if (fn.size() > 0 and fn.back() == '/') {
    return is_directory_cache(fn);
  }
  const std::string dir = dirname(fn);
  if (not is_directory_cache(dir)) {
    return false;
  }
  return does_file_exist(fn);
}

int qremove(const std::string& path)
{
  displayln(0, ssprintf("qremove: '%s'", path.c_str()));
  return std::remove(path.c_str());
}

int qremove_all(const std::string& path)
{
  clear_is_directory_cache();
  return qremove_all_aux(remove_trailing_slashes(path));
}

int qmkdir(const std::string& path, const mode_t mode)
{
  TIMER("qmkdir");
  clear_is_directory_cache();
  mkdir(path.c_str(), mode);
  return check_dir(path, mode);
}

int qmkdir_p(const std::string& path_, const mode_t mode)
// return 0 if successful
{
  TIMER("qmkdir_p");
  clear_is_directory_cache();
  std::string path = remove_trailing_slashes(path_);
  if (is_directory(path)) {
    return 0;
  }
  std::vector<std::string> paths;
  while (true) {
    if (0 == mkdir(path.c_str(), mode)) {
      break;
    } else {
      paths.push_back(path);
      path = dirname(path);
      if (does_file_exist(path)) {
        // qwarn(fname + ssprintf(": '%s' failed.", path_.c_str()));
        break;
      }
    }
  }
  for (long i = paths.size() - 1; i >= 0; i -= 1) {
    if (not(0 == mkdir(paths[i].c_str(), mode))) {
      // qwarn(fname + ssprintf(": '%s' failed.", path_.c_str()));
      continue;
    }
  }
  if (not is_directory(path)) {
    qwarn(fname + ssprintf(": '%s' failed.", path_.c_str()));
    return 1;
  }
  return 0;
}

int qmkdir_info(const std::string& path, const mode_t mode)
{
  TIMER("qmkdir_info");
  clear_is_directory_cache();
  if (0 == get_id_node()) {
    return qmkdir(path, mode);
  } else {
    return 0;
  }
}

int qmkdir_p_info(const std::string& path, const mode_t mode)
{
  TIMER("qmkdir_info");
  clear_is_directory_cache();
  if (0 == get_id_node()) {
    return qmkdir_p(path, mode);
  } else {
    return 0;
  }
}

}  // namespace qlat
