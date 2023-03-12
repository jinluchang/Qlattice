#include <qlat-utils/qutils-io.h>

namespace qlat
{  //

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

int qremove(const std::string& path)
{
  displayln(0, ssprintf("qremove: '%s'", path.c_str()));
  return std::remove(path.c_str());
}

int qremove_all(const std::string& path)
{
  return qremove_all_aux(remove_trailing_slashes(path));
}

int qmkdir(const std::string& path, const mode_t mode)
{
  TIMER("qmkdir");
  mkdir(path.c_str(), mode);
  return check_dir(path, mode);
}

int qmkdir_p(const std::string& path_, const mode_t mode)
// return 0 if successful
{
  TIMER("qmkdir_p");
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
  if (0 == get_id_node()) {
    return qmkdir(path, mode);
  } else {
    return 0;
  }
}

int qmkdir_p_info(const std::string& path, const mode_t mode)
{
  TIMER("qmkdir_info");
  if (0 == get_id_node()) {
    return qmkdir_p(path, mode);
  } else {
    return 0;
  }
}

}  // namespace qlat
