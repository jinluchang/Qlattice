#include <qlat-utils/qar.h>
#include <qlat-utils/utils-io.h>

namespace qlat
{  //

static std::vector<std::string> qls_aux(const std::string& path,
                                        const bool is_sort = true)
{
  std::vector<std::string> contents;
  if (not is_directory(path)) {
    return contents;
  }
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

static std::vector<std::string> qls_all_aux(
    const std::string& path, const bool is_folder_before_files = false,
    const bool is_sort = true)
// list all files and folder in path (not including it self)
{
  std::vector<std::string> all_contents;
  if (not is_directory(path)) {
    return all_contents;
  }
  const std::vector<std::string> contents = qls_aux(path, is_sort);
  for (Long i = 0; i < (Long)contents.size(); ++i) {
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

static Int qremove_aux(const std::string& path)
{
  TIMER("qremove_aux")
  displayln_info(0, ssprintf("qremove_aux: '%s' (id_node=%d).", path.c_str(),
                             get_id_node()));
  if (does_file_exist(path)) {
    const Int ret = std::remove(path.c_str());
    if (ret != 0) {
      qwarn(fname + ssprintf(": '%s' failed ret='%d'", path.c_str(), ret));
    }
    return ret;
  } else {
    displayln_info(
        0, ssprintf("qremove_aux: '%s' file does not exist. (id_node=%d).",
                    path.c_str(), get_id_node()));
    return 0;
  }
}

static Int qremove_all_aux(const std::string& path)
{
  if (not is_directory(path)) {
    return qremove_aux(path);
  } else {
    Int ret = 0;
    const std::vector<std::string> paths = qls_aux(path);
    for (Long i = 0; i < (Long)paths.size(); ++i) {
      ret += qremove_all_aux(paths[i]);
    }
    return ret + qremove_aux(path);
  }
}

// --------------------------------

std::string basename(const std::string& fn)
// try to follow libgen.h version see man 3 basename
{
  Long cur = fn.size() - 1;
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
    const Long pos_stop = cur + 1;
    // skip last component
    while (cur >= 0 and fn[cur] != '/') {
      cur -= 1;
    }
    return std::string(fn, cur + 1, pos_stop);
  }
  Qassert(false);
  return std::string();
}

std::string dirname(const std::string& fn)
// try to follow libgen.h version see man 3 dirname
{
  Long cur = fn.size() - 1;
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
  Qassert(false);
  return std::string();
}

std::vector<std::string> all_dirname_vec(const std::string& fn)
// use dirname
// return [ dirname(fn), dirname(dirname(fn)), ..., ] until results does not
// change.
{
  std::vector<std::string> ret;
  std::string dn = dirname(fn);
  ret.push_back(dn);
  while (true) {
    dn = dirname(dn);
    if (dn == ret.back()) {
      break;
    }
    ret.push_back(dn);
  }
  return ret;
}

std::string remove_trailing_slashes(const std::string& fn)
// remove trailing slashes (but won't remove the first slash)
// e.g.
// remove_trailing_slashes("/home/") = "/home"
// remove_trailing_slashes("//") = "/"
{
  Long cur = fn.size() - 1;
  while (cur > 0 and fn[cur] == '/') {
    cur -= 1;
  }
  return std::string(fn, 0, cur + 1);
}

// --------------------------------

bool does_file_exist(const std::string& fn)
{
  TIMER("does_file_exist")
  struct stat sb;
  return 0 == stat(fn.c_str(), &sb);
}

bool is_directory(const std::string& fn)
{
  TIMER("is_directory")
  displayln_info(
      0, fname + ssprintf(": '%s' (id_node=%d).", fn.c_str(), get_id_node()));
  struct stat sb;
  if (0 != stat(fn.c_str(), &sb)) {
    displayln_info(0,
                   fname + ssprintf(": '%s' file does not exist. (id_node=%d).",
                                    fn.c_str(), get_id_node()));
    return false;
  }
  const bool ret = S_ISDIR(sb.st_mode);
  displayln_info(0, fname + ssprintf(": '%s' file exists. ret=%d (id_node=%d).",
                                     fn.c_str(), ret, get_id_node()));
  return ret;
}

bool is_regular_file(const std::string& fn)
{
  TIMER("is_regular_file")
  struct stat sb;
  if (0 != stat(fn.c_str(), &sb)) {
    return false;
  }
  return S_ISREG(sb.st_mode);
}

void remove_entry_directory_cache(const std::string& dir_)
// if operation may be done related to dir_ remove relevant entries.
// (1) remove entries start with dir.
// (2) remove entries that are parents of dir but value being false.
{
  TIMER("remove_entry_directory_cache");
  const std::string dir = remove_trailing_slashes(dir_) + "/";
  Cache<std::string, bool>& cache = get_is_directory_cache();
  std::vector<std::string> key_vec;
  // find entries start with dir
  auto& m = cache.m;
  for (auto it = m.begin(); it != m.end(); ++it) {
    const std::string& key = it->first;
    if (key.substr(0, dir.size()) == dir) {
      key_vec.push_back(key);
    }
  }
  // find entries that are parents of dir but value being false
  for (Long cur = 0; cur < (Long)dir.size(); ++cur) {
    if (dir[cur] != '/') {
      continue;
    }
    const std::string subdir = dir.substr(0, cur + 1);
    if (cache.has(subdir)) {
      if (cache[subdir] == false) {
        key_vec.push_back(subdir);
      }
    }
  }
  // remove entries
  for (Long i = 0; i < (Long)key_vec.size(); ++i) {
    cache.erase(key_vec[i]);
  }
}

void add_entry_directory_cache(const std::string& dir_, bool is_directory)
// if is_directory: mark directory and all subdirectories exist.
{
  TIMER("add_entry_directory_cache");
  const std::string dir = remove_trailing_slashes(dir_) + "/";
  Cache<std::string, bool>& cache = get_is_directory_cache();
  if (is_directory) {
    for (Long cur = 0; cur < (Long)dir.size(); ++cur) {
      if (dir[cur] != '/') {
        continue;
      }
      const std::string subdir = dir.substr(0, cur + 1);
      cache[subdir] = true;
    }
  } else {
    remove_entry_directory_cache(dir);
    cache[dir] = false;
  }
}

bool is_directory_cache(const std::string& dir_)
// Check the existence of directory from the root directory of dir.
// If some intermediate step is false, return false.
// Cache all intermediate results.
{
  TIMER("is_directory_cache");
  const std::string dir = remove_trailing_slashes(dir_) + "/";
  Cache<std::string, bool>& cache = get_is_directory_cache();
  for (Long cur = 0; cur < (Long)dir.size(); ++cur) {
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

// --------------------------------

void flush() { fflush(stdout); }

Int qtruncate(const std::string& path, const Long offset)
// return true if successful.
{
  const Int ret = truncate(path.c_str(), offset);
  return ret;
}

std::vector<std::string> qls(const std::string& path, const bool is_sort)
{
  TIMER("qls");
  return qls_aux(remove_trailing_slashes(path), is_sort);
}

std::vector<std::string> qls_all(const std::string& path,
                                 const bool is_folder_before_files,
                                 const bool is_sort)
// list files before its folder
{
  TIMER("qls_all");
  return qls_all_aux(remove_trailing_slashes(path), is_folder_before_files,
                     is_sort);
}

// --------------------------------

Int qrename(const std::string& old_path, const std::string& new_path)
{
  TIMER_VERBOSE("qrename");
  displayln_info(0,
                 fname + ssprintf(": '%s' '%s' (id_node=%d).", old_path.c_str(),
                                  new_path.c_str(), get_id_node()));
  remove_entry_directory_cache(old_path);
  remove_entry_directory_cache(new_path);
  update_qar_cache_due_to_change_of_directory(old_path);
  update_qar_cache_due_to_change_of_qar_file(old_path);
  update_qar_cache_due_to_change_of_qar_file(new_path);
  const Int ret = rename(old_path.c_str(), new_path.c_str());
  if (ret != 0) {
    qwarn(fname + ssprintf(": '%s' '%s' ret=%d", old_path.c_str(),
                           new_path.c_str(), ret));
  }
  return ret;
}

Int qremove(const std::string& path)
{
  TIMER_VERBOSE("qremove")
  remove_entry_directory_cache(path);
  update_qar_cache_due_to_change_of_directory(path);
  update_qar_cache_due_to_change_of_qar_file(path);
  const Int ret = qremove_aux(path);
  return ret;
}

Int qremove_all(const std::string& path)
{
  TIMER_VERBOSE("qremove_all")
  remove_entry_directory_cache(path);
  update_qar_cache_due_to_change_of_directory(path);
  update_qar_cache_due_to_change_of_qar_file(path);
  const Int ret = qremove_all_aux(remove_trailing_slashes(path));
  return ret;
}

Int qmkdir(const std::string& path, const mode_t mode)
{
  TIMER_VERBOSE("qmkdir");
  displayln_info(
      0, fname + ssprintf(": '%s' (id_node=%d).", path.c_str(), get_id_node()));
  remove_entry_directory_cache(path);
  const Int ret = mkdir(path.c_str(), mode);
  if (ret != 0) {
    qwarn(fname + ssprintf(": qmkdir failed '%s' ret=%d", path.c_str(), ret));
  }
  return ret;
}

Int qmkdir_p(const std::string& dir_, const mode_t mode)
// return 0 if successful
{
  TIMER_VERBOSE("qmkdir_p");
  const std::string dir = remove_trailing_slashes(dir_) + "/";
  if (is_directory_cache(dir)) {
    return 0;
  }
  for (Long cur = 0; cur < (Long)dir.size(); ++cur) {
    if (dir[cur] != '/') {
      continue;
    }
    const std::string subdir = dir.substr(0, cur + 1);
    if (not is_directory_cache(subdir)) {
      const Int ret = qmkdir(subdir, mode);
      if (ret != 0) {
        qwarn(fname +
              ssprintf(": '%s' failed at '%s'.", dir_.c_str(), subdir.c_str()));
        return ret;
      }
    }
  }
  if (not is_directory_cache(dir)) {
    qwarn(fname + ssprintf(": '%s' failed.", dir_.c_str()));
    return 1;
  }
  return 0;
}

bool check_dir(const std::string& path, const mode_t mode)
// try to create if does not exist.
{
  TIMER_VERBOSE("check_dir");
  if (is_directory(path)) {
    add_entry_directory_cache(path, true);
    return true;
  }
  const Int max_attempts = 8;
  for (Int i = 0; i < max_attempts; ++i) {
    qmkdir(path, mode);
    ssleep(0.1);
    if (is_directory(path)) {
      add_entry_directory_cache(path, true);
      return true;
    }
  }
  return false;
}

// --------------------------------

Int qrename_info(const std::string& old_path, const std::string& new_path)
{
  TIMER("qrename_info");
  if (0 == get_id_node()) {
    return qrename(old_path, new_path);
  } else {
    remove_entry_directory_cache(new_path);
    remove_entry_directory_cache(old_path);
    update_qar_cache_due_to_change_of_directory(old_path);
    update_qar_cache_due_to_change_of_qar_file(old_path);
    update_qar_cache_due_to_change_of_qar_file(new_path);
    return 0;
  }
}

Int qremove_info(const std::string& path)
{
  TIMER("qremove_info");
  if (0 == get_id_node()) {
    return qremove(path);
  } else {
    remove_entry_directory_cache(path);
    update_qar_cache_due_to_change_of_directory(path);
    update_qar_cache_due_to_change_of_qar_file(path);
    return 0;
  }
}

Int qremove_all_info(const std::string& path)
{
  TIMER("qremove_all_info");
  if (0 == get_id_node()) {
    return qremove_all(path);
  } else {
    remove_entry_directory_cache(path);
    update_qar_cache_due_to_change_of_directory(path);
    update_qar_cache_due_to_change_of_qar_file(path);
    return 0;
  }
}

Int qmkdir_info(const std::string& path, const mode_t mode)
{
  TIMER("qmkdir_info");
  if (0 == get_id_node()) {
    return qmkdir(path, mode);
  } else {
    remove_entry_directory_cache(path);
    return 0;
  }
}

Int qmkdir_p_info(const std::string& path, const mode_t mode)
{
  TIMER("qmkdir_info");
  if (0 == get_id_node()) {
    return qmkdir_p(path, mode);
  } else {
    remove_entry_directory_cache(path);
    return 0;
  }
}

// --------------------------------

std::vector<std::string> qls_sync_node(const std::string& path,
                                       const bool is_sort)
{
  TIMER("qls_sync_node");
  std::vector<std::string> ret;
  if (0 == get_id_node()) {
    ret = qls(path, is_sort);
  }
  Int bret = bcast_val(ret, 0);
  Qassert(bret == 0);
  return ret;
}

std::vector<std::string> qls_all_sync_node(const std::string& path,
                                           const bool is_folder_before_files,
                                           const bool is_sort)
{
  TIMER("qls_all_sync_node");
  std::vector<std::string> ret;
  if (0 == get_id_node()) {
    ret = qls_all(path, is_folder_before_files, is_sort);
  }
  Int bret = bcast_val(ret, 0);
  Qassert(bret == 0);
  return ret;
}

bool does_file_exist_sync_node(const std::string& fn)
{
  Long nfile = 0;
  if (0 == get_id_node()) {
    if (does_file_exist(fn)) {
      nfile = 1;
    }
  }
  glb_sum_val(nfile);
  return 0 != nfile;
}

bool is_directory_sync_node(const std::string& fn)
{
  Long nfile = 0;
  if (0 == get_id_node()) {
    if (is_directory(fn)) {
      nfile = 1;
    }
  }
  glb_sum_val(nfile);
  return 0 != nfile;
}

bool is_regular_file_sync_node(const std::string& fn)
{
  Long nfile = 0;
  if (0 == get_id_node()) {
    if (is_regular_file(fn)) {
      nfile = 1;
    }
  }
  glb_sum_val(nfile);
  return 0 != nfile;
}

bool does_file_exist_cache_sync_node(const std::string& fn)
{
  Long nfile = 0;
  if (0 == get_id_node()) {
    if (does_file_exist_cache(fn)) {
      nfile = 1;
    }
  }
  glb_sum_val(nfile);
  return 0 != nfile;
}

bool is_directory_cache_sync_node(const std::string& fn)
{
  Long nfile = 0;
  if (0 == get_id_node()) {
    if (is_directory_cache(fn)) {
      nfile = 1;
    }
  }
  glb_sum_val(nfile);
  return 0 != nfile;
}

bool is_regular_file_cache_sync_node(const std::string& fn)
{
  Long nfile = 0;
  if (0 == get_id_node()) {
    if (is_regular_file_cache(fn)) {
      nfile = 1;
    }
  }
  glb_sum_val(nfile);
  return 0 != nfile;
}

Int qmkdir_sync_node(const std::string& path, const mode_t mode)
{
  TIMER("qmkdir_sync_node");
  Long ret = 0;
  if (0 == get_id_node()) {
    ret = qmkdir(path, mode);
  } else {
    remove_entry_directory_cache(path);
  }
  glb_sum_val(ret);
  return ret;
}

Int qmkdir_p_sync_node(const std::string& path, const mode_t mode)
{
  TIMER("qmkdir_p_sync_node");
  Long ret = 0;
  if (0 == get_id_node()) {
    ret = qmkdir_p(path, mode);
  } else {
    remove_entry_directory_cache(path);
  }
  glb_sum_val(ret);
  return ret;
}

Int qremove_sync_node(const std::string& path)
{
  TIMER("qremove_sync_node");
  Long ret = 0;
  if (0 == get_id_node()) {
    ret = qremove(path);
  } else {
    remove_entry_directory_cache(path);
    update_qar_cache_due_to_change_of_directory(path);
    update_qar_cache_due_to_change_of_qar_file(path);
  }
  glb_sum_val(ret);
  return ret;
}

Int qremove_all_sync_node(const std::string& path)
{
  TIMER("qremove_all_sync_node");
  Long ret = 0;
  if (0 == get_id_node()) {
    ret = qremove_all(path);
  } else {
    remove_entry_directory_cache(path);
    update_qar_cache_due_to_change_of_directory(path);
    update_qar_cache_due_to_change_of_qar_file(path);
  }
  glb_sum_val(ret);
  return ret;
}

}  // namespace qlat
