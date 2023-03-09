#include <qlat-utils/qar-cache.h>

namespace qlat
{  //

void register_file(const QarFileVol& qar, const std::string& fn,
                   const QarSegmentInfo& qsinfo)
{
  if (not has(qar.p->qsinfo_map, fn)) {
    qar.p->fn_list.push_back(fn);
    qar.p->qsinfo_map[fn] = qsinfo;
    std::string dir = dirname(fn);
    while (dir != ".") {
      if (has(qar.p->directories, dir)) {
        break;
      } else {
        qar.p->directories.insert(dir);
        dir = dirname(dir);
      }
    }
    if (qar.p->max_offset < qsinfo.offset_end) {
      qar.p->max_offset = qsinfo.offset_end;
    }
  } else {
    qassert(qar.p->qsinfo_map[fn].offset == qsinfo.offset);
  }
}

bool does_regular_file_exist_qar(const std::string& path)
// interface function
// Note: should only check file, not directory.
{
  TIMER("does_regular_file_exist_qar");
  const std::string key = get_qar_read_cache_key(path);
  if (key == path) {
    return true;
  } else if (key == "") {
    return false;
  }
  qassert(key == path.substr(0, key.size()));
  const std::string fn = path.substr(key.size());
  QarFile& qar = get_qar_read_cache()[key];
  return has_regular_file(qar, fn);
}

bool does_file_exist_qar(const std::string& path)
// interface function
{
  TIMER("does_file_exist_qar");
  const std::string key = get_qar_read_cache_key(path);
  if (key == path) {
    return true;
  } else if (key == "") {
    return false;
  }
  qassert(key == path.substr(0, key.size()));
  const std::string fn = path.substr(key.size());
  QarFile& qar = get_qar_read_cache()[key];
  return has(qar, fn);
}

QFile qfopen(const std::string& path, const std::string& mode)
// interface function
// qfile.null() == true if qopen failed.
// Will open files in qar for read
// Will create directories needed for write / append
{
  TIMER("qfopen(path,mode)");
  if (mode == "r") {
    const std::string key = get_qar_read_cache_key(path);
    if (key == "") {
      return QFile();
    } else if (key == path) {
      return QFile(path, mode);
    } else {
      qassert(key == path.substr(0, key.size()));
      const std::string fn = path.substr(key.size());
      QarFile& qar = get_qar_read_cache()[key];
      QFile qfile;
      read(qar, fn, qfile);
      return qfile;
    }
  } else if (mode == "w" or mode == "a") {
    const std::string path_dir = dirname(path);
    qmkdir_p(path_dir);
    return QFile(path, mode);
  } else {
    qassert(false);
  }
  return QFile();
}

std::string qgetline(const QFile& qfile)
// interface function
// read an entire line including the final '\n' char.
{
  qassert(not qfile.null());
  const int code = qfseek(qfile, qfile.p->pos, SEEK_SET);
  qassert(code == 0);
  char* lineptr = NULL;
  size_t n = 0;
  const long size = getline(&lineptr, &n, qfile.get_fp());
  qfile.p->is_eof = feof(qfile.get_fp()) != 0;
  if (size > 0) {
    std::string ret;
    const long pos = ftell(qfile.get_fp()) - qfile.p->offset_start;
    if (pos < 0) {
      qerr(
          ssprintf("qgetline: '%s' initial_pos='%ld' final_pos='%ld' "
                   "getline_size='%ld'.",
                   qfile.path().c_str(), qfile.p->pos, pos, size));
    }
    qassert(pos >= 0);
    if (qfile.p->offset_end != -1 and
        qfile.p->offset_start + pos > qfile.p->offset_end) {
      qfseek(qfile, 0, SEEK_END);
      qfile.p->is_eof = true;
      const long size_truncate =
          size - (qfile.p->offset_start + pos - qfile.p->offset_end);
      qassert(size_truncate >= 0);
      ret = std::string(lineptr, size_truncate);
    } else {
      qfile.p->pos = pos;
      ret = std::string(lineptr, size);
    }
    std::free(lineptr);
    return ret;
  } else {
    std::free(lineptr);
    return std::string();
  }
}

std::vector<std::string> list(const QarFile& qar)
// interface function
{
  if (qar.null()) {
    return std::vector<std::string>();
  }
  std::vector<std::string> fn_list;
  for (unsigned long i = 0; i < qar.size(); ++i) {
    const QarFileVol& qar_v = qar[i];
    qassert(not qar_v.null());
    qassert(qar_v.mode() == "r");
    vector_append(fn_list, list(qar_v));
  }
  return fn_list;
}

void save_qar_index(const QarFile& qar, const std::string& fn)
// interface function
{
  if (qar.null()) {
    return;
  }
  std::vector<std::string> lines;
  for (long i = 0; i < (long)qar.size(); ++i) {
    const QarFileVol& qar_v = qar[i];
    qassert(not qar_v.null());
    qassert(qar_v.mode() == "r");
    read_through(qar_v);
    const std::vector<std::string>& fn_list = qar_v.p->fn_list;
    const std::map<std::string, QarSegmentInfo> qsinfo_map =
        qar_v.p->qsinfo_map;
    for (long j = 0; j < (long)fn_list.size(); ++j) {
      const std::string& fn = fn_list[j];
      const QarSegmentInfo& qsinfo = qsinfo_map.at(fn);
      qassert(qsinfo.fn_len == (long)fn.size());
      const std::string line1 =
          ssprintf("QAR-FILE %ld %ld %ld\n", i, j, fn.size());
      const std::string line2 = fn + "\n";
      const std::string line3 = ssprintf(
          "%ld %ld %ld %ld %ld %ld %ld %ld\n", qsinfo.offset, qsinfo.offset_fn,
          qsinfo.offset_info, qsinfo.offset_data, qsinfo.offset_end,
          qsinfo.fn_len, qsinfo.info_len, qsinfo.data_len);
      lines.push_back(line1);
      lines.push_back(line2);
      lines.push_back(line3);
    }
  }
  qtouch(fn, lines);
}

void save_qar_index_info(const QarFile& qar, const std::string& fn)
{
  if (0 == get_id_node()) {
    save_qar_index(qar, fn);
  }
}

void parse_qar_index(const QarFile& qar, const std::string& qar_index_content)
{
}

std::string qcat(const std::string& path)
{
  TIMER("qcat");
  QFile qfile = qfopen(path, "r");
  if (qfile.null()) {
    return "";
  }
  qfseek(qfile, 0, SEEK_END);
  const long length = qftell(qfile);
  qfseek(qfile, 0, SEEK_SET);
  std::string ret(length, 0);
  const long length_actual = qfread(&ret[0], 1, length, qfile);
  qassert(length == length_actual);
  qfclose(qfile);
  return ret;
}

int qar_create(const std::string& path_qar, const std::string& path_folder_,
               const bool is_remove_folder_after)
// interface function
// return 0 if successful.
{
  TIMER_VERBOSE_FLOPS("qar_create");
  const std::string path_folder = remove_trailing_slashes(path_folder_);
  displayln(
      0, fname + ssprintf(
                     ": '%s' '%s' %s.", path_qar.c_str(), path_folder.c_str(),
                     is_remove_folder_after ? "remove folder" : "keep folder"));
  if (not is_directory(path_folder)) {
    qwarn(fname + ssprintf(": '%s' '%s' no folder.", path_qar.c_str(),
                           path_folder.c_str()));
    return 1;
  }
  if (does_file_exist(path_qar)) {
    qwarn(fname + ssprintf(": '%s' '%s' qar already exist.", path_qar.c_str(),
                           path_folder.c_str()));
    return 2;
  }
  if (does_file_exist(path_qar + ".acc")) {
    qwarn(fname + ssprintf(": '%s' '%s' qar.acc already exist.",
                           path_qar.c_str(), path_folder.c_str()));
    return 3;
  }
  const std::string path_prefix = path_folder + "/";
  const long path_prefix_len = path_prefix.size();  // including the final "/"
  const std::vector<std::string> contents = qls_all(path_folder);
  std::vector<std::string> reg_files;
  for (long i = 0; i < (long)contents.size(); ++i) {
    const std::string path = contents[i];
    qassert(path.substr(0, path_prefix_len) == path_prefix);
    if (not is_directory(path)) {
      if (not is_regular_file(path)) {
        qwarn(fname + ssprintf(": '%s' '%s' '%s' not regular file.",
                               path_qar.c_str(), path_folder.c_str(),
                               path.c_str()));
        return 4;
      }
      reg_files.push_back(path);
    }
  }
  long num_vol;
  {
    long iv = 0;
    QarFileVol qar;
    for (long i = 0; i < (long)reg_files.size(); ++i) {
      const std::string path = reg_files[i];
      const std::string fn = path.substr(path_prefix_len);
      QFile qfile_in(path, "r");
      qassert(not qfile_in.null());
      if (not qar.null()) {
        const long total_size_of_current_vol = qftell(qar.qfile());
        const long data_size = qfile_remaining_size(qfile_in);
        if (get_qar_multi_vol_max_size() >= 0 and
            total_size_of_current_vol + data_size >
                get_qar_multi_vol_max_size()) {
          qar.close();
          iv += 1;
        }
      }
      const std::string path_qar_v = path_qar + qar_file_multi_vol_suffix(iv);
      displayln(1, fname + ssprintf(": '%s' '%s'/'%s' %ld/%ld.",
                                    path_qar_v.c_str(), path_prefix.c_str(),
                                    fn.c_str(), i + 1, reg_files.size()));
      if (qar.null()) {
        qar.init(path_qar_v + ".acc", "w");
        qassert(not qar.null());
      }
      timer.flops += write_from_qfile(qar, fn, "", qfile_in);
    }
    qar.close();
    num_vol = iv + 1;
  }
  for (long iv = 0; iv < num_vol; ++iv) {
    const std::string path_qar_v = path_qar + qar_file_multi_vol_suffix(iv);
    qrename(path_qar_v + ".acc", path_qar_v);
  }
  if (is_remove_folder_after) {
    for (long iv = 0; iv < num_vol; ++iv) {
      const std::string path_qar_v = path_qar + qar_file_multi_vol_suffix(iv);
      qassert(does_file_exist(path_qar_v));
    }
    qremove_all(path_folder);
  }
  return 0;
}

int qar_extract(const std::string& path_qar, const std::string& path_folder_,
                const bool is_remove_qar_after)
// interface function
// return 0 if successful.
{
  TIMER_VERBOSE_FLOPS("qar_extract");
  const std::string path_folder = remove_trailing_slashes(path_folder_);
  displayln(
      0,
      fname + ssprintf(": '%s' '%s' %s.", path_qar.c_str(), path_folder.c_str(),
                       is_remove_qar_after ? "remove qar " : "keep qar"));
  if (not does_regular_file_exist_qar(path_qar)) {
    qwarn(fname + ssprintf(": '%s' '%s' qar does not exist.", path_qar.c_str(),
                           path_folder.c_str()));
    return 1;
  }
  if (does_file_exist(path_folder)) {
    qwarn(fname + ssprintf(": '%s' '%s' folder exist.", path_qar.c_str(),
                           path_folder.c_str()));
    return 2;
  }
  if (does_file_exist(path_folder + ".acc")) {
    qwarn(fname + ssprintf(": '%s' '%s' folder.acc already exist.",
                           path_qar.c_str(), path_folder.c_str()));
    return 3;
  }
  QarFile qar(path_qar, "r");
  const std::vector<std::string> contents = list(qar);
  std::set<std::string> dirs;
  qmkdir_p(path_folder + ".acc");
  dirs.insert(".");
  for (long i = 0; i < (long)contents.size(); ++i) {
    const std::string& fn = contents[i];
    const std::string dn = dirname(fn);
    if (not has(dirs, dn)) {
      const int code = qmkdir_p(path_folder + ".acc/" + dn);
      qassert(code == 0);
      dirs.insert(dn);
    }
    QFile qfile_in;
    const bool b = read(qar, fn, qfile_in);
    qassert(b);
    QFile qfile_out(path_folder + ".acc/" + fn, "w");
    qassert(not qfile_out.null());
    timer.flops += write_from_qfile(qfile_out, qfile_in);
  }
  const long num_vol = qar.size();
  qar.close();
  qrename(path_folder + ".acc", path_folder);
  if (is_remove_qar_after) {
    qassert(is_directory(path_folder));
    for (long iv = 0; iv < num_vol; ++iv) {
      const std::string path_qar_v = path_qar + qar_file_multi_vol_suffix(iv);
      qremove(path_qar_v);
    }
  }
  return 0;
}

int qcopy_file(const std::string& path_src, const std::string& path_dst)
// interface function
// return 0 if successful.
{
  TIMER_VERBOSE_FLOPS("qcopy_file");
  displayln(0,
            fname + ssprintf(": '%s' %s.", path_src.c_str(), path_dst.c_str()));
  if (not does_regular_file_exist_qar(path_src)) {
    qwarn(fname + ssprintf(": '%s' does not exist.", path_src.c_str()));
    return 1;
  }
  if (does_file_exist(path_dst)) {
    qwarn(fname + ssprintf(": '%s' already exist.", path_dst.c_str()));
    return 2;
  }
  if (does_file_exist(path_dst + ".acc")) {
    qwarn(fname +
          ssprintf(": '%s' already exist.", (path_dst + ".acc").c_str()));
    return 3;
  }
  QFile qfile_in = qfopen(path_src, "r");
  qassert(not qfile_in.null());
  QFile qfile_out = qfopen(path_dst + ".acc", "w");
  qassert(not qfile_out.null());
  timer.flops += write_from_qfile(qfile_out, qfile_in);
  qfclose(qfile_out);
  qfclose(qfile_in);
  qrename(path_dst + ".acc", path_dst);
  return 0;
}

std::vector<std::string> list_qar(const std::string& path)
{
  TIMER_VERBOSE("list_qar");
  QarFile qar(path, "r");
  return list(qar);
}

bool has_regular_file(const QarFile& qar, const std::string& fn)
// interface function
{
  for (unsigned long i = 0; i < qar.size(); ++i) {
    const QarFileVol& qar_v = qar[i];
    qassert(not qar_v.null());
    qassert(qar_v.mode() == "r");
    if (has_regular_file(qar_v, fn)) {
      return true;
    }
  }
  return false;
}

bool has(const QarFile& qar, const std::string& fn)
// interface function
{
  for (unsigned long i = 0; i < qar.size(); ++i) {
    const QarFileVol& qar_v = qar[i];
    qassert(not qar_v.null());
    qassert(qar_v.mode() == "r");
    if (has(qar_v, fn)) {
      return true;
    }
  }
  return false;
}

bool read(const QarFile& qar, const std::string& fn, QFile& qfile_in)
// interface function
{
  for (unsigned long i = 0; i < qar.size(); ++i) {
    const QarFileVol& qar_v = qar[i];
    qassert(not qar_v.null());
    qassert(qar_v.mode() == "r");
    if (read(qar_v, fn, qfile_in)) {
      return true;
    }
  }
  return false;
}

}  // namespace qlat
