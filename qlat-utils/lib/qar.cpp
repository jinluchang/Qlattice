#include <qlat-utils/qar-cache.h>

namespace qlat
{  //

static void check_all_files_crc32_aux(
    std::vector<std::pair<std::string, crc32_t> >& acc, const std::string& path)
{
  if (not is_directory(path)) {
    acc.push_back(check_file_crc32(path));
  } else {
    const std::vector<std::string> paths = qls_aux(path);
    for (long i = 0; i < (long)paths.size(); ++i) {
      check_all_files_crc32_aux(acc, paths[i]);
    }
  }
}

// ----------------------------------------------------

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
  TIMER("save_qar_index");
  std::vector<std::string> lines;
  lines.push_back(qar_idx_header);
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
          ssprintf("QAR-FILE-IDX %ld %ld %ld\n", i, j, fn.size());
      const std::string line2 = fn + "\n";
      const std::string line3 = ssprintf(
          "%ld %ld %ld %ld %ld %ld %ld %ld\n", qsinfo.offset, qsinfo.offset_fn,
          qsinfo.offset_info, qsinfo.offset_data, qsinfo.offset_end,
          qsinfo.fn_len, qsinfo.info_len, qsinfo.data_len);
      lines.push_back(line1);
      lines.push_back(line2);
      lines.push_back(line3);
      lines.push_back("\n");
    }
  }
  qtouch(fn, lines);
}

void save_qar_index_info(const QarFile& qar, const std::string& fn)
// interface function
{
  if (0 == get_id_node()) {
    save_qar_index(qar, fn);
  }
}

void parse_qar_index(const QarFile& qar, const std::string& qar_index_content)
// interface function
{
  if (qar.null()) {
    qwarn("parse_qar_index: qar is null.");
    return;
  }
  if (qar_index_content == "") {
    return;
  }
  const long header_len = qar_idx_header.size();
  if (0 != qar_index_content.compare(0, header_len, qar_idx_header)) {
    qwarn("parse_qar_index: not qar-idx file format.");
    return;
  }
  long cur = header_len;
  while (cur < (long)qar_index_content.size()) {
    long i, j, fn_len;
    std::string fn;
    QarSegmentInfo qsinfo;
    long cur1 = 0;
    long cur3 = 0;
    std::string line1;
    std::string line3;
    // line1: QAR-FILE-IDX i j fn_len
    if (not parse_line(line1, cur, qar_index_content)) {
      qwarn("parse_qar_index: not qar-idx file format.");
      return;
    }
    if (not parse_literal(cur1, line1, "QAR-FILE-IDX ")) {
      qwarn(ssprintf(
          "parse_qar_index: not qar-idx file format. cur1=%ld line1='%s'.", cur1,
          line1.c_str()));
      return;
    }
    if (not parse_long(i, cur1, line1)) {
      qwarn("parse_qar_index: not qar-idx file format.");
      return;
    }
    if (not parse_literal(cur1, line1, ' ')) {
      qwarn("parse_qar_index: not qar-idx file format.");
      return;
    }
    if (not parse_long(j, cur1, line1)) {
      qwarn("parse_qar_index: not qar-idx file format.");
      return;
    }
    if (not parse_literal(cur1, line1, ' ')) {
      qwarn("parse_qar_index: not qar-idx file format.");
      return;
    }
    if (not parse_long(fn_len, cur1, line1)) {
      qwarn("parse_qar_index: not qar-idx file format.");
      return;
    }
    if (not parse_literal(cur1, line1, '\n')) {
      qwarn(ssprintf(
          "parse_qar_index: not qar-idx file format. cur1=%ld line1='%s'.",
          cur1, line1.c_str()));
      return;
    }
    if (not parse_end(cur1, line1)) {
      qwarn("parse_qar_index: not qar-idx file format.");
      return;
    }
    // line2: fn
    if (not parse_len(fn, cur, qar_index_content, fn_len)) {
      qwarn("parse_qar_index: not qar-idx file format.");
      return;
    }
    if (not parse_literal(cur, qar_index_content, '\n')) {
      qwarn("parse_qar_index: not qar-idx file format.");
      return;
    }
    // line3: offset offset_fn offset_info offset_data offset_end fn_len
    // info_len data_len
    if (not parse_line(line3, cur, qar_index_content)) {
      qwarn("parse_qar_index: not qar-idx file format.");
      return;
    }
    if (not parse_long(qsinfo.offset, cur3, line3)) {
      qwarn("parse_qar_index: not qar-idx file format.");
      return;
    }
    if (not parse_literal(cur3, line3, ' ')) {
      qwarn("parse_qar_index: not qar-idx file format.");
      return;
    }
    if (not parse_long(qsinfo.offset_fn, cur3, line3)) {
      qwarn("parse_qar_index: not qar-idx file format.");
      return;
    }
    if (not parse_literal(cur3, line3, ' ')) {
      qwarn("parse_qar_index: not qar-idx file format.");
      return;
    }
    if (not parse_long(qsinfo.offset_info, cur3, line3)) {
      qwarn("parse_qar_index: not qar-idx file format.");
      return;
    }
    if (not parse_literal(cur3, line3, ' ')) {
      qwarn("parse_qar_index: not qar-idx file format.");
      return;
    }
    if (not parse_long(qsinfo.offset_data, cur3, line3)) {
      qwarn("parse_qar_index: not qar-idx file format.");
      return;
    }
    if (not parse_literal(cur3, line3, ' ')) {
      qwarn("parse_qar_index: not qar-idx file format.");
      return;
    }
    if (not parse_long(qsinfo.offset_end, cur3, line3)) {
      qwarn("parse_qar_index: not qar-idx file format.");
      return;
    }
    if (not parse_literal(cur3, line3, ' ')) {
      qwarn("parse_qar_index: not qar-idx file format.");
      return;
    }
    if (not parse_long(qsinfo.fn_len, cur3, line3)) {
      qwarn("parse_qar_index: not qar-idx file format.");
      return;
    }
    if (not parse_literal(cur3, line3, ' ')) {
      qwarn("parse_qar_index: not qar-idx file format.");
      return;
    }
    if (not parse_long(qsinfo.info_len, cur3, line3)) {
      qwarn("parse_qar_index: not qar-idx file format.");
      return;
    }
    if (not parse_literal(cur3, line3, ' ')) {
      qwarn("parse_qar_index: not qar-idx file format.");
      return;
    }
    if (not parse_long(qsinfo.data_len, cur3, line3)) {
      qwarn("parse_qar_index: not qar-idx file format.");
      return;
    }
    if (not parse_literal(cur3, line3, '\n')) {
      qwarn("parse_qar_index: not qar-idx file format.");
      return;
    }
    if (not parse_end(cur3, line3)) {
      qwarn("parse_qar_index: not qar-idx file format.");
      return;
    }
    // line4:
    if (not parse_literal(cur, qar_index_content, '\n')) {
      qwarn("parse_qar_index: not qar-idx file format.");
      return;
    }
    // register_file
    qassert(0 <= i);
    qassert(i < (long)qar.size());
    register_file(qar[i], fn, qsinfo);
  }
}

void load_qar_index(const QarFile& qar, const std::string& fn)
// interface function
{
  TIMER_VERBOSE("load_qar_index");
  const std::string qar_index_content = qcat(fn);
  parse_qar_index(qar, qar_index_content);
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

void qar_build_index(const std::string& path_qar)
{
  TIMER_VERBOSE("qar_build_index");
  displayln(fname + ssprintf(": '%s'.", path_qar.c_str()));
  QarFile qar(path_qar, "r");
  save_qar_index(qar, path_qar + ".idx");
  qar.close();
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
  qar_build_index(path_qar);
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
    if (does_file_exist(path_qar + ".idx")) {
      qremove(path_qar + ".idx");
    }
    for (long iv = 0; iv < num_vol; ++iv) {
      const std::string path_qar_v = path_qar + qar_file_multi_vol_suffix(iv);
      if (does_file_exist(path_qar_v)) {
        qremove(path_qar_v);
      }
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

std::string mk_new_qar_read_cache_key(const QarFile& qar,
                                      const std::string& key,
                                      const std::string& path)
// (1) Find the first new qar file in qar that match the prefix of path and
// register the new qar file in qar_read_cache.
// (2) If qar not found, return key.
// (3) If path exists in the qar, return the new key of the new qar.
// (4) If not found, repeat the procedure for the new qar.
{
  Cache<std::string, QarFile>& cache = get_qar_read_cache();
  std::string path_dir = remove_trailing_slashes(path);
  const std::string pathd = path_dir + "/";
  while (true) {
    if (path_dir == "/" or path_dir == ".") {
      return key;
    }
    if (has_regular_file(qar, path_dir + ".qar")) {
      const std::string key_new = path_dir + "/";
      qassert(pathd.substr(0, key_new.size()) == key_new);
      QarFile& qar_new = cache[key + key_new];
      if (qar_new.null()) {
        qar_new.init(key + path_dir + ".qar", "r");
      }
      qassert(not qar_new.null());
      const std::string path_new = pathd.substr(key_new.size());
      if (has(qar_new, path_new)) {
        return key + key_new;
      } else {
        return mk_new_qar_read_cache_key(qar_new, key + key_new, path_new);
      }
    }
    path_dir = dirname(path_dir);
  }
  qassert(false);
  return "";
}

std::string mk_new_qar_read_cache_key(const std::string& path)
// (1) Find first qar file that match the prefix of path and register the qar
// file in qar_read_cache.
// (2) If qar not found, return "".
// (2) If path exists in the qar, return the key of qar.
// (4) If not found, find qar within qar recursively, return the key of the
// closest qar.
{
  Cache<std::string, QarFile>& cache = get_qar_read_cache();
  std::string path_dir = remove_trailing_slashes(path);
  const std::string pathd = path_dir + "/";
  while (true) {
    if (path_dir == "/" or path_dir == ".") {
      return "";
    }
    if (does_file_exist(path_dir + ".qar")) {
      const std::string key = path_dir + "/";
      qassert(pathd.substr(0, key.size()) == key);
      qassert(not cache.has(key));
      QarFile& qar = cache[key];
      qar.init(path_dir + ".qar", "r");
      qassert(not qar.null());
      const std::string path_new = pathd.substr(key.size());
      if (has(qar, path_new)) {
        return key;
      } else {
        return mk_new_qar_read_cache_key(qar, key, path_new);
      }
    }
    path_dir = dirname(path_dir);
  }
  qassert(false);
  return "";
}

std::string get_qar_read_cache_key(const std::string& path)
// return key of get_qar_read_cache() that may contain path
// return empty string if no cached key is found.
// Note: key should end with '/'.
// Steps:
// (1) Search in Cache. If found matching key, try to find within this qar file
// recursively. Return the key of the closest match.
// (2) If not found, check if path exists. If exists, return path.
// (3) If does not exist, try to find qar file yet to be in cache recursively.
// Return values:
// valid key: valid key for a qar found. (qar may not actually contain path).
// "": no key is found and path does not exist.
// path: path exist.
{
  TIMER("get_qar_read_cache_key");
  Cache<std::string, QarFile>& cache = get_qar_read_cache();
  for (auto it = cache.m.cbegin(); it != cache.m.cend(); ++it) {
    const std::string& key = it->first;
    if (key == path.substr(0, key.size())) {
      const QarFile& qar = cache[key];
      const std::string path_new = path.substr(key.size());
      if (has(qar, path_new)) {
        return key;
      } else {
        return mk_new_qar_read_cache_key(qar, key, path_new);
      }
    }
  }
  if (does_file_exist(path)) {
    return path;
  }
  return mk_new_qar_read_cache_key(path);
}

int qtouch(const std::string& path)
// return 0 if success
{
  TIMER("qtouch");
  QFile qfile = qfopen(path, "w");
  if (qfile.null()) {
    return 1;
  }
  qfclose(qfile);
  return 0;
}

int qtouch(const std::string& path, const std::string& content)
{
  TIMER("qtouch");
  QFile qfile = qfopen(path + ".partial", "w");
  if (qfile.null()) {
    return 1;
  }
  const long total_bytes = qwrite_data(content, qfile);
  qassert(total_bytes == long(content.size()));
  qfclose(qfile);
  return qrename(path + ".partial", path);
}

int qtouch(const std::string& path, const std::vector<std::string>& content)
{
  TIMER("qtouch");
  QFile qfile = qfopen(path + ".partial", "w");
  if (qfile.null()) {
    return 1;
  }
  long total_bytes = 0;
  long total_bytes_expect = 0;
  for (long i = 0; i < (long)content.size(); ++i) {
    total_bytes_expect += content[i].size();
    total_bytes += qwrite_data(content[i], qfile);
  }
  qassert(total_bytes == total_bytes_expect);
  qfclose(qfile);
  return qrename(path + ".partial", path);
}

int qappend(const std::string& path, const std::string& content)
{
  TIMER("qappend");
  QFile qfile = qfopen(path, "a");
  if (qfile.null()) {
    return 1;
  }
  const long total_bytes = qwrite_data(content, qfile);
  qassert(total_bytes == long(content.size()));
  qfclose(qfile);
  return 0;
}

DataTable qload_datatable_serial(QFile& qfile)
{
  TIMER("qload_datatable_serial(qfile)");
  DataTable ret;
  while (not qfeof(qfile)) {
    const std::string line = qgetline(qfile);
    if (line.length() > 0 && line[0] != '#') {
      const std::vector<double> xs = read_doubles(line);
      if (xs.size() > 0) {
        ret.push_back(xs);
      }
    }
  }
  return ret;
}

DataTable qload_datatable_par(QFile& qfile)
{
  TIMER("qload_datatable_qar(qfile)");
  const size_t line_buf_size = 1024;
  DataTable ret;
  std::vector<std::string> lines;
  DataTable xss;
  while (not qfeof(qfile)) {
    lines.clear();
    for (size_t i = 0; i < line_buf_size; ++i) {
      lines.push_back(qgetline(qfile));
      if (qfeof(qfile)) {
        break;
      }
    }
    xss.resize(lines.size());
#pragma omp parallel for
    for (size_t i = 0; i < lines.size(); ++i) {
      const std::string& line = lines[i];
      if (line.length() > 0 && line[0] != '#') {
        xss[i] = read_doubles(line);
      } else {
        clear(xss[i]);
      }
    }
    for (size_t i = 0; i < xss.size(); ++i) {
      if (xss[i].size() > 0) {
        ret.push_back(xss[i]);
      }
    }
  }
  return ret;
}

DataTable qload_datatable_serial(const std::string& path)
{
  TIMER("qload_datatable_serial(path)");
  if (not does_regular_file_exist_qar(path)) {
    return DataTable();
  }
  QFile qfile = qfopen(path, "r");
  qassert(not qfile.null());
  DataTable ret = qload_datatable_serial(qfile);
  qfclose(qfile);
  return ret;
}

DataTable qload_datatable_par(const std::string& path)
{
  TIMER("qload_datatable_par(path)");
  if (not does_regular_file_exist_qar(path)) {
    return DataTable();
  }
  QFile qfile = qfopen(path, "r");
  qassert(not qfile.null());
  DataTable ret = qload_datatable_par(qfile);
  qfclose(qfile);
  return ret;
}

DataTable qload_datatable(const std::string& path, const bool is_par)
{
  if (is_par) {
    return qload_datatable_par(path);
  } else {
    return qload_datatable_serial(path);
  }
}

crc32_t compute_crc32(QFile& qfile)
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

crc32_t compute_crc32(const std::string& path)
{
  QFile qfile = qfopen(path, "r");
  const crc32_t ret = compute_crc32(qfile);
  qfclose(qfile);
  return ret;
}

std::vector<std::pair<std::string, crc32_t> > check_all_files_crc32(
    const std::string& path)
{
  TIMER_VERBOSE("check_all_files_crc32");
  std::vector<std::pair<std::string, crc32_t> > ret;
  check_all_files_crc32_aux(ret, remove_trailing_slashes(path));
  return ret;
}

void check_all_files_crc32_info(const std::string& path)
// interface function
{
  TIMER_VERBOSE("check_all_files_crc32_info");
  if (0 == get_id_node()) {
    displayln(fname + ssprintf(": start checking path='%s'", path.c_str()));
    std::vector<std::pair<std::string, crc32_t> > fcrcs;
    fcrcs = check_all_files_crc32(path);
    displayln(fname + ssprintf(": summary for path='%s'", path.c_str()));
    display(show_files_crc32(fcrcs));
  }
}

}  // namespace qlat
