#include <qlat-utils/lib/utils-io.h>
#include <qlat-utils/qar-cache.h>
#include <qlat-utils/qutils-io.h>

namespace qlat
{  //

void flush() { fflush(get_output_file()); }

int cc_qtouch(const std::string& path) { return qtouch(path); }

int cc_qtouch_info(const std::string& path) { return qtouch_info(path); }

int cc_qtouch(const std::string& path, const std::string& content)
{
  return qtouch(path, content);
}

int cc_qtouch_info(const std::string& path, const std::string& content)
{
  return qtouch_info(path, content);
}

int cc_qappend(const std::string& path, const std::string& content)
{
  return qappend(path, content);
}

int cc_qappend_info(const std::string& path, const std::string& content)
{
  return qappend_info(path, content);
}

int cc_qrename(const std::string& old_path, const std::string& new_path)
{
  return qrename(old_path, new_path);
}

int cc_qrename_info(const std::string& old_path, const std::string& new_path)
{
  return qrename_info(old_path, new_path);
}

std::vector<std::string> cc_qls(const std::string& path, const bool is_sort)
{
  return qls(path, is_sort);
}

std::vector<std::string> cc_qls_all(const std::string& path,
                                    const bool is_folder_before_files,
                                    const bool is_sort)
{
  return qls_all(path, is_folder_before_files, is_sort);
}

bool cc_is_directory(const std::string& fn) { return is_directory(fn); }

bool cc_is_regular_file(const std::string& fn) { return is_regular_file(fn); }

int cc_qremove(const std::string& path) { return qremove(path); }

int cc_qremove_all(const std::string& path) { return qremove_all(path); }

mode_t& cc_default_dir_mode() { return default_dir_mode(); }

int cc_qmkdir(const std::string& path, const mode_t mode)
{
  return qmkdir(path, mode);
}

int cc_qmkdir_p(const std::string& path, const mode_t mode)
{
  return qmkdir_p(path, mode);
}

int cc_qmkdir_info(const std::string& path, const mode_t mode)
{
  return qmkdir_info(path, mode);
}

int cc_qmkdir_p_info(const std::string& path, const mode_t mode)
{
  return qmkdir_p_info(path, mode);
}

bool cc_does_file_exist(const std::string& fn) { return does_file_exist(fn); }

crc32_t cc_compute_crc32(const std::string& path)
{
  return compute_crc32(path);
}

void cc_check_all_files_crc32_info(const std::string& path)
{
  return check_all_files_crc32_info(path);
}

DataTable cc_qload_datatable(const std::string& path, const bool is_par)
{
  return qload_datatable(path, is_par);
}

}  // namespace qlat
