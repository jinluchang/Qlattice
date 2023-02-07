#pragma once

#include <qlat-utils/types.h>

#include <string>
#include <vector>

namespace qlat
{  //

void flush();

int cc_qtouch(const std::string& path);

int cc_qtouch_info(const std::string& path);

int cc_qtouch(const std::string& path, const std::string& content);

int cc_qtouch_info(const std::string& path, const std::string& content);

int cc_qappend(const std::string& path, const std::string& content);

int cc_qappend_info(const std::string& path, const std::string& content);

int cc_qrename(const std::string& old_path, const std::string& new_path);

int cc_qrename_info(const std::string& old_path, const std::string& new_path);

std::vector<std::string> cc_qls(const std::string& path,
                                const bool is_sort = true);

std::vector<std::string> cc_qls_all(const std::string& path,
                                    const bool is_folder_before_files = false,
                                    const bool is_sort = true);

bool cc_is_directory(const std::string& fn);

bool cc_is_regular_file(const std::string& fn);

int cc_qremove(const std::string& path);

int cc_qremove_all(const std::string& path);

mode_t& cc_default_dir_mode();

int cc_qmkdir(const std::string& path, const mode_t mode = cc_default_dir_mode());

int cc_qmkdir_p(const std::string& path,
                const mode_t mode = cc_default_dir_mode());

int cc_qmkdir_info(const std::string& path,
                   const mode_t mode = cc_default_dir_mode());

int cc_qmkdir_p_info(const std::string& path,
                     const mode_t mode = cc_default_dir_mode());

bool cc_does_file_exist(const std::string& fn);

crc32_t cc_compute_crc32(const std::string& path);

void cc_check_all_files_crc32_info(const std::string& path);

DataTable cc_qload_datatable(const std::string& path,
                             const bool is_par = false);

}  // namespace qlat
