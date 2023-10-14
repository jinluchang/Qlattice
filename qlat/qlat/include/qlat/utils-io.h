// vim: set ts=2 sw=2 expandtab:

#pragma once

#include <qlat/config.h>
#include <qlat/mpi.h>
#include <sys/types.h>
#include <unistd.h>

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

namespace qlat
{  //

inline double get_remaining_time()
{
  return get_time_limit() - get_actual_total_time();
}

bool obtain_lock(const std::string& path);

void release_lock();

void qquit(const std::string& msg);

void check_time_limit(const double budget = get_default_budget());

void check_stop(const std::string& fn = "stop.txt");

bool does_file_exist_sync_node(const std::string& fn);

bool does_regular_file_exist_qar_sync_node(const std::string& fn);

bool does_file_exist_qar_sync_node(const std::string& fn);

bool is_directory_sync_node(const std::string& fn);

bool is_regular_file_sync_node(const std::string& fn);

int qmkdir_sync_node(const std::string& path,
                     const mode_t mode = default_dir_mode());

int mkdir_lock(const std::string& path, const mode_t mode = default_dir_mode());

int mkdir_lock_all_node(const std::string& path,
                        const mode_t mode = default_dir_mode());

int rmdir_lock(const std::string& path);

int rmdir_lock_all_node(const std::string& path);

std::vector<std::string> qls_sync_node(const std::string& path);

std::vector<std::string> qls_all_sync_node(
    const std::string& path, const bool is_folder_before_files = false);

int qremove_info(const std::string& path);

int qremove_all_info(const std::string& path);

int qar_create_info(const std::string& path_qar,
                    const std::string& path_folder_,
                    const bool is_remove_folder_after = false);

int qar_extract_info(const std::string& path_qar,
                     const std::string& path_folder_,
                     const bool is_remove_qar_after = false);

int qcopy_file_info(const std::string& path_src, const std::string& path_dst);

std::string qcat_sync_node(const std::string& path);

void load_qar_index_sync_node(const QarFile& qar, const std::string& fn);

DataTable qload_datatable_sync_node(const std::string& path,
                                    const bool is_par = false);

LatData lat_data_load_info(const std::string& path);

void check_sigterm();

bool check_status();

bool obtain_lock_all_node(const std::string& path);

void release_lock_all_node();

}  // namespace qlat
