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

bool obtain_lock(const std::string& path);

void release_lock();

void qquit(const std::string& msg);

void check_time_limit(const double budget = get_time_budget());

void check_stop(const std::string& fn = "stop.txt");

int mkdir_lock(const std::string& path, const mode_t mode = default_dir_mode());

int mkdir_lock_all_node(const std::string& path,
                        const mode_t mode = default_dir_mode());

int rmdir_lock(const std::string& path);

int rmdir_lock_all_node(const std::string& path);

void check_sigterm();

bool check_status();

bool obtain_lock_all_node(const std::string& path);

void release_lock_all_node();

// ------------------------

std::vector<std::string> qls_sync_node(const std::string& path);

std::vector<std::string> qls_all_sync_node(
    const std::string& path, const bool is_folder_before_files = false);

std::string qcat_sync_node(const std::string& path);

DataTable qload_datatable_sync_node(const std::string& path,
                                    const bool is_par = false);

LatData lat_data_load_sync_node(const std::string& path);

void load_qar_index_sync_node(const QarFile& qar, const std::string& fn);

}  // namespace qlat
