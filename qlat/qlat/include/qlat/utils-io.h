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

Int mkdir_lock(const std::string& path, const mode_t mode = default_dir_mode());

Int mkdir_lock_all_node(const std::string& path,
                        const mode_t mode = default_dir_mode());

Int rmdir_lock(const std::string& path);

Int rmdir_lock_all_node(const std::string& path);

void check_sigterm();

bool check_status();

bool obtain_lock_all_node(const std::string& path);

void release_lock_all_node();

}  // namespace qlat
