// vim: set ts=2 sw=2 expandtab:

#pragma once

#include <qlat/config.h>
#include <qlat/mpi.h>

#include <show.h>

#include <array>
#include <vector>
#include <iostream>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdarg>
#include <fstream>

#include <sys/stat.h>
#include <unistd.h>

QLAT_START_NAMESPACE

inline bool truncate(const std::string &evilFile) {
  std::ofstream evil;
  evil.open(evilFile.c_str());
  bool does_exist = evil.good();
  if(does_exist) evil.close();
  return does_exist;
}

inline bool does_file_exist(const std::string& fn) {
  struct stat sb;
  return 0 == stat(fn.c_str(), &sb);
}

inline bool does_file_exist_sync_node(const std::string& fn) {
  long nfile = 0;
  if (0 == get_id_node()) {
    if (does_file_exist(fn)) {
      nfile = 1;
    }
  }
  glb_sum(nfile);
  return 0 != nfile;
}

inline mode_t& default_dir_mode()
{
  static mode_t mode = 0775;
  return mode;
}

inline int ssleep(const double seconds) {
  return usleep((useconds_t)(seconds * 1.0e6));
}

inline int check_dir(const std::string& path, const mode_t mode = default_dir_mode())
{
  TIMER("check_dir");
  int ret = 0;
  while (!does_file_exist(path)) {
    ret = mkdir(path.c_str(), mode);
    ssleep(0.001);
  }
  return ret;
}

inline int check_mkdir(const std::string& path, const mode_t mode = default_dir_mode())
{
  TIMER("check_mkdir");
  mkdir(path.c_str(), mode);
  return check_dir(path, mode);
}

inline int check_mkdir_sync_node(const std::string& path, const mode_t mode = default_dir_mode())
{
  TIMER("check_mkdir_sync_node");
  if (0 == get_id_node()) {
    check_mkdir(path, mode);
  }
  sync_node();
  return check_dir(path, mode);
}

inline int mkdir_sync_node(const std::string& path, const mode_t mode = default_dir_mode())
{
  TIMER("mkdir_sync_node");
  long ret = 0;
  if (0 == get_id_node()) {
    ret = mkdir(path.c_str(), mode);
  }
  glb_sum(ret);
  return ret;
}

inline int rmdir_sync_node(const std::string& path)
{
  TIMER("rmdir_sync_node");
  long ret = 0;
  if (0 == get_id_node()) {
    ret = rmdir(path.c_str());
  }
  glb_sum(ret);
  return ret;
}

QLAT_END_NAMESPACE
