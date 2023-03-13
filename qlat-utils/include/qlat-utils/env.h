#pragma once

#include <qlat-utils/timer.h>
#include <dirent.h>

namespace qlat
{  //

API inline long& get_qar_multi_vol_max_size()
// qlat parameter
// size in bytes
{
  static long size = get_env_long_default("q_qar_multi_vol_max_size",
                                          500L * 1000L * 1000L * 1000L);
  return size;
}

API inline mode_t& default_dir_mode()
// qlat parameter
{
  static mode_t mode = 0775;
  return mode;
}

API inline long& write_from_qfile_chunk_size()
// qlat parameter
// size in bytes
{
  static long size = 512 * 1024;
  return size;
}

}  // namespace qlat
