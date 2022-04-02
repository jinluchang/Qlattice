#pragma once

#include <stdint.h>
#include <qutils/qutils-vec.h>
#include <qutils/qutils-io.h>

#include <cassert>
#include <string>
#include <vector>

namespace qlat
{  //

const std::string lat_data_header = "#!/usr/bin/env qar-glimpse\n";

struct QarFile {
  std::string path;
  std::string mode;  // can be "r", "a", "w"
  FILE* fp;
  //
  QarWriter()
  {
    fp = NULL;
    init();
  }
  ~QarWriter() { close(); }
  //
  void init()
  {
    close();
    path = "";
    qassert(fp == NULL);
  }
  void init(const std::string& path_, const std::string& mode_)
  {
    close();
    path = path_;
    mode = mode_;
    qassert(fp == NULL);
    displayln(
        ssprintf("QarFile: open '%s' with '%s'.", path.c_str(), mode.c_str()));
    fp = qopen(path, mode);
    qassert(NULL != fp);
  }
  //
  void close()
  {
    const bool is_need_close = fp != NULL;
    if (is_need_close) {
      displayln(ssprintf("QarFile: close '%s' with '%s'.", path.c_str(),
                         mode.c_str()));
    }
    qclose(fp);
    qassert(fp == NULL);
  }
};

}  // namespace qlat
