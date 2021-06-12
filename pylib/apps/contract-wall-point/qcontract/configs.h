#pragma once

#include <qlat/qlat.h>

namespace qlat
{  //

inline std::vector<int> get_data_trajs(const std::string& job_tag)
{
  std::vector<int> ret;
  if (job_tag == "64I") {
    for (int i = 500; i < 3000; i += 10) {
      ret.push_back(i);
    }
  } else if (job_tag == "48I" or job_tag == "32Dfine") {
    for (int i = 500; i < 3000; i += 10) {
      ret.push_back(i);
    }
  } else if (job_tag == "test-4nt16") {
    for (int i = 1000; i < 2000; i += 100) {
      ret.push_back(i);
    }
  } else if (job_tag == "24D") {
    for (int i = 1010; i < 3000; i += 10) {
      ret.push_back(i);
    }
  } else if (job_tag == "24DH") {
    for (int i = 210; i < 2000; i += 10) {
      ret.push_back(i);
    }
  } else if (job_tag == "32D") {
    for (int i = 690; i < 3000; i += 10) {
      ret.push_back(i);
    }
  } else {
    qassert(false);
  }
  return ret;
}

}  // namespace qlat
