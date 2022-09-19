#pragma once

#include <qlat/qlat.h>

namespace qlat
{  //

inline std::vector<int> get_todo_trajs(const std::string& job_tag)
{
  std::vector<int> ret;
  if (job_tag == "64I-0.000678") {
    for (int i = 500; i < 3000; i += 10) {
      ret.push_back(i);
    }
  } else if (job_tag == "48I-0.00078" or job_tag == "32Dfine-0.0001") {
    for (int i = 500; i < 3000; i += 10) {
      ret.push_back(i);
    }
  } else if (job_tag == "24D-0.00107") {
    for (int i = 1010; i < 3000; i += 10) {
      ret.push_back(i);
    }
  } else if (job_tag == "24D-0.0174") {
    for (int i = 210; i < 2000; i += 10) {
      ret.push_back(i);
    }
  } else if (job_tag == "32D-0.00107") {
    for (int i = 690; i < 3000; i += 10) {
      ret.push_back(i);
    }
  } else if (job_tag == "16I-0.01") {
    qassert(false);
  } else {
    qassert(false);
  }
  return ret;
}

}  // namespace qlat
