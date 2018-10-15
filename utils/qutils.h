#pragma once

#include <vector>

namespace qutils {
  
template <class M>
void clear(std::vector<M>& vec)
{
  std::vector<M> empty;
  swap(empty, vec);
}

}
