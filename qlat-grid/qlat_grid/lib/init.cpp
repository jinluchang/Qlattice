#include <qlat-grid/qlat-grid.h>

namespace qlat
{  //

void begin_with_grid(const std::vector<std::string>& sargs,
                     const std::vector<Coordinate>& node_size_list)
{
  // make cargs
  std::vector<const char*> cargs(sargs.size() + 1);
  for (long i = 0; i < (long)sargs.size(); ++i) {
    cargs[i] = sargs[i].c_str();
  }
  cargs.back() = NULL;
  //
  int argc = (int)sargs.size();
  char** argv = (char**)&cargs[0];
  //
  grid_begin(&argc, &argv, node_size_list);
}

void end_with_grid(const bool is_preserving_cache)
{
  grid_end(is_preserving_cache);
}

}  // namespace qlat
