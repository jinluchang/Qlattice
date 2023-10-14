#include <qlat-cps/qlat-cps.h>

namespace qlat
{  //

void begin_with_cps(const std::vector<std::string>& sargs,
                    const Coordinate& total_site)
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
  cps_begin(&argc, &argv, total_site);
}

void end_with_cps(const bool is_preserving_cache)
{
  cps_end(is_preserving_cache);
}

}  // namespace qlat
