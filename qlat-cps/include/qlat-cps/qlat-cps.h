#pragma once

#include <qlat/cps.h>

namespace qlat
{

void begin_with_cps(const std::vector<std::string>& sargs,
                    const Coordinate& total_site);

void end_with_cps(const bool is_preserving_cache = false);

}  // namespace qlat
