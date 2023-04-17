#pragma once

#include <qlat/cps.h>

namespace qlat
{

void begin_with_cps(const std::vector<std::string>& sargs,
                    const Coordinate& total_site);

void end_with_cps(const bool is_preserving_cache = false);

void save_cps_prop_double(const Field<WilsonMatrix>& prop, const std::string& path);

void load_cps_prop_double(Field<WilsonMatrix>& prop, const std::string& path);

}  // namespace qlat
