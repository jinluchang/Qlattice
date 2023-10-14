#pragma once

#include <qlat/grid.h>

namespace qlat
{

void begin_with_grid(
    const std::vector<std::string>& sargs = std::vector<std::string>(),
    const std::vector<Coordinate>& node_size_list = std::vector<Coordinate>());

void end_with_grid(const bool is_preserving_cache = false);

void save_grid_prop_float(const Field<WilsonMatrix>& prop, const std::string& path);

void save_grid_prop_double(const Field<WilsonMatrix>& prop, const std::string& path);

void load_grid_prop_float(Field<WilsonMatrix>& prop, const std::string& path);

void load_grid_prop_double(Field<WilsonMatrix>& prop, const std::string& path);

}  // namespace qlat
