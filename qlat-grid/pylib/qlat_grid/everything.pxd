from qlat.everything cimport *

cdef extern from "qlat-grid/qlat-grid.h" namespace "qlat":

    void begin_with_grid(const std_vector[std_string]& sargs,
                         const std_vector[Coordinate]& node_size_list)

    void end_with_grid(const bool is_preserving_cache)

    void save_grid_prop_float(const Field[WilsonMatrix]& prop, const std_string& path)

    void save_grid_prop_double(const Field[WilsonMatrix]& prop, const std_string& path)

    void load_grid_prop_float(Field[WilsonMatrix]& prop, const std_string& path)

    void load_grid_prop_double(Field[WilsonMatrix]& prop, const std_string& path)
