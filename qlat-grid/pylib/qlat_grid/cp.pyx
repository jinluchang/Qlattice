from . cimport everything as c

def end_with_grid(is_preserving_cache = False):
    c.grid_end(is_preserving_cache)
