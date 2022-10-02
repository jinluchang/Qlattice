import cqlat_grid_io as c

import sys

def begin_with_grid(argv = None, size_node_list = None):
    if size_node_list is None:
        size_node_list = []
    if argv is None:
        argv = sys.argv
    return c.begin_with_grid(argv, size_node_list)

def end_with_grid(is_preserving_cache = False):
    return c.end_with_grid(is_preserving_cache)
