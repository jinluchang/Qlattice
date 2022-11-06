import os
import qlat_grid.init
import qlat as q

def get_qlat_grid_include():
    return os.path.join(os.path.dirname(qlat_grid.init.__file__), 'include')

def get_include_list():
    return [ get_qlat_grid_include(), ] + q.get_include_list()
