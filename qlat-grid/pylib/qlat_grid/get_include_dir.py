import os
import qlat_grid.init
import qlat as q

def get_qlat_grid_dir():
    return os.path.dirname(qlat_grid.__init__.__file__)

def get_qlat_grid_include():
    return os.path.join(get_qlat_grid_dir(), 'include')

def get_dir_list():
    return [ get_qlat_grid_dir(), ] + q.get_dir_list()

def get_include_list():
    return [ os.path.join(p, 'include') for p in get_dir_list() ]

def get_lib_list():
    return [ os.path.join(p, 'lib') for p in get_dir_list() ]
