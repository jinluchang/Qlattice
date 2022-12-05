import qlat_utils.lat_io
import os

def get_qlat_utils_dir():
    return os.path.dirname(qlat_utils.lat_io.__file__)

def get_qlat_utils_include():
    return os.path.join(get_qlat_utils_dir(), 'include')

def get_dir_list():
    return [ get_qlat_utils_dir(), ]

def get_include_list():
    return [ os.path.join(p, 'include') for p in get_dir_list() ]

def get_lib_list():
    return [ os.path.join(p, 'lib') for p in get_dir_list() ]
