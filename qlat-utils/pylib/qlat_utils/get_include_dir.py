import qlat_utils.lat_io
import os

def get_qlat_utils_include():
    return os.path.join(os.path.dirname(qlat_utils.lat_io.__file__), 'include')

def get_include_list():
    return [ get_qlat_utils_include(), ]
