import os
import qlat.field
import qlat_utils as qu

def get_qlat_dir():
    return os.path.dirname(qlat.field.__file__)

def get_qlat_include():
    return os.path.join(get_qlat_dir(), 'include')

def get_dir_list():
    return [ get_qlat_dir(), ] + qu.get_dir_list()

def get_include_list():
    return [ os.path.join(p, 'include') for p in get_dir_list() ]

def get_lib_list():
    return [ os.path.join(p, 'lib') for p in get_dir_list() ]
