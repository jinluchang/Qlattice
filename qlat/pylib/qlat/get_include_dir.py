import os
from . import c
from qlat_utils import qls, qls_all, get_dir_list as qu_get_dir_list

def get_qlat_dir():
    return os.path.dirname(c.__file__)

def get_qlat_include():
    return os.path.join(get_qlat_dir(), 'include')

def get_dir_list():
    return [ get_qlat_dir(), ] + qu_get_dir_list()

def get_include_list():
    return [ os.path.join(p, 'include') for p in get_dir_list() ]

def get_lib_list():
    return [ os.path.join(p, 'lib') for p in get_dir_list() ]

def get_pxd_list():
    l = []
    for d in get_dir_list():
        fn_list = qls(d)
        for fn in fn_list:
            if fn.endswith(".pxd"):
                l.append(fn)
    return l

def get_header_list():
    l = []
    for d in get_include_list():
        fn_list = qls_all(d)
        for fn in fn_list:
            if fn.endswith(".h"):
                l.append(fn)
    return l
