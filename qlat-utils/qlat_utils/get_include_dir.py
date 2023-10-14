import os

def get_qlat_utils_dir():
    return os.path.dirname(__file__)

def get_qlat_utils_include():
    return os.path.join(get_qlat_utils_dir(), 'include')

def get_dir_list():
    return [ get_qlat_utils_dir(), ]

def get_include_list():
    return [ os.path.join(p, 'include') for p in get_dir_list() ]

def get_lib_list():
    return [ os.path.join(p, 'lib') for p in get_dir_list() ]

def get_new_ld_library_path():
    ld_lib_path = os.getenv('LD_LIBRARY_PATH')
    path_list = ld_lib_path.split(':')
    new_path_list = []
    for p in get_lib_list() + path_list:
        if p not in new_path_list:
            new_path_list.append(p)
    return ':'.join(new_path_list)

def get_pxd_list():
    from .utils_io import qls, qls_all
    l = []
    for d in get_dir_list():
        fn_list = qls(d)
        for fn in fn_list:
            if fn.endswith(".pxd"):
                l.append(fn)
    return l

def get_header_list():
    from .utils_io import qls, qls_all
    l = []
    for d in get_include_list():
        fn_list = qls_all(d)
        for fn in fn_list:
            if fn.endswith(".h"):
                l.append(fn)
    return l
