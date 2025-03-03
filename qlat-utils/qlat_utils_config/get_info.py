import os
from .config import *

def qls(directory):
    """List full paths of all files (excluding directories) in the given directory."""
    regular_files = []
    with os.scandir(directory) as entries:
        for entry in entries:
            if entry.is_file():
                # Append the full path of the file
                regular_files.append(entry.path)
    return regular_files

def qls_all(directory):
    """Recursively list all files in a directory."""
    file_list = []
    for root, dirs, files in os.walk(directory):
        for filename in files:
            filepath = os.path.join(root, filename)
            # Check if it's a regular file and not a symbolic link
            if os.path.isfile(filepath):
                file_list.append(filepath)
    return file_list

def get_qlat_utils_dir():
    return os.path.join(os.path.dirname(os.path.dirname(__file__)), 'qlat_utils')

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
    if ld_lib_path is None:
        path_list = []
    else:
        path_list = ld_lib_path.split(':')
    new_path_list = []
    for p in get_lib_list() + path_list:
        if p not in new_path_list:
            new_path_list.append(p)
    return ':'.join(new_path_list)

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

def get_eigen_type():
    return eigen_type
