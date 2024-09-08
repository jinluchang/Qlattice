# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from . cimport everything as cc

### -------------------------------------------------------------------

def basename(const cc.std_string& fn):
    return cc.basename(fn)

def dirname(const cc.std_string& fn):
    return cc.dirname(fn)

def all_dirname_vec(const cc.std_string& fn):
    return cc.all_dirname_vec(fn)

def remove_trailing_slashes(const cc.std_string& fn):
    return cc.remove_trailing_slashes(fn)

def qls(const cc.std_string& path, const cc.Bool is_sort=True):
    cdef list l = cc.qls(path, is_sort)
    return [ <str>v for v in l ]

def qls_all(const cc.std_string& path,
        const cc.Bool is_folder_before_files=False,
        const cc.Bool is_sort=True):
    cdef list l = cc.qls_all(path, is_folder_before_files, is_sort)
    return [ <str>v for v in l ]

def does_file_exist(const cc.std_string& path):
    return cc.does_file_exist(path)

def is_directory(const cc.std_string& path):
    return cc.is_directory(path)

def is_regular_file(const cc.std_string& path):
    return cc.is_regular_file(path)

def qmkdir(const cc.std_string& path):
    return cc.qmkdir(path)

def qmkdir_p(const cc.std_string& path):
    return cc.qmkdir_p(path)

def qrename(const cc.std_string& old_path, const cc.std_string& new_path):
    return cc.qrename(old_path, new_path)

def qremove(const cc.std_string& path):
    return cc.qremove(path)

def qremove_all(const cc.std_string& path):
    return cc.qremove_all(path)

def qtruncate(const cc.std_string& path, const cc.Long offset=0):
    return cc.qtruncate(path, offset)

### -------------------------------------------------------------------

def clear_is_directory_cache():
    return cc.clear_is_directory_cache()

def remove_entry_directory_cache(const cc.std_string& path):
    return cc.remove_entry_directory_cache(path)

def is_directory_cache(const cc.std_string& path):
    return cc.is_directory_cache(path)

def is_regular_file_cache(const cc.std_string& path):
    return cc.is_regular_file_cache(path)

def does_file_exist_cache(const cc.std_string& path):
    return cc.does_file_exist_cache(path)

### -------------------------------------------------------------------

def qmkdir_info(const cc.std_string& path):
    return cc.qmkdir_info(path)

def qmkdir_p_info(const cc.std_string& path):
    return cc.qmkdir_p_info(path)

def qrename_info(const cc.std_string& old_path, const cc.std_string& new_path):
    return cc.qrename_info(old_path, new_path)

def qremove_info(const cc.std_string& path):
    return cc.qremove_info(path)

def qremove_all_info(const cc.std_string& path):
    return cc.qremove_all_info(path)

### -------------------------------------------------------------------

def qls_sync_node(const cc.std_string& path, const cc.Bool is_sort=True):
    cdef list l = cc.qls_sync_node(path, is_sort)
    return [ <str>v for v in l ]

def qls_all_sync_node(const cc.std_string& path,
        const cc.Bool is_folder_before_files=False,
        const cc.Bool is_sort=True):
    cdef list l = cc.qls_all_sync_node(path, is_folder_before_files, is_sort)
    return [ <str>v for v in l ]

def does_file_exist_sync_node(path):
    return cc.does_file_exist_sync_node(path)

def is_directory_sync_node(path):
    return cc.is_directory_sync_node(path)

def is_regular_file_sync_node(path):
    return cc.is_regular_file_sync_node(path)

def qmkdir_sync_node(path):
    return cc.qmkdir_sync_node(path)

def qmkdir_p_sync_node(path):
    return cc.qmkdir_p_sync_node(path)

def qremove_sync_node(path):
    return cc.qremove_sync_node(path)

def qremove_all_sync_node(path):
    return cc.qremove_all_sync_node(path)

### -------------------------------------------------------------------

def displayln_malloc_stats():
    cc.displayln_malloc_stats()

### -------------------------------------------------------------------

def get_all_caches_info():
    cdef list l = cc.get_all_caches_info()
    return [ str(v) for v in l ]

def clear_all_caches():
    cc.clear_all_caches()

def clear_mem_cache():
    cc.clear_mem_cache()

### -------------------------------------------------------------------

def get_eigen_type():
    return cc.get_eigen_type()
