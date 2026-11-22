# cython: binding=True, embedsignature=True

"""
Module ``qlat_utils.cutils``
=============================\n
Low-level C++ filesystem utilities exposed to Python: path manipulation, directory listing, file/directory creation and removal, cached existence checks, MPI-synchronized variants, and cache management.\n
Documentation: ``docs/qlat-utils/qlat_cutils.md``\n
.. note:: Update the documentation when updating this source file.
"""

from . cimport everything as cc

### -------------------------------------------------------------------

def basename(fn):
    cdef bytes fn_b = cc.as_bytes(fn)
    return cc.as_str_std(cc.basename(fn_b))

def dirname(fn):
    cdef bytes fn_b = cc.as_bytes(fn)
    return cc.as_str_std(cc.dirname(fn_b))

def all_dirname_vec(fn):
    cdef bytes fn_b = cc.as_bytes(fn)
    return [cc.as_str_std(v) for v in cc.all_dirname_vec(fn_b)]

def remove_trailing_slashes(fn):
    cdef bytes fn_b = cc.as_bytes(fn)
    return cc.as_str_std(cc.remove_trailing_slashes(fn_b))

def qls(path, const cc.Bool is_sort=True):
    cdef bytes path_b = cc.as_bytes(path)
    return [cc.as_str_std(v) for v in cc.qls(path_b, is_sort)]

def qls_all(path,
        const cc.Bool is_folder_before_files=False,
        const cc.Bool is_sort=True):
    cdef bytes path_b = cc.as_bytes(path)
    return [cc.as_str_std(v) for v in cc.qls_all(path_b, is_folder_before_files, is_sort)]

def does_file_exist(path):
    cdef bytes path_b = cc.as_bytes(path)
    return cc.does_file_exist(path_b)

def is_directory(path):
    cdef bytes path_b = cc.as_bytes(path)
    return cc.is_directory(path_b)

def is_regular_file(path):
    cdef bytes path_b = cc.as_bytes(path)
    return cc.is_regular_file(path_b)

def qmkdir(path):
    cdef bytes path_b = cc.as_bytes(path)
    return cc.qmkdir(path_b)

def qmkdir_p(path):
    cdef bytes path_b = cc.as_bytes(path)
    return cc.qmkdir_p(path_b)

def qrename(old_path, new_path):
    cdef bytes old_path_b = cc.as_bytes(old_path)
    cdef bytes new_path_b = cc.as_bytes(new_path)
    return cc.qrename(old_path_b, new_path_b)

def qremove(path):
    cdef bytes path_b = cc.as_bytes(path)
    return cc.qremove(path_b)

def qremove_all(path):
    cdef bytes path_b = cc.as_bytes(path)
    return cc.qremove_all(path_b)

def qtruncate(path, const cc.Long offset=0):
    cdef bytes path_b = cc.as_bytes(path)
    return cc.qtruncate(path_b, offset)

### -------------------------------------------------------------------

def clear_is_directory_cache():
    return cc.clear_is_directory_cache()

def remove_entry_directory_cache(path):
    cdef bytes path_b = cc.as_bytes(path)
    return cc.remove_entry_directory_cache(path_b)

def is_directory_cache(path):
    cdef bytes path_b = cc.as_bytes(path)
    return cc.is_directory_cache(path_b)

def is_regular_file_cache(path):
    cdef bytes path_b = cc.as_bytes(path)
    return cc.is_regular_file_cache(path_b)

def does_file_exist_cache(path):
    cdef bytes path_b = cc.as_bytes(path)
    return cc.does_file_exist_cache(path_b)

### -------------------------------------------------------------------

def qmkdir_info(path):
    cdef bytes path_b = cc.as_bytes(path)
    return cc.qmkdir_info(path_b)

def qmkdir_p_info(path):
    cdef bytes path_b = cc.as_bytes(path)
    return cc.qmkdir_p_info(path_b)

def qrename_info(old_path, new_path):
    cdef bytes old_path_b = cc.as_bytes(old_path)
    cdef bytes new_path_b = cc.as_bytes(new_path)
    return cc.qrename_info(old_path_b, new_path_b)

def qremove_info(path):
    cdef bytes path_b = cc.as_bytes(path)
    return cc.qremove_info(path_b)

def qremove_all_info(path):
    cdef bytes path_b = cc.as_bytes(path)
    return cc.qremove_all_info(path_b)

### -------------------------------------------------------------------

def qls_sync_node(path, const cc.Bool is_sort=True):
    cdef bytes path_b = cc.as_bytes(path)
    return [cc.as_str_std(v) for v in cc.qls_sync_node(path_b, is_sort)]

def qls_all_sync_node(path,
        const cc.Bool is_folder_before_files=False,
        const cc.Bool is_sort=True):
    cdef bytes path_b = cc.as_bytes(path)
    return [cc.as_str_std(v) for v in cc.qls_all_sync_node(path_b, is_folder_before_files, is_sort)]

def does_file_exist_sync_node(path):
    cdef bytes path_b = cc.as_bytes(path)
    return cc.does_file_exist_sync_node(path_b)

def is_directory_sync_node(path):
    cdef bytes path_b = cc.as_bytes(path)
    return cc.is_directory_sync_node(path_b)

def is_regular_file_sync_node(path):
    cdef bytes path_b = cc.as_bytes(path)
    return cc.is_regular_file_sync_node(path_b)

def qmkdir_sync_node(path):
    cdef bytes path_b = cc.as_bytes(path)
    return cc.qmkdir_sync_node(path_b)

def qmkdir_p_sync_node(path):
    cdef bytes path_b = cc.as_bytes(path)
    return cc.qmkdir_p_sync_node(path_b)

def qremove_sync_node(path):
    cdef bytes path_b = cc.as_bytes(path)
    return cc.qremove_sync_node(path_b)

def qremove_all_sync_node(path):
    cdef bytes path_b = cc.as_bytes(path)
    return cc.qremove_all_sync_node(path_b)

### -------------------------------------------------------------------

def displayln_malloc_stats():
    cc.displayln_malloc_stats()

### -------------------------------------------------------------------

def get_all_caches_info():
    return [cc.as_str_std(v) for v in cc.get_all_caches_info()]

def clear_all_caches():
    cc.clear_all_caches()

def clear_mem_cache():
    cc.clear_mem_cache()

### -------------------------------------------------------------------

def get_eigen_type():
    return cc.as_str_std(cc.get_eigen_type())

### -------------------------------------------------------------------
