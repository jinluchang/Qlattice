# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from . cimport everything as cc

from .timer import timer

### ----------------------------------------------------------

def show_all_qfile():
    return cc.show_all_qfile()

def get_qar_multi_vol_max_size():
    """
    Parameter controls the size of a single `qar` file in number of bytes. Note, `qar` never splits a single file into multiple `qar` volume.
    """
    return cc.get_qar_multi_vol_max_size()

def set_qar_multi_vol_max_size(size=None):
    """
    Parameter controls the size of a single `qar` file in number of bytes. Note, `qar` never splits a single file into multiple `qar` volume.
    """
    if size is None:
        size = cc.get_qar_multi_vol_max_size_default()
    assert isinstance(size, int)
    cdef cc.Long* p_size = &cc.get_qar_multi_vol_max_size()
    p_size[0] = size
    assert cc.get_qar_multi_vol_max_size() == size

def does_regular_file_exist_qar(const cc.std_string& path):
    return cc.does_regular_file_exist_qar(path)

def does_file_exist_qar(const cc.std_string& path):
    return cc.does_file_exist_qar(path)

def qcat(const cc.std_string& path):
    """Return contents of file as `str`"""
    return <str>cc.qcat(path)

def qcat_bytes(const cc.std_string& path):
    """Return contents of file as `bytes`"""
    return <bytes>cc.qcat(path)

@timer
def qar_build_index(const cc.std_string& path_qar):
    """
    create "path_qar.idx" file
    """
    cc.qar_build_index(path_qar)

@timer
def qar_create(const cc.std_string& path_qar, const cc.std_string& path_folder,
               *, const cc.bool is_remove_folder_after=False):
    return cc.qar_create(path_qar, path_folder, is_remove_folder_after)

@timer
def qar_extract(const cc.std_string& path_qar, const cc.std_string& path_folder,
                *, const cc.bool is_remove_qar_after=False):
    return cc.qar_extract(path_qar, path_folder, is_remove_qar_after)

@timer
def qcopy_file(const cc.std_string& path_src, const cc.std_string& path_dst):
    return cc.qcopy_file(path_src, path_dst)

@timer
def list_qar(const cc.std_string& path_qar):
    cdef list l = cc.list_qar(path_qar)
    return [ str(fn) for fn in l ]

### ----------------------------------------------------------

def does_regular_file_exist_qar_sync_node(const cc.std_string& path):
    return cc.does_regular_file_exist_qar_sync_node(path)

def does_file_exist_qar_sync_node(const cc.std_string& path):
    return cc.does_file_exist_qar_sync_node(path)

### ----------------------------------------------------------

@timer
def qar_build_index_info(const cc.std_string& path_qar):
    """
    create "path_qar.idx" file (only on node 0)
    """
    cc.qar_build_index_info(path_qar)

@timer
def qar_create_info(const cc.std_string& path_qar, const cc.std_string& path_folder,
               *, const cc.bool is_remove_folder_after=False):
    return cc.qar_create_info(path_qar, path_folder, is_remove_folder_after)

@timer
def qar_extract_info(const cc.std_string& path_qar, const cc.std_string& path_folder,
                *, const cc.bool is_remove_qar_after=False):
    return cc.qar_extract_info(path_qar, path_folder, is_remove_qar_after)

@timer
def qcopy_file_info(const cc.std_string& path_src, const cc.std_string& path_dst):
    return cc.qcopy_file_info(path_src, path_dst)

### ----------------------------------------------------------
