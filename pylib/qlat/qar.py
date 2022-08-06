import cqlat as c

from cqlat import get_qar_multi_vol_max_size

from qlat.timer import *

@timer
def qar_create(path_qar, path_folder, *, is_remove_folder_after = False):
    return c.qar_create(path_qar, path_folder, is_remove_folder_after)

@timer
def qar_create_info(path_qar, path_folder, *, is_remove_folder_after = False):
    return c.qar_create_info(path_qar, path_folder, is_remove_folder_after)

@timer
def qar_extract(path_qar, path_folder, *, is_remove_qar_after = False):
    return c.qar_extract(path_qar, path_folder, is_remove_qar_after)

@timer
def qar_extract_info(path_qar, path_folder, *, is_remove_qar_after = False):
    return c.qar_extract_info(path_qar, path_folder, is_remove_qar_after)

@timer
def qcopy_file(path_src, path_dst):
    return c.qcopy_file(path_src, path_dst)

@timer
def qcopy_file_info(path_src, path_dst):
    return c.qcopy_file_info(path_src, path_dst)

@timer
def list_qar(path_qar):
    return c.list_qar(path_qar)
