import cqlat_utils as cu

from cqlat_utils import get_qar_multi_vol_max_size
from cqlat_utils import does_file_exist_qar
from cqlat_utils import does_file_or_directory_exist_qar
from cqlat_utils import qcat
from cqlat_utils import qcat_bytes

from qlat_utils.timer import *

@timer
def qar_create(path_qar, path_folder, *, is_remove_folder_after = False):
    return cu.qar_create(path_qar, path_folder, is_remove_folder_after)

@timer
def qar_extract(path_qar, path_folder, *, is_remove_qar_after = False):
    return cu.qar_extract(path_qar, path_folder, is_remove_qar_after)

@timer
def qcopy_file(path_src, path_dst):
    return cu.qcopy_file(path_src, path_dst)

@timer
def list_qar(path_qar):
    return cu.list_qar(path_qar)

def does_file_exist_qar_sync_node(path):
    if get_num_node() != 1:
        import cqlat as c
        return c.does_file_exist_qar_sync_node(path)
    return does_file_exist_qar(path)

def does_file_or_directory_exist_qar_sync_node(path):
    if get_num_node() != 1:
        import cqlat as c
        return c.does_file_or_directory_exist_qar_sync_node(path)
    return does_file_or_directory_exist_qar(path)

def qcat_sync_node(path):
    if get_num_node() != 1:
        import cqlat as c
        return c.qcat_sync_node(path)
    return qcat(path)

def qcat_bytes_sync_node(path):
    if get_num_node() != 1:
        import cqlat as c
        return c.qcat_bytes_sync_node(path)
    return qcat_bytes(path)

@timer
def qar_create_info(path_qar, path_folder, *, is_remove_folder_after = False):
    if get_num_node() != 1:
        import cqlat as c
        return c.qar_create_info(path_qar, path_folder, is_remove_folder_after)
    return qar_create(path_qar, path_folder, is_remove_folder_after = is_remove_folder_after)

@timer
def qar_extract_info(path_qar, path_folder, *, is_remove_qar_after = False):
    if get_num_node() != 1:
        import cqlat as c
        return c.qar_extract_info(path_qar, path_folder, is_remove_qar_after)
    return qar_extract(path_qar, path_folder, is_remove_qar_after = is_remove_qar_after)

@timer
def qcopy_file_info(path_src, path_dst):
    if get_num_node() != 1:
        import cqlat as c
        return c.qcopy_file_info(path_src, path_dst)
    return qcopy_file(path_src, path_dst)
