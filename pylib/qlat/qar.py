import cqlat as c

from qlat.timer import *

@timer
def qar_create(path_qar, path_folder, *, is_remove_folder_after = False):
    return c.qar_create(path_qar, path_folder, is_remove_folder_after)

@timer
def qar_create_info(path_qar, path_folder, *, is_remove_folder_after = False):
    return c.qar_create_info(path_qar, path_folder, is_remove_folder_after)

@timer
def qar_extract(path_qar, path_folder, *, is_remove_folder_after = False):
    return c.qar_extract(path_qar, path_folder, is_remove_folder_after)

@timer
def qar_extract_info(path_qar, path_folder, *, is_remove_qar_after = False):
    return c.qar_extract_info(path_qar, path_folder, is_remove_qar_after)

@timer
def list_qar(path_qar):
    return c.list_qar(path_qar)
