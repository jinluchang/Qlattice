import cqlat_utils as cu

import os
import pickle

from cqlat_utils import qremove
from cqlat_utils import qremove_all
from cqlat_utils import qmkdir, qmkdir_info
from cqlat_utils import qrename_info
from cqlat_utils import does_file_exist
from cqlat_utils import is_directory
from cqlat_utils import is_regular_file
from cqlat_utils import qrename, qrename_info
from cqlat_utils import qls
from cqlat_utils import qls_all

from qlat_utils.qar import *

@timer
def qmkdirs(path):
    os.makedirs(path, exist_ok=True)

@timer
def qmkdirs_info(path):
    if get_id_node() == 0:
        displayln(f"qmkdirs_info: '{path}'.")
        qmkdirs(path)

@timer
def mk_dirs(path):
    os.makedirs(path, exist_ok=True)

@timer
def mk_dirs_info(path):
    if get_id_node() == 0:
        displayln(f"mk_dirs_info: '{path}'.")
        mk_dirs(path)

@timer
def mk_file_dirs(fn):
    path = os.path.dirname(fn)
    if path != "":
        os.makedirs(path, exist_ok=True)

@timer
def mk_file_dirs_info(path):
    if get_id_node() == 0:
        displayln(f"mk_file_dirs_info: '{path}'.")
        mk_file_dirs(path)

def qtouch(path, content = None):
    # mk_file_dirs(path)
    if content is None:
        return cu.qtouch(path)
    else:
        return cu.qtouch(path, content)

def qtouch_info(path, content = None):
    # mk_file_dirs_info(path)
    if content is None:
        return cu.qtouch_info(path)
    else:
        return cu.qtouch_info(path, content)

def qappend(path, content = None):
    # mk_file_dirs(path)
    if content is None:
        return cu.qappend(path)
    else:
        return cu.qappend(path, content)

def qappend_info(path, content = None):
    # mk_file_dirs_info(path)
    if content is None:
        return cu.qappend_info(path)
    else:
        return cu.qappend_info(path, content)

@timer
def compute_crc32(path):
    return cu.compute_crc32(path)

def qload_datatable(path, is_par = False):
    return cu.qload_datatable(path, is_par)

@timer
def save_pickle_obj(obj, path):
    # only save from node 0
    # mk_file_dirs_info(path)
    if get_id_node() == 0:
        qtouch(path, pickle.dumps(obj))

@timer
def check_all_files_crc32_info(path):
    return cu.check_all_files_crc32_info(path)

@timer
def load_pickle_obj(path, default_value = None):
    # all the nodes read the same data
    if does_file_exist_qar_sync_node(path):
        obj = pickle.loads(qcat_bytes_sync_node(path))
        return obj
    else:
        return default_value

@timer
def pickle_cache_call(func, path):
    if not does_file_exist_qar_sync_node(path):
        obj = func()
        save_pickle_obj(obj, path)
    else:
        obj = load_pickle_obj(path)
    return obj

def qremove_info(path):
    if get_num_node() != 1:
        import cqlat as c
        return c.qremove_info(path)
    return qremove(path)

def qremove_all_info(path):
    if get_num_node() != 1:
        import cqlat as c
        return c.qremove_all_info(path)
    return qremove_all(path)

def qmkdir_sync_node(path):
    if get_num_node() != 1:
        import cqlat as c
        return c.qmkdir_sync_node(path)
    return qmkdir(path)

def does_file_exist_sync_node(path):
    if get_num_node() != 1:
        import cqlat as c
        return c.does_file_exist_sync_node(path)
    return does_file_exist(path)

def is_directory_sync_node(path):
    if get_num_node() != 1:
        import cqlat as c
        return c.is_directory_sync_node(path)
    return is_directory(path)

def is_regular_file_sync_node(path):
    if get_num_node() != 1:
        import cqlat as c
        return c.is_regular_file_sync_node(path)
    return is_regular_file(path)

def qls_sync_node(path):
    if get_num_node() != 1:
        import cqlat as c
        return c.qls_sync_node(path)
    return qls(path)

def qls_all_sync_node(path):
    if get_num_node() != 1:
        import cqlat as c
        return c.qls_all_sync_node(path)
    return qls_all(path)

def qload_datatable_sync_node(path, is_par = False):
    if get_num_node() != 1:
        import cqlat as c
        return c.qload_datatable_sync_node(path, is_par)
    return qload_datatable(path, is_par)
