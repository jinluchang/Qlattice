import qlat_utils.c as c

import os
import pickle

from .c import qremove, qremove_all
from .c import qmkdir, qmkdir_info
from .c import is_directory
from .c import is_regular_file
from .c import does_file_exist
from .c import is_directory_cache
from .c import is_regular_file_cache
from .c import does_file_exist_cache
from .c import clear_is_directory_cache
from .c import qrename, qrename_info
from .c import qls
from .c import qls_all
from .c import compute_crc32
from .c import qload_datatable
from .c import check_all_files_crc32_info

from .qar_utils import *

from .c import *

@timer
def qmkdirs(path):
    os.makedirs(path, exist_ok=True)
    clear_is_directory_cache()

@timer
def qmkdirs_info(path):
    clear_is_directory_cache()
    if get_id_node() == 0:
        displayln(f"qmkdirs_info: '{path}'.")
        qmkdirs(path)

@timer
def mk_dirs(path):
    clear_is_directory_cache()
    os.makedirs(path, exist_ok=True)

@timer
def mk_dirs_info(path):
    clear_is_directory_cache()
    if get_id_node() == 0:
        displayln(f"mk_dirs_info: '{path}'.")
        mk_dirs(path)

@timer
def mk_file_dirs(fn):
    clear_is_directory_cache()
    path = os.path.dirname(fn)
    if path != "":
        os.makedirs(path, exist_ok=True)

@timer
def mk_file_dirs_info(path):
    clear_is_directory_cache()
    if get_id_node() == 0:
        displayln(f"mk_file_dirs_info: '{path}'.")
        mk_file_dirs(path)

def qtouch(path, content = None):
    # mk_file_dirs(path)
    if content is None:
        return c.qtouch(path)
    else:
        return c.qtouch(path, content)

def qtouch_info(path, content = None):
    # mk_file_dirs_info(path)
    if content is None:
        return c.qtouch_info(path)
    else:
        return c.qtouch_info(path, content)

def qappend(path, content = None):
    # mk_file_dirs(path)
    if content is None:
        return c.qappend(path)
    else:
        return c.qappend(path, content)

def qappend_info(path, content = None):
    # mk_file_dirs_info(path)
    if content is None:
        return c.qappend_info(path)
    else:
        return c.qappend_info(path, content)

@timer
def save_pickle_obj(obj, path, *, is_sync_node=True):
    """
    only save from node 0 when is_sync_node
    mk_file_dirs_info(path)
    """
    if not is_sync_node or get_id_node() == 0:
        qtouch(path, pickle.dumps(obj))

@timer
def load_pickle_obj(path, default_value=None, *, is_sync_node=True):
    """
    all the nodes read the same data
    if is_sync_node:
        one node read and broadcast to other nodes
    else:
        all nodes individually read the data
    """
    if is_sync_node:
        b = does_file_exist_qar_sync_node(path)
    else:
        b = does_file_exist_qar(path)
    if b:
        if is_sync_node:
            obj = pickle.loads(qcat_bytes_sync_node(path))
        else:
            obj = pickle.loads(qcat_bytes(path))
        return obj
    else:
        return default_value

@timer
def pickle_cache_call(func, path, *, is_sync_node=True):
    """
    all the nodes compute or load the same data
    """
    if is_sync_node:
        b = does_file_exist_qar_sync_node(path)
    else:
        b = does_file_exist_qar(path)
    if not b:
        obj = func()
        save_pickle_obj(obj, path, is_sync_node=is_sync_node)
    else:
        obj = load_pickle_obj(path, is_sync_node=is_sync_node)
    return obj

def qremove_info(path):
    if get_num_node() != 1:
        import qlat.c as c
        return c.qremove_info(path)
    return qremove(path)

def qremove_all_info(path):
    if get_num_node() != 1:
        import qlat.c as c
        return c.qremove_all_info(path)
    return qremove_all(path)

def qmkdir_sync_node(path):
    clear_is_directory_cache()
    if get_num_node() != 1:
        import qlat.c as c
        return c.qmkdir_sync_node(path)
    return qmkdir(path)

def does_file_exist_sync_node(path):
    if get_num_node() != 1:
        import qlat.c as c
        return c.does_file_exist_sync_node(path)
    return does_file_exist(path)

def is_directory_sync_node(path):
    if get_num_node() != 1:
        import qlat.c as c
        return c.is_directory_sync_node(path)
    return is_directory(path)

def is_regular_file_sync_node(path):
    if get_num_node() != 1:
        import qlat.c as c
        return c.is_regular_file_sync_node(path)
    return is_regular_file(path)

def qls_sync_node(path):
    if get_num_node() != 1:
        import qlat.c as c
        return c.qls_sync_node(path)
    return qls(path)

def qls_all_sync_node(path):
    if get_num_node() != 1:
        import qlat.c as c
        return c.qls_all_sync_node(path)
    return qls_all(path)

def qload_datatable_sync_node(path, is_par = False):
    if get_num_node() != 1:
        import qlat.c as c
        return c.qload_datatable_sync_node(path, is_par)
    return qload_datatable(path, is_par)
