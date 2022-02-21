import cqlat as c

import os
import pickle

from cqlat import qremove, qremove_info
from cqlat import qremove_all, qremove_all_info
from cqlat import qmkdir, qmkdir_info, qmkdir_sync_node
from cqlat import does_file_exist, does_file_exist_sync_node
from cqlat import is_directory, is_directory_sync_node
from cqlat import qrename, qrename_info
from cqlat import qcat, qcat_sync_node
from cqlat import qls, qls_sync_node
from cqlat import qload_datatable, qload_datatable_sync_node
from cqlat import check_time_limit, check_stop

from qlat.mpi import *

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
    mk_file_dirs(path)
    if content is None:
        return c.qtouch(path)
    else:
        return c.qtouch(path, content)

def qtouch_info(path, content = None):
    mk_file_dirs_info(path)
    if content is None:
        return c.qtouch_info(path)
    else:
        return c.qtouch_info(path, content)

def qappend(path, content = None):
    mk_file_dirs(path)
    if content is None:
        return c.qappend(path)
    else:
        return c.qappend(path, content)

def qappend_info(path, content = None):
    mk_file_dirs_info(path)
    if content is None:
        return c.qappend_info(path)
    else:
        return c.qappend_info(path, content)

def obtain_lock(path):
    mk_file_dirs_info(path)
    return c.obtain_lock(path)

def release_lock():
    return c.release_lock()

@timer
def save_pickle_obj(obj, path):
    # only save from node 0
    mk_file_dirs_info(path)
    if get_id_node() == 0:
        with open(path + ".partial", "wb") as f:
            pickle.dump(obj, f)
        qrename(path + ".partial", path)

@timer
def load_pickle_obj(path, default_value = None):
    # all the nodes read the same data
    if does_file_exist_sync_node(path):
        with open(path, "rb") as f:
            obj = pickle.load(f)
        return obj
    else:
        return default_value

@timer
def pickle_cache_call(func, path):
    if not does_file_exist_sync_node(path):
        obj = func()
        save_pickle_obj(obj, path)
    else:
        obj = load_pickle_obj(path)
    return obj
