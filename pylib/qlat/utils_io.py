import cqlat as c

import os
import pickle

from cqlat import qremove, qremove_info
from cqlat import qremove_all, qremove_all_info
from cqlat import qmkdir, qmkdir_info, qmkdir_sync_node
from cqlat import does_file_exist, does_file_exist_sync_node
from cqlat import does_file_exist_qar, does_file_exist_qar_sync_node
from cqlat import does_file_or_directory_exist_qar, does_file_or_directory_exist_qar_sync_node
from cqlat import is_directory, is_directory_sync_node
from cqlat import is_regular_file, is_regular_file_sync_node
from cqlat import qrename, qrename_info
from cqlat import qcat, qcat_sync_node
from cqlat import qcat_bytes, qcat_bytes_sync_node
from cqlat import qls, qls_sync_node
from cqlat import qls_all, qls_all_sync_node
from cqlat import qload_datatable, qload_datatable_sync_node
from cqlat import check_time_limit, check_stop
from cqlat import get_time_limit, get_default_budget

from qlat.timer import *
from qlat.mpi import *
from qlat.cache import *

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
def save_pickle_obj(obj, path):
    # only save from node 0
    # mk_file_dirs_info(path)
    if get_id_node() == 0:
        qtouch(path, pickle.dumps(obj))

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
    if not does_file_exist_sync_node(path):
        obj = func()
        save_pickle_obj(obj, path)
    else:
        obj = load_pickle_obj(path)
    return obj

@timer
def compute_crc32(path):
    return c.compute_crc32(path)

@timer
def check_all_files_crc32_info(path):
    return c.check_all_files_crc32_info(path)

def obtain_lock(path):
    mk_file_dirs_info(path)
    return c.obtain_lock(path)

def release_lock():
    return c.release_lock()

def qquit(msg):
    # clean python cache and then call c.qquit(msg) (which clear all the C++ level cache and then quit)
    clean_cache()
    return c.qquit(msg)

def get_remaining_time():
  return get_time_limit() - get_actual_total_time()
