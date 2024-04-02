import os
import pickle
import hashlib
import functools

from .c import *

@timer
def qmkdirs(path):
    os.makedirs(path, exist_ok=True)
    remove_entry_directory_cache(path)

@timer
def qmkdirs_info(path):
    remove_entry_directory_cache(path)
    if get_id_node() == 0:
        displayln(f"qmkdirs_info: '{path}'.")
        qmkdirs(path)

@timer
def mk_dirs(path):
    remove_entry_directory_cache(path)
    os.makedirs(path, exist_ok=True)

@timer
def mk_dirs_info(path):
    remove_entry_directory_cache(path)
    if get_id_node() == 0:
        displayln(f"mk_dirs_info: '{path}'.")
        mk_dirs(path)

@timer
def mk_file_dirs(fn):
    remove_entry_directory_cache(fn)
    path = os.path.dirname(fn)
    if path != "":
        os.makedirs(path, exist_ok=True)

@timer
def mk_file_dirs_info(path):
    remove_entry_directory_cache(path)
    if get_id_node() == 0:
        displayln(f"mk_file_dirs_info: '{path}'.")
        mk_file_dirs(path)

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

def hash_sha256(s):
    """
    compute sha256 of str (or bytes) `s`.
    """
    m = hashlib.sha256()
    if isinstance(s, str):
        s = s.encode('utf8')
    m.update(s)
    return m.hexdigest()

def pickle_cache(path, is_sync_node=True):
    """
    `path` is the directory to cache results
    sha256 hash based on pickle.dumps of the input parameters
    """
    def dec(func):
        @functools.wraps(func)
        def f(*args, **kwargs):
            func_args = (func.__name__, args, kwargs)
            func_args_str = pickle.dumps(func_args)
            key = hash_sha256(func_args_str)
            fn = f"{path}/{key}.pickle"
            c_res = load_pickle_obj(fn, is_sync_node=is_sync_node)
            if c_res is not None:
                c_func_args, c_ret = c_res
                return c_ret
            ret = func(*args, **kwargs)
            res = (func_args, ret,)
            save_pickle_obj(res, fn, is_sync_node=is_sync_node)
            return ret
        return f
    return dec

class SetDisplayMethod:

    def __init__(self):
        set_display_method("py_stdout")
        # displayln_info(0, f"set_display_method('py_stdout')")

    def __del__(self):
        displayln_info(0, f"set_display_method()")
        set_display_method()

###
