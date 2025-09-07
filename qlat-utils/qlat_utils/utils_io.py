import os
import pickle
import hashlib
import functools
import numpy as np

from .c import *
from .lru_cache import *
from .json import json_dumps, json_loads

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
def save_json_obj(obj, path, *, indent=None, is_sync_node=True):
    """
    only save from node 0 when is_sync_node
    mk_file_dirs_info(path)
    """
    if not is_sync_node or get_id_node() == 0:
        qtouch(path, json_dumps(obj, indent=indent))

@timer
def load_json_obj(path, default_value=None, *, is_sync_node=True):
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
            obj = json_loads(qcat_bytes_sync_node(path))
        else:
            obj = json_loads(qcat_bytes(path))
        return obj
    else:
        return default_value

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
    all the nodes compute or load the same data if is_sync_node
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
    Compute sha256 of str (or bytes) `s`.
    Also support:
    Custom object that has `hash_sha256` attribute.
    Python `tuple` and `list`.
    Numpy `ndarray`.
    """
    if isinstance(s, (str, bytes,)):
        m = hashlib.sha256()
        if isinstance(s, str):
            s = s.encode('utf8')
        m.update(s)
        return m.hexdigest()
    elif hasattr(s, "hash_sha256"):
        return s.hash_sha256()
    elif isinstance(s, (tuple, list,)):
        if isinstance(s, tuple):
            m = hashlib.sha256(b'tuple:')
        elif isinstance(s, list):
            m = hashlib.sha256(b'list:')
        else:
            assert False
        for v in s:
            m.update(hash_sha256(v).encode('utf8'))
        return m.hexdigest()
    elif isinstance(s, np.ndarray):
        m = hashlib.sha256(b'np.ndarray:')
        m.update(repr(s.tolist()).encode('utf8'))
        return m.hexdigest()
    elif hasattr(s, "__repr__"):
        m = hashlib.sha256(b'repr:')
        m.update(repr(s).encode('utf8'))
        return m.hexdigest()
    else:
        assert False

def pickle_cache(path, is_sync_node=True):
    """
    `path` is the directory to cache results
    sha256 hash based on pickle.dumps of the input parameters
    """
    def dec(func):
        @functools.wraps(func)
        def f(*args, **kwargs):
            func_args = (func.__qualname__, args, kwargs)
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

def cache_call(
        *,
        maxsize=128,
        get_state=None,
        is_hash_args=True,
        path=None,
        is_sync_node=True,
        cache=None,
        ):
    """
    get_state() => object to be used as extra key of cache
    #
    `maxsize` can be `0`, where the LRUCache is effectively turned off.
    #
    if is_hash_args:
        Pickle all the keys and use hash as the key (default)
    else:
        Use `(func.__qualname__, args, state)` directly as `key`.
        Note that `kwargs` has to be empty in this case.
    if path is None:
        Only cache using LRUCache (default)
    else:
        Also cache the results in f"{path}/{key}.pickle".
        if is_sync_node:
            Only read/write to f"{path}/{key}.pickle" from process 0 (broadcast to all nodes)
        else:
            All the processes independently do the calculation and read/write to f"{path}/{key}.pickle"
            Use this if (1) `path` or `key` is different for different processes;
                     or (2) this function is only called from a certain process.
    if cache is None:
        cache = LRUCache(maxsize)
    else:
        The input `cache` will be used. This cache may be shared for other purpose
    #
    # Usage example:
    #
    @cache_call(maxsize=128, get_state=q.get_jk_state)
    def func(x):
        return x**2
    #
    block_size = 10
    block_size_dict = { "48I": 10, }
    @cache_call(maxsize=128, get_state=lambda: (block_size, block_size_dict,), is_hash_args=True)
    def func(x):
        return x**2
    #
    func.cache # returns the underlying lru_cache object.
    #
    func.clear() # clears the underlying lru_cache cache object.
    #
    """
    if cache is None:
        cache = LRUCache(maxsize)
    def dec(f):
        @functools.wraps(f)
        def func(*args, **kwargs):
            if get_state is None:
                state = None
            else:
                state = get_state()
            func_args = (f.__qualname__, args, kwargs, state,)
            if is_hash_args:
                func_args_str = pickle.dumps(func_args)
                key = hash_sha256(func_args_str)
            else:
                assert kwargs == dict()
                key = (f.__qualname__, args, state,)
            if key in cache:
                c_res = cache[key]
                c_func_args, c_ret = c_res
                return c_ret
            if path is not None:
                fn = f"{path}/{key}.pickle"
                c_res = load_pickle_obj(fn, is_sync_node=is_sync_node)
                if c_res is not None:
                    cache[key] = c_res
                    c_func_args, c_ret = c_res
                    return c_ret
            ret = f(*args, **kwargs)
            res = (func_args, ret,)
            cache[key] = res
            if path is not None:
                save_pickle_obj(res, fn, is_sync_node=is_sync_node)
            return ret
        func.cache = cache
        func.clear = lambda: cache.clear()
        return func
    return dec

class SetDisplayMethod:

    def __init__(self):
        set_display_method("py_stdout")
        # displayln_info(0, f"set_display_method('py_stdout')")

    def __del__(self):
        displayln_info(0, f"set_display_method()")
        set_display_method()

###
