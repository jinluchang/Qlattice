from qlat_utils.cache import *
from qlat_utils.rng_state import *
from qlat_utils.utils import *

import sys
import os
import multiprocessing as mp
import gc

pool_function = None

def call_pool_function(*args, **kwargs):
    assert pool_function is not None
    return pool_function(*args, **kwargs)

def process_initialization():
    verbose_level(-1)
    timer_reset(0)
    # gc.unfreeze()
    # clean_cache()
    # gc.collect()
    # gc.freeze()
    # clear_all_caches()

def get_q_num_mp_processes():
    global get_q_num_mp_processes
    s = getenv("q_num_mp_processes", "q_num_threads", "OMP_NUM_THREADS", default = "0")
    v = int(s)
    get_q_num_mp_processes = lambda : v
    return v

def set_q_num_mp_processes(v):
    global get_q_num_mp_processes
    get_q_num_mp_processes = lambda : v

def get_q_verbose_parallel_map():
    global get_q_verbose_parallel_map
    s = getenv("q_verbose_parallel_map", default = "1")
    v = int(s)
    get_q_verbose_parallel_map = lambda : v
    return v

def set_q_verbose_parallel_map(v):
    global get_q_verbose_parallel_map
    get_q_verbose_parallel_map = lambda : v

@timer
def parallel_map(func, iterable,
        _to_be_removed_ = None,
        *,
        n_proc = None,
        chunksize = 1,
        process_initialization = process_initialization,
        verbose = None):
    # iterable = [ i1, i2, ... ]
    # v1 = func(i1)
    # v2 = func(i2)
    # ...
    # return [ v1, v2, ... ]
    #####
    if isinstance(func, int):
        # assuming it is called with parallel_map(n_proc, func, iterable, ...)
        return parallel_map(iterable, _to_be_removed_,
                n_proc = func,
                sum_function = sum_function,
                chunksize = chunksize,
                process_initialization = process_initialization,
                verbose = verbose)
    assert _to_be_removed_ is None
    #####
    if n_proc is None:
        n_proc = get_q_num_mp_processes()
    if verbose is None:
        verbose = get_q_verbose_parallel_map()
    if verbose > 0:
        displayln_info(f"parallel_map(n_proc={n_proc})")
    if n_proc == 0:
        return list(map(func, iterable))
    assert n_proc >= 1
    global pool_function
    assert pool_function is None
    try:
        pool_function = func
        gc.collect()
        gc.freeze()
        with mp.Pool(n_proc, process_initialization, []) as p:
            if verbose > 0:
                p.apply(show_memory_usage)
            res = p.map(call_pool_function, iterable, chunksize = chunksize)
            if verbose > 0:
                p.apply(show_memory_usage)
                p.apply(timer_display)
    finally:
        gc.unfreeze()
        gc.collect()
        pool_function = None
    return res

@timer
def parallel_map_sum(func, iterable,
        _to_be_removed_ = None,
        *,
        n_proc = None,
        sum_function = None,
        sum_start = None,
        chunksize = 1,
        process_initialization = process_initialization,
        verbose = None):
    # iterable = [ i1, i2, ... ]
    # v1 = func(i1)
    # v2 = func(i2)
    # ...
    # return sum_function([ v1, v2, ... ])
    #####
    if isinstance(func, int):
        # assuming it is called with parallel_map_sum(n_proc, func, iterable, ...)
        return parallel_map_sum(iterable, _to_be_removed_,
                n_proc = func,
                sum_function = sum_function,
                sum_start = sum_start,
                chunksize = chunksize,
                process_initialization = process_initialization,
                verbose = verbose)
    assert _to_be_removed_ is None
    #####
    if n_proc is None:
        n_proc = get_q_num_mp_processes()
    if verbose is None:
        verbose = get_q_verbose_parallel_map()
    if verbose > 0:
        displayln_info(f"parallel_map_sum(n_proc={n_proc})")
    if sum_function is None:
        sum_function = sum
    if n_proc == 0:
        return sum_function(map(func, iterable))
    assert n_proc >= 1
    global pool_function
    assert pool_function is None
    try:
        pool_function = func
        gc.collect()
        gc.freeze()
        with mp.Pool(n_proc, process_initialization, []) as p:
            if verbose > 0:
                p.apply(show_memory_usage)
            res = p.imap(call_pool_function, iterable, chunksize = chunksize)
            if sum_start is None:
                ret = sum_function(res)
            else:
                ret = sum_function(res, sum_start)
            if verbose > 0:
                p.apply(show_memory_usage)
                p.apply(timer_display)
    finally:
        gc.unfreeze()
        gc.collect()
        pool_function = None
    return ret

def sum_list(res, start = None):
    # res = [ [ va1, vb1, ... ], [ va2, vb2, ... ], ... ]
    # return [ sum([va1, va2, ...]), sum([vb1, vb2, ...]), ... ]
    if start is not None:
        ret = list(start)
    else:
        ret = None
    for r in res:
        if ret is None:
            ret = list(r)
            continue
        for i, v in enumerate(r):
            ret[i] += v
    return ret
