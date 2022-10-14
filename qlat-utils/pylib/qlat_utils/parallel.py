import cqlat_utils as cu

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
    v = os.getenv("q_num_mp_processes")
    if v is None:
        v = os.getenv("q_num_threads")
    if v is None:
        v = os.getenv("OMP_NUM_THREADS")
    if v is None:
        v = 2
    displayln_info(0, f"q_num_mp_processes='{v}'")
    return int(v)

q_num_mp_processes = get_q_num_mp_processes()

def get_q_is_verbose_parallel_map():
    v = os.getenv("q_is_verbose_parallel_map")
    if v is None:
        v = True
    displayln_info(0, f"q_is_verbose_parallel_map='{v}'")
    return bool(v)

q_is_verbose_parallel_map = get_q_is_verbose_parallel_map()

@timer
def parallel_map(func, iterable,
        _to_be_removed_ = None,
        *,
        n_proc = None,
        chunksize = 1,
        process_initialization = process_initialization,
        is_verbose = None):
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
                is_verbose = is_verbose)
    assert _to_be_removed_ is None
    #####
    if n_proc is None:
        n_proc = q_num_mp_processes
    if is_verbose is None:
        is_verbose = q_is_verbose_parallel_map
    if is_verbose:
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
            if is_verbose:
                p.apply(show_memory_usage)
            res = p.map(call_pool_function, iterable, chunksize = chunksize)
            if is_verbose:
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
        is_verbose = False):
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
                is_verbose = is_verbose)
    assert _to_be_removed_ is None
    #####
    if n_proc is None:
        n_proc = q_num_mp_processes
    if is_verbose is None:
        is_verbose = q_is_verbose_parallel_map
    if is_verbose:
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
            if is_verbose:
                p.apply(show_memory_usage)
            res = p.imap(call_pool_function, iterable, chunksize = chunksize)
            if sum_start is None:
                ret = sum_function(res)
            else:
                ret = sum_function(res, sum_start)
            if is_verbose:
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
