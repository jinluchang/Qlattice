import qlat.cqlat as c

from qlat_utils import *

from qlat.field import *
from qlat.coordinate import *

import numpy as np
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

@timer
def parallel_map(q_mp_proc, func, iterable,
        *,
        chunksize = 1,
        process_initialization = process_initialization,
        is_verbose = False):
    if is_verbose:
        displayln_info(f"parallel_map(q_mp_proc={q_mp_proc})")
    if q_mp_proc == 0:
        return list(map(func, iterable))
    assert q_mp_proc >= 1
    global pool_function
    assert pool_function is None
    try:
        pool_function = func
        gc.collect()
        gc.freeze()
        with mp.Pool(q_mp_proc, process_initialization, []) as p:
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
def parallel_map_sum(q_mp_proc, func, iterable,
        *,
        sum_function = None,
        sum_start = 0,
        chunksize = 1,
        process_initialization = process_initialization,
        is_verbose = False):
    # iterable = [ i1, i2, ... ]
    # va1, vb1, ... = func(i1)
    # return [ sum([va1, va2, ...]), sum([vb1, vb2, ...]), ... ]
    if is_verbose:
        displayln_info(f"parallel_map_sum(q_mp_proc={q_mp_proc})")
    if sum_function is None:
        sum_function = sum
    if q_mp_proc == 0:
        return sum_function(map(func, iterable))
    assert q_mp_proc >= 1
    global pool_function
    assert pool_function is None
    try:
        pool_function = func
        gc.collect()
        gc.freeze()
        with mp.Pool(q_mp_proc, process_initialization, []) as p:
            if is_verbose:
                p.apply(show_memory_usage)
            res = p.imap(call_pool_function, iterable, chunksize = chunksize)
            if sum_start == 0:
                ret = sum_function(res)
            else:
                # ret = sum_function(res, start = sum_start)
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

def glb_sum_list(ret):
    # deprecated (use glb_sum instead)
    # ret = [ va, vb, ... ]
    # return [ glb_sum(va), glb_sum(vb), ... ]
    return glb_sum(ret)

def get_q_mp_proc():
    v = os.getenv("q_mp_proc")
    if v is None:
        return 0
    return int(v)
