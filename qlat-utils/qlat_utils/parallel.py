from .cache import *
from .rng_state import *
from .utils import *

import sys
import os
import multiprocessing as mp
import gc

pool_function = None

def call_pool_function(*args, **kwargs):
    assert pool_function is not None
    return pool_function(*args, **kwargs)

@timer
def gc_collect():
    gc.collect()

@timer
def gc_freeze():
    gc.freeze()

@timer
def gc_unfreeze():
    gc.unfreeze()

def process_initialization():
    set_verbose_level(-2)
    timer_reset(0)
    # gc_unfreeze()
    # clean_cache()
    # gc_collect()
    # gc_freeze()
    # clear_all_caches()

q_num_mp_processes = None

def get_q_num_mp_processes():
    global q_num_mp_processes
    if q_num_mp_processes is not None:
        return q_num_mp_processes
    s = getenv("q_num_mp_processes", "q_num_threads", "OMP_NUM_THREADS", default="2")
    v = int(s)
    q_num_mp_processes = v
    return v

def set_q_num_mp_processes(v):
    global q_num_mp_processes
    q_num_mp_processes = v

q_verbose_parallel_map = None

def get_q_verbose_parallel_map():
    global q_verbose_parallel_map
    if q_verbose_parallel_map is not None:
        return q_verbose_parallel_map
    s = getenv("q_verbose_parallel_map", default="2")
    v = int(s)
    q_verbose_parallel_map = v
    return v

def set_q_verbose_parallel_map(v):
    global q_verbose_parallel_map
    q_verbose_parallel_map = v

@timer
def parallel_map(func, iterable,
        *,
        n_proc=None,
        chunksize=1,
        process_initialization=process_initialization,
        verbose=None):
    """
    iterable = [ i1, i2, ... ]
    v1 = func(i1)
    v2 = func(i2)
    ...
    return [ v1, v2, ... ]
    """
    if n_proc is None:
        n_proc = get_q_num_mp_processes()
    if verbose is None:
        verbose = get_q_verbose_parallel_map()
    if verbose > 0:
        displayln_info(f"parallel_map(n_proc={n_proc})")
    if n_proc == 0:
        res = map(func, iterable)
        ret = []
        for idx, v in enumerate(res):
            if verbose > 0 and idx % chunksize == 0:
                displayln_info(-2, f"parallel_map: idx={idx} done")
            ret.append(v)
        return ret
    assert n_proc >= 1
    global pool_function
    assert pool_function is None
    try:
        pool_function = func
        gc_collect()
        gc_freeze()
        with mp.Pool(n_proc, process_initialization, []) as p:
            if verbose > 0:
                p.apply(show_memory_usage)
            timer = Timer("parallel_map.p.imap")
            try:
                timer.start()
                res = p.imap(call_pool_function, iterable, chunksize=chunksize)
            finally:
                timer.stop()
            timer = Timer("parallel_map.mk_list")
            try:
                timer.start()
                ret = []
                for idx, v in enumerate(res):
                    if verbose > 0 and idx % chunksize == 0:
                        displayln_info(-2, f"parallel_map: idx={idx} done")
                    ret.append(v)
            finally:
                timer.stop()
            if verbose > 0:
                p.apply(show_memory_usage)
                if verbose > 1:
                    p.apply(timer_display)
    finally:
        gc_unfreeze()
        gc_collect()
        pool_function = None
    return ret

@timer
def parallel_map_sum(func, iterable,
        *,
        n_proc=None,
        sum_function=None,
        sum_start=None,
        chunksize=1,
        process_initialization=process_initialization,
        verbose=None):
    """
    iterable = [ i1, i2, ... ]
    v1 = func(i1)
    v2 = func(i2)
    ...
    return sum_function([ v1, v2, ... ])
    """
    if n_proc is None:
        n_proc = get_q_num_mp_processes()
    if verbose is None:
        verbose = get_q_verbose_parallel_map()
    if verbose > 0:
        displayln_info(f"parallel_map_sum(n_proc={n_proc})")
    if sum_function is None:
        sum_function = sum
    if n_proc == 0:
        if sum_start is None:
            return sum_function(map(func, iterable))
        else:
            return sum_function(map(func, iterable), sum_start)
    assert n_proc >= 1
    global pool_function
    assert pool_function is None
    try:
        pool_function = func
        gc_collect()
        gc_freeze()
        with mp.Pool(n_proc, process_initialization, []) as p:
            if verbose > 0:
                p.apply(show_memory_usage)
            timer = Timer("parallel_map_sum.p.imap")
            try:
                timer.start()
                res = p.imap(call_pool_function, iterable, chunksize=chunksize)
            finally:
                timer.stop()
            timer = Timer("parallel_map_sum.sum_function")
            try:
                timer.start()
                if sum_start is None:
                    ret = sum_function(res)
                else:
                    ret = sum_function(res, sum_start)
            finally:
                timer.stop()
            if verbose > 0:
                p.apply(show_memory_usage)
                if verbose > 1:
                    p.apply(timer_display)
    finally:
        gc_unfreeze()
        gc_collect()
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
