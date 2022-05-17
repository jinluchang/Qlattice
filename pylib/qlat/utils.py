import cqlat as c

from qlat.field import *
from qlat.lat_io import *
from qlat.rng_state import *
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
    gc.unfreeze()
    clean_cache()
    gc.collect()
    gc.freeze()
    clear_all_caches()

@timer
def parallel_map(q_mp_proc, func, iterable):
    displayln_info(f"parallel_map(q_mp_proc={q_mp_proc})")
    if q_mp_proc == 0:
        return list(map(func, iterable))
    assert q_mp_proc >= 1
    global pool_function
    assert pool_function is None
    pool_function = func
    gc.collect()
    gc.freeze()
    with mp.Pool(q_mp_proc, process_initialization, []) as p:
        p.apply(malloc_stats)
        res = p.map(call_pool_function, iterable, chunksize = 1)
        p.apply(malloc_stats)
        p.apply(timer_display)
    gc.unfreeze()
    gc.collect()
    pool_function = None
    return res

@timer
def parallel_map_sum(q_mp_proc, func, iterable, *, sum_function = None, sum_initial = None, chunksize = 1):
    # iterable = [ i1, i2, ... ]
    # va1, vb1, ... = func(i1)
    # return [ sum([va1, va2, ...]), sum([vb1, vb2, ...]), ... ]
    displayln_info(f"parallel_map(q_mp_proc={q_mp_proc})")
    if sum_function is None:
        sum_function = lambda x: sum_list(x, sum_initial = sum_initial)
    if q_mp_proc == 0:
        return sum_function(map(func, iterable))
    assert q_mp_proc >= 1
    global pool_function
    assert pool_function is None
    pool_function = func
    gc.collect()
    gc.freeze()
    with mp.Pool(q_mp_proc, process_initialization, []) as p:
        p.apply(malloc_stats)
        res = p.imap(call_pool_function, iterable, chunksize = chunksize)
        ret = sum_function(res)
        p.apply(malloc_stats)
        p.apply(timer_display)
    gc.unfreeze()
    gc.collect()
    pool_function = None
    return ret

def sum_list(res, *, sum_initial = None):
    # res = [ [ va1, vb1, ... ], [ va2, vb2, ... ], ... ]
    # return [ sum([va1, va2, ...]), sum([vb1, vb2, ...]), ... ]
    if sum_initial is not None:
        ret = list(sum_initial)
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
    # ret = [ va, vb, ... ]
    # return [ glb_sum(va), glb_sum(vb), ... ]
    return [ glb_sum(r) for r in ret ]

def get_q_mp_proc():
    v = os.getenv("q_mp_proc")
    if v is None:
        return 0
    return int(v)

def lazy_call(f, *args, **kwargs):
    is_thunk = True
    ret = None
    def get():
        nonlocal ret, is_thunk
        if is_thunk:
            ret = f(*args, **kwargs)
            is_thunk = False
        return ret
    return get

def sqr(x):
    return x * x

def set_zero(x):
    x.set_zero()

def set_unit(x, coef = 1.0):
    x.set_unit(coef)

def qnorm(x):
    # qnorm(2) == 4
    if isinstance(x, np.ndarray):
        return np.abs(np.vdot(x, x))
    elif isinstance(x, (int,float)):
        return x * x
    elif isinstance(x, complex):
        return x.real * x.real + x.imag * x.imag
    else:
        return x.qnorm()
    assert False

def show(x):
    return x.show()

def unitarize(x):
    x.unitarize()

def random_permute(l, rs):
    # Do not change ``l''.
    # Return a new permuted list.
    assert isinstance(l, list)
    assert isinstance(rs, RngState)
    return c.random_permute(l, rs)

def get_all_caches_info():
    return c.get_all_caches_info()

def clear_all_caches():
    return c.clear_all_caches();

def malloc_stats():
    return c.malloc_stats()
