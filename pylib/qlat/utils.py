import cqlat as c

from qlat.field import *
from qlat.lat_io import *
from qlat.rng_state import *
from qlat.coordinate import *

import numpy as np
import multiprocessing as mp

pool_function = None

def call_pool_function(*args, **kwargs):
    assert pool_function is not None
    return pool_function(*args, **kwargs)

def parallel_map(n_processes, func, iterable):
    if n_processes == 0:
        return list(map(func, iterable))
    global pool_function
    pool_function = func
    with mp.Pool(n_processes) as p:
        res = p.map(call_pool_function, iterable)
    pool_function = None
    return res

def get_n_processes():
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
