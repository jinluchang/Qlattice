import qlat_utils.c as cu

from qlat_utils.cache import *
from qlat_utils.rng_state import *

import math
import sys
import os

def getenv(*names, default = None):
    assert len(names) > 0
    for name in names:
        val = os.getenv(name)
        if val is not None:
            displayln_info(0, f"{name}='{val}'")
            return val
    val = default
    displayln_info(0, f"{names[0]}='{val}' (default)")
    return val

def get_arg(option, default = None):
    argv = sys.argv
    i_max = len(argv) - 1
    for i, arg in enumerate(argv):
        if arg == option:
            if i == i_max:
                return ""
            else:
                return argv[i + 1]
    return default

def random_permute(l, rs):
    # Do not change ``l''.
    # Return a new permuted list.
    assert isinstance(l, list)
    assert isinstance(rs, RngState)
    return cu.random_permute(l, rs)

def show_memory_usage():
    try:
        import psutil
        rss = psutil.Process().memory_info().rss / (1024 * 1024 * 1024)
        displayln_info(f"show_memory_usage: rss = {rss:.6f} GB")
        # displayln_info_malloc_stats()
    except:
        displayln_info(f"show_memory_usage: no psutil.")

def displayln_malloc_stats():
    return cu.displayln_malloc_stats()

def displayln_info_malloc_stats():
    if get_id_node() == 0:
        return displayln_malloc_stats()

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

def show(x):
    return x.show()

def unitarize(x):
    x.unitarize()

def rel_mod(x, size):
    x = x % size
    assert x >= 0
    if 2 * x >= size:
        return x - size
    else:
        return x

def rel_mod_sym(x, size):
    x = x % size
    assert x >= 0
    if 2 * x > size:
        return x - size
    elif 2 * x < size:
        return x
    else:
        assert 2 * x == size
        return 0

def c_sqr(x):
    return sum([ sqr(v) for v in x ])

def c_rel_mod(x, size):
    l = len(size)
    assert l == len(x)
    return [ rel_mod(x[i], size[i]) for i in range(l) ]

def c_rel_mod_sqr(x, size):
    l = len(size)
    assert l == len(x)
    return sum([ sqr(rel_mod(x[i], size[i])) for i in range(l) ])

def phat_sqr(q, size):
    l = len(size)
    assert l == len(q)
    return 4 * sum([ sqr(math.sin(math.pi * (q[i] % size[i]) / size[i])) for i in range(l) ])
