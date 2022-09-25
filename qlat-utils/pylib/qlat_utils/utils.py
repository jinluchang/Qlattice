import cqlat_utils as cu

from qlat_utils.cache import *
from qlat_utils.rng_state import *

import sys

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
    import psutil
    rss = psutil.Process().memory_info().rss / (1024 * 1024 * 1024)
    displayln_info(f"show_memory_usage: rss = {rss:.6f} GB")
    malloc_stats_info()

def get_all_caches_info():
    return cu.get_all_caches_info()

def clear_all_caches():
    # clean python level cache and then C++ level cache
    clean_cache()
    cu.clear_all_caches()

def malloc_stats():
    return cu.malloc_stats()

def malloc_stats_info():
    if get_id_node() == 0:
        return malloc_stats()

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
