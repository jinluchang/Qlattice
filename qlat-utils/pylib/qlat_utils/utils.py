import qlat_utils.c as c

from qlat_utils.cache import *
from qlat_utils.rng_state import *

from qlat_utils.c import random_permute, displayln_malloc_stats

import math
import sys
import os
import numpy as np

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

def show_memory_usage():
    try:
        import psutil
        rss = psutil.Process().memory_info().rss / (1024 * 1024 * 1024)
        displayln_info(f"show_memory_usage: rss = {rss:.6f} GB")
        # displayln_info_malloc_stats()
    except:
        displayln_info(f"show_memory_usage: no psutil.")

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
    """
    Return ``x % size`` or ``x % size - size``
    """
    x = x % size
    assert x >= 0
    if 2 * x >= size:
        return x - size
    else:
        return x

def rel_mod_sym(x, size):
    """
    Return ``x % size`` or ``x % size - size`` or ``0``
    """
    x = x % size
    assert x >= 0
    if 2 * x > size:
        return x - size
    elif 2 * x < size:
        return x
    else:
        assert 2 * x == size
        return 0

def rel_mod_arr(x, size):
    """
    Return ``x % size`` or ``x % size - size`` where ``x`` and ``size`` are np.array of same shape
    """
    assert x.shape == size.shape
    ans = x % size
    assert np.all(ans >= 0)
    mask = 2 * ans >= size
    ans[mask] = ans[mask] - size[mask]
    return ans

def rel_mod_sym_arr(x, size):
    """
    Return ``x % size`` or ``x % size - size`` or ``0`` where ``x`` and ``size`` are np.array of same shape
    """
    assert x.shape == size.shape
    ans = x % size
    assert np.all(ans >= 0)
    mask1 = 2 * ans > size
    mask2 = 2 * ans == size
    ans[mask1] = ans[mask1] - size[mask1]
    ans[mask2] = 0
    return ans

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

def get_r_sq(x_rel):
    """
    get spatial distance square as int
    """
    return sum([ x * x for x in x_rel[:3] ])

def get_r_limit(total_site):
    """
    Return the limit for spatial ``r`` as float.\n
    :params total_site: must be Coordinate type
    """
    return math.sqrt(sum([ (l / 2)**2 for l in total_site.list()[:3] ]))

def mk_r_sq_list_3d(r_sq_limit):
    r_limit = int(math.sqrt(r_sq_limit))
    r_sq_set = set()
    for x in range(0, r_limit + 1):
        for y in range(0, x + 1):
            for z in range(0, y + 1):
                r_sq = x**2 + y**2 + z**2
                if r_sq > r_sq_limit:
                    continue
                r_sq_set.add(r_sq)
    return sorted(list(r_sq_set))

def mk_r_sq_list(r_sq_limit, dimension = '3D'):
    if dimension == '4D':
        # Lagrange's four-square theorem
        # https://en.wikipedia.org/wiki/Lagrange%27s_four-square_theorem
        return list(range(0, r_sq_limit))
    elif dimension == '3D':
        return mk_r_sq_list_3d(r_sq_limit)
    else:
        raise Exception(f"mk_r_sq_list: dimension='{dimension}' not recognized.")

def mk_r_list(r_limit, *, r_all_limit = 28.0, r_scaling_factor = 5.0, dimension = '3D'):
    """
    Make a list of `r` values from `0` up to `r_limit`.\n
    Parameters
    ----------
    r_limit: the limit for the generated `r` list.
    r_scaling_factor: After `r_all_limit`, include `r` with integer values divide `r_scaling_factor`
    r_all_limit: include all possible `r` values up to (include) this limit.
    dimension: '3D' or '4D'
    """
    r_list = [ math.sqrt(r_sq) for r_sq in mk_r_sq_list(int(min(r_limit, r_all_limit)**2 + 0.5)) ]
    r_second_start = r_all_limit
    if r_list:
        r_second_start = min(r_second_start, r_list[-1])
    r_second_start_idx = int(r_second_start * r_scaling_factor + 1.5)
    r_second_stop_idx = math.ceil(r_limit * r_scaling_factor + 1.5)
    for i in range(r_second_start_idx, r_second_stop_idx):
        r = i / r_scaling_factor
        r_list.append(r)
    return r_list

def mk_interp_tuple(x, x0, x1, x_idx):
    """
    Returns `(x_idx_low, x_idx_high, coef_low, coef_high,)`\n
    `x_idx` corresponds to `x0`
    `x_idx + 1` corresponds to `x1`
    """
    assert x0 <= x and x <= x1
    x_idx_low = x_idx
    x_idx_high = x_idx_low + 1
    x_interval = x1 - x0
    coef_low = (x1 - x) / x_interval
    coef_high = (x - x0) / x_interval
    return (x_idx_low, x_idx_high, coef_low, coef_high,)

def mk_r_sq_interp_idx_coef_list(r_list):
    """
    Return a list of tuples:\n
    ``r_sq_interp_idx_coef_list = [ (r_idx_low, r_idx_high, coef_low, coef_high,), ... ]``
    where:
    ``r_sq_interp_idx_coef_list[r_sq] = (r_idx_low, r_idx_high, coef_low, coef_high,)``
    `r_sq` ranges from `0` to `int(r_list[-1]**2 + 1.5)`
    """
    r_list_len = len(r_list)
    assert r_list_len >= 2
    assert r_list[0] == 0.0
    r_sq_list = list(range(0, int(r_list[-1]**2 + 1.5)))
    r_idx = 0
    r_sq_interp_idx_coef_list = []
    for r_sq in r_sq_list:
        r = math.sqrt(r_sq)
        while True:
            if r_idx + 1 >= r_list_len:
                r_sq_interp_idx_coef_list.append((r_idx, r_idx, 1.0, 0.0,))
                break
            r0 = r_list[r_idx]
            r1 = r_list[r_idx + 1]
            if r0 <= r and r < r1:
                r_sq_interp_idx_coef_list.append(mk_interp_tuple(r, r0, r1, r_idx))
                break
            r_idx += 1
    return r_sq_interp_idx_coef_list
