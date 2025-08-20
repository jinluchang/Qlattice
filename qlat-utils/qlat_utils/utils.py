from .timer import *
from .cache import *
from .c import *
from .json import *

import math
import sys
import os
import numpy as np
import inspect
import importlib
import importlib.util

def getenv(*names, default=None):
    assert len(names) > 0
    for name in names:
        val = os.getenv(name)
        if val is not None:
            displayln_info(0, f"{name}='{val}'")
            return val
    val = default
    displayln_info(0, f"{names[0]}='{val}' (default)")
    return val

def get_arg(option, default=None, *, argv=None, is_removing_from_argv=False):
    """
    Get the `arg` of the option when it first appears.
    Remove the option and its arg if `is_removing_from_argv`.
    """
    if argv is None:
        argv = sys.argv
    i_max = len(argv) - 1
    for i in range(len(argv)):
        if argv[i] == option:
            if i == i_max:
                if is_removing_from_argv:
                    argv.pop(i)
                return ""
            else:
                arg = argv[i + 1]
                if is_removing_from_argv:
                    argv.pop(i)
                    argv.pop(i)
                return arg
    return default

def get_arg_list(option, *, argv=None, is_removing_from_argv=False):
    """
    Get all the `arg` of the option when it appears as `arg_list`, it may appear multiple times.
    Remove the options and the args if `is_removing_from_argv`.
    """
    if argv is None:
        argv = sys.argv
    arg_list = []
    i = 0
    while i < len(argv):
        i_max = len(argv) - 1
        if argv[i] == option:
            if i == i_max:
                if is_removing_from_argv:
                    argv.pop(i)
                else:
                    i += 1
                arg_list.append("")
            else:
                arg = argv[i + 1]
                if is_removing_from_argv:
                    argv.pop(i)
                    argv.pop(i)
                else:
                    i += 2
                arg_list.append(arg)
        else:
            i += 1
    return arg_list

def get_option(option, *, argv=None, is_removing_from_argv=False):
    """
    Return if `option` in `argv`.
    Remove the option if `is_removing_from_argv`
    """
    if argv is None:
        argv = sys.argv
    if option in argv:
        if is_removing_from_argv:
            argv.remove(option)
        return True
    else:
        return False

def get_all_arg_list(option, default=None, *, argv=None, is_removing_from_argv=False):
    """
    Get all the following args as `arg_list` of the option when it first appears.
    Remove the option and the remaining args if `is_removing_from_argv`.
    """
    if argv is None:
        argv = sys.argv
    i_max = len(argv) - 1
    for i in range(len(argv)):
        if argv[i] == option:
            arg_list = argv[i + 1:]
            if is_removing_from_argv:
                argv[i:] = []
            return arg_list
    return default

is_test_state = get_option("--test")

def is_test():
    return is_test_state

def show_memory_usage():
    try:
        import psutil
        rss = psutil.Process().memory_info().rss / (1024 * 1024 * 1024)
        displayln_info(f"show_memory_usage: rss = {rss:.6f} GB")
        # displayln_info_malloc_stats()
    except:
        displayln_info(f"show_memory_usage: no psutil.")

def import_file(module_name, file_path):
    """
    return the imported module
    """
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module

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

@timer
def get_fname():
    """
    Return the function name of the current function ``fname``
    """
    f = inspect.currentframe().f_back
    return f.f_code.co_name

def sqr(x):
    return x * x

def set_zero(x):
    x.set_zero()

def set_unit(x, coef=1.0):
    x.set_unit(coef)

def show(x):
    return x.show()

def unitarize(x):
    x.unitarize()

def get_chunk_list(total_list, *, chunk_size=None, chunk_number=None, rng_state=None):
    """
    Split ``total_list`` into ``chunk_number`` chunks or chunks with ``chunk_size``.
    One (and only one) of ``chunk_size`` and ``chunk_number`` should not be ``None``.
    #
    Returns a list of chunks.
    Number of chunks is less or equal to ``chunk_number``.
    Chunk sizes are less or equal to ``chunk_size``.
    if rng_state is None:
        Do not randomly permute the list
    """
    assert chunk_size is None or chunk_number is None
    assert chunk_size is not None or chunk_number is not None
    chunk_list = []
    if rng_state is not None:
        assert isinstance(rng_state, RngState)
        total_list = random_permute(total_list, rng_state)
    total = len(total_list)
    if chunk_size is not None:
        assert isinstance(chunk_size, int)
        assert chunk_size >= 1
        chunk_number = (total - 1) // chunk_size + 1
    elif chunk_number is not None:
        assert isinstance(chunk_number, int)
        assert chunk_number >= 1
        chunk_size = (total - 1) // chunk_number + 1
    for i in range(chunk_number):
        start = min(i * chunk_size, total);
        stop = min(start + chunk_size, total);
        if stop > start:
            chunk_list.append(total_list[start:stop])
    return chunk_list

def parse_grid_coordinate_str(x_str):
    x_str_list = x_str.split(".")
    x_list = [ int(s) for s in x_str_list ]
    x = Coordinate(x_list)
    return x

def mk_epsilon_array():
    arr = np.zeros((4, 4, 4, 4,), dtype = np.int8)
    def setv(i, j, k, l):
        setw(i, j, k, l, 1)
        setw(i, j, l, k, -1)
    def setw(i, j, k, l, val):
        arr[i, j, k, l] = val;
        arr[j, k, l, i] = -val;
        arr[k, l, i, j] = val;
        arr[l, i, j, k] = -val;
    setv(0, 1, 2, 3)
    setv(0, 2, 3, 1);
    setv(0, 3, 1, 2);
    return arr

epsilon_array = mk_epsilon_array()

def epsilon_tensor(i, j, k, l=3):
    """
    epsilon_tensor(0, 1, 2, 3) == 1
    """
    return epsilon_array[i, j, k, l].item()

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
    # assert x.shape == size.shape
    ans = x % size
    assert np.all(ans >= 0)
    mask = 2 * ans >= size
    ans[mask] = (ans - size)[mask]
    return ans

def rel_mod_sym_arr(x, size):
    """
    Return ``x % size`` or ``x % size - size`` or ``0`` where ``x`` and ``size`` are np.array of same shape
    """
    # assert x.shape == size.shape
    ans = x % size
    assert np.all(ans >= 0)
    mask1 = 2 * ans > size
    mask2 = 2 * ans == size
    ans[mask1] = (ans - size)[mask1]
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
    return math.sqrt(sum([ (l / 2)**2 for l in total_site.to_list()[:3] ]))

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

def mk_r_sq_list(r_sq_limit, dimension='3D'):
    if dimension == '4D':
        # Lagrange's four-square theorem
        # https://en.wikipedia.org/wiki/Lagrange%27s_four-square_theorem
        return list(range(0, r_sq_limit))
    elif dimension == '3D':
        return mk_r_sq_list_3d(r_sq_limit)
    else:
        raise Exception(f"mk_r_sq_list: dimension='{dimension}' not recognized.")

def mk_r_list(r_limit, *, r_all_limit=28.0, r_scaling_factor=5.0, dimension='3D'):
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

@timer
def get_data_sig_arr(x, rs, sig_len):
    """
    Return a signature (an array of floating point number, real or complex) of data viewed as a 1-D array of numbers.\n
    Call `get_data_sig` several times internally.
    Result only depends on the value of the data, not the structure.
    ``x`` can be an instance of ``LatData``, ``np.ndarray``, etc.
    """
    assert isinstance(rs, RngState)
    assert isinstance(sig_len, int)
    sig_list = []
    for i in range(sig_len):
        sig = get_data_sig(x, rs.split(f"{i}"))
        sig_list.append(sig)
    sig_arr = np.array(sig_list)
    return sig_arr

global_json_results = [] # Default value for param `json_results` in functions `check_log_json` and `json_results_append`

def json_results_append(*args, json_results=None):
    if json_results is None:
        json_results = global_json_results
    displayln_info(0, r"//------------------------------------------------------------\\")
    displayln_info(0, *args)
    displayln_info(0, r"\\------------------------------------------------------------//")
    json_results.append(args)

@timer_verbose
def check_log_json(script_file, *, json_results=None, check_eps=1e-5):
    """
    q.check_log_json(__file__, json_results=json_results, check_eps=check_eps)
    #
    json_results = [ (name, value, check_eps,), (name, value,), (name,), ... ]
    #
    default:
    check_eps=1e-5
    json_results=q.global_json_results
    #
    """
    fname = get_fname()
    if json_results is None:
        json_results = global_json_results
    if 0 == get_id_node():
        json_fn_name = os.path.splitext(script_file)[0] + ".log.json"
        qtouch(json_fn_name + ".new", json_dumps(json_results, indent=1))
        if does_file_exist_qar(json_fn_name):
            json_results_load = json_loads(qcat(json_fn_name))
            for i, (p, pl,) in enumerate(zip(json_results, json_results_load)):
                if len(p) != len(pl):
                    displayln(-1, f"CHECK: {i} {p} load:{pl}")
                    displayln(-1, f"CHECK: ERROR: JSON results length does not match.")
                    continue
                if len(p) == 1:
                    eps = check_eps
                    epsl = check_eps
                    v = 0
                    vl = 0
                    n, = p
                    nl, = pl
                elif len(p) == 2:
                    eps = check_eps
                    epsl = check_eps
                    n, v = p
                    nl, vl = pl
                elif len(p) == 3:
                    n, v, eps = p
                    nl, vl, epsl = pl
                else:
                    displayln(-1, f"CHECK: {i} {p} load:{pl}")
                    displayln(-1, f"CHECK: ERROR: JSON results length not 2 or 3.")
                    continue
                if n != nl:
                    displayln(-1, f"CHECK: {i} {p} load:{pl}")
                    displayln(-1, f"CHECK: ERROR: JSON results item does not match.")
                    continue
                if eps != epsl:
                    displayln(-1, f"CHECK: {i} {p} load:{pl}")
                    displayln(-1, f"CHECK: ERROR: JSON results eps does not match.")
                    continue
                actual_eps = 0.0
                v_norm = np.linalg.norm(v)
                vl_norm = np.linalg.norm(vl)
                diff_norm = np.linalg.norm(v - vl)
                if (v_norm + vl_norm) > 0:
                    actual_eps = 2 * diff_norm / (v_norm + vl_norm)
                if actual_eps > eps:
                    displayln(-1, f"CHECK: {i} '{n}' actual: {v} ; load: {vl} .")
                    displayln(-1, f"CHECK: target eps: {eps} ; actual eps: {actual_eps} .")
                    displayln(-1, f"CHECK: ERROR: JSON results value does not match.")
                elif actual_eps != 0.0:
                    displayln(-1, f"INFO: {fname}: {i} '{n}'")
                    displayln(-1, f"INFO: {fname}: target eps: {eps} ; actual eps: {actual_eps} .")
            if len(json_results) != len(json_results_load):
                displayln(-1, f"CHECK: len(json_results)={len(json_results)} load:{len(json_results_load)}")
                displayln(-1, f"CHECK: ERROR: JSON results len does not match.")
