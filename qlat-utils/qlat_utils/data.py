from .rng_state import *
from .timer import *

class q:
    from .utils import (
            get_fname,
            )

import math
import copy
import functools
import numpy as np

alpha_qed = 1.0 / 137.035999084
fminv_gev = 0.197326979 # hbar * c / (1e-15 m * 1e9 electron charge * 1 volt)

float_types = (float, np.float32, np.float64,)
complex_types = (complex, np.complex64, np.complex128,)
int_types = (int, np.int32, np.int64,)

try:
    float_types = float_types + (np.float128,)
    complex_types = complex_types + (np.complex256,)
except:
    pass

real_types = float_types + int_types
number_types = real_types + complex_types

class use_kwargs:

    """
    self.default_kwargs
    self.keys
    """

    def __init__(self, default_kwargs, keys = None):
        self.default_kwargs = default_kwargs
        self.keys = None

    def __call__(self, func):
        @functools.wraps(func)
        def f(*args, **kwargs):
            if "is_default_kwargs_applied" not in kwargs:
                d = self.default_kwargs.copy()
                d.update(kwargs)
                kwargs = d
            if self.keys is not None:
                kwargs = { k: kwargs[k] for k in self.keys }
            return func(*args, **kwargs)
        return f

###

def interp_i_arr(data_x_arr, x_arr):
    r"""
    return i_arr
    q.interp(data_x_arr, i_arr) \approx x_arr
    #
    `x_arr` can be either an 1-D array-like object or a single number.
    #
    e.g.:
    data(x)
    data_arr[:] = data(data_x_arr)
    q.interp(data_arr, i_arr) \approx data(x_arr)
    """
    data_i_arr = np.arange(len(data_x_arr))
    i_arr = np.interp(x_arr, data_x_arr, data_i_arr)
    return i_arr

def interp(data_arr, i_arr, axis=-1):
    """
    return approximately data_arr[..., i_arr] if axis=-1
    #
    Note that `i_arr` can be non-integer.
    #
    `i_arr` can be either an 1-D array-like object or a single number.
    """
    v_arr = np.asarray(data_arr)
    v_arr = np.moveaxis(v_arr, axis, 0)
    i_arr = np.asarray(i_arr)
    shape = i_arr.shape
    if shape == ():
        i = i_arr.item()
        size = len(v_arr)
        i1 = math.floor(i)
        assert i1 >= 0
        i2 = i1 + 1
        if i2 >= size:
            return v_arr[size - 1]
        elif i1 < 0:
            return v_arr[0]
        v1 = v_arr[i1]
        v2 = v_arr[i2]
        a1 = i2 - i
        a2 = i - i1
        return a1 * v1 + a2 * v2
    elif shape == (len(i_arr),):
        iv_arr = np.array([ interp(v_arr, i, 0) for i in i_arr ], v_arr.dtype)
        iv_arr = np.moveaxis(iv_arr, 0, axis)
        return iv_arr
    else:
        fname = q.get_fname()
        raise Exception(f"{fname}: i_arr={i_arr}")

def interp_x(data_arr, data_x_arr, x_arr, axis=-1):
    """
    return interpolated_data_arr
    #
    `x_arr` can be either an 1-D array-like object or a single number.
    #
    `data_x_arr` is the x values for `data_arr`
    `x_arr` is the x values for `interpolated_data_arr`
    ``
    data_x_arr.shape == (data_arr.shape[axis],)
    ``
    If len(x_arr)
    ``
    interpolated_data_arr.shape[axis] == len(x_arr)
    len(data_arr.shape) == len(interpolated_data_arr.shape)
    ``
    """
    assert data_x_arr.shape == (data_arr.shape[axis],)
    i_arr = interp_i_arr(data_x_arr, x_arr)
    interpolated_data_arr = interp(data_arr, i_arr, axis)
    return interpolated_data_arr

def get_threshold_idx(arr, threshold):
    """
    return x
    #
    ``
    q.interp(arr, [ x, ]) = np.array([ threshold, ])
    arr.shape == (len(arr),)
    ``
    """
    i1 = 0
    i2 = len(arr) - 1
    v1 = arr[i1]
    v2 = arr[i2]
    if v1 >= v2:
        i1, i2 = i2, i1
        v1, v2 = v2, v1
    while True:
        assert v2 >= v1
        if v1 <= threshold and threshold <= v2:
            if i2 - i1 == 1:
                d_v = v2 - v1
                d_i = i2 - i1
                i3 = i1 + (threshold - v1) / d_v * d_i
                return i3
            i3 = (i1 + i2) // 2
            v3 = arr[i3]
            if threshold <= v3:
                i2 = i3
                v2 = v3
                continue
            elif v3 <= threshold:
                i1 = i3
                v1 = v3
                continue
            else:
                assert False
        elif threshold <= v1:
            return i1
        elif v2 <= threshold:
            return i2
        else:
            assert False
    assert False

def get_threshold_i_arr(data_arr, threshold_arr, axis=-1):
    r"""
    return i_arr
    #
    let `shape` = `np.moveaxis(data_arr, axis, -1)[..., 0].shape`
    #
    threshold_arr = np.broadcast_to(threshold_arr, shape)
    #
    such that
    for index in np.ndindex(shape):
        q.interp(data_arr[index], i_arr[index]) \approx threshold_arr[index]
    """
    v_arr = np.asarray(data_arr)
    threshold_arr = np.asarray(threshold_arr)
    v_arr = np.moveaxis(v_arr, axis, -1)
    shape = v_arr[..., 0].shape
    threshold_arr = np.broadcast_to(threshold_arr, shape)
    i_arr = np.zeros(shape, dtype=np.float64)
    for index in np.ndindex(shape):
        t = threshold_arr[index]
        arr = v_arr[index]
        i_arr[index] = get_threshold_idx(arr, t)
    return i_arr

def get_threshold_x_arr(data_arr, data_x_arr, threshold_arr, axis=-1):
    r"""
    return x_arr
    #
    ``
    data_x_arr.shape == (data_arr.shape[axis],)
    ``
    #
    let `shape` = `np.moveaxis(data_arr, axis, -1)[..., 0].shape`
    #
    threshold_arr = np.broadcast_to(threshold_arr, shape)
    #
    such that
    for index in np.ndindex(shape):
        q.interp_x(data_arr[index], data_x_arr, x_arr[index]) \approx threshold_arr[index]
    """
    assert data_x_arr.shape == (data_arr.shape[axis],)
    i_arr = get_threshold_i_arr(data_arr, threshold_arr, axis)
    x_arr = np.zeros(i_arr.shape, dtype=np.float64)
    x_arr.ravel()[:] = interp(data_x_arr, i_arr.ravel())
    return x_arr

def partial_sum_list(x, *, is_half_last = False):
    """Modify in-place, preserve length"""
    s = 0
    for i, v in enumerate(x):
        sp = s
        s += v
        if is_half_last:
            x[i] = (s + sp) / 2
        else:
            x[i] = s

def partial_sum(x, *, is_half_last = False):
    """Modify in-place, preserve length"""
    shape = x.shape
    if len(shape) == 0:
        return
    elif len(shape) == 1:
        partial_sum_list(x, is_half_last = is_half_last)
    elif len(shape) == 2:
        for v in x:
            partial_sum_list(v, is_half_last = is_half_last)
    else:
        assert False

def check_zero(x):
    if isinstance(x, real_types) and 0 == x:
        return True
    return False

def qnorm(x):
    """
    qnorm(2) == 4
    """
    if isinstance(x, np.ndarray):
        return np.abs(np.vdot(x, x))
    elif isinstance(x, real_types):
        return x * x
    elif isinstance(x, complex_types):
        return x.real * x.real + x.imag * x.imag
    elif isinstance(x, (list, tuple,)):
        return sum([ qnorm(x_i) for x_i in x ])
    else:
        return x.qnorm()
    assert False

class Data:

    def __init__(self, val):
        """
        # supported value types:
        # numeric
        # numpy.array
        # q.LatData
        # list
        """
        if isinstance(val, Data):
            self.val = val.val
            assert not isinstance(self.val, Data)
        else:
            self.val = val

    def __str__(self):
        return f"Data({self.val})"

    def get_val(self):
        return self.val

    def __copy__(self):
        return Data(copy.copy(self.val))

    def __deepcopy__(self, memo):
        return Data(copy.deepcopy(self.val, memo))

    def __add__(self, other):
        if isinstance(other, Data):
            if check_zero(self.val):
                return other
            elif check_zero(other.val):
                return self
            elif isinstance(self.val, list) and isinstance(other.val, list):
                assert len(self.val) == len(other.val)
                return Data([ v1 + v2 for v1, v2 in zip(self.val, other.val) ])
            elif isinstance(self.val, list):
                return Data([ v + other.val for v in self.val ])
            elif isinstance(other.val, list):
                return Data([ self.val + v for v in other.val ])
            else:
                return Data(self.val + other.val)
        else:
            return self + Data(other)

    def __radd__(self, other):
        if isinstance(other, Data):
            assert False
            return None
        else:
            return Data(other) + self

    def __mul__(self, other):
        if isinstance(other, Data):
            if check_zero(self.val) or check_zero(other.val):
                return Data(0)
            elif isinstance(self.val, list) and isinstance(other.val, list):
                return Data([ v1 * v2 for v1, v2 in zip(self.val, other.val) ])
            elif isinstance(self.val, list):
                return Data([ v * other.val for v in self.val ])
            elif isinstance(other.val, list):
                return Data([ self.val * v for v in other.val ])
            return Data(self.val * other.val)
        else:
            return self * Data(other)

    def __rmul__(self, other):
        if isinstance(other, Data):
            assert False
            return None
        else:
            return Data(other) * self

    def __neg__(self):
        if check_zero(self.val):
            return Data(0)
        elif isinstance(self.val, list):
            return Data([ -v for v in self.val ])
        else:
            return Data(-self.val)

    def __pos__(self):
        return self

    def __sub__(self, other):
        if isinstance(other, Data):
            if check_zero(self.val):
                return Data(-other.val)
            elif check_zero(other.val):
                return self
            elif isinstance(self.val, list) and isinstance(other.val, list):
                return Data([ v1 - v2 for v1, v2 in zip(self.val, other.val) ])
            elif isinstance(self.val, list):
                return Data([ v - other.val for v in self.val ])
            elif isinstance(other.val, list):
                return Data([ self.val - v for v in other.val ])
            else:
                return Data(self.val - other.val)
        else:
            return self - Data(other)

    def __rsub__(self, other):
        if isinstance(other, Data):
            assert False
            return None
        else:
            return Data(other) - self

    def qnorm(self):
        return qnorm(self.val)

    def glb_sum(self):
        from qlat.mpi import glb_sum
        return Data(glb_sum(self.val))

###

def average(data_list):
    n = len(data_list)
    v = sum(data_list)
    return 1/n * v

def average_ignore_nan(value_arr_list):
    if len(value_arr_list) == 0:
        return None
    shape = value_arr_list[0].shape
    dtype = value_arr_list[0].dtype
    count_arr = np.zeros(shape, dtype=np.int64)
    sum_arr = np.zeros(shape, dtype=dtype)
    for v_arr in value_arr_list:
        assert v_arr.shape == shape
        assert v_arr.dtype == dtype
        sel = ~np.isnan(v_arr)
        count_arr[sel] += 1
        sum_arr[sel] += v_arr[sel]
    avg_arr = np.zeros(shape, dtype=dtype)
    sel = count_arr > 0
    avg_arr[sel] = sum_arr[sel] / count_arr[sel]
    avg_arr[~sel] = np.nan
    return avg_arr

def block_data(data_list, block_size, is_overlapping=True):
    """
    return the list of block averages
    the blocks may overlap if is_overlapping == True
    """
    if block_size == 1:
        return data_list
    assert block_size >= 1
    size = len(data_list)
    if block_size >= size:
        return [ average(data_list), ]
    blocks = []
    start = 0
    stop = block_size
    while stop <= size:
        b = average(data_list[start:stop])
        blocks.append(b)
        if is_overlapping:
            start += 1
            stop += 1
        else:
            start += block_size
            stop += block_size
    return blocks

def avg_err(data_list, eps=1, *, block_size=1):
    avg = average(data_list)
    blocks = block_data(data_list, block_size)
    diff_sqr = average([ fsqr(d - avg) for d in blocks ])
    fac = abs(eps) * math.sqrt(block_size / (len(data_list) - 1))
    err = fac * fsqrt(diff_sqr)
    return (avg, err,)

def jackknife(data_list, eps=1):
    r"""
    Return jk[i] = avg - \frac{eps}{N} (v[i] - avg)
    normal jackknife uses eps=1, scale the fluctuation by eps
    """
    is_np_arr = isinstance(data_list, np.ndarray)
    data_list_real = [ d for d in data_list if d is not None ]
    n = len(data_list_real)
    fac = eps / n
    avg = average(data_list_real)
    jks = [ avg, ]
    for data in data_list:
        if data is None:
            jks.append(avg)
        else:
            jks.append(avg - fac * (data - avg))
    if is_np_arr:
        jks = np.array(jks, dtype=data_list.dtype)
    return jks

def fsqr(data):
    if isinstance(data, real_types):
        return data * data
    elif isinstance(data, complex_types):
        r = data.real
        i = data.imag
        return complex(r * r, i * i)
    elif isinstance(data, Data):
        return Data(fsqr(data.val))
    else:
        # Assuming np.ndarray like object
        if data.dtype in real_types:
            return np.square(data)
        elif data.dtype in complex_types:
            return np.square(data.real) + 1j * np.square(data.imag)
        else:
            raise Exception(f"fsqr data={data} type not supported")

def fsqrt(data):
    if isinstance(data, real_types):
        return math.sqrt(data)
    elif isinstance(data, complex_types):
        r = data.real
        i = data.imag
        return complex(math.sqrt(r), math.sqrt(i))
    elif isinstance(data, Data):
        return Data(fsqrt(data.val))
    else:
        # Assuming np.ndarray like object
        if data.dtype in real_types:
            return np.sqrt(data)
        elif data.dtype in complex_types:
            return np.sqrt(data.real) + 1j * np.sqrt(data.imag)
        else:
            raise Exception(f"fsqr data={data} type not supported")

def err_sum(*vs):
    """
    e.g.: `q.err_sum(1.4, 2.1, 1.0)` ==> `2.7147743920996454`
    """
    err_sqr = sum([ fsqr(v) for v in vs])
    err = fsqrt(err_sqr)
    return err

def jk_avg(jk_list):
    is_np_arr = isinstance(jk_list, np.ndarray)
    val = jk_list[0]
    if is_np_arr and val.size == 1:
        return val.item()
    else:
        return val

def jk_err(jk_list, eps=1, *, block_size=1):
    r"""
    Return \frac{1}{eps} \sqrt{ \sum_{i=1}^N (jk[i] - jk_avg)^2 } when block_size=1
    Note: len(jk_list) = N + 1
    Same eps as the eps used in the 'jackknife' function
    """
    is_np_arr = isinstance(jk_list, np.ndarray)
    avg = jk_avg(jk_list)
    blocks = block_data(jk_list[1:], block_size)
    diff_sqr = average([ fsqr(jk - avg) for jk in blocks ])
    fac = math.sqrt(block_size * (len(jk_list) - 1)) / abs(eps)
    val = fac * fsqrt(diff_sqr)
    if is_np_arr and val.size == 1:
        return val.item()
    else:
        return val

def jk_avg_err(jk_list, eps=1, *, block_size=1):
    return jk_avg(jk_list), jk_err(jk_list, eps, block_size=block_size)

def merge_jk_idx(*jk_idx_list):
    for jk_idx in jk_idx_list:
        assert jk_idx[0] == "avg"
    return [ "avg", ] + [ idx for jk_idx in jk_idx_list for idx in jk_idx[1:] ]

def rejk_list(jk_list, jk_idx_list, all_jk_idx):
    """
    super jackknife
    """
    assert jk_idx_list[0] == "avg"
    assert all_jk_idx[0] == "avg"
    assert len(jk_idx_list) == len(jk_list)
    assert len(jk_idx_list) <= len(all_jk_idx)
    is_np_arr = isinstance(jk_list, np.ndarray)
    jk_avg = jk_list[0]
    size_new = len(all_jk_idx)
    i_new = 0
    jk_list_new = []
    for i, idx in enumerate(jk_idx_list):
        while all_jk_idx[i_new] != idx:
            jk_list_new.append(jk_avg)
            i_new += 1
            assert i_new < size_new
        jk_list_new.append(jk_list[i])
        i_new += 1
    while i_new < size_new:
        jk_list_new.append(jk_avg)
        i_new += 1
    assert i_new == size_new
    assert size_new == len(jk_list_new)
    if is_np_arr:
        jk_list_new = np.array(jk_list_new, dtype=jk_list.dtype)
    return jk_list_new

# ----------

@timer
def rjk_jk_list(jk_list, jk_idx_list, n_rand_sample, rng_state, jk_blocking_func=None, is_normalizing_rand_sample=True, is_use_old_rand_alg=False):
    r"""
    return rjk_list
    len(rjk_list) == 1 + n_rand_sample
    distribution of rjk_list should be similar as the distribution of avg
    r_{i,j} ~ N(0, 1)
    avg = jk_list[0]
    len(jk_list) = n + 1
    rjk_list[i] = avg + \sum_{j=1}^{n} r_{i,j} (jk_list[j] - avg)
    #
    jk_blocking_func(jk_idx) => blocked jk_idx
    """
    assert jk_idx_list[0] == "avg"
    assert isinstance(n_rand_sample, int_types)
    assert n_rand_sample >= 0
    assert isinstance(rng_state, RngState)
    is_np_arr = isinstance(jk_list, np.ndarray)
    rs = rng_state
    n = len(jk_list) - 1
    if jk_blocking_func is None:
        blocked_jk_idx_list = jk_idx_list
    else:
        blocked_jk_idx_list = [ jk_blocking_func(idx) for idx in jk_idx_list[:] ]
    assert len(blocked_jk_idx_list[1:]) == n
    r_arr = np.empty((n_rand_sample, n,), dtype=np.float64)
    if is_use_old_rand_alg == "v1":
        assert not is_normalizing_rand_sample
        for i in range(n_rand_sample):
            rsi = rs.split(str(i))
            r = [ rsi.split(str(idx)).g_rand_gen() for idx in blocked_jk_idx_list[1:] ]
            for j in range(n):
                r_arr[i, j] = r[j]
    else:
        assert is_use_old_rand_alg == False
        r_arr_dict = dict()
        for jk_idx in blocked_jk_idx_list[1:]:
            jk_idx_str = str(jk_idx)
            if jk_idx_str in r_arr_dict:
                continue
            rsi = rs.split(str(jk_idx))
            garr = rsi.g_rand_arr(n_rand_sample)
            if is_normalizing_rand_sample:
                garr_qnorm = qnorm(garr) # garr_qnorm \approx n_rand_sample
                garr = garr * np.sqrt(n_rand_sample / garr_qnorm)
                assert abs(qnorm(garr) / n_rand_sample - 1) < 1e-8
            r_arr_dict[jk_idx_str] = garr
        for j, jk_idx in enumerate(blocked_jk_idx_list[1:]):
            jk_idx_str = str(jk_idx)
            r_arr[:, j] = r_arr_dict[jk_idx_str]
    avg = jk_list[0]
    if is_np_arr:
        jk_arr = jk_list
        jk_diff = jk_arr[1:] - avg
        rjk_arr = np.empty((1 + n_rand_sample, *avg.shape,), dtype=jk_arr.dtype)
        rjk_arr[:] = avg
        for j in range(n):
            for i in range(n_rand_sample):
                rjk_arr[i + 1] += r_arr[i, j] * jk_diff[j]
        return rjk_arr
    else:
        rjk_list = [ avg, ]
        jk_diff = [ jk_list[j] - avg for j in range(1, n + 1) ]
        for i in range(n_rand_sample):
            rjk_list.append(avg + sum([ r_arr[i, j] * jk_diff[j] for j in range(n) ]))
        return rjk_list

@timer
def rjk_mk_jk_val(rs_tag, val, err, n_rand_sample, rng_state, is_normalizing_rand_sample=True):
    """
    return rjk_list
    n = n_rand_sample
    len(rjk_list) == 1 + n
    rjk_list[i] = val + err * r[i] for i in 1..n
    where r[i] ~ N(0, 1)
    """
    assert n_rand_sample >= 0
    assert isinstance(rng_state, RngState)
    assert isinstance(val, real_types)
    assert isinstance(err, real_types)
    rs = rng_state.split(str(rs_tag))
    rjk_arr = np.zeros((n_rand_sample + 1,), dtype=np.float64)
    rjk_arr[0] = val
    r_arr = rs.g_rand_arr((n_rand_sample,))
    if is_normalizing_rand_sample:
        r_arr_qnorm = qnorm(r_arr)
        r_arr = r_arr * np.sqrt(n_rand_sample / r_arr_qnorm)
        assert abs(qnorm(r_arr) / n_rand_sample - 1) < 1e-8
    rjk_arr[1:] = val + r_arr * err
    return rjk_arr

def rjackknife(data_list, jk_idx_list, n_rand_sample, rng_state, *, eps=1):
    jk_list = jackknife(data_list, eps)
    return rjk_jk_list(jk_list, jk_idx_list, n_rand_sample, rng_state)

def rjk_avg(rjk_list):
    return jk_avg(rjk_list)

def rjk_err(rjk_list, eps=1):
    """Return \\frac{1}{eps} \\sqrt{ \\frac{1}{N} \\sum_{i=1}^N (jk[i] - jk_avg)^2 }
    Note: len(jk_list) = N + 1
    Same eps as the eps used in the 'jackknife' function"""
    n = len(rjk_list) - 1
    return jk_err(rjk_list, abs(eps) * np.sqrt(n))

def rjk_avg_err(rjk_list, eps=1):
    return rjk_avg(rjk_list), rjk_err(rjk_list, eps)

# ----------

default_g_jk_kwargs = dict()

def mk_g_jk_kwargs():
    """
    Return the predefined `default_g_jk_kwargs`.
    """
    g_jk_kwargs = dict()
    #
    g_jk_kwargs["jk_type"] = "rjk"  # choices: "rjk", "super"
    g_jk_kwargs["eps"] = 1
    #
    # for jk_type = "super"
    g_jk_kwargs["all_jk_idx"] = None
    g_jk_kwargs["get_all_jk_idx"] = None
    #
    # for jk_type = "rjk"
    g_jk_kwargs["n_rand_sample"] = 1024
    g_jk_kwargs["rng_state"] = RngState("rejk")
    g_jk_kwargs["is_normalizing_rand_sample"] = True
    #
    # Is only needed to reproduce old results
    # Possible choice: "v1" (also need default_g_jk_kwargs["is_normalizing_rand_sample"] == False)
    g_jk_kwargs["is_use_old_rand_alg"] = False
    #
    # these parameters are used in jk_blocking_func_default
    g_jk_kwargs["block_size"] = 1
    g_jk_kwargs["block_size_dict"] = {
            "job_tag": 1,
            }
    g_jk_kwargs["all_jk_idx_set"] = set()
    #
    # jk_blocking_func(jk_idx) => blocked jk_idx
    g_jk_kwargs["jk_blocking_func"] = jk_blocking_func_default
    #
    return g_jk_kwargs

@use_kwargs(default_g_jk_kwargs)
def get_jk_state(
        *,
        jk_type,
        eps,
        n_rand_sample,
        is_normalizing_rand_sample,
        is_use_old_rand_alg,
        block_size,
        block_size_dict,
        **_kwargs,
        ):
    """
    Currently only useful if we set
    #
    q.default_g_jk_kwargs["jk_type"] = "rjk" # this not yet the default
    #
    and
    #
    q.default_g_jk_kwargs["jk_blocking_func"] = jk_blocking_func_default
    #
    """
    assert jk_type == "rjk"
    return (
            jk_type,
            eps,
            n_rand_sample,
            is_normalizing_rand_sample,
            is_use_old_rand_alg,
            block_size,
            block_size_dict,
            )

@use_kwargs(default_g_jk_kwargs)
def jk_blocking_func_default(
        jk_idx,
        *,
        block_size,
        block_size_dict,
        all_jk_idx_set,
        **_kwargs,
        ):
    """
    block_size_for_this_job_tag = block_size_dict.get(job_tag, block_size)
    #
    use default_g_jk_kwargs for
    block_size, block_size_dict, all_jk_idx_set
    """
    block_size = default_g_jk_kwargs["block_size"]
    block_size_dict = default_g_jk_kwargs["block_size_dict"]
    all_jk_idx_set = default_g_jk_kwargs["all_jk_idx_set"]
    if block_size_dict is None:
        block_size_dict = dict()
    if all_jk_idx_set is not None:
        all_jk_idx_set.add(jk_idx)
    if isinstance(jk_idx, int_types):
        traj = jk_idx
        return traj // block_size
    elif isinstance(jk_idx, tuple) and len(jk_idx) == 2 and isinstance(jk_idx[1], int_types):
        job_tag, traj = jk_idx
        assert isinstance(job_tag, str)
        assert isinstance(traj, int_types)
        block_size_for_this_job_tag = block_size_dict.get(job_tag, block_size)
        assert isinstance(block_size_for_this_job_tag, int_types)
        return (job_tag, traj // block_size_for_this_job_tag,)
    else:
        return jk_idx
    assert False

@use_kwargs(default_g_jk_kwargs)
@timer
def g_jk(data_list, *, eps, **_kwargs):
    """
    Perform initial Jackknife for the original data set.\n
    """
    return jackknife(data_list, eps)

@use_kwargs(default_g_jk_kwargs)
@timer
def g_rejk(jk_list, jk_idx_list, *,
        jk_type,
        all_jk_idx,
        get_all_jk_idx,
        n_rand_sample,
        rng_state,
        jk_blocking_func,
        is_normalizing_rand_sample,
        is_use_old_rand_alg,
        **_kwargs,):
    """
    Perform (randomized) Super-Jackknife for the Jackknife data set.
    :jk_list: usually the Jackknife data set obtained with ``g_jk(data_list)``.
    :jk_idx_list: should be list of indices that names the ``jk_list``.
    :jk_type: [ "rjk", "super", ]``\n
    :returns: (randomized) Super-Jackknife data set.
    Note that::\n
        len(jk_list) == len(jk_idx_list)
        jk_idx_list[0] == "avg"
    Example
    """
    if jk_type == "super":
        if jk_blocking_func is not None:
            displayln_info(f"g_rejk: jk_type={jk_type} does not support jk_blocking_func={jk_blocking_func}")
        if all_jk_idx is None:
            assert get_all_jk_idx is not None
            all_jk_idx = get_all_jk_idx()
        return rejk_list(jk_list, jk_idx_list, all_jk_idx)
    elif jk_type == "rjk":
        return rjk_jk_list(jk_list, jk_idx_list, n_rand_sample, rng_state, jk_blocking_func, is_normalizing_rand_sample, is_use_old_rand_alg)
    else:
        assert False
    return None

@use_kwargs(default_g_jk_kwargs)
@timer
def g_mk_jk_val(rs_tag, val, err, *, jk_type, n_rand_sample, rng_state, is_normalizing_rand_sample, **_kwargs):
    """
    Create a jackknife sample with random numbers based on central value ``val`` and error ``err``.\n
    Need::\n
        default_g_jk_kwargs["jk_type"] = "rjk"
        default_g_jk_kwargs["n_rand_sample"] = n_rand_sample
        # e.g. n_rand_sample = 1024
        default_g_jk_kwargs["rng_state"] = rng_state
        # e.g. rng_state = RngState("rejk")
        default_g_jk_kwargs["is_normalizing_rand_sample"] = is_normalizing_rand_sample
        # e.g. is_normalizing_rand_sample = True
    """
    assert jk_type == "rjk"
    return rjk_mk_jk_val(rs_tag, val, err, n_rand_sample, rng_state, is_normalizing_rand_sample)

def g_jk_avg(jk_list, **_kwargs):
    """
    Return ``avg`` of the ``jk_list``.
    """
    if isinstance(jk_list, number_types):
        return jk_list
    return jk_avg(jk_list)

@use_kwargs(default_g_jk_kwargs)
def g_jk_err(jk_list, *, eps, jk_type, **_kwargs):
    """
    Return ``err`` of the ``jk_list``.
    """
    if isinstance(jk_list, number_types):
        return 0
    if jk_type == "super":
        return jk_err(jk_list, eps)
    elif jk_type == "rjk":
        return rjk_err(jk_list, eps)
    else:
        assert False
    return None

@timer
def g_jk_avg_err(jk_list, **kwargs):
    """
    Return ``(avg, err,)`` of the ``jk_list``.
    """
    return g_jk_avg(jk_list), g_jk_err(jk_list, **kwargs)

@timer
def g_jk_avg_err_arr(jk_list, **kwargs):
    avg, err = g_jk_avg_err(jk_list, **kwargs)
    avg_err_arr = np.stack([ avg, err, ])
    avg_err_arr = np.moveaxis(avg_err_arr, 0, -1).copy()
    return avg_err_arr

@use_kwargs(default_g_jk_kwargs)
def g_jk_size(**kwargs):
    """
    Return number of samples for the (randomized) Super-Jackknife data set.
    """
    jk_type = kwargs["jk_type"]
    all_jk_idx = kwargs["all_jk_idx"]
    get_all_jk_idx = kwargs["get_all_jk_idx"]
    n_rand_sample = kwargs["n_rand_sample"]
    # jk_type in [ "rjk", "super", ]
    if jk_type == "super":
        if all_jk_idx is None:
            assert get_all_jk_idx is not None
            all_jk_idx = get_all_jk_idx()
        return len(all_jk_idx)
    elif jk_type == "rjk":
        return 1 + n_rand_sample
    else:
        assert False
    return None

@use_kwargs(default_g_jk_kwargs)
def g_jk_blocking_func(jk_idx, *, jk_blocking_func, **_kwargs):
    """
    Return ``jk_blocking_func(jk_idx)``.
    """
    if jk_blocking_func is None:
        return jk_idx
    else:
        return jk_blocking_func(jk_idx)

@use_kwargs(default_g_jk_kwargs)
def g_jk_sample_size(job_tag, traj_list, **kwargs):
    jk_idx_list = [ (job_tag, traj,) for traj in traj_list ]
    b_jk_idx_set = set( g_jk_blocking_func(jk_idx, **kwargs) for jk_idx in jk_idx_list )
    return len(b_jk_idx_set)

default_g_jk_kwargs.update(mk_g_jk_kwargs())

class JkKwargs:

    """
    Example:
    #
    with q.JkKwargs(n_rand_sample=1024, block_size=10, block_size_dict={ "48I": 20, }):
        ...
    #
    """

    def __init__(self, **kwargs):
        self.new_kwargs = kwargs
        self.original = dict()

    def __enter__(self):
        for key, value in self.new_kwargs.items():
            self.original[key] = default_g_jk_kwargs[key]
            default_g_jk_kwargs[key] = self.new_kwargs[key]

    def __exit__(self, exc_type, exc_value, traceback):
        assert exc_type is None
        assert exc_value is None
        assert traceback is None
        for key, value in self.new_kwargs.items():
            default_g_jk_kwargs[key] = self.original[key]
        self.new_kwargs = None
        self.original = None

# ----

# ---- old funcs

def interpolate_list(data_arr, i):
    """
    Old function.
    return approximately v[i]
    Use `q.interp(data_arr, i, 0)` instead
    """
    return interp(v, i, 0)

def mk_jk_blocking_func(block_size=1, block_size_dict=None, all_jk_idx_set=None):
    """
    Recommend to use `jk_blocking_func_default` instead.
    #
    block_size_for_this_job_tag = block_size_dict.get(job_tag, block_size)
    """
    if block_size_dict is None:
        block_size_dict = dict()
    def jk_blocking_func(jk_idx):
        if all_jk_idx_set is not None:
            all_jk_idx_set.add(jk_idx)
        if isinstance(jk_idx, int_types):
            traj = jk_idx
            return traj // block_size
        elif isinstance(jk_idx, tuple) and len(jk_idx) == 2 and isinstance(jk_idx[1], int_types):
            job_tag, traj = jk_idx
            assert isinstance(job_tag, str)
            assert isinstance(traj, int_types)
            block_size_for_this_job_tag = block_size_dict.get(job_tag, block_size)
            assert isinstance(block_size_for_this_job_tag, int_types)
            return (job_tag, traj // block_size_for_this_job_tag,)
        else:
            return jk_idx
    return jk_blocking_func

def interpolate(data_arr, i_arr):
    """
    Old function. Use `q.interp(data_arr, i_arr, -1)` instead.
    #
    return approximately data_arr[..., i_arr]
    """
    vt = data_arr.transpose()
    if isinstance(i_arr, real_types):
        return interpolate_list(vt, i_arr).transpose()
    else:
        return np.array([ interpolate_list(vt, i) for i in i_arr ], data_arr.dtype).transpose()

def add_jk_idx(arr):
    """
    arr: no jk index
    return: add trivial jk index in the LAST axis
    """
    return arr.reshape(arr.shape + (1,))

def jk_transpose(arr):
    """
    arr: jk index is the 0th axis
    return: jk index is the last axis
    """
    shape = arr.shape
    ndim = len(shape)
    if ndim <= 1:
        return arr
    axes = list(range(1, ndim)) + [ 0, ]
    return arr.transpose(axes)

def jk_transpose_back(arr):
    """
    jk_transpose_back(jk_transpose(arr)) == arr
    """
    shape = arr.shape
    ndim = len(shape)
    if ndim <= 1:
        return arr
    axes = [ ndim - 1, ] + list(range(0, ndim - 1))
    return arr.transpose(axes)
