import math
import copy
import functools
import numpy as np


class q:
    from qlat_utils.utils import (
        get_fname,
    )
    from qlat_utils.timer import (
        timer,
        displayln_info,
    )
    from qlat_utils.rng_state import (
        RngState,
    )


alpha_qed = 1.0 / 137.035999084
fminv_gev = 0.197326979  # hbar * c / (1e-15 m * 1e9 electron charge * 1 volt)

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

    def __init__(self, default_kwargs, keys=None):
        """
        If ``keys`` is specified, then only the specified keys will be passed to the underlying function.
        """
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
                kwargs = {k: kwargs[k] for k in self.keys}
            return func(*args, **kwargs)
        return f

###


def interp_i_arr(data_x_arr, x_arr):
    r"""
    return ``i_arr``
    ``
    q.interp(data_x_arr, i_arr) \approx x_arr
    ``
    ``x_arr`` can be either an 1-D array-like object or a single number.
    e.g.:
    ``
    data(x)
    data_arr[:] = data(data_x_arr)
    q.interp(data_arr, i_arr) \approx data(x_arr)
    ``
    """
    data_i_arr = np.arange(len(data_x_arr))
    i_arr = np.interp(x_arr, data_x_arr, data_i_arr)
    return i_arr


def interp(data_arr, i_arr, axis=-1):
    """
    return approximately ``data_arr[..., i_arr]`` if ``axis=-1``.
    Note that ``i_arr`` can be non-integer.
    ``i_arr`` can be either an 1-D array-like object or a single number.
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
        iv_arr = np.array([interp(v_arr, i, 0) for i in i_arr], v_arr.dtype)
        iv_arr = np.moveaxis(iv_arr, 0, axis)
        return iv_arr
    else:
        fname = q.get_fname()
        raise Exception(f"{fname}: i_arr={i_arr}")


def interp_x(data_arr, data_x_arr, x_arr, axis=-1):
    """
    return ``interpolated_data_arr``
    ``x_arr`` can be either an 1-D array-like object or a single number.
    ``data_x_arr`` is the x values for ``data_arr``
    ``x_arr`` is the x values for `interpolated_data_arr`
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
    return ``x``
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
    return ``i_arr``
    #
    let ``shape`` = ``np.moveaxis(data_arr, axis, -1)[..., 0].shape``
    ``
    threshold_arr = np.broadcast_to(threshold_arr, shape)
    ``
    such that
    ``
    for index in np.ndindex(shape):
        q.interp(data_arr[index], i_arr[index]) \approx threshold_arr[index]
    ``
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
    ``
    data_x_arr.shape == (data_arr.shape[axis],)
    ``
    let `shape` = `np.moveaxis(data_arr, axis, -1)[..., 0].shape`
    ``
    threshold_arr = np.broadcast_to(threshold_arr, shape)
    ``
    such that
    ``
    for index in np.ndindex(shape):
        q.interp_x(data_arr[index], data_x_arr, x_arr[index]) \approx threshold_arr[index]
    ``
    """
    assert data_x_arr.shape == (data_arr.shape[axis],)
    i_arr = get_threshold_i_arr(data_arr, threshold_arr, axis)
    x_arr = np.zeros(i_arr.shape, dtype=np.float64)
    x_arr.ravel()[:] = interp(data_x_arr, i_arr.ravel())
    return x_arr


def partial_sum_list(x, *, is_half_last=False):
    """Modify in-place, preserve length"""
    s = 0
    for i, v in enumerate(x):
        sp = s
        s += v
        if is_half_last:
            x[i] = (s + sp) / 2
        else:
            x[i] = s


def partial_sum(x, *, is_half_last=False):
    """Modify in-place, preserve length"""
    shape = x.shape
    if len(shape) == 0:
        return
    elif len(shape) == 1:
        partial_sum_list(x, is_half_last=is_half_last)
    elif len(shape) == 2:
        for v in x:
            partial_sum_list(v, is_half_last=is_half_last)
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
        return sum([qnorm(x_i) for x_i in x])
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
                return Data([v1 + v2 for v1, v2 in zip(self.val, other.val)])
            elif isinstance(self.val, list):
                return Data([v + other.val for v in self.val])
            elif isinstance(other.val, list):
                return Data([self.val + v for v in other.val])
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
                return Data([v1 * v2 for v1, v2 in zip(self.val, other.val)])
            elif isinstance(self.val, list):
                return Data([v * other.val for v in self.val])
            elif isinstance(other.val, list):
                return Data([self.val * v for v in other.val])
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
            return Data([-v for v in self.val])
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
                return Data([v1 - v2 for v1, v2 in zip(self.val, other.val)])
            elif isinstance(self.val, list):
                return Data([v - other.val for v in self.val])
            elif isinstance(other.val, list):
                return Data([self.val - v for v in other.val])
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


def filter_np_results(val):
    if not hasattr(val, "size"):
        return val
    if val.size != 1:
        return val
    if not hasattr(val, "item"):
        return val
    return val.item()


def average(data_list):
    n = len(data_list)
    v = sum(data_list)
    return filter_np_results(1/n * v)


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
        return [average(data_list), ]
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


def avg_err(data_list, *, eps=1, block_size=1):
    """
    Compute `(avg, err)` of `data_list`.
        --
    :param data_list: list of data
    :param eps: additional scaling factor for the error
    :param block_size: blocking the list of data
    :return: (avg, err,) where avg and err have the same type as data
    :rtype: (avg, err,)
    """
    assert block_size >= 1
    avg = average(data_list)
    n = len(data_list)
    if n <= 1:
        err = abs(eps) * avg
        err = filter_np_results(err)
        return (avg, err,)
    if n < 2 * block_size:
        block_size = 1
    assert n > block_size
    blocks = block_data(data_list, block_size)
    diff_sqr = average([fsqr(d - avg) for d in blocks])
    fac = abs(eps) * math.sqrt(block_size / (n - block_size))
    err = fac * fsqrt(diff_sqr)
    err = filter_np_results(err)
    return (avg, err,)


def jackknife(data_list, *, eps=1):
    r"""
    Return jk[i] = avg - \frac{eps}{N} (v[i] - avg)
    normal jackknife uses eps=1, scale the fluctuation by eps
    """
    is_np_arr = isinstance(data_list, np.ndarray)
    data_list_real = [d for d in data_list if d is not None]
    n = len(data_list_real)
    fac = eps / n
    avg = average(data_list_real)
    jks = [avg, ]
    for data in data_list:
        if data is None:
            jks.append(avg)
        else:
            jks.append(avg - fac * (data - avg))
    if is_np_arr:
        jks = np.array(jks, dtype=data_list.dtype)
    return jks


def fsqr(data):
    """
    Separately square real and imag part in case of complex types.
        --
    :param data: real, complex, np.ndarray like objects.
    :return: squared `data`.
    :rtype: same type as `data`.
    """
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
    """
    Separately calculate the square root real and imag part in case of complex types.
        --
    :param data: real, complex, np.ndarray like objects.
    :return: squared `data`.
    :rtype: same type as `data`.
    """
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
    err_sqr = sum([fsqr(v) for v in vs])
    err = fsqrt(err_sqr)
    return err


def jk_avg(jk_arr):
    val = jk_arr[0]
    return filter_np_results(val)


def jk_err(jk_arr, *, eps=1, block_size=1):
    r"""
    Return
    $$
    \frac{1}{eps} \sqrt{ N/(N-block_size) \sum_{i=1}^N (jk[i] - jk_avg)^2 }.
    $$
    when ``block_size=1``.
    Note: ``len(jk_arr) = N + 1``.
    Same ``eps`` as the ``eps`` used in the ``jackknife`` function.
    Does not properly honor the $(N-1)$ formula in error calculation
    if there were missing data in the original ``data_list`` in the ``jackknife`` function.
    """
    assert block_size >= 1
    avg = jk_avg(jk_arr)
    n = len(jk_arr) - 1
    if n <= 1:
        fac = 1 / abs(eps)
        val = fac * avg
        val = filter_np_results(val)
        return val
    if n < 2 * block_size:
        block_size = 1
    assert n > block_size
    blocks = block_data(jk_arr[1:], block_size)
    diff_sqr = average([fsqr(jk - avg) for jk in blocks])
    fac = math.sqrt(block_size / (n - block_size)) * n / abs(eps)
    val = fac * fsqrt(diff_sqr)
    val = filter_np_results(val)
    return val


def jk_avg_err(jk_arr, *, eps=1, block_size=1):
    return jk_avg(jk_arr), jk_err(jk_arr, eps=eps, block_size=block_size)


@q.timer
def sjackknife(
    data_list,
    jk_idx_list,
    *,
    avg=None,
    is_hash_jk_idx=True,
    jk_idx_hash_size=None,
    rng_state=None,
    all_jk_idx=None,
    get_all_jk_idx=None,
    jk_blocking_func=None,
    eps=1,
    ):
    """
    Super jackknife.
    Return ``jk_arr``.
    ``len(jk_idx_list) == len(data_list)``
    ``len(jk_arr) == len(all_jk_idx)``
    ``jk_idx_list`` (after processed by ``jk_blocking_func``) should be contained in ``all_jk_idx``,
    otherwise, if ``is_hash_jk_idx`` is true, then, a hash of ``jk_idx`` will be used instead.
    Ideally, ``all_jk_idx`` should only contain distinct indices.
    However, if there are repeatations, indices appear later take precedence.
    if ``all_jk_idx`` and ``get_all_jk_idx`` are both ``None``, then a trivial
    ``all_jk_idx`` will be created based on ``jk_idx_hash_size``.
    """
    if jk_idx_hash_size is None:
        jk_idx_hash_size = 1024
    if rng_state is None:
        rng_state = q.RngState("rejk")
    rs = rng_state
    assert len(jk_idx_list) == len(data_list)
    if isinstance(data_list, np.ndarray):
        dtype = data_list.dtype
    else:
        dtype = None
    data_list_real = [d for d in data_list if d is not None]
    data_arr = np.array(data_list_real, dtype=dtype)
    if avg is None:
        avg = average(data_arr)
    dtype = data_arr.dtype
    jk_idx_list = [
        jk_idx
        for jk_idx, d in zip(jk_idx_list, data_list)
        if d is not None
    ]
    if jk_blocking_func is not None:
        jk_idx_list = [
            jk_blocking_func(0, jk_idx)
            for jk_idx in jk_idx_list
        ]
    n = len(data_arr)
    assert n == len(jk_idx_list)
    if all_jk_idx is None:
        if get_all_jk_idx is None:
            assert is_hash_jk_idx
            all_jk_idx = ["avg", ] + list(range(jk_idx_hash_size))
        else:
            all_jk_idx = get_all_jk_idx()
    assert all_jk_idx[0] == "avg"
    n_super_sample = len(all_jk_idx) - 1
    i_dict = dict()
    for i, jk_idx in enumerate(all_jk_idx):
        jk_idx_str = str(jk_idx)
        i_dict[jk_idx_str] = i
    i_arr = np.zeros(n, dtype=np.int32)
    for j in range(n):
        jk_idx = jk_idx_list[j]
        jk_idx_str = str(jk_idx)
        if jk_idx_str in i_dict:
            i = i_dict[jk_idx_str]
        else:
            assert is_hash_jk_idx
            rsi = rs.split(jk_idx_str)
            i = 1 + int(rsi.rand_gen() % n_super_sample)
        assert i > 0
        i_arr[j] = i
    count_dict = dict()
    for j in range(n):
        i = i_arr[j]
        if i in count_dict:
            count_dict[i] += 1
        else:
            count_dict[i] = 1
    jk_arr = np.empty((1 + n_super_sample, *data_arr[0].shape,), dtype=dtype)
    jk_arr[:] = avg
    data_diff = data_arr - avg
    for j in range(n):
        i = i_arr[j]
        assert i > 0
        assert i in count_dict
        if n > count_dict[i]:
            n_b = n - count_dict[i]
            fac = -eps * np.sqrt(1 / (n * n_b))
            jk_arr[i] += fac * data_diff[j]
    return jk_arr


@q.timer
def sjk_mk_jk_val(
    rs_tag,
    val,
    err,
    *,
    is_hash_jk_idx=True,
    jk_idx_hash_size=None,
    rng_state=None,
    all_jk_idx=None,
    get_all_jk_idx=None,
    eps=1,
):
    """
    return jk_arr
    n = n_rand_sample
    len(jk_arr) == 1 + n
    jk_arr[i] = val + err * r[i] for i in 1..n
    where r[i] ~ N(0, 1)
    """
    if jk_idx_hash_size is None:
        jk_idx_hash_size = 1024
    assert jk_idx_hash_size >= 0
    if rng_state is None:
        rng_state = q.RngState("rejk")
    rs = rng_state
    if all_jk_idx is None:
        if get_all_jk_idx is None:
            assert is_hash_jk_idx
            all_jk_idx = ["avg", ] + list(range(jk_idx_hash_size))
        else:
            all_jk_idx = get_all_jk_idx()
    assert all_jk_idx[0] == "avg"
    n_super_sample = len(all_jk_idx) - 1
    assert n_super_sample >= 0
    assert isinstance(rng_state, q.RngState)
    assert isinstance(val, real_types)
    assert isinstance(err, real_types)
    rs = rng_state.split(str(rs_tag))
    jk_arr = np.zeros((n_super_sample + 1,), dtype=np.float64)
    jk_arr[0] = val
    r_arr = rs.g_rand_arr((n_super_sample,))
    r_arr_qnorm = qnorm(r_arr)
    r_arr = r_arr * np.sqrt(1 / r_arr_qnorm)
    assert abs(qnorm(r_arr) - 1) < 1e-8
    jk_arr[1:] = val + eps * r_arr * err
    return jk_arr


def sjk_avg(jk_arr):
    return jk_avg(jk_arr)


def sjk_err(jk_arr, *, eps=1):
    r"""
    Return
    $$
    \frac{1}{eps} \sqrt{ \sum_{i=1}^N (jk[i] - jk_avg)^2 }.
    $$
    Note: ``len(jk_arr) = N + 1``.
    Same ``eps`` as the ``eps`` used in the ``jackknife`` function.
    """
    avg = jk_avg(jk_arr)
    n = len(jk_arr) - 1
    if n <= 1:
        fac = 1 / abs(eps)
        val = fac * avg
        val = filter_np_results(val)
        return val
    diff_sqr = average([fsqr(jk - avg) for jk in jk_arr[1:]])
    fac = math.sqrt(n) / abs(eps)
    val = fac * fsqrt(diff_sqr)
    val = filter_np_results(val)
    return val


def sjk_avg_err(jk_arr, *, eps=1):
    return sjk_avg(jk_arr), sjk_err(jk_arr, eps=eps)


# ----------


@q.timer
def mk_r_i_j_mat(
    n_rand_sample,
    jk_idx_list,
    rng_state,
    *,
    jk_blocking_func,
    is_normalizing_rand_sample,
    is_apply_rand_sample_jk_idx_blocking_shift,
    is_use_old_rand_alg,
):
    assert n_rand_sample >= 0
    rs = rng_state
    n = len(jk_idx_list)
    r_arr = np.empty((n_rand_sample, n,), dtype=np.float64)
    b_arr = np.empty((n_rand_sample, n,), dtype=np.int32)
    jk_idx_str_arr = np.empty((n_rand_sample, n,), dtype=object)
    jk_idx_str_set = set()
    if jk_blocking_func is None:
        is_apply_rand_sample_jk_idx_blocking_shift = False

    @q.timer
    def set_jk_idx():
        if is_apply_rand_sample_jk_idx_blocking_shift:
            for i in range(n_rand_sample):
                count_dict = dict()
                for j in range(n):
                    jk_idx = jk_blocking_func(i + 1, jk_idx_list[j])
                    jk_idx_str = str(jk_idx)
                    jk_idx_str_arr[i, j] = jk_idx_str
                    jk_idx_str_set.add(jk_idx_str)
                    if jk_idx_str in count_dict:
                        count_dict[jk_idx_str] += 1
                    else:
                        count_dict[jk_idx_str] = 1
                for j in range(n):
                    jk_idx_str = jk_idx_str_arr[i, j]
                    b_arr[i, j] = count_dict[jk_idx_str]
        else:
            count_dict = dict()
            for j in range(n):
                jk_idx = jk_idx_list[j]
                if jk_blocking_func is not None:
                    jk_idx = jk_blocking_func(0, jk_idx)
                jk_idx_str = str(jk_idx)
                jk_idx_str_arr[:, j] = jk_idx_str
                jk_idx_str_set.add(jk_idx_str)
                if jk_idx_str in count_dict:
                    count_dict[jk_idx_str] += 1
                else:
                    count_dict[jk_idx_str] = 1
            for j in range(n):
                jk_idx_str = jk_idx_str_arr[0, j]
                b_arr[:, j] = count_dict[jk_idx_str]

    set_jk_idx()
    if is_use_old_rand_alg == "v1":
        assert not is_normalizing_rand_sample
        for i in range(n_rand_sample):
            rsi = rs.split(str(i))
            r = [rsi.split(jk_idx_str).g_rand_gen()
                 for jk_idx_str in jk_idx_str_arr[i]]
            for j in range(n):
                r_arr[i, j] = r[j]
        return r_arr, b_arr
    assert is_use_old_rand_alg == False
    r_arr_dict = dict()

    @q.timer
    def set_r():
        for jk_idx_str in jk_idx_str_set:
            rsi = rs.split(jk_idx_str)
            garr = rsi.g_rand_arr(n_rand_sample)
            if is_normalizing_rand_sample:
                # garr_qnorm \approx n_rand_sample
                garr_qnorm = qnorm(garr)
                garr = garr * np.sqrt(n_rand_sample / garr_qnorm)
                assert abs(qnorm(garr) / n_rand_sample - 1) < 1e-8
            r_arr_dict[jk_idx_str] = garr

    set_r()

    @q.timer
    def set_r_arr():
        if is_apply_rand_sample_jk_idx_blocking_shift:
            for i in range(n_rand_sample):
                for j in range(n):
                    jk_idx_str = jk_idx_str_arr[i, j]
                    garr = r_arr_dict[jk_idx_str]
                    r_arr[i, j] = garr[i]
        else:
            for j in range(n):
                jk_idx_str = jk_idx_str_arr[0, j]
                garr = r_arr_dict[jk_idx_str]
                r_arr[:, j] = garr

    set_r_arr()
    return r_arr, b_arr


@q.timer
def rjackknife(
    data_list,
    jk_idx_list,
    *,
    avg=None,
    rng_state=None,
    n_rand_sample=None,
    jk_blocking_func=None,
    is_normalizing_rand_sample=False,
    is_apply_rand_sample_jk_idx_blocking_shift=True,
    is_use_old_rand_alg=False,
    eps=1,
):
    r"""
    Jackknife-bootstrap hybrid resampling.
    Return ``jk_arr``.
    ``len(jk_arr) == 1 + n_rand_sample``
    distribution of ``jk_arr`` should be similar as the distribution of ``avg``.
    ``r_{i,j} ~ N(0, 1)``
    ``
    if is_normalizing_rand_sample:
        n_j = \sum_i r_{i,j}^2
        r_{i,j} <- \sqrt{n_rand_sample / n_j} r_{i,j}
    data_list_real = [d for d in data_list if d is not None]
    data_arr = np.array(data_list_real, dtype=dtype)
    avg = average(data_arr)
    len(data_list_real) = n
    jk_arr[0] = avg
    jk_arr[i] = avg + \sum_{j=1}^{n} (-eps/\sqrt{n (n - b(i,j))}) r_{i,j} (data_list_real[j] - avg)
    ``
    where ``b(i,j)`` represent the ``block_size``.
    if ``jk_blocking_func`` is provided:
        ``jk_blocking_func(i, jk_idx) => blocked jk_idx``
    ``
    jk_arr[i] = avg + \sum_{j=1}^{n} r_{i,jk_block_func(j)} (jk_arr[j] - avg)
    ``
    """
    if n_rand_sample is None:
        n_rand_sample = 1024
    if rng_state is None:
        rng_state = q.RngState("rejk")
    assert len(data_list) == len(jk_idx_list)
    assert isinstance(n_rand_sample, int_types)
    assert n_rand_sample >= 0
    assert isinstance(rng_state, q.RngState)
    if isinstance(data_list, np.ndarray):
        dtype = data_list.dtype
    else:
        dtype = None
    data_list_real = [d for d in data_list if d is not None]
    data_arr = np.array(data_list_real, dtype=dtype)
    if avg is None:
        avg = average(data_arr)
    dtype = data_arr.dtype
    jk_idx_list = [
        jk_idx
        for jk_idx, d in zip(jk_idx_list, data_list)
        if d is not None
    ]
    n = len(data_arr)
    r_arr, b_arr = mk_r_i_j_mat(
        n_rand_sample,
        jk_idx_list,
        rng_state,
        jk_blocking_func=jk_blocking_func,
        is_normalizing_rand_sample=is_normalizing_rand_sample,
        is_apply_rand_sample_jk_idx_blocking_shift=is_apply_rand_sample_jk_idx_blocking_shift,
        is_use_old_rand_alg=is_use_old_rand_alg,
    )
    n_b_arr = n - b_arr
    n_b_arr[n <= b_arr] = 1
    fac_arr = -eps / np.sqrt(n * n_b_arr)
    fac_arr[n <= b_arr] = 0
    fac_r_arr = fac_arr * r_arr
    pad_shape = (1,) * len(data_arr[0].shape)
    fac_r_arr = fac_r_arr.reshape(fac_r_arr.shape + pad_shape)
    data_diff = data_arr - avg
    jk_arr = np.empty((1 + n_rand_sample, *data_arr[0].shape,), dtype=dtype)
    jk_arr[0] = avg
    jk_arr[1:] = avg + np.sum(fac_r_arr * data_diff, axis=1)
    return jk_arr


@q.timer
def rjk_mk_jk_val(
    rs_tag,
    val,
    err,
    *,
    n_rand_sample=None,
    rng_state=None,
    eps=1,
):
    """
    return jk_arr
    n = n_rand_sample
    len(jk_arr) == 1 + n
    jk_arr[i] = val + err * r[i] for i in 1..n
    where r[i] ~ N(0, 1)
    """
    if n_rand_sample is None:
        n_rand_sample = 1024
    if rng_state is None:
        rng_state = q.RngState("rejk")
    assert n_rand_sample >= 0
    assert isinstance(rng_state, q.RngState)
    assert isinstance(val, real_types)
    assert isinstance(err, real_types)
    rs = rng_state.split(str(rs_tag))
    jk_arr = np.zeros((n_rand_sample + 1,), dtype=np.float64)
    jk_arr[0] = val
    r_arr = rs.g_rand_arr((n_rand_sample,))
    r_arr_qnorm = qnorm(r_arr)
    r_arr = r_arr * np.sqrt(n_rand_sample / r_arr_qnorm)
    assert abs(qnorm(r_arr) / n_rand_sample - 1) < 1e-8
    jk_arr[1:] = val + eps * r_arr * err
    return jk_arr


def rjk_avg(jk_arr):
    return jk_avg(jk_arr)


def rjk_err(jk_arr, eps=1):
    r"""
    Return
    $$
    \frac{1}{eps} \sqrt{ 1/N \sum_{i=1}^N (jk[i] - jk_avg)^2 }.
    $$
    Note: ``
    len(jk_arr) = N + 1.
    jk_avg = jk_arr[0]
    ``
    Same ``eps`` as the ``eps`` used in the ``jackknife`` function.
    """
    avg = jk_avg(jk_arr)
    n = len(jk_arr) - 1
    if n <= 0:
        fac = 1 / abs(eps)
        val = fac * avg
        val = filter_np_results(val)
        return val
    diff_sqr = average([fsqr(jk - avg) for jk in jk_arr[1:]])
    fac = 1 / abs(eps)
    val = fac * fsqrt(diff_sqr)
    val = filter_np_results(val)
    return val


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
    # for jk_type = "rjk"
    g_jk_kwargs["n_rand_sample"] = 1024
    g_jk_kwargs["is_normalizing_rand_sample"] = False
    g_jk_kwargs["is_apply_rand_sample_jk_idx_blocking_shift"] = True
    #
    # for jk_type = "super"
    g_jk_kwargs["is_hash_jk_idx"] = True
    g_jk_kwargs["jk_idx_hash_size"] = 1024
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
    #
    # Below are items which are not touched in
    # ``get_jk_state`` or ``set_jk_state``
    #
    g_jk_kwargs["rng_state"] = q.RngState("rejk")
    #
    g_jk_kwargs["all_jk_idx"] = None
    g_jk_kwargs["get_all_jk_idx"] = None
    #
    g_jk_kwargs["all_jk_idx_set"] = set()
    #
    # jk_blocking_func(i, jk_idx) => blocked_jk_idx
    g_jk_kwargs["jk_blocking_func"] = jk_blocking_func_default
    #
    return g_jk_kwargs


def reset_default_g_jk_kwargs():
    default_g_jk_kwargs.clear()
    default_g_jk_kwargs.update(mk_g_jk_kwargs())


@use_kwargs(default_g_jk_kwargs)
def get_jk_state(
    *,
    jk_type,
    eps,
    n_rand_sample,
    is_normalizing_rand_sample,
    is_apply_rand_sample_jk_idx_blocking_shift,
    is_hash_jk_idx,
    jk_idx_hash_size,
    is_use_old_rand_alg,
    block_size,
    block_size_dict,
    **_kwargs,
):
    """
    Currently only useful if we set
    #
    q.default_g_jk_kwargs["jk_type"] = "rjk" # this is the default now
    and
    q.default_g_jk_kwargs["jk_blocking_func"] = jk_blocking_func_default
    #
    Used for `q.cache_call`.
    Example:
    q.cache_call(get_state=q.get_jk_state)
    def func(...):
        ...
    """
    assert jk_type == "rjk"
    return (
        jk_type,
        eps,
        n_rand_sample,
        is_normalizing_rand_sample,
        is_apply_rand_sample_jk_idx_blocking_shift,
        is_hash_jk_idx,
        jk_idx_hash_size,
        is_use_old_rand_alg,
        block_size,
        block_size_dict,
    )


def set_jk_state(state):
    (
        jk_type,
        eps,
        n_rand_sample,
        is_normalizing_rand_sample,
        is_apply_rand_sample_jk_idx_blocking_shift,
        is_hash_jk_idx,
        jk_idx_hash_size,
        is_use_old_rand_alg,
        block_size,
        block_size_dict,
    ) = state
    g_dict = default_g_jk_kwargs
    g_dict["jk_type"] = jk_type
    g_dict["eps"] = eps
    g_dict["n_rand_sample"] = n_rand_sample
    g_dict["is_normalizing_rand_sample"] = is_normalizing_rand_sample
    g_dict["is_apply_rand_sample_jk_idx_blocking_shift"] = is_apply_rand_sample_jk_idx_blocking_shift
    g_dict["is_hash_jk_idx"] = is_hash_jk_idx
    g_dict["jk_idx_hash_size"] = jk_idx_hash_size
    g_dict["is_use_old_rand_alg"] = is_use_old_rand_alg
    g_dict["block_size"] = block_size
    g_dict["block_size_dict"] = block_size_dict


jk_blocking_traj_shift_arr = (
    q.RngState("jk_blocking_traj_shift_arr").rand_arr(16 * 1024)
    %
    (1024 * 1024 * 1024 * 1024)
)


@use_kwargs(default_g_jk_kwargs)
def jk_blocking_func_default(
    i,
    jk_idx,
    *,
    block_size,
    block_size_dict,
    all_jk_idx_set,
    **_kwargs,
):
    """
    return ``blocked_jk_idx``.
    ``blocked_jk_idx`` should uniquely identify the block
    that configuration identified by ``jk_idx`` belongs to.
    The block scheme can be different for different J-B hybrid sample.
    The J-B hybrid sample is indexed by ``i`` (``1 <= i <= n_rand_sample``).
    ``
    block_size_for_this_job_tag = block_size_dict.get(job_tag, block_size)
    ``
    use default_g_jk_kwargs for
    block_size, block_size_dict, all_jk_idx_set
    """
    if i == 0:
        shift = 0
    else:
        assert i >= 1
        shift = int(jk_blocking_traj_shift_arr[
            (i - 1) % len(jk_blocking_traj_shift_arr)
        ])
    if block_size_dict is None:
        block_size_dict = dict()
    if all_jk_idx_set is not None:
        all_jk_idx_set.add(jk_idx)
    if isinstance(jk_idx, int_types):
        traj = jk_idx
        b_shift = shift % block_size
        return (traj + b_shift) // block_size
    elif isinstance(jk_idx, tuple) and len(jk_idx) == 2 and isinstance(jk_idx[1], int_types):
        job_tag, traj = jk_idx
        assert isinstance(job_tag, str)
        assert isinstance(traj, int_types)
        block_size_for_this_job_tag = block_size_dict.get(job_tag, block_size)
        assert isinstance(block_size_for_this_job_tag, int_types)
        b_shift = shift % block_size_for_this_job_tag
        return (job_tag, (traj + b_shift) // block_size_for_this_job_tag,)
    else:
        return jk_idx
    assert False


@use_kwargs(default_g_jk_kwargs)
@q.timer
def g_mk_jk(
    data_list,
    jk_idx_list,
    *,
    avg=None,
    jk_type,
    all_jk_idx,
    get_all_jk_idx,
    n_rand_sample,
    rng_state,
    jk_blocking_func,
    is_normalizing_rand_sample,
    is_apply_rand_sample_jk_idx_blocking_shift,
    is_use_old_rand_alg,
    is_hash_jk_idx,
    jk_idx_hash_size,
    eps,
    **_kwargs,
):
    """
    Perform (randomized) Super-Jackknife for the Jackknife data set.
        --
    :param data_list: initial un-jackknifed data.
    :param jk_idx_list: should be list of indices that names the ``jk_arr``.
    :param jk_type: ``[ "rjk", "super", ]``
    :param eps: Error scaling factor.
    :return: (randomized) Super-Jackknife data set.
    Note that::
        len(data_list) == len(jk_idx_list)
        jk_idx_list = [(job_tag, traj,) for traj in traj_list]
    We can set ``eps`` to be factor ``len(data_list)`` larger.
    """
    if jk_type == "super":
        jk_arr = sjackknife(
            data_list,
            jk_idx_list,
            avg=avg,
            is_hash_jk_idx=is_hash_jk_idx,
            jk_idx_hash_size=jk_idx_hash_size,
            rng_state=rng_state,
            all_jk_idx=all_jk_idx,
            get_all_jk_idx=get_all_jk_idx,
            jk_blocking_func=jk_blocking_func,
            eps=eps,
        )
    elif jk_type == "rjk":
        jk_arr = rjackknife(
            data_list,
            jk_idx_list,
            avg=avg,
            n_rand_sample=n_rand_sample,
            rng_state=rng_state,
            jk_blocking_func=jk_blocking_func,
            is_normalizing_rand_sample=is_normalizing_rand_sample,
            is_apply_rand_sample_jk_idx_blocking_shift=is_apply_rand_sample_jk_idx_blocking_shift,
            is_use_old_rand_alg=is_use_old_rand_alg,
            eps=eps,
        )
    else:
        assert False
    return jk_arr


@use_kwargs(default_g_jk_kwargs)
@q.timer
def g_mk_jk_val(
    rs_tag,
    val,
    err,
    *,
    jk_type,
    all_jk_idx,
    get_all_jk_idx,
    n_rand_sample,
    rng_state,
    is_hash_jk_idx,
    jk_idx_hash_size,
    eps,
    **_kwargs,
):
    """
    Create a jackknife sample with random numbers based on central value ``val`` and error ``err``.
        --
    Need:
        default_g_jk_kwargs["jk_type"] = "rjk"
        default_g_jk_kwargs["n_rand_sample"] = n_rand_sample
        # e.g. n_rand_sample = 1024
        default_g_jk_kwargs["rng_state"] = rng_state
        # e.g. rng_state = q.RngState("rejk")
    """
    if jk_type == "super":
        jk_val = sjk_mk_jk_val(
            rs_tag,
            val,
            err,
            is_hash_jk_idx=is_hash_jk_idx,
            jk_idx_hash_size=jk_idx_hash_size,
            rng_state=rng_state,
            all_jk_idx=all_jk_idx,
            get_all_jk_idx=get_all_jk_idx,
            eps=eps,
        )
    elif jk_type == "rjk":
        jk_val = rjk_mk_jk_val(
            rs_tag,
            val,
            err,
            n_rand_sample=n_rand_sample,
            rng_state=rng_state,
            eps=eps,
        )
    else:
        assert False
    return jk_val


def g_jk_avg(jk_arr, **_kwargs):
    """
    Return ``avg`` of the ``jk_arr``.
    """
    if isinstance(jk_arr, number_types):
        return jk_arr
    return jk_avg(jk_arr)


@use_kwargs(default_g_jk_kwargs)
def g_jk_err(jk_arr, *, eps, jk_type, **_kwargs):
    """
    Return ``err`` of the ``jk_arr``.
    """
    if isinstance(jk_arr, number_types):
        return 0
    if jk_type == "super":
        return sjk_err(jk_arr, eps=eps)
    elif jk_type == "rjk":
        return rjk_err(jk_arr, eps=eps)
    else:
        assert False
    return None


@q.timer
def g_jk_avg_err(jk_arr, **kwargs):
    """
    Return ``(avg, err,)`` of the ``jk_arr``.
    """
    return g_jk_avg(jk_arr), g_jk_err(jk_arr, **kwargs)


@q.timer
def g_jk_avg_err_arr(jk_arr, **kwargs):
    """
    Return ``avg_err_arr`` of the ``jk_arr``.
    ``
    avg_err_arr.shape = jk_arr[0].shape + (2,)
    ``
    """
    avg, err = g_jk_avg_err(jk_arr, **kwargs)
    avg_err_arr = np.stack([avg, err, ])
    avg_err_arr = np.moveaxis(avg_err_arr, 0, -1).copy()
    return avg_err_arr


@use_kwargs(default_g_jk_kwargs)
def g_jk_size(
    *,
    jk_type,
    all_jk_idx,
    get_all_jk_idx,
    n_rand_sample,
    is_hash_jk_idx,
    jk_idx_hash_size,
    **_kwargs,
):
    """
    Return number of samples for the (randomized) Super-Jackknife data set.
    """
    if jk_type == "super":
        if all_jk_idx is None:
            if get_all_jk_idx is None:
                assert is_hash_jk_idx
                all_jk_idx = ["avg", ] + list(range(jk_idx_hash_size))
            else:
                all_jk_idx = get_all_jk_idx()
        assert all_jk_idx[0] == "avg"
        n_super_sample = len(all_jk_idx) - 1
        assert n_super_sample >= 0
        return 1 + n_super_sample
    elif jk_type == "rjk":
        return 1 + n_rand_sample
    else:
        assert False
    return None


@use_kwargs(default_g_jk_kwargs)
def g_jk_blocking_func(
    i, jk_idx,
    *,
    jk_blocking_func,
    **_kwargs,
):
    """
    Return ``jk_blocking_func(jk_idx)``.
    """
    if jk_blocking_func is None:
        return jk_idx
    else:
        return jk_blocking_func(i, jk_idx)


@use_kwargs(default_g_jk_kwargs)
def g_jk_sample_size(
    job_tag,
    traj_list,
    **_kwargs,
):
    jk_idx_list = [(job_tag, traj,) for traj in traj_list]
    b_jk_idx_set = set(g_jk_blocking_func(0, jk_idx, **kwargs)
                       for jk_idx in jk_idx_list)
    return len(b_jk_idx_set)


reset_default_g_jk_kwargs()

# ----

default_show_val_kwargs = dict()


def mk_show_val_kwargs():
    d = dict()
    d["is_latex"] = True
    d["num_float_digit"] = None
    d["num_exp_digit"] = None
    d["exponent"] = None
    return d


default_show_val_kwargs.update(mk_show_val_kwargs())


def get_val_exp(val, exp=0):
    """
    return val, exp
    where
    `val * 10**exp` is the same as input
    """
    assert isinstance(val, (int, float))
    assert isinstance(exp, int)
    if val == 0.0:
        return 0.0, 0
    while abs(val) >= 10.0:
        val /= 10
        exp += 1
    while abs(val) < 1.0:
        val *= 10
        exp -= 1
    return val, exp


@use_kwargs(default_show_val_kwargs)
def show_val(
    val,
    *,
    is_latex,
    num_float_digit,
    num_exp_digit,
    exponent,
):
    """
    `is_latex` can be in [ None, False, True, ]
    `num_float_digit` or `num_exp_digit` can be in [ None, False, True, int, ]
    `exponent` can be in [ None, int, ]
    """
    assert isinstance(val, (int, float))
    if is_latex is None:
        is_latex = True
    if exponent is not None:
        assert isinstance(exponent, int)
        num_float_digit = False
        assert num_exp_digit is not False
        e = exponent
        v = val / 10**e
    else:
        v, e = get_val_exp(val)
    if (num_float_digit is None) and (num_exp_digit is None):
        if -2 <= e <= 4:
            num_float_digit = True
            num_exp_digit = False
        else:
            num_exp_digit = True
            num_float_digit = False
    if num_float_digit is None:
        if num_exp_digit is False:
            num_float_digit = True
        else:
            num_float_digit = False
    if num_exp_digit is None:
        if num_float_digit is False:
            num_exp_digit = True
        else:
            num_exp_digit = False
    if num_float_digit is True:
        num_float_digit = max(1, 5 - e)
    else:
        assert (num_float_digit is False) or isinstance(num_float_digit, int)
    if num_exp_digit is True:
        num_exp_digit = 5
    else:
        assert (num_exp_digit is False) or isinstance(num_exp_digit, int)
    assert not ((num_float_digit is False)
                and (num_exp_digit is False))
    if num_exp_digit is False:
        assert isinstance(num_float_digit, int)
        assert num_float_digit >= 0
        return (f"{{:.{num_float_digit}f}}").format(val)
    else:
        assert isinstance(num_exp_digit, int)
        assert num_exp_digit >= 0
        v_str = (f"{{:.{num_exp_digit}f}}").format(v)
        if is_latex:
            return f"{v_str} \\times 10^{{{e}}}"
        else:
            return f"{v_str}E{e}"


@use_kwargs(default_show_val_kwargs)
def show_val_err(
    val_err,
    *,
    is_latex,
    num_float_digit,
    num_exp_digit,
    exponent,
):
    """
    `is_latex` can be in [ None, False, True, ]
    `num_float_digit` or `num_exp_digit` can be in [ None, False, True, int, ]
    `exponent` can be in [ None, int, ]
    #
    Examples:
    print(show_val_err((1.12e16, 12), num_float_digit=1))
    print(show_val_err((1.12e16, 12e6), num_exp_digit=True))
    print(show_val_err((1.12e16, 12e6)))
    print(show_val_err((1.12e16, 12e7), exponent=10))
    print(show_val_err((1.12e16, 12e7), exponent=10, is_latex=False))
    """
    if isinstance(val_err, (int, float)):
        val = val_err
        return show_val(
            val,
            is_latex=is_latex,
            num_float_digit=num_float_digit,
            num_exp_digit=num_exp_digit,
            exponent=exponent,
        )
    val, err = val_err
    if err == 0:
        return show_val(
            val,
            is_latex=is_latex,
            num_float_digit=num_float_digit,
            num_exp_digit=num_exp_digit,
            exponent=exponent,
        )
    assert isinstance(val, (int, float))
    assert isinstance(err, (int, float))
    if is_latex is None:
        is_latex = True
    e_v, e_e = get_val_exp(err)
    if abs(e_v) <= 2.5:
        e_v *= 100
        e_e -= 2
    else:
        e_v *= 10
        e_e -= 1
    if exponent is not None:
        assert isinstance(exponent, int)
        num_float_digit = False
        assert num_exp_digit is not False
        e = exponent
        v = val / 10**e
    else:
        v, e = get_val_exp(val)
        if e_e > e:
            e = e_e
            v = val / 10**e
    if (num_float_digit is None) and (num_exp_digit is None):
        if -2 <= e <= 4:
            num_float_digit = True
            num_exp_digit = False
        else:
            num_exp_digit = True
            num_float_digit = False
    if num_float_digit is None:
        if num_exp_digit is False:
            num_float_digit = True
        else:
            num_float_digit = False
    if num_exp_digit is None:
        if num_float_digit is False:
            num_exp_digit = True
        else:
            num_exp_digit = False
    if num_float_digit is True:
        num_float_digit = max(0, -e_e)
    else:
        assert (num_float_digit is False) or isinstance(num_float_digit, int)
    if num_exp_digit is True:
        num_exp_digit = max(0, e - e_e)
    else:
        assert (num_exp_digit is False) or isinstance(num_exp_digit, int)
    assert not ((num_float_digit is False)
                and (num_exp_digit is False))
    if num_exp_digit is False:
        assert isinstance(num_float_digit, int)
        assert num_float_digit >= 0
        if abs(err) >= 1.0:
            return (f"{{0:.{num_float_digit}f}}({{1:.{num_float_digit}f}})").format(val, err)
        else:
            e_e = -num_float_digit
            e_v = err / 10**e_e
            return (f"{{0:.{num_float_digit}f}}({{1}})").format(val, round(e_v))
    else:
        assert isinstance(num_exp_digit, int)
        assert num_exp_digit >= 0
        e_e = e
        e_v = err / 10**e_e
        if abs(e_v) >= 1.0:
            v_str = (f"{{0:.{num_exp_digit}f}}({{1:.{num_exp_digit}f}})").format(
                v, e_v)
        else:
            e_e = e - num_exp_digit
            e_v = err / 10**e_e
            v_str = (f"{{0:.{num_exp_digit}f}}({{1}})").format(v, round(e_v))
        if is_latex:
            return f"{v_str} \\times 10^{{{e}}}"
        else:
            return f"{v_str}E{e}"

# ----


class NewDictValues:

    """
    Example:
    #
    with q.NewDictValues(dictionary, k1=v1, k2=v2, ...):
        ...
    #
    """

    def __init__(self, dictionary, **kwargs):
        self.dictionary = dictionary
        self.new_kwargs = kwargs
        self.original = dict()

    def __enter__(self):
        for key in self.new_kwargs.keys():
            self.original[key] = self.dictionary[key]
            self.dictionary[key] = self.new_kwargs[key]

    def __exit__(self, exc_type, exc_value, traceback):
        assert exc_type is None
        assert exc_value is None
        assert traceback is None
        for key in self.new_kwargs.keys():
            self.dictionary[key] = self.original[key]
        self.new_kwargs = None
        self.original = None

# ----


class JkKwargs(NewDictValues):

    """
    Example:
    #
    with q.JkKwargs(n_rand_sample=1024, block_size=10, block_size_dict={ "48I": 20, }):
        ...
    #
    """

    def __init__(self, **kwargs):
        super().__init__(default_g_jk_kwargs, **kwargs)

# ----


class ShowKwargs(NewDictValues):

    """
    Example:
    #
    with q.ShowKwargs(is_latex=True, exponent=-10):
        ...
    #
    """

    def __init__(self, **kwargs):
        super().__init__(default_show_val_kwargs, **kwargs)

# ----

# ---- old funcs

def merge_jk_idx(*jk_idx_list_list):
    for jk_idx_list in jk_idx_list_list:
        assert jk_idx_list[0] == "avg"
    return ["avg", ] + [jk_idx for jk_idx_list in jk_idx_list_list for jk_idx in jk_idx_list[1:]]


@q.timer
def rejk_list(jk_list, jk_idx_list, all_jk_idx):
    """
    Super jackknife
    ``jk_idx_list`` should be contained in ``all_jk_idx`` and have the same order.
    Does not properly honor the (N-1) formula in error calculation.
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


@q.timer
def rjk_jk_list(
    jk_list,
    jk_idx_list,
    n_rand_sample,
    rng_state,
    jk_blocking_func=None,
    is_normalizing_rand_sample=False,
    is_apply_rand_sample_jk_idx_blocking_shift=True,
    is_use_old_rand_alg=False,
):
    r"""
    return jk_list
    len(jk_list) == 1 + n_rand_sample
    distribution of jk_list should be similar as the distribution of avg
    r_{i,j} ~ N(0, 1)
    if is_normalizing_rand_sample:
        n_j = \sum_i r_{i,j}^2
        r_{i,j} <- \sqrt{n_rand_sample / n_j} r_{i,j}
    avg = jk_list[0]
    len(jk_list) = n + 1
    jk_list[i] = avg + \sum_{j=1}^{n} r_{i,j} (jk_list[j] - avg)
    #
    if `jk_blocking_func` is provided:
    ``
    jk_blocking_func(i, jk_idx) => blocked jk_idx
    ``
    Note that: ``1 <= i <= n_rand_sample``
    ``
    jk_list[i] = avg + \sum_{j=1}^{n} r_{i,jk_block_func(i, j)} (jk_list[j] - avg)
    ``
    """
    assert jk_idx_list[0] == "avg"
    assert isinstance(n_rand_sample, int_types)
    assert n_rand_sample >= 0
    assert isinstance(rng_state, q.RngState)
    is_np_arr = isinstance(jk_list, np.ndarray)
    n = len(jk_list) - 1
    r_arr, b_arr = mk_r_i_j_mat(
        n_rand_sample, jk_idx_list[1:], rng_state,
        jk_blocking_func=jk_blocking_func,
        is_normalizing_rand_sample=is_normalizing_rand_sample,
        is_apply_rand_sample_jk_idx_blocking_shift=is_apply_rand_sample_jk_idx_blocking_shift,
        is_use_old_rand_alg=is_use_old_rand_alg,
    )
    avg = jk_list[0]
    if is_np_arr:
        jk_arr = jk_list
        jk_diff = jk_arr[1:] - avg
        rjk_arr = np.empty((1 + n_rand_sample, *avg.shape,),
                           dtype=jk_arr.dtype)
        rjk_arr[:] = avg
        for j in range(n):
            for i in range(n_rand_sample):
                rjk_arr[i + 1] += r_arr[i, j] * jk_diff[j]
        return rjk_arr
    else:
        rjk_list = [avg, ]
        jk_diff = [jk_list[j] - avg for j in range(1, n + 1)]
        for i in range(n_rand_sample):
            rjk_list.append(
                avg + sum([r_arr[i, j] * jk_diff[j] for j in range(n)]))
        return rjk_list


@use_kwargs(default_g_jk_kwargs)
@q.timer
def g_jk(data_list, *, eps, **_kwargs):
    """
    Obsolete, call ``g_mk_jk`` instead.
        --
    Perform initial Jackknife for the original data set.\n
    """
    return jackknife(data_list, eps=eps)


@use_kwargs(default_g_jk_kwargs)
@q.timer
def g_rejk(
    jk_list, jk_idx_list, *,
    jk_type,
    all_jk_idx,
    get_all_jk_idx,
    n_rand_sample,
    rng_state,
    jk_blocking_func,
    is_normalizing_rand_sample,
    is_apply_rand_sample_jk_idx_blocking_shift,
    is_use_old_rand_alg,
    **_kwargs,
):
    """
    Obsolete, call ``g_mk_jk`` instead.
        --
    Perform (randomized) Super-Jackknife for the Jackknife data set.
        --
    :jk_list: usually the Jackknife data set obtained with ``g_jk(data_list)``.
    :jk_idx_list: should be list of indices that names the ``jk_list``.
    :jk_type: ``[ "rjk", "super", ]``
    :returns: (randomized) Super-Jackknife data set.
    Note that::
        len(jk_list) == len(jk_idx_list)
        jk_idx_list[0] == "avg"
    """
    if jk_type == "super":
        if jk_blocking_func is not None:
            q.displayln_info(
                f"g_rejk: jk_type={jk_type} does not support jk_blocking_func={jk_blocking_func}")
        if all_jk_idx is None:
            assert get_all_jk_idx is not None
            all_jk_idx = get_all_jk_idx()
        return rejk_list(
            jk_list,
            jk_idx_list,
            all_jk_idx,
        )
    elif jk_type == "rjk":
        return rjk_jk_list(
            jk_list,
            jk_idx_list,
            n_rand_sample,
            rng_state,
            jk_blocking_func,
            is_normalizing_rand_sample,
            is_apply_rand_sample_jk_idx_blocking_shift,
            is_use_old_rand_alg,
        )
    else:
        assert False
    return None


def interpolate_list(data_arr, i):
    """
    Old function.
    return approximately data_arr[i]
    Use `q.interp(data_arr, i, 0)` instead
    """
    return interp(data_arr, i, 0)


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
            block_size_for_this_job_tag = block_size_dict.get(
                job_tag, block_size)
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
        return np.array([interpolate_list(vt, i) for i in i_arr], data_arr.dtype).transpose()


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
    axes = list(range(1, ndim)) + [0, ]
    return arr.transpose(axes)


def jk_transpose_back(arr):
    """
    jk_transpose_back(jk_transpose(arr)) == arr
    """
    shape = arr.shape
    ndim = len(shape)
    if ndim <= 1:
        return arr
    axes = [ndim - 1, ] + list(range(0, ndim - 1))
    return arr.transpose(axes)
