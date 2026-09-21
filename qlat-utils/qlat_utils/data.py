"""
Module ``qlat_utils.data``
==========================\n
Data analysis utilities: interpolation, basic statistics (averaging, blocking
and error estimation), value display, and the generic helpers used by the other
modules.  The jackknife resampling functions live in
``qlat_utils.jackknife_utils``.\n
Documentation: ``docs/qlat-utils/qlat_data.md``\n
.. note:: Update the documentation when updating this source file.
"""

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

float_types = (
    float,
    np.float32,
    np.float64,
)
complex_types = (
    complex,
    np.complex64,
    np.complex128,
)
int_types = (
    int,
    np.int32,
    np.int64,
)

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
        #
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
    ``x_arr`` is the x values for ``interpolated_data_arr``
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
    let ``shape`` = ``np.moveaxis(data_arr, axis, -1)[..., 0].shape``\n
    ::\n
        threshold_arr = np.broadcast_to(threshold_arr, shape)\n
    such that\n
    ::\n
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
    return x_arr\n
    ::\n
        data_x_arr.shape == (data_arr.shape[axis],)\n
    let ``shape`` = ``np.moveaxis(data_arr, axis, -1)[..., 0].shape``\n
    ::\n
        threshold_arr = np.broadcast_to(threshold_arr, shape)\n
    such that\n
    ::\n
        for index in np.ndindex(shape):
            q.interp_x(data_arr[index], data_x_arr, x_arr[index]) \approx threshold_arr[index]
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
    elif isinstance(
        x,
        (
            list,
            tuple,
        ),
    ):
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
        #
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
    return filter_np_results(1 / n * v)

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
        return [
            average(data_list),
        ]
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
    Compute ``(avg, err)`` of ``data_list``.\n
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
        return (
            avg,
            err,
        )
    if n < 2 * block_size:
        block_size = 1
    assert n > block_size
    blocks = block_data(data_list, block_size)
    diff_sqr = average([fsqr(d - avg) for d in blocks])
    fac = abs(eps) * math.sqrt(block_size / (n - block_size))
    err = fac * fsqrt(diff_sqr)
    err = filter_np_results(err)
    return (
        avg,
        err,
    )

def fsqr(data):
    """
    Separately square real and imag part in case of complex types.\n
    :param data: real, complex, np.ndarray like objects.
    :return: squared ``data``.
    :rtype: same type as ``data``.
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
    Separately calculate the square root real and imag part in case of complex types.\n
    :param data: real, complex, np.ndarray like objects.
    :return: squared ``data``.
    :rtype: same type as ``data``.
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
    e.g.: ``q.err_sum(1.4, 2.1, 1.0)`` ==> ``2.7147743920996454``
    """
    err_sqr = sum([fsqr(v) for v in vs])
    err = fsqrt(err_sqr)
    return err

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
    ``val * 10**exp`` is the same as input
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
    ``is_latex`` can be in [ None, False, True, ]
    ``num_float_digit`` or ``num_exp_digit`` can be in [ None, False, True, int, ]
    ``exponent`` can be in [ None, int, ]
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
    assert not ((num_float_digit is False) and (num_exp_digit is False))
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
    ``is_latex`` can be in [ None, False, True, ]
    ``num_float_digit`` or ``num_exp_digit`` can be in [ None, False, True, int, ]
    ``exponent`` can be in [ None, int, ]
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
        if (e_e > e) or (v == 0.0):
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
    assert not ((num_float_digit is False) and (num_exp_digit is False))
    if num_exp_digit is False:
        assert isinstance(num_float_digit, int)
        assert num_float_digit >= 0
        if abs(err) >= 1.0:
            return (f"{{0:.{num_float_digit}f}}({{1:.{num_float_digit}f}})").format(
                val, err
            )
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
            v_str = (f"{{0:.{num_exp_digit}f}}({{1:.{num_exp_digit}f}})").format(v, e_v)
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

def interpolate_list(data_arr, i):
    """
    Old function.
    return approximately data_arr[i]
    Use ``q.interp(data_arr, i, 0)`` instead
    """
    return interp(data_arr, i, 0)

def interpolate(data_arr, i_arr):
    """
    Old function. Use ``q.interp(data_arr, i_arr, -1)`` instead.
    #
    return approximately data_arr[..., i_arr]
    """
    vt = data_arr.transpose()
    if isinstance(i_arr, real_types):
        return interpolate_list(vt, i_arr).transpose()
    else:
        return np.array(
            [interpolate_list(vt, i) for i in i_arr], data_arr.dtype
        ).transpose()
