from qlat.lat_io import *

import math
import copy
import numpy as np

def check_zero(x):
    if isinstance(x, (int, float, complex)) and 0 == x:
        return True
    return False

class Data:

    def __init__(self, val):
        # supported value types:
        # numeric
        # numpy.array
        # q.LatData
        if isinstance(val, Data):
            self.val = val.val
            assert not isinstance(self.val)
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
            else:
                return Data(self.val + other.val)
        else:
            return self + Data(other)

    __radd__ = __add__

    def __mul__(self, other):
        if isinstance(other, Data):
            if check_zero(self.val) or check_zero(other.val):
                return Data(0)
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
        return Data(-self.val)

    def __pos__(self):
        return self

    def __sub__(self, other):
        if isinstance(other, Data):
            if check_zero(self.val):
                return Data(-other.val)
            elif check_zero(other.val):
                return self
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

def average(data_list):
    n = len(data_list)
    v = sum(data_list)
    return 1/n * v

def jackknife(data_list, eps = 1):
    # normal jackknife uses eps = 1 / sqrt(len(data_list))
    n = len(data_list)
    fac = eps / math.sqrt(n)
    avg = average(data_list)
    jks = [ avg, ]
    for data in data_list:
        jks.append(avg + fac * (data - avg))
    return jks

def fsqr(data):
    if isinstance(data, (float, int,)):
        return data * data
    elif isinstance(data, complex):
        r = data.real
        i = data.imag
        return complex(r * r, i * i)
    elif isinstance(data, np.ndarray):
        if data.dtype in [ np.float, np.int, ]:
            return np.square(data)
        elif data.dtype in [ np.complex, ]:
            return np.square(data.real) + 1j * np.square(data.imag)
        else:
            raise Exception(f"fsqr data={data} type not supported")
    elif isinstance(data, Data):
        return Data(fsqr(data.val))
    else:
        raise Exception(f"fsqr data={data} type not supported")

def fsqrt(data):
    if isinstance(data, (float, int,)):
        return math.sqrt(data)
    elif isinstance(data, complex):
        r = data.real
        i = data.imag
        return complex(math.sqrt(r), math.sqrt(i))
    elif isinstance(data, np.ndarray):
        if data.dtype in [ np.float, np.int, ]:
            return np.sqrt(data)
        elif data.dtype in [ np.complex, ]:
            return np.sqrt(data.real) + 1j * np.sqrt(data.imag)
        else:
            raise Exception(f"fsqr data={data} type not supported")
    elif isinstance(data, Data):
        return Data(fsqrt(data.val))
    else:
        raise Exception(f"fsqr {data} type not supported")

def jk_avg(jk_list):
    return jk_list[0]

def jk_err(jk_list, eps = 1):
    # same eps as the eps used in the 'jackknife' function
    avg = jk_avg(jk_list)
    diff_sqr = average([ fsqr(jk - avg) for jk in jk_list[1:] ])
    return (1/eps) * fsqrt(diff_sqr)
