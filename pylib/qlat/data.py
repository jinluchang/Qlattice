from qlat.lat_io import *

import numpy as np
import math
import copy

class Data:

    def __init__(self, val):
        if isinstance(val, Data):
            self.val = val.val
            assert not isinstance(self.val)
        else:
            self.val = val

    def get_val(self):
        return self.val

    def __copy__(self):
        return Data(copy.copy(self.val))

    def __deepcopy__(self, memo):
        return Data(copy.deepcopy(self.val, memo))

    def __add__(self, other):
        if isinstance(other, Data):
            if self.val == 0:
                return other
            elif other.val == 0:
                return self
            else:
                return Data(self.val + other.val)
        else:
            if self.val == 0:
                return Data(other)
            elif other == 0:
                return self
            else:
                return Data(self.val + other)

    __radd__ = __add__

    def __mul__(self, other):
        if isinstance(other, Data):
            return Data(self.val * other.val)
        else:
            return Data(self.val * other)

    def __rmul__(self, other):
        if isinstance(other, Data):
            return Data(other.val * self.val)
        else:
            return Data(other * self.val)

    def __neg__(self):
        return Data(-self.val)

    def __pos__(self):
        return self

    def __sub__(self, other):
        if isinstance(other, Data):
            if self.val == 0:
                return other
            elif other.val == 0:
                return self
            else:
                return Data(self.val - other.val)
        else:
            if self.val == 0:
                return Data(other)
            elif other == 0:
                return self
            else:
                return Data(self.val - other)

    def __rsub__(self, other):
        if isinstance(other, Data):
            if self.val == 0:
                return -other
            elif other.val == 0:
                return self
            else:
                return Data(other.val - self.val)
        else:
            if self.val == 0:
                return Data(-self.val)
            elif other == 0:
                return self
            else:
                return Data(other - self.val)

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
    if isinstance(data, float) or isinstance(data, int):
        return data * data
    elif isinstance(data, complex):
        r = data.real
        i = data.imag
        return complex(r * r, i * i)
    else:
        raise Exception(f"fsqr {data} type not supported")

def fsqrt(data):
    if isinstance(data, float) or isinstance(data, int):
        return math.sqrt(data)
    elif isinstance(data, complex):
        r = data.real
        i = data.imag
        return complex(math.sqrt(r), math.sqrt(i))
    else:
        raise Exception(f"fsqr {data} type not supported")

def jk_avg(jk_list):
    return jk_list[0]

def jk_err(jk_list, eps = 1):
    # same eps as the eps used in the 'jackknife' function
    avg = jk_avg(jk_list)
    diff_sqr = average([ fsqr(jk - avg) for jk in jk_list[1:] ])
    return (1/eps) * fsqrt(diff_sqr)
