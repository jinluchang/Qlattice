import cqlat as c

from qlat.field import *
from qlat.lat_io import *
from qlat.rng_state import *

from cqlat import index_from_coordinate, coordinate_from_index

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
    return x.qnorm()

def show(x):
    return x.show()

def unitarize(x):
    x.unitarize()

def random_permute(l, rs):
    # Do not change ``l''.
    # Return a new permutated list.
    assert isinstance(l, list)
    assert isinstance(rs, RngState)
    return c.random_permute(l, rs)
