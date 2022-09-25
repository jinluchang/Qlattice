import qlat.cqlat as c

from qlat.cqlat import get_size_node, get_coor_node
from qlat.cqlat import sync_node

from qlat_utils import *

import numpy as np

def begin(*args):
    c.begin(*args)

def end():
    clean_cache()
    c.end()

@timer
def glb_sum_np(x):
    from qlat_utils.lat_io import LatData
    # x does NOT change
    shape = x.shape
    dtype = x.dtype
    l = list(x.flatten())
    ld = LatData()
    if dtype == np.dtype('float64'):
        ld.from_list(l, is_complex = False)
    elif dtype == np.dtype('int64'):
        ld.from_list(list(map(float, l)), is_complex = False)
    elif dtype == np.dtype('complex128'):
        ld.from_list(l, is_complex = True)
    else:
        displayln(dtype)
        assert False
    ld.glb_sum_in_place()
    return np.array(ld.to_list(), dtype = dtype).reshape(shape)

@timer
def glb_sum(x):
    # x does NOT change
    if isinstance(x, float):
        return c.glb_sum_double(x)
    elif isinstance(x, complex):
        return c.glb_sum_complex(x)
    elif isinstance(x, int):
        return c.glb_sum_long(x)
    elif isinstance(x, np.ndarray):
        return glb_sum_np(x)
    elif isinstance(x, list):
        return [ glb_sum(x_i) for x_i in x ]
    elif isinstance(x, tuple):
        return tuple([ glb_sum(x_i) for x_i in x ])
    else:
        # should not change x
        return x.glb_sum()

@timer_verbose
def show_machine():
    displayln(f"id_node: {get_id_node():4} / {get_num_node()}"
            f" ; coor_node: {str(get_coor_node()):9}"
            f" / {str(get_size_node())}")
