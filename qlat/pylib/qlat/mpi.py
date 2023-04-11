from qlat_utils import *
import qlat.c as c
from qlat.c import get_size_node, get_coor_node
from qlat.c import sync_node
from qlat.c import begin, end
from math import prod

import numpy as np

default_size_node_list = list(map(Coordinate, [
    [ 1, 1, 1, 1, ],
    [ 1, 1, 1, 2, ],
    [ 1, 1, 2, 2, ],
    [ 1, 2, 2, 2, ],
    [ 1, 2, 2, 2, ],
    [ 2, 2, 2, 2, ],
    [ 2, 2, 2, 4, ],
    [ 2, 2, 4, 4, ],
    [ 2, 4, 4, 4, ],
    [ 4, 4, 4, 4, ],
    [ 4, 4, 4, 8, ],
    [ 4, 4, 8, 8, ],
    [ 4, 8, 8, 8, ],
    [ 8, 8, 8, 8, ],
    [ 8, 8, 8, 16, ],
    [ 1, 1, 1, 3, ],
    [ 1, 1, 2, 3, ],
    [ 1, 2, 2, 3, ],
    [ 2, 2, 2, 3, ],
    [ 2, 2, 2, 6, ],
    [ 2, 2, 4, 6, ],
    [ 2, 4, 4, 6, ],
    [ 4, 4, 4, 6, ],
    [ 4, 4, 4, 12, ],
    [ 4, 4, 8, 12, ],
    [ 4, 8, 8, 12, ],
    [ 8, 8, 8, 12, ],
    ]))

def begin_with_mpi(size_node_list = None):
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    num_node = comm.size
    id_node = comm.rank
    if size_node_list is None:
        size_node_list = []
    else:
        assert isinstance(size_node_list, list)
        size_node_list = list(map(Coordinate, size_node_list))
    size_node_list = size_node_list + default_size_node_list
    size_node = None
    for size_node_check in size_node_list:
        if prod(size_node_check.list()) == num_node:
            size_node = size_node_check
            break
    if size_node is None:
        if id_node == 0:
            displayln(size_node_list)
        comm.barrier()
        raise Exception("begin_with_mpi: size_node_list not match num_node")
    c.begin(id_node, size_node)

def end_with_mpi(is_preserving_cache = False):
    c.end(is_preserving_cache)
    from mpi4py import MPI
    MPI.Finalize()

@timer
def glb_sum_np(x):
    """
    x does NOT change
    """
    from qlat_utils.lat_io import LatData
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
    """
    x does NOT change
    """
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
        # possible types: q.LatData
        return x.glb_sum()

@timer_verbose
def show_machine():
    displayln(f"id_node: {get_id_node():4} / {get_num_node()}"
            f" ; coor_node: {str(get_coor_node()):9}"
            f" / {str(get_size_node())}")
