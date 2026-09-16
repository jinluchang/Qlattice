"""
Module ``qlat.mpi_utils``
==========================\n
MPI initialisation, teardown, and data-distribution helpers for qlat.
Wraps ``mpi4py`` to set up the Cartesian node grid expected by the C++
runtime and provides ``get_mpi_chunk`` for splitting work across ranks,
together with the pure-Python ``glb_sum`` reduction dispatch.\n
Documentation: ``docs/qlat/qlat_mpi_utils.md``\n
.. note:: Update the documentation when updating this source file.
"""

class q:
    from qlat_utils import (
        timer,
        timer_verbose,
        get_num_node,
        get_id_node,
        Coordinate,
        get_chunk_list,
        displayln_info,
        displayln,
    )
    from .mpi import (
        is_initialized,
        get_size_node,
        get_coor_node,
        begin,
        end,
    )

import numpy as np

from qlat_utils import LatData

from .mpi import (
    glb_sum_complex,
    glb_sum_double,
    glb_sum_long,
)

default_size_node_list = list(
    map(
        q.Coordinate,
        [
            [
                1,
                1,
                1,
                1,
            ],
            [
                1,
                1,
                1,
                2,
            ],
            [
                1,
                1,
                2,
                2,
            ],
            [
                1,
                2,
                2,
                2,
            ],
            [
                1,
                2,
                2,
                2,
            ],
            [
                2,
                2,
                2,
                2,
            ],
            [
                2,
                2,
                2,
                4,
            ],
            [
                2,
                2,
                4,
                4,
            ],
            [
                2,
                4,
                4,
                4,
            ],
            [
                4,
                4,
                4,
                4,
            ],
            [
                4,
                4,
                4,
                8,
            ],
            [
                4,
                4,
                8,
                8,
            ],
            [
                4,
                8,
                8,
                8,
            ],
            [
                8,
                8,
                8,
                8,
            ],
            [
                8,
                8,
                8,
                16,
            ],
            [
                1,
                1,
                1,
                3,
            ],
            [
                1,
                1,
                2,
                3,
            ],
            [
                1,
                2,
                2,
                3,
            ],
            [
                2,
                2,
                2,
                3,
            ],
            [
                2,
                2,
                2,
                6,
            ],
            [
                2,
                2,
                4,
                6,
            ],
            [
                2,
                4,
                4,
                6,
            ],
            [
                4,
                4,
                4,
                6,
            ],
            [
                4,
                4,
                4,
                12,
            ],
            [
                4,
                4,
                8,
                12,
            ],
            [
                4,
                8,
                8,
                12,
            ],
            [
                8,
                8,
                8,
                12,
            ],
        ],
    )
)

comm = None

def set_comm(x):
    global comm
    comm = x

def get_comm():
    return comm

def is_initialized():
    """
    Return whether the global geometry node (`geon`) is initialized, i.e.
    whether `q.begin_with_mpi()` or `q.begin(id_node, size_node)` has been
    called.
    """
    return q.is_initialized()

def begin_with_mpi(size_node_list=None):
    global comm
    from mpi4py import MPI
    #
    comm = MPI.COMM_WORLD
    num_node = comm.size
    id_node = comm.rank
    if size_node_list is None:
        size_node_list = []
    else:
        assert isinstance(size_node_list, list)
        size_node_list = list(map(q.Coordinate, size_node_list))
    size_node_list = size_node_list + default_size_node_list
    size_node = None
    for size_node_check in size_node_list:
        if size_node_check.volume() == num_node:
            size_node = size_node_check
            break
    if size_node is None:
        if id_node == 0:
            q.displayln(size_node_list)
        comm.barrier()
        raise Exception("begin_with_mpi: size_node_list not match num_node")
    q.begin(id_node, size_node)

def end_with_mpi(is_preserving_cache=False):
    q.end(is_preserving_cache)
    from mpi4py import MPI
    #
    MPI.Finalize()

@q.timer_verbose
def show_machine():
    q.displayln(
        f"id_node: {q.get_id_node():4} / {q.get_num_node()}"
        f" ; coor_node: {str(q.get_coor_node()):9}"
        f" / {str(q.get_size_node())}"
    )

def get_mpi_chunk(total_list, *, rng_state=None):
    """
    rng_state has to be the same on all the nodes
    e.g. rng_state = q.RngState("get_mpi_chunk")
    """
    chunk_number = q.get_num_node()
    chunk_id = q.get_id_node()
    chunk_list = q.get_chunk_list(
        total_list, chunk_number=chunk_number, rng_state=rng_state
    )
    if chunk_id < len(chunk_list):
        return chunk_list[chunk_id]
    else:
        return []

### -------------------------------------------------------------------

@q.timer
def glb_sum_np(x):
    """
    x does NOT change
    """
    shape = x.shape
    dtype = x.dtype
    l = list(x.flatten())
    ld = LatData()
    if dtype == np.dtype('float64'):
        ld.from_list(l, is_complex=False)
    elif dtype == np.dtype('int64'):
        ld.from_list(list(map(float, l)), is_complex=False)
    elif dtype == np.dtype('complex128'):
        ld.from_list(l, is_complex=True)
    else:
        q.displayln(dtype)
        assert False
    ld.glb_sum_in_place()
    return np.array(ld.to_list(), dtype=dtype).reshape(shape)

@q.timer
def glb_sum(x):
    """
    x does NOT change
    """
    if isinstance(x, float):
        return glb_sum_double(x)
    elif isinstance(x, complex):
        return glb_sum_complex(x)
    elif isinstance(x, (int, np.int64)):
        return glb_sum_long(x)
    elif isinstance(x, np.ndarray):
        return glb_sum_np(x)
    elif isinstance(x, list):
        return [ glb_sum(x_i) for x_i in x ]
    elif isinstance(x, tuple):
        return tuple(glb_sum(x_i) for x_i in x)
    else:
        # possible types: q.LatData
        return x.glb_sum()

