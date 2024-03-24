from qlat_utils import *
from . import c
from .c import get_size_node, get_coor_node
from .c import sync_node
from .c import begin, end

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

def begin_with_mpi(size_node_list=None):
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
        if size_node_check.volume() == num_node:
            size_node = size_node_check
            break
    if size_node is None:
        if id_node == 0:
            displayln(size_node_list)
        comm.barrier()
        raise Exception("begin_with_mpi: size_node_list not match num_node")
    c.begin(id_node, size_node)

def end_with_mpi(is_preserving_cache=False):
    c.end(is_preserving_cache)
    from mpi4py import MPI
    MPI.Finalize()

@timer_verbose
def show_machine():
    displayln(f"id_node: {get_id_node():4} / {get_num_node()}"
            f" ; coor_node: {str(get_coor_node()):9}"
            f" / {str(get_size_node())}")

def get_mpi_chunk(total_list, *, rng_state=None):
    """
    rng_state has to be the same on all the nodes
    e.g. rng_state = q.RngState("get_mpi_chunk")
    """
    chunk_number = get_num_node()
    chunk_id = get_id_node()
    chunk_list = get_chunk_list(total_list, chunk_number=chunk_number, rng_state=rng_state)
    if chunk_id < len(chunk_list):
        return chunk_list[chunk_id]
    else:
        return []

def glb_sum_list(ret):
    displayln_info("glb_sum_list: deprecated")
    # deprecated (use glb_sum instead)
    # ret = [ va, vb, ... ]
    # return [ glb_sum(va), glb_sum(vb), ... ]
    return glb_sum(ret)
