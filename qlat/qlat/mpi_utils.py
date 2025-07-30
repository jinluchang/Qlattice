import numpy as np

class q:
    from qlat_utils import (
            timer_verbose,
            get_num_node,
            get_id_node,
            Coordinate,
            get_chunk_list,
            displayln_info,
            displayln,
            )
    from .c import (
            get_size_node,
            get_coor_node,
            begin,
            end,
            )

default_size_node_list = list(map(q.Coordinate, [
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

comm = None

def set_comm(x):
    global comm
    comm = x

def get_comm():
    return comm

def begin_with_mpi(size_node_list=None):
    global comm
    from mpi4py import MPI
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
    MPI.Finalize()

@q.timer_verbose
def show_machine():
    q.displayln(f"id_node: {q.get_id_node():4} / {q.get_num_node()}"
            f" ; coor_node: {str(q.get_coor_node()):9}"
            f" / {str(q.get_size_node())}")

def get_mpi_chunk(total_list, *, rng_state=None):
    """
    rng_state has to be the same on all the nodes
    e.g. rng_state = q.RngState("get_mpi_chunk")
    """
    chunk_number = q.get_num_node()
    chunk_id = q.get_id_node()
    chunk_list = q.get_chunk_list(total_list, chunk_number=chunk_number, rng_state=rng_state)
    if chunk_id < len(chunk_list):
        return chunk_list[chunk_id]
    else:
        return []
