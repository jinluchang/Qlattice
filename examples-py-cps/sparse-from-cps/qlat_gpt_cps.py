# import numpy as np
import qlat as q
import qlat_cps as qc
import qlat_gpt as qg
import gpt as g

def begin_with_gpt_cps(total_site):
    """
    Note: It will likely be consistent only when MPI is only along one axis. e.g. with "-qmp-geom 1 1 1 4 --mpi 1.1.1.4"
    """
    qg.begin_with_gpt()
    show_geo_total_site(total_site)
    qc.begin_with_cps(total_site)
    check_cps_gpt_id_node()

def end_with_gpt_cps():
    q.end() # end begin_with_cps but do not finalize MPI
    qg.end_with_gpt()

@q.timer_verbose
def check_cps_gpt_id_node():
    """
    Check if the `id_node` is consistent between CPS and GPT.
    Call this function after `qc.begin_with_cps(total_site)`.
    #
    Note: It will likely be consistent only when MPI is only along one axis. e.g. with "-qmp-geom 1 1 1 4 --mpi 1.1.1.4"
    """
    q.sync_node()
    error_count = 0
    id_node = q.get_id_node()
    size_node = q.get_size_node()
    coor_node = q.get_coor_node()
    grid = qg.mk_grid()
    size_node_gpt = q.Coordinate(grid.mpi)
    coor_node_gpt = q.Coordinate(grid.processor_coor)
    id_node_gpt = q.index_from_coordinate(coor_node, size_node)
    if id_node != id_node_gpt:
        q.displayln(f"ERROR: id_node={id_node} ; id_node_gpt={id_node_gpt}")
        error_count += 1
    error_count = q.glb_sum(error_count)
    if error_count > 0:
        raise Exception(f"id_node != id_node_gpt ; error_count={error_count}")
    if size_node != size_node_gpt:
        q.displayln(f"ERROR: id_node={id_node} ; size_node={size_node} ; size_node_gpt={size_node_gpt}")
        error_count += 1
    error_count = q.glb_sum(error_count)
    if error_count > 0:
        raise Exception(f"size_node != size_node_gpt ; error_count={error_count}")
    if coor_node != coor_node_gpt:
        q.displayln(f"ERROR: id_node={id_node} ; coor_node={coor_node} ; coor_node_gpt={coor_node_gpt}")
        error_count += 1
    error_count = q.glb_sum(error_count)
    if error_count > 0:
        raise Exception(f"coor_node != coor_node_gpt ; error_count={error_count}")
    id_node_gpt_rank = g.rank()
    if id_node != id_node_gpt_rank:
        q.displayln(f"WARNING: id_node={id_node} ; id_node_gpt_rank={id_node_gpt_rank}")
        error_count += 1
    error_count = q.glb_sum(error_count)
    if error_count > 0:
        q.displayln_info(f"WARNING: id_node != id_node_gpt_rank ; error_count={error_count}")
    q.sync_node()

@q.timer_verbose
def show_geo_total_site(total_site):
    geo = q.Geometry(total_site)
    q.displayln_info("-" * 80)
    q.displayln_info(f"total_site = {total_site}")
    q.displayln_info(f"geo = {geo}")
    q.displayln_info("-" * 80)