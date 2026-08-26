"""
Module ``qlat.selected_points_io``
===================================\n
Utilities for loading and saving ``SelectedPoints`` objects across MPI ranks.\n
Documentation: ``docs/qlat/qlat_selected_points_io.md``
.. note:: Update the documentation when updating this source file.
"""

from qlat_utils import timer, get_id_node, get_num_node, get_fname
from qlat_utils import sync_node
from .field_selection import PointsSelection, SelectedShufflePlan

@timer
def load_selected_points_list(cls, psel_list, path_list):
    """
    Load SelectedPoints from disk and shuffle to Local distribution.\n
    Distributed read: each file is assigned to one MPI rank (round-robin).
    That rank loads both the psel and sp from disk. Non-loading ranks
    contribute no data. The ``"l_from_g"`` shuffle plan then redistributes
    the data so every node holds only its locally-owned points.\n
    - ``cls``: a SelectedPoints type, e.g. ``q.SelectedPointsRealD``.
    - ``psel_list``: list of PointsSelection (Global) or file paths (str).
      If a string, the psel is loaded from the path on the owning node only.
      If a PointsSelection, it is used directly (must be Global type).
    - ``path_list``: list of sp file paths, each saved by ``sp.save(path)``.\n
    Returns a list of locally-distributed SelectedPoints objects.
    """
    id_node = get_id_node()
    num_node = get_num_node()
    n_paths = len(path_list)
    assert len(psel_list) == n_paths
    sync_node()
    root_list = [i % num_node for i in range(n_paths)]
    psel_empty = PointsSelection()
    # Load psel and sp from file only on the owning node
    psel_list_resolved = []
    sp_list = []
    for i, path in enumerate(path_list):
        if root_list[i] == id_node:
            if isinstance(psel_list[i], str):
                psel = PointsSelection()
                psel.load(psel_list[i], is_sync_node=False)
            else:
                psel = psel_list[i]
            sp = cls(psel)
            sp.load(path, is_sync_node=False)
        else:
            psel = psel_empty
            sp = cls()
        sp_list.append(sp)
        psel_list_resolved.append(psel)
    ssp = SelectedShufflePlan("l_from_g", psel_list_resolved, root_list)
    sp_dst_list = ssp.shuffle_sp_list(cls, sp_list)
    return sp_dst_list

@timer
def convert_selected_points_dist_type(
    sp, points_dist_type, *, sp_points_dist_type=None, root=None
):
    fname = get_fname()
    cls = type(sp)
    if sp_points_dist_type is not None:
        assert sp.points_dist_type == sp_points_dist_type
    if sp.points_dist_type == "l":
        if points_dist_type == "g":
            psel = sp.psel
            if root is None:
                root = 0
                is_bcast = True
            else:
                is_bcast = False
            geo = psel.geo
            ssp = SelectedShufflePlan("g_from_l", psel, root, geo)
            sp = ssp.shuffle_sp(cls, sp)
            if is_bcast:
                sp.bcast(root=root)
            return sp
        else:
            raise Exception(f"{fname}: {sp.points_dist_type=} {points_dist_type=}")
    elif sp.points_dist_type == "g":
        if points_dist_type == "l":
            psel = sp.psel
            if root is None:
                root = 0
            ssp = SelectedShufflePlan("l_from_g", psel, root)
            sp = ssp.shuffle_sp(cls, sp)
            return sp
        else:
            raise Exception(f"{fname}: {sp.points_dist_type=} {points_dist_type=}")
    else:
        raise Exception(f"{fname}: {sp.points_dist_type=}")
