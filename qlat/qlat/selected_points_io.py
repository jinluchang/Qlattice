"""
Module ``qlat.selected_points_io``
===================================\n
Utilities for loading and saving ``SelectedPoints`` objects across MPI ranks.\n
Documentation: ``docs/qlat/qlat_selected_points_io.md``
.. note:: Update the documentation when updating this source file.
"""

from qlat_utils import timer, get_id_node, get_num_node, displayln_info, get_fname
from qlat_utils import sync_node
from .field_selection import PointsSelection, SelectedShufflePlan


@timer
def load_selected_points_list(cls, psel_list, path_list):
    """
    Load SelectedPoints from disk and shuffle to Local distribution.

    Distributed read: each file is assigned to one MPI rank (round-robin).
    That rank loads both the psel and sp from disk. Non-loading ranks
    contribute no data. The ``"l_from_g"`` shuffle plan then redistributes
    the data so every node holds only its locally-owned points.

    - ``cls``: a SelectedPoints type, e.g. ``q.SelectedPointsRealD``.
    - ``psel_list``: list of PointsSelection (Global) or file paths (str).
      If a string, the psel is loaded from the path on the owning node only.
      If a PointsSelection, it is used directly (must be Global type).
    - ``path_list``: list of sp file paths, each saved by ``sp.save(path)``.

    Returns a list of locally-distributed SelectedPoints objects.
    """
    fname = get_fname()
    id_node = get_id_node()
    num_node = get_num_node()
    n_paths = len(path_list)
    assert len(psel_list) == n_paths
    sync_node()
    root_list = [i % num_node for i in range(n_paths)]
    # Load psel from file only on the owning node (where sp will also be loaded)
    psel_resolved = []
    psel_empty = PointsSelection()
    for i in range(n_paths):
        if root_list[i] == id_node:
            if isinstance(psel_list[i], str):
                psel_g = PointsSelection()
                psel_g.load(psel_list[i], is_sync_node=False)
                psel_resolved.append(psel_g)
            else:
                psel_resolved.append(psel_list[i])
        else:
            psel_resolved.append(psel_empty)
    ssp = SelectedShufflePlan("l_from_g", psel_resolved, root_list)
    sp_list = []
    for i, path in enumerate(path_list):
        if root_list[i] == id_node:
            sp = cls(psel_resolved[i])
            sp.load(path, is_sync_node=False)
        else:
            sp = cls()
        sp_list.append(sp)
    sp_dst_list = ssp.shuffle_sp_list(cls, sp_list)
    return sp_dst_list
