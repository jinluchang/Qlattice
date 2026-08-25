"""
Module ``qlat.selected_points_io``
===================================\n
Utilities for loading and saving ``SelectedPoints`` objects across MPI ranks.\n
Documentation: ``docs/qlat/qlat_selected_points_io.md``
.. note:: Update the documentation when updating this source file.
"""

from qlat_utils import timer, get_id_node, get_num_node, displayln_info, get_fname
from .field_selection import PointsSelection, SelectedShufflePlan


@timer
def load_selected_points_list(cls, psel_list, path_list):
    """
    Load SelectedPoints from disk and shuffle to Local distribution.

    Distributed read: each file is assigned to one MPI rank (round-robin).
    That rank reads the file from disk. Non-loading ranks create empty sp
    objects with 0 points. The ``"l_from_g"`` shuffle plan then redistributes
    the data so every node holds only its locally-owned points.

    - ``cls``: a SelectedPoints type, e.g. ``q.SelectedPointsRealD``.
    - ``psel_list``: list of PointsSelection (Global distribution), one per
      path, available on every node.
    - ``path_list``: list of file paths, each saved by ``sp.save(path)``.

    Returns a list of locally-distributed SelectedPoints objects.
    """
    fname = get_fname()
    id_node = get_id_node()
    num_node = get_num_node()
    n_paths = len(path_list)
    assert len(psel_list) == n_paths
    root_list = [i % num_node for i in range(n_paths)]
    ssp = SelectedShufflePlan("l_from_g", psel_list, root_list)
    sp_list = []
    for i, path in enumerate(path_list):
        if i % num_node == id_node:
            sp = cls(psel_list[i])
            sp.load(path, is_sync_node=False)
        else:
            sp = cls(None)
        sp_list.append(sp)
    sp_dst_list = ssp.shuffle_sp_list(cls, sp_list)
    return sp_dst_list
