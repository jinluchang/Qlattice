"""
Module ``qlat.fields_io_utils``
================================\n
Pure-Python helpers for shuffled field I/O that do not call C++ functions
directly: opening a field container by mode and checking its contents.\n
"""

import qlat_utils as q

from qlat_utils import Coordinate

from .fields_io import (
        ShuffledFieldsReader,
        ShuffledFieldsWriter,
        properly_truncate_fields,
        )

@q.timer
def open_fields(path, mode, new_size_node=None):
    """
    path can be the folder path or the 'geon-info.txt' path
    """
    if path[-14:] == "/geon-info.txt":
        path = path[:-14]
    if mode == "r":
        return ShuffledFieldsReader(path, new_size_node)
    elif mode == "w":
        assert new_size_node is not None
        return ShuffledFieldsWriter(path, new_size_node)
    elif mode == "a":
        if new_size_node is None:
            new_size_node = Coordinate()
        return ShuffledFieldsWriter(path, new_size_node, True)
    else:
        raise Exception("open_fields")

@q.timer
def check_fields(path, is_check_all=True, new_size_node=None):
    """
    return list of field that is stored successful
    """
    if path[-14:] == "/geon-info.txt":
        path = path[:-14]
    is_only_check = True
    return properly_truncate_fields(path, is_check_all, is_only_check, new_size_node)
