import qlat_utils as q

import functools
import numpy as np
from .geometry import Geometry

@functools.lru_cache(maxsize=16)
@q.timer
def mk_shift_xg_idx_arr(total_site, xg_shift):
    """
    only work without comm
    `total_site` and `xg_shift` should be tuple of ints.
    """
    total_site = np.array(total_site, dtype=np.int32)
    xg_shift = np.array(xg_shift, dtype=np.int32)
    total_volume = np.prod(total_site)
    xg_idx_arr = np.arange(total_volume).reshape(tuple(reversed(total_site))).T
    geo = Geometry(q.Coordinate(total_site))
    xg_arr = geo.xg_arr()
    xg_arr = (xg_arr + xg_shift) % total_site
    return xg_idx_arr[tuple(xg_arr.T)]

@q.timer
def smear_field(field, coef, n_steps=1):
    """
    only work without comm
    """
    if n_steps == 0:
        return field.copy()
    geo = field.geo()
    assert geo.num_node() == 1
    total_site = geo.total_site()
    xg_shift_list = [
        [ 0, 0, 0, 1],
        [ 0, 0, 1, 0],
        [ 0, 1, 0, 0],
        [ 1, 0, 0, 0],
        [ 0, 0, 0, -1],
        [ 0, 0, -1, 0],
        [ 0, -1, 0, 0],
        [ -1, 0, 0, 0],
    ]
    n_dirs = len(xg_shift_list)
    xg_idx_arr_list = [ mk_shift_xg_idx_arr(tuple(total_site), tuple(xg_shift)) for xg_shift in xg_shift_list ]
    new_field = field.copy()
    coef_dir = coef / n_dirs / (1.0 - coef)
    for k in range(n_steps):
        for xg_idx_arr in xg_idx_arr_list:
            new_field[:] += coef_dir * field[xg_idx_arr]
        new_field *= 1.0 - coef
    return new_field
