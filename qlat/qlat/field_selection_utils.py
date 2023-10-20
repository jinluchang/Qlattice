from qlat_utils import *

from .c import *

from . import c

cache_point_selection = mk_cache("point_selection")

def get_psel_tslice(total_site, *, t_dir = 3):
    """
    if t_dir = 3, then [ [0,0,0,0,], [0,0,0,1,], ..., [0,0,0,total_site[3]-1],]
    if t_dir = 2, then [ [0,0,0,0,], [0,0,1,0,], ..., [0,0,total_site[2]-1],0,]
    need total_site to set the psel.geo property
    """
    assert 0 <= t_dir and t_dir < 4
    total_site = Coordinate(total_site)
    param_tuple = (total_site[0], total_site[1], total_site[2], total_site[3], t_dir,)
    if param_tuple not in cache_point_selection:
        psel = PointsSelection(None, Geometry(total_site))
        c.set_tslice_psel(psel, total_site[t_dir], t_dir)
        cache_point_selection[param_tuple] = psel
    return cache_point_selection[param_tuple]

def is_matching_fsel(fsel1, fsel2):
    return c.is_matching_fsel(fsel1, fsel2)
