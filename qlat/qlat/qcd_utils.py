"""
Module ``qlat.qcd_utils``
==========================\n
Gauge field helpers that do not call C++ functions directly: gauge field
summaries, Wilson-line path loops, and left-expanded fields.\n
"""

class q:
    from qlat_utils import (
        timer_verbose,
        displayln_info,
    )
    from .field_types import (
        FieldColorMatrix,
    )
    from .geometry import (
        geo_resize,
    )
    from .field_utils_utils import (
        field_expanded,
        refresh_expanded_1,
    )
    from .qcd import (
        gf_wilson_line_no_comm,
    )

@q.timer_verbose
def gf_show_info(gf):
    assert gf is not None
    q.displayln_info(
        f"gf_show_info: plaq = {gf.plaq():.16F} ; link_trace = {gf.link_trace():.16F}."
    )

def gf_wilson_lines_no_comm(gf_ext, path_list):
    """
    path_list = [ path_spec, ... ]
    e.g. path_spec = [ mu, mu, nu, -mu-1, -mu-1, ]
    e.g. path_spec = ([ mu, nu, -mu-1, ], [ 2, 1, 2, ],)
    return wlf
    """
    multiplicity = len(path_list)
    geo = q.geo_resize(gf_ext.geo)
    wlf = q.FieldColorMatrix(geo, multiplicity)
    for m, p in enumerate(path_list):
        if isinstance(p, tuple) and len(p) == 2:
            path, path_n = p
            q.gf_wilson_line_no_comm(wlf, m, gf_ext, path, path_n)
        else:
            path = p
            q.gf_wilson_line_no_comm(wlf, m, gf_ext, path)
    return wlf

def mk_left_expanded_field(gf):
    """
    Return left expanded field.
    Similar to ``set_left_expanded_gauge_field`` in C++
    """
    gf1 = q.field_expanded(gf, 1, 0)
    q.refresh_expanded_1(gf1)
    return gf1
