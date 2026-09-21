"""
Module ``qlat.topology_utils``
===============================\n
Topological charge helpers that do not call C++ functions directly: global
sums of the clover-leaf and 5-loop improved topology fields.\n
"""

class q:
    from qlat_utils import (
        timer,
    )
    from .topology import (
        gf_topology_field_clf,
        gf_topology_field,
        gf_topology_terms_field,
    )

@q.timer
def gf_topology_clf(gf):
    r"""
    return top
    ininstance(top, float)
    Use the basic gf_clover_leaf_field
    NOT using 5 loop improved definition
    """
    return q.gf_topology_field_clf(gf).glb_sum()[:].item()

@q.timer
def gf_topology(gf):
    r"""
    return top
    ininstance(top, float)
    Using the 5 loop improved definition Eq. (2-7)
    https://arxiv.org/pdf/hep-lat/9701012v2.pdf
    """
    return q.gf_topology_field(gf).glb_sum()[:].item()

@q.timer
def gf_topology_terms(gf):
    r"""
    return top_terms;
    top_terms.shape == (5,)
    top_terms.dtype == np.float64
    sum of the 5 terms should equal to gf_topology
    """
    return q.gf_topology_terms_field(gf).glb_sum()[0, :]
