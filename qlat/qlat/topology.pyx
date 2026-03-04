# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from qlat_utils.all cimport *
from . cimport everything as cc
from .field_types cimport FieldRealD
from .geometry cimport Geometry
from .qcd cimport GaugeField
from .hmc cimport GaugeMomentum
from .gauge_action cimport GaugeAction

from .hmc import set_gm_force, gf_evolve
from .wilson_flow import gf_wilson_flow_step, gf_energy_density_field, gf_energy_derivative_density_field
import qlat_utils as q

@q.timer
def gf_plaq_action_density_field(GaugeField gf):
    r"""
    return paf
    paf.geo.multiplicity == 1
    \sum_P (1 - 1/3 * Re Tr U_P)
    #
    6 plaqs in the above sum (for each site).
    #
    Action = beta * total_volume * 6 * (1 - avg_plaq)
    Action = beta * total_volume * action_density
    Single instanton action = 8 * sqr(PI) / g^2
    beta = 6/g^2
    """
    paf = FieldRealD()
    cc.clf_plaq_action_density_field(paf.xx, gf.xxx().val())
    return paf

@q.timer
def gf_spatial_plaq_action_density_field(GaugeField gf):
    r"""
    return paf
    paf.geo.multiplicity == 1
    \sum_P(spatial only) (1 - 1/3 * Re Tr U_P)
    """
    paf = FieldRealD()
    cc.clf_spatial_plaq_action_density_field(paf.xx, gf.xxx().val())
    return paf

@q.timer
def gf_plaq_action_density(GaugeField gf):
    r"""
    return pa
    ininstance(pa, float)
    pa = gf_plaq_action_density_field(gf).glb_sum()[:].item() / total_volume
    """
    cdef Geometry geo = gf.geo
    cdef cc.Long total_volume = geo.total_volume
    return gf_plaq_action_density_field(gf).glb_sum()[:].item() / total_volume

@q.timer
def gf_spatial_plaq_action_density(GaugeField gf):
    r"""
    return pa
    ininstance(pa, float)
    pa = gf_spatial_plaq_action_density_field(gf).glb_sum()[:].item() / total_volume
    """
    cdef Geometry geo = gf.geo
    cdef cc.Long total_volume = geo.total_volume
    return gf_spatial_plaq_action_density_field(gf).glb_sum()[:].item() / total_volume

@q.timer
def gf_topology_field_clf(GaugeField gf):
    r"""
    return topf
    topf.geo.multiplicity == 1
    Use the basic gf_clover_leaf_field
    NOT using 5 loop improved definition
    """
    topf = FieldRealD()
    cc.clf_topology_field(topf.xx, gf.xxx().val())
    return topf

@q.timer
def gf_topology_clf(GaugeField gf):
    r"""
    return top
    ininstance(top, float)
    Use the basic gf_clover_leaf_field
    NOT using 5 loop improved definition
    """
    return gf_topology_field_clf(gf).glb_sum()[:].item()

@q.timer
def gf_topology_field(GaugeField gf):
    r"""
    return topf
    topf.geo.multiplicity == 1
    Using the 5 loop improved definition
    https://arxiv.org/pdf/hep-lat/9701012v2.pdf
    """
    topf = FieldRealD()
    cc.clf_topology_field_5(topf.xx, gf.xxx().val())
    return topf

@q.timer
def gf_topology(GaugeField gf):
    r"""
    return top
    ininstance(top, float)
    Using the 5 loop improved definition Eq. (2-7)
    https://arxiv.org/pdf/hep-lat/9701012v2.pdf
    """
    return gf_topology_field(gf).glb_sum()[:].item()

@q.timer
def gf_topology_terms_field(GaugeField gf):
    r"""
    return topf;
    topf.geo.multiplicity == 5
    sum of the 5 terms should equal to gf_topology_field
    """
    topf = FieldRealD()
    cc.clf_topology_field_5_terms(topf.xx, gf.xxx().val())
    return topf

@q.timer
def gf_topology_terms(GaugeField gf):
    r"""
    return top_terms;
    top_terms.shape == (5,)
    top_terms.dtype == np.float64
    sum of the 5 terms should equal to gf_topology
    """
    return gf_topology_terms_field(gf).glb_sum()[0, :]
