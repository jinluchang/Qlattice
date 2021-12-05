import cqlat as c

from qlat.qcd import *

def gf_topology_field(gf : GaugeField):
    geo = gf.geo()
    topf = Field("double", geo, 1)
    c.gf_topology_field(topf, gf)
    return topf

@timer
def gf_topology(gf : GaugeField):
    return gf_topology_field(gf).glb_sum()[0]

def gf_topology_terms_field(gf : GaugeField):
    geo = gf.geo()
    topf = Field("double", geo, 5)
    c.gf_topology_terms_field(topf, gf)
    return topf

@timer
def gf_topology_terms(gf : GaugeField):
    return gf_topology_terms_field(gf).glb_sum()
