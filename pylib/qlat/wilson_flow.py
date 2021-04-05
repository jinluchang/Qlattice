import cqlat as c

from qlat.qcd import *

def gf_wilson_flow_step(gf, epsilon, c1 = 0.0):
    assert isinstance(gf, GaugeField)
    assert isinstance(epsilon, float)
    assert isinstance(c1, float)
    c.gf_wilson_flow_step(gf, epsilon, c1)

def gf_energy_density(gf):
    assert isinstance(gf, GaugeField)
    c.gf_energy_density(gf)
