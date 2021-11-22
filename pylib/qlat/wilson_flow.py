import cqlat as c

from qlat.qcd import *

def gf_wilson_flow_step(gf : GaugeField, epsilon : float, c1 : float = 0.0):
    return c.gf_wilson_flow_step(gf, epsilon, c1)

def gf_energy_density(gf : GaugeField):
    return c.gf_energy_density(gf)

def gf_ape_smear(gf : GaugeField, alpha : float, steps : int = 1):
    return c.gf_ape_smear(gf, gf, alpha, steps)

def gf_hyp_smear(gf : GaugeField, alpha1 : float, alpha2 : float, alpha3 : float):
    # values in paper is 0.75 0.6 0.3
    # 10.1103/PhysRevD.64.034504 Eq(4)
    return c.gf_hyp_smear(gf, gf, alpha1, alpha2, alpha3)

@timer
def gf_wilson_flow(gf: GaugeField, flow_time : float, steps : int,
        *, c1 : float = 0.0, existing_flow_time : float = 0.0):
    epsilon = flow_time / steps
    energy_density_list = []
    for i in range(steps):
        gf_wilson_flow_step(gf, epsilon, c1)
        t = (i + 1) * epsilon + existing_flow_time
        energy_density = gf_energy_density(gf)
        energy_density_list.append(energy_density)
        displayln_info(f"gf_wilson_flow: t={t} ; E={energy_density} ; t^2 E={t*t*energy_density}")
    return energy_density_list

def gf_topology_field(gf : GaugeField):
    geo = gf.geo()
    topf = Field("double", geo, 1)
    c.gf_topology_field(topf, gf)
    return topf

@timer
def gf_topology(gf : GaugeField):
    return gf_topology_field(gf).glb_sum()[0]
