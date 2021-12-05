import cqlat as c

from qlat.qcd import *

def gf_energy_density(gf : GaugeField):
    return c.gf_energy_density(gf)

def gf_wilson_flow_step(gf : GaugeField, epsilon : float, c1 : float = 0.0):
    return c.gf_wilson_flow_step(gf, epsilon, c1)

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
