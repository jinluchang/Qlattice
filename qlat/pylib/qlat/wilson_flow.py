import qlat.cqlat as c

from qlat.qcd import *
from qlat.hmc import *

def gf_energy_density(gf : GaugeField):
    return c.gf_energy_density(gf)

def gf_wilson_flow_force(gf, c1 = 0.0):
    ga = GaugeAction(3.0, c1)
    gm_force = GaugeMomentum()
    set_gm_force(gm_force, gf, ga)
    return gm_force

def gf_wilson_flow_step(gf, epsilon, *, c1 = 0.0, wilson_flow_integrator_type = "runge-kutta"):
    if wilson_flow_integrator_type == "runge-kutta":
        return gf_wilson_flow_runge_kutta(gf, epsilon, c1 = c1)
    elif wilson_flow_integrator_type == "euler":
        return gf_wilson_flow_euler(gf, epsilon, c1 = c1)

def gf_wilson_flow_runge_kutta(gf, epsilon, *, c1 = 0.0):
    # Runge-Kutta scheme
    # http://arxiv.org/abs/1006.4518v3
    return c.gf_wilson_flow_step(gf, epsilon, c1)

def gf_wilson_flow_euler(gf, epsilon, *, c1 = 0.0):
    force = gf_wilson_flow_force(gf, c1 = c1)
    return gf_evolve(gf, force, epsilon)

@timer
def gf_wilson_flow(gf: GaugeField, flow_time : float, steps : int,
        *, c1 : float = 0.0, existing_flow_time : float = 0.0, wilson_flow_integrator_type = "runge-kutta"):
    epsilon = flow_time / steps
    energy_density_list = []
    for i in range(steps):
        gf_wilson_flow_step(gf, epsilon, c1 = c1, wilson_flow_integrator_type = wilson_flow_integrator_type)
        t = (i + 1) * epsilon + existing_flow_time
        energy_density = gf_energy_density(gf)
        energy_density_list.append(energy_density)
        displayln_info(f"gf_wilson_flow: t={t} ; E={energy_density} ; t^2 E={t*t*energy_density}")
    return energy_density_list
