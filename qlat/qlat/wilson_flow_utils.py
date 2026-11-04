"""
Module ``qlat.wilson_flow_utils``
==================================\n
Wilson flow drivers that do not call C++ functions directly; they compose the
Cython-wrapped flow steps.\n
Documentation: ``docs/qlat/qlat_wilson_flow.md``\n
.. note:: Update the documentation when updating this source file.
"""

class q:
    from qlat_utils import (
        timer,
        get_fname,
        displayln_info,
    )
    from .wilson_flow import (
        gf_wilson_flow_step,
        gf_energy_density,
    )

@q.timer
def gf_wilson_flow(
    gf,
    flow_time,
    steps,
    *,
    c1=0.0,
    existing_flow_time=0.0,
    wilson_flow_integrator_type=None,
):
    fname = q.get_fname()
    epsilon = flow_time / steps
    energy_density_list = []
    for i in range(steps):
        q.gf_wilson_flow_step(
            gf, epsilon, c1=c1, wilson_flow_integrator_type=wilson_flow_integrator_type
        )
        t = (i + 1) * epsilon + existing_flow_time
        energy_density = q.gf_energy_density(gf)
        energy_density_list.append(energy_density)
        q.displayln_info(
            f"{fname}: t={t} ; E={energy_density} ; t^2 E={t * t * energy_density}"
        )
    return energy_density_list
