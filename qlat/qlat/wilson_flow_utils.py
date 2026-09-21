"""
Module ``qlat.wilson_flow_utils``
==================================\n
Wilson flow and stout smearing drivers that do not call C++ functions
directly; they compose the Cython-wrapped flow steps.\n
"""

class q:
    from qlat_utils import (
        timer,
        get_fname,
        displayln_info,
        Coordinate,
    )
    from .wilson_flow import (
        gf_wilson_flow_step,
        gf_energy_density,
        gf_block_stout_smear,
        gf_wilson_flow_force,
    )
    from .hmc import (
        gf_evolve,
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

@q.timer
def gf_stout_smear(gf, step_size, num_step=1, *, method=None):
    """
    Apply stout smearing to gf in place for num_step steps.
    method: None or "force" uses gf_wilson_flow_force + gf_evolve directly (default, fastest).
            "stout" uses gf_block_stout_smear.
            "wilson-flow" uses gf_wilson_flow_step with euler integrator and c1=0.
    """
    fname = q.get_fname()
    if method is None:
        method = "force"
    if method == "stout":
        for step in range(num_step):
            q.gf_block_stout_smear(gf, q.Coordinate(), step_size)
    elif method == "wilson-flow":
        for step in range(num_step):
            q.gf_wilson_flow_step(
                gf, step_size, c1=0.0, wilson_flow_integrator_type="euler"
            )
    elif method == "force":
        for step in range(num_step):
            gm = q.gf_wilson_flow_force(gf, c1=0.0)
            q.gf_evolve(gf, gm, step_size)
    else:
        raise ValueError(f"{fname}: unknown method={method}")
    q.displayln_info(
        f"{fname}: method={method}, num_step={num_step}, step_size={step_size}, plaq={gf.plaq()}"
    )
