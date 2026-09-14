#!/usr/bin/env python3

# Tests for the cqlat InverterDomainWall interface (qlat/cqlat/inverter.cpp):
#     free_inverter_domain_wall,
#     get_stop_rsd_inverter_domain_wall, set_stop_rsd_inverter_domain_wall,
#     get_max_num_iter_inverter_domain_wall, set_max_num_iter_inverter_domain_wall,
#     get_max_mixed_precision_cycle_inverter_domain_wall,
#     set_max_mixed_precision_cycle_inverter_domain_wall
#
# These are reached both directly through ``qlat.c`` (``qc.<name>(...)``) and
# through the Python class ``q.InverterDomainWall`` (qlat/qlat/inverter.pyx).
#
# The setters in inverter.cpp assign the requested value verbatim: there is no
# clamping and no range assertion, so every value must round-trip exactly.
# The defaults come from ``InverterParams::init()`` in
# qlat/qlat/include/qlat/dslash.h.
#
# The solver parameters are checked functionally as well: with a single CG
# cycle / iteration the inversion is far from converged, while the default
# parameters reproduce the analytic free-field DWF inverse for unit links.

import gc

import numpy as np

import qlat as q
import qlat.c as qc

# --- defaults from InverterParams::init() (qlat/qlat/include/qlat/dslash.h)
default_stop_rsd = 1e-8
default_max_num_iter = 200
default_max_mixed_precision_cycle = 300

check_eps = 1e-10

def sol_diff_qnorm(sol, ref):
    # distance between two propagators (global norm)
    diff = sol.copy()
    diff -= ref
    return float(diff.qnorm())

# the runner uses ``mpiexec -n 2``, so the [2, 1, 1, 1] layout is selected
q.begin_with_mpi([[2, 1, 1, 1], [1, 1, 1, 1]])

geo = q.Geometry(q.Coordinate([4, 4, 4, 8]))
q.json_results_append(f"cqlat-dwf-inverter: geo.show()={geo.show()}")

# --- unit links keep the inversion analytically known, so that the CG result
#     can be compared with q.InverterDwfFreeField
gf = q.GaugeField(geo)
gf.set_unit()

fa = q.FermionAction(mass=0.1, ls=8, m5=1.0, mobius_scale=1.0)

inv = q.InverterDomainWall(gf=gf, fa=fa)

# --- the three getter defaults (Python wrapper and direct cqlat call)
assert isinstance(inv.stop_rsd(), float)
assert isinstance(inv.max_num_iter(), int)
assert isinstance(inv.max_mixed_precision_cycle(), int)
assert inv.stop_rsd() == default_stop_rsd
assert inv.max_num_iter() == default_max_num_iter
assert inv.max_mixed_precision_cycle() == default_max_mixed_precision_cycle
assert qc.get_stop_rsd_inverter_domain_wall(inv) == default_stop_rsd
assert qc.get_max_num_iter_inverter_domain_wall(inv) == default_max_num_iter
assert (
    qc.get_max_mixed_precision_cycle_inverter_domain_wall(inv)
    == default_max_mixed_precision_cycle
)
q.json_results_append(
    "cqlat-dwf-inverter: default stop_rsd", float(inv.stop_rsd()), check_eps
)
q.json_results_append(
    "cqlat-dwf-inverter: default max_num_iter", float(inv.max_num_iter()), check_eps
)
q.json_results_append(
    "cqlat-dwf-inverter: default max_mixed_precision_cycle",
    float(inv.max_mixed_precision_cycle()),
    check_eps,
)

# --- Python wrapper setters round-trip through the cqlat setters/getters
inv.set_stop_rsd(1e-3)
assert inv.stop_rsd() == 1e-3
assert qc.get_stop_rsd_inverter_domain_wall(inv) == 1e-3
inv.set_max_num_iter(3)
assert inv.max_num_iter() == 3
assert qc.get_max_num_iter_inverter_domain_wall(inv) == 3
inv.set_max_mixed_precision_cycle(7)
assert inv.max_mixed_precision_cycle() == 7
assert qc.get_max_mixed_precision_cycle_inverter_domain_wall(inv) == 7
q.json_results_append(
    "cqlat-dwf-inverter: set_stop_rsd(1e-3)", float(inv.stop_rsd()), check_eps
)
q.json_results_append(
    "cqlat-dwf-inverter: set_max_num_iter(3)", float(inv.max_num_iter()), check_eps
)
q.json_results_append(
    "cqlat-dwf-inverter: set_max_mixed_precision_cycle(7)",
    float(inv.max_mixed_precision_cycle()),
    check_eps,
)

# --- direct cqlat setters behave identically
qc.set_stop_rsd_inverter_domain_wall(inv, 1e-6)
assert inv.stop_rsd() == 1e-6
qc.set_max_num_iter_inverter_domain_wall(inv, 11)
assert inv.max_num_iter() == 11
qc.set_max_mixed_precision_cycle_inverter_domain_wall(inv, 13)
assert inv.max_mixed_precision_cycle() == 13
q.json_results_append(
    "cqlat-dwf-inverter: qc.set_stop_rsd(1e-6)", float(inv.stop_rsd()), check_eps
)
q.json_results_append(
    "cqlat-dwf-inverter: qc.set_max_num_iter(11)", float(inv.max_num_iter()), check_eps
)
q.json_results_append(
    "cqlat-dwf-inverter: qc.set_max_mixed_precision_cycle(13)",
    float(inv.max_mixed_precision_cycle()),
    check_eps,
)

# --- no clamping and no range assertion: out-of-range values are stored as-is
inv.set_stop_rsd(2.5)
assert inv.stop_rsd() == 2.5
inv.set_stop_rsd(-1.0)
assert inv.stop_rsd() == -1.0
inv.set_max_num_iter(0)
assert inv.max_num_iter() == 0
inv.set_max_num_iter(-5)
assert inv.max_num_iter() == -5
inv.set_max_mixed_precision_cycle(0)
assert inv.max_mixed_precision_cycle() == 0
inv.set_max_mixed_precision_cycle(-9)
assert inv.max_mixed_precision_cycle() == -9
q.json_results_append(
    "cqlat-dwf-inverter: setters do not clamp", float(inv.stop_rsd() == -1.0)
)
q.json_results_append(
    "cqlat-dwf-inverter: set_max_num_iter(-5) stored", float(inv.max_num_iter())
)
q.json_results_append(
    "cqlat-dwf-inverter: set_max_mixed_precision_cycle(-9) stored",
    float(inv.max_mixed_precision_cycle()),
)

# --- restore the defaults
inv.set_stop_rsd(default_stop_rsd)
inv.set_max_num_iter(default_max_num_iter)
inv.set_max_mixed_precision_cycle(default_max_mixed_precision_cycle)
assert inv.stop_rsd() == default_stop_rsd
assert inv.max_num_iter() == default_max_num_iter
assert inv.max_mixed_precision_cycle() == default_max_mixed_precision_cycle

# --- point source (deterministic; no rank-independent random data is used)
src = q.mk_point_src(geo, q.Coordinate([0, 0, 0, 0]))
qinv_free = q.InverterDwfFreeField(mass=fa.mass(), m5=fa.m5())

# --- the solver parameters are actually used: a single CG iteration cannot
#     converge, while the defaults reproduce the free-field inverse
inv.set_max_num_iter(1)
inv.set_max_mixed_precision_cycle(1)
sol_coarse = inv * src
err_coarse = sol_diff_qnorm(sol_coarse, qinv_free * src)
assert np.isfinite(err_coarse) and err_coarse > 1e-3, err_coarse
q.json_results_append(
    "cqlat-dwf-inverter: 1 iteration is not converged", float(err_coarse > 1e-3)
)

inv.set_max_num_iter(default_max_num_iter)
inv.set_max_mixed_precision_cycle(default_max_mixed_precision_cycle)

# ``InverterDwfFreeField`` is the free DWF inverse for infinite L_s
# (``free_mom_invert`` in qlat/qlat/include/qlat/qcd-prop.h), so for the finite
# ls=8 action the two agree only up to the finite-ls truncation (~2e-6), while
# the CG itself reaches a DWF-equation residual of ~4e-9.
sol = inv * src
sol_free = qinv_free * src
q.json_results_append("cqlat-dwf-inverter: sol qnorm", float(sol.qnorm()), 1e-6)
q.json_results_append(
    "cqlat-dwf-inverter: sol_free qnorm", float(sol_free.qnorm()), 1e-6
)
err_inv = sol_diff_qnorm(sol, sol_free)
assert np.isfinite(err_inv) and err_inv < 1e-5, err_inv
q.json_results_append(
    "cqlat-dwf-inverter: inversion matches free field", float(err_inv < 1e-5)
)

# --- free_inverter_domain_wall: exercised by __del__ + gc.collect()
inv_tmp = q.InverterDomainWall(gf=gf, fa=fa)
inv_tmp.set_stop_rsd(1e-4)
assert inv_tmp.stop_rsd() == 1e-4
del inv_tmp
gc.collect()
q.json_results_append("cqlat-dwf-inverter: free_inverter_domain_wall (__del__)", 1.0)

del inv
gc.collect()

q.timer_display()
if q.is_test():
    q.check_log_json(__file__)
q.end_with_mpi()
q.displayln_info("CHECK: finished successfully.")
