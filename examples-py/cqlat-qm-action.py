#!/usr/bin/env python3

# Tests for the cqlat QMAction interface (qlat/cqlat/qm-action.cpp):
#     mk_qm_action, free_qm_action, set_qm_action,
#     get_alpha_qm_action, get_beta_qm_action, get_barrier_strength_qm_action,
#     get_M_qm_action, get_L_qm_action, get_t_FV_out_qm_action,
#     get_t_FV_mid_qm_action, get_dt_qm_action,
#     V_qm_action, dV_qm_action, action_node_qm_action,
#     hmc_m_hamilton_node_qm_action, sum_sq_qm_action, hmc_set_force_qm_action,
#     hmc_field_evolve_qm_action, hmc_set_rand_momentum_qm_action
#
# These are reached through the ``q.QMAction`` Python class (``qlat.qm_action``).
# The checks below are analytic where possible:
#   * V(x, t) is compared against an independent numpy re-implementation of the
#     piecewise potential of ``qlat/qlat/include/qlat/qm-action.h`` (all the
#     time-branch boundaries t_FV_out / t_FV_mid / dtTV are exercised),
#   * dV(x, t, idx) is compared against a central finite difference of V,
#   * sum_sq / hmc_m_hamilton_node / hmc_field_evolve are compared with numpy,
#   * action_node is compared with the global action built from V,
#   * hmc_set_force is compared against a finite difference of action_node.

import numpy as np

import qlat as q
import qlat.c as qc

check_eps = 1e-10

total_site = q.Coordinate([4, 4, 4, 8])

# The t direction must NOT be split: action_node() shifts in t inside the node
# geometry, so the periodic wrap is only correct when a single node owns the
# whole time direction.
size_node_list = [
    [2, 1, 1, 1],
    [1, 1, 1, 2],
    [1, 1, 1, 1],
]

# --- the QMAction parameters used throughout the test
qm_args = dict(
    alpha=2.0,
    beta=0.9,
    V_FV_min=0.1,
    FV_offset=0.05,
    TV_offset=0.1,
    barrier_strength=0.5,
    L=1.0,  # -> vtype = GET_M
    M=0.0,
    epsilon=0.0,
    t_FV_out=3,
    t_FV_mid=2,
    dt=0.25,
    measure_offset_L=False,
    measure_offset_M=False,
)
center_bar = 1.0 / qm_args["alpha"] ** 2

class Ref:
    pass

p = Ref()
for k, v in qm_args.items():
    setattr(p, k, v)
p.center_bar = center_bar

# --- independent re-implementation of qm-action.h (GET_M branch)

def ref_order_param(x):
    return x[0] * x[0] + x[1] * x[1]

def ref_V_phi4(x):
    return p.beta * (
        (x[0] * x[0] + x[1] * x[1]) / 2.0
        + p.alpha * (x[1] * x[0] * x[0] - x[1] ** 3 / 3.0)
    )

def ref_d_order_param(x, idx):
    return 2.0 * x[idx]

def ref_V_full_xy(x):
    rtn = ref_V_phi4(x)
    return 0.0 if rtn < -p.V_FV_min else rtn + p.V_FV_min

def ref_V_bar_max():
    return p.beta / p.alpha / p.alpha / 6.0 + p.V_FV_min

def ref_V_full_op_fixed(x, op):
    norm = (ref_order_param(x) / op) ** 0.5
    vfull = ref_V_full_xy(x)
    vfull_op = ref_V_full_xy((x[0] / norm, x[1] / norm))
    return max(vfull, vfull_op)

def ref_V_FV_out(x):
    op = ref_order_param(x)
    if op > p.center_bar + p.FV_offset:
        return (
            ref_V_full_op_fixed(x, p.center_bar + p.FV_offset)
            + p.barrier_strength * (op - p.center_bar - p.FV_offset) ** 2
        )
    return ref_V_full_xy(x)

def ref_V_FV_mid(x):
    op = ref_order_param(x)
    if op > p.center_bar:
        return (
            ref_V_full_op_fixed(x, p.center_bar)
            + p.barrier_strength * (op - p.center_bar) ** 2
        )
    return ref_V_full_xy(x)

def ref_V_TV(x):
    op = ref_order_param(x)
    if op < p.center_bar + p.TV_offset:
        return (
            ref_V_full_xy(x)
            + p.barrier_strength * (op - p.center_bar - p.TV_offset) ** 2
        )
    return ref_V_full_xy(x)

def ref_V_proj(x):
    Vbar = ref_V_FV_out(x) - ref_V_full_xy(x)
    rtn = -np.log((1.0 - np.exp(-(Vbar + p.epsilon) * p.dt)) / p.dt) / p.dt
    if ref_order_param(x) < p.center_bar + p.FV_offset:
        rtn += (
            p.barrier_strength
            * (p.center_bar + p.FV_offset - ref_order_param(x)) ** 0.5
        )
    return rtn

def ref_V_max(V_D, V_N, P):
    return V_D if V_N < V_D else (1.0 - P) * V_D + P * V_N

def ref_V_FV_floored(x, P):
    v = ref_V_FV_mid(x)
    floor = p.V_FV_min + P * (ref_V_bar_max() - p.V_FV_min)
    return floor if v < floor else v

def ref_V_M(x, V_N):
    return ref_V_max(V_N, ref_V_FV_floored(x, 1.0), p.M)

def ref_V(x, t):
    # GET_M branch of QMAction::V_t_M
    if t == 0:
        return ref_V_M(x, ref_V_full_xy(x) + ref_V_proj(x))
    elif t <= p.t_FV_out:
        return ref_V_FV_out(x)
    elif t <= p.t_FV_out + p.t_FV_mid:
        return ref_V_FV_mid(x)
    elif t <= 2 * p.t_FV_out + p.t_FV_mid:
        return ref_V_FV_out(x)
    elif t == 2 * p.t_FV_out + p.t_FV_mid + 1:
        return ref_V_M(x, ref_V_full_xy(x) + ref_V_proj(x))
    return ref_V_M(x, ref_V_TV(x))

def mk_field_from_global(geo, g0, g1):
    # pickle the local part of a *global* pattern, so that the two MPI ranks
    # hold different values (a bare RngState is rank independent)
    f = q.FieldRealD(geo, 2)
    buf = np.asarray(f)
    for index in range(geo.local_volume):
        xg = geo.coordinate_g_from_l(geo.coordinate_from_index(index))
        buf[index, 0] = g0[xg[0], xg[1], xg[2], xg[3]]
        buf[index, 1] = g1[xg[0], xg[1], xg[2], xg[3]]
    return f

q.begin_with_mpi(size_node_list)

geo = q.Geometry(total_site)
q.json_results_append(f"cqlat-qm-action: geo.show()={geo.show()}")
rs = q.RngState("cqlat-qm-action")

qma = q.QMAction(
    qm_args["alpha"],
    qm_args["beta"],
    qm_args["V_FV_min"],
    qm_args["FV_offset"],
    qm_args["TV_offset"],
    qm_args["barrier_strength"],
    qm_args["L"],
    qm_args["M"],
    qm_args["epsilon"],
    qm_args["t_FV_out"],
    qm_args["t_FV_mid"],
    qm_args["dt"],
    qm_args["measure_offset_L"],
    qm_args["measure_offset_M"],
)

# --- scalar parameters are read back exactly
q.json_results_append("cqlat-qm-action: alpha", qma.alpha(), check_eps)
q.json_results_append("cqlat-qm-action: beta", qma.beta(), check_eps)
q.json_results_append(
    "cqlat-qm-action: barrier_strength", qma.barrier_strength(), check_eps
)
q.json_results_append("cqlat-qm-action: M", qma.M(), check_eps)
q.json_results_append("cqlat-qm-action: L", qma.L(), check_eps)
q.json_results_append("cqlat-qm-action: t_FV_out", float(qma.t_FV_out()), check_eps)
q.json_results_append("cqlat-qm-action: t_FV", float(qma.t_FV()), check_eps)
q.json_results_append("cqlat-qm-action: dt", qma.dt(), check_eps)
assert qma.alpha() == qm_args["alpha"]
assert qma.beta() == qm_args["beta"]
assert qma.barrier_strength() == qm_args["barrier_strength"]
assert qma.M() == qm_args["M"]
assert qma.L() == qm_args["L"]
assert qma.t_FV_out() == qm_args["t_FV_out"]
assert qma.t_FV_mid() == qm_args["t_FV_mid"]
assert qma.t_FV() == 2 * qm_args["t_FV_out"] + qm_args["t_FV_mid"]
assert qma.dt() == qm_args["dt"]

# --- set_qm_action: qma2 @= qma copies all parameters
qma2 = q.QMAction(1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 5, 5, 1.0, False, False)
qma2 @= qma
err_set = 0.0
for f in ["alpha", "beta", "barrier_strength", "M", "L", "t_FV_out", "t_FV_mid", "dt"]:
    err_set = max(err_set, abs(float(getattr(qma2, f)()) - float(getattr(qma, f)())))
assert err_set == 0.0, err_set
q.json_results_append("cqlat-qm-action: set_qm_action", float(err_set == 0.0))
del qma2
import gc

gc.collect()

# --- V(x, t): compare against the numpy reference on all time branches.
#     For |x|^2 <= center_bar + FV_offset the t == 0 (and t == t_FV+1) value
#     involves V_proj whose log term diverges, so those x are checked only for
#     t >= 1.
x_list_hi = [(1.0, 0.6), (2.0, 1.5), (-1.2, 0.8), (0.05, -0.7)]
x_list_lo = [(0.2, 0.3), (0.35, 0.3)]
t_list_all = [0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 20]
t_list_pos = [1, 2, 3, 4, 5, 6, 8, 10, 11, 20]
check_points = [(x, t) for x in x_list_hi for t in t_list_all]
check_points += [(x, t) for x in x_list_lo for t in t_list_pos]
err_V = 0.0
for x, t in check_points:
    v = qma.V(x, t)
    vl = ref_V(x, t)
    assert np.isfinite(v), (x, t, v)
    err_V = max(err_V, float(abs(v - vl)))
assert err_V < check_eps, err_V
q.json_results_append("cqlat-qm-action: V vs reference", float(err_V < check_eps))
q.json_results_append("cqlat-qm-action: V max error", err_V, check_eps)
q.json_results_append(
    "cqlat-qm-action: V((1.0,0.6), 1)", qma.V((1.0, 0.6), 1), check_eps
)
q.json_results_append(
    "cqlat-qm-action: V((0.35,0.3), 5)", qma.V((0.35, 0.3), 5), check_eps
)
q.json_results_append(
    "cqlat-qm-action: V((2.0,1.5), 20)", qma.V((2.0, 1.5), 20), check_eps
)

# --- dV(x, t, idx): central finite difference of V (both components)
err_dV = 0.0
h = 1e-6
for x, t in check_points:
    for idx in [0, 1]:
        xp = list(x)
        xm = list(x)
        xp[idx] += h
        xm[idx] -= h
        fd = (qma.V(xp, t) - qma.V(xm, t)) / (2.0 * h)
        dv0 = qma.dV(x, t) if idx == 0 else qc.dV_qm_action(qma, x[0], x[1], t, idx)
        err_dV = max(err_dV, float(abs(fd - dv0)))
assert err_dV < 1e-7, err_dV
q.json_results_append("cqlat-qm-action: dV vs finite difference", float(err_dV < 1e-7))
q.json_results_append("cqlat-qm-action: dV max error", err_dV, 1e-7)
q.json_results_append(
    "cqlat-qm-action: dV((0.2,0.3), 1, 0)", qma.dV((0.2, 0.3), 1), check_eps
)
q.json_results_append(
    "cqlat-qm-action: dV((0.2,0.3), 1, 1)",
    qc.dV_qm_action(qma, 0.2, 0.3, 1, 1),
    check_eps,
)
q.json_results_append(
    "cqlat-qm-action: dV((2.0,1.5), 20, 1)",
    qc.dV_qm_action(qma, 2.0, 1.5, 20, 1),
    check_eps,
)

# --- field data: the local values come from a global pattern so that they
#     differ between the two MPI ranks (a bare RngState is rank independent).
#     |psi|^2 > center_bar + FV_offset everywhere, so that the t == 0 part of
#     the potential (which involves V_proj) stays finite.
latt_size = tuple(int(geo.total_site[i]) for i in range(4))
n_site = int(np.prod(latt_size))
rs_g = q.RngState("cqlat-qm-action-field")
psi0_g = 1.1 + 0.4 * rs_g.split("psi0").u_rand_arr(latt_size)
psi1_g = 0.9 + 0.4 * rs_g.split("psi1").u_rand_arr(latt_size)
m0_g = rs_g.split("m0").u_rand_arr(latt_size) - 0.5
m1_g = rs_g.split("m1").u_rand_arr(latt_size) - 0.5
h0_g = 1e-6 * (rs_g.split("h0").u_rand_arr(latt_size) - 0.5)
h1_g = 1e-6 * (rs_g.split("h1").u_rand_arr(latt_size) - 0.5)
assert float(np.min(psi0_g**2 + psi1_g**2)) > center_bar + qm_args["FV_offset"]

f = mk_field_from_global(geo, psi0_g, psi1_g)
psi_arr = np.asarray(f)

# --- sum_sq and hmc_m_hamilton_node
sum_sq_ref = float((psi0_g**2 + psi1_g**2).sum())
sum_sq = q.glb_sum(qma.sum_sq(f))
assert abs(sum_sq - sum_sq_ref) < check_eps * max(1.0, abs(sum_sq_ref))
q.json_results_append("cqlat-qm-action: sum_sq", sum_sq, 1e-10)
err_ham = abs(qma.hmc_m_hamilton_node(f) - 0.5 * qma.sum_sq(f))
assert err_ham < check_eps * max(1.0, abs(sum_sq_ref))
q.json_results_append(
    "cqlat-qm-action: hmc_m_hamilton_node - sum_sq/2", float(err_ham), 1e-10
)

# --- action_node: glb_sum(action_node(f)) equals the global action built from V
dpsi_sq = (np.roll(psi0_g, -1, axis=3) - psi0_g) ** 2 + (
    np.roll(psi1_g, -1, axis=3) - psi1_g
) ** 2
V_g = np.zeros(latt_size, dtype=np.float64)
for x0 in range(latt_size[0]):
    for x1 in range(latt_size[1]):
        for x2 in range(latt_size[2]):
            for t in range(latt_size[3]):
                V_g[x0, x1, x2, t] = qma.V(
                    (psi0_g[x0, x1, x2, t], psi1_g[x0, x1, x2, t]), t
                )
S_ref = qm_args["dt"] * float(
    np.sum(qm_args["beta"] / (2.0 * qm_args["dt"] ** 2) * dpsi_sq + V_g)
    - n_site * np.log(qm_args["dt"]) / qm_args["dt"]
)
S = q.glb_sum(qma.action_node(f))
assert abs(S - S_ref) < 1e-9 * max(1.0, abs(S_ref)), (S, S_ref)
q.json_results_append("cqlat-qm-action: action_node", S, 1e-9)
q.json_results_append(
    "cqlat-qm-action: action_node vs reference", float(abs(S - S_ref)), 1e-9
)

# --- hmc_set_force: the force is the gradient of the action
force = mk_field_from_global(geo, np.zeros(latt_size), np.zeros(latt_size))
qma.hmc_set_force(force, f)
h_field = mk_field_from_global(geo, h0_g, h1_g)
f_p = mk_field_from_global(geo, psi0_g + h0_g, psi1_g + h1_g)
f_m = mk_field_from_global(geo, psi0_g - h0_g, psi1_g - h1_g)
dS = q.glb_sum(qma.action_node(f_p)) - q.glb_sum(qma.action_node(f_m))
pred = 2.0 * q.glb_sum(float(np.sum(np.asarray(force) * np.asarray(h_field))))
err_force = abs(dS - pred) / max(1.0, abs(pred))
assert err_force < 1e-8, (dS, pred, err_force)
q.json_results_append("cqlat-qm-action: hmc_set_force vs dS/dh", float(err_force), 1e-8)
q.json_results_append(
    "cqlat-qm-action: hmc_set_force sum",
    q.glb_sum(float(np.asarray(force).sum())),
    1e-8,
)

# --- hmc_field_evolve: f += m * step_size
m_field = mk_field_from_global(geo, m0_g, m1_g)
f_ev = mk_field_from_global(geo, psi0_g, psi1_g)
step_size = 0.37
qma.hmc_field_evolve(f_ev, m_field, step_size)
err_ev = float(
    np.max(np.abs(np.asarray(f_ev) - (np.asarray(f) + np.asarray(m_field) * step_size)))
)
assert err_ev < check_eps * max(1.0, float(np.max(np.abs(psi_arr)))), err_ev
q.json_results_append("cqlat-qm-action: hmc_field_evolve", float(err_ev < check_eps))
q.json_results_append(
    "cqlat-qm-action: hmc_field_evolve sum",
    q.glb_sum(float(np.asarray(f_ev).sum())),
    1e-8,
)

# --- hmc_set_rand_momentum: deterministic in the RngState, gaussian values
m1 = mk_field_from_global(geo, np.zeros(latt_size), np.zeros(latt_size))
m2 = mk_field_from_global(geo, np.zeros(latt_size), np.zeros(latt_size))
qma.hmc_set_rand_momentum(m1, q.RngState("cqlat-qm-action-rs"))
qma.hmc_set_rand_momentum(m2, q.RngState("cqlat-qm-action-rs"))
assert np.array_equal(np.asarray(m1), np.asarray(m2))
q.json_results_append("cqlat-qm-action: hmc_set_rand_momentum reproducible", 1.0)
m1_arr = np.asarray(m1)
mean_g = q.glb_sum(float(m1_arr.sum())) / (n_site * 2)
var_g = q.glb_sum(float((m1_arr**2).sum())) / (n_site * 2) - mean_g**2
q.json_results_append(
    "cqlat-qm-action: hmc_set_rand_momentum mean", float(mean_g), 1e-3
)
q.json_results_append(
    "cqlat-qm-action: hmc_set_rand_momentum std", float(var_g**0.5), 1e-3
)
assert abs(mean_g) < 0.05
assert abs(var_g**0.5 - 1.0) < 0.05

# --- free_qm_action: exercised by __del__ (with an explicit gc.collect())
qma_tmp = q.QMAction(
    qm_args["alpha"],
    qm_args["beta"],
    qm_args["V_FV_min"],
    qm_args["FV_offset"],
    qm_args["TV_offset"],
    qm_args["barrier_strength"],
    qm_args["L"],
    qm_args["M"],
    qm_args["epsilon"],
    qm_args["t_FV_out"],
    qm_args["t_FV_mid"],
    qm_args["dt"],
    qm_args["measure_offset_L"],
    qm_args["measure_offset_M"],
)
qma_tmp @= qma
del qma_tmp
gc.collect()
q.json_results_append("cqlat-qm-action: free_qm_action (via __del__)", 1.0)

del qma
gc.collect()

q.timer_display()
if q.is_test():
    q.check_log_json(__file__)
q.end_with_mpi()
q.displayln_info("CHECK: finished successfully.")
