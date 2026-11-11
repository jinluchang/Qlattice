#!/usr/bin/env python3

import qlat as q
import numpy as np
from mpi4py import MPI

def get_link_norm_sqr_arr(gf):
    # Tr(U^dag U) = sum_{i,j} |U_ij|^2 for every link, shape (local_volume, DIMN)
    arr = np.asarray(gf)
    n2_arr = np.sum(np.abs(arr) ** 2, axis=(2, 3))
    return n2_arr

def get_glb_avg_link_norm_sqr(gf):
    # Global <N^2> = (1/(V*DIMN)) sum_{x,mu} Tr(U_mu(x)^dag U_mu(x)), MPI-reduced
    n2_sum = float(np.sum(get_link_norm_sqr_arr(gf)))
    n_link = int(gf.geo.total_volume) * int(gf.multiplicity)
    return float(q.glb_sum(n2_sum)) / n_link

def get_glb_max_link_norm_sqr(gf):
    # Global max_{x,mu} Tr(U_mu(x)^dag U_mu(x)), MPI-reduced with MPI.MAX.
    # A max is exact, so the step sequence of `gf_ape_smear_no_proj` does not
    # depend on the node decomposition.
    local_max = float(np.max(get_link_norm_sqr_arr(gf)))
    comm = q.get_comm()
    return float(comm.allreduce(local_max, op=MPI.MAX))

def get_glb_max_unitarity_dev(gf):
    # Global max_{x,mu} |U_mu(x)^dag U_mu(x) - 1|, MPI-reduced with MPI.MAX
    arr = np.asarray(gf)
    udag_u = np.einsum("...ki,...kj->...ij", arr.conj(), arr)
    local_max = float(np.max(np.abs(udag_u - np.eye(3))))
    comm = q.get_comm()
    return float(comm.allreduce(local_max, op=MPI.MAX))

def get_glb_max_det_dev(gf):
    # Global max_{x,mu} |det U_mu(x) - 1|, MPI-reduced with MPI.MAX
    local_max = float(np.max(np.abs(np.linalg.det(np.asarray(gf)) - 1.0)))
    comm = q.get_comm()
    return float(comm.allreduce(local_max, op=MPI.MAX))

def get_max_rel_dev(arr, arr_ref):
    # max |arr - arr_ref| normalized by max |arr_ref|
    return float(np.max(np.abs(arr - arr_ref)) / np.max(np.abs(arr_ref)))

def replicate_no_proj_target_loop(gf, alpha, target):
    # Hand-rolled replication of the target-driven loop of
    # `gf_ape_smear_no_proj(gf, alpha, accumulated_alpha=target)`: per step
    # `M_k = global max Tr(U^dag U)`, `factor_k = 3 / M_k`,
    # `a_k = min(alpha * factor_k, target - acc)` and the raw unprojected
    # update `U' = U + (a_k / 6) staple`, with no rescaling.
    # Each step is one call with `accumulated_alpha = a_k`, and the loop stops
    # when the accumulated alpha makes no further floating-point progress.
    # Returns (field, accumulated alpha, list of a_k, list of M_k).
    gf_rep = gf.copy()
    acc = 0.0
    a_list = []
    m_list = []
    # the 1000 iteration cap only guards against a pathological stagnation
    for _ in range(1000):
        remaining = target - acc
        if not (remaining > 0.0):
            break
        m_sq = get_glb_max_link_norm_sqr(gf_rep)
        if not (m_sq > 0.0) or not np.isfinite(m_sq):
            break
        a_k = min(alpha * (3.0 / m_sq), remaining)
        if not (a_k > 0.0):
            break
        gf_rep = q.gf_ape_smear_no_proj(gf_rep, alpha, a_k)
        acc_prev = acc
        acc += a_k
        a_list.append(a_k)
        m_list.append(m_sq)
        if not (acc > acc_prev):
            break
    return gf_rep, acc, a_list, m_list

size_node_list = [
    [1, 1, 1, 1],
    [1, 1, 1, 2],
    [1, 1, 1, 4],
    [1, 1, 1, 8],
    [2, 2, 2, 2],
    [2, 2, 2, 4],
]

q.begin_with_mpi(size_node_list)

rs = q.RngState("seed")
total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
gf = q.GaugeField(geo)
gf.set_rand(rs.split("gf-init"), 0.5, 10)

q.json_results_append(f"total_site={total_site}")

# ---- gf_flow_topo with all flow types ----

flow_types = [
    "Wilson",
    "Iwasaki",
    "DBW2",
    "Freeze",
    "Shrink",
    "Localize",
    "Preserve",
]

for flow_type in flow_types:
    gf_test = gf.copy()
    for i in range(3):
        q.gf_flow_topo(gf_test, 0.1, flow_type)
    q.json_results_append(f"gf_flow_topo {flow_type} plaq", gf_test.plaq(), 1e-10)

# ---- gf_flow_topo with numeric c1 ----

for c1 in [0.0, -0.331, -1.4008]:
    gf_test = gf.copy()
    for i in range(3):
        q.gf_flow_topo(gf_test, 0.1, c1)
    q.json_results_append(f"gf_flow_topo c1={c1} plaq", gf_test.plaq(), 1e-10)

# ---- gf_flow_topo integrator types ----

for int_type in ["euler", "runge-kutta"]:
    gf_test = gf.copy()
    for i in range(3):
        q.gf_flow_topo(gf_test, 0.1, "Wilson", int_type)
    q.json_results_append(f"gf_flow_topo Wilson {int_type} plaq", gf_test.plaq(), 1e-10)

# ---- gf_wilson_flow_step ----

for int_type in ["runge-kutta", "euler"]:
    gf_test = gf.copy()
    q.gf_wilson_flow_step(gf_test, 0.1, c1=0.0, wilson_flow_integrator_type=int_type)
    q.json_results_append(f"gf_wilson_flow_step {int_type} plaq", gf_test.plaq(), 1e-10)

for c1 in [0.0, -0.331, -1.4008]:
    gf_test = gf.copy()
    q.gf_wilson_flow_step(gf_test, 0.1, c1=c1)
    q.json_results_append(f"gf_wilson_flow_step c1={c1} plaq", gf_test.plaq(), 1e-10)

# ---- gf_wilson_flow (multi-step) ----

gf_wf = gf.copy()
ed_list = q.gf_wilson_flow(gf_wf, 0.3, 3, c1=0.0)
q.json_results_append("gf_wilson_flow final plaq", gf_wf.plaq(), 1e-10)
q.json_results_append("gf_wilson_flow energy_list", np.array(ed_list), 1e-10)

# ---- gf_flow_scale ----

for int_type in ["euler", "runge-kutta"]:
    gf_test = gf.copy()
    q.gf_flow_scale(gf_test, 0.1, integrator_type=int_type)
    q.json_results_append(f"gf_flow_scale {int_type} plaq", gf_test.plaq(), 1e-10)

for is_spatial in [False, True]:
    gf_test = gf.copy()
    q.gf_flow_scale(gf_test, 0.1, is_spatial=is_spatial)
    q.json_results_append(
        f"gf_flow_scale is_spatial={is_spatial} plaq", gf_test.plaq(), 1e-10
    )

for t_dir in [0, 1, 2, 3]:
    gf_test = gf.copy()
    q.gf_flow_scale(gf_test, 0.1, t_dir=t_dir)
    q.json_results_append(f"gf_flow_scale t_dir={t_dir} plaq", gf_test.plaq(), 1e-10)

# ---- gf_ape_smear ----

for steps in [1, 3, 5]:
    gf_test = q.gf_ape_smear(gf, 0.5, steps)
    q.json_results_append(f"gf_ape_smear steps={steps} plaq", gf_test.plaq(), 1e-10)

for alpha in [0.25, 0.5, 0.75]:
    gf_test = q.gf_ape_smear(gf, alpha, 3)
    q.json_results_append(f"gf_ape_smear alpha={alpha} plaq", gf_test.plaq(), 1e-10)

# ---- gf_ape_smear_spatial ----

for steps in [1, 3, 5]:
    gf_test = q.gf_ape_smear_spatial(gf, 0.5, steps)
    q.json_results_append(
        f"gf_spatial_ape_smear steps={steps} plaq", gf_test.plaq(), 1e-10
    )

# ---- gf_hyp_smear ----

alpha_sets = [
    [0.75, 0.6, 0.3],
    [0.5, 0.5, 0.5],
    [0.9, 0.7, 0.4],
]
for a1, a2, a3 in alpha_sets:
    gf_test = q.gf_hyp_smear(gf, a1, a2, a3)
    q.json_results_append(f"gf_hyp_smear {a1},{a2},{a3} plaq", gf_test.plaq(), 1e-10)

# ---- gf_ape_smear_no_proj (v5: target accumulated alpha, additive kernel) ----
# Per step: `M_k = global max Tr(U^dag U)`, `factor_k = 3 / M_k`,
# `a_k = min(alpha * factor_k, target - acc)` and
# `U' = U + (a_k / 6) staple`, with NO rescaling of the field.
# The link term is not damped, so the field grows: `M_k > 3`, hence
# `factor_k < 1` and `a_k <= alpha`.  The per-step weight decays roughly
# geometrically, so `sum(a_k)` converges and a large `target` is never
# reached (the loop exits through its no-progress guard instead).

# `target` picks a total effective smearing (not a step count) and stays
# below 1, as every `accumulated_alpha` in the tests must; both values are
# reachable at alpha = 0.1 within a handful of steps.
for target in [0.1, 0.3]:
    gf_test = q.gf_ape_smear_no_proj(gf, 0.1, target)
    q.json_results_append(
        f"gf_ape_smear_no_proj target={target} plaq", gf_test.plaq(), 1e-10
    )

no_proj_alpha = 0.1
# 0.9 is comfortably reachable for this field at alpha = 0.1 (the accumulated
# weight converges to about 1.20, so there is ample headroom below the
# ceiling), and it takes many decaying steps rather than the nine that
# `target / alpha` would suggest.  The target stays below 1, as every
# `accumulated_alpha` in the tests must.
no_proj_target = 0.9
gf_no_proj = gf.copy()
arr_no_proj_before = np.asarray(gf_no_proj).copy()
gf_target = q.gf_ape_smear_no_proj(gf_no_proj, no_proj_alpha, no_proj_target)

# (a) target reached (the last step is truncated to the remaining target, so
# only the final rounding can leave a residue): replicate the loop in Python,
# taking `M_k` from the intermediate fields with an MPI.MAX reduction, then
# compare both the accounting and the field it produces against the single
# call.
gf_rep, acc_rep, a_list_rep, m_list_rep = replicate_no_proj_target_loop(
    gf_no_proj, no_proj_alpha, no_proj_target
)
is_target_reached = bool(
    abs(acc_rep - no_proj_target) <= 1e-13 * max(1.0, no_proj_target)
)
q.json_results_append(
    f"gf_ape_smear_no_proj target reached to 1e-13 = {is_target_reached}"
)
assert is_target_reached, acc_rep
target_loop_dev = get_max_rel_dev(
    np.asarray(gf_rep).copy(), np.asarray(gf_target).copy()
)
is_target_loop_equal = bool(target_loop_dev < 1e-13)
q.json_results_append(
    "gf_ape_smear_no_proj target call vs replicated loop"
    f" rel dev < 1e-13 = {is_target_loop_equal}"
)
assert is_target_loop_equal, target_loop_dev

# (b) unitary input: a unitary input has `M_0 = 3`, so `factor_0 = 1` and the
# first step is the plain unprojected step of weight `a_0 = min(alpha, target)`;
# in particular `target = alpha` is one plain unprojected update of weight
# `alpha`.  For `n > 1` the bookkeeping relation `target = alpha * n` is NOT a
# trajectory equivalence: the loop recomputes `M_k` every step, and with the
# additive kernel the field grows, so `factor_k < 1` and the target is reached
# in more, smaller steps (see the step count below).
gf_uni = gf.copy()
gf_uni.unitarize()
m_uni = get_glb_max_link_norm_sqr(gf_uni)
q.json_results_append(
    "gf_ape_smear_no_proj unitary input max Tr(U^dag U)", m_uni, 1e-10
)
is_unitary_max = bool(abs(m_uni - 3.0) <= 1e-13)
q.json_results_append(
    f"gf_ape_smear_no_proj unitary input M_0 = 3 to 1e-13 = {is_unitary_max}"
)
assert is_unitary_max, m_uni
gf_uni_rep, _, a_list_uni, _ = replicate_no_proj_target_loop(
    gf_uni, no_proj_alpha, no_proj_alpha
)
assert len(a_list_uni) >= 1, a_list_uni
is_one_plain_step = bool(
    abs(a_list_uni[0] - no_proj_alpha) < 1e-13 and sum(a_list_uni[1:]) < 1e-13
)
q.json_results_append(
    "gf_ape_smear_no_proj unitary target=alpha is one plain update"
    f" of weight alpha = {is_one_plain_step}"
)
assert is_one_plain_step, a_list_uni
gf_uni_one = q.gf_ape_smear_no_proj(gf_uni, no_proj_alpha, no_proj_alpha)
uni_dev = get_max_rel_dev(np.asarray(gf_uni_rep).copy(), np.asarray(gf_uni_one).copy())
is_uni_match = bool(uni_dev < 1e-13)
q.json_results_append(
    "gf_ape_smear_no_proj unitary target=alpha call vs replica"
    f" rel dev < 1e-13 = {is_uni_match}"
)
assert is_uni_match, uni_dev

# (c) accounting: the replica derives every step as
# `a_k = min(alpha * 3 / M_k, remaining)` and truncates the last one to the
# remaining target; (a) shows that the implementation follows exactly that
# sequence.  The `a_k` are recorded as a bookkeeping regression.  This is a
# bookkeeping equivalence only: each call truncates its own last step, so two
# successive calls are not interchangeable with one call of the summed target.
q.json_results_append(
    "gf_ape_smear_no_proj target a_k sequence", np.array(a_list_rep), 1e-10
)

# (d) `target = 0` and `alpha = 0` are bit-exact identities, and the input
# field is never modified.
gf_target_zero = q.gf_ape_smear_no_proj(gf_no_proj, no_proj_alpha, 0.0)
is_target_zero_identity = bool(
    np.array_equal(np.asarray(gf_target_zero), arr_no_proj_before)
)
q.json_results_append(
    f"gf_ape_smear_no_proj target=0 bitwise identity = {is_target_zero_identity}"
)
assert is_target_zero_identity
gf_alpha_zero = q.gf_ape_smear_no_proj(gf_no_proj, 0.0, no_proj_target)
is_alpha_zero_identity = bool(
    np.array_equal(np.asarray(gf_alpha_zero), arr_no_proj_before)
)
q.json_results_append(
    f"gf_ape_smear_no_proj alpha=0 bitwise identity = {is_alpha_zero_identity}"
)
assert is_alpha_zero_identity
assert np.array_equal(np.asarray(gf_no_proj), arr_no_proj_before)

# (e) informational: step count, accumulating `a_k`, the range of `M_k`, and
# the invariant that the additive kernel implies.  Because the field grows,
# `M_k > 3`, so `factor_k < 1` and no step can request more than `alpha`; the
# target is reached in more, smaller steps than `target / alpha`.
q.json_results_append(
    f"gf_ape_smear_no_proj step count = {len(a_list_rep)}"
    f" (vs target/alpha = {no_proj_target / no_proj_alpha})"
)
q.json_results_append("gf_ape_smear_no_proj max a_k", max(a_list_rep), 1e-10)
is_a_capped = bool(max(a_list_rep) <= no_proj_alpha * (1.0 + 1e-12))
q.json_results_append(
    f"gf_ape_smear_no_proj max a_k <= alpha (1 + 1e-12) = {is_a_capped}"
)
assert is_a_capped, (max(a_list_rep), no_proj_alpha)
q.json_results_append("gf_ape_smear_no_proj M_k min", min(m_list_rep), 1e-10)
q.json_results_append("gf_ape_smear_no_proj M_k max", max(m_list_rep), 1e-10)
is_m_grows = bool(m_list_rep[-1] > m_list_rep[0])
q.json_results_append(f"gf_ape_smear_no_proj M_k grows (last > first) = {is_m_grows}")
assert is_m_grows, (m_list_rep[0], m_list_rep[-1])
q.json_results_append(
    "gf_ape_smear_no_proj global avg <N^2> after target",
    get_glb_avg_link_norm_sqr(gf_target),
    1e-10,
)

# (g) the additive kernel itself: on a constant unit field the staple is
# exactly `6 U`, so one plain step is `U' = U + (a / 6) * 6 U = (1 + a) U`.
# `alpha = target = 0.75` is a deliberately special exact-arithmetic value,
# not the test's default alpha of 0.1: `a / 6 = 1 / 8` is exact, so the result
# is exactly `1.75 U`; the old `(1 - a) U + (a / 6) staple` kernel would leave
# the unit field unchanged.  With alpha = 0.1 the identity would be `1.1 U`,
# but `(0.1 / 6) * 6` is not exactly 0.1 in double precision, so an exact
# comparison would be fragile.
gf_unit = q.GaugeField(geo)
gf_unit.set_unit()
m_unit = get_glb_max_link_norm_sqr(gf_unit)
q.json_results_append("gf_ape_smear_no_proj unit field max Tr(U^dag U)", m_unit, 1e-13)
assert abs(m_unit - 3.0) <= 1e-13, m_unit
gf_unit_step = q.gf_ape_smear_no_proj(gf_unit, 0.75, 0.75)
arr_unit_step = np.asarray(gf_unit_step)
expected_unit_step = np.zeros_like(arr_unit_step)
for mu in range(4):
    expected_unit_step[..., mu, :, :] = 1.75 * np.eye(3)
kernel_dev = float(np.max(np.abs(arr_unit_step - expected_unit_step)))
is_additive_kernel = bool(kernel_dev < 1e-13)
q.json_results_append(
    "gf_ape_smear_no_proj unit field one step is 1.75 U"
    f" (special exact alpha = 0.75, additive kernel) = {is_additive_kernel}"
)
assert is_additive_kernel, kernel_dev

# ---- gf_ape_smear_no_proj_rescaling (additive kernel + per-step rescale) ----
# No SU(3) projection: each step is the additive kernel `U' = U + (alpha / 6) S`
# (the same step kernel as `gf_ape_smear_no_proj`, without its `3 / M`
# pre-scaling), followed by the global scalar rescale
# `U -> sqrt(3 / <N^2>) U` with
# `<N^2> = qnorm(gf) / (V * DIMN)` (MPI-reduced), so that `<N^2> = 3` (the
# unitary value) after every step.  Unlike `gf_ape_smear`, the field is
# therefore not unitary after a step.

rescale_alpha = 0.5
gf_rescale_sig = np.asarray(gf).copy()

# (a) one unprojected step: `<N^2>` is restored to 3, but the field is far
# from SU(3) (`|U^dag U - 1|` and `|det U - 1|` are O(1)); `gf_ape_smear`
# instead projects back to SU(3) after every step.
gf_rescale_1 = q.gf_ape_smear_no_proj_rescaling(gf, rescale_alpha, 1)
avg_rescale_1 = get_glb_avg_link_norm_sqr(gf_rescale_1)
q.json_results_append(
    "gf_ape_smear_no_proj_rescaling 1 step global avg <N^2>", avg_rescale_1, 1e-10
)
is_avg_3 = bool(abs(avg_rescale_1 - 3.0) <= 1e-13)
q.json_results_append(
    f"gf_ape_smear_no_proj_rescaling 1 step avg <N^2> = 3 to 1e-13 = {is_avg_3}"
)
assert is_avg_3, avg_rescale_1
max_ut_dev_1 = get_glb_max_unitarity_dev(gf_rescale_1)
max_det_dev_1 = get_glb_max_det_dev(gf_rescale_1)
q.json_results_append(
    "gf_ape_smear_no_proj_rescaling 1 step global max |U^dag U - 1|",
    max_ut_dev_1,
    1e-10,
)
q.json_results_append(
    "gf_ape_smear_no_proj_rescaling 1 step global max |det U - 1|",
    max_det_dev_1,
    1e-10,
)
is_not_unitary = bool(max_ut_dev_1 > 0.5 and max_det_dev_1 > 0.5)
q.json_results_append(
    "gf_ape_smear_no_proj_rescaling 1 step max |U^dag U - 1| > 0.5"
    f" and max |det U - 1| > 0.5 = {is_not_unitary}"
)
assert is_not_unitary, (max_ut_dev_1, max_det_dev_1)
q.json_results_append(
    "gf_ape_smear_no_proj_rescaling 1 step plaq", gf_rescale_1.plaq(), 1e-10
)
assert abs(gf_rescale_1.plaq() - 0.30241595774022034) <= 1e-10
max_ut_dev_ape = get_glb_max_unitarity_dev(q.gf_ape_smear(gf, rescale_alpha, 1))
is_ape_unitary = bool(max_ut_dev_ape < 1e-13)
q.json_results_append(
    f"gf_ape_smear 1 step is unitary (max |U^dag U - 1| < 1e-13) = {is_ape_unitary}"
)
assert is_ape_unitary, max_ut_dev_ape

# (b) the rescale is applied after every step: four `num_step = 1` iterations
# keep the global `<N^2>` at 3, and a pure scalar rescale of a low field
# (`alpha = 0`, `num_step = 1`) multiplies it by `sqrt(3 / <N^2>)`.
gf_rescale_walk = gf.copy()
rescale_avg_list = []
for _ in range(4):
    gf_rescale_walk = q.gf_ape_smear_no_proj_rescaling(
        gf_rescale_walk, rescale_alpha, 1
    )
    rescale_avg_list.append(get_glb_avg_link_norm_sqr(gf_rescale_walk))
q.json_results_append(
    "gf_ape_smear_no_proj_rescaling 4 x num_step=1 global avg <N^2>",
    np.array(rescale_avg_list),
    1e-10,
)
is_avg_kept = bool(max(abs(v - 3.0) for v in rescale_avg_list) < 1e-13)
q.json_results_append(
    f"gf_ape_smear_no_proj_rescaling per-step avg <N^2> = 3 to 1e-13 = {is_avg_kept}"
)
assert is_avg_kept, rescale_avg_list
gf_rescale_low = gf.copy()
gf_rescale_low.unitarize()
gf_rescale_low *= 0.5
arr_rescale_low = np.asarray(gf_rescale_low).copy()
gf_rescale_restore = q.gf_ape_smear_no_proj_rescaling(gf_rescale_low, 0.0, 1)
restore_dev = get_max_rel_dev(
    np.asarray(gf_rescale_restore).copy(), 2.0 * arr_rescale_low
)
is_restore_scaled = bool(restore_dev < 1e-12)
q.json_results_append(
    "gf_ape_smear_no_proj_rescaling alpha=0 on low field is 2 x field"
    f" (max rel dev < 1e-12) = {is_restore_scaled}"
)
assert is_restore_scaled, restore_dev
avg_restore = get_glb_avg_link_norm_sqr(gf_rescale_restore)
q.json_results_append(
    "gf_ape_smear_no_proj_rescaling low field restored global avg <N^2>",
    avg_restore,
    1e-10,
)
assert abs(avg_restore - 3.0) <= 1e-13, avg_restore

# (c) agreement and edges: `num_step = 0` / `-1` return bit-identical copies,
# `alpha = 0` on a unit field is exactly the identity (for both this function
# and `gf_ape_smear`), `num_step = 3` is three `num_step = 1` calls, and the
# input field is never modified.
is_zero_copy = bool(
    np.array_equal(
        np.asarray(q.gf_ape_smear_no_proj_rescaling(gf, rescale_alpha, 0)),
        gf_rescale_sig,
    )
)
q.json_results_append(
    f"gf_ape_smear_no_proj_rescaling num_step=0 bitwise copy = {is_zero_copy}"
)
assert is_zero_copy
is_neg_copy = bool(
    np.array_equal(
        np.asarray(q.gf_ape_smear_no_proj_rescaling(gf, rescale_alpha, -1)),
        gf_rescale_sig,
    )
)
q.json_results_append(
    f"gf_ape_smear_no_proj_rescaling num_step=-1 bitwise copy = {is_neg_copy}"
)
assert is_neg_copy
gf_rescale_unit = q.GaugeField(geo)
gf_rescale_unit.set_unit()
arr_rescale_unit = np.asarray(gf_rescale_unit).copy()
is_rescale_unit_identity = bool(
    np.array_equal(
        np.asarray(q.gf_ape_smear_no_proj_rescaling(gf_rescale_unit, 0.0, 1)),
        arr_rescale_unit,
    )
)
q.json_results_append(
    "gf_ape_smear_no_proj_rescaling alpha=0 on unit field bitwise identity"
    f" = {is_rescale_unit_identity}"
)
assert is_rescale_unit_identity
is_ape_unit_identity = bool(
    np.array_equal(
        np.asarray(q.gf_ape_smear(gf_rescale_unit, 0.0, 1)), arr_rescale_unit
    )
)
q.json_results_append(
    f"gf_ape_smear alpha=0 on unit field bitwise identity = {is_ape_unit_identity}"
)
assert is_ape_unit_identity
gf_rescale_3 = q.gf_ape_smear_no_proj_rescaling(gf, rescale_alpha, 3)
gf_rescale_loop = gf.copy()
for _ in range(3):
    gf_rescale_loop = q.gf_ape_smear_no_proj_rescaling(
        gf_rescale_loop, rescale_alpha, 1
    )
rescale_loop_dev = get_max_rel_dev(
    np.asarray(gf_rescale_3).copy(), np.asarray(gf_rescale_loop).copy()
)
is_rescale_loop_equal = bool(rescale_loop_dev < 1e-10)
q.json_results_append(
    "gf_ape_smear_no_proj_rescaling num_step=3 vs three num_step=1 calls"
    f" rel dev < 1e-10 = {is_rescale_loop_equal}"
)
assert is_rescale_loop_equal, rescale_loop_dev
q.json_results_append(
    "gf_ape_smear_no_proj_rescaling 3 steps plaq", gf_rescale_3.plaq(), 1e-10
)
assert abs(gf_rescale_3.plaq() - 0.6407780936466481) <= 1e-10
is_rescale_input_kept = bool(np.array_equal(np.asarray(gf), gf_rescale_sig))
q.json_results_append(
    "gf_ape_smear_no_proj_rescaling input field never modified"
    f" = {is_rescale_input_kept}"
)
assert is_rescale_input_kept

# ---- gf_stout_smear ----

gf_sig = q.get_data_sig_arr(gf, q.RngState(), 3)
for num_step in [1, 3, 5]:
    gf_test = q.gf_stout_smear(gf, 0.1, num_step)
    q.json_results_append(
        f"gf_stout_smear steps={num_step} plaq", gf_test.plaq(), 1e-10
    )
assert np.all(q.get_data_sig_arr(gf, q.RngState(), 3) == gf_sig)

# ---- gf_stout_smear_block ----

gf_bs = q.gf_stout_smear_block(gf, q.Coordinate(), 0.1)
q.json_results_append("gf_block_stout_smear plaq", gf_bs.plaq(), 1e-10)

gf_bs2 = q.gf_stout_smear_block(gf, q.Coordinate([2, 2, 2, 2]), 0.1)
q.json_results_append("gf_block_stout_smear block=[2,2,2,2] plaq", gf_bs2.plaq(), 1e-10)

# ---- gf_stout_smear_local ----

gf_ls = q.gf_stout_smear_local(gf, q.Coordinate(), 0.1)
q.json_results_append("gf_local_stout_smear plaq", gf_ls.plaq(), 1e-10)

assert np.all(q.get_data_sig_arr(gf, q.RngState(), 3) == gf_sig)

# ---- field_smear (Gaussian smearing) ----

f1 = q.FieldComplexD(geo, 3)
f1.set_rand(rs.split("f1"))
for radius in [1.0, 2.0, 3.0]:
    sf = q.field_smear(f1, radius)
    q.json_results_append(
        f"smear_field r={radius} sig",
        q.get_data_sig(sf, rs.split(f"sf{radius}")),
        1e-10,
    )

sf_sp = q.field_smear(f1, 2.0, spatial_only=True)
q.json_results_append(
    "smear_field spatial sig", q.get_data_sig(sf_sp, rs.split("sfsp")), 1e-10
)

# ---- field_smear_step (density field smearing) ----

f_real = q.FieldRealD(geo, 1)
f_real.set_rand(rs.split("f_real"))
for coef in [0.5, 0.8]:
    sf_step = q.field_smear_step(f_real, coef, 3)
    q.json_results_append(
        f"smear_field_step coef={coef} sig",
        q.get_data_sig(sf_step, rs.split(f"sfs{coef}")),
        1e-10,
    )

# ---- sphere_sum_field ----

for radius in [1.0, 2.0, 3.0]:
    ss = q.sphere_sum_field(f1, radius)
    q.json_results_append(
        f"sphere_sum_field r={radius} sig",
        q.get_data_sig(ss, rs.split(f"ss{radius}")),
        1e-10,
    )

ss_sp = q.sphere_sum_field(f1, 2.0, is_only_spatial=True)
q.json_results_append(
    "sphere_sum_field spatial sig", q.get_data_sig(ss_sp, rs.split("sssp")), 1e-10
)

# ---- prop_smear ----

gf_ape = q.gf_ape_smear_spatial(gf, 0.5, 3)
gf1 = q.mk_left_expanded_field(gf_ape)
prop = q.Prop(geo)
prop.set_rand(rs.split("prop-init"))
mom = q.CoordinateD([0.0, 0.1, -0.2, 0.2])
coef = 0.9375
step = 10

mode_name_dict = {0: "standard", 1: "qlat"}
for mode in [0, 1]:
    sp = q.prop_smear(prop, gf1, coef, step, mom, mode=mode_name_dict[mode])
    q.json_results_append(
        f"prop_smear mode={mode} sig",
        q.get_data_sig_arr(sp, rs.split(f"psm{mode}"), 3),
        1e-10,
    )

sp_no_mom = q.prop_smear(prop, gf1, coef, step)
q.json_results_append(
    "prop_smear no_mom sig", q.get_data_sig_arr(sp_no_mom, rs.split("psnm"), 3), 1e-10
)

sp_time = q.prop_smear(prop, gf1, coef, step, mom, smear_in_time_dir=True)
q.json_results_append(
    "prop_smear time_dir sig", q.get_data_sig_arr(sp_time, rs.split("psmt"), 3), 1e-10
)

# ---- prop_spatial_smear with Prop ----

sp_spatial = q.prop_spatial_smear(prop, gf_ape, coef, step, mom)
q.json_results_append(
    "prop_spatial_smear Prop sig",
    q.get_data_sig_arr(sp_spatial, rs.split("pss"), 3),
    1e-10,
)

sp_spatial_nomom = q.prop_spatial_smear(prop, gf_ape, coef, step)
q.json_results_append(
    "prop_spatial_smear Prop no_mom sig",
    q.get_data_sig_arr(sp_spatial_nomom, rs.split("pssnm"), 3),
    1e-10,
)

# ---- prop_spatial_smear with FermionField4d ----

ff = q.FermionField4d(geo)
ff.set_rand(rs.split("ff-init"))

sp_ff = q.prop_spatial_smear(ff, gf_ape, coef, step, mom)
q.json_results_append(
    "prop_spatial_smear FF4d sig", q.get_data_sig_arr(sp_ff, rs.split("pssf"), 3), 1e-10
)

sp_ff_nomom = q.prop_spatial_smear(ff, gf_ape, coef, step)
q.json_results_append(
    "prop_spatial_smear FF4d no_mom sig",
    q.get_data_sig_arr(sp_ff_nomom, rs.split("pssfnm"), 3),
    1e-10,
)

# ---- prop_spatial_smear with list of Prop ----

prop2 = q.Prop(geo)
prop2.set_rand(rs.split("prop2-init"))
sp_list_prop = q.prop_spatial_smear([prop, prop2], gf_ape, coef, step, mom)
q.json_results_append(f"prop_spatial_smear [Prop,Prop] len={len(sp_list_prop)}")
q.json_results_append(
    "prop_spatial_smear [Prop,Prop] sig",
    q.get_data_sig_arr(sp_list_prop[0], rs.split("pssp0"), 3),
    1e-10,
)

# ---- prop_spatial_smear with list of FermionField4d ----

ff2 = q.FermionField4d(geo)
ff2.set_rand(rs.split("ff2-init"))
sp_list_ff = q.prop_spatial_smear([ff, ff2], gf_ape, coef, step, mom)
q.json_results_append(f"prop_spatial_smear [FF4d,FF4d] len={len(sp_list_ff)}")
q.json_results_append(
    "prop_spatial_smear [FF4d,FF4d] sig",
    q.get_data_sig_arr(sp_list_ff[0], rs.split("pssf0"), 3),
    1e-10,
)

q.timer_display()
if q.is_test():
    q.check_log_json(__file__, check_eps=1e-0)
q.end_with_mpi()
q.displayln_info("CHECK: finished successfully.")
