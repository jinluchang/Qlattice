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

def get_max_rel_dev(arr, arr_ref):
    # max |arr - arr_ref| normalized by max |arr_ref|
    return float(np.max(np.abs(arr - arr_ref)) / np.max(np.abs(arr_ref)))

def replicate_no_proj_target_loop(gf, alpha, target):
    # Hand-rolled replication of the target-driven loop of
    # `gf_ape_smear_no_proj(gf, alpha, accumulated_alpha=target)`: per step
    # `M_k = global max Tr(U^dag U)`, `factor_k = 3 / M_k`,
    # `a_k = min(alpha * factor_k, target - acc)` and the raw unprojected
    # update `U' = (1 - a_k) U + (a_k / 6) staple`, with no rescaling.
    # Each step is one call with `accumulated_alpha = a_k`.
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
        acc += a_k
        a_list.append(a_k)
        m_list.append(m_sq)
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

# ---- gf_ape_smear_no_proj (v4: target accumulated alpha) ----
# Per step: `M_k = global max Tr(U^dag U)`, `factor_k = 3 / M_k`,
# `a_k = min(alpha * factor_k, target - acc)` and
# `U' = (1 - a_k) U + (a_k / 6) staple`, with NO rescaling of the field.
# `acc` sums `a_k` until it equals `target` exactly (the last step is
# truncated); `a_k > 1` is allowed, so the map can grow the field.

# `target = alpha * n` where the old call asked for n plain steps.
for steps in [1, 3]:
    gf_test = q.gf_ape_smear_no_proj(gf, 0.5, 0.5 * steps)
    q.json_results_append(
        f"gf_ape_smear_no_proj target={0.5 * steps} plaq", gf_test.plaq(), 1e-10
    )

no_proj_alpha = 0.5
# 2.5 needs three steps: [0.5, 1.06..., 0.94...], so a_k > 1 and the
# truncated last step are both exercised.
no_proj_target = 2.5
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

# (b) unitary input, in the only form that holds: a unitary input has
# `M_0 = 3`, so `factor_0 = 1` and the first step is the plain unprojected step
# of weight `a_0 = min(alpha, target)`; in particular `target = alpha` is one
# plain unprojected update of weight `alpha`.  For `n > 1` the bookkeeping
# relation `target = alpha * n` is NOT a trajectory equivalence: the loop
# recomputes `M_k` every step, so once the field shrinks `factor_k > 1` and the
# target is reached in fewer, larger steps (see the step count below).
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

# (e) informational: step count, accumulating `a_k` and the range of `M_k` the
# loop went through; `max Tr(U^dag U)` is no longer pinned to 3, and the target
# is reached in fewer steps than `target / alpha`.
q.json_results_append(
    f"gf_ape_smear_no_proj step count = {len(a_list_rep)}"
    f" (< target/alpha = {no_proj_target / no_proj_alpha})"
)
q.json_results_append("gf_ape_smear_no_proj max a_k", max(a_list_rep), 1e-10)
q.json_results_append(
    f"gf_ape_smear_no_proj a_k > 1 branch exercised = {bool(max(a_list_rep) > 1.0)}"
)
q.json_results_append("gf_ape_smear_no_proj M_k min", min(m_list_rep), 1e-10)
q.json_results_append("gf_ape_smear_no_proj M_k max", max(m_list_rep), 1e-10)
q.json_results_append(
    "gf_ape_smear_no_proj global avg <N^2> after target",
    get_glb_avg_link_norm_sqr(gf_target),
    1e-10,
)

# (f) `a_k > 1` is allowed: on a strongly suppressed input
# `factor_k = 3 / M_k` is huge, so the single truncated step uses
# `a_1 = target > 1` and the map grows the field (`1 - a_1 = -2`).
gf_suppressed = gf_no_proj.copy()
gf_suppressed *= 0.1
m_suppressed_before = get_glb_max_link_norm_sqr(gf_suppressed)
gf_suppressed_smeared = q.gf_ape_smear_no_proj(gf_suppressed, no_proj_alpha, 3.0)
m_suppressed_after = get_glb_max_link_norm_sqr(gf_suppressed_smeared)
q.json_results_append(
    "gf_ape_smear_no_proj large a_k grows max Tr(U^dag U)",
    float(m_suppressed_after / m_suppressed_before),
    1e-10,
)
assert m_suppressed_after > m_suppressed_before

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
