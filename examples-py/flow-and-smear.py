#!/usr/bin/env python3

import qlat as q
import numpy as np

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

def get_max_rel_dev(arr, arr_ref):
    # max |arr - arr_ref| normalized by max |arr_ref|
    return float(np.max(np.abs(arr - arr_ref)) / np.max(np.abs(arr_ref)))

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

# ---- gf_ape_smear_no_proj ----

for steps in [1, 3]:
    gf_test = q.gf_ape_smear_no_proj(gf, 0.5, steps)
    q.json_results_append(
        f"gf_simple_no_proj_smear steps={steps} plaq", gf_test.plaq(), 1e-10
    )

# ---- gf_ape_smear_no_proj: unconditional per-step rescale to <N^2> = 3 ----

gf_np_before = gf.copy()
arr_before = np.asarray(gf_np_before).copy()

# (a) after EVERY num_step=1 application the global average must be <N^2> = 3
gf_iter = gf_np_before.copy()
max_step_dev = 0.0
for k in range(1, 5):
    gf_iter = q.gf_ape_smear_no_proj(gf_iter, 0.5, 1)
    max_step_dev = max(max_step_dev, abs(get_glb_avg_link_norm_sqr(gf_iter) - 3.0))
q.json_results_append(
    "gf_ape_smear_no_proj per-step <N^2> max deviation from 3", max_step_dev
)
assert max_step_dev < 1e-13

# (b) one step is a uniform rescale: N^2_after / N^2_before == s^2 on every link.
# alpha = 0 keeps the raw update equal to the input, isolating the global
# rescale from the non-uniform staple update.  The input is multiplied by 0.5
# first (so s^2 = 4 and a missing rescale cannot hide: the random field itself
# already has <N^2> ~ 3, i.e. s ~ 1, which would make this check vacuous).
gf_low = gf_np_before.copy()
gf_low *= 0.5
avg_n2_low = get_glb_avg_link_norm_sqr(gf_low)
arr_low = np.asarray(gf_low).copy()
s_sq = 3.0 / avg_n2_low
s = float(np.sqrt(s_sq))
gf_scaled = q.gf_ape_smear_no_proj(gf_low, 0.0, 1)
n2_low = get_link_norm_sqr_arr(gf_low)
n2_scaled = get_link_norm_sqr_arr(gf_scaled)
max_ratio_dev = float(np.max(np.abs(n2_scaled / n2_low / s_sq - 1.0)))
q.json_results_append(
    "gf_ape_smear_no_proj single-step link-wise rescale uniformity", max_ratio_dev
)
assert max_ratio_dev < 1e-12

# (c) the rescale is real and non-negative: the output is s times the input
arr_scaled = np.asarray(gf_scaled).copy()
max_real_dev = get_max_rel_dev(arr_scaled, s * arr_low)
q.json_results_append(
    "gf_ape_smear_no_proj single-step rescale is a real scalar", max_real_dev
)
assert max_real_dev < 1e-12

# (d) num_step=0 returns the input field unchanged
gf_zero = q.gf_ape_smear_no_proj(gf_np_before, 0.5, 0)
is_zero_unchanged = bool(np.array_equal(np.asarray(gf_zero), arr_before))
q.json_results_append(
    f"gf_ape_smear_no_proj num_step=0 unchanged = {is_zero_unchanged}"
)
assert is_zero_unchanged

# (e) the rescale is applied after EVERY iteration, not only at the end:
# num_step=3 must equal three num_step=1 calls (the staple is quadratic, so
# scaling between iterations is observable)
gf_ns3 = q.gf_ape_smear_no_proj(gf_np_before, 0.5, 3)
gf_rep = gf_np_before.copy()
for _ in range(3):
    gf_rep = q.gf_ape_smear_no_proj(gf_rep, 0.5, 1)
max_per_step_dev = get_max_rel_dev(np.asarray(gf_ns3).copy(), np.asarray(gf_rep).copy())
q.json_results_append(
    "gf_ape_smear_no_proj num_step=3 vs three num_step=1", max_per_step_dev
)
assert max_per_step_dev < 1e-13
assert np.array_equal(np.asarray(gf_np_before), arr_before)

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
