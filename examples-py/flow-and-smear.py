#!/usr/bin/env python3

import qlat as q
import numpy as np

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
    q.json_results_append(f"gf_flow_scale is_spatial={is_spatial} plaq", gf_test.plaq(), 1e-10)

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

# ---- gf_spatial_ape_smear ----

for steps in [1, 3, 5]:
    gf_test = q.gf_spatial_ape_smear(gf, 0.5, steps)
    q.json_results_append(f"gf_spatial_ape_smear steps={steps} plaq", gf_test.plaq(), 1e-10)

# ---- gf_hyp_smear ----

alpha_sets = [
    [0.75, 0.6, 0.3],
    [0.5, 0.5, 0.5],
    [0.9, 0.7, 0.4],
]
for a1, a2, a3 in alpha_sets:
    gf_test = q.gf_hyp_smear(gf, a1, a2, a3)
    q.json_results_append(f"gf_hyp_smear {a1},{a2},{a3} plaq", gf_test.plaq(), 1e-10)

# ---- gf_stout_smear ----

for method in ["force", "stout", "wilson-flow"]:
    gf_test = gf.copy()
    q.gf_stout_smear(gf_test, 0.1, 3, method=method)
    q.json_results_append(f"gf_stout_smear {method} plaq", gf_test.plaq(), 1e-10)

for num_step in [1, 3, 5]:
    gf_test = gf.copy()
    q.gf_stout_smear(gf_test, 0.1, num_step)
    q.json_results_append(f"gf_stout_smear steps={num_step} plaq", gf_test.plaq(), 1e-10)

# ---- gf_block_stout_smear ----

gf_bs = gf.copy()
q.gf_block_stout_smear(gf_bs, q.Coordinate(), 0.1)
q.json_results_append("gf_block_stout_smear plaq", gf_bs.plaq(), 1e-10)

gf_bs2 = gf.copy()
q.gf_block_stout_smear(gf_bs2, q.Coordinate([2, 2, 2, 2]), 0.1)
q.json_results_append("gf_block_stout_smear block=[2,2,2,2] plaq", gf_bs2.plaq(), 1e-10)

# ---- gf_local_stout_smear ----

gf_ls = gf.copy()
q.gf_local_stout_smear(gf_ls, q.Coordinate(), 0.1)
q.json_results_append("gf_local_stout_smear plaq", gf_ls.plaq(), 1e-10)

# ---- smear_field (Gaussian smearing) ----

f1 = q.FieldComplexD(geo, 3)
f1.set_rand(rs.split("f1"))
for radius in [1.0, 2.0, 3.0]:
    sf = q.smear_field(f1, radius)
    q.json_results_append(f"smear_field r={radius} sig", q.get_data_sig(sf, rs.split(f"sf{radius}")), 1e-10)

sf_sp = q.smear_field(f1, 2.0, is_only_spatial=True)
q.json_results_append("smear_field spatial sig", q.get_data_sig(sf_sp, rs.split("sfsp")), 1e-10)

# ---- smear_field_step (density field smearing) ----

f_real = q.FieldRealD(geo, 1)
f_real.set_rand(rs.split("f_real"))
for coef in [0.5, 0.8]:
    sf_step = q.smear_field_step(f_real, coef, 3)
    q.json_results_append(f"smear_field_step coef={coef} sig", q.get_data_sig(sf_step, rs.split(f"sfs{coef}")), 1e-10)

# ---- sphere_sum_field ----

for radius in [1.0, 2.0, 3.0]:
    ss = q.sphere_sum_field(f1, radius)
    q.json_results_append(f"sphere_sum_field r={radius} sig", q.get_data_sig(ss, rs.split(f"ss{radius}")), 1e-10)

ss_sp = q.sphere_sum_field(f1, 2.0, is_only_spatial=True)
q.json_results_append("sphere_sum_field spatial sig", q.get_data_sig(ss_sp, rs.split("sssp")), 1e-10)

# ---- prop_smear ----

gf_ape = q.gf_spatial_ape_smear(gf, 0.5, 3)
gf1 = q.mk_left_expanded_field(gf_ape)
prop = q.Prop(geo)
prop.set_rand(rs.split("prop-init"))
mom = q.CoordinateD([0.0, 0.1, -0.2, 0.2])
coef = 0.9375
step = 10

for mode in [0, 1]:
    sp = q.prop_smear(prop, gf1, coef, step, mom, mode_smear=mode)
    q.json_results_append(f"prop_smear mode={mode} sig", q.get_data_sig_arr(sp, rs.split(f"psm{mode}"), 3), 1e-10)

sp_no_mom = q.prop_smear(prop, gf1, coef, step)
q.json_results_append("prop_smear no_mom sig", q.get_data_sig_arr(sp_no_mom, rs.split("psnm"), 3), 1e-10)

sp_time = q.prop_smear(prop, gf1, coef, step, mom, smear_in_time_dir=True)
q.json_results_append("prop_smear time_dir sig", q.get_data_sig_arr(sp_time, rs.split("psmt"), 3), 1e-10)

# ---- prop_spatial_smear with Prop ----

sp_spatial = q.prop_spatial_smear(prop, gf_ape, coef, step, mom)
q.json_results_append("prop_spatial_smear Prop sig", q.get_data_sig_arr(sp_spatial, rs.split("pss"), 3), 1e-10)

sp_spatial_nomom = q.prop_spatial_smear(prop, gf_ape, coef, step)
q.json_results_append("prop_spatial_smear Prop no_mom sig", q.get_data_sig_arr(sp_spatial_nomom, rs.split("pssnm"), 3), 1e-10)

# ---- prop_spatial_smear with FermionField4d ----

ff = q.FermionField4d(geo)
ff.set_rand(rs.split("ff-init"))

sp_ff = q.prop_spatial_smear(ff, gf_ape, coef, step, mom)
q.json_results_append("prop_spatial_smear FF4d sig", q.get_data_sig_arr(sp_ff, rs.split("pssf"), 3), 1e-10)

sp_ff_nomom = q.prop_spatial_smear(ff, gf_ape, coef, step)
q.json_results_append("prop_spatial_smear FF4d no_mom sig", q.get_data_sig_arr(sp_ff_nomom, rs.split("pssfnm"), 3), 1e-10)

# ---- prop_spatial_smear with list of Prop ----

prop2 = q.Prop(geo)
prop2.set_rand(rs.split("prop2-init"))
sp_list_prop = q.prop_spatial_smear([prop, prop2], gf_ape, coef, step, mom)
q.json_results_append("prop_spatial_smear [Prop,Prop] len", len(sp_list_prop), 1e-10)
q.json_results_append("prop_spatial_smear [Prop,Prop] sig", q.get_data_sig_arr(sp_list_prop[0], rs.split("pssp0"), 3), 1e-10)

# ---- prop_spatial_smear with list of FermionField4d ----

ff2 = q.FermionField4d(geo)
ff2.set_rand(rs.split("ff2-init"))
sp_list_ff = q.prop_spatial_smear([ff, ff2], gf_ape, coef, step, mom)
q.json_results_append("prop_spatial_smear [FF4d,FF4d] len", len(sp_list_ff), 1e-10)
q.json_results_append("prop_spatial_smear [FF4d,FF4d] sig", q.get_data_sig_arr(sp_list_ff[0], rs.split("pssf0"), 3), 1e-10)

q.timer_display()
q.check_log_json(__file__, check_eps=1e-0)
q.end_with_mpi()
q.displayln_info("CHECK: finished successfully.")
