#!/usr/bin/env python3

"""
Test time-direction truncation utilities:
- mk_field_truncated: generic field truncation
- mk_gf_truncated: symmetric gauge field truncation with divisor padding
- mk_gt_truncated: gauge transform truncation
- mk_gf_truncated_evolve: asymmetric truncation with unity-link padding and masked HMC
- mk_selected_points_truncated: point selection truncation
- gf_evolve_masked / gm_evolve_fg_pure_gauge_masked / run_hmc_evolve_pure_gauge_masked:
    masked HMC evolution that only updates links where mask is True
"""

import qlat as q
import numpy as np

### ------

q.begin_with_mpi()

total_site = q.Coordinate([4, 4, 4, 16])
geo = q.Geometry(total_site)
q.displayln_info(f"geo.show()={geo.show()}")
rs = q.RngState("seed")

gf = q.GaugeField(geo)
gf.set_rand(rs.split("gf-init"), 0.3, 1)

plaq = gf.plaq()
q.json_results_append("plaq", plaq, 1e-8)

ga = q.GaugeAction(5.5, 0.0)
rs_evolve = rs.split("mk_gf_truncated_evolve")

### ------

# Test 1: basic truncation, no padding
t_center = 8
t_left = 3
t_right = 4
t_pad = 0
gf_trunc, t_start, t_size_trunc = q.mk_gf_truncated_evolve(
    gf,
    t_center=t_center,
    t_left=t_left,
    t_right=t_right,
    t_pad=t_pad,
    ga=ga,
    rs=rs_evolve,
    num_traj=1,
    n_step=6,
    md_time=1.0,
)
q.json_results_append(f"test1 t_start={t_start}")
q.json_results_append(f"test1 t_size_trunc={t_size_trunc}")
q.json_results_append(
    "test1 gf_trunc sig", q.get_data_sig(gf_trunc, rs.split("test1-sig")), 1e-8
)
assert t_size_trunc == 8  # 3 + 4 + 1
assert t_start == (t_center - t_left) % 16

# Test 2: truncation with padding
t_pad = 4
gf_trunc, t_start, t_size_trunc = q.mk_gf_truncated_evolve(
    gf,
    t_center=t_center,
    t_left=t_left,
    t_right=t_right,
    t_pad=t_pad,
    ga=ga,
    rs=rs_evolve,
    num_traj=1,
    n_step=6,
    md_time=1.0,
)
q.json_results_append(f"test2 t_start={t_start}")
q.json_results_append(f"test2 t_size_trunc={t_size_trunc}")
q.json_results_append(
    "test2 gf_trunc sig", q.get_data_sig(gf_trunc, rs.split("test2-sig")), 1e-8
)
assert t_size_trunc == 12  # 3 + 4 + 1 + 4

# Test 3: verify padded slices are unity and valid slices are not
t_size_trunc_valid = t_left + t_right + 1  # 8
geo_trunc = gf_trunc.geo
gf_arr = np.asarray(gf_trunc)
xg_arr = q.mk_xg_field(geo_trunc)[:]
eye3 = np.eye(3, dtype=gf_arr.dtype)
n_pad_ok = 0
n_pad_total = 0
n_valid_ok = 0
n_valid_total = 0
for index in range(geo_trunc.local_volume):
    xg = xg_arr[index]
    if xg[3] >= t_size_trunc_valid:
        n_pad_total += 1
        if np.allclose(gf_arr[index], eye3):
            n_pad_ok += 1
    else:
        n_valid_total += 1
        if not np.allclose(gf_arr[index], eye3):
            n_valid_ok += 1
q.json_results_append(f"test3 pad_ok={n_pad_ok}")
q.json_results_append(f"test3 pad_total={n_pad_total}")
q.json_results_append(f"test3 valid_not_identity={n_valid_ok}")
q.json_results_append(f"test3 valid_total={n_valid_total}")
assert n_pad_ok == 0, (
    f"padded sites should have evolved (not identity): {n_pad_ok}/{n_pad_total}"
)
if n_valid_total > 0:
    assert n_valid_ok > 0, "valid sites should not all be identity"

# Test 4: verify valid region data matches original gauge field via plaquette
plaq_trunc = gf_trunc.plaq()
q.json_results_append("test4 plaq_trunc", plaq_trunc, 1e-8)
assert 0.0 < plaq_trunc < 1.0, f"truncated plaq should be physical: {plaq_trunc}"

# Test 5: wrap-around at boundary
t_center_wrap = 2
t_left = 3
t_right = 2
t_pad = 0
gf_trunc, t_start, t_size_trunc = q.mk_gf_truncated_evolve(
    gf,
    t_center=t_center_wrap,
    t_left=t_left,
    t_right=t_right,
    t_pad=t_pad,
    ga=ga,
    rs=rs_evolve,
    num_traj=1,
    n_step=6,
    md_time=1.0,
)
q.json_results_append(f"test5 t_start={t_start}")
q.json_results_append(f"test5 t_size_trunc={t_size_trunc}")
q.json_results_append(
    "test5 gf_trunc sig", q.get_data_sig(gf_trunc, rs.split("test5-sig")), 1e-8
)
assert t_start == 15  # (2 - 3) % 16 = 15
assert t_size_trunc == 6  # 3 + 2 + 1

# Test 6: verify unity padding with wrap-around
t_pad = 2
gf_trunc, t_start, t_size_trunc = q.mk_gf_truncated_evolve(
    gf,
    t_center=t_center_wrap,
    t_left=t_left,
    t_right=t_right,
    t_pad=t_pad,
    ga=ga,
    rs=rs_evolve,
    num_traj=1,
    n_step=6,
    md_time=1.0,
)
assert t_size_trunc == 8  # 3 + 2 + 1 + 2, already divisible by 2
t_size_trunc_valid = t_left + t_right + 1  # 6
geo_trunc = gf_trunc.geo
gf_arr = np.asarray(gf_trunc)
xg_arr = q.mk_xg_field(geo_trunc)[:]
eye3 = np.eye(3, dtype=gf_arr.dtype)
n_pad_ok = 0
n_pad_total = 0
for index in range(geo_trunc.local_volume):
    xg = xg_arr[index]
    if xg[3] >= t_size_trunc_valid:
        n_pad_total += 1
        if np.allclose(gf_arr[index], eye3):
            n_pad_ok += 1
assert n_pad_ok == 0, (
    f"padded sites should have evolved (not identity, wrap): {n_pad_ok}/{n_pad_total}"
)
q.json_results_append(f"test6 t_start={t_start}")
q.json_results_append(f"test6 t_size_trunc={t_size_trunc}")
q.json_results_append(
    "test6 gf_trunc sig", q.get_data_sig(gf_trunc, rs.split("test6-sig")), 1e-8
)
q.json_results_append(f"test6 pad_ok={n_pad_ok}")
q.json_results_append(f"test6 pad_total={n_pad_total}")

# Test 7: mk_gf_truncated symmetric, divisor=2 for MPI time-rank compatibility
t_center = 8
t_half = 3
gf_trunc, t_start, t_size_trunc = q.mk_gf_truncated(
    gf, t_center=t_center, t_half=t_half, t_size_divisor=2
)
q.json_results_append(f"test7 t_start={t_start}")
q.json_results_append(f"test7 t_size_trunc={t_size_trunc}")
q.json_results_append(
    "test7 gf_trunc sig", q.get_data_sig(gf_trunc, rs.split("test7-sig")), 1e-8
)
assert t_size_trunc == 8  # ceil((2 * 3 + 1) / 2) * 2
assert t_start == (t_center - t_half) % 16
plaq_trunc = gf_trunc.plaq()
q.json_results_append("test7 plaq_trunc", plaq_trunc, 1e-8)
assert 0.0 < plaq_trunc < 1.0

# Test 8: mk_gf_truncated with divisor padding
t_size_divisor = 4
gf_trunc, t_start, t_size_trunc = q.mk_gf_truncated(
    gf, t_center=t_center, t_half=t_half, t_size_divisor=t_size_divisor
)
q.json_results_append(f"test8 t_start={t_start}")
q.json_results_append(f"test8 t_size_trunc={t_size_trunc}")
q.json_results_append(
    "test8 gf_trunc sig", q.get_data_sig(gf_trunc, rs.split("test8-sig")), 1e-8
)
assert t_size_trunc == 8  # ceil(7 / 4) * 4

# Test 9: mk_gf_truncated padding slices have zero links, boundary has zero temporal link
geo_trunc = gf_trunc.geo
gf_arr = np.asarray(gf_trunc)
xg_arr = q.mk_xg_field(geo_trunc)[:]
zero3 = np.zeros((3, 3), dtype=gf_arr.dtype)
n_pad_ok = 0
n_pad_total = 0
n_boundary_ok = 0
n_boundary_total = 0
for index in range(geo_trunc.local_volume):
    xg = xg_arr[index]
    if xg[3] >= 2 * t_half + 1:
        # Fully padded slice: all links zero
        n_pad_total += 1
        if np.allclose(gf_arr[index], 0):
            n_pad_ok += 1
    elif xg[3] == 2 * t_half:
        # Boundary slice: only temporal link zero
        n_boundary_total += 1
        if np.allclose(gf_arr[index, 3, :, :], 0):
            n_boundary_ok += 1
q.json_results_append(f"test9 pad_zero_ok={n_pad_ok}")
q.json_results_append(f"test9 pad_zero_total={n_pad_total}")
q.json_results_append(f"test9 boundary_zero_ok={n_boundary_ok}")
q.json_results_append(f"test9 boundary_zero_total={n_boundary_total}")
assert n_pad_ok == n_pad_total, (
    f"padded slices should have zero links: {n_pad_ok}/{n_pad_total}"
)
assert n_boundary_ok == n_boundary_total, (
    f"boundary slice temporal links should be zero: {n_boundary_ok}/{n_boundary_total}"
)

# Test 10: mk_gt_truncated
gt = q.GaugeTransform(geo)
gt.set_rand(rs.split("gt-init"))
gt_trunc, t_start, t_size_trunc = q.mk_gt_truncated(
    gt, t_center=t_center, t_half=t_half, t_size_divisor=t_size_divisor
)
q.json_results_append(f"test10 t_start={t_start}")
q.json_results_append(f"test10 t_size_trunc={t_size_trunc}")
assert t_size_trunc == 8
assert gt_trunc.geo.total_site[3] == t_size_trunc

# Test 11: mk_gt_truncated data matches original for valid region
gt_trunc_arr = np.asarray(gt_trunc)
gt_arr = np.asarray(gt)
xg_trunc_arr = q.mk_xg_field(gt_trunc.geo)[:]
xg_full_arr = q.mk_xg_field(geo)[:]
n_match = 0
n_total = 0
for index in range(gt_trunc.geo.local_volume):
    xg = xg_trunc_arr[index]
    t_local = xg[3]
    if t_local >= 2 * t_half:
        continue
    t_orig = (t_local + (t_center - t_half)) % 16
    for full_idx in range(geo.local_volume):
        xg_full = xg_full_arr[full_idx]
        if xg_full[3] == t_orig and all(xg_full[i] == xg[i] for i in range(3)):
            n_total += 1
            if np.allclose(gt_trunc_arr[index], gt_arr[full_idx]):
                n_match += 1
            break
q.json_results_append(f"test11 data_match={n_match}")
q.json_results_append(f"test11 data_total={n_total}")

# Test 12: mk_selected_points_truncated
prop = q.Prop(geo)
prop.set_rand(rs.split("prop-init"))
ps = prop.glb_sum_tslice()
q.json_results_append("test12 ps qnorm", ps.qnorm(), 1e-8)
n_total_ps = len(ps)
ps_trunc = q.mk_selected_points_truncated(ps, idx_start=2, idx_end=6)
q.json_results_append("test12 ps_trunc qnorm", ps_trunc.qnorm(), 1e-8)
assert len(ps_trunc) == 4
q.json_results_append(f"test12 n_total={n_total_ps}")
q.json_results_append(f"test12 n_trunc={len(ps_trunc)}")

### ------

# Test 13: gf_evolve_masked — verify ~mask links unchanged
rs13 = rs.split("test13")
gm13 = q.GaugeMomentum(geo)
gm13.set_rand(rs13.split("set_rand_gauge_momentum"), 1.0)
gf13 = q.GaugeField(geo)
gf13 @= gf
mask13 = np.ones((geo.local_volume, 4), dtype=bool)
mask13[:, 3] = False  # freeze temporal links
gf13_copy = q.GaugeField(geo)
gf13_copy @= gf13
q.gf_evolve_masked(gf13, gm13, 0.1, mask13)
gf13_arr = np.asarray(gf13)
gf13_copy_arr = np.asarray(gf13_copy)
assert np.allclose(gf13_arr[~mask13], gf13_copy_arr[~mask13]), (
    "non-masked (frozen) links should be unchanged"
)
assert not np.allclose(gf13_arr[mask13], gf13_copy_arr[mask13]), (
    "masked links should have changed"
)
q.json_results_append(f"test13 gf_evolve_masked={1}")

# Test 14: run_hmc_evolve_pure_gauge_masked — freeze temporal links
rs14 = rs.split("test14")
gf14 = q.GaugeField(geo)
gf14 @= gf
gm14 = q.GaugeMomentum(geo)
gm14.set_rand(rs14.split("set_rand_gauge_momentum"), 1.0)
mask14 = np.ones((geo.local_volume, 4), dtype=bool)
mask14[:, 3] = False
gf14_copy = q.GaugeField(geo)
gf14_copy @= gf14
delta_h14 = q.run_hmc_evolve_pure_gauge_masked(gm14, gf14, ga, n_step=6, mask=mask14)
q.json_results_append("test14 delta_h", delta_h14, 1e-8)
gf14_arr = np.asarray(gf14)
gf14_copy_arr = np.asarray(gf14_copy)
assert np.allclose(gf14_arr[:, 3], gf14_copy_arr[:, 3]), (
    "frozen temporal links should be unchanged"
)
assert abs(delta_h14) > 1e-8, "delta_h should be non-zero"

# Test 15: mask all-True — full evolution
rs15 = rs.split("test15")
gm15 = q.GaugeMomentum(geo)
gm15.set_rand(rs15.split("set_rand_gauge_momentum"), 1.0)
gf15 = q.GaugeField(geo)
gf15 @= gf
mask15 = np.ones((geo.local_volume, 4), dtype=bool)
delta_h15 = q.run_hmc_evolve_pure_gauge_masked(gm15, gf15, ga, n_step=6, mask=mask15)
q.json_results_append("test15 delta_h", delta_h15, 1e-8)
assert abs(delta_h15) > 1e-8, "all-True mask should give non-zero delta_h"

# Test 16: mask all-False — no evolution
rs16 = rs.split("test16")
gm16 = q.GaugeMomentum(geo)
gm16.set_rand(rs16.split("set_rand_gauge_momentum"), 1.0)
gf16 = q.GaugeField(geo)
gf16 @= gf
gf16_copy = q.GaugeField(geo)
gf16_copy @= gf16
mask16 = np.zeros((geo.local_volume, 4), dtype=bool)
delta_h16 = q.run_hmc_evolve_pure_gauge_masked(gm16, gf16, ga, n_step=6, mask=mask16)
q.json_results_append("test16 delta_h", delta_h16, 1e-8)
assert abs(delta_h16) < 1e-8, f"all-False mask should give zero delta_h: {delta_h16}"
gf16_arr = np.asarray(gf16)
gf16_copy_arr = np.asarray(gf16_copy)
assert np.allclose(gf16_arr, gf16_copy_arr), (
    "all-False mask should leave field unchanged"
)

# Test 17: mk_gf_truncated_evolve with HMC
t_center = 8
t_left = 3
t_right = 4
t_pad = 4
rs17 = rs.split("test17")
gf_trunc17, t_start, t_size_trunc = q.mk_gf_truncated_evolve(
    gf,
    t_center=t_center,
    t_left=t_left,
    t_right=t_right,
    t_pad=t_pad,
    ga=ga,
    rs=rs17,
    num_traj=1,
    n_step=6,
    md_time=1.0,
)
geo_trunc17 = gf_trunc17.geo
t_size_trunc_valid17 = t_left + t_right + 1
gf_arr17 = np.asarray(gf_trunc17)
xg_arr17 = q.mk_xg_field(geo_trunc17)[:]
eye3 = np.eye(3, dtype=gf_arr17.dtype)
n_pad_ok = 0
n_pad_total = 0
for index in range(geo_trunc17.local_volume):
    xg = xg_arr17[index]
    if xg[3] >= t_size_trunc_valid17:
        n_pad_total += 1
        if np.allclose(gf_arr17[index], eye3):
            n_pad_ok += 1
q.json_results_append(f"test17 pad_ok={n_pad_ok}")
q.json_results_append(f"test17 pad_total={n_pad_total}")
q.json_results_append(
    "test17 gf_trunc sig", q.get_data_sig(gf_trunc17, rs.split("test17-sig")), 1e-8
)
assert n_pad_ok == 0, (
    f"padded sites should have evolved (not identity): {n_pad_ok}/{n_pad_total}"
)

# Test 18: run_hmc_pure_gauge_masked — freeze temporal links
rs18 = rs.split("test18")
gf18 = q.GaugeField(geo)
gf18 @= gf
mask18 = np.ones((geo.local_volume, 4), dtype=bool)
mask18[:, 3] = False  # freeze temporal links
gf18_copy = q.GaugeField(geo)
gf18_copy @= gf18
delta_h18 = q.run_hmc_pure_gauge_masked(gf18, ga, 18, rs18, mask18, n_step=6)
q.json_results_append("test18 delta_h", delta_h18, 1e-8)
gf18_arr = np.asarray(gf18)
gf18_copy_arr = np.asarray(gf18_copy)
assert np.allclose(gf18_arr[:, 3], gf18_copy_arr[:, 3]), (
    "frozen temporal links should be unchanged"
)
assert abs(delta_h18) > 1e-8, "delta_h should be non-zero"

# Test 19: run_hmc_pure_gauge_masked — mask all-False
rs19 = rs.split("test19")
gf19 = q.GaugeField(geo)
gf19 @= gf
mask19 = np.zeros((geo.local_volume, 4), dtype=bool)
gf19_copy = q.GaugeField(geo)
gf19_copy @= gf19
delta_h19 = q.run_hmc_pure_gauge_masked(gf19, ga, 19, rs19, mask19, n_step=6)
q.json_results_append("test19 delta_h", delta_h19, 1e-8)
assert abs(delta_h19) < 1e-8, f"all-False mask should give zero delta_h: {delta_h19}"
gf19_arr = np.asarray(gf19)
gf19_copy_arr = np.asarray(gf19_copy)
assert np.allclose(gf19_arr, gf19_copy_arr), (
    "all-False mask should leave field unchanged"
)

# Test 20: run_hmc_pure_gauge_masked — mask all-True
rs20 = rs.split("test20")
gf20 = q.GaugeField(geo)
gf20 @= gf
mask20 = np.ones((geo.local_volume, 4), dtype=bool)
delta_h20 = q.run_hmc_pure_gauge_masked(gf20, ga, 20, rs20, mask20, n_step=6)
q.json_results_append("test20 delta_h", delta_h20, 1e-8)
assert abs(delta_h20) > 1e-8, "all-True mask should give non-zero delta_h"

# Test 21: delta_h decreases as n_step increases (same initial field, same RNG → same momentum)
rs21 = rs.split("test21")
mask21 = np.ones((geo.local_volume, 4), dtype=bool)
gf21_low = q.GaugeField(geo)
gf21_low @= gf
gf21_high = q.GaugeField(geo)
gf21_high @= gf
dh21_low = q.run_hmc_pure_gauge_masked(gf21_low, ga, 21, rs21, mask21, n_step=6)
dh21_high = q.run_hmc_pure_gauge_masked(gf21_high, ga, 21, rs21, mask21, n_step=24)
q.json_results_append("test21 dh_low_n6", dh21_low, 1e-8)
q.json_results_append("test21 dh_high_n24", dh21_high, 1e-8)
assert abs(dh21_high) < abs(dh21_low), (
    f"higher n_step should give smaller |delta_h|: "
    f"|dh_low(n6)|={abs(dh21_low):.4g} vs |dh_high(n24)|={abs(dh21_high):.4g}"
)

### ------

q.timer_display()
if q.is_test():
    q.check_log_json(__file__)
q.end_with_mpi()
q.displayln_info("CHECK: finished successfully.")
