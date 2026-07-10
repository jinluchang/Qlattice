#!/usr/bin/env python3

"""
Test time-direction truncation utilities:
- mk_field_truncated: generic field truncation
- mk_gf_truncated: symmetric gauge field truncation with divisor padding
- mk_gt_truncated: gauge transform truncation
- mk_gf_truncated_evolve: asymmetric truncation with unity-link padding
- mk_selected_points_truncated: point selection truncation
- gf_evolve_masked / gm_evolve_fg_pure_gauge_masked / run_hmc_evolve_pure_gauge_masked:
    masked HMC evolution that only updates links where mask is True
"""

import math

import qlat as q
import numpy as np

### ------

@q.timer
def mk_field_truncated(field, *, t_start, t_end):
    geo = field.geo
    total_site = geo.total_site
    t_size = total_site[3]
    t_size_trunc = (t_end - t_start) % t_size
    if t_size_trunc == 0:
        t_size_trunc = t_size
    assert t_size_trunc <= t_size
    t_size_trunc_valid = t_size_trunc
    size_node_t = geo.size_node[3]
    t_size_trunc = ((t_size_trunc + size_node_t - 1) // size_node_t) * size_node_t
    total_site_trunc = q.Coordinate(
        [total_site[0], total_site[1], total_site[2], t_size_trunc]
    )
    geo_trunc = q.Geometry(total_site_trunc)
    field_trunc = type(field)(geo_trunc, field.multiplicity)
    t_offset = t_start % t_size
    node_site_t = geo.node_site[3]
    n_shift = geo.size_node[3]
    shift = q.Coordinate([0, 0, 0, node_site_t])
    node_site_trunc = geo_trunc.node_site
    coor_node = geo.coor_node
    current_field = field.copy()
    current_xg_field = q.mk_xg_field(geo)
    field_trunc_arr = np.asarray(field_trunc)
    for i_shift in range(n_shift):
        current_field_arr = np.asarray(current_field)
        xg_arr = current_xg_field[:]
        for index in range(current_field.geo.local_volume):
            xg = xg_arr[index]
            t_orig = xg[3] % t_size
            t_local = (t_orig - t_offset) % t_size
            if t_local >= t_size_trunc_valid:
                continue
            xl = [xg[i] - coor_node[i] * node_site_trunc[i] for i in range(3)]
            xl3 = t_local - coor_node[3] * node_site_trunc[3]
            if xl3 < 0 or xl3 >= node_site_trunc[3]:
                continue
            index_trunc = xl[0] + node_site_trunc[0] * (
                xl[1] + node_site_trunc[1] * (xl[2] + node_site_trunc[2] * xl3)
            )
            field_trunc_arr[index_trunc] = current_field_arr[index]
        if i_shift < n_shift - 1:
            current_field = current_field.shift(shift)
            current_xg_field = current_xg_field.shift(shift)
    return field_trunc

@q.timer
def mk_gf_truncated(gf, *, t_center, t_half, t_size_divisor=1):
    total_site = gf.geo.total_site
    t_size = total_site[3]
    t_size_trunc = 2 * t_half + 1
    r = t_size_trunc % t_size_divisor
    if r != 0:
        t_size_trunc += t_size_divisor - r
    assert t_size_trunc <= t_size
    t_start = (t_center - t_half) % t_size
    t_end = (t_start + t_size_trunc) % t_size
    gf_trunc = mk_field_truncated(gf, t_start=t_start, t_end=t_end)
    geo_trunc = gf_trunc.geo
    tslice_target_list = list(range(2 * t_half, t_size_trunc))
    gf_arr = np.asarray(gf_trunc)
    xg_arr = q.mk_xg_field(geo_trunc)[:]
    for index in range(geo_trunc.local_volume):
        xg = xg_arr[index]
        xg_t = xg[3]
        if xg_t in tslice_target_list:
            gf_arr[index, 3, :, :] = 0
    return gf_trunc, t_start, t_size_trunc

@q.timer
def mk_gt_truncated(gt, *, t_center, t_half, t_size_divisor=1):
    total_site = gt.geo.total_site
    t_size = total_site[3]
    t_size_trunc = 2 * t_half + 1
    r = t_size_trunc % t_size_divisor
    if r != 0:
        t_size_trunc += t_size_divisor - r
    assert t_size_trunc <= t_size
    t_start = (t_center - t_half) % t_size
    t_end = (t_start + t_size_trunc) % t_size
    gt_trunc = mk_field_truncated(gt, t_start=t_start, t_end=t_end)
    return gt_trunc, t_start, t_size_trunc

@q.timer(is_verbose=True)
def mk_gf_truncated_evolve(gf, *, t_center, t_left, t_right, t_pad):
    total_site = gf.geo.total_site
    t_size = total_site[3]
    t_size_trunc_valid = t_left + t_right + 1
    t_size_trunc = t_size_trunc_valid + t_pad
    size_node_t = gf.geo.size_node[3]
    t_size_trunc = ((t_size_trunc + size_node_t - 1) // size_node_t) * size_node_t
    assert t_size_trunc <= t_size
    t_start = (t_center - t_left) % t_size
    t_end = (t_start + t_size_trunc) % t_size
    gf_trunc = mk_field_truncated(gf, t_start=t_start, t_end=t_end)
    geo_trunc = gf_trunc.geo
    gf_arr = np.asarray(gf_trunc)
    assert gf_arr.shape == (geo_trunc.local_volume, 4, 3, 3)
    xg_arr = q.mk_xg_field(geo_trunc)[:]
    eye3 = np.eye(3, dtype=gf_arr.dtype)
    mask = np.zeros((geo_trunc.local_volume, 4), dtype=bool)
    mask[xg_arr[:, 3] >= t_size_trunc_valid, :] = True
    mask[xg_arr[:, 3] == t_size_trunc_valid - 1, 3] = True
    gf_arr[mask] = eye3
    return gf_trunc, t_start, t_size_trunc

@q.timer
def mk_selected_points_truncated(sp, *, idx_start, idx_end):
    n_keep = idx_end - idx_start
    total_site = sp.psel.total_site
    psel_sub = q.PointsSelection(total_site, n_keep)
    np.asarray(psel_sub)[:] = np.asarray(sp.psel)[idx_start:idx_end]
    sp_trunc = type(sp)(psel_sub)
    sp_trunc @= sp
    return sp_trunc

### ------
### Masked HMC evolution functions
### ------

@q.timer
def gf_evolve_masked(gf, gm, dt, mask):
    gf_saved = q.GaugeField(gf.geo)
    gf_saved @= gf
    q.gf_evolve(gf, gm, dt)
    gf_arr = np.asarray(gf)
    gf_saved_arr = np.asarray(gf_saved)
    gf_arr[~mask] = gf_saved_arr[~mask]

@q.timer
def gf_unitarize_masked(gf, mask):
    gf_saved = q.GaugeField(gf.geo)
    gf_saved @= gf
    gf.unitarize()
    gf_arr = np.asarray(gf)
    gf_saved_arr = np.asarray(gf_saved)
    gf_arr[~mask] = gf_saved_arr[~mask]

@q.timer
def gm_evolve_fg_pure_gauge_masked(gm, gf_init, ga, fg_dt, dt, mask):
    geo = gf_init.geo
    gf = q.GaugeField(geo)
    gf @= gf_init
    gm_force = q.GaugeMomentum(geo)
    q.set_gm_force(gm_force, gf, ga)
    gf_evolve_masked(gf, gm_force, fg_dt, mask)
    q.set_gm_force(gm_force, gf, ga)
    gm_force *= dt
    gm_arr = np.asarray(gm)
    gm_force_arr = np.asarray(gm_force)
    gm_arr[mask] += gm_force_arr[mask]

@q.timer(is_timer_fork=True)
def run_hmc_evolve_pure_gauge_masked(gm, gf, ga, n_step, mask, md_time=1.0):
    energy = q.gm_hamilton_node(gm) + q.gf_hamilton_node(gf, ga)
    dt = md_time / n_step
    lam = 0.5 * (1.0 - 1.0 / math.sqrt(3.0))
    theta = (2.0 - math.sqrt(3.0)) / 48.0
    ttheta = theta * dt * dt * dt
    gf_evolve_masked(gf, gm, lam * dt, mask)
    for i in range(n_step):
        gm_evolve_fg_pure_gauge_masked(gm, gf, ga, 4.0 * ttheta / dt, 0.5 * dt, mask)
        gf_evolve_masked(gf, gm, (1.0 - 2.0 * lam) * dt, mask)
        gm_evolve_fg_pure_gauge_masked(gm, gf, ga, 4.0 * ttheta / dt, 0.5 * dt, mask)
        if i < n_step - 1:
            gf_evolve_masked(gf, gm, 2.0 * lam * dt, mask)
        else:
            gf_evolve_masked(gf, gm, lam * dt, mask)
    gf_unitarize_masked(gf, mask)
    delta_h = q.gm_hamilton_node(gm) + q.gf_hamilton_node(gf, ga) - energy
    delta_h = q.glb_sum_double(delta_h)
    return delta_h

@q.timer(is_timer_fork=True)
def run_hmc_pure_gauge_masked(
    gf,
    ga,
    traj,
    rs,
    mask,
    *,
    n_step=6,
    md_time=1.0,
):
    fname = q.get_fname()
    rs = rs.split(f"{traj}")
    geo = gf.geo
    gf_orig = q.GaugeField(geo)
    gf_orig @= gf
    gf0 = q.GaugeField(geo)
    gf0 @= gf
    gm = q.GaugeMomentum(geo)
    gm.set_rand(rs.split("set_rand_gauge_momentum"), 1.0)
    gm_arr = np.asarray(gm)
    gm_arr[~mask] = 0.0
    delta_h = run_hmc_evolve_pure_gauge_masked(gm, gf0, ga, n_step, mask, md_time)
    q.metropolis_accept(delta_h, traj, rs.split("metropolis_accept"))
    gf0_arr = np.asarray(gf0)
    gf_orig_arr = np.asarray(gf_orig)
    assert (gf0_arr[~mask] == gf_orig_arr[~mask]).all(), (
        "non-masked gf links should be unchanged"
    )
    q.displayln_info(f"{fname}: update gf (traj={traj})")
    gf @= gf0
    return delta_h

### ------

q.begin_with_mpi()

total_site = q.Coordinate([4, 4, 4, 16])
geo = q.Geometry(total_site)
q.displayln_info(f"geo.show()={geo.show()}")
rs = q.RngState("seed")

gf = q.GaugeField(geo)
gf.set_rand(rs.split("gf-init"), 0.3, 1)

plaq = gf.plaq()
q.json_results_append("plaq", plaq, 1e-12)

### ------

# Test 1: basic truncation, no padding
t_center = 8
t_left = 3
t_right = 4
t_pad = 0
gf_trunc, t_start, t_size_trunc = mk_gf_truncated_evolve(
    gf, t_center=t_center, t_left=t_left, t_right=t_right, t_pad=t_pad
)
q.json_results_append("test1 t_start", t_start)
q.json_results_append("test1 t_size_trunc", t_size_trunc)
assert t_size_trunc == 8  # 3 + 4 + 1
assert t_start == (t_center - t_left) % 16

# Test 2: truncation with padding
t_pad = 4
gf_trunc, t_start, t_size_trunc = mk_gf_truncated_evolve(
    gf, t_center=t_center, t_left=t_left, t_right=t_right, t_pad=t_pad
)
q.json_results_append("test2 t_start", t_start)
q.json_results_append("test2 t_size_trunc", t_size_trunc)
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
q.json_results_append("test3 pad_ok", n_pad_ok)
q.json_results_append("test3 pad_total", n_pad_total)
q.json_results_append("test3 valid_not_identity", n_valid_ok)
q.json_results_append("test3 valid_total", n_valid_total)
assert n_pad_ok == n_pad_total, (
    f"padded sites should be unity: {n_pad_ok}/{n_pad_total}"
)
if n_valid_total > 0:
    assert n_valid_ok > 0, "valid sites should not all be identity"

# Test 4: verify valid region data matches original gauge field via plaquette
plaq_trunc = gf_trunc.plaq()
q.json_results_append("test4 plaq_trunc", plaq_trunc, 1e-12)
assert 0.0 < plaq_trunc < 1.0, f"truncated plaq should be physical: {plaq_trunc}"

# Test 5: wrap-around at boundary
t_center_wrap = 2
t_left = 3
t_right = 2
t_pad = 0
gf_trunc, t_start, t_size_trunc = mk_gf_truncated_evolve(
    gf, t_center=t_center_wrap, t_left=t_left, t_right=t_right, t_pad=t_pad
)
q.json_results_append("test5 t_start", t_start)
q.json_results_append("test5 t_size_trunc", t_size_trunc)
assert t_start == 15  # (2 - 3) % 16 = 15
assert t_size_trunc == 6  # 3 + 2 + 1

# Test 6: verify unity padding with wrap-around
t_pad = 2
gf_trunc, t_start, t_size_trunc = mk_gf_truncated_evolve(
    gf, t_center=t_center_wrap, t_left=t_left, t_right=t_right, t_pad=t_pad
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
assert n_pad_ok == n_pad_total, (
    f"padded sites should be unity (wrap): {n_pad_ok}/{n_pad_total}"
)
q.json_results_append("test6 t_start", t_start)
q.json_results_append("test6 t_size_trunc", t_size_trunc)
q.json_results_append("test6 pad_ok", n_pad_ok)
q.json_results_append("test6 pad_total", n_pad_total)

# Test 7: mk_gf_truncated symmetric, no divisor
t_center = 8
t_half = 3
gf_trunc, t_start, t_size_trunc = mk_gf_truncated(gf, t_center=t_center, t_half=t_half)
q.json_results_append("test7 t_start", t_start)
q.json_results_append("test7 t_size_trunc", t_size_trunc)
assert t_size_trunc == 7  # 2 * 3 + 1
assert t_start == (t_center - t_half) % 16
plaq_trunc = gf_trunc.plaq()
q.json_results_append("test7 plaq_trunc", plaq_trunc, 1e-12)
assert 0.0 < plaq_trunc < 1.0

# Test 8: mk_gf_truncated with divisor padding
t_size_divisor = 4
gf_trunc, t_start, t_size_trunc = mk_gf_truncated(
    gf, t_center=t_center, t_half=t_half, t_size_divisor=t_size_divisor
)
q.json_results_append("test8 t_start", t_start)
q.json_results_append("test8 t_size_trunc", t_size_trunc)
assert t_size_trunc == 8  # ceil(7 / 4) * 4

# Test 9: mk_gf_truncated padding slices have zero temporal links
geo_trunc = gf_trunc.geo
gf_arr = np.asarray(gf_trunc)
xg_arr = q.mk_xg_field(geo_trunc)[:]
n_pad_ok = 0
n_pad_total = 0
for index in range(geo_trunc.local_volume):
    xg = xg_arr[index]
    if xg[3] >= 2 * t_half:
        n_pad_total += 1
        if np.allclose(gf_arr[index, 3, :, :], 0):
            n_pad_ok += 1
q.json_results_append("test9 pad_zero_temporal_ok", n_pad_ok)
q.json_results_append("test9 pad_zero_temporal_total", n_pad_total)
assert n_pad_ok == n_pad_total, (
    f"padded slices should have zero temporal links: {n_pad_ok}/{n_pad_total}"
)

# Test 10: mk_gt_truncated
gt = q.GaugeTransform(geo)
gt.set_rand(rs.split("gt-init"))
gt_trunc, t_start, t_size_trunc = mk_gt_truncated(
    gt, t_center=t_center, t_half=t_half, t_size_divisor=t_size_divisor
)
q.json_results_append("test10 t_start", t_start)
q.json_results_append("test10 t_size_trunc", t_size_trunc)
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
q.json_results_append("test11 data_match", n_match)
q.json_results_append("test11 data_total", n_total)

# Test 12: mk_selected_points_truncated
prop = q.Prop(geo)
prop.set_rand(rs.split("prop-init"))
ps = prop.glb_sum_tslice()
q.json_results_append("test12 ps qnorm", ps.qnorm(), 1e-12)
n_total_ps = len(ps)
ps_trunc = mk_selected_points_truncated(ps, idx_start=2, idx_end=6)
q.json_results_append("test12 ps_trunc qnorm", ps_trunc.qnorm(), 1e-12)
assert len(ps_trunc) == 4
q.json_results_append("test12 n_total", n_total_ps)
q.json_results_append("test12 n_trunc", len(ps_trunc))

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
gf_evolve_masked(gf13, gm13, 0.1, mask13)
gf13_arr = np.asarray(gf13)
gf13_copy_arr = np.asarray(gf13_copy)
assert np.allclose(gf13_arr[~mask13], gf13_copy_arr[~mask13]), (
    "non-masked (frozen) links should be unchanged"
)
assert not np.allclose(gf13_arr[mask13], gf13_copy_arr[mask13]), (
    "masked links should have changed"
)
q.json_results_append("test13 gf_evolve_masked", 1)

# Test 14: run_hmc_evolve_pure_gauge_masked — freeze temporal links
ga = q.GaugeAction(5.5, 0.0)
rs14 = rs.split("test14")
gf14 = q.GaugeField(geo)
gf14 @= gf
gm14 = q.GaugeMomentum(geo)
gm14.set_rand(rs14.split("set_rand_gauge_momentum"), 1.0)
mask14 = np.ones((geo.local_volume, 4), dtype=bool)
mask14[:, 3] = False
gf14_copy = q.GaugeField(geo)
gf14_copy @= gf14
delta_h14 = run_hmc_evolve_pure_gauge_masked(gm14, gf14, ga, n_step=6, mask=mask14)
q.json_results_append("test14 delta_h", delta_h14, 1e-12)
gf14_arr = np.asarray(gf14)
gf14_copy_arr = np.asarray(gf14_copy)
assert np.allclose(gf14_arr[:, 3], gf14_copy_arr[:, 3]), (
    "frozen temporal links should be unchanged"
)
assert abs(delta_h14) > 1e-12, "delta_h should be non-zero"

# Test 15: mask all-True — full evolution
rs15 = rs.split("test15")
gm15 = q.GaugeMomentum(geo)
gm15.set_rand(rs15.split("set_rand_gauge_momentum"), 1.0)
gf15 = q.GaugeField(geo)
gf15 @= gf
mask15 = np.ones((geo.local_volume, 4), dtype=bool)
delta_h15 = run_hmc_evolve_pure_gauge_masked(gm15, gf15, ga, n_step=6, mask=mask15)
q.json_results_append("test15 delta_h", delta_h15, 1e-12)
assert abs(delta_h15) > 1e-12, "all-True mask should give non-zero delta_h"

# Test 16: mask all-False — no evolution
rs16 = rs.split("test16")
gm16 = q.GaugeMomentum(geo)
gm16.set_rand(rs16.split("set_rand_gauge_momentum"), 1.0)
gf16 = q.GaugeField(geo)
gf16 @= gf
gf16_copy = q.GaugeField(geo)
gf16_copy @= gf16
mask16 = np.zeros((geo.local_volume, 4), dtype=bool)
delta_h16 = run_hmc_evolve_pure_gauge_masked(gm16, gf16, ga, n_step=6, mask=mask16)
q.json_results_append("test16 delta_h", delta_h16, 1e-12)
assert abs(delta_h16) < 1e-12, f"all-False mask should give zero delta_h: {delta_h16}"
gf16_arr = np.asarray(gf16)
gf16_copy_arr = np.asarray(gf16_copy)
assert np.allclose(gf16_arr, gf16_copy_arr), (
    "all-False mask should leave field unchanged"
)

# Test 17: masked evolution on truncated geometry
rs17 = rs.split("test17")
t_center = 8
t_left = 3
t_right = 4
t_pad = 4
gf_trunc17, t_start, t_size_trunc = mk_gf_truncated_evolve(
    gf, t_center=t_center, t_left=t_left, t_right=t_right, t_pad=t_pad
)
geo_trunc17 = gf_trunc17.geo
t_size_trunc_valid17 = t_left + t_right + 1
xg_arr17 = q.mk_xg_field(geo_trunc17)[:]
mask17 = np.ones((geo_trunc17.local_volume, 4), dtype=bool)
mask17[xg_arr17[:, 3] >= t_size_trunc_valid17, :] = False
mask17[xg_arr17[:, 3] == t_size_trunc_valid17 - 1, 3] = False
gm17 = q.GaugeMomentum(geo_trunc17)
gm17.set_rand(rs17.split("set_rand_gauge_momentum"), 1.0)
gf17 = q.GaugeField(geo_trunc17)
gf17 @= gf_trunc17
gf17_copy = q.GaugeField(geo_trunc17)
gf17_copy @= gf17
delta_h17 = run_hmc_evolve_pure_gauge_masked(gm17, gf17, ga, n_step=6, mask=mask17)
q.json_results_append("test17 delta_h", delta_h17, 1e-12)
gf17_arr = np.asarray(gf17)
gf17_copy_arr = np.asarray(gf17_copy)
assert np.allclose(gf17_arr[~mask17], gf17_copy_arr[~mask17]), (
    "frozen links should be unchanged"
)
assert not np.allclose(gf17_arr[mask17], gf17_copy_arr[mask17]), (
    "active links should have changed"
)

# Test 18: run_hmc_pure_gauge_masked — freeze temporal links
rs18 = rs.split("test18")
gf18 = q.GaugeField(geo)
gf18 @= gf
mask18 = np.ones((geo.local_volume, 4), dtype=bool)
mask18[:, 3] = False  # freeze temporal links
gf18_copy = q.GaugeField(geo)
gf18_copy @= gf18
delta_h18 = run_hmc_pure_gauge_masked(gf18, ga, 18, rs18, mask18, n_step=6)
q.json_results_append("test18 delta_h", delta_h18, 1e-12)
gf18_arr = np.asarray(gf18)
gf18_copy_arr = np.asarray(gf18_copy)
assert np.allclose(gf18_arr[:, 3], gf18_copy_arr[:, 3]), (
    "frozen temporal links should be unchanged"
)
assert abs(delta_h18) > 1e-12, "delta_h should be non-zero"

# Test 19: run_hmc_pure_gauge_masked — mask all-False
rs19 = rs.split("test19")
gf19 = q.GaugeField(geo)
gf19 @= gf
mask19 = np.zeros((geo.local_volume, 4), dtype=bool)
gf19_copy = q.GaugeField(geo)
gf19_copy @= gf19
delta_h19 = run_hmc_pure_gauge_masked(gf19, ga, 19, rs19, mask19, n_step=6)
q.json_results_append("test19 delta_h", delta_h19, 1e-12)
assert abs(delta_h19) < 1e-12, f"all-False mask should give zero delta_h: {delta_h19}"
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
delta_h20 = run_hmc_pure_gauge_masked(gf20, ga, 20, rs20, mask20, n_step=6)
q.json_results_append("test20 delta_h", delta_h20, 1e-12)
assert abs(delta_h20) > 1e-12, "all-True mask should give non-zero delta_h"

# Test 21: delta_h decreases as n_step increases (same initial field, same RNG → same momentum)
rs21 = rs.split("test21")
mask21 = np.ones((geo.local_volume, 4), dtype=bool)
gf21_low = q.GaugeField(geo)
gf21_low @= gf
gf21_high = q.GaugeField(geo)
gf21_high @= gf
dh21_low = run_hmc_pure_gauge_masked(gf21_low, ga, 21, rs21, mask21, n_step=6)
dh21_high = run_hmc_pure_gauge_masked(gf21_high, ga, 21, rs21, mask21, n_step=24)
q.json_results_append("test21 dh_low_n6", dh21_low, 1e-12)
q.json_results_append("test21 dh_high_n24", dh21_high, 1e-12)
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
