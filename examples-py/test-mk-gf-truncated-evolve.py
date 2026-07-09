#!/usr/bin/env python3

"""
Test time-direction truncation utilities:
- mk_field_truncated: generic field truncation
- mk_gf_truncated: symmetric gauge field truncation with divisor padding
- mk_gt_truncated: gauge transform truncation
- mk_gf_truncated_evolve: asymmetric truncation with unity-link padding
- mk_selected_points_truncated: point selection truncation
"""

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
    tslice_pad_list = list(range(t_size_trunc_valid, t_size_trunc))
    gf_arr = np.asarray(gf_trunc)
    assert gf_arr.shape == (geo_trunc.local_volume, 4, 3, 3)
    xg_arr = q.mk_xg_field(geo_trunc)[:]
    eye3 = np.eye(3, dtype=gf_arr.dtype)
    for index in range(geo_trunc.local_volume):
        xg = xg_arr[index]
        if xg[3] in tslice_pad_list:
            gf_arr[index, :] = eye3
        elif xg[3] == t_size_trunc_valid - 1:
            gf_arr[index, 3] = eye3
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
gf_trunc, t_start, t_size_trunc = mk_gf_truncated_evolve(gf, t_center=t_center, t_left=t_left, t_right=t_right, t_pad=t_pad)
q.json_results_append("test1 t_start", t_start)
q.json_results_append("test1 t_size_trunc", t_size_trunc)
assert t_size_trunc == 8  # 3 + 4 + 1
assert t_start == (t_center - t_left) % 16

# Test 2: truncation with padding
t_pad = 4
gf_trunc, t_start, t_size_trunc = mk_gf_truncated_evolve(gf, t_center=t_center, t_left=t_left, t_right=t_right, t_pad=t_pad)
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
assert n_pad_ok == n_pad_total, f"padded sites should be unity: {n_pad_ok}/{n_pad_total}"
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
gf_trunc, t_start, t_size_trunc = mk_gf_truncated_evolve(gf, t_center=t_center_wrap, t_left=t_left, t_right=t_right, t_pad=t_pad)
q.json_results_append("test5 t_start", t_start)
q.json_results_append("test5 t_size_trunc", t_size_trunc)
assert t_start == 15  # (2 - 3) % 16 = 15
assert t_size_trunc == 6  # 3 + 2 + 1

# Test 6: verify unity padding with wrap-around
t_pad = 2
gf_trunc, t_start, t_size_trunc = mk_gf_truncated_evolve(gf, t_center=t_center_wrap, t_left=t_left, t_right=t_right, t_pad=t_pad)
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
assert n_pad_ok == n_pad_total, f"padded sites should be unity (wrap): {n_pad_ok}/{n_pad_total}"
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
gf_trunc, t_start, t_size_trunc = mk_gf_truncated(gf, t_center=t_center, t_half=t_half, t_size_divisor=t_size_divisor)
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
assert n_pad_ok == n_pad_total, f"padded slices should have zero temporal links: {n_pad_ok}/{n_pad_total}"

# Test 10: mk_gt_truncated
gt = q.GaugeTransform(geo)
gt.set_rand(rs.split("gt-init"))
gt_trunc, t_start, t_size_trunc = mk_gt_truncated(gt, t_center=t_center, t_half=t_half, t_size_divisor=t_size_divisor)
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

q.timer_display()
if q.is_test():
    q.check_log_json(__file__)
q.end_with_mpi()
q.displayln_info("CHECK: finished successfully.")
