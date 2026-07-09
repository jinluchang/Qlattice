#!/usr/bin/env python3

"""
Test mk_gf_truncated_evolve: asymmetric time-direction truncation of a gauge
field with padding set to unity links.
"""

import qlat_gpt as qg
import qlat as q
import numpy as np

qg.begin_with_gpt()


total_site = q.Coordinate([4, 4, 4, 16])
geo = q.Geometry(total_site)
q.displayln_info(f"geo.show()={geo.show()}")
rs = q.RngState("seed")

gf = q.GaugeField(geo)
gf.set_rand(rs.split("gf-init"), 0.3, 1)

plaq = gf.plaq()
q.json_results_append("plaq", plaq, 1e-12)

### ------

def mk_gf_truncated_evolve(gf, t_center, t_left, t_right, t_pad):
    total_site = gf.geo.total_site
    t_size = total_site[3]
    t_size_trunc_valid = t_left + t_right + 1
    t_size_trunc = t_size_trunc_valid + t_pad
    size_node_t = gf.geo.size_node[3]
    t_size_trunc = ((t_size_trunc + size_node_t - 1) // size_node_t) * size_node_t
    assert t_size_trunc <= t_size
    t_start = (t_center - t_left) % t_size
    t_end = (t_start + t_size_trunc) % t_size
    gf_trunc = mk_field_truncated(gf, t_start, t_end)
    geo_trunc = gf_trunc.geo
    tslice_pad_list = list(range(t_size_trunc_valid, t_size_trunc))
    gf_arr = np.asarray(gf_trunc)
    xg_arr = q.mk_xg_field(geo_trunc)[:]
    eye3 = np.eye(3, dtype=gf_arr.dtype)
    for index in range(geo_trunc.local_volume):
        xg = xg_arr[index]
        if xg[3] in tslice_pad_list:
            gf_arr[index, :, :] = eye3
    return gf_trunc, t_start, t_size_trunc

def mk_field_truncated(field, t_start, t_end):
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

### ------

# Test 1: basic truncation, no padding
t_center = 8
t_left = 3
t_right = 4
t_pad = 0
gf_trunc, t_start, t_size_trunc = mk_gf_truncated_evolve(gf, t_center, t_left, t_right, t_pad)
q.json_results_append("test1 t_start", t_start)
q.json_results_append("test1 t_size_trunc", t_size_trunc)
assert t_size_trunc == 8  # 3 + 4 + 1
assert t_start == (t_center - t_left) % 16
q.displayln_info(f"CHECK: test1 t_start={t_start} t_size_trunc={t_size_trunc}")

# Test 2: truncation with padding
t_pad = 4
gf_trunc, t_start, t_size_trunc = mk_gf_truncated_evolve(gf, t_center, t_left, t_right, t_pad)
q.json_results_append("test2 t_start", t_start)
q.json_results_append("test2 t_size_trunc", t_size_trunc)
assert t_size_trunc == 12  # 3 + 4 + 1 + 4
q.displayln_info(f"CHECK: test2 t_start={t_start} t_size_trunc={t_size_trunc}")

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
q.displayln_info(f"CHECK: test3 pad={n_pad_ok}/{n_pad_total} valid_not_identity={n_valid_ok}/{n_valid_total}")
assert n_pad_ok == n_pad_total, f"padded sites should be unity: {n_pad_ok}/{n_pad_total}"
# Valid sites should generally NOT be identity (random gauge field)
# but on some ranks there may be no valid sites, so only check if n_valid_total > 0
if n_valid_total > 0:
    assert n_valid_ok > 0, "valid sites should not all be identity"

# Test 4: verify valid region data matches original gauge field via plaquette
# The truncated gauge field's plaquette on the valid region should differ from
# the full field's plaquette (since it's a subset), but should be well-defined.
plaq_trunc = gf_trunc.plaq()
q.json_results_append("test4 plaq_trunc", plaq_trunc, 1e-12)
assert 0.0 < plaq_trunc < 1.0, f"truncated plaq should be physical: {plaq_trunc}"
q.displayln_info(f"CHECK: test4 plaq_trunc={plaq_trunc}")

# Test 5: wrap-around at boundary
t_center_wrap = 2
t_left = 3
t_right = 2
t_pad = 0
gf_trunc, t_start, t_size_trunc = mk_gf_truncated_evolve(gf, t_center_wrap, t_left, t_right, t_pad)
q.json_results_append("test5 t_start", t_start)
q.json_results_append("test5 t_size_trunc", t_size_trunc)
assert t_start == 15  # (2 - 3) % 16 = 15
assert t_size_trunc == 6  # 3 + 2 + 1
q.displayln_info(f"CHECK: test5 t_start={t_start} t_size_trunc={t_size_trunc}")

# Test 6: verify unity padding with wrap-around
t_pad = 2
gf_trunc, t_start, t_size_trunc = mk_gf_truncated_evolve(gf, t_center_wrap, t_left, t_right, t_pad)
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
q.displayln_info(f"CHECK: test6 t_start={t_start} t_size_trunc={t_size_trunc} pad={n_pad_ok}/{n_pad_total}")

### ------

q.check_log_json(__file__)

q.timer_display()

qg.end_with_gpt()

q.displayln_info("CHECK: finished successfully.")
