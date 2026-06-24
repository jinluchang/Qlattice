#!/usr/bin/env python3

import qlat as q

from qlat_scripts.v1.load_data import (
    dict_flavor_inv_type,
    f_get_elem_wm,
    mk_field_norm_sqrt,
    check_cache_assign,
    get_prop_lookup_snk_src,
    get_prop_norm_lookup_snk_src,
    load_prop_wsrc_psel,
)

q.begin_with_mpi()

q.json_results_append("test load_data documentation examples")

q.json_results_append("dict_flavor_inv_type")

assert dict_flavor_inv_type["l"] == 0
assert dict_flavor_inv_type["u"] == 0
assert dict_flavor_inv_type["d"] == 0
assert dict_flavor_inv_type["s"] == 1
assert dict_flavor_inv_type["c"] == 2
q.json_results_append(f"dict_flavor_inv_type: {dict_flavor_inv_type}")

q.json_results_append("f_get_elem_wm with zero field")

result = f_get_elem_wm(0, 0)
assert result == 0
q.json_results_append(f"f_get_elem_wm(0, 0) = {result}")

q.json_results_append("mk_field_norm_sqrt with zero field")

result_norm = mk_field_norm_sqrt(0)
assert result_norm == 0
q.json_results_append(f"mk_field_norm_sqrt(0) = {result_norm}")

q.json_results_append("check_cache_assign")

cache = {}
check_cache_assign(cache, "key1", 42)
assert cache["key1"] == 42
check_cache_assign(cache, "key1", 42)
assert cache["key1"] == 42
cache["key2"] = "existing"
check_cache_assign(cache, "key2", "existing")
assert cache["key2"] == "existing"
q.json_results_append("check_cache_assign passed")

q.json_results_append("get_prop_lookup_snk_src direct lookup")

prop_lookup_cache = {}

def mock_get(pos_snk):
    return f"val_at_{pos_snk}"

prop_lookup_cache[("l", 0, "wall", "wall")] = mock_get
result = get_prop_lookup_snk_src(prop_lookup_cache, "l", ("wall", 3), ("wall", 0))
assert result == "val_at_3"
q.json_results_append(f"direct lookup: {result}")

q.json_results_append("get_prop_lookup_snk_src g5_herm fallback")

prop_lookup_cache2 = {}

def mock_get2(pos_snk):
    return f"val_at_{pos_snk}"

prop_lookup_cache2[("l", 7, "point-snk", "wall")] = mock_get2
result2 = get_prop_lookup_snk_src(
    prop_lookup_cache2, "l", ("point-snk", 7), ("wall", 5)
)
assert result2 == ("g5_herm", "val_at_5"), f"got {result2!r}"
q.json_results_append(f"g5_herm fallback: {result2}")

q.json_results_append("get_prop_norm_lookup_snk_src direct lookup")

prop_norm_lookup_cache = {}

def mock_norm_get(pos_snk):
    return 2.0

prop_norm_lookup_cache[("s", 0, "wall", "wall")] = mock_norm_get
result_norm_lu = get_prop_norm_lookup_snk_src(
    prop_norm_lookup_cache, "s", ("wall", 4), ("wall", 0)
)
assert abs(result_norm_lu - 2.0) < 1e-14
q.json_results_append(f"norm direct lookup: {result_norm_lu}")

q.json_results_append("get_prop_norm_lookup_snk_src g5_herm fallback")

prop_norm_lookup_cache2 = {}

def mock_norm_get2(pos_snk):
    return 0.9

prop_norm_lookup_cache2[("s", 5, "point-snk", "wall")] = mock_norm_get2
result_norm2 = get_prop_norm_lookup_snk_src(
    prop_norm_lookup_cache2, "s", ("point-snk", 5), ("wall", 3)
)
assert abs(result_norm2 - 0.9) < 1e-14, f"got {result_norm2!r}"
q.json_results_append(f"norm g5_herm fallback: {result_norm2}")

q.json_results_append("load_prop_wsrc_psel returns None when no data")

job_tag = "test-4nt8"
traj = 1000

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
gt = q.GaugeTransform(geo)
gt.set_unit()
psel = q.PointsSelection()
fsel = q.FieldSelection()
fsel.set_uniform(geo, 1)

result_load = load_prop_wsrc_psel(job_tag, traj, "l", psel=psel, fsel=fsel, gt=gt)
assert result_load is None
q.json_results_append(f"load_prop_wsrc_psel with no data: {result_load}")

q.check_log_json(__file__, check_eps=1e-10)

q.timer_display()

q.end_with_mpi()

q.displayln_info("CHECK: finished successfully.")
