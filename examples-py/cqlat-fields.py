#!/usr/bin/env python3

# Tests for the field helpers that the qlat Python API still relies on
# through the cqlat extension, covering:
#   qlat/cqlat/field-utils.cpp:
#       make_field_expand_comm_plan, set_marks_field_all,
#       refresh_expanded_field, refresh_expanded_1_field,
#       merge_fields_ms_field
#   qlat/cqlat/field-double.cpp:
#       multiply_double_field, invert_double_field
#   qlat/cqlat/field.cpp:
#       get_mview_field
#
# The checks use real content wherever possible:
#   * the expanded halo is filled with a *global* periodic pattern, so the two
#     MPI ranks hold different data and the halo values can only come from the
#     neighbouring rank; refresh_expanded / refresh_expanded_1 / the explicit
#     CommPlan are all compared against that pattern,
#   * set_marks_field_all is compared against a numpy reimplementation of the
#     marking rule (mark expanded offsets which are not local),
#   * Field.shift(is_reflect=True) is compared against a numpy reflection
#     new_f[xg] == f[mod(-xg, total_site)],
#   * merge_fields_ms_field / q.sqrt_field / multiply_double_field /
#     invert_double_field are compared against numpy reimplementations of the
#     elementwise operations,
#   * get_mview_field is exercised by writing through the returned writable
#     memoryview and reading the field back,
#   * Field.write_direct / read_direct are checked with a round trip (including
#     a Coordinate new_size_node) and the returned byte counts.
#
# ``CommMarks`` is re-exported as ``q.CommMarks`` (it used to be missing from
# the ``__all__`` of ``qlat/qlat/c.py.in``).

import gc
import os

import numpy as np

import qlat as q
import qlat.c as qc

check_eps = 1e-10

# ---------------------------------------------------------------------------
# helpers (defined before q.begin_with_mpi; they must not build a Geometry)
# ---------------------------------------------------------------------------

def mk_field_from_global(geo, g, mult, is_complex=False):
    """
    Field holding the local part of the global pattern ``g``.
    ``g`` has shape ``tuple(geo.total_site) + (mult,)``.
    A global pattern is required so that the two MPI ranks hold different data
    (a bare RngState gives identical values on all ranks).
    """
    if is_complex:
        f = q.FieldComplexD(geo, mult)
    else:
        f = q.FieldRealD(geo, mult)
    buf = np.asarray(f)
    for index in range(geo.local_volume):
        xg = geo.coordinate_g_from_l(geo.coordinate_from_index(index))
        for m in range(mult):
            buf[index, m] = g[xg[0], xg[1], xg[2], xg[3], m]
    return f

def flat_xg(geo, is_expanded=False):
    """
    Global coordinate (mod total_site) of every flat local field offset.
    Valid for eo == 0 geometries.  With ``is_expanded=True`` the expanded data
    layout ``node_site + expansion_left + expansion_right`` is used.
    """
    assert int(geo.eo) == 0
    node_site = [int(geo.node_site[i]) for i in range(4)]
    expansion_left = [int(geo.expansion_left[i]) for i in range(4)]
    expansion_right = [int(geo.expansion_right[i]) for i in range(4)]
    coor_node = [int(geo.coor_node[i]) for i in range(4)]
    total_site = [int(geo.total_site[i]) for i in range(4)]
    if is_expanded:
        size = [node_site[i] + expansion_left[i] + expansion_right[i] for i in range(4)]
        n_points = int(geo.local_volume_expanded)
    else:
        expansion_left = [0, 0, 0, 0]
        size = node_site
        n_points = int(geo.local_volume)
    coords = np.zeros((n_points, 4), dtype=np.int64)
    index = np.arange(n_points, dtype=np.int64)
    x_list = []
    for i in range(4):
        x_list.append(index % size[i])
        index = index // size[i]
    for i in range(4):
        xg = x_list[i] - expansion_left[i] + coor_node[i] * node_site[i]
        coords[:, i] = xg % total_site[i]
    return coords

def pattern_at(g, coords):
    return g[coords[:, 0], coords[:, 1], coords[:, 2], coords[:, 3]]

# ---------------------------------------------------------------------------

# The runner uses ``mpiexec -n 2`` so the [2, 1, 1, 1] layout is selected.
# Only direction 0 is split, therefore the halo communication test below only
# expands direction 0 (a single-node direction is aliased in the field buffer).
q.begin_with_mpi([[2, 1, 1, 1], [1, 1, 1, 1]])

geo = q.Geometry(q.Coordinate([4, 4, 4, 8]))
latt_size = tuple(int(geo.total_site[i]) for i in range(4))
q.json_results_append(f"cqlat-fields: geo.show()={geo.show()}")

rs = q.RngState("cqlat-fields")
coords_l = flat_xg(geo)

# --- field expansion, communication plan and halo refresh ------------------
# The expanded field is filled with a sentinel, the local part is assigned,
# and then the halo is refreshed from the neighbouring rank.  The halo global
# coordinates are outside this node's local volume, so the comparison against
# the global pattern can only pass if the MPI communication happened.

g_expand0 = rs.split("expand-0").u_rand_arr(latt_size)
g_expand1 = rs.split("expand-1").u_rand_arr(latt_size)
g_expand = np.stack([g_expand0, g_expand1], axis=-1)
f_expand = mk_field_from_global(geo, g_expand, 2)

expansion_left = [1, 0, 0, 0]
expansion_right = [1, 0, 0, 0]
geo_e = q.geo_resize(geo, expansion_left, expansion_right)
coords_e = flat_xg(geo_e, is_expanded=True)
g_expand_e = pattern_at(g_expand, coords_e)
n_halo = int(geo_e.local_volume_expanded) - int(geo_e.local_volume)

sentinel = -1.2345e30

f_e = q.field_expanded(f_expand, expansion_left, expansion_right)
np.asarray(f_e)[:] = sentinel
f_e @= f_expand
arr_e = np.asarray(f_e)
is_halo = arr_e[:, 0] < sentinel * 0.5
assert int(is_halo.sum()) == n_halo, (int(is_halo.sum()), n_halo)
assert np.array_equal(arr_e[~is_halo], g_expand_e[~is_halo])
q.json_results_append("cqlat-fields: field_expanded n_halo", float(n_halo), 1e-12)

cm = q.CommMarks(geo_e, 2)
qc.set_marks_field_all(cm, geo_e, 2, "cqlat-fields")
marks = np.asarray(cm)
assert marks.shape == (int(geo_e.local_volume_expanded), 2), marks.shape
assert marks.dtype == np.int8, marks.dtype
expected_marks = np.zeros(marks.shape, dtype=np.int8)
expected_marks[is_halo, :] = 1
assert np.array_equal(marks, expected_marks)
q.json_results_append(
    "cqlat-fields: set_marks_field_all n_marks", q.glb_sum(float(marks.sum())), 1e-12
)

cp = q.make_field_expand_comm_plan(cm)
q.refresh_expanded(f_e, cp)
assert np.array_equal(np.asarray(f_e), g_expand_e)
q.json_results_append(
    "cqlat-fields: refresh_expanded_field(f, comm_plan) max err", 0.0, 1e-12
)
q.json_results_append(
    "cqlat-fields: refresh_expanded_field(f, comm_plan) halo sum",
    q.glb_sum(float(np.asarray(f_e)[is_halo].sum())),
    1e-12,
)

f_e2 = q.field_expanded(f_expand, expansion_left, expansion_right)
q.refresh_expanded(f_e2)
assert np.array_equal(np.asarray(f_e2), g_expand_e)
q.json_results_append("cqlat-fields: refresh_expanded_field(f) max err", 0.0, 1e-12)

f_e3 = q.field_expanded(f_expand, expansion_left, expansion_right)
q.refresh_expanded_1(f_e3)
assert np.array_equal(np.asarray(f_e3), g_expand_e)
assert np.array_equal(np.asarray(f_e3), np.asarray(f_e))
q.json_results_append("cqlat-fields: refresh_expanded_1_field(f) max err", 0.0, 1e-12)
q.json_results_append(
    "cqlat-fields: refresh_expanded_1 == refresh_expanded (1D halo)", 1.0, 1e-12
)

del f_e3
gc.collect()

# --- reflection: new_f[xg] == f[mod(-xg, total_site)] -----------------------

g_reflect = rs.split("reflect").u_rand_arr(latt_size)
f_reflect = mk_field_from_global(geo, g_reflect[..., None], 1)
reflect_idx = tuple((-np.arange(n)) % n for n in latt_size)
g_reflect_r = g_reflect[np.ix_(*reflect_idx)]

f_refl = f_reflect.shift(is_reflect=True)
exp_refl = pattern_at(g_reflect_r, coords_l)
err_refl = float(np.max(np.abs(np.asarray(f_refl)[:, 0] - exp_refl)))
assert err_refl == 0.0, err_refl

f_refl2 = f_refl.shift(is_reflect=True)
assert np.array_equal(np.asarray(f_refl2), np.asarray(f_reflect))

q.json_results_append("cqlat-fields: reflect_field max err", err_refl, 1e-12)
q.json_results_append("cqlat-fields: reflect_field twice is identity", 1.0, 1e-12)
q.json_results_append("cqlat-fields: reflect_field qnorm", f_refl.qnorm(), 1e-12)

# --- merge_fields_ms_field -------------------------------------------------
# f.get_elem(x, m) == fms[m][0].get_elem(x, fms[m][1])

g_merge0 = rs.split("merge-0").u_rand_arr(latt_size)
g_merge1 = rs.split("merge-1").u_rand_arr(latt_size)
g_merge2 = rs.split("merge-2").u_rand_arr(latt_size)
f_merge0 = mk_field_from_global(geo, np.stack([g_merge0, g_merge1], axis=-1), 2)
f_merge1 = mk_field_from_global(geo, g_merge2[..., None], 1)

f_merge = q.FieldRealD(geo, 4)
q.merge_fields_ms(f_merge, [(f_merge0, 1), (f_merge1, 0), (f_merge0, 0), (f_merge1, 0)])
exp_merge = np.stack(
    [
        np.asarray(f_merge0)[:, 1],
        np.asarray(f_merge1)[:, 0],
        np.asarray(f_merge0)[:, 0],
        np.asarray(f_merge1)[:, 0],
    ],
    axis=1,
)
err_merge = float(np.max(np.abs(np.asarray(f_merge) - exp_merge)))
assert err_merge == 0.0, err_merge
q.json_results_append("cqlat-fields: merge_fields_ms_field max err", err_merge, 1e-12)
q.json_results_append(
    "cqlat-fields: merge_fields_ms_field sum",
    q.glb_sum(float(np.asarray(f_merge).sum())),
    1e-12,
)

# --- q.sqrt_field: f = sqrt(f1), elementwise -------------------------------

g_sqrt0 = 0.25 + rs.split("sqrt-0").u_rand_arr(latt_size)
g_sqrt1 = 0.50 + rs.split("sqrt-1").u_rand_arr(latt_size)
f_sqrt_in = mk_field_from_global(geo, np.stack([g_sqrt0, g_sqrt1], axis=-1), 2)

f_sqrt = q.sqrt_field(f_sqrt_in)
assert np.asarray(f_sqrt).shape == np.asarray(f_sqrt_in).shape
assert f_sqrt.multiplicity == f_sqrt_in.multiplicity
err_sqrt = float(np.max(np.abs(np.asarray(f_sqrt) - np.sqrt(np.asarray(f_sqrt_in)))))
assert err_sqrt == 0.0, err_sqrt
q.json_results_append("cqlat-fields: set_sqrt_field max err", err_sqrt, 1e-12)
q.json_results_append(
    "cqlat-fields: set_sqrt_field sum",
    q.glb_sum(float(np.asarray(f_sqrt).sum())),
    1e-12,
)

# --- multiply_double_field: f *= factor, elementwise -----------------------

g_mul_a = np.stack(
    [
        rs.split("mul-a0").u_rand_arr(latt_size),
        rs.split("mul-a1").u_rand_arr(latt_size),
    ],
    axis=-1,
)
g_mul_b = np.stack(
    [
        rs.split("mul-b0").u_rand_arr(latt_size),
        rs.split("mul-b1").u_rand_arr(latt_size),
    ],
    axis=-1,
)
f_mul_a = mk_field_from_global(geo, g_mul_a, 2)
f_mul_b = mk_field_from_global(geo, g_mul_b, 2)
arr_mul_a = np.asarray(f_mul_a).copy()
arr_mul_b = np.asarray(f_mul_b).copy()
q.field_double.multiply_double(f_mul_a, f_mul_b)
err_mul = float(np.max(np.abs(np.asarray(f_mul_a) - arr_mul_a * arr_mul_b)))
assert err_mul == 0.0, err_mul
q.json_results_append("cqlat-fields: multiply_double_field max err", err_mul, 1e-12)
q.json_results_append(
    "cqlat-fields: multiply_double_field sum",
    q.glb_sum(float(np.asarray(f_mul_a).sum())),
    1e-12,
)

# --- invert_double_field: f = 1 / f, elementwise ---------------------------

f_inv = mk_field_from_global(geo, g_mul_a + 1.0, 2)
arr_inv = np.asarray(f_inv).copy()
q.field_double.invert_double(f_inv)
err_inv = float(np.max(np.abs(np.asarray(f_inv) - 1.0 / arr_inv)))
assert err_inv == 0.0, err_inv
q.field_double.invert_double(f_inv)
err_inv2 = float(np.max(np.abs(np.asarray(f_inv) - arr_inv)))
assert err_inv2 < 1e-12, err_inv2
q.json_results_append("cqlat-fields: invert_double_field max err", err_inv, 1e-12)
q.json_results_append(
    "cqlat-fields: invert_double_field twice max err", err_inv2, 1e-12
)

# --- complex field *= complex factor ---------------------------------------

g_mulc = (0.2 + rs.split("mulc-re").u_rand_arr(latt_size)) + 1j * (
    0.3 + rs.split("mulc-im").u_rand_arr(latt_size)
)
f_mulc = mk_field_from_global(geo, g_mulc[..., None], 1, is_complex=True)
arr_mulc = np.asarray(f_mulc).copy()
z_mulc = 0.7 - 0.4j

f_mulc_1 = f_mulc.copy()
f_mulc_1 *= z_mulc
err_mulc = float(np.max(np.abs(np.asarray(f_mulc_1) - arr_mulc * z_mulc)))
assert err_mulc < 1e-14, err_mulc
q.json_results_append(
    "cqlat-fields: complex field *= factor max err", err_mulc, check_eps
)
q.json_results_append(
    "cqlat-fields: complex field *= factor qnorm", f_mulc_1.qnorm(), 1e-12
)

# --- a real field rejects a genuinely complex factor -----------------------
#     (earlier this reached ``operator*=(RealD&, const ComplexD&)`` whose body
#     was ``assert(false)``: a hard abort, or a silent no-op with NDEBUG)
f_rej = mk_field_from_global(geo, g_mulc.real[..., None], 1)
arr_rej = np.asarray(f_rej).copy()
try:
    f_rej *= z_mulc
except ValueError:
    rejected = 1.0
    q.json_results_append("cqlat-fields: real field *= complex factor rejected", 1.0)
else:
    rejected = 0.0
    q.json_results_append("cqlat-fields: real field *= complex factor rejected", 0.0)
assert rejected == 1.0
assert np.array_equal(np.asarray(f_rej), arr_rej), "rejected factor must not modify"
# a real valued complex factor is still accepted for a real field
f_rej *= complex(2.0, 0.0)
err_rej = float(np.max(np.abs(np.asarray(f_rej) - arr_rej * 2.0)))
assert err_rej == 0.0, err_rej
q.json_results_append(
    "cqlat-fields: real field *= complex(2.0, 0.0) max err", err_rej, check_eps
)

# --- get_mview_field: writable memoryview of the field data ----------------

g_mview = np.stack(
    [
        rs.split("mview-0").u_rand_arr(latt_size),
        rs.split("mview-1").u_rand_arr(latt_size),
    ],
    axis=-1,
)
f_mview = mk_field_from_global(geo, g_mview, 2)
target_mview = 3.0 + np.asarray(f_mview).copy()

mv = f_mview.mview()
raw_mv = np.asarray(mv)
assert mv.readonly is False
assert raw_mv.flags["WRITEABLE"]
assert not raw_mv.flags["OWNDATA"]
assert np.shares_memory(raw_mv, np.asarray(f_mview))
vals_mview = raw_mv.view(np.float64)
assert vals_mview.size == int(geo.local_volume) * 2, vals_mview.size
vals_mview[:] = target_mview.reshape(-1)
back_mview = np.asarray(f_mview)
assert np.array_equal(back_mview, target_mview)
q.json_results_append(
    "cqlat-fields: get_mview_field nbytes", float(raw_mv.nbytes), 1e-12
)
q.json_results_append("cqlat-fields: get_mview_field write/read max err", 0.0, 1e-12)
q.json_results_append("cqlat-fields: get_mview_field qnorm", f_mview.qnorm(), 1e-12)

# --- Field.write_direct / read_direct (qlat distributed field IO) ----------

g_save = np.stack(
    [
        rs.split("save-0").u_rand_arr(latt_size),
        rs.split("save-1").u_rand_arr(latt_size),
    ],
    axis=-1,
)
f_save = mk_field_from_global(geo, g_save, 2)
arr_save = np.asarray(f_save).copy()
n_bytes_expected = float(int(geo.total_volume) * 2 * 8)

path_plain = os.path.join(".", "cqlat-fields-save.field")
if os.path.isfile(path_plain):
    os.remove(path_plain)
n_bytes = f_save.write_direct(path_plain)
assert float(n_bytes) == n_bytes_expected, (n_bytes, n_bytes_expected)
f_load = q.FieldRealD()
n_read = f_load.read_direct(path_plain)
assert float(n_read) == n_bytes_expected, (n_read, n_bytes_expected)
assert f_load.geo == f_save.geo
assert f_load.multiplicity == f_save.multiplicity
err_load = float(np.max(np.abs(np.asarray(f_load) - arr_save)))
assert err_load == 0.0, err_load
q.json_results_append("cqlat-fields: write_direct bytes", float(n_bytes), 1e-12)
q.json_results_append("cqlat-fields: read_direct bytes", float(n_read), 1e-12)
q.json_results_append("cqlat-fields: write_direct/read_direct max err", err_load, 1e-12)

# a non-trivial new_size_node repartitions the field on disk and must load back
path_nsn = os.path.join(".", "cqlat-fields-save-nsn.field")
n_bytes_nsn = f_save.write_direct(path_nsn, q.Coordinate([1, 1, 1, 8]))
assert float(n_bytes_nsn) == n_bytes_expected, (n_bytes_nsn, n_bytes_expected)
f_load_nsn = q.FieldRealD()
n_read_nsn = f_load_nsn.read_direct(path_nsn)
assert float(n_read_nsn) == n_bytes_expected, (n_read_nsn, n_bytes_expected)
assert f_load_nsn.geo == f_save.geo
err_load_nsn = float(np.max(np.abs(np.asarray(f_load_nsn) - arr_save)))
assert err_load_nsn == 0.0, err_load_nsn
q.json_results_append(
    "cqlat-fields: write_direct(new_size_node) bytes", float(n_bytes_nsn), 1e-12
)
q.json_results_append(
    "cqlat-fields: write_direct(new_size_node)/read_direct max err",
    err_load_nsn,
    1e-12,
)

del f_load, f_load_nsn
gc.collect()

q.timer_display()
if q.is_test():
    q.check_log_json(__file__)
q.end_with_mpi()
q.displayln_info("CHECK: finished successfully.")
