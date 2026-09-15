#!/usr/bin/env python3

# Tests for the cqlat SelectedField interface (qlat/cqlat/selected-field.cpp):
#     set_add_sfield, set_mul_double_sfield, acc_field_sfield,
#     glb_sum_tslice_long_sfield
#
# The three arithmetic exports are reached through the public Python
# operators of ``SelectedFieldBase`` / ``FieldBase``:
#     sf1 += sf2          -> set_add_sfield
#     sf *= 2.0           -> set_mul_double_sfield
#     field += sf         -> acc_field_sfield
# and the time-slice global sum through ``SelectedFieldLong.glb_sum_tslice``:
#     sf_long.glb_sum_tslice(t_dir=3)  -> glb_sum_tslice_long_sfield
#
# All results are checked against numpy on the local part of the selection
# (``fsel.to_psel_local().xg_arr``) and the field data is built from *global*
# arrays so that the two MPI ranks hold different random values.

import gc
import os

import numpy as np

import qlat as q

check_eps = 1e-10

total_site = q.Coordinate([4, 4, 4, 8])

# The runner uses ``mpiexec -n 2`` so the [2, 1, 1, 1] layout is selected.
size_node_list = [
    [2, 1, 1, 1],
    [1, 1, 1, 1],
]

def mk_field_real_d_from_global(geo, a_g):
    # pick up the local part of a *global* pattern, so that the two MPI ranks
    # hold different values (a bare RngState is rank independent)
    f = q.FieldRealD(geo, 1)
    buf = np.asarray(f)
    for index in range(geo.local_volume):
        xg = geo.coordinate_g_from_l(geo.coordinate_from_index(index))
        buf[index, 0] = a_g[xg[0], xg[1], xg[2], xg[3]]
    return f

def mk_local_field_values(geo, a_g):
    # local values of a global array, in flat local site order
    buf = np.empty(geo.local_volume, dtype=np.float64)
    for index in range(geo.local_volume):
        xg = geo.coordinate_g_from_l(geo.coordinate_from_index(index))
        buf[index] = a_g[xg[0], xg[1], xg[2], xg[3]]
    return buf

def glb_tslice_sums(xg_arr, t_dir, values, t_size):
    # reference for glb_sum_tslice: sum ``values`` over the selected points,
    # grouped by the coordinate ``t_dir``
    rtn = np.zeros(t_size, dtype=np.int64)
    for i in range(xg_arr.shape[0]):
        rtn[int(xg_arr[i, t_dir])] += int(values[i])
    return rtn

q.begin_with_mpi(size_node_list)

geo = q.Geometry(total_site)
latt_size = tuple(int(geo.total_site[i]) for i in range(4))
q.json_results_append(f"cqlat-selected-field: geo.show()={geo.show()}")
rs = q.RngState("cqlat-selected-field")

# --- a random field selection with 16 points per time slice
fsel = q.FieldSelection()
fsel.set_rand(total_site, 16, rs.split("fsel"))
n_elems = int(fsel.n_elems)
n_elems_glb = float(q.glb_sum(float(n_elems)))
q.json_results_append("cqlat-selected-field: fsel glb n_elems", n_elems_glb)
assert n_elems_glb > 0.0
assert fsel.geo.total_site == total_site

# local global-coordinates of the selected elements, in fsel.indices order
psel_l = fsel.to_psel_local()
xg_arr = np.array(psel_l.xg_arr, dtype=np.int64)
assert xg_arr.shape == (n_elems, 4)
assert psel_l.points_dist_type == "l"

# --- rank-distinct random data coming from global arrays
rs_g = q.RngState("cqlat-selected-field-global")
a0_g = rs_g.split("a0").u_rand_arr(latt_size)
a1_g = rs_g.split("a1").u_rand_arr(latt_size)

f_a = mk_field_real_d_from_global(geo, a0_g)
f_b = mk_field_real_d_from_global(geo, a1_g)

sf_a = q.SelectedField(q.ElemTypeRealD, fsel, 1)
sf_a @= f_a
sf_b = q.SelectedField(q.ElemTypeRealD, fsel, 1)
sf_b @= f_b

sel = (xg_arr[:, 0], xg_arr[:, 1], xg_arr[:, 2], xg_arr[:, 3])
ref_a = a0_g[sel]
ref_b = a1_g[sel]

arr_a = np.asarray(sf_a)
assert arr_a.shape == (n_elems, 1)
assert np.array_equal(arr_a[:, 0], ref_a)
q.json_results_append(
    "cqlat-selected-field: selected data sig",
    q.get_data_sig_arr(np.asarray(sf_a), q.RngState("cqlat-selected-field-sig"), 3),
    check_eps,
)

# --- set_add_sfield: sf1 += sf2
sf_c = sf_a.copy()
sf_c += sf_b
ref_sum = ref_a + ref_b
err_add = float(np.max(np.abs(np.asarray(sf_c)[:, 0] - ref_sum)))
q.json_results_append(
    "cqlat-selected-field: set_add_sfield max err", err_add, check_eps
)
assert err_add < 1e-15
q.json_results_append(
    "cqlat-selected-field: set_add_sfield sum", float(ref_sum.sum()), check_eps
)

# --- set_mul_double_sfield: sf *= factor
sf_c *= 2.0
ref_mul = 2.0 * ref_sum
err_mul = float(np.max(np.abs(np.asarray(sf_c)[:, 0] - ref_mul)))
q.json_results_append(
    "cqlat-selected-field: set_mul_double_sfield max err", err_mul, check_eps
)
assert err_mul == 0.0
sf_c *= -1.5
ref_mul2 = -1.5 * ref_mul
err_mul2 = float(np.max(np.abs(np.asarray(sf_c)[:, 0] - ref_mul2)))
q.json_results_append(
    "cqlat-selected-field: set_mul_double_sfield (-1.5) max err", err_mul2, check_eps
)
assert err_mul2 < 1e-15
q.json_results_append(
    "cqlat-selected-field: set_mul_double_sfield sum",
    float(np.asarray(sf_c)[:, 0].sum()),
)

# --- acc_field_sfield: field += sf
f_c = mk_field_real_d_from_global(geo, a1_g)
f_c += sf_c
exp_full = a1_g.copy()
exp_full[sel] += ref_mul2
exp_local = mk_local_field_values(geo, exp_full)
got_local = np.asarray(f_c)[:, 0]
err_acc = float(np.max(np.abs(got_local - exp_local)))
q.json_results_append(
    "cqlat-selected-field: acc_field_sfield max err", err_acc, check_eps
)
assert err_acc == 0.0
# only the selected sites may change
n_changed = int(np.sum(got_local != mk_local_field_values(geo, a1_g)))
n_changed = int(q.glb_sum(float(n_changed)))
q.json_results_append(f"cqlat-selected-field: acc_field_sfield n changed sites = {n_changed}")
assert n_changed == int(n_elems_glb)
q.json_results_append(
    "cqlat-selected-field: acc_field_sfield sum", float(got_local.sum()), check_eps
)

# --- Field.__isub__ with a SelectedField
#     (regression: the subtraction used to be applied and then an unconditional
#     ``assert False`` raised, so the caller got a half-applied result plus a
#     bogus AssertionError)
f_d = mk_field_real_d_from_global(geo, a1_g)
f_d -= sf_c
exp_sub_g = a1_g.copy()
exp_sub_g[sel] -= np.asarray(sf_c)[:, 0]
err_sub = float(
    np.max(np.abs(np.asarray(f_d)[:, 0] - mk_local_field_values(geo, exp_sub_g)))
)
assert err_sub == 0.0, err_sub
q.json_results_append(
    "cqlat-selected-field: field -= selected_field max err", err_sub, check_eps
)
q.json_results_append(
    "cqlat-selected-field: field -= selected_field sum",
    float(np.asarray(f_d)[:, 0].sum()),
    check_eps,
)

# --- acc_field_spfield: Field += / -= SelectedPoints
#     (the cqlat function needs the geometry and the points selection, which
#     Field.__iadd__ / __isub__ take from the SelectedPoints object)
sp_a = q.SelectedPoints(q.ElemTypeRealD, psel_l, 1)
assert np.asarray(sp_a).shape == (n_elems, 1)
np.asarray(sp_a)[:, 0] = ref_a

f_e = mk_field_real_d_from_global(geo, a1_g)
f_e += sp_a
exp_sp_g = a1_g.copy()
exp_sp_g[sel] += ref_a
err_sp_add = float(
    np.max(np.abs(np.asarray(f_e)[:, 0] - mk_local_field_values(geo, exp_sp_g)))
)
assert err_sp_add == 0.0, err_sp_add
q.json_results_append(
    "cqlat-selected-field: acc_field_spfield (+=) max err", err_sp_add, check_eps
)
q.json_results_append(
    "cqlat-selected-field: acc_field_spfield (+=) sum",
    float(np.asarray(f_e)[:, 0].sum()),
    check_eps,
)

f_f = mk_field_real_d_from_global(geo, a1_g)
f_f -= sp_a
exp_sp_sub_g = a1_g.copy()
exp_sp_sub_g[sel] -= ref_a
err_sp_sub = float(
    np.max(np.abs(np.asarray(f_f)[:, 0] - mk_local_field_values(geo, exp_sp_sub_g)))
)
assert err_sp_sub == 0.0, err_sp_sub
q.json_results_append(
    "cqlat-selected-field: acc_field_spfield (-=) max err", err_sp_sub, check_eps
)
q.json_results_append(
    "cqlat-selected-field: acc_field_spfield (-=) sum",
    float(np.asarray(f_f)[:, 0].sum()),
    check_eps,
)

# from an empty Field: acc_field_spfield initializes it
f_g = q.FieldRealD()
f_g += sp_a
exp_sp_empty_g = np.zeros(latt_size, dtype=np.float64)
exp_sp_empty_g[sel] = ref_a
err_sp_empty = float(
    np.max(np.abs(np.asarray(f_g)[:, 0] - mk_local_field_values(geo, exp_sp_empty_g)))
)
assert err_sp_empty == 0.0, err_sp_empty
q.json_results_append(
    "cqlat-selected-field: acc_field_spfield (empty field) max err",
    err_sp_empty,
    check_eps,
)
q.json_results_append(
    "cqlat-selected-field: acc_field_spfield (empty field) sum",
    float(np.asarray(f_g)[:, 0].sum()),
    check_eps,
)

# --- start from an empty Field: acc_field_sfield initializes it
f_empty = q.FieldRealD()
sf_d = sf_b.copy()
sf_d *= 2.0
f_empty += sf_d
exp_empty_g = np.zeros(latt_size, dtype=np.float64)
exp_empty_g[sel] = 2.0 * ref_b
assert np.array_equal(
    np.asarray(f_empty)[:, 0], mk_local_field_values(geo, exp_empty_g)
)
q.json_results_append(
    "cqlat-selected-field: acc_field_sfield init from empty sum",
    float(np.asarray(f_empty)[:, 0].sum()),
    check_eps,
)

# --- glb_sum_tslice_long_sfield: SelectedFieldLong.glb_sum_tslice
sf_long = q.SelectedField(q.ElemTypeLong, fsel, 1)
v_g = np.floor(1000.0 * rs_g.split("v-long").u_rand_arr(latt_size)).astype(np.int64)
arr_long = np.asarray(sf_long)
assert arr_long.shape == (n_elems, 1)
for idx in range(n_elems):
    xg = xg_arr[idx]
    arr_long[idx, 0] = v_g[xg[0], xg[1], xg[2], xg[3]]

# the global selection gathers the contributions of every MPI rank
xg_all = np.array(fsel.to_psel().xg_arr, dtype=np.int64)
n_all = float(q.glb_sum(float(n_elems)))
assert xg_all.shape[0] == int(n_all)
v_sel_all = v_g[xg_all[:, 0], xg_all[:, 1], xg_all[:, 2], xg_all[:, 3]]

for t_dir in [3, 0]:
    t_size = int(total_site[t_dir])
    sp = sf_long.glb_sum_tslice(t_dir=t_dir)
    got_t = np.array(np.asarray(sp)[:, 0], dtype=np.int64)
    assert got_t.shape == (t_size,)
    exp_t = glb_tslice_sums(xg_all, t_dir, v_sel_all, t_size)
    err_tslice = float(np.max(np.abs(got_t - exp_t)))
    q.json_results_append(
        f"cqlat-selected-field: glb_sum_tslice_long_sfield t_dir={t_dir} max err",
        err_tslice,
        check_eps,
    )
    assert err_tslice == 0.0
    q.json_results_append(
        f"cqlat-selected-field: glb_sum_tslice_long_sfield t_dir={t_dir} total",
        float(got_t.sum()),
    )
    q.json_results_append(
        f"cqlat-selected-field: glb_sum_tslice_long_sfield t_dir={t_dir} sig",
        q.get_data_sig_arr(
            got_t, q.RngState(f"cqlat-selected-field-tslice-{t_dir}"), 3
        ),
        check_eps,
    )

del sf_a, sf_b, sf_c, sf_d, sf_long, sp_a, f_a, f_b, f_c, f_d, f_e, f_f, f_g, f_empty
gc.collect()
q.json_results_append("cqlat-selected-field: gc.collect() done")
q.json_results_append(f"cqlat-selected-field: cwd={os.path.basename(os.getcwd())}")

q.timer_display()
if q.is_test():
    q.check_log_json(__file__)
q.end_with_mpi()
q.displayln_info("CHECK: finished successfully.")
