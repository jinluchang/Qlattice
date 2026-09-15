#!/usr/bin/env python3

# Tests for the cqlat propagator random-U1 interface
# (qlat/qlat/propagator.pyx):
#     set_rand_u1_src_psel, set_rand_u1_sol_psel,
#     set_rand_u1_src_fsel, set_rand_u1_sol_fsel,
#     flip_tpbc_with_tslice_sp_prop, flip_tpbc_with_tslice_s_prop
#
# They are reached through the public Python wrappers:
#     q.mk_rand_u1_src(sel, rs)          -> set_rand_u1_src_psel / _fsel
#     q.get_rand_u1_sol(prop_sol, fu1, sel)
#                                        -> set_rand_u1_sol_psel / _fsel
#     q.flip_tpbc_with_tslice(prop, tslice)
#                                        -> flip_tpbc_with_tslice_sp_prop
#                                           flip_tpbc_with_tslice_s_prop
#
# ``set_rand_u1_src_*`` sets ``prop(x) = fu1(x) * 1`` and stores the random
# phase in ``fu1``; ``set_rand_u1_sol_*`` multiplies the propagator at the
# selected points by ``conj(fu1)``.  Since ``|fu1| == 1``, the resulting
# "self loop" is the identity WilsonMatrix, which is checked here, together
# with the explicit ``sol = prop * conj(fu1)`` relation.
#
# ``flip_tpbc_with_tslice`` negates the propagator on the half of the time
# direction away from ``tslice_flip_tpbc``; the flipped range is recomputed
# independently below from the global time coordinate of every selected
# point.

import gc
import os

import numpy as np

import qlat as q

check_eps = 1e-10

total_site = q.Coordinate([4, 4, 4, 8])
t_size = int(total_site[3])
t_size_half = t_size // 2

# The runner uses ``mpiexec -n 2`` so the [2, 1, 1, 1] layout is selected.
size_node_list = [
    [2, 1, 1, 1],
    [1, 1, 1, 1],
]

def flip_t_range(tslice_flip_tpbc):
    # reference for set_t_range_flip_tpbc_with_tslice (qcd-prop.cpp)
    assert 0 <= tslice_flip_tpbc < t_size
    if tslice_flip_tpbc + t_size_half < t_size:
        return (tslice_flip_tpbc + t_size_half, t_size)
    if tslice_flip_tpbc - t_size_half >= 0:
        return (0, tslice_flip_tpbc - t_size_half)
    raise AssertionError(tslice_flip_tpbc)

def flip_signs(t_arr, tslice_flip_tpbc):
    t_start, t_stop = flip_t_range(tslice_flip_tpbc)
    return np.where((t_arr >= t_start) & (t_arr < t_stop), -1.0, 1.0)

def check_flip_prop(prop, t_arr, arr0, tag):
    # flip twice with the same tslice exactly restores the original data, and
    # a single flip is a sign flip on the expected time range
    for tslice_flip_tpbc in [7, 0, 5]:
        signs = flip_signs(t_arr, tslice_flip_tpbc)
        n_flip = float(q.glb_sum(float(np.sum(signs < 0.0))))
        q.json_results_append(f"{tag} tslice={tslice_flip_tpbc} glb n flipped", n_flip)
        assert n_flip > 0.0
        q.flip_tpbc_with_tslice(prop, tslice_flip_tpbc)
        err = float(
            np.max(np.abs(np.asarray(prop) - arr0 * signs[:, None, None, None]))
        )
        q.json_results_append(
            f"{tag} tslice={tslice_flip_tpbc} sign flip max err = {err < check_eps}"
        )
        assert err == 0.0
        q.flip_tpbc_with_tslice(prop, tslice_flip_tpbc)
        assert np.array_equal(np.asarray(prop), arr0)
    # a negative tslice is a documented no-op
    q.flip_tpbc_with_tslice(prop, -1)
    assert np.array_equal(np.asarray(prop), arr0)
    q.json_results_append(f"{tag} tslice=-1 no-op")

q.begin_with_mpi(size_node_list)

geo = q.Geometry(total_site)
q.json_results_append(f"cqlat-propagator: geo.show()={geo.show()}")
rs = q.RngState("cqlat-propagator")
eye_wm = np.eye(12, dtype=np.complex128)

# --- random selections: a FieldSelection and a PointsSelection
fsel = q.FieldSelection()
fsel.set_rand(total_site, 16, rs.split("fsel"))
n_elems = int(fsel.n_elems)
n_elems_glb = float(q.glb_sum(float(n_elems)))
q.json_results_append("cqlat-propagator: fsel glb n_elems", n_elems_glb)
assert n_elems_glb > 0.0
psel_l = fsel.to_psel_local()
xg_f = np.array(psel_l.xg_arr, dtype=np.int64)

psel = q.PointsSelection()
psel.set_rand(total_site, 32, rs.split("psel"))
n_points = int(psel.n_points)
q.json_results_append(f"cqlat-propagator: psel n_points = {n_points}")
xg_p = np.array(psel.xg_arr, dtype=np.int64)

# --- set_rand_u1_src_fsel / set_rand_u1_sol_fsel
prop_src_f, fu1_f = q.mk_rand_u1_src(fsel, rs.split("u1-fsel"))
q.json_results_append(f"cqlat-propagator: fsel prop_src shape = {np.asarray(prop_src_f).shape}")
q.json_results_append(f"cqlat-propagator: fsel fu1 shape = {np.asarray(fu1_f).shape}")

arr_fu1_f = np.asarray(fu1_f)[:, 0]
loc_idx_f = np.array(
    [
        geo.index_from_coordinate(geo.coordinate_l_from_g(q.Coordinate(list(xg))))
        for xg in xg_f
    ],
    dtype=np.int64,
)
nz_f = np.nonzero(arr_fu1_f)[0]
assert np.array_equal(np.sort(nz_f), np.sort(loc_idx_f))
assert np.max(np.abs(np.abs(arr_fu1_f[loc_idx_f]) - 1.0)) < 1e-14
q.json_results_append(
    "cqlat-propagator: fsel fu1 phases", float(np.sum(np.abs(arr_fu1_f)))
)

s_src_f = q.SelProp(fsel)
s_src_f @= prop_src_f
s_fu1_f = q.SelectedField(q.ElemTypeComplexD, fsel, 1)
s_fu1_f @= fu1_f
src_f = np.asarray(s_src_f)
fu1_sel_f = np.asarray(s_fu1_f)
assert src_f.shape == (n_elems, 1, 12, 12)
err_src_f = float(
    np.max(np.abs(src_f - fu1_sel_f[:, :, None, None] * eye_wm[None, None, :, :]))
)
q.json_results_append(
    f"cqlat-propagator: fsel src = fu1 * 1 max err = {err_src_f < check_eps}"
)
assert err_src_f == 0.0

sol_f = q.get_rand_u1_sol(prop_src_f, fu1_f, fsel)
sol_arr_f = np.asarray(sol_f)
assert sol_arr_f.shape == (n_elems, 1, 12, 12)
err_sol_f = float(
    np.max(np.abs(sol_arr_f - src_f * np.conj(fu1_sel_f)[:, :, None, None]))
)
err_sol_id_f = float(np.max(np.abs(sol_arr_f - eye_wm[None, None, :, :])))
q.json_results_append(
    f"cqlat-propagator: fsel sol = src * conj(fu1) max err = {err_sol_f < check_eps}"
)
q.json_results_append(
    f"cqlat-propagator: fsel sol = 1 max err = {err_sol_id_f < check_eps}"
)
assert err_sol_f < 1e-14
assert err_sol_id_f < 1e-14
q.json_results_append(
    "cqlat-propagator: fsel sol sig",
    q.get_data_sig_arr(sol_arr_f, q.RngState("cqlat-propagator-fsel-sol"), 3),
    check_eps,
)

# --- set_rand_u1_src_psel / set_rand_u1_sol_psel
prop_src_p, fu1_p = q.mk_rand_u1_src(psel, rs.split("u1-psel"))
q.json_results_append(f"cqlat-propagator: psel prop_src shape = {np.asarray(prop_src_p).shape}")
arr_fu1_p = np.asarray(fu1_p)[:, 0]
# only the points owned by this rank are stored in the local fu1 field
loc_idx_p = np.array(
    [
        geo.index_from_coordinate(xl)
        for xg in xg_p
        for xl in [geo.coordinate_l_from_g(q.Coordinate(list(xg)))]
        if geo.is_local(xl)
    ],
    dtype=np.int64,
)
nz_p = np.nonzero(arr_fu1_p)[0]
assert np.array_equal(np.sort(nz_p), np.sort(loc_idx_p))
assert np.max(np.abs(np.abs(arr_fu1_p[loc_idx_p]) - 1.0)) < 1e-14
q.json_results_append(
    "cqlat-propagator: psel fu1 phases", float(np.sum(np.abs(arr_fu1_p)))
)

sp_src_p = q.PselProp(psel)
sp_src_p @= prop_src_p
sp_fu1_p = q.SelectedPoints(q.ElemTypeComplexD, psel, 1)
sp_fu1_p @= fu1_p
src_p = np.asarray(sp_src_p)
fu1_sel_p = np.asarray(sp_fu1_p)
assert src_p.shape == (n_points, 1, 12, 12)
err_src_p = float(
    np.max(np.abs(src_p - fu1_sel_p[:, :, None, None] * eye_wm[None, None, :, :]))
)
q.json_results_append(
    f"cqlat-propagator: psel src = fu1 * 1 max err = {err_src_p < check_eps}"
)
assert err_src_p == 0.0

sol_p = q.get_rand_u1_sol(prop_src_p, fu1_p, psel)
sol_arr_p = np.asarray(sol_p)
assert sol_arr_p.shape == (n_points, 1, 12, 12)
err_sol_p = float(
    np.max(np.abs(sol_arr_p - src_p * np.conj(fu1_sel_p)[:, :, None, None]))
)
err_sol_id_p = float(np.max(np.abs(sol_arr_p - eye_wm[None, None, :, :])))
q.json_results_append(
    f"cqlat-propagator: psel sol = src * conj(fu1) max err = {err_sol_p < check_eps}"
)
q.json_results_append(
    f"cqlat-propagator: psel sol = 1 max err = {err_sol_id_p < check_eps}"
)
assert err_sol_p < 1e-14
assert err_sol_id_p < 1e-14
q.json_results_append(
    "cqlat-propagator: psel sol sig",
    q.get_data_sig_arr(sol_arr_p, q.RngState("cqlat-propagator-psel-sol"), 3),
    check_eps,
)

# --- flip_tpbc_with_tslice_s_prop (SelProp)
s_prop = q.SelProp(fsel)
s_prop @= prop_src_f
s_flip0 = np.asarray(s_prop).copy()
check_flip_prop(s_prop, xg_f[:, 3], s_flip0, "cqlat-propagator: s_prop")

# --- flip_tpbc_with_tslice_sp_prop (PselProp)
sp_prop = q.PselProp(psel)
sp_prop @= prop_src_p
sp_flip0 = np.asarray(sp_prop).copy()
check_flip_prop(sp_prop, xg_p[:, 3], sp_flip0, "cqlat-propagator: sp_prop")

del prop_src_f, prop_src_p, fu1_f, fu1_p, s_prop, sp_prop, sol_f, sol_p
gc.collect()
q.json_results_append("cqlat-propagator: gc.collect() done")
q.json_results_append(f"cqlat-propagator: cwd={os.path.basename(os.getcwd())}")

q.timer_display()
if q.is_test():
    q.check_log_json(__file__)
q.end_with_mpi()
q.displayln_info("CHECK: finished successfully.")
