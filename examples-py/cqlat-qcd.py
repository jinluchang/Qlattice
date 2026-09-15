#!/usr/bin/env python3

# Tests for the cqlat QCD interface (qlat/cqlat/qcd.cpp):
#     gf_wilson_line_no_comm, gf_twist_boundary_at_boundary,
#     save_gauge_transform_cps, load_gauge_transform_cps
#
# They are reached through the public Python wrappers:
#     q.gf_wilson_line_no_comm(wlf, m, gf_ext, path, path_n=None)
#     q.gf_wilson_lines_no_comm(gf_ext, path_list)
#     q.gf_twist_boundary_at_boundary(gf, lmom, mu) / gf.twist_boundary_at_boundary
#     q.GaugeTransform.save_cps(path) / load_cps(path)
#     q.GaugeTransform.inv()   (the former cqlat ``gt_invert`` export is gone;
#                               inv() calls the Cython ``cc.gt_invert`` binding)
#
# ``gf_wilson_line_no_comm`` is checked against an independent numpy
# re-implementation of the path-ordered product of links (all paths used here
# stay inside the spatial directions that are not split across the two MPI
# ranks, so the re-implementation is purely local).
# ``gf_twist_boundary_at_boundary`` is checked link by link: only the links in
# direction ``mu`` on the ``mu``-boundary time slice may change, and the
# average plaquette must be unaffected.
# The CPS gauge-transform round trip is checked to be bit exact, and
# ``gt.inv()`` to produce the pointwise adjoint with ``gt * gt.inv() == 1``.

import gc
import os

import numpy as np

import qlat as q
import qlat.c as qc

check_eps = 1e-9

total_site = q.Coordinate([4, 4, 4, 8])

# The runner uses ``mpiexec -n 2`` so the [2, 1, 1, 1] layout is selected.
size_node_list = [
    [2, 1, 1, 1],
    [1, 1, 1, 1],
]

def ref_wilson_line_arr(geo, gf_arr, latt_size, path, path_n=None):
    # independent numpy re-implementation of qlat::gf_wilson_line_no_comm
    # (qcd-acc.h); only valid for paths that stay local
    if path_n is None:
        path_n = [1] * len(path)
    assert len(path_n) == len(path)
    rtn = np.zeros((geo.local_volume, 3, 3), dtype=np.complex128)
    for index in range(geo.local_volume):
        xg = [int(v) for v in geo.coordinate_g_from_l(geo.coordinate_from_index(index))]
        m = np.eye(3, dtype=np.complex128)
        for d, n in zip(path, path_n):
            for _ in range(n):
                if d >= 0:
                    xl = geo.coordinate_l_from_g(q.Coordinate(xg))
                    assert geo.is_local(xl)
                    m = m @ gf_arr[geo.index_from_coordinate(xl), d]
                    xg[d] = (xg[d] + 1) % latt_size[d]
                else:
                    mu = -d - 1
                    xg[mu] = (xg[mu] - 1) % latt_size[mu]
                    xl = geo.coordinate_l_from_g(q.Coordinate(xg))
                    assert geo.is_local(xl)
                    m = m @ np.conj(gf_arr[geo.index_from_coordinate(xl), mu].T)
        rtn[index] = m
    return rtn

def mk_boundary_mask(geo, mu, len_mu):
    rtn = np.zeros(geo.local_volume, dtype=bool)
    for index in range(geo.local_volume):
        xg = geo.coordinate_g_from_l(geo.coordinate_from_index(index))
        if xg[mu] == len_mu - 1:
            rtn[index] = True
    return rtn

q.begin_with_mpi(size_node_list)

geo = q.Geometry(total_site)
latt_size = tuple(int(geo.total_site[i]) for i in range(4))
q.json_results_append(f"cqlat-qcd: geo.show()={geo.show()}")
rs = q.RngState("cqlat-qcd")
eye_cm = np.eye(3, dtype=np.complex128)

gf = q.GaugeField(geo)
gf.set_rand(rs.split("gf"), 0.5, 10)
q.json_results_append("cqlat-qcd: gf plaq", gf.plaq(), check_eps)
q.json_results_append("cqlat-qcd: gf link_trace", gf.link_trace(), check_eps)
gf_arr = np.asarray(gf)
assert gf_arr.shape == (geo.local_volume, 4, 3, 3)

# --- gf_wilson_line_no_comm: closed spatial loops (x0 is the split direction)
path_a = [1, 2, -2, -3]
path_b = [2, 3, -3, -2]
path_n = [1, 2, -2, -3]
path_n_n = [3, 1, 3, 1]

wlf_a = q.FieldColorMatrix(geo, 1)
q.gf_wilson_line_no_comm(wlf_a, 0, gf, path_a)
ref_a = ref_wilson_line_arr(geo, gf_arr, latt_size, path_a)
err_a = float(np.max(np.abs(np.asarray(wlf_a)[:, 0] - ref_a)))
q.json_results_append(
    f"cqlat-qcd: gf_wilson_line_no_comm path_a max err = {err_a < check_eps}"
)
assert err_a < 1e-14
q.json_results_append(
    "cqlat-qcd: gf_wilson_line_no_comm path_a sum",
    float(np.sum(np.asarray(wlf_a)[:, 0])),
    check_eps,
)

wlf_b = q.FieldColorMatrix(geo, 1)
q.gf_wilson_line_no_comm(wlf_b, 0, gf, path_b)
ref_b = ref_wilson_line_arr(geo, gf_arr, latt_size, path_b)
err_b = float(np.max(np.abs(np.asarray(wlf_b)[:, 0] - ref_b)))
q.json_results_append(
    f"cqlat-qcd: gf_wilson_line_no_comm path_b max err = {err_b < check_eps}"
)
assert err_b < 1e-14

wlf_n = q.FieldColorMatrix(geo, 1)
q.gf_wilson_line_no_comm(wlf_n, 0, gf, path_n, path_n_n)
ref_n = ref_wilson_line_arr(geo, gf_arr, latt_size, path_n, path_n_n)
err_n = float(np.max(np.abs(np.asarray(wlf_n)[:, 0] - ref_n)))
q.json_results_append(
    f"cqlat-qcd: gf_wilson_line_no_comm path_n max err = {err_n < check_eps}"
)
assert err_n < 1e-14

# --- gf_wilson_lines_no_comm: the multi-path wrapper must agree with the
#     single-path calls above
wlf_all = q.gf_wilson_lines_no_comm(gf, [path_a, path_b, (path_n, path_n_n)])
wlf_all_arr = np.asarray(wlf_all)
assert wlf_all_arr.shape == (geo.local_volume, 3, 3, 3)
err_all = max(
    float(np.max(np.abs(wlf_all_arr[:, 0] - np.asarray(wlf_a)[:, 0]))),
    float(np.max(np.abs(wlf_all_arr[:, 1] - np.asarray(wlf_b)[:, 0]))),
    float(np.max(np.abs(wlf_all_arr[:, 2] - np.asarray(wlf_n)[:, 0]))),
)
q.json_results_append(
    f"cqlat-qcd: gf_wilson_lines_no_comm vs single max err = {err_all < check_eps}"
)
assert err_all == 0.0
q.json_results_append(
    "cqlat-qcd: gf_wilson_lines_no_comm sig",
    q.get_data_sig_arr(wlf_all, q.RngState("cqlat-qcd-wlf"), 3),
    check_eps,
)

# --- gf_twist_boundary_at_boundary: only the boundary links in direction mu
#     change, and the plaquette is unchanged
for mu, lmom in [(3, -0.5), (0, 0.25)]:
    gf_t = gf.copy()
    q.gf_twist_boundary_at_boundary(gf_t, lmom, mu)
    gf_t_m = gf.copy()
    gf_t_m.twist_boundary_at_boundary(lmom, mu)
    assert np.array_equal(np.asarray(gf_t), np.asarray(gf_t_m))
    gf_t_d = gf.copy()
    qc.gf_twist_boundary_at_boundary(gf_t_d, lmom, mu)
    assert np.array_equal(np.asarray(gf_t), np.asarray(gf_t_d))
    arr_base = np.asarray(gf)
    arr_t = np.asarray(gf_t)
    phase = np.exp(1j * 2.0 * np.pi * lmom)
    mask = mk_boundary_mask(geo, mu, latt_size[mu])
    n_mask = float(q.glb_sum(float(np.sum(mask))))
    q.json_results_append(
        f"cqlat-qcd: gf_twist_boundary_at_boundary mu={mu} glb n links", n_mask
    )
    assert n_mask > 0.0
    exp_t = arr_base.copy()
    exp_t[mask, mu] *= phase
    err_t = float(np.max(np.abs(arr_t - exp_t)))
    q.json_results_append(
        f"cqlat-qcd: gf_twist_boundary_at_boundary mu={mu} max err"
        f" = {err_t < check_eps}"
    )
    assert err_t < 1e-14
    q.json_results_append(
        f"cqlat-qcd: gf_twist_boundary_at_boundary mu={mu} plaq", gf_t.plaq(), check_eps
    )
    err_plaq = abs(gf_t.plaq() - gf.plaq())
    q.json_results_append(
        f"cqlat-qcd: gf_twist_boundary_at_boundary mu={mu} plaq diff"
        f" = {err_plaq < check_eps}"
    )
    assert err_plaq < 1e-12
    # the number of changed links equals the number of boundary links
    n_changed = int(np.sum(np.any(arr_t != arr_base, axis=(2, 3))))
    n_changed = int(q.glb_sum(float(n_changed)))
    q.json_results_append(
        f"cqlat-qcd: gf_twist_boundary_at_boundary mu={mu} n changed = {n_changed}"
    )
    assert float(n_changed) == n_mask

# --- GaugeTransform: set_rand + unitarize, CPS round trip, inversion
gt = q.GaugeTransform(geo)
gt.set_rand(rs.split("gt"), 0.5, 10)
gt.unitarize()
gt_arr = np.asarray(gt)
assert gt_arr.shape == (geo.local_volume, 1, 3, 3)
err_unitary = float(
    np.max(
        np.abs(gt_arr @ np.conj(gt_arr.transpose(0, 1, 3, 2)))
        - eye_cm[None, None, :, :]
    )
)
q.json_results_append(f"cqlat-qcd: gt unitary max err = {err_unitary < check_eps}")
q.json_results_append(
    "cqlat-qcd: gt sig",
    q.get_data_sig_arr(gt, q.RngState("cqlat-qcd-gt"), 3),
    check_eps,
)

path_cps = "cqlat-qcd-gt-cps.tmp"
n_bytes = int(gt.save_cps(path_cps))
gt_load = q.GaugeTransform()
n_bytes_load = int(gt_load.load_cps(path_cps))
q.json_results_append(f"cqlat-qcd: save_gauge_transform_cps bytes = {n_bytes}")
q.json_results_append(f"cqlat-qcd: load_gauge_transform_cps bytes = {n_bytes_load}")
q.json_results_append(
    f"cqlat-qcd: save_gauge_transform_cps file size = {os.path.getsize(path_cps)}"
)
assert n_bytes == n_bytes_load
assert n_bytes > 0
assert np.array_equal(np.asarray(gt), np.asarray(gt_load))
q.json_results_append(
    "cqlat-qcd: cps round trip sig",
    q.get_data_sig_arr(gt_load, q.RngState("cqlat-qcd-gt-load"), 3),
    check_eps,
)

gt_inv = gt.inv()
gt_inv_arr = np.asarray(gt_inv)
err_inv = float(np.max(np.abs(gt_inv_arr - np.conj(gt_arr.transpose(0, 1, 3, 2)))))
q.json_results_append(f"cqlat-qcd: gt_invert vs adjoint max err = {err_inv < check_eps}")
assert err_inv == 0.0
q.json_results_append(
    "cqlat-qcd: gt_invert sig",
    q.get_data_sig_arr(gt_inv, q.RngState("cqlat-qcd-gt-inv"), 3),
    check_eps,
)

gt_prod = gt * gt_inv
err_prod = float(np.max(np.abs(np.asarray(gt_prod) - eye_cm[None, None, :, :])))
q.json_results_append(f"cqlat-qcd: gt * gt.inv() vs 1 max err = {err_prod < check_eps}")
assert err_prod < 1e-14
gt_prod2 = gt_inv * gt
err_prod2 = float(np.max(np.abs(np.asarray(gt_prod2) - eye_cm[None, None, :, :])))
q.json_results_append(f"cqlat-qcd: gt.inv() * gt vs 1 max err = {err_prod2 < check_eps}")
assert err_prod2 < 1e-14

del gf, gf_t, gf_t_m, gf_t_d, gt, gt_load, gt_inv, gt_prod, gt_prod2
del wlf_a, wlf_b, wlf_n, wlf_all
gc.collect()
q.json_results_append("cqlat-qcd: gc.collect() done")
q.json_results_append(f"cqlat-qcd: cwd={os.path.basename(os.getcwd())}")

q.timer_display()
if q.is_test():
    q.check_log_json(__file__)
q.end_with_mpi()
q.displayln_info("CHECK: finished successfully.")
