#!/usr/bin/env python3

# Tests for the cqlat ftHMC interface:
#     qlat/cqlat/fthmc.cpp:     add_flow_flow_info, free_flow_info,
#                               set_gm_force_flowed_no_det
#     qlat/cqlat/hmc-stats.cpp: get_gm_force_magnitudes
#
# They are reached through
#     q.FlowInfo().add_flow(...)         -> c.add_flow_flow_info
#     q.FlowInfo.__del__ (via gc.collect) -> c.free_flow_info
#     q.c.set_gm_force_flowed_no_det      (no Python wrapper exists)
#     q.get_gm_force_magnitudes           (qlat.hmc_stats wrapper)
#
# The checks are:
#   * an empty FlowInfo propagates the force unchanged,
#   * a trivial flow (epsilon = 0) reproduces both the unflowed gauge field
#     (``q.gf_flow``) and the unflowed force,
#   * a genuine flow gives a finite, anti-Hermitian force that changes the
#     gauge field,
#   * ``get_gm_force_magnitudes`` returns the norm family
#     (index 0 = mean magnitude, index 1 = L2 norm, ..., index n_elems-1 =
#     global maximum magnitude, see hmc-stats.h) and the values are compared
#     against a numpy reference built from the global random basis
#     coefficients (``Tr[T_a T_b] = -2 delta_ab`` implies
#     ``-Tr[g^2] = 2 sum_a b_a^2`` for the anti-Hermitian force g).

import gc

import numpy as np

import qlat as q
import qlat.c as qc

check_eps = 1e-10

def rand_real_arr(rs, tag, shape):
    """Global random array; differs between the two MPI ranks because it is
    indexed by the *global* coordinate."""
    return rs.split(tag).u_rand_arr(shape) - 0.5

def rand_complex_arr(rs, tag, shape):
    return rand_real_arr(rs, tag + "-re", shape) + 1j * rand_real_arr(
        rs, tag + "-im", shape
    )

def fill_field_from_global(field, geo, arr_g):
    """copy the local part of a global array into a lattice field"""
    arr = np.asarray(field)
    for index in range(geo.local_volume):
        xg = geo.coordinate_g_from_l(geo.coordinate_from_index(index))
        arr[index] = arr_g[xg[0], xg[1], xg[2], xg[3]]

def max_diff(a, b):
    return float(np.max(np.abs(np.asarray(a) - np.asarray(b))))

def gm_local_magnitudes(gm):
    """sqrt(-Tr[g^2]) for the anti-Hermitian force g at every local site/dir"""
    arr = np.asarray(gm)
    return np.sqrt(np.real(-np.trace(arr @ arr, axis1=-2, axis2=-1)))

size_node_list = [
    [2, 1, 1, 1],
    [1, 1, 1, 1],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
q.json_results_append(f"cqlat-fthmc: geo.show()={geo.show()}")
latt_size = tuple(int(geo.total_site[i]) for i in range(4))
total_volume = int(np.prod(latt_size))

rs = q.RngState("cqlat-fthmc")

# --- a random gauge field (rank dependent), projected onto U(3)
gf_arr_g = rand_complex_arr(rs, "gf", latt_size + (4, 3, 3))
gf0 = q.GaugeField(geo)
fill_field_from_global(gf0, geo, gf_arr_g)
gf0.unitarize()
q.json_results_append("cqlat-fthmc: gf0 avg plaq", q.gf_avg_plaq(gf0), check_eps)

# --- a random anti-Hermitian gauge momentum built from a global basis array
basis_arr_g = rand_real_arr(rs, "basis", latt_size + (4 * 8,))
basis = q.FieldRealD(geo, 4 * 8)
fill_field_from_global(basis, geo, basis_arr_g)
gm_pre = q.GaugeMomentum(geo)
q.set_anti_hermitian_matrix_from_basis(gm_pre, basis)
gm_pre_arr = np.asarray(gm_pre)
assert np.max(np.abs(gm_pre_arr + gm_pre_arr.conj().swapaxes(-1, -2))) == 0.0
q.json_results_append("cqlat-fthmc: gm_pre qnorm", q.qnorm(gm_pre), check_eps)

# --- zero and trivial flows built with add_flow_flow_info
fi_empty = q.FlowInfo()
assert fi_empty.show() == ""
q.json_results_append("cqlat-fthmc: empty FlowInfo step count")

fi_trivial = q.FlowInfo()
for mu in range(4):
    for eo in [1, 2]:
        fi_trivial.add_flow(eo, mu, 0.0, 1)
n_steps_trivial = len(fi_trivial.show().strip().split("\n"))
assert n_steps_trivial == 8
q.json_results_append(f"cqlat-fthmc: trivial FlowInfo step count = {n_steps_trivial}")

# --- set_gm_force_flowed_no_det: an empty flow is the identity
gm_empty = q.GaugeMomentum(geo)
qc.set_gm_force_flowed_no_det(gm_empty, gm_pre, gf0, fi_empty)
err_empty = max_diff(gm_empty, gm_pre)
assert err_empty == 0.0, err_empty
q.json_results_append(
    f"cqlat-fthmc: set_gm_force_flowed_no_det (empty flow) = {err_empty < check_eps}"
)

# --- set_gm_force_flowed_no_det: a trivial (epsilon = 0) flow is the identity
gm_trivial = q.GaugeMomentum(geo)
qc.set_gm_force_flowed_no_det(gm_trivial, gm_pre, gf0, fi_trivial)
err_trivial = max_diff(gm_trivial, gm_pre)
assert err_trivial < check_eps, err_trivial
q.json_results_append(
    f"cqlat-fthmc: set_gm_force_flowed_no_det (trivial flow)"
    f" = {err_trivial < check_eps}"
)

# --- q.gf_flow with the trivial flow leaves the gauge field unchanged
gf_trivial = q.GaugeField(geo)
q.gf_flow(gf_trivial, gf0, fi_trivial)
err_gf_trivial = max_diff(gf_trivial, gf0)
assert err_gf_trivial == 0.0, err_gf_trivial
q.json_results_append(f"cqlat-fthmc: gf_flow (trivial flow) = {err_gf_trivial < check_eps}")

# --- a genuine flow
flow_eps = 0.05
fi_flow = q.FlowInfo()
for mu in range(4):
    for eo in [1, 2]:
        fi_flow.add_flow(eo, mu, flow_eps, 1)
n_steps_flow = len(fi_flow.show().strip().split("\n"))
assert n_steps_flow == 8
q.json_results_append(f"cqlat-fthmc: flow FlowInfo step count = {n_steps_flow}")

gm_flow = q.GaugeMomentum(geo)
qc.set_gm_force_flowed_no_det(gm_flow, gm_pre, gf0, fi_flow)
gm_flow_arr = np.asarray(gm_flow)
assert np.all(np.isfinite(gm_flow_arr))
err_antiherm = float(np.max(np.abs(gm_flow_arr + gm_flow_arr.conj().swapaxes(-1, -2))))
assert err_antiherm < check_eps, err_antiherm
q.json_results_append(
    f"cqlat-fthmc: flowed force is anti-Hermitian = {err_antiherm < check_eps}"
)
q.json_results_append(f"cqlat-fthmc: flowed force is finite = {np.all(np.isfinite(gm_flow_arr))}")
q.json_results_append(
    f"cqlat-fthmc: flowed force differs from input = {max_diff(gm_flow, gm_pre) > 0.0}"
)

gf_flowed = q.GaugeField(geo)
q.gf_flow(gf_flowed, gf0, fi_flow)
q.json_results_append(
    f"cqlat-fthmc: gf_flow changes the gauge field = {max_diff(gf_flowed, gf0) > 0.0}"
)
q.json_results_append(
    "cqlat-fthmc: flowed gf avg plaq", q.gf_avg_plaq(gf_flowed), check_eps
)

# --- get_gm_force_magnitudes (via the qlat.hmc_stats wrapper)
n_elems = 5
mag_pre = q.get_gm_force_magnitudes(gm_pre, n_elems)
assert len(mag_pre) == n_elems
assert all(float(v) >= 0.0 for v in mag_pre)
q.json_results_append(
    f"cqlat-fthmc: get_gm_force_magnitudes all non-negative = "
    f"{all(float(v) >= 0.0 for v in mag_pre)}"
)
q.json_results_append(
    "cqlat-fthmc: get_gm_force_magnitudes(gm_pre)",
    np.array(mag_pre, dtype=np.float64),
    check_eps,
)

# --- n_elems < 3 must be rejected cleanly ---------------------------------
#     (n_elems == 2 used to pass qassert(n_elems >= 2) and then index a
#     multiplicity-1 vector inside an OpenMP loop -> uncatchable SIGABRT)
try:
    q.get_gm_force_magnitudes(gm_pre, 2)
except RuntimeError:
    n2_rejected = True
else:
    n2_rejected = False
assert n2_rejected
q.json_results_append(
    f"cqlat-fthmc: get_gm_force_magnitudes n_elems=2 rejected = {n2_rejected}"
)

# --- independent reference from the global basis coefficients:
#     gm(x, mu) = sum_a basis[x, mu * 8 + a] * T_a with Tr[T_a T_b] = -2 delta
#     => magnitude = sqrt(2 * sum_a basis^2)
basis_site = basis_arr_g.reshape(total_volume, 4, 8)
ref_mag = np.sqrt(2.0 * np.sum(basis_site**2, axis=-1))
ref_mean = float(ref_mag.mean())
ref_l2 = float(np.sqrt(np.mean(ref_mag**2)))
ref_max = float(ref_mag.max())
assert abs(float(mag_pre[0]) - ref_mean) < check_eps * max(1.0, ref_mean)
assert abs(float(mag_pre[1]) - ref_l2) < check_eps * max(1.0, ref_l2)
assert abs(float(mag_pre[n_elems - 1]) - ref_max) < check_eps * max(1.0, ref_max)
q.json_results_append(
    f"cqlat-fthmc: magnitude 0 vs numpy mean"
    f" = {abs(float(mag_pre[0]) - ref_mean) < check_eps}"
)
q.json_results_append(
    f"cqlat-fthmc: magnitude 1 vs numpy l2"
    f" = {abs(float(mag_pre[1]) - ref_l2) < check_eps}"
)
q.json_results_append(
    f"cqlat-fthmc: magnitude n_elems-1 vs numpy max"
    f" = {abs(float(mag_pre[n_elems - 1]) - ref_max) < check_eps}"
)
q.json_results_append(
    f"cqlat-fthmc: magnitude n_elems-1 equals largest magnitude = "
    f"{abs(float(mag_pre[n_elems - 1]) - ref_max) < check_eps * max(1.0, ref_max)}"
)

# --- the cqlat export directly
mag_pre_direct = qc.get_gm_force_magnitudes(gm_pre, n_elems)
assert len(mag_pre_direct) == n_elems
q.json_results_append(
    f"cqlat-fthmc: c.get_gm_force_magnitudes direct max diff"
    f" = {float(np.max(np.abs(np.array(mag_pre_direct) - np.array(mag_pre)))) < check_eps}"
)

# --- magnitudes of the flowed force (checked against the field itself)
mag_flow = q.get_gm_force_magnitudes(gm_flow, n_elems)
assert len(mag_flow) == n_elems
assert all(float(v) >= 0.0 for v in mag_flow)
mag_flow_local = gm_local_magnitudes(gm_flow)
mean_flow_ref = q.glb_sum(float(np.sum(mag_flow_local))) / (4.0 * total_volume)
max_flow_ref_local = float(np.max(mag_flow_local))
assert abs(float(mag_flow[0]) - mean_flow_ref) < check_eps * max(1.0, mean_flow_ref)
assert float(mag_flow[n_elems - 1]) >= max_flow_ref_local
assert float(mag_flow[n_elems - 1]) <= q.glb_sum(max_flow_ref_local) + check_eps
q.json_results_append(
    "cqlat-fthmc: get_gm_force_magnitudes(gm_flow)",
    np.array(mag_flow, dtype=np.float64),
    check_eps,
)
q.json_results_append(
    f"cqlat-fthmc: flowed magnitude 0 vs numpy mean"
    f" = {abs(float(mag_flow[0]) - mean_flow_ref) < check_eps}"
)

# --- free_flow_info is exercised by FlowInfo.__del__ (explicit gc.collect())
fi_tmp = q.FlowInfo()
fi_tmp.add_flow(1, 0, 0.1, 1)
assert len(fi_tmp.show().strip().split("\n")) == 1
del fi_tmp
gc.collect()
q.json_results_append("cqlat-fthmc: free_flow_info (FlowInfo.__del__)")

del fi_empty, fi_trivial, fi_flow, basis, gm_pre, gm_empty, gm_trivial, gm_flow
gc.collect()

q.timer_display()
if q.is_test():
    q.check_log_json(__file__)
q.end_with_mpi()
q.displayln_info("CHECK: finished successfully.")
