#!/usr/bin/env python3

# Tests for the cqlat contraction interface:
#     qlat/cqlat/contraction-field.cpp: contract_chvp_16_field
#     qlat/cqlat/contraction-hvp.cpp:   contract_chvp3_sfield
#     qlat/cqlat/contraction-pion.cpp:  contract_pion_sfield
#
# They are reached through the Python wrappers
#     q.contract_chvp_16       -> c.contract_chvp_16_field
#     q.contract_chvp3_field   -> c.contract_chvp3_sfield
#     q.contract_pion_field    -> cc.contract_pion (dense Prop) or
#                                 c.contract_pion_sfield (SelProp)
# and the cqlat exports are also invoked directly as ``qc.<name>``.
#
# The checks are analytic: each result is compared against a numpy
# re-implementation of the corresponding C++ header.
#   * contract_chvp_16: chvp(x, mu*4+nu) ==
#       tr(g5_herm(prop2(x)) * gamma[mu] * prop1(x) * gamma[nu])
#     (``q.get_gamma_matrix`` returns CPS's convention gamma matrices, which
#     are exactly ``SpinMatrixConstants::get_cps_gammas()``; the 12x12
#     WilsonMatrix structure is ``gamma kron unit_color``).
#   * contract_chvp3_field: sums the same traces for mu = 0, 1, 2 over the
#     selected sites of each time slice, divides by the field-selection
#     probability, and does a global sum.
#   * contract_pion_field: sums qnorm(prop(x)) over the sites of each time
#     slice (the dense and the selected propagator paths).

import gc

import numpy as np

import qlat as q
import qlat.c as qc

check_eps = 1e-10

def rand_complex_arr(rs, tag, shape):
    """Global random complex array; differs between the values of the two MPI
    ranks because it is indexed by the *global* coordinate."""
    return (rs.split(tag + "-re").u_rand_arr(shape) - 0.5) + 1j * (
        rs.split(tag + "-im").u_rand_arr(shape) - 0.5
    )

def mk_gammas():
    """CPS convention gamma matrices as 12x12 (spin kron color) matrices."""
    eye3 = np.eye(3, dtype=np.complex128)
    G5 = np.kron(np.asarray(q.get_gamma_matrix(5), dtype=np.complex128), eye3)
    Gs = [
        np.kron(np.asarray(q.get_gamma_matrix(mu), dtype=np.complex128), eye3)
        for mu in range(4)
    ]
    return G5, Gs

def g5_herm_np(m, G5):
    """numpy version of ``g5_herm``: gamma5 * adjoint(m) * gamma5"""
    return G5 @ m.conj().swapaxes(-1, -2) @ G5

def fill_prop_from_global(prop, geo, arr_g):
    """copy the local part of a global (spin-color matrix valued) array"""
    arr = np.asarray(prop)
    for index in range(geo.local_volume):
        xg = geo.coordinate_g_from_l(geo.coordinate_from_index(index))
        arr[index, 0] = arr_g[xg[0], xg[1], xg[2], xg[3]]

def max_diff(a, b):
    return float(np.max(np.abs(np.asarray(a) - np.asarray(b))))

size_node_list = [
    [2, 1, 1, 1],
    [1, 1, 1, 1],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
q.json_results_append(f"cqlat-contraction: geo.show()={geo.show()}")
latt_size = tuple(int(geo.total_site[i]) for i in range(4))
t_size = latt_size[3]
t_slice_src = 2

rs = q.RngState("cqlat-contraction")

# --- rank dependent propagator data (global arrays, local parts copied in)
arr1_g = rand_complex_arr(rs, "prop1", latt_size + (12, 12))
arr2_g = rand_complex_arr(rs, "prop2", latt_size + (12, 12))

prop1 = q.Prop(geo)
fill_prop_from_global(prop1, geo, arr1_g)
prop2 = q.Prop(geo)
fill_prop_from_global(prop2, geo, arr2_g)

G5, Gs = mk_gammas()
prop1_arr = np.asarray(prop1)[:, 0]
prop2_arr = np.asarray(prop2)[:, 0]
local_volume = geo.local_volume

# --- contract_chvp_16_field (via q.contract_chvp_16)
chvp = q.contract_chvp_16(prop1, prop2)
chvp_arr = np.asarray(chvp).reshape(local_volume, 16)
assert np.asarray(chvp).shape == (local_volume, 16)
assert chvp.multiplicity == 16
ref_chvp = np.zeros((local_volume, 16), dtype=np.complex128)
g5h_prop2 = g5_herm_np(prop2_arr, G5)
for mu in range(4):
    for nu in range(4):
        wm = g5h_prop2 @ Gs[mu] @ prop1_arr @ Gs[nu]
        ref_chvp[:, mu * 4 + nu] = np.trace(wm, axis1=-2, axis2=-1)
err_chvp = float(np.max(np.abs(chvp_arr - ref_chvp)))
assert err_chvp < check_eps, err_chvp
q.json_results_append(f"cqlat-contraction: contract_chvp_16 vs numpy = {err_chvp < check_eps}")
q.json_results_append(
    "cqlat-contraction: contract_chvp_16 max error", err_chvp, check_eps
)
q.json_results_append(
    "cqlat-contraction: contract_chvp_16 sum(mu=0,nu=0)",
    q.glb_sum(complex(np.sum(chvp_arr[:, 0]))),
    check_eps,
)

# --- contract_chvp_16_field invoked directly as a cqlat export
chvp_direct = q.FieldComplexD(geo, 16)
qc.contract_chvp_16_field(chvp_direct, prop1, prop2)
err_chvp_direct = max_diff(chvp_direct, chvp)
assert err_chvp_direct == 0.0, err_chvp_direct
q.json_results_append(
    f"cqlat-contraction: c.contract_chvp_16_field direct vs wrapper = {err_chvp_direct == 0.0}"
)

# --- contract_chvp_16 is linear in each propagator
prop1x2 = q.Prop(geo)
fill_prop_from_global(prop1x2, geo, 2.0 * arr1_g)
chvp_x2 = np.asarray(q.contract_chvp_16(prop1x2, prop2)).reshape(local_volume, 16)
err_lin = float(np.max(np.abs(chvp_x2 - 2.0 * chvp_arr)))
assert err_lin < check_eps, err_lin
q.json_results_append(
    "cqlat-contraction: contract_chvp_16 linearity", err_lin, check_eps
)

# --- field selection and selected propagators
fsel = q.FieldSelection()
fsel.set_rand(geo.total_site, 16, rs.split("fsel"))
fsel_prob = q.glb_sum(fsel.n_elems) / geo.total_volume
assert fsel_prob > 0.0
q.json_results_append("cqlat-contraction: fsel prob", float(fsel_prob), check_eps)

sp1 = q.SelProp(fsel)
sp1 @= prop1
sp2 = q.SelProp(fsel)
sp2 @= prop2

sp1_arr = np.asarray(sp1)[:, 0]
sp2_arr = np.asarray(sp2)[:, 0]

# --- contract_chvp3_sfield (via q.contract_chvp3_field)
ld_hvp = q.contract_chvp3_field(sp1, sp2, t_slice_src)
hvp_arr = np.asarray(ld_hvp)
assert hvp_arr.shape == (t_size, 3)
ref_hvp = np.zeros((t_size, 3), dtype=np.complex128)
g5h_sp2 = g5_herm_np(sp2_arr, G5)
for idx in range(fsel.n_elems):
    xg = fsel.coordinate_from_idx(idx)
    tsep = (xg[3] - t_slice_src) % t_size
    for mu in range(3):
        wm = g5h_sp2[idx] @ Gs[mu] @ sp1_arr[idx] @ Gs[mu]
        ref_hvp[tsep, mu] += np.trace(wm)
ref_hvp = q.glb_sum(ref_hvp) / fsel_prob
err_hvp = float(np.max(np.abs(hvp_arr - ref_hvp)))
assert err_hvp < check_eps, err_hvp
q.json_results_append(f"cqlat-contraction: contract_chvp3_field vs numpy = {err_hvp < check_eps}")
q.json_results_append(
    "cqlat-contraction: contract_chvp3_field max error", err_hvp, check_eps
)
q.json_results_append(
    "cqlat-contraction: contract_chvp3_field tsep=0 mu=0",
    float(np.real(hvp_arr[0, 0])),
    check_eps,
)

# --- contract_chvp3_sfield invoked directly as a cqlat export
ld_hvp_direct = q.LatData()
qc.contract_chvp3_sfield(ld_hvp_direct, sp1, sp2, t_slice_src)
err_hvp_direct = float(np.max(np.abs(np.asarray(ld_hvp_direct) - hvp_arr)))
assert err_hvp_direct == 0.0, err_hvp_direct
q.json_results_append(
    f"cqlat-contraction: c.contract_chvp3_sfield direct vs wrapper = {err_hvp_direct == 0.0}"
)

# --- contract_pion_field, dense Prop path (cc.contract_pion)
ld_pion_prop = q.contract_pion_field(prop1, t_slice_src)
pion_prop_arr = np.asarray(ld_pion_prop)
assert pion_prop_arr.shape == (t_size,)
ref_pion_prop = np.zeros(t_size, dtype=np.complex128)
for index in range(local_volume):
    xg = geo.coordinate_g_from_l(geo.coordinate_from_index(index))
    tsep = (xg[3] - t_slice_src) % t_size
    ref_pion_prop[tsep] += np.sum(np.abs(prop1_arr[index]) ** 2)
ref_pion_prop = q.glb_sum(ref_pion_prop)
err_pion_prop = float(np.max(np.abs(pion_prop_arr - ref_pion_prop)))
assert err_pion_prop < check_eps, err_pion_prop
q.json_results_append(
    f"cqlat-contraction: contract_pion_field(Prop) vs numpy = {err_pion_prop < check_eps}"
)
q.json_results_append(
    "cqlat-contraction: contract_pion_field(Prop) max error", err_pion_prop, check_eps
)
q.json_results_append(
    "cqlat-contraction: contract_pion_field(Prop) tsep=0",
    float(np.real(pion_prop_arr[0])),
    check_eps,
)

# --- contract_pion_field, SelProp path (c.contract_pion_sfield)
ld_pion_sel = q.contract_pion_field(sp1, t_slice_src)
pion_sel_arr = np.asarray(ld_pion_sel)
assert pion_sel_arr.shape == (t_size,)
ref_pion_sel = np.zeros(t_size, dtype=np.complex128)
for idx in range(fsel.n_elems):
    xg = fsel.coordinate_from_idx(idx)
    tsep = (xg[3] - t_slice_src) % t_size
    ref_pion_sel[tsep] += np.sum(np.abs(sp1_arr[idx]) ** 2)
ref_pion_sel = q.glb_sum(ref_pion_sel) / fsel_prob
err_pion_sel = float(np.max(np.abs(pion_sel_arr - ref_pion_sel)))
assert err_pion_sel < check_eps, err_pion_sel
q.json_results_append(
    f"cqlat-contraction: contract_pion_field(SelProp) vs numpy = {err_pion_sel < check_eps}"
)
q.json_results_append(
    "cqlat-contraction: contract_pion_field(SelProp) max error", err_pion_sel, check_eps
)
q.json_results_append(
    "cqlat-contraction: contract_pion_field(SelProp) tsep=0",
    float(np.real(pion_sel_arr[0])),
    check_eps,
)

# --- contract_pion_sfield invoked directly as a cqlat export
ld_pion_sel_direct = q.LatData()
qc.contract_pion_sfield(ld_pion_sel_direct, sp1, t_slice_src)
err_pion_sel_direct = float(
    np.max(np.abs(np.asarray(ld_pion_sel_direct) - pion_sel_arr))
)
assert err_pion_sel_direct == 0.0, err_pion_sel_direct
q.json_results_append(
    f"cqlat-contraction: c.contract_pion_sfield direct vs wrapper = {err_pion_sel_direct == 0.0}"
)

del prop1, prop2, prop1x2, sp1, sp2, fsel
gc.collect()

q.timer_display()
if q.is_test():
    q.check_log_json(__file__)
q.end_with_mpi()
q.displayln_info("CHECK: finished successfully.")
