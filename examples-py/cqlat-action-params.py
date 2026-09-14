#!/usr/bin/env python3

# Tests for the cqlat ScalarAction and FermionAction interfaces:
#   qlat/cqlat/scalar-action.cpp:
#     free_scalar_action, set_scalar_action, get_m_sq_scalar_action,
#     get_lmbd_scalar_action, get_alpha_scalar_action,
#     hmc_estimate_mass_scalar_action, to_mass_factor_scalar_action
#   qlat/cqlat/fermion-action.cpp:
#     free_fermion_action, set_fermion_action, get_ls_fermion_action,
#     get_omega_fermion_action, get_mobius_scale_fermion_action
#
# These are reached through the ``q.ScalarAction`` (``qlat.scalar_action``) and
# ``q.FermionAction`` (``qlat.fermion_action``) Python classes.
# ``hmc_estimate_mass_scalar_action`` and ``to_mass_factor_scalar_action`` are
# checked against independent numpy re-implementations of
# ``qlat/qlat/include/qlat/scalar-action.h``.

import numpy as np

import qlat as q
import qlat.c as qc

check_eps = 1e-10

total_site = q.Coordinate([4, 4, 4, 8])

size_node_list = [
    [2, 1, 1, 1],
    [1, 1, 1, 1],
]

# --- ScalarAction parameters
sa_args = dict(m_sq=0.7, lmbd=1.1, alpha=0.03)

# --- FermionAction parameters
fa_args = dict(mass=0.12, ls=6, m5=1.3, mobius_scale=1.7)
zm_omega = [
    complex(0.80, 0.30),
    complex(1.10, -0.20),
    complex(0.95, 0.45),
]

def ref_to_mass_factor(v):
    # ScalarAction::to_mass_factor
    return (1.0 + 2.0 * np.arcsin(v) / np.pi) ** 2

def ref_hmc_estimate_mass(fld, frc, orig_loc, phi0):
    # ScalarAction::hmc_estimate_mass; the phi0 shift applies at the global
    # origin for the m == 0 component only (orig_loc marks the local origin)
    den = np.abs(fld) ** 2
    den[orig_loc, 0] = np.abs(fld[orig_loc, 0] - phi0) ** 2
    return 4.0 / np.pi**2 * np.sqrt(np.abs(frc) ** 2 / den)

def mk_complex_field_from_global(geo, g_re, g_im):
    f = q.FieldComplexD(geo, g_re.shape[-1])
    buf = np.asarray(f)
    for index in range(geo.local_volume):
        xg = geo.coordinate_g_from_l(geo.coordinate_from_index(index))
        buf[index, :] = (
            g_re[xg[0], xg[1], xg[2], xg[3], :]
            + 1j * g_im[xg[0], xg[1], xg[2], xg[3], :]
        )
    return f

q.begin_with_mpi(size_node_list)

geo = q.Geometry(total_site)
q.json_results_append(f"cqlat-action-params: geo.show()={geo.show()}")
latt_size = tuple(int(geo.total_site[i]) for i in range(4))
n_site = int(np.prod(latt_size))
orig_loc = np.zeros(geo.local_volume, dtype=bool)
for index in range(geo.local_volume):
    if geo.coordinate_g_from_l(geo.coordinate_from_index(index)) == q.Coordinate(
        [0, 0, 0, 0]
    ):
        orig_loc[index] = True

# =====================================================================
# ScalarAction
# =====================================================================
sa = q.ScalarAction(sa_args["m_sq"], sa_args["lmbd"], sa_args["alpha"])
assert sa.m_sq() == sa_args["m_sq"]
assert sa.lmbd() == sa_args["lmbd"]
assert sa.alpha() == sa_args["alpha"]
q.json_results_append(
    "cqlat-action-params: get_m_sq_scalar_action", sa.m_sq(), check_eps
)
q.json_results_append(
    "cqlat-action-params: get_lmbd_scalar_action", sa.lmbd(), check_eps
)
q.json_results_append(
    "cqlat-action-params: get_alpha_scalar_action", sa.alpha(), check_eps
)

# --- set_scalar_action: sa2 @= sa
sa2 = q.ScalarAction(1.0, 2.0, 0.5)
sa2 @= sa
assert sa2.m_sq() == sa.m_sq()
assert sa2.lmbd() == sa.lmbd()
assert sa2.alpha() == sa.alpha()
err_set = (
    abs(sa2.m_sq() - sa.m_sq())
    + abs(sa2.lmbd() - sa.lmbd())
    + abs(sa2.alpha() - sa.alpha())
)
q.json_results_append("cqlat-action-params: set_scalar_action", float(err_set), 0.0)
del sa2
import gc

gc.collect()

# --- to_mass_factor_scalar_action: v -> (1 + 2 asin(v) / pi)^2, in place
rs = q.RngState("cqlat-action-params")
sin_domega = q.FieldRealD(geo, 2)
sd_arr = 2.0 * rs.split("sd").u_rand_arr((geo.local_volume, 2)) - 1.0  # in [-1, 1)
np.asarray(sin_domega)[:] = sd_arr
sa.to_mass_factor(sin_domega)
sd_ref = ref_to_mass_factor(sd_arr)
err_mf = float(np.max(np.abs(np.asarray(sin_domega) - sd_ref)))
assert err_mf < check_eps, err_mf
q.json_results_append(
    "cqlat-action-params: to_mass_factor max error", err_mf, check_eps
)
q.json_results_append(
    "cqlat-action-params: to_mass_factor sum",
    q.glb_sum(float(np.asarray(sin_domega).sum())),
    1e-10,
)

# --- hmc_estimate_mass_scalar_action
rs_m = q.RngState("cqlat-action-params-mass")
f_re = rs_m.split("fr").u_rand_arr(latt_size + (2,))
f_im = rs_m.split("fi").u_rand_arr(latt_size + (2,))
g_re = rs_m.split("gr").u_rand_arr(latt_size + (2,))
g_im = rs_m.split("gi").u_rand_arr(latt_size + (2,))
# keep the origin component (which enters with the phi0 shift) away from zero
f_re[0, 0, 0, 0, 0] += 1.0
field_ft = mk_complex_field_from_global(geo, f_re, f_im)
force_ft = mk_complex_field_from_global(geo, g_re, g_im)
masses = q.FieldRealD(geo, 2)
np.asarray(masses)[:] = -1.0
phi0 = 0.3
sa.hmc_estimate_mass(masses, field_ft, force_ft, phi0)
mas_ref = ref_hmc_estimate_mass(
    np.asarray(field_ft), np.asarray(force_ft), orig_loc, phi0
)
err_hem = float(np.max(np.abs(np.asarray(masses) - mas_ref)))
assert err_hem < check_eps, err_hem
q.json_results_append(
    "cqlat-action-params: hmc_estimate_mass max error", err_hem, check_eps
)
mas0_g = np.where(orig_loc, np.asarray(masses)[:, 0], 0.0)
q.json_results_append(
    "cqlat-action-params: hmc_estimate_mass origin value",
    q.glb_sum(float(mas0_g.sum())),
    1e-10,
)
q.json_results_append(
    "cqlat-action-params: hmc_estimate_mass sum",
    q.glb_sum(float(np.asarray(masses).sum())),
    1e-10,
)
# the same through the raw cqlat export must agree with the Python wrapper
masses_direct = q.FieldRealD(geo, 2)
np.asarray(masses_direct)[:] = -1.0
qc.hmc_estimate_mass_scalar_action(sa, masses_direct, field_ft, force_ft, phi0)
err_hem_direct = float(np.max(np.abs(np.asarray(masses_direct) - np.asarray(masses))))
assert err_hem_direct == 0.0, err_hem_direct
q.json_results_append(
    "cqlat-action-params: hmc_estimate_mass wrapper vs direct cqlat",
    err_hem_direct,
    0.0,
)

# --- set_complex_from_double / set_double_from_complex
#     (ScalarAction methods that delegate to qlat.field_double; they used to
#     call the non-existent cqlat functions set_*_from_*_scalar_action and
#     raised AttributeError)
sf_real = q.FieldRealD(geo, 2)
sf_real_arr = rs.split("scfd").u_rand_arr((geo.local_volume, 2))
np.asarray(sf_real)[:] = sf_real_arr
cf_conv = q.FieldComplexD(geo, 2)
np.asarray(cf_conv)[:] = 0.0
sa.set_complex_from_double(cf_conv, sf_real)
err_cfd = float(np.max(np.abs(np.asarray(cf_conv) - sf_real_arr)))
assert err_cfd == 0.0, err_cfd
q.json_results_append(
    "cqlat-action-params: set_complex_from_double max err", err_cfd, check_eps
)
q.json_results_append(
    "cqlat-action-params: set_complex_from_double sum",
    q.glb_sum(float(np.asarray(cf_conv).real.sum())),
    1e-10,
)
# a genuinely complex source: only the real part may be written back
cf_conv2 = q.FieldComplexD(geo, 2)
cf_conv2_arr = sf_real_arr + 1j * rs.split("scfd-im").u_rand_arr((geo.local_volume, 2))
np.asarray(cf_conv2)[:] = cf_conv2_arr
sf_back = q.FieldRealD(geo, 2)
np.asarray(sf_back)[:] = -1.0
sa.set_double_from_complex(sf_back, cf_conv2)
err_dfc = float(np.max(np.abs(np.asarray(sf_back) - cf_conv2_arr.real)))
assert err_dfc == 0.0, err_dfc
q.json_results_append(
    "cqlat-action-params: set_double_from_complex max err", err_dfc, check_eps
)
q.json_results_append(
    "cqlat-action-params: set_double_from_complex sum",
    q.glb_sum(float(np.asarray(sf_back).sum())),
    1e-10,
)

# --- free_scalar_action, exercised by __del__
del sa
gc.collect()
q.json_results_append("cqlat-action-params: free_scalar_action (via __del__)", 1.0)

# =====================================================================
# FermionAction
# =====================================================================
fa = q.FermionAction(
    mass=fa_args["mass"],
    ls=fa_args["ls"],
    m5=fa_args["m5"],
    mobius_scale=fa_args["mobius_scale"],
)
assert fa.mass() == fa_args["mass"]
assert fa.ls() == fa_args["ls"]
assert fa.m5() == fa_args["m5"]
q.json_results_append("cqlat-action-params: get_ls_fermion_action", float(fa.ls()), 0.0)
q.json_results_append(
    "cqlat-action-params: get_mobius_scale_fermion_action (mobius)",
    fa.mobius_scale(),
    check_eps,
)
# a Mobius action has no omega list
assert fa.omega() is None
q.json_results_append(
    "cqlat-action-params: get_omega_fermion_action (mobius) is None", 1.0
)

fa2 = q.FermionAction(mass=0.5, ls=4, m5=1.0)
fa2 @= fa
assert fa2.mass() == fa.mass()
assert fa2.ls() == fa.ls()
assert fa2.m5() == fa.m5()
assert fa2.mobius_scale() == fa.mobius_scale()
q.json_results_append("cqlat-action-params: set_fermion_action", 1.0)
del fa2
gc.collect()

fa_z = q.FermionAction(
    mass=fa_args["mass"], ls=len(zm_omega), m5=fa_args["m5"], omega=zm_omega
)
assert fa_z.ls() == len(zm_omega)
omega_arr = np.asarray(fa_z.omega())
err_omega = float(np.max(np.abs(omega_arr - np.array(zm_omega))))
assert err_omega < check_eps, (omega_arr, zm_omega)
q.json_results_append(
    "cqlat-action-params: get_omega_fermion_action max error", err_omega, check_eps
)
q.json_results_append(
    "cqlat-action-params: get_omega_fermion_action real sum",
    float(omega_arr.real.sum()),
    1e-10,
)
q.json_results_append(
    "cqlat-action-params: get_omega_fermion_action imag sum",
    float(omega_arr.imag.sum()),
    1e-10,
)
assert fa_z.mobius_scale() == 0.0
q.json_results_append(
    "cqlat-action-params: get_mobius_scale_fermion_action (zmobius)",
    fa_z.mobius_scale(),
    0.0,
)

# --- free_fermion_action, exercised by __del__
del fa_z
gc.collect()
del fa
gc.collect()
q.json_results_append("cqlat-action-params: free_fermion_action (via __del__)", 1.0)

q.timer_display()
if q.is_test():
    q.check_log_json(__file__)
q.end_with_mpi()
q.displayln_info("CHECK: finished successfully.")
