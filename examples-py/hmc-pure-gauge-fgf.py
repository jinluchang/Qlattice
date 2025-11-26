#!/usr/bin/env python3

"""
Flexible gauge fixed HMC
"""

import sys
import math
import numpy as np

import qlat as q

from qlat_scripts.v1 import *

load_path_list[:] = [
        "results",
        ]

def mk_evolve_campostrini_from_leapfrog(evolve_leapfrog):
    """
    return evolve_campostrini
    """
    fac = 1.0 / (2.0 - math.cbrt(2.0))
    fac2 = 1.0 - 2 * fac
    @q.timer
    def evolve_campostrini(dt):
        evolve_leapfrog(fac * dt)
        evolve_leapfrog(fac2 * dt)
        evolve_leapfrog(fac * dt)
    return evolve_campostrini

def mk_evolve_force_gradient(qq_evolve, pp_evolve_fg):
    """
    return evolve_force_gradient
    """
    lam = 0.5 * (1.0 - 1.0 / math.sqrt(3.0))
    theta = (2.0 - math.sqrt(3.0)) / 48.0
    @q.timer
    def evolve_force_gradient(dt):
        ttheta = theta * dt * dt
        qq_evolve(lam * dt)
        pp_evolve_fg(4.0 * ttheta, 0.5 * dt)
        qq_evolve((1.0 - 2.0 * lam) * dt)
        pp_evolve_fg(4.0 * ttheta, 0.5 * dt)
        qq_evolve(lam * dt)
    return evolve_force_gradient

def mk_evolve_leapfrog(qq_evolve_1, qq_evolve_2, pp_evolve):
    """
    return evolve_leapfrog
    """
    @q.timer
    def evolve_leapfrog(dt):
        qq_evolve_1(0.5 * dt)
        pp_evolve(dt)
        qq_evolve_2(0.5 * dt)
    return evolve_leapfrog

@q.timer
def mk_mass_mats_for_gm(job_tag):
    r"""
    return sqrt_mass_matrix, mass_inv_matrix
    """
    block_site = q.Coordinate(get_param(job_tag, "hmc", "gauge_fixing", "block_site"))
    mat_dim = block_site.volume() * 4
    shape = (mat_dim, mat_dim,)
    sqrt_mass = get_param(job_tag, "hmc", "fourier_acceleration", "sqrt_mass")
    sqrt_mass_matrix = sqrt_mass * np.eye(mat_dim, dtype=np.float64)
    assert sqrt_mass_matrix.shape == shape
    sqrt_mass_inv_matrix = np.linalg.inv(sqrt_mass_matrix)
    mass_inv_matrix = sqrt_mass_inv_matrix @ sqrt_mass_inv_matrix
    return sqrt_mass_matrix, mass_inv_matrix

@q.timer
def mk_mass_mats_for_af(job_tag):
    r"""
    return sqrt_af_mass_matrix, af_mass_inv_matrix
    """
    block_site = q.Coordinate(get_param(job_tag, "hmc", "gauge_fixing", "block_site"))
    mat_dim = block_site.volume() * 4
    shape = (mat_dim, mat_dim,)
    sqrt_af_mass = get_param(job_tag, "hmc", "fourier_acceleration", "sqrt_af_mass")
    sqrt_af_mass_matrix = sqrt_af_mass * np.eye(mat_dim, dtype=np.float64)
    assert sqrt_af_mass_matrix.shape == shape
    sqrt_af_mass_inv_matrix = np.linalg.inv(sqrt_af_mass_matrix)
    af_mass_inv_matrix = sqrt_af_mass_inv_matrix @ sqrt_af_mass_inv_matrix
    return sqrt_af_mass_matrix, af_mass_inv_matrix

@q.timer
def mk_gm_v_from_gm(job_tag, geo, mass_inv_matrix):
    r"""
    return gm_v_from_gm
    or
    return af_v_from_af
    Usage:
    gm_v_from_gm = mk_gm_v_from_gm(job_tag, geo, mass_inv_matrix)
    af_v_from_af = mk_gm_v_from_gm(job_tag, geo, af_mass_inv_matrix)
    """
    total_site = geo.total_site
    block_site = q.Coordinate(get_param(job_tag, "hmc", "gauge_fixing", "block_site"))
    mat_dim = block_site.volume() * 4
    new_size_node = total_site // block_site
    assert isinstance(mass_inv_matrix, np.ndarray)
    assert mass_inv_matrix.shape == (mat_dim, mat_dim,)
    assert mass_inv_matrix.dtype == np.float64
    @q.timer
    def gm_v_from_gm(gm):
        assert gm.geo == geo
        f_basis = q.FieldRealD()
        q.set_basis_from_anti_hermitian_matrix(f_basis, gm)
        f_basis_list = q.shuffle_field(f_basis, new_size_node)
        for f_basis_local in f_basis_list:
            local_arr = np.asarray(f_basis_local).reshape((mat_dim, 8,), copy=False)
            assert local_arr.shape == (mat_dim, 8,)
            local_arr[:] = mass_inv_matrix @ local_arr[:]
        q.shuffle_field_back(f_basis, f_basis_list, new_size_node)
        gm_v = q.GaugeMomentum(geo)
        q.set_anti_hermitian_matrix_from_basis(gm_v, f_basis)
        return gm_v
    return gm_v_from_gm

@q.timer
def mk_gm_set_rand(job_tag, geo, sqrt_mass_matrix):
    r"""
    return gm_set_rand
    or
    return af_set_rand
    Usage:
    gm_set_rand = mk_gm_set_rand(job_tag, geo, sqrt_mass_matrix)
    af_set_rand = mk_gm_set_rand(job_tag, geo, sqrt_af_mass_matrix)
    """
    total_site = geo.total_site
    block_site = q.Coordinate(get_param(job_tag, "hmc", "gauge_fixing", "block_site"))
    mat_dim = block_site.volume() * 4
    new_size_node = total_site // block_site
    assert isinstance(sqrt_mass_matrix, np.ndarray)
    assert sqrt_mass_matrix.shape == (mat_dim, mat_dim,)
    assert sqrt_mass_matrix.dtype == np.float64
    @q.timer
    def gm_set_rand(gm, rs):
        assert gm.geo == geo
        assert gm.multiplicity == 4
        gm.set_rand(rs, 1.0)
        f_basis = q.FieldRealD()
        q.set_basis_from_anti_hermitian_matrix(f_basis, gm)
        f_basis_list = q.shuffle_field(f_basis, new_size_node)
        for f_basis_local in f_basis_list:
            local_arr = np.asarray(f_basis_local).reshape((mat_dim, 8,), copy=False)
            assert local_arr.shape == (mat_dim, 8,)
            local_arr[:] = sqrt_mass_matrix @ local_arr[:]
        q.shuffle_field_back(f_basis, f_basis_list, new_size_node)
        q.set_anti_hermitian_matrix_from_basis(gm, f_basis)
    return gm_set_rand

@q.timer
def gm_gm_v_hamilton_node(gm, gm_v):
    geo = gm.geo
    assert geo.is_only_local
    assert gm_v.geo == geo
    assert gm.multiplicity == 4
    assert gm_v.multiplicity == 4
    f_basis = q.FieldRealD()
    f_basis_v = q.FieldRealD()
    q.set_basis_from_anti_hermitian_matrix(f_basis, gm)
    q.set_basis_from_anti_hermitian_matrix(f_basis_v, gm_v)
    energy = np.sum(f_basis[:] * f_basis_v[:]).item()
    return energy

@q.timer
def gm_blocked_correlation(gm1, gm2, block_site):
    """
    return corr
    `corr` is global averaged correlation matrix.
    isinstance(corr, np.ndarray)
    corr.shape == (mat_dim, mat_dim,)
    corr.dtype == np.float64
    mat_dim = block_site.volume() * 4
    """
    geo = gm1.geo
    assert geo == gm2.geo
    assert gm1.multiplicity == 4
    assert gm2.multiplicity == 4
    total_site = geo.total_site
    mat_dim = block_site.volume() * 4
    new_size_node = total_site // block_site
    num_block = new_size_node.volume()
    f1_basis = q.FieldRealD()
    f2_basis = q.FieldRealD()
    q.set_basis_from_anti_hermitian_matrix(f1_basis, gm1)
    q.set_basis_from_anti_hermitian_matrix(f2_basis, gm2)
    f1_basis_list = q.shuffle_field(f1_basis, new_size_node)
    f2_basis_list = q.shuffle_field(f2_basis, new_size_node)
    assert len(f1_basis_list) == len(f2_basis_list)
    corr = np.zeros((mat_dim, mat_dim,), dtype=np.float64)
    for f1_basis_local, f2_basis_local in zip(f1_basis_list, f2_basis_list):
        local_arr1 = np.asarray(f1_basis_local).reshape((mat_dim, 8,), copy=False)
        local_arr2 = np.asarray(f2_basis_local).reshape((mat_dim, 8,), copy=False)
        assert local_arr1.shape == (mat_dim, 8,)
        assert local_arr2.shape == (mat_dim, 8,)
        corr += np.sum(local_arr1[:, None, :] * local_arr2[None, :, :], axis=-1)
    corr = q.glb_sum(corr)
    corr /= num_block
    return corr

runtime_info = dict()

@q.timer
def mk_acc_runtime_info(job_tag, ga, get_gm_force, af_v_from_af):
    """
    return acc_runtime_info
    acc_runtime_info = mk_acc_runtime_info(job_tag, ga, get_gm_force, af_v_from_af)
    """
    block_site = q.Coordinate(get_param(job_tag, "hmc", "gauge_fixing", "block_site"))
    runtime_info.clear()
    runtime_info["list"] = []
    step = 0
    gm_init = None
    gm_v_init = None
    gf_init = None
    af_init = None
    gm_force_init = None
    gm_force_qcd_init = None
    gm_force_gauge_fixing_init = None
    energy_init = None
    @q.timer
    def acc_runtime_info(time, gm, gm_v, gf, af):
        nonlocal step, gm_init, gm_v_init, gf_init, af_init
        gm_force = get_gm_force(gf, af)
        gm_force_qcd = get_gm_force(gf, af, tag="qcd")
        gm_force_gauge_fixing = get_gm_force(gf, af, tag="gauge_fixing")
        energy = (
            gm_gm_v_hamilton_node(gm, gm_v)
            + gm_gm_v_hamilton_node(af, af_v_from_af(af))
            + q.gf_hamilton_node(gf, ga)
            )
        nonlocal gm_force_init, gm_force_qcd_init, gm_force_gauge_fixing_init, energy_init
        if step == 0:
            assert time == 0.0
            gm_init = gm.copy()
            gm_v_init = gm_v.copy()
            gf_init = gf.copy()
            af_init = af.copy()
            gm_force_init = gm_force.copy()
            gm_force_qcd_init = gm_force_qcd.copy()
            gm_force_gauge_fixing_init = gm_force_gauge_fixing.copy()
            energy_init = energy
        info = dict()
        info["time"] = time
        info["plaq"] = gf.plaq()
        info["link_trace"] = gf.link_trace()
        info["delta_h"] = q.glb_sum(energy - energy_init)
        info["gm_corr"] = gm_blocked_correlation(gm_init, gm, block_site)
        info["gm_force_corr"] = gm_blocked_correlation(gm_force_init, gm_force, block_site)
        # info["gm_force_qcd_corr"] = gm_blocked_correlation(gm_force_qcd_init, gm_force_qcd, block_site)
        # info["gm_force_gauge_fixing_corr"] = gm_blocked_correlation(gm_force_gauge_fixing_init, gm_force_gauge_fixing, block_site)
        info["gm_gm_force_corr"] = gm_blocked_correlation(gm_init, gm_force, block_site)
        info["gm_force_gm_corr"] = gm_blocked_correlation(gm_force_init, gm, block_site)
        runtime_info["list"].append(info)
        step += 1
    return acc_runtime_info

@q.timer
def mk_fgf(job_tag, rs):
    r"""
    return fgf, fgf_g
    gt_inv = fgf(gf)
    And:
    fgf_g(gf, gm_v) * diff_eps
    \approx
    fgf(gf2) - fgf(gf1)
    where
    gf2 = gf.copy()
    q.gf_evolve(gf2, gm_v, 0.5 * diff_eps)
    gf1 = gf.copy()
    q.gf_evolve(gf1, gm_v, -0.5 * diff_eps)
    #
    Assume `rs` is already the rng for this traj.
    """
    rs_f_dir = rs.split("seed-gt_block_tree_gauge-rs_f_dir")
    f_dir_list = None
    block_site = q.Coordinate(get_param(job_tag, "hmc", "gauge_fixing", "block_site"))
    new_size_node = q.Coordinate(get_param(job_tag, "hmc", "gauge_fixing", "new_size_node"))
    stout_smear_step_size = get_param(job_tag, "hmc", "gauge_fixing", "stout_smear_step_size")
    num_smear_step = get_param(job_tag, "hmc", "gauge_fixing", "num_smear_step")
    diff_eps = get_param(job_tag, "hmc", "diff_eps")
    @q.timer
    def fgf(gf):
        nonlocal f_dir_list
        gt_inv, f_dir_list = q.gt_block_tree_gauge(
            gf,
            block_site=block_site,
            new_size_node=new_size_node,
            stout_smear_step_size=stout_smear_step_size,
            num_smear_step=num_smear_step,
            f_dir_list=f_dir_list,
            rs_f_dir=rs_f_dir,
            )
        return gt_inv
    @q.timer
    def fgf_g(gf, gm_v):
        gf2 = gf.copy()
        gf1 = gf.copy()
        q.gf_evolve(gf2, gm_v, 0.5 * diff_eps)
        q.gf_evolve(gf1, gm_v, -0.5 * diff_eps)
        gt = fgf(gf2)
        gt -= fgf(gf1)
        gt *= 1.0 / diff_eps
        q.set_tr_less_anti_herm_matrix(gt)
        return gt
    return fgf, fgf_g

@q.timer
def mk_fgf_gf_evolve(job_tag, geo, fgf, fgf_g):
    """
    Evolve `gf` in place.
    `gf` should be initially in the gauge fixed state.
    gt = fgf(gf)
    gf_gfixed = gt.inv() * gf
    """
    implicity_integrator_eps = get_param(job_tag, "hmc", "implicity_integrator_eps")
    implicity_integrator_max_iter = get_param(job_tag, "hmc", "implicity_integrator_max_iter")
    gt_unit = q.GaugeTransform(geo)
    gt_unit.set_unit()
    gt_norm = q.qnorm(gt_unit)
    gf_norm = gt_norm * 4
    @q.timer
    def gf_evolve(gf, af, gm_v, dt, *, tag=None, is_initial_gauge_fixed=False):
        # Update `gf` and `af` in place.
        #
        # Always return gauge field in a gauge fixed and unitarized state.
        #
        # tag in [ None, "fix_midpoint", "fix_endpoint", "fix_start", "fix_stop", "no_fix", ]
        #
        # Note that `af` represents $(1/\sqrt{2}) A_a(x,\mu) T_a$
        #
        # Preparation
        fname = q.get_fname()
        #
        if tag is None:
            tag = "fix_midpoint"
        #
        gf_init = gf
        af_init = af
        #
        if not is_initial_gauge_fixed:
            gf = fgf(gf).inv() * gf
            gf.unitarize()
        #
        egm_p = None
        egm_p_dt = None
        #
        @q.timer
        def evolve_fix_stop(dt_step):
            nonlocal gf, af, egm_p, egm_p_dt
            #
            if egm_p_dt != dt_step:
                egm_p = q.field_color_matrix_exp(gm_v, dt_step)
                egm_p_dt = dt_step
            #
            gf0 = gf
            gf = q.GaugeField()
            # Find proper initial gauge transformation so that the midpoint evolution is exactly gauge fixed
            gt = q.GaugeTransform(geo)
            gt.set_unit()
            for i in range(implicity_integrator_max_iter):
                gf.swap(q.field_color_matrix_mul(egm_p, gt * gf0))
                gt_new = fgf(gf)
                gt = gt_new.inv() * gt
                gt.unitarize()
                gt_new -= gt_unit
                gf_eps = np.sqrt(q.qnorm(gt_new) / gt_norm).item()
                # q.displayln_info(1, f"{fname}: iter={i} gf_eps: {gf_eps} (target: {implicity_integrator_eps})")
                if gf_eps < implicity_integrator_eps:
                    break
            q.displayln_info(0, f"{fname}: iter={i} gf_eps: {gf_eps} (target: {implicity_integrator_eps})")
            #
            # Evolve a step
            gf.swap(q.field_color_matrix_mul(egm_p, gt * gf0))
            gf = fgf(gf).inv() * gf
            gf.unitarize()
            #
            gf_fixed = gf
            #
            dg = fgf_g(gf_fixed, af)
            prod = q.field_color_matrix_mul(gm_v, dg)
            prod -= q.field_color_matrix_mul(dg, gm_v)
            q.set_tr_less_anti_herm_matrix(prod)
            prod *= dt_step
            af += prod
            #
            af_fixed = af
        #
        @q.timer
        def evolve_fix_start(dt_step):
            nonlocal gf, af, egm_p, egm_p_dt
            #
            if egm_p_dt != dt_step:
                egm_p = q.field_color_matrix_exp(gm_v, dt_step)
                egm_p_dt = dt_step
            #
            gf_fixed = gf
            gf = q.GaugeField()
            #
            # Evolve a step
            gf.swap(q.field_color_matrix_mul(egm_p, gf_fixed))
            gf = fgf(gf).inv() * gf
            gf.unitarize()
            #
            af_fixed = af.copy()
            #
            for i in range(implicity_integrator_max_iter):
                dg = fgf_g(gf_fixed, af)
                prod = q.field_color_matrix_mul(gm_v, dg)
                prod -= q.field_color_matrix_mul(dg, gm_v)
                q.set_tr_less_anti_herm_matrix(prod)
                prod *= dt_step
                af_diff = af
                af = af_fixed.copy()
                af += prod
                af_diff -= af
                af_eps = np.sqrt(q.qnorm(af_diff) / gf_norm).item()
                if af_eps < implicity_integrator_eps:
                    break
            q.displayln_info(0, f"{fname}: iter={i} af_eps: {af_eps} (target: {implicity_integrator_eps})")
        #
        @q.timer
        def evolve_no_fix(dt_step):
            nonlocal gf, af, egm_p, egm_p_dt
            #
            if egm_p_dt != dt_step:
                egm_p = q.field_color_matrix_exp(gm_v, dt_step)
                egm_p_dt = dt_step
            #
            gf0 = gf
            gf = q.GaugeField()
            #
            gf.swap(q.field_color_matrix_mul(egm_p, gf0))
            gf = fgf(gf).inv() * gf
            gf.unitarize()
            #
            dg = fgf_g(gf0, af)
            prod = q.field_color_matrix_mul(gm_v, dg)
            prod -= q.field_color_matrix_mul(dg, gm_v)
            q.set_tr_less_anti_herm_matrix(prod)
            prod *= dt_step
            af += prod
        #
        @q.timer
        def evolve_fix_endpoint(dt_step):
            evolve_fix_start(0.5 * dt_step)
            evolve_fix_stop(0.5 * dt_step)
        #
        @q.timer
        def evolve_fix_midpoint(dt_step):
            evolve_fix_stop(0.5 * dt_step)
            evolve_fix_start(0.5 * dt_step)
        #
        evolve_fix_endpoint_4th = mk_evolve_campostrini_from_leapfrog(evolve_fix_endpoint)
        evolve_fix_midpoint_4th = mk_evolve_campostrini_from_leapfrog(evolve_fix_midpoint)
        #
        if tag == "no_fix":
            evolve_no_fix(dt)
        elif tag == "fix_start":
            evolve_fix_start(dt)
        elif tag == "fix_stop":
            evolve_fix_stop(dt)
        elif tag == "fix_midpoint":
            evolve_fix_midpoint(dt)
        elif tag == "fix_endpoint":
            evolve_fix_endpoint(dt)
        elif tag == "fix_midpoint_4th":
            evolve_fix_midpoint_4th(dt)
        elif tag == "fix_endpoint_4th":
            evolve_fix_endpoint_4th(dt)
        else:
            raise Exception(f"{fname}: tag={tag}")
        #
        # Update input object
        #
        gf_init.swap(gf)
        af_init.swap(af)
    return gf_evolve

@q.timer
def mk_fgf_get_gm_force(job_tag, geo, fgf_g, af_v_from_af):
    """
    return `get_gm_force` for `gf` and `af`.
    `gf` should be in the gauge fixed state.
    """
    beta = get_param(job_tag, "hmc", "beta")
    c1 = get_param(job_tag, "hmc", "c1")
    ga = q.GaugeAction(beta, c1)
    total_volume = geo.total_volume
    @q.timer
    def get_gm_force(gf, af, tag=None):
        #
        fname = q.get_fname()
        #
        if tag is None:
            tag = "all"
        #
        assert tag in ["all", "qcd", "gauge_fixing", ]
        #
        if tag in [ "qcd", "all", ]:
            # Set QCD gauge action force
            gm_force_qcd = q.GaugeMomentum()
            q.set_gm_force(gm_force_qcd, gf, ga)
            force_size_qcd = math.sqrt(q.qnorm(gm_force_qcd) / (total_volume * 4))
            q.displayln_info(0, f"force_size_qcd={force_size_qcd:.5f}")
        #
        if tag in [ "gauge_fixing", "all", ]:
            # Add force due to gauge fixing
            geo = gf.geo
            gm_force_gauge_fixing = q.GaugeMomentum(geo)
            gm_force_gauge_fixing.set_zero()
            dg = fgf_g(gf, af)
            af_v = af_v_from_af(af)
            prod = q.field_color_matrix_mul(af_v, dg)
            prod -= q.field_color_matrix_mul(dg, af_v)
            q.set_tr_less_anti_herm_matrix(prod)
            gm_force_gauge_fixing += prod
            force_size_gauge_fixing = math.sqrt(q.qnorm(gm_force_gauge_fixing) / (total_volume * 4))
            q.displayln_info(0, f"force_size_gauge_fixing={force_size_gauge_fixing:.5f}")
        #
        if tag == "all":
            geo = gf.geo
            gm_force = q.GaugeMomentum(geo)
            gm_force.set_zero()
            gm_force += gm_force_qcd
            gm_force += gm_force_gauge_fixing
            force_size_total = math.sqrt(q.qnorm(gm_force) / (total_volume * 4))
            q.displayln_info(0, f"force_size_total={force_size_total:.5f}")
            cos_alpha_qcd_gauge_fixing = (
                (force_size_total**2 - force_size_qcd**2 - force_size_gauge_fixing**2)
                / (force_size_qcd * force_size_gauge_fixing))
            q.displayln_info(0, f"cos_alpha_qcd_gauge_fixing={cos_alpha_qcd_gauge_fixing:.5f}")
        elif tag == "qcd":
            gm_force = gm_force_qcd
        elif tag == "gauge_fixing":
            gm_force = gm_force_gauge_fixing
        else:
            raise Exception(f"{fname}: tag={tag}")
        #
        return gm_force
    return get_gm_force

@q.timer
def mk_fgf_gm_evolve_fg(get_gm_force, gf_evolve, gm_v_from_gm):
    """
    return `gm_evolve_fg` for `gf` and `af`.
    `gf` should be in the gauge fixed state.
    """
    @q.timer
    def gm_evolve_fg(gm, gm_v, gf, af, fg_dt, dt):
        """
        Modify `gm` and `gm_v` in place.
        """
        gm_force = get_gm_force(gf, af)
        if fg_dt != 0.0:
            gm_force_v = gm_v_from_gm(gm_force)
            gf_g = gf.copy()
            af_g = af.copy()
            gf_evolve(gf_g, af_g, gm_force_v, fg_dt, tag="no_fix", is_initial_gauge_fixed=True)
            gm_force = get_gm_force(gf_g, af_g)
        gm_force *= dt
        gm += gm_force
        gm_v.swap(gm_v_from_gm(gm))
    return gm_evolve_fg

@q.timer(is_timer_fork=True)
def run_hmc_evolve_pure_gauge(
        gm, gm_v, gf, af,
        *,
        ga,
        fgf,
        gf_evolve,
        gm_evolve_fg,
        af_v_from_af,
        acc_runtime_info,
        gf_integrator_tag,
        n_step,
        md_time,
        ):
    fname = q.get_fname()
    @q.timer
    def qq_evolve_4th(dt):
        gf_evolve(gf, af, gm_v, dt, tag="fix_endpoint_4th", is_initial_gauge_fixed=True)
    @q.timer
    def pp_evolve_fg(fg_dt, dt):
        gm_evolve_fg(gm, gm_v, gf, af, fg_dt, dt)
    @q.timer
    def qq_evolve_1(dt):
        gf_evolve(gf, af, gm_v, dt, tag="fix_start", is_initial_gauge_fixed=True)
    @q.timer
    def qq_evolve_2(dt):
        gf_evolve(gf, af, gm_v, dt, tag="fix_stop", is_initial_gauge_fixed=True)
    @q.timer
    def pp_evolve(dt):
        pp_evolve_fg(0.0, dt)
    evolve_leapfrog = mk_evolve_leapfrog(qq_evolve_1, qq_evolve_2, pp_evolve)
    evolve_campostrini = mk_evolve_campostrini_from_leapfrog(evolve_leapfrog)
    evolve_force_gradient = mk_evolve_force_gradient(qq_evolve_4th, pp_evolve_fg)
    if gf_integrator_tag == "leapfrog":
        evolve = evolve_leapfrog
    elif gf_integrator_tag == "campostrini":
        evolve = evolve_campostrini
    elif gf_integrator_tag == "force_gradient":
        evolve = evolve_force_gradient
    else:
        raise Exception(f"{fname}: gf_integrator_tag={gf_integrator_tag}")
    energy = (
        gm_gm_v_hamilton_node(gm, gm_v)
        + gm_gm_v_hamilton_node(af, af_v_from_af(af))
        + q.gf_hamilton_node(gf, ga)
        )
    gf @= fgf(gf).inv() * gf
    gf.unitarize()
    dt = md_time / n_step
    time = 0.0
    if acc_runtime_info is not None:
        acc_runtime_info(time, gm, gm_v, gf, af)
    for i in range(n_step):
        evolve(dt)
        time += dt
        if acc_runtime_info is not None:
            acc_runtime_info(time, gm, gm_v, gf, af)
    gf.unitarize()
    delta_h = (
        gm_gm_v_hamilton_node(gm, gm_v)
        + gm_gm_v_hamilton_node(af, af_v_from_af(af))
        + q.gf_hamilton_node(gf, ga)
        - energy
        )
    delta_h = q.glb_sum(delta_h)
    return delta_h

@q.timer(is_timer_fork=True)
def run_hmc_pure_gauge(job_tag, gf, traj, rs):
    fname = q.get_fname()
    gf_in = gf
    gf = gf.copy()
    rs = rs.split(f"{traj}")
    geo = gf.geo
    sqrt_mass_matrix, mass_inv_matrix = mk_mass_mats_for_gm(job_tag)
    sqrt_af_mass_matrix, af_mass_inv_matrix = mk_mass_mats_for_af(job_tag)
    gm_set_rand = mk_gm_set_rand(job_tag, geo, sqrt_mass_matrix)
    af_set_rand = mk_gm_set_rand(job_tag, geo, sqrt_af_mass_matrix)
    gm_v_from_gm = mk_gm_v_from_gm(job_tag, geo, mass_inv_matrix)
    af_v_from_af = mk_gm_v_from_gm(job_tag, geo, af_mass_inv_matrix)
    fgf, fgf_g, = mk_fgf(job_tag, rs)
    gf_evolve = mk_fgf_gf_evolve(job_tag, geo, fgf, fgf_g)
    get_gm_force = mk_fgf_get_gm_force(job_tag, geo, fgf_g, af_v_from_af)
    gm_evolve_fg = mk_fgf_gm_evolve_fg(get_gm_force, gf_evolve, gm_v_from_gm)
    gf_integrator_tag = get_param(job_tag, "hmc", "gf_integrator_tag")
    md_time = get_param(job_tag, "hmc", "md_time")
    n_step = get_param(job_tag, "hmc", "n_step")
    beta = get_param(job_tag, "hmc", "beta")
    c1 = get_param(job_tag, "hmc", "c1")
    ga = q.GaugeAction(beta, c1)
    acc_runtime_info = mk_acc_runtime_info(job_tag, ga, get_gm_force, af_v_from_af)
    max_traj_always_accept = get_param(job_tag, "hmc", "max_traj_always_accept")
    max_traj_reverse_test= get_param(job_tag, "hmc", "max_traj_reverse_test")
    is_always_accept = traj < max_traj_always_accept
    is_reverse_test = traj < max_traj_reverse_test
    total_site = geo.total_site
    xg_field_shift = rs.split("xg_field_shift").c_rand_gen(total_site)
    gm = q.GaugeMomentum(geo)
    gm_set_rand(gm, rs.split("set_rand_gauge_momentum"))
    gm_v = q.GaugeMomentum()
    gm_v.swap(gm_v_from_gm(gm))
    af = q.GaugeMomentum(geo)
    af_set_rand(af, rs.split("set_rand_af"))
    gf = gf.shift(xg_field_shift)
    gf0 = gf.copy()
    delta_h = run_hmc_evolve_pure_gauge(
        gm, gm_v, gf, af, ga=ga,
        fgf=fgf, gf_evolve=gf_evolve, gm_evolve_fg=gm_evolve_fg,
        acc_runtime_info=acc_runtime_info,
        af_v_from_af=af_v_from_af,
        gf_integrator_tag=gf_integrator_tag,
        n_step=n_step, md_time=md_time,
        )
    gf = gf.shift(-xg_field_shift)
    if is_reverse_test:
        gf_r = gf.copy()
        gm_r = gm.copy()
        gm_v_r = gm_v.copy()
        af_r = af.copy()
        gf_r = gf_r.shift(xg_field_shift)
        delta_h_rev = run_hmc_evolve_pure_gauge(
            gm_r, gm_v_r, gf_r, af_r, ga=ga,
            fgf=fgf, gf_evolve=gf_evolve, gm_evolve_fg=gm_evolve_fg,
            acc_runtime_info=None,
            af_v_from_af=af_v_from_af,
            gf_integrator_tag=gf_integrator_tag,
            n_step=n_step, md_time=-md_time,
            )
        gf_r -= fgf(gf0).inv() * gf0
        gf_r = gf_r.shift(-xg_field_shift)
        q.displayln_info(f"{fname}: reversed delta_diff: {delta_h + delta_h_rev} / {delta_h}")
        gf_diff_norm = q.qnorm(gf_r)
        gf_norm = q.qnorm(gf_in)
        q.displayln_info(f"{fname}: reversed gf_diff: {gf_diff_norm} / {gf_norm}")
        assert gf_diff_norm <= 1e-16 * gf_norm
    flag, accept_prob = q.metropolis_accept(delta_h, traj, rs.split("metropolis_accept"))
    if flag or is_always_accept:
        q.displayln_info(f"{fname}: update gf (traj={traj})")
        gf_in @= gf
    return flag, delta_h

@q.timer(is_timer_fork=True)
def run_topo_info(job_tag, traj, gf):
    info_path = get_save_path(f"{job_tag}/topo-measure-wilson-flow/traj-{traj}")
    flow_time = 6
    flow_n_step = 80
    smear_info_list = [
            [ 1.0 / flow_n_step, flow_n_step, 0.0, "runge-kutta", ],
            ] * flow_time
    energy_derivative_info = [ 1.0 / flow_n_step, 0.0, "runge-kutta", ]
    topo_list, energy_list, = q.smear_measure_topo(
            gf.copy(),
            smear_info_list=smear_info_list,
            energy_derivative_info=energy_derivative_info,
            info_path=info_path,
            density_field_path=info_path,
            )

@q.timer(is_timer_fork=True)
def run_hmc(job_tag):
    fname = q.get_fname()
    rs = q.RngState(f"run_hmc-{get_job_seed(job_tag)}")
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    max_traj = get_param(job_tag, "hmc", "max_traj")
    save_traj_interval = get_param(job_tag, "hmc", "save_traj_interval")
    is_saving_topo_info = get_param(job_tag, "hmc", "is_saving_topo_info")
    geo = q.Geometry(total_site)
    gf = q.GaugeField(geo)
    traj_load = None
    if get_load_path(f"{job_tag}/configs") is not None:
        for traj in range(max_traj):
            fn = get_load_path(f"{job_tag}/configs/ckpoint_lat.{traj}")
            if fn is not None:
                traj_load = traj
    if traj_load is None:
        traj = 0
        gf.set_rand(rs.split("init"), 0.1, 2)
        gf.save(get_save_path(f"{job_tag}/configs/ckpoint_lat.{traj}"))
    else:
        traj = traj_load
        gf.load(get_load_path(f"{job_tag}/configs/ckpoint_lat.{traj}"))
    for traj in range(traj, max_traj):
        traj += 1
        flag, delta_h = run_hmc_pure_gauge(
                job_tag,
                gf, traj, rs.split("run_hmc_pure_gauge"),
                )
        plaq = gf.plaq()
        info = dict()
        info["traj"] = traj
        info["plaq"] = plaq
        info["flag"] = flag
        info["delta_h"] = delta_h
        q.qtouch_info(get_save_path(f"{job_tag}/configs/ckpoint_lat_info.{traj}.txt"), pformat(info))
        q.json_results_append(f"{fname}: {traj} plaq", plaq, 1e-10)
        q.json_results_append(f"{fname}: {traj} delta_h", delta_h, 1e-4)
        q.save_pickle_obj(runtime_info, get_save_path(f"{job_tag}/runtime_info/traj-{traj}.pickle"))
        if traj % save_traj_interval == 0:
            gf.save(get_save_path(f"{job_tag}/configs/ckpoint_lat.{traj}"))
            if is_saving_topo_info:
                run_topo_info(job_tag, traj, gf)

# ----

job_tag = "test-4nt8"
set_param(job_tag, "total_site")((4, 4, 4, 8,))
set_param(job_tag, "hmc", "max_traj")(4)
set_param(job_tag, "hmc", "max_traj_always_accept")(3)
set_param(job_tag, "hmc", "max_traj_reverse_test")(2)
set_param(job_tag, "hmc", "md_time")(0.6)
set_param(job_tag, "hmc", "n_step")(2)
set_param(job_tag, "hmc", "beta")(2.13)
set_param(job_tag, "hmc", "c1")(-0.331)
set_param(job_tag, "hmc", "diff_eps")(1e-5)
set_param(job_tag, "hmc", "gf_integrator_tag")("force_gradient")
set_param(job_tag, "hmc", "implicity_integrator_eps")(1e-11)
set_param(job_tag, "hmc", "implicity_integrator_max_iter")(50)
set_param(job_tag, "hmc", "gauge_fixing", "block_site")((4, 4, 4, 4,))
set_param(job_tag, "hmc", "gauge_fixing", "new_size_node")((1, 1, 1, 2,))
set_param(job_tag, "hmc", "gauge_fixing", "stout_smear_step_size")(0.125)
set_param(job_tag, "hmc", "gauge_fixing", "num_smear_step")(4)
set_param(job_tag, "hmc", "fourier_acceleration", "sqrt_mass")(2.0)
set_param(job_tag, "hmc", "fourier_acceleration", "sqrt_af_mass")(2.0)
set_param(job_tag, "hmc", "save_traj_interval")(4)
set_param(job_tag, "hmc", "is_saving_topo_info")(True)

job_tag = "16I_b2p8_fgf_md4"
set_param(job_tag, "total_site")((16, 16, 16, 32,))
set_param(job_tag, "a_inv_gev")(2.646) # 2003 lattice spacing 0309017.pdf
set_param(job_tag, "hmc", "max_traj")(5000)
set_param(job_tag, "hmc", "max_traj_always_accept")(0)
set_param(job_tag, "hmc", "max_traj_reverse_test")(0)
set_param(job_tag, "hmc", "md_time")(4.0)
set_param(job_tag, "hmc", "n_step")(32 * 4)
set_param(job_tag, "hmc", "beta")(2.80)
set_param(job_tag, "hmc", "c1")(-0.331)
set_param(job_tag, "hmc", "diff_eps")(1e-5)
set_param(job_tag, "hmc", "gf_integrator_tag")("force_gradient")
set_param(job_tag, "hmc", "implicity_integrator_eps")(1e-11)
set_param(job_tag, "hmc", "implicity_integrator_max_iter")(50)
set_param(job_tag, "hmc", "gauge_fixing", "block_site")((4, 4, 4, 4,))
set_param(job_tag, "hmc", "gauge_fixing", "new_size_node")((1, 1, 1, 2,))
set_param(job_tag, "hmc", "gauge_fixing", "stout_smear_step_size")(0.125)
set_param(job_tag, "hmc", "gauge_fixing", "num_smear_step")(6)
set_param(job_tag, "hmc", "fourier_acceleration", "sqrt_mass")(1.0)
set_param(job_tag, "hmc", "fourier_acceleration", "sqrt_af_mass")(2.0)
set_param(job_tag, "hmc", "save_traj_interval")(1)
set_param(job_tag, "hmc", "is_saving_topo_info")(True)

job_tag = "32I_b2p8_fgf_md4"
set_param(job_tag, "total_site")((32, 32, 32, 64,))
set_param(job_tag, "a_inv_gev")(2.646) # 2003 lattice spacing 0309017.pdf
set_param(job_tag, "hmc", "max_traj")(5000)
set_param(job_tag, "hmc", "max_traj_always_accept")(0)
set_param(job_tag, "hmc", "max_traj_reverse_test")(0)
set_param(job_tag, "hmc", "md_time")(4.0)
set_param(job_tag, "hmc", "n_step")(32 * 4)
set_param(job_tag, "hmc", "beta")(2.80)
set_param(job_tag, "hmc", "c1")(-0.331)
set_param(job_tag, "hmc", "diff_eps")(1e-5)
set_param(job_tag, "hmc", "gf_integrator_tag")("force_gradient")
set_param(job_tag, "hmc", "implicity_integrator_eps")(1e-11)
set_param(job_tag, "hmc", "implicity_integrator_max_iter")(50)
set_param(job_tag, "hmc", "gauge_fixing", "block_site")((4, 4, 4, 4,))
set_param(job_tag, "hmc", "gauge_fixing", "new_size_node")((1, 1, 1, 2,))
set_param(job_tag, "hmc", "gauge_fixing", "stout_smear_step_size")(0.125)
set_param(job_tag, "hmc", "gauge_fixing", "num_smear_step")(6)
set_param(job_tag, "hmc", "fourier_acceleration", "sqrt_mass")(1.0)
set_param(job_tag, "hmc", "fourier_acceleration", "sqrt_af_mass")(2.0)
set_param(job_tag, "hmc", "save_traj_interval")(1)
set_param(job_tag, "hmc", "is_saving_topo_info")(True)

# ----

size_node_list = [
        [1, 1, 1, 1],
        [1, 1, 1, 2],
        [1, 1, 1, 3],
        [1, 1, 1, 4],
        [1, 1, 1, 6],
        [1, 1, 1, 8],
        [1, 2, 2, 4],
        [2, 2, 2, 4],
        [2, 2, 2, 4],
        ]

if __name__ == "__main__":

    q.begin_with_mpi(size_node_list)

    ##################### CMD options #####################

    job_tag_list = q.get_arg("--job_tag_list", default="").split(",")

    #######################################################

    job_tag_list_default = [
            "test-4nt8",
            ]

    if job_tag_list == [ "", ]:
        job_tag_list = job_tag_list_default

    for job_tag in job_tag_list:
        run_params(job_tag)
        run_hmc(job_tag)

    q.timer_display()

    q.check_log_json(__file__)
    q.end_with_mpi()
    q.displayln_info(f"CHECK: finished successfully.")

# ----
