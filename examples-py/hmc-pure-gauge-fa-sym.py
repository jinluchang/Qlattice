#!/usr/bin/env python3

json_results = []

import sys
import math as m
import numpy as np

import qlat as q

from qlat_scripts.v1 import *

from qlat import (
        FieldRealD,
        GaugeField,
        GaugeMomentum,
        metropolis_accept,
        set_gm_force,
        set_gm_force_dual,
        gm_hamilton_node_fa,
        gf_hamilton_node,
        project_gauge_transform,
        )

load_path_list[:] = [
        "results",
        ]

# ----

@q.timer
def gf_evolve(gf, gm, gm_dual, mf, mf_dual, dt):
    q.gf_evolve_fa(gf, gm, mf, dt)
    q.gf_evolve_fa_dual(gf, gm_dual, mf_dual, dt)

@q.timer
def gm_evolve_fg(gm, gm_dual, gf_init, mf, mf_dual, ga, fg_dt, dt):
    fname = q.get_fname()
    geo = gf_init.geo
    gf = GaugeField(geo)
    gf @= gf_init
    gm_force = GaugeMomentum(geo)
    gm_force_dual = GaugeMomentum(geo)
    set_gm_force(gm_force, gf, ga)
    set_gm_force_dual(gm_force_dual, gf, gm_force)
    project_gauge_transform(gm_force, gm_force_dual, mf, mf_dual)
    gf_evolve(gf, gm_force, gm_force_dual, mf, mf_dual, fg_dt)
    set_gm_force(gm_force, gf, ga)
    set_gm_force_dual(gm_force_dual, gf, gm_force)
    gm_force *= dt
    gm_force_dual *= dt
    gm += gm_force
    gm_dual += gm_force_dual
    project_gauge_transform(gm, gm_dual, mf, mf_dual)

@q.timer_verbose
def run_hmc_evolve(gm, gm_dual, gf, mf, mf_dual, ga, rs, n_step, md_time=1.0):
    energy = gm_hamilton_node_fa(gm, mf) + gm_hamilton_node_fa(gm_dual, mf_dual) + gf_hamilton_node(gf, ga)
    dt = md_time / n_step
    lam = 0.5 * (1.0 - 1.0 / math.sqrt(3.0));
    theta = (2.0 - math.sqrt(3.0)) / 48.0;
    ttheta = theta * dt * dt * dt;
    gf_evolve(gf, gm, gm_dual, mf, mf_dual, lam * dt)
    for i in range(n_step):
        gm_evolve_fg(gm, gm_dual, gf, mf, mf_dual, ga, 4.0 * ttheta / dt, 0.5 * dt);
        gf_evolve(gf, gm, gm_dual, mf, mf_dual, (1.0 - 2.0 * lam) * dt);
        gm_evolve_fg(gm, gm_dual, gf, mf, mf_dual, ga, 4.0 * ttheta / dt, 0.5 * dt);
        if i < n_step - 1:
            gf_evolve(gf, gm, gm_dual, mf, mf_dual, 2.0 * lam * dt);
        else:
            gf_evolve(gf, gm, gm_dual, mf, mf_dual, lam * dt);
    gf.unitarize()
    delta_h = gm_hamilton_node_fa(gm, mf) + gm_hamilton_node_fa(gm_dual, mf_dual) + gf_hamilton_node(gf, ga) - energy;
    delta_h = q.glb_sum(delta_h)
    return delta_h

@q.timer_verbose
def run_hmc_traj(
        gf,
        gm, gm_dual,
        mf, mf_dual,
        ga, traj,
        rs,
        *,
        is_reverse_test=False,
        n_step=6,
        md_time=1.0,
        is_always_accept=False,
        ):
    fname = q.get_fname()
    rs = rs.split(f"{traj}")
    geo = gf.geo
    gf0 = GaugeField(geo)
    gf0 @= gf
    gm0 = GaugeMomentum(geo)
    gm_dual0 = GaugeMomentum(geo)
    gm0 @= gm
    gm_dual0 @= gm_dual
    change_qnorm = project_gauge_transform(gm0, gm_dual0, mf, mf_dual)
    if change_qnorm > 1e-8:
        q.displayln_info(f"{fname}: project_gauge_transform: change_qnorm={change_qnorm}")
    delta_h = run_hmc_evolve(gm0, gm_dual0, gf0, mf, mf_dual, ga, rs, n_step, md_time)
    if is_reverse_test:
        gm_r = GaugeMomentum(geo)
        gm_dual_r = GaugeMomentum(geo)
        gm_r @= gm0
        gm_dual_r @= gm_dual0
        gf0_r = GaugeField(geo)
        gf0_r @= gf0
        delta_h_rev = run_hmc_evolve(gm_r, gm_dual_r, gf0_r, mf, mf_dual, ga, rs, n_step, -md_time)
        gf0_r -= gf
        gm_r -= gm
        gm_dual_r -= gm_dual
        q.displayln_info(f"{fname}: reversed delta_diff: {delta_h + delta_h_rev} / {delta_h}")
        gm_diff_norm = q.qnorm(gm_r)
        gm_dual_diff_norm = q.qnorm(gm_dual_r)
        gm_norm = q.qnorm(gm)
        gm_dual_norm = q.qnorm(gm_dual)
        gf_diff_norm = q.qnorm(gf0_r)
        gf_norm = q.qnorm(gf)
        q.displayln_info(f"{fname}: reversed gf_diff: {gf_diff_norm} / {gf_norm}")
        q.displayln_info(f"{fname}: reversed gm_diff: {gm_diff_norm} / {gm_norm}")
        q.displayln_info(f"{fname}: reversed gm_dual_diff: {gm_dual_diff_norm} / {gm_dual_norm}")
        assert gf_diff_norm <= 1e-12 * gf_norm
        assert gm_diff_norm <= 1e-12 * gm_norm
        assert gm_dual_diff_norm <= 1e-12 * gm_dual_norm
    flag, accept_prob = metropolis_accept(delta_h, traj, rs.split("metropolis_accept"))
    q.displayln_info(f"{fname}: delta_h={delta_h}, flag={flag}, accept_prob={accept_prob}")
    if flag or is_always_accept:
        q.displayln_info(f"{fname}: update gf (traj={traj})")
        gf @= gf0
        gm @= gm0
        gm_dual @= gm_dual0
    return flag, delta_h

@q.timer_verbose
def run_hmc_mass_mom_refresh(
        job_tag, traj,
        is_force_refresh,
        rs,
        mf, mf_dual,
        gm, gm_dual,
        sel_list,
        sel_dual_list,
        ):
    fname = q.get_fname()
    rs = rs.split(f"{traj}")
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    geo = q.Geometry(total_site)
    complete_refresh_interval = get_param(job_tag, "hmc", "fa", "complete_refresh_interval", default=1)
    mass_type = get_param(job_tag, "hmc", "fa", "mass_type")
    mass_list = get_param(job_tag, "hmc", "fa", "mass_list")
    interval_list = get_param(job_tag, "hmc", "fa", "interval_list")
    is_refresh = traj % complete_refresh_interval == 0
    if is_force_refresh:
        if not is_refresh:
            q.displayln_info(f"{fname}: Force complete refresh of mass and momentum. {job_tag} {traj}")
        is_refresh = True
    if is_refresh:
        # Completely refresh all mass and momentum (and selection)
        if mass_type is None:
            q.set_unit(mf)
            q.set_unit(mf_dual)
            sel = np.zeros(mf[:].shape, dtype=bool)
            sel[:] = True
            sel_list[:] = [
                    sel,
                    ]
            sel_dual_list[:] = [
                    sel,
                    ]
        elif mass_type == "random":
            mf.set_rand(rs.split("fa_mass"), 4.0, 1.0)
            mf_dual.set_rand(rs.split("fa_mass_dual"), 4.0, 1.0)
            sel = mf[:] < 2.0
            sel_dual = mf_dual[:] < 2.0
            sel_list[:] = [
                    sel,
                    ~sel,
                    ]
            sel_dual_list[:] = [
                    sel_dual,
                    ~sel_dual,
                    ]
            assert isinstance(interval_list, list)
            assert len(sel_list) == len(interval_list)
        elif mass_type == "grid-2":
            xg_arr = geo.xg_arr()
            c_offset = rs.c_rand_gen(total_site).to_numpy()
            xg_rel_arr = xg_arr - c_offset
            count = np.zeros(xg_arr.shape, np.int32)
            for m in range(4):
                for n in range(4):
                    if n == m:
                        continue
                    count[xg_rel_arr[:, n] % 2 == 0, m] += 1
            assert np.all(count >= 0)
            assert np.all(count <= 3)
            sel_list[:] = [
                count == 0,
                count == 1,
                count == 2,
                count == 3,
            ]
            sel_dual_list[:] = sel_list
            assert isinstance(interval_list, list)
            assert len(sel_list) == len(interval_list)
            assert isinstance(mass_list, list)
            assert len(mass_list) == len(interval_list)
            for idx, mass in enumerate(mass_list):
                mf[sel_list[idx]] = mass
                mf_dual[sel_dual_list[idx]] = mass
        else:
            raise Exception(f"{fname}: mass_type={mass_type}")
        gm.set_rand_fa(mf, rs.split("set_rand_gauge_momentum"))
        gm_dual.set_rand_fa(mf_dual, rs.split("set_rand_gauge_momentum_dual"))
    else:
        # Only refresh some selection momentum (do not change mass and selection)
        # IMPORTANT: Need to first add the gauge transform component back
        gm1 = GaugeMomentum(geo)
        gm_dual1 = GaugeMomentum(geo)
        gm1.set_rand_fa(mf, rs.split("set_rand_gauge_momentum_g"))
        gm_dual1.set_rand_fa(mf_dual, rs.split("set_rand_gauge_momentum_dual_g"))
        gm2 = gm1.copy()
        gm_dual2 = gm_dual1.copy()
        project_gauge_transform(gm2, gm_dual2, mf, mf_dual)
        gm2 -= gm1
        gm_dual2 -= gm_dual1
        gm -= gm2
        gm_dual -= gm_dual2
        # Then update the selected momentum
        gm1.set_rand_fa(mf, rs.split("set_rand_gauge_momentum"))
        gm_dual1.set_rand_fa(mf_dual, rs.split("set_rand_gauge_momentum_dual"))
        assert len(sel_list) == len(interval_list)
        assert len(sel_dual_list) == len(interval_list)
        for sel, interval in zip(sel_list, interval_list):
            if traj % interval == 0:
                gm[sel] = gm1[sel]
        for sel_dual, interval in zip(sel_dual_list, interval_list):
            if traj % interval == 0:
                gm_dual[sel_dual] = gm_dual1[sel_dual]
    # Project out the gauge transform movement
    project_gauge_transform(gm, gm_dual, mf, mf_dual)

@q.timer_verbose
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

@q.timer_verbose
def run_hmc(job_tag):
    fname = q.get_fname()
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    max_traj = get_param(job_tag, "hmc", "max_traj")
    max_traj_always_accept = get_param(job_tag, "hmc", "max_traj_always_accept")
    max_traj_reverse_test= get_param(job_tag, "hmc", "max_traj_reverse_test")
    save_traj_interval = get_param(job_tag, "hmc", "save_traj_interval")
    is_saving_topo_info = get_param(job_tag, "hmc", "is_saving_topo_info")
    md_time = get_param(job_tag, "hmc", "md_time")
    n_step = get_param(job_tag, "hmc", "n_step")
    beta = get_param(job_tag, "hmc", "beta")
    c1 = get_param(job_tag, "hmc", "c1")
    ga = q.GaugeAction(beta, c1)
    geo = q.Geometry(total_site)
    rs = q.RngState(f"run_hmc-{job_tag}")
    gf = q.GaugeField(geo)
    mf = FieldRealD(geo, 4)
    mf_dual = FieldRealD(geo, 4)
    gm = GaugeMomentum(geo)
    gm_dual = GaugeMomentum(geo)
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
    is_force_refresh = True
    sel_list = []
    sel_dual_list = []
    for traj in range(traj, max_traj):
        run_hmc_mass_mom_refresh(
                job_tag, traj,
                is_force_refresh,
                rs.split("run_hmc_mass_mom_refresh"),
                mf, mf_dual,
                gm, gm_dual,
                sel_list,
                sel_dual_list,
                )
        is_force_refresh = False
        traj += 1
        is_always_accept = traj < max_traj_always_accept
        is_reverse_test = traj < max_traj_reverse_test
        flag, delta_h = run_hmc_traj(
                gf,
                gm, gm_dual,
                mf, mf_dual,
                ga, traj,
                rs.split("run_hmc_traj"),
                n_step=n_step,
                md_time=md_time,
                is_always_accept=is_always_accept,
                is_reverse_test=is_reverse_test,
                )
        plaq = gf.plaq()
        info = dict()
        info["traj"] = traj
        info["plaq"] = plaq
        info["flag"] = flag
        info["delta_h"] = delta_h
        q.qtouch_info(get_save_path(f"{job_tag}/configs/ckpoint_lat_info.{traj}.txt"), pformat(info))
        json_results.append((f"{fname}: {traj} plaq", plaq,))
        if traj % save_traj_interval == 0:
            gf.save(get_save_path(f"{job_tag}/configs/ckpoint_lat.{traj}"))
            if is_saving_topo_info:
                run_topo_info(job_tag, traj, gf)
        q.timer_display()

# ----

job_tag = "test0-4nt8"
set_param(job_tag, "total_site")((4, 4, 4, 8,))
set_param(job_tag, "hmc", "max_traj")(8)
set_param(job_tag, "hmc", "max_traj_always_accept")(4)
set_param(job_tag, "hmc", "max_traj_reverse_test")(2)
set_param(job_tag, "hmc", "md_time")(1.0)
set_param(job_tag, "hmc", "n_step")(10)
set_param(job_tag, "hmc", "beta")(2.13)
set_param(job_tag, "hmc", "c1")(-0.331)
set_param(job_tag, "hmc", "save_traj_interval")(4)
set_param(job_tag, "hmc", "is_saving_topo_info")(True)

job_tag = "test1-4nt8"
set_param(job_tag, "total_site")((4, 4, 4, 8,))
set_param(job_tag, "hmc", "max_traj")(8)
set_param(job_tag, "hmc", "max_traj_always_accept")(4)
set_param(job_tag, "hmc", "max_traj_reverse_test")(2)
set_param(job_tag, "hmc", "md_time")(1.0)
set_param(job_tag, "hmc", "n_step")(10)
set_param(job_tag, "hmc", "beta")(2.13)
set_param(job_tag, "hmc", "c1")(-0.331)
set_param(job_tag, "hmc", "save_traj_interval")(4)
set_param(job_tag, "hmc", "is_saving_topo_info")(True)
set_param(job_tag, "hmc", "fa", "complete_refresh_interval")(2)
set_param(job_tag, "hmc", "fa", "mass_type")("random")
set_param(job_tag, "hmc", "fa", "interval_list")([ 1, 2, ])

job_tag = "test2-4nt8"
set_param(job_tag, "total_site")((4, 4, 4, 8,))
set_param(job_tag, "hmc", "max_traj")(8)
set_param(job_tag, "hmc", "max_traj_always_accept")(4)
set_param(job_tag, "hmc", "max_traj_reverse_test")(2)
set_param(job_tag, "hmc", "md_time")(1.0)
set_param(job_tag, "hmc", "n_step")(10)
set_param(job_tag, "hmc", "beta")(2.13)
set_param(job_tag, "hmc", "c1")(-0.331)
set_param(job_tag, "hmc", "save_traj_interval")(4)
set_param(job_tag, "hmc", "is_saving_topo_info")(True)
set_param(job_tag, "hmc", "fa", "complete_refresh_interval")(2)
set_param(job_tag, "hmc", "fa", "mass_type")("grid-2")
set_param(job_tag, "hmc", "fa", "mass_list")([ 1.0, 1.25, 1.5, 4.0, ])
set_param(job_tag, "hmc", "fa", "interval_list")([ 1, 1, 1, 2, ])

job_tag = "32I_b2p8_fa_sym"
set_param(job_tag, "total_site")((32, 32, 32, 64,))
set_param(job_tag, "a_inv_gev")(2.646) # 2003 lattice spacing 0309017.pdf
set_param(job_tag, "hmc", "max_traj")(20000)
set_param(job_tag, "hmc", "max_traj_always_accept")(100)
set_param(job_tag, "hmc", "max_traj_reverse_test")(2)
set_param(job_tag, "hmc", "md_time")(1.0)
set_param(job_tag, "hmc", "n_step")(32)
set_param(job_tag, "hmc", "beta")(2.80)
set_param(job_tag, "hmc", "c1")(-0.331)
set_param(job_tag, "hmc", "save_traj_interval")(10)
set_param(job_tag, "hmc", "is_saving_topo_info")(True)

job_tag = "32I_b2p8_fa_sym_md2"
set_param(job_tag, "total_site")((32, 32, 32, 64,))
set_param(job_tag, "a_inv_gev")(2.646) # 2003 lattice spacing 0309017.pdf
set_param(job_tag, "hmc", "max_traj")(20000)
set_param(job_tag, "hmc", "max_traj_always_accept")(100)
set_param(job_tag, "hmc", "max_traj_reverse_test")(2)
set_param(job_tag, "hmc", "md_time")(2.0)
set_param(job_tag, "hmc", "n_step")(32 * 2)
set_param(job_tag, "hmc", "beta")(2.80)
set_param(job_tag, "hmc", "c1")(-0.331)
set_param(job_tag, "hmc", "save_traj_interval")(5)
set_param(job_tag, "hmc", "is_saving_topo_info")(True)

job_tag = "32I_b2p8_fa_sym_md3"
set_param(job_tag, "total_site")((32, 32, 32, 64,))
set_param(job_tag, "a_inv_gev")(2.646) # 2003 lattice spacing 0309017.pdf
set_param(job_tag, "hmc", "max_traj")(20000)
set_param(job_tag, "hmc", "max_traj_always_accept")(100)
set_param(job_tag, "hmc", "max_traj_reverse_test")(2)
set_param(job_tag, "hmc", "md_time")(3.0)
set_param(job_tag, "hmc", "n_step")(32 * 3)
set_param(job_tag, "hmc", "beta")(2.80)
set_param(job_tag, "hmc", "c1")(-0.331)
set_param(job_tag, "hmc", "save_traj_interval")(4)
set_param(job_tag, "hmc", "is_saving_topo_info")(True)

job_tag = "32I_b2p8_fa_sym_md4"
set_param(job_tag, "total_site")((32, 32, 32, 64,))
set_param(job_tag, "a_inv_gev")(2.646) # 2003 lattice spacing 0309017.pdf
set_param(job_tag, "hmc", "max_traj")(20000)
set_param(job_tag, "hmc", "max_traj_always_accept")(100)
set_param(job_tag, "hmc", "max_traj_reverse_test")(2)
set_param(job_tag, "hmc", "md_time")(4.0)
set_param(job_tag, "hmc", "n_step")(32 * 4)
set_param(job_tag, "hmc", "beta")(2.80)
set_param(job_tag, "hmc", "c1")(-0.331)
set_param(job_tag, "hmc", "save_traj_interval")(3)
set_param(job_tag, "hmc", "is_saving_topo_info")(True)

job_tag = "32I_b2p8_fa_sym_md5"
set_param(job_tag, "total_site")((32, 32, 32, 64,))
set_param(job_tag, "a_inv_gev")(2.646) # 2003 lattice spacing 0309017.pdf
set_param(job_tag, "hmc", "max_traj")(20000)
set_param(job_tag, "hmc", "max_traj_always_accept")(100)
set_param(job_tag, "hmc", "max_traj_reverse_test")(2)
set_param(job_tag, "hmc", "md_time")(5.0)
set_param(job_tag, "hmc", "n_step")(32 * 5)
set_param(job_tag, "hmc", "beta")(2.80)
set_param(job_tag, "hmc", "c1")(-0.331)
set_param(job_tag, "hmc", "save_traj_interval")(2)
set_param(job_tag, "hmc", "is_saving_topo_info")(True)

job_tag = "32I_b2p8_fa_sym_v1"
set_param(job_tag, "total_site")((32, 32, 32, 64,))
set_param(job_tag, "a_inv_gev")(2.646) # 2003 lattice spacing 0309017.pdf
set_param(job_tag, "hmc", "max_traj")(20000)
set_param(job_tag, "hmc", "max_traj_always_accept")(100)
set_param(job_tag, "hmc", "max_traj_reverse_test")(2)
set_param(job_tag, "hmc", "md_time")(1.0)
set_param(job_tag, "hmc", "n_step")(32)
set_param(job_tag, "hmc", "beta")(2.80)
set_param(job_tag, "hmc", "c1")(-0.331)
set_param(job_tag, "hmc", "save_traj_interval")(10)
set_param(job_tag, "hmc", "is_saving_topo_info")(True)
set_param(job_tag, "hmc", "fa", "mass_type")("grid-2")
set_param(job_tag, "hmc", "fa", "mass_list")([ 1.0, 1.25, 1.5, 4.0, ])
set_param(job_tag, "hmc", "fa", "interval_list")([ 1, 1, 1, 2, ])

job_tag = "32I_b2p8_fa_sym_v2"
set_param(job_tag, "total_site")((32, 32, 32, 64,))
set_param(job_tag, "a_inv_gev")(2.646) # 2003 lattice spacing 0309017.pdf
set_param(job_tag, "hmc", "max_traj")(20000)
set_param(job_tag, "hmc", "max_traj_always_accept")(100)
set_param(job_tag, "hmc", "max_traj_reverse_test")(2)
set_param(job_tag, "hmc", "md_time")(2.0)
set_param(job_tag, "hmc", "n_step")(32)
set_param(job_tag, "hmc", "beta")(2.80)
set_param(job_tag, "hmc", "c1")(-0.331)
set_param(job_tag, "hmc", "save_traj_interval")(10)
set_param(job_tag, "hmc", "is_saving_topo_info")(True)
set_param(job_tag, "hmc", "fa", "mass_type")("grid-2")
set_param(job_tag, "hmc", "fa", "mass_list")([ 1.0, 1.25, 1.5, 4.0, ])
set_param(job_tag, "hmc", "fa", "interval_list")([ 1, 1, 1, 2, ])

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

    job_tags = q.get_arg("--job_tags", default="").split(",")

    #######################################################

    job_tags_default = [
            "test0-4nt8",
            "test1-4nt8",
            "test2-4nt8",
            ]

    if job_tags == [ "", ]:
        job_tags = job_tags_default

    for job_tag in job_tags:
        run_params(job_tag)
        run_hmc(job_tag)

    q.check_log_json(__file__, json_results)

    q.timer_display()

    q.end_with_mpi()

    q.displayln_info(f"CHECK: finished successfully.")

# ----
