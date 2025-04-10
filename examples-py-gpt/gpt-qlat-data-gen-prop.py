#!/usr/bin/env python3

json_results = []
check_eps = 5e-5

def json_results_append(*args):
    q.displayln_info(r"//------------------------------------------------------------\\")
    q.displayln_info(-1, *args)
    q.displayln_info(r"\\------------------------------------------------------------//")
    json_results.append(args)

from auto_contractor import (
        mk_meson,
        mk_fac,
        mk_vec5_mu,
        show_dagger,
        contract_simplify_compile,
        cache_compiled_cexpr,
        benchmark_eval_cexpr,
        get_expr_names,
        eval_cexpr,
        )

from qlat_scripts.v1 import (
        load_path_list,
        get_param,
        set_param,
        get_load_path,
        get_save_path,
        check_job,
        run_params,
        run_gf,
        run_gt,
        run_eig,
        run_f_weight_uniform,
        run_f_rand_01,
        run_fsel_prob,
        run_psel_prob,
        run_fsel_from_fsel_prob,
        run_psel_from_psel_prob,
        )

import functools
import math
import os
import time
import importlib
import sys

import qlat as q
import qlat_gpt as qg
import qlat_scripts.v1 as qs
import numpy as np

is_cython = False

### ------

load_path_list[:] = [
        "results",
        "qcddata",
        "/lustre20/volatile/qcdqedta/qcddata",
        "/lustre20/volatile/decay0n2b/qcddata",
        "/lustre20/volatile/pqpdf/ljin/qcddata",
        "/lustre20/volatile/decay0n2b/qcddata/qcddata4",
        "/lustre20/volatile/decay0n2b/qcddata/qcddata3",
        "/lustre20/volatile/decay0n2b/qcddata/qcddata1",
        "/data1/qcddata4",
        "/data1/qcddata3",
        "/data2/qcddata3-prop",
        "/data1/qcddata1",
        ]

### ------

@q.timer_verbose
def run_prop_rand_vol_u1_src(
        job_tag, traj,
        *,
        inv_type,
        get_gf,
        get_psel,
        get_fsel,
        get_eig=None,
        ):
    """
    fsel should contain psel
    likely strange quark mass is included as well.
    """
    fname = q.get_fname()
    inv_type_name_list = [ "light", "strange", ] + [ f"charm-{idx}" for idx in range(8) ]
    quark_flavor_list = get_param(job_tag, "quark_flavor_list")
    quark_flavor = quark_flavor_list[inv_type]
    path_s = f"{job_tag}/prop-rand-vol-u1-{quark_flavor}/traj-{traj}"
    path_sp = f"{job_tag}/psel-prop-rand-vol-u1-{quark_flavor}/traj-{traj}"
    if get_load_path(path_s + "/geon-info.txt") is not None:
        assert get_load_path(path_sp + "/checkpoint.txt") is not None
        return
    gf = get_gf()
    geo = gf.geo
    fsel = get_fsel()
    psel = get_psel()
    total_site = geo.total_site
    if get_eig is None:
        eig = None
    else:
        eig = get_eig()
    gt = None
    num_prop_rand_vol_u1 = get_param(job_tag, "num_prop_rand_vol_u1")
    prob_acc_1_rand_vol_u1 = get_param(job_tag, "prob_acc_1_rand_vol_u1")
    prob_acc_2_rand_vol_u1 = get_param(job_tag, "prob_acc_2_rand_vol_u1")
    sfw = q.open_fields(get_save_path(path_s + ".acc"), "a", q.Coordinate([ 2, 2, 2, 4, ]))
    qar_sp = q.open_qar_info(get_save_path(path_sp + ".qar"), "a")
    rs_rand_u1 = q.RngState(f"seed {job_tag} {traj}").split(f"compute_prop_rand_vol_u1_src(rand_u1)")
    rs_ama = q.RngState(f"seed {job_tag} {traj}").split(f"compute_prop_rand_u1(ama)")
    @q.timer
    def compute_and_save(idx_rand_vol_u1, inv_acc):
        tag = f"idx_rand_vol_u1={idx_rand_vol_u1} ; type={inv_type} ; accuracy={inv_acc}"
        if tag in sfw:
            assert f"{tag} ; fu1" in sfw
            if q.get_id_node() == 0:
                assert f"{tag}.lat" in qar_sp
            return
        if q.get_id_node() == 0:
            if f"{tag}.lat" in qar_sp:
                q.displayln_info(-1, f"WARNING: {fname}: '{tag}.lat' already exist in {qar_sp.path()}")
        if f"{tag} ; fu1" in sfw:
            q.displayln_info(-1, f"WARNING: {fname}: '{tag} ; fu1' already exist in {sfw.path()}")
        rsi = rs_rand_u1.split(str(idx_rand_vol_u1))
        prop_src, fu1 = q.mk_rand_vol_u1_src(geo, rsi)
        inv = qs.get_inv(gf, job_tag, inv_type, inv_acc, gt=gt, eig=eig)
        prop_sol = inv * prop_src
        s_prop = q.SelProp(fsel)
        ps_prop = q.PselProp(psel)
        s_prop @= prop_sol
        ps_prop @= prop_sol
        qar_sp.write(f"{tag}.lat", "", ps_prop.save_str(), skip_if_exist=True)
        fu1.save_float_from_double(sfw, f"{tag} ; fu1", skip_if_exist=True)
        s_prop.save_float_from_double(sfw, tag, skip_if_exist=True)
        qar_sp.flush()
        sfw.flush()
    for idx_rand_vol_u1 in range(num_prop_rand_vol_u1):
        r = rs_ama.split(str(idx_rand_vol_u1)).u_rand_gen()
        inv_acc = 0
        assert 0 <= r and r <= 1
        compute_and_save(idx_rand_vol_u1, inv_acc)
        inv_acc = 1
        if r <= prob_acc_1_rand_vol_u1:
            compute_and_save(idx_rand_vol_u1, inv_acc)
        inv_acc = 2
        if r <= prob_acc_2_rand_vol_u1:
            compute_and_save(idx_rand_vol_u1, inv_acc)
    sfw.close()
    qar_sp.write("checkpoint.txt", "", "", skip_if_exist=True)
    qar_sp.flush()
    qar_sp.close()
    q.clean_cache(q.cache_inv)

### ------

def get_param_quark_mass_list(job_tag):
    quark_mass_list = get_param(job_tag, "quark_mass_list")
    return quark_mass_list

@q.timer_verbose
def run_quark_mass_list(job_tag, traj):
    """
    return quark_mass_list
    """
    quark_mass_list = get_param_quark_mass_list(job_tag)
    for inv_type, quark_mass in enumerate(quark_mass_list):
        assert quark_mass == get_param(job_tag, "fermion_params", inv_type, 0, "mass")
    obj = quark_mass_list
    fn = f"{job_tag}/params/traj-{traj}/quark_mass_list.json"
    path = get_load_path(fn)
    if path is None:
        path = get_save_path(fn)
        q.save_json_obj(obj, path)
    else:
        obj_load = q.load_json_obj(path)
        assert obj_load is not None
        assert obj_load == obj
    return obj

### ------

@q.timer(is_timer_fork=True)
def run_job(job_tag, traj):
    #
    traj_gf = traj
    #
    if job_tag[:5] == "test-":
        traj_gf = 1000
    #
    fns_produce = [
            f"{job_tag}/auto-contract/traj-{traj}/checkpoint.txt",
            ]
    fns_need = [
            (f"{job_tag}/configs/ckpoint_lat.{traj_gf}", f"{job_tag}/configs/ckpoint_lat.IEEE64BIG.{traj_gf}",),
            f"{job_tag}/eig/traj-{traj}/metadata.txt",
            f"{job_tag}/gauge-transform/traj-{traj_gf}.field",
            f"{job_tag}/point-selection/traj-{traj}.txt",
            f"{job_tag}/field-selection/traj-{traj}.field",
            ]
    #
    if job_tag[:5] == "test-":
        fns_need = []
    #
    if not check_job(job_tag, traj, fns_produce, fns_need):
        return
    #
    get_gf = run_gf(job_tag, traj_gf)
    get_gt = run_gt(job_tag, traj_gf, get_gf)
    get_eig_light = run_eig(job_tag, traj_gf, get_gf)
    #
    get_f_weight = run_f_weight_uniform(job_tag, traj)
    get_f_rand_01 = run_f_rand_01(job_tag, traj)
    get_fsel_prob = run_fsel_prob(job_tag, traj, get_f_rand_01=get_f_rand_01, get_f_weight=get_f_weight)
    get_psel_prob = run_psel_prob(job_tag, traj, get_f_rand_01=get_f_rand_01, get_f_weight=get_f_weight)
    get_fsel = run_fsel_from_fsel_prob(get_fsel_prob)
    get_psel = run_psel_from_psel_prob(get_psel_prob)
    #
    quark_mass_list = run_quark_mass_list(job_tag, traj)
    #
    fn_checkpoint = f"{job_tag}/auto-contract/traj-{traj}/checkpoint.txt"
    if get_load_path(fn_checkpoint) is not None:
        return
    if not q.obtain_lock(f"locks/{job_tag}-{traj}-auto-contract"):
        return
    for inv_type, quark_mass in enumerate(quark_mass_list):
        if inv_type == 0:
            get_eig = get_eig_light
        else:
            get_eig = None
        run_prop_rand_vol_u1_src(
                job_tag, traj,
                inv_type=inv_type,
                get_gf=get_gf,
                get_psel=get_psel,
                get_fsel=get_fsel,
                get_eig=get_eig,
                )
    q.qtouch_info(get_save_path(fn_checkpoint))
    q.release_lock()

### ------

def get_all_cexpr():
    ...

### ------

set_param("16IH2", "trajs")(list(range(1000, 4200, 100)))
set_param("16IH2", "quark_mass_list")([ 0.01, 0.04, 0.02963, 0.05358, 0.07945, 0.10852, ])
for inv_type in range(1, len(get_param("16IH2", "quark_mass_list"))):
    mass = get_param("16IH2", "quark_mass_list")[inv_type]
    for inv_acc in [ 0, 1, 2, ]:
        set_param("16IH2", 'fermion_params', inv_type, inv_acc)(
            {
                'M5': 1.8,
                'boundary_phases': [1.0, 1.0, 1.0, 1.0],
                'b': 1.0,
                'c': 0.0,
                'Ls': 16,
                'mass': mass,
            }
        )
        set_param("16IH2", f"cg_params-{inv_type}-{inv_acc}", "maxiter")(300)
    set_param("16IH2", f"cg_params-{inv_type}-0", "maxcycle")(1)
    set_param("16IH2", f"cg_params-{inv_type}-1", "maxcycle")(2)
    set_param("16IH2", f"cg_params-{inv_type}-2", "maxcycle")(50)
set_param("16IH2", f"cg_params-0-0", "maxiter")(200)
set_param("16IH2", f"cg_params-0-0", "maxcycle")(1)
set_param("16IH2", f"cg_params-0-1", "maxiter")(200)
set_param("16IH2", f"cg_params-0-1", "maxcycle")(2)
set_param("16IH2", f"cg_params-0-2", "maxiter")(300)
set_param("16IH2", f"cg_params-0-2", "maxcycle")(50)
if 1 in get_param("16IH2", 'lanc_params'):
    get_param("16IH2", 'lanc_params').pop(1)
if 1 in get_param("16IH2", 'clanc_params'):
    get_param("16IH2", 'clanc_params').pop(1)
set_param("16IH2", "num_prop_rand_vol_u1")(256)
set_param("16IH2", "prob_acc_1_rand_vol_u1")(1 / 32)
set_param("16IH2", "prob_acc_2_rand_vol_u1")(1 / 128)
set_param("16IH2", "measurement", "auto_contractor_chunk_size")(128)

set_param("24D", "trajs")(list(range(1000, 5100, 10)))
set_param("24D", "quark_mass_list")([ 0.00107, 0.085, 0.07819, 0.13207, 0.19829, ])
for inv_type in range(1, len(get_param("24D", "quark_mass_list"))):
    mass = get_param("24D", "quark_mass_list")[inv_type]
    for inv_acc in [ 0, 1, 2, ]:
        set_param("24D", 'fermion_params', inv_type, inv_acc)(
            {
                'M5': 1.8,
                'boundary_phases': [1.0, 1.0, 1.0, 1.0],
                'b': 2.5,
                'c': 1.5,
                'mass': mass,
                'Ls': 24,
            }
        )
        set_param("24D", f"cg_params-{inv_type}-{inv_acc}", "maxiter")(300)
    set_param("24D", f"cg_params-{inv_type}-0", "maxcycle")(1)
    set_param("24D", f"cg_params-{inv_type}-1", "maxcycle")(2)
    set_param("24D", f"cg_params-{inv_type}-2", "maxcycle")(50)
set_param("24D", f"cg_params-0-0", "maxiter")(200)
set_param("24D", f"cg_params-0-0", "maxcycle")(1)
set_param("24D", f"cg_params-0-1", "maxiter")(200)
set_param("24D", f"cg_params-0-1", "maxcycle")(2)
set_param("24D", f"cg_params-0-2", "maxiter")(300)
set_param("24D", f"cg_params-0-2", "maxcycle")(50)
if 1 in get_param("24D", 'lanc_params'):
    get_param("24D", 'lanc_params').pop(1)
if 1 in get_param("24D", 'clanc_params'):
    get_param("24D", 'clanc_params').pop(1)
set_param("24D", "num_prop_rand_vol_u1")(256)
set_param("24D", "prob_acc_1_rand_vol_u1")(1 / 32)
set_param("24D", "prob_acc_2_rand_vol_u1")(1 / 128)
set_param("24D", "measurement", "auto_contractor_chunk_size")(128)

set_param("32Dfine", "trajs")(list(range(1000, 5100, 10)))
set_param("32Dfine", "quark_mass_list")([ 0.0001, 0.045, 0.04635, 0.07794, 0.11333, 0.15327, ])
for inv_type in range(1, len(get_param("32Dfine", "quark_mass_list"))):
    mass = get_param("32Dfine", "quark_mass_list")[inv_type]
    for inv_acc in [ 0, 1, 2, ]:
        set_param("32Dfine", 'fermion_params', inv_type, inv_acc)(
            {
                'M5': 1.8,
                'boundary_phases': [1.0, 1.0, 1.0, 1.0],
                'b': 1.8333333333333333,
                'c': 0.8333333333333333,
                'mass': mass,
                'Ls': 12,
            }
        )
        set_param("32Dfine", f"cg_params-{inv_type}-{inv_acc}", "maxiter")(300)
    set_param("32Dfine", f"cg_params-{inv_type}-0", "maxcycle")(1)
    set_param("32Dfine", f"cg_params-{inv_type}-1", "maxcycle")(2)
    set_param("32Dfine", f"cg_params-{inv_type}-2", "maxcycle")(50)
set_param("32Dfine", f"cg_params-0-0", "maxiter")(200)
set_param("32Dfine", f"cg_params-0-0", "maxcycle")(1)
set_param("32Dfine", f"cg_params-0-1", "maxiter")(200)
set_param("32Dfine", f"cg_params-0-1", "maxcycle")(2)
set_param("32Dfine", f"cg_params-0-2", "maxiter")(300)
set_param("32Dfine", f"cg_params-0-2", "maxcycle")(50)
if 1 in get_param("32Dfine", 'lanc_params'):
    get_param("32Dfine", 'lanc_params').pop(1)
if 1 in get_param("32Dfine", 'clanc_params'):
    get_param("32Dfine", 'clanc_params').pop(1)
set_param("32Dfine", "num_prop_rand_vol_u1")(256)
set_param("32Dfine", "prob_acc_1_rand_vol_u1")(1 / 32)
set_param("32Dfine", "prob_acc_2_rand_vol_u1")(1 / 128)
set_param("32Dfine", "measurement", "auto_contractor_chunk_size")(128)

# ----

job_tag = "test-4nt8-checker"
#
set_param(job_tag, "trajs")([ 1000, ])
#
set_param(job_tag, "total_site")([ 4, 4, 4, 8, ])
set_param(job_tag, "load_config_params", "twist_boundary_at_boundary")([ 0.0, 0.0, 0.0, -0.5, ])
#
set_param(job_tag, "mk_sample_gauge_field", "rand_n_step")(2)
set_param(job_tag, "mk_sample_gauge_field", "flow_n_step")(8)
set_param(job_tag, "mk_sample_gauge_field", "hmc_n_traj")(5)
#
set_param(job_tag, "quark_mass_list")([ 0.01, 0.04, 0.1, 0.2, ])
set_param(job_tag, "quark_flavor_list")([ "light", "strange", "charm-1", "charm-2", ])
set_param(job_tag, "fermion_params", 0, 0)({ 'Ls': 8, 'M5': 1.8, 'b': 1.5, 'c': 0.5, 'boundary_phases': [1.0, 1.0, 1.0, 1.0], })
for inv_type, mass in enumerate(get_param(job_tag, "quark_mass_list")):
    set_param(job_tag, "fermion_params", inv_type, 0)(get_param(job_tag, "fermion_params", 0, 0).copy())
    set_param(job_tag, "fermion_params", inv_type, 0, "mass")(mass)
    for inv_acc in [ 0, 1, 2, ]:
        set_param(job_tag, "fermion_params", inv_type, inv_acc)(get_param(job_tag, "fermion_params", inv_type, 0).copy())
        set_param(job_tag, f"cg_params-{inv_type}-{inv_acc}", "maxiter")(10)
        set_param(job_tag, f"cg_params-{inv_type}-{inv_acc}", "maxcycle")(1 + inv_acc)
#
set_param(job_tag, "lanc_params", 0, 0, "cheby_params")({ "low": 0.5, "high": 5.5, "order": 40, })
set_param(job_tag, "lanc_params", 0, 0, "irl_params")({ "Nstop": 100, "Nk": 150, "Nm": 200, "resid": 1e-8, "betastp": 0.0, "maxiter": 20, "Nminres": 0, })
set_param(job_tag, "lanc_params", 0, 0, "pit_params")({ 'eps': 0.01, 'maxiter': 500, 'real': True, })
set_param(job_tag, "lanc_params", 0, 0, "fermion_params")(get_param(job_tag, "fermion_params", inv_type, 0).copy())
#
set_param(job_tag, "clanc_params", 0, 0, "nbasis")(100)
set_param(job_tag, "clanc_params", 0, 0, "block")([ 4, 4, 2, 2, ])
set_param(job_tag, "clanc_params", 0, 0, "cheby_params")({ "low": 0.5, "high": 5.5, "order": 40, })
set_param(job_tag, "clanc_params", 0, 0, "save_params")({ "nsingle": 100, "mpi": [ 1, 1, 1, 4, ], })
set_param(job_tag, "clanc_params", 0, 0, "irl_params")({ "Nstop": 100, "Nk": 150, "Nm": 200, "resid": 1e-8, "betastp": 0.0, "maxiter": 20, "Nminres": 0, })
set_param(job_tag, "clanc_params", 0, 0, "smoother_params")({'eps': 1e-08, 'maxiter': 20})
#
set_param(job_tag, "field-selection-fsel-rate")(0.1)
set_param(job_tag, "field-selection-psel-rate")(0.01)
set_param(job_tag, "field-selection-fsel-psrc-prop-norm-threshold")(0.05)
#
set_param(job_tag, "prob_exact_wsrc")(0.20)
#
set_param(job_tag, "prob_acc_1_psrc")(0.25)
set_param(job_tag, "prob_acc_2_psrc")(0.10)
#
set_param(job_tag, "num_prop_rand_vol_u1")(4)
set_param(job_tag, "prob_acc_1_rand_vol_u1")(0.25)
set_param(job_tag, "prob_acc_2_rand_vol_u1")(0.10)
#
set_param(job_tag, "measurement", "auto_contractor_chunk_size")(2)

# ----

def gracefully_finish():
    q.displayln_info("Begin to gracefully_finish.")
    q.timer_display()
    qg.end_with_gpt()
    q.displayln_info("CHECK: finished successfully.")
    exit()

if __name__ == "__main__":

    qg.begin_with_gpt()

    ##################### CMD options #####################

    job_tags = q.get_arg("--job_tags", default="").split(",")

    #######################################################

    job_tags_default = [
            "test-4nt8-checker",
            ]

    if job_tags == [ "", ]:
        job_tags = job_tags_default
    else:
        is_cython = True

    q.check_time_limit()

    get_all_cexpr()

    for job_tag in job_tags:
        run_params(job_tag)
        for traj in get_param(job_tag, "trajs"):
            q.check_time_limit()
            run_job(job_tag, traj)
            q.clean_cache()
            if q.obtained_lock_history_list:
                json_results_append(f"q.obtained_lock_history_list={q.obtained_lock_history_list}")
                if job_tag[:5] != "test-":
                    gracefully_finish()

    q.check_log_json(__file__, json_results, check_eps=check_eps)

    gracefully_finish()

# ----
