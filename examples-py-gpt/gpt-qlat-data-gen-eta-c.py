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
        run_params,
        check_job,
        run_gf,
        run_gt,
        get_load_path,
        get_save_path,
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

def mk_eta_c(p:str, is_dagger=False):
    """
    cbar g5 c  #dag: same
    """
    return mk_meson("c", "c", p, is_dagger) + f"eta_c({p}){show_dagger(is_dagger)}"

def mk_j5_eta_c_mu(p:str, mu, is_dagger=False):
    return mk_vec5_mu("c", "c", p, mu, is_dagger) + f"j5_eta_c_mu({p},{mu}){show_dagger(is_dagger)}"

@q.timer
def get_cexpr_eta_c_corr():
    fn_base = "cache/auto_contract_cexpr/get_cexpr_eta_c_corr"
    def calc_cexpr():
        diagram_type_dict = dict()
        diagram_type_dict[((('x_1', 'x_2'), 1), (('x_2', 'x_1'), 1))] = 'Type1'
        diagram_type_dict[((('x_1', 'x_1'), 1), (('x_2', 'x_2'), 1))] = None
        exprs = [
                mk_fac(1) + f"1",
                mk_eta_c("x_2", True)    * mk_eta_c("x_1")                + f"eta_c^dag(0) * eta_c(-tsep)",
                mk_j5_eta_c_mu("x_2", 3) * mk_eta_c("x_1")                + f"j5_eta_c_t(0) * eta_c(-tsep)",
                mk_eta_c("x_2", True)    * mk_j5_eta_c_mu("x_1", 3, True) + f"eta_c^dag(0) * j5_eta_c_t^dag(-tsep)",
                mk_j5_eta_c_mu("x_2", 3) * mk_j5_eta_c_mu("x_1", 3, True) + f"j5_eta_c_t(0) * j5_eta_c_t^dag(-tsep)",
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit=True, diagram_type_dict=diagram_type_dict)
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

### ------

@q.timer_verbose
def auto_contract_eta_c_corr(job_tag, traj, get_get_prop, charm_mass_idx, tslice_list):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract/traj-{traj}/eta_c/charm_mass_idx-{charm_mass_idx}.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_eta_c_corr()
    expr_names = get_expr_names(cexpr)
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    t_size = total_site[3]
    geo = q.Geometry(total_site)
    get_prop = get_get_prop()
    def load_data():
        t_t_list = q.get_mpi_chunk(
                [ (t_src, t_snk,) for t_snk in range(total_site[3]) for t_src in tslice_list ],
                rng_state=None)
        for t_src, t_snk in t_t_list:
            yield t_src, t_snk
    @q.timer
    def feval(args):
        t_src, t_snk = args
        t = (t_snk - t_src) % total_site[3]
        pd = {
                "x_2" : ("wall", t_snk,),
                "x_1" : ("wall", t_src,),
                "size" : total_site,
                }
        val = eval_cexpr(cexpr, positions_dict=pd, get_prop=get_prop)
        return val, t
    def sum_function(val_list):
        values = np.zeros((total_site[3], len(expr_names),), dtype=complex)
        for val, t in val_list:
            values[t] += val
        return q.glb_sum(values.transpose(1, 0))
    auto_contractor_chunk_size = get_param(job_tag, "measurement", "auto_contractor_chunk_size", default=128)
    with q.TimerFork(max_call_times_for_always_show_info=0):
        res_sum = q.parallel_map_sum(feval, load_data(), sum_function=sum_function, chunksize=auto_contractor_chunk_size)
        q.displayln_info(f"{fname}: timer_display for parallel_map_sum")
    res_sum *= 1.0 / len(tslice_list)
    assert q.qnorm(res_sum[0] - 1.0) < 1e-10
    ld = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t_sep", t_size, [ str(q.rel_mod(t, t_size)) for t in range(t_size) ], ],
        ])
    ld.from_numpy(res_sum)
    ld.save(get_save_path(fn))
    json_results_append(f"{fname}: ld sig", q.get_data_sig(ld, q.RngState()))
    for i, en in enumerate(expr_names):
        json_results_append(f"{fname}: ld '{en}' sig", q.get_data_sig(ld[i], q.RngState()))

@q.timer_verbose
def auto_contract_eta_c_corr_psnk(job_tag, traj, get_get_prop, charm_mass_idx, tslice_list):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract/traj-{traj}/eta_c-psnk/charm_mass_idx-{charm_mass_idx}.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_eta_c_corr()
    expr_names = get_expr_names(cexpr)
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    t_size = total_site[3]
    geo = q.Geometry(total_site)
    xg_arr = geo.xg_arr()
    xg_arr_list = [ xg_arr[xg_arr[:, 3] == t] for t in range(t_size) ]
    spatial_volume = geo.total_volume / t_size
    get_prop = get_get_prop()
    def load_data():
        t_t_list = q.get_mpi_chunk(
                [ (t_src, t_snk,) for t_snk in range(total_site[3]) for t_src in tslice_list ],
                rng_state=None)
        for t_src, t_snk in t_t_list:
            yield t_src, t_snk
    @q.timer
    def feval(args):
        t_src, t_snk = args
        t = (t_snk - t_src) % total_site[3]
        xg_arr_t_snk = xg_arr_list[t_snk]
        val = 0
        for xg_snk in xg_arr_t_snk:
            xg_snk = tuple(xg_snk)
            pd = {
                    "x_2" : ("point-snk", xg_snk,),
                    "x_1" : ("wall", t_src,),
                    "size" : total_site,
                    }
            val += eval_cexpr(cexpr, positions_dict=pd, get_prop=get_prop)
        return val, t
    def sum_function(val_list):
        values = np.zeros((total_site[3], len(expr_names),), dtype=complex)
        for val, t in val_list:
            values[t] += val
        return q.glb_sum(values.transpose(1, 0))
    auto_contractor_chunk_size = get_param(job_tag, "measurement", "auto_contractor_chunk_size", default=128)
    with q.TimerFork(max_call_times_for_always_show_info=0):
        res_sum = q.parallel_map_sum(feval, load_data(), sum_function=sum_function, chunksize=auto_contractor_chunk_size)
        q.displayln_info(f"{fname}: timer_display for parallel_map_sum")
    res_sum *= 1.0 / len(tslice_list) / spatial_volume
    assert q.qnorm(res_sum[0] - 1.0) < 1e-10
    ld = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t_sep", t_size, [ str(q.rel_mod(t, t_size)) for t in range(t_size) ], ],
        ])
    ld.from_numpy(res_sum)
    ld.save(get_save_path(fn))
    json_results_append(f"{fname}: ld sig", q.get_data_sig(ld, q.RngState()))
    for i, en in enumerate(expr_names):
        json_results_append(f"{fname}: ld '{en}' sig", q.get_data_sig(ld[i], q.RngState()))

### ------

@q.timer_verbose
def mk_get_prop(prop_dict):
    """
    return get_prop
    #
    key = (flavor, pos_src, type_src, type_snk,)
    prop_dict[key](pos_snk) => wm or AmaVal(wm)
    #
    e.g. prop_dict = q.mk_cache("prop-dict")
    """
    def get_prop(flavor, p_snk, p_src):
        assert isinstance(p_snk, tuple) and isinstance(p_src, tuple)
        type_snk, pos_snk = p_snk
        type_src, pos_src = p_src
        key = (flavor, pos_src, type_src, type_snk,)
        get = prop_dict.get(key)
        if get is not None:
            return get(p_snk)
        key = (flavor, pos_snk, type_snk, type_src,)
        get = prop_dict.get(key)
        if get is not None:
            return ("g5_herm", get(p_src),)
        fname = q.get_fname()
        raise Exception(f"{fname}: {flavor} {p_snk} {p_src}")
    return get_prop

@q.timer_verbose
def run_get_prop_wsrc_charm(job_tag, traj, *, get_gf, get_gt, charm_mass, tslice_list):
    """
    return get_get_prop
    """
    @q.timer_verbose
    def get_get_prop():
        """
        return get_prop
        """
        gf = get_gf()
        gt = get_gt()
        gt_inv = gt.inv()
        geo = gf.geo
        total_site = geo.total_site
        inv_type = 2
        inv_acc = 2
        eig = None
        mass_initial = get_param(job_tag, "fermion_params", inv_type, inv_acc, "mass")
        set_param(job_tag, "fermion_params", inv_type, inv_acc, "mass")(charm_mass)
        inv = qs.get_inv(gf, job_tag, inv_type, inv_acc, gt=gt, eig=eig)
        prop_dict = dict()
        for tslice in tslice_list:
            src = q.mk_wall_src(geo, tslice)
            sol = inv * src
            prop = sol
            ps_prop_ws = prop.glb_sum_tslice()
            prop_ps = gt_inv * prop
            def get_prop_point_snk(p_snk):
                type_snk, pos_snk = p_snk
                assert type_snk == "point-snk"
                assert isinstance(pos_snk, tuple)
                pos_snk = q.Coordinate(pos_snk)
                index = geo.index_from_g_coordinate(pos_snk)
                return prop_ps.get_elem_wm(index)
            def get_prop_wall_snk(p_snk):
                type_snk, pos_snk = p_snk
                assert type_snk == "wall"
                assert isinstance(pos_snk, int)
                return ps_prop_ws.get_elem_wm(pos_snk)
            key = ("c", tslice, "wall", "point-snk")
            prop_dict[key] = get_prop_point_snk
            key = ("c", tslice, "wall", "wall")
            prop_dict[key] = get_prop_wall_snk
        q.clean_cache(q.cache_inv)
        set_param(job_tag, "fermion_params", inv_type, inv_acc, "mass")(mass_initial)
        prop_cache = q.mk_cache("prop_cache", f"{job_tag}", f"{traj}")
        prop_cache["prop-dict"] = prop_dict
        return mk_get_prop(prop_dict)
    return q.lazy_call(get_get_prop)

### ------

@q.timer(is_timer_fork=True)
def run_eta_c_corr(job_tag, traj, get_gf, get_gt):
    charm_quark_mass_list = get_param_charm_mass_list(job_tag)
    tslice_list = get_param_charm_wall_src_tslice_list(job_tag, traj)
    for charm_mass_idx, charm_mass in enumerate(charm_quark_mass_list):
        get_get_prop = run_get_prop_wsrc_charm(
                job_tag, traj,
                get_gf=get_gf,
                get_gt=get_gt,
                charm_mass=charm_mass,
                tslice_list=tslice_list,
                )
        auto_contract_eta_c_corr(job_tag, traj, get_get_prop, charm_mass_idx, tslice_list)
        auto_contract_eta_c_corr_psnk(job_tag, traj, get_get_prop, charm_mass_idx, tslice_list)

### ------

def get_param_charm_mass_list(job_tag):
    charm_quark_mass_list = get_param(job_tag, "measurement", "charm_quark_mass_list")
    return charm_quark_mass_list

@q.cache_call()
@q.timer_verbose
def get_param_charm_wall_src_tslice_list(job_tag, traj):
    num_charm_wall_src = get_param(job_tag, "measurement", "num_charm_wall_src")
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    t_size = total_site[3]
    rs = q.RngState(f"{job_tag}-{traj}-get_param_charm_wall_src_tslice_list")
    t_start = rs.u_rand_gen() * t_size
    t_sep = t_size / num_charm_wall_src
    tslice_list = [ round(t_start + i * t_sep) % t_size for i in range(num_charm_wall_src) ]
    return tslice_list

@q.timer_verbose
def run_charm_wall_src_prop_params(job_tag, traj):
    charm_wall_src_tslice_list = get_param_charm_wall_src_tslice_list(job_tag, traj)
    charm_quark_mass_list = get_param_charm_mass_list(job_tag)
    obj = dict()
    obj["charm_wall_src_tslice_list"] = charm_wall_src_tslice_list
    obj["charm_quark_mass_list"] = charm_quark_mass_list
    fn = f"{job_tag}/params/traj-{traj}/charm_wall_src_prop.json"
    path = get_load_path(fn)
    if path is None:
        path = get_save_path(fn)
        q.save_json_obj(obj, path)
    else:
        obj_load = q.load_json_obj(path)
        assert obj_load is not None
        assert obj_load == obj

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
            f"{job_tag}/gauge-transform/traj-{traj_gf}.field",
            (f"{job_tag}/configs/ckpoint_lat.{traj_gf}", f"{job_tag}/configs/ckpoint_lat.IEEE64BIG.{traj_gf}",),
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
    #
    run_charm_wall_src_prop_params(job_tag, traj)
    #
    fn_checkpoint = f"{job_tag}/auto-contract/traj-{traj}/checkpoint.txt"
    if get_load_path(fn_checkpoint) is not None:
        return
    if not q.obtain_lock(f"locks/{job_tag}-{traj}-auto-contract"):
        return
    run_eta_c_corr(job_tag, traj, get_gf, get_gt)
    q.qtouch_info(get_save_path(fn_checkpoint))
    q.release_lock()

### ------

def get_all_cexpr():
    benchmark_eval_cexpr(get_cexpr_eta_c_corr())

### ------

set_param("16IH2", "trajs")(list(range(1000, 4020, 10)))
set_param("16IH2", "measurement", "auto_contractor_chunk_size")(128)
set_param("16IH2", "measurement", "charm_quark_mass_list")([ 0.04, 0.05, 0.1, 0.15, 0.2, ])
set_param("16IH2", "measurement", "num_charm_wall_src")(2)
set_param("16IH2", f"cg_params-2-2", "maxiter")(200)
set_param("16IH2", f"cg_params-2-2", "maxcycle")(50)

set_param("24D", "trajs")(list(range(1000, 5100, 10)))
set_param("24D", "measurement", "charm_quark_mass_list")([ 0.0850, 0.05, 0.1, 0.15, 0.2, ])
set_param("24D", "measurement", "num_charm_wall_src")(2)
set_param("24D", f"cg_params-2-2", "maxiter")(200)
set_param("24D", f"cg_params-2-2", "maxcycle")(50)

set_param("32Dfine", "trajs")(list(range(1000, 2600, 10)))
set_param("32Dfine", "measurement", "charm_quark_mass_list")([ 0.045, 0.05, 0.1, 0.15, 0.2, ])
set_param("32Dfine", "measurement", "num_charm_wall_src")(2)
set_param("32Dfine", f"cg_params-2-2", "maxiter")(200)
set_param("32Dfine", f"cg_params-2-2", "maxcycle")(50)

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
set_param(job_tag, "fermion_params", 0, 0)({ 'Ls': 8, 'M5': 1.8, 'b': 1.5, 'c': 0.5, 'boundary_phases': [1.0, 1.0, 1.0, 1.0], })
for inv_type in [ 1, 2, ]:
    set_param(job_tag, "fermion_params", inv_type, 0)(get_param(job_tag, "fermion_params", 0, 0).copy())
set_param(job_tag, "fermion_params", 0, 0, "mass")(0.01)
set_param(job_tag, "fermion_params", 1, 0, "mass")(0.04)
set_param(job_tag, "fermion_params", 2, 0, "mass")(0.10)
for inv_type in [ 0, 1, 2, ]:
    for inv_acc in [ 1, 2, ]:
        set_param(job_tag, "fermion_params", inv_type, inv_acc)(get_param(job_tag, "fermion_params", inv_type, 0).copy())
#
set_param(job_tag, "lanc_params", 0, 0, "cheby_params")({ "low": 0.5, "high": 5.5, "order": 40, })
set_param(job_tag, "lanc_params", 0, 0, "irl_params")({ "Nstop": 100, "Nk": 150, "Nm": 200, "resid": 1e-8, "betastp": 0.0, "maxiter": 20, "Nminres": 0, })
set_param(job_tag, "lanc_params", 0, 0, "pit_params")({ 'eps': 0.01, 'maxiter': 500, 'real': True, })
set_param(job_tag, "lanc_params", 1, 0)(get_param(job_tag, "lanc_params", 0, 0).copy())
#
for inv_type in [ 0, 1, ]:
    set_param(job_tag, "lanc_params", inv_type, 0, "fermion_params")(get_param(job_tag, "fermion_params", inv_type, 0).copy())
#
set_param(job_tag, "clanc_params", 0, 0, "nbasis")(100)
set_param(job_tag, "clanc_params", 0, 0, "block")([ 4, 4, 2, 2, ])
set_param(job_tag, "clanc_params", 0, 0, "cheby_params")({ "low": 0.5, "high": 5.5, "order": 40, })
set_param(job_tag, "clanc_params", 0, 0, "save_params")({ "nsingle": 100, "mpi": [ 1, 1, 1, 4, ], })
set_param(job_tag, "clanc_params", 0, 0, "irl_params")({ "Nstop": 100, "Nk": 150, "Nm": 200, "resid": 1e-8, "betastp": 0.0, "maxiter": 20, "Nminres": 0, })
set_param(job_tag, "clanc_params", 0, 0, "smoother_params")({'eps': 1e-08, 'maxiter': 20})
set_param(job_tag, "clanc_params", 1, 0)(get_param(job_tag, "clanc_params", 0, 0).copy())
#
for inv_type in [ 0, 1, 2, ]:
    for inv_acc in [ 0, 1, 2, ]:
        set_param(job_tag, f"cg_params-{inv_type}-{inv_acc}", "maxiter")(10)
        set_param(job_tag, f"cg_params-{inv_type}-{inv_acc}", "maxcycle")(inv_acc)
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
set_param(job_tag, "measurement", "auto_contractor_chunk_size")(2)
set_param(job_tag, "measurement", "charm_quark_mass_list")([ 0.1, 0.2, ])
set_param(job_tag, "measurement", "num_charm_wall_src")(2)

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
