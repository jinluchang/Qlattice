#!/usr/bin/env python3

import functools
import math
import os
import time
import importlib
import sys

import qlat_gpt as qg
import qlat as q
import qlat_scripts.v1 as qs
import numpy as np

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
        get_job_seed,
        run_params,
        check_job,
        run_gf,
        run_gt,
        get_load_path,
        get_save_path,
        is_test,
        )

### ------

load_path_list[:] = [
        "results",
        "qcddata",
        "/lustre/orion/lgt119/proj-shared/ljin/qcddata4",
        "/lustre/orion/lgt119/proj-shared/ljin/qcddata5",
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

is_cython = not is_test()

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
def auto_contract_eta_c_corr(job_tag, traj, get_get_prop, inv_type, tslice_list):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract-eta-c/traj-{traj}/eta_c/inv_type-{inv_type}.lat"
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
    q.json_results_append(f"{fname}: ld sig", q.get_data_sig(ld, q.RngState()))
    for i, en in enumerate(expr_names):
        q.json_results_append(f"{fname}: ld '{en}' sig", q.get_data_sig(ld[i], q.RngState()))

@q.timer_verbose
def auto_contract_eta_c_corr_psnk(job_tag, traj, get_get_prop, inv_type, tslice_list):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract-eta-c/traj-{traj}/eta_c-psnk/inv_type-{inv_type}.lat"
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
        t_t_list = [ (t_src, t_snk,) for t_snk in range(total_site[3]) for t_src in tslice_list ]
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
    with q.TimerFork(max_call_times_for_always_show_info=0):
        res_sum = q.parallel_map_sum(feval, load_data(), sum_function=sum_function, chunksize=1)
        q.displayln_info(f"{fname}: timer_display for parallel_map_sum")
    res_sum *= 1.0 / len(tslice_list) / spatial_volume
    assert q.qnorm(res_sum[0] - 1.0) < 1e-10
    ld = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t_sep", t_size, [ str(q.rel_mod(t, t_size)) for t in range(t_size) ], ],
        ])
    ld.from_numpy(res_sum)
    ld.save(get_save_path(fn))
    q.json_results_append(f"{fname}: ld sig", q.get_data_sig(ld, q.RngState()))
    for i, en in enumerate(expr_names):
        q.json_results_append(f"{fname}: ld '{en}' sig", q.get_data_sig(ld[i], q.RngState()))

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

def mk_get_prop_point_snk(prop_ps, geo):
    def get(p_snk):
        type_snk, pos_snk = p_snk
        assert type_snk == "point-snk"
        assert isinstance(pos_snk, tuple)
        pos_snk = q.Coordinate(pos_snk)
        index = geo.index_from_g_coordinate(pos_snk)
        return prop_ps.get_elem_wm(index)
    return get

def mk_get_prop_wall_snk(ps_prop_ws):
    def get(p_snk):
        type_snk, pos_snk = p_snk
        assert type_snk == "wall"
        assert isinstance(pos_snk, int)
        return ps_prop_ws.get_elem_wm(pos_snk)
    return get

@q.timer_verbose
def run_get_prop_wsrc_charm(job_tag, traj, *, get_gf, get_gt, inv_type, tslice_list):
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
        inv_acc = 2
        eig = None
        inv = qs.get_inv(gf, job_tag, inv_type, inv_acc, gt=gt, eig=eig)
        prop_dict = dict()
        for tslice in tslice_list:
            src = q.mk_wall_src(geo, tslice)
            sol = inv * src
            prop = sol
            ps_prop_ws = prop.glb_sum_tslice()
            prop_ps = gt_inv * prop
            key = ("c", tslice, "wall", "point-snk")
            prop_dict[key] = mk_get_prop_point_snk(prop_ps, geo)
            key = ("c", tslice, "wall", "wall")
            prop_dict[key] = mk_get_prop_wall_snk(ps_prop_ws)
        q.clean_cache(q.cache_inv)
        prop_cache = q.mk_cache("prop_cache", f"{job_tag}", f"{traj}")
        prop_cache["prop-dict"] = prop_dict
        return mk_get_prop(prop_dict)
    return q.lazy_call(get_get_prop)

### ------

@q.timer(is_timer_fork=True)
def run_eta_c_corr(job_tag, traj, get_gf, get_gt):
    fname = q.get_fname()
    fn_checkpoint = f"{job_tag}/auto-contract-eta-c/traj-{traj}/checkpoint.txt"
    if get_load_path(fn_checkpoint) is not None:
        return
    if not q.obtain_lock(f"locks/{job_tag}-{traj}-{fname}"):
        return
    charm_wall_src_prop_params = run_charm_wall_src_prop_params(job_tag, traj)
    tslice_list = get_param_charm_wall_src_tslice_list(job_tag, traj)
    for inv_type in charm_wall_src_prop_params["charm_quark_inv_type_list"]:
        get_get_prop = run_get_prop_wsrc_charm(
                job_tag, traj,
                get_gf=get_gf,
                get_gt=get_gt,
                inv_type=inv_type,
                tslice_list=tslice_list,
                )
        auto_contract_eta_c_corr(job_tag, traj, get_get_prop, inv_type, tslice_list)
        auto_contract_eta_c_corr_psnk(job_tag, traj, get_get_prop, inv_type, tslice_list)
    q.qtouch_info(get_save_path(fn_checkpoint))
    q.release_lock()

### ------

def get_param_charm_mass_list(job_tag):
    charm_quark_mass_list = get_param(job_tag, "quark_mass_list")[2:]
    return charm_quark_mass_list

@q.cache_call()
@q.timer_verbose
def get_param_charm_wall_src_tslice_list(job_tag, traj):
    num_charm_wall_src = get_param(job_tag, "measurement", "num_charm_wall_src")
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    t_size = total_site[3]
    rs = q.RngState(f"{get_job_seed(job_tag)}-{traj}-get_param_charm_wall_src_tslice_list")
    t_start = rs.u_rand_gen() * t_size
    t_sep = t_size / num_charm_wall_src
    tslice_list = [ round(t_start + i * t_sep) % t_size for i in range(num_charm_wall_src) ]
    return tslice_list

@q.timer_verbose
def run_charm_wall_src_prop_params(job_tag, traj):
    assert len(get_param(job_tag, "quark_flavor_list")) == len(get_param(job_tag, "quark_mass_list"))
    inv_type_ref = 1
    n_flavor = len(get_param(job_tag, "quark_flavor_list"))
    inv_type_list = list(range(n_flavor))
    charm_wall_src_tslice_list = get_param_charm_wall_src_tslice_list(job_tag, traj)
    obj = dict()
    obj["charm_quark_inv_type_list"] = inv_type_list[inv_type_ref + 1:]
    obj["charm_quark_flavor_list"] = get_param(job_tag, "quark_flavor_list")[inv_type_ref + 1:]
    obj["charm_quark_mass_list"] = get_param(job_tag, "quark_mass_list")[inv_type_ref + 1:]
    obj["charm_wall_src_tslice_list"] = charm_wall_src_tslice_list
    fn = f"{job_tag}/params-eta-c/traj-{traj}/charm_wall_src_prop.json"
    path = get_load_path(fn)
    if path is None:
        path = get_save_path(fn)
        q.save_json_obj(obj, path, indent=2, is_sync_node=True)
    else:
        obj_load = q.load_json_obj(path, is_sync_node=True)
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
            f"{job_tag}/auto-contract-eta-c/traj-{traj}/checkpoint.txt",
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
    run_eta_c_corr(job_tag, traj, get_gf, get_gt)

### ------

def get_all_cexpr():
    benchmark_eval_cexpr(get_cexpr_eta_c_corr())

### ------

job_tag = "16IH2"
set_param(job_tag, "traj_list")(list(range(1000, 4020, 100)))
set_param(job_tag, "measurement", "auto_contractor_chunk_size")(128)
set_param(job_tag, "measurement", "num_charm_wall_src")(2)
set_param(job_tag, "quark_flavor_list")([ "light", "strange", ] + [ f"charm-{idx+1}" for idx in range(17) ])
set_param(job_tag, "quark_mass_list")([ 0.01, 0.04, 0.02963, 0.05358, 0.07945, 0.10852, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.6, 0.7, 0.8, 0.9, 1.0, ])
for inv_type, mass in list(enumerate(get_param(job_tag, "quark_mass_list")))[2:]:
    set_param(job_tag, "fermion_params", inv_type, 0)(get_param(job_tag, "fermion_params", 0, 0).copy())
    set_param(job_tag, "fermion_params", inv_type, 0, "mass")(mass)
    for inv_acc in [ 0, 1, 2, ]:
        set_param(job_tag, "fermion_params", inv_type, inv_acc)(get_param(job_tag, "fermion_params", inv_type, 0).copy())
        set_param(job_tag, f"cg_params-{inv_type}-{inv_acc}", "maxiter")(200)
        set_param(job_tag, f"cg_params-{inv_type}-{inv_acc}", "maxcycle")(50)

job_tag = "24D"
set_param(job_tag, "traj_list")(list(range(1000, 5100, 80)))
set_param(job_tag, "measurement", "auto_contractor_chunk_size")(128)
set_param(job_tag, "measurement", "num_charm_wall_src")(2)
set_param(job_tag, "quark_flavor_list")([ "light", "strange", ] + [ f"charm-{idx+1}" for idx in range(14) ])
set_param(job_tag, "quark_mass_list")([ 0.00107, 0.0850, 0.07819, 0.13207, 0.19829, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.6, 0.7, 0.8, 0.9, 1.0, ])
for inv_type, mass in list(enumerate(get_param(job_tag, "quark_mass_list")))[2:]:
    set_param(job_tag, "fermion_params", inv_type, 0)(get_param(job_tag, "fermion_params", 0, 2).copy())
    set_param(job_tag, "fermion_params", inv_type, 0, "mass")(mass)
    for inv_acc in [ 0, 1, 2, ]:
        set_param(job_tag, "fermion_params", inv_type, inv_acc)(get_param(job_tag, "fermion_params", inv_type, 0).copy())
        set_param(job_tag, f"cg_params-{inv_type}-{inv_acc}", "maxiter")(200)
        set_param(job_tag, f"cg_params-{inv_type}-{inv_acc}", "maxcycle")(50)

job_tag = "32Dfine"
set_param(job_tag, "traj_list")(list(range(520, 2600, 40)))
set_param(job_tag, "measurement", "auto_contractor_chunk_size")(128)
set_param(job_tag, "measurement", "num_charm_wall_src")(2)
set_param(job_tag, "quark_flavor_list")([ "light", "strange", ] + [ f"charm-{idx+1}" for idx in range(16) ])
set_param(job_tag, "quark_mass_list")([ 0.0001, 0.045, 0.04635, 0.07794, 0.11333, 0.15327, 0.2, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.6, 0.7, 0.8, 0.9, 1.0, ])
for inv_type, mass in list(enumerate(get_param(job_tag, "quark_mass_list")))[2:]:
    set_param(job_tag, "fermion_params", inv_type, 0)(get_param(job_tag, "fermion_params", 0, 2).copy())
    set_param(job_tag, "fermion_params", inv_type, 0, "mass")(mass)
    for inv_acc in [ 0, 1, 2, ]:
        set_param(job_tag, "fermion_params", inv_type, inv_acc)(get_param(job_tag, "fermion_params", inv_type, 0).copy())
        set_param(job_tag, f"cg_params-{inv_type}-{inv_acc}", "maxiter")(200)
        set_param(job_tag, f"cg_params-{inv_type}-{inv_acc}", "maxcycle")(50)

job_tag = "64I"
set_param(job_tag, "traj_list")(list(range(1200, 3680, 80)))
set_param(job_tag, "measurement", "auto_contractor_chunk_size")(128)
set_param(job_tag, "measurement", "num_charm_wall_src")(2)
set_param(job_tag, "quark_flavor_list")([ "light", "strange", ] + [ f"charm-{idx+1}" for idx in range(5) ])
set_param(job_tag, "quark_mass_list")([ 0.000678, 0.02661, 0.0611417 , 0.08234643, 0.17112621, 0.29854376, 0.33262794, ])
for inv_type, mass in list(enumerate(get_param(job_tag, "quark_mass_list")))[2:]:
    set_param(job_tag, "fermion_params", inv_type, 0)(get_param(job_tag, "fermion_params", 0, 2).copy())
    set_param(job_tag, "fermion_params", inv_type, 0, "mass")(mass)
    for inv_acc in [ 0, 1, 2, ]:
        set_param(job_tag, "fermion_params", inv_type, inv_acc)(get_param(job_tag, "fermion_params", inv_type, 0).copy())
        set_param(job_tag, f"cg_params-{inv_type}-{inv_acc}", "maxiter")(200)
        set_param(job_tag, f"cg_params-{inv_type}-{inv_acc}", "maxcycle")(50)

# ----

job_tag = "test-4nt8-checker"
#
set_param(job_tag, "traj_list")([ 1000, ])
#
set_param(job_tag, "seed")("test-4nt8-checker")
set_param(job_tag, "total_site")([ 4, 4, 4, 8, ])
set_param(job_tag, "load_config_params", "twist_boundary_at_boundary")([ 0.0, 0.0, 0.0, -0.5, ])
#
set_param(job_tag, "mk_sample_gauge_field", "rand_n_step")(2)
set_param(job_tag, "mk_sample_gauge_field", "flow_n_step")(8)
set_param(job_tag, "mk_sample_gauge_field", "hmc_n_traj")(5)
#
set_param(job_tag, "quark_flavor_list")([ "light", "strange", "charm-1", "charm-2", ])
set_param(job_tag, "quark_mass_list")([ 0.01, 0.04, 0.1, 0.2, ])
set_param(job_tag, "fermion_params", 0, 0)({ 'Ls': 8, 'M5': 1.8, 'b': 1.5, 'c': 0.5, 'boundary_phases': [1.0, 1.0, 1.0, 1.0], })
for inv_type, mass in enumerate(get_param(job_tag, "quark_mass_list")):
    set_param(job_tag, "fermion_params", inv_type, 0)(get_param(job_tag, "fermion_params", 0, 0).copy())
    set_param(job_tag, "fermion_params", inv_type, 0, "mass")(mass)
    for inv_acc in [ 0, 1, 2, ]:
        set_param(job_tag, "fermion_params", inv_type, inv_acc)(get_param(job_tag, "fermion_params", inv_type, 0).copy())
        set_param(job_tag, f"cg_params-{inv_type}-{inv_acc}", "maxiter")(10)
        set_param(job_tag, f"cg_params-{inv_type}-{inv_acc}", "maxcycle")(1 + inv_acc)
#
set_param(job_tag, "measurement", "auto_contractor_chunk_size")(2)
set_param(job_tag, "measurement", "num_charm_wall_src")(2)

# ----

##################### CMD options #####################

job_tag_list_default = [
        "test-4nt8-checker",
        ]
job_tag_list_str_default = ",".join(job_tag_list_default)
job_tag_list = q.get_arg("--job_tag_list", default=job_tag_list_str_default).split(",")

#######################################################

def gracefully_finish():
    q.displayln_info("Begin to gracefully_finish.")
    q.timer_display()
    if is_test():
        q.json_results_append(f"q.obtained_lock_history_list={q.obtained_lock_history_list}")
        q.check_log_json(__file__)
    qg.end_with_gpt()
    q.displayln_info("CHECK: finished successfully.")
    exit()

def try_gracefully_finish():
    """
    Call `gracefully_finish` if not test and if some work is done (q.obtained_lock_history_list != [])
    """
    if (not is_test()) and (len(q.obtained_lock_history_list) > 0):
        gracefully_finish()

if __name__ == "__main__":

    qg.begin_with_gpt()
    q.check_time_limit()
    get_all_cexpr()

    job_tag_traj_list = []
    for job_tag in job_tag_list:
        run_params(job_tag)
        traj_list = get_param(job_tag, "traj_list")
        for traj in traj_list:
            job_tag_traj_list.append((job_tag, traj,))
    if not is_test():
        job_tag_traj_list = q.random_permute(job_tag_traj_list, q.RngState(f"{q.get_time()}"))
        job_tag_traj_list = q.get_comm().bcast(job_tag_traj_list)
    for job_tag, traj in job_tag_traj_list:
        q.check_time_limit()
        run_job(job_tag, traj)
        q.clean_cache()
        try_gracefully_finish()

    gracefully_finish()

# ----
