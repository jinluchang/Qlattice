#!/usr/bin/env python3

json_results = []
check_eps = 5e-5

def json_results_append(*args):
    q.displayln_info(r"//------------------------------------------------------------\\")
    q.displayln_info(-1, *args)
    q.displayln_info(r"\\------------------------------------------------------------//")
    json_results.append(args)

from auto_contractor.operators import *

import functools
import math
import os
import time
import importlib
import sys

import qlat_gpt as qg

from qlat_scripts.v1 import *

is_cython = False

# ----

load_path_list[:] = [
        "results",
        "qcddata",
        "/lustre20/volatile/qcdqedta/qcddata",
        "/lustre20/volatile/decay0n2b/qcddata",
        "/lustre20/volatile/pqpdf/ljin/qcddata",
        "/data1/qcddata3",
        "/data2/qcddata3-prop",
        ]

# ----

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
            return get(pos_snk)
        key = (flavor, pos_snk, type_snk, type_src,)
        get = prop_dict.get(key)
        if get is not None:
            return ("g5_herm", get(pos_src),)
        fname = q.get_fname()
        raise Exception(f"{fname}: {flavor} {p_snk} {p_src}")
    return get_prop

def run_get_prop_wsrc_charm(job_tag, traj, *, get_gf, get_gt, charm_mass, tslice_list):
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
    inv = ru.get_inv(gf, job_tag, inv_type, inv_acc, gt=gt, eig=eig)
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
            assert isinstance(pos_snk, q.Coordinate)
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
    return mk_get_prop(prop_dict)

@q.timer_verbose
def auto_contract_eta_c_corr(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract/traj-{traj}/eta_c.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_eta_c_corr()
    expr_names = get_expr_names(cexpr)
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    t_size = total_site[3]
    get_prop = get_get_prop()
    psel_prob = get_psel_prob()
    fsel_prob = get_fsel_prob()
    psel = psel_prob.psel
    fsel = fsel_prob.fsel
    if not fsel.is_containing(psel):
        q.displayln_info(-1, f"WARNING: fsel is not containing psel. The probability weighting may be wrong.")
    fsel_n_elems = fsel.n_elems
    fsel_prob_arr = fsel_prob[:].ravel()
    psel_prob_arr = psel_prob[:].ravel()
    xg_fsel_arr = fsel.to_psel_local()[:]
    geo = q.Geometry(total_site)
    total_volume = geo.total_volume
    def load_data():
        t_t_list = q.get_mpi_chunk(
                [ (t_src, t_snk,) for t_snk in range(total_site[3]) for t_src in range(total_site[3]) ],
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
    q.timer_fork(0)
    res_sum = q.parallel_map_sum(feval, load_data(), sum_function=sum_function, chunksize=auto_contractor_chunk_size)
    q.displayln_info(f"{fname}: timer_display for parallel_map_sum")
    q.timer_display()
    q.timer_merge()
    res_sum *= 1.0 / total_site[3]
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
def run_job(job_tag, traj):
    #
    traj_gf = traj
    if job_tag[:5] == "test-":
        # ADJUST ME
        traj_gf = 1000
        #
    #
    fns_produce = [
            f"{job_tag}/auto-contract/traj-{traj}/checkpoint.txt",
            #
            (f"{job_tag}/psel-prop-wsrc-charm/traj-{traj}.qar", f"{job_tag}/psel-prop-wsrc-light/traj-{traj}/checkpoint.txt",),
            ]
    fns_need = [
            # f"{job_tag}/gauge-transform/traj-{traj}.field",
            # f"{job_tag}/point-selection/traj-{traj}.txt",
            # f"{job_tag}/field-selection/traj-{traj}.field",
            # f"{job_tag}/wall-src-info-light/traj-{traj}.txt",
            # f"{job_tag}/wall-src-info-strange/traj-{traj}.txt",
            # (f"{job_tag}/configs/ckpoint_lat.{traj}", f"{job_tag}/configs/ckpoint_lat.IEEE64BIG.{traj}",),
            ]
    if not check_job(job_tag, traj, fns_produce, fns_need):
        return
    #
    get_gf = run_gf(job_tag, traj_gf)
    get_gt = run_gt(job_tag, traj_gf, get_gf)
    #
    fn_checkpoint = f"{job_tag}/auto-contract/traj-{traj}/checkpoint.txt"
    if get_load_path(fn_checkpoint) is None:
        if q.obtain_lock(f"locks/{job_tag}-{traj}-auto-contract"):
            get_prop = get_get_prop()
            if get_prop is not None:
                q.timer_fork()
                # ADJUST ME
                auto_contract_eta_c_corr(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob)
                #
                q.qtouch_info(get_save_path(fn_checkpoint))
                q.displayln_info("timer_display for runjob")
                q.timer_display()
                q.timer_merge()
            q.release_lock()
    #
    q.clean_cache()
    if q.obtained_lock_history_list:
        q.timer_display()

@q.timer_verbose
def run_job_contraction(job_tag, traj):
    #
    traj_gf = traj
    if job_tag[:5] == "test-":
        # ADJUST ME
        traj_gf = 1000
        #
    #
    fns_produce = [
            f"{job_tag}/auto-contract/traj-{traj}/checkpoint.txt",
            #
            ]
    fns_need = [
            (f"{job_tag}/prop-psrc-light/traj-{traj}.qar", f"{job_tag}/prop-psrc-light/traj-{traj}/geon-info.txt",),
            (f"{job_tag}/psel-prop-psrc-light/traj-{traj}.qar", f"{job_tag}/psel-prop-psrc-light/traj-{traj}/checkpoint.txt",),
            (f"{job_tag}/prop-psrc-strange/traj-{traj}.qar", f"{job_tag}/prop-psrc-strange/traj-{traj}/geon-info.txt",),
            (f"{job_tag}/psel-prop-psrc-strange/traj-{traj}.qar", f"{job_tag}/psel-prop-psrc-strange/traj-{traj}/checkpoint.txt",),
            #
            (f"{job_tag}/prop-wsrc-light/traj-{traj}.qar", f"{job_tag}/prop-wsrc-light/traj-{traj}/geon-info.txt",),
            (f"{job_tag}/psel-prop-wsrc-light/traj-{traj}.qar", f"{job_tag}/psel-prop-wsrc-light/traj-{traj}/checkpoint.txt",),
            (f"{job_tag}/prop-wsrc-strange/traj-{traj}.qar", f"{job_tag}/prop-wsrc-strange/traj-{traj}/geon-info.txt",),
            (f"{job_tag}/psel-prop-wsrc-strange/traj-{traj}.qar", f"{job_tag}/psel-prop-wsrc-strange/traj-{traj}/checkpoint.txt",),
            f"{job_tag}/gauge-transform/traj-{traj_gf}.field",
            f"{job_tag}/point-selection/traj-{traj}.txt",
            f"{job_tag}/field-selection/traj-{traj}.field",
            # f"{job_tag}/wall-src-info-light/traj-{traj}.txt",
            # f"{job_tag}/wall-src-info-strange/traj-{traj}.txt",
            # (f"{job_tag}/configs/ckpoint_lat.{traj}", f"{job_tag}/configs/ckpoint_lat.IEEE64BIG.{traj}",),
            ]
    if not check_job(job_tag, traj, fns_produce, fns_need):
        return
    #
    get_gf = None
    get_gt = run_gt(job_tag, traj_gf, get_gf)
    #
    get_fsel_prob = run_fsel_prob(job_tag, traj, get_f_rand_01=None, get_f_weight=None)
    get_psel_prob = run_psel_prob(job_tag, traj, get_f_rand_01=None, get_f_weight=None)
    get_fsel = run_fsel_from_fsel_prob(get_fsel_prob)
    get_psel = run_psel_from_psel_prob(get_psel_prob)
    #
    get_get_prop = run_get_prop(job_tag, traj,
            get_gf = get_gf,
            get_gt = get_gt,
            get_psel = get_psel,
            get_fsel = get_fsel,
            prop_types = [
                "wsrc psel s",
                "wsrc psel l",
                "wsrc fsel s",
                "wsrc fsel l",
                "psrc psel s",
                "psrc psel l",
                "psrc fsel s",
                "psrc fsel l",
                # "rand_u1 fsel c",
                # "rand_u1 fsel s",
                # "rand_u1 fsel l",
                ],
            )
    #
    run_r_list(job_tag)
    #
    fn_checkpoint = f"{job_tag}/auto-contract/traj-{traj}/checkpoint.txt"
    if get_load_path(fn_checkpoint) is None:
        if q.obtain_lock(f"locks/{job_tag}-{traj}-auto-contract"):
            get_prop = get_get_prop()
            if get_prop is not None:
                q.timer_fork()
                # ADJUST ME
                auto_contract_eta_c_corr(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob)
                #
                q.qtouch_info(get_save_path(fn_checkpoint))
                q.displayln_info("timer_display for runjob")
                q.timer_display()
                q.timer_merge()
            q.release_lock()
    q.clean_cache()
    if q.obtained_lock_history_list:
        q.timer_display()

### ------

def get_all_cexpr():
    benchmark_eval_cexpr(get_cexpr_eta_c_corr())

### ------

set_param("test-4nt8", "trajs")(list(range(1000, 1001)))
set_param("test-4nt8", "measurement", "meson_tensor_t_sep")(1)
set_param("test-4nt8", "measurement", "auto_contractor_chunk_size")(2)

set_param("test-4nt8", "mk_sample_gauge_field", "rand_n_step")(2)
set_param("test-4nt8", "mk_sample_gauge_field", "flow_n_step")(8)
set_param("test-4nt8", "mk_sample_gauge_field", "hmc_n_traj")(1)
set_param("test-4nt8", "lanc_params", 0, 0, "cheby_params")({ "low": 0.5, "high": 5.5, "order": 40, })
set_param("test-4nt8", "lanc_params", 0, 0, "irl_params")({ "Nstop": 100, "Nk": 150, "Nm": 200, "resid": 1e-8, "betastp": 0.0, "maxiter": 20, "Nminres": 0, })
set_param("test-4nt8", "clanc_params", 0, 0, "nbasis")(100)
set_param("test-4nt8", "clanc_params", 0, 0, "block")([ 4, 4, 2, 2, ])
set_param("test-4nt8", "clanc_params", 0, 0, "cheby_params")({ "low": 0.5, "high": 5.5, "order": 40, })
set_param("test-4nt8", "clanc_params", 0, 0, "save_params")({ "nsingle": 100, "mpi": [ 1, 1, 1, 4, ], })
set_param("test-4nt8", "clanc_params", 0, 0, "irl_params")({ "Nstop": 100, "Nk": 150, "Nm": 200, "resid": 1e-8, "betastp": 0.0, "maxiter": 20, "Nminres": 0, })
set_param("test-4nt8", "clanc_params", 1, 0)(get_param("test-4nt8", "clanc_params", 0, 0).copy())
set_param("test-4nt8", "lanc_params", 1, 0)(get_param("test-4nt8", "lanc_params", 0, 0).copy())
set_param("test-4nt8", "lanc_params", 1, 0, "fermion_params")(get_param("test-4nt8", "fermion_params", 1, 0).copy())
set_param("test-4nt8", "cg_params-0-0", "maxiter")(5)
set_param("test-4nt8", "cg_params-0-1", "maxiter")(5)
set_param("test-4nt8", "cg_params-0-2", "maxiter")(5)
set_param("test-4nt8", "cg_params-1-0", "maxiter")(5)
set_param("test-4nt8", "cg_params-1-1", "maxiter")(5)
set_param("test-4nt8", "cg_params-1-2", "maxiter")(5)
set_param("test-4nt8", "cg_params-0-0", "maxcycle")(1)
set_param("test-4nt8", "cg_params-0-1", "maxcycle")(2)
set_param("test-4nt8", "cg_params-0-2", "maxcycle")(3)
set_param("test-4nt8", "cg_params-1-0", "maxcycle")(1)
set_param("test-4nt8", "cg_params-1-1", "maxcycle")(2)
set_param("test-4nt8", "cg_params-1-2", "maxcycle")(3)
set_param("test-4nt8", "fermion_params", 0, 2, "Ls")(8)
set_param("test-4nt8", "fermion_params", 1, 2, "Ls")(8)
set_param("test-4nt8", "fermion_params", 2, 2, "Ls")(8)

set_param("16IH2", "trajs")(list(range(1000, 4020, 10)))
set_param("16IH2", "measurement", "auto_contractor_chunk_size")(128)
set_param("16IH2", "measurement", "meson_tensor_t_sep")(2)
set_param("16IH2", "measurement", "pipi_op_t_sep")(2)
set_param("16IH2", "measurement", "pipi_corr_t_sep_list")(list(range(1, 10)))
set_param("16IH2", "measurement", "pipi_tensor_t_sep_list")([ 1, 2, ])
set_param("16IH2", "measurement", "pipi_tensor_t_max")(6)
set_param("16IH2", "measurement", "pipi_tensor_r_max")(16)

set_param("24D", "trajs")([ 2430, 2550, 2590, 2610, 2630, 2940, 2960, ])
set_param("24D", "measurement", "meson_tensor_t_sep")(8)
set_param("24D", "measurement", "auto_contractor_chunk_size")(128)

set_param("48I", "trajs")(list(range(1000, 2000, 20)))
set_param("48I", "measurement", "meson_tensor_t_sep")(12)
set_param("48I", "measurement", "auto_contractor_chunk_size")(128)

set_param("64I", "trajs")(list(range(1200, 3000, 40)))
set_param("64I", "measurement", "meson_tensor_t_sep")(18)
set_param("64I", "measurement", "auto_contractor_chunk_size")(128)

# ----

job_tag = "test-4nt8-checker"

set_param(job_tag, "trajs")([ 1000, ])

set_param(job_tag, "total_site")([ 4, 4, 4, 8, ])
set_param(job_tag, "load_config_params", "twist_boundary_at_boundary")([ 0.0, 0.0, 0.0, -0.5, ])

set_param(job_tag, "mk_sample_gauge_field", "rand_n_step")(2)
set_param(job_tag, "mk_sample_gauge_field", "flow_n_step")(8)
set_param(job_tag, "mk_sample_gauge_field", "hmc_n_traj")(5)

set_param(job_tag, "fermion_params", 0, 0)({ 'Ls': 8, 'M5': 1.8, 'b': 1.5, 'boundary_phases': [1.0, 1.0, 1.0, 1.0], 'c': 0.5, })
for inv_type in [ 1, 2, ]:
    set_param(job_tag, "fermion_params", inv_type, 0)(get_param(job_tag, "fermion_params", 0, 0).copy())
set_param(job_tag, "fermion_params", 0, 0, "mass")(0.01)
set_param(job_tag, "fermion_params", 1, 0, "mass")(0.04)
set_param(job_tag, "fermion_params", 2, 0, "mass")(0.04)
for inv_type in [ 0, 1, 2, ]:
    for inv_acc in [ 1, 2, ]:
        set_param(job_tag, "fermion_params", inv_type, inv_acc)(get_param(job_tag, "fermion_params", inv_type, 0).copy())

set_param(job_tag, "lanc_params", 0, 0, "cheby_params")({ "low": 0.5, "high": 5.5, "order": 40, })
set_param(job_tag, "lanc_params", 0, 0, "irl_params")({ "Nstop": 100, "Nk": 150, "Nm": 200, "resid": 1e-8, "betastp": 0.0, "maxiter": 20, "Nminres": 0, })
set_param(job_tag, "lanc_params", 0, 0, "pit_params")({ 'eps': 0.01, 'maxiter': 500, 'real': True, })
set_param(job_tag, "lanc_params", 1, 0)(get_param(job_tag, "lanc_params", 0, 0).copy())

for inv_type in [ 0, 1, ]:
    set_param(job_tag, "lanc_params", inv_type, 0, "fermion_params")(get_param(job_tag, "fermion_params", inv_type, 0).copy())

set_param(job_tag, "clanc_params", 0, 0, "nbasis")(100)
set_param(job_tag, "clanc_params", 0, 0, "block")([ 4, 4, 2, 2, ])
set_param(job_tag, "clanc_params", 0, 0, "cheby_params")({ "low": 0.5, "high": 5.5, "order": 40, })
set_param(job_tag, "clanc_params", 0, 0, "save_params")({ "nsingle": 100, "mpi": [ 1, 1, 1, 4, ], })
set_param(job_tag, "clanc_params", 0, 0, "irl_params")({ "Nstop": 100, "Nk": 150, "Nm": 200, "resid": 1e-8, "betastp": 0.0, "maxiter": 20, "Nminres": 0, })
set_param(job_tag, "clanc_params", 0, 0, "smoother_params")({'eps': 1e-08, 'maxiter': 20})
set_param(job_tag, "clanc_params", 1, 0)(get_param(job_tag, "clanc_params", 0, 0).copy())

for inv_type in [ 0, 1, 2, ]:
    for inv_acc in [ 0, 1, 2, ]:
        set_param(job_tag, f"cg_params-{inv_type}-{inv_acc}", "maxiter")(10)
        set_param(job_tag, f"cg_params-{inv_type}-{inv_acc}", "maxcycle")(inv_acc)

set_param(job_tag, "field-selection-fsel-rate")(0.1)
set_param(job_tag, "field-selection-psel-rate")(0.01)
set_param(job_tag, "field-selection-fsel-psrc-prop-norm-threshold")(0.05)

set_param(job_tag, "prob_exact_wsrc")(0.20)

set_param(job_tag, "prob_acc_1_psrc")(0.25)
set_param(job_tag, "prob_acc_2_psrc")(0.10)

set_param(job_tag, "measurement", "auto_contractor_chunk_size")(2)
set_param(job_tag, "measurement", "meson_tensor_t_sep")(1)
set_param(job_tag, "measurement", "pipi_op_t_sep")(1)
set_param(job_tag, "measurement", "pipi_corr_t_sep_list")([ 1, 2, 3, 4, ])
set_param(job_tag, "measurement", "pipi_tensor_t_sep_list")([ 1, 2, 3, ])
set_param(job_tag, "measurement", "pipi_tensor_t_max")(3)
set_param(job_tag, "measurement", "pipi_tensor_r_max")(4)

# ----

def gracefully_finish():
    q.timer_display()
    qg.end_with_gpt()
    q.displayln_info("CHECK: finished successfully.")
    exit()

if __name__ == "__main__":

    qg.begin_with_gpt()

    ##################### CMD options #####################

    job_tags = q.get_arg("--job_tags", default="").split(",")

    is_performing_inversion = not q.get_option("--no-inversion")

    is_performing_contraction = not q.get_option("--no-contraction")

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
            if is_performing_inversion:
                q.check_time_limit()
                run_job_inversion(job_tag, traj)
                if q.obtained_lock_history_list:
                    json_results_append(f"q.obtained_lock_history_list={q.obtained_lock_history_list}")
                    if job_tag[:5] != "test-":
                        gracefully_finish()
        for traj in get_param(job_tag, "trajs"):
            if is_performing_contraction:
                q.check_time_limit()
                run_job_contraction(job_tag, traj)
                if q.obtained_lock_history_list:
                    json_results_append(f"q.obtained_lock_history_list={q.obtained_lock_history_list}")
                    if job_tag[:5] != "test-":
                        gracefully_finish()

    q.check_log_json(__file__, json_results, check_eps=check_eps)

    gracefully_finish()

# ----
