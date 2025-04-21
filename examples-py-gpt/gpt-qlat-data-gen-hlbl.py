#!/usr/bin/env python3

json_results = []

from auto_contractor.operators import *

import functools
import math
import os
import time
import importlib
import sys
import pprint
from copy import deepcopy

import qlat_gpt as qg

from qlat_scripts.v1 import *

is_cython = False

# ----

load_path_list[:] = [
        "results",
        "results-props",
        "/lustre/orion/lgt119/proj-shared/ljin/qcddata4",
        "/lustre/orion/lgt119/proj-shared/ljin/qcddata5",
        "/lustre/orion/lgt119/proj-shared/ljin/hlbl-muon-line-data/hlbl-muon-line",
        ]

# ----

@q.timer
def get_cexpr_meson_corr():
    fn_base = "cache/auto_contract_cexpr/get_cexpr_meson_corr"
    def calc_cexpr():
        diagram_type_dict = dict()
        diagram_type_dict[((('x_1', 'x_2'), 1), (('x_2', 'x_1'), 1))] = 'Type1'
        diagram_type_dict[((('x_1', 'x_1'), 1), (('x_2', 'x_2'), 1))] = None
        exprs = [
                mk_fac(1) + f"1",
                mk_pi_p("x_2", True)    * mk_pi_p("x_1")             + f"pi+^dag(0) * pi+(-tsep)",
                mk_k_p("x_2", True)     * mk_k_p("x_1")              + f"K+^dag(0) * K+(-tsep)",
                mk_sw5("x_2")           * mk_pi_p("x_1")             + f"sw5(0) * pi+(-tsep)",
                mk_sw5("x_2")           * mk_k_p("x_1")              + f"sw5(0) * K+(-tsep)",
                mk_sw5("x_2")           * mk_j5pi_mu("x_1", 3, True) + f"sw5(0) * j5pi_t^dag(-tsep)",
                mk_sw5("x_2")           * mk_j5k_mu("x_1", 3, True)  + f"sw5(0) * j5k_t^dag(-tsep)",
                mk_jw_a_mu("x_2", 3)    * mk_pi_p("x_1")             + f"jw_a_t(0) * pi+(-tsep)",
                mk_jw_a_mu("x_2", 3)    * mk_k_p("x_1")              + f"jw_a_t(0) * K+(-tsep)",
                mk_jw_a_mu("x_2", 3)    * mk_j5pi_mu("x_1", 3, True) + f"jw_a_t(0) * j5pi_t^dag(-tsep)",
                mk_jw_a_mu("x_2", 3)    * mk_j5k_mu("x_1", 3, True)  + f"jw_a_t(0) * j5k_t^dag(-tsep)",
                mk_a0_p("x_2", True)    * mk_a0_p("x_1")             + f"a0+^dag(0) * a0+(-tsep)",
                mk_kappa_p("x_2", True) * mk_kappa_p("x_1")          + f"kappa+^dag(0) * kappa+(-tsep)",
                mk_j_mu("x_2", 3)       * mk_j_mu("x_1", 3)          + f"j_t(0) * j_t(-tsep)",
                sum([ mk_j_mu("x_2", mu) * mk_j_mu("x_1", mu) for mu in range(4) ])
                + f"j_mu(0) * j_mu(-tsep)",
                sum([ mk_jw_a_mu("x_2", mu) * mk_j5pi_mu("x_1", mu, True) for mu in range(4) ])
                + f"jw_a_mu(0) * j5pi_mu^dag(-tsep)",
                sum([ mk_jw_a_mu("x_2", mu) * mk_j5k_mu("x_1", mu, True) for mu in range(4) ])
                + f"jw_a_mu(0) * j5k_mu^dag(-tsep)",
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit=True, diagram_type_dict=diagram_type_dict)
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

@q.timer_verbose
def auto_contract_meson_corr(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_corr.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_corr()
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
    json_results.append((f"{fname}: {traj} ld sig", q.get_data_sig(ld, q.RngState()),))

@q.timer_verbose
def auto_contract_meson_corr_psnk(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_corr_psnk.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_corr()
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
        for t_src in range(total_site[3]):
            for idx in range(fsel_n_elems):
                yield t_src, idx
    @q.timer
    def feval(args):
        t_src, idx = args
        xg_snk = tuple(xg_fsel_arr[idx])
        prob_snk = fsel_prob_arr[idx]
        t = (xg_snk[3] - t_src) % total_site[3]
        pd = {
                "x_2" : ("point-snk", xg_snk,),
                "x_1" : ("wall", t_src,),
                "size" : total_site,
                }
        val = eval_cexpr(cexpr, positions_dict=pd, get_prop=get_prop)
        return val / prob_snk, t
    def sum_function(val_list):
        values = np.zeros((total_site[3], len(expr_names),), dtype=complex)
        for val, t in val_list:
            values[t] += val
        return values.transpose(1, 0)
    auto_contractor_chunk_size = get_param(job_tag, "measurement", "auto_contractor_chunk_size", default=128)
    q.timer_fork(0)
    res_sum = q.glb_sum(
            q.parallel_map_sum(feval, load_data(), sum_function=sum_function, chunksize=auto_contractor_chunk_size))
    q.displayln_info("timer_display for auto_contract_meson_corr_psnk")
    q.timer_display()
    q.timer_merge()
    res_sum *= 1.0 / total_volume
    q.displayln_info(0, res_sum[0])
    ld = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t_sep", t_size, [ str(q.rel_mod(t, t_size)) for t in range(t_size) ], ],
        ])
    ld.from_numpy(res_sum)
    ld.save(get_save_path(fn))
    json_results.append((f"{fname}: {traj} ld sig", q.get_data_sig(ld, q.RngState()),))

@q.timer_verbose
def auto_contract_meson_corr_psrc(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_corr_psrc.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_corr()
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
    xg_psel_arr = psel[:]
    geo = q.Geometry(total_site)
    total_volume = geo.total_volume
    def load_data():
        x_t_list = q.get_mpi_chunk(
                [ (pidx, t_snk,) for t_snk in range(total_site[3]) for pidx in range(len(xg_psel_arr)) ],
                rng_state=None)
        for pidx, t_snk in x_t_list:
            yield pidx, t_snk
    @q.timer
    def feval(args):
        pidx, t_snk = args
        xg_src = tuple(xg_psel_arr[pidx])
        prob_src = psel_prob_arr[pidx]
        t = (xg_src[3] - t_snk) % total_site[3]
        pd = {
                "x_2" : ("point", xg_src,),
                "x_1" : ("wall", t_snk,),
                "size" : total_site,
                }
        val = eval_cexpr(cexpr, positions_dict=pd, get_prop=get_prop)
        return val / prob_src, t
    def sum_function(val_list):
        values = np.zeros((total_site[3], len(expr_names),), dtype=complex)
        for val, t in val_list:
            values[t] += val
        return values.transpose(1, 0)
    auto_contractor_chunk_size = get_param(job_tag, "measurement", "auto_contractor_chunk_size", default=128)
    q.timer_fork(0)
    res_sum = q.glb_sum(
            q.parallel_map_sum(feval, load_data(), sum_function=sum_function, chunksize=auto_contractor_chunk_size))
    q.displayln_info("timer_display for auto_contract_meson_corr_psrc")
    q.timer_display()
    q.timer_merge()
    res_sum *= 1.0 / total_volume
    q.displayln_info(0, res_sum[0])
    ld = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t_sep", t_size, [ str(q.rel_mod(t, t_size)) for t in range(t_size) ], ],
        ])
    ld.from_numpy(res_sum)
    ld.save(get_save_path(fn))
    json_results.append((f"{fname}: {traj} ld sig", q.get_data_sig(ld, q.RngState()),))

@q.timer_verbose
def auto_contract_meson_corr_psnk_psrc(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_corr_psnk_psrc.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_corr()
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
    xg_psel_arr = psel[:]
    geo = q.Geometry(total_site)
    total_volume = geo.total_volume
    r_list = get_r_list(job_tag)
    r_sq_interp_idx_coef_list = get_r_sq_interp_idx_coef_list(job_tag)
    def load_data():
        for pidx in range(len(xg_psel_arr)):
            yield pidx
    @q.timer
    def feval(args):
        pidx = args
        xg_src = tuple(xg_psel_arr[pidx])
        prob_src = psel_prob_arr[pidx]
        res_list = []
        for idx in range(len(xg_fsel_arr)):
            xg_snk = tuple(xg_fsel_arr[idx])
            if xg_snk == xg_src:
                prob_snk = 1.0
            else:
                prob_snk = fsel_prob_arr[idx]
            prob = prob_src * prob_snk
            x_rel = [ q.rel_mod(xg_snk[mu] - xg_src[mu], total_site[mu]) for mu in range(4) ]
            r_sq = q.get_r_sq(x_rel)
            r_idx_low, r_idx_high, coef_low, coef_high = r_sq_interp_idx_coef_list[r_sq]
            t = (xg_snk[3] - xg_src[3]) % total_site[3]
            pd = {
                    "x_2": ("point-snk", xg_snk,),
                    "x_1": ("point", xg_src,),
                    "size": total_site,
                    }
            val = eval_cexpr(cexpr, positions_dict=pd, get_prop=get_prop)
            res_list.append((val / prob, t, r_idx_low, r_idx_high, coef_low, coef_high))
        return res_list
    def sum_function(val_list):
        values = np.zeros((total_site[3], len(r_list), len(expr_names),), dtype=complex)
        for idx, res_list in enumerate(val_list):
            for val, t, r_idx_low, r_idx_high, coef_low, coef_high in res_list:
                values[t, r_idx_low] += coef_low * val
                values[t, r_idx_high] += coef_high * val
            q.displayln_info(f"{fname}: {idx+1}/{len(xg_psel_arr)}")
        return values.transpose(2, 0, 1)
    q.timer_fork(0)
    res_sum = q.glb_sum(
            q.parallel_map_sum(feval, load_data(), sum_function=sum_function, chunksize=1))
    q.displayln_info("timer_display for auto_contract_meson_corr_psnk_psrc")
    q.timer_display()
    q.timer_merge()
    res_sum *= 1.0 / (total_volume**2 / total_site[3])
    q.displayln_info(res_sum[0].sum(1))
    ld = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t_sep", t_size, [ str(q.rel_mod(t, t_size)) for t in range(t_size) ], ],
        [ "r", len(r_list), [ f"{r:.5f}" for r in r_list ], ],
        ])
    ld.from_numpy(res_sum)
    ld.save(get_save_path(fn))
    json_results.append((f"{fname}: {traj} ld sig", q.get_data_sig(ld, q.RngState()),))

# ----

def get_all_cexpr():
    benchmark_eval_cexpr(get_cexpr_meson_corr())

# ----

@q.timer
def run_job_global_hvp_average(job_tag, *, inv_type):
    """
    get_glb_hvp_avg_light = run_job_global_hvp_average(job_tag, inv_type=0)
    get_glb_hvp_avg_strange = run_job_global_hvp_average(job_tag, inv_type=1)
    #
    get_glb_hvp_avg = run_job_global_hvp_average(job_tag, inv_type=inv_type)
    glb_hvp_avg_trajs, glb_hvp_avg = get_glb_hvp_avg()
    """
    fname = q.get_fname()
    inv_type_name_list = [ "light", "strange", ]
    inv_type_name = inv_type_name_list[inv_type]
    fn = f"{job_tag}/hvp-average/hvp_average_{inv_type_name}.field"
    fn_trajs = f"{job_tag}/hvp-average/hvp_average_{inv_type_name}.trajs.txt"
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    geo = q.Geometry(total_site)
    @q.lazy_call
    @q.timer_verbose
    def load_glb_hvp_avg():
        trajs = [ int(x) for x in q.qcat_sync_node(get_load_path(fn_trajs)).splitlines() ]
        assert len(trajs) > 0
        hvp_average = q.FieldComplexD(geo, 16)
        total_bytes = hvp_average.load_double_from_float(get_load_path(fn))
        assert total_bytes > 0
        return trajs, hvp_average
    ret = load_glb_hvp_avg
    if get_load_path(fn) is not None:
        assert get_load_path(fn_trajs) is not None
        return ret
    #
    q.check_stop()
    q.check_time_limit()
    if not q.obtain_lock(f"locks/{job_tag}-{fname}-{inv_type_name}"):
        return None
    q.timer_fork()
    #
    @q.timer_verbose
    def compute_glb_hvp_average():
        trajs = get_param(job_tag, "trajs")
        trajs_hvp_avg = []
        for traj in trajs:
            fns_produce = [
                    fn,
                    ]
            fns_need = [
                    f"{job_tag}/hvp-average/traj-{traj}/hvp_average_{inv_type_name}.field",
                    f"{job_tag}/wall-src-info-light/traj-{traj}.txt",
                    f"{job_tag}/wall-src-info-strange/traj-{traj}.txt",
                    f"{job_tag}/point-selection/traj-{traj}.txt",
                    f"{job_tag}/field-selection/traj-{traj}.field",
                    f"{job_tag}/field-selection-weight/traj-{traj}/f-rand-01.field",
                    f"{job_tag}/field-selection-weight/traj-{traj}/weight.field",
                    f"{job_tag}/field-selection-weight/traj-{traj}/fsel-prob.sfield",
                    f"{job_tag}/field-selection-weight/traj-{traj}/psel-prob.lat",
                    ]
            if check_job(job_tag, traj, fns_produce, fns_need):
                q.displayln_info(0, f"{fname}: {job_tag} {traj} {inv_type_name} OK")
                trajs_hvp_avg.append(traj)
            else:
                q.displayln_info(0, f"{fname}: {job_tag} {traj} {inv_type_name} not yet finished")
        #
        q.displayln_info(0, f"{fname}: {job_tag} {traj} {inv_type_name} num_trajs={len(trajs_hvp_avg)} start.")
        assert len(trajs_hvp_avg) > 0
        #
        hvp_average = q.FieldComplexD(geo, 16)
        hvp_average.set_zero()
        for traj in trajs_hvp_avg:
            get_wi = run_wi(job_tag, traj)
            get_f_weight = run_f_weight_from_wsrc_prop_full(job_tag, traj)
            get_f_rand_01 = run_f_rand_01(job_tag, traj)
            get_psel_prob = run_psel_prob(job_tag, traj, get_f_rand_01=get_f_rand_01, get_f_weight=get_f_weight)
            get_hvp_average = run_hvp_average(
                    job_tag, traj,
                    inv_type=inv_type,
                    get_psel_prob=get_psel_prob,
                    )
            hvp_average += get_hvp_average()
        hvp_average *= 1 / len(trajs_hvp_avg)
        q.qtouch_info(get_save_path(fn_trajs), [ f"{traj}\n" for traj in trajs_hvp_avg ])
        hvp_average.save_float_from_double(get_save_path(fn))
        q.displayln_info(-1, f"{fname}: {job_tag} {inv_type_name} num_trajs={len(trajs_hvp_avg)} finished.")
    #
    compute_glb_hvp_average()
    #
    q.release_lock()
    q.clean_cache()
    q.timer_display()
    q.timer_merge()
    return ret

@q.timer_verbose
def run_job_global_hvp_average_for_subtract(job_tag, traj, *, inv_type, get_glb_hvp_avg, get_hvp_average):
    """
    get_glb_hvp_avg_for_sub = run_job_global_hvp_average_for_subtract(job_tag, traj, inv_type=inv_type, get_glb_hvp_avg=get_glb_hvp_avg, get_hvp_average=get_hvp_average)
    glb_hvp_avg_for_sub = get_glb_hvp_avg_for_sub()
    #
    Get global hvp average excluding data from this traj.
    Suitable for use in subtraction.
    """
    fname = q.get_fname()
    inv_type_name_list = [ "light", "strange", ]
    inv_type_name = inv_type_name_list[inv_type]
    fn = f"{job_tag}/hlbl/glb-hvp-avg-for-sub/traj-{traj}/hvp_average_{inv_type_name}.field"
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    geo = q.Geometry(total_site)
    @q.lazy_call
    def get_glb_hvp_avg_for_sub():
        glb_hvp_avg_for_sub = q.FieldComplexD(geo, 16)
        glb_hvp_avg_for_sub.load_double_from_float(get_load_path(fn))
        return glb_hvp_avg_for_sub
    ret = get_glb_hvp_avg_for_sub
    if get_load_path(fn) is not None:
        return ret
    if get_glb_hvp_avg is None:
        q.displayln_info(-1, f"{fname}: get_glb_hvp_avg is None.")
        return None
    if get_hvp_average is None:
        q.displayln_info(-1, f"{fname}: get_hvp_average is None.")
        return None
    if not q.obtain_lock(f"locks/{job_tag}-{traj}-{fname}-{inv_type_name}"):
        return
    q.timer_fork()
    glb_hvp_avg_trajs, glb_hvp_avg = get_glb_hvp_avg()
    glb_hvp_avg_for_sub = glb_hvp_avg.copy()
    if traj not in glb_hvp_avg_trajs:
        q.displayln_info(-1, f"WARNING: {fname} glb_hvp_avg traj={traj} not in glb_hvp_avg_trajs.")
    else:
        num_trajs = len(glb_hvp_avg_trajs)
        assert num_trajs > 0
        if num_trajs == 1:
            q.displayln_info(-1, f"WARNING: {fname} glb_hvp_avg num_trajs={num_trajs}")
        else:
            hvp_average = get_hvp_average()
            glb_hvp_avg_for_sub *= num_trajs
            glb_hvp_avg_for_sub -= hvp_average
            glb_hvp_avg_for_sub *= 1.0 / (num_trajs - 1.0)
    glb_hvp_avg_for_sub.save_float_from_double(get_save_path(fn))
    q.release_lock()
    q.timer_display()
    q.timer_merge()
    return ret

# ----

def get_muon_mass(job_tag):
    a_inv_gev = get_param(job_tag, "a_inv_gev")
    muon_mass_mev = 105.6583745
    return muon_mass_mev / 1e3 / a_inv_gev

def get_r_sq_limit(job_tag):
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    return q.sqr(total_site[0])

def get_r_coordinate(xg, total_site):
    xg_rel = q.smod_coordinate(xg, total_site)
    r_sq = xg_rel.sqr()
    return np.sqrt(r_sq)

def show_lslt(labels, lslt, *, label=None):
    nlabel, nshort, nlong = lslt.shape
    assert nlabel == len(labels)
    if label is None:
        return pprint.pformat([ (labels[i], list(lslt[i,-1,:].real),) for i in range(nlabel) ])
    else:
        for i in range(nlabel):
            if label == labels[i]:
                return pprint.pformat((labels[i], list(lslt[i,-1,:].real),))
    return None

@q.timer
def get_psrc_prop(job_tag, traj, xg, inv_type, inv_acc, *, sfr, fsel):
    tag = mk_psrc_tag(xg, inv_type, inv_acc)
    if tag not in sfr:
        return None
    s_prop = q.SelProp(fsel)
    s_prop.load_double_from_float(sfr, tag)
    return s_prop

# ----

@q.timer_verbose
def load_or_compute_muon_line_interpolation():
    """
    Actually load (compute if data does not exist) muon_line_interpolation.
    """
    fname = q.get_fname()
    fn = "huge-data-muon-line-interpolation-data"
    path_qar = get_load_path(f"{fn}.qar")
    if path_qar is None:
        # make a quick sample data
        q.displayln_info(-1, f"WARNING: {fname}: '{fn}' does not exist. Create a test sample data instead.")
        path = get_save_path(fn)
        epsabs = 1e-8
        epsrel = 1e-3
        mineval = 1024 * 1024
        maxeval = 1024 * 1024 * 1024
        eps0 = (epsabs * 1e2, epsrel * 1e2, mineval // 1024, maxeval // 1024,)
        q.compute_save_muonline_interpolation(
                f"{path}/{0:010d}",
                [ 3, 2, 2, 2, 2, ],
                eps0)
        q.compute_save_muonline_interpolation(
                f"{path}/{1:010d}",
                [ 2, 2, 2, 2, 2, ],
                eps0)
        q.qar_create_info(f"{path}.qar", path, is_remove_folder_after=True)
        q.set_muon_line_m_extra_weights([ [ 4/3, -1/3, ], [ 1.0, ], ])
        q.load_multiple_muonline_interpolations(path, [ 0, 1, ])
    else:
        # load actual data
        assert path_qar[-4:] == ".qar"
        path = path_qar[:-4]
        q.load_multiple_muonline_interpolations(path, [ 1, 3, 5, ])
    q.sync_node();
    weights = q.get_muon_line_m_extra_weights()
    num_muon_line_interps = q.get_number_of_muon_line_interpolations()
    q.displayln_info(-1, f"{fname}: extrapolation weights={weights} loaded num_muon_line_interps={num_muon_line_interps}")

is_loaded_multiple_muonline_interpolations = False

def force_load_muon_line_interpolation():
    """
    Ensure muon line interpolation data is loaded.
    Only load the data if not loaded.
    """
    global is_loaded_multiple_muonline_interpolations
    if is_loaded_multiple_muonline_interpolations:
        return
    load_or_compute_muon_line_interpolation()
    is_loaded_multiple_muonline_interpolations = True

# ----

def get_prob_func(job_tag, inv_type, r_sq_limit, r_sq):
    if job_tag == "48I" and inv_type == 0:
        if r_sq > r_sq_limit:
            prob = 0.0
        elif r_sq == 0:
            prob = 1.0 / 128.0
        elif r_sq <= 8 * 8:
            prob = 1.0
        else:
            prob = (8.0 / np.sqrt(r_sq))**3
    elif job_tag == "64I" and inv_type == 0:
        if r_sq > r_sq_limit:
            prob = 0.0
        elif r_sq == 0:
            prob = 1.0 / 128.0
        elif r_sq <= 20 * 20:
            prob = 1.0
        else:
            prob = (20.0 / np.sqrt(r_sq))**3
    elif job_tag == "24D" and inv_type == 0:
        if r_sq > r_sq_limit:
            prob = 0.0
        elif r_sq == 0:
            prob = 1.0 / 128.0
        elif r_sq <= 10 * 10:
            prob = 1.0
        else:
            prob = (10.0 / np.sqrt(r_sq))**3
    elif job_tag == "test-4nt8":
        if r_sq > r_sq_limit:
            prob = 0.0
        elif r_sq == 0:
            prob = 1.0 / 16.0
        elif r_sq <= 1.5 * 1.5:
            prob = 1.0 / 3.0
        else:
            prob = 1.0 / 3.0 * (1.5 / np.sqrt(r_sq))**3
    else:
        if r_sq > r_sq_limit:
            prob = 0.0
        elif r_sq <= 4 * 4:
            prob = 1.0
        elif r_sq <= 8 * 8:
            prob = 0.5
        else:
            prob = (6.0 / np.sqrt(r_sq))**3
    assert 0 <= prob
    assert prob <= 1
    return prob

def mk_hlbl_four_get_prob(job_tag, inv_type):
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    r_sq_limit = get_r_sq_limit(job_tag)
    factor = 1.0
    if inv_type == 0:
        factor *= get_param(job_tag, "hlbl_four_prob_scaling_factor")
    elif inv_type == 1:
        factor *= get_param(job_tag, "hlbl_four_prob_scaling_factor_strange")
    else:
        assert False
    def get_prob(xg):
        xg_rel = q.smod_coordinate(xg, total_site)
        r_sq = xg_rel.sqr()
        prob = factor * get_prob_func(job_tag, inv_type, r_sq_limit, r_sq)
        return prob
    return get_prob

@q.timer
def get_total_prob(total_site, get_prob):
    """
    get_prob(xg) = prob
    """
    geo = q.Geometry(total_site)
    f_prob = q.FieldRealD(geo)
    f_prob_v = f_prob[:]
    local_volume = geo.local_volume()
    xg_arr = geo.xg_arr()
    assert len(xg_arr) == local_volume
    for index in range(local_volume):
        xg = q.Coordinate(xg_arr[index])
        prob = get_prob(xg)
        assert isinstance(prob, float)
        f_prob_v[index] = prob
    total_prob = f_prob.glb_sum()
    total_prob = total_prob[0, 0]
    return total_prob

@q.timer
def get_hlbl_four_total_prob(job_tag, inv_type):
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    get_prob = mk_hlbl_four_get_prob(job_tag, inv_type)
    total_prob = get_total_prob(total_site, get_prob)
    return total_prob

@q.timer
def mk_hlbl_four_point_pairs(job_tag, traj, *, inv_type, get_psel_prob):
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    geo = q.Geometry(total_site)
    total_volume = geo.total_volume
    inv_type_name_list = [ "light", "strange", ]
    inv_type_name = inv_type_name_list[inv_type]
    #
    psel_prob = get_psel_prob()
    psel = psel_prob.psel
    xg_arr = psel.xg_arr
    n_xg_arr = len(xg_arr)
    #
    get_prob = mk_hlbl_four_get_prob(job_tag, inv_type)
    # total_prob = get_total_prob(total_site, get_prob)
    #
    point_pairs = []
    rs = q.RngState(f"seed {job_tag} {traj}").split(f"mk_hlbl_four_point_pairs")
    for i in range(n_xg_arr):
        xg_x = q.Coordinate(xg_arr[i])
        for j in range(i + 1):
            xg_y = q.Coordinate(xg_arr[j])
            xg_diff = xg_x - xg_y
            prob = get_prob(xg_diff)
            prob_accept = prob
            if prob_accept == 0:
                continue
            rsi = rs.split(f"{j} {i}")
            r = rsi.u_rand_gen()
            if r <= prob_accept:
                dict_val = dict()
                dict_val["idx_xg_x"] = i
                dict_val["idx_xg_y"] = j
                dict_val["xg_x"] = xg_x
                dict_val["xg_y"] = xg_y
                dict_val["r"] = get_r_coordinate(xg_diff, total_site)
                dict_val["prob_accept"] = prob_accept
                if xg_x == xg_y:
                    dict_val["weight_pair"] = 1.0 / prob_accept
                else:
                    # NOTE: need to account for contribution for j > i which is not included in the loop.
                    dict_val["weight_pair"] = 2.0 / prob_accept
                point_pairs.append(dict_val)
    point_pairs = q.random_permute(point_pairs, q.RngState(f"mk_hlbl_four_point_pairs {job_tag} {traj} {inv_type}"))
    q.displayln_info(f"mk_hlbl_four_point_pairs: {job_tag} {traj} {inv_type_name} len(point_pairs)={len(point_pairs)}")
    return point_pairs

@q.timer
def run_hlbl_four_point_pairs_info(job_tag, traj, *, inv_type, get_psel_prob):
    """
    return get_point_pairs
    #
    point_paris = get_point_pairs()
    """
    fname = q.get_fname()
    inv_type_name_list = [ "light", "strange", ]
    inv_type_name = inv_type_name_list[inv_type]
    fn = f"{job_tag}/hlbl/clbl-{inv_type_name}/traj-{traj}/point-pairs.pickle"
    @q.lazy_call
    def load_point_pairs():
        return q.load_pickle_obj(get_load_path(fn))
    ret = load_point_pairs
    if get_load_path(fn) is not None:
        return ret
    if not q.obtain_lock(f"locks/{job_tag}-{traj}-{fname}-{inv_type_name}"):
        return None
    q.timer_fork()
    point_pairs = mk_hlbl_four_point_pairs(
            job_tag, traj,
            inv_type=inv_type,
            get_psel_prob=get_psel_prob,
            )
    q.save_pickle_obj(point_pairs, get_save_path(fn))
    q.release_lock()
    q.timer_display()
    q.timer_merge()
    return ret

def get_hlbl_clbl_info_ref_tags(job_tag):
    return [ "ref-far", "ref-close", "ref-center", ]

@q.timer
def contract_hlbl_four_labels(job_tag):
    tags = get_hlbl_clbl_info_ref_tags(job_tag)
    return q.contract_four_pair_labels(tags)

@q.timer_verbose
def contract_hlbl_four_ama(
        job_tag,
        *,
        inv_type,
        idx_xg_x,
        idx_xg_y,
        psel_prob,
        weight_pair,
        prob_pair,
        psel_d_prob_xy,
        ama_sc_xy,
        ama_sc_yx,
        ama_cm_xy,
        ama_cm_yx,
        ):
    """
    get_prop(xg) => sprop_ama
    """
    psel = psel_prob.psel
    psel_d = psel_d_prob_xy.psel
    xg_x = psel.coordinate_from_idx(idx_xg_x)
    xg_y = psel.coordinate_from_idx(idx_xg_y)
    muon_mass = get_muon_mass(job_tag)
    coef = complex(weight_pair / prob_pair)
    force_load_muon_line_interpolation()
    # q.displayln_info(f"INFO: contract_hlbl_four_ama: {psel_d.geo.total_site}")
    smf_d = q.mk_m_z_field_tag(psel_d, xg_x, xg_y, a=muon_mass, tag=0)
    tags = get_hlbl_clbl_info_ref_tags(job_tag)
    r_sq_limit = get_r_sq_limit(job_tag)
    zz_vv = get_param(job_tag, "zz_vv")
    def f(sc_xy, sc_yx, cm_xy, cm_yx):
        return q.contract_four_pair_no_glb_sum(
                coef,
                psel_prob,
                psel_d_prob_xy,
                idx_xg_x,
                idx_xg_y,
                smf_d,
                sc_xy,
                sc_yx,
                cm_xy,
                cm_yx,
                inv_type,
                tags,
                r_sq_limit,
                muon_mass,
                zz_vv,
                )
    ama_val = ama_apply(f, ama_sc_xy, ama_sc_yx, ama_cm_xy, ama_cm_yx)
    return (ama_extract(ama_val, is_sloppy=False), ama_extract(ama_val, is_sloppy=True),)

@q.timer
def run_hlbl_four_chunk(job_tag, traj, *, inv_type, get_psel_prob, get_fsel_prob, get_point_pairs, prop_cache, id_chunk_list, num_chunk):
    fname = q.get_fname()
    inv_type_name_list = [ "light", "strange", ]
    inv_type_name = inv_type_name_list[inv_type]
    def mk_fn(id_chunk):
        fn = f"{job_tag}/hlbl/clbl-pairs-data-{inv_type_name}/traj-{traj}/chunk_{id_chunk}_{num_chunk}.pickle"
        return fn
    is_all_done = True
    for id_chunk in id_chunk_list:
        fn = mk_fn(id_chunk)
        if get_load_path(fn) is None:
            is_all_done = False
            break
    if is_all_done:
        q.displayln_info(f"{fname}: {job_tag} {traj} {inv_type_name} {id_chunk_list} {num_chunk} all done.")
        return
    q.check_stop()
    q.check_time_limit()
    assert len(id_chunk_list) > 0
    if not q.obtain_lock(f"locks/{job_tag}-{traj}-{fname}-{inv_type_name}-{id_chunk_list[0]}-{num_chunk}"):
        return
    q.timer_fork()
    #
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    geo = q.Geometry(total_site)
    total_volume = geo.total_volume
    #
    hlbl_four_contract_sparse_ratio = get_param(job_tag, "hlbl_four_contract_sparse_ratio")
    #
    psel_prob = get_psel_prob()
    psel = psel_prob.psel
    fsel_prob = get_fsel_prob()
    fsel = fsel_prob.fsel
    #
    ssp = q.SelectedShufflePlan(q.PointsSelection(fsel), q.RngState(f"{job_tag}-{traj}-hlbl-four-fsel-permute"))
    psel_d_prob = q.SelectedPointsRealD(fsel_prob, ssp)
    #
    point_pairs = get_point_pairs()
    all_point_pairs_chunk_list = q.get_chunk_list(point_pairs, chunk_number=num_chunk)
    point_pairs_chunk_list = []
    for id_chunk in id_chunk_list:
        if id_chunk < len(all_point_pairs_chunk_list):
            point_pairs_chunk = all_point_pairs_chunk_list[id_chunk]
        else:
            point_pairs_chunk = []
        point_pairs_chunk_list.append(point_pairs_chunk)
    #
    rel_acc_list = [ 0, 1, 2, ]
    prob_list = [ get_param(job_tag, f"prob_acc_{inv_acc}_psrc") for inv_acc in rel_acc_list ]
    #
    sfr = q.open_fields(get_load_path(f"{job_tag}/prop-psrc-{inv_type_name}/traj-{traj}/geon-info.txt"), "r")
    #
    @q.timer
    def get_prop_cache(xg, inv_acc):
        key = (xg.to_tuple(), inv_type, inv_acc,)
        if key not in prop_cache:
            sprop = get_psrc_prop(job_tag, traj, xg, inv_type, inv_acc, sfr=sfr, fsel=fsel)
            if sprop is None:
                psprop = None
            else:
                psprop = q.PselProp(sprop, ssp)
            prop_cache[key] = psprop
        return prop_cache[key]
    #
    @q.timer
    def get_prop(xg):
        val_list = [ get_prop_cache(xg, inv_acc) for inv_acc in rel_acc_list ]
        return mk_ama_val(val_list[0], xg.to_tuple(), val_list, rel_acc_list, prob_list)
    #
    @q.timer
    def mk_ama_current(ama_sprop1, ama_sprop2):
        ama_val = ama_apply2(q.mk_local_current_from_props, ama_sprop1, ama_sprop2)
        return ama_val
    #
    @q.timer
    def mk_psel_d_prob_xy(idx_xg_x, idx_xg_y):
        prob_pair, psel_d_prob_xy = q.mk_psel_d_prob_xy(psel_prob, psel_d_prob, idx_xg_x, idx_xg_y)
        return prob_pair, psel_d_prob_xy
    #
    @q.timer
    def mk_ama_cm(ama_current, psel_d_prob_xy):
        def f(current):
            cm = q.CurrentMoments(current, psel_d_prob_xy)
            return cm
        return ama_apply1(f, ama_current)
    #
    @q.timer
    def ama_cm_glb_sum(ama_cm):
        def f(cm):
            return cm.glb_sum()
        return ama_apply1(f, ama_cm)
    #
    q.displayln_info(f"{fname}: len(prop_cache)={len(prop_cache)}")
    labels = contract_hlbl_four_labels(job_tag)
    for id_chunk, point_pairs_chunk in zip(id_chunk_list, point_pairs_chunk_list):
        if get_load_path(mk_fn(id_chunk)) is not None:
            q.displayln_info(f"{fname}: {job_tag} {traj} {inv_type_name} {id_chunk}/{num_chunk} done.")
            continue
        q.check_stop()
        q.check_time_limit()
        pairs_data = []
        len_chunk = len(point_pairs_chunk)
        len_info_str = f"len_chunk={len_chunk} ; id_chunk={id_chunk} ; num_chunk={num_chunk}"
        for idx, pp in enumerate(point_pairs_chunk):
            q.displayln_info(f"{fname}: load prop ; idx={idx} ; {len_info_str}")
            xg_x = pp["xg_x"]
            xg_y = pp["xg_y"]
            get_prop(xg_x)
            get_prop(xg_y)
        pairs_int_results = dict()
        for idx, pp in enumerate(point_pairs_chunk):
            q.displayln_info(f"{fname}: make intermediate results; idx={idx} ; {len_info_str}")
            idx_xg_x = pp["idx_xg_x"]
            idx_xg_y = pp["idx_xg_y"]
            xg_x = pp["xg_x"]
            xg_y = pp["xg_y"]
            assert xg_x == psel.coordinate_from_idx(idx_xg_x)
            assert xg_y == psel.coordinate_from_idx(idx_xg_y)
            key = (idx_xg_x, idx_xg_y,)
            prob_pair, psel_d_prob_xy = mk_psel_d_prob_xy(idx_xg_x, idx_xg_y)
            ama_sprop_x = get_prop(xg_x)
            ama_sprop_y = get_prop(xg_y)
            ama_sc_xy = mk_ama_current(ama_sprop_y, ama_sprop_x)
            ama_cm_xy = mk_ama_cm(ama_sc_xy, psel_d_prob_xy)
            sp_norm = q.qnorm_field(ama_extract(ama_sc_xy, is_sloppy=False))
            ama_sc_yx = mk_ama_current(ama_sprop_x, ama_sprop_y)
            ama_cm_yx = mk_ama_cm(ama_sc_yx, psel_d_prob_xy)
            sp_norm += q.qnorm_field(ama_extract(ama_sc_yx, is_sloppy=False))
            sp_norm = q.sqrt_field(sp_norm)
            int_results = dict()
            int_results['prob_pair'] = prob_pair
            int_results['psel_d_prob_xy'] = psel_d_prob_xy
            int_results['ama_sc_xy'] = ama_sc_xy
            int_results['ama_sc_yx'] = ama_sc_yx
            int_results['ama_cm_xy'] = ama_cm_xy
            int_results['ama_cm_yx'] = ama_cm_yx
            int_results['sp_norm'] = sp_norm
            pairs_int_results[key] = int_results
        @q.timer_verbose
        def sync_node_after_hlbl_four_current():
            q.sync_node()
        sync_node_after_hlbl_four_current()
        for idx, pp in enumerate(point_pairs_chunk):
            q.displayln_info(f"{fname}: cm glb sum; idx={idx} ; {len_info_str}")
            idx_xg_x = pp["idx_xg_x"]
            idx_xg_y = pp["idx_xg_y"]
            xg_x = pp["xg_x"]
            xg_y = pp["xg_y"]
            prob_accept = pp["prob_accept"]
            weight_pair = pp["weight_pair"]
            r = pp["r"]
            key = (idx_xg_x, idx_xg_y,)
            int_results = pairs_int_results[key]
            prob_pair = int_results['prob_pair']
            psel_d_prob_xy = int_results['psel_d_prob_xy']
            ama_sc_xy = int_results['ama_sc_xy']
            ama_sc_yx = int_results['ama_sc_yx']
            ama_cm_xy = int_results['ama_cm_xy']
            ama_cm_yx = int_results['ama_cm_yx']
            ama_cm_xy = ama_cm_glb_sum(ama_cm_xy)
            ama_cm_yx = ama_cm_glb_sum(ama_cm_yx)
            int_results['ama_cm_xy'] = ama_cm_xy
            int_results['ama_cm_yx'] = ama_cm_yx
            @q.timer_verbose
            def hlbl_four_contract_sparse():
                sf_pair_f_rand_01 = q.SelectedFieldRealD(fsel, 1)
                sf_pair_f_rand_01.set_rand(q.RngState(f"{job_tag} {traj} {inv_type} {xg_x.to_tuple()} {xg_y.to_tuple()}"), 1.0, 0.0)
                sp_pair_f_rand_01 = q.SelectedPointsRealD(sf_pair_f_rand_01, ssp)
                assert len(sp_pair_f_rand_01) == len(psel_d_prob_xy)
                psel_d = psel_d_prob_xy.psel
                sp_norm = int_results['sp_norm'].copy()
                sp_norm *= weight_pair / prob_pair
                sp_norm[:] = sp_norm[:] / psel_d_prob_xy[:]
                glb_avg = q.glb_sum(sp_norm[:].sum()) / q.glb_sum(len(sp_norm))
                q.displayln_info(f"{fname}: {inv_type_name} ; r={r} ; weight_pair={weight_pair} ; prob_pair={prob_pair} ; sp_norm_avg={glb_avg}")
                # q.displayln_info(f"INFO: {fname}: sp_norm=\n{sp_norm[:, 0]}")
                sp_norm[:] = np.minimum(1.0, sp_norm[:] / glb_avg / hlbl_four_contract_sparse_ratio)
                selection = sp_pair_f_rand_01[:, 0] <= sp_norm[:, 0]
                psel_d_sel = q.PointsSelection(psel_d.total_site, psel_d[selection, :], psel_d.points_dist_type)
                assert psel_d_sel.points_dist_type == psel_d.points_dist_type
                psel_d_sel_prob_xy = q.SelectedPointsRealD(psel_d_sel, 1)
                psel_d_sel_prob_xy @= psel_d_prob_xy
                psel_d_sel_prob_xy[:] *= sp_norm[selection, :]
                n_avail = q.glb_sum(len(psel_d))
                n_sel = q.glb_sum(len(psel_d_sel))
                sel_ratio = n_sel / n_avail
                q.displayln_info(f"{fname}: {inv_type_name} ; n_avail={n_avail} ; n_sel={n_sel} ; ratio={sel_ratio}")
                def sel_sc(sc):
                    multiplicity = sc.multiplicity
                    assert multiplicity == 4
                    sc_sel = q.SelectedPointsWilsonMatrix(psel_d_sel, multiplicity)
                    sc_sel @= sc
                    return sc_sel
                ama_sel_sc_xy = ama_apply1(sel_sc, ama_sc_xy)
                ama_sel_sc_yx = ama_apply1(sel_sc, ama_sc_yx)
                int_results['psel_d_prob_xy'] = psel_d_sel_prob_xy
                int_results['ama_sc_xy'] = ama_sel_sc_xy
                int_results['ama_sc_yx'] = ama_sel_sc_yx
            hlbl_four_contract_sparse()
        for idx, pp in enumerate(point_pairs_chunk):
            q.displayln_info(f"{fname}: contract ; idx={idx} ; {len_info_str}")
            idx_xg_x = pp["idx_xg_x"]
            idx_xg_y = pp["idx_xg_y"]
            xg_x = pp["xg_x"]
            xg_y = pp["xg_y"]
            r = pp["r"]
            prob_accept = pp["prob_accept"]
            weight_pair = pp["weight_pair"]
            key = (idx_xg_x, idx_xg_y,)
            int_results = pairs_int_results[key]
            pairs_int_results[key] = None
            prob_pair = int_results['prob_pair']
            psel_d_prob_xy = int_results['psel_d_prob_xy']
            ama_sc_xy = int_results['ama_sc_xy']
            ama_sc_yx = int_results['ama_sc_yx']
            ama_cm_xy = int_results['ama_cm_xy']
            ama_cm_yx = int_results['ama_cm_yx']
            lslt, lslt_sloppy = contract_hlbl_four_ama(
                    job_tag,
                    inv_type=inv_type,
                    idx_xg_x=idx_xg_x,
                    idx_xg_y=idx_xg_y,
                    psel_prob=psel_prob,
                    weight_pair=weight_pair,
                    prob_pair=prob_pair,
                    psel_d_prob_xy=psel_d_prob_xy,
                    ama_sc_xy=ama_sc_xy,
                    ama_sc_yx=ama_sc_yx,
                    ama_cm_xy=ama_cm_xy,
                    ama_cm_yx=ama_cm_yx,
                    )
            dict_val = dict()
            dict_val["lslt"] = lslt
            dict_val["lslt_sloppy"] = lslt_sloppy
            dict_val["xg_x"] = xg_x
            dict_val["xg_y"] = xg_y
            dict_val["r"] = r
            dict_val["prob_accept"] = prob_accept
            dict_val["weight_pair"] = weight_pair
            pairs_data.append(dict_val)
        @q.timer_verbose
        def sync_node_after_hlbl_four_contract():
            q.sync_node()
        sync_node_after_hlbl_four_contract()
        for idx, pp in enumerate(point_pairs_chunk):
            q.displayln_info(f"{fname}: collect results ; idx={idx} ; {len_info_str}")
            dict_val = pairs_data[idx]
            lslt = q.glb_sum(dict_val["lslt"])
            lslt_sloppy = q.glb_sum(dict_val["lslt_sloppy"])
            dict_val["lslt"] = lslt
            dict_val["lslt_sloppy"] = lslt_sloppy
            xg_x = dict_val["xg_x"]
            xg_y = dict_val["xg_y"]
            r = dict_val["r"] = r
            prob_accept = dict_val["prob_accept"]
            weight_pair = dict_val["weight_pair"]
            info_str = f"{job_tag} {traj} {inv_type_name} {id_chunk}/{num_chunk} {idx}/{len(point_pairs_chunk)} {xg_x} {xg_y}"
            q.displayln_info(
                    f"{fname}: {info_str}\n",
                    f"r={r} weight_pair={weight_pair}\n",
                    show_lslt(labels, lslt * len(point_pairs), label="ref-far proj-all"),
                    "\nsloppy:\n",
                    show_lslt(labels, lslt_sloppy * len(point_pairs), label="ref-far proj-all"))
            json_results.append((
                f"{fname}: {info_str} lslt",
                q.get_data_sig(lslt, q.RngState()),
                15e-2,
                ))
            json_results.append((
                f"{fname}: {info_str} lslt_sloppy",
                q.get_data_sig(lslt_sloppy, q.RngState()),
                15e-2,
                ))
        if len(point_pairs_chunk) != len(pairs_data):
            raise Exception(f"len(point_pairs_chunk)={len(point_pairs_chunk)} len(pairs_data)={len(pairs_data)}")
        if len(pairs_data) > 0:
            q.displayln_info(f"{fname}: {job_tag} {traj} {inv_type_name} {id_chunk}/{num_chunk}\n",
                    show_lslt(labels, sum([ d["lslt"] for d in pairs_data ]) / len(point_pairs_chunk) * len(point_pairs)))
        q.save_pickle_obj(pairs_data, get_save_path(mk_fn(id_chunk)))
    sfr.close()
    q.displayln_info(f"{fname}: len(prop_cache)={len(prop_cache)}")
    q.release_lock()
    q.timer_display()
    q.timer_merge()

@q.timer
def run_hlbl_four(job_tag, traj, *, inv_type, get_psel_prob, get_fsel_prob, get_point_pairs):
    fname = q.get_fname()
    inv_type_name_list = [ "light", "strange", ]
    inv_type_name = inv_type_name_list[inv_type]
    if get_point_pairs is None:
        q.displayln_info(f"{fname}: {job_tag} {traj} {inv_type_name} get_point_pairs is None")
        return
    fn_s = f"{job_tag}/hlbl/clbl-{inv_type_name}/traj-{traj}/results-brief.pickle"
    if get_load_path(fn_s) is not None:
        return
    num_chunk = get_param(job_tag, "hlbl_four_num_chunk")
    prop_cache = q.mk_cache(f"{fname}-prop_cache", job_tag, f"{traj}")
    id_chunk_list_list = q.get_chunk_list(list(range(num_chunk)), chunk_number=4)
    for id_chunk_list in id_chunk_list_list:
        run_hlbl_four_chunk(
                job_tag, traj,
                inv_type=inv_type,
                get_psel_prob=get_psel_prob,
                get_fsel_prob=get_fsel_prob,
                get_point_pairs=get_point_pairs,
                prop_cache=prop_cache,
                id_chunk_list=id_chunk_list,
                num_chunk=num_chunk,
                )
    q.clean_cache(prop_cache)
    fn_chunk_list = []
    for id_chunk in range(num_chunk):
        fn_chunk = f"{job_tag}/hlbl/clbl-pairs-data-{inv_type_name}/traj-{traj}/chunk_{id_chunk}_{num_chunk}.pickle"
        if get_load_path(fn_chunk) is None:
            q.displayln_info(f"{fname}: {job_tag} {traj} {inv_type_name} not finished (missing '{fn_chunk})'")
            return
        else:
            fn_chunk_list.append(fn_chunk)
    assert len(fn_chunk_list) == num_chunk
    if get_load_path(fn_s) is not None:
        return
    if not q.obtain_lock(f"locks/{job_tag}-{traj}-{fname}-{inv_type_name}"):
        return
    q.sync_node()
    q.timer_fork()
    labels = contract_hlbl_four_labels(job_tag)
    pairs_data_n_pairs = 0
    pairs_data_lslt_sum = 0
    pairs_data_lslt_sloppy_sum = 0
    for fn_chunk in fn_chunk_list:
        path_chunk = get_load_path(fn_chunk)
        assert path_chunk is not None
        if q.get_id_node() == 0:
            pairs_data_chunk = q.load_pickle_obj(path_chunk, is_sync_node=False)
            pairs_data_n_pairs += len(pairs_data_chunk)
            pairs_data_lslt_sum += sum([ d["lslt"] for d in pairs_data_chunk ])
            pairs_data_lslt_sloppy_sum += sum([ d["lslt_sloppy"] for d in pairs_data_chunk ])
    q.sync_node()
    q.displayln_info(-1, f"{fname}: {job_tag} {traj} {inv_type_name} load all chunk done.")
    pairs_data_n_pairs_ref = len(get_point_pairs())
    if q.get_id_node() == 0:
        if pairs_data_n_pairs != pairs_data_n_pairs_ref:
            q.displayln_info(-1, f"ERROR: {fname}: pairs_data_n_pairs={pairs_data_n_pairs} len(get_point_pairs())={pairs_data_n_pairs_ref}")
            raise Exception(f"pairs_data_n_pairs={pairs_data_n_pairs} len(get_point_pairs())={pairs_data_n_pairs_ref}")
    q.sync_node()
    if q.get_id_node() == 0:
        results = dict()
        results["labels"] = labels
        results["lslt_sum"] = pairs_data_lslt_sum
        results["lslt_sloppy_sum"] = pairs_data_lslt_sloppy_sum
        results["n_pairs"] = pairs_data_n_pairs
        q.save_pickle_obj(results, get_save_path(fn_s))
        q.displayln_info(-1, f"{fname}: {job_tag} {traj} {inv_type_name}\n", show_lslt(labels, results["lslt_sum"]))
        json_results.append((
            f"{fname}: {job_tag} {traj} {inv_type_name} lslt_sum",
            q.get_data_sig(results["lslt_sum"], q.RngState()),
            1e-2,
            ))
        json_results.append((
            f"{fname}: {job_tag} {traj} {inv_type_name} lslt_sloppy_sum",
            q.get_data_sig(results["lslt_sloppy_sum"], q.RngState()),
            1e-2,
            ))
        json_results.append((
            f"{fname}: {job_tag} {traj} {inv_type_name} lslt_sum[labels.index('ref-far proj-all'), -1, -1]",
            results["lslt_sum"][labels.index('ref-far proj-all'), -1, -1],
            3e-2,
            ))
        json_results.append((
            f"{fname}: {job_tag} {traj} {inv_type_name} lslt_sloppy_sum[labels.index('ref-far proj-all'), -1, -1]",
            results["lslt_sloppy_sum"][labels.index('ref-far proj-all'), -1, -1],
            3e-2,
            ))
    q.sync_node()
    q.displayln_info(f"{fname}: {job_tag} {traj} {inv_type_name} done.")
    q.release_lock()
    q.timer_display()
    q.timer_merge()
    return [ f"{fname}: {job_tag} {traj} {inv_type_name} done.", ]

# ----

@q.timer
def run_hvp_sum_tslice_accs(job_tag, traj, *, inv_type, get_psel):
    """
    (1) mu is the sink polarization and nu is the src polarization
    (2) hvp field is simply the trace of the products of gamma matrix and propagators.
        It does not include the any minus sign (e.g. The minus sign due to the loop).
    #
    accs stand for information of many accuracy levels are included
    #
    return get_hvp_sum_tslice_accs
    #
    hvp_sum_tslice_accs = get_hvp_sum_tslice_accs()
    nparr = hvp_sum_tslice_accs[tag]
    nparr[t_dir, tslice, mu, nu] is complex
    """
    fname = q.get_fname()
    inv_type_name_list = [ "light", "strange", ]
    inv_type_name = inv_type_name_list[inv_type]
    if get_psel is None:
        q.displayln_info(-1, f"{fname}: {job_tag} {traj} {inv_type_name} get_psel is None.")
        return None
    path = f"{job_tag}/hvp-sum-tslice-psrc-{inv_type_name}/traj-{traj}"
    if get_load_path(f"{path}/checkpoint.txt") is None:
        q.displayln_info(-1, f"{fname}: '{path}' data is not ready.")
        return None
    lpath = get_load_path(f"{path}")
    assert lpath is not None
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    t_size = total_site[3]
    #
    @q.lazy_call
    @q.timer_verbose
    def get_hvp_sum_tslice_accs():
        hvp_sum_tslice_accs = dict()
        psel = get_psel()
        rel_acc_list = [ 0, 1, 2, ]
        for xg in psel:
            for inv_acc in rel_acc_list:
                tag = mk_psrc_tag(xg, inv_type, inv_acc)
                fn = f"{lpath}/{tag}.lat"
                if not q.does_file_exist_qar_sync_node(fn):
                    continue
                ld = q.load_lat_data(fn)
                arr = ld.to_numpy()
                assert arr.shape == (4, t_size, 4, 4,)
                hvp_sum_tslice_accs[tag] = arr
        q.displayln_info(f"{fname}: {job_tag} {traj} {inv_type_name} loaded.")
        return hvp_sum_tslice_accs
    return get_hvp_sum_tslice_accs

@q.timer
def run_hvp_sum_tslice(job_tag, traj, *, inv_type, get_psel, get_hvp_sum_tslice_accs):
    """
    return the ama corrected: get_hvp_sum_tslice()
    hvp_sum_tslice[idx, t_dir, tslice, mu, nu] is complex
    hvp_sum_tslice.shape == (len(psel), 4, t_size, 4, 4)
    #
    get_hvp_sum_tslice_accs = run_hvp_sum_tslice_accs(job_tag, traj, inv_type=inv_type, get_psel=get_psel)
    """
    fname = q.get_fname()
    if get_hvp_sum_tslice_accs is None:
        q.displayln_info(-1, f"{fname}: get_hvp_sum_tslice_accs is None")
        return None
    if get_psel is None:
        q.displayln_info(-1, f"{fname}: get_psel is None")
        return None
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    t_size = total_site[3]
    @q.timer_verbose
    @q.lazy_call
    def get_hvp_sum_tslice():
        psel = get_psel()
        hvp_sum_tslice_accs = get_hvp_sum_tslice_accs()
        hvp_sum_tslice = np.zeros((len(psel), 4, t_size, 4, 4,), dtype=np.complex128)
        rel_acc_list = [ 0, 1, 2, ]
        prob_list = [ get_param(job_tag, f"prob_acc_{inv_acc}_psrc") for inv_acc in rel_acc_list ]
        for idx, xg in enumerate(psel):
            val_list = [ hvp_sum_tslice_accs.get(mk_psrc_tag(xg, inv_type, inv_acc)) for inv_acc in rel_acc_list ]
            ama_val = mk_ama_val(val_list[0], xg.to_tuple(), val_list, rel_acc_list, prob_list)
            val = ama_extract(ama_val)
            hvp_sum_tslice[idx] = val
        q.displayln_info(-1, f"{fname}: {job_tag} {traj} {inv_type} hvp_sum_tslice.shape={hvp_sum_tslice.shape}")
        return hvp_sum_tslice
    return get_hvp_sum_tslice

def get_edl_from_hvp_sum_tslice(xg, total_site, hvp_sum_tslice):
    """
    edl[k, nu] is complex
    include -1 from fermion loop (not yet included in this hvp)
    include 1/2 ii for the magnetic moment projection
    include - ii for the current couple to internal photon
    #
    hvp_sum_tslice[t_dir, tslice, mu, nu]
    #
    inner_products[i, mu, nu]
    #
    edl[k, nu] = (-1 * 0.5j * -1j) * sum_{i, j} epsilon_{i, j, k} * (t_arr - xg[i]) * hvp_sum_tslice[i, t_arr, j, nu]
    """
    assert isinstance(xg, q.Coordinate)
    edl = np.zeros((3, 4,), dtype=np.complex128)
    inner_products = np.zeros((3, 4, 4), dtype=np.complex128)
    t_size = hvp_sum_tslice.shape[1]
    t_arr = np.arange(t_size)
    assert hvp_sum_tslice.shape == (4, t_size, 4, 4,)
    for t_dir in range(3):
        x_t = xg[t_dir]
        size = total_site[t_dir]
        xrel = q.rel_mod_sym_arr(t_arr - x_t, size)[:, None, None]
        inner_products[t_dir] = (xrel * hvp_sum_tslice[t_dir]).sum(axis=0)
    edl[0] = inner_products[1, 2] - inner_products[2, 1]
    edl[1] = inner_products[2, 0] - inner_products[0, 2]
    edl[2] = inner_products[0, 1] - inner_products[1, 0]
    edl = (-1 * 0.5j * -1j) * edl
    return edl

@q.timer
def run_edl(job_tag, traj, *, inv_type, get_psel, get_hvp_sum_tslice):
    """
    edl stand for external (photon) disconnected loop
    return hvp_edl
    isinstance(hvp_edl, q.SelectedPointsComplexD)
    where hvp_edl[idx, k * 4 + nu] is complex
    hvp_edl[:].shape == (len(psel), 3 * 4,)
    #
    include -1 from fermion loop (not yet included in this hvp)
    include 1/2 ii for the magnetic moment projection
    include - ii for the current couple to internal photon
    include 1 or 1/5 as the hvp_type_charge_factor
    """
    fname = q.get_fname()
    inv_type_name_list = [ "light", "strange", ]
    inv_type_name = inv_type_name_list[inv_type]
    fn = f"{job_tag}/hlbl/edl/traj-{traj}/edl-{inv_type_name}.lat"
    @q.timer_verbose
    @q.lazy_call
    def get_edl():
        psel = get_psel()
        hvp_edl = q.SelectedPointsComplexD(psel)
        hvp_edl.load(get_load_path(fn))
        return hvp_edl
    ret = get_edl
    if get_load_path(fn) is not None:
        return ret
    if get_hvp_sum_tslice is None:
        q.displayln_info(-1, f"{fname} get_hvp_sum_tslice is None.")
        return None
    if get_psel is None:
        q.displayln_info(-1, f"{fname} get_psel is None.")
        return None
    if not q.obtain_lock(f"locks/{job_tag}-{traj}-{fname}-{inv_type_name}"):
        return None
    q.timer_fork()
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    psel = get_psel()
    hvp_sum_tslice = get_hvp_sum_tslice()
    hvp_type_charge_factor_list = [ 1.0, 1.0 / 5.0, ]
    hvp_type_charge_factor = hvp_type_charge_factor_list[inv_type]
    hvp_edl = q.SelectedPointsComplexD(psel, 3 * 4)
    hvp_edl.set_zero()
    hvp_edl_view = hvp_edl[:].reshape(len(psel), 3, 4)
    for idx, xg in enumerate(psel):
        hvp_edl_view[idx] = (
                hvp_type_charge_factor
                * get_edl_from_hvp_sum_tslice(xg, total_site, hvp_sum_tslice[idx])
                )
    hvp_edl.save(get_save_path(fn))
    json_results.append((
        f"{fname}: {job_tag} {traj} {inv_type_name} edl",
        q.get_data_sig(hvp_edl[:], q.RngState()),
        1e-3,
        ))
    q.release_lock()
    q.timer_display()
    q.timer_merge()
    return ret

@q.timer_verbose
def run_check_hvp_avg(job_tag, traj, *, inv_type, get_psel_prob, get_hvp_sum_tslice, get_hvp_average):
    fname = q.get_fname()
    inv_type_name_list = [ "light", "strange", ]
    inv_type_name = inv_type_name_list[inv_type]
    fn = f"{job_tag}/hlbl/check-hvp-avg-{inv_type_name}/traj-{traj}.txt"
    if get_load_path(fn) is not None:
        return
    if get_psel_prob is None:
        return
    if get_hvp_sum_tslice is None:
        return
    if get_hvp_average is None:
        return
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    total_site_arr = total_site.to_numpy()
    t_size = total_site[3]
    geo = q.Geometry(total_site)
    psel_prob = get_psel_prob()
    psel = psel_prob.psel
    xg_arr = psel.xg_arr
    hvp_average = get_hvp_average()
    hvp_sum_tslice = get_hvp_sum_tslice()
    idx_arr = np.arange(len(psel))
    t_arr = np.arange(t_size)
    t_dir_arr = np.arange(4)
    tslice_arr = (t_arr + xg_arr[:, :, None]) % total_site_arr[:, None]
    tslice_arr_sel = t_arr < total_site_arr[:, None]
    tslice_arr[:, ~tslice_arr_sel] = np.broadcast_to(np.arange(t_size), (4, t_size,))[~tslice_arr_sel]
    assert tslice_arr.shape == (len(psel), 4, t_size,)
    hvp_sum_tslice = hvp_sum_tslice[idx_arr[:, None, None], t_dir_arr[None, :, None], tslice_arr]
    hvp_sum_tslice_avg1 = np.sum(hvp_sum_tslice / psel_prob[:, None, None, None], axis=0) / geo.total_volume
    hvp_sum_tslice_avg2 = calc_hvp_sum_tslice(hvp_average)[:]
    norm_diff_ratio = (
            2 * np.sqrt(q.qnorm(hvp_sum_tslice_avg1 - hvp_sum_tslice_avg2))
            / (np.sqrt(q.qnorm(hvp_sum_tslice_avg1)) + np.sqrt(q.qnorm(hvp_sum_tslice_avg2)))
            )
    q.displayln_info(-1, f"{fname}: {job_tag} {traj} {inv_type_name} {norm_diff_ratio}")
    assert norm_diff_ratio < 1e-6
    json_results.append((
        f"{fname}: {job_tag} {traj} {inv_type_name} hvp_sum_tslice",
        q.get_data_sig(hvp_sum_tslice, q.RngState()),
        1e-4,
        ))
    json_results.append((
        f"{fname}: {job_tag} {traj} {inv_type_name} hvp_sum_tslice_avg1",
        q.get_data_sig(hvp_sum_tslice_avg1, q.RngState()),
        1e-4,
        ))
    json_results.append((
        f"{fname}: {job_tag} {traj} {inv_type_name} hvp_average",
        q.get_data_sig(hvp_average, q.RngState()),
        1e-4,
        ))
    json_results.append((
        f"{fname}: {job_tag} {traj} {inv_type_name} hvp_sum_tslice_avg2",
        q.get_data_sig(hvp_sum_tslice_avg2, q.RngState()),
        1e-4,
        ))
    q.qtouch_info(get_save_path(fn))

@q.timer_verbose
def run_hlbl_sub_hvp_sfield(
        job_tag,
        traj,
        *,
        inv_type,
        get_psel_prob,
        get_glb_hvp_avg_for_sub,
        get_f_rand_01,
        ):
    """
    Save the glb_avg subtracted, importance sparsed, charge factor multiplied, hvp sparse fields.
    """
    fname = q.get_fname()
    inv_type_name_list = [ "light", "strange", ]
    inv_type_name = inv_type_name_list[inv_type]
    fn = f"{job_tag}/hlbl/sub-hvp-{inv_type_name}/traj-{traj}"
    if get_load_path(fn) is not None:
        return
    if get_glb_hvp_avg_for_sub is None:
        return
    if get_psel_prob is None:
        return
    if get_f_rand_01 is None:
        return
    q.check_stop()
    q.check_time_limit()
    if not q.obtain_lock(f"locks/{job_tag}-{traj}-{fname}-{inv_type_name}"):
        return
    q.timer_fork()
    #
    sfw = q.open_fields(get_save_path(fn + ".acc"), "w", q.Coordinate([ 1, 2, 2, 2, ]))
    # qar_sp = q.open_qar_info(get_save_path(fn), "a")
    #
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    geo = q.Geometry(total_site)
    total_volume = geo.total_volume
    #
    hvp_sel_threshold = get_param(job_tag, "hlbl_two_plus_two_num_hvp_sel_threshold")
    #
    psel_prob = get_psel_prob()
    psel = psel_prob.psel
    #
    # q.FieldRealD(geo, 1) with uniform [0, 1] random number used to sample points.
    f_rand_01 = get_f_rand_01()
    #
    glb_hvp_avg = get_glb_hvp_avg_for_sub()
    #
    hvp_type_charge_factor_list = [ 1.0, 1.0 / 5.0, ]
    hvp_type_charge_factor = hvp_type_charge_factor_list[inv_type]
    #
    xg_arr = psel.xg_arr
    #
    idx_xg_list = list(range(len(psel)))
    #
    sfr = q.open_fields(get_load_path(f"{job_tag}/hvp-psrc-{inv_type_name}/traj-{traj}/geon-info.txt"), "r")
    tags = sfr.list()
    rel_acc_list = [ 0, 1, 2, ]
    prob_list = [ get_param(job_tag, f"prob_acc_{inv_acc}_psrc") for inv_acc in rel_acc_list ]
    for idx in idx_xg_list:
        xg = q.Coordinate(xg_arr[idx])
        val_list = []
        for inv_acc in rel_acc_list:
            tag = mk_psrc_tag(xg, inv_type, inv_acc)
            if tag not in tags:
                val_list.append(None)
            else:
                chvp_16 = q.FieldComplexD(geo, 16)
                chvp_16.load_double_from_float(sfr, tag)
                val_list.append(chvp_16)
        assert val_list[0] is not None
        ama_val = mk_ama_val(val_list[0], xg.to_tuple(), val_list, rel_acc_list, prob_list)
        hvp = ama_extract(ama_val)
        hvp_avg = glb_hvp_avg.copy().shift(xg)
        hvp -= hvp_avg
        hvp_sel_prob = q.sqrt_field(q.qnorm_field(hvp))
        hvp_sel_prob *= 1.0 / hvp_sel_threshold
        hvp_sel_prob[:] = np.minimum(1.0, hvp_sel_prob[:])
        ps_sel = f_rand_01[:, 0] <= hvp_sel_prob[:, 0]
        fsel_ps = q.FieldSelection(geo)
        fsel_ps[ps_sel] = 0
        fsel_ps.update()
        num_fsel_ps = q.glb_sum(fsel_ps.n_elems)
        tag = mk_psrc_tag(xg, inv_type, inv_acc="ama")
        q.displayln_info(0, f"{fname}: {tag} ; num_fsel_ps={num_fsel_ps} ratio={num_fsel_ps/total_volume}")
        assert num_fsel_ps > 0
        fsel_ps_prob = q.SelectedFieldRealD(fsel_ps, 1)
        fsel_ps_prob @= hvp_sel_prob
        fsel_ps_prob.save_double(sfw, f"{tag} ; fsel-prob", skip_if_exist=True)
        s_hvp = q.SelectedFieldComplexD(fsel_ps, 16)
        s_hvp @= hvp
        s_hvp.save_float_from_double(sfw, tag, skip_if_exist=True)
        sfw.flush()
    sfw.close()
    sfr.close()
    q.qrename_info(get_save_path(fn + ".acc"), get_save_path(fn))
    #
    q.release_lock()
    q.timer_display()
    q.timer_merge()

@q.timer_verbose
def run_hlbl_two_plus_two_chunk(
        job_tag,
        traj,
        *,
        inv_type,
        inv_type_e,
        get_psel_prob,
        get_edl_light,
        get_edl_strange,
        get_glb_hvp_avg_for_sub_light,
        get_glb_hvp_avg_for_sub_strange,
        get_f_rand_01,
        id_chunk,
        num_chunk,
        ):
    fname = q.get_fname()
    inv_type_name_list = [ "light", "strange", ]
    inv_type_name = inv_type_name_list[inv_type]
    inv_type_e_name = inv_type_name_list[inv_type_e]
    fn = f"{job_tag}/hlbl/dlbl-points-data-{inv_type_name}-{inv_type_e_name}/traj-{traj}/chunk_{id_chunk}_{num_chunk}.pickle"
    if get_load_path(fn) is not None:
        return
    sub_hvp_fn = f"{job_tag}/hlbl/sub-hvp-{inv_type_name}/traj-{traj}/geon-info.txt"
    if get_load_path(sub_hvp_fn) is None:
        return
    q.check_stop()
    q.check_time_limit()
    if not q.obtain_lock(f"locks/{job_tag}-{traj}-{fname}-{inv_type_name}-{inv_type_e_name}-{id_chunk}_{num_chunk}"):
        return
    q.timer_fork()
    #
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    geo = q.Geometry(total_site)
    #
    r_sq_limit = get_r_sq_limit(job_tag)
    muon_mass = get_muon_mass(job_tag)
    zz_vv = get_param(job_tag, "zz_vv")
    #
    psel_prob = get_psel_prob()
    psel = psel_prob.psel
    #
    # edl already have charge_factor (and many other factors) multiplied
    if inv_type_e == 0:
        edl = get_edl_light()
    elif inv_type_e == 1:
        edl = get_edl_strange()
    else:
        assert False
    #
    # q.FieldRealD(geo, 1) with uniform [0, 1] random number used to sample points.
    f_rand_01 = get_f_rand_01()
    #
    if inv_type == 0:
        glb_hvp_avg = get_glb_hvp_avg_for_sub_light()
    elif inv_type == 1:
        glb_hvp_avg = get_glb_hvp_avg_for_sub_strange()
    else:
        assert False
    #
    hvp_type_charge_factor_list = [ 1.0, 1.0 / 5.0, ]
    hvp_type_charge_factor = hvp_type_charge_factor_list[inv_type]
    #
    xg_arr = psel.xg_arr
    #
    idx_xg_list = list(range(len(psel)))
    idx_xg_list_chunk_list = q.get_chunk_list(
            idx_xg_list,
            chunk_number=num_chunk,
            rng_state=q.RngState(f"xg_list permute for hlbl_two_plus_two"),
            )
    if id_chunk < len(idx_xg_list_chunk_list):
        idx_xg_list_chunk = idx_xg_list_chunk_list[id_chunk]
    else:
        idx_xg_list_chunk = []
    #
    hvp_list = []
    #
    sfr = q.open_fields(get_load_path(sub_hvp_fn), "r")
    tags = sfr.list()
    for idx, idx_xg_x in enumerate(idx_xg_list_chunk):
        q.check_stop()
        q.check_time_limit()
        xg = q.Coordinate(xg_arr[idx_xg_x])
        tag = mk_psrc_tag(xg, inv_type, inv_acc="ama")
        if tag not in tags:
            raise Exception(f"{fname}: idx_xg_x={idx_xg_x} '{tag}' {tags}")
        if f"{tag} ; fsel-prob" not in tags:
            raise Exception(f"{fname}: idx_xg_x={idx_xg_x} '{tag} ; fsel-prob' {tags}")
        fsel_ps_prob = q.SelectedFieldRealD(None)
        fsel_ps_prob.load_double(sfr, f"{tag} ; fsel-prob")
        fsel_ps = fsel_ps_prob.fsel
        s_hvp = q.SelectedFieldComplexD(fsel_ps, 16)
        s_hvp.load_double_from_float(sfr, tag)
        ssp = q.SelectedShufflePlan(q.PointsSelection(fsel_ps), q.RngState(f"psel_ps-permute-{idx_xg_x}"))
        psel_lps_prob = q.SelectedPointsRealD(fsel_ps_prob, ssp)
        lps_hvp = q.SelectedPointsComplexD(s_hvp, ssp)
        lps_hvp *= hvp_type_charge_factor
        num = len(psel_lps_prob.psel)
        tot_num = q.glb_sum(num)
        hvp_list.append((psel_lps_prob, lps_hvp,))
        q.displayln_info(f"{fname}: idx={idx} ; idx_xg_x={idx_xg_x} ; xg={xg} ; tot_num={tot_num} ; num={num} .")
    sfr.close()
    assert len(hvp_list) == len(idx_xg_list_chunk)
    #
    info_str = f"{fname}: {job_tag} {traj} {inv_type_name} {inv_type_e_name} id_chunk/num_chunk={id_chunk}/{num_chunk}"
    #
    force_load_muon_line_interpolation()
    #
    points_data = []
    for idx, idx_xg_x in enumerate(idx_xg_list_chunk):
        q.displayln_info(0,
                f"{info_str} idx/chunk_size={idx}/{len(idx_xg_list_chunk)}")
        xg_x = q.Coordinate(xg_arr[idx_xg_x])
        psel_lps_prob, lps_hvp = hvp_list[idx]
        n_points_selected, n_points_computed, lslt = q.contract_two_plus_two_pair_no_glb_sum(
                complex(1.0),
                psel_prob,
                psel_lps_prob,
                idx_xg_x,
                lps_hvp,
                edl,
                r_sq_limit,
                muon_mass,
                zz_vv,
                )
        dict_val = dict()
        dict_val["idx_xg_x"] = idx_xg_x
        dict_val["xg_x"] = xg_x
        dict_val["lslt"] = lslt
        dict_val["n_points_selected"] = n_points_selected
        dict_val["n_points_computed"] = n_points_computed
        points_data.append(dict_val)
    @q.timer_verbose
    def sync_node_after_hlbl_two_plus_two_contract():
        q.sync_node()
    sync_node_after_hlbl_two_plus_two_contract()
    for d in points_data:
        d["lslt"] = q.glb_sum(d["lslt"])
        d["n_points_selected"] = q.glb_sum(d["n_points_selected"])
        d["n_points_computed"] = q.glb_sum(d["n_points_computed"])
        json_results.append((
            f"{fname}: {info_str} idx_xg_x={d['idx_xg_x']} {d['xg_x']} lslt",
            q.get_data_sig(d["lslt"], q.RngState()),
            20e-2,
            ))
    q.save_pickle_obj(points_data, get_save_path(fn))
    if len(points_data) > 0:
        labels = q.contract_two_plus_two_pair_labels()
        q.displayln_info(-1,
                f"{info_str}\n",
                show_lslt(labels, sum([ d["lslt"] for d in points_data ]) / len(points_data) * len(psel)))
        n_points_selected = sum([ d["n_points_selected"] for d in points_data ]) / len(points_data)
        n_points_computed = sum([ d["n_points_computed"] for d in points_data ]) / len(points_data)
        q.displayln_info(-1,
                f"{info_str}\n avg n_points_selected={n_points_selected} avg n_points_computed={n_points_computed}")
    q.displayln_info(0, f"{info_str} done.")
    q.release_lock()
    q.timer_display()
    q.timer_merge()

@q.timer_verbose
def run_hlbl_two_plus_two(
        job_tag,
        traj,
        *,
        inv_type,
        inv_type_e,
        get_psel_prob,
        get_edl_light,
        get_edl_strange,
        get_glb_hvp_avg_for_sub_light,
        get_glb_hvp_avg_for_sub_strange,
        get_f_rand_01,
        ):
    """
    inv_type: inv_type of the internal quark loop
    inv_type_e: inv_type of the external quark loop (with external photon)
    """
    fname = q.get_fname()
    inv_type_name_list = [ "light", "strange", ]
    inv_type_name = inv_type_name_list[inv_type]
    inv_type_e_name = inv_type_name_list[inv_type_e]
    fn_s = f"{job_tag}/hlbl/dlbl-{inv_type_name}-{inv_type_e_name}/traj-{traj}/results-brief.pickle"
    if get_load_path(fn_s) is not None:
        return
    if get_psel_prob is None:
        return
    if get_glb_hvp_avg_for_sub_light is None:
        return
    if get_glb_hvp_avg_for_sub_strange is None:
        return
    if get_edl_light is None:
        return
    if get_edl_strange is None:
        return
    if get_f_rand_01 is None:
        return
    sub_hvp_fn = f"{job_tag}/hlbl/sub-hvp-{inv_type_name}/traj-{traj}/geon-info.txt"
    if get_load_path(sub_hvp_fn) is None:
        return
    num_chunk = get_param(job_tag, "hlbl_two_plus_two_num_chunk")
    for id_chunk in range(num_chunk):
        run_hlbl_two_plus_two_chunk(
                job_tag,
                traj,
                inv_type=inv_type,
                inv_type_e=inv_type_e,
                get_psel_prob=get_psel_prob,
                get_edl_light=get_edl_light,
                get_edl_strange=get_edl_strange,
                get_glb_hvp_avg_for_sub_light=get_glb_hvp_avg_for_sub_light,
                get_glb_hvp_avg_for_sub_strange=get_glb_hvp_avg_for_sub_strange,
                get_f_rand_01=get_f_rand_01,
                id_chunk=id_chunk,
                num_chunk=num_chunk,
                )
    q.check_stop()
    q.check_time_limit()
    info_str = f"{fname}: {job_tag} {traj} {inv_type_name} {inv_type_e_name}"
    fn_chunk_list = []
    for id_chunk in range(num_chunk):
        fn_chunk = f"{job_tag}/hlbl/dlbl-points-data-{inv_type_name}-{inv_type_e_name}/traj-{traj}/chunk_{id_chunk}_{num_chunk}.pickle"
        if get_load_path(fn_chunk) is None:
            q.displayln_info(-1, f"{info_str} not finished (missing '{fn_chunk})'")
            return
        else:
            fn_chunk_list.append(fn_chunk)
    assert len(fn_chunk_list) == num_chunk
    if get_load_path(fn_s) is not None:
        return
    if not q.obtain_lock(f"locks/{job_tag}-{traj}-{fname}-{inv_type_name}-{inv_type_e_name}"):
        return
    q.timer_fork()
    labels = q.contract_two_plus_two_pair_labels()
    points_data = []
    for fn_chunk in fn_chunk_list:
        path_chunk = get_load_path(fn_chunk)
        if q.get_id_node() == 0:
            points_data += q.load_pickle_obj(path_chunk, is_sync_node=False)
    psel = get_psel_prob().psel
    if q.get_id_node() == 0:
        assert len(points_data) == len(psel)
        results = dict()
        results["labels"] = labels
        results["lslt_sum"] = sum([ d["lslt"] for d in points_data ])
        results["n_points_selected"] = sum([ d["n_points_selected"] for d in points_data ]) / len(points_data)
        results["n_points_computed"] = sum([ d["n_points_computed"] for d in points_data ]) / len(points_data)
        results["n_points"] = len(points_data)
        q.save_pickle_obj(results, get_save_path(fn_s))
        if results["n_points"] > 0:
            q.displayln_info(0, f"{info_str}\n",
                    show_lslt(labels, results["lslt_sum"]))
            n_points_selected = results["n_points_selected"]
            n_points_computed = results["n_points_computed"]
            q.displayln_info(-1,
                    f"{info_str}\n avg n_points_selected={n_points_selected} avg n_points_computed={n_points_computed}")
        json_results.append((
            f"{info_str} lslt_sum",
            q.get_data_sig(results["lslt_sum"], q.RngState()),
            3e-2,
            ))
        json_results.append((
            f"{info_str} lslt_sum[labels.index('sub proj-all'), -1, -1]",
            results["lslt_sum"][labels.index('sub proj-all'), -1, -1],
            3e-2,
            ))
    q.sync_node()
    q.displayln_info(0, f"{info_str} done.")
    q.release_lock()
    q.timer_display()
    q.timer_merge()
    return [ f"{fname}: {job_tag} {traj} {inv_type_name} {inv_type_e_name} done.", ]

# ----

@q.timer_verbose
def run_job(job_tag, traj):
    fname = q.get_fname()
    #
    traj_gf = traj
    if job_tag[:5] == "test-":
        # ADJUST ME
        traj_gf = 1000
        #
    #
    fns_produce = [
            # f"{job_tag}/auto-contract/traj-{traj}/checkpoint.txt",
            #
            # (f"{job_tag}/prop-rand-u1-light/traj-{traj}.qar", f"{job_tag}/prop-rand-u1-light/traj-{traj}/geon-info.txt",),
            # (f"{job_tag}/prop-rand-u1-strange/traj-{traj}.qar", f"{job_tag}/prop-rand-u1-strange/traj-{traj}/geon-info.txt",),
            # (f"{job_tag}/prop-rand-u1-charm/traj-{traj}.qar", f"{job_tag}/prop-rand-u1-charm/traj-{traj}/geon-info.txt",),
            #
            (f"{job_tag}/prop-psrc-light/traj-{traj}.qar", f"{job_tag}/prop-psrc-light/traj-{traj}/geon-info.txt",),
            (f"{job_tag}/prop-psrc-strange/traj-{traj}.qar", f"{job_tag}/prop-psrc-strange/traj-{traj}/geon-info.txt",),
            #
            (f"{job_tag}/prop-wsrc-light/traj-{traj}.qar", f"{job_tag}/prop-wsrc-light/traj-{traj}/geon-info.txt",),
            (f"{job_tag}/prop-wsrc-strange/traj-{traj}.qar", f"{job_tag}/prop-wsrc-strange/traj-{traj}/geon-info.txt",),
            #
            # (f"{job_tag}/prop-smear-light/traj-{traj}.qar", f"{job_tag}/prop-smear-light/traj-{traj}/geon-info.txt",),
            # (f"{job_tag}/prop-smear-strange/traj-{traj}.qar", f"{job_tag}/prop-smear-strange/traj-{traj}/geon-info.txt",),
            #
            (f"{job_tag}/psel-prop-psrc-light/traj-{traj}.qar", f"{job_tag}/psel-prop-psrc-light/traj-{traj}/checkpoint.txt",),
            (f"{job_tag}/psel-prop-psrc-strange/traj-{traj}.qar", f"{job_tag}/psel-prop-psrc-strange/traj-{traj}/checkpoint.txt",),
            #
            (f"{job_tag}/psel-prop-wsrc-light/traj-{traj}.qar", f"{job_tag}/psel-prop-wsrc-light/traj-{traj}/checkpoint.txt",),
            (f"{job_tag}/psel-prop-wsrc-strange/traj-{traj}.qar", f"{job_tag}/psel-prop-wsrc-strange/traj-{traj}/checkpoint.txt",),
            #
            # (f"{job_tag}/psel-prop-smear-light/traj-{traj}.qar", f"{job_tag}/psel-prop-smear-light/traj-{traj}/checkpoint.txt",),
            # (f"{job_tag}/psel-prop-smear-strange/traj-{traj}.qar", f"{job_tag}/psel-prop-smear-strange/traj-{traj}/checkpoint.txt",),
            #
            f"{job_tag}/hvp-average/traj-{traj}/hvp_average_light.field",
            f"{job_tag}/hvp-average/traj-{traj}/hvp_average_strange.field",
            ]
    if job_tag[:5] == "test-":
        has_eig = True
        fns_need = []
    else:
        has_eig = get_load_path(f"{job_tag}/eig/traj-{traj}/metadata.txt") is not None
        fns_need = [
                (f"{job_tag}/configs/ckpoint_lat.{traj_gf}", f"{job_tag}/configs/ckpoint_lat.IEEE64BIG.{traj_gf}",),
                # f"{job_tag}/eig/traj-{traj_gf}/metadata.txt",
                # f"{job_tag}/gauge-transform/traj-{traj_gf}.field",
                # f"{job_tag}/point-selection/traj-{traj}.txt",
                # f"{job_tag}/field-selection/traj-{traj}.field",
                # f"{job_tag}/wall-src-info-light/traj-{traj}.txt",
                # f"{job_tag}/wall-src-info-strange/traj-{traj}.txt",
                ]
    #
    if not check_job(job_tag, traj, fns_produce, fns_need):
        return
    #
    run_ret_list = []
    def add_to_run_ret_list(v):
        nonlocal run_ret_list
        if v is None:
            return
        if isinstance(v, list):
            run_ret_list += v
        else:
            run_ret_list.append(v)
    #
    get_gf = run_gf(job_tag, traj_gf)
    get_gt = run_gt(job_tag, traj_gf, get_gf)
    get_gf_ape = run_gf_ape(job_tag, get_gf)
    #
    get_wi = run_wi(job_tag, traj)
    #
    if has_eig:
        get_eig_light = run_eig(job_tag, traj_gf, get_gf)
    get_eig_strange = run_eig_strange(job_tag, traj_gf, get_gf)
    #
    def run_wsrc_full():
        if has_eig:
            get_eig = get_eig_light
            # run_get_inverter(job_tag, traj, inv_type=0, get_gf=get_gf, get_gt=get_gt, get_eig=get_eig)
            v = run_prop_wsrc_full(job_tag, traj, inv_type=0, get_gf=get_gf, get_eig=get_eig, get_gt=get_gt, get_wi=get_wi)
            add_to_run_ret_list(v)
            q.clean_cache(q.cache_inv)
        #
        get_eig = get_eig_strange
        # run_get_inverter(job_tag, traj, inv_type=1, get_gf=get_gf, get_gt=get_gt, get_eig=get_eig)
        v = run_prop_wsrc_full(job_tag, traj, inv_type=1, get_gf=get_gf, get_eig=get_eig, get_gt=get_gt, get_wi=get_wi)
        add_to_run_ret_list(v)
        q.clean_cache(q.cache_inv)
    #
    run_wsrc_full()
    #
    get_f_weight = run_f_weight_from_wsrc_prop_full(job_tag, traj)
    get_f_rand_01 = run_f_rand_01(job_tag, traj)
    get_fsel_prob = run_fsel_prob(job_tag, traj, get_f_rand_01=get_f_rand_01, get_f_weight=get_f_weight)
    get_psel_prob = run_psel_prob(job_tag, traj, get_f_rand_01=get_f_rand_01, get_f_weight=get_f_weight)
    get_fsel = run_fsel_from_fsel_prob(get_fsel_prob)
    get_psel = run_psel_from_psel_prob(get_psel_prob)
    #
    if get_fsel is None:
        q.clean_cache()
        return
    #
    get_fselc = run_fselc(job_tag, traj, get_fsel, get_psel)
    #
    v = run_prop_wsrc_sparse(job_tag, traj, inv_type=0, get_gt=get_gt, get_psel=get_psel, get_fsel=get_fsel, get_wi=get_wi)
    add_to_run_ret_list(v)
    v = run_prop_wsrc_sparse(job_tag, traj, inv_type=1, get_gt=get_gt, get_psel=get_psel, get_fsel=get_fsel, get_wi=get_wi)
    add_to_run_ret_list(v)
    #
    get_psel_smear = run_psel_smear(job_tag, traj)
    #
    def run_with_eig():
        if has_eig:
            get_eig = get_eig_light
            # run_get_inverter(job_tag, traj, inv_type=0, get_gf=get_gf, get_eig=get_eig)
            # v = run_prop_rand_u1(job_tag, traj, inv_type=0, get_gf=get_gf, get_fsel=get_fsel, get_eig=get_eig)
            # add_to_run_ret_list(v)
            v = run_prop_psrc(job_tag, traj, inv_type=0, get_gf=get_gf, get_eig=get_eig, get_gt=get_gt, get_psel=get_psel, get_fsel=get_fselc, get_f_rand_01=get_f_rand_01)
            add_to_run_ret_list(v)
            # v = run_prop_smear(job_tag, traj, inv_type=0, get_gf=get_gf, get_gf_ape=get_gf_ape, get_eig=get_eig, get_gt=get_gt, get_psel=get_psel, get_fsel=get_fselc, get_psel_smear=get_psel_smear)
            # add_to_run_ret_list(v)
            q.clean_cache(q.cache_inv)
    #
    def run_with_eig_strange():
        get_eig = get_eig_strange
        # run_get_inverter(job_tag, traj, inv_type=1, get_gf=get_gf, get_eig=get_eig)
        # v = run_prop_rand_u1(job_tag, traj, inv_type=1, get_gf=get_gf, get_fsel=get_fsel, get_eig=get_eig)
        # add_to_run_ret_list(v)
        v = run_prop_psrc(job_tag, traj, inv_type=1, get_gf=get_gf, get_eig=get_eig, get_gt=get_gt, get_psel=get_psel, get_fsel=get_fselc, get_f_rand_01=get_f_rand_01)
        add_to_run_ret_list(v)
        # v = run_prop_smear(job_tag, traj, inv_type=1, get_gf=get_gf, get_gf_ape=get_gf_ape, get_eig=get_eig, get_gt=get_gt, get_psel=get_psel, get_fsel=get_fselc, get_psel_smear=get_psel_smear)
        # add_to_run_ret_list(v)
        q.clean_cache(q.cache_inv)
    #
    def run_charm():
        # run_get_inverter(job_tag, traj, inv_type=2, get_gf=get_gf)
        # v = run_prop_rand_u1(job_tag, traj, inv_type=2, get_gf=get_gf, get_fsel=get_fsel)
        # add_to_run_ret_list(v)
        q.clean_cache(q.cache_inv)
    #
    run_with_eig()
    run_with_eig_strange()
    run_charm()
    #
    get_hvp_average_light = run_hvp_average(job_tag, traj, inv_type=0, get_psel_prob=get_psel_prob)
    get_hvp_average_strange = run_hvp_average(job_tag, traj, inv_type=1, get_psel_prob=get_psel_prob)
    #
    q.clean_cache()
    #
    q.sync_node()
    q.displayln_info(f"{fname}: run_ret_list={run_ret_list}")
    if job_tag[:5] != "test-":
        if run_ret_list:
            q.qquit(f"{fname} {job_tag} {traj} (partly) done.")

@q.timer_verbose
def run_job_contract(job_tag, traj):
    fname = q.get_fname()
    #
    traj_gf = traj
    if job_tag[:5] == "test-":
        # ADJUST ME
        traj_gf = 1000
        #
    #
    fn_checkpoint_auto_contract = f"{job_tag}/auto-contract/traj-{traj}/checkpoint.txt"
    #
    fns_produce = [
            fn_checkpoint_auto_contract,
            f"{job_tag}/hlbl/dlbl-light-light/traj-{traj}/results-brief.pickle",
            f"{job_tag}/hlbl/dlbl-light-strange/traj-{traj}/results-brief.pickle",
            f"{job_tag}/hlbl/dlbl-strange-light/traj-{traj}/results-brief.pickle",
            f"{job_tag}/hlbl/dlbl-strange-strange/traj-{traj}/results-brief.pickle",
            f"{job_tag}/hlbl/clbl-light/traj-{traj}/results-brief.pickle",
            f"{job_tag}/hlbl/clbl-strange/traj-{traj}/results-brief.pickle",
            ]
    fns_need = [
            #
            (f"{job_tag}/configs/ckpoint_lat.{traj_gf}", f"{job_tag}/configs/ckpoint_lat.IEEE64BIG.{traj_gf}",),
            f"{job_tag}/gauge-transform/traj-{traj_gf}.field",
            f"{job_tag}/wall-src-info-light/traj-{traj}.txt",
            f"{job_tag}/wall-src-info-strange/traj-{traj}.txt",
            f"{job_tag}/point-selection/traj-{traj}.txt",
            f"{job_tag}/field-selection/traj-{traj}.field",
            f"{job_tag}/field-selection-weight/traj-{traj}/f-rand-01.field",
            f"{job_tag}/field-selection-weight/traj-{traj}/weight.field",
            f"{job_tag}/field-selection-weight/traj-{traj}/fsel-prob.sfield",
            f"{job_tag}/field-selection-weight/traj-{traj}/psel-prob.lat",
            #
            f"{job_tag}/hvp-average/hvp_average_light.field",
            f"{job_tag}/hvp-average/hvp_average_light.trajs.txt",
            f"{job_tag}/hvp-average/hvp_average_strange.field",
            f"{job_tag}/hvp-average/hvp_average_strange.trajs.txt",
            #
            f"{job_tag}/hvp-average/traj-{traj}/hvp_average_light.field",
            f"{job_tag}/hvp-average/traj-{traj}/hvp_average_strange.field",
            #
            # (f"{job_tag}/prop-rand-u1-light/traj-{traj}.qar", f"{job_tag}/prop-rand-u1-light/traj-{traj}/geon-info.txt",),
            # (f"{job_tag}/prop-rand-u1-strange/traj-{traj}.qar", f"{job_tag}/prop-rand-u1-strange/traj-{traj}/geon-info.txt",),
            # (f"{job_tag}/prop-rand-u1-charm/traj-{traj}.qar", f"{job_tag}/prop-rand-u1-charm/traj-{traj}/geon-info.txt",),
            #
            (f"{job_tag}/prop-psrc-light/traj-{traj}.qar", f"{job_tag}/prop-psrc-light/traj-{traj}/geon-info.txt",),
            (f"{job_tag}/prop-psrc-strange/traj-{traj}.qar", f"{job_tag}/prop-psrc-strange/traj-{traj}/geon-info.txt",),
            #
            (f"{job_tag}/prop-wsrc-light/traj-{traj}.qar", f"{job_tag}/prop-wsrc-light/traj-{traj}/geon-info.txt",),
            (f"{job_tag}/prop-wsrc-strange/traj-{traj}.qar", f"{job_tag}/prop-wsrc-strange/traj-{traj}/geon-info.txt",),
            #
            # (f"{job_tag}/prop-smear-light/traj-{traj}.qar", f"{job_tag}/prop-smear-light/traj-{traj}/geon-info.txt",),
            # (f"{job_tag}/prop-smear-strange/traj-{traj}.qar", f"{job_tag}/prop-smear-strange/traj-{traj}/geon-info.txt",),
            #
            (f"{job_tag}/psel-prop-psrc-light/traj-{traj}.qar.idx", f"{job_tag}/psel-prop-psrc-light/traj-{traj}/checkpoint.txt",),
            (f"{job_tag}/psel-prop-psrc-strange/traj-{traj}.qar.idx", f"{job_tag}/psel-prop-psrc-strange/traj-{traj}/checkpoint.txt",),
            #
            (f"{job_tag}/psel-prop-wsrc-light/traj-{traj}.qar.idx", f"{job_tag}/psel-prop-wsrc-light/traj-{traj}/checkpoint.txt",),
            (f"{job_tag}/psel-prop-wsrc-strange/traj-{traj}.qar.idx", f"{job_tag}/psel-prop-wsrc-strange/traj-{traj}/checkpoint.txt",),
            #
            # (f"{job_tag}/psel-prop-smear-light/traj-{traj}.qar.idx", f"{job_tag}/psel-prop-smear-light/traj-{traj}/checkpoint.txt",),
            # (f"{job_tag}/psel-prop-smear-strange/traj-{traj}.qar.idx", f"{job_tag}/psel-prop-smear-strange/traj-{traj}/checkpoint.txt",),
            ]
    #
    if not check_job(job_tag, traj, fns_produce, fns_need):
        return
    #
    run_ret_list = []
    def add_to_run_ret_list(v):
        nonlocal run_ret_list
        if v is None:
            return
        if isinstance(v, list):
            run_ret_list += v
        else:
            run_ret_list.append(v)
    #
    run_r_list(job_tag)
    #
    get_gf = run_gf(job_tag, traj_gf)
    get_gt = run_gt(job_tag, traj_gf, get_gf)
    get_gf_ape = run_gf_ape(job_tag, get_gf)
    #
    get_wi = run_wi(job_tag, traj)
    #
    get_f_weight = run_f_weight_from_wsrc_prop_full(job_tag, traj)
    get_f_rand_01 = run_f_rand_01(job_tag, traj)
    get_psel_prob = run_psel_prob(job_tag, traj, get_f_rand_01=get_f_rand_01, get_f_weight=get_f_weight)
    get_fsel_prob = run_fsel_prob(job_tag, traj, get_f_rand_01=get_f_rand_01, get_f_weight=get_f_weight)
    get_fsel = run_fsel_from_fsel_prob(get_fsel_prob)
    get_psel = run_psel_from_psel_prob(get_psel_prob)
    #
    get_psel_smear = run_psel_smear(job_tag, traj)
    #
    for inv_type in [ 0, 1, ]:
        get_glb_hvp_avg = run_job_global_hvp_average(job_tag, inv_type=inv_type)
        get_hvp_average = run_hvp_average(job_tag, traj, inv_type=inv_type, get_psel_prob=get_psel_prob)
        get_glb_hvp_avg_for_sub = run_job_global_hvp_average_for_subtract(job_tag, traj, inv_type=inv_type, get_glb_hvp_avg=get_glb_hvp_avg, get_hvp_average=get_hvp_average)
        get_hvp_sum_tslice_accs = run_hvp_sum_tslice_accs(job_tag, traj, inv_type=inv_type, get_psel=get_psel)
        get_hvp_sum_tslice = run_hvp_sum_tslice(job_tag, traj, inv_type=inv_type, get_psel=get_psel, get_hvp_sum_tslice_accs=get_hvp_sum_tslice_accs)
        get_edl = run_edl(job_tag, traj, inv_type=inv_type, get_psel=get_psel, get_hvp_sum_tslice=get_hvp_sum_tslice)
        run_check_hvp_avg(job_tag, traj, inv_type=inv_type, get_psel_prob=get_psel_prob, get_hvp_sum_tslice=get_hvp_sum_tslice, get_hvp_average=get_hvp_average)
        run_hlbl_sub_hvp_sfield(job_tag, traj, inv_type=inv_type, get_psel_prob=get_psel_prob, get_glb_hvp_avg_for_sub=get_glb_hvp_avg_for_sub, get_f_rand_01=get_f_rand_01)
        if inv_type == 0:
            get_edl_light = get_edl
            get_hvp_average_light = get_hvp_average
            get_glb_hvp_avg_for_sub_light = get_glb_hvp_avg_for_sub
        elif inv_type == 1:
            get_edl_strange = get_edl
            get_hvp_average_strange = get_hvp_average
            get_glb_hvp_avg_for_sub_strange = get_glb_hvp_avg_for_sub
        else:
            raise Exception(f"{fname}: inv_type={inv_type} wrong.")
    #
    if is_performing_hlbl_contraction:
        #
        if job_tag[:5] == "test-":
            force_load_muon_line_interpolation()
        #
        for inv_type in [ 0, 1, ]:
            get_point_pairs = run_hlbl_four_point_pairs_info(
                    job_tag,
                    traj,
                    inv_type=inv_type,
                    get_psel_prob=get_psel_prob,
                    )
            v = run_hlbl_four(
                    job_tag,
                    traj,
                    inv_type=inv_type,
                    get_psel_prob=get_psel_prob,
                    get_fsel_prob=get_fsel_prob,
                    get_point_pairs=get_point_pairs,
                    )
            add_to_run_ret_list(v)
        #
        for inv_type in [ 0, 1, ]:
            for inv_type_e in [ 0, 1, ]:
                v = run_hlbl_two_plus_two(
                        job_tag,
                        traj,
                        inv_type=inv_type,
                        inv_type_e=inv_type_e,
                        get_psel_prob=get_psel_prob,
                        get_edl_light=get_edl_light,
                        get_edl_strange=get_edl_strange,
                        get_glb_hvp_avg_for_sub_light=get_glb_hvp_avg_for_sub_light,
                        get_glb_hvp_avg_for_sub_strange=get_glb_hvp_avg_for_sub_strange,
                        get_f_rand_01=get_f_rand_01,
                        )
                add_to_run_ret_list(v)
    #
    get_get_prop = run_get_prop(job_tag, traj,
            get_gf=get_gf,
            get_gt=get_gt,
            get_psel=get_psel,
            get_fsel=get_fsel,
            get_psel_smear=get_psel_smear,
            prop_types=[
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
    if is_performing_auto_contraction and q.obtain_lock(f"locks/{job_tag}-{traj}-{fname}"):
        q.timer_fork()
        get_prop = get_get_prop()
        if get_prop is not None:
            # ADJUST ME
            auto_contract_meson_corr(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob)
            auto_contract_meson_corr_psnk(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob)
            auto_contract_meson_corr_psrc(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob)
            auto_contract_meson_corr_psnk_psrc(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob)
            #
            json_results.append((f"get_hvp_average_light: {traj}", q.get_data_sig(get_hvp_average_light(), q.RngState()),))
            json_results.append((f"get_hvp_average_strange: {traj}:", q.get_data_sig(get_hvp_average_strange(), q.RngState()),))
            json_results.append((f"get_glb_hvp_avg_for_sub_light: {traj}:", q.get_data_sig(get_glb_hvp_avg_for_sub_light(), q.RngState()),))
            json_results.append((f"get_glb_hvp_avg_for_sub_strange: {traj}:", q.get_data_sig(get_glb_hvp_avg_for_sub_strange(), q.RngState()),))
            #
            q.qtouch_info(get_save_path(fn_checkpoint_auto_contract))
            q.displayln_info("timer_display for runjob")
        q.release_lock()
        q.timer_display()
        q.timer_merge()
        v = [ f"{fname} {job_tag} {traj} done", ]
        add_to_run_ret_list(v)
    #
    q.clean_cache()
    #
    q.sync_node()
    q.displayln_info(f"{fname}: run_ret_list={run_ret_list}")
    if job_tag[:5] != "test-":
        if run_ret_list:
            q.qquit(f"{fname} {job_tag} {traj} (partly) done.")

# ----

set_param("test-4nt8", "mk_sample_gauge_field", "rand_n_step")(2)
set_param("test-4nt8", "mk_sample_gauge_field", "flow_n_step")(8)
set_param("test-4nt8", "mk_sample_gauge_field", "hmc_n_traj")(1)
set_param("test-4nt8", "fermion_params", 0, 2, "Ls")(8)
set_param("test-4nt8", "fermion_params", 1, 2, "Ls")(8)
set_param("test-4nt8", "fermion_params", 2, 2, "Ls")(8)
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
set_param("test-4nt8", "cg_params-0-2", "pv_maxiter", value=5)
set_param("test-4nt8", "cg_params-1-2", "pv_maxiter", value=5)
set_param("test-4nt8", "a_inv_gev")(1.73)
set_param("test-4nt8", "zz_vv")(0.71)

set_param("24D", "lanc_params", 1, value=None)
set_param("24D", "clanc_params", 1, value=None)
set_param("24D", 'fermion_params', 0, 2, value=deepcopy(get_param("24D", 'fermion_params', 0, 0)))
set_param("24D", 'fermion_params', 1, 0, value=deepcopy(get_param("24D", 'fermion_params', 1, 2)))
set_param("24D", 'fermion_params', 1, 1, value=deepcopy(get_param("24D", 'fermion_params', 1, 2)))

set_param("test-4nt8", "trajs", value=[ 1000, 2000, ])
set_param("test-8nt16", "trajs", value=[ 1000, 2000, ])
set_param("24D", "trajs", value=list(range(2000, 3000, 10)))
set_param("48I", "trajs", value=list(range(975, 2185, 10)) + list(range(1102, 1502, 10)))
set_param("64I", "trajs", value=list(range(1200, 3680, 40)))

set_param("test-4nt8", "hlbl_four_prob_scaling_factor", value=1.0)
set_param("test-8nt16", "hlbl_four_prob_scaling_factor", value=1.0)
set_param("24D", "hlbl_four_prob_scaling_factor", value=1.0)
set_param("48I", "hlbl_four_prob_scaling_factor", value=1.0)
set_param("64I", "hlbl_four_prob_scaling_factor", value=1.0)

set_param("test-4nt8", "hlbl_four_prob_scaling_factor_strange", value=1.0)
set_param("test-8nt16", "hlbl_four_prob_scaling_factor_strange", value=1.0)
set_param("24D", "hlbl_four_prob_scaling_factor_strange", value=1.0)
set_param("48I", "hlbl_four_prob_scaling_factor_strange", value=1.0)
set_param("64I", "hlbl_four_prob_scaling_factor_strange", value=1.0)

set_param("test-4nt8", "hlbl_four_num_chunk", value=3)
set_param("test-8nt16", "hlbl_four_num_chunk", value=6)
set_param("24D", "hlbl_four_num_chunk", value=512)
set_param("48I", "hlbl_four_num_chunk", value=1024)
set_param("64I", "hlbl_four_num_chunk", value=2048)

# larger value means less computation
set_param("test-4nt8", "hlbl_four_contract_sparse_ratio", value=2.0)
set_param("test-8nt16", "hlbl_four_contract_sparse_ratio", value=2.0)
set_param("24D", "hlbl_four_contract_sparse_ratio", value=10.0)
set_param("48I", "hlbl_four_contract_sparse_ratio", value=20.0)
set_param("64I", "hlbl_four_contract_sparse_ratio", value=20.0)

set_param("test-4nt8", "hlbl_two_plus_two_num_hvp_sel_threshold", value=5e-3)
set_param("test-8nt16", "hlbl_two_plus_two_num_hvp_sel_threshold", value=5e-3)
set_param("24D", "hlbl_two_plus_two_num_hvp_sel_threshold", value=5e-4)
set_param("48I", "hlbl_two_plus_two_num_hvp_sel_threshold", value=5e-5)
set_param("64I", "hlbl_two_plus_two_num_hvp_sel_threshold", value=5e-5)

set_param("test-4nt8", "hlbl_two_plus_two_num_chunk", value=3)
set_param("test-8nt16", "hlbl_two_plus_two_num_chunk", value=6)
set_param("24D", "hlbl_two_plus_two_num_chunk", value=4)
set_param("48I", "hlbl_two_plus_two_num_chunk", value=8)
set_param("64I", "hlbl_two_plus_two_num_chunk", value=8)

# ----

if __name__ == "__main__":

    qg.begin_with_gpt()

    ##################### CMD options #####################

    job_tags = q.get_arg("--job_tags", default="").split(",")

    is_performing_inversion = q.get_arg("--no-inversion", default=None) is None

    is_performing_contraction = q.get_arg("--no-contract", default=None) is None

    is_performing_hlbl_contraction = q.get_arg("--no-hlbl-contract", default=None) is None

    is_performing_auto_contraction = q.get_arg("--no-auto-contract", default=None) is None

    #######################################################

    job_tags_default = [
            "test-4nt8",
            # "test-8nt16",
            # "24D",
            # "64I",
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
                run_job(job_tag, traj)
        if is_performing_contraction:
            for inv_type in [ 0, 1, ]:
                q.check_time_limit()
                run_job_global_hvp_average(job_tag, inv_type=inv_type)
        for traj in get_param(job_tag, "trajs"):
            if is_performing_contraction:
                q.check_time_limit()
                run_job_contract(job_tag, traj)

    q.check_log_json(__file__, json_results)

    q.timer_display()

    qg.end_with_gpt()

    q.displayln_info("CHECK: finished successfully.")

# ----
