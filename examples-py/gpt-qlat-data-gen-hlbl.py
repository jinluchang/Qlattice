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
    fsel_n_elems = fsel.n_elems()
    fsel_prob_arr = fsel_prob[:].ravel()
    psel_prob_arr = psel_prob[:].ravel()
    xg_fsel_arr = fsel.to_psel_local()[:]
    geo = q.Geometry(total_site, 1)
    total_volume = geo.total_volume()
    def load_data():
        t_t_list = get_mpi_chunk(
                [ (t_src, t_snk,) for t_snk in range(total_site[3]) for t_src in range(total_site[3]) ],
                rng_state = None)
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
    fsel_n_elems = fsel.n_elems()
    fsel_prob_arr = fsel_prob[:].ravel()
    psel_prob_arr = psel_prob[:].ravel()
    xg_fsel_arr = fsel.to_psel_local()[:]
    geo = q.Geometry(total_site, 1)
    total_volume = geo.total_volume()
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
    fsel_n_elems = fsel.n_elems()
    fsel_prob_arr = fsel_prob[:].ravel()
    psel_prob_arr = psel_prob[:].ravel()
    xg_fsel_arr = fsel.to_psel_local()[:]
    xg_psel_arr = psel[:]
    geo = q.Geometry(total_site, 1)
    total_volume = geo.total_volume()
    def load_data():
        x_t_list = get_mpi_chunk(
                [ (pidx, t_snk,) for t_snk in range(total_site[3]) for pidx in range(len(xg_psel_arr)) ],
                rng_state = None)
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
    fsel_n_elems = fsel.n_elems()
    fsel_prob_arr = fsel_prob[:].ravel()
    psel_prob_arr = psel_prob[:].ravel()
    xg_fsel_arr = fsel.to_psel_local()[:]
    xg_psel_arr = psel[:]
    geo = q.Geometry(total_site, 1)
    total_volume = geo.total_volume()
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
    geo = q.Geometry(total_site, 1)
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
            get_f_weight = run_f_weight_from_wsrc_prop_full(job_tag, traj, get_wi=get_wi)
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
    return ret

@q.timer_verbose
def run_job_global_hvp_average_for_subtract(job_tag, traj, *, get_glb_hvp_avg, get_hvp_average):
    """
    get_glb_hvp_avg_for_sub = run_job_global_hvp_average_for_subtract(job_tag, traj, get_glb_hvp_avg=get_glb_hvp_avg, get_hvp_average=get_hvp_average)
    glb_hvp_avg_for_sub = get_glb_hvp_avg_for_sub()
    #
    Get global hvp average excluding data from this traj.
    Suitable for use in subtraction.
    """
    fname = q.get_fname()
    if get_glb_hvp_avg is None:
        q.displayln_info(-1, f"{fname}: get_glb_hvp_avg is None.")
        return None
    if get_hvp_average is None:
        q.displayln_info(-1, f"{fname}: get_hvp_average is None.")
        return None
    @q.lazy_call
    @q.timer_verbose
    def get_glb_hvp_avg_for_sub():
        glb_hvp_avg_trajs, glb_hvp_avg = get_glb_hvp_avg()
        glb_hvp_avg_for_sub = glb_hvp_avg.copy()
        if traj not in glb_hvp_avg_trajs:
            return glb_hvp_avg_for_sub
        num_trajs = len(glb_hvp_avg_trajs)
        assert num_trajs > 0
        if num_trajs == 1:
            q.displayln_info(-1, f"WARNING: {fname} glb_hvp_avg num_trajs={num_trajs}")
            return glb_hvp_avg_for_sub
        hvp_average = get_hvp_average()
        glb_hvp_avg_for_sub *= num_trajs
        glb_hvp_avg_for_sub -= hvp_average
        glb_hvp_avg_for_sub *= 1.0 / (num_trajs - 1.0)
        return glb_hvp_avg_for_sub
    ret = get_glb_hvp_avg_for_sub
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
    if not sfr.has(tag):
        return None
    s_prop = q.SelProp(fsel)
    s_prop.load_double_from_float(sfr, tag)
    return s_prop

def get_list_chunk(x_list, id_chunk, num_chunk):
    x_list_size = len(x_list)
    chunk_size = (x_list_size - 1) // num_chunk + 1
    assert chunk_size >= 0
    chunk_start = chunk_size * id_chunk
    chunk_stop = chunk_start + chunk_size
    if id_chunk == num_chunk - 1:
        assert chunk_stop >= x_list_size
    x_list_chunk = x_list[chunk_start:chunk_stop]
    return x_list_chunk

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
        elif r_sq <= 12 * 12:
            prob = 1.0
        else:
            prob = (12.0 / np.sqrt(r_sq))**3
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
    geo = q.Geometry(total_site, 1)
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
    geo = q.Geometry(total_site, 1)
    total_volume = geo.total_volume()
    inv_type_name_list = [ "light", "strange", ]
    inv_type_name = inv_type_name_list[inv_type]
    #
    psel_prob = get_psel_prob()
    psel = psel_prob.psel
    xg_arr = psel.xg_arr()
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
                dict_val = {}
                dict_val["idx_xg_x"] = i
                dict_val["idx_xg_y"] = j
                dict_val["xg_x"] = xg_x
                dict_val["xg_y"] = xg_y
                dict_val["r"] = get_r_coordinate(xg_diff, total_site)
                dict_val["prob_accept"] = prob_accept
                dict_val["weight_pair"] = 1.0 / prob_accept
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
    point_pairs = mk_hlbl_four_point_pairs(
            job_tag, traj,
            inv_type=inv_type,
            get_psel_prob=get_psel_prob,
            )
    q.save_pickle_obj(point_pairs, get_save_path(fn))
    q.release_lock()
    return ret

def get_hlbl_clbl_info_ref_tags(job_tag):
    return [ "ref-far", "ref-close", "ref-center", ]

@q.timer
def contract_hlbl_four_labels(job_tag):
    tags = get_hlbl_clbl_info_ref_tags(job_tag)
    return q.contract_four_pair_labels(tags)

@q.timer_verbose
def contract_hlbl_four_ama(job_tag, *, inv_type, get_prop, idx_xg_x, idx_xg_y, psel_prob, fsel_prob, weight_pair):
    """
    get_prop(xg) => sprop_ama
    """
    psel = psel_prob.psel
    fsel = fsel_prob.fsel
    xg_x = psel.coordinate_from_idx(idx_xg_x)
    xg_y = psel.coordinate_from_idx(idx_xg_y)
    muon_mass = get_muon_mass(job_tag)
    coef = complex(weight_pair)
    force_load_muon_line_interpolation()
    smf_d = q.mk_m_z_field_tag(fsel, xg_x, xg_y, a=muon_mass, tag=0)
    tags = get_hlbl_clbl_info_ref_tags(job_tag)
    r_sq_limit = get_r_sq_limit(job_tag)
    zz_vv = get_param(job_tag, "zz_vv")
    def f(sprop_x, sprop_y):
        return q.contract_four_pair(
                coef,
                psel_prob,
                fsel_prob,
                idx_xg_x,
                idx_xg_y,
                smf_d,
                sprop_x,
                sprop_y,
                inv_type,
                tags,
                r_sq_limit,
                muon_mass,
                zz_vv,
                )
    ama_val = ama_apply2(f, get_prop(xg_x), get_prop(xg_y))
    return (ama_extract(ama_val, is_sloppy=False), ama_extract(ama_val, is_sloppy=True),)

@q.timer
def run_hlbl_four_chunk(job_tag, traj, *, inv_type, get_psel_prob, get_fsel_prob, get_point_pairs, id_chunk, num_chunk):
    fname = q.get_fname()
    inv_type_name_list = [ "light", "strange", ]
    inv_type_name = inv_type_name_list[inv_type]
    fn = f"{job_tag}/hlbl/clbl-{inv_type_name}/traj-{traj}/chunk_{id_chunk}_{num_chunk}.pickle"
    if get_load_path(fn) is not None:
        return
    q.check_stop()
    q.check_time_limit()
    if not q.obtain_lock(f"locks/{job_tag}-{traj}-{fname}-{inv_type_name}-{id_chunk}-{num_chunk}"):
        return
    #
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    geo = q.Geometry(total_site, 1)
    total_volume = geo.total_volume()
    #
    psel_prob = get_psel_prob()
    psel = psel_prob.psel
    fsel_prob = get_fsel_prob()
    fsel = fsel_prob.fsel
    #
    point_pairs = get_point_pairs()
    point_pairs_chunk = get_list_chunk(point_pairs, id_chunk, num_chunk)
    #
    rel_acc_list = [ 0, 1, 2, ]
    prob_list = [ get_param(job_tag, f"prob_acc_{inv_acc}_psrc") for inv_acc in rel_acc_list ]
    #
    sfr = q.open_fields(get_load_path(f"{job_tag}/prop-psrc-{inv_type_name}/traj-{traj}/geon-info.txt"), "r")
    #
    @functools.lru_cache(maxsize=3000)
    @q.timer
    def get_prop_cache(xg, inv_acc):
        """
        xg must be tuple
        """
        return get_psrc_prop(job_tag, traj, q.Coordinate(xg), inv_type, inv_acc, sfr=sfr, fsel=fsel)
    #
    @q.timer
    def get_prop(xg):
        val_list = [ get_prop_cache(xg.to_tuple(), inv_acc) for inv_acc in rel_acc_list ]
        return mk_ama_val(val_list[0], xg.to_tuple(), val_list, rel_acc_list, prob_list)
    #
    q.displayln_info(f"get_prop_cache info {get_prop_cache.cache_info()}")
    labels = contract_hlbl_four_labels(job_tag)
    pairs_data = []
    for idx, pp in enumerate(point_pairs_chunk):
        idx_xg_x = pp["idx_xg_x"]
        idx_xg_y = pp["idx_xg_y"]
        xg_x = pp["xg_x"]
        xg_y = pp["xg_y"]
        assert xg_x == psel.coordinate_from_idx(idx_xg_x)
        assert xg_y == psel.coordinate_from_idx(idx_xg_y)
        r = pp["r"]
        prob_accept = pp["prob_accept"]
        weight_pair = pp["weight_pair"]
        lslt, lslt_sloppy = contract_hlbl_four_ama(
                job_tag,
                inv_type=inv_type,
                get_prop=get_prop,
                idx_xg_x=idx_xg_x,
                idx_xg_y=idx_xg_y,
                psel_prob=psel_prob,
                fsel_prob=fsel_prob,
                weight_pair=weight_pair,
                )
        info_tag = f"{job_tag} {traj} {inv_type_name} {id_chunk}/{num_chunk} {idx}/{len(point_pairs_chunk)} {xg_x} {xg_y}"
        q.displayln_info(
                f"{fname}: {info_tag}\n",
                f"r={r} weight_pair={weight_pair}\n",
                show_lslt(labels, lslt * len(point_pairs), label="ref-far proj-all"),
                "\nsloppy:\n",
                show_lslt(labels, lslt_sloppy * len(point_pairs), label="ref-far proj-all"))
        json_results.append((
            f"{fname}: {info_tag} lslt",
            q.get_data_sig(lslt, q.RngState()),
            5e-3,
            ))
        json_results.append((
            f"{fname}: {info_tag} lslt_sloppy",
            q.get_data_sig(lslt_sloppy, q.RngState()),
            5e-3,
            ))
        dict_val = {}
        dict_val["lslt"] = lslt
        dict_val["lslt_sloppy"] = lslt_sloppy
        dict_val["xg_x"] = xg_x
        dict_val["xg_y"] = xg_y
        dict_val["r"] = r
        dict_val["prob_accept"] = prob_accept
        dict_val["weight_pair"] = weight_pair
        pairs_data.append(dict_val)
    sfr.close()
    if len(point_pairs_chunk) != len(pairs_data):
        raise Exception(f"len(point_pairs_chunk)={len(point_pairs_chunk)} len(pairs_data)={len(pairs_data)}")
    if len(pairs_data) > 0:
        q.displayln_info(f"{fname}: {job_tag} {traj} {inv_type_name} {id_chunk}/{num_chunk}\n",
                show_lslt(labels, sum([ d["lslt"] for d in pairs_data ]) / len(point_pairs_chunk) * len(point_pairs)))
    q.save_pickle_obj(pairs_data, get_save_path(fn))
    q.release_lock()
    q.displayln_info(f"get_prop_cache info {get_prop_cache.cache_info()}")

@q.timer
def run_hlbl_four(job_tag, traj, *, inv_type, get_psel_prob, get_fsel_prob, get_point_pairs):
    fname = q.get_fname()
    inv_type_name_list = [ "light", "strange", ]
    inv_type_name = inv_type_name_list[inv_type]
    if get_point_pairs is None:
        q.displayln_info(f"{fname}: {job_tag} {traj} {inv_type_name} get_point_pairs is None")
        return
    fn = f"{job_tag}/hlbl/clbl-{inv_type_name}/traj-{traj}/results.pickle"
    if get_load_path(fn) is not None:
        return
    num_chunk = get_param(job_tag, "hlbl_four_num_chunk")
    for id_chunk in range(num_chunk):
        run_hlbl_four_chunk(
                job_tag, traj,
                inv_type=inv_type,
                get_psel_prob=get_psel_prob,
                get_fsel_prob=get_fsel_prob,
                get_point_pairs=get_point_pairs,
                id_chunk=id_chunk,
                num_chunk=num_chunk,
                )
    fn_chunk_list = []
    for id_chunk in range(num_chunk):
        fn_chunk = f"{job_tag}/hlbl/clbl-{inv_type_name}/traj-{traj}/chunk_{id_chunk}_{num_chunk}.pickle"
        if get_load_path(fn_chunk) is None:
            q.displayln_info(f"{fname}: {job_tag} {traj} {inv_type_name} not finished (missing '{fn_chunk})'")
            return
        else:
            fn_chunk_list.append(fn_chunk)
    assert len(fn_chunk_list) == num_chunk
    if get_load_path(fn) is not None:
        return
    if not q.obtain_lock(f"locks/{job_tag}-{traj}-{fname}-{inv_type_name}"):
        return
    labels = contract_hlbl_four_labels(job_tag)
    pairs_data = []
    for fn_chunk in fn_chunk_list:
        pairs_data += q.load_pickle_obj(get_load_path(fn_chunk))
    if len(pairs_data) != len(get_point_pairs()):
        raise Exception(f"len(pairs_data)={len(pairs_data)} len(get_point_pairs())=len(get_point_pairs())")
    results = {}
    results["labels"] = labels
    results["pairs_data"] = pairs_data
    results["n_pairs"] = len(pairs_data)
    results["lslt_sum"] = sum([ d["lslt"] for d in pairs_data ])
    results["lslt_sloppy_sum"] = sum([ d["lslt_sloppy"] for d in pairs_data ])
    q.displayln_info(-1, f"{fname}: {job_tag} {traj} {inv_type_name}\n", show_lslt(labels, results["lslt_sum"]))
    json_results.append((
        f"{fname}: {job_tag} {traj} {inv_type_name} lslt_sum",
        q.get_data_sig(results["lslt_sum"], q.RngState()),
        5e-3,
        ))
    json_results.append((
        f"{fname}: {job_tag} {traj} {inv_type_name} lslt_sloppy_sum",
        q.get_data_sig(results["lslt_sloppy_sum"], q.RngState()),
        5e-3,
        ))
    json_results.append((
        f"{fname}: {job_tag} {traj} {inv_type_name} lslt_sum[labels.index('ref-far proj-all'), -1, -1]",
        results["lslt_sum"][labels.index('ref-far proj-all'), -1, -1],
        5e-3,
        ))
    json_results.append((
        f"{fname}: {job_tag} {traj} {inv_type_name} lslt_sloppy_sum[labels.index('ref-far proj-all'), -1, -1]",
        results["lslt_sloppy_sum"][labels.index('ref-far proj-all'), -1, -1],
        5e-3,
        ))
    q.save_pickle_obj(results, get_save_path(fn))
    for fn_chunk in fn_chunk_list:
        q.qremove_info(get_load_path(fn_chunk))
    q.displayln_info(f"{fname}: {job_tag} {traj} {inv_type_name} done.")
    q.release_lock()

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
            f"{job_tag}/auto-contract/traj-{traj}/checkpoint.txt",
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
    run_ret_list = []
    #
    def add_to_run_ret_list(v):
        nonlocal run_ret_list
        if v is None:
            return
        if isinstance(v, list):
            run_ret_list += v
        else:
            run_ret_list.append(v)
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
    get_f_weight = run_f_weight_from_wsrc_prop_full(job_tag, traj, get_wi=get_wi)
    get_f_rand_01 = run_f_rand_01(job_tag, traj)
    get_fsel_prob = run_fsel_prob(job_tag, traj, get_f_rand_01=get_f_rand_01, get_f_weight=get_f_weight)
    get_fsel = lambda : get_fsel_prob().fsel
    get_psel_prob = run_psel_prob(job_tag, traj, get_f_rand_01=get_f_rand_01, get_f_weight=get_f_weight)
    get_psel = lambda : get_psel_prob().psel
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
    run_r_list(job_tag)
    #
    get_gf = run_gf(job_tag, traj_gf)
    get_gt = run_gt(job_tag, traj_gf, get_gf)
    get_gf_ape = run_gf_ape(job_tag, get_gf)
    #
    get_wi = run_wi(job_tag, traj)
    #
    get_f_weight = run_f_weight_from_wsrc_prop_full(job_tag, traj, get_wi=get_wi)
    get_f_rand_01 = run_f_rand_01(job_tag, traj)
    get_psel_prob = run_psel_prob(job_tag, traj, get_f_rand_01=get_f_rand_01, get_f_weight=get_f_weight)
    get_psel = lambda : get_psel_prob().psel
    get_fsel_prob = run_fsel_prob(job_tag, traj, get_f_rand_01=get_f_rand_01, get_f_weight=get_f_weight)
    get_fsel = lambda : get_fsel_prob().fsel
    #
    get_psel_smear = run_psel_smear(job_tag, traj)
    #
    get_glb_hvp_avg_light = run_job_global_hvp_average(job_tag, inv_type=0)
    get_glb_hvp_avg_strange = run_job_global_hvp_average(job_tag, inv_type=1)
    get_hvp_average_light = run_hvp_average(job_tag, traj, inv_type=0, get_psel_prob=get_psel_prob)
    get_hvp_average_strange = run_hvp_average(job_tag, traj, inv_type=1, get_psel_prob=get_psel_prob)
    #
    get_glb_hvp_avg_for_sub_light = run_job_global_hvp_average_for_subtract(job_tag, traj, get_glb_hvp_avg=get_glb_hvp_avg_light, get_hvp_average=get_hvp_average_light)
    get_glb_hvp_avg_for_sub_strange = run_job_global_hvp_average_for_subtract(job_tag, traj, get_glb_hvp_avg=get_glb_hvp_avg_strange, get_hvp_average=get_hvp_average_strange)
    #
    if job_tag[:5] == "test-":
        force_load_muon_line_interpolation()
    #
    get_point_pairs_light = run_hlbl_four_point_pairs_info(job_tag, traj, inv_type=0, get_psel_prob=get_psel_prob)
    get_point_pairs_strange = run_hlbl_four_point_pairs_info(job_tag, traj, inv_type=1, get_psel_prob=get_psel_prob)
    #
    run_hlbl_four(job_tag, traj, inv_type=0, get_psel_prob=get_psel_prob, get_fsel_prob=get_fsel_prob, get_point_pairs=get_point_pairs_light)
    run_hlbl_four(job_tag, traj, inv_type=1, get_psel_prob=get_psel_prob, get_fsel_prob=get_fsel_prob, get_point_pairs=get_point_pairs_strange)
    #
    get_get_prop = run_get_prop(job_tag, traj,
            get_gf=get_gf,
            get_gt=get_gt,
            get_psel=get_psel,
            get_fsel=get_fsel,
            get_psel_smear=get_psel_smear,
            get_wi=get_wi,
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
    if q.obtain_lock(f"locks/{job_tag}-{traj}-{fname}"):
        get_prop = get_get_prop()
        if get_prop is not None:
            q.timer_fork()
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
            q.timer_display()
            q.timer_merge()
        q.release_lock()
    q.clean_cache()

# ----

set_param("test-4nt8", "mk_sample_gauge_field", "rand_n_step", value=2)
set_param("test-4nt8", "mk_sample_gauge_field", "flow_n_step", value=8)
set_param("test-4nt8", "mk_sample_gauge_field", "hmc_n_traj", value=1)
set_param("test-4nt8", "fermion_params", 0, 2, "Ls", value=8)
set_param("test-4nt8", "fermion_params", 1, 2, "Ls", value=8)
set_param("test-4nt8", "fermion_params", 2, 2, "Ls", value=8)
set_param("test-4nt8", "cg_params-0-0", "maxiter", value=5)
set_param("test-4nt8", "cg_params-0-1", "maxiter", value=5)
set_param("test-4nt8", "cg_params-0-2", "maxiter", value=5)
set_param("test-4nt8", "cg_params-1-0", "maxiter", value=5)
set_param("test-4nt8", "cg_params-1-1", "maxiter", value=5)
set_param("test-4nt8", "cg_params-1-2", "maxiter", value=5)
set_param("test-4nt8", "cg_params-0-0", "maxcycle", value=1)
set_param("test-4nt8", "cg_params-0-1", "maxcycle", value=2)
set_param("test-4nt8", "cg_params-0-2", "maxcycle", value=3)
set_param("test-4nt8", "cg_params-1-0", "maxcycle", value=1)
set_param("test-4nt8", "cg_params-1-1", "maxcycle", value=2)
set_param("test-4nt8", "cg_params-1-2", "maxcycle", value=3)
set_param("test-4nt8", "cg_params-0-2", "pv_maxiter", value=5)
set_param("test-4nt8", "cg_params-1-2", "pv_maxiter", value=5)
set_param("test-4nt8", "a_inv_gev", value=1.73)
set_param("test-4nt8", "zz_vv", value=0.71)

set_param("24D", "lanc_params", 1, value=None)
set_param("24D", "clanc_params", 1, value=None)
set_param("24D", 'fermion_params', 0, 2, value=deepcopy(get_param("24D", 'fermion_params', 0, 0)))
set_param("24D", 'fermion_params', 1, 0, value=deepcopy(get_param("24D", 'fermion_params', 1, 2)))
set_param("24D", 'fermion_params', 1, 1, value=deepcopy(get_param("24D", 'fermion_params', 1, 2)))

tag = "trajs"
set_param("test-4nt8", tag, value=[ 1000, 2000, ])
set_param("test-8nt16", tag, value=[ 1000, 2000, ])
set_param("24D", tag, value=list(range(2000, 3000, 10)))
set_param("48I", tag, value=list(range(975, 2185, 10)) + list(range(1102, 1502, 10)))
set_param("64I", tag, value=list(range(1200, 3680, 40)))

tag = "hlbl_four_prob_scaling_factor"
set_param("test-4nt8", tag, value=1.0)
set_param("test-8nt16", tag, value=1.0)
set_param("24D", tag, value=1.0)
set_param("48I", tag, value=1.0)
set_param("64I", tag, value=1.0)

tag = "hlbl_four_prob_scaling_factor_strange"
set_param("test-4nt8", tag, value=1.0)
set_param("test-8nt16", tag, value=1.0)
set_param("24D", tag, value=1.0)
set_param("48I", tag, value=1.0)
set_param("64I", tag, value=1.0)

tag = "hlbl_four_num_chunk"
set_param("test-4nt8", tag, value=3)
set_param("test-8nt16", tag, value=6)
set_param("24D", tag, value=8)
set_param("48I", tag, value=8)
set_param("64I", tag, value=8)

# ----

if __name__ == "__main__":

    qg.begin_with_gpt()

    ##################### CMD options #####################

    job_tags = q.get_arg("--job_tags", default="").split(",")

    is_performing_inversion = q.get_arg("--no-inversion", default=None) is None

    is_performing_contraction = q.get_arg("--no-contract", default=None) is None

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
            q.check_time_limit()
            if is_performing_inversion:
                run_job(job_tag, traj)
        if is_performing_contraction:
            run_job_global_hvp_average(job_tag, inv_type=0)
            run_job_global_hvp_average(job_tag, inv_type=1)
        for traj in get_param(job_tag, "trajs"):
            q.check_time_limit()
            if is_performing_contraction:
                q.check_time_limit()
                run_job_contract(job_tag, traj)

    q.check_log_json(__file__, json_results)

    q.timer_display()

    qg.end_with_gpt()

    q.displayln_info("CHECK: finished successfully.")

# ----
