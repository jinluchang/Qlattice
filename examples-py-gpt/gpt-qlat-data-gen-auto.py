#!/usr/bin/env python3

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
    q.json_results_append(f"{fname}: ld sig", q.get_data_sig(ld, q.RngState()))
    for i, en in enumerate(expr_names):
        q.json_results_append(f"{fname}: ld '{en}' sig", q.get_data_sig(ld[i], q.RngState()))

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
    q.displayln_info(f"{fname}: timer_display for parallel_map_sum")
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
    q.json_results_append(f"{fname}: ld sig", q.get_data_sig(ld, q.RngState()))
    for i, en in enumerate(expr_names):
        q.json_results_append(f"{fname}: ld '{en}' sig", q.get_data_sig(ld[i], q.RngState()))

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
    q.displayln_info(f"{fname}: timer_display for parallel_map_sum")
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
    q.json_results_append(f"{fname}: ld sig", q.get_data_sig(ld, q.RngState()))
    for i, en in enumerate(expr_names):
        q.json_results_append(f"{fname}: ld '{en}' sig", q.get_data_sig(ld[i], q.RngState()))

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
    q.displayln_info(f"{fname}: timer_display for parallel_map_sum")
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
    q.json_results_append(f"{fname}: ld sig", q.get_data_sig(ld, q.RngState()))
    for i, en in enumerate(expr_names):
        q.json_results_append(f"{fname}: ld '{en}' sig", q.get_data_sig(ld[i], q.RngState()))

# ----

@q.timer
def get_cexpr_meson_jt():
    fn_base = "cache/auto_contract_cexpr/get_cexpr_meson_jt"
    def calc_cexpr():
        diagram_type_dict = dict()
        diagram_type_dict[((('t_1', 't_2'), 1), (('t_2', 'x'), 1), (('x', 't_1'), 1))] = 'Type1'
        diagram_type_dict[((('t_1', 't_2'), 1), (('t_2', 't_1'), 1), (('x', 'x'), 1))] = None
        diagram_type_dict[((('t_1p', 't_2p'), 1), (('t_2p', 'x'), 1), (('x', 't_1p'), 1))] = 'Type2'
        diagram_type_dict[((('t_1p', 't_2p'), 1), (('t_2p', 't_1p'), 1), (('x', 'x'), 1))] = None
        mm_list = [
                mk_pi_p("t_2", True) * mk_pi_p("t_1") + "pi+^dag(+tsep) * pi+(-tsep)",
                mk_k_p("t_2", True) * mk_k_p("t_1") + "K+^dag(+tsep) * K+(-tsep)",
                mk_pi_p("t_2p", True) * mk_pi_p("t_1p") + "pi+^dag(T/2+tsep) * pi+(T/2-tsep)",
                mk_k_p("t_2p", True) * mk_k_p("t_1p") + "K+^dag(T/2+tsep) * K+(T/2-tsep)",
                ]
        op_list = [
                mk_vec_mu("u", "u", "x", 3) + "ubar_gt_u(0)",
                mk_vec_mu("s", "s", "x", 3) + "sbar_gt_s(0)",
                ]
        exprs = [
                mk_expr(1) + f"1",
                op_list[0] * mm_list[0],
                op_list[0] * mm_list[1],
                -op_list[1] * mm_list[1],
                op_list[0] * mm_list[2],
                op_list[0] * mm_list[3],
                -op_list[1] * mm_list[3],
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit=True, diagram_type_dict=diagram_type_dict)
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

@q.timer_verbose
def auto_contract_meson_jt(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_jt.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_jt()
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
    tsep = get_param(job_tag, "meson_tensor_tsep")
    geo = q.Geometry(total_site)
    total_volume = geo.total_volume
    def load_data():
        for idx in range(len(xg_fsel_arr)):
            yield idx
    @q.timer
    def feval(args):
        idx = args
        xg_snk = tuple(xg_fsel_arr[idx])
        prob_snk = fsel_prob_arr[idx]
        t = xg_snk[3]
        t_1 = (t - tsep) % total_site[3]
        t_2 = (t + tsep) % total_site[3]
        t_1p = (t_1 + total_site[3] // 2) % total_site[3]
        t_2p = (t_2 + total_site[3] // 2) % total_site[3]
        pd = {
                "x" : ("point-snk", xg_snk,),
                "t_1" : ("wall", t_1),
                "t_2" : ("wall", t_2),
                "t_1p" : ("wall", t_1p),
                "t_2p" : ("wall", t_2p),
                "size" : total_site,
                }
        val = eval_cexpr(cexpr, positions_dict=pd, get_prop=get_prop)
        return val / prob_snk
    def sum_function(val_list):
        values = np.zeros(len(expr_names), dtype=complex)
        for val in val_list:
            values += val
        return values
    auto_contractor_chunk_size = get_param(job_tag, "measurement", "auto_contractor_chunk_size", default=128)
    q.timer_fork(0)
    res_sum = q.glb_sum(
            q.parallel_map_sum(feval, load_data(), sum_function=sum_function, chunksize=auto_contractor_chunk_size))
    q.displayln_info(f"{fname}: timer_display for parallel_map_sum")
    q.timer_display()
    q.timer_merge()
    res_sum *= 1.0 / total_volume
    q.displayln_info(0, res_sum[0])
    ld_sum = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        ])
    ld_sum.from_numpy(res_sum)
    ld_sum.save(get_save_path(fn))
    q.json_results_append(f"{fname}: ld_sum sig", q.get_data_sig(ld_sum, q.RngState()))
    for i, en in enumerate(expr_names):
        q.json_results_append(f"{fname}: ld_sum '{en}' sig", q.get_data_sig(ld_sum[i], q.RngState()))

# ----

@q.timer
def get_cexpr_meson_m():
    fn_base = "cache/auto_contract_cexpr/get_cexpr_meson_m"
    def calc_cexpr():
        diagram_type_dict = dict()
        diagram_type_dict[((('t_1', 't_2'), 1), (('t_2', 't_1'), 1), (('x', 'x'), 1))] = None
        diagram_type_dict[((('t_1', 't_2'), 1), (('t_2', 'x'), 1), (('x', 't_1'), 1))] = 'Type1'
        mm_list = [
                mk_pi_0("t_2", True) * mk_pi_0("t_1") + "pi0^dag(+tsep) * pi0(-tsep)",
                mk_pi_p("t_2", True) * mk_pi_p("t_1") + "pi+^dag(+tsep) * pi+(-tsep)",
                mk_k_0("t_2", True) * mk_k_0("t_1") + "K0^dag(+tsep) * K0(-tsep)",
                mk_k_p("t_2", True) * mk_k_p("t_1") + "K+^dag(+tsep) * K+(-tsep)",
                ]
        m_list = [
                mk_m("u", "x") + "ubar_u(0)",
                mk_m("d", "x") + "dbar_d(0)",
                mk_m("s", "x") + "sbar_s(0)",
                ]
        exprs = [ mk_expr(1) + f"1", ]
        exprs += [ m * mm for mm in mm_list for m in m_list ]
        cexpr = contract_simplify_compile(
                *exprs,
                is_isospin_symmetric_limit=True,
                diagram_type_dict=diagram_type_dict)
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

@q.timer_verbose
def auto_contract_meson_m(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_m.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_m()
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
    tsep = get_param(job_tag, "meson_tensor_tsep")
    geo = q.Geometry(total_site)
    total_volume = geo.total_volume
    def load_data():
        for idx in range(len(xg_fsel_arr)):
            yield idx
    @q.timer
    def feval(args):
        idx = args
        xg_snk = tuple(xg_fsel_arr[idx])
        prob_snk = fsel_prob_arr[idx]
        t = xg_snk[3]
        t_2 = (t + tsep) % total_site[3]
        t_1 = (t - tsep) % total_site[3]
        pd = {
                "x" : ("point-snk", xg_snk,),
                "t_1" : ("wall", t_1,),
                "t_2" : ("wall", t_2,),
                "size" : total_site,
                }
        val = eval_cexpr(cexpr, positions_dict=pd, get_prop=get_prop)
        return val / prob_snk
    def sum_function(val_list):
        values = np.zeros(len(expr_names), dtype=complex)
        for val in val_list:
            values += val
        return values
    auto_contractor_chunk_size = get_param(job_tag, "measurement", "auto_contractor_chunk_size", default=128)
    q.timer_fork(0)
    res_sum = q.glb_sum(
            q.parallel_map_sum(feval, load_data(), sum_function=sum_function, chunksize=auto_contractor_chunk_size))
    q.displayln_info(f"{fname}: timer_display for parallel_map_sum")
    q.timer_display()
    q.timer_merge()
    total_volume = geo.total_volume
    res_sum *= 1.0 / total_volume
    q.displayln_info(0, res_sum[0])
    ld_sum = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        ])
    ld_sum.from_numpy(res_sum)
    ld_sum.save(get_save_path(fn))
    q.json_results_append(f"{fname}: ld_sum sig", q.get_data_sig(ld_sum, q.RngState()))
    for i, en in enumerate(expr_names):
        q.json_results_append(f"{fname}: ld_sum '{en}' sig", q.get_data_sig(ld_sum[i], q.RngState()))

# ----

@q.timer
def get_cexpr_meson_jj():
    fn_base = "cache/auto_contract_cexpr/get_cexpr_meson_jj"
    def calc_cexpr():
        diagram_type_dict = dict()
        diagram_type_dict[((('t_1', 'x_1'), 1), (('t_2', 'x_2'), 1), (('x_1', 't_1'), 1), (('x_2', 't_2'), 1))] = 'Type1'
        diagram_type_dict[((('t_1', 'x_1'), 1), (('t_2', 'x_2'), 1), (('x_1', 't_2'), 1), (('x_2', 't_1'), 1))] = 'Type2'
        diagram_type_dict[((('t_1', 't_2'), 1), (('t_2', 'x_1'), 1), (('x_1', 'x_2'), 1), (('x_2', 't_1'), 1))] = 'Type3'
        diagram_type_dict[((('t_1', 't_2'), 1), (('t_2', 'x_1'), 1), (('x_1', 't_1'), 1), (('x_2', 'x_2'), 1))] = None
        diagram_type_dict[((('t_1', 't_2'), 1), (('t_2', 't_1'), 1), (('x_1', 'x_1'), 1), (('x_2', 'x_2'), 1))] = None
        diagram_type_dict[((('t_1', 't_2'), 1), (('t_2', 't_1'), 1), (('x_1', 'x_2'), 1), (('x_2', 'x_1'), 1))] = None
        diagram_type_dict[((('t_1', 'x_1'), 1), (('x_1', 'x_2'), 1), (('x_2', 't_1'), 1))] = 'TypeD1'
        diagram_type_dict[((('t_2', 'x_1'), 1), (('x_1', 'x_2'), 1), (('x_2', 't_2'), 1))] = 'TypeD2'
        diagram_type_dict[((('t_1', 'x_1'), 1), (('x_1', 't_1'), 1), (('x_2', 'x_2'), 1))] = None
        diagram_type_dict[((('t_2', 'x_1'), 1), (('x_1', 't_2'), 1), (('x_2', 'x_2'), 1))] = None
        jj_list = [
                sum([ mk_j_mu("x_2", mu) * mk_j_mu("x_1", mu) for mu in range(4) ])
                + "j_mu(x) * j_mu(0)",
                #
                mk_j_mu("x_2", 3) * mk_j_mu("x_1", 3)
                + "j_t(x) * j_t(0)",
                #
                sum([ mk_j_mu("x_2", mu) * mk_j_mu("x_1", mu) for mu in range(3) ])
                + "j_i(x) * j_i(0)",
                #
                sum([
                    mk_fac(f"rel_mod_sym(x_2[1][{mu}] - x_1[1][{mu}], size[{mu}])")
                    * mk_fac(f"rel_mod_sym(x_2[1][{nu}] - x_1[1][{nu}], size[{nu}])")
                    * mk_j_mu("x_2", mu) * mk_j_mu("x_1", nu)
                    for mu in range(3) for nu in range(3) ])
                + "x[i] * x[j] * j_i(x) * j_j(0)",
                #
                sum([
                    mk_fac(f"rel_mod_sym(x_2[1][{mu}] - x_1[1][{mu}], size[{mu}])")
                    * mk_j_mu("x_2", mu) * mk_j_mu("x_1", 3)
                    for mu in range(3) ])
                + "x[i] * j_i(x) * j_t(0)",
                #
                sum([
                    mk_fac(f"rel_mod_sym(x_1[1][{mu}] - x_2[1][{mu}], size[{mu}])")
                    * mk_j_mu("x_1", mu) * mk_j_mu("x_2", 3)
                    for mu in range(3) ])
                + "-x[i] * j_i(-x) * j_t(0)",
                ]
        assert len(jj_list) == 6
        m2_list = [
                mk_m("u", "x_2") + "ubar_u(x)",
                mk_m("d", "x_2") + "dbar_d(x)",
                mk_m("s", "x_2") + "sbar_s(x)",
                ]
        assert len(m2_list) == 3
        m1_list = [
                mk_m("u", "x_1") + "ubar_u(0)",
                mk_m("d", "x_1") + "dbar_d(0)",
                mk_m("s", "x_1") + "sbar_s(0)",
                ]
        assert len(m1_list) == 3
        m1m2_list = [ m2 * m1 for m2 in m2_list for m1 in m1_list ]
        assert len(m1m2_list) == 9
        op_list = jj_list + m1m2_list
        assert len(op_list) == 15
        mm_list = [
                mk_pi_0("t_2", True) * mk_pi_0("t_1")
                + "pi0^dag(x[t]+tsep) * pi0(-tsep)",
                mk_sym(1)/2 * (mk_pi_p("t_2", True) * mk_pi_p("t_1") + mk_pi_m("t_2", True) * mk_pi_m("t_1"))
                + "pi+^dag(x[t]+tsep) * pi+(-tsep)",
                mk_sym(1)/2 * (mk_k_0("t_2", True) * mk_k_0("t_1") + mk_k_0_bar("t_2", True) * mk_k_0_bar("t_1"))
                + "K0^dag(x[t]+tsep) * K0(-tsep)",
                mk_sym(1)/2 * (mk_k_p("t_2", True) * mk_k_p("t_1") + mk_k_m("t_2", True) * mk_k_m("t_1"))
                + "K+^dag(x[t]+tsep) * K+(-tsep)",
                ]
        assert len(mm_list) == 4
        exprs_self_energy = [ op * mm for mm in mm_list for op in op_list ]
        assert len(exprs_self_energy) == 60
        #
        op_l_ope_list = [
                sum([ mk_vec_mu("d", "d", "x_2", mu) * mk_vec_mu("d", "d", "x_1", mu) for mu in range(4) ])
                + "jd_mu(x) * jd_mu(0)",
                mk_vec_mu("d", "d", "x_2", 3) * mk_vec_mu("d", "d", "x_1", 3)
                + "jd_t(x) * jd_t(0)",
                ]
        assert len(op_l_ope_list) == 2
        op_s_ope_list = [
                sum([ mk_vec_mu("s", "s", "x_2", mu) * mk_vec_mu("s", "s", "x_1", mu) for mu in range(4) ])
                + "js_mu(x) * js_mu(0)",
                mk_vec_mu("s", "s", "x_2", 3) * mk_vec_mu("s", "s", "x_1", 3)
                + "js_t(x) * js_t(0)",
                ]
        assert len(op_s_ope_list) == 2
        mm_l_ope_list = [
                mk_sym(1)/2 * (mk_pi_p("t_2", True) * mk_pi_p("t_1") + mk_pi_m("t_2", True) * mk_pi_m("t_1"))
                + "pi+^dag(x[t]+tsep) * pi+(-tsep)",
                ]
        assert len(mm_l_ope_list) == 1
        mm_s_ope_list = [
                mk_sym(1)/2 * (mk_k_p("t_2", True) * mk_k_p("t_1") + mk_k_m("t_2", True) * mk_k_m("t_1"))
                + "K+^dag(x[t]+tsep) * K+(-tsep)",
                ]
        assert len(mm_s_ope_list) == 1
        exprs_ope = [ op * mm for mm in mm_l_ope_list for op in op_l_ope_list ] + [ op * mm for mm in mm_s_ope_list for op in op_s_ope_list ]
        assert len(exprs_ope) == 4
        #
        jwj_list = [
                mk_jw_a_mu("x_1", 3) * mk_j_mu("x_2", 3)
                + "jw_a_t(0) * j_t(x)",
                #
                mk_sw5("x_1") * mk_j_mu("x_2", 3)
                + "sw5(0) * j_t(x)",
                #
                sum([
                    mk_jw_a_mu("x_1", mu) * mk_j_mu("x_2", mu)
                    for mu in range(3) ])
                + "jw_a_i(0) * j_i(x)",
                #
                sum([
                    mk_fac(f"rel_mod_sym(x_2[1][{mu}] - x_1[1][{mu}], size[{mu}])")
                    * mk_jw_a_mu("x_1", 3) * mk_j_mu("x_2", mu)
                    for mu in range(3) ])
                + "x[i] * jw_a_t(0) * j_i(x)",
                #
                sum([
                    mk_fac(f"rel_mod_sym(x_2[1][{mu}] - x_1[1][{mu}], size[{mu}])")
                    * mk_sw5("x_1") * mk_j_mu("x_2", mu)
                    for mu in range(3) ])
                + "x[i] * sw5(0) * j_i(x)",
                #
                sum([
                    mk_fac(f"rel_mod_sym(x_2[1][{mu}] - x_1[1][{mu}], size[{mu}])")
                    * mk_jw_a_mu("x_1", mu) * mk_j_mu("x_2", 3)
                    for mu in range(3) ])
                + "x[i] * jw_a_i(0) * j_t(x)",
                #
                sum([
                    mk_fac(f"rel_mod_sym(x_2[1][{mu}] - x_1[1][{mu}], size[{mu}])")
                    * mk_fac(f"rel_mod_sym(x_2[1][{nu}] - x_1[1][{nu}], size[{nu}])")
                    * mk_jw_a_mu("x_1", mu) * mk_j_mu("x_2", nu)
                    for mu in range(3) for nu in range(3) ])
                + "x[i] * x[j] * jw_a_i(0) * j_j(x)",
                #
                sum([
                    q.epsilon_tensor(mu, nu, rho)
                    * mk_fac(f"rel_mod_sym(x_2[1][{mu}] - x_1[1][{mu}], size[{mu}])")
                    * mk_jw_v_mu("x_1", nu) * mk_j_mu("x_2", rho)
                    for mu in range(3) for nu in range(3) for rho in range(3) ])
                + "e(i,j,k) * x[i] * jw_v_j(0) * j_k(x)",
                ]
        assert len(jwj_list) == 8
        jjw_list = [
                mk_jw_a_mu("x_2", 3) * mk_j_mu("x_1", 3)
                + "jw_a_t(0) * j_t(-x)",
                #
                mk_sw5("x_2") * mk_j_mu("x_1", 3)
                + "sw5(0) * j_t(-x)",
                #
                sum([
                    mk_jw_a_mu("x_2", mu) * mk_j_mu("x_1", mu)
                    for mu in range(3) ])
                + "jw_a_i(0) * j_i(-x)",
                #
                sum([
                    mk_fac(f"rel_mod_sym(x_1[1][{mu}] - x_2[1][{mu}], size[{mu}])")
                    * mk_jw_a_mu("x_2", 3) * mk_j_mu("x_1", mu)
                    for mu in range(3) ])
                + "-x[i] * jw_a_t(0) * j_i(-x)",
                #
                sum([
                    mk_fac(f"rel_mod_sym(x_1[1][{mu}] - x_2[1][{mu}], size[{mu}])")
                    * mk_sw5("x_2") * mk_j_mu("x_1", mu)
                    for mu in range(3) ])
                + "-x[i] * sw5(0) * j_i(-x)",
                #
                sum([
                    mk_fac(f"rel_mod_sym(x_1[1][{mu}] - x_2[1][{mu}], size[{mu}])")
                    * mk_jw_a_mu("x_2", mu) * mk_j_mu("x_1", 3)
                    for mu in range(3) ])
                + "-x[i] * jw_a_i(0) * j_t(-x)",
                #
                sum([
                    mk_fac(f"rel_mod_sym(x_1[1][{mu}] - x_2[1][{mu}], size[{mu}])")
                    * mk_fac(f"rel_mod_sym(x_1[1][{nu}] - x_2[1][{nu}], size[{nu}])")
                    * mk_jw_a_mu("x_2", mu) * mk_j_mu("x_1", nu)
                    for mu in range(3) for nu in range(3) ])
                + "-x[i] * -x[j] * jw_a_i(0) * j_j(-x)",
                #
                sum([
                    q.epsilon_tensor(mu, nu, rho)
                    * mk_fac(f"rel_mod_sym(x_1[1][{mu}] - x_2[1][{mu}], size[{mu}])")
                    * mk_jw_v_mu("x_2", nu) * mk_j_mu("x_1", rho)
                    for mu in range(3) for nu in range(3) for rho in range(3) ])
                + "e(i,j,k) * -x[i] * jw_v_j(0) * j_k(-x)",
                ]
        assert len(jjw_list) == 8
        md_list = [
                mk_pi_p("t_1") + "pi+(-tsep)",
                mk_k_p("t_1") + "K+(-tsep)",
                mk_pi_p("t_2") + "pi+(x[t]+tsep)",
                mk_k_p("t_2") + "K+(x[t]+tsep)",
                ]
        assert len(md_list) == 4
        exprs_decay1 = [ jwj * md for md in md_list for jwj in jwj_list ]
        assert len(exprs_decay1) == 32
        exprs_decay2 = [ jjw * md for md in md_list for jjw in jjw_list ]
        assert len(exprs_decay2) == 32
        #
        jwm_list = [
                mk_jw_a_mu("x_1", 3) * mk_m("u", "x_2") + "jw_a_t(0) ubar_u(x)",
                mk_jw_a_mu("x_1", 3) * mk_m("d", "x_2") + "jw_a_t(0) dbar_d(x)",
                mk_jw_a_mu("x_1", 3) * mk_m("s", "x_2") + "jw_a_t(0) sbar_s(x)",
                mk_jw_a_mu("x_2", 3) * mk_m("u", "x_1") + "jw_a_t(0) ubar_u(-x)",
                mk_jw_a_mu("x_2", 3) * mk_m("d", "x_1") + "jw_a_t(0) dbar_d(-x)",
                mk_jw_a_mu("x_2", 3) * mk_m("s", "x_1") + "jw_a_t(0) sbar_s(-x)",
                mk_sw5("x_1") * mk_m("u", "x_2") + "sw5(0) ubar_u(x)",
                mk_sw5("x_1") * mk_m("d", "x_2") + "sw5(0) dbar_d(x)",
                mk_sw5("x_1") * mk_m("s", "x_2") + "sw5(0) sbar_s(x)",
                mk_sw5("x_2") * mk_m("u", "x_1") + "sw5(0) ubar_u(-x)",
                mk_sw5("x_2") * mk_m("d", "x_1") + "sw5(0) dbar_d(-x)",
                mk_sw5("x_2") * mk_m("s", "x_1") + "sw5(0) sbar_s(-x)",
                ]
        assert len(jwm_list) == 12
        exprs_decay_m = [ jwm * md for md in md_list for jwm in jwm_list ]
        assert len(exprs_decay_m) == 48
        #
        jj_d_list = [
                sum([
                    q.epsilon_tensor(mu, nu, rho)
                    * mk_fac(f"rel_mod_sym(x_2[1][{mu}] - x_1[1][{mu}], size[{mu}])")
                    * mk_j_mu("x_2", nu) * mk_j_mu("x_1", rho)
                    for mu in range(3) for nu in range(3) for rho in range(3) ])
                + "e(i,j,k) * x[i] * j_j(x) * j_k(0)",
                ]
        assert len(jj_d_list) == 1
        pi0d_list = [
                mk_pi_0("t_1") + "pi0(-tsep)",
                mk_pi_0("t_2") + "pi0(x[t]+tsep)",
                ]
        assert len(pi0d_list) == 2
        exprs_pi0_decay = [ jj_d * pi0d for pi0d in pi0d_list for jj_d in jj_d_list ]
        assert len(exprs_pi0_decay) == 2
        #
        exprs = [ mk_expr(1) + f"1", ]
        exprs += exprs_self_energy + exprs_ope + exprs_decay1 + exprs_decay2 + exprs_decay_m + exprs_pi0_decay
        assert len(exprs) == 179
        cexpr = contract_simplify_compile(
                *exprs,
                is_isospin_symmetric_limit=True,
                diagram_type_dict=diagram_type_dict)
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

@q.timer_verbose
def auto_contract_meson_jj(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_jj.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_jj()
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
    tsep = get_param(job_tag, "meson_tensor_tsep")
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
            x_rel_t = x_rel[3]
            x_2_t = xg_src[3]
            x_1_t = x_2_t + x_rel_t
            t_2 = (max(x_1_t, x_2_t) + tsep) % total_site[3]
            t_1 = (min(x_1_t, x_2_t) - tsep) % total_site[3]
            pd = {
                    "x_2" : ("point-snk", xg_snk,),
                    "x_1" : ("point", xg_src,),
                    "t_2" : ("wall", t_2),
                    "t_1" : ("wall", t_1),
                    "size" : total_site,
                    }
            t = x_rel_t % t_size
            val = eval_cexpr(cexpr, positions_dict=pd, get_prop=get_prop)
            res_list.append((val / prob, t, r_idx_low, r_idx_high, coef_low, coef_high,))
        return res_list
    def sum_function(val_list):
        values = np.zeros((t_size, len(r_list), len(expr_names),), dtype=complex)
        for idx, res_list in enumerate(val_list):
            for val, t, r_idx_low, r_idx_high, coef_low, coef_high in res_list:
                values[t, r_idx_low] += coef_low * val
                values[t, r_idx_high] += coef_high * val
            q.displayln_info(f"{fname}: {idx+1}/{len(xg_psel_arr)}")
        return q.glb_sum(values.transpose(2, 0, 1))
    q.timer_fork(0)
    res_sum = q.parallel_map_sum(feval, load_data(), sum_function=sum_function, chunksize=1)
    q.displayln_info(f"{fname}: timer_display for parallel_map_sum")
    q.timer_display()
    q.timer_merge()
    res_sum *= 1.0 / total_volume
    ld_sum = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t", t_size, [ str(q.rel_mod(t, t_size)) for t in range(t_size) ], ],
        [ "r", len(r_list), [ f"{r:.5f}" for r in r_list ], ],
        ])
    ld_sum.from_numpy(res_sum)
    ld_sum.save(get_save_path(fn))
    q.json_results_append(f"{fname}: ld_sum sig", q.get_data_sig(ld_sum, q.RngState()))
    for i, en in enumerate(expr_names):
        q.json_results_append(f"{fname}: ld_sum '{en}' sig", q.get_data_sig(ld_sum[i], q.RngState()))

# ----

@q.timer
def get_cexpr_meson_jwjj():
    fn_base = "cache/auto_contract_cexpr/get_cexpr_meson_jwjj"
    def calc_cexpr():
        diagram_type_dict = dict()
        diagram_type_dict[((('t_1', 'x_1'), 1), (('w', 'x_2'), 1), (('x_1', 'w'), 1), (('x_2', 't_1'), 1))] = 'Type1'
        diagram_type_dict[((('t_1', 'w'), 1), (('w', 'x_1'), 1), (('x_1', 'x_2'), 1), (('x_2', 't_1'), 1))] = 'Type2'
        diagram_type_dict[((('t_1', 'x_1'), 1), (('w', 't_1'), 1), (('x_1', 'x_2'), 1), (('x_2', 'w'), 1))] = 'Type2'
        diagram_type_dict[((('t_1', 'w'), 1), (('w', 'x_1'), 1), (('x_1', 't_1'), 1), (('x_2', 'x_2'), 1))] = None
        diagram_type_dict[((('t_1', 'x_1'), 1), (('w', 't_1'), 1), (('x_1', 'w'), 1), (('x_2', 'x_2'), 1))] = None
        diagram_type_dict[((('t_1', 'w'), 1), (('w', 't_1'), 1), (('x_1', 'x_1'), 1), (('x_2', 'x_2'), 1))] = None
        diagram_type_dict[((('t_1', 'w'), 1), (('w', 't_1'), 1), (('x_1', 'x_2'), 1), (('x_2', 'x_1'), 1))] = None
        diagram_type_dict[((('t_2', 'w'), 1), (('w', 'x_1'), 1), (('x_1', 't_2'), 1), (('x_2', 'x_2'), 1))] = None
        diagram_type_dict[((('t_2', 'x_1'), 1), (('w', 'x_2'), 1), (('x_1', 'w'), 1), (('x_2', 't_2'), 1))] = 'Type1'
        diagram_type_dict[((('t_2', 'w'), 1), (('w', 'x_1'), 1), (('x_1', 'x_2'), 1), (('x_2', 't_2'), 1))] = 'Type2'
        diagram_type_dict[((('t_2', 'x_1'), 1), (('w', 't_2'), 1), (('x_1', 'x_2'), 1), (('x_2', 'w'), 1))] = 'Type2'
        diagram_type_dict[((('t_2', 'x_1'), 1), (('w', 't_2'), 1), (('x_1', 'w'), 1), (('x_2', 'x_2'), 1))] = None
        diagram_type_dict[((('t_2', 'w'), 1), (('w', 't_2'), 1), (('x_1', 'x_1'), 1), (('x_2', 'x_2'), 1))] = None
        diagram_type_dict[((('t_2', 'w'), 1), (('w', 't_2'), 1), (('x_1', 'x_2'), 1), (('x_2', 'x_1'), 1))] = None
        jj_list = [
                (sum([ mk_j_mu("x_2", mu) * mk_j_mu("x_1", mu) for mu in range(4) ])
                    + "j_mu(x) * j_mu(y)"),
                (mk_j_mu("x_2", 3) * mk_j_mu("x_1", 3)
                    + "j_t(x) * j_t(y)"),
                ]
        assert len(jj_list) == 2
        jm_list = [
                mk_jw_a_mu("w", 3) * mk_pi_p("t_1") + "jw_a_t(0) * pi+(-tsep)",
                mk_jw_a_mu("w", 3) * mk_k_p("t_1") + "jw_a_t(0) * K+(-tsep)",
                mk_jw_a_mu("w", 3) * mk_pi_p("t_2") + "jw_a_t(0) * pi+(tsep)",
                mk_jw_a_mu("w", 3) * mk_k_p("t_2") + "jw_a_t(0) * K+(tsep)",
                mk_sw5("w") * mk_pi_p("t_1") + "sw5(0) * pi+(-tsep)",
                mk_sw5("w") * mk_k_p("t_1") + "sw5(0) * K+(-tsep)",
                mk_sw5("w") * mk_pi_p("t_2") + "sw5(0) * pi+(tsep)",
                mk_sw5("w") * mk_k_p("t_2") + "sw5(0) * K+(tsep)",
                ]
        assert len(jm_list) == 8
        exprs = [ mk_expr(1) + f"1", ]
        exprs += [ jj * jm for jm in jm_list for jj in jj_list ]
        assert len(exprs) == 17
        cexpr = contract_simplify_compile(
                *exprs,
                is_isospin_symmetric_limit=True,
                diagram_type_dict=diagram_type_dict)
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

@q.timer_verbose
def auto_contract_meson_jwjj(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_jwjj.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_jwjj()
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
    fsel_prob_inv_avg = np.average(1.0 / fsel_prob_arr)
    psel_prob_inv_avg = np.average(1.0 / psel_prob_arr)
    xg_fsel_arr = fsel.to_psel_local()[:]
    xg_psel_arr = psel[:]
    tsep = get_param(job_tag, "meson_tensor_tsep")
    geo = q.Geometry(total_site)
    total_volume = geo.total_volume
    r_list = get_r_list(job_tag)
    r_sq_interp_idx_coef_list = get_r_sq_interp_idx_coef_list(job_tag)
    n_elems = len(xg_fsel_arr)
    n_points = len(xg_psel_arr)
    n_pairs = n_points * (n_points - 1) // 2 + n_points
    #
    threshold = get_param(job_tag, "meson_jwjj_threshold")
    u_rand_prob = q.SelectedFieldRealD(fsel, 1)
    u_rand_prob.set_rand(q.RngState(f"auto_contract_meson_jwjj,{job_tag},{traj}"), 1.0, 0.0)
    u_rand_prob_arr = np.asarray(u_rand_prob).ravel()
    fn_meson_corr = f"{job_tag}/auto-contract/traj-{traj}/meson_corr_psnk.lat"
    if get_load_path(fn_meson_corr) is None:
        q.displayln_info(f"{fname}: '{fn_meson_corr}' does not exist. Skipping.")
        return
    ld_meson_corr = q.load_lat_data(get_load_path(fn_meson_corr))
    meson_corr_arr = ld_meson_corr.to_numpy()
    def get_prop_norm_sqrt(*args):
        is_sloppy = True
        return abs(ama_extract(get_prop(*args, is_norm_sqrt=True), is_sloppy=is_sloppy))
    def load_psrc_psrc_prop_norm_sqrt(flavor, i):
        xg1_src = tuple(xg_psel_arr[i])
        x_1 = ("point", xg1_src,)
        v_list = []
        for j in range(n_points):
            xg2_src = tuple(xg_psel_arr[j])
            x_2 = ("point", xg2_src,)
            v = get_prop_norm_sqrt(flavor, x_1, x_2)
            v_list.append(v)
        return np.array(v_list)
    psrc_psrc_prop_norm_sqrt = np.array([
        q.parallel_map(lambda i: load_psrc_psrc_prop_norm_sqrt(flavor, i), range(n_points), verbose=1)
        for flavor in [ "l", "s", ]
        ], dtype = float)
    def load_wsrc_psrc_prop_norm_sqrt(flavor, t):
        ts = ("wall", t,)
        v_list = []
        for j in range(n_points):
            xg2_src = tuple(xg_psel_arr[j])
            x_2 = ("point", xg2_src,)
            v = get_prop_norm_sqrt(flavor, x_2, ts)
            v_list.append(v)
        return np.array(v_list)
    wsrc_psrc_prop_norm_sqrt = np.array([
        q.parallel_map(lambda t: load_wsrc_psrc_prop_norm_sqrt(flavor, t), range(t_size), verbose=1)
        for flavor in [ "l", "s", ]
        ], dtype = float)
    def load_wsrc_psnk_prop_norm_sqrt(flavor, t):
        ts = ("wall", t,)
        v_list = []
        for j in range(n_elems):
            xg2_snk = tuple(xg_fsel_arr[j])
            x_2 = ("point-snk", xg2_snk,)
            v = get_prop_norm_sqrt(flavor, x_2, ts)
            v_list.append(v)
        return np.array(v_list)
    wsrc_psnk_prop_norm_sqrt = np.array([
        q.parallel_map(lambda t: load_wsrc_psnk_prop_norm_sqrt(flavor, t), range(t_size), verbose=1)
        for flavor in [ "l", "s", ]
        ], dtype = float)
    def load_psrc_psnk_prop_norm_sqrt(flavor, i):
        xg1_src = tuple(xg_psel_arr[i])
        x_1 = ("point", xg1_src,)
        v_list = []
        for j in range(n_elems):
            xg2_snk = tuple(xg_fsel_arr[j])
            x_2 = ("point-snk", xg2_snk,)
            v = get_prop_norm_sqrt(flavor, x_2, x_1)
            v_list.append(v)
        return np.array(v_list)
    psrc_psnk_prop_norm_sqrt = np.array([
        q.parallel_map(lambda i: load_psrc_psnk_prop_norm_sqrt(flavor, i), range(n_points), verbose=1)
        for flavor in [ "l", "s", ]
        ], dtype = float)
    def get_estimate(idx_snk, idx1, idx2, t_1, t_2):
        flavor_l = 0
        flavor_s = 1
        prob_1 = psel_prob_arr[idx1] * psel_prob_inv_avg
        prob_2 = psel_prob_arr[idx2] * psel_prob_inv_avg
        prob_w = fsel_prob_arr[idx_snk] * fsel_prob_inv_avg
        xg_snk_t = xg_fsel_arr[idx_snk, 3]
        corr1 = np.abs(meson_corr_arr[1, (xg_snk_t - t_1) % t_size])
        corr2 = np.abs(meson_corr_arr[1, (xg_snk_t - t_2) % t_size])
        p1t1 = wsrc_psrc_prop_norm_sqrt[flavor_l, t_1, idx1]
        p2t1 = wsrc_psrc_prop_norm_sqrt[flavor_l, t_1, idx2]
        wt1 = wsrc_psnk_prop_norm_sqrt[flavor_l, t_1, idx_snk]
        wt2 = wsrc_psnk_prop_norm_sqrt[flavor_l, t_2, idx_snk]
        p1t2 = wsrc_psrc_prop_norm_sqrt[flavor_l, t_2, idx1]
        p2t2 = wsrc_psrc_prop_norm_sqrt[flavor_l, t_2, idx2]
        p1p2 = psrc_psrc_prop_norm_sqrt[flavor_l, idx1, idx2]
        wp1 = psrc_psnk_prop_norm_sqrt[flavor_l, idx1, idx_snk]
        wp2 = psrc_psnk_prop_norm_sqrt[flavor_l, idx2, idx_snk]
        value = 0
        value += 2 * (p1t1 * p2t1 / corr1 + p1t2 * p2t2 / corr2) * (wp1 * wp2)
        value += 5 * (p1t1 * wt1 / corr1 + p1t2 * wt2 / corr2) * (p1p2 * wp2)
        value += 5 * (p2t1 * wt1 / corr1 + p2t2 * wt2 / corr2) * (p1p2 * wp1)
        value /= prob_1 * prob_w * prob_2
        assert np.all(value > 0)
        return value
    @q.timer
    def get_weight(idx_snk, idx1, idx2, t_1, t_2):
        """
        return weight for this point (1 / prob or zero)
        """
        idx_snk = np.asarray(idx_snk).ravel()
        est = get_estimate(idx_snk, idx1, idx2, t_1, t_2)
        assert est.shape == idx_snk.shape
        prob = est / threshold
        assert prob.shape == idx_snk.shape
        rand = u_rand_prob_arr[idx_snk]
        assert rand.shape == idx_snk.shape
        weight = 1.0 / prob
        weight[prob >= 1] = 1.0
        weight[rand >= prob] = 0.0
        return weight
    #
    def load_data():
        idx_pair = 0
        for idx1, xg1_src in enumerate(xg_psel_arr):
            xg1_src = tuple(xg1_src.tolist())
            for idx2, xg2_src in enumerate(xg_psel_arr):
                xg2_src = tuple(xg2_src.tolist())
                if idx2 > idx1:
                    continue
                idx_pair += 1
                yield idx1, idx2
    @q.timer
    def feval(args):
        idx1, idx2 = args
        xg1_src = tuple(xg_psel_arr[idx1])
        xg2_src = tuple(xg_psel_arr[idx2])
        #
        # Adaptive sampling method:
        prob_1 = psel_prob_arr[idx1]
        prob_2 = psel_prob_arr[idx2]
        if idx1 != idx2:
            prob = total_volume * prob_1 * prob_2
        else:
            prob = total_volume * prob_1
        weight_base = 1.0 / prob / total_volume
        #
        xg1_src_t = xg1_src[3]
        xg2_src_t = xg2_src[3]
        x_rel = [ q.rel_mod(xg2_src[mu] - xg1_src[mu], total_site[mu]) for mu in range(4) ]
        r_sq = q.get_r_sq(x_rel)
        idx_snk_arr = np.arange(n_elems)
        xg_t_arr = xg_fsel_arr[idx_snk_arr, 3]
        t_size_arr = np.broadcast_to(t_size, xg_t_arr.shape)
        xg1_xg_t_arr = q.rel_mod_arr(xg1_src_t - xg_t_arr, t_size_arr)
        xg2_xg_t_arr = q.rel_mod_arr(xg2_src_t - xg_t_arr, t_size_arr)
        t_1_arr = (np.minimum(0, np.minimum(xg1_xg_t_arr, xg2_xg_t_arr)) + xg_t_arr - tsep) % t_size
        t_2_arr = (np.maximum(0, np.maximum(xg1_xg_t_arr, xg2_xg_t_arr)) + xg_t_arr + tsep) % t_size
        weight_arr = weight_base * get_weight(idx_snk_arr, idx1, idx2, t_1_arr, t_2_arr)
        weight_arr[np.abs(xg2_xg_t_arr - xg1_xg_t_arr) >= t_size_arr // 2] = 0.0
        results = []
        for idx_snk in idx_snk_arr[weight_arr > 0]:
            xg_snk = tuple(xg_fsel_arr[idx_snk])
            prob_snk = fsel_prob_arr[idx_snk]
            weight = weight_arr[idx_snk]
            #
            # Adaptive sampling method:
            if xg_snk != xg1_src and xg_snk != xg2_src:
                weight = weight / prob_snk
            #
            xg_t = xg_t_arr[idx_snk]
            xg1_xg_t = xg1_xg_t_arr[idx_snk]
            xg2_xg_t = xg2_xg_t_arr[idx_snk]
            t_1 = t_1_arr[idx_snk]
            t_2 = t_2_arr[idx_snk]
            pd = {
                    "w" : ("point-snk", xg_snk,),
                    "x_1" : ("point", xg1_src,),
                    "x_2" : ("point", xg2_src,),
                    "t_1" : ("wall", t_1,),
                    "t_2" : ("wall", t_2,),
                    "size" : total_site,
                    }
            t1 = xg1_xg_t
            t2 = xg2_xg_t
            val = eval_cexpr(cexpr, positions_dict=pd, get_prop=get_prop)
            r_idx_low, r_idx_high, coef_low, coef_high = r_sq_interp_idx_coef_list[r_sq]
            results.append((weight * val, t1, t2, r_idx_low, r_idx_high, coef_low, coef_high,))
        return idx1, idx2, results
    def sum_function(val_list):
        n_total = 0
        n_selected = 0
        idx_pair = 0
        values = np.zeros((t_size, t_size, len(r_list), len(expr_names),), dtype=complex)
        for idx1, idx2, results in val_list:
            idx_pair += 1
            n_total += n_elems
            xg1_src = tuple(xg_psel_arr[idx1])
            xg2_src = tuple(xg_psel_arr[idx2])
            for val, t1, t2, r_idx_low, r_idx_high, coef_low, coef_high in results:
                n_selected += 1
                values[t1, t2, r_idx_low] += coef_low * val
                values[t1, t2, r_idx_high] += coef_high * val
            if idx_pair % (n_pairs // 1000 + 100) == 0:
                q.displayln_info(1, f"{fname}: {idx_pair}/{n_pairs} {xg1_src} {xg2_src} {len(results)}/{n_elems} n_total={n_total} n_selected={n_selected} ratio={n_selected/n_total}")
        q.displayln_info(1, f"{fname}: Final: n_total={n_total} n_selected={n_selected} ratio={n_selected/n_total}")
        return values.transpose(3, 0, 1, 2)
    q.timer_fork(0)
    res_sum = q.glb_sum(
            q.parallel_map_sum(feval, load_data(), sum_function=sum_function, chunksize=1))
    q.displayln_info(f"{fname}: timer_display")
    q.timer_display()
    q.timer_merge()
    res_sum *= 2.0 / (total_volume / t_size) # factor of 2 account for only calculating half of the pairs
    ld_sum = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t1", t_size, [ str(q.rel_mod(t, t_size)) for t in range(t_size) ], ],
        [ "t2", t_size, [ str(q.rel_mod(t, t_size)) for t in range(t_size) ], ],
        [ "r", len(r_list), [ f"{r:.5f}" for r in r_list ], ],
        ])
    ld_sum.from_numpy(res_sum)
    ld_sum.save(get_save_path(fn))
    q.json_results_append(f"{fname}: ld_sum sig", q.get_data_sig(ld_sum, q.RngState()))
    for i, en in enumerate(expr_names):
        q.json_results_append(f"{fname}: ld_sum '{en}' sig", q.get_data_sig(ld_sum[i], q.RngState()))

@q.timer_verbose
def auto_contract_meson_jwjj2(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_jwjj2.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_jwjj()
    expr_names = get_expr_names(cexpr)
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    t_size = total_site[3]
    point_distribution = load_point_distribution(job_tag)
    total_site_array = np.array(total_site.to_list())
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
    fsel_prob_inv_avg = np.average(1.0 / fsel_prob_arr)
    psel_prob_inv_avg = np.average(1.0 / psel_prob_arr)
    xg_fsel_arr = fsel.to_psel_local()[:]
    xg_psel_arr = psel[:]
    n_points = len(xg_psel_arr)
    tsep = get_param(job_tag, "meson_tensor_tsep")
    geo = q.Geometry(total_site)
    total_volume = geo.total_volume
    r_list = get_r_list(job_tag)
    r_sq_interp_idx_coef_list = get_r_sq_interp_idx_coef_list(job_tag)
    n_elems = len(xg_fsel_arr)
    n_pairs = n_points * n_points
    total_site_arr = np.array(total_site.to_list())
    total_site_arr = np.broadcast_to(total_site_arr, (n_elems, 4,))
    #
    threshold = get_param(job_tag, "meson_jwjj_threshold")
    u_rand_prob = q.SelectedFieldRealD(fsel, 1)
    u_rand_prob.set_rand(q.RngState(f"auto_contract_meson_jwjj2,{job_tag},{traj}"), 1.0, 0.0)
    u_rand_prob_arr = np.asarray(u_rand_prob).ravel()
    fn_meson_corr = f"{job_tag}/auto-contract/traj-{traj}/meson_corr_psnk.lat"
    if get_load_path(fn_meson_corr) is None:
        q.displayln_info(f"{fname}: '{fn_meson_corr}' does not exist. Skipping.")
        return
    ld_meson_corr = q.load_lat_data(get_load_path(fn_meson_corr))
    meson_corr_arr = ld_meson_corr.to_numpy()
    def get_prop_norm_sqrt(*args):
        is_sloppy = True
        return abs(ama_extract(get_prop(*args, is_norm_sqrt=True), is_sloppy=is_sloppy))
    def load_psrc_psrc_prop_norm_sqrt(flavor, i):
        xg1_src = tuple(xg_psel_arr[i])
        x_1 = ("point", xg1_src,)
        v_list = []
        for j in range(n_points):
            xg2_src = tuple(xg_psel_arr[j])
            x_2 = ("point", xg2_src,)
            v = get_prop_norm_sqrt(flavor, x_1, x_2)
            v_list.append(v)
        return np.array(v_list)
    psrc_psrc_prop_norm_sqrt = np.array([
        q.parallel_map(lambda i: load_psrc_psrc_prop_norm_sqrt(flavor, i), range(n_points), verbose=1)
        for flavor in [ "l", "s", ]
        ], dtype = float)
    def load_wsrc_psrc_prop_norm_sqrt(flavor, t):
        ts = ("wall", t,)
        v_list = []
        for j in range(n_points):
            xg2_src = tuple(xg_psel_arr[j])
            x_2 = ("point", xg2_src,)
            v = get_prop_norm_sqrt(flavor, x_2, ts)
            v_list.append(v)
        return np.array(v_list)
    wsrc_psrc_prop_norm_sqrt = np.array([
        q.parallel_map(lambda t: load_wsrc_psrc_prop_norm_sqrt(flavor, t), range(t_size), verbose=1)
        for flavor in [ "l", "s", ]
        ], dtype = float)
    def load_wsrc_psnk_prop_norm_sqrt(flavor, t):
        ts = ("wall", t,)
        v_list = []
        for j in range(n_elems):
            xg2_snk = tuple(xg_fsel_arr[j])
            x_2 = ("point-snk", xg2_snk,)
            v = get_prop_norm_sqrt(flavor, x_2, ts)
            v_list.append(v)
        return np.array(v_list)
    wsrc_psnk_prop_norm_sqrt = np.array([
        q.parallel_map(lambda t: load_wsrc_psnk_prop_norm_sqrt(flavor, t), range(t_size), verbose=1)
        for flavor in [ "l", "s", ]
        ], dtype = float)
    def load_psrc_psnk_prop_norm_sqrt(flavor, i):
        xg1_src = tuple(xg_psel_arr[i])
        x_1 = ("point", xg1_src,)
        v_list = []
        for j in range(n_elems):
            xg2_snk = tuple(xg_fsel_arr[j])
            x_2 = ("point-snk", xg2_snk,)
            v = get_prop_norm_sqrt(flavor, x_2, x_1)
            v_list.append(v)
        return np.array(v_list)
    psrc_psnk_prop_norm_sqrt = np.array([
        q.parallel_map(lambda i: load_psrc_psnk_prop_norm_sqrt(flavor, i), range(n_points), verbose=1)
        for flavor in [ "l", "s", ]
        ], dtype = float)
    def get_estimate(idx_w, idx_1, idx_2, xg_w_t, t_1, t_2):
        flavor_l = 0
        flavor_s = 1
        prob_1 = psel_prob_arr[idx_1] * psel_prob_inv_avg
        prob_w = psel_prob_arr[idx_w] * psel_prob_inv_avg
        prob_2 = fsel_prob_arr[idx_2] * fsel_prob_inv_avg
        corr1 = np.abs(meson_corr_arr[1, (xg_w_t - t_1) % t_size])
        corr2 = np.abs(meson_corr_arr[1, (xg_w_t - t_2) % t_size])
        p1t1 = wsrc_psrc_prop_norm_sqrt[flavor_l, t_1, idx_1]
        p2t1 = wsrc_psnk_prop_norm_sqrt[flavor_l, t_1, idx_2]
        wt1 = wsrc_psrc_prop_norm_sqrt[flavor_l, t_1, idx_w]
        wt2 = wsrc_psrc_prop_norm_sqrt[flavor_l, t_2, idx_w]
        p1t2 = wsrc_psrc_prop_norm_sqrt[flavor_l, t_2, idx_1]
        p2t2 = wsrc_psnk_prop_norm_sqrt[flavor_l, t_2, idx_2]
        p1p2 = psrc_psnk_prop_norm_sqrt[flavor_l, idx_1, idx_2]
        wp1 = psrc_psrc_prop_norm_sqrt[flavor_l, idx_1, idx_w]
        wp2 = psrc_psnk_prop_norm_sqrt[flavor_l, idx_w, idx_2]
        value = 0
        value += 2 * (p1t1 * p2t1 / corr1 + p1t2 * p2t2 / corr2) * (wp1 * wp2)
        value += 5 * (p1t1 * wt1 / corr1 + p1t2 * wt2 / corr2) * (p1p2 * wp2)
        value += 5 * (p2t1 * wt1 / corr1 + p2t2 * wt2 / corr2) * (p1p2 * wp1)
        value /= prob_1 * prob_w * prob_2
        assert np.all(value > 0)
        return value
    @q.timer
    def get_weight(idx_w, idx_1, idx_2, xg_w_t, t_1, t_2):
        """
        return weight for point ``idx_2`` (1 / prob or zero)
        """
        est = get_estimate(idx_w, idx_1, idx_2, xg_w_t, t_1, t_2)
        assert est.shape == idx_2.shape
        prob = est / threshold
        assert prob.shape == idx_2.shape
        rand = u_rand_prob_arr[idx_2]
        assert rand.shape == idx_2.shape
        weight = 1.0 / prob
        weight[prob >= 1] = 1.0
        weight[rand >= prob] = 0.0
        return weight
    #
    def load_data():
        idx_pair = 0
        for idx_1 in range(n_points):
            for idx_w in range(n_points):
                idx_pair += 1
                yield idx_1, idx_w
    @q.timer
    def feval(args):
        idx_1, idx_w = args
        idx_2_arr = np.arange(n_elems)
        xg_1 = tuple(xg_psel_arr[idx_1])
        xg_w = tuple(xg_psel_arr[idx_w])
        #
        # Old method to obtain the prob with previous sampling strategy:
        #
        # xg_rel_array = np.array(xg_1) - np.array(xg_w)
        # prob = get_point_xrel_prob(xg_rel_array, total_site_array, point_distribution, n_points)
        #
        # Adaptive sampling method:
        prob_1 = psel_prob_arr[idx_1]
        prob_w = psel_prob_arr[idx_w]
        if idx_1 != idx_w:
            prob = total_volume * prob_1 * prob_w
        else:
            prob = total_volume * prob_1
        #
        weight_base = 1.0 / prob / total_volume
        xg_1_arr = np.broadcast_to(np.array(xg_1), total_site_arr.shape)
        xg_2_arr = xg_fsel_arr[idx_2_arr]
        xg_w_arr = np.broadcast_to(np.array(xg_w), total_site_arr.shape)
        t_size_arr = total_site_arr[:, 3]
        xg_1_t_arr = xg_1_arr[:, 3]
        xg_2_t_arr = xg_2_arr[:, 3]
        xg_w_t_arr = xg_w_arr[:, 3]
        x_rel_arr = q.rel_mod_arr(xg_2_arr - xg_1_arr, total_site_arr)
        r_sq_arr = q.sqr(x_rel_arr[:, :3]).sum(-1)
        xg_1_xg_t_arr = q.rel_mod_arr(xg_1_t_arr - xg_w_t_arr, t_size_arr)
        xg_2_xg_t_arr = q.rel_mod_arr(xg_2_t_arr - xg_w_t_arr, t_size_arr)
        t_1_arr = (np.minimum(0, np.minimum(xg_1_xg_t_arr, xg_2_xg_t_arr)) + xg_w_t_arr - tsep) % t_size
        t_2_arr = (np.maximum(0, np.maximum(xg_1_xg_t_arr, xg_2_xg_t_arr)) + xg_w_t_arr + tsep) % t_size
        weight_arr = weight_base * get_weight(idx_w, idx_1, idx_2_arr, xg_w_t_arr, t_1_arr, t_2_arr)
        weight_arr[np.abs(xg_2_xg_t_arr - xg_1_xg_t_arr) >= t_size_arr // 2] = 0.0
        results = []
        for idx_2 in idx_2_arr[weight_arr > 0]:
            xg_2 = tuple(xg_fsel_arr[idx_2])
            weight = weight_arr[idx_2]
            #
            # Old method to obtain the prob with previous sampling strategy:
            # prob_2 = fsel_prob_arr[idx_2]
            #
            # Adaptive sampling method:
            if xg_2 != xg_1 and xg_2 != xg_w:
                prob_2 = fsel_prob_arr[idx_2]
            else:
                prob_2 = 1
            #
            r_sq = r_sq_arr[idx_2]
            xg_1_xg_t = xg_1_xg_t_arr[idx_2]
            xg_2_xg_t = xg_2_xg_t_arr[idx_2]
            t_1 = t_1_arr[idx_2]
            t_2 = t_2_arr[idx_2]
            pd = {
                    "w" : ("point", xg_w,),
                    "x_1" : ("point", xg_1,),
                    "x_2" : ("point-snk", xg_2,),
                    "t_1" : ("wall", t_1,),
                    "t_2" : ("wall", t_2,),
                    "size" : total_site.to_list(),
                    }
            t_1 = xg_1_xg_t
            t_2 = xg_2_xg_t
            val = eval_cexpr(cexpr, positions_dict=pd, get_prop=get_prop)
            r_idx_low, r_idx_high, coef_low, coef_high = r_sq_interp_idx_coef_list[r_sq]
            results.append((weight * val / prob_2, t_1, t_2, r_idx_low, r_idx_high, coef_low, coef_high,))
        return idx_1, idx_w, results
    def sum_function(val_list):
        n_total = 0
        n_selected = 0
        idx_pair = 0
        values = np.zeros((t_size, t_size, len(r_list), len(expr_names),), dtype = complex)
        for idx_1, idx_w, results in val_list:
            idx_pair += 1
            n_total += n_elems
            for val, t_1, t_2, r_idx_low, r_idx_high, coef_low, coef_high in results:
                n_selected += 1
                values[t_1, t_2, r_idx_low] += coef_low * val
                values[t_1, t_2, r_idx_high] += coef_high * val
            if idx_pair % (n_pairs // 1000 + 100) == 0:
                xg_1_src = tuple(xg_psel_arr[idx_1])
                xg_w_src = tuple(xg_psel_arr[idx_w])
                q.displayln_info(1, f"{fname}: {idx_pair}/{n_pairs} {xg_1_src} {xg_w_src} {len(results)}/{n_elems} n_total={n_total} n_selected={n_selected} ratio={n_selected/n_total}")
        q.displayln_info(1, f"{fname}: Final: n_total={n_total} n_selected={n_selected} ratio={n_selected/n_total}")
        n_total = q.glb_sum(n_total)
        n_selected = q.glb_sum(n_selected)
        q.displayln_info(1, f"{fname}: Final(glb_sum): n_total={n_total} n_selected={n_selected} ratio={n_selected/n_total}")
        return values.transpose(3, 0, 1, 2)
    q.timer_fork(0)
    res_sum = q.parallel_map_sum(feval, load_data(), sum_function = sum_function, chunksize = 1)
    res_sum = q.glb_sum(res_sum)
    q.displayln_info(f"{fname}: timer_display")
    q.timer_display()
    q.timer_merge()
    res_sum *= 1.0 / (total_volume / t_size)
    ld_sum = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t1", t_size, [ str(q.rel_mod(t, t_size)) for t in range(t_size) ], ],
        [ "t2", t_size, [ str(q.rel_mod(t, t_size)) for t in range(t_size) ], ],
        [ "r", len(r_list), [ f"{r:.5f}" for r in r_list ], ],
        ])
    ld_sum.from_numpy(res_sum)
    ld_sum.save(get_save_path(fn))
    q.json_results_append(f"{fname}: ld_sum sig", q.get_data_sig(ld_sum, q.RngState()))
    for i, en in enumerate(expr_names):
        q.json_results_append(f"{fname}: ld_sum '{en}' sig", q.get_data_sig(ld_sum[i], q.RngState()))

### ------

@q.timer
def get_cexpr_pi0_gg():
    fn_base = "cache/auto_contract_cexpr/get_cexpr_pi0_gg"
    def calc_cexpr():
        diagram_type_dict = dict()
        diagram_type_dict[()] = '1'
        diagram_type_dict[((('t_1', 'x_1'), 1), (('x_1', 't_1'), 1), (('x_2', 'x_2'), 1))] = 'TypeD'
        diagram_type_dict[((('t_1', 'x_1'), 1), (('x_1', 'x_2'), 1), (('x_2', 't_1'), 1))] = 'TypeC'
        diagram_type_dict[((('t_2', 'x_1'), 1), (('x_1', 't_2'), 1), (('x_2', 'x_2'), 1))] = 'TypeD'
        diagram_type_dict[((('t_2', 'x_1'), 1), (('x_1', 'x_2'), 1), (('x_2', 't_2'), 1))] = 'TypeC'
        #
        jj_d_list = [
                sum([
                    q.epsilon_tensor(mu, nu, rho)
                    * mk_fac(f"rel_mod_sym(x_2[1][{mu}] - x_1[1][{mu}], size[{mu}])")
                    * mk_j_mu("x_2", nu) * mk_j_mu("x_1", rho)
                    for mu in range(3) for nu in range(3) for rho in range(3) ])
                + "e(i,j,k) * x[i] * j_j(x) * j_k(0)",
                ]
        assert len(jj_d_list) == 1
        #
        pi0d_list = [
                mk_pi_0("t_1") + "pi0(-tsep)",
                mk_pi_0("t_2") + "pi0(x[t]+tsep)",
                ]
        assert len(pi0d_list) == 2
        #
        exprs_list_pi0_decay = [
                (jj_d * pi0d, None, "TypeC", "TypeD",)
                for pi0d in pi0d_list for jj_d in jj_d_list
                ]
        assert len(exprs_list_pi0_decay) == 2
        exprs_pi0_decay = [ e for e in exprs_list_pi0_decay ]
        #
        exprs = [
                mk_expr(1) + f"1",
                ]
        exprs += exprs_pi0_decay
        #
        cexpr = contract_simplify_compile(
                *exprs,
                is_isospin_symmetric_limit=True,
                diagram_type_dict=diagram_type_dict)
        return cexpr
        #
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

@q.timer
def get_cexpr_tadpole_current():
    """
    Intend to calculating < J_mu(x_1) >
    !!! Results needs to be multiplied by (m_s - m_l) and sum over "x_2" over the entire space-time volume !!!
    After summing over "x_2", the first and second expr should be same (and so is the following exprs). (We can average these two.)
    1/3 charge factor is already multiplied.
    """
    fn_base = "cache/auto_contract_cexpr/get_cexpr_tadpole_current"
    def calc_cexpr():
        diagram_type_dict = dict()
        diagram_type_dict[()] = 'T1'
        diagram_type_dict[((('x_1', 'x_2'), 1), (('x_2', 'x_1'), 1))] = 'TypeC'
        diagram_type_dict[((('x_1', 'x_1'), 1), (('x_2', 'x_2'), 1))] = None
        #
        jj_list = []
        jj_list += [
                1/3 * mk_vec_mu("l", "s", "x_1", mu) * mk_scalar("s", "l", "x_2")
                for mu in range(4)
                ]
        jj_list += [
                1/3 * mk_vec_mu("s", "l", "x_1", mu) * mk_scalar("l", "s", "x_2")
                for mu in range(4)
                ]
        jj_list += [
                1/3 * mk_vec_mu("l", "l", "x_1", mu) * mk_scalar("l", "l", "x_2")
                for mu in range(4)
                ]
        jj_list += [
                1/3 * mk_vec_mu("s", "s", "x_1", mu) * mk_scalar("s", "s", "x_2")
                for mu in range(4)
                ]
        #
        exprs = [
                mk_expr(1) + f"1",
                ]
        exprs += jj_list
        #
        cexpr = contract_simplify_compile(
                *exprs,
                is_isospin_symmetric_limit=True,
                diagram_type_dict=diagram_type_dict)
        return cexpr
        #
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

@q.timer
def get_cexpr_pi0_current():
    """
    Intend to calculating < J_mu(x_1) pi0(t_1) >
    """
    fn_base = "cache/auto_contract_cexpr/get_cexpr_pi0_current"
    def calc_cexpr():
        diagram_type_dict = dict()
        diagram_type_dict[()] = 'T0'
        diagram_type_dict[((('t_1', 'x_1'), 1), (('x_1', 't_1'), 1))] = 'TypeC'
        #
        jj_d_list = [
            mk_j_mu("x_1", mu)
            for mu in range(4)
        ]
        assert len(jj_d_list) == 4
        #
        pi0d_list = [
                mk_pi_0("t_1") + "pi0(t_1)",
                ]
        assert len(pi0d_list) == 1
        #
        exprs_list_pi0_decay = [
            jj_d * pi0d
            for pi0d in pi0d_list for jj_d in jj_d_list
        ]
        #
        exprs = [
                mk_expr(1) + f"1",
                ]
        exprs += exprs_list_pi0_decay
        #
        cexpr = contract_simplify_compile(
                *exprs,
                is_isospin_symmetric_limit=True,
                diagram_type_dict=diagram_type_dict)
        return cexpr
        #
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

@q.timer_verbose
def auto_contract_tadpole_current(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob):
    fname = q.get_fname()
    fn = f"{job_tag}/tadpole-current/traj-{traj}/tadpole-current.sfield"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_tadpole_current()
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
    sf_tadpole_current = q.SelectedFieldComplexD(fsel, len(expr_names))
    q.set_zero(sf_tadpole_current)
    sf_tadpole_current_arr = sf_tadpole_current[:]
    def load_data():
        for fidx in range(len(xg_fsel_arr)):
            yield fidx
    @q.timer
    def feval(args):
        fidx = args
        xg_snk = tuple(xg_fsel_arr[fidx])
        prob_snk = fsel_prob_arr[fidx]
        val_sum = 0
        for pidx in range(len(xg_psel_arr)):
            xg_src = tuple(xg_psel_arr[pidx])
            prob_src = psel_prob_arr[pidx]
            if xg_snk == xg_src:
                prob_src = prob_src / prob_snk
            prob = prob_src * prob_snk
            pd = {
                    "x_1": ("point-snk", xg_snk,),
                    "x_2": ("point", xg_src,),
                    "size": total_site,
                    }
            val = eval_cexpr(cexpr, positions_dict=pd, get_prop=get_prop)
            val_sum = val_sum + val / prob
        res = (fidx, val_sum,)
        return res
    def sum_function(val_list):
        for idx, res in enumerate(val_list):
            fidx, val_sum, = res
            sf_tadpole_current_arr[fidx] = val_sum
            if (idx + 1) % (len(xg_fsel_arr) // 128 + 16) == 0:
                q.displayln_info(f"{fname}: {idx+1}/{len(xg_fsel_arr)}")
        return None
    auto_contractor_chunk_size = get_param(job_tag, "measurement", "auto_contractor_chunk_size", default=128)
    q.timer_fork(0)
    q.parallel_map_sum(feval, load_data(), sum_function=sum_function, chunksize=auto_contractor_chunk_size)
    q.displayln_info(f"{fname}: timer_display for parallel_map_sum")
    q.timer_display()
    q.timer_merge()
    sf_tadpole_current.save_double(get_save_path(fn))
    q.json_results_append(f"{fname}: sf_tadpole_current sig", q.get_data_sig(sf_tadpole_current, q.RngState()))
    for i, en in enumerate(expr_names):
        q.json_results_append(f"{fname}: sf_tadpole_current '{en}' sig", q.glb_sum(q.get_data_sig(sf_tadpole_current[:, i], q.RngState())))

@q.timer_verbose
def auto_contract_pi0_current(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob):
    fname = q.get_fname()
    fn = f"{job_tag}/pi0-current/traj-{traj}/pi0-current"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_pi0_current()
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
    sf_pi0_current_list = []
    sf_pi0_current_arr_list = []
    for t_src in range(t_size):
        sf = q.SelectedFieldComplexD(fsel, len(expr_names))
        q.set_zero(sf)
        sf_pi0_current_list.append(sf)
        sf_pi0_current_arr_list.append(sf[:])
    def load_data():
        for fidx in range(len(xg_fsel_arr)):
            yield fidx
    @q.timer
    def feval(args):
        fidx = args
        xg_snk = tuple(xg_fsel_arr[fidx])
        prob_snk = fsel_prob_arr[fidx]
        val_list = []
        for t_src in range(t_size):
            prob = prob_snk
            pd = {
                    "x_1": ("point-snk", xg_snk,),
                    "t_1": ("wall", t_src,),
                    "size": total_site,
                    }
            val = eval_cexpr(cexpr, positions_dict=pd, get_prop=get_prop)
            val_list.append(val / prob)
        res = (fidx, val_list,)
        return res
    def sum_function(val_list):
        for idx, res in enumerate(val_list):
            fidx, val_list, = res
            for t_src in range(t_size):
                sf_pi0_current_arr_list[t_src][fidx] = val_list[t_src]
            if (idx + 1) % (len(xg_fsel_arr) // 128 + 16) == 0:
                q.displayln_info(f"{fname}: {idx+1}/{len(xg_fsel_arr)}")
        return None
    auto_contractor_chunk_size = get_param(job_tag, "measurement", "auto_contractor_chunk_size", default=128)
    q.timer_fork(0)
    q.parallel_map_sum(feval, load_data(), sum_function=sum_function, chunksize=auto_contractor_chunk_size)
    q.displayln_info(f"{fname}: timer_display for parallel_map_sum")
    q.timer_display()
    q.timer_merge()
    sfw = q.open_fields(get_save_path(fn + ".acc"), "w", q.Coordinate([ 2, 2, 2, 4, ]))
    for t_src in range(t_size):
        sf = sf_pi0_current_list[t_src]
        tag = f"sf_pi0_current ; t_src={t_src}"
        sf.save_double(sfw, tag)
    sfw.close()
    q.qrename_info(get_save_path(fn + ".acc"), get_save_path(fn))
    sig_arr = np.zeros((t_size, len(expr_names),), dtype=np.complex128)
    for t_src in range(t_size):
        sf = sf_pi0_current_list[t_src]
        q.json_results_append(f"{fname}: sf_pi0_current t_src={t_src} sig", q.get_data_sig(sf, q.RngState()))
        for i, en in enumerate(expr_names):
            sig_arr[t_src, i] = q.glb_sum(q.get_data_sig(sf[:, i], q.RngState(f"t_src={t_src}")))
    for i, en in enumerate(expr_names):
        q.json_results_append(f"{fname}: sf_pi0_current '{en}' sig", sig_arr[:, i].sum())

@q.timer_verbose
def auto_contract_pi0_gg_disc(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract/traj-{traj}/pi0-gg-disc.lat"
    if get_load_path(fn) is not None:
        return
    fn_tadpole_current = get_load_path(f"{job_tag}/tadpole-current/traj-{traj}/tadpole-current.sfield")
    fn_pi0_current = get_load_path(f"{job_tag}/pi0-current/traj-{traj}/pi0-current")
    if fn_tadpole_current is None or fn_pi0_current is None:
        return
    tadpole_current_expr_names = get_expr_names(get_cexpr_tadpole_current())
    pi0_current_expr_names = get_expr_names(get_cexpr_pi0_current())
    assert len(tadpole_current_expr_names) == 1 + 4 * 4
    assert len(pi0_current_expr_names) == 1 + 4
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    t_size = total_site[3]
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
    tsep = get_param(job_tag, "meson_tensor_tsep")
    m_l = get_param(job_tag, "m_l")
    m_h = get_param(job_tag, "m_h")
    geo = q.Geometry(total_site)
    total_volume = geo.total_volume
    r_list = get_r_list(job_tag)
    r_sq_interp_idx_coef_list = get_r_sq_interp_idx_coef_list(job_tag)
    sf_tadpole_current = q.SelectedFieldComplexD(fsel)
    sf_tadpole_current.load_double(fn_tadpole_current)
    assert sf_tadpole_current.multiplicity == len(tadpole_current_expr_names)
    sfr = q.open_fields(fn_pi0_current, "r")
    sf_pi0_current_list = []
    for t_src in range(t_size):
        sf = q.SelectedFieldComplexD(fsel)
        tag = f"sf_pi0_current ; t_src={t_src}"
        sf.load_double(sfr, tag)
        sf_pi0_current_list.append(sf)
        assert sf.multiplicity == len(pi0_current_expr_names)
    sfr.close()
    q.displayln_info(f"{len(sf_pi0_current_list)}")
    expr_names = [
        "< 1 [disc x] [con 0] >",
        "< 1 [disc 0] [con x] >",
    ]
    for tadpole in [ "ls", "sl", "ll", "ss", ]:
        expr_names += [
            f"< e(i,j,k) * x[i] * disc[ j_j(x) ] * j_k(0) * pi0(-tsep) > (tadpole {tadpole})",
            f"< e(i,j,k) * x[i] * j_j(x) * disc[ j_k(0) ] * pi0(-tsep) > (tadpole {tadpole})",
            f"< e(i,j,k) * x[i] * disc[ j_j(x) ] * j_k(0) * pi0(x[t]+tsep) > (tadpole {tadpole})",
            f"< e(i,j,k) * x[i] * j_j(x) * disc[ j_k(0) ] * pi0(x[t]+tsep) > (tadpole {tadpole})",
        ]
    # Some final modification for sf_tadpole_current
    sf_tadpole_current[:, 0] *= 1 / total_volume
    sf_tadpole_current[:, 1:] *= m_h - m_l
    # compute zero distance counts sum
    zero_dis_counts_sum = q.glb_sum((fsel_prob_arr * sf_tadpole_current[:, 0] * sf_pi0_current_list[0][:, 0]).sum())
    # Convert tadpole to field
    f_tadpole_current = q.FieldComplexD(geo, sf_tadpole_current.multiplicity)
    q.set_zero(f_tadpole_current)
    f_tadpole_current @= sf_tadpole_current
    # Convert pi0 current list to field list
    f_pi0_current_list = []
    for t_src in range(t_size):
        sf = sf_pi0_current_list[t_src]
        f = q.FieldComplexD(geo, sf.multiplicity)
        q.set_zero(f)
        f @= sf
        f_pi0_current_list.append(f)
    # make new pi0 current field list based on pi0_current_sep
    xg_arr = geo.xg_arr()
    t_arr = xg_arr[:, 3]
    f_pi0_current_sep_list = []
    for pi0_current_sep in range(t_size):
        f = q.FieldComplexD(geo, sf.multiplicity)
        q.set_zero(f)
        for t_src in range(t_size):
            f_src = f_pi0_current_list[t_src]
            sel = (t_arr - t_src) % t_size == pi0_current_sep
            f[sel] = f_src[sel]
        f_pi0_current_sep_list.append(f)
    # perform all needed convolution
    ff_idx = 0
    idx1 = [ 0, ] # for con 0 (f_pi0_current_sep_list)
    idx2 = [ 0, ] # for disc x (f_tadpole_current)
    info = [ (ff_idx,), ] # for (i, j, k,) tuple
    ff_idx += 1
    for i in range(3):
        for j in range(3):
            if i == j:
                continue
            for k in range(3):
                if k == i or k == j:
                    continue
                for tadpole_idx in range(4):
                    idx1.append(1 + k)
                    idx2.append(1 + tadpole_idx * 4 + j)
                    info.append((ff_idx, i, j, k, tadpole_idx,))
                    ff_idx += 1
    idx1 = np.array(idx1, dtype=np.int32)
    idx2 = np.array(idx2, dtype=np.int32)
    info = np.array(info, dtype=object)
    assert ff_idx == len(idx1)
    assert ff_idx == len(idx2)
    assert ff_idx == len(info)
    ff_list = [
        q.field_convolution(f_pi0_current_sep_list[pi0_current_sep], f_tadpole_current, idx1, idx2)
        for pi0_current_sep in range(t_size)
    ]
    pos_t_expr_name_idx_arr = np.arange(0, len(expr_names), 2, dtype=np.int32)
    neg_t_expr_name_idx_arr = np.arange(1, len(expr_names), 2, dtype=np.int32)
    values = np.zeros((t_size, len(r_list), len(expr_names),), dtype=np.complex128)
    for xg_idx, xg_rel in enumerate(geo.xg_arr()):
        x_rel = q.rel_mod_arr(xg_rel, total_site.to_numpy())
        r_sq = q.get_r_sq(x_rel)
        r_idx_low, r_idx_high, coef_low, coef_high = r_sq_interp_idx_coef_list[r_sq]
        x_rel_t = x_rel[3]
        t = x_rel_t % t_size
        if x_rel_t >= 0:
            pi0_current_sep1 = tsep
            pi0_current_sep2 = -tsep - abs(x_rel_t)
        else:
            pi0_current_sep1 = tsep + abs(x_rel_t)
            pi0_current_sep2 = -tsep
        val = np.zeros(len(expr_names) // 2, dtype=np.complex128)
        val[0] = ff_list[0][xg_idx, 0]
        for ff_idx, i, j, k, tadpole_idx in info[1:]:
            x_i = q.rel_mod_sym(x_rel[i], total_site[i]) # x[i]
            eijk = q.epsilon_tensor(i, j, k) # e(i,j,k)
            jjpi0_t1 = ff_list[pi0_current_sep1][xg_idx, ff_idx] # disc[ j_j(x) ] * j_k(0) * pi0(-tsep)
            jjpi0_t2 = ff_list[pi0_current_sep2][xg_idx, ff_idx] # disc[ j_j(x) ] * j_k(0) * pi0(x[t]+tsep)
            val[1 + tadpole_idx * 2] += eijk * x_i * jjpi0_t1
            val[1 + tadpole_idx * 2 + 1] += eijk * x_i * jjpi0_t2
        if np.all(x_rel == 0):
            val[0] = zero_dis_counts_sum
        if abs(x_rel_t) == t_size // 2:
            continue
        values[t, r_idx_low, pos_t_expr_name_idx_arr] += coef_low * val
        values[t, r_idx_high, pos_t_expr_name_idx_arr] += coef_high * val
        values[-t, r_idx_low, neg_t_expr_name_idx_arr] += coef_low * val
        values[-t, r_idx_high, neg_t_expr_name_idx_arr] += coef_high * val
    res_sum = q.glb_sum(values.transpose(2, 0, 1))
    res_sum *= 1.0 / total_volume
    ld_sum = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t", t_size, [ str(q.rel_mod(t, t_size)) for t in range(t_size) ], ],
        [ "r", len(r_list), [ f"{r:.5f}" for r in r_list ], ],
        ])
    ld_sum.from_numpy(res_sum)
    ld_sum.save(get_save_path(fn))
    q.json_results_append(f"{fname}: ld_sum sig", q.get_data_sig(ld_sum, q.RngState()))
    for i, en in enumerate(expr_names):
        q.json_results_append(f"{fname}: ld_sum '{en}' sig", q.get_data_sig(ld_sum[i], q.RngState()))

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
            (f"{job_tag}/prop-rand-u1-light/traj-{traj}.qar", f"{job_tag}/prop-rand-u1-light/traj-{traj}/geon-info.txt",),
            (f"{job_tag}/prop-rand-u1-strange/traj-{traj}.qar", f"{job_tag}/prop-rand-u1-strange/traj-{traj}/geon-info.txt",),
            (f"{job_tag}/prop-rand-u1-charm/traj-{traj}.qar", f"{job_tag}/prop-rand-u1-charm/traj-{traj}/geon-info.txt",),
            #
            (f"{job_tag}/prop-psrc-light/traj-{traj}.qar", f"{job_tag}/prop-psrc-light/traj-{traj}/geon-info.txt",),
            (f"{job_tag}/psel-prop-psrc-light/traj-{traj}.qar", f"{job_tag}/psel-prop-psrc-light/traj-{traj}/checkpoint.txt",),
            (f"{job_tag}/prop-psrc-strange/traj-{traj}.qar", f"{job_tag}/prop-psrc-strange/traj-{traj}/geon-info.txt",),
            (f"{job_tag}/psel-prop-psrc-strange/traj-{traj}.qar", f"{job_tag}/psel-prop-psrc-strange/traj-{traj}/checkpoint.txt",),
            #
            (f"{job_tag}/prop-smear-light/traj-{traj}.qar", f"{job_tag}/prop-smear-light/traj-{traj}/geon-info.txt",),
            (f"{job_tag}/psel-prop-smear-light/traj-{traj}.qar", f"{job_tag}/psel-prop-smear-light/traj-{traj}/checkpoint.txt",),
            (f"{job_tag}/prop-smear-strange/traj-{traj}.qar", f"{job_tag}/prop-smear-strange/traj-{traj}/geon-info.txt",),
            (f"{job_tag}/psel-prop-smear-strange/traj-{traj}.qar", f"{job_tag}/psel-prop-smear-strange/traj-{traj}/checkpoint.txt",),
            #
            (f"{job_tag}/prop-wsrc-light/traj-{traj}.qar", f"{job_tag}/prop-wsrc-light/traj-{traj}/geon-info.txt",),
            (f"{job_tag}/psel-prop-wsrc-light/traj-{traj}.qar", f"{job_tag}/psel-prop-wsrc-light/traj-{traj}/checkpoint.txt",),
            (f"{job_tag}/prop-wsrc-strange/traj-{traj}.qar", f"{job_tag}/prop-wsrc-strange/traj-{traj}/geon-info.txt",),
            (f"{job_tag}/psel-prop-wsrc-strange/traj-{traj}.qar", f"{job_tag}/psel-prop-wsrc-strange/traj-{traj}/checkpoint.txt",),
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
    get_gf_ape = run_gf_ape(job_tag, get_gf)
    #
    get_wi = run_wi(job_tag, traj)
    #
    get_eig_light = run_eig(job_tag, traj_gf, get_gf)
    get_eig_strange = run_eig_strange(job_tag, traj_gf, get_gf)
    #
    def run_wsrc_full():
        get_eig = get_eig_light
        # run_get_inverter(job_tag, traj, inv_type=0, get_gf=get_gf, get_gt=get_gt, get_eig=get_eig)
        run_prop_wsrc_full(job_tag, traj, inv_type=0, get_gf=get_gf, get_eig=get_eig, get_gt=get_gt, get_wi=get_wi)
        #
        get_eig = get_eig_strange
        # run_get_inverter(job_tag, traj, inv_type=1, get_gf=get_gf, get_gt=get_gt, get_eig=get_eig)
        run_prop_wsrc_full(job_tag, traj, inv_type=1, get_gf=get_gf, get_eig=get_eig, get_gt=get_gt, get_wi=get_wi)
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
    get_fselc = run_fselc(job_tag, traj, get_fsel, get_psel)
    #
    get_eig = get_eig_light
    run_prop_wsrc_sparse(job_tag, traj, inv_type=0, get_gf=get_gf, get_eig=get_eig, get_gt=get_gt, get_psel=get_psel, get_fsel=get_fsel, get_wi=get_wi)
    get_eig = get_eig_strange
    run_prop_wsrc_sparse(job_tag, traj, inv_type=1, get_gf=get_gf, get_eig=get_eig, get_gt=get_gt, get_psel=get_psel, get_fsel=get_fsel, get_wi=get_wi)
    #
    get_psel_smear = run_psel_smear(job_tag, traj)
    #
    def run_with_eig():
        get_eig = get_eig_light
        # run_get_inverter(job_tag, traj, inv_type=0, get_gf=get_gf, get_eig=get_eig)
        # run_prop_wsrc(job_tag, traj, inv_type=0, get_gf=get_gf, get_eig=get_eig, get_gt=get_gt, get_psel=get_psel, get_fsel=get_fselc, get_wi=get_wi)
        run_prop_rand_u1(job_tag, traj, inv_type=0, get_gf=get_gf, get_fsel=get_fsel, get_eig=get_eig)
        run_prop_psrc(job_tag, traj, inv_type=0, get_gf=get_gf, get_eig=get_eig, get_gt=get_gt, get_psel=get_psel, get_fsel=get_fselc, get_f_rand_01=get_f_rand_01)
        run_prop_smear(job_tag, traj, inv_type=0, get_gf=get_gf, get_gf_ape=get_gf_ape, get_eig=get_eig, get_gt=get_gt, get_psel=get_psel, get_fsel=get_fselc, get_psel_smear=get_psel_smear)
        q.clean_cache(q.cache_inv)
    #
    def run_with_eig_strange():
        get_eig = get_eig_strange
        # run_get_inverter(job_tag, traj, inv_type=1, get_gf=get_gf, get_eig=get_eig)
        # run_prop_wsrc(job_tag, traj, inv_type=1, get_gf=get_gf, get_eig=get_eig, get_gt=get_gt, get_psel=get_psel, get_fsel=get_fselc, get_wi=get_wi)
        run_prop_rand_u1(job_tag, traj, inv_type=1, get_gf=get_gf, get_fsel=get_fsel, get_eig=get_eig)
        run_prop_psrc(job_tag, traj, inv_type=1, get_gf=get_gf, get_eig=get_eig, get_gt=get_gt, get_psel=get_psel, get_fsel=get_fselc, get_f_rand_01=get_f_rand_01)
        run_prop_smear(job_tag, traj, inv_type=1, get_gf=get_gf, get_gf_ape=get_gf_ape, get_eig=get_eig, get_gt=get_gt, get_psel=get_psel, get_fsel=get_fselc, get_psel_smear=get_psel_smear)
        q.clean_cache(q.cache_inv)
    #
    def run_charm():
        # run_get_inverter(job_tag, traj, inv_type=2, get_gf=get_gf)
        run_prop_rand_u1(job_tag, traj, inv_type=2, get_gf=get_gf, get_fsel=get_fsel)
        q.clean_cache(q.cache_inv)
    #
    run_with_eig()
    run_with_eig_strange()
    run_charm()
    #
    q.clean_cache()
    if q.obtained_lock_history_list:
        q.timer_display()

@q.timer_verbose
def run_job_contract(job_tag, traj):
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
    get_wi = run_wi(job_tag, traj)
    #
    get_f_weight = run_f_weight_from_wsrc_prop_full(job_tag, traj)
    get_f_rand_01 = run_f_rand_01(job_tag, traj)
    get_fsel_prob = run_fsel_prob(job_tag, traj, get_f_rand_01=get_f_rand_01, get_f_weight=get_f_weight)
    get_psel_prob = run_psel_prob(job_tag, traj, get_f_rand_01=get_f_rand_01, get_f_weight=get_f_weight)
    get_fsel = run_fsel_from_fsel_prob(get_fsel_prob)
    get_psel = run_psel_from_psel_prob(get_psel_prob)
    #
    get_psel_smear = run_psel_smear(job_tag, traj)
    #
    get_get_prop = run_get_prop(job_tag, traj,
            get_gf = get_gf,
            get_gt = get_gt,
            get_psel = get_psel,
            get_fsel = get_fsel,
            get_psel_smear = get_psel_smear,
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
                auto_contract_meson_corr_psnk(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob)
                auto_contract_meson_jwjj2(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob)
                auto_contract_meson_jwjj(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob)
                auto_contract_meson_jj(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob)
                auto_contract_meson_jt(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob)
                auto_contract_meson_m(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob)
                auto_contract_meson_corr(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob)
                auto_contract_meson_corr_psrc(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob)
                auto_contract_meson_corr_psnk_psrc(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob)
                auto_contract_tadpole_current(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob)
                auto_contract_pi0_current(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob)
                auto_contract_pi0_gg_disc(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob)
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
    benchmark_eval_cexpr(get_cexpr_meson_corr())
    benchmark_eval_cexpr(get_cexpr_meson_m())
    benchmark_eval_cexpr(get_cexpr_meson_jt())
    benchmark_eval_cexpr(get_cexpr_meson_jj())
    benchmark_eval_cexpr(get_cexpr_meson_jwjj())
    benchmark_eval_cexpr(get_cexpr_pi0_gg())
    benchmark_eval_cexpr(get_cexpr_tadpole_current())
    benchmark_eval_cexpr(get_cexpr_pi0_current())

### ------

set_param("test-4nt8", "trajs")(list(range(1000, 1001)))
set_param("test-4nt8", "meson_tensor_tsep")(1)
set_param("test-4nt8", "meson_jwjj_threshold")(0.1)
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

set_param("24D", "trajs")([ 2430, 2550, 2590, 2610, 2630, 2940, 2960, ])
set_param("24D", "meson_tensor_tsep")(8)
set_param("24D", "meson_jwjj_threshold")(0.02)
set_param("24D", "measurement", "auto_contractor_chunk_size")(128)

set_param("48I", "trajs")(list(range(1000, 2000, 20)))
set_param("48I", "meson_tensor_tsep")(12)
set_param("48I", "meson_jwjj_threshold")(0.01)
set_param("48I", "measurement", "auto_contractor_chunk_size")(128)

set_param("64I", "trajs")(list(range(1200, 3680, 20)))
set_param("64I", "meson_tensor_tsep")(18)
set_param("64I", "meson_jwjj_threshold")(0.0005)
set_param("64I", "measurement", "auto_contractor_chunk_size")(128)

# ----

job_tag = "test-4nt16-checker"

set_param(job_tag, "trajs")(list(range(1000, 1032)))

set_param(job_tag, "total_site")([ 4, 4, 4, 16, ])
set_param(job_tag, "load_config_params", "twist_boundary_at_boundary")([ 0.0, 0.0, 0.0, -0.5, ])

set_param(job_tag, "mk_sample_gauge_field", "rand_n_step")(2)
set_param(job_tag, "mk_sample_gauge_field", "flow_n_step")(8)
set_param(job_tag, "mk_sample_gauge_field", "hmc_n_traj")(5)

set_param(job_tag, "fermion_params", 0, 0)({ 'Ls': 8, 'M5': 1.8, 'b': 1.5, 'boundary_phases': [1.0, 1.0, 1.0, 1.0], 'c': 0.5, })
for inv_type in [ 1, 2, ]:
    set_param(job_tag, "fermion_params", inv_type, 0)(get_param(job_tag, "fermion_params", 0, 0).copy())
set_param(job_tag, "fermion_params", 0, 0, "mass")(0.01)
set_param(job_tag, "fermion_params", 1, 0, "mass")(0.04)
set_param(job_tag, "fermion_params", 2, 0, "mass")(0.2)
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
    set_param(job_tag, f"cg_params-{inv_type}-2", "maxiter")(500)
    set_param(job_tag, f"cg_params-{inv_type}-2", "maxcycle")(50)

set_param(job_tag, "m_l")(get_param(job_tag, "fermion_params", 0, 0, "mass"))
set_param(job_tag, "m_h")(get_param(job_tag, "fermion_params", 1, 0, "mass"))

set_param(job_tag, "meson_tensor_tsep")(3)

set_param(job_tag, "gf_ape_smear_coef")(0.5)
set_param(job_tag, "gf_ape_smear_step")(30)

set_param(job_tag, "measurement", "auto_contract_meson_corr_wf", "sample_num")(32)
set_param(job_tag, "measurement", "auto_contract_meson_corr_wf", "sample_size")(2)
set_param(job_tag, "measurement", "auto_contract_meson_corr_wf", "t_sep_range")(6)
set_param(job_tag, "measurement", "auto_contract_meson_meson_i0_j0_corr_wf", "sample_num")(32)
set_param(job_tag, "measurement", "auto_contract_meson_meson_i0_j0_corr_wf", "sample_size")(2)
set_param(job_tag, "measurement", "auto_contract_meson_meson_i0_j0_corr_wf", "t_sep_range")(6)

set_param(job_tag, "field-selection-fsel-rate")(20.0)
set_param(job_tag, "field-selection-psel-rate")(10.0)
set_param(job_tag, "field-selection-fsel-psrc-prop-norm-threshold")(0.05)

set_param(job_tag, "prob_exact_wsrc")(0.20)

set_param(job_tag, "prob_acc_1_psrc")(0.25)
set_param(job_tag, "prob_acc_2_psrc")(0.10)

set_param(job_tag, "n_per_tslice_smear")(4)
set_param(job_tag, "prop_smear_coef")(0.9375)
set_param(job_tag, "prop_smear_step")(10)
set_param(job_tag, "prob_acc_1_smear")(0.25)
set_param(job_tag, "prob_acc_2_smear")(0.10)

set_param(job_tag, "n_rand_u1_fsel")(16)
set_param(job_tag, "prob_acc_1_rand_u1")(0.25)
set_param(job_tag, "prob_acc_2_rand_u1")(0.10)

set_param(job_tag, "meson_jwjj_threshold")(2.5)

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

    is_performing_contraction = not q.get_option("--no-contract")

    #######################################################

    job_tags_default = [
            "test-4nt8",
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
                if q.obtained_lock_history_list:
                    q.json_results_append(f"q.obtained_lock_history_list={q.obtained_lock_history_list}")
                    if job_tag[:5] != "test-":
                        gracefully_finish()
        for traj in get_param(job_tag, "trajs"):
            if is_performing_contraction:
                q.check_time_limit()
                run_job_contract(job_tag, traj)
                if q.obtained_lock_history_list:
                    q.json_results_append(f"q.obtained_lock_history_list={q.obtained_lock_history_list}")
                    if job_tag[:5] != "test-":
                        gracefully_finish()

    q.check_log_json(__file__, check_eps=5e-5)

    gracefully_finish()

# ----
