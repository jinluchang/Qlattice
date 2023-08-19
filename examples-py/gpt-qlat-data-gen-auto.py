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
        ]

# ----

def mk_jw_v_mu(p, mu):
    return mk_jpi_mu(p, mu) + mk_jk_mu(p, mu)

def mk_jw_a_mu(p, mu):
    return mk_j5pi_mu(p, mu) + mk_j5k_mu(p, mu)

def mk_sw5(p):
    return mk_pi_p(p, is_dagger = True) + mk_k_p(p, is_dagger = True)

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
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True, diagram_type_dict = diagram_type_dict)
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

@q.timer_verbose
def auto_contract_meson_corr(job_tag, traj, get_get_prop, get_psel, get_fsel):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_corr.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_corr()
    expr_names = get_cexpr_names(cexpr)
    total_site = dict_params[job_tag]["total_site"]
    t_size = total_site[3]
    get_prop = get_get_prop()
    psel = get_psel()
    fsel, fselc = get_fsel()
    xg_fsel_list = np.array(fsel.to_psel_local().to_list())
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
        val = eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop)
        return val, t
    def sum_function(val_list):
        values = np.zeros((total_site[3], len(expr_names),), dtype = complex)
        for val, t in val_list:
            values[t] += val
        return q.glb_sum(values.transpose(1, 0))
    auto_contractor_chunk_size = get_param(job_tag, "measurement", "auto_contractor_chunk_size", default = 128)
    q.timer_fork(0)
    res_sum = q.parallel_map_sum(feval, load_data(), sum_function = sum_function, chunksize = auto_contractor_chunk_size)
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
    q.displayln_info(f"CHECK: {fname}: ld sig: {q.get_double_sig(ld, q.RngState()):.5E}")

@q.timer_verbose
def auto_contract_meson_corr_psnk(job_tag, traj, get_get_prop, get_psel, get_fsel):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_corr_psnk.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_corr()
    expr_names = get_cexpr_names(cexpr)
    total_site = q.Coordinate(dict_params[job_tag]["total_site"])
    t_size = total_site[3]
    get_prop = get_get_prop()
    fsel, fselc = get_fsel()
    xg_fsel_list = np.array(fsel.to_psel_local().to_list())
    geo = q.Geometry(total_site, 1)
    total_volume = geo.total_volume()
    def load_data():
        for t_src in range(total_site[3]):
            for xg_snk in xg_fsel_list:
                yield t_src, xg_snk
    @q.timer
    def feval(args):
        t_src, xg_snk = args
        xg_snk = tuple(xg_snk.tolist())
        t = (xg_snk[3] - t_src) % total_site[3]
        pd = {
                "x_2" : ("point-snk", xg_snk,),
                "x_1" : ("wall", t_src,),
                "size" : total_site,
                }
        val = eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop)
        return val, t
    def sum_function(val_list):
        values = np.zeros((total_site[3], len(expr_names),), dtype = complex)
        for val, t in val_list:
            values[t] += val
        return values.transpose(1, 0)
    auto_contractor_chunk_size = get_param(job_tag, "measurement", "auto_contractor_chunk_size", default = 128)
    q.timer_fork(0)
    res_sum = q.glb_sum(
            q.parallel_map_sum(feval, load_data(), sum_function = sum_function, chunksize = auto_contractor_chunk_size))
    q.displayln_info("timer_display for auto_contract_meson_corr_psnk")
    q.timer_display()
    q.timer_merge()
    res_sum *= 1.0 / (total_volume * fsel.prob())
    assert q.qnorm(res_sum[0] - 1.0) < 1e-10
    ld = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t_sep", t_size, [ str(q.rel_mod(t, t_size)) for t in range(t_size) ], ],
        ])
    ld.from_numpy(res_sum)
    ld.save(get_save_path(fn))
    sig_msg_list = [
            f"CHECK: {fname}: ld sig: {q.get_double_sig(ld, q.RngState()):.5E}",
            ]
    for msg in sig_msg_list:
        q.displayln_info(msg)

@q.timer_verbose
def auto_contract_meson_corr_psrc(job_tag, traj, get_get_prop, get_psel, get_fsel):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_corr_psrc.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_corr()
    expr_names = get_cexpr_names(cexpr)
    total_site = q.Coordinate(dict_params[job_tag]["total_site"])
    t_size = total_site[3]
    get_prop = get_get_prop()
    psel = get_psel()
    fsel, fselc = get_fsel()
    xg_fsel_list = np.array(fsel.to_psel_local().to_list())
    xg_psel_list = np.array(psel.to_list())
    geo = q.Geometry(total_site, 1)
    total_volume = geo.total_volume()
    def load_data():
        x_t_list = get_mpi_chunk(
                [ (tuple(xg_src.tolist()), t_snk,) for t_snk in range(total_site[3]) for xg_src in xg_psel_list ],
                rng_state = None)
        for xg_src, t_snk in x_t_list:
            yield xg_src, t_snk
    @q.timer
    def feval(args):
        xg_src, t_snk = args
        t = (xg_src[3] - t_snk) % total_site[3]
        pd = {
                "x_2" : ("point", xg_src,),
                "x_1" : ("wall", t_snk,),
                "size" : total_site,
                }
        val = eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop)
        return val, t
    def sum_function(val_list):
        values = np.zeros((total_site[3], len(expr_names),), dtype = complex)
        for val, t in val_list:
            values[t] += val
        return values.transpose(1, 0)
    auto_contractor_chunk_size = get_param(job_tag, "measurement", "auto_contractor_chunk_size", default = 128)
    q.timer_fork(0)
    res_sum = q.glb_sum(
            q.parallel_map_sum(feval, load_data(), sum_function = sum_function, chunksize = auto_contractor_chunk_size))
    q.displayln_info("timer_display for auto_contract_meson_corr_psrc")
    q.timer_display()
    q.timer_merge()
    res_sum *= 1.0 / len(xg_psel_list)
    assert q.qnorm(res_sum[0] - 1.0) < 1e-10
    ld = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t_sep", t_size, [ str(q.rel_mod(t, t_size)) for t in range(t_size) ], ],
        ])
    ld.from_numpy(res_sum)
    ld.save(get_save_path(fn))
    sig_msg_list = [
            f"CHECK: {fname}: ld sig: {q.get_double_sig(ld, q.RngState()):.5E}",
            ]
    for msg in sig_msg_list:
        q.displayln_info(msg)

@q.timer_verbose
def auto_contract_meson_corr_psnk_psrc(job_tag, traj, get_get_prop, get_psel, get_fsel):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_corr_psnk_psrc.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_corr()
    expr_names = get_cexpr_names(cexpr)
    total_site = q.Coordinate(dict_params[job_tag]["total_site"])
    t_size = total_site[3]
    get_prop = get_get_prop()
    psel = get_psel()
    fsel, fselc = get_fsel()
    xg_fsel_list = np.array(fsel.to_psel_local().to_list())
    xg_psel_list = np.array(psel.to_list())
    geo = q.Geometry(total_site, 1)
    total_volume = geo.total_volume()
    r_list = get_r_list(job_tag)
    r_sq_interp_idx_coef_list = get_r_sq_interp_idx_coef_list(job_tag)
    def load_data():
        for xg_src in xg_psel_list:
            yield xg_src
    @q.timer
    def feval(args):
        xg_src = args
        xg_src = tuple(xg_src.tolist())
        res_list = []
        for xg_snk in xg_fsel_list:
            xg_snk = tuple(xg_snk.tolist())
            x_rel = [ q.rel_mod(xg_snk[mu] - xg_src[mu], total_site[mu]) for mu in range(4) ]
            r_sq = q.get_r_sq(x_rel)
            r_idx_low, r_idx_high, coef_low, coef_high = r_sq_interp_idx_coef_list[r_sq]
            t = (xg_snk[3] - xg_src[3]) % total_site[3]
            pd = {
                    "x_2": ("point-snk", xg_snk,),
                    "x_1": ("point", xg_src,),
                    "size": total_site,
                    }
            val = eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop)
            res_list.append((val, t, r_idx_low, r_idx_high, coef_low, coef_high))
        return res_list
    def sum_function(val_list):
        values = np.zeros((total_site[3], len(r_list), len(expr_names),), dtype = complex)
        for idx, res_list in enumerate(val_list):
            for val, t, r_idx_low, r_idx_high, coef_low, coef_high in res_list:
                values[t, r_idx_low] += coef_low * val
                values[t, r_idx_high] += coef_high * val
            q.displayln_info(f"{fname}: {idx+1}/{len(xg_psel_list)}")
        return values.transpose(2, 0, 1)
    q.timer_fork(0)
    res_sum = q.glb_sum(
            q.parallel_map_sum(feval, load_data(), sum_function = sum_function, chunksize = 1))
    q.displayln_info("timer_display for auto_contract_meson_corr_psnk_psrc")
    q.timer_display()
    q.timer_merge()
    res_sum *= 1.0 / (len(xg_psel_list) * total_volume * fsel.prob() / total_site[3])
    assert q.qnorm(res_sum[0].sum(1) - 1.0) < 1e-10
    ld = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t_sep", t_size, [ str(q.rel_mod(t, t_size)) for t in range(t_size) ], ],
        [ "r", len(r_list), [ f"{r:.5f}" for r in r_list ], ],
        ])
    ld.from_numpy(res_sum)
    ld.save(get_save_path(fn))
    q.displayln_info(f"CHECK: {fname}: ld sig: {q.get_double_sig(ld, q.RngState()):.5E}")

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
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True, diagram_type_dict = diagram_type_dict)
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

@q.timer_verbose
def auto_contract_meson_jt(job_tag, traj, get_get_prop, get_psel, get_fsel):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_jt.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_jt()
    expr_names = get_cexpr_names(cexpr)
    total_site = q.Coordinate(dict_params[job_tag]["total_site"])
    t_size = total_site[3]
    get_prop = get_get_prop()
    psel = get_psel()
    fsel, fselc = get_fsel()
    xg_fsel_list = np.array(fsel.to_psel_local().to_list())
    xg_psel_list = np.array(psel.to_list())
    tsep = dict_params[job_tag]["meson_tensor_tsep"]
    geo = q.Geometry(total_site, 1)
    total_volume = geo.total_volume()
    def load_data():
        for xg_snk in xg_fsel_list:
            yield xg_snk
    @q.timer
    def feval(args):
        xg_snk = args
        xg_snk = tuple(xg_snk.tolist())
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
        val = eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop)
        return val
    def sum_function(val_list):
        values = np.zeros(len(expr_names), dtype = complex)
        for val in val_list:
            values += val
        return values
    auto_contractor_chunk_size = get_param(job_tag, "measurement", "auto_contractor_chunk_size", default = 128)
    q.timer_fork(0)
    res_sum = q.glb_sum(
            q.parallel_map_sum(feval, load_data(), sum_function = sum_function, chunksize = auto_contractor_chunk_size))
    q.displayln_info("timer_display for auto_contract_meson_jt")
    q.timer_display()
    q.timer_merge()
    res_sum *= 1.0 / (total_volume * fsel.prob())
    assert q.qnorm(res_sum[0] - 1.0) < 1e-10
    ld_sum = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        ])
    ld_sum.from_numpy(res_sum)
    ld_sum.save(get_save_path(fn))
    sig_msg_list = [
            f"CHECK: {fname}: ld_sum sig: {q.get_double_sig(ld_sum, q.RngState()):.5E}",
            ]
    for msg in sig_msg_list:
        q.displayln_info(msg)

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
                is_isospin_symmetric_limit = True,
                diagram_type_dict = diagram_type_dict)
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

@q.timer_verbose
def auto_contract_meson_m(job_tag, traj, get_get_prop, get_psel, get_fsel):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_m.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_m()
    expr_names = get_cexpr_names(cexpr)
    total_site = q.Coordinate(dict_params[job_tag]["total_site"])
    t_size = total_site[3]
    get_prop = get_get_prop()
    psel = get_psel()
    fsel, fselc = get_fsel()
    xg_fsel_list = np.array(fsel.to_psel_local().to_list())
    xg_psel_list = np.array(psel.to_list())
    tsep = dict_params[job_tag]["meson_tensor_tsep"]
    geo = q.Geometry(total_site, 1)
    total_volume = geo.total_volume()
    def load_data():
        for xg_snk in xg_fsel_list:
            yield xg_snk
    @q.timer
    def feval(args):
        xg_snk = args
        xg_snk = tuple(xg_snk.tolist())
        t = xg_snk[3]
        t_2 = (t + tsep) % total_site[3]
        t_1 = (t - tsep) % total_site[3]
        pd = {
                "x" : ("point-snk", xg_snk,),
                "t_1" : ("wall", t_1,),
                "t_2" : ("wall", t_2,),
                "size" : total_site,
                }
        val = eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop)
        return val
    def sum_function(val_list):
        values = np.zeros(len(expr_names), dtype = complex)
        for val in val_list:
            values += val
        return values
    auto_contractor_chunk_size = get_param(job_tag, "measurement", "auto_contractor_chunk_size", default = 128)
    q.timer_fork(0)
    res_sum = q.glb_sum(
            q.parallel_map_sum(feval, load_data(), sum_function = sum_function, chunksize = auto_contractor_chunk_size))
    q.displayln_info("timer_display for auto_contract_meson_m")
    q.timer_display()
    q.timer_merge()
    res_sum *= 1.0 / (total_volume * fsel.prob())
    assert q.qnorm(res_sum[0] - 1.0) < 1e-10
    ld_sum = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        ])
    ld_sum.from_numpy(res_sum)
    ld_sum.save(get_save_path(fn))
    sig_msg_list = [
            f"CHECK: {fname}: ld_sum sig: {q.get_double_sig(ld_sum, q.RngState()):.5E}",
            ]
    for msg in sig_msg_list:
        q.displayln_info(msg)

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
                is_isospin_symmetric_limit = True,
                diagram_type_dict = diagram_type_dict)
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

@q.timer_verbose
def auto_contract_meson_jj(job_tag, traj, get_get_prop, get_psel, get_fsel):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_jj.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_jj()
    expr_names = get_cexpr_names(cexpr)
    total_site = q.Coordinate(dict_params[job_tag]["total_site"])
    t_size = total_site[3]
    get_prop = get_get_prop()
    psel = get_psel()
    fsel, fselc = get_fsel()
    xg_fsel_list = np.asarray(fsel.to_psel_local())
    xg_psel_list = np.asarray(psel)
    tsep = dict_params[job_tag]["meson_tensor_tsep"]
    geo = q.Geometry(total_site, 1)
    r_list = get_r_list(job_tag)
    r_sq_interp_idx_coef_list = get_r_sq_interp_idx_coef_list(job_tag)
    def load_data():
        for idx, xg_src in enumerate(xg_psel_list):
            xg_src = tuple(xg_src.tolist())
            yield xg_src
    @q.timer
    def feval(args):
        xg_src = args
        res_list = []
        for xg_snk in xg_fsel_list:
            xg_snk = tuple(xg_snk.tolist())
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
            val = eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop)
            res_list.append((val, t, r_idx_low, r_idx_high, coef_low, coef_high,))
        return res_list
    def sum_function(val_list):
        values = np.zeros((t_size, len(r_list), len(expr_names),), dtype = complex)
        for idx, res_list in enumerate(val_list):
            for val, t, r_idx_low, r_idx_high, coef_low, coef_high in res_list:
                values[t, r_idx_low] += coef_low * val
                values[t, r_idx_high] += coef_high * val
            q.displayln_info(f"{fname}: {idx+1}/{len(xg_psel_list)}")
        return q.glb_sum(values.transpose(2, 0, 1))
    q.timer_fork(0)
    res_sum = q.parallel_map_sum(feval, load_data(), sum_function = sum_function, chunksize = 1)
    q.displayln_info(f"{fname}: timer_display for parallel_map_sum")
    q.timer_display()
    q.timer_merge()
    res_sum *= 1.0 / (len(xg_psel_list) * fsel.prob())
    ld_sum = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t", t_size, [ str(q.rel_mod(t, t_size)) for t in range(t_size) ], ],
        [ "r", len(r_list), [ f"{r:.5f}" for r in r_list ], ],
        ])
    ld_sum.from_numpy(res_sum)
    ld_sum.save(get_save_path(fn))
    q.displayln_info(f"CHECK: {fname}: ld_sum sig: {q.get_double_sig(ld_sum, q.RngState()):.5E}")

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
                is_isospin_symmetric_limit = True,
                diagram_type_dict = diagram_type_dict)
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

@q.timer_verbose
def auto_contract_meson_jwjj(job_tag, traj, get_get_prop, get_psel, get_fsel):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_jwjj.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_jwjj()
    expr_names = get_cexpr_names(cexpr)
    total_site = q.Coordinate(dict_params[job_tag]["total_site"])
    t_size = total_site[3]
    get_prop = get_get_prop()
    psel = get_psel()
    fsel, fselc = get_fsel()
    xg_fsel_list = np.array(fsel.to_psel_local().to_list())
    xg_psel_list = np.array(psel.to_list())
    tsep = dict_params[job_tag]["meson_tensor_tsep"]
    geo = q.Geometry(total_site, 1)
    r_list = get_r_list(job_tag)
    r_sq_interp_idx_coef_list = get_r_sq_interp_idx_coef_list(job_tag)
    n_elems = len(xg_fsel_list)
    n_points = len(xg_psel_list)
    n_pairs = n_points * (n_points - 1) // 2 + n_points
    #
    threshold = dict_params[job_tag]["meson_jwjj_threshold"]
    u_rand_prob = q.SelectedField(q.ElemTypeDouble, fsel, 1)
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
        return abs(ama_extract(get_prop(*args, is_norm_sqrt = True), is_sloppy = is_sloppy))
    def load_psrc_psrc_prop_norm_sqrt(flavor, i):
        xg1_src = tuple(xg_psel_list[i])
        x_1 = ("point", xg1_src,)
        v_list = []
        for j in range(n_points):
            xg2_src = tuple(xg_psel_list[j])
            x_2 = ("point", xg2_src,)
            v = get_prop_norm_sqrt(flavor, x_1, x_2)
            v_list.append(v)
        return np.array(v_list)
    psrc_psrc_prop_norm_sqrt = np.array([
        q.parallel_map(lambda i: load_psrc_psrc_prop_norm_sqrt(flavor, i), range(n_points))
        for flavor in [ "l", "s", ]
        ], dtype = float)
    def load_wsrc_psrc_prop_norm_sqrt(flavor, t):
        ts = ("wall", t,)
        v_list = []
        for j in range(n_points):
            xg2_src = tuple(xg_psel_list[j])
            x_2 = ("point", xg2_src,)
            v = get_prop_norm_sqrt(flavor, x_2, ts)
            v_list.append(v)
        return np.array(v_list)
    wsrc_psrc_prop_norm_sqrt = np.array([
        q.parallel_map(lambda t: load_wsrc_psrc_prop_norm_sqrt(flavor, t), range(t_size))
        for flavor in [ "l", "s", ]
        ], dtype = float)
    def load_wsrc_psnk_prop_norm_sqrt(flavor, t):
        ts = ("wall", t,)
        v_list = []
        for j in range(n_elems):
            xg2_snk = tuple(xg_fsel_list[j])
            x_2 = ("point-snk", xg2_snk,)
            v = get_prop_norm_sqrt(flavor, x_2, ts)
            v_list.append(v)
        return np.array(v_list)
    wsrc_psnk_prop_norm_sqrt = np.array([
        q.parallel_map(lambda t: load_wsrc_psnk_prop_norm_sqrt(flavor, t), range(t_size))
        for flavor in [ "l", "s", ]
        ], dtype = float)
    def load_psrc_psnk_prop_norm_sqrt(flavor, i):
        xg1_src = tuple(xg_psel_list[i])
        x_1 = ("point", xg1_src,)
        v_list = []
        for j in range(n_elems):
            xg2_snk = tuple(xg_fsel_list[j])
            x_2 = ("point-snk", xg2_snk,)
            v = get_prop_norm_sqrt(flavor, x_2, x_1)
            v_list.append(v)
        return np.array(v_list)
    psrc_psnk_prop_norm_sqrt = np.array([
        q.parallel_map(lambda i: load_psrc_psnk_prop_norm_sqrt(flavor, i), range(n_points))
        for flavor in [ "l", "s", ]
        ], dtype = float)
    def get_estimate(idx_snk, idx1, idx2, t_1, t_2):
        flavor_l = 0
        flavor_s = 1
        xg_snk_t = xg_fsel_list[idx_snk, 3]
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
        assert np.all(value > 0)
        return value
    @q.timer
    def get_weight(idx_snk, idx1, idx2, t_1, t_2):
        # return weight for this point (1 / prob or zero)
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
        for idx1, xg1_src in enumerate(xg_psel_list):
            xg1_src = tuple(xg1_src.tolist())
            for idx2, xg2_src in enumerate(xg_psel_list):
                xg2_src = tuple(xg2_src.tolist())
                if idx2 > idx1:
                    continue
                idx_pair += 1
                yield idx1, idx2
    @q.timer
    def feval(args):
        idx1, idx2 = args
        xg1_src = tuple(xg_psel_list[idx1])
        xg2_src = tuple(xg_psel_list[idx2])
        xg1_src_t = xg1_src[3]
        xg2_src_t = xg2_src[3]
        x_rel = [ q.rel_mod(xg2_src[mu] - xg1_src[mu], total_site[mu]) for mu in range(4) ]
        r_sq = q.get_r_sq(x_rel)
        idx_snk_arr = np.arange(n_elems)
        xg_t_arr = xg_fsel_list[idx_snk_arr, 3]
        t_size_arr = np.broadcast_to(t_size, xg_t_arr.shape)
        xg1_xg_t_arr = q.rel_mod_arr(xg1_src_t - xg_t_arr, t_size_arr)
        xg2_xg_t_arr = q.rel_mod_arr(xg2_src_t - xg_t_arr, t_size_arr)
        t_1_arr = (np.minimum(0, np.minimum(xg1_xg_t_arr, xg2_xg_t_arr)) + xg_t_arr - tsep) % t_size
        t_2_arr = (np.maximum(0, np.maximum(xg1_xg_t_arr, xg2_xg_t_arr)) + xg_t_arr + tsep) % t_size
        weight_arr = get_weight(idx_snk_arr, idx1, idx2, t_1_arr, t_2_arr)
        weight_arr[np.abs(xg2_xg_t_arr - xg1_xg_t_arr) >= t_size_arr // 2] = 0.0
        results = []
        for idx_snk in idx_snk_arr[weight_arr > 0]:
            weight = weight_arr[idx_snk]
            xg_snk = tuple(xg_fsel_list[idx_snk])
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
                    "size" : total_site.to_list(),
                    }
            t1 = xg1_xg_t
            t2 = xg2_xg_t
            val = eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop)
            r_idx_low, r_idx_high, coef_low, coef_high = r_sq_interp_idx_coef_list[r_sq]
            results.append((weight, val, t1, t2, r_idx_low, r_idx_high, coef_low, coef_high,))
        return idx1, idx2, results
    def sum_function(val_list):
        n_total = 0
        n_selected = 0
        idx_pair = 0
        values = np.zeros((t_size, t_size, len(r_list), len(expr_names),), dtype = complex)
        for idx1, idx2, results in val_list:
            idx_pair += 1
            n_total += n_elems
            xg1_src = tuple(xg_psel_list[idx1])
            xg2_src = tuple(xg_psel_list[idx2])
            for weight, val, t1, t2, r_idx_low, r_idx_high, coef_low, coef_high in results:
                n_selected += 1
                values[t1, t2, r_idx_low] += coef_low * weight * val
                values[t1, t2, r_idx_high] += coef_high * weight * val
            if idx_pair % (n_pairs // 1000 + 100) == 0:
                q.displayln_info(1, f"{fname}: {idx_pair}/{n_pairs} {xg1_src} {xg2_src} {len(results)}/{n_elems} n_total={n_total} n_selected={n_selected} ratio={n_selected/n_total}")
        q.displayln_info(1, f"{fname}: Final: n_total={n_total} n_selected={n_selected} ratio={n_selected/n_total}")
        return values.transpose(3, 0, 1, 2)
    q.timer_fork(0)
    res_sum = q.glb_sum(
            q.parallel_map_sum(feval, load_data(), sum_function = sum_function, chunksize = 1))
    q.displayln_info("{fname}: timer_display")
    q.timer_display()
    q.timer_merge()
    total_volume = geo.total_volume()
    res_sum *= 1.0 / (total_volume * fsel.prob() / t_size)
    ld_sum = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t1", t_size, [ str(q.rel_mod(t, t_size)) for t in range(t_size) ], ],
        [ "t2", t_size, [ str(q.rel_mod(t, t_size)) for t in range(t_size) ], ],
        [ "r", len(r_list), [ f"{r:.5f}" for r in r_list ], ],
        ])
    ld_sum.from_numpy(res_sum)
    ld_sum.save(get_save_path(fn))
    sig_msg_list = [
            f"CHECK: {fname}: ld_sum sig: {q.get_double_sig(ld_sum, q.RngState()):.5E}",
            ]
    for msg in sig_msg_list:
        q.displayln_info(msg)

### ------

@q.timer_verbose
def run_job(job_tag, traj):
    fns_produce = [
            f"{job_tag}/auto-contract/traj-{traj}/checkpoint.txt",
            #
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
    traj_gf = traj
    if job_tag[:5] == "test-":
        # ADJUST ME
        traj_gf = 1000
        #
    #
    get_gf = run_gf(job_tag, traj_gf)
    get_gt = run_gt(job_tag, traj_gf, get_gf)
    get_gf_ape = run_gf_ape(job_tag, get_gf)
    #
    get_psel = run_psel(job_tag, traj)
    get_fsel = run_fsel(job_tag, traj, get_psel)
    #
    get_wi = run_wi(job_tag, traj)
    get_psel_smear = run_psel_smear(job_tag, traj)
    #
    def run_with_eig():
        get_eig = run_eig(job_tag, traj_gf, get_gf)
        run_get_inverter(job_tag, traj, inv_type = 0, get_gf = get_gf, get_eig = get_eig)
        run_get_inverter(job_tag, traj, inv_type = 0, get_gf = get_gf, get_gt = get_gt, get_eig = get_eig)
        run_prop_wsrc(job_tag, traj, inv_type = 0, get_gf = get_gf, get_eig = get_eig, get_gt = get_gt, get_psel = get_psel, get_fsel = get_fsel, get_wi = get_wi)
        run_prop_rand_u1(job_tag, traj, inv_type = 0, get_gf = get_gf, get_fsel = get_fsel, get_eig = get_eig)
        run_prop_psrc(job_tag, traj, inv_type = 0, get_gf = get_gf, get_eig = get_eig, get_gt = get_gt, get_psel = get_psel, get_fsel = get_fsel)
        run_prop_smear(job_tag, traj, inv_type = 0, get_gf = get_gf, get_gf_ape = get_gf_ape, get_eig = get_eig, get_gt = get_gt, get_psel = get_psel, get_fsel = get_fsel, get_psel_smear = get_psel_smear)
        q.clean_cache(q.cache_inv)
    #
    def run_with_eig_strange():
        get_eig = run_eig_strange(job_tag, traj_gf, get_gf)
        run_get_inverter(job_tag, traj, inv_type = 1, get_gf = get_gf, get_eig = get_eig)
        run_get_inverter(job_tag, traj, inv_type = 1, get_gf = get_gf, get_gt = get_gt, get_eig = get_eig)
        run_prop_wsrc(job_tag, traj, inv_type = 1, get_gf = get_gf, get_eig = get_eig, get_gt = get_gt, get_psel = get_psel, get_fsel = get_fsel, get_wi = get_wi)
        run_prop_rand_u1(job_tag, traj, inv_type = 1, get_gf = get_gf, get_fsel = get_fsel, get_eig = get_eig)
        run_prop_psrc(job_tag, traj, inv_type = 1, get_gf = get_gf, get_eig = get_eig, get_gt = get_gt, get_psel = get_psel, get_fsel = get_fsel)
        run_prop_smear(job_tag, traj, inv_type = 1, get_gf = get_gf, get_gf_ape = get_gf_ape, get_eig = get_eig, get_gt = get_gt, get_psel = get_psel, get_fsel = get_fsel, get_psel_smear = get_psel_smear)
        q.clean_cache(q.cache_inv)
    #
    def run_charm():
        run_get_inverter(job_tag, traj, inv_type = 2, get_gf = get_gf)
        run_prop_rand_u1(job_tag, traj, inv_type = 2, get_gf = get_gf, get_fsel = get_fsel)
        q.clean_cache(q.cache_inv)
    #
    run_with_eig()
    run_with_eig_strange()
    run_charm()
    #
    get_get_prop = run_get_prop(job_tag, traj,
            get_gf = get_gf,
            get_gt = get_gt,
            get_psel = get_psel,
            get_fsel = get_fsel,
            get_psel_smear = get_psel_smear,
            get_wi = get_wi,
            prop_types = [
                "wsrc psel s",
                "wsrc psel l",
                "wsrc fsel s",
                "wsrc fsel l",
                "psrc psel s",
                "psrc psel l",
                "psrc fsel s",
                "psrc fsel l",
                "rand_u1 fsel c",
                "rand_u1 fsel s",
                "rand_u1 fsel l",
                ],
            )
    #
    run_r_list(job_tag)
    #
    fn_checkpoint = f"{job_tag}/auto-contract/traj-{traj}/checkpoint.txt"
    if get_load_path(fn_checkpoint) is None:
        if q.obtain_lock(f"locks/{job_tag}-{traj}-auto-contract"):
            q.timer_fork()
            # ADJUST ME
            auto_contract_meson_corr_psnk(job_tag, traj, get_get_prop, get_psel, get_fsel)
            auto_contract_meson_jwjj(job_tag, traj, get_get_prop, get_psel, get_fsel)
            auto_contract_meson_jj(job_tag, traj, get_get_prop, get_psel, get_fsel)
            auto_contract_meson_jt(job_tag, traj, get_get_prop, get_psel, get_fsel)
            auto_contract_meson_m(job_tag, traj, get_get_prop, get_psel, get_fsel)
            auto_contract_meson_corr(job_tag, traj, get_get_prop, get_psel, get_fsel)
            auto_contract_meson_corr_psrc(job_tag, traj, get_get_prop, get_psel, get_fsel)
            auto_contract_meson_corr_psnk_psrc(job_tag, traj, get_get_prop, get_psel, get_fsel)
            #
            q.qtouch_info(get_save_path(fn_checkpoint))
            q.release_lock()
            q.displayln_info("timer_display for runjob")
            q.timer_display()
            q.timer_merge()
            #
            q.clean_cache()

def get_all_cexpr():
    benchmark_eval_cexpr(get_cexpr_meson_corr())
    benchmark_eval_cexpr(get_cexpr_meson_m())
    benchmark_eval_cexpr(get_cexpr_meson_jt())
    benchmark_eval_cexpr(get_cexpr_meson_jj())
    benchmark_eval_cexpr(get_cexpr_meson_jwjj())

size_node_list = [
        [1, 1, 1, 1],
        [1, 1, 2, 1],
        [1, 2, 2, 1],
        [2, 2, 2, 1],
        [2, 2, 4, 1],
        [2, 4, 4, 1],
        [4, 4, 4, 1],
        [4, 4, 8, 1],
        [4, 8, 8, 1],
        [8, 8, 8, 1],
        ]

set_param("test-4nt8", "mk_sample_gauge_field", "rand_n_step", value = 2)
set_param("test-4nt8", "mk_sample_gauge_field", "flow_n_step", value = 8)
set_param("test-4nt8", "mk_sample_gauge_field", "hmc_n_traj", value = 1)
set_param("test-4nt8", "trajs", value = [ 1000, ])

qg.begin_with_gpt()

q.qremove_all_info("results")

job_tags = [
        "test-4nt8",
        ]

q.check_time_limit()

get_all_cexpr()

for job_tag in job_tags:
    q.displayln_info(pprint.pformat(get_param(job_tag)))
    q.displayln_info("CHECK: ", get_param(job_tag))
    for traj in get_param(job_tag, "trajs"):
        run_job(job_tag, traj)

q.timer_display()

q.displayln_info("CHECK: finished successfully.")

qg.end_with_gpt()
