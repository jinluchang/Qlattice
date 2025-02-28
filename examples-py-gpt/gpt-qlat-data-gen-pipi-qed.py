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
                mk_j5pi_mu("x_2", 3)    * mk_pi_p("x_1")             + f"j5pi_t(0) * pi+(-tsep)",
                mk_pi_p("x_2", True)    * mk_j5pi_mu("x_1", 3, True) + f"pi+^dag(0) * j5pi_t^dag(-tsep)",
                mk_j5pi_mu("x_2", 3)    * mk_j5pi_mu("x_1", 3, True) + f"j5pi_t(0) * j5pi_t^dag(-tsep)",
                mk_k_p("x_2", True)     * mk_k_p("x_1")              + f"K+^dag(0) * K+(-tsep)",
                mk_j5k_mu("x_2", 3)     * mk_k_p("x_1")              + f"j5k_t(0) * K+(-tsep)",
                mk_k_p("x_2", True)     * mk_j5k_mu("x_1", 3, True)  + f"K+^dag(0) * j5k_t^dag(-tsep)",
                mk_j5k_mu("x_2", 3)     * mk_j5k_mu("x_1", 3, True)  + f"j5k_t(0) * j5k_t^dag(-tsep)",
                #
                mk_a0_p("x_2", True)    * mk_a0_p("x_1")             + f"a0+^dag(0) * a0+(-tsep)",
                mk_kappa_p("x_2", True) * mk_kappa_p("x_1")          + f"kappa+^dag(0) * kappa+(-tsep)",
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
    json_results_append(f"{fname}: ld sig", q.get_data_sig(ld, q.RngState()))
    for i, en in enumerate(expr_names):
        json_results_append(f"{fname}: ld '{en}' sig", q.get_data_sig(ld[i], q.RngState()))

# ----

@q.timer
def get_cexpr_pipi_corr():
    fn_base = "cache/auto_contract_cexpr/get_cexpr_pipi_corr"
    def calc_cexpr():
        diagram_type_dict = dict()
        exprs = [
                mk_fac(1) + f"1",
                mk_pipi_i22("snk_1", "snk_2", True) * mk_pipi_i22("src_1", "src_2") + f"pipi_i22+^dag(0) * pipi_i22(-tsep)",
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit=True, diagram_type_dict=diagram_type_dict)
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

@q.timer_verbose
def auto_contract_pipi_corr(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract/traj-{traj}/pipi_corr.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_pipi_corr()
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
    pipi_op_t_sep = get_param(job_tag, "measurement", "pipi_op_t_sep")
    pipi_corr_t_sep_list = get_param(job_tag, "measurement", "pipi_corr_t_sep_list")
    def load_data():
        t_t_list = q.get_mpi_chunk(
                [ (t_src, t_sep_idx,)
                 for t_src in range(total_site[3])
                 for t_sep_idx in range(len(pipi_corr_t_sep_list))
                 ],
                rng_state=None)
        for t_src, t_sep_idx in t_t_list:
            t_sep = pipi_corr_t_sep_list[t_sep_idx]
            t_snk = (t_src + t_sep) % t_size
            yield t_snk, t_src, t_sep_idx
    @q.timer
    def feval(args):
        t_snk, t_src, t_sep_idx, = args
        pd = {
                "snk_1" : ("wall", t_snk,),
                "snk_2" : ("wall", (t_snk + pipi_op_t_sep) % t_size,),
                "src_1" : ("wall", t_src,),
                "src_2" : ("wall", (t_src - pipi_op_t_sep) % t_size,),
                "size" : total_site,
                }
        val = eval_cexpr(cexpr, positions_dict=pd, get_prop=get_prop)
        return val, t_sep_idx
    def sum_function(val_list):
        values = np.zeros((len(pipi_corr_t_sep_list), len(expr_names),), dtype=complex)
        for val, t_sep_idx in val_list:
            values[t_sep_idx] += val
        return q.glb_sum(values.transpose(1, 0))
    auto_contractor_chunk_size = get_param(job_tag, "measurement", "auto_contractor_chunk_size", default=128)
    q.timer_fork(0)
    res_sum = q.parallel_map_sum(feval, load_data(), sum_function=sum_function, chunksize=auto_contractor_chunk_size)
    q.displayln_info(f"{fname}: timer_display for parallel_map_sum")
    q.timer_display()
    q.timer_merge()
    res_sum *= 1.0 / t_size
    assert q.qnorm(res_sum[0] - 1.0) < 1e-10
    ld = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t_sep", len(pipi_corr_t_sep_list), pipi_corr_t_sep_list, ],
        ])
    ld.from_numpy(res_sum)
    ld.save(get_save_path(fn))
    json_results_append(f"{fname}: ld sig", q.get_data_sig(ld, q.RngState()))
    for i, en in enumerate(expr_names):
        json_results_append(f"{fname}: ld '{en}' sig", q.get_data_sig(ld[i], q.RngState()))

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
        exprs = [ mk_expr(1) + f"1", ]
        exprs += exprs_self_energy + exprs_ope
        assert len(exprs) == 65
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
    t_sep = get_param(job_tag, "measurement", "meson_tensor_t_sep")
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
            t_2 = (max(x_1_t, x_2_t) + t_sep) % total_site[3]
            t_1 = (min(x_1_t, x_2_t) - t_sep) % total_site[3]
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
    json_results_append(f"{fname}: ld_sum sig", q.get_data_sig(ld_sum, q.RngState()))
    for i, en in enumerate(expr_names):
        json_results_append(f"{fname}: ld_sum '{en}' sig", q.get_data_sig(ld_sum[i], q.RngState()))

# ----

@q.timer
def get_cexpr_pipi_jj():
    fn_base = "cache/auto_contract_cexpr/get_cexpr_pipi_jj"
    def calc_cexpr():
        diagram_type_dict = dict()
        diagram_type_dict[((('snk_1', 'src_1'), 1), (('snk_2', 'src_2'), 1), (('src_1', 'snk_2'), 1), (('src_2', 'x_1'), 1), (('x_1', 'x_2'), 1), (('x_2', 'snk_1'), 1))] = 'Type1'
        diagram_type_dict[((('snk_1', 'src_1'), 1), (('snk_2', 'src_2'), 1), (('src_1', 'snk_1'), 1), (('src_2', 'x_1'), 1), (('x_1', 'x_2'), 1), (('x_2', 'snk_2'), 1))] = 'Type2'
        diagram_type_dict[((('snk_1', 'src_1'), 1), (('snk_2', 'x_1'), 1), (('src_1', 'snk_2'), 1), (('src_2', 'snk_1'), 1), (('x_1', 'x_2'), 1), (('x_2', 'src_2'), 1))] = 'Type3'
        diagram_type_dict[((('snk_1', 'src_1'), 1), (('snk_2', 'x_1'), 1), (('src_1', 'snk_1'), 1), (('src_2', 'snk_2'), 1), (('x_1', 'x_2'), 1), (('x_2', 'src_2'), 1))] = 'Type4'
        diagram_type_dict[((('snk_1', 'src_1'), 1), (('snk_2', 'x_1'), 1), (('src_1', 'snk_2'), 1), (('src_2', 'x_2'), 1), (('x_1', 'src_2'), 1), (('x_2', 'snk_1'), 1))] = 'Type5'
        diagram_type_dict[((('snk_1', 'src_1'), 1), (('snk_2', 'src_2'), 1), (('src_1', 'x_1'), 1), (('src_2', 'x_2'), 1), (('x_1', 'snk_1'), 1), (('x_2', 'snk_2'), 1))] = 'Type6'
        diagram_type_dict[((('snk_1', 'src_1'), 1), (('snk_2', 'x_1'), 1), (('src_1', 'x_2'), 1), (('src_2', 'snk_2'), 1), (('x_1', 'src_2'), 1), (('x_2', 'snk_1'), 1))] = 'Type7'
        diagram_type_dict[((('snk_1', 'src_1'), 1), (('snk_2', 'src_2'), 1), (('src_1', 'x_1'), 1), (('src_2', 'x_2'), 1), (('x_1', 'snk_2'), 1), (('x_2', 'snk_1'), 1))] = 'Type8'
        diagram_type_dict[((('snk_1', 'src_1'), 1), (('snk_2', 'x_1'), 1), (('src_1', 'x_2'), 1), (('src_2', 'snk_1'), 1), (('x_1', 'src_2'), 1), (('x_2', 'snk_2'), 1))] = 'Type9'
        diagram_type_dict[((('snk_1', 'src_1'), 1), (('snk_2', 'x_1'), 1), (('src_1', 'snk_1'), 1), (('src_2', 'x_2'), 1), (('x_1', 'src_2'), 1), (('x_2', 'snk_2'), 1))] = 'Type10'
        diagram_type_dict[((('snk_1', 'x_1'), 1), (('snk_2', 'x_2'), 1), (('src_1', 'snk_1'), 1), (('src_2', 'snk_2'), 1), (('x_1', 'src_1'), 1), (('x_2', 'src_2'), 1))] = 'Type11'
        diagram_type_dict[((('snk_1', 'x_1'), 1), (('snk_2', 'x_2'), 1), (('src_1', 'snk_1'), 1), (('src_2', 'snk_2'), 1), (('x_1', 'src_2'), 1), (('x_2', 'src_1'), 1))] = 'Type12'
        diagram_type_dict[((('snk_1', 'src_1'), 1), (('snk_2', 'src_2'), 1), (('src_1', 'snk_2'), 1), (('src_2', 'x_1'), 1), (('x_1', 'snk_1'), 1), (('x_2', 'x_2'), 1))] = None
        diagram_type_dict[((('snk_1', 'src_1'), 1), (('snk_2', 'src_2'), 1), (('src_1', 'snk_1'), 1), (('src_2', 'x_1'), 1), (('x_1', 'snk_2'), 1), (('x_2', 'x_2'), 1))] = None
        diagram_type_dict[((('snk_1', 'src_1'), 1), (('snk_2', 'x_1'), 1), (('src_1', 'snk_2'), 1), (('src_2', 'snk_1'), 1), (('x_1', 'src_2'), 1), (('x_2', 'x_2'), 1))] = None
        diagram_type_dict[((('snk_1', 'src_1'), 1), (('snk_2', 'x_1'), 1), (('src_1', 'snk_1'), 1), (('src_2', 'snk_2'), 1), (('x_1', 'src_2'), 1), (('x_2', 'x_2'), 1))] = None
        diagram_type_dict[((('snk_1', 'src_1'), 1), (('snk_2', 'src_2'), 1), (('src_1', 'snk_1'), 1), (('src_2', 'snk_2'), 1), (('x_1', 'x_1'), 1), (('x_2', 'x_2'), 1))] = None
        diagram_type_dict[((('snk_1', 'src_1'), 1), (('snk_2', 'src_2'), 1), (('src_1', 'snk_2'), 1), (('src_2', 'snk_1'), 1), (('x_1', 'x_1'), 1), (('x_2', 'x_2'), 1))] = None
        diagram_type_dict[((('snk_1', 'src_1'), 1), (('snk_2', 'src_2'), 1), (('src_1', 'snk_1'), 1), (('src_2', 'snk_2'), 1), (('x_1', 'x_2'), 1), (('x_2', 'x_1'), 1))] = None
        diagram_type_dict[((('snk_1', 'src_1'), 1), (('snk_2', 'src_2'), 1), (('src_1', 'snk_2'), 1), (('src_2', 'snk_1'), 1), (('x_1', 'x_2'), 1), (('x_2', 'x_1'), 1))] = None
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
                mk_pipi_i22("snk_1", "snk_2", True) * mk_pipi_i22("src_1", "src_2") + f"pipi_i22+^dag(0) * pipi_i22(-tsep)",
                ]
        exprs_self_energy = [ op * mm for mm in mm_list for op in op_list ]
        assert len(exprs_self_energy) == 15
        exprs = [ mk_fac(1) + f"1", ]
        exprs += exprs_self_energy
        assert len(exprs) == 16
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit=True, diagram_type_dict=diagram_type_dict)
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

@q.timer_verbose
def auto_contract_pipi_jj(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract/traj-{traj}/pipi_jj.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_pipi_jj()
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
    pipi_op_t_sep = get_param(job_tag, "measurement", "pipi_op_t_sep")
    pipi_tensor_t_sep_list = get_param(job_tag, "measurement", "pipi_tensor_t_sep_list")
    pipi_tensor_t_max = get_param(job_tag, "measurement", "pipi_tensor_t_max")
    pipi_tensor_r_max = get_param(job_tag, "measurement", "pipi_tensor_r_max")
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
            x_rel_t = x_rel[3]
            r_sq = q.get_r_sq(x_rel)
            r_norm = np.sqrt(r_sq)
            t_norm = abs(x_rel_t)
            if t_norm > pipi_tensor_t_max:
                continue
            if r_norm > pipi_tensor_r_max:
                continue
            r_idx_low, r_idx_high, coef_low, coef_high = r_sq_interp_idx_coef_list[r_sq]
            x_2_t = xg_src[3]
            x_1_t = x_2_t + x_rel_t
            val = np.zeros((len(pipi_tensor_t_sep_list), len(pipi_tensor_t_sep_list), len(expr_names),), dtype=complex)
            for t_sep_1_idx, t_sep_1 in enumerate(pipi_tensor_t_sep_list):
                for t_sep_2_idx, t_sep_2 in enumerate(pipi_tensor_t_sep_list):
                    t_snk = (max(x_1_t, x_2_t) + t_sep_2) % total_site[3]
                    t_src = (min(x_1_t, x_2_t) - t_sep_1) % total_site[3]
                    pd = {
                            "x_2" : ("point-snk", xg_snk,),
                            "x_1" : ("point", xg_src,),
                            "snk_1" : ("wall", t_snk,),
                            "snk_2" : ("wall", (t_snk + pipi_op_t_sep) % t_size,),
                            "src_1" : ("wall", t_src,),
                            "src_2" : ("wall", (t_src - pipi_op_t_sep) % t_size,),
                            "size" : total_site,
                            }
                    t = x_rel_t % t_size
                    val[t_sep_1_idx, t_sep_2_idx] = eval_cexpr(cexpr, positions_dict=pd, get_prop=get_prop)
            res_list.append((val / prob, t, r_idx_low, r_idx_high, coef_low, coef_high,))
        return res_list
    def sum_function(val_list):
        values = np.zeros((
            t_size,
            len(r_list),
            len(pipi_tensor_t_sep_list),
            len(pipi_tensor_t_sep_list),
            len(expr_names),
            ), dtype=complex)
        for idx, res_list in enumerate(val_list):
            for val, t, r_idx_low, r_idx_high, coef_low, coef_high in res_list:
                values[t, r_idx_low] += coef_low * val
                values[t, r_idx_high] += coef_high * val
            q.displayln_info(f"{fname}: {idx+1}/{len(xg_psel_arr)}")
        return q.glb_sum(values.transpose(4, 2, 3, 0, 1))
    q.timer_fork(0)
    res_sum = q.parallel_map_sum(feval, load_data(), sum_function=sum_function, chunksize=1)
    q.displayln_info(f"{fname}: timer_display for parallel_map_sum")
    q.timer_display()
    q.timer_merge()
    res_sum *= 1.0 / total_volume
    ld_sum = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t_sep_1", len(pipi_tensor_t_sep_list), pipi_tensor_t_sep_list, ],
        [ "t_sep_2", len(pipi_tensor_t_sep_list), pipi_tensor_t_sep_list, ],
        [ "t", t_size, [ str(q.rel_mod(t, t_size)) for t in range(t_size) ], ],
        [ "r", len(r_list), [ f"{r:.5f}" for r in r_list ], ],
        ])
    ld_sum.from_numpy(res_sum)
    ld_sum.save(get_save_path(fn))
    json_results_append(f"{fname}: ld_sum sig", q.get_data_sig(ld_sum, q.RngState()))
    for i, en in enumerate(expr_names):
        json_results_append(f"{fname}: ld_sum '{en}' sig", q.get_data_sig(ld_sum[i], q.RngState()))

# ----

### ------

@q.timer_verbose
def run_job_inversion(job_tag, traj):
    #
    traj_gf = traj
    if job_tag[:5] == "test-":
        # ADJUST ME
        traj_gf = 1000
        #
    #
    fns_produce = [
            (f"{job_tag}/prop-psrc-light/traj-{traj}.qar", f"{job_tag}/prop-psrc-light/traj-{traj}/geon-info.txt",),
            (f"{job_tag}/psel-prop-psrc-light/traj-{traj}.qar", f"{job_tag}/psel-prop-psrc-light/traj-{traj}/checkpoint.txt",),
            (f"{job_tag}/prop-psrc-strange/traj-{traj}.qar", f"{job_tag}/prop-psrc-strange/traj-{traj}/geon-info.txt",),
            (f"{job_tag}/psel-prop-psrc-strange/traj-{traj}.qar", f"{job_tag}/psel-prop-psrc-strange/traj-{traj}/checkpoint.txt",),
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
    #
    get_wi = run_wi(job_tag, traj)
    #
    def run_wsrc_full():
        get_eig = run_eig(job_tag, traj_gf, get_gf)
        # run_get_inverter(job_tag, traj, inv_type=0, get_gf=get_gf, get_gt=get_gt, get_eig=get_eig)
        run_prop_wsrc_full(job_tag, traj, inv_type=0, get_gf=get_gf, get_eig=get_eig, get_gt=get_gt, get_wi=get_wi)
        #
        get_eig = run_eig_strange(job_tag, traj_gf, get_gf)
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
    run_prop_wsrc_sparse(job_tag, traj, inv_type=0, get_gt=get_gt, get_psel=get_psel, get_fsel=get_fsel, get_wi=get_wi)
    run_prop_wsrc_sparse(job_tag, traj, inv_type=1, get_gt=get_gt, get_psel=get_psel, get_fsel=get_fsel, get_wi=get_wi)
    #
    def run_with_eig():
        get_eig = run_eig(job_tag, traj_gf, get_gf)
        # run_get_inverter(job_tag, traj, inv_type=0, get_gf=get_gf, get_eig=get_eig)
        # run_prop_wsrc(job_tag, traj, inv_type=0, get_gf=get_gf, get_eig=get_eig, get_gt=get_gt, get_psel=get_psel, get_fsel=get_fselc, get_wi=get_wi)
        run_prop_psrc(job_tag, traj, inv_type=0, get_gf=get_gf, get_eig=get_eig, get_gt=get_gt, get_psel=get_psel, get_fsel=get_fselc, get_f_rand_01=get_f_rand_01)
        q.clean_cache(q.cache_inv)
    #
    def run_with_eig_strange():
        get_eig = run_eig_strange(job_tag, traj_gf, get_gf)
        # run_get_inverter(job_tag, traj, inv_type=1, get_gf=get_gf, get_eig=get_eig)
        # run_prop_wsrc(job_tag, traj, inv_type=1, get_gf=get_gf, get_eig=get_eig, get_gt=get_gt, get_psel=get_psel, get_fsel=get_fselc, get_wi=get_wi)
        run_prop_psrc(job_tag, traj, inv_type=1, get_gf=get_gf, get_eig=get_eig, get_gt=get_gt, get_psel=get_psel, get_fsel=get_fselc, get_f_rand_01=get_f_rand_01)
        q.clean_cache(q.cache_inv)
    #
    def run_charm():
        # run_get_inverter(job_tag, traj, inv_type=2, get_gf=get_gf)
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
                auto_contract_meson_corr(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob)
                auto_contract_meson_jj(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob)
                auto_contract_pipi_corr(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob)
                auto_contract_pipi_jj(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob)
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
    benchmark_eval_cexpr(get_cexpr_meson_jj())
    benchmark_eval_cexpr(get_cexpr_pipi_corr())
    benchmark_eval_cexpr(get_cexpr_pipi_jj())

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
set_param("16IH2", "measurement", "meson_tensor_t_sep")(4)
set_param("16IH2", "measurement", "pipi_op_t_sep")(2)
set_param("16IH2", "measurement", "pipi_corr_t_sep_list")(list(range(1, 10)))
set_param("16IH2", "measurement", "pipi_tensor_t_sep_list")([ 1, 2, 3, 4, ])
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
                run_job_contract(job_tag, traj)
                if q.obtained_lock_history_list:
                    json_results_append(f"q.obtained_lock_history_list={q.obtained_lock_history_list}")
                    if job_tag[:5] != "test-":
                        gracefully_finish()

    q.check_log_json(__file__, json_results, check_eps=check_eps)

    gracefully_finish()

# ----
