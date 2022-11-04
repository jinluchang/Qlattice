#!/usr/bin/env python3

from auto_contractor.operators import *

import functools
import math
import os
import time
import importlib
import sys

from jobs import *
from load_data import *
from params import *

# ----

load_path_list[:] = [
        "results",
        "../qcddata",
        os.path.join(os.getenv("HOME"), "qcddata"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-gf-gt/results"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-sel/results"),
        # os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-fsel-self-loop/results"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-selected-data/results"),
        # os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-wsrc-prop/results"),
        # os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-psrc-prop/results"),
        # os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-smear-prop/results"),
        "/sdcc/u/jluchang/qcdqedta/luchang/data-gen/fill-wsnk-prop/results",
        ]

# ----

def r_scaling_factor():
    return 5.0

def get_r(x_rel):
    fac = r_scaling_factor()
    return fac * math.sqrt(q.c_sqr(x_rel))

def get_r_limit(total_site):
    return math.ceil(get_r([ total_site[i] // 2 for i in range(4) ])) + 1

def get_interp_idx_coef(x, limit = None):
    # return x_idx_low, x_idx_high, coef_low, coef_high
    x_idx_low = math.floor(x)
    x_idx_high = x_idx_low + 1
    if limit is not None:
        assert x_idx_high < limit
    coef_low = x_idx_high - x
    coef_high = x - x_idx_low
    return x_idx_low, x_idx_high, coef_low, coef_high

# ----

@q.timer
def get_cexpr_meson_corr():
    fn_base = "cache/auto_contract_cexpr/get_cexpr_meson_corr"
    def calc_cexpr():
        diagram_type_dict = dict()
        diagram_type_dict[((('t_1', 't_2'), 1), (('t_2', 't_1'), 1))] = 'Type1'
        exprs = [
                mk_pi_p("t_2", True) * mk_pi_p("t_1") + f"pi+^dag(t_sep) * pi+(0)",
                mk_k_p("t_2", True) * mk_k_p("t_1") + f"K+^dag(t_sep) * K+(0)",
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True, diagram_type_dict = diagram_type_dict)
        q.qtouch_info(fn_base + ".info.txt", display_cexpr(cexpr))
        cexpr.optimize()
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, fn_base)

@q.timer_verbose
def auto_contract_meson_corr(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_corr.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_corr()
    expr_names = get_cexpr_names(cexpr)
    total_site = rup.get_total_site(job_tag)
    fsel, fselc = get_fsel()
    xg_fsel_list = np.array(fsel.to_psel_local().to_list())
    geo = q.Geometry(total_site, 1)
    total_volume = geo.total_volume()
    def load_data():
        t_t_list = get_mpi_chunk(
                [ (t_src, t_snk,) for t_snk in range(total_site[3]) for t_src in range(total_site[3]) ],
                rng_state = q.RngState("get_mpi_chunk"))
        for t_src, t_snk in t_t_list:
            t = (t_snk - t_src) % total_site[3]
            pd = {
                    "t_2" : ("wall", t_snk,),
                    "t_1" : ("wall", t_src,),
                    "size" : total_site,
                    }
            yield pd, t
    @q.timer
    def feval(args):
        pd, t = args
        val = eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop)
        return val, t
    def sum_function(val_list):
        counts = np.zeros(total_site[3], dtype = complex)
        values = np.zeros((len(expr_names), total_site[3],), dtype = complex)
        for val, t in val_list:
            counts[t] += 1
            values[:, t] += val
        return counts, values
    q.timer_fork(0)
    res_count, res_sum = q.glb_sum(
            q.parallel_map_sum(feval, load_data(), sum_function = sum_function, chunksize = 16))
    q.displayln_info("timer_display for auto_contract_meson_corr")
    q.timer_display()
    q.timer_merge()
    res_count *= 1.0 / total_site[3]
    res_sum *= 1.0 / total_site[3]
    assert q.qnorm(res_count - 1.0) < 1e-10
    ld = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t_sep", total_site[3], ],
        ])
    ld.from_numpy(res_sum)
    # q.displayln_info(ld.show())
    ld.save(get_save_path(fn))

@q.timer_verbose
def auto_contract_meson_corr_psnk(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_corr_psnk.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_corr()
    expr_names = get_cexpr_names(cexpr)
    total_site = rup.get_total_site(job_tag)
    fsel, fselc = get_fsel()
    xg_fsel_list = np.array(fsel.to_psel_local().to_list())
    geo = q.Geometry(total_site, 1)
    total_volume = geo.total_volume()
    def load_data():
        for t_src in range(total_site[3]):
            for xg_snk in xg_fsel_list:
                xg_snk = tuple(xg_snk.tolist())
                t = (xg_snk[3] - t_src) % total_site[3]
                pd = {
                        "t_2" : ("point-snk", xg_snk,),
                        "t_1" : ("wall", t_src,),
                        "size" : total_site,
                        }
                yield pd, t
    @q.timer
    def feval(args):
        pd, t = args
        val = eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop)
        return val, t
    def sum_function(val_list):
        counts = np.zeros(total_site[3], dtype = complex)
        values = np.zeros((len(expr_names), total_site[3],), dtype = complex)
        for val, t in val_list:
            counts[t] += 1
            values[:, t] += val
        return counts, values
    q.timer_fork(0)
    res_count, res_sum = q.glb_sum(
            q.parallel_map_sum(feval, load_data(), sum_function = sum_function, chunksize = 16))
    q.displayln_info("timer_display for auto_contract_meson_corr_psnk")
    q.timer_display()
    q.timer_merge()
    res_count *= 1.0 / (total_volume * fsel.prob())
    res_sum *= 1.0 / (total_volume * fsel.prob())
    assert q.qnorm(res_count - 1.0) < 1e-10
    ld = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t_sep", total_site[3], ],
        ])
    ld.from_numpy(res_sum)
    # q.displayln_info(ld.show())
    ld.save(get_save_path(fn))

@q.timer_verbose
def auto_contract_meson_corr_psrc(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_corr_psrc.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_corr()
    expr_names = get_cexpr_names(cexpr)
    total_site = rup.get_total_site(job_tag)
    psel = get_psel()
    fsel, fselc = get_fsel()
    xg_fsel_list = np.array(fsel.to_psel_local().to_list())
    xg_psel_list = np.array(psel.to_list())
    geo = q.Geometry(total_site, 1)
    total_volume = geo.total_volume()
    def load_data():
        x_t_list = get_mpi_chunk(
                [ (tuple(xg_src.tolist()), t_snk,) for t_snk in range(total_site[3]) for xg_src in xg_psel_list ],
                rng_state = q.RngState("get_mpi_chunk"))
        for xg_src, t_snk in x_t_list:
            t = (xg_src[3] - t_snk) % total_site[3]
            pd = {
                    "t_2" : ("point", xg_src,),
                    "t_1" : ("wall", t_snk,),
                    "size" : total_site,
                    }
            yield pd, t
    @q.timer
    def feval(args):
        pd, t = args
        val = eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop)
        return val, t
    def sum_function(val_list):
        counts = np.zeros(total_site[3], dtype = complex)
        values = np.zeros((len(expr_names), total_site[3],), dtype = complex)
        for val, t in val_list:
            counts[t] += 1
            values[:, t] += val
        return counts, values
    q.timer_fork(0)
    res_count, res_sum = q.glb_sum(
            q.parallel_map_sum(feval, load_data(), sum_function = sum_function, chunksize = 16))
    q.displayln_info("timer_display for auto_contract_meson_corr_psrc")
    q.timer_display()
    q.timer_merge()
    res_count *= 1.0 / len(xg_psel_list)
    res_sum *= 1.0 / len(xg_psel_list)
    assert q.qnorm(res_count - 1.0) < 1e-10
    ld = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t_sep", total_site[3], ],
        ])
    ld.from_numpy(res_sum)
    # q.displayln_info(ld.show())
    ld.save(get_save_path(fn))

@q.timer_verbose
def auto_contract_meson_corr_psnk_psrc(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_corr_psnk_psrc.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_corr()
    expr_names = get_cexpr_names(cexpr)
    total_site = rup.get_total_site(job_tag)
    psel = get_psel()
    fsel, fselc = get_fsel()
    xg_fsel_list = np.array(fsel.to_psel_local().to_list())
    xg_psel_list = np.array(psel.to_list())
    geo = q.Geometry(total_site, 1)
    total_volume = geo.total_volume()
    def load_data():
        for idx, xg_src in enumerate(xg_psel_list):
            xg_src = tuple(xg_src.tolist())
            q.displayln_info(f"auto_contract_meson_corr_psnk_psrc: {idx+1}/{len(xg_psel_list)} {xg_src}")
            for xg_snk in xg_fsel_list:
                xg_snk = tuple(xg_snk.tolist())
                t = (xg_snk[3] - xg_src[3]) % total_site[3]
                pd = {
                        "t_2" : ("point-snk", xg_snk,),
                        "t_1" : ("point", xg_src,),
                        "size" : total_site,
                        }
                yield pd, t
    @q.timer
    def feval(args):
        pd, t = args
        val = eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop)
        return val, t
    def sum_function(val_list):
        counts = np.zeros(total_site[3], dtype = complex)
        values = np.zeros((len(expr_names), total_site[3],), dtype = complex)
        for val, t in val_list:
            counts[t] += 1
            values[:, t] += val
        return counts, values
    q.timer_fork(0)
    res_count, res_sum = q.glb_sum(
            q.parallel_map_sum(feval, load_data(), sum_function = sum_function, chunksize = 16))
    q.displayln_info("timer_display for auto_contract_meson_corr_psnk_psrc")
    q.timer_display()
    q.timer_merge()
    res_count *= 1.0 / (len(xg_psel_list) * total_volume * fsel.prob() / total_site[3])
    res_sum *= 1.0 / (len(xg_psel_list) * total_volume * fsel.prob() / total_site[3])
    assert q.qnorm(res_count - 1.0) < 1e-10
    ld = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t_sep", total_site[3], ],
        ])
    ld.from_numpy(res_sum)
    # q.displayln_info(ld.show())
    ld.save(get_save_path(fn))

# ----

@q.timer
def get_cexpr_meson_f_corr():
    fn_base = "cache/auto_contract_cexpr/get_cexpr_meson_f_corr"
    def calc_cexpr():
        diagram_type_dict = dict()
        diagram_type_dict[((('t_1', 'x_1'), 1), (('x_1', 't_1'), 1))] = 'Type1'
        exprs = [
                mk_j5pi_mu("x_1", 3) * mk_pi_p("t_1") + "A_pi(0) * pi+(-t_sep)",
                mk_j5k_mu("x_1", 3)  * mk_k_p("t_1") + "A_K(0) * K+(-t_sep)",
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True, diagram_type_dict = diagram_type_dict)
        q.qtouch_info(fn_base + ".info.txt", display_cexpr(cexpr))
        cexpr.optimize()
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, fn_base)

@q.timer_verbose
def auto_contract_meson_f_corr(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_f_corr.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_f_corr()
    expr_names = get_cexpr_names(cexpr)
    total_site = rup.get_total_site(job_tag)
    fsel, fselc = get_fsel()
    xg_fsel_list = np.array(fsel.to_psel_local().to_list())
    geo = q.Geometry(total_site, 1)
    total_volume = geo.total_volume()
    def load_data():
        for t_src in range(total_site[3]):
            for xg_snk in xg_fsel_list:
                xg_snk = tuple(xg_snk.tolist())
                t = (xg_snk[3] - t_src) % total_site[3]
                pd = {
                        "x_1" : ("point-snk", xg_snk,),
                        "t_1" : ("wall", t_src,),
                        "size" : total_site,
                        }
                yield pd, t
    @q.timer
    def feval(args):
        pd, t = args
        val = eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop)
        return val, t
    def sum_function(val_list):
        counts = np.zeros(total_site[3], dtype = complex)
        values = np.zeros((len(expr_names), total_site[3],), dtype = complex)
        for val, t in val_list:
            counts[t] += 1
            values[:, t] += val
        return counts, values
    q.timer_fork(0)
    res_count, res_sum = q.glb_sum(
            q.parallel_map_sum(feval, load_data(), sum_function = sum_function, chunksize = 16))
    q.displayln_info("timer_display for auto_contract_meson_f_corr")
    q.timer_display()
    q.timer_merge()
    res_count *= 1.0 / (total_volume * fsel.prob())
    res_sum *= 1.0 / (total_volume * fsel.prob())
    assert q.qnorm(res_count - 1.0) < 1e-10
    ld = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t_sep", total_site[3], ],
        ])
    ld.from_numpy(res_sum)
    # q.displayln_info(ld.show())
    ld.save(get_save_path(fn))

@q.timer_verbose
def auto_contract_meson_f_corr_psrc(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_f_corr_psrc.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_f_corr()
    expr_names = get_cexpr_names(cexpr)
    total_site = rup.get_total_site(job_tag)
    psel = get_psel()
    fsel, fselc = get_fsel()
    xg_fsel_list = np.array(fsel.to_psel_local().to_list())
    xg_psel_list = np.array(psel.to_list())
    geo = q.Geometry(total_site, 1)
    total_volume = geo.total_volume()
    def load_data():
        x_t_list = get_mpi_chunk(
                [ (tuple(xg_src.tolist()), t_snk,) for t_snk in range(total_site[3]) for xg_src in xg_psel_list ],
                rng_state = q.RngState("get_mpi_chunk"))
        for xg_src, t_snk in x_t_list:
            t = (xg_src[3] - t_snk) % total_site[3]
            pd = {
                    "x_1" : ("point", xg_src,),
                    "t_1" : ("wall", t_snk,),
                    "size" : total_site,
                    }
            yield pd, t
    @q.timer
    def feval(args):
        pd, t = args
        val = eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop)
        return val, t
    def sum_function(val_list):
        counts = np.zeros(total_site[3], dtype = complex)
        values = np.zeros((len(expr_names), total_site[3],), dtype = complex)
        for val, t in val_list:
            counts[t] += 1
            values[:, t] += val
        return counts, values
    q.timer_fork(0)
    res_count, res_sum = q.glb_sum(
            q.parallel_map_sum(feval, load_data(), sum_function = sum_function, chunksize = 16))
    q.displayln_info("timer_display for auto_contract_meson_f_corr_psrc")
    q.timer_display()
    q.timer_merge()
    res_count *= 1.0 / len(xg_psel_list)
    res_sum *= 1.0 / len(xg_psel_list)
    assert q.qnorm(res_count - 1.0) < 1e-10
    ld = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t_sep", total_site[3], ],
        ])
    ld.from_numpy(res_sum)
    # q.displayln_info(ld.show())
    ld.save(get_save_path(fn))

@q.timer_verbose
def auto_contract_meson_f_corr_psnk_psrc(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_f_corr_psnk_psrc.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_f_corr()
    expr_names = get_cexpr_names(cexpr)
    total_site = rup.get_total_site(job_tag)
    psel = get_psel()
    fsel, fselc = get_fsel()
    xg_fsel_list = np.array(fsel.to_psel_local().to_list())
    xg_psel_list = np.array(psel.to_list())
    geo = q.Geometry(total_site, 1)
    total_volume = geo.total_volume()
    def load_data():
        for idx, xg_src in enumerate(xg_psel_list):
            xg_src = tuple(xg_src.tolist())
            q.displayln_info(f"auto_contract_meson_f_corr_psnk_psrc: {idx+1}/{len(xg_psel_list)} {xg_src}")
            for xg_snk in xg_fsel_list:
                xg_snk = tuple(xg_snk.tolist())
                t = (xg_snk[3] - xg_src[3]) % total_site[3]
                pd = {
                        "x_1" : ("point-snk", xg_snk,),
                        "t_1" : ("point", xg_src,),
                        "size" : total_site,
                        }
                yield pd, t
    @q.timer
    def feval(args):
        pd, t = args
        val = eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop)
        return val, t
    def sum_function(val_list):
        counts = np.zeros(total_site[3], dtype = complex)
        values = np.zeros((len(expr_names), total_site[3],), dtype = complex)
        for val, t in val_list:
            counts[t] += 1
            values[:, t] += val
        return counts, values
    q.timer_fork(0)
    res_count, res_sum = q.glb_sum(
            q.parallel_map_sum(feval, load_data(), sum_function = sum_function, chunksize = 16))
    q.displayln_info("timer_display for auto_contract_meson_f_corr_psnk_psrc")
    q.timer_display()
    q.timer_merge()
    res_count *= 1.0 / (len(xg_psel_list) * total_volume * fsel.prob() / total_site[3])
    res_sum *= 1.0 / (len(xg_psel_list) * total_volume * fsel.prob() / total_site[3])
    assert q.qnorm(res_count - 1.0) < 1e-10
    ld = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t_sep", total_site[3], ],
        ])
    ld.from_numpy(res_sum)
    # q.displayln_info(ld.show())
    ld.save(get_save_path(fn))

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
        exprs = [
                mk_pi_p("t_2", True) * mk_vec_mu("u", "u", "x", 3) * mk_pi_p("t_1") + "pi+^dag(tsep) * ubar gt u * pi+(-tsep)",
                mk_k_p("t_2", True) * mk_vec_mu("u", "u", "x", 3) * mk_k_p("t_1") + "K+^dag(tsep) * ubar gt u * K+(-tsep)",
                -mk_k_p("t_2", True) * mk_vec_mu("s", "s", "x", 3) * mk_k_p("t_1") + "- K+^dag(tsep) * sbar gt s * K+(-tsep)",
                mk_pi_p("t_2p", True) * mk_vec_mu("u", "u", "x", 3) * mk_pi_p("t_1p") + "pi+^dag(T/2+tsep) * ubar gt u * pi+(T/2-tsep)",
                mk_k_p("t_2p", True) * mk_vec_mu("u", "u", "x", 3) * mk_k_p("t_1p") + "K+^dag(T/2+tsep) * ubar gt u * K+(T/2-tsep)",
                -mk_k_p("t_2p", True) * mk_vec_mu("s", "s", "x", 3) * mk_k_p("t_1p") + "- K+^dag(T/2+tsep) * sbar gt s * K+(T/2-tsep)",
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True, diagram_type_dict = diagram_type_dict)
        q.qtouch_info(fn_base + ".info.txt", display_cexpr(cexpr))
        cexpr.optimize()
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, fn_base)

@q.timer_verbose
def auto_contract_meson_jt(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_jt.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_jt()
    expr_names = get_cexpr_names(cexpr)
    total_site = rup.get_total_site(job_tag)
    psel = get_psel()
    fsel, fselc = get_fsel()
    xg_fsel_list = np.array(fsel.to_psel_local().to_list())
    xg_psel_list = np.array(psel.to_list())
    tsep = rup.dict_params[job_tag]["meson_tensor_tsep"]
    geo = q.Geometry(total_site, 1)
    total_volume = geo.total_volume()
    def load_data():
        for xg_snk in xg_fsel_list:
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
            yield pd
    @q.timer
    def feval(args):
        pd = args
        val = eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop)
        return val
    def sum_function(val_list):
        counts = 0.0
        values = np.zeros(len(expr_names), dtype = complex)
        for val in val_list:
            counts += 1.0
            values += val
        return counts, values
    q.timer_fork(0)
    res_count, res_sum = q.glb_sum(
            q.parallel_map_sum(feval, load_data(), sum_function = sum_function, chunksize = 16))
    q.displayln_info("timer_display for auto_contract_meson_jt")
    q.timer_display()
    q.timer_merge()
    res_count *= 1.0 / (total_volume * fsel.prob())
    res_sum *= 1.0 / (total_volume * fsel.prob())
    assert q.qnorm(res_count - 1.0) < 1e-10
    ld_sum = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        ])
    ld_sum.from_numpy(res_sum)
    # q.displayln_info(ld_sum.show())
    ld_sum.save(get_save_path(fn))

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
        exprs = [ m * mm for mm in mm_list for m in m_list ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True, diagram_type_dict = diagram_type_dict)
        q.qtouch_info(fn_base + ".info.txt", display_cexpr(cexpr))
        cexpr.optimize()
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, fn_base)

@q.timer_verbose
def auto_contract_meson_m(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_m.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_m()
    expr_names = get_cexpr_names(cexpr)
    total_site = rup.get_total_site(job_tag)
    psel = get_psel()
    fsel, fselc = get_fsel()
    xg_fsel_list = np.array(fsel.to_psel_local().to_list())
    xg_psel_list = np.array(psel.to_list())
    tsep = rup.dict_params[job_tag]["meson_tensor_tsep"]
    geo = q.Geometry(total_site, 1)
    total_volume = geo.total_volume()
    def load_data():
        for xg_snk in xg_fsel_list:
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
            yield pd
    @q.timer
    def feval(args):
        pd = args
        val = eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop)
        return val
    def sum_function(val_list):
        counts = 0.0
        values = np.zeros(len(expr_names), dtype = complex)
        for val in val_list:
            counts += 1.0
            values += val
        return counts, values
    q.timer_fork(0)
    res_count, res_sum = q.glb_sum(
            q.parallel_map_sum(feval, load_data(), sum_function = sum_function, chunksize = 16))
    q.displayln_info("timer_display for auto_contract_meson_m")
    q.timer_display()
    q.timer_merge()
    res_count *= 1.0 / (total_volume * fsel.prob())
    res_sum *= 1.0 / (total_volume * fsel.prob())
    assert q.qnorm(res_count - 1.0) < 1e-10
    ld_sum = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        ])
    ld_sum.from_numpy(res_sum)
    # q.displayln_info(ld_sum.show())
    ld_sum.save(get_save_path(fn))

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
                (sum([ mk_j_mu("x_2", mu) * mk_j_mu("x_1", mu) for mu in range(4) ])
                 + "j_mu(x) * j_mu(0)"),
                (mk_j_mu("x_2", 3) * mk_j_mu("x_1", 3)
                 + "j_t(x) * j_t(0)"),
                (sum([ mk_j_mu("x_2", mu) * mk_j_mu("x_1", mu) for mu in range(3) ])
                 + "j_i(x) * j_i(0)"),
                (sum([
                    mk_fac(f"rel_mod_sym(x_2[{mu}] - x_1[{mu}], size[{mu}])")
                    * mk_fac(f"rel_mod_sym(x_2[{nu}] - x_1[{nu}], size[{nu}])")
                    * mk_j_mu("x_2", mu) * mk_j_mu("x_1", nu)
                    for mu in range(3) for nu in range(3) ])
                 + "x[i] * x[j] * j_i(x) * j_j(0)"),
                (sum([
                    mk_fac(f"rel_mod_sym(x_2[{mu}] - x_1[{mu}], size[{mu}])")
                    * mk_j_mu("x_2", mu) * mk_j_mu("x_1", 3)
                    for mu in range(3) ])
                 + "x[i] * j_i(x) * j_t(0)"),
                (sum([
                    mk_fac(f"rel_mod_sym(x_1[{mu}] - x_2[{mu}], size[{mu}])")
                    * mk_j_mu("x_1", mu) * mk_j_mu("x_2", 3)
                    for mu in range(3) ])
                 + "-x[i] * j_i(-x) * j_t(0)"),
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
        def mk_jw_v_mu(p, mu):
            return mk_jpi_mu(p, mu) + mk_jk_mu(p, mu) + f"jw_v({p},{mu})"
        def mk_jw_a_mu(p, mu):
            return mk_j5pi_mu(p, mu) + mk_j5k_mu(p, mu) + f"jw_a({p},{mu})"
        jwj_list = [
                mk_jw_a_mu("x_1", 3) * mk_j_mu("x_2", 3)
                + "jw_a_t(0) * j_t(x)",
                sum([
                    mk_jw_a_mu("x_1", mu) * mk_j_mu("x_2", mu)
                    for mu in range(3) ])
                + "jw_a_i(0) * j_i(x)",
                sum([
                    mk_fac(f"rel_mod_sym(x_2[{mu}] - x_1[{mu}], size[{mu}])")
                    * mk_jw_a_mu("x_1", 3) * mk_j_mu("x_2", mu)
                    for mu in range(3) ])
                + "x[i] * jw_a_t(0) * j_i(x)",
                sum([
                    mk_fac(f"rel_mod_sym(x_2[{mu}] - x_1[{mu}], size[{mu}])")
                    * mk_jw_a_mu("x_1", mu) * mk_j_mu("x_2", 3)
                    for mu in range(3) ])
                + "x[i] * jw_a_i(0) * j_t(x)",
                sum([
                    mk_fac(f"rel_mod_sym(x_2[{mu}] - x_1[{mu}], size[{mu}])")
                    * mk_fac(f"rel_mod_sym(x_2[{nu}] - x_1[{nu}], size[{nu}])")
                    * mk_jw_a_mu("x_1", mu) * mk_j_mu("x_2", nu)
                    for mu in range(3) for nu in range(3) ])
                + "x[i] * x[j] * jw_a_i(0) * j_j(x)",
                sum([
                    q.epsilon_tensor(mu, nu, rho)
                    * mk_fac(f"rel_mod_sym(x_2[{mu}] - x_1[{mu}], size[{mu}])")
                    * mk_jw_v_mu("x_1", nu) * mk_j_mu("x_2", rho)
                    for mu in range(3) for nu in range(3) for rho in range(3) ])
                + "e(i,j,k) * x[i] * jw_v_j(0) * j_k(x)",
                ]
        assert len(jwj_list) == 6
        jjw_list = [
                mk_jw_a_mu("x_2", 3) * mk_j_mu("x_1", 3)
                + "jw_a_t(0) * j_t(-x)",
                sum([
                    mk_jw_a_mu("x_2", mu) * mk_j_mu("x_1", mu)
                    for mu in range(3) ])
                + "jw_a_i(0) * j_i(-x)",
                sum([
                    mk_fac(f"rel_mod_sym(x_1[{mu}] - x_2[{mu}], size[{mu}])")
                    * mk_jw_a_mu("x_2", 3) * mk_j_mu("x_1", mu)
                    for mu in range(3) ])
                + "-x[i] * jw_a_t(0) * j_i(-x)",
                sum([
                    mk_fac(f"rel_mod_sym(x_1[{mu}] - x_2[{mu}], size[{mu}])")
                    * mk_jw_a_mu("x_2", mu) * mk_j_mu("x_1", 3)
                    for mu in range(3) ])
                + "-x[i] * jw_a_i(0) * j_t(-x)",
                sum([
                    mk_fac(f"rel_mod_sym(x_1[{mu}] - x_2[{mu}], size[{mu}])")
                    * mk_fac(f"rel_mod_sym(x_1[{nu}] - x_2[{nu}], size[{nu}])")
                    * mk_jw_a_mu("x_2", mu) * mk_j_mu("x_1", nu)
                    for mu in range(3) for nu in range(3) ])
                + "-x[i] * -x[j] * jw_a_i(0) * j_j(-x)",
                sum([
                    q.epsilon_tensor(mu, nu, rho)
                    * mk_fac(f"rel_mod_sym(x_1[{mu}] - x_2[{mu}], size[{mu}])")
                    * mk_jw_v_mu("x_2", nu) * mk_j_mu("x_1", rho)
                    for mu in range(3) for nu in range(3) for rho in range(3) ])
                + "e(i,j,k) * -x[i] * jw_v_j(0) * j_k(-x)",
                ]
        assert len(jjw_list) == 6
        md_list = [
                mk_pi_p("t_1") + "pi+(-tsep)",
                mk_k_p("t_1") + "K+(-tsep)",
                mk_pi_p("t_2") + "pi+(x[t]+tsep)",
                mk_k_p("t_2") + "K+(x[t]+tsep)",
                ]
        assert len(md_list) == 4
        exprs_decay1 = [ jwj * md for md in md_list for jwj in jwj_list ]
        assert len(exprs_decay1) == 24
        exprs_decay2 = [ jjw * md for md in md_list for jjw in jjw_list ]
        assert len(exprs_decay2) == 24
        #
        jwm_list = [
                mk_jw_a_mu("x_1", 3) * mk_m("u", "x_2") + "jw_a_t(0) ubar_u(x)",
                mk_jw_a_mu("x_1", 3) * mk_m("d", "x_2") + "jw_a_t(0) dbar_d(x)",
                mk_jw_a_mu("x_1", 3) * mk_m("s", "x_2") + "jw_a_t(0) sbar_s(x)",
                mk_jw_a_mu("x_2", 3) * mk_m("u", "x_1") + "jw_a_t(0) ubar_u(-x)",
                mk_jw_a_mu("x_2", 3) * mk_m("d", "x_1") + "jw_a_t(0) dbar_d(-x)",
                mk_jw_a_mu("x_2", 3) * mk_m("s", "x_1") + "jw_a_t(0) sbar_s(-x)",
                ]
        assert len(jwm_list) == 6
        exprs_decay_m = [ jwm * md for md in md_list for jwm in jwm_list ]
        assert len(exprs_decay_m) == 24 # 4 * 6
        #
        exprs = exprs_self_energy + exprs_decay1 + exprs_decay2 + exprs_decay_m
        assert len(exprs) == 132
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True, diagram_type_dict = diagram_type_dict)
        q.qtouch_info(fn_base + ".info.txt", display_cexpr(cexpr))
        cexpr.optimize()
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, fn_base)

@q.timer_verbose
def auto_contract_meson_jj(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_jj.lat"
    fn_counts = f"{job_tag}/auto-contract/traj-{traj}/meson_jj_counts.lat"
    if get_load_path(fn) is not None:
        assert get_load_path(fn_counts) is not None
        return
    cexpr = get_cexpr_meson_jj()
    expr_names = get_cexpr_names(cexpr)
    total_site = rup.get_total_site(job_tag)
    psel = get_psel()
    fsel, fselc = get_fsel()
    xg_fsel_list = np.array(fsel.to_psel_local().to_list())
    xg_psel_list = np.array(psel.to_list())
    tsep = rup.dict_params[job_tag]["meson_tensor_tsep"]
    geo = q.Geometry(total_site, 1)
    t_size = total_site[3]
    r_limit = get_r_limit(total_site)
    def load_data():
        for idx, xg_src in enumerate(xg_psel_list):
            xg_src = tuple(xg_src.tolist())
            q.displayln_info(f"auto_contract_meson_jj: {idx+1}/{len(xg_psel_list)} {xg_src}")
            for xg_snk in xg_fsel_list:
                xg_snk = tuple(xg_snk.tolist())
                x_rel = [ q.rel_mod(xg_snk[mu] - xg_src[mu], total_site[mu]) for mu in range(4) ]
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
                r = get_r(x_rel)
                yield pd, t, r
    @q.timer
    def feval(args):
        pd, t, r = args
        val = eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop)
        return val, t, r
    def sum_function(val_list):
        counts = np.zeros((t_size, r_limit,), dtype = complex)
        values = np.zeros((t_size, r_limit, len(expr_names),), dtype = complex)
        for val, t, r in val_list:
            r_idx_low, r_idx_high, coef_low, coef_high = get_interp_idx_coef(r, r_limit)
            counts[t, r_idx_low] += coef_low
            counts[t, r_idx_high] += coef_high
            values[t, r_idx_low] += coef_low * val
            values[t, r_idx_high] += coef_high * val
        return counts, values
    q.timer_fork(0)
    res_count, res_sum = q.glb_sum(
            q.parallel_map_sum(feval, load_data(), sum_function = sum_function, chunksize = 16))
    q.displayln_info("timer_display for auto_contract_meson_jj")
    q.timer_display()
    q.timer_merge()
    res_count *= 1.0 / (len(xg_psel_list) * fsel.prob())
    res_sum *= 1.0 / (len(xg_psel_list) * fsel.prob())
    ld_count = q.mk_lat_data([
        [ "t", t_size, [ str(q.rel_mod(t, t_size)) for t in range(t_size) ], ],
        [ "r", r_limit, [ f"{r / r_scaling_factor():.2f}" for r in range(r_limit) ], ],
        ])
    ld_sum = q.mk_lat_data([
        [ "t", t_size, [ str(q.rel_mod(t, t_size)) for t in range(t_size) ], ],
        [ "r", r_limit, [ f"{r / r_scaling_factor():.2f}" for r in range(r_limit) ], ],
        [ "expr_name", len(expr_names), expr_names, ],
        ])
    ld_count.from_numpy(res_count)
    ld_sum.from_numpy(res_sum)
    ld_count.save(get_save_path(fn_counts))
    ld_sum.save(get_save_path(fn))
    # q.displayln_info(ld_count.show())
    # q.displayln_info(ld_sum.show())

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
        jj_list = [
                (sum([ mk_j_mu("x_2", mu) * mk_j_mu("x_1", mu) for mu in range(4) ])
                 + "j_mu(x) * j_mu(y)"),
                (mk_j_mu("x_2", 3) * mk_j_mu("x_1", 3)
                 + "j_t(x) * j_t(y)"),
                (sum([ mk_j_mu("x_2", mu) * mk_j_mu("x_1", mu) for mu in range(3) ])
                 + "j_i(x) * j_i(y)"),
                ]
        assert len(jj_list) == 3
        m2_list = [
                mk_m("u", "x_2") + "ubar_u(x)",
                mk_m("d", "x_2") + "dbar_d(x)",
                mk_m("s", "x_2") + "sbar_s(x)",
                ]
        assert len(m2_list) == 3
        m1_list = [
                mk_m("u", "x_1") + "ubar_u(y)",
                mk_m("d", "x_1") + "dbar_d(y)",
                mk_m("s", "x_1") + "sbar_s(y)",
                ]
        assert len(m1_list) == 3
        m1m2_list = [ m2 * m1 for m2 in m2_list for m1 in m1_list ]
        assert len(m1m2_list) == 9
        op_list = jj_list + m1m2_list
        assert len(op_list) == 12
        def mk_jw_v_mu(p, mu):
            return mk_jpi_mu(p, mu) + mk_jk_mu(p, mu)
        def mk_jw_a_mu(p, mu):
            return mk_j5pi_mu(p, mu) + mk_j5k_mu(p, mu)
        jm_list = [
                mk_jw_a_mu("w", 3) * mk_pi_p("t_1") + "jw_a_t(0) * pi+(-tsep)",
                mk_jw_a_mu("w", 3) * mk_k_p("t_1") + "jw_a_t(0) * K+(-tsep)",
                mk_jw_a_mu("w", 3) * mk_pi_p("t_2") + "jw_a_t(0) * pi+(tsep)",
                mk_jw_a_mu("w", 3) * mk_k_p("t_2") + "jw_a_t(0) * K+(tsep)",
                ]
        assert len(jm_list) == 4
        exprs = [ op * jm for jm in jm_list for op in op_list ]
        assert len(exprs) == 48
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True, diagram_type_dict = diagram_type_dict)
        q.qtouch_info(fn_base + ".info.txt", display_cexpr(cexpr))
        cexpr.optimize()
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, fn_base)

@q.timer_verbose
def auto_contract_meson_jwjj(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_jwjj.lat"
    fn_counts = f"{job_tag}/auto-contract/traj-{traj}/meson_jwjj.lat"
    if get_load_path(fn) is not None:
        assert get_load_path(fn_counts) is not None
        return
    cexpr = get_cexpr_meson_jwjj()
    expr_names = get_cexpr_names(cexpr)
    total_site = rup.get_total_site(job_tag)
    psel = get_psel()
    fsel, fselc = get_fsel()
    xg_fsel_list = np.array(fsel.to_psel_local().to_list())
    xg_psel_list = np.array(psel.to_list())
    tsep = rup.dict_params[job_tag]["meson_tensor_tsep"]
    geo = q.Geometry(total_site, 1)
    t_size = total_site[3]
    r_limit = get_r_limit(total_site)
    n_points = len(xg_psel_list)
    n_pairs = n_points * (n_points - 1) // 2 + n_points
    #
    threshold = rup.dict_params[job_tag]["meson_jwjj_threshold"]
    u_rand_prob = q.SelectedField("double", fsel, 1)
    u_rand_prob.set_rand(q.RngState(f"auto_contract_meson_jwjj,{job_tag},{traj}"), 1.0, 0.0)
    fn_meson_corr = f"{job_tag}/auto-contract/traj-{traj}/meson_corr_psnk.lat"
    if get_load_path(fn_meson_corr) is None:
        q.displayln_info("auto_contract_meson_jwjj: '{fn_meson_corr}' does not exist. Skipping.")
        return
    ld_meson_corr = q.load_lat_data(get_load_path(fn_meson_corr))
    arr_meson_corr = ld_meson_corr.to_numpy()
    def get_corr(t_op, t_meson, meson_type):
        # meson_type = 0 (for pion)
        # meson_type = 1 (for kaon)
        t = (t_op - t_meson) % t_size
        return abs(arr_meson_corr[meson_type, t].item())
    def get_prop_norm_sqrt(*args):
        is_sloppy = False
        return ama_extract(get_prop(*args, is_norm_sqrt = True), is_sloppy = is_sloppy)
    def get_estimate(xg_snk, xg1_src, xg2_src, t_1, t_2):
        corr1 = get_corr(xg_snk[3], t_1, 0)
        corr1s = get_corr(xg_snk[3], t_1, 1)
        corr2 = get_corr(xg_snk[3], t_2, 0)
        corr2s = get_corr(xg_snk[3], t_2, 1)
        w = ("point-snk", xg_snk,)
        x_1 = ("point", xg1_src,)
        x_2 = ("point", xg2_src,)
        t_1s = ("wall", t_1,)
        t_2s = ("wall", t_2,)
        p1t1 = get_prop_norm_sqrt("l", x_1, t_1s)
        p2t1 = get_prop_norm_sqrt("l", x_2, t_1s)
        p1t2 = get_prop_norm_sqrt("l", x_1, t_2s)
        p2t2 = get_prop_norm_sqrt("l", x_2, t_2s)
        wt1 = get_prop_norm_sqrt("l", w, t_1s)
        wt2 = get_prop_norm_sqrt("l", w, t_2s)
        p1t1s = get_prop_norm_sqrt("s", x_1, t_1s)
        p2t1s = get_prop_norm_sqrt("s", x_2, t_1s)
        p1t2s = get_prop_norm_sqrt("s", x_1, t_2s)
        p2t2s = get_prop_norm_sqrt("s", x_2, t_2s)
        wt1s = get_prop_norm_sqrt("s", w, t_1s)
        wt2s = get_prop_norm_sqrt("s", w, t_2s)
        wp1 = get_prop_norm_sqrt("l", x_1, w)
        wp2 = get_prop_norm_sqrt("l", x_2, w)
        wp1s = get_prop_norm_sqrt("s", x_1, w)
        wp2s = get_prop_norm_sqrt("s", x_2, w)
        p1p2 = get_prop_norm_sqrt("l", x_1, x_2)
        p1p2s = get_prop_norm_sqrt("s", x_1, x_2)
        values = [
                2 * (p1t1 * p2t1 / corr1 + p1t2 * p2t2 / corr2) * (wp1 * wp2),
                2 * (p1t1s * p2t1 / corr1s + p1t2s * p2t2 / corr2s) * (wp1s * wp2),
                2 * (p1t1 * p2t1s / corr1s + p1t2 * p2t2s / corr2s) * (wp1 * wp2s),
                5 * (p1t1 * wt1 / corr1 + p1t2 * wt2 / corr2) * (p1p2 * wp2),
                4 * (p1t1 * wt1s / corr1s + p1t2 * wt2s / corr2s) * (p1p2 * wp2),
                1 * (p1t1s * wt1 / corr1s + p1t2s * wt2 / corr2s) * (p1p2s * wp2s),
                5 * (p2t1 * wt1 / corr1 + p2t2 * wt2 / corr2) * (p1p2 * wp1),
                4 * (p2t1 * wt1s / corr1s + p2t2 * wt2s / corr2s) * (p1p2 * wp1),
                1 * (p2t1s * wt1 / corr1s + p2t2s * wt2 / corr2s) * (p1p2s * wp1s),
                ]
        q.displayln_info(6, f"get_estimate: {values}")
        return sum(values)
    def get_weight(idx_snk, xg_snk, xg1_src, xg2_src, t_1, t_2):
        # return weight for this point (1 / prob or zero)
        est = get_estimate(xg_snk, xg1_src, xg2_src, t_1, t_2)
        prob = est / threshold
        q.displayln_info(5, f"get_weight: {prob} {est} {threshold}")
        if prob >= 1:
            return 1
        rand = u_rand_prob.get_elem(idx_snk).item()
        if rand < prob:
            return 1.0 / prob
        else:
            return 0
    #
    def load_data():
        idx_pair = 0
        n_total = 0
        n_selected = 0
        for idx1, xg1_src in enumerate(xg_psel_list):
            xg1_src = tuple(xg1_src.tolist())
            xg1_src_t = xg1_src[3]
            for idx2, xg2_src in enumerate(xg_psel_list):
                xg2_src = tuple(xg2_src.tolist())
                xg2_src_t = xg2_src[3]
                x_rel = [ q.rel_mod(xg2_src[mu] - xg1_src[mu], total_site[mu]) for mu in range(4) ]
                if idx2 > idx1:
                    continue
                idx_pair += 1
                q.displayln_info(1, f"auto_contract_meson_jwjj: {idx_pair}/{n_pairs} {xg1_src} {xg2_src}")
                for idx_snk, xg_snk in enumerate(xg_fsel_list):
                    xg_snk = tuple(xg_snk.tolist())
                    xg_t = xg_snk[3]
                    xg1_xg_t = q.rel_mod(xg1_src_t - xg_t, t_size)
                    xg2_xg_t = q.rel_mod(xg2_src_t - xg_t, t_size)
                    t_1 = (min(0, xg1_xg_t, xg2_xg_t) + xg_t - tsep) % total_site[3]
                    t_2 = (max(0, xg1_xg_t, xg2_xg_t) + xg_t + tsep) % total_site[3]
                    weight = get_weight(idx_snk, xg_snk, xg1_src, xg2_src, t_1, t_2)
                    n_total += 1
                    if weight == 0:
                        continue
                    n_selected += 1
                    pd = {
                            "w" : ("point-snk", xg_snk,),
                            "x_1" : ("point", xg1_src,),
                            "x_2" : ("point", xg2_src,),
                            "t_1" : ("wall", t_1,),
                            "t_2" : ("wall", t_2,),
                            "size" : total_site,
                            }
                    t1 = t_1 % t_size
                    t2 = t_2 % t_size
                    r = get_r(x_rel)
                    yield weight, pd, t1, t2, r
                q.displayln_info(1, f"auto_contract_meson_jwjj: {idx_pair}/{n_pairs} n_total={n_total} n_selected={n_selected} ratio={n_selected/n_total}")
    @q.timer
    def feval(args):
        weight, pd, t1, t2, r = args
        val = eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop)
        return weight, val, t1, t2, r
    def sum_function(val_list):
        counts = np.zeros((t_size, t_size, r_limit,), dtype = complex)
        values = np.zeros((t_size, t_size, r_limit, len(expr_names),), dtype = complex)
        for weight, val, t1, t2, r in val_list:
            r_idx_low, r_idx_high, coef_low, coef_high = get_interp_idx_coef(r, r_limit)
            counts[t1, t2, r_idx_low] += coef_low * weight
            counts[t1, t2, r_idx_high] += coef_high * weight
            values[t1, t2, r_idx_low] += coef_low * weight * val
            values[t1, t2, r_idx_high] += coef_high * weight * val
        return counts, values
    q.timer_fork(0)
    res_count, res_sum = q.glb_sum(
            q.parallel_map_sum(feval, load_data(), sum_function = sum_function, chunksize = 16))
    q.displayln_info("timer_display for auto_contract_meson_jj")
    q.timer_display()
    q.timer_merge()
    res_count *= 1.0 / (len(xg_psel_list) * fsel.prob())
    res_sum *= 1.0 / (len(xg_psel_list) * fsel.prob())
    ld_count = q.mk_lat_data([
        [ "t1", t_size, [ str(q.rel_mod(t, t_size)) for t in range(t_size) ], ],
        [ "t2", t_size, [ str(q.rel_mod(t, t_size)) for t in range(t_size) ], ],
        [ "r", r_limit, [ f"{r / r_scaling_factor():.2f}" for r in range(r_limit) ], ],
        ])
    ld_sum = q.mk_lat_data([
        [ "t1", t_size, [ str(q.rel_mod(t, t_size)) for t in range(t_size) ], ],
        [ "t2", t_size, [ str(q.rel_mod(t, t_size)) for t in range(t_size) ], ],
        [ "r", r_limit, [ f"{r / r_scaling_factor():.2f}" for r in range(r_limit) ], ],
        [ "expr_name", len(expr_names), expr_names, ],
        ])
    ld_count.from_numpy(res_count)
    ld_sum.from_numpy(res_sum)
    ld_count.save(get_save_path(fn_counts))
    ld_sum.save(get_save_path(fn))
    # q.displayln_info(ld_count.show())
    # q.displayln_info(ld_sum.show())

### ------

@q.timer_verbose
def run_job(job_tag, traj):
    fns_produce = [
            f"{job_tag}/auto-contract/traj-{traj}/checkpoint.txt",
            ]
    fns_need = [
            # f"{job_tag}/prop-rand-u1-light/traj-{traj}/geon-info.txt",
            # f"{job_tag}/prop-rand-u1-strange/traj-{traj}/geon-info.txt",
            # f"{job_tag}/prop-rand-u1-charm/traj-{traj}/geon-info.txt",
            f"{job_tag}/prop-psrc-light/traj-{traj}/geon-info.txt",
            f"{job_tag}/prop-psrc-strange/traj-{traj}/geon-info.txt",
            f"{job_tag}/prop-wsrc-light/traj-{traj}/geon-info.txt",
            f"{job_tag}/prop-wsrc-strange/traj-{traj}/geon-info.txt",
            f"{job_tag}/psel-prop-wsrc-light/traj-{traj}/checkpoint.txt",
            f"{job_tag}/psel-prop-wsrc-strange/traj-{traj}/checkpoint.txt",
            f"{job_tag}/psel-prop-psrc-light/traj-{traj}/checkpoint.txt",
            f"{job_tag}/psel-prop-psrc-strange/traj-{traj}/checkpoint.txt",
            f"{job_tag}/gauge-transform/traj-{traj}.field",
            f"{job_tag}/point-selection/traj-{traj}.txt",
            f"{job_tag}/field-selection/traj-{traj}.field",
            f"{job_tag}/wall-src-info-light/traj-{traj}.txt",
            f"{job_tag}/wall-src-info-strange/traj-{traj}.txt",
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
    #
    get_psel = run_psel(job_tag, traj)
    get_fsel = run_fsel(job_tag, traj, get_psel)
    #
    get_wi = run_wi(job_tag, traj)
    get_psel_smear = run_psel_smear(job_tag, traj)
    #
    get_get_prop = run_get_prop(job_tag, traj,
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
                ],
            )
    #
    fn_checkpoint = f"{job_tag}/auto-contract/traj-{traj}/checkpoint.txt"
    if get_load_path(fn_checkpoint) is None:
        if q.obtain_lock(f"locks/{job_tag}-{traj}-auto-contract"):
            q.timer_fork()
            get_prop = get_get_prop()
            # ADJUST ME
            auto_contract_meson_corr_psnk(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contract_meson_jwjj(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contract_meson_jj(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contract_meson_jt(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contract_meson_m(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contract_meson_corr(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contract_meson_corr_psrc(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contract_meson_corr_psnk_psrc(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contract_meson_f_corr(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contract_meson_f_corr_psrc(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contract_meson_f_corr_psnk_psrc(job_tag, traj, get_prop, get_psel, get_fsel)
            #
            q.qtouch_info(get_save_path(fn_checkpoint))
            q.release_lock()
            q.displayln_info("timer_display for runjob")
            q.timer_display()
            q.timer_merge()
    #
    q.clean_cache()
    q.clear_all_caches()

def get_all_cexpr():
    benchmark_eval_cexpr(get_cexpr_meson_corr())
    benchmark_eval_cexpr(get_cexpr_meson_f_corr())
    benchmark_eval_cexpr(get_cexpr_meson_m())
    benchmark_eval_cexpr(get_cexpr_meson_jt())
    benchmark_eval_cexpr(get_cexpr_meson_jj())
    benchmark_eval_cexpr(get_cexpr_meson_jwjj())

def test():
    # ADJUST ME
    q.qremove_all_info("locks")
    q.qremove_all_info("cache")
    q.qremove_all_info("results")
    get_all_cexpr()
    run_job("test-4nt8", 1000)
    # run_job("test-4nt16", 1000)
    # run_job("16IH2", 1000)

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

q.begin(sys.argv, size_node_list)

# ADJUST ME
test()

# ADJUST ME
job_tags = [
        # "test-4nt8", "test-4nt16",
        # "32IH1",
        # "32IH2",
        # "24IH1",
        # "24IH2",
        # "24IH3",
        # "64I",
        # "48I",
        # "32D",
        # "32Dfine",
        # "24DH",
        # "24D",
        # "16IH2",
        # "32IfineH",
        # "32IcoarseH1",
        ]

q.check_time_limit()

for job_tag in job_tags:
    if job_tag == "48I":
        if q.get_num_node() != 4 * 32:
            continue
    elif job_tag == "64I":
        if q.get_num_node() != 4 * 64:
            continue
    elif q.get_num_node() > 4 * 16:
        continue
    q.displayln_info(pprint.pformat(rup.dict_params[job_tag]))
    for traj in rup.dict_params[job_tag]["trajs"]:
        run_job(job_tag, traj)

q.timer_display()

q.end()
