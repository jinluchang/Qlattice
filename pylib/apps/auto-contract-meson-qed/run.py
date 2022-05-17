#!/usr/bin/env python3

from auto_contractor.operators import *

import functools
import math
import os
import time
import importlib
import sys
import gc

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

def rel_mod(x, size):
    x = x % size
    assert x >= 0
    if 2 * x >= size:
        return x - size
    else:
        return x

def rel_mod_sym(x, size):
    x = x % size
    assert x >= 0
    if 2 * x > size:
        return x - size
    elif 2 * x < size:
        return x
    else:
        assert 2 * x == size
        return 0

# ----

@q.timer
def get_cexpr_meson_corr():
    def calc_cexpr():
        t_1, t_2 = ['t_1', 't_2']
        terms = [
                tr(gamma_5*S_l(t_1,t_2)*gamma_5*S_l(t_2,t_1)), # term_Type1_0001
                tr(gamma_5*S_l(t_1,t_2)*gamma_5*S_s(t_2,t_1)), # term_Type1_0002
                tr(gamma_5*S_l(t_2,t_1)*gamma_5*S_s(t_1,t_2)), # term_Type1_0003
                ]
        cexpr = contract_simplify_compile(*terms, is_isospin_symmetric_limit = True)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.collect_op()
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, "cache/auto_contract_cexpr/meson_corr_cexpr")

@q.timer_verbose
def auto_contract_meson_corr(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"auto-contract/{job_tag}/traj={traj}/meson_corr.lat"
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
                    }
            props = eval_cexpr_get_props(cexpr, positions_dict = pd, get_prop = get_prop)
            yield props, t
    @q.timer
    def feval(args):
        props, t = args
        val = eval_cexpr_eval(cexpr, props = props)
        return val, t
    def sum_function(val_list):
        counts = np.zeros(total_site[3], dtype = complex)
        values = np.zeros((len(expr_names), total_site[3],), dtype = complex)
        for val, t in val_list:
            counts[t] += 1
            values[:, t] += val
        return counts, values
    gc.collect()
    q.timer_fork(0)
    res_count, res_sum = q.glb_sum_list(
            q.parallel_map_sum(q.get_q_mp_proc(), feval, load_data(), sum_function = sum_function, chunksize = 16))
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
    q.displayln_info(ld.show())
    ld.save(get_save_path(fn))

@q.timer_verbose
def auto_contract_meson_corr_psnk(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"auto-contract/{job_tag}/traj={traj}/meson_corr_psnk.lat"
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
                xg_snk = xg_snk.tolist()
                t = (xg_snk[3] - t_src) % total_site[3]
                pd = {
                        "t_2" : ("point-snk", xg_snk,),
                        "t_1" : ("wall", t_src,),
                        }
                props = eval_cexpr_get_props(cexpr, positions_dict = pd, get_prop = get_prop)
                yield props, t
    @q.timer
    def feval(args):
        props, t = args
        val = eval_cexpr_eval(cexpr, props = props)
        return val, t
    def sum_function(val_list):
        counts = np.zeros(total_site[3], dtype = complex)
        values = np.zeros((len(expr_names), total_site[3],), dtype = complex)
        for val, t in val_list:
            counts[t] += 1
            values[:, t] += val
        return counts, values
    gc.collect()
    q.timer_fork(0)
    res_count, res_sum = q.glb_sum_list(
            q.parallel_map_sum(q.get_q_mp_proc(), feval, load_data(), sum_function = sum_function, chunksize = 16))
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
    q.displayln_info(ld.show())
    ld.save(get_save_path(fn))

@q.timer_verbose
def auto_contract_meson_corr_psrc(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"auto-contract/{job_tag}/traj={traj}/meson_corr_psrc.lat"
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
                [ (xg_src.tolist(), t_snk,) for t_snk in range(total_site[3]) for xg_src in xg_psel_list ],
                rng_state = q.RngState("get_mpi_chunk"))
        for xg_src, t_snk in x_t_list:
            t = (xg_src[3] - t_snk) % total_site[3]
            pd = {
                    "t_2" : ("point", xg_src,),
                    "t_1" : ("wall", t_snk,),
                    }
            props = eval_cexpr_get_props(cexpr, positions_dict = pd, get_prop = get_prop)
            yield props, t
    @q.timer
    def feval(args):
        props, t = args
        val = eval_cexpr_eval(cexpr, props = props)
        return val, t
    def sum_function(val_list):
        counts = np.zeros(total_site[3], dtype = complex)
        values = np.zeros((len(expr_names), total_site[3],), dtype = complex)
        for val, t in val_list:
            counts[t] += 1
            values[:, t] += val
        return counts, values
    gc.collect()
    q.timer_fork(0)
    res_count, res_sum = q.glb_sum_list(
            q.parallel_map_sum(q.get_q_mp_proc(), feval, load_data(), sum_function = sum_function, chunksize = 16))
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
    q.displayln_info(ld.show())
    ld.save(get_save_path(fn))

@q.timer_verbose
def auto_contract_meson_corr_psnk_psrc(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"auto-contract/{job_tag}/traj={traj}/meson_corr_psnk_psrc.lat"
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
            xg_src = xg_src.tolist()
            q.displayln_info(f"auto_contract_meson_corr_psnk_psrc: {idx+1}/{len(xg_psel_list)} {xg_src}")
            for xg_snk in xg_fsel_list:
                xg_snk = xg_snk.tolist()
                t = (xg_snk[3] - xg_src[3]) % total_site[3]
                pd = {
                        "t_2" : ("point-snk", xg_snk,),
                        "t_1" : ("point", xg_src,),
                        }
                props = eval_cexpr_get_props(cexpr, positions_dict = pd, get_prop = get_prop)
                yield props, t
    @q.timer
    def feval(args):
        props, t = args
        val = eval_cexpr_eval(cexpr, props = props)
        return val, t
    def sum_function(val_list):
        counts = np.zeros(total_site[3], dtype = complex)
        values = np.zeros((len(expr_names), total_site[3],), dtype = complex)
        for val, t in val_list:
            counts[t] += 1
            values[:, t] += val
        return counts, values
    gc.collect()
    q.timer_fork(0)
    res_count, res_sum = q.glb_sum_list(
            q.parallel_map_sum(q.get_q_mp_proc(), feval, load_data(), sum_function = sum_function, chunksize = 16))
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
    q.displayln_info(ld.show())
    ld.save(get_save_path(fn))

# ----

@q.timer
def get_cexpr_meson_f_corr():
    def calc_cexpr():
        t_1, x_2 = ['t_1', 'x_2']
        terms = [
                tr(gamma_t*gamma_5*S_l(x_2,t_1)*gamma_5*S_l(t_1,x_2)), # term_Type1_0001
                tr(gamma_t*gamma_5*S_l(x_2,t_1)*gamma_5*S_s(t_1,x_2)), # term_Type1_0002
                tr(gamma_t*gamma_5*S_s(x_2,t_1)*gamma_5*S_l(t_1,x_2)), # term_Type1_0003
                ]
        cexpr = contract_simplify_compile(*terms, is_isospin_symmetric_limit = True)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.collect_op()
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, "cache/auto_contract_cexpr/meson_f_corr_cexpr")

@q.timer_verbose
def auto_contract_meson_f_corr_psnk(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"auto-contract/{job_tag}/traj={traj}/meson_f_corr_psnk.lat"
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
                xg_snk = xg_snk.tolist()
                t = (xg_snk[3] - t_src) % total_site[3]
                pd = {
                        "x_2" : ("point-snk", xg_snk,),
                        "t_1" : ("wall", t_src,),
                        }
                props = eval_cexpr_get_props(cexpr, positions_dict = pd, get_prop = get_prop)
                yield props, t
    @q.timer
    def feval(args):
        props, t = args
        val = eval_cexpr_eval(cexpr, props = props)
        return val, t
    def sum_function(val_list):
        counts = np.zeros(total_site[3], dtype = complex)
        values = np.zeros((len(expr_names), total_site[3],), dtype = complex)
        for val, t in val_list:
            counts[t] += 1
            values[:, t] += val
        return counts, values
    gc.collect()
    q.timer_fork(0)
    res_count, res_sum = q.glb_sum_list(
            q.parallel_map_sum(q.get_q_mp_proc(), feval, load_data(), sum_function = sum_function, chunksize = 16))
    q.displayln_info("timer_display for auto_contract_meson_f_corr_psnk")
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
    q.displayln_info(ld.show())
    ld.save(get_save_path(fn))

@q.timer_verbose
def auto_contract_meson_f_corr_psrc(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"auto-contract/{job_tag}/traj={traj}/meson_f_corr_psrc.lat"
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
                [ (xg_src.tolist(), t_snk,) for t_snk in range(total_site[3]) for xg_src in xg_psel_list ],
                rng_state = q.RngState("get_mpi_chunk"))
        for xg_src, t_snk in x_t_list:
            t = (xg_src[3] - t_snk) % total_site[3]
            pd = {
                    "x_2" : ("point", xg_src,),
                    "t_1" : ("wall", t_snk,),
                    }
            props = eval_cexpr_get_props(cexpr, positions_dict = pd, get_prop = get_prop)
            yield props, t
    @q.timer
    def feval(args):
        props, t = args
        val = eval_cexpr_eval(cexpr, props = props)
        return val, t
    def sum_function(val_list):
        counts = np.zeros(total_site[3], dtype = complex)
        values = np.zeros((len(expr_names), total_site[3],), dtype = complex)
        for val, t in val_list:
            counts[t] += 1
            values[:, t] += val
        return counts, values
    gc.collect()
    q.timer_fork(0)
    res_count, res_sum = q.glb_sum_list(
            q.parallel_map_sum(q.get_q_mp_proc(), feval, load_data(), sum_function = sum_function, chunksize = 16))
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
    q.displayln_info(ld.show())
    ld.save(get_save_path(fn))

@q.timer_verbose
def auto_contract_meson_f_corr_psnk_psrc(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"auto-contract/{job_tag}/traj={traj}/meson_f_corr_psnk_psrc.lat"
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
            xg_src = xg_src.tolist()
            q.displayln_info(f"auto_contract_meson_f_corr_psnk_psrc: {idx+1}/{len(xg_psel_list)} {xg_src}")
            for xg_snk in xg_fsel_list:
                xg_snk = xg_snk.tolist()
                t = (xg_snk[3] - xg_src[3]) % total_site[3]
                pd = {
                        "x_2" : ("point-snk", xg_snk,),
                        "t_1" : ("point", xg_src,),
                        }
                props = eval_cexpr_get_props(cexpr, positions_dict = pd, get_prop = get_prop)
                yield props, t
    @q.timer
    def feval(args):
        props, t = args
        val = eval_cexpr_eval(cexpr, props = props)
        return val, t
    def sum_function(val_list):
        counts = np.zeros(total_site[3], dtype = complex)
        values = np.zeros((len(expr_names), total_site[3],), dtype = complex)
        for val, t in val_list:
            counts[t] += 1
            values[:, t] += val
        return counts, values
    gc.collect()
    q.timer_fork(0)
    res_count, res_sum = q.glb_sum_list(
            q.parallel_map_sum(q.get_q_mp_proc(), feval, load_data(), sum_function = sum_function, chunksize = 16))
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
    q.displayln_info(ld.show())
    ld.save(get_save_path(fn))

# ----

@q.timer
def get_cexpr_meson_m():
    def calc_cexpr():
        t_1, t_2, x_1 = ['t_1', 't_2', 'x_1']
        exprs = [
                tr(gamma_5*S_l(t_1,t_2)*gamma_5*S_l(t_2,x_1)*S_l(x_1,t_1)), # term_Type1_0001
                tr(gamma_5*S_l(t_1,x_1)*S_l(x_1,t_2)*gamma_5*S_l(t_2,t_1)), # term_Type1_0002
                tr(gamma_5*S_l(t_1,x_1)*S_l(x_1,t_2)*gamma_5*S_s(t_2,t_1)), # term_Type1_0003
                tr(gamma_5*S_l(t_2,x_1)*S_l(x_1,t_1)*gamma_5*S_s(t_1,t_2)), # term_Type1_0004
                tr(gamma_5*S_l(t_1,t_2)*gamma_5*S_s(t_2,x_1)*S_s(x_1,t_1)), # term_Type1_0005
                tr(gamma_5*S_l(t_2,t_1)*gamma_5*S_s(t_1,x_1)*S_s(x_1,t_2)), # term_Type1_0006
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.collect_op()
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, "cache/auto_contract_cexpr/meson_m_cexpr")

@q.timer_verbose
def auto_contract_meson_m(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"auto-contract/{job_tag}/traj={traj}/meson_m.lat"
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
            xg_snk = xg_snk.tolist()
            t = xg_snk[3]
            t_1 = (t + tsep) % total_site[3]
            t_2 = (t - tsep) % total_site[3]
            pd = {
                    "x_1" : ("point-snk", xg_snk,),
                    "t_1" : ("wall", t_1),
                    "t_2" : ("wall", t_2),
                    }
            props = eval_cexpr_get_props(cexpr, positions_dict = pd, get_prop = get_prop)
            yield props
    @q.timer
    def feval(args):
        props = args
        val = eval_cexpr_eval(cexpr, props = props)
        return val
    def sum_function(val_list):
        counts = 0.0
        values = np.zeros(len(expr_names), dtype = complex)
        for val in val_list:
            counts += 1.0
            values += val
        return counts, values
    gc.collect()
    q.timer_fork(0)
    res_count, res_sum = q.glb_sum_list(
            q.parallel_map_sum(q.get_q_mp_proc(), feval, load_data(), sum_function = sum_function, chunksize = 16))
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
    ld_sum.save(get_save_path(fn))
    q.displayln_info(ld_sum.show())

# ----

@q.timer
def get_cexpr_meson_jt():
    def calc_cexpr():
        t_1, t_1p, t_2, t_2p, x_1 = ['t_1', 't_1p', 't_2', 't_2p', 'x_1']
        terms = [
                tr(gamma_t*S_l(x_1,t_1)*gamma_5*S_l(t_1,t_2)*gamma_5*S_l(t_2,x_1)), # term_Type1_0001
                tr(gamma_t*S_l(x_1,t_2)*gamma_5*S_l(t_2,t_1)*gamma_5*S_l(t_1,x_1)), # term_Type1_0002
                tr(gamma_t*S_l(x_1,t_1)*gamma_5*S_s(t_1,t_2)*gamma_5*S_l(t_2,x_1)), # term_Type1_0003
                tr(gamma_t*S_l(x_1,t_2)*gamma_5*S_s(t_2,t_1)*gamma_5*S_l(t_1,x_1)), # term_Type1_0004
                tr(gamma_t*S_s(x_1,t_1)*gamma_5*S_l(t_1,t_2)*gamma_5*S_s(t_2,x_1)), # term_Type1_0005
                tr(gamma_t*S_s(x_1,t_2)*gamma_5*S_l(t_2,t_1)*gamma_5*S_s(t_1,x_1)), # term_Type1_0006
                tr(gamma_t*S_l(x_1,t_1p)*gamma_5*S_l(t_1p,t_2p)*gamma_5*S_l(t_2p,x_1)), # term_Type2_0001
                tr(gamma_t*S_l(x_1,t_2p)*gamma_5*S_l(t_2p,t_1p)*gamma_5*S_l(t_1p,x_1)), # term_Type2_0002
                tr(gamma_t*S_l(x_1,t_1p)*gamma_5*S_s(t_1p,t_2p)*gamma_5*S_l(t_2p,x_1)), # term_Type2_0003
                tr(gamma_t*S_l(x_1,t_2p)*gamma_5*S_s(t_2p,t_1p)*gamma_5*S_l(t_1p,x_1)), # term_Type2_0004
                tr(gamma_t*S_s(x_1,t_1p)*gamma_5*S_l(t_1p,t_2p)*gamma_5*S_s(t_2p,x_1)), # term_Type2_0005
                tr(gamma_t*S_s(x_1,t_2p)*gamma_5*S_l(t_2p,t_1p)*gamma_5*S_s(t_1p,x_1)), # term_Type2_0006
                ]
        cexpr = contract_simplify_compile(*terms, is_isospin_symmetric_limit = True)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.collect_op()
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, "cache/auto_contract_cexpr/meson_jt_cexpr")

@q.timer_verbose
def auto_contract_meson_jt(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"auto-contract/{job_tag}/traj={traj}/meson_jt.lat"
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
            xg_snk = xg_snk.tolist()
            t = xg_snk[3]
            t_1 = (t + tsep) % total_site[3]
            t_2 = (t - tsep) % total_site[3]
            t_1p = (t_1 + total_site[3] // 2) % total_site[3]
            t_2p = (t_2 + total_site[3] // 2) % total_site[3]
            pd = {
                    "x_1" : ("point-snk", xg_snk,),
                    "t_1" : ("wall", t_1),
                    "t_2" : ("wall", t_2),
                    "t_1p" : ("wall", t_1p),
                    "t_2p" : ("wall", t_2p),
                    }
            props = eval_cexpr_get_props(cexpr, positions_dict = pd, get_prop = get_prop)
            yield props
    @q.timer
    def feval(args):
        props = args
        val = eval_cexpr_eval(cexpr, props = props)
        return val
    def sum_function(val_list):
        counts = 0.0
        values = np.zeros(len(expr_names), dtype = complex)
        for val in val_list:
            counts += 1.0
            values += val
        return counts, values
    gc.collect()
    q.timer_fork(0)
    res_count, res_sum = q.glb_sum_list(
            q.parallel_map_sum(q.get_q_mp_proc(), feval, load_data(), sum_function = sum_function, chunksize = 16))
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
    ld_sum.save(get_save_path(fn))
    q.displayln_info(ld_sum.show())

# ----

@q.timer
def get_cexpr_meson_jj():
    def calc_cexpr():
        t_1, t_2, x_1, x_2 = ['t_1', 't_2', 'x_1', 'x_2']
        def mk_terms(mu, nu):
            return [
                    tr(gamma(mu)*S_l(x_1,x_2)*gamma(nu)*S_l(x_2,x_1)), # term_Type0_0001
                    tr(gamma(mu)*S_s(x_1,x_2)*gamma(nu)*S_s(x_2,x_1)), # term_Type0_0002
                    tr(gamma(mu)*S_l(x_1,t_1)*gamma_5*S_l(t_1,x_1))*tr(gamma(nu)*S_l(x_2,t_2)*gamma_5*S_l(t_2,x_2)), # term_Type1_0001
                    tr(gamma(mu)*S_l(x_1,t_2)*gamma_5*S_l(t_2,x_1))*tr(gamma(nu)*S_l(x_2,t_1)*gamma_5*S_l(t_1,x_2)), # term_Type1_0002
                    tr(gamma(mu)*S_l(x_1,t_1)*gamma_5*S_l(t_1,x_2)*gamma(nu)*S_l(x_2,t_2)*gamma_5*S_l(t_2,x_1)), # term_Type2_0001
                    tr(gamma(mu)*S_l(x_1,t_2)*gamma_5*S_l(t_2,x_2)*gamma(nu)*S_l(x_2,t_1)*gamma_5*S_l(t_1,x_1)), # term_Type2_0002
                    tr(gamma(mu)*S_l(x_1,t_1)*gamma_5*S_s(t_1,x_2)*gamma(nu)*S_s(x_2,t_2)*gamma_5*S_l(t_2,x_1)), # term_Type2_0003
                    tr(gamma(mu)*S_l(x_1,t_2)*gamma_5*S_s(t_2,x_2)*gamma(nu)*S_s(x_2,t_1)*gamma_5*S_l(t_1,x_1)), # term_Type2_0004
                    tr(gamma(mu)*S_s(x_1,t_1)*gamma_5*S_l(t_1,x_2)*gamma(nu)*S_l(x_2,t_2)*gamma_5*S_s(t_2,x_1)), # term_Type2_0005
                    tr(gamma(mu)*S_s(x_1,t_2)*gamma_5*S_l(t_2,x_2)*gamma(nu)*S_l(x_2,t_1)*gamma_5*S_s(t_1,x_1)), # term_Type2_0006
                    tr(gamma(mu)*S_l(x_1,t_1)*gamma_5*S_l(t_1,t_2)*gamma_5*S_l(t_2,x_2)*gamma(nu)*S_l(x_2,x_1)), # term_Type3_0001
                    tr(gamma(mu)*S_l(x_1,t_2)*gamma_5*S_l(t_2,t_1)*gamma_5*S_l(t_1,x_2)*gamma(nu)*S_l(x_2,x_1)), # term_Type3_0002
                    tr(gamma(mu)*S_l(x_1,x_2)*gamma(nu)*S_l(x_2,t_1)*gamma_5*S_l(t_1,t_2)*gamma_5*S_l(t_2,x_1)), # term_Type3_0003
                    tr(gamma(mu)*S_l(x_1,x_2)*gamma(nu)*S_l(x_2,t_2)*gamma_5*S_l(t_2,t_1)*gamma_5*S_l(t_1,x_1)), # term_Type3_0004
                    tr(gamma(mu)*S_l(x_1,t_1)*gamma_5*S_s(t_1,t_2)*gamma_5*S_l(t_2,x_2)*gamma(nu)*S_l(x_2,x_1)), # term_Type3_0005
                    tr(gamma(mu)*S_l(x_1,t_2)*gamma_5*S_s(t_2,t_1)*gamma_5*S_l(t_1,x_2)*gamma(nu)*S_l(x_2,x_1)), # term_Type3_0006
                    tr(gamma(mu)*S_l(x_1,x_2)*gamma(nu)*S_l(x_2,t_1)*gamma_5*S_s(t_1,t_2)*gamma_5*S_l(t_2,x_1)), # term_Type3_0007
                    tr(gamma(mu)*S_l(x_1,x_2)*gamma(nu)*S_l(x_2,t_2)*gamma_5*S_s(t_2,t_1)*gamma_5*S_l(t_1,x_1)), # term_Type3_0008
                    tr(gamma(mu)*S_s(x_1,t_1)*gamma_5*S_l(t_1,t_2)*gamma_5*S_s(t_2,x_2)*gamma(nu)*S_s(x_2,x_1)), # term_Type3_0009
                    tr(gamma(mu)*S_s(x_1,t_2)*gamma_5*S_l(t_2,t_1)*gamma_5*S_s(t_1,x_2)*gamma(nu)*S_s(x_2,x_1)), # term_Type3_0010
                    tr(gamma(mu)*S_s(x_1,x_2)*gamma(nu)*S_s(x_2,t_1)*gamma_5*S_l(t_1,t_2)*gamma_5*S_s(t_2,x_1)), # term_Type3_0011
                    tr(gamma(mu)*S_s(x_1,x_2)*gamma(nu)*S_s(x_2,t_2)*gamma_5*S_l(t_2,t_1)*gamma_5*S_s(t_1,x_1)), # term_Type3_0012
                    tr(gamma(mu)*S_l(x_1,x_2)*gamma(nu)*S_l(x_2,x_1))*tr(gamma_5*S_l(t_1,t_2)*gamma_5*S_l(t_2,t_1)), # term_Type4_0001
                    tr(gamma(mu)*S_s(x_1,x_2)*gamma(nu)*S_s(x_2,x_1))*tr(gamma_5*S_l(t_1,t_2)*gamma_5*S_l(t_2,t_1)), # term_Type4_0002
                    tr(gamma(mu)*S_l(x_1,x_2)*gamma(nu)*S_l(x_2,x_1))*tr(gamma_5*S_l(t_1,t_2)*gamma_5*S_s(t_2,t_1)), # term_Type4_0003
                    tr(gamma(mu)*S_l(x_1,x_2)*gamma(nu)*S_l(x_2,x_1))*tr(gamma_5*S_l(t_2,t_1)*gamma_5*S_s(t_1,t_2)), # term_Type4_0004
                    tr(gamma(mu)*S_s(x_1,x_2)*gamma(nu)*S_s(x_2,x_1))*tr(gamma_5*S_l(t_1,t_2)*gamma_5*S_s(t_2,t_1)), # term_Type4_0005
                    tr(gamma(mu)*S_s(x_1,x_2)*gamma(nu)*S_s(x_2,x_1))*tr(gamma_5*S_l(t_2,t_1)*gamma_5*S_s(t_1,t_2)), # term_Type4_0006
                    ]
        n_tensor = 28
        terms_mu_nu = [ [ mk_terms(mu, nu) for nu in range(4) ] for mu in range(4) ]
        for t_nu in terms_mu_nu:
            for t in t_nu:
                assert n_tensor == len(t)
        terms = []
        terms += [ terms_mu_nu[mu][nu][i] for i in range(n_tensor) for mu in range(4) for nu in range(4) ]
        # Also calculate the meson corr
        terms += [
                tr(gamma_5*S_l(t_1,t_2)*gamma_5*S_l(t_2,t_1)), # term_Type1_0001
                tr(gamma_5*S_l(t_1,t_2)*gamma_5*S_s(t_2,t_1)), # term_Type1_0002
                tr(gamma_5*S_l(t_2,t_1)*gamma_5*S_s(t_1,t_2)), # term_Type1_0003
                ]
        cexpr = contract_simplify_compile(*terms, is_isospin_symmetric_limit = True)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.collect_op()
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, "cache/auto_contract_cexpr/meson_jj_cexpr")

def r_scaling_factor():
    return 5.0

def get_r(x_rel):
    fac = r_scaling_factor()
    return fac * math.sqrt(x_rel[0] * x_rel[0] + x_rel[1] * x_rel[1] + x_rel[2] * x_rel[2])

def get_r_limit(total_site):
    return math.ceil(get_r([ total_site[i] // 2 for i in range(4) ])) + 1

def jj_proj_mm(res_arr, x_rel_sym):
    # res_arr is 3-D np.array with n*4*4 elements
    # res_arr[idx_tensor, mu, nu]
    return np.array([ np.trace(res_arr[idx_tensor]) for idx_tensor in range(len(res_arr)) ])

def jj_proj_tt(res_arr, x_rel_sym):
    # res_arr is 3-D np.array with n*4*4 elements
    # res_arr[idx_tensor, mu, nu]
    return np.array([ res_arr[idx_tensor, 3, 3] for idx_tensor in range(len(res_arr)) ])

def jj_proj_ii(res_arr, x_rel_sym):
    # res_arr is 3-D np.array with n*4*4 elements
    # res_arr[idx_tensor, mu, nu]
    return np.array([ sum([ res_arr[idx_tensor, i, i] for i in range(3) ]) for idx_tensor in range(len(res_arr)) ])

def jj_proj_xx(res_arr, x_rel_sym):
    # res_arr is 3-D np.array with n*4*4 elements
    # res_arr[idx_tensor, mu, nu]
    return np.array(
            [ sum(
                [ x_rel_sym[i] * x_rel_sym[j] * res_arr[idx_tensor, i, j]
                    for i in range(3)
                    for j in range(3) ])
                for idx_tensor in range(len(res_arr)) ])

all_jj_projections = [
        jj_proj_mm,
        jj_proj_tt,
        jj_proj_ii,
        jj_proj_xx,
        ]

all_jj_projection_names = [ "mm", "tt", "ii", "xx", ]

def get_interp_idx_coef(x, limit = None):
    # return x_idx_low, x_idx_high, coef_low, coef_high
    x_idx_low = math.floor(x)
    x_idx_high = x_idx_low + 1
    if limit is not None:
        assert x_idx_high < limit
    coef_low = x_idx_high - x
    coef_high = x - x_idx_low
    return x_idx_low, x_idx_high, coef_low, coef_high

@q.timer
def proj_meson_jj(res_arr, res_meson_corr, x_rel, total_site):
    r = get_r(x_rel)
    r_limit = get_r_limit(total_site)
    r_idx_low, r_idx_high, coef_low, coef_high = get_interp_idx_coef(r, r_limit)
    x_rel_sym = [ rel_mod_sym(x_rel[mu], total_site[mu]) for mu in range(4) ]
    counts = np.array([ coef_low, coef_high, ]) # r
    v = np.array([ proj(res_arr, x_rel_sym) for proj in all_jj_projections ]) # idx_proj, idx_tensor
    values = np.array([ coef_low * v, coef_high * v, ]).transpose(2, 0, 1) # idx_tensor, r, idx_proj
    v = np.array([ res_meson_corr[idx_meson] for idx_meson in range(3) ]) # r, idx_meson
    values_meson_corr = np.array([ coef_low * v, coef_high * v, ]).transpose(1, 0) # idx_meson, r
    return counts, values, values_meson_corr

def accumulate_proj_meson_jj(counts, values, values_meson_corr, proj_acc, x_rel, total_site):
    # counts[t, r_idx]
    # values[idx_tensor, t, r_idx, idx_proj,]
    # values_meson_corr[idx_meson, t, r_idx]
    # 0 <= idx_proj < len(all_jj_projections): mm, tt, ii, xx
    # 0 <= idx_tensor < n_tensor. ( n_tensor = 28 )
    # 0 <= idx_meson < 3
    # 0 <= t < total_site[3]
    # 0 <= r < r_limit (scale by factor of 5.0)
    (n_tensor, t_size, r_limit, n_proj,) = values.shape
    assert (t_size, r_limit,) == counts.shape
    assert t_size == total_site[3]
    assert r_limit == get_r_limit(total_site)
    assert len(all_jj_projections) == n_proj
    counts_acc, values_acc, values_meson_corr_acc = proj_acc
    r = get_r(x_rel)
    r_idx_low, r_idx_high, coef_low, coef_high = get_interp_idx_coef(r, r_limit)
    t = x_rel[3] % total_site[3]
    counts[t, r_idx_low:r_idx_high + 1] += counts_acc
    values[:, t, r_idx_low:r_idx_high + 1, :] += values_acc
    values_meson_corr[:, t, r_idx_low:r_idx_high + 1] += values_meson_corr_acc

@q.timer
def accumulate_meson_jj(counts, values, values_meson_corr, res_arr, res_meson_corr, x_rel, total_site):
    proj_acc = proj_meson_jj(res_arr, res_meson_corr, x_rel, total_site)
    accumulate_proj_meson_jj(counts, values, values_meson_corr, proj_acc, x_rel, total_site)

@q.timer_verbose
def auto_contract_meson_jj(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"auto-contract/{job_tag}/traj={traj}/meson_jj.lat"
    fn_counts = f"auto-contract/{job_tag}/traj={traj}/meson_jj_counts.lat"
    fn_meson_corr = f"auto-contract/{job_tag}/traj={traj}/meson_jj_meson_corr.lat"
    if get_load_path(fn) is not None:
        assert get_load_path(fn_counts) is not None
        assert get_load_path(fn_meson_corr) is not None
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
    n_proj = len(all_jj_projections)
    assert n_proj == len(all_jj_projection_names)
    n_tensor = (len(expr_names) - 3) // 16
    assert n_tensor * 16 + 3 == len(expr_names)
    def load_data():
        for idx, xg_src in enumerate(xg_psel_list):
            xg_src = xg_src.tolist()
            q.displayln_info(f"auto_contract_meson_jj: {idx+1}/{len(xg_psel_list)} {xg_src}")
            for xg_snk in xg_fsel_list:
                xg_snk = xg_snk.tolist()
                x_rel = [ rel_mod(xg_snk[mu] - xg_src[mu], total_site[mu]) for mu in range(4) ]
                x_rel_t = x_rel[3]
                x_2_t = xg_src[3]
                x_1_t = x_2_t + x_rel_t
                t_1 = (max(x_1_t, x_2_t) + tsep) % total_site[3]
                t_2 = (min(x_1_t, x_2_t) - tsep) % total_site[3]
                pd = {
                        "x_1" : ("point-snk", xg_snk,),
                        "x_2" : ("point", xg_src,),
                        "t_1" : ("wall", t_1),
                        "t_2" : ("wall", t_2),
                        }
                props = eval_cexpr_get_props(cexpr, positions_dict = pd, get_prop = get_prop)
                yield props, x_rel, total_site
    @q.timer
    def feval(args):
        props, x_rel, total_site = args
        val = eval_cexpr_eval(cexpr, props = props)
        assert val.shape[0] == 16 * n_tensor + 3
        res_arr = val[:-3].reshape((n_tensor, 4, 4))
        res_meson_corr = val[-3:]
        proj_acc = proj_meson_jj(res_arr, res_meson_corr, x_rel, total_site)
        return proj_acc, x_rel
    def sum_function(val_list):
        counts = np.zeros((t_size, r_limit,), dtype = complex)
        values = np.zeros((n_tensor, t_size, r_limit, n_proj,), dtype = complex)
        values_meson_corr = np.zeros((3, t_size, r_limit), dtype = complex)
        for proj_acc, x_rel in val_list:
            accumulate_proj_meson_jj(counts, values, values_meson_corr, proj_acc, x_rel, total_site)
        return counts, values, values_meson_corr
    gc.collect()
    q.timer_fork(0)
    res_count, res_sum, res_meson_corr_sum = q.glb_sum_list(
            q.parallel_map_sum(q.get_q_mp_proc(), feval, load_data(), sum_function = sum_function, chunksize = 16))
    q.displayln_info("timer_display for auto_contract_meson_jj")
    q.timer_display()
    q.timer_merge()
    res_count *= 1.0 / (len(xg_psel_list) * fsel.prob())
    res_sum *= 1.0 / (len(xg_psel_list) * fsel.prob())
    res_meson_corr_sum *= 1.0 / (len(xg_psel_list) * fsel.prob())
    ld_count = q.mk_lat_data([
        [ "t", t_size, [ str(rel_mod(t, t_size)) for t in range(t_size) ], ],
        [ "r", r_limit, ],
        ])
    ld_sum = q.mk_lat_data([
        [ "idx_tensor", n_tensor, ],
        [ "t", t_size, [ str(rel_mod(t, t_size)) for t in range(t_size) ], ],
        [ "r", r_limit, ],
        [ "idx_proj", len(all_jj_projection_names), all_jj_projection_names, ],
        ])
    ld_meson_corr_sum = q.mk_lat_data([
        [ "idx_meson", 3, [ "pi", "kp", "km", ], ],
        [ "t", t_size, [ str(rel_mod(t, t_size)) for t in range(t_size) ], ],
        [ "r", r_limit, ],
        ])
    ld_count.from_numpy(res_count)
    ld_sum.from_numpy(res_sum)
    ld_meson_corr_sum.from_numpy(res_meson_corr_sum)
    ld_count.save(get_save_path(fn_counts))
    ld_sum.save(get_save_path(fn))
    ld_meson_corr_sum.save(get_save_path(fn_meson_corr))
    # q.displayln_info(ld_count.show())
    # q.displayln_info(ld_sum.show())
    # q.displayln_info(ld_meson_corr_sum.show())

### ------

@q.timer_verbose
def run_job(job_tag, traj):
    fns_produce = [
            f"auto-contract/{job_tag}/traj={traj}/checkpoint.txt",
            ]
    fns_need = [
            # (f"configs/{job_tag}/ckpoint_lat.{traj}", f"configs/{job_tag}/ckpoint_lat.IEEE64BIG.{traj}",),
            f"point-selection/{job_tag}/traj={traj}.txt",
            f"field-selection/{job_tag}/traj={traj}.field",
            f"gauge-transform/{job_tag}/traj={traj}.field",
            f"wall-src-info-light/{job_tag}/traj={traj}.txt",
            f"wall-src-info-strange/{job_tag}/traj={traj}.txt",
            f"psel-prop-wsrc-light/{job_tag}/traj={traj}/checkpoint.txt",
            f"psel-prop-wsrc-strange/{job_tag}/traj={traj}/checkpoint.txt",
            f"psel-prop-psrc-light/{job_tag}/traj={traj}/checkpoint.txt",
            f"psel-prop-psrc-strange/{job_tag}/traj={traj}/checkpoint.txt",
            f"prop-wsrc-light/{job_tag}/traj={traj}/geon-info.txt",
            f"prop-wsrc-strange/{job_tag}/traj={traj}/geon-info.txt",
            f"prop-psrc-light/{job_tag}/traj={traj}/geon-info.txt",
            f"prop-psrc-strange/{job_tag}/traj={traj}/geon-info.txt",
            # f"prop-rand-u1-light/{job_tag}/traj={traj}/geon-info.txt",
            # f"prop-rand-u1-strange/{job_tag}/traj={traj}/geon-info.txt",
            # f"prop-rand-u1-charm/{job_tag}/traj={traj}/geon-info.txt",
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
            )
    #
    fn_checkpoint = f"auto-contract/{job_tag}/traj={traj}/checkpoint.txt"
    if get_load_path(fn_checkpoint) is None:
        if q.obtain_lock(f"locks/{job_tag}-{traj}-auto-contract"):
            q.timer_fork()
            get_prop = get_get_prop()
            # ADJUST ME
            auto_contract_meson_corr(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contract_meson_corr_psnk(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contract_meson_corr_psrc(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contract_meson_corr_psnk_psrc(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contract_meson_f_corr_psnk(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contract_meson_f_corr_psrc(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contract_meson_f_corr_psnk_psrc(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contract_meson_m(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contract_meson_jt(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contract_meson_jj(job_tag, traj, get_prop, get_psel, get_fsel)
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
    benchmark_eval_cexpr(get_cexpr_meson_f_corr())
    benchmark_eval_cexpr(get_cexpr_meson_m())
    benchmark_eval_cexpr(get_cexpr_meson_jt())
    benchmark_eval_cexpr(get_cexpr_meson_jj())

def test():
    q.qremove_all_info("locks")
    q.qremove_all_info("results")
    run_job("test-4nt8", 1000)
    # run_job("test-4nt16", 1000)
    # run_job("16IH2", 1000)

size_node_list = [
        [1, 1, 1, 1],
        [1, 1, 1, 2],
        [1, 1, 1, 4],
        [1, 1, 1, 8],
        [2, 2, 2, 2],
        [2, 2, 2, 4],
        [2, 2, 2, 8],
        [2, 2, 4, 8],
        [2, 4, 4, 8],
        [4, 4, 4, 8],
        [4, 4, 4, 16],
        [4, 4, 8, 16],
        [4, 8, 8, 16],
        [8, 8, 8, 16],
        ]

q.begin(sys.argv, size_node_list)

# ADJUST ME
q.qremove_all_info("cache")
get_all_cexpr()
test()

# ADJUST ME
job_tags = [
        # "test-4nt8", "test-4nt16",
        # "32IH1",
        # "32IH2",
        # "24IH1",
        # "24IH2",
        # "24IH3",
        # "32D",
        # "32Dfine",
        # "24D",
        # "24DH",
        # "16IH2",
        # "32IfineH",
        # "32IcoarseH1",
        # "48I",
        ]

q.check_time_limit()

for job_tag in job_tags:
    q.displayln_info(pprint.pformat(rup.dict_params[job_tag]))
    for traj in rup.dict_params[job_tag]["trajs"]:
        run_job(job_tag, traj)

q.timer_display()

q.end()
