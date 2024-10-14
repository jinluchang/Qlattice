#!/usr/bin/env python3

json_results = []

from auto_contractor.operators import *

import functools
import math
import os
import time
import importlib
import sys
import math
import cmath

import qlat_gpt as qg

from qlat_scripts.v1 import *

is_cython = False

# ----

load_path_list[:] = [
        "results",
        "qcddata",
        "qcddata-1",
        "qcddata-2",
        "qcddata-3",
        "qcddata-4",
        "qcddata-5",
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
                mk_pi_p("x_2", True)    * mk_j5pi_mu("x_1", 3, True) + f"pi+^dag(0) * j5pi_t^dag(-tsep)",
                mk_k_p("x_2", True)     * mk_j5k_mu("x_1", 3, True)  + f"K+^dag(0) * j5k_t^dag(-tsep)",
                mk_j5pi_mu("x_2", 3)    * mk_pi_p("x_1")             + f"j5pi_t(0) * pi+(-tsep)",
                mk_j5k_mu("x_2", 3)     * mk_k_p("x_1")              + f"j5k_t(0) * K+(-tsep)",
                mk_j5pi_mu("x_2", 3)    * mk_j5pi_mu("x_1", 3, True) + f"j5pi_t(0) * j5pi_t^dag(-tsep)",
                mk_j5k_mu("x_2", 3)     * mk_j5k_mu("x_1", 3, True)  + f"j5k_t(0) * j5k_t^dag(-tsep)",
                mk_a0_p("x_2", True)    * mk_a0_p("x_1")             + f"a0+^dag(0) * a0+(-tsep)",
                mk_kappa_p("x_2", True) * mk_kappa_p("x_1")          + f"kappa+^dag(0) * kappa+(-tsep)",
                mk_j_mu("x_2", 3)       * mk_j_mu("x_1", 3)          + f"j_t(0) * j_t(-tsep)",
                sum([ mk_j_mu("x_2", mu) * mk_j_mu("x_1", mu) for mu in range(4) ])
                + f"j_mu(0) * j_mu(-tsep)",
                sum([ mk_j5pi_mu("x_2", mu) * mk_j5pi_mu("x_1", mu, True) for mu in range(4) ])
                + f"j5pi_mu(0) * j5pi_mu^dag(-tsep)",
                sum([ mk_j5k_mu("x_2", mu) * mk_j5k_mu("x_1", mu, True) for mu in range(4) ])
                + f"j5k_mu(0) * j5k_mu^dag(-tsep)",
                ]
        cexpr = contract_simplify_compile(
                *exprs,
                is_isospin_symmetric_limit=True,
                diagram_type_dict=diagram_type_dict,
                )
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

@q.timer_verbose
def auto_contract_meson_corr(job_tag, traj, get_get_prop):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_corr.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_corr()
    expr_names = get_expr_names(cexpr)
    get_prop = get_get_prop()
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    t_size = total_site[3]
    geo = q.Geometry(total_site)
    total_volume = geo.total_volume
    def load_data():
        t_t_list = q.get_mpi_chunk(
                [ (t_src, t_snk,) for t_snk in range(t_size) for t_src in range(t_size) ],
                rng_state=None)
        for t_src, t_snk in t_t_list:
            yield t_src, t_snk
    @q.timer
    def feval(args):
        t_src, t_snk = args
        t = (t_snk - t_src) % t_size
        pd = {
                "x_2": ("wall", t_snk,),
                "x_1": ("wall", t_src,),
                "size": total_site,
                }
        val = eval_cexpr(cexpr, positions_dict=pd, get_prop=get_prop)
        return val, t
    def sum_function(val_list):
        values = np.zeros((t_size, len(expr_names),), dtype=complex)
        for val, t in val_list:
            values[t] += val
        return q.glb_sum(values.transpose(1, 0))
    chunk_size = get_param(job_tag, "measurement", "auto_contractor_chunk_size", default=128)
    q.timer_fork(0)
    res_sum = q.parallel_map_sum(feval, load_data(), sum_function=sum_function, chunksize=chunk_size)
    q.displayln_info(f"{fname}: timer_display for parallel_map_sum")
    q.timer_display()
    q.timer_merge()
    res_sum *= 1.0 / t_size
    assert q.qnorm(res_sum[0] - 1.0) < 1e-10
    ld = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t_sep", t_size, [ str(q.rel_mod(t, t_size)) for t in range(t_size) ], ],
        ])
    ld.from_numpy(res_sum)
    ld.save(get_save_path(fn))
    json_results.append((f"{fname}: ld sig", q.get_data_sig(ld, q.RngState()),))
    for i, en in enumerate(expr_names):
        json_results.append((f"{fname}: ld '{en}' sig", q.get_data_sig(ld[i], q.RngState()),))

@q.timer_verbose
def auto_contract_meson_corr_psnk(job_tag, traj, get_get_prop):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_corr_psnk.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_corr()
    expr_names = get_expr_names(cexpr)
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    t_size = total_site[3]
    geo = q.Geometry(total_site)
    total_volume = geo.total_volume
    get_prop = get_get_prop()
    xg_list = get_all_points(total_site)
    xg_local_list = [ q.Coordinate(xg) for xg in geo.xg_arr() ]
    def load_data():
        for t_src in range(t_size):
            for xg_snk in xg_local_list:
                yield t_src, xg_snk
    @q.timer
    def feval(args):
        t_src, xg_snk = args
        t = (xg_snk[3] - t_src) % t_size
        pd = {
                "x_2" : ("point", xg_snk,),
                "x_1" : ("wall", t_src,),
                "size" : total_site,
                }
        val = eval_cexpr(cexpr, positions_dict=pd, get_prop=get_prop)
        return val, t
    def sum_function(val_list):
        values = np.zeros((t_size, len(expr_names),), dtype=complex)
        for val, t in val_list:
            values[t] += val
        return q.glb_sum(values.transpose(1, 0))
    chunk_size = get_param(job_tag, "measurement", "auto_contractor_chunk_size", default=128)
    q.timer_fork(0)
    res_sum = q.parallel_map_sum(feval, load_data(), sum_function=sum_function, chunksize=chunk_size)
    q.displayln_info(f"{fname}: timer_display for parallel_map_sum")
    q.timer_display()
    q.timer_merge()
    res_sum *= 1.0 / total_volume
    assert q.qnorm(res_sum[0] - 1.0) < 1e-10
    ld = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t_sep", t_size, [ str(q.rel_mod(t, t_size)) for t in range(t_size) ], ],
        ])
    ld.from_numpy(res_sum)
    ld.save(get_save_path(fn))
    json_results.append((f"{fname}: ld sig", q.get_data_sig(ld, q.RngState()),))
    for i, en in enumerate(expr_names):
        json_results.append((f"{fname}: ld '{en}' sig", q.get_data_sig(ld[i], q.RngState()),))

@q.timer_verbose
def auto_contract_meson_corr_psnk_psrc(job_tag, traj, get_get_prop):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_corr_psnk_psrc.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_corr()
    expr_names = get_expr_names(cexpr)
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    t_size = total_site[3]
    geo = q.Geometry(total_site)
    total_volume = geo.total_volume
    get_prop = get_get_prop()
    xg_list = get_all_points(total_site)
    xg_local_list = [ q.Coordinate(xg) for xg in geo.xg_arr() ]
    r_list = get_r_list(job_tag)
    r_sq_interp_idx_coef_list = get_r_sq_interp_idx_coef_list(job_tag)
    def load_data():
        for xg_src in xg_local_list:
            yield xg_src
    @q.timer
    def feval(args):
        xg_src = args
        res_list = []
        for xg_snk in xg_list:
            x_rel = q.smod_coordinate(xg_snk - xg_src, total_site)
            r_sq = x_rel.r_sqr()
            r_idx_low, r_idx_high, coef_low, coef_high = r_sq_interp_idx_coef_list[r_sq]
            t = x_rel[3]
            pd = {
                    "x_2": ("point", xg_snk,),
                    "x_1": ("point", xg_src,),
                    "size": total_site,
                    }
            val = eval_cexpr(cexpr, positions_dict=pd, get_prop=get_prop)
            res_list.append((val, t, r_idx_low, r_idx_high, coef_low, coef_high,))
        return res_list
    def sum_function(val_list):
        values = np.zeros((t_size, len(r_list), len(expr_names),), dtype=complex)
        for idx, res_list in enumerate(val_list):
            for val, t, r_idx_low, r_idx_high, coef_low, coef_high in res_list:
                values[t, r_idx_low] += coef_low * val
                values[t, r_idx_high] += coef_high * val
            q.displayln_info(f"{fname}: {idx+1}/{len(xg_local_list)}")
        return q.glb_sum(values.transpose(2, 0, 1))
    q.timer_fork(0)
    res_sum = q.parallel_map_sum(feval, load_data(), sum_function=sum_function, chunksize=1)
    q.displayln_info(f"{fname}: timer_display for parallel_map_sum")
    q.timer_display()
    q.timer_merge()
    res_sum *= 1.0 / (total_volume * total_volume / t_size)
    assert q.qnorm(res_sum[0].sum(1) - 1.0) < 1e-10
    ld = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t_sep", t_size, [ str(q.rel_mod(t, t_size)) for t in range(t_size) ], ],
        [ "r", len(r_list), [ f"{r:.5f}" for r in r_list ], ],
        ])
    ld.from_numpy(res_sum)
    ld.save(get_save_path(fn))
    json_results.append((f"{fname}: ld sig", q.get_data_sig(ld, q.RngState()),))
    for i, en in enumerate(expr_names):
        json_results.append((f"{fname}: ld '{en}' sig", q.get_data_sig(ld[i], q.RngState()),))

@q.timer_verbose
def auto_contract_meson_corr_psnk_psrc_rand(job_tag, traj, get_get_prop):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_corr_psnk_psrc_rand.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_corr()
    expr_names = get_expr_names(cexpr)
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    t_size = total_site[3]
    geo = q.Geometry(total_site)
    total_volume = geo.total_volume
    get_prop = get_get_prop()
    xg_list = get_all_points(total_site)
    xg_local_list = [ q.Coordinate(xg) for xg in geo.xg_arr() ]
    r_list = get_r_list(job_tag)
    r_sq_interp_idx_coef_list = get_r_sq_interp_idx_coef_list(job_tag)
    sample_num = get_param(job_tag, "measurement", "auto_contract_meson_corr_psnk_psrc_rand", "sample_num", default=128)
    sample_size = get_param(job_tag, "measurement", "auto_contract_meson_corr_psnk_psrc_rand", "sample_size", default=128)
    rs = q.RngState(f"{job_tag}-{traj}-{fname}")
    mpi_chunk = q.get_mpi_chunk(list(range(sample_num)))
    def load_data():
        for idx in mpi_chunk:
            yield idx
    @q.timer
    def feval(args):
        idx = args
        res_list = []
        for idx2 in range(sample_size):
            rsi = rs.split(f"{idx}-{idx2}")
            xg_src = rsi.c_rand_gen(total_site)
            xg_snk = rsi.c_rand_gen(total_site)
            x_rel = q.smod_coordinate(xg_snk - xg_src, total_site)
            r_sq = x_rel.r_sqr()
            r_idx_low, r_idx_high, coef_low, coef_high = r_sq_interp_idx_coef_list[r_sq]
            t = x_rel[3]
            pd = {
                    "x_2": ("point", xg_snk,),
                    "x_1": ("point", xg_src,),
                    "size": total_site,
                    }
            val = eval_cexpr(cexpr, positions_dict=pd, get_prop=get_prop)
            res_list.append((val, t, r_idx_low, r_idx_high, coef_low, coef_high,))
        return res_list
    def sum_function(val_list):
        values = np.zeros((t_size, len(r_list), len(expr_names),), dtype=complex)
        for idx, res_list in enumerate(val_list):
            for val, t, r_idx_low, r_idx_high, coef_low, coef_high in res_list:
                values[t, r_idx_low] += coef_low * val
                values[t, r_idx_high] += coef_high * val
            q.displayln_info(f"{fname}: {idx+1}/{len(mpi_chunk)}")
        return q.glb_sum(values.transpose(2, 0, 1))
    q.timer_fork(0)
    res_sum = q.parallel_map_sum(feval, load_data(), sum_function=sum_function, chunksize=1)
    q.displayln_info(f"{fname}: timer_display for parallel_map_sum")
    q.timer_display()
    q.timer_merge()
    res_sum *= 1.0 / (sample_num * sample_size / t_size)
    # assert q.qnorm(res_sum[0].sum(1) - 1.0) < 1e-10
    ld = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t_sep", t_size, [ str(q.rel_mod(t, t_size)) for t in range(t_size) ], ],
        [ "r", len(r_list), [ f"{r:.5f}" for r in r_list ], ],
        ])
    ld.from_numpy(res_sum)
    ld.save(get_save_path(fn))
    json_results.append((f"{fname}: ld sig", q.get_data_sig(ld, q.RngState()),))
    for i, en in enumerate(expr_names):
        json_results.append((f"{fname}: ld '{en}' sig", q.get_data_sig(ld[i], q.RngState()),))

# ----

all_mom_list_dict = {
        0: [
            q.Coordinate([ 0, 0, 0, 0, ]),
            ],
        1: [
            q.Coordinate([ 0, 0, 1, 0, ]),
            q.Coordinate([ 0, 1, 0, 0, ]),
            q.Coordinate([ 1, 0, 0, 0, ]),
            q.Coordinate([ 0, 0, -1, 0, ]),
            q.Coordinate([ 0, -1, 0, 0, ]),
            q.Coordinate([ -1, 0, 0, 0, ]),
            ],
        2: [
            q.Coordinate([ 0, 1, 1, 0, ]),
            q.Coordinate([ 1, 0, 1, 0, ]),
            q.Coordinate([ 1, 1, 0, 0, ]),
            q.Coordinate([ 0, -1, 1, 0, ]),
            q.Coordinate([ -1, 0, 1, 0, ]),
            q.Coordinate([ -1, 1, 0, 0, ]),
            q.Coordinate([ 0, 1, -1, 0, ]),
            q.Coordinate([ 1, 0, -1, 0, ]),
            q.Coordinate([ 1, -1, 0, 0, ]),
            q.Coordinate([ 0, -1, -1, 0, ]),
            q.Coordinate([ -1, 0, -1, 0, ]),
            q.Coordinate([ -1, -1, 0, 0, ]),
            ],
        3: [
            q.Coordinate([ 1, 1, 1, 0, ]),
            q.Coordinate([ -1, 1, 1, 0, ]),
            q.Coordinate([ 1, -1, 1, 0, ]),
            q.Coordinate([ -1, -1, 1, 0, ]),
            q.Coordinate([ 1, 1, -1, 0, ]),
            q.Coordinate([ -1, 1, -1, 0, ]),
            q.Coordinate([ 1, -1, -1, 0, ]),
            q.Coordinate([ -1, -1, -1, 0, ]),
            ],
        4: [
            q.Coordinate([ 0, 0, 2, 0, ]),
            q.Coordinate([ 0, 2, 0, 0, ]),
            q.Coordinate([ 2, 0, 0, 0, ]),
            q.Coordinate([ 0, 0, -2, 0, ]),
            q.Coordinate([ 0, -2, 0, 0, ]),
            q.Coordinate([ -2, 0, 0, 0, ]),
            ],
        }

def get_mom_list(mom_idx):
    return all_mom_list_dict[mom_idx]

# ----

def wave_function(p1, p2, radius, size):
    p1_tag, c1 = p1
    p2_tag, c2 = p2
    c12 = q.smod_coordinate(c1 - c2, size)
    # assert c12[3] == 0
    c12_r_sqr = c12.r_sqr()
    dis = math.sqrt(c12_r_sqr)
    vol = size[0] * size[1] * size[2]
    wf = vol**2 * math.exp(-dis / radius)
    return wf

def momentum_factor(mom, p, size, is_dagger=False):
    p_tag, c = p
    assert mom[3] == 0
    phase = mom[0] * c[0] / size[0] + mom[1] * c[1] / size[1] + mom[2] * c[2] / size[2]
    phase = phase * 2.0 * math.pi
    if not is_dagger:
        mf = cmath.rect(1.0, phase)
    else:
        mf = cmath.rect(1.0, -phase)
    return mf

def mk_meson_wf(f1, f2, p1, p2, radius, mom, is_dagger=False):
    """
    i q1bar g5 q2 #dag: i q2bar g5 q1
    return the actual dagger of the operator
    """
    s1 = new_spin_index()
    s2 = new_spin_index()
    c = new_color_index()
    g5 = G(5, s1, s2)
    wf = mk_fac(f"wave_function({p1},{p2},{radius},size)")
    if not is_dagger:
        q1b = Qb(f1, p1, s1, c)
        q2v = Qv(f2, p2, s2, c)
        mf = mk_fac(f"momentum_factor({mom},{p2},size)")
        return sympy.I * wf * mf * q1b * g5 * q2v + f"(i {f1}bar g5 {f2})({p1},{p2})"
    else:
        q2b = Qb(f2, p2, s1, c)
        q1v = Qv(f1, p1, s2, c)
        mf = mk_fac(f"momentum_factor(-{mom},{p2},size)")
        return sympy.I * wf * mf * q2b * g5 * q1v + f"(i {f2}bar g5 {f1})({p2},{p1},{radius},{mom})"

def mk_scalar_meson_wf(f1, f2, p1, p2, radius, mom, is_dagger=False):
    """
    i q1bar g5 q2 #dag: i q2bar g5 q1
    return the actual dagger of the operator
    """
    s = new_spin_index()
    c = new_color_index()
    wf = mk_fac(f"wave_function({p1},{p2},{radius},size)")
    if not is_dagger:
        q1b = Qb(f1, p1, s, c)
        q2v = Qv(f2, p2, s, c)
        mf = mk_fac(f"momentum_factor({mom},{p2},size)")
        return wf * mf * q1b * q2v + f"({f1}bar {f2})({p1},{p2})"
    else:
        q2b = Qb(f2, p2, s, c)
        q1v = Qv(f1, p1, s, c)
        mf = mk_fac(f"momentum_factor(-{mom},{p2},size)")
        return wf * mf * q2b * q1v + f"({f2}bar {f1})({p2},{p1},{radius},{mom})"

# ----

def mk_pi_0_wf(p1, p2, mom, is_dagger=False):
    """
    i/sqrt(2) * (ubar g5 u - dbar g5 d)  #dag: same
    """
    radius = "r_pi"
    return 1 / sympy.sqrt(2) * (
            mk_meson_wf("u", "u", p1, p2, radius, mom, is_dagger)
            - mk_meson_wf("d", "d", p1, p2, radius, mom, is_dagger)
            ) + f"pi0({p1},{p2},{radius},{mom}){show_dagger(is_dagger)}"

def mk_pi_p_wf(p1, p2, mom, is_dagger=False):
    """
    i ubar g5 d  #dag: i dbar g5 u
    """
    radius = "r_pi"
    return (mk_meson_wf("u", "d", p1, p2, radius, mom, is_dagger)
            + f"pi+({p1},{p2},{radius},{mom}){show_dagger(is_dagger)}")

def mk_pi_m_wf(p1, p2, mom, is_dagger=False):
    """
    -i dbar g5 u  #dag: -i ubar g5 d
    """
    radius = "r_pi"
    return (-mk_meson_wf("d", "u", p1, p2, radius, mom, is_dagger)
            + f"pi-({p1},{p2},{radius},{mom}){show_dagger(is_dagger)}")

def mk_sigma_wf(p1, p2, mom, is_dagger=False):
    """
    1/sqrt(2) * (ubar u + dbar d)
    """
    radius = "r_sigma"
    return 1 / sympy.sqrt(2) * (
            mk_scalar_meson_wf("u", "u", p1, p2, radius, mom, is_dagger)
            + mk_scalar_meson_wf("d", "d", p1, p2, radius, mom, is_dagger)
            ) + f"sigma({p1},{p2},{radius},{mom}){show_dagger(is_dagger)}"

def mk_kk_p_wf(p1, p2, mom, is_dagger=False):
    """
    i ubar g5 s  #dag:  i sbar g5 u
    """
    radius = "r_kk"
    return (mk_meson_wf("u", "s", p1, p2, radius, mom, is_dagger)
            + f"K+({p1},{p2},{radius},{mom}){show_dagger(is_dagger)}")

def mk_kk_m_wf(p1, p2, mom, is_dagger=False):
    """
    -i sbar g5 u  #dag: -i ubar g5 s
    """
    radius = "r_kk"
    return (-mk_meson_wf("s", "u", p1, p2, radius, mom, is_dagger)
            + f"K-({p1},{p2},{radius},{mom}){show_dagger(is_dagger)}")

def mk_kk_0_wf(p1, p2, mom, is_dagger=False):
    """
    i dbar g5 s  #dag: i sbar g5 d
    """
    radius = "r_kk"
    return (mk_meson_wf("d", "s", p1, p2, radius, mom, is_dagger)
            + f"K0({p1},{p2},{radius},{mom}){show_dagger(is_dagger)}")

def mk_kk_0_bar_wf(p1, p2, mom, is_dagger=False):
    """
    -i sbar g5 d  #dag: -i dbar g5 s
    """
    radius = "r_kk"
    return (-mk_meson_wf("s", "d", p1, p2, radius, mom, is_dagger)
            + f"K0b({p1},{p2},{radius},{mom}){show_dagger(is_dagger)}")

# ----

def mk_pipi_i0_j0_wf(p11, p12, p21, p22, mom, is_dagger=False):
    """
    note that two pion may not at the same time slice
    mom is for individual pion
    total momentum is zero
    """
    def mk_op_mom(mom_pi):
        mom1 = mom_pi
        mom2 = -mom_pi
        return 1 / sympy.sqrt(3) * (
                - mk_pi_0_wf(p11, p12, mom1, is_dagger) * mk_pi_0_wf(p21, p22, mom2, is_dagger)
                + mk_pi_m_wf(p11, p12, mom1, is_dagger) * mk_pi_p_wf(p21, p22, mom2, is_dagger)
                + mk_pi_p_wf(p11, p12, mom1, is_dagger) * mk_pi_m_wf(p21, p22, mom2, is_dagger)
                )
    mom_pi_list = get_mom_list(mom)
    n_op = 0
    expr = 0
    for mom_pi in mom_pi_list:
        n_op += 1
        expr = expr + mk_op_mom(mom_pi)
    return 1 / sympy.sqrt(n_op) * expr + f"pipi_I0({p11},{p12},{p21},{p22},mom={mom}){show_dagger(is_dagger)}"

def mk_kkkk_i0_j0_wf(p11, p12, p21, p22, mom, is_dagger=False, *, is_sym=True):
    """
    note that two kaon may not at the same time slice
    symmetrization is used by default (different from the library version)
    mom is for individual kaon
    total momentum is zero
    """
    if is_sym:
        return 1 / sympy.sqrt(2) * (
                mk_kkkk_i0_j0_wf(p11, p12, p21, p22, mom, is_dagger, is_sym=False)
                + mk_kkkk_i0_j0_wf(p21, p22, p11, p12, mom, is_dagger, is_sym=False)
                ) + f"KK_I0({p11},{p12},{p21},{p22},mom={mom},sym){show_dagger(is_dagger)}"
    def mk_op_mom(mom_kk):
        mom1 = mom_kk
        mom2 = -mom_kk
        return 1 / sympy.sqrt(2) * (
                mk_kk_0_wf(p11, p12, mom1, is_dagger) * mk_kk_0_bar_wf(p21, p22, mom2, is_dagger)
                + mk_kk_p_wf(p11, p12, mom1, is_dagger) * mk_kk_m_wf(p21, p22, mom2, is_dagger)
                )
    mom_kk_list = get_mom_list(mom)
    n_op = 0
    expr = 0
    for mom_kk in mom_kk_list:
        n_op += 1
        expr = expr + mk_op_mom(mom_kk)
    return 1 / sympy.sqrt(n_op) * expr + f"KK_I0({p11},{p12},{p21},{p22},mom={mom}){show_dagger(is_dagger)}"

# ----

@q.timer
def get_cexpr_meson_corr_wf():
    fname = q.get_fname()
    fn_base = f"cache/auto_contract_cexpr/{fname}"
    def calc_cexpr():
        diagram_type_dict = dict()
        def get_mom_avg_expr_list(f):
            """
            f(mom) -> expr
            """
            expr_list = []
            for mom_idx, mom_list in all_mom_list_dict.items():
                expr = 0
                fac = 1 / mk_sym(len(mom_list))
                desp = ""
                for mom in mom_list:
                    e = f(mom)
                    if desp == "":
                        desp = e.description
                    expr += fac * e
                expr = expr + desp
                expr_list.append(expr)
            return expr_list
        exprs = [ mk_fac(1) + f"1", ]
        exprs += get_mom_avg_expr_list(
                lambda mom:
                (-mk_pi_m_wf("x2_1", "x2_2", -mom) * mk_pi_p_wf("x1_1", "x1_2", mom))
                )
        exprs += get_mom_avg_expr_list(
                lambda mom:
                (-mk_kk_m_wf("x2_1", "x2_2", -mom) * mk_kk_p_wf("x1_1", "x1_2", mom))
                )
        exprs += get_mom_avg_expr_list(
                lambda mom:
                (mk_sigma_wf("x2_1", "x2_2", -mom) * mk_sigma_wf("x1_1", "x1_2", mom))
                )
        cexpr = contract_simplify_compile(
                *exprs,
                is_isospin_symmetric_limit=True,
                diagram_type_dict=diagram_type_dict,
                )
        return cexpr
    base_positions_dict = {
            "wave_function": wave_function,
            "momentum_factor": momentum_factor,
            "r_pi": 1.5,
            "r_sigma": 1.5,
            "r_kk": 2.0,
            "Coordinate": q.Coordinate,
            }
    return cache_compiled_cexpr(calc_cexpr, fn_base, base_positions_dict=base_positions_dict, is_cython=is_cython)

@q.timer_verbose
def auto_contract_meson_corr_wf(job_tag, traj, get_get_prop):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract/traj-{traj}/{fname}.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_corr_wf()
    expr_names = get_expr_names(cexpr)
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    t_size = total_site[3]
    geo = q.Geometry(total_site)
    total_volume = geo.total_volume
    get_prop = get_get_prop()
    xg_list = get_all_points(total_site)
    xg_local_list = [ q.Coordinate(xg) for xg in geo.xg_arr() ]
    sample_num = get_param(job_tag, "measurement", fname,
                           "sample_num", default=512)
    sample_size = get_param(job_tag, "measurement", fname,
                            "sample_size", default=128)
    t_sep_range = get_param(job_tag, "measurement", fname,
                            "t_sep_range", default=17)
    t_sep_range = min(t_size, t_sep_range)
    rs = q.RngState(f"{job_tag}-{traj}-{fname}")
    mpi_chunk = q.get_mpi_chunk(list(range(sample_num)))
    def load_data():
        for idx in mpi_chunk:
            yield idx
    @q.timer
    def feval_single(idx, idx2, t_src, t_sep):
        rsi = rs.split(f"{idx}-{idx2}-{t_src}-{t_sep}")
        x1_1 = rsi.c_rand_gen(total_site)
        x1_2 = rsi.c_rand_gen(total_site)
        x2_1 = rsi.c_rand_gen(total_site)
        x2_2 = rsi.c_rand_gen(total_site)
        x1_1[3] = t_src
        x2_1[3] = (x1_1[3] + t_sep) % t_size
        x1_2[3] = x1_1[3]
        x2_2[3] = x2_1[3]
        pd = {
                "x1_1": ("point", x1_1,),
                "x1_2": ("point", x1_2,),
                "x2_1": ("point", x2_1,),
                "x2_2": ("point", x2_2,),
                "size": total_site,
                }
        val = eval_cexpr(cexpr, positions_dict=pd, get_prop=get_prop)
        return val
    @q.timer
    def feval(args):
        idx = args
        values = np.zeros((t_sep_range, len(expr_names),), dtype=complex)
        for idx2 in range(sample_size):
            for t_src in range(t_size):
                for t_sep in range(t_sep_range):
                    val = feval_single(idx, idx2, t_src, t_sep)
                    values[t_sep] += val
        return values
    def sum_function(val_list):
        values = np.zeros((t_sep_range, len(expr_names),), dtype=complex)
        for idx, val in enumerate(val_list):
            values += val
            q.displayln_info(f"{fname}: {idx+1}/{len(mpi_chunk)}")
        return q.glb_sum(values.transpose(1, 0,))
    q.timer_fork(0)
    res_sum = q.parallel_map_sum(feval, load_data(), sum_function=sum_function, chunksize=1)
    q.displayln_info(f"{fname}: timer_display for parallel_map_sum")
    q.timer_display()
    q.timer_merge()
    res_sum *= 1.0 / (sample_num * sample_size * t_size)
    assert q.qnorm(res_sum[0] - 1.0) < 1e-10
    ld = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t_sep", t_sep_range, [ str(t) for t in range(t_sep_range) ], ],
        ])
    ld.from_numpy(res_sum)
    ld.save(get_save_path(fn))
    json_results.append((f"{fname}: ld sig", q.get_data_sig(ld, q.RngState()),))
    for i, en in enumerate(expr_names):
        json_results.append((f"{fname}: ld '{en}' sig", q.get_data_sig(ld[i], q.RngState()),))

# ----

def time_selector_8(p):
    p_tag, c = p
    t = c[3]
    if t % 8 == 0:
        return 8
    else:
        return 0

@q.timer
def get_cexpr_meson_meson_i0_j0_corr_wf():
    fname = q.get_fname()
    fn_base = f"cache/auto_contract_cexpr/{fname}"
    def calc_cexpr():
        diagram_type_dict = dict()
        #
        diagram_type_dict[()] = 'TypeV'
        diagram_type_dict[((('x1_1_1', 'x1_1_2'), 1),)] = 'TypeV'
        diagram_type_dict[((('x2_1_1', 'x2_1_2'), 1),)] = 'TypeV'
        diagram_type_dict[((('x1_1_1', 'x1_2_3'), 1), (('x1_2_4', 'x1_1_2'), 1))] = 'TypeV'
        diagram_type_dict[((('x1_1_1', 'x1_2_3p'), 1), (('x1_2_4p', 'x1_1_2'), 1))] = 'TypeV'
        diagram_type_dict[((('x2_1_1', 'x2_2_3'), 1), (('x2_2_4', 'x2_1_2'), 1))] = 'TypeV'
        diagram_type_dict[((('x2_1_1', 'x2_2_3p'), 1), (('x2_2_4p', 'x2_1_2'), 1))] = 'TypeV'
        diagram_type_dict[((('x1_1_1', 'x1_1_2'), 1), (('x2_1_1', 'x2_1_2'), 1))] = 'TypeV'
        diagram_type_dict[((('x1_1_1', 'x1_2_3p'), 1), (('x1_2_4p', 'x1_1_2'), 1), (('x2_1_1', 'x2_1_2'), 1))] = 'TypeV'
        diagram_type_dict[((('x1_1_1', 'x1_1_2'), 1), (('x2_1_1', 'x2_2_3'), 1), (('x2_2_4', 'x2_1_2'), 1))] = 'TypeV'
        diagram_type_dict[((('x1_1_1', 'x1_1_2'), 1), (('x2_1_1', 'x2_2_3p'), 1), (('x2_2_4p', 'x2_1_2'), 1))] = 'TypeV'
        diagram_type_dict[((('x1_1_1', 'x1_2_3'), 1), (('x1_2_4', 'x1_1_2'), 1), (('x2_1_1', 'x2_1_2'), 1))] = 'TypeV'
        diagram_type_dict[((('x1_1_1', 'x1_2_3'), 1), (('x1_2_4', 'x1_1_2'), 1), (('x2_1_1', 'x2_2_3'), 1), (('x2_2_4', 'x2_1_2'), 1))] = 'TypeV'
        diagram_type_dict[((('x1_1_1', 'x1_2_3'), 1), (('x1_2_4', 'x1_1_2'), 1), (('x2_1_1', 'x2_2_3p'), 1), (('x2_2_4p', 'x2_1_2'), 1))] = 'TypeV'
        diagram_type_dict[((('x1_1_1', 'x1_2_3p'), 1), (('x1_2_4p', 'x1_1_2'), 1), (('x2_1_1', 'x2_2_3'), 1), (('x2_2_4', 'x2_1_2'), 1))] = 'TypeV'
        diagram_type_dict[((('x1_1_1', 'x1_2_3p'), 1), (('x1_2_4p', 'x1_1_2'), 1), (('x2_1_1', 'x2_2_3p'), 1), (('x2_2_4p', 'x2_1_2'), 1))] = 'TypeV'
        #
        diagram_type_dict[((('x1_1_1', 'x2_2_3'), 1), (('x2_1_1', 'x1_1_2'), 1), (('x2_2_4', 'x2_1_2'), 1))] = 'TypeR'
        diagram_type_dict[((('x1_1_1', 'x2_1_1'), 1), (('x2_1_2', 'x2_2_3'), 1), (('x2_2_4', 'x1_1_2'), 1))] = 'TypeR'
        diagram_type_dict[((('x1_1_1', 'x2_2_3p'), 1), (('x2_1_1', 'x1_1_2'), 1), (('x2_2_4p', 'x2_1_2'), 1))] = 'TypeR'
        diagram_type_dict[((('x1_1_1', 'x2_1_1'), 1), (('x2_1_2', 'x2_2_3p'), 1), (('x2_2_4p', 'x1_1_2'), 1))] = 'TypeR'
        diagram_type_dict[((('x1_1_1', 'x1_2_3'), 1), (('x1_2_4', 'x2_1_1'), 1), (('x2_1_2', 'x1_1_2'), 1))] = 'TypeR'
        diagram_type_dict[((('x1_1_1', 'x2_1_1'), 1), (('x1_2_3', 'x1_1_2'), 1), (('x2_1_2', 'x1_2_4'), 1))] = 'TypeR'
        diagram_type_dict[((('x1_1_1', 'x2_1_1'), 1), (('x1_2_3p', 'x1_1_2'), 1), (('x2_1_2', 'x1_2_4p'), 1))] = 'TypeR'
        diagram_type_dict[((('x1_1_1', 'x1_2_3p'), 1), (('x1_2_4p', 'x2_1_1'), 1), (('x2_1_2', 'x1_1_2'), 1))] = 'TypeR'
        diagram_type_dict[((('x1_1_1', 'x1_2_3'), 1), (('x1_2_4', 'x2_1_1'), 1), (('x2_1_2', 'x2_2_3'), 1), (('x2_2_4', 'x1_1_2'), 1))] = 'TypeR'
        diagram_type_dict[((('x1_1_1', 'x1_2_3'), 1), (('x1_2_4', 'x2_2_3'), 1), (('x2_1_1', 'x1_1_2'), 1), (('x2_2_4', 'x2_1_2'), 1))] = 'TypeR'
        diagram_type_dict[((('x1_1_1', 'x2_1_1'), 1), (('x1_2_3', 'x1_1_2'), 1), (('x2_1_2', 'x2_2_3'), 1), (('x2_2_4', 'x1_2_4'), 1))] = 'TypeR'
        diagram_type_dict[((('x1_1_1', 'x2_2_3'), 1), (('x1_2_3', 'x1_1_2'), 1), (('x2_1_1', 'x1_2_4'), 1), (('x2_2_4', 'x2_1_2'), 1))] = 'TypeR'
        diagram_type_dict[((('x1_1_1', 'x1_2_3'), 1), (('x1_2_4', 'x2_1_1'), 1), (('x2_1_2', 'x2_2_3p'), 1), (('x2_2_4p', 'x1_1_2'), 1))] = 'TypeR'
        diagram_type_dict[((('x1_1_1', 'x1_2_3'), 1), (('x1_2_4', 'x2_2_3p'), 1), (('x2_1_1', 'x1_1_2'), 1), (('x2_2_4p', 'x2_1_2'), 1))] = 'TypeR'
        diagram_type_dict[((('x1_1_1', 'x2_1_1'), 1), (('x1_2_3', 'x1_1_2'), 1), (('x2_1_2', 'x2_2_3p'), 1), (('x2_2_4p', 'x1_2_4'), 1))] = 'TypeR'
        diagram_type_dict[((('x1_1_1', 'x2_2_3p'), 1), (('x1_2_3', 'x1_1_2'), 1), (('x2_1_1', 'x1_2_4'), 1), (('x2_2_4p', 'x2_1_2'), 1))] = 'TypeR'
        diagram_type_dict[((('x1_1_1', 'x2_1_1'), 1), (('x1_2_3p', 'x1_1_2'), 1), (('x2_1_2', 'x2_2_3'), 1), (('x2_2_4', 'x1_2_4p'), 1))] = 'TypeR'
        diagram_type_dict[((('x1_1_1', 'x2_2_3'), 1), (('x1_2_3p', 'x1_1_2'), 1), (('x2_1_1', 'x1_2_4p'), 1), (('x2_2_4', 'x2_1_2'), 1))] = 'TypeR'
        diagram_type_dict[((('x1_1_1', 'x1_2_3p'), 1), (('x1_2_4p', 'x2_1_1'), 1), (('x2_1_2', 'x2_2_3'), 1), (('x2_2_4', 'x1_1_2'), 1))] = 'TypeR'
        diagram_type_dict[((('x1_1_1', 'x1_2_3p'), 1), (('x1_2_4p', 'x2_2_3'), 1), (('x2_1_1', 'x1_1_2'), 1), (('x2_2_4', 'x2_1_2'), 1))] = 'TypeR'
        diagram_type_dict[((('x1_1_1', 'x1_2_3p'), 1), (('x1_2_4p', 'x2_1_1'), 1), (('x2_1_2', 'x2_2_3p'), 1), (('x2_2_4p', 'x1_1_2'), 1))] = 'TypeR'
        diagram_type_dict[((('x1_1_1', 'x1_2_3p'), 1), (('x1_2_4p', 'x2_2_3p'), 1), (('x2_1_1', 'x1_1_2'), 1), (('x2_2_4p', 'x2_1_2'), 1))] = 'TypeR'
        diagram_type_dict[((('x1_1_1', 'x2_1_1'), 1), (('x1_2_3p', 'x1_1_2'), 1), (('x2_1_2', 'x2_2_3p'), 1), (('x2_2_4p', 'x1_2_4p'), 1))] = 'TypeR'
        diagram_type_dict[((('x1_1_1', 'x2_2_3p'), 1), (('x1_2_3p', 'x1_1_2'), 1), (('x2_1_1', 'x1_2_4p'), 1), (('x2_2_4p', 'x2_1_2'), 1))] = 'TypeR'
        #
        diagram_type_dict[((('x1_1_1', 'x2_1_1'), 1), (('x1_2_3', 'x2_2_3'), 1), (('x2_1_2', 'x1_2_4'), 1), (('x2_2_4', 'x1_1_2'), 1))] = 'TypeC'
        diagram_type_dict[((('x1_1_1', 'x2_2_3'), 1), (('x1_2_3', 'x2_1_1'), 1), (('x2_1_2', 'x1_1_2'), 1), (('x2_2_4', 'x1_2_4'), 1))] = 'TypeC'
        #
        diagram_type_dict[((('x1_1_1', 'x2_1_1'), 1), (('x2_1_2', 'x1_1_2'), 1))] = 'TypeD'
        diagram_type_dict[((('x1_1_1', 'x2_1_1'), 1), (('x1_2_3', 'x2_2_3'), 1), (('x2_1_2', 'x1_1_2'), 1), (('x2_2_4', 'x1_2_4'), 1))] = 'TypeD'
        diagram_type_dict[((('x1_1_1', 'x2_1_1'), 1), (('x1_2_3p', 'x2_2_3p'), 1), (('x2_1_2', 'x1_1_2'), 1), (('x2_2_4p', 'x1_2_4p'), 1))] = 'TypeD'
        diagram_type_dict[((('x1_1_1', 'x2_2_3'), 1), (('x1_2_3', 'x2_1_1'), 1), (('x2_1_2', 'x1_2_4'), 1), (('x2_2_4', 'x1_1_2'), 1))] = 'TypeD'
        diagram_type_dict[((('x1_1_1', 'x2_2_3p'), 1), (('x1_2_3p', 'x2_1_1'), 1), (('x2_1_2', 'x1_2_4p'), 1), (('x2_2_4p', 'x1_1_2'), 1))] = 'TypeD'
        #
        mom0 = q.Coordinate([ 0, 0, 0, 0, ])
        src_op_list = [
                mk_sigma_wf("x1_1_1", "x1_1_2", mom0) + f"sigma(x1_1_1,x1_1_2)",
                mk_pipi_i0_j0_wf("x1_1_1", "x1_1_2", "x1_2_3",  "x1_2_4", 0),
                mk_pipi_i0_j0_wf("x1_1_1", "x1_1_2", "x1_2_3",  "x1_2_4", 1),
                mk_pipi_i0_j0_wf("x1_1_1", "x1_1_2", "x1_2_3",  "x1_2_4", 2),
                mk_pipi_i0_j0_wf("x1_1_1", "x1_1_2", "x1_2_3",  "x1_2_4", 3),
                mk_kkkk_i0_j0_wf("x1_1_1", "x1_1_2", "x1_2_3p", "x1_2_4p", 0),
                mk_kkkk_i0_j0_wf("x1_1_1", "x1_1_2", "x1_2_3p", "x1_2_4p", 1),
                ]
        snk_op_list = [
                mk_sigma_wf("x2_1_1", "x2_1_2", mom0) + f"sigma(x2_1_1,x2_1_2)",
                mk_pipi_i0_j0_wf("x2_1_1", "x2_1_2", "x2_2_3",  "x2_2_4", 0),
                mk_pipi_i0_j0_wf("x2_1_1", "x2_1_2", "x2_2_3",  "x2_2_4", 1),
                mk_pipi_i0_j0_wf("x2_1_1", "x2_1_2", "x2_2_3",  "x2_2_4", 2),
                mk_pipi_i0_j0_wf("x2_1_1", "x2_1_2", "x2_2_3",  "x2_2_4", 3),
                mk_kkkk_i0_j0_wf("x2_1_1", "x2_1_2", "x2_2_3p", "x2_2_4p", 0),
                mk_kkkk_i0_j0_wf("x2_1_1", "x2_1_2", "x2_2_3p", "x2_2_4p", 1),
                ]
        t_sel = mk_fac("time_selector_8(x1_1_1)") + f"time_selector_8(x1_1_1)"
        exprs = [
                mk_fac(1) + f"1",
                t_sel,
                ]
        for src_op in src_op_list:
            exprs.append((mk_fac(1) + f"src_op_only") * src_op)
        for snk_op in snk_op_list:
            exprs.append((mk_fac(1) + f"snk_op_only") * snk_op)
        for src_op in src_op_list:
            for snk_op in snk_op_list:
                exprs.append((snk_op * src_op, None, "TypeV", "TypeR", "TypeC", "TypeD",))
                exprs.append((t_sel * snk_op * src_op, "TypeR", "TypeC",))
        cexpr = contract_simplify_compile(
                *exprs,
                is_isospin_symmetric_limit=True,
                diagram_type_dict=diagram_type_dict,
                )
        return cexpr
    base_positions_dict = {
            "wave_function": wave_function,
            "momentum_factor": momentum_factor,
            "time_selector_8": time_selector_8,
            "r_pi": 1.5,
            "r_sigma": 1.5,
            "r_kk": 2.0,
            "Coordinate": q.Coordinate,
            }
    return cache_compiled_cexpr(calc_cexpr, fn_base, base_positions_dict=base_positions_dict, is_cython=is_cython)

@q.timer_verbose
def auto_contract_meson_meson_i0_j0_corr_wf(job_tag, traj, get_get_prop):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract/traj-{traj}/{fname}.lat"
    fn_src_op = f"{job_tag}/auto-contract/traj-{traj}/{fname}_src_op_only.lat"
    fn_snk_op = f"{job_tag}/auto-contract/traj-{traj}/{fname}_snk_op_only.lat"
    if all([ get_load_path(f) is not None for f in [ fn, fn_src_op, fn_snk_op, ] ]):
        return
    cexpr = get_cexpr_meson_meson_i0_j0_corr_wf()
    expr_names = np.array(get_expr_names(cexpr))
    idx_arr_of_src_op_only = np.array([ 0, ] + [ idx for idx, name in enumerate(expr_names) if " src_op_only " in name ], int)
    idx_arr_of_snk_op_only = np.array([ 0, ] + [ idx for idx, name in enumerate(expr_names) if " snk_op_only " in name ], int)
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    t_size = total_site[3]
    geo = q.Geometry(total_site)
    total_volume = geo.total_volume
    get_prop = get_get_prop()
    xg_list = get_all_points(total_site)
    xg_local_list = [ q.Coordinate(xg) for xg in geo.xg_arr() ]
    sample_num = get_param(job_tag, "measurement", fname,
                           "sample_num", default=512)
    sample_size = get_param(job_tag, "measurement", fname,
                            "sample_size", default=128)
    rs = q.RngState(f"{job_tag}-{traj}-{fname}")
    mpi_chunk = q.get_mpi_chunk(list(range(sample_num)))
    t_sep_range = get_param(job_tag, "measurement", fname,
                            "t_sep_range", default=17)
    t_sep_range = min(t_size, t_sep_range)
    def load_data():
        for idx in mpi_chunk:
            yield idx
    @q.timer
    def feval_single(idx, idx2, t_src, t_sep):
        rsi = rs.split(f"{idx}-{idx2}-{t_src}-{t_sep}")
        x1_1_1 = rsi.c_rand_gen(total_site)
        x1_1_2 = rsi.c_rand_gen(total_site)
        x1_2_3 = rsi.c_rand_gen(total_site)
        x1_2_4 = rsi.c_rand_gen(total_site)
        x1_2_3p = x1_2_3.copy()
        x1_2_4p = x1_2_4.copy()
        x2_1_1 = rsi.c_rand_gen(total_site)
        x2_1_2 = rsi.c_rand_gen(total_site)
        x2_2_3 = rsi.c_rand_gen(total_site)
        x2_2_4 = rsi.c_rand_gen(total_site)
        x2_2_3p = x2_2_3.copy()
        x2_2_4p = x2_2_4.copy()
        x1_1_1[3] = t_src
        x1_1_2[3] = x1_1_1[3]
        x1_2_3[3] = (x1_1_1[3] - 3) % t_size
        x1_2_3p[3] = (x1_1_1[3] - 2) % t_size
        x1_2_4[3] = x1_2_3[3]
        x1_2_4p[3] = x1_2_3p[3]
        x2_1_1[3] = (x1_1_1[3] + t_sep) % t_size
        x2_1_2[3] = x2_1_1[3]
        x2_2_3[3] = (x2_1_1[3] + 3) % t_size
        x2_2_3p[3] = (x2_1_1[3] + 2) % t_size
        x2_2_4[3] = x2_2_3[3]
        x2_2_4p[3] = x2_2_3p[3]
        pd = {
                "x1_1_1":  ("point", x1_1_1,),
                "x1_1_2":  ("point", x1_1_2,),
                "x1_2_3":  ("point", x1_2_3,),
                "x1_2_3p": ("point", x1_2_3p,),
                "x1_2_4":  ("point", x1_2_4,),
                "x1_2_4p": ("point", x1_2_4p,),
                "x2_1_1":  ("point", x2_1_1,),
                "x2_1_2":  ("point", x2_1_2,),
                "x2_2_3":  ("point", x2_2_3,),
                "x2_2_3p": ("point", x2_2_3p,),
                "x2_2_4":  ("point", x2_2_4,),
                "x2_2_4p": ("point", x2_2_4p,),
                "size": total_site,
                }
        val = eval_cexpr(cexpr, positions_dict=pd, get_prop=get_prop)
        return val
    @q.timer
    def feval(args):
        idx = args
        values = np.zeros((t_sep_range, len(expr_names),), dtype=complex)
        values_src_op = np.zeros((t_size, len(expr_names[idx_arr_of_src_op_only]),), dtype=complex)
        values_snk_op = np.zeros((t_size, len(expr_names[idx_arr_of_snk_op_only]),), dtype=complex)
        for idx2 in range(sample_size):
            for t_src in range(t_size):
                for t_sep in range(t_sep_range):
                    t_snk = (t_src + t_sep) % t_size
                    val = feval_single(idx, idx2, t_src, t_sep)
                    values[t_sep] += val
                    values_src_op[t_src] += val[idx_arr_of_src_op_only]
                    values_snk_op[t_snk] += val[idx_arr_of_snk_op_only]
        return values, values_src_op, values_snk_op,
    def sum_function(val_list):
        values = np.zeros((t_sep_range, len(expr_names),), dtype=complex)
        values_src_op = np.zeros((t_size, len(expr_names[idx_arr_of_src_op_only]),), dtype=complex)
        values_snk_op = np.zeros((t_size, len(expr_names[idx_arr_of_snk_op_only]),), dtype=complex)
        for idx, (val, val_src_op, val_snk_op,) in enumerate(val_list):
            values += val
            values_src_op += val_src_op
            values_snk_op += val_snk_op
            q.displayln_info(f"{fname}: {idx+1}/{len(mpi_chunk)}")
        return (
                q.glb_sum(values.transpose(1, 0,)),
                q.glb_sum(values_src_op.transpose(1, 0,)),
                q.glb_sum(values_snk_op.transpose(1, 0,)),
                )
    q.timer_fork(0)
    res_sum, res_sum_src_op, res_sum_snk_op = q.parallel_map_sum(feval, load_data(), sum_function=sum_function, chunksize=1)
    q.displayln_info(f"{fname}: timer_display for parallel_map_sum")
    q.timer_display()
    q.timer_merge()
    res_sum *= 1.0 / (sample_num * sample_size * t_size)
    res_sum_src_op *= 1.0 / (sample_num * sample_size * t_sep_range)
    res_sum_snk_op *= 1.0 / (sample_num * sample_size * t_sep_range)
    assert q.qnorm(res_sum[0] - 1.0) < 1e-10
    assert q.qnorm(res_sum_src_op[0] - 1.0) < 1e-10
    assert q.qnorm(res_sum_snk_op[0] - 1.0) < 1e-10
    ld = q.mk_lat_data([
        [ "expr_name", len(expr_names), list(expr_names), ],
        [ "t_sep", t_sep_range, [ str(t) for t in range(t_sep_range) ], ],
        ])
    ld.from_numpy(res_sum)
    ld.save(get_save_path(fn))
    json_results.append((f"{fname}: ld sig", q.get_data_sig(ld, q.RngState()),))
    for i, en in enumerate(expr_names):
        json_results.append((f"{fname}: ld '{en}' sig", q.get_data_sig(ld[i], q.RngState()),))
    ld_src_op = q.mk_lat_data([
        [ "expr_name", len(expr_names[idx_arr_of_src_op_only]), list(expr_names[idx_arr_of_src_op_only]), ],
        [ "t_op", t_size, [ str(t) for t in range(t_size) ], ],
        ])
    ld_src_op.from_numpy(res_sum_src_op)
    ld_src_op.save(get_save_path(fn_src_op))
    json_results.append((f"{fname}: ld_src_op sig", q.get_data_sig(ld_src_op, q.RngState()),))
    for i, en in enumerate(expr_names[idx_arr_of_src_op_only]):
        json_results.append((f"{fname}: ld '{en}' sig", q.get_data_sig(ld_src_op[i], q.RngState()),))
    ld_snk_op = q.mk_lat_data([
        [ "expr_name", len(expr_names[idx_arr_of_snk_op_only]), list(expr_names[idx_arr_of_snk_op_only]), ],
        [ "t_op", t_size, [ str(t) for t in range(t_size) ], ],
        ])
    ld_snk_op.from_numpy(res_sum_snk_op)
    ld_snk_op.save(get_save_path(fn_snk_op))
    json_results.append((f"{fname}: ld_snk_op sig", q.get_data_sig(ld_snk_op, q.RngState()),))
    for i, en in enumerate(expr_names[idx_arr_of_snk_op_only]):
        json_results.append((f"{fname}: ld '{en}' sig", q.get_data_sig(ld_snk_op[i], q.RngState()),))

# ----

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
                [
                    jj_d * pi0d,
                    (jj_d * pi0d, "TypeC"),
                    (jj_d * pi0d, "TypeD"),
                ]
                for pi0d in pi0d_list for jj_d in jj_d_list
                ]
        assert len(exprs_list_pi0_decay) == 2
        exprs_pi0_decay = [ e for el in exprs_list_pi0_decay for e in el ]
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

@q.timer_verbose
def auto_contract_pi0_gg(job_tag, traj, get_get_prop):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract/traj-{traj}/pi0-gg.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_pi0_gg()
    expr_names = get_expr_names(cexpr)
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    t_size = total_site[3]
    geo = q.Geometry(total_site)
    total_volume = geo.total_volume
    get_prop = get_get_prop()
    xg_list = get_all_points(total_site)
    xg_local_list = [ q.Coordinate(xg) for xg in geo.xg_arr() ]
    r_list = get_r_list(job_tag)
    r_sq_interp_idx_coef_list = get_r_sq_interp_idx_coef_list(job_tag)
    tsep = get_param(job_tag, "meson_tensor_tsep")
    def load_data():
        for xg_snk in xg_local_list:
            yield xg_snk
    @q.timer
    def feval(args):
        xg_snk = args
        res_list = []
        for xg_src in xg_list:
            x_rel = [ q.rel_mod(xg_snk[mu] - xg_src[mu], total_site[mu]) for mu in range(4) ]
            r_sq = q.get_r_sq(x_rel)
            r_idx_low, r_idx_high, coef_low, coef_high = r_sq_interp_idx_coef_list[r_sq]
            x_rel_t = x_rel[3]
            x_2_t = xg_src[3]
            x_1_t = x_2_t + x_rel_t
            t_2 = (max(x_1_t, x_2_t) + tsep) % total_site[3]
            t_1 = (min(x_1_t, x_2_t) - tsep) % total_site[3]
            pd = {
                    "x_2" : ("point", xg_snk,),
                    "x_1" : ("point", xg_src,),
                    "t_2" : ("wall", t_2),
                    "t_1" : ("wall", t_1),
                    "size" : total_site,
                    }
            t = x_rel_t % t_size
            val = eval_cexpr(cexpr, positions_dict=pd, get_prop=get_prop)
            res_list.append((val, t, r_idx_low, r_idx_high, coef_low, coef_high,))
        return res_list
    def sum_function(val_list):
        values = np.zeros((t_size, len(r_list), len(expr_names),), dtype=complex)
        for idx, res_list in enumerate(val_list):
            for val, t, r_idx_low, r_idx_high, coef_low, coef_high in res_list:
                values[t, r_idx_low] += coef_low * val
                values[t, r_idx_high] += coef_high * val
            q.displayln_info(f"{fname}: {idx+1}/{len(xg_local_list)}")
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
    json_results.append((f"{fname}: ld_sum sig", q.get_data_sig(ld_sum, q.RngState()),))
    for i, en in enumerate(expr_names):
        json_results.append((f"{fname}: ld_sum '{en}' sig", q.get_data_sig(ld_sum[i], q.RngState()),))

# ----

@q.timer_verbose
def run_job(job_tag, traj):
    #
    # fix gauge field in checking
    # all the props (including psrc props) are generated with Coulomb gauge fixing
    #
    traj_gf = 1000
    #
    fns_props = [
            (f"{job_tag}/prop-psrc-light/traj-{traj_gf}.qar", f"{job_tag}/prop-psrc-light/traj-{traj_gf}/geon-info.txt",),
            (f"{job_tag}/prop-psrc-strange/traj-{traj_gf}.qar", f"{job_tag}/prop-psrc-strange/traj-{traj_gf}/geon-info.txt",),
            (f"{job_tag}/prop-wsrc-light/traj-{traj_gf}.qar", f"{job_tag}/prop-wsrc-light/traj-{traj_gf}/geon-info.txt",),
            (f"{job_tag}/prop-wsrc-strange/traj-{traj_gf}.qar", f"{job_tag}/prop-wsrc-strange/traj-{traj_gf}/geon-info.txt",),
            (f"{job_tag}/psel-prop-wsrc-light/traj-{traj_gf}.qar", f"{job_tag}/psel-prop-wsrc-light/traj-{traj_gf}/checkpoint.txt",),
            (f"{job_tag}/psel-prop-wsrc-strange/traj-{traj_gf}.qar", f"{job_tag}/psel-prop-wsrc-strange/traj-{traj_gf}/checkpoint.txt",),
            ]
    #
    fns_produce = fns_props + [
            f"{job_tag}/auto-contract/traj-{traj}/checkpoint.txt",
            #
            (f"{job_tag}/configs/ckpoint_lat.{traj_gf}", f"{job_tag}/configs/ckpoint_lat.IEEE64BIG.{traj_gf}",),
            #
            f"{job_tag}/gauge-transform/traj-{traj_gf}.field",
            #
            f"{job_tag}/eig/traj-{traj_gf}",
            f"{job_tag}/eig-strange/traj-{traj_gf}",
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
    def run_with_eig():
        get_eig = run_eig(job_tag, traj_gf, get_gf)
        if get_eig is None:
            return
        # run_get_inverter_checker(job_tag, traj_gf, inv_type=0, get_gf=get_gf, get_gt=get_gt, get_eig=get_eig)
        run_prop_psrc_checker(job_tag, traj_gf, inv_type=0, get_gf=get_gf, get_eig=get_eig, get_gt=get_gt)
        run_prop_wsrc_checker(job_tag, traj_gf, inv_type=0, get_gf=get_gf, get_eig=get_eig, get_gt=get_gt)
        q.clean_cache(q.cache_inv)
    #
    def run_with_eig_strange():
        get_eig = run_eig_strange(job_tag, traj_gf, get_gf)
        if get_eig is None:
            return
        # run_get_inverter_checker(job_tag, traj_gf, inv_type=1, get_gf=get_gf, get_gt=get_gt, get_eig=get_eig)
        run_prop_psrc_checker(job_tag, traj_gf, inv_type=1, get_gf=get_gf, get_eig=get_eig, get_gt=get_gt)
        run_prop_wsrc_checker(job_tag, traj_gf, inv_type=1, get_gf=get_gf, get_eig=get_eig, get_gt=get_gt)
        q.clean_cache(q.cache_inv)
    #
    run_with_eig()
    run_with_eig_strange()
    #
    run_r_list(job_tag)
    get_get_prop = run_get_prop_checker(job_tag, traj_gf, get_gf=get_gf, get_gt=get_gt)
    #
    fn_checkpoint = f"{job_tag}/auto-contract/traj-{traj}/checkpoint.txt"
    if get_load_path(fn_checkpoint) is None and get_get_prop is not None:
        if q.obtain_lock(f"locks/{job_tag}-{traj}-auto-contract"):
            q.timer_fork()
            # ADJUST ME
            auto_contract_pi0_gg(job_tag, traj, get_get_prop)
            auto_contract_meson_meson_i0_j0_corr_wf(job_tag, traj, get_get_prop)
            auto_contract_meson_corr_wf(job_tag, traj, get_get_prop)
            auto_contract_meson_corr_psnk_psrc_rand(job_tag, traj, get_get_prop)
            #
            auto_contract_meson_corr(job_tag, traj_gf, get_get_prop)
            auto_contract_meson_corr_psnk(job_tag, traj_gf, get_get_prop)
            auto_contract_meson_corr_psnk_psrc(job_tag, traj_gf, get_get_prop)
            #
            q.qtouch_info(get_save_path(fn_checkpoint))
            q.release_lock()
            q.displayln_info("timer_display for runjob")
            q.timer_display()
            q.timer_merge()
            # q.clean_cache()

@q.timer_verbose
def get_all_cexpr():
    benchmark_eval_cexpr(get_cexpr_pi0_gg())
    benchmark_eval_cexpr(get_cexpr_meson_corr())
    benchmark_eval_cexpr(get_cexpr_meson_corr_wf())
    benchmark_eval_cexpr(get_cexpr_meson_meson_i0_j0_corr_wf())

# ----

set_param("test-4nt8", "trajs")(list(range(1000, 1001)))
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
set_param("test-4nt8", "cg_params-0-2", "maxiter")(5)
set_param("test-4nt8", "cg_params-0-2", "maxcycle")(1)
set_param("test-4nt8", "cg_params-1-2", "maxiter")(5)
set_param("test-4nt8", "cg_params-1-2", "maxcycle")(1)
set_param("test-4nt8", "fermion_params", 0, 2, "Ls")(8)
set_param("test-4nt8", "fermion_params", 1, 2, "Ls")(8)
set_param("test-4nt8", "fermion_params", 2, 2, "Ls")(8)

set_param("test-4nt16", "trajs")(list(range(1000, 1010)))
set_param("test-4nt16", "mk_sample_gauge_field", "rand_n_step")(2)
set_param("test-4nt16", "mk_sample_gauge_field", "flow_n_step")(8)
set_param("test-4nt16", "mk_sample_gauge_field", "hmc_n_traj")(1)
set_param("test-4nt16", "lanc_params", 0, 0, "cheby_params")({ "low": 0.45, "high": 5.5, "order": 20, })
set_param("test-4nt16", "lanc_params", 0, 0, "irl_params")({ "Nstop": 1000, "Nk": 1100, "Nm": 1300, "resid": 1e-8, "betastp": 0.0, "maxiter": 20, "Nminres": 0, })
set_param("test-4nt16", "clanc_params", 0, 0, "nbasis")(1000)
set_param("test-4nt16", "clanc_params", 0, 0, "block")([ 4, 4, 2, 2, ])
set_param("test-4nt16", "clanc_params", 0, 0, "cheby_params")({ "low": 0.45, "high": 5.5, "order": 20, })
set_param("test-4nt16", "clanc_params", 0, 0, "save_params")({ "nsingle": 100, "mpi": [ 1, 1, 1, 4, ], })
set_param("test-4nt16", "clanc_params", 0, 0, "irl_params")({ "Nstop": 1000, "Nk": 1100, "Nm": 1300, "resid": 1e-8, "betastp": 0.0, "maxiter": 20, "Nminres": 0, })
set_param("test-4nt16", "clanc_params", 1, 0)(get_param("test-4nt16", "clanc_params", 0, 0).copy())
set_param("test-4nt16", "lanc_params", 1, 0)(get_param("test-4nt16", "lanc_params", 0, 0).copy())
set_param("test-4nt16", "lanc_params", 1, 0, "fermion_params")(get_param("test-4nt16", "fermion_params", 1, 0).copy())

set_param("test-4nt64", "trajs")(list(range(1000, 10000)))
set_param("test-4nt64", "lanc_params", 0, 0, "cheby_params")({ "low": 0.22, "high": 5.5, "order": 30, })
set_param("test-4nt64", "lanc_params", 0, 0, "irl_params")({ "Nstop": 1000, "Nk": 1100, "Nm": 1300, "resid": 1e-8, "betastp": 0.0, "maxiter": 20, "Nminres": 0, })
set_param("test-4nt64", "clanc_params", 0, 0, "nbasis")(1000)
set_param("test-4nt64", "clanc_params", 0, 0, "block")([ 4, 4, 2, 8, ])
set_param("test-4nt64", "clanc_params", 0, 0, "cheby_params")({ "low": 0.22, "high": 5.5, "order": 30, })
set_param("test-4nt64", "clanc_params", 0, 0, "save_params")({ "nsingle": 100, "mpi": [ 1, 1, 1, 4, ], })
set_param("test-4nt64", "clanc_params", 0, 0, "irl_params")({ "Nstop": 1000, "Nk": 1100, "Nm": 1300, "resid": 1e-8, "betastp": 0.0, "maxiter": 20, "Nminres": 0, })
set_param("test-4nt64", "clanc_params", 1, 0)(get_param("test-4nt64", "clanc_params", 0, 0).copy())
set_param("test-4nt64", "lanc_params", 1, 0)(get_param("test-4nt64", "lanc_params", 0, 0).copy())
set_param("test-4nt64", "lanc_params", 1, 0, "fermion_params")(get_param("test-4nt64", "fermion_params", 1, 0).copy())

set_param("test-4nt8", "measurement", "auto_contract_meson_corr_wf", "sample_num")(32)
set_param("test-4nt8", "measurement", "auto_contract_meson_corr_wf", "sample_size")(2)
set_param("test-4nt8", "measurement", "auto_contract_meson_corr_wf", "t_sep_range")(4)
set_param("test-4nt8", "measurement", "auto_contract_meson_meson_i0_j0_corr_wf", "sample_num")(32)
set_param("test-4nt8", "measurement", "auto_contract_meson_meson_i0_j0_corr_wf", "sample_size")(2)
set_param("test-4nt8", "measurement", "auto_contract_meson_meson_i0_j0_corr_wf", "t_sep_range")(4)

set_param("test-4nt16", "measurement", "auto_contract_meson_corr_wf", "sample_num")(32)
set_param("test-4nt16", "measurement", "auto_contract_meson_corr_wf", "sample_size")(2)
set_param("test-4nt16", "measurement", "auto_contract_meson_corr_wf", "t_sep_range")(6)
set_param("test-4nt16", "measurement", "auto_contract_meson_meson_i0_j0_corr_wf", "sample_num")(32)
set_param("test-4nt16", "measurement", "auto_contract_meson_meson_i0_j0_corr_wf", "sample_size")(2)
set_param("test-4nt16", "measurement", "auto_contract_meson_meson_i0_j0_corr_wf", "t_sep_range")(6)

set_param("test-4nt64", "measurement", "auto_contract_meson_corr_wf", "sample_num")(512)
set_param("test-4nt64", "measurement", "auto_contract_meson_corr_wf", "sample_size")(128)
set_param("test-4nt64", "measurement", "auto_contract_meson_corr_wf", "t_sep_range")(17)
set_param("test-4nt64", "measurement", "auto_contract_meson_meson_i0_j0_corr_wf", "sample_num")(512)
set_param("test-4nt64", "measurement", "auto_contract_meson_meson_i0_j0_corr_wf", "sample_size")(128)
set_param("test-4nt64", "measurement", "auto_contract_meson_meson_i0_j0_corr_wf", "t_sep_range")(17)

tag = "meson_tensor_tsep"
set_param("test-4nt8", tag)(1)
set_param("test-4nt16", tag)(2)
set_param("test-4nt64", tag)(8)

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

    is_performing_inversion = q.get_arg("--no-inversion", default=None) is None

    is_performing_contraction = q.get_arg("--no-contract", default=None) is None

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
            q.check_time_limit()
            run_job(job_tag, traj)

    q.check_log_json(__file__, json_results, check_eps=5e-5)

    gracefully_finish()

# ----
