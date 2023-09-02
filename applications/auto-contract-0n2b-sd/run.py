#!/usr/bin/env python3

"""

@article{Detmold:2022jwu,
    author = "Detmold, William and Jay, William I. and Murphy, David J. and Oare, Patrick R. and Shanahan, Phiala E.",
    title = "{Neutrinoless Double Beta Decay from Lattice QCD: The Short-Distance $\pi^-\rightarrow\pi^+ e^- e^-$ Amplitude}",
    eprint = "2208.05322",
    archivePrefix = "arXiv",
    primaryClass = "hep-lat",
    reportNumber = "MIT-CTP/5414",
    month = "8",
    year = "2022"
}

Operators' definition follow Eq. (13) of the above reference.

"""

from auto_contractor.operators import *

import functools
import math
import os
import time
import importlib
import sys

from qlat_scripts.v1.jobs import *
from qlat_scripts.v1.load_data import *
from qlat_scripts.v1.params import *

# ----

load_path_list[:] = [
        "results",
        "qcddata",
        "qcddata-1",
        "qcddata-2",
        "qcddata-3",
        "qcddata-4",
        "qcddata-5",
        "../qcddata",
        "../qcddata-1",
        "../qcddata-2",
        "../qcddata-3",
        "../qcddata-4",
        "../qcddata-5",
        "../mk-gf-gt/results",
        "../mk-sel/results",
        "../mk-selected-data/results",
        "../mk-wsrc-prop/results",
        os.path.join(os.getenv("HOME"), "qcddata"),
        os.path.join(os.getenv("HOME"), "qcddata-1"),
        os.path.join(os.getenv("HOME"), "qcddata-2"),
        os.path.join(os.getenv("HOME"), "qcddata-3"),
        os.path.join(os.getenv("HOME"), "qcddata-4"),
        os.path.join(os.getenv("HOME"), "qcddata-5"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-gf-gt/results"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-sel/results"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-selected-data/results"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-wsrc-prop/results"),
        ]

auto_contractor_chunk_size = 128

# ----

def mk_jw_v_mu(p, mu):
    return mk_jpi_mu(p, mu) + mk_jk_mu(p, mu)

def mk_jw_a_mu(p, mu):
    return mk_j5pi_mu(p, mu) + mk_j5k_mu(p, mu)

def mk_sw5(p):
    return mk_pi_p(p, is_dagger = True) + mk_k_p(p, is_dagger = True)

def mk_ss(p, f1, f2):
    s = new_spin_index()
    c = new_color_index()
    v1 = Qb(f1, p, s, c) * Qv(f2, p, s, c)
    s = new_spin_index()
    c = new_color_index()
    v2 = Qb(f1, p, s, c) * Qv(f2, p, s, c)
    ret = v1 * v2
    return ret + f"SS({p},{f1},{f2})"

def mk_pp(p, f1, f2):
    s1 = new_spin_index()
    s2 = new_spin_index()
    c = new_color_index()
    v1 = Qb(f1, p, s1, c) * G(5, s1, s2) * Qv(f2, p, s2, c)
    s1 = new_spin_index()
    s2 = new_spin_index()
    c = new_color_index()
    v2 = Qb(f1, p, s1, c) * G(5, s1, s2) * Qv(f2, p, s2, c)
    ret = v1 * v2
    return ret + f"PP({p},{f1},{f2})"

def mk_vv(p, f1, f2):
    s1 = new_spin_index()
    s2 = new_spin_index()
    c = new_color_index()
    v1_list = []
    for mu in range(4):
        v1 = Qb(f1, p, s1, c) * G(mu, s1, s2) * Qv(f2, p, s2, c)
        v1_list.append(v1)
    s1 = new_spin_index()
    s2 = new_spin_index()
    c = new_color_index()
    v2_list = []
    for mu in range(4):
        v2 = Qb(f1, p, s1, c) * G(mu, s1, s2) * Qv(f2, p, s2, c)
        v2_list.append(v2)
    ret = 0
    for mu in range(4):
        ret = ret + v1_list[mu] * v2_list[mu]
    return ret + f"VV({p},{f1},{f2})"

def mk_aa(p, f1, f2):
    s1 = new_spin_index()
    s2 = new_spin_index()
    s3 = new_spin_index()
    c = new_color_index()
    v1_list = []
    for mu in range(4):
        v1 = Qb(f1, p, s1, c) * G(mu, s1, s2) * G(5, s2, s3) * Qv(f2, p, s3, c)
        v1_list.append(v1)
    s1 = new_spin_index()
    s2 = new_spin_index()
    s3 = new_spin_index()
    c = new_color_index()
    v2_list = []
    for mu in range(4):
        v2 = Qb(f1, p, s1, c) * G(mu, s1, s2) * G(5, s2, s3) * Qv(f2, p, s3, c)
        v2_list.append(v2)
    ret = 0
    for mu in range(4):
        ret = ret + v1_list[mu] * v2_list[mu]
    return ret + f"AA({p},{f1},{f2})"

def mk_tt(p, f1, f2):
    mu_nu_list = []
    for mu in range(4):
        for nu in range (4):
            if mu < nu:
                mu_nu_list.append((mu, nu,))
    s1 = new_spin_index()
    s2 = new_spin_index()
    s3 = new_spin_index()
    c = new_color_index()
    v1_list = []
    for mu, nu in mu_nu_list:
        v1 = Qb(f1, p, s1, c) * G(mu, s1, s2) * G(nu, s2, s3) * Qv(f2, p, s3, c)
        v1_list.append(v1)
    s1 = new_spin_index()
    s2 = new_spin_index()
    s3 = new_spin_index()
    c = new_color_index()
    v2_list = []
    for mu, nu in mu_nu_list:
        v2 = Qb(f1, p, s1, c) * G(mu, s1, s2) * G(nu, s2, s3) * Qv(f2, p, s3, c)
        v2_list.append(v2)
    ret = 0
    for idx in range(len(mu_nu_list)):
        ret = ret + v1_list[idx] * v2_list[idx]
    return ret + f"TT({p},{f1},{f2})"

def mk_0n2b_sd_op_list(p, f1, f2):
    ss = mk_ss(p, f1, f2)
    pp = mk_pp(p, f1, f2)
    vv = mk_vv(p, f1, f2)
    aa = mk_aa(p, f1, f2)
    tt = mk_tt(p, f1, f2)
    return [ ss, pp, vv, aa, tt, ]

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
    return cache_compiled_cexpr(calc_cexpr, fn_base)

@q.timer_verbose
def auto_contract_meson_corr(job_tag, traj, get_prop, get_psel, get_fsel):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_corr.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_corr()
    expr_names = get_expr_names(cexpr)
    total_site = rup.get_total_site(job_tag)
    t_size = total_site[3]
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
        return values.transpose(1, 0)
    q.timer_fork(0)
    res_sum = q.glb_sum(
            q.parallel_map_sum(feval, load_data(), sum_function = sum_function, chunksize = auto_contractor_chunk_size))
    q.displayln_info("timer_display for auto_contract_meson_corr")
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
    sig_msg_list = [
            f"CHECK: {fname}: ld sig: {q.get_double_sig(ld, q.RngState()):.14E}",
            ]
    for msg in sig_msg_list:
        q.displayln_info(msg)

@q.timer_verbose
def auto_contract_meson_corr_psnk(job_tag, traj, get_prop, get_psel, get_fsel):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_corr_psnk.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_corr()
    expr_names = get_expr_names(cexpr)
    total_site = rup.get_total_site(job_tag)
    t_size = total_site[3]
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
            f"CHECK: {fname}: ld sig: {q.get_double_sig(ld, q.RngState()):.14E}",
            ]
    for msg in sig_msg_list:
        q.displayln_info(msg)

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
    return cache_compiled_cexpr(calc_cexpr, fn_base)

@q.timer_verbose
def auto_contract_meson_jt(job_tag, traj, get_prop, get_psel, get_fsel):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_jt.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_jt()
    expr_names = get_expr_names(cexpr)
    total_site = rup.get_total_site(job_tag)
    t_size = total_site[3]
    psel = get_psel()
    fsel, fselc = get_fsel()
    xg_fsel_list = np.array(fsel.to_psel_local().to_list())
    xg_psel_list = np.array(psel.to_list())
    tsep = rup.dict_params[job_tag]["meson_tensor_tsep"]
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
            f"CHECK: {fname}: ld_sum sig: {q.get_double_sig(ld_sum, q.RngState()):.14E}",
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
    return cache_compiled_cexpr(calc_cexpr, fn_base)

@q.timer_verbose
def auto_contract_meson_m(job_tag, traj, get_prop, get_psel, get_fsel):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_m.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_m()
    expr_names = get_expr_names(cexpr)
    total_site = rup.get_total_site(job_tag)
    t_size = total_site[3]
    psel = get_psel()
    fsel, fselc = get_fsel()
    xg_fsel_list = np.array(fsel.to_psel_local().to_list())
    xg_psel_list = np.array(psel.to_list())
    tsep = rup.dict_params[job_tag]["meson_tensor_tsep"]
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
            f"CHECK: {fname}: ld_sum sig: {q.get_double_sig(ld_sum, q.RngState()):.14E}",
            ]
    for msg in sig_msg_list:
        q.displayln_info(msg)

### ------

@q.timer
def get_cexpr_meson_0n2b_sd():
    fn_base = "cache/auto_contract_cexpr/get_cexpr_meson_0n2b_sd"
    def calc_cexpr():
        diagram_type_dict = dict()
        diagram_type_dict[((('t_1', 'x'), 1), (('t_2', 'x'), 1), (('x', 't_1'), 1), (('x', 't_2'), 1))] = 'Type1'
        pi_mm = mk_pi_p("t_2", True) * mk_pi_m("t_1") + "pi+^dag(+tsep) * pi-(-tsep)"
        k_mm = mk_k_p("t_2", True) * mk_k_m("t_1") + "K+^dag(+tsep) * K-(-tsep)"
        pi_op_list = mk_0n2b_sd_op_list("x", "u", "d")
        k_op_list = mk_0n2b_sd_op_list("x", "u", "s")
        exprs = [ mk_expr(1) + f"1", ]
        exprs += [ op * pi_mm for op in pi_op_list ]
        exprs += [ op * k_mm for op in k_op_list ]
        cexpr = contract_simplify_compile(
                *exprs,
                is_isospin_symmetric_limit = True,
                diagram_type_dict = diagram_type_dict)
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, fn_base)

@q.timer_verbose
def auto_contract_meson_0n2b(job_tag, traj, get_prop, get_psel, get_fsel):
    fname = q.get_fname()
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_0n2b.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_0n2b_sd()
    expr_names = get_expr_names(cexpr)
    total_site = rup.get_total_site(job_tag)
    t_size = total_site[3]
    psel = get_psel()
    fsel, fselc = get_fsel()
    xg_fsel_list = np.array(fsel.to_psel_local().to_list())
    xg_fsel_list_chunk_list = q.get_chunk_list(xg_fsel_list, chunk_size = auto_contractor_chunk_size)
    xg_psel_list = np.array(psel.to_list())
    tsep_list = rup.dict_params[job_tag]["meson_tsep_list"]
    geo = q.Geometry(total_site, 1)
    total_volume = geo.total_volume()
    def load_data():
        for tsep1_idx, tsep1 in enumerate(tsep_list):
            for tsep2_idx, tsep2 in enumerate(tsep_list):
                for xg_snk_list in xg_fsel_list_chunk_list:
                    yield tsep1_idx, tsep2_idx, xg_snk_list
    @q.timer
    def feval(args):
        tsep1_idx, tsep2_idx, xg_snk_list = args
        tsep1 = tsep_list[tsep1_idx]
        tsep2 = tsep_list[tsep2_idx]
        val = 0
        for xg_snk in xg_snk_list:
            xg_snk = tuple(xg_snk.tolist())
            t = xg_snk[3]
            t_2 = (t + tsep2) % total_site[3]
            t_1 = (t - tsep1) % total_site[3]
            pd = {
                    "x" : ("point-snk", xg_snk,),
                    "t_1" : ("wall", t_1,),
                    "t_2" : ("wall", t_2,),
                    "size" : total_site,
                    }
            val = val + eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop)
        return tsep1_idx, tsep2_idx, val, len(xg_snk_list)
    def sum_function(val_list):
        values = np.zeros((len(tsep_list), len(tsep_list), len(expr_names),), dtype = complex)
        total_eval = 0
        total_eval_print = 0
        total_eval_print_step = len(xg_fsel_list)
        for tsep1_idx, tsep2_idx, val, n_eval in val_list:
            values[tsep1_idx, tsep2_idx] += val
            total_eval += n_eval
            if total_eval >= total_eval_print:
                q.displayln_info(f"{fname}: tsep_idx={tsep1_idx},{tsep2_idx}/{len(tsep_list)} total_eval={total_eval}/{len(xg_fsel_list) * len(tsep_list)**2}")
                total_eval_print += total_eval_print_step
        return values.transpose(2, 0, 1)
    q.timer_fork(0)
    res_sum = q.glb_sum(
            q.parallel_map_sum(feval, load_data(), sum_function = sum_function, chunksize = 1))
    q.displayln_info("timer_display for auto_contract_meson_m")
    q.timer_display()
    q.timer_merge()
    res_sum *= 1.0 / (total_volume * fsel.prob())
    assert q.qnorm(res_sum[0] - 1.0) < 1e-10
    ld_sum = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "tsep1", len(tsep_list), [ f"{t}" for t in tsep_list ], ],
        [ "tsep2", len(tsep_list), [ f"{t}" for t in tsep_list ], ],
        ])
    ld_sum.from_numpy(res_sum)
    ld_sum.save(get_save_path(fn))
    sig_msg_list = [
            f"CHECK: {fname}: ld_sum sig: {q.get_double_sig(ld_sum, q.RngState()):.14E}",
            ]
    for msg in sig_msg_list:
        q.displayln_info(msg)

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
            # f"{job_tag}/prop-psrc-light/traj-{traj}/geon-info.txt",
            # f"{job_tag}/prop-psrc-strange/traj-{traj}/geon-info.txt",
            f"{job_tag}/prop-wsrc-light/traj-{traj}/geon-info.txt",
            f"{job_tag}/prop-wsrc-strange/traj-{traj}/geon-info.txt",
            (f"{job_tag}/psel-prop-wsrc-light/traj-{traj}.qar", f"{job_tag}/psel-prop-wsrc-light/traj-{traj}/checkpoint.txt",),
            (f"{job_tag}/psel-prop-wsrc-strange/traj-{traj}.qar", f"{job_tag}/psel-prop-wsrc-strange/traj-{traj}/checkpoint.txt",),
            # (f"{job_tag}/psel-prop-psrc-light/traj-{traj}.qar", f"{job_tag}/psel-prop-psrc-light/traj-{traj}/checkpoint.txt",),
            # (f"{job_tag}/psel-prop-psrc-strange/traj-{traj}.qar", f"{job_tag}/psel-prop-psrc-strange/traj-{traj}/checkpoint.txt",),
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
    get_get_prop = run_get_prop(
            job_tag, traj,
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
                ],
            )
    #
    run_r_list(job_tag)
    #
    fn_checkpoint = f"{job_tag}/auto-contract/traj-{traj}/checkpoint.txt"
    if get_load_path(fn_checkpoint) is None:
        if q.obtain_lock(f"locks/{job_tag}-{traj}-auto-contract"):
            q.timer_fork()
            get_prop = get_get_prop()
            # ADJUST ME
            auto_contract_meson_0n2b(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contract_meson_corr(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contract_meson_corr_psnk(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contract_meson_jt(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contract_meson_m(job_tag, traj, get_prop, get_psel, get_fsel)
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
    benchmark_eval_cexpr(get_cexpr_meson_0n2b_sd())

def test():
    # ADJUST ME
    assert q.get_num_node() <= 4
    q.qremove_all_info("results/test-4nt8")
    q.qremove_info("results")
    assert not q.does_file_exist_sync_node("results")
    q.qremove_all_info("locks")
    # q.qremove_all_info("cache")
    # get_all_cexpr()
    run_job("test-4nt8", 1000)
    run_job("test-4nt16", 1000)
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

q.begin_with_mpi(size_node_list)

# ADJUST ME
# test()
get_all_cexpr()

# ADJUST ME
job_tags = [
        "test-4nt8", "test-4nt16",
        # "32IH1",
        # "32IH2",
        # "24IH1",
        # "24IH2",
        # "24IH3",
        # "64I",
        "24D",
        "48I",
        # "32D",
        # "32Dfine",
        # "24DH",
        # "16IH2",
        # "32IfineH",
        # "32IcoarseH1",
        ]

q.check_time_limit()

for job_tag in job_tags:
    if job_tag.startswith("test-"):
        if q.get_num_node() > 4:
            continue
    else:
        if q.get_num_node() <= 4:
            continue
    q.displayln_info(pprint.pformat(rup.dict_params[job_tag]))
    for traj in rup.dict_params[job_tag]["trajs"]:
        run_job(job_tag, traj)

q.clear_all_caches()
q.timer_display()

q.displayln_info("CHECK: finished successfully.")

q.end_with_mpi()
