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
                mk_pi_p("x_2")          * mk_j5pi_mu("x_1", 3, True) + f"pi+^dag(0) * j5pi_t^dag(-tsep)",
                mk_k_p("x_2")           * mk_j5k_mu("x_1", 3, True)  + f"K+^dag(0) * j5k_t^dag(-tsep)",
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
    total_site = get_param(job_tag, "total_site")
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
    total_site = get_param(job_tag, "total_site")
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
                "x_2" : ("point", xg_snk,),
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
    total_site = get_param(job_tag, "total_site")
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
    total_site = get_param(job_tag, "total_site")
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
                    "x_2": ("point", xg_snk,),
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
def get_all_points(total_site):
    n_points = total_site[0] * total_site[1] * total_site[2] * total_site[3]
    xg_list = []
    for index in range(n_points):
        xg = q.Coordinate()
        xg.from_index(index, total_site)
        xg_list.append(xg)
    return xg_list

@q.timer
def get_all_points_psel(total_site):
    geo = q.Geometry(total_site, 1)
    xg_list = get_all_points(total_site)
    psel = q.PointsSelection([ xg.to_list() for xg in xg_list ], geo)
    return psel

# ----

@q.timer
def run_get_inverter_checker(job_tag, traj, *, inv_type, get_gf, get_gt = None, get_eig = None):
    if None in [ get_gf, ]:
        return
    if get_gt is None:
        get_gt = lambda: None
    if get_eig is None:
        get_eig = lambda: None
    gf = get_gf()
    gt = get_gt()
    eig = get_eig()
    inv_acc = 2
    ru.get_inv(gf, job_tag, inv_type, inv_acc, gt = gt, eig = eig)

@q.timer_verbose
def compute_prop_1_checker(inv, src, *, tag, sfw, path_sp):
    fn_spw = os.path.join(path_sp, f"{tag} ; wsnk.lat")
    sol = inv * src
    sol.save_float_from_double(sfw, tag)
    sfw.flush()
    sol_ws = sol.glb_sum_tslice()
    sol_ws.save(get_save_path(fn_spw))
    return sol

@q.timer
def compute_prop_wsrc_checker(job_tag, tslice, inv_type, inv_acc, *,
                              idx, gf, gt, sfw, path_sp, eig, finished_tags):
    tag = f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc}"
    if tag in finished_tags:
        return None
    q.displayln_info(f"compute_prop_wsrc: idx={idx} tslice={tslice}", job_tag, inv_type, inv_acc)
    inv = ru.get_inv(gf, job_tag, inv_type, inv_acc, gt = gt, eig = eig)
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    geo = q.Geometry(total_site, 1)
    src = q.mk_wall_src(geo, tslice)
    prop = compute_prop_1_checker(inv, src, tag = tag, sfw = sfw, path_sp = path_sp)

@q.timer_verbose
def compute_prop_wsrc_all_checker(job_tag, traj, *,
                                  inv_type, gf, gt, eig):
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    inv_type_names = [ "light", "strange", ]
    inv_type_name = inv_type_names[inv_type]
    path_s = f"{job_tag}/prop-wsrc-{inv_type_name}/traj-{traj}"
    path_sp = f"{job_tag}/psel-prop-wsrc-{inv_type_name}/traj-{traj}"
    finished_tags = q.properly_truncate_fields(get_save_path(path_s + ".acc"))
    sfw = q.open_fields(get_save_path(path_s + ".acc"), "a", [ 1, 1, 1, 4, ])
    inv_acc = 2
    for idx, tslice in enumerate(range(total_site[3])):
        compute_prop_wsrc_checker(job_tag, tslice, inv_type, inv_acc = 2,
                                  idx = idx, gf = gf, gt = gt, sfw = sfw, path_sp = path_sp,
                                  eig = eig, finished_tags = finished_tags)
    sfw.close()
    q.qtouch_info(get_save_path(os.path.join(path_sp, "checkpoint.txt")))
    q.qrename_info(get_save_path(path_s + ".acc"), get_save_path(path_s))
    q.qar_create_info(get_save_path(path_sp + ".qar"), get_save_path(path_sp), is_remove_folder_after = True)
    q.qar_create_info(get_save_path(path_s + ".qar"), get_save_path(path_s), is_remove_folder_after = True)

@q.timer
def run_prop_wsrc_checker(job_tag, traj, *, inv_type, get_gf, get_eig, get_gt):
    if None in [ get_gf, get_gt, ]:
        return
    if get_eig is None:
        if inv_type == 0:
            return
        get_eig = lambda: None
    inv_type_names = [ "light", "strange", ]
    inv_type_name = inv_type_names[inv_type]
    if get_load_path(f"{job_tag}/prop-wsrc-{inv_type_name}/traj-{traj}/geon-info.txt") is not None:
        return
    if q.obtain_lock(f"locks/{job_tag}-{traj}-wsrc-{inv_type_name}"):
        gf = get_gf()
        gt = get_gt()
        eig = get_eig()
        compute_prop_wsrc_all_checker(job_tag, traj,
                                      inv_type = inv_type, gf = gf, gt = gt, eig = eig)
        q.release_lock()

# ----

@q.timer_verbose
def compute_prop_2_checker(inv, src, *, tag, sfw):
    sol = inv * src
    sol.save_float_from_double(sfw, tag)
    sfw.flush()
    return sol

@q.timer
def compute_prop_psrc_checker(job_tag, xg_src, inv_type, inv_acc, *,
                              idx, gf, gt, sfw, eig, finished_tags):
    xg = xg_src.to_list()
    xg_str = f"({xg[0]},{xg[1]},{xg[2]},{xg[3]})"
    tag = f"xg={xg_str} ; type={inv_type} ; accuracy={inv_acc}"
    if tag in finished_tags:
        return None
    q.displayln_info(f"compute_prop_psrc: {job_tag} idx={idx} tag='{tag}'")
    inv = ru.get_inv(gf, job_tag, inv_type, inv_acc, gt = gt, eig = eig)
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    geo = q.Geometry(total_site, 1)
    src = q.mk_point_src(geo, xg)
    prop = compute_prop_2_checker(inv, src, tag = tag, sfw = sfw)

@q.timer_verbose
def compute_prop_psrc_all_checker(job_tag, traj, *,
                                  inv_type, gf, gt, eig):
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    inv_type_names = [ "light", "strange", ]
    inv_type_name = inv_type_names[inv_type]
    path_s = f"{job_tag}/prop-psrc-{inv_type_name}/traj-{traj}"
    finished_tags = q.properly_truncate_fields(get_save_path(path_s + ".acc"))
    sfw = q.open_fields(get_save_path(path_s + ".acc"), "a", [ 1, 1, 1, 4, ])
    def comp(idx, xg_src, inv_acc):
        compute_prop_psrc_checker(job_tag, xg_src, inv_type, inv_acc,
                idx = idx, gf = gf, gt = gt, sfw = sfw,
                eig = eig, finished_tags = finished_tags)
    for idx, xg_src in enumerate(get_all_points(total_site)):
        comp(idx, xg_src, inv_acc = 2)
    sfw.close()
    q.qrename_info(get_save_path(path_s + ".acc"), get_save_path(path_s))
    q.qar_create_info(get_save_path(path_s + ".qar"), get_save_path(path_s), is_remove_folder_after = True)

@q.timer
def run_prop_psrc_checker(job_tag, traj, *, inv_type, get_gf, get_eig, get_gt):
    if None in [ get_gf, get_gt, ]:
        return
    if get_eig is None:
        if inv_type == 0:
            return
        get_eig = lambda: None
    inv_type_names = [ "light", "strange", ]
    inv_type_name = inv_type_names[inv_type]
    if get_load_path(f"{job_tag}/prop-psrc-{inv_type_name}/traj-{traj}/geon-info.txt") is not None:
        return
    if q.obtain_lock(f"locks/{job_tag}-{traj}-psrc-{inv_type_name}"):
        gf = get_gf()
        gt = get_gt()
        eig = get_eig()
        compute_prop_psrc_all_checker(job_tag, traj,
                                      inv_type = inv_type, gf = gf, gt = gt,
                                      eig = eig)
        q.release_lock()

# ----

@q.timer_verbose
def load_prop_psrc(job_tag, traj, inv_type):
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    inv_type_names = [ "light", "strange", ]
    inv_type_name = inv_type_names[inv_type]
    inv_tag_list = [ "l", "s", ]
    inv_tag = inv_tag_list[inv_type]
    inv_acc = 2
    path_s = f"{job_tag}/prop-psrc-{inv_type_name}/traj-{traj}/geon-info.txt"
    psel = get_all_points_psel(total_site)
    prop_list = []
    xg_list = [ q.Coordinate(xg) for xg in psel.to_list() ]
    sfr = q.open_fields(get_load_path(path_s), "r")
    for xg_src in xg_list:
        xg_idx = xg_src.to_index(total_site)
        xg = xg_src.to_list()
        xg_str = f"({xg[0]},{xg[1]},{xg[2]},{xg[3]})"
        tag = f"xg={xg_str} ; type={inv_type} ; accuracy={inv_acc}"
        prop = q.Prop()
        prop.load_double_from_float(sfr, tag)
        sp_prop = q.PselProp(psel)
        sp_prop @= prop
        prop_list.append(sp_prop)
    sfr.close()
    cache = q.mk_cache(f"prop_cache", f"{job_tag}", f"{traj}", inv_tag)
    cache["psnk-psrc"] = prop_list

@q.timer_verbose
def load_prop_wsrc(job_tag, traj, inv_type):
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    inv_type_names = [ "light", "strange", ]
    inv_type_name = inv_type_names[inv_type]
    inv_tag_list = [ "l", "s", ]
    inv_tag = inv_tag_list[inv_type]
    inv_acc = 2
    path_s = f"{job_tag}/prop-wsrc-{inv_type_name}/traj-{traj}/geon-info.txt"
    path_sp = f"{job_tag}/psel-prop-wsrc-{inv_type_name}/traj-{traj}/"
    psel = get_all_points_psel(total_site)
    psel_ts = q.get_psel_tslice(total_site)
    prop_list = []
    prop2_list = []
    tslice_list = list(range(total_site[3]))
    sfr = q.open_fields(get_load_path(path_s), "r")
    for tslice in tslice_list:
        tag = f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc}"
        prop = q.Prop()
        prop.load_double_from_float(sfr, tag)
        sp_prop = q.PselProp(psel)
        sp_prop @= prop
        prop_list.append(sp_prop)
        fn_spw = os.path.join(path_sp, f"{tag} ; wsnk.lat")
        sp_prop2 = q.PselProp(psel_ts)
        sp_prop2.load(get_load_path(fn_sp))
        prop2_list.append(sp_prop2)
    sfr.close()
    cache = q.mk_cache(f"prop_cache", f"{job_tag}", f"{traj}", inv_tag)
    cache["psnk-wsrc"] = prop_list
    cache["wsnk-wsrc"] = prop2_list

@q.timer_verbose
def run_get_prop_checker(job_tag, traj, *,
                         get_gf,
                         get_gt):
    traj_gf = traj
    fns_props = [
            (f"{job_tag}/prop-psrc-light/traj-{traj_gf}.qar", f"{job_tag}/prop-psrc-light/traj-{traj_gf}/geon-info.txt",),
            (f"{job_tag}/prop-psrc-strange/traj-{traj_gf}.qar", f"{job_tag}/prop-psrc-strange/traj-{traj_gf}/geon-info.txt",),
            (f"{job_tag}/prop-wsrc-light/traj-{traj_gf}.qar", f"{job_tag}/prop-wsrc-light/traj-{traj_gf}/geon-info.txt",),
            (f"{job_tag}/prop-wsrc-strange/traj-{traj_gf}.qar", f"{job_tag}/prop-wsrc-strange/traj-{traj_gf}/geon-info.txt",),
            (f"{job_tag}/psel-prop-wsrc-light/traj-{traj_gf}.qar", f"{job_tag}/psel-prop-wsrc-light/traj-{traj_gf}/checkpoint.txt",),
            (f"{job_tag}/psel-prop-wsrc-strange/traj-{traj_gf}.qar", f"{job_tag}/psel-prop-wsrc-strange/traj-{traj_gf}/checkpoint.txt",),
            ]
    for fn in fns_props:
        if get_load_path(fn) is None:
            return None
    @q.timer_verbose
    def mk_get_prop():
        q.timer_fork()
        total_site = q.Coordinate(get_param(job_tag, "total_site"))
        gf = get_gf()
        gt = get_gt()
        #
        load_prop_psrc(job_tag, traj, inv_type = 0)
        load_prop_psrc(job_tag, traj, inv_type = 1)
        load_prop_wsrc(job_tag, traj, inv_type = 0)
        load_prop_wsrc(job_tag, traj, inv_type = 1)
        #
        prop_cache = q.mk_cache(f"prop_cache", f"{job_tag}", f"{traj}")
        def get_prop(flavor, p_snk, p_src):
            cache = prop_cache[flavor]
            p_snk_tag, p_snk_xg = p_snk
            p_src_tag, p_src_xg = p_src
            if p_snk_tag == "point" and p_src_tag == "point":
                prop_list = cache["psnk-psrc"]
                p_src_idx = p_src_xg.to_index(total_site)
                p_snk_idx = p_snk_xg.to_index(total_site)
                return prop_list[p_src_idx].get_elem_wm(p_snk_idx)
            elif p_snk_tag == "point" and p_src_tag == "wall":
                prop_list = cache["psnk-wsrc"]
                p_snk_idx = p_snk_xg.to_index(total_site)
                return prop_list[p_src_xg].get_elem_wm(p_snk_idx)
            elif p_snk_tag == "wall" and p_src_tag == "point":
                prop_list = cache["psnk-wsrc"]
                p_src_idx = p_src_xg.to_index(total_site)
                return wilson_matrix_g5_herm(prop_list[p_snk_xg].get_elem_wm(p_src_idx))
            elif p_snk_tag == "wall" and p_src_tag == "wall":
                prop_list = cache["wsnk-wsrc"]
                return prop_list[p_src_xg].get_elem_wm(p_snk_xg)
            else:
                raise Exception(f"get_prop: f={flavor} snk={p_snk} src={p_src}")
        #
        q.timer_display()
        q.timer_merge()
        return get_prop
    return q.lazy_call(mk_get_prop)

# ----

@q.timer_verbose
def run_job(job_tag, traj):
    #
    traj_gf = 1000 # fix gauge field in checking
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
        # run_get_inverter_checker(job_tag, traj_gf, inv_type = 0, get_gf = get_gf, get_gt = get_gt, get_eig = get_eig)
        run_prop_psrc_checker(job_tag, traj_gf, inv_type = 0, get_gf = get_gf, get_eig = get_eig, get_gt = get_gt)
        run_prop_wsrc_checker(job_tag, traj_gf, inv_type = 0, get_gf = get_gf, get_eig = get_eig, get_gt = get_gt)
        q.clean_cache(q.cache_inv)
    #
    def run_with_eig_strange():
        get_eig = run_eig_strange(job_tag, traj_gf, get_gf)
        if get_eig is None:
            return
        # run_get_inverter_checker(job_tag, traj_gf, inv_type = 1, get_gf = get_gf, get_gt = get_gt, get_eig = get_eig)
        run_prop_psrc_checker(job_tag, traj_gf, inv_type = 1, get_gf = get_gf, get_eig = get_eig, get_gt = get_gt)
        run_prop_wsrc_checker(job_tag, traj_gf, inv_type = 1, get_gf = get_gf, get_eig = get_eig, get_gt = get_gt)
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
            get_prop = get_get_prop()
            # auto_contract_meson_corr(job_tag, traj, get_get_prop, get_psel, get_fsel)
            # auto_contract_meson_corr_psnk(job_tag, traj, get_get_prop, get_psel, get_fsel)
            # auto_contract_meson_corr_psrc(job_tag, traj, get_get_prop, get_psel, get_fsel)
            # auto_contract_meson_corr_psnk_psrc(job_tag, traj, get_get_prop, get_psel, get_fsel)
            #
            # q.qtouch_info(get_save_path(fn_checkpoint))
            q.release_lock()
            q.displayln_info("timer_display for runjob")
            q.timer_display()
            q.timer_merge()
            #
            q.clean_cache()

def get_all_cexpr():
    benchmark_eval_cexpr(get_cexpr_meson_corr())

set_param("test-4nt8", "trajs", value=[ 1000, 1010, ])
set_param("test-4nt8", "mk_sample_gauge_field", "rand_n_step", value=2)
set_param("test-4nt8", "mk_sample_gauge_field", "flow_n_step", value=8)
set_param("test-4nt8", "mk_sample_gauge_field", "hmc_n_traj", value=1)
set_param("test-4nt8", "fermion_params", 0, 2, "Ls", value=8)
set_param("test-4nt8", "fermion_params", 1, 2, "Ls", value=8)
set_param("test-4nt8", "fermion_params", 2, 2, "Ls", value=8)
set_param("test-4nt8", "lanc_params", 0, 0, "cheby_params", value={ "low": 0.8, "high": 5.5, "order": 20, })
set_param("test-4nt8", "lanc_params", 0, 0, "irl_params", value={ "Nstop": 1000, "Nk": 1100, "Nm": 1300, "resid": 1e-8, "betastp": 0.0, "maxiter": 20, "Nminres": 0, })
set_param("test-4nt8", "clanc_params", 0, 0, "nbasis", value=1000)
set_param("test-4nt8", "clanc_params", 0, 0, "block", value=[ 4, 4, 2, 2, ])
set_param("test-4nt8", "clanc_params", 0, 0, "cheby_params", value={ "low": 0.8, "high": 5.5, "order": 20, })
set_param("test-4nt8", "clanc_params", 0, 0, "save_params", value={ "nsingle": 100, "mpi": [ 1, 1, 1, 4, ], })
set_param("test-4nt8", "clanc_params", 0, 0, "irl_params", value={ "Nstop": 1000, "Nk": 1100, "Nm": 1300, "resid": 1e-8, "betastp": 0.0, "maxiter": 20, "Nminres": 0, })
set_param("test-4nt8", "clanc_params", 1, 0, value=get_param("test-4nt8", "clanc_params", 0, 0).copy())
set_param("test-4nt8", "lanc_params", 1, 0, value=get_param("test-4nt8", "lanc_params", 0, 0).copy())
set_param("test-4nt8", "lanc_params", 1, 0, "fermion_params", value=get_param("test-4nt8", "fermion_params", 1, 0).copy())

qg.begin_with_gpt()

job_tags = [
        "test-4nt8",
        # "test-4nt64",
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
