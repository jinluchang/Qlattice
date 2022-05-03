#!/usr/bin/env python3

from auto_contractor.operators import *
import rbc_ukqcd as ru
import qlat_gpt as qg

import functools
import math
import os

from jobs import *
from load_data import *
from params import *

load_path_list[:] = [
        "results",
        "../qcddata",
        os.path.join(os.getenv("HOME"), "qcddata"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-gf-gt/results"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-sel/results"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-fsel-self-loop/results"),
        # os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-selected-data/results"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-wsrc-prop/results"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-psrc-prop/results"),
        # os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-smear-prop/results"),
        ]

def rel_mod(x, size):
    x = (x + 2 * size) % size
    assert x >= 0
    if 2 * x >= size:
        return x - size
    else:
        return x

@q.timer
def get_cexpr_vev():
    def calc_cexpr():
        s = new_spin_index()
        c = new_color_index()
        p = "x"
        exprs = [
                Qb("u", "x", s, c) * Qv("u", "x", s, c) + "u_bar*u",
                Qb("s", "x", s, c) * Qv("s", "x", s, c) + "s_bar*s",
                Qb("c", "x", s, c) * Qv("c", "x", s, c) + "c_bar*c",
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.collect_op()
        return cexpr
    cexpr = q.pickle_cache_call(calc_cexpr, f"cache/auto_contractor_cexpr/vev-cexpr.pickle")
    q.displayln_info(display_cexpr_raw(cexpr))
    return cexpr

@q.timer_verbose
def auto_contractor_vev(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"auto-contract/{job_tag}/traj={traj}/vev.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_vev()
    expr_names = get_cexpr_names(cexpr)
    fsel, fselc = get_fsel()
    xg_fsel_list = fsel.to_psel_local().to_list()
    def feval(x):
        pd = {
                "x" : ("point-snk", x,),
                }
        res = eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop, is_only_total = "total")
        return res
    res_list = q.parallel_map(q.get_n_processes(), feval, xg_fsel_list)
    res_sum = q.glb_sum(sum(res_list))
    res_count = q.glb_sum(len(res_list))
    res_avg = res_sum / res_count
    ld = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        ])
    ld.from_numpy(res_avg)
    q.displayln_info(ld.show())
    ld.save(get_save_path(fn))

@q.timer
def get_cexpr_meson_f_corr():
    def calc_cexpr():
        exprs = [
                mk_pi_p("x2", True) * mk_pi_p("x1") + "(pi     * pi)",
                -1j * mk_j5pi_mu("x2", 3) * mk_pi_p("x1") + "(a_pi/i * pi)",
                mk_k_p("x2", True)  * mk_k_p("x1")  + "(k      * k )",
                -1j * mk_j5k_mu("x2", 3)  * mk_k_p("x1")  + "(a_k/i  * k )",
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.collect_op()
        return cexpr
    cexpr = q.pickle_cache_call(calc_cexpr, f"cache/auto_contractor_cexpr/meson_f_corr-cexpr.pickle")
    q.displayln_info(display_cexpr_raw(cexpr))
    return cexpr

@q.timer_verbose
def auto_contractor_meson_f_corr(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"auto-contract/{job_tag}/traj={traj}/meson_f_corr.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_f_corr()
    expr_names = get_cexpr_names(cexpr)
    total_site = ru.get_total_site(job_tag)
    fsel, fselc = get_fsel()
    xg_fsel_list = fsel.to_psel_local().to_list()
    def feval(x):
        l = []
        for t in range(total_site[3]):
            pd = {
                    "x2" : ("point-snk", x,),
                    "x1" : ("wall", (x[3] - t) % total_site[3],),
                    }
            res = eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop, is_only_total = "total")
            l.append(res)
        return np.array(l).transpose()
    res_list = q.parallel_map(q.get_n_processes(), feval, xg_fsel_list)
    res_sum = q.glb_sum(sum(res_list))
    res_count = q.glb_sum(len(res_list))
    res_avg = res_sum / res_count
    ld = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t_sep", total_site[3], ],
        ])
    ld.from_numpy(res_avg)
    q.displayln_info(ld.show())
    ld.save(get_save_path(fn))

@q.timer
def get_cexpr_hvp():
    def calc_cexpr():
        exprs = [
                mk_expr(1) + "1",
                -mk_jl_mu("x2", 0) * mk_jl_mu("x1", 0) + "(-jl_0 * jl_0)",
                -mk_jl_mu("x2", 1) * mk_jl_mu("x1", 1) + "(-jl_1 * jl_1)",
                -mk_jl_mu("x2", 2) * mk_jl_mu("x1", 2) + "(-jl_2 * jl_2)",
                -mk_jl_mu("x2", 3) * mk_jl_mu("x1", 3) + "(-jl_3 * jl_3)",
                -mk_jl_mu("x2", 5) * mk_jl_mu("x1", 5) + "(-jl_5 * jl_5)",
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.collect_op()
        return cexpr
    cexpr = q.pickle_cache_call(calc_cexpr, f"cache/auto_contractor_cexpr/hvp-cexpr.pickle")
    q.displayln_info(display_cexpr_raw(cexpr))
    return cexpr

@q.timer_verbose
def auto_contractor_hvp(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"auto-contract/{job_tag}/traj={traj}/hvp.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_hvp()
    expr_names = get_cexpr_names(cexpr)
    total_site = ru.get_total_site(job_tag)
    psel = get_psel()
    fsel, fselc = get_fsel()
    xg_fsel_list = fsel.to_psel_local().to_list()
    xg_psel_list = psel.to_list()
    vol = total_site[0] * total_site[1] * total_site[2]
    def feval(xg_src):
        counts = np.zeros(total_site[3], dtype = int)
        values = np.zeros((total_site[3], len(expr_names)), dtype = complex)
        for xg_snk in xg_fsel_list:
            pd = {
                    "x2" : ("point-snk", xg_snk,),
                    "x1" : ("point", xg_src,),
                    }
            res = eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop, is_only_total = "total")
            tsep = (xg_snk[3] - xg_src[3]) % total_site[3]
            counts[tsep] += 1
            values[tsep] += res
        counts = counts # counts[tsep]
        values = values.transpose() # values[expr_idx, tsep]
        return counts, values
    counts_list, values_list = zip(*q.parallel_map(q.get_n_processes(), feval, xg_psel_list))
    res_count = q.glb_sum(sum(counts_list))
    res_sum = q.glb_sum(sum(values_list))
    res_avg = res_sum * (vol / res_count)
    ld = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t_sep", total_site[3], ],
        ])
    ld.from_numpy(res_avg)
    q.displayln_info(ld.show())
    ld.save(get_save_path(fn))

@q.timer_verbose
def auto_contractor_hvp_field(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"auto-contract/{job_tag}/traj={traj}/hvp.field"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_hvp()
    expr_names = get_cexpr_names(cexpr)
    total_site = ru.get_total_site(job_tag)
    psel = get_psel()
    fsel, fselc = get_fsel()
    xg_fsel_list = fsel.to_psel_local().to_list()
    xg_psel_list = psel.to_list()
    geo = q.Geometry(total_site, 1)
    field = q.Field("Complex", geo, len(expr_names))
    field.set_zero()
    for idx, xg_src in enumerate(xg_psel_list):
        def feval(xg_snk):
            pd = {
                    "x2" : ("point-snk", xg_snk,),
                    "x1" : ("point", xg_src,),
                    }
            res = eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop, is_only_total = "total")
            return res
        values_list = q.parallel_map(q.get_n_processes(), feval, xg_fsel_list)
        assert len(values_list) == fsel.n_elems()
        values = q.SelectedField("Complex", fsel, len(expr_names))
        for idx, v in enumerate(values_list):
            values[idx] = v.tobytes()
        shift = [ -x for x in xg_src ]
        values_shifted = values.field_shift(shift)
        field += values_shifted
    # scale the value appropriately
    field *= 1.0 / (len(xg_psel_list) * fsel.prob())
    field.save_float_from_double(get_save_path(fn))
    q.displayln_info(field.glb_sum_tslice().to_numpy())

@q.timer
def get_cexpr_meson_v_v_meson():
    def calc_cexpr():
        x_1 = "x_1"
        x_2 = "x_2"
        xj_1 = "xj_1"
        xj_2 = "xj_2"
        exprs = [ mk_expr(1) + "1", ]
        exprs += [
                tr(gamma_5*S_l(x_1,x_2)*gamma_5*S_l(x_2,x_1))
                + f"pion-corr"
                ]
        exprs += [
                tr(gamma_va(mu)*S_l(xj_1,xj_2)*gamma_va(nu)*S_l(xj_2,xj_1))
                + f"hvp-{mu}-{nu}"
                for mu in range(4) for nu in range(4) ]
        exprs += [
                tr(gamma_va(mu)*S_l(xj_1,x_1)*gamma_5*S_l(x_1,xj_1))*tr(gamma_va(nu)*S_l(xj_2,x_2)*gamma_5*S_l(x_2,xj_2))
                + f"type1-{mu}-{nu}"
                for mu in range(4) for nu in range(4) ]
        exprs += [
                tr(gamma_va(mu)*S_l(xj_1,x_2)*gamma_5*S_l(x_2,xj_2)*gamma_va(nu)*S_l(xj_2,x_1)*gamma_5*S_l(x_1,xj_1))
                + f"type2-{mu}-{nu}"
                for mu in range(4) for nu in range(4) ]
        exprs += [
                tr(gamma_va(mu)*S_l(xj_1,xj_2)*gamma_va(nu)*S_l(xj_2,x_1)*gamma_5*S_l(x_1,x_2)*gamma_5*S_l(x_2,xj_1))
                + f"type3-{mu}-{nu}"
                for mu in range(4) for nu in range(4) ]
        exprs += [
                tr(gamma_va(mu)*S_l(xj_1,xj_2)*gamma_va(nu)*S_l(xj_2,xj_1))*tr(gamma_5*S_l(x_1,x_2)*gamma_5*S_l(x_2,x_1))
                + f"type4-{mu}-{nu}"
                for mu in range(4) for nu in range(4) ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.collect_op()
        return cexpr
    cexpr = q.pickle_cache_call(calc_cexpr, f"cache/auto_contractor_cexpr/meson_v_v_meson-cexpr.pickle")
    q.displayln_info(display_cexpr_raw(cexpr))
    return cexpr

@q.timer_verbose
def auto_contractor_meson_v_v_meson_field(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"auto-contract/{job_tag}/traj={traj}/meson_v_v_meson.field"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_v_v_meson()
    expr_names = get_cexpr_names(cexpr)
    total_site = ru.get_total_site(job_tag)
    psel = get_psel()
    fsel, fselc = get_fsel()
    xg_fsel_list = fsel.to_psel_local().to_list()
    xg_psel_list = psel.to_list()
    tsep = rup.dict_params[job_tag]["meson_tensor_tsep"]
    geo = q.Geometry(total_site, 1)
    field = q.Field("Complex", geo, len(expr_names))
    field.set_zero()
    for idx, xg_src in enumerate(xg_psel_list):
        def feval(xg_snk):
            xj_1 = ("point-snk", xg_snk,)
            xj_2 = ("point", xg_src,)
            tj_1 = xj_1[1][3]
            tj_2 = xj_2[1][3]
            tj_1 = tj_2 + rel_mod(tj_1 - tj_2, total_site[3])
            t_1 = (max(tj_1, tj_2) + tsep) % total_site[3]
            t_2 = (min(tj_1, tj_2) - tsep) % total_site[3]
            pd = {
                    "xj_1" : xj_1,
                    "xj_2" : xj_2,
                    "x_1" : ("wall", t_1),
                    "x_2" : ("wall", t_2),
                    }
            res = eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop, is_only_total = "total")
            return res
        values_list = q.parallel_map(q.get_n_processes(), feval, xg_fsel_list)
        assert len(values_list) == fsel.n_elems()
        values = q.SelectedField("Complex", fsel, len(expr_names))
        for idx, v in enumerate(values_list):
            values[idx] = v.tobytes()
        shift = [ -x for x in xg_src ]
        values_shifted = values.field_shift(shift)
        field += values_shifted
    field_r = q.Field("Complex", geo, len(expr_names))
    field_r.set_zero()
    for idx, xg_src in enumerate(xg_psel_list):
        def feval(xg_snk):
            xj_1 = ("point", xg_src,)
            xj_2 = ("point-snk", xg_snk,)
            tj_1 = xj_1[1][3]
            tj_2 = xj_2[1][3]
            tj_1 = tj_2 + rel_mod(tj_1 - tj_2, total_site[3])
            t_1 = (max(tj_1, tj_2) + tsep) % total_site[3]
            t_2 = (min(tj_1, tj_2) - tsep) % total_site[3]
            pd = {
                    "xj_1" : xj_1,
                    "xj_2" : xj_2,
                    "x_1" : ("wall", t_1),
                    "x_2" : ("wall", t_2),
                    }
            res = eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop, is_only_total = "total")
            return res
        values_list = q.parallel_map(q.get_n_processes(), feval, xg_fsel_list)
        assert len(values_list) == fsel.n_elems()
        values = q.SelectedField("Complex", fsel, len(expr_names))
        for idx, v in enumerate(values_list):
            values[idx] = v.tobytes()
        shift = [ -x for x in xg_src ]
        values_shifted = values.field_shift(shift)
        field_r += values_shifted
    field_r.reflect()
    field += field_r
    field *= 0.5
    # scale the value appropriately
    field *= 1 / (len(xg_psel_list) * fsel.prob())
    field.save_float_from_double(get_save_path(fn))
    ld = q.mk_lat_data([
        [ "t_sep", total_site[3], ],
        [ "expr_name", len(expr_names), expr_names, ],
        ])
    ld.from_list(field.glb_sum_tslice().to_lat_data().to_list())
    q.displayln_info(ld.show())

# ----

@q.timer
def get_cexpr_meson_corr():
    def calc_cexpr():
        exprs = [
                mk_pi_p("x_2", True) * mk_pi_p("x_1") + "pi * pi",
                mk_k_p("x_2", True)  * mk_k_p("x_1")  + "k  * k ",
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.collect_op()
        return cexpr
    cexpr = q.pickle_cache_call(calc_cexpr, f"cache/auto_contractor_cexpr/meson_corr-cexpr.pickle")
    q.displayln_info(display_cexpr_raw(cexpr))
    return cexpr

@q.timer_verbose
def auto_contractor_meson_corr(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"auto-contract/{job_tag}/traj={traj}/meson_corr.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_corr()
    expr_names = get_cexpr_names(cexpr)
    total_site = ru.get_total_site(job_tag)
    fsel, fselc = get_fsel()
    xg_fsel_list = fsel.to_psel_local().to_list()
    geo = q.Geometry(total_site, 1)
    total_volume = geo.total_volume()
    def fempty():
        values = np.zeros((len(expr_names), total_site[3],), dtype = complex)
        return 0, values
    @q.timer
    def feval(t_snk):
        counts, values = fempty()
        counts += 1
        values = values.transpose()
        for t_src in range(total_site[3]):
            t = (t_snk - t_src) % total_site[3]
            pd = {
                    "x_2" : ("wall", t_snk,),
                    "x_1" : ("wall", t_src,),
                    }
            values[t] += eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop, is_only_total = "total")
        values = values.transpose() # res[expr_name, t_sep]
        return counts, values
    t_snk_list = get_mpi_chunk(list(range(total_site[3])))
    counts_list, values_list = zip(fempty(), *q.parallel_map(q.get_n_processes(), feval, t_snk_list))
    res_count = q.glb_sum(sum(counts_list))
    res_sum = q.glb_sum(sum(values_list))
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
def auto_contractor_meson_corr_psnk(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"auto-contract/{job_tag}/traj={traj}/meson_corr_psnk.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_corr()
    expr_names = get_cexpr_names(cexpr)
    total_site = ru.get_total_site(job_tag)
    fsel, fselc = get_fsel()
    xg_fsel_list = fsel.to_psel_local().to_list()
    geo = q.Geometry(total_site, 1)
    total_volume = geo.total_volume()
    @q.timer
    def feval(xg_snk):
        res = np.zeros((total_site[3], len(expr_names),), dtype = complex)
        for t_src in range(total_site[3]):
            t = (xg_snk[3] - t_src) % total_site[3]
            pd = {
                    "x_2" : ("point-snk", xg_snk,),
                    "x_1" : ("wall", t_src,),
                    }
            res[t] += eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop, is_only_total = "total")
        res = res.transpose() # res[expr_name, t_sep]
        return 1.0, res
    counts_list, values_list = zip(*q.parallel_map(q.get_n_processes(), feval, xg_fsel_list))
    res_count = q.glb_sum(sum(counts_list))
    res_sum = q.glb_sum(sum(values_list))
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
def auto_contractor_meson_corr_psrc(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"auto-contract/{job_tag}/traj={traj}/meson_corr_psrc.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_corr()
    expr_names = get_cexpr_names(cexpr)
    total_site = ru.get_total_site(job_tag)
    psel = get_psel()
    fsel, fselc = get_fsel()
    xg_fsel_list = fsel.to_psel_local().to_list()
    xg_psel_list = psel.to_list()
    geo = q.Geometry(total_site, 1)
    total_volume = geo.total_volume()
    def fempty():
        values = np.zeros((len(expr_names), total_site[3],), dtype = complex)
        return 0, values
    @q.timer
    def feval(xg_src):
        counts, values = fempty()
        counts += 1
        values = values.transpose()
        for t_snk in range(total_site[3]):
            t = (t_snk - xg_src[3]) % total_site[3]
            pd = {
                    "x_2" : ("wall", t_snk,),
                    "x_1" : ("point", xg_src,),
                    }
            values[t] += eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop, is_only_total = "total")
        values = values.transpose() # values[expr_name, t_sep]
        return counts, values
    xg_src_list = get_mpi_chunk(xg_psel_list, rng_state = q.RngState("get_mpi_chunk"))
    counts_list, values_list = zip(fempty(), *q.parallel_map(q.get_n_processes(), feval, xg_src_list))
    res_count = q.glb_sum(sum(counts_list))
    res_sum = q.glb_sum(sum(values_list))
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
def auto_contractor_meson_corr_psnk_psrc(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"auto-contract/{job_tag}/traj={traj}/meson_corr_psnk_psrc.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_corr()
    expr_names = get_cexpr_names(cexpr)
    total_site = ru.get_total_site(job_tag)
    psel = get_psel()
    fsel, fselc = get_fsel()
    xg_fsel_list = fsel.to_psel_local().to_list()
    xg_psel_list = psel.to_list()
    geo = q.Geometry(total_site, 1)
    total_volume = geo.total_volume()
    @q.timer
    def feval(xg_src):
        counts = np.zeros(total_site[3], dtype = complex)
        values = np.zeros((total_site[3], len(expr_names),), dtype = complex)
        for xg_snk in xg_fsel_list:
            t = (xg_snk[3] - xg_src[3]) % total_site[3]
            pd = {
                    "x_2" : ("point-snk", xg_snk,),
                    "x_1" : ("point", xg_src,),
                    }
            counts[t] += 1
            values[t] += eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop, is_only_total = "total")
        values = values.transpose() # values[expr_name, t_sep]
        return counts, values
    counts_list, values_list = zip(*q.parallel_map(q.get_n_processes(), feval, xg_psel_list))
    res_count = q.glb_sum(sum(counts_list))
    res_sum = q.glb_sum(sum(values_list))
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
                tr(gamma_5*S_l(t_1,t_2)*gamma_5*S_s(t_2,x_1)*S_s(x_1,t_1)), # term_Type1_0004
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.collect_op()
        return cexpr
    cexpr = q.pickle_cache_call(calc_cexpr, f"cache/auto_contractor_cexpr/meson_m-cexpr.pickle")
    q.displayln_info(display_cexpr_raw(cexpr))
    return cexpr

@q.timer_verbose
def auto_contractor_meson_m(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"auto-contract/{job_tag}/traj={traj}/meson_m.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_m()
    expr_names = get_cexpr_names(cexpr)
    total_site = ru.get_total_site(job_tag)
    psel = get_psel()
    fsel, fselc = get_fsel()
    xg_fsel_list = fsel.to_psel_local().to_list()
    xg_psel_list = psel.to_list()
    tsep = rup.dict_params[job_tag]["meson_tensor_tsep"]
    geo = q.Geometry(total_site, 1)
    total_volume = geo.total_volume()
    @q.timer
    def feval(xg_snk):
        t = xg_snk[3]
        t_1 = (t + tsep) % total_site[3]
        t_2 = (t - tsep) % total_site[3]
        pd = {
                "x_1" : ("point-snk", xg_snk,),
                "t_1" : ("wall", t_1),
                "t_2" : ("wall", t_2),
                }
        values = eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop, is_only_total = "total")
        return 1.0, values
    counts_list, values_list = zip(*q.parallel_map(q.get_n_processes(), feval, xg_fsel_list))
    res_count = q.glb_sum(sum(counts_list))
    res_sum = q.glb_sum(sum(values_list))
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

def r_scaling_factor():
    return 5.0

def get_r(x_rel):
    fac = r_scaling_factor()
    return fac * math.sqrt(x_rel[0] * x_rel[0] + x_rel[1] * x_rel[1] + x_rel[2] * x_rel[2])

def get_r_limit(total_site):
    return math.ceil(get_r([ total_site[i] // 2 for i in range(4) ])) + 1

def jj_proj_mm(res_arr, x_rel):
    # res_arr is 3-D np.array with n*4*4 elements
    # res_arr[idx_tensor, mu, nu]
    return np.array([ np.trace(res_arr[idx_tensor]) for idx_tensor in range(len(res_arr)) ])

def jj_proj_tt(res_arr, x_rel):
    # res_arr is 3-D np.array with n*4*4 elements
    # res_arr[idx_tensor, mu, nu]
    return np.array([ res_arr[idx_tensor, 3, 3] for idx_tensor in range(len(res_arr)) ])

def jj_proj_ii(res_arr, x_rel):
    # res_arr is 3-D np.array with n*4*4 elements
    # res_arr[idx_tensor, mu, nu]
    return np.array([ sum([ res_arr[idx_tensor, i, i] for i in range(3) ]) for idx_tensor in range(len(res_arr)) ])

def jj_proj_xx(res_arr, x_rel):
    # res_arr is 3-D np.array with n*4*4 elements
    # res_arr[idx_tensor, mu, nu]
    return np.array(
            [ sum(
                [ x_rel[i] * x_rel[j] * res_arr[idx_tensor, i, j]
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

@q.timer
def accumulate_meson_jj(counts, values, res_arr, x_rel, total_site):
    # counts[t, r_idx]
    # values[idx_proj, idx_tensor, t, r_idx]
    # 0 <= idx_proj < len(all_jj_projections): mm, tt, ii, xx
    # 0 <= idx_tensor < n_tensor. ( n_tensor = 14 )
    # 0 <= t < total_site[3]
    # 0 <= r < r_limit (scale by factor of 5.0)
    (n_proj, n_tensor, t_size, r_limit,) = values.shape
    assert (t_size, r_limit,) == counts.shape
    assert t_size == total_site[3]
    assert r_limit == get_r_limit(total_site)
    assert len(all_jj_projections) == n_proj
    t = x_rel[3] % total_site[3]
    r = get_r(x_rel)
    r_idx_low = math.floor(r)
    r_idx_high = math.ceil(r)
    if r_idx_high == r_idx_low:
        r_idx_high += 1
    assert r_idx_high < r_limit
    coef_low = r_idx_high - r
    coef_high = r - r_idx_low
    counts[t, r_idx_low] += coef_low
    counts[t, r_idx_high] += coef_high
    for idx_proj, proj in enumerate(all_jj_projections):
        v = proj(res_arr, x_rel)
        v_low = coef_low * v
        v_high = coef_high * v
        for idx_tensor in range(n_tensor):
            values[idx_proj, idx_tensor, t, r_idx_low, ] += v_low[idx_tensor]
            values[idx_proj, idx_tensor, t, r_idx_high, ] += v_high[idx_tensor]

@q.timer
def get_cexpr_meson_jj():
    def calc_cexpr():
        t_1, t_2, x_1, x_2 = ['t_1', 't_2', 'x_1', 'x_2']
        def mk_terms(mu, nu):
            return [
                    tr(gamma(mu)*S_l(x_1,t_1)*gamma_5*S_l(t_1,x_1))*tr(gamma(nu)*S_l(x_2,t_2)*gamma_5*S_l(t_2,x_2)), # term_Type1_0001
                    tr(gamma(mu)*S_l(x_1,t_2)*gamma_5*S_l(t_2,x_1))*tr(gamma(nu)*S_l(x_2,t_1)*gamma_5*S_l(t_1,x_2)), # term_Type1_0002
                    tr(gamma(mu)*S_l(x_1,t_1)*gamma_5*S_l(t_1,x_2)*gamma(nu)*S_l(x_2,t_2)*gamma_5*S_l(t_2,x_1)), # term_Type2_0001
                    tr(gamma(mu)*S_l(x_1,t_2)*gamma_5*S_l(t_2,x_2)*gamma(nu)*S_l(x_2,t_1)*gamma_5*S_l(t_1,x_1)), # term_Type2_0002
                    tr(gamma(mu)*S_l(x_1,t_2)*gamma_5*S_s(t_2,x_2)*gamma(nu)*S_s(x_2,t_1)*gamma_5*S_l(t_1,x_1)), # term_Type2_0003
                    tr(gamma(mu)*S_s(x_1,t_1)*gamma_5*S_l(t_1,x_2)*gamma(nu)*S_l(x_2,t_2)*gamma_5*S_s(t_2,x_1)), # term_Type2_0004
                    tr(gamma(mu)*S_l(x_1,t_1)*gamma_5*S_l(t_1,t_2)*gamma_5*S_l(t_2,x_2)*gamma(nu)*S_l(x_2,x_1)), # term_Type3_0001
                    tr(gamma(mu)*S_l(x_1,t_2)*gamma_5*S_l(t_2,t_1)*gamma_5*S_l(t_1,x_2)*gamma(nu)*S_l(x_2,x_1)), # term_Type3_0002
                    tr(gamma(mu)*S_l(x_1,x_2)*gamma(nu)*S_l(x_2,t_1)*gamma_5*S_l(t_1,t_2)*gamma_5*S_l(t_2,x_1)), # term_Type3_0003
                    tr(gamma(mu)*S_l(x_1,x_2)*gamma(nu)*S_l(x_2,t_2)*gamma_5*S_l(t_2,t_1)*gamma_5*S_l(t_1,x_1)), # term_Type3_0004
                    tr(gamma(mu)*S_l(x_1,t_2)*gamma_5*S_s(t_2,t_1)*gamma_5*S_l(t_1,x_2)*gamma(nu)*S_l(x_2,x_1)), # term_Type3_0005
                    tr(gamma(mu)*S_l(x_1,x_2)*gamma(nu)*S_l(x_2,t_2)*gamma_5*S_s(t_2,t_1)*gamma_5*S_l(t_1,x_1)), # term_Type3_0006
                    tr(gamma(mu)*S_s(x_1,t_1)*gamma_5*S_l(t_1,t_2)*gamma_5*S_s(t_2,x_2)*gamma(nu)*S_s(x_2,x_1)), # term_Type3_0007
                    tr(gamma(mu)*S_s(x_1,x_2)*gamma(nu)*S_s(x_2,t_1)*gamma_5*S_l(t_1,t_2)*gamma_5*S_s(t_2,x_1)), # term_Type3_0008
                    ]
        terms_mu_nu = [ [ mk_terms(mu, nu) for nu in range(4) ] for mu in range(4) ]
        n_tensor = 14
        for t_nu in terms_mu_nu:
            for t in t_nu:
                assert n_tensor == len(t)
        exprs = [ terms_mu_nu[mu][nu][i] for i in range(n_tensor) for mu in range(4) for nu in range(4) ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.collect_op()
        return cexpr
    cexpr = q.pickle_cache_call(calc_cexpr, f"cache/auto_contractor_cexpr/meson_jj-cexpr.pickle")
    q.displayln_info(display_cexpr_raw(cexpr))
    return cexpr

@q.timer_verbose
def auto_contractor_meson_jj(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"auto-contract/{job_tag}/traj={traj}/meson_jj.lat"
    fn_counts = f"auto-contract/{job_tag}/traj={traj}/meson_jj_counts.lat"
    if (get_load_path(fn_counts) is not None) and (get_load_path(fn) is not None):
        return
    cexpr = get_cexpr_meson_jj()
    expr_names = get_cexpr_names(cexpr)
    total_site = ru.get_total_site(job_tag)
    psel = get_psel()
    fsel, fselc = get_fsel()
    xg_fsel_list = fsel.to_psel_local().to_list()
    xg_psel_list = psel.to_list()
    tsep = rup.dict_params[job_tag]["meson_tensor_tsep"]
    geo = q.Geometry(total_site, 1)
    t_size = total_site[3]
    r_limit = get_r_limit(total_site)
    n_proj = len(all_jj_projections)
    assert n_proj == len(all_jj_projection_names)
    n_tensor = len(expr_names) // 16
    assert n_tensor * 16 == len(expr_names)
    @q.timer_verbose
    def feval(xg_src):
        counts = np.zeros((t_size, r_limit,), dtype = complex)
        values = np.zeros((n_proj, n_tensor, t_size, r_limit,), dtype = complex)
        for xg_snk in xg_fsel_list:
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
            res = eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop, is_only_total = "total")
            assert res.shape[0] == 16 * n_tensor
            res_arr = res.reshape((n_tensor, 4, 4))
            accumulate_meson_jj(counts, values, res_arr, x_rel, total_site)
        return counts, values
    counts_list, values_list = zip(*q.parallel_map(q.get_n_processes(), feval, xg_psel_list))
    res_count = q.glb_sum(sum(counts_list))
    res_sum = q.glb_sum(sum(values_list))
    res_count *= 1.0 / (len(xg_psel_list) * fsel.prob())
    res_sum *= 1.0 / (len(xg_psel_list) * fsel.prob())
    ld_count = q.mk_lat_data([
        [ "t", t_size, [ str(rel_mod(t, t_size)) for t in range(t_size) ], ],
        [ "r", r_limit, ],
        ])
    ld_sum = q.mk_lat_data([
        [ "idx_proj", len(all_jj_projection_names), all_jj_projection_names, ],
        [ "idx_tensor", n_tensor, ],
        [ "t", t_size, [ str(rel_mod(t, t_size)) for t in range(t_size) ], ],
        [ "r", r_limit, ],
        ])
    ld_count.from_numpy(res_count)
    ld_sum.from_numpy(res_sum)
    ld_count.save(get_save_path(fn_counts))
    ld_sum.save(get_save_path(fn))
    q.displayln_info(ld_count.show())
    q.displayln_info(ld_sum.show())

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
            f"prop-rand-u1-light/{job_tag}/traj={traj}/geon-info.txt",
            f"prop-rand-u1-strange/{job_tag}/traj={traj}/geon-info.txt",
            f"prop-rand-u1-charm/{job_tag}/traj={traj}/geon-info.txt",
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
            get_prop = get_get_prop()
            # ADJUST ME
            # auto_contractor_vev(job_tag, traj, get_prop, get_psel, get_fsel)
            # auto_contractor_meson_f_corr(job_tag, traj, get_prop, get_psel, get_fsel)
            # auto_contractor_hvp(job_tag, traj, get_prop, get_psel, get_fsel)
            # auto_contractor_hvp_field(job_tag, traj, get_prop, get_psel, get_fsel)
            # auto_contractor_meson_v_v_meson_field(job_tag, traj, get_prop, get_psel, get_fsel)
            #
            auto_contractor_meson_corr(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contractor_meson_corr_psnk(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contractor_meson_corr_psrc(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contractor_meson_corr_psnk_psrc(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contractor_meson_m(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contractor_meson_jj(job_tag, traj, get_prop, get_psel, get_fsel)
            #
            q.qtouch_info(get_save_path(fn_checkpoint))
            q.release_lock()
    #
    q.clean_cache()
    q.timer_display()

qg.begin_with_gpt()

q.qremove_all_info("locks")
q.qremove_all_info("results")
q.qremove_all_info("cache")

# run_job("16IH2", 1000)
run_job("test-4nt8", 1000)

qg.end_with_gpt()
