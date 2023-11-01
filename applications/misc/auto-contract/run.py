#!/usr/bin/env python3

from auto_contractor.operators import *

import functools
import math
import os
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
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-fsel-self-loop/results"),
        # os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-selected-data/results"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-wsrc-prop/results"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-psrc-prop/results"),
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
    return cache_compiled_cexpr(calc_cexpr, f"cache/auto_contract_cexpr/vev-cexpr")

@q.timer_verbose
def auto_contract_vev(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"{job_tag}/auto-contract/traj-{traj}/vev.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_vev()
    expr_names = get_expr_names(cexpr)
    fsel, fselc = get_fsel()
    xg_fsel_list = list(map(tuple, fsel.to_psel_local().to_list()))
    @q.timer
    def feval(x):
        pd = {
                "x" : ("point-snk", x,),
                }
        res = eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop)
        return q.Data([ 1, res, ])
    q.timer_fork(0)
    res_count, res_sum = q.glb_sum(
            q.parallel_map_sum(feval, xg_fsel_list, chunksize = 16)).get_val()
    q.displayln_info("timer_display for auto_contract_vev")
    q.timer_display()
    q.timer_merge()
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
    return cache_compiled_cexpr(calc_cexpr, f"cache/auto_contract_cexpr/meson_f_corr-cexpr")

@q.timer_verbose
def auto_contract_meson_f_corr(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_f_corr.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_f_corr()
    expr_names = get_expr_names(cexpr)
    total_site = rup.get_total_site(job_tag)
    fsel, fselc = get_fsel()
    xg_fsel_list = list(map(tuple, fsel.to_psel_local().to_list()))
    @q.timer
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
    res_sum = q.glb_sum(q.parallel_map_sum(feval, xg_fsel_list))
    res_count = q.glb_sum(len(xg_fsel_list))
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
    return cache_compiled_cexpr(calc_cexpr, f"cache/auto_contract_cexpr/hvp-cexpr")

@q.timer_verbose
def auto_contract_hvp(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"{job_tag}/auto-contract/traj-{traj}/hvp.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_hvp()
    expr_names = get_expr_names(cexpr)
    total_site = rup.get_total_site(job_tag)
    psel = get_psel()
    fsel, fselc = get_fsel()
    xg_fsel_list = list(map(tuple, fsel.to_psel_local().to_list()))
    xg_psel_list = list(map(tuple, psel.to_list()))
    vol = total_site[0] * total_site[1] * total_site[2]
    @q.timer
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
        return q.Data([ counts, values, ])
    res_count, res_sum = q.glb_sum(q.parallel_map_sum(feval, xg_psel_list)).get_val()
    res_avg = res_sum * (vol / res_count)
    ld = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t_sep", total_site[3], ],
        ])
    ld.from_numpy(res_avg)
    q.displayln_info(ld.show())
    ld.save(get_save_path(fn))

@q.timer_verbose
def auto_contract_hvp_field(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"{job_tag}/auto-contract/traj-{traj}/hvp.field"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_hvp()
    expr_names = get_expr_names(cexpr)
    total_site = rup.get_total_site(job_tag)
    psel = get_psel()
    fsel, fselc = get_fsel()
    xg_fsel_list = list(map(tuple, fsel.to_psel_local().to_list()))
    xg_psel_list = list(map(tuple, psel.to_list()))
    geo = q.Geometry(total_site, 1)
    field = q.Field("ComplexD", geo, len(expr_names))
    field.set_zero()
    for idx, xg_src in enumerate(xg_psel_list):
        @q.timer
        def feval(xg_snk):
            pd = {
                    "x2" : ("point-snk", xg_snk,),
                    "x1" : ("point", xg_src,),
                    }
            res = eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop, is_only_total = "total")
            return res
        def sum_function(val_list, start):
            values = start
            count = 0
            for idx, v in enumerate(val_list):
                count += 1
                values[idx] = v
            assert count == fsel.n_elems()
            shift = [ -x for x in xg_src ]
            values_shifted = values.shift(shift)
            return values_shifted
        values_shifted = q.parallel_map_sum(feval, xg_fsel_list,
                sum_function = sum_function, sum_start = q.SelectedField("ComplexD", fsel, len(expr_names)))
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
    return cache_compiled_cexpr(calc_cexpr, f"cache/auto_contract_cexpr/meson_v_v_meson-cexpr")

@q.timer_verbose
def auto_contract_meson_v_v_meson_field(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_v_v_meson.field"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_v_v_meson()
    expr_names = get_expr_names(cexpr)
    total_site = rup.get_total_site(job_tag)
    psel = get_psel()
    fsel, fselc = get_fsel()
    xg_fsel_list = list(map(tuple, fsel.to_psel_local().to_list()))
    xg_psel_list = list(map(tuple, psel.to_list()))
    tsep = rup.dict_params[job_tag]["meson_tensor_tsep"]
    geo = q.Geometry(total_site, 1)
    field = q.Field("ComplexD", geo, len(expr_names))
    field.set_zero()
    for idx, xg_src in enumerate(xg_psel_list):
        @q.timer
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
        def sum_function(val_list, start):
            values = start
            count = 0
            for idx, v in enumerate(val_list):
                count += 1
                values[idx] = v.tobytes()
            assert count == fsel.n_elems()
            shift = [ -x for x in xg_src ]
            values_shifted = values.shift(shift)
            return values_shifted
        values_shifted = q.parallel_map_sum(feval, xg_fsel_list,
                sum_function = sum_function, sum_start = q.SelectedField("ComplexD", fsel, len(expr_names)))
        field += values_shifted
    field_r = q.Field("ComplexD", geo, len(expr_names))
    field_r.set_zero()
    for idx, xg_src in enumerate(xg_psel_list):
        @q.timer
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
        def sum_function(val_list, start):
            values = start
            count = 0
            for idx, v in enumerate(val_list):
                count += 1
                values[idx] = v.tobytes()
            assert count == fsel.n_elems()
            shift = [ -x for x in xg_src ]
            values_shifted = values.shift(shift)
            return values_shifted
        values_shifted = q.parallel_map_sum(feval, xg_fsel_list,
                sum_function = sum_function, sum_start = q.SelectedField("ComplexD", fsel, len(expr_names)))
        field_r += values_shifted
    field_r = field_r.reflect()
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
    return cache_compiled_cexpr(calc_cexpr, f"cache/auto_contract_cexpr/meson_corr-cexpr")

@q.timer_verbose
def auto_contract_meson_corr(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_corr.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_corr()
    expr_names = get_expr_names(cexpr)
    total_site = rup.get_total_site(job_tag)
    fsel, fselc = get_fsel()
    xg_fsel_list = list(map(tuple, fsel.to_psel_local().to_list()))
    geo = q.Geometry(total_site, 1)
    total_volume = geo.total_volume()
    def fempty():
        values = np.zeros((len(expr_names), total_site[3],), dtype = complex)
        return q.Data([ 0, values, ])
    @q.timer
    def feval(t_snk):
        counts, values = fempty().get_val()
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
        return q.Data([ counts, values, ])
    t_snk_list = get_mpi_chunk(list(range(total_site[3])))
    res_count, res_sum = q.glb_sum(q.parallel_map_sum(feval, t_snk_list, sum_start = fempty())).get_val()
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
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_corr_psnk.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_corr()
    expr_names = get_expr_names(cexpr)
    total_site = rup.get_total_site(job_tag)
    fsel, fselc = get_fsel()
    xg_fsel_list = list(map(tuple, fsel.to_psel_local().to_list()))
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
        return res
    res_count = q.glb_sum(len(xg_fsel_list))
    res_sum = q.glb_sum(q.parallel_map_sum(feval, xg_fsel_list))
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
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_corr_psrc.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_corr()
    expr_names = get_expr_names(cexpr)
    total_site = rup.get_total_site(job_tag)
    psel = get_psel()
    fsel, fselc = get_fsel()
    xg_fsel_list = list(map(tuple, fsel.to_psel_local().to_list()))
    xg_psel_list = list(map(tuple, psel.to_list()))
    geo = q.Geometry(total_site, 1)
    total_volume = geo.total_volume()
    def fempty():
        counts = 0
        values = np.zeros((len(expr_names), total_site[3],), dtype = complex)
        return q.Data([ counts, values, ])
    @q.timer
    def feval(xg_src):
        counts, values = fempty().get_val()
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
        return q.Data([ counts, values, ])
    xg_src_list = get_mpi_chunk(xg_psel_list, rng_state = q.RngState("get_mpi_chunk"))
    res_count, res_sum = q.glb_sum(q.parallel_map_sum(feval, xg_src_list, sum_start = fempty())).get_val()
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
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_corr_psnk_psrc.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_corr()
    expr_names = get_expr_names(cexpr)
    total_site = rup.get_total_site(job_tag)
    psel = get_psel()
    fsel, fselc = get_fsel()
    xg_fsel_list = list(map(tuple, fsel.to_psel_local().to_list()))
    xg_psel_list = list(map(tuple, psel.to_list()))
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
        return q.Data([ counts, values, ])
    res_count, res_sum = q.glb_sum(q.parallel_map_sum(feval, xg_psel_list)).get_val()
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

### ------

@q.timer_verbose
def run_job(job_tag, traj):
    fns_produce = [
            f"{job_tag}/auto-contract/traj-{traj}/checkpoint.txt",
            ]
    fns_need = [
            # (f"{job_tag}/configs/ckpoint_lat.{traj}", f"{job_tag}/configs/ckpoint_lat.IEEE64BIG.{traj}",),
            f"{job_tag}/point-selection/traj-{traj}.txt",
            f"{job_tag}/field-selection/traj-{traj}.field",
            f"{job_tag}/gauge-transform/traj-{traj}.field",
            f"{job_tag}/wall-src-info-light/traj-{traj}.txt",
            f"{job_tag}/wall-src-info-strange/traj-{traj}.txt",
            f"{job_tag}/psel-prop-wsrc-light/traj-{traj}/checkpoint.txt",
            f"{job_tag}/psel-prop-wsrc-strange/traj-{traj}/checkpoint.txt",
            f"{job_tag}/psel-prop-psrc-light/traj-{traj}/checkpoint.txt",
            f"{job_tag}/psel-prop-psrc-strange/traj-{traj}/checkpoint.txt",
            f"{job_tag}/prop-wsrc-light/traj-{traj}/geon-info.txt",
            f"{job_tag}/prop-wsrc-strange/traj-{traj}/geon-info.txt",
            f"{job_tag}/prop-psrc-light/traj-{traj}/geon-info.txt",
            f"{job_tag}/prop-psrc-strange/traj-{traj}/geon-info.txt",
            f"{job_tag}/prop-rand-u1-light/traj-{traj}/geon-info.txt",
            f"{job_tag}/prop-rand-u1-strange/traj-{traj}/geon-info.txt",
            f"{job_tag}/prop-rand-u1-charm/traj-{traj}/geon-info.txt",
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
    fn_checkpoint = f"{job_tag}/auto-contract/traj-{traj}/checkpoint.txt"
    if get_load_path(fn_checkpoint) is None:
        if q.obtain_lock(f"locks/{job_tag}-{traj}-auto-contract"):
            get_prop = get_get_prop()
            # ADJUST ME
            auto_contract_vev(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contract_meson_f_corr(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contract_hvp(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contract_hvp_field(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contract_meson_v_v_meson_field(job_tag, traj, get_prop, get_psel, get_fsel)
            #
            auto_contract_meson_corr(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contract_meson_corr_psnk(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contract_meson_corr_psrc(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contract_meson_corr_psnk_psrc(job_tag, traj, get_prop, get_psel, get_fsel)
            #
            q.qtouch_info(get_save_path(fn_checkpoint))
            q.release_lock()
    #
    q.clean_cache()
    q.timer_display()

def get_all_cexpr():
    benchmark_eval_cexpr(get_cexpr_vev())
    benchmark_eval_cexpr(get_cexpr_meson_f_corr())
    benchmark_eval_cexpr(get_cexpr_hvp())
    benchmark_eval_cexpr(get_cexpr_meson_v_v_meson())
    benchmark_eval_cexpr(get_cexpr_meson_corr())

def test():
    q.qremove_all_info("locks")
    q.qremove_all_info("results")
    run_job("test-4nt8", 1000)
    # run_job("16IH2", 1000)

size_node_list = [
        [ 1, 1, 1, 1, ],
        [ 1, 1, 1, 2, ],
        [ 1, 1, 1, 4, ],
        [ 1, 1, 1, 8, ],
        ]

q.begin_with_mpi(sys.argv, size_node_list)

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

q.end_with_mpi()
