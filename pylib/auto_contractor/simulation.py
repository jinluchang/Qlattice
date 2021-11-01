#    Qlattice (https://github.com/waterret/qlattice)
#
#    Copyright (C) 2021
#
#    Author: Luchang Jin (ljin.luchang@gmail.com)
#    Author: Masaaki Tomii
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to the Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

import qlat_gpt as qg
import qlat as q
import gpt as g
import rbc_ukqcd as ru
import rbc_ukqcd_params as rup
import numpy as np
import pprint
import os

from auto_contractor.eval import *
from auto_contractor.operators import *

def get_save_path(fn):
    return os.path.join("results", fn)

def get_load_path(fn):
    if fn is None:
        return None
    path_list = [ "results" ]
    for path in path_list:
        p = os.path.join(path, fn)
        if q.does_file_exist_sync_node(p):
            return p
    return None

@q.timer
def compute_prop(inv, src, *, tag, sfw):
    sol = inv * src
    sol.save_double(sfw, tag)
    return sol

@q.timer
def compute_prop_psrc(gf, xg, job_tag, inv_type, inv_acc, *, idx, sfw, eig, finished_tags):
    tag = f"xg=({xg[0]},{xg[1]},{xg[2]},{xg[3]}) ; type={inv_type} ; accuracy={inv_acc}"
    if tag in finished_tags:
        return None
    q.check_stop()
    q.check_time_limit()
    q.displayln_info(f"compute_prop_psrc: idx={idx} xg={xg}", job_tag, inv_type, inv_acc)
    inv = ru.get_inv(gf, job_tag, inv_type, inv_acc, eig = eig)
    total_site = ru.get_total_site(job_tag)
    geo = q.Geometry(total_site, 1)
    src = q.mk_point_src(geo, xg)
    compute_prop(inv, src, tag = tag, sfw = sfw)

def get_all_points(total_site, *, tslice = None):
    all_points = []
    if tslice is None:
        for t in range(total_site[3]):
            all_points += get_all_points(total_site, tslice = t)
        return all_points
    t = tslice
    for x in range(total_site[0]):
        for y in range(total_site[1]):
            for z in range(total_site[2]):
                all_points.append([x, y, z, t,])
    return all_points

@q.timer
def compute_prop_psrc_all(gf, job_tag, inv_type, *, path_s, eig):
    finished_tags = q.properly_truncate_fields(get_save_path(path_s + ".acc"))
    sfw = q.open_fields(get_save_path(path_s + ".acc"), "a", [ 1, 1, 1, 2 ])
    total_site = ru.get_total_site(job_tag)
    inv_acc = 2
    for idx, xg in enumerate(get_all_points(total_site)):
        compute_prop_psrc(gf, xg, job_tag, inv_type, inv_acc,
                idx = idx, sfw = sfw, eig = eig,
                finished_tags = finished_tags)
    q.clean_cache(q.cache_inv)
    sfw.close()
    q.qrename_info(get_save_path(path_s + ".acc"), get_save_path(path_s))

@q.timer
def compute_prop_wsrc(gf, gt, tslice, job_tag, inv_type, inv_acc, *, idx, sfw, path_sp, eig, finished_tags):
    tag = f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc}"
    if tag in finished_tags:
        return None
    q.check_stop()
    q.check_time_limit()
    q.displayln_info(f"compute_prop_wsrc: idx={idx} tslice={tslice}", job_tag, inv_type, inv_acc)
    inv = ru.get_inv(gf, job_tag, inv_type, inv_acc, gt = gt, eig = eig)
    total_site = ru.get_total_site(job_tag)
    geo = q.Geometry(total_site, 1)
    src = q.mk_wall_src(geo, tslice)
    prop = compute_prop(inv, src, tag = tag, sfw = sfw)
    fn_spw = os.path.join(path_sp, f"{tag} ; wsnk.lat")
    prop.glb_sum_tslice().save(get_save_path(fn_spw))

def get_all_walls(time_vol):
    all_walls = []
    for t in range(time_vol):
        all_walls.append(t)
    return all_walls

@q.timer
def compute_prop_wsrc_all(gf, gt, job_tag, inv_type, *, path_s, path_sp, eig):
    finished_tags = q.properly_truncate_fields(get_save_path(path_s + ".acc"))
    sfw = q.open_fields(get_save_path(path_s + ".acc"), "a", [ 1, 1, 1, 2 ])
    total_site = ru.get_total_site(job_tag)
    inv_acc = 2
    for idx, tslice in enumerate(get_all_walls(total_site[3])):
        compute_prop_wsrc(gf, gt, tslice, job_tag, inv_type, inv_acc,
                idx = idx, sfw = sfw, path_sp = path_sp, eig = eig,
                finished_tags = finished_tags)
    q.clean_cache(q.cache_inv)
    sfw.close()
    q.qrename_info(get_save_path(path_s + ".acc"), get_save_path(path_s))

@q.timer
def load_prop_psrc_all(job_tag, traj, flavor : str, path_s : str, psel):
    cache = q.mk_cache(f"prop_cache-{job_tag}-{traj}", flavor)
    total_site = ru.get_total_site(job_tag)
    if flavor in ["l", "u", "d",]:
        inv_type = 0
    elif flavor in ["s",]:
        inv_type = 1
    elif flavor in ["c",]:
        inv_type = 2
    else:
        assert False
    inv_acc = 2
    sfr = q.open_fields(get_load_path(path_s), "r")
    for idx, xg in enumerate(get_all_points(total_site)):
        q.displayln_info(f"load_prop_psrc_all: idx={idx} xg={xg} flavor={flavor} path_s={path_s}")
        prop = q.Prop()
        tag = f"xg=({xg[0]},{xg[1]},{xg[2]},{xg[3]}) ; type={inv_type} ; accuracy={inv_acc}"
        prop.load_double(sfr, tag)
        # convert to GPT/Grid prop mspincolor order
        prop_msc = q.convert_mspincolor_from_wm(prop)
        sp_prop_msc = q.PselProp(psel)
        sp_prop_msc @= prop_msc
        cache[f"xg=({xg[0]},{xg[1]},{xg[2]},{xg[3]})"] = sp_prop_msc
    sfr.close()

@q.timer
def load_prop_wsrc_all(job_tag, traj, flavor : str, path_s : str, path_sp : str, psel, gt):
    cache = q.mk_cache(f"prop_cache-{job_tag}-{traj}", flavor)
    total_site = ru.get_total_site(job_tag)
    psel_ts = q.get_psel_tslice(total_site)
    if flavor in ["l", "u", "d",]:
        inv_type = 0
    elif flavor in ["s",]:
        inv_type = 1
    elif flavor in ["c",]:
        inv_type = 2
    else:
        assert False
    inv_acc = 2
    gt_inv = gt.inv()
    sfr = q.open_fields(get_load_path(path_s), "r")
    for idx, tslice in enumerate(get_all_walls(total_site[3])):
        q.displayln_info(f"load_prop_wsrc_all: idx={idx} tslice={tslice} flavor={flavor} path_s={path_s} path_sp={path_sp}")
        prop = q.Prop()
        tag = f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc}"
        prop.load_double(sfr, tag)
        prop = gt_inv * prop
        # convert to GPT/Grid prop mspincolor order
        prop_msc = q.convert_mspincolor_from_wm(prop)
        sp_prop_msc = q.PselProp(psel)
        sp_prop_msc @= prop_msc
        cache[f"tslice={tslice}"] = sp_prop_msc
        # load wsnk prop
        fn_spw = os.path.join(path_sp, f"{tag} ; wsnk.lat")
        spw_prop = q.PselProp(psel_ts)
        spw_prop.load(get_load_path(fn_spw))
        spw_prop_msc = q.convert_mspincolor_from_wm(spw_prop)
        cache[f"tslice={tslice} ; wsnk"] = spw_prop_msc
    sfr.close()

@q.timer
def get_prop_psrc(prop_cache, flavor : str, xg_src):
    # prop_cache[flavor][src_p] = prop
    # call load_prop_psrc_all(flavor, path_s) first
    return prop_cache[flavor][f"xg=({xg_src[0]},{xg_src[1]},{xg_src[2]},{xg_src[3]})"]

@q.timer
def get_prop_wsrc(prop_cache, flavor : str, t_src):
    # prop_cache[flavor][src_p] = prop
    # call load_prop_wsrc_all(flavor, path_s) first
    return prop_cache[flavor][f"tslice={t_src}"]

@q.timer
def get_prop_wsnk_wsrc(prop_cache, flavor : str, t_snk, t_src):
    # prop_cache[flavor][src_p] = prop
    # call load_prop_wsrc_all(flavor, path_s) first
    return prop_cache[flavor][f"tslice={t_src} ; wsnk"].get_elem(t_snk)

@q.timer
def get_prop_psnk_psrc(prop_cache, flavor : str, p_snk, p_src, *, psel_pos_dict):
    if isinstance(p_snk, list) and isinstance(p_src, list):
        assert 4 == len(p_snk)
        assert 4 == len(p_src)
        xg_snk = p_snk
        xg_src = p_src
        msc = get_prop_psrc(prop_cache, flavor, xg_src).get_elem(psel_pos_dict[tuple(xg_snk)])
    else:
        assert isinstance(p_snk, tuple) and isinstance(p_src, tuple)
        assert 2 == len(p_snk)
        assert 2 == len(p_src)
        type_snk, pos_snk = p_snk
        type_src, pos_src = p_src
        if type_snk == "wall" and type_src == "wall":
            msc = get_prop_wsnk_wsrc(prop_cache, flavor, pos_snk, pos_src)
        elif type_snk == "point" and type_src == "wall":
            msc = get_prop_wsrc(prop_cache, flavor, pos_src).get_elem(psel_pos_dict[tuple(pos_snk)])
        elif type_snk == "wall" and type_src == "point":
            msc = ascontiguoustensor(
                    ascontiguoustensor(
                        g.gamma[5]
                        * adj_msc(
                            as_mspincolor(
                                get_prop_wsrc(prop_cache, flavor, pos_snk).get_elem(psel_pos_dict[tuple(pos_src)]))))
                    * g.gamma[5])
        elif type_snk == "point" and type_src == "point":
            msc = get_prop_psrc(prop_cache, flavor, pos_src).get_elem(psel_pos_dict[tuple(pos_snk)])
        else:
            raise Exception("get_prop_psnk_psrc unknown p_snk={p_snk} p_src={p_src}")
    return as_mspincolor(msc)

def mk_get_prop(prop_cache, *, psel_pos_dict):
    def get_prop(flavor, p_snk, p_src):
        return get_prop_psnk_psrc(prop_cache, flavor, p_snk, p_src, psel_pos_dict = psel_pos_dict)
    return get_prop

@q.timer
def auto_contractor_simple_test(job_tag, traj, get_prop):
    q.displayln_info(g.gamma[5] * get_prop("l", [1, 2, 3, 2], [1, 2, 3, 2]))
    q.displayln_info(g.trace(g.gamma[5] * get_prop("l", [1, 2, 3, 2], [1, 2, 3, 2])))
    q.displayln_info(g.gamma[5] * get_prop("l", [1, 2, 3, 2], [1, 2, 3, 2]) - get_prop("l", [1, 2, 3, 2], [1, 2, 3, 2]))
    q.displayln_info(g.gamma[5] * adj_msc(get_prop("l", [1, 2, 3, 2], [1, 2, 3, 2])))
    q.displayln_info(g.norm2(g.gamma[5] * adj_msc(get_prop("l", [1, 2, 3, 2], [1, 2, 3, 4])) * g.gamma[5] - get_prop("l", [1, 2, 3, 4], [1, 2, 3, 2])))
    expr = (1
            * Qb("d", "x1", "s1", "c1")
            * G(5, "s1", "s2")
            * Qv("u", "x1", "s2", "c1")
            * Qb("u", "x2", "s3", "c2")
            * G(5, "s3", "s4")
            * Qv("d", "x2", "s4", "c2"))
    expr = contract_expr(expr)
    expr.simplify(is_isospin_symmetric_limit = True)
    cexpr = mk_cexpr(expr)
    cexpr.collect_op()
    q.displayln_info(cexpr)
    positions_dict = {}
    positions_dict["x1"] = [1, 2, 3, 4]
    positions_dict["x2"] = [1, 2, 3, 2]
    val = eval_cexpr(cexpr, positions_dict = positions_dict, get_prop = get_prop, is_only_total = "total")
    q.displayln_info("eval_cexpr: ", val)
    q.displayln_info("gpt_direct: ",
            -g.trace(
                g.gamma[5]
                * get_prop("l", [1, 2, 3, 4], [1, 2, 3, 2])
                * g.gamma[5]
                * get_prop("l", [1, 2, 3, 2], [1, 2, 3, 4])
                ))

@q.timer
def get_cexpr_vev():
    def calc_cexpr():
        s = new_spin_index()
        c = new_color_index()
        p = "x"
        exprs = [
                Qb("u", "x", s, c) * Qv("u", "x", s, c) + "u_bar*u",
                Qb("d", "x", s, c) * Qv("d", "x", s, c) + "d_bar*d",
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

@q.timer
def auto_contractor_vev(job_tag, traj, get_prop, num_trials):
    total_site = ru.get_total_site(job_tag)
    vol = total_site[0] * total_site[1] * total_site[2]
    cexpr = get_cexpr_vev()
    rng_state = q.RngState("seed")
    def positions_dict_maker(idx):
        rs = rng_state.split(str(idx))
        x = rs.c_rand_gen(total_site)
        pd = {
                "x" : ("point", x,),
                }
        facs = [ 1.0, ]
        return pd, facs
    names_fac = [ "rest", ]
    trial_indices = get_mpi_chunk(range(num_trials))
    results_list = eval_cexpr_simulation(
            cexpr,
            positions_dict_maker = positions_dict_maker,
            trial_indices = trial_indices,
            get_prop = get_prop,
            is_only_total = "total")
    for name_fac, results in zip(names_fac, results_list):
        q.displayln_info(f"{name_fac} :")
        for k, v in results.items():
            q.displayln_info(f"{name_fac} {k}:\n  {v}")

@q.timer
def auto_contractor_meson_corr(job_tag, traj, get_prop, num_trials):
    total_site = ru.get_total_site(job_tag)
    vol = total_site[0] * total_site[1] * total_site[2]
    exprs = [
            vol**2 * mk_pi_p("x2", True) * mk_pi_p("x1"),
            vol**2 * mk_k_p("x2", True) * mk_k_p("x1"),
            vol**2 * mk_k_m("x2", True) * mk_k_m("x1"),
            vol * mk_pi_p("x2", True) * mk_pi_p("t1"),
            vol * mk_pi_p("t2", True) * mk_pi_p("x1"),
            mk_pi_p("t2", True) * mk_pi_p("t1"),
            ]
    cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True)
    q.displayln_info(display_cexpr(cexpr))
    cexpr.collect_op()
    q.displayln_info(display_cexpr_raw(cexpr))
    rng_state = q.RngState("seed")
    def positions_dict_maker(idx):
        rs = rng_state.split(str(idx))
        t2 = 5
        x1 = rs.c_rand_gen(total_site)
        x2 = rs.c_rand_gen(total_site)
        x2[3] = (x1[3] + t2) % total_site[3]
        pd = {
                "x1" : ("point", x1,),
                "x2" : ("point", x2,),
                "t1" : ("wall", x1[3],),
                "t2" : ("wall", x2[3],),
                }
        lmom = [ 2 * math.pi / total_site[i] for i in range(3) ]
        facs = [
                1.0,
                sum([ cmath.rect(1.0, (x2[i] - x1[i]) * lmom[i]).real for i in range(3) ]) / 3.0,
                sum([ cmath.rect(1.0, (x2[i] - x1[i]) * 2 * lmom[i]).real for i in range(3) ]) / 3.0,
                ]
        return pd, facs
    names_fac = ["rest", "mom1", "mom2",]
    trial_indices = get_mpi_chunk(range(num_trials))
    results_list = eval_cexpr_simulation(cexpr, positions_dict_maker = positions_dict_maker, trial_indices = trial_indices, get_prop = get_prop, is_only_total = "typed_total")
    for name_fac, results in zip(names_fac, results_list):
        q.displayln_info(f"{name_fac} :")
        for k, v in results.items():
            q.displayln_info(f"{name_fac} {k}:\n  {v}")

@q.timer
def get_cexpr_test_corr(vol):
    def calc_cexpr():
        exprs = [
                vol**2 * mk_k0pi0("x2_1", "x2_2", True) * mk_k0pi0("x1_1", "x1_2"),
                vol**2 * mk_k0pi0("x2_1", "x2_2", True) * mk_kmpip("x1_1", "x1_2"),
                vol**2 * mk_kmpip("x2_1", "x2_2", True) * mk_kmpip("x1_1", "x1_2"),
                vol**2 * mk_kppim("x2_1", "x2_2", True) * mk_kppim("x1_1", "x1_2"),
                vol**2 * mk_k0barpi0("x2_1", "x2_2", True) * mk_k0barpi0("x1_1", "x1_2"),
                vol**2 * mk_kpi_0_i1half("x2_1", "x2_2", True) * mk_kpi_0_i1half("x1_1", "x1_2"),
                vol**2 * mk_kpi_0_i1half("x2_1", "x2_2", True) * mk_kpi_m_i1half("x1_1", "x1_2"),
                vol**2 * mk_kpi_m_i1half("x2_1", "x2_2", True) * mk_kpi_m_i1half("x1_1", "x1_2"),
                vol**2 * mk_kpi_0_i1half("x2_1", "x2_2", True) * mk_kpi_0_i3halves("x1_1", "x1_2"),
                vol**2 * mk_kpi_m_i1half("x2_1", "x2_2", True) * mk_kpi_m1_i3halves("x1_1", "x1_2"),
                vol**2 * mk_kpi_0_i3halves("x2_1", "x2_2", True) * mk_kpi_0_i3halves("x1_1", "x1_2"),
                vol**2 * mk_k_0("x2_1", True) * mk_kpi_0_i1half("x1_1", "x1_2"),
                vol**2 * mk_k_0("x2_1", True) * mk_kpi_0_i3halves("x1_1", "x1_2"),
                vol**2 * mk_k_m("x2_1", True) * mk_kpi_m_i1half("x1_1", "x1_2"),
                vol**2 * mk_k_m("x2_1", True) * mk_kpi_m1_i3halves("x1_1", "x1_2"),
                vol**2 * mk_kpi_p_i3halves("x2_1", "x2_2", True) * mk_kpi_0_i3halves("x1_1", "x1_2"),
                vol**2 * mk_kpi_m1_i3halves("x2_1", "x2_2", True) * mk_kpi_0_i3halves("x1_1", "x1_2"),
                vol**2 * mk_kpi_m1_i3halves("x2_1", "x2_2", True) * mk_kpi_m1_i3halves("x1_1", "x1_2"),
                vol**2 * mk_kpi_m2_i3halves("x2_1", "x2_2", True) * mk_kpi_0_i3halves("x1_1", "x1_2"),
                vol**2 * mk_kpi_m2_i3halves("x2_1", "x2_2", True) * mk_kpi_m2_i3halves("x1_1", "x1_2"),
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.collect_op()
        return cexpr
    cexpr = q.pickle_cache_call(calc_cexpr, f"cache/auto_contractor_cexpr/test_corr-cexpr.{vol}.pickle")
    q.displayln_info(display_cexpr_raw(cexpr))
    return cexpr

@q.timer
def auto_contractor_test_corr(job_tag, traj, get_prop, num_trials):
    total_site = ru.get_total_site(job_tag)
    vol = total_site[0] * total_site[1] * total_site[2]
    cexpr = get_cexpr_test_corr(vol)
    rng_state = q.RngState("seed")
    def positions_dict_maker(idx):
        rs = rng_state.split(str(idx))
        t1_2 = 2
        t2_1 = 5
        t2_2 = 7
        x1_1 = rs.c_rand_gen(total_site)
        x1_2 = rs.c_rand_gen(total_site)
        x2_1 = rs.c_rand_gen(total_site)
        x2_2 = rs.c_rand_gen(total_site)
        x1_2[3] = (x1_1[3] + t1_2) % total_site[3]
        x2_1[3] = (x1_1[3] + t2_1) % total_site[3]
        x2_2[3] = (x1_1[3] + t2_2) % total_site[3]
        pd = {
                "x1_1" : ("point", x1_1,),
                "x1_2" : ("point", x1_2,),
                "x2_1" : ("point", x2_1,),
                "x2_2" : ("point", x2_2,),
                }
        lmom = [ 2 * math.pi / total_site[i] for i in range(3) ]
        fac1 = sum([ 2 * cmath.rect(1.0, (x1_1[i] - x1_2[i]) * lmom[i]).real for i in range(3) ]) / np.sqrt(6)
        fac2 = sum([ 2 * cmath.rect(1.0, (x2_1[i] - x2_2[i]) * lmom[i]).real for i in range(3) ]) / np.sqrt(6)
        facs = [1.0, fac1, fac2, fac1 * fac2,]
        return pd, facs
    names_fac = ["rest-rest", "rest-moving", "moving-rest", "moving-moving",]
    trial_indices = get_mpi_chunk(range(num_trials))
    results_list = eval_cexpr_simulation(cexpr, positions_dict_maker = positions_dict_maker, trial_indices = trial_indices, get_prop = get_prop, is_only_total = "typed_total")
    for name_fac, results in zip(names_fac, results_list):
        q.displayln_info(f"{name_fac} :")
        for k, v in results.items():
            q.displayln_info(f"{name_fac} {k}:\n  {v}")

@q.timer
def get_cexpr_pipi_corr(vol):
    def calc_cexpr():
        exprs = [
                vol**4 * mk_pipi_i22("x2_1", "x2_2", True) * mk_pipi_i22("x1_1", "x1_2"),
                vol**4 * mk_pipi_i11("x2_1", "x2_2", True) * mk_pipi_i11("x1_1", "x1_2"),
                vol**4 * mk_pipi_i0("x2_1", "x2_2", True) * mk_pipi_i0("x1_1", "x1_2"),
                vol**2 * mk_pipi_i0("x2_1", "x2_2", True),
                vol**2 * mk_pipi_i0("x1_1", "x1_2"),
                vol**2 * mk_sigma("x2_1", True) * mk_sigma("x1_1"),
                vol**2 * mk_sigma("x2_1", True) * mk_sigma("x1_2"),
                vol**2 * mk_sigma("x2_2", True) * mk_sigma("x1_1"),
                vol**2 * mk_sigma("x2_2", True) * mk_sigma("x1_2"),
                vol * mk_sigma("x1_1"),
                vol * mk_sigma("x1_2"),
                vol * mk_sigma("x2_1", True),
                vol * mk_sigma("x2_2", True),
                vol**3 * mk_sigma("x2_1", True) * mk_pipi_i0("x1_1", "x1_2"),
                vol**3 * mk_sigma("x2_2", True) * mk_pipi_i0("x1_1", "x1_2"),
                vol**3 * mk_pipi_i0("x2_1", "x2_2", True) * mk_sigma("x1_1"),
                vol**3 * mk_pipi_i0("x2_1", "x2_2", True) * mk_sigma("x1_2"),
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.collect_op()
        return cexpr
    cexpr = q.pickle_cache_call(calc_cexpr, f"cache/auto_contractor_cexpr/pipi_corr-cexpr.{vol}.pickle")
    q.displayln_info(display_cexpr_raw(cexpr))
    return cexpr

@q.timer
def auto_contractor_pipi_corr(job_tag, traj, get_prop, num_trials):
    total_site = ru.get_total_site(job_tag)
    vol = total_site[0] * total_site[1] * total_site[2]
    cexpr = get_cexpr_pipi_corr(vol)
    rng_state = q.RngState("seed")
    def positions_dict_maker(idx):
        rs = rng_state.split(str(idx))
        t1_2 = 2
        t2_1 = 5
        t2_2 = 7
        x1_1 = rs.c_rand_gen(total_site)
        x1_2 = rs.c_rand_gen(total_site)
        x2_1 = rs.c_rand_gen(total_site)
        x2_2 = rs.c_rand_gen(total_site)
        x1_2[3] = (x1_1[3] + t1_2) % total_site[3]
        x2_1[3] = (x1_1[3] + t2_1) % total_site[3]
        x2_2[3] = (x1_1[3] + t2_2) % total_site[3]
        pd = {
                "x1_1" : ("point", x1_1,),
                "x1_2" : ("point", x1_2,),
                "x2_1" : ("point", x2_1,),
                "x2_2" : ("point", x2_2,),
                }
        lmom = [ 2 * math.pi / total_site[i] for i in range(3) ]
        fac1 = sum([ 2 * cmath.rect(1.0, (x1_1[i] - x1_2[i]) * lmom[i]).real for i in range(3) ]) / np.sqrt(6)
        fac2 = sum([ 2 * cmath.rect(1.0, (x2_1[i] - x2_2[i]) * lmom[i]).real for i in range(3) ]) / np.sqrt(6)
        facs = [1.0, fac1, fac2, fac1 * fac2,]
        return pd, facs
    names_fac = ["rest-rest", "rest-moving", "moving-rest", "moving-moving",]
    trial_indices = get_mpi_chunk(range(num_trials))
    results_list = eval_cexpr_simulation(cexpr, positions_dict_maker = positions_dict_maker, trial_indices = trial_indices, get_prop = get_prop, is_only_total = "typed_total")
    for name_fac, results in zip(names_fac, results_list):
        q.displayln_info(f"{name_fac} :")
        for k, v in results.items():
            q.displayln_info(f"{name_fac} {k}:\n  {v}")

@q.timer
def get_cexpr_kpipi(vol):
    @q.timer
    def calc_cexpr():
        total_site = ru.get_total_site(job_tag)
        vol = total_site[0] * total_site[1] * total_site[2]
        exprs_odd_ops = [
                vol * mk_Q1("x", "odd") + "Q1(o)",
                vol * mk_Q2("x", "odd") + "Q2(o)",
                vol * mk_Q3("x", "odd") + "Q3(o)",
                vol * mk_Q4("x", "odd") + "Q4(o)",
                vol * mk_Q5("x", "odd") + "Q5(o)",
                vol * mk_Q6("x", "odd") + "Q6(o)",
                vol * mk_Q7("x", "odd") + "Q7(o)",
                vol * mk_Q8("x", "odd") + "Q8(o)",
                vol * mk_Q9("x", "odd") + "Q9(o)",
                vol * mk_Q10("x", "odd") + "Q10(o)",
                vol * mk_Qsub("x", "odd") + "Qs(o)",
                ]
        exprs_even_ops = [
                vol * mk_Q1("x", "even") + "Q1(e)",
                vol * mk_Q2("x", "even") + "Q2(e)",
                vol * mk_Q3("x", "even") + "Q3(e)",
                vol * mk_Q4("x", "even") + "Q4(e)",
                vol * mk_Q5("x", "even") + "Q5(e)",
                vol * mk_Q6("x", "even") + "Q6(e)",
                vol * mk_Q7("x", "even") + "Q7(e)",
                vol * mk_Q8("x", "even") + "Q8(e)",
                vol * mk_Q9("x", "even") + "Q9(e)",
                vol * mk_Q10("x", "even") + "Q10(e)",
                vol * mk_Qsub("x", "even") + "Qs(e)",
                ]
        exprs_ops = exprs_odd_ops + exprs_even_ops
        exprs_k = [
                vol * mk_k_0("x2") + "K0",
                ]
        exprs_pipi = [
                vol**2 * mk_pipi_i0("x1_1", "x1_2", True) + "pipi_I0",
                vol**2 * mk_pipi_i20("x1_1", "x1_2", True) + "pipi_I2",
                vol * mk_sigma("x1_1", True) + "sigma_1",
                vol * mk_sigma("x1_2", True) + "sigma_2",
                vol * mk_pi_0("x1_1", True) + "pi0_1",
                vol * mk_pi_0("x1_2", True) + "pi0_2",
                mk_expr(1) + "1",
                ]
        exprs = []
        for expr_k in exprs_k:
            for expr_pipi in exprs_pipi:
                for expr_op in exprs_ops:
                    exprs.append(expr_pipi * expr_op * expr_k)
        diagram_type_dict = dict()
        diagram_type_dict[((('x', 'x'), 1), (('x', 'x1_1'), 1), (('x1_1', 'x1_2'), 1), (('x1_2', 'x2'), 1), (('x2', 'x'), 1))] = "Type3"
        diagram_type_dict[((('x', 'x'), 1), (('x', 'x2'), 1), (('x1_1', 'x1_2'), 1), (('x1_2', 'x1_1'), 1), (('x2', 'x'), 1))] = "Type4"
        diagram_type_dict[((('x', 'x1_1'), 1), (('x', 'x1_2'), 1), (('x1_1', 'x'), 1), (('x1_2', 'x2'), 1), (('x2', 'x'), 1))] = "Type1"
        diagram_type_dict[((('x', 'x1_1'), 1), (('x', 'x2'), 1), (('x1_1', 'x1_2'), 1), (('x1_2', 'x'), 1), (('x2', 'x'), 1))] = "Type2"
        diagram_type_dict[((('x', 'x1_1'), 1), (('x1_1', 'x1_2'), 1), (('x1_2', 'x2'), 1), (('x2', 'x'), 1))] = "Type3"
        diagram_type_dict[((('x', 'x2'), 1), (('x1_1', 'x1_2'), 1), (('x1_2', 'x1_1'), 1), (('x2', 'x'), 1))] = "Type4"
        diagram_type_dict[((('x', 'x'), 1), (('x', 'x1_1'), 1), (('x1_1', 'x2'), 1), (('x2', 'x'), 1))] = "Type3"
        diagram_type_dict[((('x', 'x'), 1), (('x', 'x2'), 1), (('x1_1', 'x1_1'), 1), (('x2', 'x'), 1))] = "Type4"
        diagram_type_dict[((('x', 'x1_1'), 1), (('x', 'x2'), 1), (('x1_1', 'x'), 1), (('x2', 'x'), 1))] = "Type2"
        diagram_type_dict[((('x', 'x1_1'), 1), (('x1_1', 'x2'), 1), (('x2', 'x'), 1))] = "Type3"
        diagram_type_dict[((('x', 'x2'), 1), (('x1_1', 'x1_1'), 1), (('x2', 'x'), 1))] = "Type4"
        diagram_type_dict[((('x', 'x'), 1), (('x', 'x1_2'), 1), (('x1_2', 'x2'), 1), (('x2', 'x'), 1))] = "Type3"
        diagram_type_dict[((('x', 'x'), 1), (('x', 'x2'), 1), (('x1_2', 'x1_2'), 1), (('x2', 'x'), 1))] = "Type4"
        diagram_type_dict[((('x', 'x1_2'), 1), (('x', 'x2'), 1), (('x1_2', 'x'), 1), (('x2', 'x'), 1))] = "Type2"
        diagram_type_dict[((('x', 'x1_2'), 1), (('x1_2', 'x2'), 1), (('x2', 'x'), 1))] = "Type3"
        diagram_type_dict[((('x', 'x2'), 1), (('x1_2', 'x1_2'), 1), (('x2', 'x'), 1))] = "Type4"
        diagram_type_dict[((('x', 'x'), 1), (('x', 'x2'), 1), (('x2', 'x'), 1))] = "Type4"
        diagram_type_dict[((('x', 'x2'), 1), (('x2', 'x'), 1))] = "Type4"
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True, diagram_type_dict = diagram_type_dict)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.collect_op()
        return cexpr
    cexpr = q.pickle_cache_call(calc_cexpr, f"cache/auto_contractor_cexpr/kpipi-cexpr.{vol}.pickle")
    q.displayln_info(display_cexpr_raw(cexpr))
    return cexpr

@q.timer
def auto_contractor_kpipi_corr(job_tag, traj, get_prop, num_trials):
    total_site = ru.get_total_site(job_tag)
    vol = total_site[0] * total_site[1] * total_site[2]
    cexpr = get_cexpr_kpipi(vol)
    rng_state = q.RngState("seed")
    def positions_dict_maker(idx):
        rs = rng_state.split(str(idx))
        t1_1 = 5
        t1_2 = 7
        t = 3
        x1_1 = rs.c_rand_gen(total_site)
        x1_2 = rs.c_rand_gen(total_site)
        x = rs.c_rand_gen(total_site)
        x2 = rs.c_rand_gen(total_site)
        x1_1[3] = (x2[3] + t1_1) % total_site[3]
        x1_2[3] = (x2[3] + t1_2) % total_site[3]
        x[3] = (x2[3] + t) % total_site[3]
        pd = {
                "x1_1" : ("point", x1_1,),
                "x1_2" : ("point", x1_2,),
                "x" : ("point", x,),
                "x2" : ("point", x2,),
                }
        lmom = [ 2 * math.pi / total_site[i] for i in range(3) ]
        fac1 = sum([ 2 * cmath.rect(1.0, (x1_1[i] - x1_2[i]) * lmom[i]).real for i in range(3) ]) / np.sqrt(6)
        facs = [1.0, fac1]
        return pd, facs
    names_fac = [ "rest", "moving", ]
    trial_indices = get_mpi_chunk(range(num_trials))
    results_list = eval_cexpr_simulation(cexpr, positions_dict_maker = positions_dict_maker, trial_indices = trial_indices, get_prop = get_prop, is_only_total = "typed_total")
    q.qremove_all_info("analysis/kpipi")
    def mk_fn(info):
        def f(c):
            if c in "()<>/* ":
                return "_"
            else:
                return c
        info = "".join(map(f, info))
        while True:
            fn = info.replace("__", "_")
            if fn == info:
                break
            info = fn
        if fn[-1] == "_":
            fn = fn[:-1]
        return fn
    for name_fac, results in zip(names_fac, results_list):
        q.displayln_info(f"{name_fac} :")
        for k, v in results.items():
            q.displayln_info(f"{name_fac} {k}:\n  {v}")
            fn = "analysis/kpipi/" + mk_fn(f"{name_fac} {k}") + ".txt"
            [ a, e, ] = v
            q.qtouch(fn, f"0 {a.real} {a.imag} {e.real} {e.imag}\n")

@q.timer
def auto_contractor_kpipi_corr_81oprs(job_tag, traj, get_prop, num_trials):
    total_site = ru.get_total_site(job_tag)
    vol = total_site[0] * total_site[1] * total_site[2]
    exprs_odd_ops = [
            vol * mk_Q0_b81("x", "odd") + "Q0_b81(o)",
            vol * mk_Q1_b81("x", "odd") + "Q1_b81(o)",
            vol * mk_Q2_b81("x", "odd") + "Q2_b81(o)",
            vol * mk_Q3_b81("x", "odd") + "Q3_b81(o)",
            vol * mk_Q4_b81("x", "odd") + "Q4_b81(o)",
            vol * mk_Q5_b81("x", "odd") + "Q5_b81(o)",
            vol * mk_Q6_b81("x", "odd") + "Q6_b81(o)",
            vol * mk_Q7_b81("x", "odd") + "Q7_b81(o)",
            vol * mk_Q8_b81("x", "odd") + "Q8_b81(o)",
            ]
    exprs_even_ops = [
            vol * mk_Q0_b81("x", "even") + "Q0_b81(e)",
            vol * mk_Q1_b81("x", "even") + "Q1_b81(e)",
            vol * mk_Q2_b81("x", "even") + "Q2_b81(e)",
            vol * mk_Q3_b81("x", "even") + "Q3_b81(e)",
            vol * mk_Q4_b81("x", "even") + "Q4_b81(e)",
            vol * mk_Q5_b81("x", "even") + "Q5_b81(e)",
            vol * mk_Q6_b81("x", "even") + "Q6_b81(e)",
            vol * mk_Q7_b81("x", "even") + "Q7_b81(e)",
            vol * mk_Q8_b81("x", "even") + "Q8_b81(e)",
            ]
    exprs_ops = exprs_odd_ops + exprs_even_ops
    exprs_k = [
            vol * mk_k_0("x2") + "K0",
            ]
    exprs_pipi = [
            #vol**2 * mk_pipi_i0("x1_1", "x1_2", True) + "pipi_I0",
            vol**2 * mk_pipi_i20("x1_1", "x1_2", True) + "pipi_I2",
            #vol * mk_sigma("x1_1", True) + "sigma_1",
            #vol * mk_sigma("x1_2", True) + "sigma_2",
            #vol * mk_pi_0("x1_1", True) + "pi0_1",
            #vol * mk_pi_0("x1_2", True) + "pi0_2",
            #mk_expr(1) + "1",
            ]
    exprs = []
    for expr_k in exprs_k:
        for expr_pipi in exprs_pipi:
            for expr_op in exprs_ops:
                exprs.append(expr_pipi * expr_op * expr_k)
    diagram_type_dict = dict()
    cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True, diagram_type_dict = diagram_type_dict)
    q.displayln_info(display_cexpr(cexpr))
    cexpr.collect_op()
    q.displayln_info(display_cexpr_raw(cexpr))
    rng_state = q.RngState("seed")
    def positions_dict_maker(idx):
        rs = rng_state.split(str(idx))
        t1_1 = 5
        t1_2 = 7
        t = 3
        x1_1 = rs.c_rand_gen(total_site)
        x1_2 = rs.c_rand_gen(total_site)
        x = rs.c_rand_gen(total_site)
        x2 = rs.c_rand_gen(total_site)
        x1_1[3] = (x2[3] + t1_1) % total_site[3]
        x1_2[3] = (x2[3] + t1_2) % total_site[3]
        x[3] = (x2[3] + t) % total_site[3]
        pd = {
                "x1_1" : ("point", x1_1,),
                "x1_2" : ("point", x1_2,),
                "x" : ("point", x,),
                "x2" : ("point", x2,),
                }
        lmom = [ 2 * math.pi / total_site[i] for i in range(3) ]
        fac1 = sum([ 2 * cmath.rect(1.0, (x1_1[i] - x1_2[i]) * lmom[i]).real for i in range(3) ]) / np.sqrt(6)
        facs = [1.0, fac1]
        return pd, facs
    names_fac = [ "rest", "moving", ]
    trial_indices = get_mpi_chunk(range(num_trials))
    results_list = eval_cexpr_simulation(cexpr, positions_dict_maker = positions_dict_maker, trial_indices = trial_indices, get_prop = get_prop, is_only_total = "typed_total")
    q.qmkdir_info("analysis")
    q.qremove_all_info("analysis/kpipi_b81")
    q.qmkdir_info("analysis/kpipi_b81")
    def mk_fn(info):
        def f(c):
            if c in "()<>/* ":
                return "_"
            else:
                return c
        info = "".join(map(f, info))
        while True:
            fn = info.replace("__", "_")
            if fn == info:
                break
            info = fn
        if fn[-1] == "_":
            fn = fn[:-1]
        return fn
    for name_fac, results in zip(names_fac, results_list):
        q.displayln_info(f"{name_fac} :")
        for k, v in results.items():
            q.displayln_info(f"{name_fac} {k}:\n  {v}")
            fn = "analysis/kpipi_b81/" + mk_fn(f"{name_fac} {k}") + ".txt"
            [ a, e, ] = v
            q.qtouch(fn, f"0 {a.real} {a.imag} {e.real} {e.imag}\n")

@q.timer
def get_cexpr_3f4f(vol):
    @q.timer
    def calc_cexpr():
        total_site = ru.get_total_site(job_tag)
        vol = total_site[0] * total_site[1] * total_site[2]
        exprs_odd_ops = [
            vol * mk_Q0_b81("x", "odd") + "Q0_b81(o)",
            vol * mk_Q1_b81("x", "odd") + "Q1_b81(o)",
            vol * mk_Q2_b81("x", "odd") + "Q2_b81(o)",
            vol * mk_Q3_b81("x", "odd") + "Q3_b81(o)",
            vol * mk_Q4_b81("x", "odd") + "Q4_b81(o)",
            vol * mk_Q5_b81("x", "odd") + "Q5_b81(o)",
            vol * mk_Q6_b81("x", "odd") + "Q6_b81(o)",
            vol * mk_Q7_b81("x", "odd") + "Q7_b81(o)",
            vol * mk_Q8_b81("x", "odd") + "Q8_b81(o)",
        ]
        exprs_even_ops = [
            vol * mk_Q0_b81("x", "even") + "Q0_b81(e)",
            vol * mk_Q1_b81("x", "even") + "Q1_b81(e)",
            vol * mk_Q2_b81("x", "even") + "Q2_b81(e)",
            vol * mk_Q3_b81("x", "even") + "Q3_b81(e)",
            vol * mk_Q4_b81("x", "even") + "Q4_b81(e)",
            vol * mk_Q5_b81("x", "even") + "Q5_b81(e)",
            vol * mk_Q6_b81("x", "even") + "Q6_b81(e)",
            vol * mk_Q7_b81("x", "even") + "Q7_b81(e)",
            vol * mk_Q8_b81("x", "even") + "Q8_b81(e)",
        ]
        exprs_ops = exprs_odd_ops + exprs_even_ops
        exprs_src = [
            vol * mk_k_0("t2_1") + "K0",
            vol * mk_kpi_0_i1half("t2_1", "t2_2") + "Kpi_0_I1half",
            vol * mk_kpi_0_i3halves("t2_1", "t2_2") + "Kpi_0_I3halves",
        ]
        exprs_snk = [
            vol**2 * mk_pipi_i20("t1_1", "t1_2", True) + "pipi_I2",
            vol**2 * mk_pipi_i0("t1_1", "t1_2", True) + "pipi_I0",
            vol * mk_sigma("t1_1", True) + "sigma",
            vol * mk_pi_0("t1_1", True) + "pi0",
            mk_expr(1) + "1",
        ]
        #    exprs_src_snk = [
        #        vol**2 * mk_
        #        ]
        exprs = []
        for expr_src in exprs_src:
            for expr_snk in exprs_snk:
                for expr_op in exprs_ops:
                    exprs.append(expr_snk * expr_op * expr_src)
        diagram_type_dict = dict()
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True, diagram_type_dict = diagram_type_dict)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.collect_op()
        return cexpr
    cexpr = q.pickle_cache_call(calc_cexpr, f"cache/auto_contractor_cexpr/3f4f-cexpr.{vol}.pickle")
    q.displayln_info(display_cexpr(cexpr))
    return cexpr

@q.timer
def auto_contractor_3f4f_matching(job_tag, traj, get_prop, num_trials):
    total_site = ru.get_total_site(job_tag)
    vol = total_site[0] * total_site[1] * total_site[2]
    cexpr = get_cexpr_3f4f(vol)
    rng_state = q.RngState("seed")
    src_snk_seps = [8,10,12,14,16]
    tsep_src = -1
    tsep_snk = 1
    q.qmkdir_info("analysis")
    q.qremove_all_info("analysis/3f4f_b81")
    q.qmkdir_info("analysis/3f4f_b81")
    def mk_key(info):
        def f(c):
            if c in "()<>/* ":
                return "_"
            else:
                return c
        info = "".join(map(f, info))
        while True:
            fn = info.replace("__", "_")
            if fn == info:
                break
            info = fn
        if fn[-1] == "_":
            fn = fn[:-1]
        return fn
    #
    for tsnk_tsrc in src_snk_seps:
        max_top_tsrc = int(tsnk_tsrc/2)
        min_top_tsrc = int(tsnk_tsrc/2)
        fn = f"analysis/3f4f_b81/tsnk_tsrc{tsnk_tsrc}_top_tsrc{min_top_tsrc}-{max_top_tsrc}.bin"
        with open(fn, mode='wb') as f:
            for top_tsrc in range(min_top_tsrc,max_top_tsrc+1):
                def positions_dict_maker(idx):
                    rs = rng_state.split(str(idx))
                    tsrc1_top = - top_tsrc
                    tsrc2_top = tsep_src  + tsrc1_top
                    tsnk1_top = tsnk_tsrc + tsrc1_top
                    tsnk2_top = tsep_snk  + tsnk1_top
                    dt1_2 = 4
                    x = rs.c_rand_gen(total_site)
                    t2_1 = ( tsrc1_top + x[3] + total_site[3] ) % total_site[3]
                    t2_2 = ( tsrc2_top + x[3] + total_site[3] ) % total_site[3]
                    t1_1 = ( tsnk1_top + x[3] + total_site[3] ) % total_site[3]
                    t1_2 = ( tsnk2_top + x[3] + total_site[3] ) % total_site[3]
                    pd = {
                        "t1_1" : ("wall", t1_1,),
                        "t1_2" : ("wall", t1_2,),
                        "x" : ("point", x,),
                        "t2_1" : ("wall", t2_1,),
                        "t2_2" : ("wall", t2_2,),
                    }
                    facs = [1.0]
                    return pd, facs
                trial_indices = get_mpi_chunk(range(num_trials))
                results_list = eval_cexpr_simulation(cexpr, positions_dict_maker = positions_dict_maker, trial_indices = trial_indices, get_prop = get_prop, is_only_total = "total")
                #for results in results_list:
                results = results_list[0]
                for k, v in results.items():
                    if v[1].real == 0:
                        ratio_real = None
                    else:
                        ratio_real = v[0].real / v[1].real
                    if v[1].imag == 0:
                        ratio_imag = None
                    else:
                        ratio_imag = v[0].imag / v[1].imag
                    q.displayln_info(f"{k}:\n  {v}, ({ratio_real}, {ratio_imag})")
                    ###
                    [ a, e, ] = v
                    f.write(a)
                    f.write(e)
                if top_tsrc == min_top_tsrc and tsnk_tsrc == src_snk_seps[0]:
                    metafn = f"analysis/3f4f_b81/meta.txt"
                    with open(metafn, mode='w') as metaf:
                        for k, v in results.items():
                            key = mk_key(f"{k}")
                            metaf.write(f"{key}\n")

@q.timer_verbose
def run_gf(job_tag, traj):
    path_gf = get_load_path(f"configs/{job_tag}/ckpoint_lat.{traj}")
    if path_gf is None:
        if job_tag[:5] == "test-":
            gf = ru.mk_sample_gauge_field(job_tag, f"{traj}")
            path_gf = get_save_path(f"configs/{job_tag}/ckpoint_lat.{traj}")
            # gf.save(path_gf)
            qg.save_gauge_field(gf, path_gf)
        else:
            assert False
    get_gf = ru.load_config_lazy(job_tag, path_gf)
    return get_gf

@q.timer_verbose
def run_gt(job_tag, traj, get_gf):
    if None in [ get_gf, ]:
        return None
    path_gt = get_load_path(f"gauge-transform/{job_tag}/traj={traj}.field")
    if path_gt is None:
        if q.obtain_lock(f"locks/{job_tag}-{traj}-gauge_fix_coulomb"):
            gf = get_gf()
            gt = qg.gauge_fix_coulomb(gf)
            gt.save_double(get_save_path(f"gauge-transform/{job_tag}/traj={traj}.field"))
            q.release_lock()
            return lambda : gt
        else:
            return None
    else:
        @q.timer_verbose
        def load_gt():
            gt = q.GaugeTransform()
            gt.load_double(path_gt)
            # ADJUST ME
            # qg.check_gauge_fix_coulomb(get_gf(), gt)
            #
            return gt
        get_gt = q.lazy_call(load_gt)
    return get_gt

@q.timer
def compute_eig(gf, job_tag, inv_type = 0, inv_acc = 0, *, path = None):
    # return a function ``get_eig''
    # ``get_eig()'' return the ``eig''
    load_eig = ru.load_eig_lazy(get_load_path(path), job_tag)
    if load_eig is not None:
        return load_eig
    # evec, evals = ru.mk_eig(gf, job_tag, inv_type, inv_acc)
    basis, cevec, smoothed_evals = ru.mk_ceig(gf, job_tag, inv_type, inv_acc)
    eig = [ basis, cevec, smoothed_evals ]
    ru.save_ceig(get_save_path(path + ".partial"), eig, job_tag, inv_type, inv_acc);
    q.qrename_info(get_save_path(path + ".partial"), get_save_path(path))
    test_eig(gf, eig, job_tag, inv_type)
    def get_eig():
        return eig
    return get_eig

@q.timer
def test_eig(gf, eig, job_tag, inv_type):
    geo = gf.geo()
    src = q.FermionField4d(geo)
    q.displayln_info(f"src norm {src.qnorm()}")
    src.set_rand(q.RngState("test_eig:{id(inv)}"))
    sol_ref = ru.get_inv(gf, job_tag, inv_type, inv_acc = 2, eig = eig, eps = 1e-10, mpi_split = False, timer = False) * src
    q.displayln_info(f"sol_ref norm {sol_ref.qnorm()} with eig")
    for inv_acc in [0, 1, 2]:
        sol = ru.get_inv(gf, job_tag, inv_type, inv_acc, eig = eig, mpi_split = False, timer = False) * src
        sol -= sol_ref
        q.displayln_info(f"sol diff norm {sol.qnorm()} inv_acc={inv_acc} with eig")
        sol = ru.get_inv(gf, job_tag, inv_type, inv_acc, mpi_split = False, timer = False) * src
        sol -= sol_ref
        q.displayln_info(f"sol diff norm {sol.qnorm()} inv_acc={inv_acc} without eig")

@q.timer_verbose
def run_eig(job_tag, traj, get_gf):
    if None in [ get_gf, ]:
        return None
    get_eig = ru.load_eig_lazy(get_load_path(f"eig/{job_tag}/traj={traj}"), job_tag)
    if get_eig is None and get_gf is not None:
        if q.obtain_lock(f"locks/{job_tag}-{traj}-run-eig"):
            q.qmkdir_info(get_save_path(f"eig"))
            q.qmkdir_info(get_save_path(f"eig/{job_tag}"))
            get_eig = compute_eig(get_gf(), job_tag, inv_type = 0, path = f"eig/{job_tag}/traj={traj}")
            q.release_lock()
    return get_eig

@q.timer_verbose
def run_prop(job_tag, traj, get_gf, get_gt, get_eig):
    for inv_type in [0, 1, 2,]:
        if get_load_path(f"prop-wsrc-{inv_type}/{job_tag}/traj={traj}") is None:
            if inv_type == 0 and get_eig is None:
                continue
            if q.obtain_lock(f"locks/{job_tag}-{traj}-compute-prop-wsrc-all-{inv_type}"):
                q.qmkdir_info(get_save_path(f"prop-wsrc-{inv_type}"))
                q.qmkdir_info(get_save_path(f"prop-wsrc-{inv_type}/{job_tag}"))
                if inv_type == 0:
                    eig = get_eig()
                else:
                    eig = None
                compute_prop_wsrc_all(get_gf(), get_gt(), job_tag, inv_type,
                        path_s = f"prop-wsrc-{inv_type}/{job_tag}/traj={traj}",
                        path_sp = f"psel-prop-wsrc-{inv_type}/{job_tag}/traj={traj}",
                        eig = eig)
                q.release_lock()
        #
        if get_load_path(f"prop-psrc-{inv_type}/{job_tag}/traj={traj}") is None:
            if inv_type == 0 and get_eig is None:
                continue
            if q.obtain_lock(f"locks/{job_tag}-{traj}-compute-prop-psrc-all-{inv_type}"):
                q.qmkdir_info(get_save_path(f"prop-psrc-{inv_type}"))
                q.qmkdir_info(get_save_path(f"prop-psrc-{inv_type}/{job_tag}"))
                if inv_type == 0:
                    eig = get_eig()
                else:
                    eig = None
                compute_prop_psrc_all(get_gf(), job_tag, inv_type,
                        path_s = f"prop-psrc-{inv_type}/{job_tag}/traj={traj}",
                        eig = eig)
                q.release_lock()
    #
    path_prop_list = \
            [ get_load_path(f"prop-psrc-{inv_type}/{job_tag}/traj={traj}") for inv_type in [0, 1, 2,] ] + \
            [ get_load_path(f"prop-wsrc-{inv_type}/{job_tag}/traj={traj}") for inv_type in [0, 1, 2,] ]
    if all(map(lambda x : x is not None, path_prop_list)):
        total_site = ru.get_total_site(job_tag)
        psel = q.PointSelection(get_all_points(total_site))
        psel_pos_dict = dict([ (tuple(pos), i) for i, pos in enumerate(psel.to_list()) ])
        load_prop_wsrc_all(job_tag, traj, "l",
                f"prop-wsrc-0/{job_tag}/traj={traj}",
                f"psel-prop-wsrc-0/{job_tag}/traj={traj}",
                psel, get_gt())
        load_prop_wsrc_all(job_tag, traj, "s",
                f"prop-wsrc-1/{job_tag}/traj={traj}",
                f"psel-prop-wsrc-1/{job_tag}/traj={traj}",
                psel, get_gt())
        load_prop_wsrc_all(job_tag, traj, "c",
                f"prop-wsrc-2/{job_tag}/traj={traj}",
                f"psel-prop-wsrc-2/{job_tag}/traj={traj}",
                psel, get_gt())
        load_prop_psrc_all(job_tag, traj, "l", f"prop-psrc-0/{job_tag}/traj={traj}", psel)
        load_prop_psrc_all(job_tag, traj, "s", f"prop-psrc-1/{job_tag}/traj={traj}", psel)
        load_prop_psrc_all(job_tag, traj, "c", f"prop-psrc-2/{job_tag}/traj={traj}", psel)
        prop_cache = q.mk_cache(f"prop_cache-{job_tag}-{traj}")
        get_prop = mk_get_prop(prop_cache, psel_pos_dict = psel_pos_dict)
        return get_prop
    return None

@q.timer
def run_job(job_tag, traj):
    q.check_stop()
    q.check_time_limit()
    #
    get_gf = run_gf(job_tag, traj)
    get_gt = run_gt(job_tag, traj, get_gf)
    get_eig = run_eig(job_tag, traj, get_gf)
    #
    get_prop = run_prop(job_tag, traj, get_gf, get_gt, get_eig)
    #
    if get_prop is not None:
        # auto_contractor_simple_test(job_tag, traj, get_prop)
        num_trials = 100
        # auto_contractor_vev(job_tag, traj, get_prop, num_trials)
        # auto_contractor_test_corr(job_tag, traj, get_prop, num_trials)
        # auto_contractor_meson_corr(job_tag, traj, get_prop, num_trials)
        # auto_contractor_pipi_corr(job_tag, traj, get_prop, num_trials)
        # auto_contractor_kpipi_corr(job_tag, traj, get_prop, num_trials)
        # auto_contractor_kpipi_corr_81oprs(job_tag, traj, get_prop, num_trials)
        auto_contractor_3f4f_matching(job_tag, traj, get_prop, num_trials)
    #
    q.clean_cache()

if __name__ == "__main__":
    qg.begin_with_gpt()
    job_tag = "test-4nt16"
    q.displayln_info(pprint.pformat(rup.dict_params[job_tag]))
    traj = 1000
    run_job(job_tag, traj)
    q.timer_display()
    qg.end_with_gpt()
