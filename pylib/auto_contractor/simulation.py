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
def mk_sample_gauge_field(job_tag, fn):
    rs = q.RngState(f"seed {job_tag} {fn}").split("mk_sample_gauge_field")
    total_site = ru.get_total_site(job_tag)
    geo = q.Geometry(total_site, 1)
    gf = q.GaugeField(geo)
    gf.set_rand(rs, sigma = 0.25, n_step = 12)
    for i in range(7):
        q.gf_wilson_flow_step(gf, 0.05)
    gf.unitarize()
    return gf

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
    ru.save_ceig(get_save_path(path), eig, job_tag, inv_type, inv_acc);
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
    for inv_acc in [0, 1,]:
        sol = ru.get_inv(gf, job_tag, inv_type, inv_acc, eig = eig, mpi_split = False, timer = False) * src
        sol -= sol_ref
        q.displayln_info(f"sol diff norm {sol.qnorm()} inv_acc={inv_acc} with eig")
        sol = ru.get_inv(gf, job_tag, inv_type, inv_acc, mpi_split = False, timer = False) * src
        sol -= sol_ref
        q.displayln_info(f"sol diff norm {sol.qnorm()} inv_acc={inv_acc} without eig")

@q.timer
def compute_prop(inv, src, *, tag, sfw):
    sol = inv * src
    sol.save_double(sfw, tag)

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
    finished_tags = q.properly_truncate_fields_sync_node(get_save_path(path_s + ".acc"))
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

def adj_msc(x):
    x = g.adj(x)
    return g.tensor(np.ascontiguousarray(x.array), x.otype)

@q.timer
def load_prop_psrc_all(job_tag, traj, flavor : str, path_s : str):
    cache = q.mk_cache(f"prop_cache-{job_tag}-{traj}", flavor)
    total_site = ru.get_total_site(job_tag)
    get_all_points(total_site)
    if flavor in ["l", "u", "d",]:
        inv_type = 0
    elif flavor in ["s",]:
        inv_type = 1
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
        prop_msc = q.Prop()
        q.convert_mspincolor_from_wm_prop(prop_msc, prop)
        cache[f"xg=({xg[0]},{xg[1]},{xg[2]},{xg[3]})"] = prop_msc
    sfr.close()

@q.timer
def auto_contractor_simple_test(job_tag, traj):
    prop_cache = q.mk_cache(f"prop_cache-{job_tag}-{traj}")
    q.displayln_info(g.gamma[5] * get_prop_psnk_psrc(prop_cache, "l", [1, 2, 3, 2], [1, 2, 3, 2]))
    q.displayln_info(g.trace(g.gamma[5] * get_prop_psnk_psrc(prop_cache, "l", [1, 2, 3, 2], [1, 2, 3, 2])))
    q.displayln_info(g.gamma[5] * get_prop_psnk_psrc(prop_cache, "l", [1, 2, 3, 2], [1, 2, 3, 2]) - get_prop_psnk_psrc(prop_cache, "l", [1, 2, 3, 2], [1, 2, 3, 2]))
    q.displayln_info(g.gamma[5] * adj_msc(get_prop_psnk_psrc(prop_cache, "l", [1, 2, 3, 2], [1, 2, 3, 2])))
    q.displayln_info(g.norm2(g.gamma[5] * adj_msc(get_prop_psnk_psrc(prop_cache, "l", [1, 2, 3, 2], [1, 2, 3, 4])) * g.gamma[5] - get_prop_psnk_psrc(prop_cache, "l", [1, 2, 3, 4], [1, 2, 3, 2])))
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
    val = eval_cexpr(cexpr, positions_dict = positions_dict, prop_cache = prop_cache)
    q.displayln_info("eval_cexpr: ", val)
    q.displayln_info("gpt_direct: ",
            -g.trace(
                g.gamma[5]
                * get_prop_psnk_psrc(prop_cache, "l", [1, 2, 3, 4], [1, 2, 3, 2])
                * g.gamma[5]
                * get_prop_psnk_psrc(prop_cache, "l", [1, 2, 3, 2], [1, 2, 3, 4])
                ))

@q.timer
def auto_contractor_meson_corr(job_tag, traj, num_trials):
    total_site = ru.get_total_site(job_tag)
    vol = total_site[0] * total_site[1] * total_site[2]
    exprs = [
            vol**2 * mk_pi_p("x2", True) * mk_pi_p("x1"),
            vol**2 * mk_k_p("x2", True) * mk_k_p("x1"),
            ]
    names = [
            "pi pi",
            "k k",
            ]
    cexpr = contract_simplify_round_compile(*exprs, is_isospin_symmetric_limit = True)
    q.displayln_info(display_cexpr(cexpr))
    cexpr.collect_op()
    q.displayln_info(display_cexpr(cexpr))
    def positions_dict_maker(idx, rs, total_site):
        t2 = 5
        x1 = rs.c_rand_gen(total_site)
        x2 = rs.c_rand_gen(total_site)
        x2[3] = (x1[3] + t2) % total_site[3]
        pd = {
                "x1" : x1,
                "x2" : x2,
                }
        lmom = [ 2 * math.pi / total_site[i] for i in range(3) ]
        facs = [
                1.0,
                sum([ cmath.rect(1.0, (x2[i] - x1[i]) * lmom[i]).real for i in range(3) ]) / 3.0,
                sum([ cmath.rect(1.0, (x2[i] - x1[i]) * 2 * lmom[i]).real for i in range(3) ]) / 3.0,
                ]
        return pd, facs
    names_fac = ["rest", "mom1", "mom2",]
    rng_state = q.RngState("seed")
    trial_indices = range(num_trials)
    prop_cache = q.mk_cache(f"prop_cache-{job_tag}-{traj}")
    results_list = eval_cexpr_simulation(cexpr, positions_dict_maker = positions_dict_maker, rng_state = rng_state, trial_indices = trial_indices, total_site = total_site, prop_cache = prop_cache, is_only_total = True)
    for name_fac, results in zip(names_fac, results_list):
        q.displayln_info(f"{name_fac} :")
        for n, (k, v) in zip(names, results.items()):
            q.displayln_info(f"{k:>10} {n:>40} : {v} ")

@q.timer
def auto_contractor_pipi_corr(job_tag, traj, num_trials):
    total_site = ru.get_total_site(job_tag)
    vol = total_site[0] * total_site[1] * total_site[2]
    exprs = [
            vol**2 * mk_pipi_i0("x21", "x22", True),
            vol**2 * mk_pipi_i0("x11", "x12"),
            vol * mk_sigma("x11"),
            vol * mk_sigma("x12"),
            vol * mk_sigma("x21", True),
            vol * mk_sigma("x22", True),
            vol**4 * mk_pipi_i0("x21", "x22", True) * mk_pipi_i0("x11", "x12"),
            vol**3 * mk_pipi_i0("x21", "x22", True) * mk_sigma("x11"),
            vol**3 * mk_pipi_i0("x21", "x22", True) * mk_sigma("x12"),
            vol**3 * mk_sigma("x21", True) * mk_pipi_i0("x11", "x12"),
            vol**3 * mk_sigma("x22", True) * mk_pipi_i0("x11", "x12"),
            vol**4 * mk_pipi_i11("x21", "x22", True) * mk_pipi_i11("x11", "x12"),
            vol**4 * mk_pipi_i22("x21", "x22", True) * mk_pipi_i22("x11", "x12"),
            ]
    names = [
            "pipi_i0_2",
            "pipi_i0_1",
            "sigma_11 ",
            "sigma_12 ",
            "sigma_21 ",
            "sigma_22 ",
            "pipi_i0_2   pipi_i0_1 ",
            "pipi_i0_2   sigma_11  ",
            "pipi_i0_2   sigma_12  ",
            "sigma_21    pipi_i0_1 ",
            "sigma_22    pipi_i0_1 ",
            "pipi_i11_2  pipi_i11_1",
            "pipi_i22_2  pipi_i22_1",
            ]
    cexpr = contract_simplify_round_compile(*exprs, is_isospin_symmetric_limit = True)
    q.displayln_info(display_cexpr(cexpr))
    cexpr.collect_op()
    q.displayln_info(display_cexpr(cexpr))
    def positions_dict_maker(idx, rs, total_site):
        t12 = 2
        t21 = 5
        t22 = 7
        x11 = rs.c_rand_gen(total_site)
        x12 = rs.c_rand_gen(total_site)
        x21 = rs.c_rand_gen(total_site)
        x22 = rs.c_rand_gen(total_site)
        x12[3] = (x11[3] + t12) % total_site[3]
        x21[3] = (x11[3] + t21) % total_site[3]
        x22[3] = (x11[3] + t22) % total_site[3]
        pd = {
                "x11" : x11,
                "x12" : x12,
                "x21" : x21,
                "x22" : x22,
                }
        lmom = [ 2 * math.pi / total_site[i] for i in range(3) ]
        fac1 = sum([ cmath.rect(1.0, (x11[i] - x12[i]) * lmom[i]).real for i in range(3) ])
        fac2 = sum([ cmath.rect(1.0, (x21[i] - x22[i]) * lmom[i]).real for i in range(3) ])
        facs = [1.0, fac1, fac2, fac1 * fac2,]
        return pd, facs
    names_fac = ["rest-rest", "rest-moving", "moving-rest", "moving-moving",]
    rng_state = q.RngState("seed")
    trial_indices = range(num_trials)
    total_site = ru.get_total_site(job_tag)
    prop_cache = q.mk_cache(f"prop_cache-{job_tag}-{traj}")
    results_list = eval_cexpr_simulation(cexpr, positions_dict_maker = positions_dict_maker, rng_state = rng_state, trial_indices = trial_indices, total_site = total_site, prop_cache = prop_cache, is_only_total = True)
    for name_fac, results in zip(names_fac, results_list):
        q.displayln_info(f"{name_fac} :")
        for n, (k, v) in zip(names, results.items()):
            q.displayln_info(f"{k:>10} {n:>40} : {v} ")

@q.timer
def auto_contractor_kpipi_corr(job_tag, traj, num_trials):
    total_site = ru.get_total_site(job_tag)
    vol = total_site[0] * total_site[1] * total_site[2]
    exprs_odd_ops = [
            vol * mk_Q1("x", "odd"),
            vol * mk_Q2("x", "odd"),
            vol * mk_Q3("x", "odd"),
            vol * mk_Q4("x", "odd"),
            vol * mk_Q5("x", "odd"),
            vol * mk_Q6("x", "odd"),
            vol * mk_Q7("x", "odd"),
            vol * mk_Q8("x", "odd"),
            vol * mk_Q9("x", "odd"),
            vol * mk_Q10("x", "odd"),
            vol * mk_Qsub("x", "odd"),
            ]
    names_odd_ops = [ f"Q{i+1}(o)" for i in range(10) ] + ["Qs(o)",]
    exprs_even_ops = [
            vol * mk_Q1("x", "even"),
            vol * mk_Q2("x", "even"),
            vol * mk_Q3("x", "even"),
            vol * mk_Q4("x", "even"),
            vol * mk_Q5("x", "even"),
            vol * mk_Q6("x", "even"),
            vol * mk_Q7("x", "even"),
            vol * mk_Q8("x", "even"),
            vol * mk_Q9("x", "even"),
            vol * mk_Q10("x", "even"),
            vol * mk_Qsub("x", "even"),
            ]
    names_even_ops = [ f"Q{i+1}(e)" for i in range(10) ] + ["Qs(e)",]
    exprs_ops = exprs_odd_ops + exprs_even_ops
    names_ops = names_odd_ops + names_even_ops
    exprs_k = [
            vol * mk_k_0("x2"),
            ]
    names_k = [ "k" ]
    exprs_pipi = [
            vol**2 * mk_pipi_i0("x11", "x12", True),
            vol**2 * mk_pipi_i20("x11", "x12", True),
            vol * mk_sigma("x11", True),
            vol * mk_sigma("x12", True),
            vol * mk_pi_0("x11", True),
            vol * mk_pi_0("x12", True),
            1,
            ]
    names_pipi= [ "pipi_i0 ", "pipi_i20", "sigma_11", "sigma_12", "pi_11", "pi_12", "",]
    exprs = []
    names = []
    for name_k, expr_k in zip(names_k, exprs_k):
        for name_pipi, expr_pipi in zip(names_pipi, exprs_pipi):
            for name_op, expr_op in zip(names_ops, exprs_ops):
                names.append(f"{name_pipi} {name_op:6} {name_k}")
                exprs.append(expr_pipi * expr_op * expr_k)
    cexpr = contract_simplify_round_compile(*exprs, is_isospin_symmetric_limit = True)
    q.displayln_info(display_cexpr(cexpr))
    cexpr.collect_op()
    q.displayln_info(display_cexpr(cexpr))
    def positions_dict_maker(idx, rs, total_site):
        t11 = 5
        t12 = 7
        t = 3
        x11 = rs.c_rand_gen(total_site)
        x12 = rs.c_rand_gen(total_site)
        x = rs.c_rand_gen(total_site)
        x2 = rs.c_rand_gen(total_site)
        x11[3] = (x2[3] + t11) % total_site[3]
        x12[3] = (x2[3] + t12) % total_site[3]
        x[3] = (x2[3] + t) % total_site[3]
        pd = {
                "x11" : x11,
                "x12" : x12,
                "x" : x,
                "x2" : x2,
                }
        lmom = [ 2 * math.pi / total_site[i] for i in range(3) ]
        fac1 = sum([ cmath.rect(1.0, (x11[i] - x12[i]) * lmom[i]).real for i in range(3) ])
        facs = [1.0, fac1]
        return pd, facs
    names_fac = ["rest", "moving",]
    rng_state = q.RngState("seed")
    trial_indices = range(num_trials)
    prop_cache = q.mk_cache(f"prop_cache-{job_tag}-{traj}")
    results_list = eval_cexpr_simulation(cexpr, positions_dict_maker = positions_dict_maker, rng_state = rng_state, trial_indices = trial_indices, total_site = total_site, prop_cache = prop_cache, is_only_total = True)
    for name_fac, results in zip(names_fac, results_list):
        q.displayln_info(f"{name_fac} :")
        for n, (k, v) in zip(names, results.items()):
            q.displayln_info(f"{k:>10} {n:>40} : {v} ")

@q.timer
def run_job(job_tag, traj):
    q.check_stop()
    q.check_time_limit()
    q.displayln_info(rup.dict_params[job_tag])
    #
    q.qmkdir_info(f"locks")
    q.qmkdir_info(get_save_path(f""))
    #
    total_site = ru.get_total_site(job_tag)
    geo = q.Geometry(total_site, 1)
    q.displayln_info("geo.show() =", geo.show())
    #
    path_gf = get_load_path(f"configs/{job_tag}/ckpoint_lat.{traj}")
    if path_gf is None:
        if job_tag[:5] == "test-":
            gf = ru.mk_sample_gauge_field(job_tag, f"{traj}")
            q.qmkdir_info(get_save_path(f"configs"))
            q.qmkdir_info(get_save_path(f"configs/{job_tag}"))
            path_gf = get_save_path(f"configs/{job_tag}/ckpoint_lat.{traj}")
            gf.save(path_gf)
        else:
            assert False
    gf = ru.load_config(job_tag, path_gf)
    gf.show_info()
    #
    get_eig = None
    if q.obtain_lock(f"locks/{job_tag}-{traj}-compute-eig"):
        q.qmkdir_info(get_save_path(f"eig"))
        q.qmkdir_info(get_save_path(f"eig/{job_tag}"))
        get_eig = compute_eig(gf, job_tag, inv_type = 0, path = f"eig/{job_tag}/traj={traj}")
        q.release_lock()
    #
    for inv_type in [0, 1,]:
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
                compute_prop_psrc_all(gf, job_tag, inv_type, path_s = f"prop-psrc-{inv_type}/{job_tag}/traj={traj}", eig = eig)
                q.release_lock()
    #
    path_prop_list = [ get_load_path(f"prop-psrc-{inv_type}/{job_tag}/traj={traj}") for inv_type in [0, 1,] ]
    if all(map(lambda x : x is not None, path_prop_list)):
        load_prop_psrc_all(job_tag, traj, "l", f"prop-psrc-0/{job_tag}/traj={traj}")
        load_prop_psrc_all(job_tag, traj, "s", f"prop-psrc-1/{job_tag}/traj={traj}")
        # auto_contractor_simple_test(job_tag, traj)
        num_trials = 100
        auto_contractor_meson_corr(job_tag, traj, num_trials)
        auto_contractor_pipi_corr(job_tag, traj, num_trials)
        auto_contractor_kpipi_corr(job_tag, traj, num_trials)
    #
    q.clean_cache()

if __name__ == "__main__":
    qg.begin_with_gpt()
    job_tag = "test-4nt16"
    traj = 1000
    rup.dict_params[job_tag]["fermion_params"][0][2] = rup.dict_params[job_tag]["fermion_params"][0][0]
    rup.dict_params[job_tag]["fermion_params"][1][2] = rup.dict_params[job_tag]["fermion_params"][1][0]
    rup.dict_params[job_tag]["load_config_params"]["twist_boundary_at_boundary"] = [0.0, 0.0, 0.0, -0.5,]
    run_job(job_tag, traj)
    q.timer_display()
    qg.end_with_gpt()
