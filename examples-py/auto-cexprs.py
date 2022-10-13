#!/usr/bin/env python3

import qlat as q

from auto_contractor.operators import *
from auto_contractor.eval import *

import sys

@q.timer
def get_cexpr_zeros():
    def calc_cexpr():
        exprs = [
                mk_expr(0),
                mk_k_p("t_2", True)     * mk_k_0("t_1")     + "k+^dag  * k0    ",
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.optimize()
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, f"cache/auto_contract_cexpr/get_cexpr_zeros")

@q.timer
def get_cexpr_meson_corr():
    def calc_cexpr():
        exprs = [
                mk_pi_p("t_2", True)    * mk_pi_p("t_1")    + "pi^dag * pi   ",
                mk_k_p("t_2", True)     * mk_k_p("t_1")     + "k^dag  * k    ",
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.optimize()
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, f"cache/auto_contract_cexpr/get_cexpr_meson_corr")

@q.timer
def get_cexpr_meson_f_corr():
    def calc_cexpr():
        exprs = [
                mk_j5pi_mu("x_2", 3)    * mk_pi_p("t_1")    + "a_pi   * pi   ",
                mk_j5k_mu("x_2", 3)     * mk_k_p("t_1")     + "a_k    * k    ",
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.optimize()
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, f"cache/auto_contract_cexpr/get_cexpr_meson_f_corr")

def mk_bk_vv_aa(p : str):
    s = 0
    for mu in range(4):
        v1 = mk_vec_mu("s", "d", p, mu) - mk_vec5_mu("s", "d", p, mu)
        v2 = mk_vec_mu("s", "d", p, mu) - mk_vec5_mu("s", "d", p, mu)
        s = s + v1 * v2
    return s + f"Ok_{{VV+AA}}"

def mk_bpi_vv_aa(p : str):
    s = 0
    for mu in range(4):
        v1 = mk_vec_mu("u", "d", p, mu) - mk_vec5_mu("u", "d", p, mu)
        v2 = mk_vec_mu("u", "d", p, mu) - mk_vec5_mu("u", "d", p, mu)
        s = s + v1 * v2
    return s + f"Opi_{{VV+AA}}"

def mk_bkpi1_vv_aa(p : str):
    s = 0
    for mu in range(4):
        v1 = mk_vec_mu("d", "u", p, mu) - mk_vec5_mu("d", "u", p, mu)
        v2 = mk_vec_mu("u'", "s", p, mu) - mk_vec5_mu("u'", "s", p, mu)
        s = s + v1 * v2
    return s + f"Okpi_{{VV+AA}}"

def mk_bkpi2_vv_aa(p : str):
    s = 0
    for mu in range(4):
        v1 = mk_vec_mu("u'", "u", p, mu) - mk_vec5_mu("u'", "u", p, mu)
        v2 = mk_vec_mu("d", "s", p, mu) - mk_vec5_mu("d", "s", p, mu)
        s = s + v1 * v2
    return s + f"Okpi_{{VV+AA}}"

@q.timer
def get_cexpr_meson_bk_bpi_corr():
    def calc_cexpr():
        exprs = [
                mk_meson("d", "s", "t_2") * mk_bk_vv_aa("x") * mk_meson("d", "s", "t_1")
                + "(i dbar g5 s)  * (sbar  gmu (1-g5) d)(sbar  gmu (1-g5) d) * (i dbar g5 s)",
                mk_meson("d", "u", "t_2") * mk_bpi_vv_aa("x") * mk_meson("d", "u", "t_1")
                + "(i dbar g5 u)  * (ubar  gmu (1-g5) d)(ubar  gmu (1-g5) d) * (i dbar g5 u)",
                mk_meson("u", "u'", "t_2") * mk_bkpi1_vv_aa("x") * mk_meson("s", "d", "t_1")
                + "(i ubar g5 u') * (dbar  gmu (1-g5) u)(u'bar gmu (1-g5) s) * (i sbar g5 d)",
                mk_meson("u", "u'", "t_2") * mk_bkpi2_vv_aa("x") * mk_meson("s", "d", "t_1")
                + "(i ubar g5 u') * (u'bar gmu (1-g5) u)(dbar  gmu (1-g5) s) * (i sbar g5 d)",
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.optimize()
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, f"cache/auto_contract_cexpr/get_cexpr_meson_bk_bpi_corr")

@q.timer
def get_cexpr_meson_corr2():
    def calc_cexpr():
        exprs = [
                mk_pi_p("x_2", True)    * mk_pi_p("x_1")    + "pi     * pi   ",
                mk_k_p("x_2", True)     * mk_k_p("x_1")     + "k      * k    ",
                mk_a0_p("x_2", True)    * mk_a0_p("x_1")    + "a0     * a0   ",
                mk_kappa_p("x_2", True) * mk_kappa_p("x_1") + "kappa  * kappa",
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.optimize()
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, f"cache/auto_contract_cexpr/get_cexpr_meson_corr2")

@q.timer
def get_cexpr_meson_corr3():
    def calc_cexpr():
        exprs = [
                mk_pi_0("t_1", True) * mk_pi_0("t_2"),
                sympy.simplify(1)/2 * (
                    mk_k_0("t_1", True) * mk_k_0("t_2")
                    + mk_k_0_bar("t_1", True) * mk_k_0_bar("t_2")
                    ),
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.optimize()
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, f"cache/auto_contract_cexpr/get_cexpr_meson_corr3")

@q.timer
def get_cexpr_meson_f_corr2():
    def calc_cexpr():
        exprs = [
                mk_j5pi_mu("x_2", 3) * mk_pi_p("t_1") + "(a_pi * pi)",
                sympy.simplify(1)/2 * (
                    mk_j5k_mu("x_2", 3)  * mk_k_p("t_1")
                    + mk_j5km_mu("x_2", 3)  * mk_k_m("t_1")
                    ) + "(a_k  * k )",
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.optimize()
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, f"cache/auto_contract_cexpr/get_cexpr_meson_f_corr2")

@q.timer
def get_cexpr_meson_quark_mass():
    def calc_cexpr():
        exprs = [
                mk_pi_0("t_1", True) * mk_m("u", "x_1") * mk_pi_0("t_2"),
                mk_pi_0("t_1", True) * mk_m("d", "x_1") * mk_pi_0("t_2"),
                sympy.simplify(1)/2 * (
                    mk_pi_p("t_1", True) * mk_m("u", "x_1") * mk_pi_p("t_2")
                    + mk_pi_m("t_1", True) * mk_m("u", "x_1") * mk_pi_m("t_2")
                    ),
                sympy.simplify(1)/2 * (
                    mk_pi_p("t_1", True) * mk_m("d", "x_1") * mk_pi_p("t_2")
                    + mk_pi_m("t_1", True) * mk_m("d", "x_1") * mk_pi_m("t_2")
                    ),
                sympy.simplify(1)/2 * (
                    mk_k_0("t_1", True) * mk_m("d", "x_1") * mk_k_0("t_2")
                    + mk_k_0_bar("t_1", True) * mk_m("d", "x_1") * mk_k_0_bar("t_2")
                    ),
                sympy.simplify(1)/2 * (
                    mk_k_0("t_1", True) * mk_m("s", "x_1") * mk_k_0("t_2")
                    + mk_k_0_bar("t_1", True) * mk_m("s", "x_1") * mk_k_0_bar("t_2")
                    ),
                sympy.simplify(1)/2 * (
                    mk_k_p("t_1", True) * mk_m("u", "x_1") * mk_k_p("t_2")
                    + mk_k_m("t_1", True) * mk_m("u", "x_1") * mk_k_m("t_2")
                    ),
                sympy.simplify(1)/2 * (
                    mk_k_p("t_1", True) * mk_m("s", "x_1") * mk_k_p("t_2")
                    + mk_k_m("t_1", True) * mk_m("s", "x_1") * mk_k_m("t_2")
                    ),
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.optimize()
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, f"cache/auto_contract_cexpr/get_cexpr_meson_quark_mass")

@q.timer
def get_cexpr_meson_jt_zv():
    def calc_cexpr():
        exprs = [
                sympy.simplify(1)/2 * (
                    mk_pi_p("t_1", True) * mk_vec_mu("u", "u", "x_1", 3) * mk_pi_p("t_2")
                    - mk_pi_m("t_1", True) * mk_vec_mu("u", "u", "x_1", 3) * mk_pi_m("t_2")
                    ),
                sympy.simplify(1)/2 * (
                    mk_k_p("t_1", True) * mk_vec_mu("u", "u", "x_1", 3) * mk_k_p("t_2")
                    - mk_k_m("t_1", True) * mk_vec_mu("u", "u", "x_1", 3) * mk_k_m("t_2")
                    ),
                sympy.simplify(1)/2 * (
                    mk_k_m("t_1", True) * mk_vec_mu("s", "s", "x_1", 3) * mk_k_m("t_2")
                    - mk_k_p("t_1", True) * mk_vec_mu("s", "s", "x_1", 3) * mk_k_p("t_2")
                    ),
                sympy.simplify(1)/2 * (
                    mk_pi_p("t_1p", True) * mk_vec_mu("u", "u", "x_1", 3) * mk_pi_p("t_2p")
                    - mk_pi_m("t_1p", True) * mk_vec_mu("u", "u", "x_1", 3) * mk_pi_m("t_2p")
                    ),
                sympy.simplify(1)/2 * (
                    mk_k_p("t_1p", True) * mk_vec_mu("u", "u", "x_1", 3) * mk_k_p("t_2p")
                    - mk_k_m("t_1p", True) * mk_vec_mu("u", "u", "x_1", 3) * mk_k_m("t_2p")
                    ),
                sympy.simplify(1)/2 * (
                    mk_k_m("t_1p", True) * mk_vec_mu("s", "s", "x_1", 3) * mk_k_m("t_2p")
                    - mk_k_p("t_1p", True) * mk_vec_mu("s", "s", "x_1", 3) * mk_k_p("t_2p")
                    ),
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.optimize()
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, f"cache/auto_contract_cexpr/get_cexpr_meson_jt_zv")

@q.timer
def get_cexpr_meson_jj_mm():
    def calc_cexpr():
        exprs = [
                sum([
                    mk_pi_0("t_1", True) * mk_j_mu("x_1", mu) * mk_j_mu("x_2", mu) * mk_pi_0("t_2")
                    for mu in range(4) ])
                + f"pi0 j_mu j_mu pi0",
                sum([
                    sympy.simplify(1)/2 * (
                        mk_pi_p("t_1", True) * mk_j_mu("x_1", mu) * mk_j_mu("x_2", mu) * mk_pi_p("t_2")
                        + mk_pi_m("t_1", True) * mk_j_mu("x_1", mu) * mk_j_mu("x_2", mu) * mk_pi_m("t_2"))
                    for mu in range(4) ])
                + f"pi+ j_mu j_mu pi+",
                sum([
                    sympy.simplify(1)/2 * (
                        mk_k_0("t_1", True) * mk_j_mu("x_1", mu) * mk_j_mu("x_2", mu) * mk_k_0("t_2")
                        + mk_k_0_bar("t_1", True) * mk_j_mu("x_1", mu) * mk_j_mu("x_2", mu) * mk_k_0_bar("t_2")
                        )
                    for mu in range(4) ])
                + f"K0 j_mu j_mu K0",
                sum([
                    sympy.simplify(1)/2 * (
                        mk_k_p("t_1", True) * mk_j_mu("x_1", mu) * mk_j_mu("x_2", mu) * mk_k_p("t_2")
                        + mk_k_m("t_1", True) * mk_j_mu("x_1", mu) * mk_j_mu("x_2", mu) * mk_k_m("t_2"))
                    for mu in range(4) ])
                + f"K+ j_mu j_mu K+",
                sum([
                    mk_j_mu("x_1", mu) * mk_j_mu("x_2", mu)
                    for mu in range(4) ])
                + f"j_mu j_mu",
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.optimize()
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, f"cache/auto_contract_cexpr/get_cexpr_meson_jj_mm")

def get_all_cexpr():
    cexprs = [
            get_cexpr_zeros(),
            get_cexpr_meson_corr(),
            get_cexpr_meson_f_corr(),
            get_cexpr_meson_bk_bpi_corr(),
            get_cexpr_meson_corr2(),
            get_cexpr_meson_corr3(),
            get_cexpr_meson_f_corr2(),
            get_cexpr_meson_quark_mass(),
            get_cexpr_meson_jt_zv(),
            get_cexpr_meson_jj_mm(),
            ]
    check_list = []
    check_ama_list = []
    for cexpr in cexprs:
        check, check_ama = benchmark_eval_cexpr(cexpr)
        check_list.append(check)
        check_ama_list.append(check_ama)
    for cexpr, check, check_ama in zip(cexprs, check_list, check_ama_list):
        names = get_cexpr_names(cexpr)
        for name in names:
            name_str = name.replace('\n', '  ')
            q.displayln_info(f"CHECK: {name_str}")
        q.displayln_info(f"CHECK: {benchmark_show_check(check)}")
        q.displayln_info(f"CHECK: {benchmark_show_check(check_ama)}")

size_node_list = [
        [1, 1, 1, 1],
        [1, 1, 1, 2],
        [1, 1, 1, 4],
        [1, 1, 1, 8],
        [2, 2, 2, 2],
        [2, 2, 2, 4],
        ]

q.begin(sys.argv, size_node_list)

q.qremove_all_info("cache")

get_all_cexpr()

q.timer_display()

q.displayln_info(f"CHECK: finished successfully.")

q.end()