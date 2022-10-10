#!/usr/bin/env python3

import qlat as q

from auto_contractor.operators import *
from auto_contractor.eval import *

import sys

@q.timer
def get_cexpr_meson_corr():
    def calc_cexpr():
        exprs = [
                mk_pi_p("t_2", True)    * mk_pi_p("t_1")    + "pi^dag * pi   ",
                mk_k_p("t_2", True)     * mk_k_p("t_1")     + "k^dag  * k    ",
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.collect_op()
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, f"cache/auto_contract_cexpr/meson_corr-cexpr")

@q.timer
def get_cexpr_meson_f_corr():
    def calc_cexpr():
        exprs = [
                mk_j5pi_mu("x_2", 3)    * mk_pi_p("t_1")    + "a_pi   * pi   ",
                mk_j5k_mu("x_2", 3)     * mk_k_p("t_1")     + "a_k    * k    ",
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.collect_op()
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, f"cache/auto_contract_cexpr/meson_f_corr-cexpr")

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
        cexpr.collect_op()
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, f"cache/auto_contract_cexpr/meson_bk_bpi-cexpr")

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
        cexpr.collect_op()
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, f"cache/auto_contract_cexpr/meson_corr2-cexpr")

def get_all_cexpr():
    cexprs = [
            get_cexpr_meson_corr(),
            get_cexpr_meson_f_corr(),
            get_cexpr_meson_bk_bpi_corr(),
            get_cexpr_meson_corr2(),
            ]
    check_list = []
    for cexpr in cexprs:
        check = benchmark_eval_cexpr(cexpr)
        check_list.append(check)
    for cexpr, check in zip(cexprs, check_list):
        names = get_cexpr_names(cexpr)
        for name in names:
            name_str = name.replace('\n', '  ')
            q.displayln_info(f"CHECK: {name_str}")
        q.displayln_info(f"CHECK: {benchmark_show_check(check)}")

size_node_list = [
        [1, 1, 1, 1],
        [1, 1, 1, 2],
        [1, 1, 1, 4],
        [1, 1, 1, 8],
        [2, 2, 2, 2],
        [2, 2, 2, 4],
        ]

q.begin(sys.argv, size_node_list)

get_all_cexpr()

q.qremove_all_info("cache")

q.timer_display()

q.displayln_info(f"CHECK: finished successfully.")

q.end()
