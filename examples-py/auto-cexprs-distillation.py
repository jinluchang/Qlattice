#!/usr/bin/env python3

import qlat as q

from auto_contractor.operators import *
from auto_contractor.eval import *

import sys

is_cython = False
is_distillation = True

@q.timer
def get_cexpr_zeros():
    fn_base = f"cache/auto_contract_cexpr/get_cexpr_zeros"
    def calc_cexpr():
        exprs = [
                mk_expr(0),
                mk_k_p("t_2", True)     * mk_k_0("t_1")     + "k+^dag  * k0    ",
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit=True)
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython, is_distillation=is_distillation)

@q.timer
def get_cexpr_meson_corr():
    fn_base = f"cache/auto_contract_cexpr/get_cexpr_meson_corr"
    def calc_cexpr():
        exprs = [
                mk_pi_p("t_2", True)    * mk_pi_p("t_1")    + "pi^dag * pi   ",
                mk_k_p("t_2", True)     * mk_k_p("t_1")     + "k^dag  * k    ",
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit=True)
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython, is_distillation=is_distillation)

@q.timer
def get_cexpr_meson_f_corr():
    fn_base = f"cache/auto_contract_cexpr/get_cexpr_meson_f_corr"
    def calc_cexpr():
        exprs = [
                mk_j5pi_mu("x_2", 3)    * mk_pi_p("t_1")    + "a_pi   * pi   ",
                mk_j5k_mu("x_2", 3)     * mk_k_p("t_1")     + "a_k    * k    ",
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit=True)
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython, is_distillation=is_distillation)

@q.timer
def get_cexpr_meson_corr2():
    fn_base = f"cache/auto_contract_cexpr/get_cexpr_meson_corr2"
    def calc_cexpr():
        exprs = [
                mk_pi_p("x_2", True)    * mk_pi_p("x_1")    + "pi     * pi   ",
                mk_k_p("x_2", True)     * mk_k_p("x_1")     + "k      * k    ",
                mk_a0_p("x_2", True)    * mk_a0_p("x_1")    + "a0     * a0   ",
                mk_kappa_p("x_2", True) * mk_kappa_p("x_1") + "kappa  * kappa",
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit=True)
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython, is_distillation=is_distillation)

@q.timer
def get_cexpr_meson_corr3():
    fn_base = f"cache/auto_contract_cexpr/get_cexpr_meson_corr3"
    def calc_cexpr():
        exprs = [
                mk_pi_0("t_1", True) * mk_pi_0("t_2"),
                sympy.simplify(1)/2 * (
                    mk_k_0("t_1", True) * mk_k_0("t_2")
                    + mk_k_0_bar("t_1", True) * mk_k_0_bar("t_2")
                    ),
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit=True)
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython, is_distillation=is_distillation)

@q.timer
def get_cexpr_meson_f_corr2():
    fn_base = f"cache/auto_contract_cexpr/get_cexpr_meson_f_corr2"
    def calc_cexpr():
        exprs = [
                mk_j5pi_mu("x_2", 3) * mk_pi_p("t_1") + "(a_pi * pi)",
                sympy.simplify(1)/2 * (
                    mk_j5k_mu("x_2", 3)  * mk_k_p("t_1")
                    + mk_j5km_mu("x_2", 3)  * mk_k_m("t_1")
                    ) + "(a_k  * k )",
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit=True)
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython, is_distillation=is_distillation)

@q.timer
def get_cexpr_meson_quark_mass():
    fn_base = f"cache/auto_contract_cexpr/get_cexpr_meson_quark_mass"
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
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit=True)
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython, is_distillation=is_distillation)

def get_all_cexpr():
    cexprs = [
            lambda : get_cexpr_zeros(),
            lambda : get_cexpr_meson_corr(),
            lambda : get_cexpr_meson_f_corr(),
            lambda : get_cexpr_meson_corr2(),
            lambda : get_cexpr_meson_corr3(),
            lambda : get_cexpr_meson_f_corr2(),
            lambda : get_cexpr_meson_quark_mass(),
            ]
    for cexpr in cexprs:
        cexpr = cexpr()
        check, check_ama = benchmark_eval_cexpr(cexpr)
        names = get_expr_names(cexpr)
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

q.begin_with_mpi(size_node_list)

q.qremove_all_info("cache")

get_all_cexpr()

q.timer_display()

q.end_with_mpi()

q.displayln_info(f"CHECK: finished successfully.")
