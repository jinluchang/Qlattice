#!/usr/bin/env python3

import qlat as q

from auto_contractor.operators import *
from auto_contractor.eval import *

import sys

is_cython = False

@q.timer
def get_cexpr_meson_uu_corr():
    fn_base = f"cache/auto_contract_cexpr/get_cexpr_meson_uu_corr"
    def calc_cexpr():
        exprs = [
                mk_pi_p("t_2", True) * mk_vec_uu_mu("u", "u", "x1", "x2", 3) * mk_pi_p("t_1"),
                mk_k_p("t_2", True)  * mk_vec_uu_mu("u", "u", "x1", "x2", 3) * mk_k_p("t_1"),
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit=True)
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

def get_all_cexpr():
    cexprs = [
            lambda : get_cexpr_meson_uu_corr(),
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
