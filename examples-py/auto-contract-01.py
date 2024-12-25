#!/usr/bin/env python3

json_results = []
check_eps = 1e-14

def json_results_append(*args):
    json_results.append(args)

import sys
import numpy as np
import qlat as q

import auto_contractor as qac

size_node_list = [
        [1, 1, 1, 1],
        [1, 1, 1, 2],
        [1, 1, 1, 4],
        [1, 1, 1, 8],
        [2, 2, 2, 2],
        [2, 2, 2, 4],
        ]

q.begin_with_mpi(size_node_list)

f = "u" # flavor, can be "u", "d", "s", "c"

p1 = "x1"
s1 = "s1"
c1 = "c1"

p2 = "x2"
s2 = "s2"
c2 = "c2"

q1v = qac.Qv(f, p1, s1, c1)
q1b = qac.Qb(f, p1, s1, c1)
q2v = qac.Qv(f, p2, s2, c2)
q2b = qac.Qb(f, p2, s2, c2)

exprs = [
    qac.mk_expr(1),
    q1v * q2b,
    q1v * q2b + q2b * q1v,
]

for expr in exprs:
    json_results_append(str(expr))

for expr in qac.contract_simplify(*exprs):
    json_results_append(str(expr))

cexpr = qac.contract_simplify_compile(*exprs)

json_results_append(
    qac.display_cexpr(cexpr)
)

for v in cexpr.list():
    json_results_append(str(v))

cexpr_opt = cexpr.copy()
cexpr_opt.optimize()

json_results_append(
    qac.display_cexpr(cexpr_opt)
)

for v in cexpr_opt.list():
    json_results_append(str(v))

for is_cython in [ False, True, ]:
    json_results_append(
        qac.cexpr_code_gen_py(cexpr_opt, is_cython=is_cython)
    )
    @q.timer
    def get_cexpr_test():
        fn_base = f"cache/auto_contract_cexpr/get_cexpr_test"
        def calc_cexpr():
            cexpr = qac.contract_simplify_compile(*exprs, is_isospin_symmetric_limit=True)
            return cexpr
        return qac.cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

    cexpr = get_cexpr_test()
    check, check_ama = qac.benchmark_eval_cexpr(cexpr)
    json_results_append(f"get_cexpr_test benchmark_eval_cexpr check get_data_sig is_cython={is_cython}", q.get_data_sig(np.array(check, dtype=np.complex128), q.RngState()))
    json_results_append(f"get_cexpr_test benchmark_eval_cexpr check_ama get_data_sig is_cython={is_cython}", q.get_data_sig(np.array(check_ama, dtype=np.complex128), q.RngState()))

q.check_log_json(__file__, json_results, check_eps=check_eps)

q.timer_display()

q.end_with_mpi()

q.displayln_info(f"CHECK: finished successfully.")
