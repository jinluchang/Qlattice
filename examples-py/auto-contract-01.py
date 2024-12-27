#!/usr/bin/env python3

json_results = []
check_eps = 1e-14

def json_results_append(*args):
    q.displayln_info(r"//------------------------------------------------------------\\")
    q.displayln_info(-1, *args)
    q.displayln_info(r"\\------------------------------------------------------------//")
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

cexpr = qac.contract_simplify_compile(*exprs, is_isospin_symmetric_limit=True)

json_results_append(
    qac.display_cexpr(cexpr)
)

for v in cexpr.list():
    json_results_append(str(v))

json_results_append(f"diagram_type_dict = qac.get_diagram_type_dict(cexpr)")
for k, v in qac.get_diagram_type_dict(cexpr).items():
    json_results_append(f"diagram_type_dict[{k!r}] = {v!r}")

json_results_append(f"qac.get_expr_names(cexpr)")
for name in qac.get_expr_names(cexpr):
    json_results_append(name)

cexpr_opt = cexpr.copy()
cexpr_opt.optimize()

json_results_append(
    qac.display_cexpr(cexpr_opt)
)

for v in cexpr_opt.list():
    json_results_append(str(v))

json_results_append(f"diagram_type_dict = qac.get_diagram_type_dict(cexpr_opt)")
for k, v in qac.get_diagram_type_dict(cexpr_opt).items():
    json_results_append(f"diagram_type_dict[{k!r}] = {v!r}")

json_results_append(f"qac.get_expr_names(cexpr_opt)")
for name in qac.get_expr_names(cexpr_opt):
    json_results_append(name)

@q.timer
def get_cexpr_test(is_cython=False):
    fn_base = f"cache/auto_contract_cexpr/get_cexpr_test"
    def calc_cexpr():
        cexpr = qac.contract_simplify_compile(*exprs, is_isospin_symmetric_limit=True)
        return cexpr
    return qac.cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

for is_cython in [ False, True, ]:
    json_results_append(
        qac.cexpr_code_gen_py(cexpr_opt, is_cython=is_cython)
    )
    ccexpr = get_cexpr_test(is_cython=is_cython)
    json_results_append(f"diagram_type_dict = qac.get_diagram_type_dict(ccexpr)")
    for k, v in qac.get_diagram_type_dict(ccexpr).items():
        json_results_append(f"diagram_type_dict[{k!r}] = {v!r}")
    json_results_append(f"qac.get_expr_names(ccexpr)")
    for name in qac.get_expr_names(ccexpr):
        json_results_append(name)
    check, check_ama = qac.benchmark_eval_cexpr(ccexpr)
    json_results_append(f"get_cexpr_test benchmark_eval_cexpr check get_data_sig is_cython={is_cython}", q.get_data_sig(np.array(check, dtype=np.complex128), q.RngState()))
    json_results_append(f"get_cexpr_test benchmark_eval_cexpr check_ama get_data_sig is_cython={is_cython}", q.get_data_sig(np.array(check_ama, dtype=np.complex128), q.RngState()))

def get_prop(flavor, p1, p2):
    tag = f"get_prop({flavor!r},{p1!r},{p2!r})"
    q.displayln_info(f"Call {tag}")
    rs = q.RngState(tag)
    wm = qac.make_rand_spin_color_matrix(rs)
    assert isinstance(wm, q.WilsonMatrix)
    return wm

wm = get_prop(
    "l",
    ("point-snk", q.Coordinate([ 1, 2, 3, 4, ]),),
    ("point", q.Coordinate([ 3, 7, 2, 1, ]),),
)

json_results_append(f"get_prop wm: {q.get_data_sig(wm, q.RngState())}")
q.displayln_info(-1, f"{wm}")

pd = {
    "x1": ("point-snk", q.Coordinate([ 1, 2, 3, 4, ]),),
    "x2": ("point", q.Coordinate([ 3, 7, 2, 1, ]),),
}

res = qac.eval_cexpr(ccexpr=get_cexpr_test(), positions_dict=pd, get_prop=get_prop)

for idx, v in enumerate(res):
    json_results_append(f"eval_cexpr res[{idx}]: {q.get_data_sig(v, q.RngState())}")
    q.displayln_info(-1, f"{v}")

q.check_log_json(__file__, json_results, check_eps=check_eps)

q.timer_display()

q.end_with_mpi()

q.displayln_info(f"CHECK: finished successfully.")
