#!/usr/bin/env python3

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

f1 = "u"
f2 = "d"
f3 = "s"

p1 = "x1"
p2 = "x2"
p3 = "x3"

diagram_type_dict = dict()
diagram_type_dict[()] = 'Type0'
diagram_type_dict[((('x1', 'x2'), 1), (('x2', 'x1'), 1))] = 'Type1'
diagram_type_dict[((('x1', 'x2'), 1), (('x2', 'x3'), 1), (('x3', 'x1'), 1))] = 'Type2'
diagram_type_dict[((('x1', 'x3'), 1), (('x2', 'x1'), 1), (('x3', 'x2'), 1))] = 'Type3'
diagram_type_dict[((('x1', 'x2'), 1), (('x2', 'x1'), 1), (('x3', 'x3'), 1))] = 'Type4'

exprs = [
    qac.mk_expr(1),
    qac.mk_pi_p(p2, is_dagger=True) * qac.mk_pi_p(p1),
    qac.mk_pi_p(p2, is_dagger=True) * qac.mk_pi_p(p1) * qac.mk_scalar(f1, f1, p3),
    qac.mk_pi_p(p2, is_dagger=True) * qac.mk_pi_p(p1) * qac.mk_scalar(f2, f2, p3),
    qac.mk_pi_p(p2, is_dagger=True) * qac.mk_pi_p(p1) * qac.mk_scalar(f3, f3, p3),
    (qac.mk_pi_p(p2, is_dagger=True) * qac.mk_pi_p(p1) * qac.mk_vec_mu(f1, f1, p3, 3), None, "Type4", [ "Type2", "Type3", ]),
    (qac.mk_pi_p(p2, is_dagger=True) * qac.mk_pi_p(p1) * qac.mk_vec_mu(f2, f2, p3, 3), None, "Type4", [ "Type2", "Type3", ]),
    (qac.mk_pi_p(p2, is_dagger=True) * qac.mk_pi_p(p1) * qac.mk_vec_mu(f3, f3, p3, 3), None, "Type4", [ "Type2", "Type3", ]),
]

for expr in exprs:
    q.json_results_append(str(expr))

for expr in qac.contract_simplify(*exprs):
    q.json_results_append(str(expr))

cexpr = qac.contract_simplify_compile(*exprs, is_isospin_symmetric_limit=True, diagram_type_dict=diagram_type_dict)

q.json_results_append(
    qac.display_cexpr(cexpr)
)

for v in cexpr.list():
    q.json_results_append(str(v))

q.json_results_append(f"diagram_type_dict = qac.get_diagram_type_dict(cexpr)")
for k, v in qac.get_diagram_type_dict(cexpr).items():
    q.json_results_append(f"diagram_type_dict[{k!r}] = {v!r}")

q.json_results_append(f"qac.get_expr_names(cexpr)")
for name in qac.get_expr_names(cexpr):
    q.json_results_append(name)

cexpr_opt = cexpr.copy()
cexpr_opt.optimize()

q.json_results_append(
    qac.display_cexpr(cexpr_opt)
)

for v in cexpr_opt.list():
    q.json_results_append(str(v))

q.json_results_append(f"diagram_type_dict = qac.get_diagram_type_dict(cexpr_opt)")
for k, v in qac.get_diagram_type_dict(cexpr_opt).items():
    q.json_results_append(f"diagram_type_dict[{k!r}] = {v!r}")

q.json_results_append(f"qac.get_expr_names(cexpr_opt)")
for name in qac.get_expr_names(cexpr_opt):
    q.json_results_append(name)

@q.timer
def get_cexpr_test(is_cython=False):
    fn_base = f"cache/auto_contract_cexpr/get_cexpr_test"
    def calc_cexpr():
        cexpr = qac.contract_simplify_compile(*exprs, is_isospin_symmetric_limit=True, diagram_type_dict=diagram_type_dict)
        return cexpr
    return qac.cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

for is_cython in [ False, True, ]:
    q.json_results_append(
        qac.cexpr_code_gen_py(cexpr_opt, is_cython=is_cython)
    )
    ccexpr = get_cexpr_test(is_cython=is_cython)
    q.json_results_append(f"diagram_type_dict = qac.get_diagram_type_dict(ccexpr)")
    for k, v in qac.get_diagram_type_dict(ccexpr).items():
        q.json_results_append(f"diagram_type_dict[{k!r}] = {v!r}")
    q.json_results_append(f"qac.get_expr_names(ccexpr)")
    for name in qac.get_expr_names(ccexpr):
        q.json_results_append(name)
    check, check_ama = qac.benchmark_eval_cexpr(ccexpr)
    q.json_results_append(f"get_cexpr_test benchmark_eval_cexpr check get_data_sig is_cython={is_cython}", q.get_data_sig(np.array(check, dtype=np.complex128), q.RngState()))
    q.json_results_append(f"get_cexpr_test benchmark_eval_cexpr check_ama get_data_sig is_cython={is_cython}", q.get_data_sig(np.array(check_ama, dtype=np.complex128), q.RngState()))

def get_prop(flavor, p1, p2):
    tag = f"get_prop({flavor!r},{p1!r},{p2!r})"
    q.displayln_info(f"Call {tag}")
    rs = q.RngState(tag)
    wm = qac.make_rand_spin_color_matrix(rs)
    assert isinstance(wm, q.WilsonMatrix)
    return wm

# S_l(x1,x2)
wm1 = get_prop(
    "l",
    ("point", q.Coordinate([ 1, 2, 3, 4, ]),),
    ("point", q.Coordinate([ 3, 7, 2, 1, ]),),
)

# S_l(x2,x1)
wm2 = get_prop(
    "l",
    ("point", q.Coordinate([ 3, 7, 2, 1, ]),),
    ("point", q.Coordinate([ 1, 2, 3, 4, ]),),
)

# tr(gamma_5*S_l(x1,x2)*gamma_5*S_l(x2,x1))

c_pi = q.mat_tr_wm_wm(
    q.mat_mul_sm_wm(q.get_gamma_matrix(5), wm1),
    q.mat_mul_sm_wm(q.get_gamma_matrix(5), wm2),
)

q.json_results_append(f"get_prop c_pi", q.get_data_sig(c_pi, q.RngState()))
q.displayln_info(-1, f"{c_pi}")

pd = {
    "x1": ("point", q.Coordinate([ 1, 2, 3, 4, ]),),
    "x2": ("point", q.Coordinate([ 3, 7, 2, 1, ]),),
    "x3": ("point-snk", q.Coordinate([ 3, 7, 2, 1, ]),),
}

res = qac.eval_cexpr(ccexpr=get_cexpr_test(), positions_dict=pd, get_prop=get_prop)

for idx, v in enumerate(res):
    q.json_results_append(f"eval_cexpr res[{idx}]", q.get_data_sig(v, q.RngState()))
    q.displayln_info(-1, f"{v}")

q.check_log_json(__file__, check_eps=1e-14)

q.timer_display()

q.end_with_mpi()

q.displayln_info(f"CHECK: finished successfully.")
