#!/usr/bin/env python3

import sys
import sympy
import math
import cmath
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

diagram_type_dict = dict()

jj_op = sum([
    qac.mk_fac(f"rel_mod_sym(xx_1[1][{mu}]-xx_2[1][{mu}],size[{mu}])")
    * qac.mk_fac(f"rel_mod_sym(xx_1[1][{nu}]-xx_2[1][{nu}],size[{nu}])")
    * qac.mk_j_mu("xx_1", mu) * qac.mk_j_mu("xx_2", nu)
    for mu in range(3) for nu in range(3)
    ]) + f"(xx_1-xx_2)_i * (xx_1-xx_2)_j * j_i(xx_1) * j_j(xx_2)"

exprs = [
        qac.mk_expr(1) + f"1",
        qac.mk_proton("x_2", "u", is_dagger=True) * qac.mk_proton("x_1", "u"),
        qac.mk_m("u", "xx") * qac.mk_proton("x_2", "u", is_dagger=True) * qac.mk_proton("x_1", "u"),
        qac.mk_m("d", "xx") * qac.mk_proton("x_2", "u", is_dagger=True) * qac.mk_proton("x_1", "u"),
        qac.mk_m("s", "xx") * qac.mk_proton("x_2", "u", is_dagger=True) * qac.mk_proton("x_1", "u"),
        jj_op * qac.mk_proton("x_2", "u", is_dagger=True) * qac.mk_proton("x_1", "u"),
        ]

for expr in exprs:
    q.json_results_append(str(expr)[:256])

for expr in qac.contract_simplify(*exprs):
    q.json_results_append(str(expr)[:256])

cexpr = qac.contract_simplify_compile(*exprs, is_isospin_symmetric_limit=True, diagram_type_dict=diagram_type_dict)

q.json_results_append(
        qac.display_cexpr(cexpr)[:256]
        )

for v in cexpr.list():
    q.json_results_append(str(v)[:256])

q.json_results_append(f"diagram_type_dict = qac.get_diagram_type_dict(cexpr)")
for k, v in qac.get_diagram_type_dict(cexpr).items():
    q.json_results_append(f"diagram_type_dict[{k!r}] = {v!r}")

q.json_results_append(f"qac.get_expr_names(cexpr)")
for name in qac.get_expr_names(cexpr):
    q.json_results_append(name)

cexpr_opt = cexpr.copy()
cexpr_opt.optimize()

q.json_results_append(
        qac.display_cexpr(cexpr_opt)[:256]
        )

for v in cexpr_opt.list():
    q.json_results_append(str(v)[:256])

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
    base_positions_dict = dict()
    return qac.cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython, base_positions_dict=base_positions_dict)

for is_cython in [ False, True, ]:
    q.json_results_append(
            qac.cexpr_code_gen_py(cexpr_opt, is_cython=is_cython)[:256]
            )
    ccexpr = get_cexpr_test(is_cython=is_cython)
    q.json_results_append(f"diagram_type_dict = qac.get_diagram_type_dict(ccexpr)")
    for k, v in qac.get_diagram_type_dict(ccexpr).items():
        q.json_results_append(f"diagram_type_dict[{k!r}] = {v!r}")
    q.json_results_append(f"qac.get_expr_names(ccexpr)")
    for name in qac.get_expr_names(ccexpr):
        q.json_results_append(name)
    base_positions_dict = dict()
    check, check_ama = qac.benchmark_eval_cexpr(ccexpr, base_positions_dict=base_positions_dict)
    q.json_results_append(f"get_cexpr_test benchmark_eval_cexpr check get_data_sig is_cython={is_cython}", q.get_data_sig(np.array(check, dtype=np.complex128), q.RngState()))
    q.json_results_append(f"get_cexpr_test benchmark_eval_cexpr check_ama get_data_sig is_cython={is_cython}", q.get_data_sig(np.array(check_ama, dtype=np.complex128), q.RngState()))

def get_prop(flavor, p1, p2):
    tag = f"get_prop({flavor!r},{p1!r},{p2!r})"
    q.displayln_info(f"Call {tag}")
    rs = q.RngState(tag)
    wm = qac.make_rand_spin_color_matrix(rs)
    assert isinstance(wm, q.WilsonMatrix)
    return wm

pd = {
        "x_1": ("point", q.Coordinate([ 1, 2, 3, 1, ]),),
        "x_2": ("point", q.Coordinate([ 0, 2, 1, 3, ]),),
        "xx_1": ("point", q.Coordinate([ 2, 0, 2, 2, ]),),
        "xx_2": ("point", q.Coordinate([ 1, 3, 2, 2, ]),),
        "xx": ("point", q.Coordinate([ 2, 0, 1, 2, ]),),
        "size": q.Coordinate([ 4, 4, 4, 8, ]),
        }

res = qac.eval_cexpr(ccexpr=get_cexpr_test(), positions_dict=pd, get_prop=get_prop)

for idx, v in enumerate(res):
    q.json_results_append(f"eval_cexpr res[{idx}]", q.get_data_sig(v, q.RngState()))
    q.displayln_info(-1, f"{v}")

q.check_log_json(__file__, check_eps=1e-14)

q.timer_display()

q.end_with_mpi()

q.displayln_info(f"CHECK: finished successfully.")
