#!/usr/bin/env python3

json_results = []
check_eps = 1e-14

def json_results_append(*args):
    q.displayln_info(0, r"//------------------------------------------------------------\\")
    q.displayln_info(0, *args)
    q.displayln_info(0, r"\\------------------------------------------------------------//")
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

diagram_type_dict = dict()
diagram_type_dict[()] = '1'
diagram_type_dict[((('t_1', 'x_1'), 1), (('x_1', 't_1'), 1), (('x_2', 'x_2'), 1))] = 'TypeD'
diagram_type_dict[((('t_1', 'x_1'), 1), (('x_1', 'x_2'), 1), (('x_2', 't_1'), 1))] = 'TypeC'
diagram_type_dict[((('t_2', 'x_1'), 1), (('x_1', 't_2'), 1), (('x_2', 'x_2'), 1))] = 'TypeD'
diagram_type_dict[((('t_2', 'x_1'), 1), (('x_1', 'x_2'), 1), (('x_2', 't_2'), 1))] = 'TypeC'

jj_d_list = [
        sum([
            q.epsilon_tensor(mu, nu, rho)
            * qac.mk_fac(f"rel_mod_sym(x_2[1][{mu}] - x_1[1][{mu}], size[{mu}])")
            * qac.mk_j_mu("x_2", nu) * qac.mk_j_mu("x_1", rho)
            for mu in range(3) for nu in range(3) for rho in range(3) ])
        + "e(i,j,k) * x[i] * j_j(x) * j_k(0)",
        ]

pi0d_list = [
        qac.mk_pi_0("t_1") + "pi0(-tsep)",
        qac.mk_pi_0("t_2") + "pi0(x[t]+tsep)",
        ]

exprs_list_pi0_decay = [
        (jj_d * pi0d, None, "TypeC", "TypeD",)
        for pi0d in pi0d_list for jj_d in jj_d_list
        ]

exprs = [
        qac.mk_expr(1) + f"1",
        ]
exprs += exprs_list_pi0_decay

for expr in exprs:
    json_results_append(str(expr)[:256])

for expr in qac.contract_simplify(*exprs):
    json_results_append(str(expr)[:256])

cexpr = qac.contract_simplify_compile(*exprs, is_isospin_symmetric_limit=True, diagram_type_dict=diagram_type_dict)

json_results_append(
        qac.display_cexpr(cexpr)[:256]
        )

for v in cexpr.list():
    json_results_append(str(v)[:256])

json_results_append(f"diagram_type_dict = qac.get_diagram_type_dict(cexpr)")
for k, v in qac.get_diagram_type_dict(cexpr).items():
    json_results_append(f"diagram_type_dict[{k!r}] = {v!r}")

json_results_append(f"qac.get_expr_names(cexpr)")
for name in qac.get_expr_names(cexpr):
    json_results_append(name)

cexpr_opt = cexpr.copy()
cexpr_opt.optimize()

json_results_append(
        qac.display_cexpr(cexpr_opt)[:256]
        )

for v in cexpr_opt.list():
    json_results_append(str(v)[:256])

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
        cexpr = qac.contract_simplify_compile(*exprs, is_isospin_symmetric_limit=True, diagram_type_dict=diagram_type_dict)
        return cexpr
    return qac.cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython, base_positions_dict=None)

for is_cython in [ False, True, ]:
    json_results_append(
        qac.cexpr_code_gen_py(cexpr_opt, is_cython=is_cython)[:256]
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

pd = {
    "x_1": ("point", q.Coordinate([ 1, 2, 3, 2, ]),),
    "x_2": ("point-snk", q.Coordinate([ 3, 1, 2, 1, ]),),
    "t_1": ("wall", 7,),
    "t_2": ("wall", 4,),
    "size": q.Coordinate([ 4, 4, 4, 8, ]),
}

res = qac.eval_cexpr(ccexpr=get_cexpr_test(), positions_dict=pd, get_prop=get_prop)

for idx, v in enumerate(res):
    json_results_append(f"eval_cexpr res[{idx}]", q.get_data_sig(v, q.RngState()))
    q.displayln_info(-1, f"{v}")

q.check_log_json(__file__, json_results, check_eps=check_eps)

q.timer_display()

q.end_with_mpi()

q.displayln_info(f"CHECK: finished successfully.")
