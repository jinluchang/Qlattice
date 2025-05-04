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

def wave_function(p1, p2, radius, size):
    p1_tag, c1 = p1
    p2_tag, c2 = p2
    c12 = q.smod_coordinate(c1 - c2, size)
    # assert c12[3] == 0
    c12_r_sqr = c12.r_sqr()
    dis = math.sqrt(c12_r_sqr)
    vol = size[0] * size[1] * size[2]
    wf = vol**2 * math.exp(-dis / radius)
    return wf

def momentum_factor(mom, p, size, is_dagger=False):
    p_tag, c = p
    assert mom[3] == 0
    phase = mom[0] * c[0] / size[0] + mom[1] * c[1] / size[1] + mom[2] * c[2] / size[2]
    phase = phase * 2.0 * math.pi
    if not is_dagger:
        mf = cmath.rect(1.0, phase)
    else:
        mf = cmath.rect(1.0, -phase)
    return mf

def mk_meson_wf(f1, f2, p1, p2, radius, mom, is_dagger=False):
    """
    i q1bar g5 q2 #dag: i q2bar g5 q1
    return the actual dagger of the operator
    """
    s1 = qac.new_spin_index()
    s2 = qac.new_spin_index()
    c = qac.new_color_index()
    g5 = qac.G(5, s1, s2)
    wf = qac.mk_fac(f"wave_function({p1},{p2},{radius},size)")
    if not is_dagger:
        q1b = qac.Qb(f1, p1, s1, c)
        q2v = qac.Qv(f2, p2, s2, c)
        mf = qac.mk_fac(f"momentum_factor({mom},{p2},size)")
        return sympy.I * wf * mf * q1b * g5 * q2v + f"(i {f1}bar g5 {f2})({p1},{p2})"
    else:
        q2b = qac.Qb(f2, p2, s1, c)
        q1v = qac.Qv(f1, p1, s2, c)
        mf = qac.mk_fac(f"momentum_factor(-{mom},{p2},size)")
        return sympy.I * wf * mf * q2b * g5 * q1v + f"(i {f2}bar g5 {f1})({p2},{p1},{radius},{mom})"

def mk_pi_0_wf(p1, p2, mom, is_dagger=False):
    """
    i/sqrt(2) * (ubar g5 u - dbar g5 d)  #dag: same
    """
    radius = "r_pi"
    return 1 / sympy.sqrt(2) * (
            mk_meson_wf("u", "u", p1, p2, radius, mom, is_dagger)
            - mk_meson_wf("d", "d", p1, p2, radius, mom, is_dagger)
            ) + f"pi0({p1},{p2},{radius},{mom}){qac.show_dagger(is_dagger)}"

diagram_type_dict = dict()
diagram_type_dict[()] = 'Type0'
diagram_type_dict[((('x12', 'x22'), 1), (('x21', 'x11'), 1))] = 'Type1'

exprs = [
        qac.mk_expr(1) + f"1",
        mk_pi_0_wf("x21", "x22", "mom", is_dagger=True) * mk_pi_0_wf("x11", "x12", "mom"),
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
    base_positions_dict = dict()
    base_positions_dict['wave_function'] = wave_function
    base_positions_dict['momentum_factor'] = momentum_factor
    base_positions_dict['r_pi'] = 1.5
    return qac.cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython, base_positions_dict=base_positions_dict)

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
    base_positions_dict = dict()
    base_positions_dict['mom'] = q.CoordinateD([ 0, 0, 1, 0, ])
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
    "x11": ("point", q.Coordinate([ 1, 2, 3, 1, ]),),
    "x12": ("point", q.Coordinate([ 3, 1, 2, 1, ]),),
    "x21": ("point", q.Coordinate([ 2, 2, 3, 4, ]),),
    "x22": ("point", q.Coordinate([ 3, 2, 0, 4, ]),),
    "mom": q.CoordinateD([ 0, 0, 1, 0, ]),
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
