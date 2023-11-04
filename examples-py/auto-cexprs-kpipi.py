#!/usr/bin/env python3

import qlat as q

from auto_contractor.operators import *
from auto_contractor.eval import *

import sys

is_cython = False

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

all_mom_list_dict = {
        0: [
            q.Coordinate([ 0, 0, 0, 0, ]),
            ],
        1: [
            q.Coordinate([ 0, 0, 1, 0, ]),
            q.Coordinate([ 0, 1, 0, 0, ]),
            q.Coordinate([ 1, 0, 0, 0, ]),
            q.Coordinate([ 0, 0, -1, 0, ]),
            q.Coordinate([ 0, -1, 0, 0, ]),
            q.Coordinate([ -1, 0, 0, 0, ]),
            ],
        2: [
            q.Coordinate([ 0, 1, 1, 0, ]),
            q.Coordinate([ 1, 0, 1, 0, ]),
            q.Coordinate([ 1, 1, 0, 0, ]),
            q.Coordinate([ 0, -1, 1, 0, ]),
            q.Coordinate([ -1, 0, 1, 0, ]),
            q.Coordinate([ -1, 1, 0, 0, ]),
            q.Coordinate([ 0, 1, -1, 0, ]),
            q.Coordinate([ 1, 0, -1, 0, ]),
            q.Coordinate([ 1, -1, 0, 0, ]),
            q.Coordinate([ 0, -1, -1, 0, ]),
            q.Coordinate([ -1, 0, -1, 0, ]),
            q.Coordinate([ -1, -1, 0, 0, ]),
            ],
        3: [
            q.Coordinate([ 1, 1, 1, 0, ]),
            q.Coordinate([ -1, 1, 1, 0, ]),
            q.Coordinate([ 1, -1, 1, 0, ]),
            q.Coordinate([ -1, -1, 1, 0, ]),
            q.Coordinate([ 1, 1, -1, 0, ]),
            q.Coordinate([ -1, 1, -1, 0, ]),
            q.Coordinate([ 1, -1, -1, 0, ]),
            q.Coordinate([ -1, -1, -1, 0, ]),
            ],
        4: [
            q.Coordinate([ 0, 0, 2, 0, ]),
            q.Coordinate([ 0, 2, 0, 0, ]),
            q.Coordinate([ 2, 0, 0, 0, ]),
            q.Coordinate([ 0, 0, -2, 0, ]),
            q.Coordinate([ 0, -2, 0, 0, ]),
            q.Coordinate([ -2, 0, 0, 0, ]),
            ],
        }

def mk_mom_fac_jj0(mom, p1, p2):
    """
    mom in [ 0, 1, 2, 3, 4, ]
    """
    mom_list = all_mom_list_dict[mom]
    fac = 0
    for mom in mom_list:
        mom1 = mom
        mom2 = -mom
        fac1 = mk_fac(f"momentum_factor({mom1},{p1},size)")
        fac2 = mk_fac(f"momentum_factor({mom2},{p2},size)")
        fac = fac + fac1 * fac2
    fac = 1 / sympy.sqrt(len(mom_list)) * fac
    return fac + f"mom_fac_jj0(mom={mom})"

@q.timer
def get_cexpr_kpipi():
    fn_base = f"cache/auto_contract_cexpr/get_cexpr_kpipi"
    def calc_cexpr():
        vol = 1
        exprs_odd_ops = [
                vol * mk_Q1("x", "odd") + "Q1(o)",
                vol * mk_Q2("x", "odd") + "Q2(o)",
                vol * mk_Q3("x", "odd") + "Q3(o)",
                vol * mk_Q4("x", "odd") + "Q4(o)",
                vol * mk_Q5("x", "odd") + "Q5(o)",
                vol * mk_Q6("x", "odd") + "Q6(o)",
                vol * mk_Q7("x", "odd") + "Q7(o)",
                vol * mk_Q8("x", "odd") + "Q8(o)",
                vol * mk_Q9("x", "odd") + "Q9(o)",
                vol * mk_Q10("x", "odd") + "Q10(o)",
                vol * mk_Qsub("x", "odd") + "Qs(o)",
                ]
        exprs_even_ops = [
                vol * mk_Q1("x", "even") + "Q1(e)",
                vol * mk_Q2("x", "even") + "Q2(e)",
                vol * mk_Q3("x", "even") + "Q3(e)",
                vol * mk_Q4("x", "even") + "Q4(e)",
                vol * mk_Q5("x", "even") + "Q5(e)",
                vol * mk_Q6("x", "even") + "Q6(e)",
                vol * mk_Q7("x", "even") + "Q7(e)",
                vol * mk_Q8("x", "even") + "Q8(e)",
                vol * mk_Q9("x", "even") + "Q9(e)",
                vol * mk_Q10("x", "even") + "Q10(e)",
                vol * mk_Qsub("x", "even") + "Qs(e)",
                ]
        exprs_ops = exprs_odd_ops + exprs_even_ops
        exprs_k = [
                vol * mk_k_0("x2") + "K0",
                ]
        exprs_pipi = [
                vol**2 * mk_mom_fac_jj0(0, "x1_1", "x1_2") * mk_pipi_i0("x1_1", "x1_2", True) + "pipi_I0",
                vol**2 * mk_mom_fac_jj0(0, "x1_1", "x1_2") * mk_pipi_i20("x1_1", "x1_2", True) + "pipi_I2",
                vol**2 * mk_mom_fac_jj0(1, "x1_1", "x1_2") * mk_pipi_i0("x1_1", "x1_2", True) + "pipi_I0_mom1",
                vol**2 * mk_mom_fac_jj0(1, "x1_1", "x1_2") * mk_pipi_i20("x1_1", "x1_2", True) + "pipi_I2_mom1",
                vol * mk_sigma("x1_1", True) + "sigma_1",
                vol * mk_sigma("x1_2", True) + "sigma_2",
                vol * mk_pi_0("x1_1", True) + "pi0_1",
                vol * mk_pi_0("x1_2", True) + "pi0_2",
                mk_expr(1) + "1",
                ]
        exprs = []
        for expr_k in exprs_k:
            for expr_pipi in exprs_pipi:
                for expr_op in exprs_ops:
                    exprs.append(expr_pipi * expr_op * expr_k)
        diagram_type_dict = dict()
        diagram_type_dict[((('x', 'x'), 1), (('x', 'x1_1'), 1), (('x1_1', 'x1_2'), 1), (('x1_2', 'x2'), 1), (('x2', 'x'), 1))] = "Type3"
        diagram_type_dict[((('x', 'x'), 1), (('x', 'x2'), 1), (('x1_1', 'x1_2'), 1), (('x1_2', 'x1_1'), 1), (('x2', 'x'), 1))] = "Type4"
        diagram_type_dict[((('x', 'x1_1'), 1), (('x', 'x1_2'), 1), (('x1_1', 'x'), 1), (('x1_2', 'x2'), 1), (('x2', 'x'), 1))] = "Type1"
        diagram_type_dict[((('x', 'x1_1'), 1), (('x', 'x2'), 1), (('x1_1', 'x1_2'), 1), (('x1_2', 'x'), 1), (('x2', 'x'), 1))] = "Type2"
        diagram_type_dict[((('x', 'x1_1'), 1), (('x1_1', 'x1_2'), 1), (('x1_2', 'x2'), 1), (('x2', 'x'), 1))] = "Type3"
        diagram_type_dict[((('x', 'x2'), 1), (('x1_1', 'x1_2'), 1), (('x1_2', 'x1_1'), 1), (('x2', 'x'), 1))] = "Type4"
        diagram_type_dict[((('x', 'x'), 1), (('x', 'x1_1'), 1), (('x1_1', 'x2'), 1), (('x2', 'x'), 1))] = "Type3"
        diagram_type_dict[((('x', 'x'), 1), (('x', 'x2'), 1), (('x1_1', 'x1_1'), 1), (('x2', 'x'), 1))] = "Type4"
        diagram_type_dict[((('x', 'x1_1'), 1), (('x', 'x2'), 1), (('x1_1', 'x'), 1), (('x2', 'x'), 1))] = "Type2"
        diagram_type_dict[((('x', 'x1_1'), 1), (('x1_1', 'x2'), 1), (('x2', 'x'), 1))] = "Type3"
        diagram_type_dict[((('x', 'x2'), 1), (('x1_1', 'x1_1'), 1), (('x2', 'x'), 1))] = "Type4"
        diagram_type_dict[((('x', 'x'), 1), (('x', 'x1_2'), 1), (('x1_2', 'x2'), 1), (('x2', 'x'), 1))] = "Type3"
        diagram_type_dict[((('x', 'x'), 1), (('x', 'x2'), 1), (('x1_2', 'x1_2'), 1), (('x2', 'x'), 1))] = "Type4"
        diagram_type_dict[((('x', 'x1_2'), 1), (('x', 'x2'), 1), (('x1_2', 'x'), 1), (('x2', 'x'), 1))] = "Type2"
        diagram_type_dict[((('x', 'x1_2'), 1), (('x1_2', 'x2'), 1), (('x2', 'x'), 1))] = "Type3"
        diagram_type_dict[((('x', 'x2'), 1), (('x1_2', 'x1_2'), 1), (('x2', 'x'), 1))] = "Type4"
        diagram_type_dict[((('x', 'x'), 1), (('x', 'x2'), 1), (('x2', 'x'), 1))] = "Type4"
        diagram_type_dict[((('x', 'x2'), 1), (('x2', 'x'), 1))] = "Type4"
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit=True, diagram_type_dict=diagram_type_dict)
        return cexpr
    base_positions_dict = {}
    base_positions_dict["momentum_factor"] = momentum_factor
    base_positions_dict["Coordinate"] = q.Coordinate
    return cache_compiled_cexpr(calc_cexpr, fn_base, base_positions_dict=base_positions_dict, is_cython=is_cython)

def get_all_cexpr():
    cexprs = [
            lambda : get_cexpr_kpipi(),
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
