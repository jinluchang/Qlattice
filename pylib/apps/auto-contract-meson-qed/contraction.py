#!/usr/bin/env python3

from auto_contractor.operators import *

def test_meson_corr():
    diagram_type_dict = dict()
    diagram_type_dict[((('t_1', 't_2'), 1), (('t_2', 't_1'), 1))] = 'Type1'
    exprs = [
            mk_pi_0("t_1", True) * mk_pi_0("t_2"),
            mk_k_0("t_1", True) * mk_k_0("t_2"),
            ]
    cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True, diagram_type_dict = diagram_type_dict)
    print("meson_corr")
    print(display_cexpr(cexpr))

def test_meson_f_corr():
    diagram_type_dict = dict()
    diagram_type_dict[((('t_1', 'x_2'), 1), (('x_2', 't_1'), 1))] = 'Type1'
    exprs = [
            mk_j5pi_mu("x_2", 3) * mk_pi_p("t_1") + "(a_pi * pi)",
            mk_j5k_mu("x_2", 3)  * mk_k_p("t_1")  + "(a_k  * k )",
            ]
    cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True, diagram_type_dict = diagram_type_dict)
    print("meson_f_corr")
    print(display_cexpr(cexpr))

def test_meson_m():
    diagram_type_dict = dict()
    diagram_type_dict[((('t_1', 't_2'), 1), (('t_2', 'x_1'), 1), (('x_1', 't_1'), 1))] = 'Type1'
    diagram_type_dict[((('t_1', 't_2'), 1), (('t_2', 't_1'), 1), (('x_1', 'x_1'), 1))] = None
    exprs = [
            mk_pi_0("t_1", True) * mk_m("u", "x_1") * mk_pi_0("t_2"),
            mk_pi_0("t_1", True) * mk_m("d", "x_1") * mk_pi_0("t_2"),
            mk_pi_p("t_1", True) * mk_m("u", "x_1") * mk_pi_p("t_2"),
            mk_pi_p("t_1", True) * mk_m("d", "x_1") * mk_pi_p("t_2"),
            mk_k_0("t_1", True) * mk_m("d", "x_1") * mk_k_0("t_2"),
            mk_k_0("t_1", True) * mk_m("s", "x_1") * mk_k_0("t_2"),
            mk_k_p("t_1", True) * mk_m("u", "x_1") * mk_k_p("t_2"),
            mk_k_p("t_1", True) * mk_m("s", "x_1") * mk_k_p("t_2"),
            ]
    cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True, diagram_type_dict = diagram_type_dict)
    print("meson_m")
    print(display_cexpr(cexpr))

def test_meson_jt():
    diagram_type_dict = dict()
    diagram_type_dict[((('t_1', 't_2'), 1), (('t_2', 'x_1'), 1), (('x_1', 't_1'), 1))] = 'Type1'
    diagram_type_dict[((('t_1', 't_2'), 1), (('t_2', 't_1'), 1), (('x_1', 'x_1'), 1))] = None
    diagram_type_dict[((('t_1p', 't_2p'), 1), (('t_2p', 'x_1'), 1), (('x_1', 't_1p'), 1))] = 'Type2'
    diagram_type_dict[((('t_1p', 't_2p'), 1), (('t_2p', 't_1p'), 1), (('x_1', 'x_1'), 1))] = None
    exprs = [
            mk_pi_p("t_1", True) * mk_vec_mu("u", "u", "x_1", 3) * mk_pi_p("t_2"),
            mk_k_p("t_1", True) * mk_vec_mu("u", "u", "x_1", 3) * mk_k_p("t_2"),
            mk_k_m("t_1", True) * mk_vec_mu("s", "s", "x_1", 3) * mk_k_m("t_2"),
            mk_pi_p("t_1p", True) * mk_vec_mu("u", "u", "x_1", 3) * mk_pi_p("t_2p"),
            mk_k_p("t_1p", True) * mk_vec_mu("u", "u", "x_1", 3) * mk_k_p("t_2p"),
            mk_k_m("t_1p", True) * mk_vec_mu("s", "s", "x_1", 3) * mk_k_m("t_2p"),
            ]
    cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True, diagram_type_dict = diagram_type_dict)
    print("meson_jt")
    print(display_cexpr(cexpr))

def test_meson_jj():
    diagram_type_dict = dict()
    diagram_type_dict[((('x_1', 'x_2'), 1), (('x_2', 'x_1'), 1))] = 'Type0'
    diagram_type_dict[((('x_1', 'x_1'), 1), (('x_2', 'x_2'), 1))] = None
    diagram_type_dict[((('t_1', 'x_1'), 1), (('t_2', 'x_2'), 1), (('x_1', 't_1'), 1), (('x_2', 't_2'), 1))] = 'Type1'
    diagram_type_dict[((('t_1', 'x_1'), 1), (('t_2', 'x_2'), 1), (('x_1', 't_2'), 1), (('x_2', 't_1'), 1))] = 'Type2'
    diagram_type_dict[((('t_1', 't_2'), 1), (('t_2', 'x_1'), 1), (('x_1', 'x_2'), 1), (('x_2', 't_1'), 1))] = 'Type3'
    diagram_type_dict[((('t_1', 't_2'), 1), (('t_2', 't_1'), 1), (('x_1', 'x_2'), 1), (('x_2', 'x_1'), 1))] = 'Type4'
    diagram_type_dict[((('t_1', 't_2'), 1), (('t_2', 'x_1'), 1), (('x_1', 't_1'), 1), (('x_2', 'x_2'), 1))] = None
    diagram_type_dict[((('t_1', 't_2'), 1), (('t_2', 't_1'), 1), (('x_1', 'x_1'), 1), (('x_2', 'x_2'), 1))] = None
    exprs = [
            mk_pi_0("t_1", True) * mk_j_mu("x_1", "mu") * mk_j_mu("x_2", "nu") * mk_pi_0("t_2"),
            mk_pi_p("t_1", True) * mk_j_mu("x_1", "mu") * mk_j_mu("x_2", "nu") * mk_pi_p("t_2"),
            mk_k_0("t_1", True) * mk_j_mu("x_1", "mu") * mk_j_mu("x_2", "nu") * mk_k_0("t_2"),
            mk_k_p("t_1", True) * mk_j_mu("x_1", "mu") * mk_j_mu("x_2", "nu") * mk_k_p("t_2"),
            mk_j_mu("x_1", "mu") * mk_j_mu("x_2", "nu"),
            ]
    cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True, diagram_type_dict = diagram_type_dict)
    print("meson_jj")
    print(display_cexpr(cexpr))

def test_meson_fj():
    diagram_type_dict = dict()
    diagram_type_dict[((('t', 'x_1'), 1), (('x_1', 'x_2'), 1), (('x_2', 't'), 1))] = 'Type1'
    diagram_type_dict[((('t', 'x_1'), 1), (('x_1', 't'), 1), (('x_2', 'x_2'), 1))] = None
    exprs = [
            mk_j5pi_mu("x_1", "mu") * mk_j_mu("x_2", "nu") * mk_pi_p("t"),
            mk_j5k_mu("x_1", "mu") * mk_j_mu("x_2", "nu") * mk_k_p("t"),
            mk_jpi_mu("x_1", "mu") * mk_j_mu("x_2", "nu") * mk_pi_p("t"),
            mk_jk_mu("x_1", "mu") * mk_j_mu("x_2", "nu") * mk_k_p("t"),
            ]
    cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True, diagram_type_dict = diagram_type_dict)
    print("meson_fj")
    print(display_cexpr(cexpr))

test_meson_corr()
print()
test_meson_f_corr()
print()
test_meson_m()
print()
test_meson_jt()
print()
test_meson_jj()
print()
test_meson_fj()
