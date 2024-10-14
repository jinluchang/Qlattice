#!/usr/bin/env python3

import qlat as q

from auto_contractor.operators import *
from auto_contractor.eval import *

import sys

is_cython = False

def mk_bk_vv_aa(p : str):
    s = 0
    for mu in range(4):
        v1 = mk_vec_mu("s", "d", p, mu) - mk_vec5_mu("s", "d", p, mu)
        v2 = mk_vec_mu("s", "d", p, mu) - mk_vec5_mu("s", "d", p, mu)
        s = s + v1 * v2
    return s + f"(sbar gmu (1-g5) d)(sbar gmu (1-g5) d)({p})"

def mk_bpi_vv_aa(p : str):
    s = 0
    for mu in range(4):
        v1 = mk_vec_mu("u", "d", p, mu) - mk_vec5_mu("u", "d", p, mu)
        v2 = mk_vec_mu("u", "d", p, mu) - mk_vec5_mu("u", "d", p, mu)
        s = s + v1 * v2
    return s + f"(ubar gmu (1-g5) d)(ubar gmu (1-g5) d)({p})"

def mk_bkp_vv_aa(p : str):
    s = 0
    for mu in range(4):
        v1 = mk_vec_mu("s", "d", p, mu) + mk_vec5_mu("s", "d", p, mu)
        v2 = mk_vec_mu("s", "d", p, mu) + mk_vec5_mu("s", "d", p, mu)
        s = s + v1 * v2
    return s + f"(sbar gmu (1+g5) d)(sbar gmu (1+g5) d)({p})"

def mk_bpip_vv_aa(p : str):
    s = 0
    for mu in range(4):
        v1 = mk_vec_mu("u", "d", p, mu) + mk_vec5_mu("u", "d", p, mu)
        v2 = mk_vec_mu("u", "d", p, mu) + mk_vec5_mu("u", "d", p, mu)
        s = s + v1 * v2
    return s + f"(ubar gmu (1+g5) d)(ubar gmu (1+g5) d)({p})"

def mk_bkpi1_vv_aa(p : str):
    s = 0
    for mu in range(4):
        v1 = mk_vec_mu("d", "u", p, mu) - mk_vec5_mu("d", "u", p, mu)
        v2 = mk_vec_mu("u'", "s", p, mu) - mk_vec5_mu("u'", "s", p, mu)
        s = s + v1 * v2
    return s + f"(dbar gmu (1-g5) u)(u'bar gmu (1-g5) s)({p})"

def mk_bkpi2_vv_aa(p : str):
    s = 0
    for mu in range(4):
        v1 = mk_vec_mu("u'", "u", p, mu) - mk_vec5_mu("u'", "u", p, mu)
        v2 = mk_vec_mu("d", "s", p, mu) - mk_vec5_mu("d", "s", p, mu)
        s = s + v1 * v2
    return s + f"(u'bar gmu (1-g5) u)(dbar gmu (1-g5) s)({p})"

def mk_bkpi1p_vv_aa(p : str):
    s = 0
    for mu in range(4):
        v1 = mk_vec_mu("d", "u", p, mu) + mk_vec5_mu("d", "u", p, mu)
        v2 = mk_vec_mu("u'", "s", p, mu) + mk_vec5_mu("u'", "s", p, mu)
        s = s + v1 * v2
    return s + f"(dbar gmu (1+g5) u)(u'bar gmu (1+g5) s)({p})"

def mk_bkpi2p_vv_aa(p : str):
    s = 0
    for mu in range(4):
        v1 = mk_vec_mu("u'", "u", p, mu) + mk_vec5_mu("u'", "u", p, mu)
        v2 = mk_vec_mu("d", "s", p, mu) + mk_vec5_mu("d", "s", p, mu)
        s = s + v1 * v2
    return s + f"(u'bar gmu (1+g5) u)(dbar gmu (1+g5) s)({p})"

def mk_bkpi3_vv_aa(p : str):
    s = 0
    for mu in range(4):
        v1 = mk_vec_mu("s", "u", p, mu) - mk_vec5_mu("s", "u", p, mu)
        v2 = mk_vec_mu("u'", "d", p, mu) - mk_vec5_mu("u'", "d", p, mu)
        s = s + v1 * v2
    return s + f"(sbar gmu (1-g5) u)(u'bar gmu (1-g5) d)({p})"

def mk_bkpi4_vv_aa(p : str):
    s = 0
    for mu in range(4):
        v1 = mk_vec_mu("u'", "u", p, mu) - mk_vec5_mu("u'", "u", p, mu)
        v2 = mk_vec_mu("s", "d", p, mu) - mk_vec5_mu("s", "d", p, mu)
        s = s + v1 * v2
    return s + f"(u'bar gmu (1-g5) u)(sbar gmu (1-g5) d)({p})"

def mk_bkpi3p_vv_aa(p : str):
    s = 0
    for mu in range(4):
        v1 = mk_vec_mu("s", "u", p, mu) + mk_vec5_mu("s", "u", p, mu)
        v2 = mk_vec_mu("u'", "d", p, mu) + mk_vec5_mu("u'", "d", p, mu)
        s = s + v1 * v2
    return s + f"(sbar gmu (1+g5) u)(u'bar gmu (1+g5) d)({p})"

def mk_bkpi4p_vv_aa(p : str):
    s = 0
    for mu in range(4):
        v1 = mk_vec_mu("u'", "u", p, mu) + mk_vec5_mu("u'", "u", p, mu)
        v2 = mk_vec_mu("s", "d", p, mu) + mk_vec5_mu("s", "d", p, mu)
        s = s + v1 * v2
    return s + f"(u'bar gmu (1+g5) u)(sbar gmu (1+g5) d)({p})"

@q.timer
def get_cexpr_meson_bk_bpi_corr():
    fn_base = f"cache/auto_contract_cexpr/get_cexpr_meson_bk_bpi_corr"
    def calc_cexpr():
        m_ds_2 = mk_meson("d", "s", "t_2")
        m_du_2 = mk_meson("d", "u", "t_2")
        m_uup_2 = mk_meson("u", "u'", "t_2")
        m_ds_1 = mk_meson("d", "s", "t_1")
        m_du_1 = mk_meson("d", "u", "t_1")
        m_sd_1 = mk_meson("s", "d", "t_1")
        bk = m_ds_2 * mk_bk_vv_aa("x") * m_ds_1
        bpi = m_du_2 * mk_bpi_vv_aa("x") * m_du_1
        bkp = m_ds_2 * mk_bkp_vv_aa("x") * m_ds_1
        bpip = m_du_2 * mk_bpip_vv_aa("x") * m_du_1
        bkpi1 = m_uup_2 * mk_bkpi1_vv_aa("x") * m_sd_1
        bkpi2 = m_uup_2 * mk_bkpi2_vv_aa("x") * m_sd_1
        bkpi3 = m_uup_2 * mk_bkpi3_vv_aa("x") * m_ds_1
        bkpi4 = m_uup_2 * mk_bkpi4_vv_aa("x") * m_ds_1
        bkpi1p = m_uup_2 * mk_bkpi1p_vv_aa("x") * m_sd_1
        bkpi2p = m_uup_2 * mk_bkpi2p_vv_aa("x") * m_sd_1
        bkpi3p = m_uup_2 * mk_bkpi3p_vv_aa("x") * m_ds_1
        bkpi4p = m_uup_2 * mk_bkpi4p_vv_aa("x") * m_ds_1
        exprs = [
                bk, bkp,
                bpi, bpip,
                bkpi1, bkpi1p,
                bkpi2, bkpi2p,
                bkpi3, bkpi3p,
                bkpi4, bkpi4p,
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit=True)
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

@q.timer
def get_cexpr_meson_jt_zv():
    fn_base = f"cache/auto_contract_cexpr/get_cexpr_meson_jt_zv"
    def calc_cexpr():
        diagram_type_dict = dict()
        diagram_type_dict[((('t_1', 't_2'), 1), (('t_2', 'x_1'), 1), (('x_1', 't_1'), 1))] = 'Type1'
        diagram_type_dict[((('t_1', 't_2'), 1), (('t_2', 't_1'), 1), (('x_1', 'x_1'), 1))] = None
        diagram_type_dict[((('t_1p', 't_2p'), 1), (('t_2p', 'x_1'), 1), (('x_1', 't_1p'), 1))] = 'Type2'
        diagram_type_dict[((('t_1p', 't_2p'), 1), (('t_2p', 't_1p'), 1), (('x_1', 'x_1'), 1))] = None
        exprs = [
                mk_sym(1)/2 * (
                    mk_pi_p("t_1", True) * mk_vec_mu("u", "u", "x_1", 3) * mk_pi_p("t_2")
                    - mk_pi_m("t_1", True) * mk_vec_mu("u", "u", "x_1", 3) * mk_pi_m("t_2")
                    ),
                mk_sym(1)/2 * (
                    mk_k_p("t_1", True) * mk_vec_mu("u", "u", "x_1", 3) * mk_k_p("t_2")
                    - mk_k_m("t_1", True) * mk_vec_mu("u", "u", "x_1", 3) * mk_k_m("t_2")
                    ),
                mk_sym(1)/2 * (
                    mk_k_m("t_1", True) * mk_vec_mu("s", "s", "x_1", 3) * mk_k_m("t_2")
                    - mk_k_p("t_1", True) * mk_vec_mu("s", "s", "x_1", 3) * mk_k_p("t_2")
                    ),
                mk_sym(1)/2 * (
                    mk_pi_p("t_1p", True) * mk_vec_mu("u", "u", "x_1", 3) * mk_pi_p("t_2p")
                    - mk_pi_m("t_1p", True) * mk_vec_mu("u", "u", "x_1", 3) * mk_pi_m("t_2p")
                    ),
                mk_sym(1)/2 * (
                    mk_k_p("t_1p", True) * mk_vec_mu("u", "u", "x_1", 3) * mk_k_p("t_2p")
                    - mk_k_m("t_1p", True) * mk_vec_mu("u", "u", "x_1", 3) * mk_k_m("t_2p")
                    ),
                mk_sym(1)/2 * (
                    mk_k_m("t_1p", True) * mk_vec_mu("s", "s", "x_1", 3) * mk_k_m("t_2p")
                    - mk_k_p("t_1p", True) * mk_vec_mu("s", "s", "x_1", 3) * mk_k_p("t_2p")
                    ),
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit=True, diagram_type_dict=diagram_type_dict)
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

@q.timer
def get_cexpr_meson_jj_mm():
    fn_base = f"cache/auto_contract_cexpr/get_cexpr_meson_jj_mm"
    def calc_cexpr():
        jj_op = sum([
            mk_j_mu("x_1", mu) * mk_j_mu("x_2", mu)
            for mu in range(4)
            ])
        exprs = [ jj_op * mk_pi_0("t_1", True) * mk_pi_0("t_2")
                + f"pi0 j_mu j_mu pi0",
                mk_sym(1)/2 * jj_op
                * (mk_pi_p("t_1", True) * mk_pi_p("t_2") + mk_pi_m("t_1", True) * mk_pi_m("t_2"))
                + f"pi+ j_mu j_mu pi+",
                mk_sym(1)/2 * jj_op
                * (mk_k_0("t_1", True) * mk_k_0("t_2") + mk_k_0_bar("t_1", True) * mk_k_0_bar("t_2"))
                + f"K0 j_mu j_mu K0",
                mk_sym(1)/2 * jj_op
                * (mk_k_p("t_1", True) * mk_k_p("t_2") + mk_k_m("t_1", True) * mk_k_m("t_2"))
                + f"K+ j_mu j_mu K+",
                jj_op
                + f"j_mu j_mu",
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit=True)
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

@q.timer
def get_cexpr_meson_jj_xx():
    fn_base = f"cache/auto_contract_cexpr/get_cexpr_meson_jj_xx"
    def calc_cexpr():
        jj_op = sum([
            mk_fac(f"rel_mod_sym(x_1[1][{mu}]-x_2[1][{mu}],size[{mu}])")
            * mk_fac(f"rel_mod_sym(x_1[1][{nu}]-x_2[1][{nu}],size[{nu}])")
            * mk_j_mu("x_1", mu) * mk_j_mu("x_2", nu)
            for mu in range(3) for nu in range(3)
            ])
        exprs = [ jj_op * mk_pi_0("t_1", True) * mk_pi_0("t_2")
                + f"x[a] x[b] pi0 j_a j_b pi0",
                mk_sym(1)/2 * jj_op
                * (mk_pi_p("t_1", True) * mk_pi_p("t_2") + mk_pi_m("t_1", True) * mk_pi_m("t_2"))
                + f"x[a] x[b] pi+ j_a j_b pi+",
                mk_sym(1)/2 * jj_op
                * (mk_k_0("t_1", True) * mk_k_0("t_2") + mk_k_0_bar("t_1", True) * mk_k_0_bar("t_2"))
                + f"x[a] x[b] K0 j_a j_b K0",
                mk_sym(1)/2 * jj_op
                * (mk_k_p("t_1", True) * mk_k_p("t_2") + mk_k_m("t_1", True) * mk_k_m("t_2"))
                + f"x[a] x[b] K+ j_a j_b K+",
                jj_op
                + f"x[a] x[b] j_a j_b",
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit=True)
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

@q.timer
def get_cexpr_meson_jj_mm_types():
    fn_base = f"cache/auto_contract_cexpr/get_cexpr_meson_jj_mm_types"
    def calc_cexpr():
        diagram_type_dict = dict()
        diagram_type_dict[((('t_1', 't_2'), 1), (('t_2', 'x_1'), 1), (('x_1', 't_1'), 1), (('x_2', 'x_2'), 1))] = None
        diagram_type_dict[((('t_1', 't_2'), 1), (('t_2', 't_1'), 1), (('x_1', 'x_1'), 1), (('x_2', 'x_2'), 1))] = None
        diagram_type_dict[((('t_1', 'x_1'), 1), (('t_2', 'x_2'), 1), (('x_1', 't_1'), 1), (('x_2', 't_2'), 1))] = 'Type1'
        diagram_type_dict[((('t_1', 'x_1'), 1), (('t_2', 'x_2'), 1), (('x_1', 't_2'), 1), (('x_2', 't_1'), 1))] = 'Type2'
        diagram_type_dict[((('t_1', 't_2'), 1), (('t_2', 'x_1'), 1), (('x_1', 'x_2'), 1), (('x_2', 't_1'), 1))] = 'Type3'
        diagram_type_dict[((('t_1', 't_2'), 1), (('t_2', 't_1'), 1), (('x_1', 'x_2'), 1), (('x_2', 'x_1'), 1))] = 'Type4'
        diagram_type_dict[((('x_1', 'x_1'), 1), (('x_2', 'x_2'), 1))] = None
        diagram_type_dict[((('x_1', 'x_2'), 1), (('x_2', 'x_1'), 1))] = 'Type4'
        jj_op = sum([
            mk_j_mu("x_1", mu) * mk_j_mu("x_2", mu)
            for mu in range(4)
            ])
        exprs = [ jj_op * mk_pi_0("t_1", True) * mk_pi_0("t_2")
                + f"pi0 j_mu j_mu pi0",
                mk_sym(1)/2 * jj_op
                * (mk_pi_p("t_1", True) * mk_pi_p("t_2") + mk_pi_m("t_1", True) * mk_pi_m("t_2"))
                + f"pi+ j_mu j_mu pi+",
                mk_sym(1)/2 * jj_op
                * (mk_k_0("t_1", True) * mk_k_0("t_2") + mk_k_0_bar("t_1", True) * mk_k_0_bar("t_2"))
                + f"K0 j_mu j_mu K0",
                mk_sym(1)/2 * jj_op
                * (mk_k_p("t_1", True) * mk_k_p("t_2") + mk_k_m("t_1", True) * mk_k_m("t_2"))
                + f"K+ j_mu j_mu K+",
                jj_op
                + f"j_mu j_mu",
                ]
        exprs = contract_simplify(*exprs, is_isospin_symmetric_limit=True, diagram_type_dict=diagram_type_dict)
        typed_exprs = []
        for expr in exprs:
            typed_exprs.append((expr, None, [ 'Type2', 'Type3', ], 'Type1', 'Type2', 'Type3', 'Type4',))
        cexpr = contract_simplify_compile(*typed_exprs, is_isospin_symmetric_limit=True, diagram_type_dict=diagram_type_dict)
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

def get_all_cexpr():
    cexprs = [
            lambda : get_cexpr_meson_bk_bpi_corr(),
            lambda : get_cexpr_meson_jt_zv(),
            lambda : get_cexpr_meson_jj_mm(),
            lambda : get_cexpr_meson_jj_xx(),
            lambda : get_cexpr_meson_jj_mm_types(),
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
