#!/usr/bin/env python3

from auto_contractor.operators import *

import functools
import math
import os
import sys

from jobs import *
from load_data import *
from params import *

# ----

load_path_list[:] = [
        "results",
        "qcddata",
        "qcddata-1",
        "qcddata-2",
        "qcddata-3",
        "qcddata-4",
        "qcddata-5",
        # "../qcddata",
        # os.path.join(os.getenv("HOME"), "qcddata"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-gf-gt/results"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-sel/results"),
        # os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-fsel-self-loop/results"),
        # os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-selected-data/results"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-wsrc-prop/results"),
        # os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-psrc-prop/results"),
        # os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-smear-prop/results"),
        # "/sdcc/u/jluchang/qcdqedta/luchang/data-gen/fill-wsnk-prop/results",
        # "/sdcc/u/jluchang/qcdqedta/summit-oakforest-data-cache",
        ]

# ----

@q.timer
def get_cexpr_meson_corr():
    def calc_cexpr():
        exprs = [
                mk_pi_p("t_2", True)    * mk_pi_p("t_1")    + "pi^dag * pi   ",
                mk_k_p("t_2", True)     * mk_k_p("t_1")     + "k^dag  * k    ",
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True)
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, f"cache/auto_contract_cexpr/get_cexpr_meson_corr")

@q.timer
def get_cexpr_meson_f_corr():
    def calc_cexpr():
        exprs = [
                mk_j5pi_mu("x_2", 3)    * mk_pi_p("t_1")    + "a_pi   * pi   ",
                mk_j5k_mu("x_2", 3)     * mk_k_p("t_1")     + "a_k    * k    ",
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True)
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, f"cache/auto_contract_cexpr/get_cexpr_meson_f_corr")

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
    def calc_cexpr():
        t_1, t_2, x = ['t_1', 't_2', 'x']
        terms = [
          tr(gamma_x*gamma_5*S_l(x,t_1)*gamma_5*S_s(t_1,x)*gamma_x*gamma_5*S_l(x,t_2)*gamma_5*S_s(t_2,x)), # term_ADT00_0000
          tr(gamma_x*gamma_5*S_l(x,t_1)*gamma_5*S_s(t_1,x)*gamma_x*S_l(x,t_2)*gamma_5*S_s(t_2,x)), # term_ADT00_0001
          tr(gamma_x*gamma_5*S_l(x,t_1)*gamma_5*S_s(t_1,x))*tr(gamma_x*gamma_5*S_l(x,t_2)*gamma_5*S_s(t_2,x)), # term_ADT00_0002
          tr(gamma_x*gamma_5*S_l(x,t_1)*gamma_5*S_s(t_1,x))*tr(gamma_x*S_l(x,t_2)*gamma_5*S_s(t_2,x)), # term_ADT00_0003
          tr(gamma_x*gamma_5*S_l(x,t_2)*gamma_5*S_s(t_2,x)*gamma_x*S_l(x,t_1)*gamma_5*S_s(t_1,x)), # term_ADT00_0004
          tr(gamma_x*gamma_5*S_l(x,t_2)*gamma_5*S_s(t_2,x))*tr(gamma_x*S_l(x,t_1)*gamma_5*S_s(t_1,x)), # term_ADT00_0005
          tr(gamma_x*S_l(x,t_1)*gamma_5*S_s(t_1,x)*gamma_x*S_l(x,t_2)*gamma_5*S_s(t_2,x)), # term_ADT00_0006
          tr(gamma_x*S_l(x,t_1)*gamma_5*S_s(t_1,x))*tr(gamma_x*S_l(x,t_2)*gamma_5*S_s(t_2,x)), # term_ADT00_0007
          tr(gamma_y*gamma_5*S_l(x,t_1)*gamma_5*S_s(t_1,x)*gamma_y*gamma_5*S_l(x,t_2)*gamma_5*S_s(t_2,x)), # term_ADT00_0008
          tr(gamma_y*gamma_5*S_l(x,t_1)*gamma_5*S_s(t_1,x)*gamma_y*S_l(x,t_2)*gamma_5*S_s(t_2,x)), # term_ADT00_0009
          tr(gamma_y*gamma_5*S_l(x,t_1)*gamma_5*S_s(t_1,x))*tr(gamma_y*gamma_5*S_l(x,t_2)*gamma_5*S_s(t_2,x)), # term_ADT00_0010
          tr(gamma_y*gamma_5*S_l(x,t_1)*gamma_5*S_s(t_1,x))*tr(gamma_y*S_l(x,t_2)*gamma_5*S_s(t_2,x)), # term_ADT00_0011
          tr(gamma_y*gamma_5*S_l(x,t_2)*gamma_5*S_s(t_2,x)*gamma_y*S_l(x,t_1)*gamma_5*S_s(t_1,x)), # term_ADT00_0012
          tr(gamma_y*gamma_5*S_l(x,t_2)*gamma_5*S_s(t_2,x))*tr(gamma_y*S_l(x,t_1)*gamma_5*S_s(t_1,x)), # term_ADT00_0013
          tr(gamma_y*S_l(x,t_1)*gamma_5*S_s(t_1,x)*gamma_y*S_l(x,t_2)*gamma_5*S_s(t_2,x)), # term_ADT00_0014
          tr(gamma_y*S_l(x,t_1)*gamma_5*S_s(t_1,x))*tr(gamma_y*S_l(x,t_2)*gamma_5*S_s(t_2,x)), # term_ADT00_0015
          tr(gamma_z*gamma_5*S_l(x,t_1)*gamma_5*S_s(t_1,x)*gamma_z*gamma_5*S_l(x,t_2)*gamma_5*S_s(t_2,x)), # term_ADT00_0016
          tr(gamma_z*gamma_5*S_l(x,t_1)*gamma_5*S_s(t_1,x)*gamma_z*S_l(x,t_2)*gamma_5*S_s(t_2,x)), # term_ADT00_0017
          tr(gamma_z*gamma_5*S_l(x,t_1)*gamma_5*S_s(t_1,x))*tr(gamma_z*gamma_5*S_l(x,t_2)*gamma_5*S_s(t_2,x)), # term_ADT00_0018
          tr(gamma_z*gamma_5*S_l(x,t_1)*gamma_5*S_s(t_1,x))*tr(gamma_z*S_l(x,t_2)*gamma_5*S_s(t_2,x)), # term_ADT00_0019
          tr(gamma_z*gamma_5*S_l(x,t_2)*gamma_5*S_s(t_2,x)*gamma_z*S_l(x,t_1)*gamma_5*S_s(t_1,x)), # term_ADT00_0020
          tr(gamma_z*gamma_5*S_l(x,t_2)*gamma_5*S_s(t_2,x))*tr(gamma_z*S_l(x,t_1)*gamma_5*S_s(t_1,x)), # term_ADT00_0021
          tr(gamma_z*S_l(x,t_1)*gamma_5*S_s(t_1,x)*gamma_z*S_l(x,t_2)*gamma_5*S_s(t_2,x)), # term_ADT00_0022
          tr(gamma_z*S_l(x,t_1)*gamma_5*S_s(t_1,x))*tr(gamma_z*S_l(x,t_2)*gamma_5*S_s(t_2,x)), # term_ADT00_0023
          tr(gamma_t*gamma_5*S_l(x,t_1)*gamma_5*S_s(t_1,x)*gamma_t*gamma_5*S_l(x,t_2)*gamma_5*S_s(t_2,x)), # term_ADT00_0024
          tr(gamma_t*gamma_5*S_l(x,t_1)*gamma_5*S_s(t_1,x)*gamma_t*S_l(x,t_2)*gamma_5*S_s(t_2,x)), # term_ADT00_0025
          tr(gamma_t*gamma_5*S_l(x,t_1)*gamma_5*S_s(t_1,x))*tr(gamma_t*gamma_5*S_l(x,t_2)*gamma_5*S_s(t_2,x)), # term_ADT00_0026
          tr(gamma_t*gamma_5*S_l(x,t_1)*gamma_5*S_s(t_1,x))*tr(gamma_t*S_l(x,t_2)*gamma_5*S_s(t_2,x)), # term_ADT00_0027
          tr(gamma_t*gamma_5*S_l(x,t_2)*gamma_5*S_s(t_2,x)*gamma_t*S_l(x,t_1)*gamma_5*S_s(t_1,x)), # term_ADT00_0028
          tr(gamma_t*gamma_5*S_l(x,t_2)*gamma_5*S_s(t_2,x))*tr(gamma_t*S_l(x,t_1)*gamma_5*S_s(t_1,x)), # term_ADT00_0029
          tr(gamma_t*S_l(x,t_1)*gamma_5*S_s(t_1,x)*gamma_t*S_l(x,t_2)*gamma_5*S_s(t_2,x)), # term_ADT00_0030
          tr(gamma_t*S_l(x,t_1)*gamma_5*S_s(t_1,x))*tr(gamma_t*S_l(x,t_2)*gamma_5*S_s(t_2,x)), # term_ADT00_0031
          tr(gamma_x*gamma_5*S_l(x,t_1)*gamma_5*S_l(t_1,x)*gamma_x*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0032
          tr(gamma_x*gamma_5*S_l(x,t_1)*gamma_5*S_l(t_1,x)*gamma_x*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0033
          tr(gamma_x*gamma_5*S_l(x,t_1)*gamma_5*S_l(t_1,x))*tr(gamma_x*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0034
          tr(gamma_x*gamma_5*S_l(x,t_1)*gamma_5*S_l(t_1,x))*tr(gamma_x*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0035
          tr(gamma_x*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x)*gamma_x*S_l(x,t_1)*gamma_5*S_l(t_1,x)), # term_ADT00_0036
          tr(gamma_x*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x))*tr(gamma_x*S_l(x,t_1)*gamma_5*S_l(t_1,x)), # term_ADT00_0037
          tr(gamma_x*S_l(x,t_1)*gamma_5*S_l(t_1,x)*gamma_x*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0038
          tr(gamma_x*S_l(x,t_1)*gamma_5*S_l(t_1,x))*tr(gamma_x*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0039
          tr(gamma_y*gamma_5*S_l(x,t_1)*gamma_5*S_l(t_1,x)*gamma_y*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0040
          tr(gamma_y*gamma_5*S_l(x,t_1)*gamma_5*S_l(t_1,x)*gamma_y*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0041
          tr(gamma_y*gamma_5*S_l(x,t_1)*gamma_5*S_l(t_1,x))*tr(gamma_y*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0042
          tr(gamma_y*gamma_5*S_l(x,t_1)*gamma_5*S_l(t_1,x))*tr(gamma_y*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0043
          tr(gamma_y*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x)*gamma_y*S_l(x,t_1)*gamma_5*S_l(t_1,x)), # term_ADT00_0044
          tr(gamma_y*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x))*tr(gamma_y*S_l(x,t_1)*gamma_5*S_l(t_1,x)), # term_ADT00_0045
          tr(gamma_y*S_l(x,t_1)*gamma_5*S_l(t_1,x)*gamma_y*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0046
          tr(gamma_y*S_l(x,t_1)*gamma_5*S_l(t_1,x))*tr(gamma_y*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0047
          tr(gamma_z*gamma_5*S_l(x,t_1)*gamma_5*S_l(t_1,x)*gamma_z*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0048
          tr(gamma_z*gamma_5*S_l(x,t_1)*gamma_5*S_l(t_1,x)*gamma_z*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0049
          tr(gamma_z*gamma_5*S_l(x,t_1)*gamma_5*S_l(t_1,x))*tr(gamma_z*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0050
          tr(gamma_z*gamma_5*S_l(x,t_1)*gamma_5*S_l(t_1,x))*tr(gamma_z*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0051
          tr(gamma_z*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x)*gamma_z*S_l(x,t_1)*gamma_5*S_l(t_1,x)), # term_ADT00_0052
          tr(gamma_z*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x))*tr(gamma_z*S_l(x,t_1)*gamma_5*S_l(t_1,x)), # term_ADT00_0053
          tr(gamma_z*S_l(x,t_1)*gamma_5*S_l(t_1,x)*gamma_z*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0054
          tr(gamma_z*S_l(x,t_1)*gamma_5*S_l(t_1,x))*tr(gamma_z*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0055
          tr(gamma_t*gamma_5*S_l(x,t_1)*gamma_5*S_l(t_1,x)*gamma_t*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0056
          tr(gamma_t*gamma_5*S_l(x,t_1)*gamma_5*S_l(t_1,x)*gamma_t*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0057
          tr(gamma_t*gamma_5*S_l(x,t_1)*gamma_5*S_l(t_1,x))*tr(gamma_t*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0058
          tr(gamma_t*gamma_5*S_l(x,t_1)*gamma_5*S_l(t_1,x))*tr(gamma_t*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0059
          tr(gamma_t*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x)*gamma_t*S_l(x,t_1)*gamma_5*S_l(t_1,x)), # term_ADT00_0060
          tr(gamma_t*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x))*tr(gamma_t*S_l(x,t_1)*gamma_5*S_l(t_1,x)), # term_ADT00_0061
          tr(gamma_t*S_l(x,t_1)*gamma_5*S_l(t_1,x)*gamma_t*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0062
          tr(gamma_t*S_l(x,t_1)*gamma_5*S_l(t_1,x))*tr(gamma_t*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0063
          tr(gamma_x*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x)*gamma_x*gamma_5*S_s(x,t_1)*gamma_5*S_l(t_1,x)), # term_ADT00_0064
          tr(gamma_x*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x)*gamma_x*S_s(x,t_1)*gamma_5*S_l(t_1,x)), # term_ADT00_0065
          tr(gamma_x*gamma_5*S_s(x,t_1)*gamma_5*S_l(t_1,x)*gamma_x*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0066
          tr(gamma_x*S_l(x,t_2)*gamma_5*S_l(t_2,x)*gamma_x*S_s(x,t_1)*gamma_5*S_l(t_1,x)), # term_ADT00_0067
          tr(gamma_y*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x)*gamma_y*gamma_5*S_s(x,t_1)*gamma_5*S_l(t_1,x)), # term_ADT00_0068
          tr(gamma_y*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x)*gamma_y*S_s(x,t_1)*gamma_5*S_l(t_1,x)), # term_ADT00_0069
          tr(gamma_y*gamma_5*S_s(x,t_1)*gamma_5*S_l(t_1,x)*gamma_y*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0070
          tr(gamma_y*S_l(x,t_2)*gamma_5*S_l(t_2,x)*gamma_y*S_s(x,t_1)*gamma_5*S_l(t_1,x)), # term_ADT00_0071
          tr(gamma_z*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x)*gamma_z*gamma_5*S_s(x,t_1)*gamma_5*S_l(t_1,x)), # term_ADT00_0072
          tr(gamma_z*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x)*gamma_z*S_s(x,t_1)*gamma_5*S_l(t_1,x)), # term_ADT00_0073
          tr(gamma_z*gamma_5*S_s(x,t_1)*gamma_5*S_l(t_1,x)*gamma_z*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0074
          tr(gamma_z*S_l(x,t_2)*gamma_5*S_l(t_2,x)*gamma_z*S_s(x,t_1)*gamma_5*S_l(t_1,x)), # term_ADT00_0075
          tr(gamma_t*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x)*gamma_t*gamma_5*S_s(x,t_1)*gamma_5*S_l(t_1,x)), # term_ADT00_0076
          tr(gamma_t*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x)*gamma_t*S_s(x,t_1)*gamma_5*S_l(t_1,x)), # term_ADT00_0077
          tr(gamma_t*gamma_5*S_s(x,t_1)*gamma_5*S_l(t_1,x)*gamma_t*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0078
          tr(gamma_t*S_l(x,t_2)*gamma_5*S_l(t_2,x)*gamma_t*S_s(x,t_1)*gamma_5*S_l(t_1,x)), # term_ADT00_0079
          tr(gamma_x*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x))*tr(gamma_x*gamma_5*S_s(x,t_1)*gamma_5*S_l(t_1,x)), # term_ADT00_0080
          tr(gamma_x*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x))*tr(gamma_x*S_s(x,t_1)*gamma_5*S_l(t_1,x)), # term_ADT00_0081
          tr(gamma_x*gamma_5*S_s(x,t_1)*gamma_5*S_l(t_1,x))*tr(gamma_x*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0082
          tr(gamma_x*S_l(x,t_2)*gamma_5*S_l(t_2,x))*tr(gamma_x*S_s(x,t_1)*gamma_5*S_l(t_1,x)), # term_ADT00_0083
          tr(gamma_y*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x))*tr(gamma_y*gamma_5*S_s(x,t_1)*gamma_5*S_l(t_1,x)), # term_ADT00_0084
          tr(gamma_y*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x))*tr(gamma_y*S_s(x,t_1)*gamma_5*S_l(t_1,x)), # term_ADT00_0085
          tr(gamma_y*gamma_5*S_s(x,t_1)*gamma_5*S_l(t_1,x))*tr(gamma_y*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0086
          tr(gamma_y*S_l(x,t_2)*gamma_5*S_l(t_2,x))*tr(gamma_y*S_s(x,t_1)*gamma_5*S_l(t_1,x)), # term_ADT00_0087
          tr(gamma_z*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x))*tr(gamma_z*gamma_5*S_s(x,t_1)*gamma_5*S_l(t_1,x)), # term_ADT00_0088
          tr(gamma_z*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x))*tr(gamma_z*S_s(x,t_1)*gamma_5*S_l(t_1,x)), # term_ADT00_0089
          tr(gamma_z*gamma_5*S_s(x,t_1)*gamma_5*S_l(t_1,x))*tr(gamma_z*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0090
          tr(gamma_z*S_l(x,t_2)*gamma_5*S_l(t_2,x))*tr(gamma_z*S_s(x,t_1)*gamma_5*S_l(t_1,x)), # term_ADT00_0091
          tr(gamma_t*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x))*tr(gamma_t*gamma_5*S_s(x,t_1)*gamma_5*S_l(t_1,x)), # term_ADT00_0092
          tr(gamma_t*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x))*tr(gamma_t*S_s(x,t_1)*gamma_5*S_l(t_1,x)), # term_ADT00_0093
          tr(gamma_t*gamma_5*S_s(x,t_1)*gamma_5*S_l(t_1,x))*tr(gamma_t*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0094
          tr(gamma_t*S_l(x,t_2)*gamma_5*S_l(t_2,x))*tr(gamma_t*S_s(x,t_1)*gamma_5*S_l(t_1,x)), # term_ADT00_0095
          tr(gamma_x*gamma_5*S_l(x,t_1)*gamma_5*S_s(t_1,x)*gamma_x*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0096
          tr(gamma_x*gamma_5*S_l(x,t_1)*gamma_5*S_s(t_1,x)*gamma_x*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0097
          tr(gamma_x*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x)*gamma_x*S_l(x,t_1)*gamma_5*S_s(t_1,x)), # term_ADT00_0098
          tr(gamma_x*S_l(x,t_1)*gamma_5*S_s(t_1,x)*gamma_x*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0099
          tr(gamma_y*gamma_5*S_l(x,t_1)*gamma_5*S_s(t_1,x)*gamma_y*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0100
          tr(gamma_y*gamma_5*S_l(x,t_1)*gamma_5*S_s(t_1,x)*gamma_y*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0101
          tr(gamma_y*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x)*gamma_y*S_l(x,t_1)*gamma_5*S_s(t_1,x)), # term_ADT00_0102
          tr(gamma_y*S_l(x,t_1)*gamma_5*S_s(t_1,x)*gamma_y*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0103
          tr(gamma_z*gamma_5*S_l(x,t_1)*gamma_5*S_s(t_1,x)*gamma_z*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0104
          tr(gamma_z*gamma_5*S_l(x,t_1)*gamma_5*S_s(t_1,x)*gamma_z*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0105
          tr(gamma_z*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x)*gamma_z*S_l(x,t_1)*gamma_5*S_s(t_1,x)), # term_ADT00_0106
          tr(gamma_z*S_l(x,t_1)*gamma_5*S_s(t_1,x)*gamma_z*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0107
          tr(gamma_t*gamma_5*S_l(x,t_1)*gamma_5*S_s(t_1,x)*gamma_t*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0108
          tr(gamma_t*gamma_5*S_l(x,t_1)*gamma_5*S_s(t_1,x)*gamma_t*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0109
          tr(gamma_t*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x)*gamma_t*S_l(x,t_1)*gamma_5*S_s(t_1,x)), # term_ADT00_0110
          tr(gamma_t*S_l(x,t_1)*gamma_5*S_s(t_1,x)*gamma_t*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0111
          tr(gamma_x*gamma_5*S_l(x,t_1)*gamma_5*S_s(t_1,x))*tr(gamma_x*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0112
          tr(gamma_x*gamma_5*S_l(x,t_1)*gamma_5*S_s(t_1,x))*tr(gamma_x*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0113
          tr(gamma_x*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x))*tr(gamma_x*S_l(x,t_1)*gamma_5*S_s(t_1,x)), # term_ADT00_0114
          tr(gamma_x*S_l(x,t_1)*gamma_5*S_s(t_1,x))*tr(gamma_x*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0115
          tr(gamma_y*gamma_5*S_l(x,t_1)*gamma_5*S_s(t_1,x))*tr(gamma_y*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0116
          tr(gamma_y*gamma_5*S_l(x,t_1)*gamma_5*S_s(t_1,x))*tr(gamma_y*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0117
          tr(gamma_y*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x))*tr(gamma_y*S_l(x,t_1)*gamma_5*S_s(t_1,x)), # term_ADT00_0118
          tr(gamma_y*S_l(x,t_1)*gamma_5*S_s(t_1,x))*tr(gamma_y*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0119
          tr(gamma_z*gamma_5*S_l(x,t_1)*gamma_5*S_s(t_1,x))*tr(gamma_z*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0120
          tr(gamma_z*gamma_5*S_l(x,t_1)*gamma_5*S_s(t_1,x))*tr(gamma_z*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0121
          tr(gamma_z*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x))*tr(gamma_z*S_l(x,t_1)*gamma_5*S_s(t_1,x)), # term_ADT00_0122
          tr(gamma_z*S_l(x,t_1)*gamma_5*S_s(t_1,x))*tr(gamma_z*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0123
          tr(gamma_t*gamma_5*S_l(x,t_1)*gamma_5*S_s(t_1,x))*tr(gamma_t*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0124
          tr(gamma_t*gamma_5*S_l(x,t_1)*gamma_5*S_s(t_1,x))*tr(gamma_t*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0125
          tr(gamma_t*gamma_5*S_l(x,t_2)*gamma_5*S_l(t_2,x))*tr(gamma_t*S_l(x,t_1)*gamma_5*S_s(t_1,x)), # term_ADT00_0126
          tr(gamma_t*S_l(x,t_1)*gamma_5*S_s(t_1,x))*tr(gamma_t*S_l(x,t_2)*gamma_5*S_l(t_2,x)), # term_ADT00_0127
        ]
        bk = mk_meson("d", "s", "t_2") * mk_bk_vv_aa("x") * mk_meson("d", "s", "t_1")
        bpi = mk_meson("d", "u", "t_2") * mk_bpi_vv_aa("x") * mk_meson("d", "u", "t_1")
        bkp = mk_meson("d", "s", "t_2") * mk_bkp_vv_aa("x") * mk_meson("d", "s", "t_1")
        bpip = mk_meson("d", "u", "t_2") * mk_bpip_vv_aa("x") * mk_meson("d", "u", "t_1")
        bkpi1 = mk_meson("u", "u'", "t_2") * mk_bkpi1_vv_aa("x") * mk_meson("s", "d", "t_1")
        bkpi2 = mk_meson("u", "u'", "t_2") * mk_bkpi2_vv_aa("x") * mk_meson("s", "d", "t_1")
        bkpi1p = mk_meson("u", "u'", "t_2") * mk_bkpi1p_vv_aa("x") * mk_meson("s", "d", "t_1")
        bkpi2p = mk_meson("u", "u'", "t_2") * mk_bkpi2p_vv_aa("x") * mk_meson("s", "d", "t_1")
        bkpi3 = mk_meson("u", "u'", "t_2") * mk_bkpi3_vv_aa("x") * mk_meson("d", "s", "t_1")
        bkpi4 = mk_meson("u", "u'", "t_2") * mk_bkpi4_vv_aa("x") * mk_meson("d", "s", "t_1")
        bkpi3p = mk_meson("u", "u'", "t_2") * mk_bkpi3p_vv_aa("x") * mk_meson("d", "s", "t_1")
        bkpi4p = mk_meson("u", "u'", "t_2") * mk_bkpi4p_vv_aa("x") * mk_meson("d", "s", "t_1")
        exprs = terms + [
                bk, bkp,
                bpi, bpip,
                bkpi1, bkpi2,
                bkpi1p, bkpi2p,
                bkpi3, bkpi4,
                bkpi3p, bkpi4p,
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True)
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, f"cache/auto_contract_cexpr/get_cexpr_meson_bk_bpi_corr")

@q.timer_verbose
def auto_contract_meson_corr(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_corr.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_corr()
    expr_names = get_cexpr_names(cexpr)
    total_site = rup.get_total_site(job_tag)
    psel = get_psel()
    fsel, fselc = get_fsel()
    xg_fsel_list = np.array(fsel.to_psel_local().to_list())
    xg_psel_list = np.array(psel.to_list())
    geo = q.Geometry(total_site, 1)
    total_volume = geo.total_volume()
    def load_data():
        t_t_list = get_mpi_chunk(
                [ (t_src, t_snk,) for t_snk in range(total_site[3]) for t_src in range(total_site[3]) ],
                rng_state = q.RngState("get_mpi_chunk"))
        for t_src, t_snk in t_t_list:
            t = (t_snk - t_src) % total_site[3]
            pd = {
                    "t_2" : ("wall", t_snk,),
                    "t_1" : ("wall", t_src,),
                    }
            yield pd, t
    @q.timer
    def feval(args):
        pd, t = args
        val = eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop)
        return val, t
    def sum_function(val_list):
        counts = np.zeros(total_site[3], dtype = complex)
        values = np.zeros((len(expr_names), total_site[3],), dtype = complex)
        for val, t in val_list:
            counts[t] += 1
            values[:, t] += val
        return counts, values
    q.timer_fork(0)
    res_count, res_sum = q.glb_sum(
            q.parallel_map_sum(feval, load_data(), sum_function = sum_function, chunksize = 16))
    q.displayln_info("timer_display for auto_contract_meson_corr")
    q.timer_display()
    q.timer_merge()
    res_count *= 1.0 / total_site[3]
    res_sum *= 1.0 / total_site[3]
    assert q.qnorm(res_count - 1.0) < 1e-10
    ld = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t_sep", total_site[3], ],
        ])
    ld.from_numpy(res_sum)
    # q.displayln_info(ld.show())
    ld.save(get_save_path(fn))

@q.timer_verbose
def auto_contract_meson_f_corr(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_f_corr.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_f_corr()
    expr_names = get_cexpr_names(cexpr)
    total_site = rup.get_total_site(job_tag)
    psel = get_psel()
    fsel, fselc = get_fsel()
    xg_fsel_list = np.array(fsel.to_psel_local().to_list())
    xg_psel_list = np.array(psel.to_list())
    geo = q.Geometry(total_site, 1)
    total_volume = geo.total_volume()
    def load_data():
        for t_src in range(total_site[3]):
            q.displayln_info(f"auto_contract_meson_f_corr: {t_src}")
            for xg_snk in xg_fsel_list:
                xg_snk = tuple(xg_snk.tolist())
                t = (xg_snk[3] - t_src) % total_site[3]
                pd = {
                        "x_2" : ("point-snk", xg_snk,),
                        "t_1" : ("wall", t_src,),
                        }
                yield pd, t
    @q.timer
    def feval(args):
        pd, t = args
        val = eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop)
        return val, t
    def sum_function(val_list):
        counts = np.zeros(total_site[3], dtype = complex)
        values = np.zeros((len(expr_names), total_site[3],), dtype = complex)
        for val, t in val_list:
            counts[t] += 1
            values[:, t] += val
        return counts, values
    q.timer_fork(0)
    res_count, res_sum = q.glb_sum(
            q.parallel_map_sum(feval, load_data(), sum_function = sum_function, chunksize = 16))
    q.displayln_info("timer_display for auto_contract_meson_f_corr")
    q.timer_display()
    q.timer_merge()
    res_count *= 1.0 / (total_volume * fsel.prob())
    res_sum *= 1.0 / (total_volume * fsel.prob())
    assert q.qnorm(res_count - 1.0) < 1e-10
    ld = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t_sep", total_site[3], ],
        ])
    ld.from_numpy(res_sum)
    # q.displayln_info(ld.show())
    ld.save(get_save_path(fn))

@q.timer_verbose
def auto_contract_meson_bk_bpi_corr(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_bk_bpi_corr.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_bk_bpi_corr()
    expr_names = get_cexpr_names(cexpr)
    total_site = rup.get_total_site(job_tag)
    psel = get_psel()
    fsel, fselc = get_fsel()
    xg_fsel_list = np.array(fsel.to_psel_local().to_list())
    xg_psel_list = np.array(psel.to_list())
    geo = q.Geometry(total_site, 1)
    total_volume = geo.total_volume()
    tsep_step = max(2, total_site[3] // 16)
    tsep_list = list(range(tsep_step, total_site[3], tsep_step))
    def load_data():
        for t_src in range(total_site[3]):
            for tt_idx, tt in enumerate(tsep_list):
                t_snk = (t_src + tt) % total_site[3]
                q.displayln_info(f"auto_contract_meson_bk_bpi_corr: {t_src} {tt}")
                for xg_snk in xg_fsel_list:
                    xg_snk = tuple(xg_snk.tolist())
                    t = (xg_snk[3] - t_src) % total_site[3]
                    pd = {
                            "t_2" : ("wall", t_snk,),
                            "x" : ("point-snk", xg_snk,),
                            "t_1" : ("wall", t_src,),
                            }
                    yield pd, tt_idx, t
    @q.timer
    def feval(args):
        pd, tt_idx, t = args
        val, val_sloppy = eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop, is_ama_and_sloppy = True)
        return val, val_sloppy, tt_idx, t
    def sum_function(val_list):
        counts = np.zeros((len(tsep_list), total_site[3],), dtype = complex)
        values = np.zeros((len(expr_names), len(tsep_list), total_site[3], 2,), dtype = complex)
        for val, val_sloppy, tt_idx, t in val_list:
            counts[tt_idx, t] += 1
            values[:, tt_idx, t, 0] += val
            values[:, tt_idx, t, 1] += val_sloppy
        return counts, values
    q.timer_fork(0)
    res_count, res_sum = q.glb_sum(
            q.parallel_map_sum(feval, load_data(), sum_function = sum_function, chunksize = 16))
    q.displayln_info("timer_display for auto_contract_meson_bk_bpi_corr")
    q.timer_display()
    q.timer_merge()
    res_count *= 1.0 / (total_volume * fsel.prob())
    res_sum *= 1.0 / (total_volume * fsel.prob())
    assert q.qnorm(res_count - 1.0) < 1e-10
    ld = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t_sep", len(tsep_list), tsep_list, ],
        [ "t_op", total_site[3], ],
        [ "ama", 2, [ "ama", "sloppy", ], ],
        ])
    ld.from_numpy(res_sum)
    # q.displayln_info(ld.show())
    ld.save(get_save_path(fn))

### ------

@q.timer_verbose
def run_job(job_tag, traj):
    fns_produce = [
            f"{job_tag}/auto-contract/traj-{traj}/checkpoint.txt",
            ]
    fns_need = [
            # (f"{job_tag}/configs/ckpoint_lat.{traj}", f"{job_tag}/configs/ckpoint_lat.IEEE64BIG.{traj}",),
            f"{job_tag}/point-selection/traj-{traj}.txt",
            f"{job_tag}/field-selection/traj-{traj}.field",
            f"{job_tag}/gauge-transform/traj-{traj}.field",
            f"{job_tag}/wall-src-info-light/traj-{traj}.txt",
            f"{job_tag}/wall-src-info-strange/traj-{traj}.txt",
            f"{job_tag}/psel-prop-wsrc-light/traj-{traj}/checkpoint.txt",
            f"{job_tag}/psel-prop-wsrc-strange/traj-{traj}/checkpoint.txt",
            # f"{job_tag}/psel-prop-psrc-light/traj-{traj}/checkpoint.txt",
            # f"{job_tag}/psel-prop-psrc-strange/traj-{traj}/checkpoint.txt",
            # f"{job_tag}/prop-wsrc-light/traj-{traj}/geon-info.txt",
            # f"{job_tag}/prop-wsrc-strange/traj-{traj}/geon-info.txt",
            # f"{job_tag}/prop-psrc-light/traj-{traj}/geon-info.txt",
            # f"{job_tag}/prop-psrc-strange/traj-{traj}/geon-info.txt",
            # f"{job_tag}/prop-rand-u1-light/traj-{traj}/geon-info.txt",
            # f"{job_tag}/prop-rand-u1-strange/traj-{traj}/geon-info.txt",
            # f"{job_tag}/prop-rand-u1-charm/traj-{traj}/geon-info.txt",
            ]
    if not check_job(job_tag, traj, fns_produce, fns_need):
        return
    #
    traj_gf = traj
    if job_tag[:5] == "test-":
        # ADJUST ME
        traj_gf = 1000
        #
    #
    get_gf = run_gf(job_tag, traj_gf)
    get_gt = run_gt(job_tag, traj_gf, get_gf)
    #
    get_psel = run_psel(job_tag, traj)
    get_fsel = run_fsel(job_tag, traj, get_psel)
    #
    get_wi = run_wi(job_tag, traj)
    get_psel_smear = run_psel_smear(job_tag, traj)
    #
    get_get_prop = run_get_prop(job_tag, traj,
            get_gt = get_gt,
            get_psel = get_psel,
            get_fsel = get_fsel,
            get_psel_smear = get_psel_smear,
            get_wi = get_wi,
            prop_types = [
                "wsrc psel s",
                "wsrc psel l",
                "wsrc fsel s",
                "wsrc fsel l",
                ],
            )
    #
    fn_checkpoint = f"{job_tag}/auto-contract/traj-{traj}/checkpoint.txt"
    if get_load_path(fn_checkpoint) is None:
        if q.obtain_lock(f"locks/{job_tag}-{traj}-auto-contract"):
            get_prop = get_get_prop()
            # ADJUST ME
            auto_contract_meson_corr(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contract_meson_f_corr(job_tag, traj, get_prop, get_psel, get_fsel)
            auto_contract_meson_bk_bpi_corr(job_tag, traj, get_prop, get_psel, get_fsel)
            #
            q.qtouch_info(get_save_path(fn_checkpoint))
            q.release_lock()
    #
    q.clean_cache()
    q.timer_display()

def get_all_cexpr():
    benchmark_eval_cexpr(get_cexpr_meson_corr())
    benchmark_eval_cexpr(get_cexpr_meson_f_corr())
    benchmark_eval_cexpr(get_cexpr_meson_bk_bpi_corr())

def test():
    # ADJUST ME
    assert q.get_num_node() <= 4
    q.qremove_all_info("results/test-4nt8")
    q.qremove_info("results")
    assert not q.does_file_exist_sync_node("results")
    q.qremove_all_info("locks")
    q.qremove_all_info("cache")
    get_all_cexpr()
    run_job("test-4nt8", 1000)
    # run_job("test-4nt16", 1000)
    # run_job("16IH2", 1000)

size_node_list = [
        [ 1, 1, 1, 1, ],
        [ 1, 1, 1, 2, ],
        [ 1, 1, 1, 4, ],
        [ 1, 1, 1, 8, ],
        ]

q.begin_with_mpi(size_node_list)

# ADJUST ME
test()

# ADJUST ME
job_tags = [
        # "test-4nt8", "test-4nt16",
        # "32IH1",
        # "32IH2",
        # "24IH1",
        # "24IH2",
        # "24IH3",
        # "24D",
        # "24DH",
        # "16IH2",
        # "32IfineH",
        # "32IcoarseH1",
        # "48I",
        ]

q.check_time_limit()

for job_tag in job_tags:
    q.displayln_info(pprint.pformat(rup.dict_params[job_tag]))
    for traj in rup.dict_params[job_tag]["trajs"]:
        run_job(job_tag, traj)

q.end_with_mpi()
