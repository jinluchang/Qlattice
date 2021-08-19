#    Qlattice (https://github.com/waterret/qlattice)
#
#    Copyright (C) 2021
#
#    Author: Luchang Jin (ljin.luchang@gmail.com)
#    Author: Masaaki Tomii
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to the Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

# operators definition follows Eq.(103,122) in PHYS. REV. D 101, 014506 (2020)

from auto_contractor.wick import *
from auto_contractor.compile import *

import math

spin_index_counter = 0

color_index_counter = 0

saved_sc_indices = []
# saved_sc_indices = None

def save_sc_indices():
    if saved_sc_indices is None:
        return
    saved_sc_indices.append([spin_index_counter, color_index_counter])

def restore_sc_indices():
    if saved_sc_indices is None:
        return
    global spin_index_counter
    global color_index_counter
    spin_index_counter, color_index_counter = saved_sc_indices.pop()

def jump_sc_indices(step = 100):
    if saved_sc_indices is None:
        return
    global spin_index_counter
    global color_index_counter
    spin_index_counter += 100
    color_index_counter += 100

def new_spin_index():
    global spin_index_counter
    spin_index_counter += 1
    return f"a_s_{spin_index_counter}"

def new_color_index():
    global color_index_counter
    color_index_counter += 1
    return f"a_c_{color_index_counter}"

def mk_meson(f1 : str, f2 : str, p : str):
    s1 = new_spin_index()
    s2 = new_spin_index()
    c = new_color_index()
    return Qb(f1, p, s1, c) * G(5, s1, s2) * Qv(f2, p, s2, c)

def mk_pi_0(p : str, is_dagger = False):
    return 1.0j / math.sqrt(2.0) * (mk_meson("u", "u", p) - mk_meson("d", "d", p))

def mk_pi_p(p : str, is_dagger = False):
    if not is_dagger:
        return 1.0j * mk_meson("u", "d", p)
    else:
        return -mk_pi_m(p)

def mk_pi_m(p : str, is_dagger = False):
    if not is_dagger:
        return -1.0j * mk_meson("d", "u", p)
    else:
        return -mk_pi_p(p)

def mk_pipi_i22(p1 : str, p2 : str, is_dagger = False):
    return mk_pi_p(p1, is_dagger) * mk_pi_p(p2, is_dagger)

def mk_pipi_i21(p1 : str, p2 : str, is_dagger = False):
    return 1.0 / math.sqrt(2.0) * (
            mk_pi_p(p1, is_dagger) * mk_pi_0(p2, is_dagger)
            + mk_pi_0(p1, is_dagger) * mk_pi_p(p2, is_dagger)
            )

def mk_pipi_i11(p1 : str, p2 : str, is_dagger = False):
    return 1.0 / math.sqrt(2.0) * (
            mk_pi_p(p1, is_dagger) * mk_pi_0(p2, is_dagger)
            - mk_pi_0(p1, is_dagger) * mk_pi_p(p2, is_dagger)
            )

def mk_pipi_i20(p1 : str, p2 : str, is_dagger = False):
    return 1.0 / math.sqrt(6.0) * (
            2.0 * mk_pi_0(p1, is_dagger) * mk_pi_0(p2, is_dagger)
            + mk_pi_m(p1, is_dagger) * mk_pi_p(p2, is_dagger)
            + mk_pi_p(p1, is_dagger) * mk_pi_m(p2, is_dagger)
            )

def mk_pipi_i10(p1 : str, p2 : str, is_dagger = False):
    return 1.0 / math.sqrt(2.0) * (
            mk_pi_p(p1, is_dagger) * mk_pi_m(p2, is_dagger)
            - mk_pi_m(p1, is_dagger) * mk_pi_p(p2, is_dagger)
            )

def mk_pipi_i0(p1 : str, p2 : str, is_dagger = False):
    return 1.0 / math.sqrt(3.0) * (
            - mk_pi_0(p1, is_dagger) * mk_pi_0(p2, is_dagger)
            + mk_pi_m(p1, is_dagger) * mk_pi_p(p2, is_dagger)
            + mk_pi_p(p1, is_dagger) * mk_pi_m(p2, is_dagger)
            )

def mk_k_p(p : str, is_dagger = False):
    if not is_dagger:
        return 1.0j * mk_meson("u", "s", p)
    else:
        return -mk_k_m(p)

def mk_k_m(p : str, is_dagger = False):
    if not is_dagger:
        return -1.0j * mk_meson("s", "u", p)
    else:
        return -mk_k_p(p)

def mk_k_0(p : str, is_dagger = False):
    if not is_dagger:
        return 1.0j * mk_meson("d", "s", p)
    else:
        return -mk_k_0_bar(p)

def mk_k_0_bar(p : str, is_dagger = False):
    if not is_dagger:
        return -1.0j * mk_meson("s", "d", p)
    else:
        return -mk_k_0(p)

def mk_kk_i11(p1 : str, p2 : str, is_dagger = False):
    return mk_k_p(p1, is_dagger) * mk_k_0_bar(p2, is_dagger)

def mk_kk_i10(p1 : str, p2 : str, is_dagger = False):
    return 1.0 / math.sqrt(2) * (
            - mk_k_0(p1, is_dagger) * mk_k_0_bar(p2, is_dagger)
            + mk_k_p(p1, is_dagger) * mk_k_m(p2, is_dagger)
            )

def mk_kk_i0(p1 : str, p2 : str, is_dagger = False):
    return 1.0 / math.sqrt(2) * (
            mk_k_0(p1, is_dagger) * mk_k_0_bar(p2, is_dagger)
            + mk_k_p(p1, is_dagger) * mk_k_m(p2, is_dagger)
            )

def mk_k0k0bar(p1 : str, p2 : str, is_dagger = False):
    return mk_k_0(p1, is_dagger) * mk_k_0_bar(p2, is_dagger)

def mk_sigma(p : str, is_dagger = False):
    s = new_spin_index()
    c = new_color_index()
    return 1.0 / math.sqrt(2.0) * (Qb("u", p, s, c) * Qv("u", p, s, c) + Qb("d", p, s, c) * Qv("d", p, s, c))

def mk_scalar(f1 : str, f2 : str, p : str):
    s = new_spin_index()
    c = new_color_index()
    return Qb(f1, p, s, c) * Qv(f2, p, s, c)

def mk_scalar5(f1 : str, f2 : str, p : str):
    s1 = new_spin_index()
    s2 = new_spin_index()
    c = new_color_index()
    return Qb(f1, p, s1, c) * G(5, s1, s2) * Qv(f2, p, s2, c)

def mk_vec_mu(f1 : str, f2 : str, p : str, mu, is_dagger = False):
    s1 = new_spin_index()
    s2 = new_spin_index()
    c = new_color_index()
    if is_dagger:
        if mu in [0, 1, 2]:
            return -Qb(f2, p, s1, c) * G(mu, s1, s2) * Qv(f1, p, s2, c)
        else:
            return Qb(f2, p, s1, c) * G(mu, s1, s2) * Qv(f1, p, s2, c)
    else:
        return Qb(f1, p, s1, c) * G(mu, s1, s2) * Qv(f2, p, s2, c)

def mk_vec5_mu(f1 : str, f2 : str, p : str, mu, is_dagger = False):
    s1 = new_spin_index()
    s2 = new_spin_index()
    s3 = new_spin_index()
    c = new_color_index()
    if is_dagger:
        if mu in [0, 1, 2]:
            return -Qb(f2, p, s1, c) * G(mu, s1, s2) * G(5, s2, s3) * Qv(f1, p, s3, c)
        else:
            return Qb(f2, p, s1, c) * G(mu, s1, s2) * G(5, s2, s3) * Qv(f1, p, s3, c)
    else:
        return Qb(f1, p, s1, c) * G(mu, s1, s2) * G(5, s2, s3) * Qv(f2, p, s3, c)

def mk_j_mu(p : str, mu, is_dagger = False):
    return 2/3 * mk_vec_mu("u", "u", p, mu, is_dagger) - 1/3 * mk_vec_mu("d", "d", p, mu, is_dagger) - 1/3 * mk_vec_mu("s", "s", p, mu, is_dagger)

def mk_jl_mu(p : str, mu, is_dagger = False):
    # jl = sqrt(2)/6 * (j0 + 3 * j10) if no s quark
    return 2/3 * mk_vec_mu("u", "u", p, mu, is_dagger) - 1/3 * mk_vec_mu("d", "d", p, mu, is_dagger)

def mk_j0_mu(p : str, mu, is_dagger = False):
    return 1.0 / math.sqrt(2.0) * (mk_vec_mu("u", "u", p, mu, is_dagger) + mk_vec_mu("d", "d", p, mu, is_dagger))

def mk_j10_mu(p : str, mu, is_dagger = False):
    return 1.0 / math.sqrt(2.0) * (mk_vec_mu("u", "u", p, mu, is_dagger) - mk_vec_mu("d", "d", p, mu, is_dagger))

def mk_j11_mu(p : str, mu, is_dagger = False):
    return mk_vec_mu("u", "d", p, mu, is_dagger)

def mk_j1n1_mu(p : str, mu, is_dagger = False):
    return -mk_vec_mu("d", "u", p, mu, is_dagger)

def mk_4qOp_VV(f1 : str, f2 : str, f3 : str, f4 : str, p, is_scalar = False, parity = None):
    if parity == "odd":
        return 0
    if is_scalar:
        return mk_4qOp_SS(f1,f2,f3,f4,p)
    sum = 0
    for mu in range(4):
        save_sc_indices()
    for mu in range(4):
        restore_sc_indices()
        sum = sum + mk_vec_mu(f1,f2,p,mu) * mk_vec_mu(f3,f4,p,mu)
    sum.simplify()
    jump_sc_indices()
    return sum

def mk_4qOp_VA(f1 : str, f2 : str, f3 : str, f4 : str, p, is_scalar = False, parity = None):
    if parity == "even":
        return 0
    if is_scalar:
        return mk_4qOp_SP(f1,f2,f3,f4,p)
    sum = 0
    for mu in range(4):
        save_sc_indices()
    for mu in range(4):
        restore_sc_indices()
        sum = sum + mk_vec_mu(f1,f2,p,mu) * mk_vec5_mu(f3,f4,p,mu)
    sum.simplify()
    jump_sc_indices()
    return sum

def mk_4qOp_AV(f1 : str, f2 : str, f3 : str, f4 : str, p, is_scalar = False, parity = None ):
    if parity == "even":
        return 0
    if is_scalar:
        return mk_4qOp_PS(f1,f2,f3,f4,p)
    sum = 0
    for mu in range(4):
        save_sc_indices()
    for mu in range(4):
        restore_sc_indices()
        sum = sum + mk_vec5_mu(f1,f2,p,mu) * mk_vec_mu(f3,f4,p,mu)
    sum.simplify()
    jump_sc_indices()
    return sum

def mk_4qOp_AA(f1 : str, f2 : str, f3 : str, f4 : str, p, is_scalar = False, parity = None ):
    if parity == "odd":
        return 0
    if is_scalar:
        return mk_4qOp_PP(f1,f2,f3,f4,p)
    sum = 0
    for mu in range(4):
        save_sc_indices()
    for mu in range(4):
        restore_sc_indices()
        sum = sum + mk_vec5_mu(f1,f2,p,mu) * mk_vec5_mu(f3,f4,p,mu)
    sum.simplify()
    jump_sc_indices()
    return sum

def mk_4qOp_SS(f1 : str, f2 : str, f3 : str, f4 : str, p ):
    return mk_scalar(f1,f2,p) * mk_scalar(f3,f4,p)

def mk_4qOp_SP(f1 : str, f2 : str, f3 : str, f4 : str, p ):
    return mk_scalar(f1,f2,p) * mk_scalar5(f3,f4,p)

def mk_4qOp_PS(f1 : str, f2 : str, f3 : str, f4 : str, p ):
    return mk_scalar5(f1,f2,p) * mk_scalar(f3,f4,p)

def mk_4qOp_PP(f1 : str, f2 : str, f3 : str, f4 : str, p ):
    return mk_scalar5(f1,f2,p) * mk_scalar5(f3,f4,p)

def rsc_call(x, *args):
    restore_sc_indices()
    return x(*args)

def mk_4qOp_LL(*args):
    for mu in range(4):
        save_sc_indices()
    expr = simplified( rsc_call(mk_4qOp_VV,*args)
                       - rsc_call(mk_4qOp_VA,*args)
                       - rsc_call(mk_4qOp_AV,*args)
                       + rsc_call(mk_4qOp_AA,*args) )
    jump_sc_indices()
    return expr

def mk_4qOp_LR(*args):
    for mu in range(4):
        save_sc_indices()
    expr = simplified( rsc_call(mk_4qOp_VV,*args)
                       + rsc_call(mk_4qOp_VA,*args)
                       - rsc_call(mk_4qOp_AV,*args)
                       - rsc_call(mk_4qOp_AA,*args) )
    jump_sc_indices()
    return expr

def mk_4qOp_RL(*args):
    for mu in range(4):
        save_sc_indices()
    expr = simplified( rsc_call(mk_4qOp_VV,*args)
                       - rsc_call(mk_4qOp_VA,*args)
                       + rsc_call(mk_4qOp_AV,*args)
                       - rsc_call(mk_4qOp_AA,*args) )
    jump_sc_indices()
    return expr

def mk_4qOp_RR(*args):
    for mu in range(4):
        save_sc_indices()
    expr = simplified( rsc_call(mk_4qOp_VV,*args)
                       + rsc_call(mk_4qOp_VA,*args)
                       + rsc_call(mk_4qOp_AV,*args)
                       + rsc_call(mk_4qOp_AA,*args) )
    jump_sc_indices()
    return expr

def mk_4qOp_LL_cmix(f1,f2,f3,f4,p,is_scalar = False, parity = None):
    assert not is_scalar
    return mk_4qOp_LL(f1,f4,f3,f2,p,is_scalar,parity)

def mk_4qOp_LR_cmix(f1,f2,f3,f4,p,is_scalar = False, parity = None):
    assert not is_scalar
    return -2 * mk_4qOp_RL(f1,f4,f3,f2,p,True,parity)

def mk_Qsub(p, parity = None):
    if parity is None:
        return mk_Qsub(p, "even") + mk_Qsub(p, "odd")
    elif parity == "even":
        return mk_scalar("s", "d", p)
    elif parity == "odd":
        return -mk_scalar5("s", "d", p)

def mk_Q1(p, parity = None):
    return mk_4qOp_LL("s","d","u","u",p,False,parity)

def mk_Q2(p, parity = None):
    return mk_4qOp_LL_cmix("s","d","u","u",p,False,parity)

def mk_Q3(p, parity = None):
    for mu in range(3):
        save_sc_indices()
    expr = simplified( rsc_call(mk_4qOp_LL,"s","d","u","u",p,False,parity)
                       + rsc_call(mk_4qOp_LL,"s","d","d","d",p,False,parity)
                       + rsc_call(mk_4qOp_LL,"s","d","s","s",p,False,parity) )
    jump_sc_indices()
    return expr

def mk_Q4(p, parity = None):
    for mu in range(3):
        save_sc_indices()
    expr = simplified( rsc_call(mk_4qOp_LL_cmix,"s","d","u","u",p,False,parity)
                       + rsc_call(mk_4qOp_LL_cmix,"s","d","d","d",p,False,parity)
                       + rsc_call(mk_4qOp_LL_cmix,"s","d","s","s",p,False,parity) )
    jump_sc_indices()
    return expr

def mk_Q5(p, parity = None):
    for mu in range(3):
        save_sc_indices()
    expr = simplified( rsc_call(mk_4qOp_LR,"s","d","u","u",p,False,parity)
                       + rsc_call(mk_4qOp_LR,"s","d","d","d",p,False,parity)
                       + rsc_call(mk_4qOp_LR,"s","d","s","s",p,False,parity) )
    jump_sc_indices()
    return expr

def mk_Q6(p, parity = None):
    for mu in range(3):
        save_sc_indices()
    expr = simplified( rsc_call(mk_4qOp_LR_cmix,"s","d","u","u",p,False,parity)
                       + rsc_call(mk_4qOp_LR_cmix,"s","d","d","d",p,False,parity)
                       + rsc_call(mk_4qOp_LR_cmix,"s","d","s","s",p,False,parity) )
    jump_sc_indices()
    return expr

def mk_Q7(p, parity = None):
    for mu in range(3):
        save_sc_indices()
    expr = simplified( rsc_call(mk_4qOp_LR,"s","d","u","u",p,False,parity)
                       - 0.5 * rsc_call(mk_4qOp_LR,"s","d","d","d",p,False,parity)
                       - 0.5 * rsc_call(mk_4qOp_LR,"s","d","s","s",p,False,parity) )
    jump_sc_indices()
    return expr

def mk_Q8(p, parity = None):
    for mu in range(3):
        save_sc_indices()
    expr = simplified( rsc_call(mk_4qOp_LR_cmix,"s","d","u","u",p,False,parity)
                       - 0.5 * rsc_call(mk_4qOp_LR_cmix,"s","d","d","d",p,False,parity)
                       - 0.5 * rsc_call(mk_4qOp_LR_cmix,"s","d","s","s",p,False,parity) )
    jump_sc_indices()
    return expr

def mk_Q9(p, parity = None):
    for mu in range(3):
        save_sc_indices()
    expr = simplified( rsc_call(mk_4qOp_LL,"s","d","u","u",p,False,parity)
                       - 0.5 * rsc_call(mk_4qOp_LL,"s","d","d","d",p,False,parity)
                       - 0.5 * rsc_call(mk_4qOp_LL,"s","d","s","s",p,False,parity) )
    jump_sc_indices()
    return expr

def mk_Q10(p, parity = None):
    for mu in range(3):
        save_sc_indices()
    expr = simplified( rsc_call(mk_4qOp_LL_cmix,"s","d","u","u",p,False,parity)
                       - 0.5 * rsc_call(mk_4qOp_LL_cmix,"s","d","d","d",p,False,parity)
                       - 0.5 * rsc_call(mk_4qOp_LL_cmix,"s","d","s","s",p,False,parity) )
    jump_sc_indices()
    return expr

def test():
    print("test")
    args = ["x", None]
    for mu in range(11):
        save_sc_indices()
    expr1 = simplified( rsc_call(mk_Q1,*args)
                       - rsc_call(mk_Q2,*args)
                       - rsc_call(mk_Q3,*args)
                       + rsc_call(mk_Q4,*args) )
    expr2 = simplified( 3 * rsc_call(mk_Q1,*args)
                       - rsc_call(mk_Q3,*args)
                       - 2 * rsc_call(mk_Q9,*args) )
    expr3 = simplified( rsc_call(mk_Q1,*args)
                       + 2 * rsc_call(mk_Q2,*args)
                       - rsc_call(mk_Q3,*args)
                       - 2 * rsc_call(mk_Q10,*args))
    jump_sc_indices()
    expr1 = mk_pipi_i0("x1_1", "x1_2", True) * expr1 * mk_k_0("x2")
    expr2 = mk_pipi_i0("x1_1", "x1_2", True) * expr2 * mk_k_0("x2")
    expr3 = mk_pipi_i0("x1_1", "x1_2", True) * expr3 * mk_k_0("x2")
    print(display_cexpr(contract_simplify_round_compile(expr1, expr2, expr3)))
    print(display_cexpr(contract_simplify_round_compile(
        mk_pipi_i0("x1_1", "x1_2", True) * mk_Q1("x") * mk_k_0("x2"),
        mk_pipi_i0("x1_1", "x1_2", True) * mk_Q2("x") * mk_k_0("x2"),
        mk_pipi_i0("x1_1", "x1_2", True) * mk_Q3("x") * mk_k_0("x2"),
        mk_pipi_i0("x1_1", "x1_2", True) * mk_Q4("x") * mk_k_0("x2"),
        mk_pipi_i0("x1_1", "x1_2", True) * mk_Q5("x") * mk_k_0("x2"),
        mk_pipi_i0("x1_1", "x1_2", True) * mk_Q6("x") * mk_k_0("x2"),
        mk_pipi_i0("x1_1", "x1_2", True) * mk_Q7("x") * mk_k_0("x2"),
        mk_pipi_i0("x1_1", "x1_2", True) * mk_Q8("x") * mk_k_0("x2"),
        mk_pipi_i0("x1_1", "x1_2", True) * mk_Q9("x") * mk_k_0("x2"),
        mk_pipi_i0("x1_1", "x1_2", True) * mk_Q10("x") * mk_k_0("x2"),
        )))

def test1():
    def A(j_p, pi_p, is_dagger = False):
        return (mk_j10_mu(j_p, 0, is_dagger) * mk_pi_p(pi_p, is_dagger)
                + mk_j11_mu(j_p, 0, is_dagger) * mk_pi_0(pi_p, is_dagger))
    def B(j_p, pi_p, is_dagger = False):
        return (2 * mk_j10_mu(j_p, 0, is_dagger) * mk_pi_0(pi_p, is_dagger)
                + mk_j1n1_mu(j_p, 0, is_dagger) * mk_pi_p(pi_p, is_dagger)
                + mk_j11_mu(j_p, 0, is_dagger) * mk_pi_m(pi_p, is_dagger))
    def C(j_p, pi_p, is_dagger = False):
        return (mk_j10_mu(j_p, 0, is_dagger) * mk_pi_p(pi_p, is_dagger)
                - mk_j11_mu(j_p, 0, is_dagger) * mk_pi_0(pi_p, is_dagger))
    def D(j_p, pi_p, is_dagger = False):
        return (mk_j1n1_mu(j_p, 0, is_dagger) * mk_pi_p(pi_p, is_dagger)
                - mk_j11_mu(j_p, 0, is_dagger) * mk_pi_m(pi_p, is_dagger))
    def E(j_p, pi_p, is_dagger = False):
        return (-mk_j10_mu(j_p, 0, is_dagger) * mk_pi_0(pi_p, is_dagger)
                + mk_j1n1_mu(j_p, 0, is_dagger) * mk_pi_p(pi_p, is_dagger)
                + mk_j11_mu(j_p, 0, is_dagger) * mk_pi_m(pi_p, is_dagger))
    def F(j_p, pi_p, is_dagger = False):
        return mk_j0_mu(j_p, 0, is_dagger) * mk_pi_p(pi_p, is_dagger)
    def G(j_p, pi_p, is_dagger = False):
        return mk_j0_mu(j_p, 0, is_dagger) * mk_pi_0(pi_p, is_dagger)
    def jpi_p(j_p, pi_p, is_dagger = False):
        return mk_jl_mu(j_p, 0, is_dagger) * mk_pi_p(pi_p, is_dagger)
    def jpi_0(j_p, pi_p, is_dagger = False):
        return mk_jl_mu(j_p, 0, is_dagger) * mk_pi_0(pi_p, is_dagger)
    expr1 = 1/4 * (A("x2_j", "x2_p", True) * A("x1_j", "x1_p") + A("x2_j", "x2_p") * A("x1_j", "x1_p", True))
    expr2 = 1/12 * (B("x2_j", "x2_p", True) * B("x1_j", "x1_p") + B("x2_j", "x2_p") * B("x1_j", "x1_p", True))
    expr3 = 1/4 * (C("x2_j", "x2_p", True) * C("x1_j", "x1_p") + C("x2_j", "x2_p") * C("x1_j", "x1_p", True))
    expr4 = 1/4 * (D("x2_j", "x2_p", True) * D("x1_j", "x1_p") + D("x2_j", "x2_p") * D("x1_j", "x1_p", True))
    expr5 = 1/6 * (E("x2_j", "x2_p", True) * E("x1_j", "x1_p") + E("x2_j", "x2_p") * E("x1_j", "x1_p", True))
    expr6 = 1/2 * (F("x2_j", "x2_p", True) * F("x1_j", "x1_p") + F("x2_j", "x2_p") * F("x1_j", "x1_p", True))
    expr7 = 1/2 * (G("x2_j", "x2_p", True) * G("x1_j", "x1_p") + G("x2_j", "x2_p") * G("x1_j", "x1_p", True))
    expr8 = 1/(4*math.sqrt(2)) * (C("x2_j", "x2_p", True) * F("x1_j", "x1_p") + F("x2_j", "x2_p", True) * C("x1_j", "x1_p")
            + C("x2_j", "x2_p") * F("x1_j", "x1_p", True) + F("x2_j", "x2_p") * C("x1_j", "x1_p", True))
    expr9 = 1/(4*math.sqrt(2)) * (D("x2_j", "x2_p", True) * G("x1_j", "x1_p") + G("x2_j", "x2_p", True) * D("x1_j", "x1_p")
            + D("x2_j", "x2_p") * G("x1_j", "x1_p", True) + G("x2_j", "x2_p") * D("x1_j", "x1_p", True))
    expr_p = 1/2*(jpi_p("x2_j", "x2_p", True) * jpi_p("x1_j", "x1_p") + jpi_p("x2_j", "x2_p") * jpi_p("x1_j", "x1_p", True))
    expr_0 = 1/2*(jpi_0("x2_j", "x2_p", True) * jpi_0("x1_j", "x1_p") + jpi_0("x2_j", "x2_p") * jpi_0("x1_j", "x1_p", True))
    exprs = [expr1, expr1 - expr2, expr3, expr3 - expr4, expr5, expr6, expr6 - expr7, expr8, expr8 - expr9,]
    print(display_cexpr(contract_simplify_round_compile(*exprs)))
    exprs1 = [
            expr_p - 1/18*(9/2*expr1 + 9/2*expr3 + expr6 + 3*math.sqrt(2)*expr8),
            expr_0 - 1/18*(6*expr2 + 3*expr5 + expr7),
            (expr_p - expr_0),
            ]
    print(display_cexpr(contract_simplify_round_compile(*exprs1)))

def test_kk():
    expr1 = mk_kk_i0("x2_1", "x2_2", True) * mk_kk_i0("x1_1", "x1_2")
    expr2 = mk_kk_i0("x2_1", "x2_2", True) * mk_pipi_i0("x1_1", "x1_2")
    expr3 = mk_kk_i0("x2_1", "x2_2", True) * mk_sigma("x1", "x1")
    expr4 = mk_k0k0bar("x2_1", "x2_2", True) * mk_k0k0bar("x1_1", "x1_2")
    expr5 = mk_k0k0bar("x2_1", "x2_2", True) * mk_pipi_i0("x1_1", "x1_2")
    expr6 = mk_k0k0bar("x2_1", "x2_2", True) * mk_sigma("x1", "x1")
    all_exprs = [
            [ expr1, expr4, ],
            [ expr2, expr5, ],
            [ expr3, expr6, ],
            ]
    names = [ "kk kk", "kk pipi", "kk sigma", ]
    for name, exprs in zip(names, all_exprs):
        print(f"\n-- {name} --")
        cexpr = contract_simplify_round_compile(*exprs)
        print(display_cexpr(cexpr))
        cexpr.collect_op()
        print(display_cexpr(cexpr))
        print(f"-- {name} --\n")

if __name__ == "__main__":
    # test1()
    # exit()
    #
    print("pi+(x1):\n", mk_pi_p("x1"))
    print("pi+(x2)^dag pi+(x1):\n", (mk_pi_p("x2", True) * mk_pi_p("x1")).round())
    print("< pi+(x2)^dag pi+(x1) >:")
    expr = mk_pi_p("x2", True) * mk_pi_p("x1")
    print(display_cexpr(contract_simplify_round_compile(expr)))
    print()
    print("pi0(x1):\n", mk_pi_0("x1"))
    print("pi0(x2)^dag pi0(x1):\n", (mk_pi_0("x2", True) * mk_pi_0("x1")).round())
    print("< pi0(x2)^dag pi0(x1) >:")
    expr = mk_pi_0("x2", True) * mk_pi_0("x1")
    print(display_cexpr(contract_simplify_round_compile(expr, is_isospin_symmetric_limit = False)))
    print(display_cexpr(contract_simplify_round_compile(expr)))
    print()
    print("< pi0(x2)^dag pi0(x1) > - < pi+(x2)^dag pi+(x1) >: ",
            simplified(contract_expr(mk_pi_0("x2", True) * mk_pi_0("x1") - mk_pi_p("x2", True) * mk_pi_p("x1")), is_isospin_symmetric_limit = True).round())
    print()
    print("< pipiI22(x2_1,x2_2)^dag pipiI22(x1_1,x1_2) >:")
    expr = mk_pipi_i22("x2_1", "x2_2", True) * mk_pipi_i22("x1_1", "x1_2")
    print(display_cexpr(contract_simplify_round_compile(expr)))
    print("< pipiI21(x2_1,x2_2)^dag pipiI21(x1_1,x1_2) - pipiI22(x2_1,x2_2)^dag pipiI22(x1_1,x1_2) >:")
    expr = mk_pipi_i21("x2_1", "x2_2", True) * mk_pipi_i21("x1_1", "x1_2") - mk_pipi_i22("x2_1", "x2_2", True) * mk_pipi_i22("x1_1", "x1_2")
    print(display_cexpr(contract_simplify_round_compile(expr)))
    print("< pipiI20(x2_1,x2_2)^dag pipiI20(x1_1,x1_2) - pipiI22(x2_1,x2_2)^dag pipiI22(x1_1,x1_2) >:")
    expr = mk_pipi_i20("x2_1", "x2_2", True) * mk_pipi_i20("x1_1", "x1_2") - mk_pipi_i22("x2_1", "x2_2", True) * mk_pipi_i22("x1_1", "x1_2")
    print(display_cexpr(contract_simplify_round_compile(expr)))
    print()
    print("< pipiI11(x2_1,x2_2)^dag pipiI11(x1_1,x1_2) >:")
    expr = mk_pipi_i11("x2_1", "x2_2", True) * mk_pipi_i11("x1_1", "x1_2")
    print(display_cexpr(contract_simplify_round_compile(expr)))
    print("< pipiI10(x2_1,x2_2)^dag pipiI10(x1_1,x1_2) - pipiI11(x2_1,x2_2)^dag pipiI11(x1_1,x1_2) >:")
    expr = mk_pipi_i10("x2_1", "x2_2", True) * mk_pipi_i10("x1_1", "x1_2") - mk_pipi_i11("x2_1", "x2_2", True) * mk_pipi_i11("x1_1", "x1_2")
    print(display_cexpr(contract_simplify_round_compile(expr)))
    print()
    print("< pipiI0(x1_1,x1_2) >:")
    expr = mk_pipi_i0("x1_1", "x1_2")
    print(display_cexpr(contract_simplify_round_compile(expr)))
    print("< pipiI0(x2_1,x2_2)^dag pipiI0(x1_1,x1_2) >:")
    expr = mk_pipi_i0("x2_1", "x2_2", True) * mk_pipi_i0("x1_1", "x1_2")
    print(display_cexpr(contract_simplify_round_compile(expr)))
    print()
    print("< sigma(x1) >:")
    expr = mk_sigma("x1")
    print(display_cexpr(contract_simplify_round_compile(expr)))
    print("< pipiI0(x2_1,x2_2)^dag sigma(x1) >:")
    expr = mk_pipi_i0("x2_1", "x2_2", True) * mk_sigma("x1")
    print(display_cexpr(contract_simplify_round_compile(expr)))
    print()
    print("< K+(x2)^dag K+(x1)>:")
    expr = mk_k_p("x2", True) * mk_k_p("x1")
    print(display_cexpr(contract_simplify_round_compile(expr)))
    print()
    print("< pi0(x2)^dag j_mu(xj_1) j_nu(xj_2) pi0(x1) >:")
    expr = mk_pi_0("x1", True) * mk_j_mu("xj_1", "mu") * mk_j_mu("xj_2", "nu") * mk_pi_0("x2")
    print(display_cexpr(contract_simplify_round_compile(expr)))
    print()
    print("< pi+(x2)^dag j_mu(xj_1) j_nu(xj_2) pi+(x1) / 2 + pi-(x2)^dag j_mu(xj_1) j_nu(xj_2) pi-(x1) / 2 >:")
    expr = (
            1/2 * mk_pi_p("x2", True) * mk_j_mu("xj_1", "mu") * mk_j_mu("xj_2", "nu") * mk_pi_p("x1")
            + 1/2 * mk_pi_m("x2", True) * mk_j_mu("xj_1", "mu") * mk_j_mu("xj_2", "nu") * mk_pi_m("x1")
            )
    print(display_cexpr(contract_simplify_round_compile(expr)))
    print()
    print("< pi+(x2)^dag j_mu(xj_1) j_nu(xj_2) pi+(x1) / 2 + pi-(x2)^dag j_mu(xj_1) j_nu(xj_2) pi-(x1) / 2 - pi0(x2)^dag j_mu(xj_1) j_nu(xj_2) pi0(x1) >:")
    expr = (
            1/2 * mk_pi_p("x2", True) * mk_j_mu("xj_1", "mu") * mk_j_mu("xj_2", "nu") * mk_pi_p("x1")
            + 1/2 * mk_pi_m("x2", True) * mk_j_mu("xj_1", "mu") * mk_j_mu("xj_2", "nu") * mk_pi_m("x1")
            - mk_pi_0("x1", True) * mk_j_mu("xj_1", "mu") * mk_j_mu("xj_2", "nu") * mk_pi_0("x2")
            )
    print(display_cexpr(contract_simplify_round_compile(expr)))
    #
    print()
    test_kk()
    #
    print()
    test()


