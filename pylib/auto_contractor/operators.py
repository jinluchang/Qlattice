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

import sympy

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

def mk_scalar(f1 : str, f2 : str, p : str, is_dagger = False):
    s = new_spin_index()
    c = new_color_index()
    if is_dagger:
        return Qb(f2, p, s, c) * Qv(f1, p, s, c)
    else:
        return Qb(f1, p, s, c) * Qv(f2, p, s, c)

def mk_scalar5(f1 : str, f2 : str, p : str, is_dagger = False):
    s1 = new_spin_index()
    s2 = new_spin_index()
    c = new_color_index()
    if is_dagger:
        return Qb(f2, p, s1, c) * G(5, s1, s2) * Qv(f1, p, s2, c)
    else:
        return Qb(f1, p, s1, c) * G(5, s1, s2) * Qv(f2, p, s2, c)

def mk_vec_mu(f1 : str, f2 : str, p : str, mu, is_dagger = False):
    s1 = new_spin_index()
    s2 = new_spin_index()
    c = new_color_index()
    if is_dagger:
        if mu in [ 0, 1, 2, 5 ]:
            return -Qb(f2, p, s1, c) * G(mu, s1, s2) * Qv(f1, p, s2, c)
        else:
            assert mu in [ 3, ]
            return Qb(f2, p, s1, c) * G(mu, s1, s2) * Qv(f1, p, s2, c)
    else:
        return Qb(f1, p, s1, c) * G(mu, s1, s2) * Qv(f2, p, s2, c)

def mk_vec5_mu(f1 : str, f2 : str, p : str, mu, is_dagger = False):
    s1 = new_spin_index()
    s2 = new_spin_index()
    s3 = new_spin_index()
    c = new_color_index()
    if is_dagger:
        if mu in [ 0, 1, 2, ]:
            return -Qb(f2, p, s1, c) * G(mu, s1, s2) * G(5, s2, s3) * Qv(f1, p, s3, c)
        else:
            assert mu in [ 3, 5, ]
            return Qb(f2, p, s1, c) * G(mu, s1, s2) * G(5, s2, s3) * Qv(f1, p, s3, c)
    else:
        return Qb(f1, p, s1, c) * G(mu, s1, s2) * G(5, s2, s3) * Qv(f2, p, s3, c)

def mk_meson(f1 : str, f2 : str, p : str, is_dagger = False):
    return mk_scalar5(f1, f2, p, is_dagger)

def show_dagger(is_dagger):
    if is_dagger:
        return "^dag"
    else:
        return ""

def mk_pi_0(p : str, is_dagger = False):
    return sympy.I / sympy.sqrt(2) * (mk_meson("u", "u", p) - mk_meson("d", "d", p)) + f"pi0({p})"

def mk_pi_p(p : str, is_dagger = False):
    if not is_dagger:
        return sympy.I * mk_meson("u", "d", p) + f"pi+({p})"
    else:
        return -mk_pi_m(p) + f"pi+({p})^dag"

def mk_pi_m(p : str, is_dagger = False):
    if not is_dagger:
        return -sympy.I * mk_meson("d", "u", p) + f"pi-({p})"
    else:
        return -mk_pi_p(p) + f"pi-({p})^dag"

def mk_k_p(p : str, is_dagger = False):
    if not is_dagger:
        return sympy.I * mk_meson("u", "s", p) + f"K+({p})"
    else:
        return -mk_k_m(p) + f"K+({p})^dag"

def mk_k_m(p : str, is_dagger = False):
    if not is_dagger:
        return -sympy.I * mk_meson("s", "u", p) + f"K-({p})"
    else:
        return -mk_k_p(p) + f"K-({p})^dag"

def mk_k_0(p : str, is_dagger = False):
    if not is_dagger:
        return sympy.I * mk_meson("d", "s", p) + f"K0({p})"
    else:
        return -mk_k_0_bar(p) + f"K0({p})^dag"

def mk_k_0_bar(p : str, is_dagger = False):
    if not is_dagger:
        return -sympy.I * mk_meson("s", "d", p) + f"K0b({p})"
    else:
        return -mk_k_0(p) + f"K0b({p})^dag"

def mk_j5pi_mu(p : str, mu, is_dagger = False):
    return mk_vec5_mu("d", "u", p, mu, is_dagger)

def mk_j5k_mu(p : str, mu, is_dagger = False):
    return mk_vec5_mu("s", "u", p, mu, is_dagger)

def mk_pipi_i22(p1 : str, p2 : str, is_dagger = False):
    return mk_pi_p(p1, is_dagger) * mk_pi_p(p2, is_dagger) + f"pipi_I22({p1},{p2}){show_dagger(is_dagger)}"

def mk_pipi_i21(p1 : str, p2 : str, is_dagger = False):
    return 1 / sympy.sqrt(2) * (
            mk_pi_p(p1, is_dagger) * mk_pi_0(p2, is_dagger)
            + mk_pi_0(p1, is_dagger) * mk_pi_p(p2, is_dagger)
            ) + f"pipi_I21({p1},{p2}){show_dagger(is_dagger)}"

def mk_pipi_i11(p1 : str, p2 : str, is_dagger = False):
    return 1 / sympy.sqrt(2) * (
            mk_pi_p(p1, is_dagger) * mk_pi_0(p2, is_dagger)
            - mk_pi_0(p1, is_dagger) * mk_pi_p(p2, is_dagger)
            ) + f"pipi_I11({p1},{p2}){show_dagger(is_dagger)}"

def mk_pipi_i20(p1 : str, p2 : str, is_dagger = False):
    return 1 / sympy.sqrt(6) * (
            2 * mk_pi_0(p1, is_dagger) * mk_pi_0(p2, is_dagger)
            + mk_pi_m(p1, is_dagger) * mk_pi_p(p2, is_dagger)
            + mk_pi_p(p1, is_dagger) * mk_pi_m(p2, is_dagger)
            ) + f"pipi_I20({p1},{p2}){show_dagger(is_dagger)}"

def mk_pipi_i10(p1 : str, p2 : str, is_dagger = False):
    return 1 / sympy.sqrt(2) * (
            mk_pi_p(p1, is_dagger) * mk_pi_m(p2, is_dagger)
            - mk_pi_m(p1, is_dagger) * mk_pi_p(p2, is_dagger)
            ) + f"pipi_I10({p1},{p2}){show_dagger(is_dagger)}"

def mk_pipi_i0(p1 : str, p2 : str, is_dagger = False):
    return 1 / sympy.sqrt(3) * (
            - mk_pi_0(p1, is_dagger) * mk_pi_0(p2, is_dagger)
            + mk_pi_m(p1, is_dagger) * mk_pi_p(p2, is_dagger)
            + mk_pi_p(p1, is_dagger) * mk_pi_m(p2, is_dagger)
            ) + f"pipi_I0({p1},{p2}){show_dagger(is_dagger)}"

def mk_k_0_star_mu(p : str, mu, is_dagger = False):
    return mk_vec_mu("d", "s", mu, is_dagger)

def mk_k_0_star_bar_mu(p : str, mu, is_dagger = False):
    return mk_vec_mu("s", "d", mu, is_dagger)

def mk_kk_i11(p1 : str, p2 : str, is_dagger = False, *, is_sym = False):
    if is_sym:
        return 1 / sympy.sqrt(2) * (mk_kk_i11(p1, p2, is_dagger) + mk_kk_i11(p2, p1, is_dagger)) + f"KK_I11({p1},{p2},sym){show_dagger(is_dagger)}"
    else:
        return mk_k_p(p1, is_dagger) * mk_k_0_bar(p2, is_dagger) + f"KK_I11({p1},{p2}){show_dagger(is_dagger)}"

def mk_kk_i10(p1 : str, p2 : str, is_dagger = False, *, is_sym = False):
    if is_sym:
        return 1 / sympy.sqrt(2) * (mk_kk_i10(p1, p2, is_dagger) + mk_kk_i10(p2, p1, is_dagger)) + f"KK_I10({p1},{p2},sym){show_dagger(is_dagger)}"
    else:
        return 1 / sympy.sqrt(2) * (
            - mk_k_0(p1, is_dagger) * mk_k_0_bar(p2, is_dagger)
            + mk_k_p(p1, is_dagger) * mk_k_m(p2, is_dagger)
            ) + f"KK_I10({p1},{p2}){show_dagger(is_dagger)}"

def mk_kk_i0(p1 : str, p2 : str, is_dagger = False, *, is_sym = False):
    if is_sym:
        return 1 / sympy.sqrt(2) * (mk_kk_i0(p1, p2, is_dagger) + mk_kk_i0(p2, p1, is_dagger)) + f"KK_I0({p1},{p2},sym){show_dagger(is_dagger)}"
    else:
        return 1 / sympy.sqrt(2) * (
            mk_k_0(p1, is_dagger) * mk_k_0_bar(p2, is_dagger)
            + mk_k_p(p1, is_dagger) * mk_k_m(p2, is_dagger)
            ) + f"KK_I0({p1},{p2}){show_dagger(is_dagger)}"

def mk_k0k0bar(p1 : str, p2 : str, is_dagger = False, *, is_sym = False):
    if is_sym:
        return 1 / sympy.sqrt(2) * (mk_k0k0bar(p1, p2, is_dagger) + mk_k0k0bar(p2, p1, is_dagger)) + f"K0K0b({p1},{p2},sym){show_dagger(is_dagger)}"
    else:
        return mk_k_0(p1, is_dagger) * mk_k_0_bar(p2, is_dagger) + f"K0K0b({p1},{p2}){show_dagger(is_dagger)}"

def mk_k0pi0(p1 : str, p2 : str, is_dagger = False, *, is_sym = False):
    if is_sym:
        return 1 / sympy.sqrt(2) * (mk_k0pi0(p1, p2, is_dagger) + mk_k0pi0(p2, p1, is_dagger)) + f"K0pi0({p1},{p2},sym){show_dagger(is_dagger)}"
    else:
        return mk_k_0(p1, is_dagger) * mk_pi_0(p2, is_dagger) + f"K0pi0({p1},{p2}){show_dagger(is_dagger)}"

def mk_k0barpi0(p1 : str, p2 : str, is_dagger = False, *, is_sym = False):
    if is_sym:
        return 1 / sympy.sqrt(2) * (mk_k0barpi0(p1, p2, is_dagger) + mk_k0barpi0(p2, p1, is_dagger)) + f"K0barpi0({p1},{p2},sym){show_dagger(is_dagger)}"
    else:
        return mk_k_0_bar(p1, is_dagger) * mk_pi_0(p2, is_dagger) + f"K0barpi0({p1},{p2}){show_dagger(is_dagger)}"

def mk_kppim(p1 : str, p2 : str, is_dagger = False, *, is_sym = False):
    if is_sym:
        return 1 / sympy.sqrt(2) * (mk_kppim(p1, p2, is_dagger) + mk_kppim(p2, p1, is_dagger)) + f"K+pi-({p1},{p2},sym){show_dagger(is_dagger)}"
    else:
        return mk_k_p(p1, is_dagger) * mk_pi_m(p2, is_dagger) + f"K+pi-({p1},{p2}){show_dagger(is_dagger)}"

def mk_kmpip(p1 : str, p2 : str, is_dagger = False, *, is_sym = False):
    if is_sym:
        return 1 / sympy.sqrt(2) * (mk_kmpip(p1, p2, is_dagger) + mk_kmpip(p2, p1, is_dagger)) + f"K-pi+({p1},{p2},sym){show_dagger(is_dagger)}"
    else:
        return mk_k_m(p1, is_dagger) * mk_pi_p(p2, is_dagger) + f"K-pi+({p1},{p2}){show_dagger(is_dagger)}"

def mk_kpi_0_i1half(p1 : str, p2: str, is_dagger = False, *, is_sym = False):# strangeness = -1
    if is_sym:
        return 1 / sympy.sqrt(2) * (mk_kpi_0_i1half(p1, p2, is_dagger) + mk_kpi_0_i1half(p2, p1, is_dagger)) + f"Kpi_0_I1half({p1},{p2},sym){show_dagger(is_dagger)}"
    else:
        return simplified( sympy.simplify(1)/sympy.sqrt(3)* mk_k_0(p1, is_dagger) * mk_pi_0(p2, is_dagger)
                           + sympy.sqrt(2)/sympy.sqrt(3)* mk_k_p(p1, is_dagger) * mk_pi_m(p2, is_dagger) ) + f"Kpi_0_I1half({p1},{p2}){show_dagger(is_dagger)}"

def mk_kpi_p_i1half(p1 : str, p2: str, is_dagger = False, *, is_sym = False):# strangeness = -1
    if is_sym:
        return 1 / sympy.sqrt(2) * (mk_kpi_m_i1half(p1, p2, is_dagger) + mk_kpi_m_i1half(p2, p1, is_dagger)) + f"Kpi_+_I1half({p1},{p2},sym){show_dagger(is_dagger)}"
    else:
        return simplified( sympy.sqrt(2)/sympy.sqrt(3)* mk_k_0(p1, is_dagger) * mk_pi_p(p2, is_dagger)
                           + sympy.simplify(1)/sympy.sqrt(3)* mk_k_p(p1, is_dagger) * mk_pi_0(p2, is_dagger) ) + f"Kpi_+_I1half({p1},{p2}){show_dagger(is_dagger)}"

def mk_kpi_m_i3halves(p1 : str, p2: str, is_dagger = False, *, is_sym = False):# strangeness = -1
    if is_sym:
        return 1 / sympy.sqrt(2) * ( mk_kpi_m_i3halves(p1, p2, is_dagger) + mk_kpi_m_i3halves(p2, p1, is_dagger) ) + f"Kpi_-_I3halves({p1},{p2},sym){show_dagger(is_dagger)}"
    else:
        return mk_k_0(p1, is_dagger) * mk_pi_m(p2, is_dagger) + f"Kpi_-_I3halves({p1},{p2}){show_dagger(is_dagger)}"

def mk_kpi_0_i3halves(p1 : str, p2: str, is_dagger = False, *, is_sym = False):# strangeness = -1
    if is_sym:
        return 1 / sympy.sqrt(2) * (mk_kpi_0_i3halves(p1, p2, is_dagger) + mk_kpi_0_i3halves(p2, p1, is_dagger)) + f"Kpi_0_I3halves({p1},{p2},sym){show_dagger(is_dagger)}"
    else:
        return simplified( - sympy.sqrt(2)/sympy.sqrt(3)* mk_k_0(p1, is_dagger) * mk_pi_0(p2, is_dagger)
                           + sympy.simplify(1)/sympy.sqrt(3)* mk_k_p(p1, is_dagger) * mk_pi_m(p2, is_dagger) ) + f"Kpi_0_I3halves({p1},{p2}){show_dagger(is_dagger)}"

def mk_kpi_p1_i3halves(p1 : str, p2: str, is_dagger = False, *, is_sym = False):# strangeness = -1
    if is_sym:
        return 1 / sympy.sqrt(2) * (mk_kpi_p1_i3halves(p1, p2, is_dagger) + mk_kpi_p1_i3halves(p2, p1, is_dagger)) + f"Kpi_+_I3halves({p1},{p2},sym){show_dagger(is_dagger)}"
    else:
        return simplified( - sympy.simplify(1)/sympy.sqrt(3)* mk_k_0(p1, is_dagger) * mk_pi_p(p2, is_dagger)
                           + sympy.sqrt(2)/sympy.sqrt(3)* mk_k_p(p1, is_dagger) * mk_pi_0(p2, is_dagger) ) + f"Kpi_+_I3halves({p1},{p2}){show_dagger(is_dagger)}"

def mk_kpi_p2_i3halves(p1 : str, p2: str, is_dagger = False, *, is_sym = False):# strangeness = -1
    if is_sym:
        return 1 / sympy.sqrt(2) * (mk_kpi_p2_i3halves(p1, p2, is_dagger) + mk_kpi_p2_i3halves(p2, p1, is_dagger)) + f"Kpi_++_I3halves({p1},{p2},sym){show_dagger(is_dagger)}"
    else:
        return mk_k_p(p1, is_dagger) * mk_pi_p(p2, is_dagger) + f"Kpi_++_I3halves({p1},{p2}){show_dagger(is_dagger)}"

def mk_sigma(p : str, is_dagger = False):
    s = new_spin_index()
    c = new_color_index()
    return 1 / sympy.sqrt(2) * (Qb("u", p, s, c) * Qv("u", p, s, c) + Qb("d", p, s, c) * Qv("d", p, s, c)) + f"sigma({p})"

def mk_j_mu(p : str, mu, is_dagger = False):
    return sympy.simplify(1)*2/3 * mk_vec_mu("u", "u", p, mu, is_dagger) - sympy.simplify(1)*1/3 * mk_vec_mu("d", "d", p, mu, is_dagger) - sympy.simplify(1)*1/3 * mk_vec_mu("s", "s", p, mu, is_dagger) + f"j_mu({p},{mu}){show_dagger(is_dagger)}"

def mk_jl_mu(p : str, mu, is_dagger = False):
    # jl = sqrt(2)/6 * (j0 + 3 * j10) if no s quark
    return sympy.simplify(1)*2/3 * mk_vec_mu("u", "u", p, mu, is_dagger) - sympy.simplify(1)*1/3 * mk_vec_mu("d", "d", p, mu, is_dagger)

def mk_j0_mu(p : str, mu, is_dagger = False):
    return sympy.simplify(1)*1/sympy.sqrt(2) * (mk_vec_mu("u", "u", p, mu, is_dagger) + mk_vec_mu("d", "d", p, mu, is_dagger))

def mk_j10_mu(p : str, mu, is_dagger = False):
    return sympy.simplify(1)*1/sympy.sqrt(2) * (mk_vec_mu("u", "u", p, mu, is_dagger) - mk_vec_mu("d", "d", p, mu, is_dagger))

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

def show_parity(parity):
    if parity is None:
        return ""
    elif parity == "even":
        return ",e"
    elif parity == "odd":
        return ",o"
    else:
        return parity

def mk_Qsub(p, parity = None):
    if parity is None:
        expr = mk_Qsub(p, "even") + mk_Qsub(p, "odd")
    elif parity == "even":
        expr = mk_scalar("s", "d", p)
    elif parity == "odd":
        expr = -mk_scalar5("s", "d", p)
    return expr + f"Qsub({p}{show_parity(parity)})"

def mk_Q1(p, parity = None):
    return mk_4qOp_LL("s","d","u","u",p,False,parity) + f"Q1({p}{show_parity(parity)})"

def mk_Q2(p, parity = None):
    return mk_4qOp_LL_cmix("s","d","u","u",p,False,parity) + f"Q2({p}{show_parity(parity)})"

def mk_Q3(p, parity = None):
    for mu in range(3):
        save_sc_indices()
    expr = simplified( rsc_call(mk_4qOp_LL,"s","d","u","u",p,False,parity)
                       + rsc_call(mk_4qOp_LL,"s","d","d","d",p,False,parity)
                       + rsc_call(mk_4qOp_LL,"s","d","s","s",p,False,parity) )
    jump_sc_indices()
    return expr + f"Q3({p}{show_parity(parity)})"

def mk_Q4(p, parity = None):
    for mu in range(3):
        save_sc_indices()
    expr = simplified( rsc_call(mk_4qOp_LL_cmix,"s","d","u","u",p,False,parity)
                       + rsc_call(mk_4qOp_LL_cmix,"s","d","d","d",p,False,parity)
                       + rsc_call(mk_4qOp_LL_cmix,"s","d","s","s",p,False,parity) )
    jump_sc_indices()
    return expr + f"Q4({p}{show_parity(parity)})"

def mk_Q5(p, parity = None):
    for mu in range(3):
        save_sc_indices()
    expr = simplified( rsc_call(mk_4qOp_LR,"s","d","u","u",p,False,parity)
                       + rsc_call(mk_4qOp_LR,"s","d","d","d",p,False,parity)
                       + rsc_call(mk_4qOp_LR,"s","d","s","s",p,False,parity) )
    jump_sc_indices()
    return expr + f"Q5({p}{show_parity(parity)})"

def mk_Q6(p, parity = None):
    for mu in range(3):
        save_sc_indices()
    expr = simplified( rsc_call(mk_4qOp_LR_cmix,"s","d","u","u",p,False,parity)
                       + rsc_call(mk_4qOp_LR_cmix,"s","d","d","d",p,False,parity)
                       + rsc_call(mk_4qOp_LR_cmix,"s","d","s","s",p,False,parity) )
    jump_sc_indices()
    return expr + f"Q6({p}{show_parity(parity)})"

def mk_Q7(p, parity = None):
    for mu in range(3):
        save_sc_indices()
    expr = simplified( rsc_call(mk_4qOp_LR,"s","d","u","u",p,False,parity)
                       - sympy.simplify(1)/2* rsc_call(mk_4qOp_LR,"s","d","d","d",p,False,parity)
                       - sympy.simplify(1)/2* rsc_call(mk_4qOp_LR,"s","d","s","s",p,False,parity) )
    jump_sc_indices()
    return expr + f"Q7({p}{show_parity(parity)})"

def mk_Q8(p, parity = None):
    for mu in range(3):
        save_sc_indices()
    expr = simplified( rsc_call(mk_4qOp_LR_cmix,"s","d","u","u",p,False,parity)
                       - sympy.simplify(1)/2* rsc_call(mk_4qOp_LR_cmix,"s","d","d","d",p,False,parity)
                       - sympy.simplify(1)/2* rsc_call(mk_4qOp_LR_cmix,"s","d","s","s",p,False,parity) )
    jump_sc_indices()
    return expr + f"Q8({p}{show_parity(parity)})"

def mk_Q9(p, parity = None):
    for mu in range(3):
        save_sc_indices()
    expr = simplified( rsc_call(mk_4qOp_LL,"s","d","u","u",p,False,parity)
                       - sympy.simplify(1)/2* rsc_call(mk_4qOp_LL,"s","d","d","d",p,False,parity)
                       - sympy.simplify(1)/2* rsc_call(mk_4qOp_LL,"s","d","s","s",p,False,parity) )
    jump_sc_indices()
    return expr + f"Q9({p}{show_parity(parity)})"

def mk_Q10(p, parity = None):
    for mu in range(3):
        save_sc_indices()
    expr = simplified( rsc_call(mk_4qOp_LL_cmix,"s","d","u","u",p,False,parity)
                       - sympy.simplify(1)/2* rsc_call(mk_4qOp_LL_cmix,"s","d","d","d",p,False,parity)
                       - sympy.simplify(1)/2* rsc_call(mk_4qOp_LL_cmix,"s","d","s","s",p,False,parity) )
    jump_sc_indices()
    return expr + f"Q10({p}{show_parity(parity)})"

# 3-flavor operators in (8,1) representation
# mk_Q1_b81
# mk_Q2_b81
# mk_Q3_b81
# mk_Q4_b81
#
# subtraction operators
# mk_Q0_b81 ( = mk_Qsub )
#
# charm-contained operators in (8,1) representation
# mk_Q5_b81
# mk_Q6_b81
# mk_Q7_b81
# mk_Q8_b81
#
# Qa^{e/o} = Aa^{e/o} Q0^{e/o} + Mai Qi^{e/o} ( i = 1, ... ,4; a = 5, ... ,8 )

def mk_Q0_b81(p, parity = None):
    return mk_Qsub(p, parity)

def mk_Q1_b81(p, parity = None):
    for mu in range (4):
        save_sc_indices()
    expr = simplified( sympy.simplify(1)/sympy.sqrt(10)* rsc_call(mk_4qOp_LL,"s","d","u","u",p,False,parity)
                       + sympy.simplify(1)/sympy.sqrt(10)* rsc_call(mk_4qOp_LL_cmix,"s","d","u","u",p,False,parity)
                       + sympy.simplify(2)/sympy.sqrt(10)* rsc_call(mk_4qOp_LL,"s","d","d","d",p,False,parity)
                       + sympy.simplify(2)/sympy.sqrt(10)* rsc_call(mk_4qOp_LL,"s","d","s","s",p,False,parity) )
    jump_sc_indices()
    return expr + f"Q1_b81({p}{show_parity(parity)})"

def mk_Q2_b81(p, parity = None):
    for mu in range (2):
        save_sc_indices()
    expr = simplified( sympy.simplify(1)/sympy.sqrt(2)* rsc_call(mk_4qOp_LL,"s","d","u","u",p,False,parity)
                       - sympy.simplify(1)/sympy.sqrt(2)* rsc_call(mk_4qOp_LL_cmix,"s","d","u","u",p,False,parity) )
    jump_sc_indices()
    return expr + f"Q2_b81({p}{show_parity(parity)})"

def mk_Q3_b81(p, parity = None):
    for mu in range (4):
        save_sc_indices()
    expr = simplified( sympy.simplify(1)/sympy.sqrt(3)* rsc_call(mk_4qOp_LR,"s","d","u","u",p,False,parity)
                       + sympy.simplify(1)/sympy.sqrt(3)* rsc_call(mk_4qOp_LR,"s","d","d","d",p,False,parity)
                       + sympy.simplify(1)/sympy.sqrt(3)* rsc_call(mk_4qOp_LR,"s","d","s","s",p,False,parity) )
    jump_sc_indices()
    return expr + f"Q3_b81({p}{show_parity(parity)})"

def mk_Q4_b81(p, parity = None):
    for mu in range (4):
        save_sc_indices()
    expr = simplified( sympy.simplify(1)/sympy.sqrt(3)* rsc_call(mk_4qOp_LR_cmix,"s","d","u","u",p,False,parity)
                       + sympy.simplify(1)/sympy.sqrt(3)* rsc_call(mk_4qOp_LR_cmix,"s","d","d","d",p,False,parity)
                       + sympy.simplify(1)/sympy.sqrt(3)* rsc_call(mk_4qOp_LR_cmix,"s","d","s","s",p,False,parity) )
    jump_sc_indices()
    return expr + f"Q4_b81({p}{show_parity(parity)})"

def mk_Q5_b81(p, parity = None):
    return mk_4qOp_LL("s","d","c","c",p,False,parity) + f"Q5_b81({p}{show_parity(parity)})"

def mk_Q6_b81(p, parity = None):
    return mk_4qOp_LL_cmix("s","d","c","c",p,False,parity) + f"Q6_b81({p}{show_parity(parity)})"

def mk_Q7_b81(p, parity = None):
    return mk_4qOp_LR("s","d","c","c",p,False,parity) + f"Q7_b81({p}{show_parity(parity)})"

def mk_Q8_b81(p, parity = None):
    return mk_4qOp_LR_cmix("s","d","c","c",p,False,parity) + f"Q8_b81({p}{show_parity(parity)})"

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
    print(display_cexpr(contract_simplify_compile(expr1, expr2, expr3)))
    # print(display_cexpr(contract_simplify_compile(
    #     mk_pipi_i0("x1_1", "x1_2", True) * mk_Q1("x") * mk_k_0("x2"),
    #     mk_pipi_i0("x1_1", "x1_2", True) * mk_Q2("x") * mk_k_0("x2"),
    #     mk_pipi_i0("x1_1", "x1_2", True) * mk_Q3("x") * mk_k_0("x2"),
    #     mk_pipi_i0("x1_1", "x1_2", True) * mk_Q4("x") * mk_k_0("x2"),
    #     mk_pipi_i0("x1_1", "x1_2", True) * mk_Q5("x") * mk_k_0("x2"),
    #     mk_pipi_i0("x1_1", "x1_2", True) * mk_Q6("x") * mk_k_0("x2"),
    #     mk_pipi_i0("x1_1", "x1_2", True) * mk_Q7("x") * mk_k_0("x2"),
    #     mk_pipi_i0("x1_1", "x1_2", True) * mk_Q8("x") * mk_k_0("x2"),
    #     mk_pipi_i0("x1_1", "x1_2", True) * mk_Q9("x") * mk_k_0("x2"),
    #     mk_pipi_i0("x1_1", "x1_2", True) * mk_Q10("x") * mk_k_0("x2"),
    #     )))

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
    expr1 = sympy.simplify(1)*1/4 * (A("x2_j", "x2_p", True) * A("x1_j", "x1_p") + A("x2_j", "x2_p") * A("x1_j", "x1_p", True))
    expr2 = sympy.simplify(1)*1/12 * (B("x2_j", "x2_p", True) * B("x1_j", "x1_p") + B("x2_j", "x2_p") * B("x1_j", "x1_p", True))
    expr3 = sympy.simplify(1)*1/4 * (C("x2_j", "x2_p", True) * C("x1_j", "x1_p") + C("x2_j", "x2_p") * C("x1_j", "x1_p", True))
    expr4 = sympy.simplify(1)*1/4 * (D("x2_j", "x2_p", True) * D("x1_j", "x1_p") + D("x2_j", "x2_p") * D("x1_j", "x1_p", True))
    expr5 = sympy.simplify(1)*1/6 * (E("x2_j", "x2_p", True) * E("x1_j", "x1_p") + E("x2_j", "x2_p") * E("x1_j", "x1_p", True))
    expr6 = sympy.simplify(1)*1/2 * (F("x2_j", "x2_p", True) * F("x1_j", "x1_p") + F("x2_j", "x2_p") * F("x1_j", "x1_p", True))
    expr7 = sympy.simplify(1)*1/2 * (G("x2_j", "x2_p", True) * G("x1_j", "x1_p") + G("x2_j", "x2_p") * G("x1_j", "x1_p", True))
    expr8 = sympy.simplify(1)*1/(4*sympy.sqrt(2)) * (C("x2_j", "x2_p", True) * F("x1_j", "x1_p") + F("x2_j", "x2_p", True) * C("x1_j", "x1_p")
            + C("x2_j", "x2_p") * F("x1_j", "x1_p", True) + F("x2_j", "x2_p") * C("x1_j", "x1_p", True))
    expr9 = sympy.simplify(1)*1/(4*sympy.sqrt(2)) * (D("x2_j", "x2_p", True) * G("x1_j", "x1_p") + G("x2_j", "x2_p", True) * D("x1_j", "x1_p")
            + D("x2_j", "x2_p") * G("x1_j", "x1_p", True) + G("x2_j", "x2_p") * D("x1_j", "x1_p", True))
    expr_p = sympy.simplify(1)*1/2*(jpi_p("x2_j", "x2_p", True) * jpi_p("x1_j", "x1_p") + jpi_p("x2_j", "x2_p") * jpi_p("x1_j", "x1_p", True))
    expr_0 = sympy.simplify(1)*1/2*(jpi_0("x2_j", "x2_p", True) * jpi_0("x1_j", "x1_p") + jpi_0("x2_j", "x2_p") * jpi_0("x1_j", "x1_p", True))
    exprs = [expr1, expr1 - expr2, expr3, expr3 - expr4, expr5, expr6, expr6 - expr7, expr8, expr8 - expr9,]
    print(display_cexpr(contract_simplify_compile(*exprs)))
    exprs1 = [
            expr_p - sympy.simplify(1)*1/18*(9/2*expr1 + 9/2*expr3 + expr6 + 3*sympy.sqrt(2)*expr8),
            expr_0 - sympy.simplify(1)*1/18*(6*expr2 + 3*expr5 + expr7),
            (expr_p - expr_0),
            ]
    print(display_cexpr(contract_simplify_compile(*exprs1)))

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
        cexpr = contract_simplify_compile(*exprs)
        print(display_cexpr(cexpr))
        cexpr.collect_op()
        print(display_cexpr(cexpr))
        print(f"-- {name} --\n")

def test_kk_sym():
    expr1 = mk_kk_i0("x2_1", "x2_2", True, is_sym = True) * mk_kk_i0("x1_1", "x1_2", is_sym = True)
    expr2 = mk_kk_i0("x2_1", "x2_2", True, is_sym = True) * mk_pipi_i0("x1_1", "x1_2")
    expr3 = mk_kk_i0("x2_1", "x2_2", True, is_sym = True) * mk_sigma("x1", "x1")
    expr4 = mk_k0k0bar("x2_1", "x2_2", True, is_sym = True) * mk_k0k0bar("x1_1", "x1_2", is_sym = True)
    expr5 = mk_k0k0bar("x2_1", "x2_2", True, is_sym = True) * mk_pipi_i0("x1_1", "x1_2")
    expr6 = mk_k0k0bar("x2_1", "x2_2", True, is_sym = True) * mk_sigma("x1", "x1")
    all_exprs = [
            [ expr1, expr4, ],
            [ expr2, expr5, ],
            [ expr3, expr6, ],
            ]
    names = [ "kk kk", "kk pipi", "kk sigma", ]
    for name, exprs in zip(names, all_exprs):
        print(f"\n-- {name} --")
        cexpr = contract_simplify_compile(*exprs)
        print(display_cexpr(cexpr))
        cexpr.collect_op()
        print(display_cexpr(cexpr))
        print(f"-- {name} --\n")

def test_kpipi():
    s1 = new_spin_index()
    s2 = new_spin_index()
    c = new_color_index()
    expr1 = Qb("s", "x", s1, c) * G("G1", s1, s2) * Qv("d", "x", s2, c) + "sb_G1_d"
    s1 = new_spin_index()
    s2 = new_spin_index()
    c = new_color_index()
    expr2 = Qb("u", "x", s1, c) * G("G2", s1, s2) * Qv("u", "x", s2, c) + "ub_G2_u"
    s1 = new_spin_index()
    s2 = new_spin_index()
    c = new_color_index()
    expr3 = Qb("d", "y", s1, c) * G("5", s1, s2) * Qv("s", "y", s2, c) + "db_g5_s"
    print(expr1)
    print(expr2)
    print(expr3)
    print(display_cexpr(contract_simplify_compile(expr1 * expr2 * expr3, is_isospin_symmetric_limit = False)))

def test_prop():
    s1 = new_spin_index()
    s2 = new_spin_index()
    c1 = new_color_index()
    c2 = new_color_index()
    expr = Qv("q", "x", s1, c1) * Qb("q", "y", s2, c2) + "q(x) qb(y)"
    print(expr)
    print(display_cexpr(contract_simplify_compile(expr, is_isospin_symmetric_limit = False)))

def test_pipi():
    expr1 = mk_pi_0("x2_1", True) * mk_pi_0("x2_2", True) * mk_pi_0("x1_1") * mk_pi_0("x1_2")
    expr2 = mk_pi_p("x2_1", True) * mk_pi_m("x2_2", True) * mk_pi_p("x1_1") * mk_pi_m("x1_2")
    expr2p = mk_pi_m("x2_1", True) * mk_pi_p("x2_2", True) * mk_pi_p("x1_1") * mk_pi_m("x1_2")
    expr = (expr2 + expr2p) - expr1
    print(expr)
    print(display_cexpr(contract_simplify_compile(expr, is_isospin_symmetric_limit = True)))

if __name__ == "__main__":
    # test_pipi()
    # test_prop()
    # test_kpipi()
    # test1()
    # exit()
    #
    print("pi+(x1):\n", mk_pi_p("x1"))
    print("pi+(x2)^dag pi+(x1):\n", mk_pi_p("x2", True) * mk_pi_p("x1"))
    print("< pi+(x2)^dag pi+(x1) >:")
    expr = mk_pi_p("x2", True) * mk_pi_p("x1")
    print(display_cexpr(contract_simplify_compile(expr)))
    print()
    print("pi0(x1):\n", mk_pi_0("x1"))
    print("pi0(x2)^dag pi0(x1):\n", mk_pi_0("x2", True) * mk_pi_0("x1"))
    print("< pi0(x2)^dag pi0(x1) >:")
    expr = mk_pi_0("x2", True) * mk_pi_0("x1")
    print(display_cexpr(contract_simplify_compile(expr, is_isospin_symmetric_limit = False)))
    print(display_cexpr(contract_simplify_compile(expr)))
    print()
    print("< pi0(x2)^dag pi0(x1) > - < pi+(x2)^dag pi+(x1) >: ",
            simplified(contract_expr(mk_pi_0("x2", True) * mk_pi_0("x1") - mk_pi_p("x2", True) * mk_pi_p("x1")), is_isospin_symmetric_limit = True))
    print()
    print("< pipiI22(x2_1,x2_2)^dag pipiI22(x1_1,x1_2) >:")
    expr = mk_pipi_i22("x2_1", "x2_2", True) * mk_pipi_i22("x1_1", "x1_2")
    print(display_cexpr(contract_simplify_compile(expr)))
    print("< pipiI21(x2_1,x2_2)^dag pipiI21(x1_1,x1_2) - pipiI22(x2_1,x2_2)^dag pipiI22(x1_1,x1_2) >:")
    expr = mk_pipi_i21("x2_1", "x2_2", True) * mk_pipi_i21("x1_1", "x1_2") - mk_pipi_i22("x2_1", "x2_2", True) * mk_pipi_i22("x1_1", "x1_2")
    print(display_cexpr(contract_simplify_compile(expr)))
    print("< pipiI20(x2_1,x2_2)^dag pipiI20(x1_1,x1_2) - pipiI22(x2_1,x2_2)^dag pipiI22(x1_1,x1_2) >:")
    expr = mk_pipi_i20("x2_1", "x2_2", True) * mk_pipi_i20("x1_1", "x1_2") - mk_pipi_i22("x2_1", "x2_2", True) * mk_pipi_i22("x1_1", "x1_2")
    print(display_cexpr(contract_simplify_compile(expr)))
    print()
    print("< pipiI11(x2_1,x2_2)^dag pipiI11(x1_1,x1_2) >:")
    expr = mk_pipi_i11("x2_1", "x2_2", True) * mk_pipi_i11("x1_1", "x1_2")
    print(display_cexpr(contract_simplify_compile(expr)))
    print("< pipiI10(x2_1,x2_2)^dag pipiI10(x1_1,x1_2) - pipiI11(x2_1,x2_2)^dag pipiI11(x1_1,x1_2) >:")
    expr = mk_pipi_i10("x2_1", "x2_2", True) * mk_pipi_i10("x1_1", "x1_2") - mk_pipi_i11("x2_1", "x2_2", True) * mk_pipi_i11("x1_1", "x1_2")
    print(display_cexpr(contract_simplify_compile(expr)))
    print()
    print("< pipiI0(x1_1,x1_2) >:")
    expr = mk_pipi_i0("x1_1", "x1_2")
    print(display_cexpr(contract_simplify_compile(expr)))
    print("< pipiI0(x2_1,x2_2)^dag pipiI0(x1_1,x1_2) >:")
    expr = mk_pipi_i0("x2_1", "x2_2", True) * mk_pipi_i0("x1_1", "x1_2")
    print(display_cexpr(contract_simplify_compile(expr)))
    print()
    print("< sigma(x1) >:")
    expr = mk_sigma("x1")
    print(display_cexpr(contract_simplify_compile(expr)))
    print("< pipiI0(x2_1,x2_2)^dag sigma(x1) >:")
    expr = mk_pipi_i0("x2_1", "x2_2", True) * mk_sigma("x1")
    print(display_cexpr(contract_simplify_compile(expr)))
    print()
    print("< K+(x2)^dag K+(x1)>:")
    expr = mk_k_p("x2", True) * mk_k_p("x1")
    print(display_cexpr(contract_simplify_compile(expr)))
    print()
    print("< pi0(x2)^dag j_mu(xj_1) j_nu(xj_2) pi0(x1) >:")
    expr = mk_pi_0("x1", True) * mk_j_mu("xj_1", "mu") * mk_j_mu("xj_2", "nu") * mk_pi_0("x2")
    print(display_cexpr(contract_simplify_compile(expr)))
    print()
    print("< pi+(x2)^dag j_mu(xj_1) j_nu(xj_2) pi+(x1) / 2 + pi-(x2)^dag j_mu(xj_1) j_nu(xj_2) pi-(x1) / 2 >:")
    expr = (
            sympy.simplify(1)/2 * mk_pi_p("x2", True) * mk_j_mu("xj_1", "mu") * mk_j_mu("xj_2", "nu") * mk_pi_p("x1")
            + sympy.simplify(1)/2 * mk_pi_m("x2", True) * mk_j_mu("xj_1", "mu") * mk_j_mu("xj_2", "nu") * mk_pi_m("x1")
            )
    print(display_cexpr(contract_simplify_compile(expr)))
    print()
    print("< pi+(x2)^dag j_mu(xj_1) j_nu(xj_2) pi+(x1) / 2 + pi-(x2)^dag j_mu(xj_1) j_nu(xj_2) pi-(x1) / 2 - pi0(x2)^dag j_mu(xj_1) j_nu(xj_2) pi0(x1) >:")
    expr = (
            sympy.simplify(1)/2 * mk_pi_p("x2", True) * mk_j_mu("xj_1", "mu") * mk_j_mu("xj_2", "nu") * mk_pi_p("x1")
            + sympy.simplify(1)/2 * mk_pi_m("x2", True) * mk_j_mu("xj_1", "mu") * mk_j_mu("xj_2", "nu") * mk_pi_m("x1")
            - mk_pi_0("x1", True) * mk_j_mu("xj_1", "mu") * mk_j_mu("xj_2", "nu") * mk_pi_0("x2")
            )
    print(display_cexpr(contract_simplify_compile(expr)))
    #
    print()
    test_kk()
    print()
    test_kk_sym()
    #
    print()
    test()
    print("< KpiI3/2_Iz1/2(x2_1,x2_2)^dag KpiI1/2_Iz1/2(x1_1,x1_2) >:")
    expr = mk_kpi_0_i3halves("x2_1", "x2_2", True) * mk_kpi_0_i1half("x1_1", "x1_2")
    print(display_cexpr(contract_simplify_compile(expr)))
    print()
    print("< KpiI3/2_Iz3/2(x2_1,x2_2)^dag KpiI3/2_Iz3/2(x1_1,x1_2) > - KpiI3/2_Iz1/2(x2_1,x2_2)^dag KpiI3/2_Iz1/2(x1_1,x1_2) >:")
    expr = mk_kpi_p2_i3halves("x2_1", "x2_2", True) * mk_kpi_p2_i3halves("x1_1", "x1_2") - mk_kpi_p1_i3halves("x2_1", "x2_2", True) * mk_kpi_p1_i3halves("x1_1", "x1_2")
    print(display_cexpr(contract_simplify_compile(expr)))
    print()
    print("< KpiI3/2_Iz3/2(x2_1,x2_2)^dag KpiI3/2_Iz3/2(x1_1,x1_2) > - KpiI3/2_Iz-1/2(x2_1,x2_2)^dag KpiI3/2_Iz-1/2(x1_1,x1_2) >:")
    expr = mk_kpi_p2_i3halves("x2_1", "x2_2", True) * mk_kpi_p2_i3halves("x1_1", "x1_2") - mk_kpi_0_i3halves("x2_1", "x2_2", True) * mk_kpi_0_i3halves("x1_1", "x1_2")
    print(display_cexpr(contract_simplify_compile(expr)))
    print()
    print("< KpiI3/2_Iz3/2(x2_1,x2_2)^dag KpiI3/2_Iz3/2(x1_1,x1_2) > - KpiI3/2_Iz-3/2(x2_1,x2_2)^dag KpiI3/2_Iz-3/2(x1_1,x1_2) >:")
    expr = mk_kpi_p2_i3halves("x2_1", "x2_2", True) * mk_kpi_p2_i3halves("x1_1", "x1_2") - mk_kpi_m_i3halves("x2_1", "x2_2", True) * mk_kpi_m_i3halves("x1_1", "x1_2")
    print(display_cexpr(contract_simplify_compile(expr)))
    print()
    print("< KpiI1/2_Iz1/2(x2_1,x2_2)^dag KpiI1/2_Iz1/2(x1_1,x1_2) > - < KpiI1/2_Iz-1/2(x2_1,x2_2)^dag KpiI1/2_Iz-1/2(x1_1,x1_2) >:")
    expr = mk_kpi_0_i1half("x2_1", "x2_2", True) * mk_kpi_0_i1half("x1_1", "x1_2") - mk_kpi_p_i1half("x2_1", "x2_2", True) * mk_kpi_p_i1half("x1_1", "x1_2")
    print(display_cexpr(contract_simplify_compile(expr)))
    print()


