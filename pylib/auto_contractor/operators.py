from auto_contract.wick import *

import math

spin_index_counter = 0

color_index_counter = 0

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
    return Qb(f1, p, s1, c) * G("5", s1, s2) * Qv(f2, p, s2, c)

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

if __name__ == "__main__":
    #
    print("pi+(x1):\n", mk_pi_p("x1"))
    print("pi+(x2)^dag pi+(x1):\n", (mk_pi_p("x2", True) * mk_pi_p("x1")).round())
    print("< pi+(x2)^dag pi+(x1) >:\n", simplified(contract_expr(mk_pi_p("x2", True) * mk_pi_p("x1")), is_isospin_symmetric_limit = True).round())
    print()
    print("pi0(x1):\n", mk_pi_0("x1"))
    print("pi0(x2)^dag pi0(x1):\n", (mk_pi_0("x2", True) * mk_pi_0("x1")).round())
    print("< pi0(x2)^dag pi0(x1) >:\n", simplified(contract_expr(mk_pi_0("x2", True) * mk_pi_0("x1")), is_isospin_symmetric_limit = False).round())
    print("< pi0(x2)^dag pi0(x1) >:\n", simplified(contract_expr(mk_pi_0("x2", True) * mk_pi_0("x1")), is_isospin_symmetric_limit = True).round())
    print()
    print("< pi0(x2)^dag pi0(x1) > - < pi+(x2)^dag pi+(x1) >: ",
            simplified(contract_expr(mk_pi_0("x2", True) * mk_pi_0("x1") - mk_pi_p("x2", True) * mk_pi_p("x1")), is_isospin_symmetric_limit = True).round())
    print()
    print("< pipiI22(x21,x22)^dag pipiI22(x11,x12) >:\n",
            simplified(contract_expr(mk_pipi_i22("x21", "x22", True) * mk_pipi_i22("x11", "x12")), is_isospin_symmetric_limit = True).round())
    print("< pipiI21(x21,x22)^dag pipiI21(x11,x12) - pipiI22(x21,x22)^dag pipiI22(x11,x12) >:",
            simplified(contract_expr(mk_pipi_i21("x21", "x22", True) * mk_pipi_i21("x11", "x12") - mk_pipi_i22("x21", "x22", True) * mk_pipi_i22("x11", "x12")), is_isospin_symmetric_limit = True).round())
    print("< pipiI20(x21,x22)^dag pipiI20(x11,x12) - pipiI22(x21,x22)^dag pipiI22(x11,x12) >:",
            simplified(contract_expr(mk_pipi_i20("x21", "x22", True) * mk_pipi_i20("x11", "x12") - mk_pipi_i22("x21", "x22", True) * mk_pipi_i22("x11", "x12")), is_isospin_symmetric_limit = True).round())
    print()
    print("< pipiI11(x21,x22)^dag pipiI11(x11,x12) >:\n",
            simplified(contract_expr(mk_pipi_i11("x21", "x22", True) * mk_pipi_i11("x11", "x12")), is_isospin_symmetric_limit = True).round())
    print("< pipiI10(x21,x22)^dag pipiI10(x11,x12) - pipiI11(x21,x22)^dag pipiI11(x11,x12) >: ",
            simplified(contract_expr(mk_pipi_i10("x21", "x22", True) * mk_pipi_i10("x11", "x12") - mk_pipi_i11("x21", "x22", True) * mk_pipi_i11("x11", "x12")), is_isospin_symmetric_limit = True).round())
    print()
    print("< pipiI0(x11,x12) >:\n", simplified(contract_expr(mk_pipi_i0("x11", "x12")), is_isospin_symmetric_limit = True).round())
    print("< pipiI0(x21,x22)^dag pipiI0(x11,x12) >:\n", simplified(contract_expr(mk_pipi_i0("x21", "x22", True) * mk_pipi_i0("x11", "x12")), is_isospin_symmetric_limit = True).round())
