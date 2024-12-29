#    Qlattice (https://github.com/jinluchang/qlattice)
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

import qlat_utils as q
import copy
import sympy
import numpy as np
from typing import Self

try:
    from . import expr_arithmetic as ea
    from .expr_arithmetic import mk_sym
except:
    import expr_arithmetic as ea
    from expr_arithmetic import mk_sym

class Op:

    """
    self.otype
    """

    def __init__(self, otype:str):
        self.otype = otype

    def is_commute(self) -> bool:
        return True

    def show(self, is_multiply = False):
        return repr(self)

    def __repr__(self) -> str:
        return f"Op({self.otype})"

    def __add__(self, other):
        return mk_expr(self) + other

    __radd__ = __add__

    def __mul__(self, other):
        return mk_expr(self) * other

    def __rmul__(self, other):
        return mk_expr(other) * self

    def __neg__(self):
        return mk_expr(-1) * mk_expr(self)

    def __pos__(self):
        return self

    def __sub__(self, other):
        return mk_expr(self) + mk_expr(-1) * other

    def __rsub__(self, other):
        return mk_expr(other) + mk_expr(-1) * self

    def sort(self) -> None:
        pass

    def isospin_symmetric_limit(self) -> None:
        pass

### ------

class Qfield(Op):

    """
    self.f
    self.p
    self.s
    self.c
    """

    def __init__(self, otype:str, flavor:str, position:str, spin:str, color:str):
        Op.__init__(self, otype)
        self.f = flavor
        self.p = position
        self.s = spin
        self.c = color

    def is_commute(self) -> bool:
        return False

    def __repr__(self) -> str:
        if self.s == "auto" and self.c == "auto":
            return f"{self.otype}({self.f!r},{self.p!r})"
        else:
            return f"{self.otype}({self.f!r},{self.p!r},{self.s!r},{self.c!r})"

    def list(self):
        return [ self.otype, self.f, self.p, self.s, self.c, ]

    def __eq__(self, other) -> bool:
        return self.list() == other.list()

### ------

class Qv(Qfield):

    """
    act as d / d Hb
    """

    def __init__(self, f, p, s, c):
        """
        interface function
        """
        Qfield.__init__(self, "Qv", f, p, s, c)

### ------

class Qb(Qfield):

    """
    act as d / d Hv
    """

    def __init__(self, f, p, s, c):
        """
        interface function
        """
        Qfield.__init__(self, "Qb", f, p, s, c)

### ------

class Hv(Qfield):

    def __init__(self, f, p, s, c):
        Qfield.__init__(self, "Hv", f, p, s, c)

### ------

class Hb(Qfield):

    def __init__(self, f, p, s, c):
        Qfield.__init__(self, "Hb", f, p, s, c)

### ------

class SHv(Qfield):

    def __init__(self, f, p, s, c):
        Qfield.__init__(self, "SHv", f, p, s, c)

### ------

class HbS(Qfield):

    def __init__(self, f, p, s, c):
        Qfield.__init__(self, "HbS", f, p, s, c)

### ------

class S(Op):

    """
    propagator
    #
    self.f
    self.p1
    self.p2
    self.s1
    self.s2
    self.c1
    self.c2
    """

    def __init__(self, flavor:str, p1:str, p2:str, s1:str="auto", s2:str="auto", c1:str="auto", c2:str="auto"):
        Op.__init__(self, "S")
        self.f = flavor
        self.p1 = p1
        self.p2 = p2
        self.s1 = s1
        self.s2 = s2
        self.c1 = c1
        self.c2 = c2

    def __repr__(self) -> str:
        if self.s1 == "auto" and self.s2 == "auto" and self.c1 == "auto" and self.c2 == "auto":
            return f"{self.otype}({self.f!r},{self.p1!r},{self.p2!r})"
        else:
            return f"{self.otype}({self.f!r},{self.p1!r},{self.p2!r},{self.s1!r},{self.s2!r},{self.c1!r},{self.c2!r})"

    def list(self):
        return [ self.otype, self.f, self.p1, self.p2, self.s1, self.s2, self.c1, self.c2, ]

    def __eq__(self, other) -> bool:
        return self.list() == other.list()

    def isospin_symmetric_limit(self) -> None:
        """
        can also use these fictitious quark field to remove some unwanted disconnected diagrams
        """
        if self.f in [ "u", "d", "u'", "d'", "u''", "d''", "u'''", "d'''", ]:
            self.f = "l"
        elif self.f in [ "s", "s'", "s''", "s'''", ]:
            self.f = "s"
        elif self.f in [ "c", "c'", "c''", "c'''", ]:
            self.f = "c"

### ------

class G(Op):

    """
    spin matrix
    #
    self.tag
    self.s1
    self.s2
    #
    tag = 0, 1, 2, 3, 5 for gamma matrices
    """

    def __init__(self, tag, s1:str="auto", s2:str="auto"):
        Op.__init__(self, "G")
        self.tag = tag
        self.s1 = s1
        self.s2 = s2

    def __repr__(self) -> str:
        if self.s1 == "auto" and self.s2 == "auto":
            return f"{self.otype}({self.tag!r})"
        else:
            return f"{self.otype}({self.tag!r},{self.s1!r},{self.s2!r})"

    def list(self):
        return [ self.otype, self.tag, self.s1, self.s2, ]

    def __eq__(self, other) -> bool:
        return self.list() == other.list()

### ------

class U(Op):

    """
    color matrix
    #
    self.tag
    self.p
    self.mu
    self.c1
    self.c2
    """

    def __init__(self, tag, p, mu, c1:str="auto", c2:str="auto"):
        Op.__init__(self, "U")
        self.tag = tag
        self.p = p
        self.mu = mu
        self.c1 = c1
        self.c2 = c2

    def __repr__(self) -> str:
        if self.c1 == "auto" and self.c2 == "auto":
            return f"{self.otype}({self.tag!r},{self.p!r},{self.mu!r})"
        else:
            return f"{self.otype}({self.tag!r},{self.p!r},{self.mu!r},{self.c1!r},{self.c2!r})"

    def list(self):
        return [ self.otype, self.tag, self.p, self.mu, self.c1, self.c2, ]

    def __eq__(self, other) -> bool:
        return self.list() == other.list()

### ------

class Tr(Op):

    """
    a collection of ops taking the trace
    #
    self.ops
    self.tag
    """

    def __init__(self, ops:list, tag=None):
        Op.__init__(self, "Tr")
        if tag is not None:
            # do not perform check if tag is set
            self.ops = ops
            self.tag = tag
            return
        for op in ops:
            assert op.is_commute()
            assert op.otype in [ "S", "G", "U", ]
        s = None
        c = None
        for op in ops + ops:
            if not check_chain_sc(op, s, c):
                raise Exception(f"ops={ops} tag={tag}")
            s, c, = update_chain_sc(op, s, c)
        if s is not None and c is not None:
            self.tag = "sc"
        elif s is not None:
            self.tag = "s"
        elif c is not None:
            self.tag = "c"
        else:
            self.tag = ""
        for i, op in enumerate(ops):
            ops[i] = copy_op_index_auto(op)
        self.ops = ops

    def __repr__(self) -> str:
        return f"{self.otype}({self.ops!r},{self.tag!r})"

    def list(self):
        return [ self.otype, self.tag, self.ops, ]

    def __eq__(self, other) -> bool:
        return self.list() == other.list()

    def sort(self):
        ops = self.ops
        if len(ops) > 1:
            self.ops = sorted([ ops[i:] + ops[:i] for i in range(len(ops)) ], key=repr)[0]

    def isospin_symmetric_limit(self) -> None:
        for op in self.ops:
            op.isospin_symmetric_limit()

### ------

class Chain(Op):

    """
    a collection of ops multiplying together but do not form a loop
    #
    self.ops
    self.tag
    """

    def __init__(self, ops:list, tag=None):
        Op.__init__(self, "Chain")
        if tag is not None:
            # do not perform check if tag is set
            self.ops = ops
            self.tag = tag
            return
        for op in ops:
            assert op.is_commute()
            assert op.otype in [ "S", "G", "U", ]
        s = None
        c = None
        for op in ops:
            if not check_chain_sc(op, s, c):
                raise Exception(f"ops={ops} tag={tag}")
            s, c, = update_chain_sc(op, s, c)
        if s is not None and c is not None:
            self.tag = "sc"
        elif s is not None:
            self.tag = "s"
        elif c is not None:
            self.tag = "c"
        else:
            self.tag = ""
        n_ops = len(ops)
        assert n_ops >= 1
        for i, op in enumerate(ops):
            if i == 0:
                if i < n_ops - 1:
                    ops[i] = copy_op_index_auto(op, is_auto_sc1=False)
            elif i == n_ops - 1:
                ops[i] = copy_op_index_auto(op, is_auto_sc2=False)
            else:
                ops[i] = copy_op_index_auto(op)
        self.ops = ops

    def __repr__(self) -> str:
        return f"{self.otype}({self.ops!r},{self.tag!r})"

    def list(self):
        return [ self.otype, self.tag, self.ops, ]

    def __eq__(self, other) -> bool:
        return self.list() == other.list()

    def sort(self):
        """
        Do not need to sort
        """

    def isospin_symmetric_limit(self) -> None:
        for op in self.ops:
            op.isospin_symmetric_limit()

### ------

def copy_op_index_auto(op:Op, is_auto_sc1=True, is_auto_sc2=True):
    if op.otype not in [ "S", "G", "U", ]:
        return op
    op = copy.copy(op)
    if op.otype in [ "S", "G", ]:
        if is_auto_sc1:
            op.s1 = "auto"
        if is_auto_sc2:
            op.s2 = "auto"
    if op.otype in [ "S", "U", ]:
        if is_auto_sc1:
            op.c1 = "auto"
        if is_auto_sc2:
            op.c2 = "auto"
    return op

def check_chain_spin_index(ops:list, s:str):
    """
    for a spin index `s`,
    return is_having_repeated_index, first_index_of_operator, second_index_of_operator
    """
    count1 = 0
    count2 = 0
    i1 = None
    i2 = None
    for i, op in enumerate(ops):
        if op.otype in [ "S", "G", ]:
            if op.s1 == s:
                i1 = i
                count1 += 1
            if op.s2 == s:
                i2 = i
                count2 += 1
    assert count1 < 2
    assert count2 < 2
    return count1 == 1 and count2 == 1, i1, i2

def check_chain_color_index(ops:list, c:str):
    """
    for a color index `s`,
    return is_having_repeated_index, first_index_of_operator, second_index_of_operator
    """
    count1 = 0
    count2 = 0
    i1 = None
    i2 = None
    for i, op in enumerate(ops):
        if op.otype in [ "S", "U", ]:
            if op.c1 == c:
                i1 = i
                count1 += 1
            if op.c2 == c:
                i2 = i
                count2 += 1
    assert count1 < 2
    assert count2 < 2
    return count1 == 1 and count2 == 1, i1, i2

def check_chain_op(ops:list, op:Op):
    """
    for a operator `op`,
    return "single" or "begin" or "middle" or "end",
    if this `op` is part of a (longer) chain of contraction
    if `op` does not have right type, return `None`.
    """
    if op.otype not in [ "S", "G", "U", ]:
        return None
    b_begin = True
    b_end = True
    if op.otype in [ "S", ]:
        if check_chain_spin_index(ops, op.s1)[0] and check_chain_color_index(ops, op.c1)[0]:
            b_begin = False
        if check_chain_spin_index(ops, op.s2)[0] and check_chain_color_index(ops, op.c2)[0]:
            b_end = False
    elif op.otype in [ "G", ]:
        if check_chain_spin_index(ops, op.s1)[0]:
            b_begin = False
        if check_chain_spin_index(ops, op.s2)[0]:
            b_end = False
    elif op.otype in [ "U", ]:
        if check_chain_color_index(ops, op.c1)[0]:
            b_begin = False
        if check_chain_color_index(ops, op.c2)[0]:
            b_end = False
    else:
        assert False
    if b_begin and b_end:
        return "single"
    elif b_begin:
        return "begin"
    elif b_end:
        return "end"
    else:
        return "middle"

def check_chain_sc(op, s, c):
    if op.otype in [ "S", "G", ]:
        if s is not None and s != op.s1:
            return False
    if op.otype in [ "S", "U", ]:
        if c is not None and c != op.c1:
            return False
    return True

def update_chain_sc(op, s, c):
    if op.otype in [ "S", "G", ]:
        s = op.s2
    if op.otype in [ "S", "U", ]:
        c = op.c2
    return s, c

def pick_chain_op(ops:list, masks:list, s, c, type_list=None):
    if type_list is None:
        type_list = [ "single", "begin", "end", "middle", ]
    for i, op in enumerate(ops):
        if masks[i]:
            continue
        if not check_chain_sc(op, s, c):
            continue
        if check_chain_op(ops, op) not in type_list:
            continue
        return i, op
    return None

def find_trace(ops:list):
    """
    return None or (Tr(tr_ops), remaining_ops,)
    """
    size = len(ops)
    for i, op in enumerate(ops):
        if check_chain_op(ops, op) not in [ "middle", ]:
            continue
        masks = [ False for op in ops ]
        tr_ops = []
        s = None
        c = None
        masks[i] = True
        tr_ops.append(op)
        s, c, = update_chain_sc(op, s, c)
        while True:
            p_op = pick_chain_op(ops, masks, s, c, [ "middle", ])
            if p_op is None:
                # trace found
                return Tr(tr_ops), [ op for i, op in enumerate(ops) if not masks[i] ]
            i2, op2 = p_op
            masks[i2] = True
            tr_ops.append(op2)
            s, c, = update_chain_sc(op2, s, c)
    return None

def find_chain(ops:list):
    """
    return None or (Chain(ch_ops), remaining_ops,)
    """
    size = len(ops)
    for i, op in enumerate(ops):
        if check_chain_op(ops, op) not in [ "single", "begin", ]:
            continue
        masks = [ False for op in ops ]
        ch_ops = []
        s = None
        c = None
        masks[i] = True
        ch_ops.append(op)
        s, c, = update_chain_sc(op, s, c)
        while True:
            p_op = pick_chain_op(ops, masks, s, c, [ "single", "begin", "middle", "end", ])
            if p_op is None:
                # chain found
                return Chain(ch_ops), [ op for i, op in enumerate(ops) if not masks[i] ]
            i2, op2 = p_op
            masks[i2] = True
            ch_ops.append(op2)
            s, c, = update_chain_sc(op2, s, c)
    return None

def collect_traces(ops:list) -> list:
    """
    First collect all the `Chain`s, then `Tr`s
    """
    chs = []
    while True:
        fc = find_chain(ops)
        if fc is None:
            break
        ch, ops = fc
        chs.append(ch)
    trs = []
    while True:
        ft = find_trace(ops)
        if ft is None:
            break
        tr, ops = ft
        trs.append(tr)
    return chs + trs + ops

### ------

class BfieldCoef:

    """
    self.coef_dict
    #
    self.coef_dict[chiral_projection] = spin_coef_dict
    spin_coef_dict[spin] = spin_coef
    #
    chiral_projection = ((1, 1,), (1, 0,), (1, 0,),)
    spin = (0, 1, 0,)
    spin_coef = 1
    """

    def __init__(self):
        self.coef_dict = dict()

    def get_spin_tensor(self, permute=None) -> np.ndarray:
        """
        return spin_tensor
        spin_tensor[s1, s2, s3] = spin_coef
        Euclidean convention, not the code convention
        """
        shape = (4, 4, 4,)
        spin_tensor = np.zeros(shape, dtype=object)
        for k, spin_coef_dict in self.coef_dict.items():
            ((c00, c01,), (c10, c11,), (c20, c21,),) = k
            for sk, spin_coef in spin_coef_dict.items():
                s1, s2, s3, = sk
                spin_tensor[s1, s2, s3] += c00 * c10 * c20 * spin_coef
                spin_tensor[s1, s2, s3 + 2] += c00 * c10 * c21 * spin_coef
                spin_tensor[s1, s2 + 2, s3] += c00 * c11 * c20 * spin_coef
                spin_tensor[s1, s2 + 2, s3 + 2] += c00 * c11 * c21 * spin_coef
                spin_tensor[s1 + 2, s2, s3] += c01 * c10 * c20 * spin_coef
                spin_tensor[s1 + 2, s2, s3 + 2] += c01 * c10 * c21 * spin_coef
                spin_tensor[s1 + 2, s2 + 2, s3] += c01 * c11 * c20 * spin_coef
                spin_tensor[s1 + 2, s2 + 2, s3 + 2] += c01 * c11 * c21 * spin_coef
        spin_tensor = q.epsilon_tensor(*permute) * spin_tensor.transpose(permute)
        return spin_tensor

    def get_spin_tensor_code(self, permute=None) -> np.ndarray:
        """
        return spin_tensor_code
        spin_tensor_code[s1, s2, s3] = spin_coef
        CPS/Grid/GPT code convention, not the Euclidean convention
        #
        psi^Eucl_0 = -1j * psi^Code_1
        psi^Eucl_1 = +1j * psi^Code_0
        psi^Eucl_2 = -1j * psi^Code_3
        psi^Eucl_3 = +1j * psi^Code_2
        """
        shape = (4, 4, 4,)
        spin_tensor_code = np.zeros(shape, dtype=object)
        sidx_arr = [
                (1, -sympy.I),
                (0, sympy.I),
                (3, -sympy.I),
                (2, sympy.I),
                ]
        spin_tensor = self.get_spin_tensor(permute=permute)
        assert spin_tensor.shape == shape
        for s1 in range(4):
            s1p, fac1, = sidx_arr[s1]
            for s2 in range(4):
                s2p, fac2, = sidx_arr[s2]
                for s3 in range(4):
                    s3p, fac3, = sidx_arr[s3]
                    spin_tensor_code[s1p, s2p, s3p] = fac1 * fac2 * fac3 * spin_tensor[s1, s2, s3]
        return spin_tensor

    def add(self, chiral_projection, spin, coef) -> Self:
        """
        return self
        Allow keep add.
        """
        if mk_sym(coef) == 0:
            return self
        ((c00, c01,), (c10, c11,), (c20, c21,),) = chiral_projection
        chiral_projection = (
                (mk_sym(c00), mk_sym(c01),),
                (mk_sym(c10), mk_sym(c11),),
                (mk_sym(c20), mk_sym(c21),),
                )
        if chiral_projection not in self.coef_dict:
            self.coef_dict[chiral_projection] = dict()
        spin_coef_dict = self.coef_dict[chiral_projection]
        s1, s2, s3, = spin
        assert isinstance(s1, int) and 0 <= s1 and s1 <= 1
        assert isinstance(s2, int) and 0 <= s2 and s2 <= 1
        assert isinstance(s3, int) and 0 <= s3 and s3 <= 1
        if spin not in spin_coef_dict:
            spin_coef_dict[spin] = coef
        else:
            v = spin_coef_dict[spin] + coef
            if v == 0:
                spin_coef_dict.pop(spin)
            else:
                spin_coef_dict[spin] = v
        return self

    def __repr__(self):
        return f"BfieldCoef()"

### ------

bfield_tag_dict = dict()

bfield_tag_dict["std-u"] = BfieldCoef()
bfield_tag_dict["std-u"].add(((1, 1,), (1, 0,), (1, 0,),), (0, 1, 0,), mk_sym(1))
bfield_tag_dict["std-u"].add(((1, 1,), (1, 0,), (1, 0,),), (0, 0, 1,), mk_sym(-1))
bfield_tag_dict["std-u"].add(((1, 1,), (0, 1,), (0, 1,),), (0, 1, 0,), mk_sym(1))
bfield_tag_dict["std-u"].add(((1, 1,), (0, 1,), (0, 1,),), (0, 0, 1,), mk_sym(-1))

bfield_tag_dict["std-d"] = BfieldCoef()
bfield_tag_dict["std-d"].add(((1, 1,), (1, 0,), (1, 0,),), (1, 1, 0,), mk_sym(1))
bfield_tag_dict["std-d"].add(((1, 1,), (1, 0,), (1, 0,),), (1, 0, 1,), mk_sym(-1))
bfield_tag_dict["std-d"].add(((1, 1,), (0, 1,), (0, 1,),), (1, 1, 0,), mk_sym(1))
bfield_tag_dict["std-d"].add(((1, 1,), (0, 1,), (0, 1,),), (1, 0, 1,), mk_sym(-1))

### ------

class Bfield(Op):

    """
    baryon tensor
    #
    self.tag
    self.s1
    self.s2
    self.s3
    self.c1
    self.c2
    self.c3
    #
    tag in bfield_tag_dict
    #
    The tensor is assumed to take the form:
        baryon(s1,s2,s3,c1,c2,c3) = spin_tensor[s1,s2,s3] * q.epsilon_tensor(c1,c2,c3)
        spin_tensor = bfield_tag_dict[self.tag].get_spin_tensor_code()
    """

    def __init__(self, tag:str, s1:str, s2:str, s3:str, c1:str, c2:str, c3:str):
        assert tag in bfield_tag_dict
        Op.__init__(self, "Bfield")
        self.tag = tag
        self.s1 = s1
        self.s2 = s2
        self.s3 = s3
        self.c1 = c1
        self.c2 = c2
        self.c3 = c3

    def __repr__(self) -> str:
        return f"{self.otype}({self.tag!r},{self.s1!r},{self.s2!r},{self.s3!r},{self.c1!r},{self.c2!r},{self.c3!r})"

    def list(self):
        return [ self.otype, self.tag, self.s1, self.s2, self.s3, self.c1, self.c2, self.c3, ]

    def __eq__(self, other) -> bool:
        return self.list() == other.list()

### ------

class BS(Op):

    """
    single baryon prop
    #
    self.tag_pair_list
    self.chain_list
    #
    tag_pair_list = [ (tag_v, permute_v, tag_b, permute_b,), ... ]
    chain_list = [ prop_0, prop_1, prop_2, ]
    tag_v or tag_b in `bfield_tag_dict`
    permute_v = (spin_color_index_that_prop_0_contract_with,
                 spin_color_index_that_prop_1_contract_with,
                 spin_color_index_that_prop_2_contract_with,)
    """

    def __init__(self, tag_pair_list:list, chain_list:list):
        assert isinstance(tag_pair_list, list)
        for v in tag_pair_list:
            (tag_v, permute_v, tag_b, permute_b,) = v
            assert isinstance(tag_v, str)
            assert isinstance(tag_b, str)
            assert tag_v in bfield_tag_dict
            assert tag_b in bfield_tag_dict
            assert isinstance(permute_v, tuple)
            assert isinstance(permute_b, tuple)
            assert len(permute_v) == 3
            assert len(permute_b) == 3
            for i in permute_v:
                assert 0 <= i and i < 3
            for i in permute_b:
                assert 0 <= i and i < 3
        assert isinstance(chain_list, list)
        assert len(chain_list) == 3
        for ch in chain_list:
            assert isinstance(ch, Chain)
            assert ch.otype == "Chain"
        Op.__init__(self, "BS")
        self.tag_pair_list = tag_pair_list
        self.prop_list = prop_list

    def sort(self):
        sp_chain_list = sorted(list(enumerate(self.chain_list)), key=lambda x: repr(x[1]))[0]
        s_chain_list = []
        i_list = []
        for sp in sp_chain_list:
            i, ch = sp
            i_list.append(i)
            s_chain_list(ch)
        def permute_permute(p):
            return tuple(p[i] for i in i_list)
        s_tag_pair_list = []
        for tag_pair in self.tag_pair_list:
            (tag_v, permute_v, tag_b, permute_b,) = tag_pair
            s_tag_pair = (tag_v, permute_permute(permute_v), tag_b, permute_permute(permute_b),)
            s_tag_pair_list.append(s_tag_pair)
        self.chain_list = s_chain_list
        self.tag_pair_list = sorted(s_tag_pair_list, key=repr)

    def isospin_symmetric_limit(self) -> None:
        for op in self.chain_list:
            op.isospin_symmetric_limit()

### ------

class Term:

    """
    self.coef
    self.c_ops
    self.a_ops
    """

    def __init__(self, c_ops, a_ops, coef = 1):
        self.coef = coef
        self.c_ops = c_ops
        self.a_ops = a_ops
        for op in c_ops:
            assert op.is_commute()
        for op in a_ops:
            assert not op.is_commute()

    def list(self):
        return [ self.coef, self.c_ops, self.a_ops, ]

    def __eq__(self, other) -> bool:
        return self.list() == other.list()

    def check_commute(self) -> bool:
        for op in self.c_ops:
            assert op.is_commute()
        for op in self.a_ops:
            assert not op.is_commute()

    def __imul__(self, factor):
        self.coef *= factor

    def __add__(self, other):
        return mk_expr(self) + other

    def __add__(self, other):
        return mk_expr(self) + other

    __radd__ = __add__

    def __mul__(self, other):
        return mk_expr(self) * other

    def __rmul__(self, other):
        return mk_expr(other) * self

    def __neg__(self):
        return mk_expr(-1) * mk_expr(self)

    def __pos__(self):
        return self

    def __sub__(self, other):
        return mk_expr(self) + mk_expr(-1) * other

    def __rsub__(self, other):
        return mk_expr(other) + mk_expr(-1) * self

    def show(self, is_multiply = False) -> str:
        return "*".join([ f"({self.coef})", ] + self.c_ops + self.a_ops)

    def __repr__(self) -> str:
        return f"Term({self.c_ops},{self.a_ops},{self.coef})"

    def sort(self) -> None:
        # only sort commutable factors
        for op in self.c_ops:
            op.sort()
        self.c_ops.sort(key=repr)

    def simplify_coef(self) -> None:
        self.coef = ea.coef_simplified(self.coef)

    def simplify_ea(self) -> None:
        self.coef = ea.simplified(self.coef)

    def collect_traces(self) -> None:
        if len(self.a_ops) == 0:
            self.c_ops = collect_traces(self.c_ops)

    def isospin_symmetric_limit(self) -> None:
        for op in self.c_ops:
            op.isospin_symmetric_limit()

### ------

def combine_two_terms(t1:Term, t2:Term, t1_sig:str, t2_sig:str):
    if t1_sig == t2_sig:
        coef = t1.coef + t2.coef
        if ea.is_zero(coef):
            return Term([], [], 0)
        else:
            return Term(t1.c_ops, t1.a_ops, coef)
    else:
        return None

class Expr:

    # self.description
    # self.terms

    def __init__(self, terms, description = None):
        self.description = description
        self.terms = terms

    @q.timer
    def copy(self):
        """
        return a deep copy of this object.
        """
        return copy.deepcopy(self)

    def __imul__(self, factor):
        for term in self.terms:
            term *= factor
        description = f"{factor} * {self.show(True)}"

    def __add__(self, other):
        # if other is str, then it is used to set the description of the result expr
        if isinstance(other, str):
            return Expr(self.terms, other)
        # otherwise return self + other
        other = mk_expr(other)
        terms = self.terms + other.terms
        return Expr(terms, f"+{self.show()} + {other.show()}")

    __radd__ = __add__

    def __mul__(self, other):
        other = mk_expr(other)
        terms = []
        for t1 in self.terms:
            for t2 in other.terms:
                coef = t1.coef * t2.coef
                if coef == 0:
                    continue
                t = Term(t1.c_ops + t2.c_ops, t1.a_ops + t2.a_ops, coef)
                terms.append(t)
        return Expr(terms, f"{self.show(True)} * {other.show(True)}")

    def __rmul__(self, other):
        return mk_expr(other) * self

    def __neg__(self):
        return mk_expr(-1) * mk_expr(self)

    def __pos__(self):
        return self

    def __sub__(self, other):
        return mk_expr(self) + mk_expr(-1) * other

    def __rsub__(self, other):
        return mk_expr(other) + mk_expr(-1) * self

    def show(self, is_multiply = False):
        if self.description is not None:
            assert isinstance(self.description, str)
            if len(self.description) >= 1 and self.description[0] == "+":
                if is_multiply:
                    return f"( {self.description[1:]} )"
                else:
                    return f"{self.description[1:]}"
            else:
                return self.description
        else:
            return repr(self)

    def __repr__(self) -> str:
        if self.description is None:
            return f"Expr({self.terms})"
        else:
            return f"Expr({self.terms},{self.description!r})"

    @q.timer
    def sort(self) -> None:
        for term in self.terms:
            term.sort()
        self.terms.sort(key=repr)

    @q.timer
    def simplify_coef(self) -> None:
        for t in self.terms:
            t.simplify_coef()
        self.drop_zeros()

    @q.timer
    def simplify_ea(self) -> None:
        for t in self.terms:
            t.simplify_ea()
        self.drop_zeros()

    @q.timer
    def combine_terms(self) -> None:
        self.terms = combine_terms_expr(self).terms

    @q.timer
    def drop_zeros(self) -> None:
        self.terms = drop_zero_terms(self).terms

    def round(self, ndigit:int=20) -> None:
        """
        interface function
        """
        sexpr = self.copy()
        for term in sexpr.terms:
            coef = term.coef
            term.coef = term.coef.evalf(ndigit)
        sexpr.terms = drop_zero_terms(sexpr).terms
        return sexpr

    @q.timer
    def collect_traces(self) -> None:
        for term in self.terms:
            term.collect_traces()

    @q.timer
    def isospin_symmetric_limit(self) -> None:
        for term in self.terms:
            term.isospin_symmetric_limit()

    @q.timer
    def simplify(self, *, is_isospin_symmetric_limit:bool=True) -> None:
        """
        interface function
        """
        if is_isospin_symmetric_limit:
            self.isospin_symmetric_limit()
        self.sort()
        self.combine_terms()
        self.collect_traces()
        self.sort()
        self.combine_terms()
        self.simplify_ea()

### ------

def mk_fac(x) -> Expr:
    """
    interface function
    Stand for "make factor", the result of this function can be used in auto contractor as a factor.
    Make an Expr obj (can be sympy expression).
    `x` can have type `str`, which will be viewed as code segment.
    The code segment can use functions and variables defined in `auto_contractor.auto_fac_funcs`, `position_dict`, `base_position_dict`.
    You can define functions in `position_dict` or `base_position_dict`.
    `position_dict` is argument in function `eval_cexpr`.
    `base_position_dict` is argument in function `cache_compiled_cexpr`.
    """
    return mk_expr(ea.mk_fac(x))

def simplified(expr:Expr, *, is_isospin_symmetric_limit:bool=True) -> Expr:
    """
    interface function
    does not change expr
    """
    sexpr = expr.copy()
    sexpr.simplify(is_isospin_symmetric_limit = is_isospin_symmetric_limit)
    return sexpr

def mk_expr(x) -> Expr:
    if isinstance(x, Op):
        if x.is_commute():
            return Expr([Term([x,], [], 1),], f"{x}")
        else:
            return Expr([Term([], [x,], 1),], f"{x}")
    elif isinstance(x, Term):
        return Expr([x,], f"x.show()")
    elif isinstance(x, Expr):
        return x
    elif isinstance(x, (int, float, complex, sympy.Basic, ea.Expr)):
        return Expr([Term([], [], x),], f"({x})")
    else:
        print(x)
        assert False

@q.timer
def combine_terms_expr(expr:Expr) -> Expr:
    if not expr.terms:
        return expr
    def get_sig(t):
        return f"{t.c_ops},{t.a_ops}"
    zero_term = Term([], [], 0)
    zero_term_sig = get_sig(zero_term)
    signatures = [ get_sig(t) for t in expr.terms ]
    terms = []
    term = expr.terms[0]
    term_sig = signatures[0]
    for t, t_sig in zip(expr.terms[1:], signatures[1:]):
        if ea.is_zero(term.coef):
            term = t
            term_sig = t_sig
        else:
            ct = combine_two_terms(t, term, t_sig, term_sig)
            if ct is None:
                terms.append(term)
                term = t
                term_sig = t_sig
            elif ea.is_zero(ct.coef):
                term = zero_term
                term_sig = zero_term_sig
            else:
                term = ct
    if not ea.is_zero(term.coef):
        terms.append(term)
    return Expr(terms, expr.description)

@q.timer
def drop_zero_terms(expr:Expr) -> Expr:
    terms = []
    for t in expr.terms:
        if not ea.is_zero(t.coef):
            terms.append(t)
    return Expr(terms, expr.description)

def op_derivative_exp(op:Op):
    if op.otype == "Qv":
        return Term([], [SHv(op.f, op.p, op.s, op.c),], -1)
    elif op.otype == "Qb":
        return Term([], [HbS(op.f, op.p, op.s, op.c),], 1)
    else:
        return None

def op_derivative_op(op:Op, op1:Op):
    if op.otype == "Qv" and op1.otype == "HbS" and op.f == op1.f:
        return Term([S(op.f, op.p, op1.p, op.s, op1.s, op.c, op1.c),], [], 1)
    elif op.otype == "Qb" and op1.otype == "SHv" and op.f == op1.f:
        return Term([S(op.f, op1.p, op.p, op1.s, op.s, op1.c, op.c),], [], 1)
    else:
        return None

def flip_sign(i:int) -> int:
    if i % 2 == 0:
        return 1
    else:
        return -1

def op_derivative_term(op:Op, term:Term) -> Expr:
    coef = term.coef
    c_ops = term.c_ops
    a_ops = term.a_ops
    terms = []
    for i in range(len(a_ops)):
        op1 = a_ops[i]
        dop1 = op_derivative_op(op, op1)
        if dop1 is not None:
            sign = flip_sign(i)
            terms.append(Term(dop1.c_ops + c_ops, a_ops[:i] + dop1.a_ops + a_ops[i+1:], sign * dop1.coef * coef))
    de = op_derivative_exp(op)
    if de is not None:
        sign = flip_sign(len(a_ops))
        terms.append(Term(de.c_ops + c_ops, a_ops + de.a_ops, sign * de.coef * coef))
    return Expr(terms)

def op_push_term(op:Op, term:Term) -> Expr:
    if op.otype == "Qv" or op.otype == "Qb":
        return op_derivative_term(op, term)
    else:
        coef = term.coef
        c_ops = term.c_ops
        a_ops = term.a_ops
        if op.is_commute():
            return Expr([Term([op,] + c_ops, a_ops, coef),])
        else:
            return Expr([Term(c_ops, [op,] + a_ops, coef),])

def op_push_expr(op:Op, expr:Expr) -> Expr:
    terms = []
    for term in expr.terms:
        terms += op_push_term(op, term).terms
    return Expr(terms)

def is_hop(op:Op) -> bool:
    if op.otype == "SHv" or op.otype == "HbS":
        return True
    if op.otype == "Hv" or op.otype == "Hb":
        return True
    return False

def has_hops(term:Term, count_limit:int = 0) -> bool:
    c = 0
    for op in term.a_ops:
        if is_hop(op):
            c += 1
            if c > count_limit:
                return True
    return False

def remove_hops(expr:Expr, count_limit:int = 0) -> Expr:
    terms = []
    for term in expr.terms:
        if not has_hops(term, count_limit):
            terms.append(term)
    return Expr(terms)

def contract_term(term:Term) -> Expr:
    coef = term.coef
    c_ops = term.c_ops
    a_ops = term.a_ops
    n_a_ops = len(a_ops)
    expr = Expr([Term(c_ops, [], coef),])
    for idx, op in enumerate(reversed(a_ops)):
        expr = op_push_expr(op, expr)
        expr = remove_hops(expr, n_a_ops - idx - 1)
    return expr

@q.timer
def contract_expr(expr: Expr) -> Expr:
    """
    interface function
    does not change expr
    """
    all_terms = []
    for term in expr.terms:
        all_terms += contract_term(term).terms
    return Expr(all_terms, f"< {expr.show()} >")

### ------

def S_l(p1, p2):
    return mk_expr(S('l', p1, p2)) + f"S_l({p1},{p2})"

def S_s(p1, p2):
    return mk_expr(S('s', p1, p2)) + f"S_s({p1},{p2})"

def S_c(p1, p2):
    return mk_expr(S('c', p1, p2)) + f"S_c({p1},{p2})"

def tr(expr):
    if isinstance(expr, Term):
        term = expr
        assert term.a_ops == []
        return Term([ Tr(term.c_ops), ], [], term.coef)
    elif isinstance(expr, Expr):
        return sum(map(tr, expr.terms)) + f"tr( {expr.show()} )"
    else:
        assert False

gamma_x = mk_expr(G(0)) + f"gamma_x"

gamma_y = mk_expr(G(1)) + f"gamma_y"

gamma_z = mk_expr(G(2)) + f"gamma_z"

gamma_t = mk_expr(G(3)) + f"gamma_t"

gamma_5 = mk_expr(G(5)) + f"gamma_5"

def gamma(tag):
    return mk_expr(G(tag)) + f"gamma({tag})"

def gamma_va(tag):
    assert isinstance(tag, int)
    if tag in [ 0, 1, 2, 3, ]:
        return mk_expr(G(tag)) + f"gamma({tag})"
    elif tag in [ 4, 5, 6, 7, ]:
        return mk_expr(G(tag - 4)) * mk_expr(G(5)) + f"gamma({tag-4})*gamma_5"
    else:
        assert False

### ------

if __name__ == "__main__":
    expr = (1
            * Qb("d", "x1", "s1", "c1")
            * G(5, "s1", "s2")
            * Qv("u", "x1", "s2", "c1")
            * Qb("u", "x2", "s3", "c2")
            * G(5, "s3", "s4")
            * Qv("d", "x2", "s4", "c2"))
    print(expr)
    c_expr = contract_expr(expr)
    c_expr.simplify(is_isospin_symmetric_limit = False)
    print(c_expr)
    c_expr_check = Expr([Term([Tr([G('5'), S('d','x2','x1'), G('5'), S('u','x1','x2')],'sc')],[],(-1+0j))]) - c_expr
    c_expr_check.simplify(is_isospin_symmetric_limit = False)
    print(c_expr_check)
    c_expr.simplify(is_isospin_symmetric_limit = True)
    print(c_expr)
    print(Qb("u", "x2", "s23", "c23") * Qv("u", "x1", "s11", "c11") - Qb("d", "x2", "s23", "c23") * Qv("d", "x1", "s11", "c11"))
    print(bfield_tag_dict['std-u'].get_spin_tensor_code())
