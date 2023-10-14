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
import io
import tokenize
import sympy

def get_variables(s):
    """
    `s` is a string
    return a list of variables contained in `s`
    """
    variables = []
    g = tokenize.tokenize(io.BytesIO(s.encode('utf-8')).readline)  # tokenize the string
    for toknum, tokval, _, _, _ in g:
        if toknum == tokenize.NAME:
            variables.append(tokval)
    variables = sorted(list(set(variables)))
    return variables

class Factor:

    """
    self.code
    self.variables
    self.otype
    #
    self.otype in [ "Expr", "Var", ]
    """

    def __init__(self, code, variables=None, otype=None):
        self.code = code
        self.otype = otype
        self.variables = variables
        if self.variables is None:
            self.variables = get_variables(self.code)
        if self.otype is None:
            if code.isidentifier():
                self.otype = "Var"
            else:
                self.otype = "Expr"

    def __repr__(self) -> str:
        return f"ea.Factor({self.code},{self.variables},{self.otype})"

    def __str__(self) -> str:
        return self.compile_py()

    def compile_py(self, var_dict=None) -> str:
        if (var_dict is None) or (self.otype != "Var") or (self.code not in var_dict):
            return f"{self.code}"
        return var_dict[self.code]

### ------

class Term:

    """
    self.coef
    self.factors
    #
    coef should be a number or a sympy expression
    """

    def __init__(self, factors, coef = 1):
        self.coef = coef
        self.factors = factors

    def __repr__(self) -> str:
        return f"ea.Term({self.factors},{self.coef})"

    def __str__(self) -> str:
        return self.compile_py()

    def sort(self) -> None:
        self.factors.sort(key = repr)

    def simplify_coef(self) -> None:
        self.coef = sympy.simplify(self.coef)

    def compile_py(self, var_dict=None) -> str:
        fs = [ f"({f.compile_py(var_dict)})" for f in self.factors ]
        if self.coef == 1:
            if self.factors == []:
                return "1"
            else:
                return '*'.join(fs)
        return '*'.join([ compile_py_complex(self.coef), ] + fs)

### ------

class Expr:

    """
    self.terms
    """

    def __init__(self, terms):
        self.terms = terms

    def __repr__(self) -> str:
        return f"ea.Expr({self.terms})"

    def __str__(self) -> str:
        return self.compile_py()

    def __add__(self, other):
        return Expr(self.terms + mk_expr(other).terms)

    __radd__ = __add__

    def __mul__(self, other):
        other = mk_expr(other)
        terms = []
        for t1 in self.terms:
            for t2 in other.terms:
                coef = t1.coef * t2.coef
                t = Term(t1.factors + t2.factors, coef)
                terms.append(t)
        return Expr(terms)

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
        for term in self.terms:
            term.sort()
        self.terms.sort(key = repr)

    def combine_terms(self) -> None:
        self.terms = combine_terms_expr(self).terms

    def drop_zeros(self) -> None:
        self.terms = drop_zero_terms(self).terms

    @q.timer
    def simplify_coef(self) -> None:
        for t in self.terms:
            t.simplify_coef()
        self.drop_zeros()

    @q.timer
    def simplify(self) -> None:
        """
        interface function
        """
        self.sort()
        self.combine_terms()
        self.drop_zeros()

    def is_zero(self) -> bool:
        return not self.terms

    def compile_py(self, var_dict=None) -> str:
        if self.terms:
            return '+'.join([ f"{t.compile_py(var_dict)}" for t in self.terms ])
        else:
            return '0'

### ------

def simplified(x):
    """
    interface function
    """
    x = copy.deepcopy(mk_expr(x))
    x.simplify()
    if not x.terms:
        x = 0
    elif len(x.terms) == 1:
        if not x.terms[0].factors:
            x = x.terms[0].coef
    return x

def coef_simplified(x):
    """
    interface function
    """
    x = copy.deepcopy(mk_expr(x))
    x.simplify_coef()
    return x

def compile_py(x, var_dict=None):
    """
    interface function
    """
    if isinstance(x, (int, float, complex, sympy.Basic,)):
        return compile_py_complex(x)
    else:
        return mk_expr(x).compile_py(var_dict)

def is_zero(x):
    """
    interface function
    """
    if isinstance(x, (int, float, complex, sympy.Basic,)):
        return x == 0
    elif isinstance(x, Expr):
        return x.is_zero()
    else:
        print(x)
        assert False

def mk_expr(x):
    """
    interface function
    """
    if isinstance(x, Expr):
        return x
    elif isinstance(x, Term):
        return Expr([ x, ])
    elif isinstance(x, Factor):
        return Expr([ Term([ x, ]), ])
    elif isinstance(x, (int, float, complex, sympy.Basic,)):
        return Expr([ Term([], x,), ])
    elif isinstance(x, str):
        # str viewed as code segment
        return Expr([ Term([ Factor(x), ]), ])
    else:
        print(x)
        assert False

def compile_py_complex(x):
    """
    interface function
    """
    v = complex(x)
    return f"{v}"

def drop_zero_terms(expr : Expr) -> Expr:
    terms = []
    for t in expr.terms:
        if t.coef != 0:
            terms.append(t)
    return Expr(terms)

def combine_two_terms(t1 : Term, t2 : Term, t1_sig : str, t2_sig : str):
    if t1_sig == t2_sig:
        coef = t1.coef + t2.coef
        if coef == 0:
            return Term([], 0)
        else:
            return Term(t1.factors, coef)
    else:
        return None

def combine_terms_expr(expr : Expr) -> Expr:
    if not expr.terms:
        return expr
    def get_sig(t):
        return f"{t.factors}"
    zero_term = Term([], 0)
    zero_term_sig = get_sig(zero_term)
    signatures = [ get_sig(t) for t in expr.terms ]
    terms = []
    term = expr.terms[0]
    term_sig = signatures[0]
    for t, t_sig in zip(expr.terms[1:], signatures[1:]):
        if term.coef == 0:
            term = t
            term_sig = t_sig
        else:
            ct = combine_two_terms(t, term, t_sig, term_sig)
            if ct is None:
                terms.append(term)
                term = t
                term_sig = t_sig
            elif ct.coef == 0:
                term = zero_term
                term_sig = zero_term_sig
            else:
                term = ct
    if term.coef != 0:
        terms.append(term)
    return Expr(terms)

### ------

def mk_sym(x):
    """
    make a sympy simplified value
    """
    return sympy.simplify(x)

def mk_fac(x):
    """
    make an Expr obj (can be sympy expression)
    """
    return mk_expr(x)

### ------

if __name__ == "__main__":
    a = mk_expr(1)
    b = mk_expr(2)
    c = mk_expr(Factor("a + b"))
    d = mk_expr(Factor("a + b"))
    c = c * b + c + d * d
    print(c)
    c.simplify()
    print(c)
    print(c.compile_py())
