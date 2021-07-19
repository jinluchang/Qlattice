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

from auto_contractor.wick import *

class Var(Op):

    def __init__(self, name : str):
        Op.__init__(self, "Var")
        self.name = name

    def __repr__(self):
        return f"{self.otype}({self.name!r})"

    def list(self):
        return [self.otype, self.name]

    def __eq__(self, other) -> bool:
        return self.list() == other.list()

def add_positions(s, x):
    if isinstance(x, Term):
        for op in x.c_ops:
            add_positions(s, op)
        for op in x.a_ops:
            add_positions(s, op)
    elif isinstance(x, Op):
        if x.otype == "S":
            s.update([x.p1, x.p2])
        elif x.otype == "Tr":
            for op in x.ops:
                add_positions(s, op)
        elif x.otype == "Qfield":
            s.add(x.p)
    elif isinstance(x, Expr):
        for t in x.terms:
            add_positions(s, t)

def add_prop_variables(variables, x):
    prop_counter = 0
    prop_dataset = {} # prop_dataset[prop_repr] = prop_var
    if isinstance(x, list):
        for i, op in enumerate(x):
            if op.otype == "S":
                prop_repr = repr(op)
                if prop_repr in prop_dataset:
                    x[i] = prop_dataset[prop_repr]
                else:
                    while True:
                        prop_counter += 1
                        name = f"S_{prop_counter}"
                        if name not in variables:
                            break
                    variables[name] = op
                    var = Var(name)
                    prop_dataset[prop_repr] = var
                    x[i] = var
            elif op.otype == "Tr":
                add_prop_variables(variables, op.ops)
    elif isinstance(x, Term):
        add_prop_variables(variables, x.c_ops)
    elif isinstance(x, Expr):
        for t in x.terms:
            add_prop_variables(variables, t)

class CExpr:

    def __init__(self, variables, named_terms, positions = None):
        self.variables = variables
        self.named_terms = named_terms
        if positions is not None:
            self.positions = positions
        else:
            s = set()
            for name, term in named_terms:
                add_positions(s, term)
            self.positions = sorted(list(s))

    def __repr__(self) -> str:
        return f"CExpr({self.variables},{self.named_terms},{self.positions})"

    def collect_prop(self):
        # interface function
        variables = dict(self.variables)
        for name, term in self.named_terms:
            add_prop_variables(variables, term)
        self.variables = sorted(variables.items())

def mk_cexpr(expr : Expr):
    # interface function
    return CExpr([], [ (f"T_{i+1}", term,) for i, term in enumerate(expr.terms) ])

def contract_simplify_round_compile(expr : Expr, is_isospin_symmetric_limit = True):
    # interface function
    expr = copy.deepcopy(expr)
    expr = contract_expr(expr)
    expr.simplify(is_isospin_symmetric_limit = is_isospin_symmetric_limit)
    cexpr = mk_cexpr(expr.round())
    return cexpr

def contract_simplify_round_compile_collect(expr : Expr, is_isospin_symmetric_limit = True):
    # interface function
    expr = copy.deepcopy(expr)
    expr = contract_expr(expr)
    expr.simplify(is_isospin_symmetric_limit = is_isospin_symmetric_limit)
    cexpr = mk_cexpr(expr.round())
    cexpr.collect_prop()
    return cexpr

def display_cexpr(cexpr : CExpr):
    # interface function
    # return a string
    lines = []
    lines.append(f"Begin CExpr")
    lines.append(f"{'Positions':>10} : {cexpr.positions}")
    for name, value in cexpr.variables:
        lines.append(f"{name:>10} : {value}")
    for name, term in cexpr.named_terms:
        lines.append(f"{name:>10} : {term}")
    lines.append(f"End CExpr")
    return "\n".join(lines)

if __name__ == "__main__":
    expr = Qb("d", "x1", "s1", "c1") * G(5, "s1", "s2") * Qv("u", "x1", "s2", "c1") * Qb("u", "x2", "s3", "c2") * G(5, "s3", "s4") * Qv("d", "x2", "s4", "c2")
    print(expr)
    expr = simplified(contract_expr(expr))
    print(expr.round())
    cexpr = mk_cexpr(expr.round())
    print(cexpr)
    cexpr.collect_prop()
    print(cexpr)
    print(display_cexpr(cexpr))
    print(CExpr([('S_1', S('d','x2','x1')), ('S_2', S('u','x1','x2'))],[('T_1', Term([Tr([G(5), Var('S_1'), G(5), Var('S_2')],'sc')],[],(-1+0j)))],['x1', 'x2']))
