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

try:
    from .wick import *
    from . import auto_fac_funcs as aff
except:
    from wick import *
    import auto_fac_funcs as aff

from itertools import permutations

import qlat as q

class Var(Op):

    def __init__(self, name:str):
        Op.__init__(self, "Var")
        self.name = name

    def __repr__(self):
        return f"{self.otype}({self.name!r})"

    def list(self):
        return [self.otype, self.name]

    def __eq__(self, other) -> bool:
        return self.list() == other.list()

### ----

def get_var_name_type(x):
    """
    types include: V_S (wilson matrix), V_G (spin matrix), V_U (color matrix), V_a (c-number)
    """
    if x.startswith("V_S_"):
        return "V_S"
    elif x.startswith("V_U_"):
        return "V_U"
    elif x.startswith("V_tr_"):
        return "V_a"
    elif x.startswith("V_chain_"):
        return "V_S"
    elif x.startswith("V_bs_"):
        return "V_a"
    elif x.startswith("V_prod_GG_"):
        return "V_G"
    elif x.startswith("V_prod_GS_"):
        return "V_S"
    elif x.startswith("V_prod_SG_"):
        return "V_S"
    elif x.startswith("V_prod_SS_"):
        return "V_S"
    elif x.startswith("V_prod_UU_"):
        return "V_U"
    elif x.startswith("V_prod_UG_"):
        return "V_S"
    elif x.startswith("V_prod_GU_"):
        return "V_S"
    elif x.startswith("V_prod_US_"):
        return "V_S"
    elif x.startswith("V_prod_SU_"):
        return "V_S"
    else:
        assert False

def get_op_type(x):
    """
    ``x`` should be an Op
    return a string which indicate the type of the op
    """
    if not isinstance(x, Op):
        return None
    if x.otype == "S":
        return "S"
    elif x.otype == "G":
        return "G"
    elif x.otype == "U":
        return "U"
    elif x.otype == "Tr":
        return "Tr"
    elif x.otype == "Chain":
        return "Chain"
    elif x.otype == "BS":
        return "BS"
    elif x.otype == "Var":
        return get_var_name_type(x.name)
    else:
        assert False
    return None

def add_positions(s, x):
    """
    ``s`` is set of str
    ``x`` is expr/sub-expr
    """
    if isinstance(x, Term):
        add_positions(s, x.coef)
        for op in x.c_ops:
            add_positions(s, op)
        for op in x.a_ops:
            add_positions(s, op)
    elif isinstance(x, Op):
        if x.otype == "S":
            if isinstance(x.p1, str):
                s.add(x.p1)
            if isinstance(x.p2, str):
                s.add(x.p2)
        elif x.otype == "G":
            if isinstance(x.tag, str):
                s.add(x.tag)
        elif x.otype == "U":
            if isinstance(x.p, str):
                s.add(x.p)
            if isinstance(x.mu, str):
                s.add(x.mu)
        elif x.otype == "Tr":
            for op in x.ops:
                add_positions(s, op)
        elif x.otype == "Chain":
            for op in x.ops:
                add_positions(s, op)
        elif x.otype == "BS":
            for op in x.chain_list:
                add_positions(s, op)
            for v, c, in x.elem_list:
                add_positions(s, c)
        elif x.otype == "Qfield":
            if isinstance(x.p, str):
                s.add(x.p)
    elif isinstance(x, Expr):
        for t in x.terms:
            add_positions(s, t)
    elif isinstance(x, ea.Expr):
        for t in x.terms:
            add_positions(s, t)
    elif isinstance(x, ea.Term):
        for f in x.factors:
            add_positions(s, f)
    elif isinstance(x, ea.Factor):
        for v in x.variables:
            s.add(v)

def get_positions(term):
    s = set()
    add_positions(s, term)
    return sorted(list(s))

@q.timer
def collect_position_in_cexpr(named_terms, named_exprs):
    s = set()
    for name, term in named_terms:
        add_positions(s, term)
    for name, expr_list in named_exprs:
        for coef, term_name in expr_list:
            add_positions(s, coef)
    positions = sorted(list(s))
    return positions

@q.timer
def find_common_prod_in_factors(variables_factor):
    """
    return [ (var_name_1, var_name_2,), ... ]
    """
    subexpr_count = {}
    for _, ea_coef in variables_factor:
        assert isinstance(ea_coef, ea.Expr)
        for t in ea_coef.terms:
            x = t.factors
            if len(x) <= 1:
                continue
            for i, f in enumerate(x[:-1]):
                assert f.otype == "Var"
                f1 = x[i+1]
                assert f1.otype == "Var"
                prod = (f.code, f1.code,)
                count = subexpr_count.get(prod, 0)
                subexpr_count[prod] = count + 1
    max_num_repeat = 0
    best_match_list = []
    for prod, num_repeat in subexpr_count.items():
        if num_repeat > max_num_repeat:
            max_num_repeat = num_repeat
            best_match_list = [ prod, ]
        elif num_repeat == max_num_repeat:
            best_match_list.append(prod)
    return best_match_list

@q.timer
def collect_common_prod_in_factors(variables_factor_intermediate, variables_factor, var_nameset, var_counter, common_prod_list, var_dataset):
    """
    common_prod_list = find_common_prod_in_factors(variables_factor)
    var = var_dataset[prod]
    var = ea.Factor(name, variables=[], otype="Var")
    """
    common_prod_set = set(common_prod_list)
    for _, ea_coef in variables_factor:
        assert isinstance(ea_coef, ea.Expr)
        for t in ea_coef.terms:
            x = t.factors
            if len(x) <= 1:
                continue
            for i in range(len(x) - 1):
                f = x[i]
                if f is None:
                    continue
                assert f.otype == "Var"
                f1 = x[i+1]
                assert f1.otype == "Var"
                prod = (f.code, f1.code,)
                if prod in common_prod_set:
                    if prod in var_dataset:
                        var = var_dataset[prod]
                    else:
                        code1, code2 = prod
                        var1 = ea.Factor(code1, variables=[], otype="Var")
                        var2 = ea.Factor(code2, variables=[], otype="Var")
                        prod_expr = ea.mk_expr(var1) * ea.mk_expr(var2)
                        while True:
                            name = f"V_factor_prod_{var_counter}"
                            var_counter += 1
                            if name not in var_nameset:
                                break
                        var_nameset.add(name)
                        variables_factor_intermediate.append((name, prod_expr,))
                        var = ea.Factor(name, variables=[], otype="Var")
                        var_dataset[prod] = var
                    x[i] = var
                    x[i+1] = None
    for _, ea_coef in variables_factor:
        assert isinstance(ea_coef, ea.Expr)
        for t in ea_coef.terms:
            x = t.factors
            x_new = []
            for v in x:
                if v is not None:
                    x_new.append(v)
            t.factors = x_new
    return var_counter

@q.timer
def find_common_sum_in_factors(variables_factor):
    """
    return [ (var_name_1, var_name_2,), ... ]
    """
    subexpr_count = {}
    for _, ea_coef in variables_factor:
        assert isinstance(ea_coef, ea.Expr)
        x = ea_coef.terms
        if len(x) < 1:
            continue
        for i, t in enumerate(x[:-1]):
            t1 = x[i+1]
            assert t.coef == 1
            assert len(t.factors) == 1
            assert t1.coef == 1
            assert len(t1.factors) == 1
            f = t.factors[0]
            f1 = t1.factors[0]
            assert f.otype == "Var"
            assert f1.otype == "Var"
            pair = (f.code, f1.code,)
            count = subexpr_count.get(pair, 0)
            subexpr_count[pair] = count + 1
    max_num_repeat = 0
    best_match_list = []
    for pair, num_repeat in subexpr_count.items():
        if num_repeat > max_num_repeat:
            max_num_repeat = num_repeat
            best_match_list = [ pair, ]
        elif num_repeat == max_num_repeat:
            best_match_list.append(pair)
    return best_match_list

@q.timer
def collect_common_sum_in_factors(variables_factor_intermediate, variables_factor, var_nameset, var_counter, common_pair_list, var_dataset):
    """
    common_pair_list = find_common_sum_in_factors(variables_factor)
    var = var_dataset[pair]
    var = ea.Factor(name, variables=[], otype="Var")
    """
    common_pair_set = set(common_pair_list)
    for _, ea_coef in variables_factor:
        assert isinstance(ea_coef, ea.Expr)
        x = ea_coef.terms
        if len(x) <= 1:
            continue
        for i in range(len(x) - 1):
            t = x[i]
            if t is None:
                continue
            t1 = x[i+1]
            assert t.coef == 1
            assert len(t.factors) == 1
            assert t1.coef == 1
            assert len(t1.factors) == 1
            f = t.factors[0]
            f1 = t1.factors[0]
            assert f.otype == "Var"
            assert f1.otype == "Var"
            pair = (f.code, f1.code,)
            if pair in common_pair_set:
                if pair in var_dataset:
                    var = var_dataset[pair]
                else:
                    code1, code2 = pair
                    var1 = ea.Factor(code1, variables=[], otype="Var")
                    var2 = ea.Factor(code2, variables=[], otype="Var")
                    pair_expr = ea.mk_expr(var1) + ea.mk_expr(var2)
                    while True:
                        name = f"V_factor_sum_{var_counter}"
                        var_counter += 1
                        if name not in var_nameset:
                            break
                    var_nameset.add(name)
                    variables_factor_intermediate.append((name, pair_expr,))
                    var = ea.Factor(name, variables=[], otype="Var")
                    var_dataset[pair] = var
                x[i].factors[0] = var
                x[i+1] = None
    for _, ea_coef in variables_factor:
        assert isinstance(ea_coef, ea.Expr)
        x = ea_coef.terms
        x_new = []
        for v in x:
            if v is not None:
                x_new.append(v)
        ea_coef.terms = x_new
    return var_counter

@q.timer
def collect_factor_in_cexpr(variables_factor, var_nameset, named_exprs, named_terms):
    """
    Make variables to all coefs of all the terms
    """
    var_counter = 0
    var_dataset = {}
    def add_variables(ea_coef):
        nonlocal var_counter
        key = repr(ea_coef)
        if key in var_dataset:
            var = var_dataset[key]
        else:
            s_ea_coef = ea.mk_expr(ea.simplified_ea(ea_coef))
            key2 = repr(s_ea_coef)
            if key2 in var_dataset:
                var = var_dataset[key2]
                var_dataset[key] = var
            else:
                while True:
                    name = f"V_factor_final_{var_counter}"
                    var_counter += 1
                    if name not in var_nameset:
                        break
                var_nameset.add(name)
                variables_factor.append((name, s_ea_coef,))
                var = ea.Factor(name, variables=[], otype="Var")
                var_dataset[key] = var
                var_dataset[key2] = var
        ea_var = ea.mk_expr(var)
        return ea_var
    for _, expr in named_exprs:
        for i, (ea_coef, term_name,) in enumerate(expr):
            ea_var = add_variables(ea_coef)
            expr[i] = (ea_var, term_name,)
    for _, term in named_terms:
        for op in term.c_ops:
            if op.otype == "BS":
                for i, (val, ea_coef,) in enumerate(op.elem_list):
                    ea_var = add_variables(ea_coef)
                    op.elem_list[i] = (val, ea_var)

@q.timer
def collect_factor_coef_in_cexpr(variables_factor_intermediate, var_nameset, variables_factor):
    """
    Add numerical coef to variables_factor_intermediate
    """
    var_counter = 0
    var_dataset = {} # var_dataset[factor_code] = factor_var
    for _, ea_coef in variables_factor:
        assert isinstance(ea_coef, ea.Expr)
        for t in ea_coef.terms:
            x = t.coef
            if x == 1 or x == -1:
                continue
            code = ea.compile_py_complex(x)
            if code in var_dataset:
                t.coef = 1
                t.factors.append(var_dataset[code])
            else:
                while True:
                    name = f"V_factor_coef_{var_counter}"
                    var_counter += 1
                    if name not in var_nameset:
                        break
                var_nameset.add(name)
                variables_factor_intermediate.append((name, ea.mk_expr(t.coef),))
                t.coef = 1
                var = ea.Factor(name, variables=[], otype="Var")
                t.factors.append(var)
                var_dataset[code] = var

@q.timer
def collect_factor_fac_in_cexpr(variables_factor_intermediate, var_nameset, variables_factor):
    """
    Add ea.Factor with `otype` "Expr" to `variables_factor_intermediate`
    """
    var_counter = 0
    var_dataset = {} # var_dataset[factor_code] = factor_var
    for _, ea_coef in variables_factor:
        assert isinstance(ea_coef, ea.Expr)
        for t in ea_coef.terms:
            x = t.factors
            for i, f in enumerate(x):
                if f.otype != "Var":
                    assert f.otype == "Expr"
                    if f.code in var_dataset:
                        x[i] = var_dataset[f.code]
                    else:
                        while True:
                            name = f"V_factor_fac_{var_counter}"
                            var_counter += 1
                            if name not in var_nameset:
                                break
                        var_nameset.add(name)
                        variables_factor_intermediate.append((name, ea.mk_expr(f),))
                        var = ea.Factor(name, f.variables)
                        x[i] = var
                        var_dataset[f.code] = var

@q.timer
def collect_factor_facs_in_cexpr(variables_factors, var_nameset, variables_factor):
    """
    Add ea.Factor with `otype` "Expr" to `variables_factors`
    """
    var_counter = 0
    var_dataset = {} # var_dataset[factor_code] = factor_var
    for _, ea_coef in variables_factor:
        assert isinstance(ea_coef, ea.Expr)
        for t in ea_coef.terms:
            x = t.factors
            if len(x) <= 1:
                continue
            key = repr(x)
            if key in var_dataset:
                var = var_dataset[key]
                t.factors = [ var, ]
            else:
                while True:
                    name = f"V_factor_facs_{var_counter}"
                    var_counter += 1
                    if name not in var_nameset:
                        break
                var_nameset.add(name)
                expr = ea.mk_expr(ea.Term(x))
                variables_factors.append((name, expr,))
                var = ea.Factor(name, variables=[], otype="Var")
                t.factors = [ var, ]
                var_dataset[key] = var

@q.timer
def collect_factor_prod_in_cexpr(variables_factor_intermediate, var_nameset, variables_factor):
    """
    Common product elimination
    """
    var_counter = 0
    var_dataset = {} # var_dataset[(code1, code2,)] = factor_var
    while True:
        prod_list = find_common_prod_in_factors(variables_factor)
        if len(prod_list) == 0:
            break
        var_counter = collect_common_prod_in_factors(variables_factor_intermediate, variables_factor, var_nameset, var_counter, prod_list, var_dataset)

@q.timer
def collect_factor_sum_in_cexpr(variables_factor_intermediate, var_nameset, variables_factor):
    """
    Common summation elimination
    """
    var_counter = 0
    var_dataset = {} # var_dataset[(code1, code2,)] = factor_var
    while True:
        pair_list = find_common_sum_in_factors(variables_factor)
        if len(pair_list) == 0:
            break
        var_counter = collect_common_sum_in_factors(variables_factor_intermediate, variables_factor, var_nameset, var_counter, pair_list, var_dataset)

@q.timer
def collect_and_optimize_factor_in_cexpr(named_exprs, named_terms):
    """
    collect the factors in all ea_coef
    collect common sub-expressions in ea_coef
    #
    variables_factor = [ (name, ea.expr,), ... ]
    variables_factor_intermediate = [ (name, ea.expr,), ... ]
    """
    variables_factor_intermediate = []
    variables_factor = []
    var_nameset = set()
    collect_factor_in_cexpr(variables_factor, var_nameset, named_exprs, named_terms)
    collect_factor_coef_in_cexpr(variables_factor_intermediate, var_nameset, variables_factor)
    collect_factor_fac_in_cexpr(variables_factor_intermediate, var_nameset, variables_factor)
    variables_factors = []
    collect_factor_facs_in_cexpr(variables_factors, var_nameset, variables_factor)
    collect_factor_prod_in_cexpr(variables_factor_intermediate, var_nameset, variables_factors)
    variables_factor_intermediate += variables_factors
    collect_factor_sum_in_cexpr(variables_factor_intermediate, var_nameset, variables_factor)
    return variables_factor_intermediate, variables_factor

@q.timer
def collect_prop_in_cexpr(named_terms):
    """
    collect the propagators
    modify the named_terms in-place and return the prop variable definitions as variables
    """
    variables_prop = []
    var_counter = 0
    var_dataset = {} # var_dataset[op_repr] = op_var
    var_nameset = set()
    def add_prop_variables(x):
        nonlocal var_counter
        if isinstance(x, list):
            for i, op in enumerate(x):
                if op.otype in [ "S", ]:
                    op_repr = repr(op)
                    if op_repr in var_dataset:
                        x[i] = var_dataset[op_repr]
                    else:
                        while True:
                            name = f"V_S_{var_counter}"
                            var_counter += 1
                            if name not in var_nameset:
                                break
                        var_nameset.add(name)
                        variables_prop.append((name, op,))
                        var = Var(name)
                        x[i] = var
                        var_dataset[op_repr] = var
                elif op.otype == "Tr":
                    add_prop_variables(op.ops)
                elif op.otype == "Chain":
                    add_prop_variables(op.ops)
                elif op.otype == "BS":
                    add_prop_variables(op.chain_list)
        elif isinstance(x, Term):
            add_prop_variables(x.c_ops)
        elif isinstance(x, Expr):
            for t in x.terms:
                add_prop_variables(t)
    for name, term in named_terms:
        add_prop_variables(term)
        term.sort()
    return variables_prop

@q.timer
def collect_color_matrix_in_cexpr(named_terms):
    """
    collect the color matrices
    modify the named_terms in-place and return the color matrix variable definitions as variables
    """
    variables_color_matrix = []
    var_counter = 0
    var_dataset = {} # var_dataset[op_repr] = op_var
    var_nameset = set()
    def add_variables(x):
        nonlocal var_counter
        if isinstance(x, list):
            for i, op in enumerate(x):
                if op.otype in [ "U", ]:
                    op_repr = repr(op)
                    if op_repr in var_dataset:
                        x[i] = var_dataset[op_repr]
                    else:
                        while True:
                            name = f"V_U_{var_counter}"
                            var_counter += 1
                            if name not in var_nameset:
                                break
                        var_nameset.add(name)
                        variables_color_matrix.append((name, op,))
                        var = Var(name)
                        x[i] = var
                        var_dataset[op_repr] = var
                elif op.otype == "Tr":
                    add_variables(op.ops)
                elif op.otype == "Chain":
                    add_variables(op.ops)
                elif op.otype == "BS":
                    add_variables(op.chain_list)
        elif isinstance(x, Term):
            add_variables(x.c_ops)
        elif isinstance(x, Expr):
            for t in x.terms:
                add_variables(t)
    for name, term in named_terms:
        add_variables(term)
        term.sort()
    return variables_color_matrix

@q.timer
def collect_chain_in_cexpr(named_terms):
    """
    collect common chains
    modify named_terms in-place and return definitions as variables_chain
    variables_chain = [ (name, value,), ... ]
    possible (name, value,) includes
    ("V_chain_0", op,) where op.otype == "Chain"
    """
    var_nameset = set()
    variables_chain = []
    var_counter = 0
    var_dataset_chain = {} # var_dataset[op_repr] = op_var
    def add_ch_varibles(x):
        nonlocal var_counter
        if isinstance(x, Term):
            add_ch_varibles(x.c_ops)
        elif isinstance(x, Op) and x.otype == "Chain":
            add_ch_varibles(x.ops) # not very necessary, do not expect Chain in Chain
        elif isinstance(x, Op) and x.otype == "BS":
            add_ch_varibles(x.chain_list)
        elif isinstance(x, list):
            for op in x:
                add_ch_varibles(op)
            for i, op in enumerate(x):
                if isinstance(op, Op) and op.otype == "Chain":
                    op_repr = repr(op)
                    if op_repr in var_dataset_chain:
                        x[i] = var_dataset_chain[op_repr]
                    else:
                        while True:
                            name = f"V_chain_{var_counter}"
                            var_counter += 1
                            if name not in var_nameset:
                                break
                        var_nameset.add(name)
                        variables_chain.append((name, op,))
                        var = Var(name)
                        x[i] = var
                        var_dataset_chain[op_repr] = var
    for name, term in named_terms:
        add_ch_varibles(term)
        term.sort()
    return variables_chain

@q.timer
def collect_tr_in_cexpr(named_terms):
    """
    collect common traces
    modify named_terms in-place and return definitions as variables_tr
    variables_tr = [ (name, value,), ... ]
    possible (name, value,) includes
    ("V_tr_0", op,) where op.otype == "Tr"
    """
    var_nameset = set()
    variables_tr = []
    var_counter = 0
    var_dataset_tr = {} # var_dataset[op_repr] = op_var
    def add_tr_varibles(x):
        nonlocal var_counter
        if isinstance(x, Term):
            add_tr_varibles(x.c_ops)
        elif isinstance(x, Op) and x.otype == "Tr":
            add_tr_varibles(x.ops) # not very necessary, do not expect Tr in Tr
        elif isinstance(x, list):
            for op in x:
                add_tr_varibles(op)
            for i, op in enumerate(x):
                if isinstance(op, Op) and op.otype == "Tr":
                    op_repr = repr(op)
                    if op_repr in var_dataset_tr:
                        x[i] = var_dataset_tr[op_repr]
                    else:
                        while True:
                            name = f"V_tr_{var_counter}"
                            var_counter += 1
                            if name not in var_nameset:
                                break
                        var_nameset.add(name)
                        variables_tr.append((name, op,))
                        var = Var(name)
                        x[i] = var
                        var_dataset_tr[op_repr] = var
    for name, term in named_terms:
        add_tr_varibles(term)
        term.sort()
    return variables_tr

@q.timer
def collect_baryon_prop_in_cexpr(named_terms):
    """
    collect common baryon_prop
    modify named_terms in-place and return definitions as variables_baryon_prop
    variables_baryon_prop = [ (name_list, value_list,), ... ]
    len(name_list) == len(value_list)
    name_list = [ name, ... ]
    value_list = [ value, ... ]
    possible (name, value,) includes
    ("V_bs_0", op,) where op.otype == "BS"
    """
    var_nameset = set()
    var_counter = 0
    var_dataset = {} # var_dataset[op_repr] = op_var
    variable_dict = {} # variable_dict[name] = op
    var_chain_list_dataset = {} # var_dataset[op_chain_list_repr] = [ name, ... ]
    def add_varibles(x):
        nonlocal var_counter
        if isinstance(x, Term):
            add_varibles(x.c_ops)
        elif isinstance(x, list):
            for op in x:
                add_varibles(op)
            for i, op in enumerate(x):
                if isinstance(op, Op) and op.otype == "BS":
                    op_repr = repr(op)
                    if op_repr in var_dataset:
                        x[i] = var_dataset[op_repr]
                    else:
                        while True:
                            name = f"V_bs_{var_counter}"
                            var_counter += 1
                            if name not in var_nameset:
                                break
                        var_nameset.add(name)
                        var = Var(name)
                        x[i] = var
                        var_dataset[op_repr] = var
                        variable_dict[name] = op
                        op_chain_list_repr = repr(op.chain_list)
                        if op_chain_list_repr not in var_chain_list_dataset:
                            var_chain_list_dataset[op_chain_list_repr] = []
                        var_chain_list_dataset[op_chain_list_repr].append(name)
    for name, term in named_terms:
        add_varibles(term)
        term.sort()
    variables = []
    for name_list in var_chain_list_dataset.values():
        value_list = [ variable_dict[name] for name in name_list ]
        variables.append((name_list, value_list,))
    return variables

@q.timer
def find_common_subexpr_in_tr(variables_tr):
    """
    return None or (op, op1,)
    """
    subexpr_count = {}
    def add(x, count_added):
        op_repr = repr(x)
        if op_repr in subexpr_count:
            c, op = subexpr_count[op_repr]
            assert x == op
            subexpr_count[op_repr] = (c + count_added, x)
        else:
            subexpr_count[op_repr] = (count_added, x)
    def find_op_pair(x:list, is_periodic=True):
        if len(x) < 2:
            return None
        for i, op in enumerate(x):
            factor = 1
            if len(x) > 2:
                factor = 2
            op_type = get_op_type(op)
            if op_type in [ "V_S", "V_G", "V_U", "S", "G", "U", ]:
                i1 = i + 1
                if i1 >= len(x):
                    if is_periodic:
                        i1 = i1 % len(x)
                    else:
                        continue
                op1 = x[i1]
                op1_type = get_op_type(op1)
                if op1_type in [ "V_S", "V_G", "V_U", "S", "G", "U", ]:
                    prod = (op, op1,)
                    if op_type in [ "V_G", "G", ] and op1_type in [ "V_G", "G", ]:
                        add(prod, 1.02 * factor)
                    elif op_type in [ "V_U", "U", ] and op1_type in [ "V_U", "U", ]:
                        add(prod, 1.02 * factor)
                    elif op_type in [ "V_U", "U", ] and op1_type in [ "V_G", "G", ]:
                        # do not multiply U with G
                        pass
                    elif op_type in [ "V_G", "G", ] and op1_type in [ "V_U", "U", ]:
                        # do not multiply U with G
                        pass
                    elif op_type in [ "V_G", "G", ] or op1_type in [ "V_G", "G", ]:
                        add(prod, 1.01 * factor)
                    elif op_type in [ "V_U", "U", ] or op1_type in [ "V_U", "U", ]:
                        add(prod, 1.01 * factor)
                    elif op_type in [ "V_S", "S", ] and op1_type in [ "V_S", "S", ]:
                        add(prod, 1 * factor)
                    else:
                        assert False
    def find(x):
        if isinstance(x, list):
            for op in x:
                find(op)
        elif isinstance(x, Op) and x.otype == "Tr" and len(x.ops) >= 2:
            find_op_pair(x.ops, is_periodic=True)
        elif isinstance(x, Op) and x.otype == "Chain" and len(x.ops) >= 2:
            find_op_pair(x.ops, is_periodic=False)
        elif isinstance(x, Op) and x.otype == "BS":
            find(x.chain_list)
        elif isinstance(x, Term):
            find(x.c_ops)
        elif isinstance(x, Expr):
            for t in x.terms:
                find(t)
    for name, tr in variables_tr:
        find(tr)
    max_num_repeat = 1.5
    best_match = None
    for num_repeat, op in subexpr_count.values():
        if num_repeat > max_num_repeat:
            max_num_repeat = num_repeat
            best_match = op
    return best_match

@q.timer
def collect_common_subexpr_in_tr(variables_tr, op_common, var):
    op_repr = repr(op_common)
    def replace(x, is_periodic=True):
        if x is None:
            return None
        elif isinstance(x, list):
            # need to represent the product of the list of operators
            for op in x:
                replace(op, is_periodic=is_periodic)
            if len(x) < 2:
                return None
            for i, op in enumerate(x):
                if isinstance(op, Op) and op.otype in [ "Var", "S", "G", "U", ]:
                    i1 = i + 1
                    if i1 >= len(x):
                        if is_periodic:
                            i1 = i1 % len(x)
                        else:
                            continue
                    op1 = x[i1]
                    if isinstance(op1, Op) and op1.otype in [ "Var", "S", "G", "U", ]:
                        prod = (op, op1,)
                        if repr(prod) == op_repr:
                            x[i1] = None
                            x[i] = var
        elif isinstance(x, Op) and x.otype == "Tr" and len(x.ops) >= 2:
            replace(x.ops, is_periodic=True)
        elif isinstance(x, Op) and x.otype == "Chain" and len(x.ops) >= 2:
            replace(x.ops, is_periodic=False)
        elif isinstance(x, Op) and x.otype == "BS":
            replace(x.chain_list, is_periodic=True)
        elif isinstance(x, Term):
            replace(x.c_ops, is_periodic=True)
        elif isinstance(x, Expr):
            for t in x.terms:
                replace(t)
    def remove_none(x):
        """
        return a None removed x
        possibly modify in-place
        """
        if isinstance(x, list):
            return [ remove_none(op) for op in x if op is not None ]
        elif isinstance(x, Op):
            if x.otype == "Tr":
                x.ops = remove_none(x.ops)
            elif x.otype == "Chain":
                x.ops = remove_none(x.ops)
            elif x.otype == "BS":
                x.chain_list = remove_none(x.chain_list)
            return x
        elif isinstance(x, Term):
            x.c_ops = remove_none(x.c_ops)
            x.a_ops = remove_none(x.a_ops)
            return x
        elif isinstance(x, Expr):
            x.terms = [ remove_none(t) for t in x.terms ]
            return x
        else:
            assert False
    for name, tr in variables_tr:
        replace(tr)
        remove_none(tr)

@q.timer
def collect_subexpr_in_cexpr(variables_tr):
    """
    collect common sub-expressions
    modify variables_tr in-place and return definitions as variables_prod
    variables_prod = [ (name, value,), ... ]
    possible (name, value,) includes
    ("V_prod_SG_0", [ op, op1, ],) where get_op_type(op) in [ "V_S", "S", ] and get_op_type(op1) in [ "V_G", "G", ]
    """
    var_nameset = set()
    var_counter_dict = {}
    var_counter_dict["V_prod_GG_"] = 0
    var_counter_dict["V_prod_GS_"] = 0
    var_counter_dict["V_prod_SG_"] = 0
    var_counter_dict["V_prod_SS_"] = 0
    var_counter_dict["V_prod_UU_"] = 0
    var_counter_dict["V_prod_UG_"] = 0
    var_counter_dict["V_prod_GU_"] = 0
    var_counter_dict["V_prod_US_"] = 0
    var_counter_dict["V_prod_SU_"] = 0
    variables_prod = []
    while True:
        subexpr = find_common_subexpr_in_tr(variables_tr)
        if subexpr is None:
            break
        op, op1 = subexpr
        op_type = get_op_type(op)
        op1_type = get_op_type(op1)
        assert op_type in [ "V_S", "V_G", "V_U", "S", "G", "U", ]
        assert op1_type in [ "V_S", "V_G", "V_U", "S", "G", "U", ]
        if op_type in [ "V_G", "G", ] and op1_type in [ "V_G", "G", ]:
            name_prefix = "V_prod_GG_"
        elif op_type in [ "V_G", "G", ] and op1_type in [ "V_S", "S", ]:
            name_prefix = "V_prod_GS_"
        elif op_type in [ "V_S", "S", ] and op1_type in [ "V_G", "G", ]:
            name_prefix = "V_prod_SG_"
        elif op_type in [ "V_S", "S", ] and op1_type in [ "V_S", "S", ]:
            name_prefix = "V_prod_SS_"
        elif op_type in [ "V_U", "U", ] and op1_type in [ "V_U", "U", ]:
            name_prefix = "V_prod_UU_"
        elif op_type in [ "V_U", "U", ] and op1_type in [ "V_G", "G", ]:
            name_prefix = "V_prod_UG_"
        elif op_type in [ "V_G", "G", ] and op1_type in [ "V_U", "U", ]:
            name_prefix = "V_prod_GU_"
        elif op_type in [ "V_U", "U", ] and op1_type in [ "V_S", "S", ]:
            name_prefix = "V_prod_US_"
        elif op_type in [ "V_S", "S", ] and op1_type in [ "V_U", "U", ]:
            name_prefix = "V_prod_SU_"
        else:
            assert False
        while True:
            name = f"{name_prefix}{var_counter_dict[name_prefix]}"
            var_counter_dict[name_prefix] += 1
            if name not in var_nameset:
                break
        var_nameset.add(name)
        variables_prod.append((name, subexpr,))
        var = Var(name)
        collect_common_subexpr_in_tr(variables_tr, subexpr, var)
    return variables_prod

class CExpr:

    """
    self.diagram_types
    self.positions
    self.variables_factor_intermediate
    self.variables_factor
    self.variables_prop
    self.variables_color_matrix
    self.variables_prod
    self.variables_chain
    self.variables_tr
    self.variables_baryon_prop
    self.named_terms
    self.named_exprs
    #
    self.named_terms[i] = (term_name, Term(c_ops, [], 1),)
    self.named_exprs[i] = (expr_name, [ (ea_coef, term_name,), ... ],)
    self.positions == sorted(list(self.positions))
    """

    def __init__(self):
        self.diagram_types = []
        self.positions = []
        self.variables_factor_intermediate = []
        self.variables_factor = []
        self.variables_prop = []
        self.variables_color_matrix = []
        self.variables_prod = []
        self.variables_chain = []
        self.variables_tr = []
        self.variables_baryon_prop = []
        self.named_terms = []
        self.named_exprs = []

    @q.timer
    def copy(self):
        """
        return a deep copy of this object.
        """
        return copy.deepcopy(self)

    @q.timer
    def optimize(self):
        """
        interface function
        #
        This is necessary step to run the `CExpr`.
        """
        self.collect_op()

    @q.timer
    def collect_op(self):
        """
        interface function
        Performing common sub-expression elimination
        Should be called after contract_simplify_compile(*exprs) or mk_cexpr(*exprs)
        The cexpr cannot be evaluated before collect_op!!!
        eval term factor
        """
        for name, term in self.named_terms:
            assert term.coef == 1
            assert term.a_ops == []
        checks = [
                self.variables_factor_intermediate == [],
                self.variables_factor == [],
                self.variables_prop == [],
                self.variables_color_matrix == [],
                self.variables_prod == [],
                self.variables_chain == [],
                self.variables_tr == [],
                ]
        if not all(checks):
            # likely collect_op is already performed
            return
        # collect ea_coef factors into variables
        self.variables_factor_intermediate, self.variables_factor = collect_and_optimize_factor_in_cexpr(self.named_exprs, self.named_terms)
        # collect prop expr into variables
        self.variables_prop = collect_prop_in_cexpr(self.named_terms)
        # collect color matrix expr into variables
        self.variables_color_matrix = collect_color_matrix_in_cexpr(self.named_terms)
        # collect chain expr into variables
        self.variables_chain = collect_chain_in_cexpr(self.named_terms)
        # collect trace expr into variables
        self.variables_tr = collect_tr_in_cexpr(self.named_terms)
        # collect baryon_prop expr into variables
        self.variables_baryon_prop = collect_baryon_prop_in_cexpr(self.named_terms)
        # collect common prod into variables
        self.variables_prod = collect_subexpr_in_cexpr(self.variables_tr)

    def get_expr_names(self):
        return [ name for name, expr in self.named_exprs ]

    def list(self):
        return [
                self.diagram_types,
                self.positions,
                self.variables_factor_intermediate,
                self.variables_factor,
                self.variables_prop,
                self.variables_color_matrix,
                self.variables_prod,
                self.variables_chain,
                self.variables_tr,
                self.variables_baryon_prop,
                self.named_terms,
                self.named_exprs,
                ]

### ----

def increase_type_dict_count(type_dict, key):
    if key in type_dict:
        type_dict[key] += 1
    else:
        type_dict[key] = 1

def drop_tag_last_subscript(tag):
    """
    Coordinates with names that are same after dropping last subscript belong to the same permutation group.
    """
    return tag.rsplit("_", 1)[0]

def mk_permuting_dicts(pos_list):
    for p_pos_list in permutations(pos_list):
        yield dict(zip(pos_list, p_pos_list))

def get_position_permutation_groups(term):
    """
    Allow permutations within the same permutation group.
    """
    pos_list = get_positions(term)
    group_dict = dict()
    for pos in pos_list:
        g_pos = drop_tag_last_subscript(pos)
        if g_pos != pos:
            if g_pos not in group_dict:
                group_dict[g_pos] = [ pos, ]
            else:
                group_dict[g_pos].append(pos)
    return group_dict.values()

def mk_combined_permuting_dicts(term):
    """
    Produce all possible permutations allowed by permuting within the permutation groups.
    Generator of dict, each dict is a map of original coordinate and the permuted coordinate.
    """
    pos_list_list = list(get_position_permutation_groups(term))
    p_dicts_list = [ list(mk_permuting_dicts(pos_list)) for pos_list in pos_list_list ]
    def loop(level):
        if level < 0:
            yield dict()
        else:
            for d1 in loop(level - 1):
                for d2 in p_dicts_list[level]:
                    d = d1.copy()
                    d.update(d2)
                    yield d
    return loop(len(p_dicts_list) - 1)

def loop_term_ops(type_dict, ops):
    for op in ops:
        if op.otype == "S":
            # diagram type elem definition
            increase_type_dict_count(type_dict, (op.p1, op.p2,))
        elif op.otype == "Tr":
            loop_term_ops(type_dict, op.ops)
        elif op.otype == "Chain":
            loop_term_ops(type_dict, op.ops)
        elif op.otype == "BS":
            loop_term_ops(type_dict, op.chain_list)

def get_term_diagram_type_info_no_permutation(term):
    type_dict = dict()
    loop_term_ops(type_dict, term.c_ops)
    return tuple(sorted(type_dict.items()))

def permute_tag(x, p_dict):
    if x in p_dict:
        return p_dict[x]
    else:
        return x

def permute_type_info_entry(x, p_dict):
    ((p1, p2,), n,) = x
    p1 = permute_tag(p1, p_dict)
    p2 = permute_tag(p2, p_dict)
    return ((p1, p2,), n,)

def permute_type_info(type_info, p_dict):
    l = [ permute_type_info_entry(e, p_dict) for e in type_info ]
    return tuple(sorted(l))

def get_term_diagram_type_info(term):
    base_type_info = get_term_diagram_type_info_no_permutation(term)
    min_type_info = base_type_info
    min_type_info_repr = str(min_type_info)
    for p_dict in mk_combined_permuting_dicts(term):
        type_info = permute_type_info(base_type_info, p_dict)
        type_info_repr = repr(type_info)
        if min_type_info is None or type_info_repr < min_type_info_repr:
            min_type_info = type_info
            min_type_info_repr = type_info_repr
    return min_type_info

def filter_diagram_type(expr, diagram_type_dict=None, included_types=None):
    """
    first: drop diagrams with diagram_type_dict[diagram_type] == None
    second:
        if included_types is None:
            return a list of a single expr with all the remaining diagrams summed together
        else:
            assert isinstance(included_types, list)
            # included_types = [ None, "Type1", [ "Type2", "Type3", ], ]
            return a list of exprs, each expr only includes the types specified in the `included_types`.
    `included_types` is a list of specs of included diagram type.
    Each spec in the list of `included_types` can be
    (1) None: means all remaining types included
    (2) a single str, with value be one diagram_type_name: means only include this type
    (3) a list of str: means include only the types listed in the list.
    """
    if diagram_type_dict is None:
        return [ expr, ]
    if included_types is None:
        included_types_list = [ None, ]
    else:
        assert isinstance(included_types, list)
        included_types_list = []
        for its in included_types:
            if its is None:
                included_types_list.append(its)
            elif isinstance(its, str):
                included_types_list.append([ its, ])
            elif isinstance(its, list):
                included_types_list.append(its)
            else:
                assert False
    expr_terms_list = [ [] for its in included_types_list ]
    for term in expr.terms:
        diagram_type = get_term_diagram_type_info(term)
        if diagram_type in diagram_type_dict:
            diagram_type_name = diagram_type_dict.get(diagram_type)
            if diagram_type_name is not None:
                for i, its in enumerate(included_types_list):
                    if (its is None) or (diagram_type_name in its):
                        expr_terms_list[i].append(term)
        else:
            for i, its in enumerate(included_types_list):
                if its is None:
                    expr_terms_list[i].append(term)
    expr_list = []
    for i, its in enumerate(included_types_list):
        if its is None:
            its_tag = ""
        else:
            its_tag = " (" + ','.join(its) + ")"
        expr_list.append(Expr(expr_terms_list[i], expr.description + its_tag))
    return expr_list

def mk_cexpr(*exprs, diagram_type_dict=None):
    """
    interface function
    exprs already finished wick contraction,
    otherwise use `contract_simplify_compile(*exprs, is_isospin_symmetric_limit, diagram_type_dict)`
    !!!if diagram_type_dict[diagram_type] == None: this diagram_type should have already be dropped!!!
    """
    if diagram_type_dict is None:
        diagram_type_dict = dict()
    descriptions = [ expr.show() for expr in exprs ]
    # build diagram_types and term names
    diagram_type_counter = 0
    diagram_type_term_dict = dict() # diagram_type_term_dict[repr_term] = diagram_type_name
    term_name_dict = dict() # term_name_dict[term_name] = term
    term_dict = dict() # term_dict[repr(term)] = term_name
    for expr in exprs:
        for term_coef in expr.terms:
            term = Term(term_coef.c_ops, term_coef.a_ops, 1)
            repr_term = repr(term)
            if repr_term in diagram_type_term_dict:
                continue
            diagram_type = get_term_diagram_type_info(term)
            if diagram_type not in diagram_type_dict:
                diagram_type_name = f"ADT{diagram_type_counter}" # ADT is short for "auto diagram type"
                diagram_type_counter += 1
                diagram_type_dict[diagram_type] = diagram_type_name
            diagram_type_name = diagram_type_dict[diagram_type]
            diagram_type_term_dict[repr_term] = diagram_type_name
            assert diagram_type_name is not None
            term_name_counter = 0
            while True:
                term_name = f"term_{diagram_type_name}_{term_name_counter}"
                term_name_counter += 1
                if term_name not in term_name_dict:
                    break
            term_name_dict[term_name] = term
            term_dict[repr_term] = term_name
    # name diagram_types
    diagram_types = []
    for diagram_type, diagram_type_name in diagram_type_dict.items():
        diagram_types.append((diagram_type_name, diagram_type,))
    # name terms
    named_terms = []
    for term_name, term in sorted(term_name_dict.items()):
        named_terms.append((term_name, term,))
    # name exprs
    named_exprs = []
    for i, expr in enumerate(exprs):
        expr_list = []
        typed_expr_list_dict = { name: [] for name, diagram_type in diagram_types }
        for j, term_coef in enumerate(expr.terms):
            coef = term_coef.coef
            term = Term(term_coef.c_ops, term_coef.a_ops, 1)
            repr_term = repr(term)
            diagram_type_name = diagram_type_term_dict[repr_term]
            assert diagram_type_name is not None
            term_name = term_dict[repr_term]
            typed_expr_list_dict[diagram_type_name].append((coef, term_name,))
            expr_list.append((coef, term_name,))
        named_exprs.append((f"{descriptions[i]}  exprs[{i}]", expr_list,))
    # positions
    positions = collect_position_in_cexpr(named_terms, named_exprs)
    # cexpr
    cexpr = CExpr()
    cexpr.diagram_types = diagram_types
    cexpr.positions = positions
    cexpr.named_terms = named_terms
    cexpr.named_exprs = named_exprs
    return cexpr

@q.timer
def contract_simplify(*exprs, is_isospin_symmetric_limit=True, diagram_type_dict=None):
    """
    interface function
    `exprs = [ expr, (expr, *included_types,), ... ]`
    #
    In case `diagram_type_dict` is not `None`, perform the following filter
    If `diagram_type_dict[diagram_type]` is `None`: term is removed.
    For `(expr, *included_types,)`, only terms (in `expr`) with `diagram_type` that is in `included_types` is kept included_types should be a list/tuple of string.
    See `filter_diagram_type` for more details.
    """
    def func(expr):
        expr = copy.deepcopy(expr)
        if isinstance(expr, tuple):
            expr, *included_types = expr
        elif isinstance(expr, Expr):
            included_types = None
        else:
            assert False
        expr = contract_expr(expr)
        expr.simplify(is_isospin_symmetric_limit=is_isospin_symmetric_limit)
        expr_list = filter_diagram_type(
                expr,
                diagram_type_dict=diagram_type_dict,
                included_types=included_types)
        return expr_list
    expr_list_list = q.parallel_map(func, exprs)
    expr_list = []
    for el in expr_list_list:
        expr_list += el
    return expr_list

@q.timer
def compile_expr(*exprs, diagram_type_dict=None):
    """
    interface function
    """
    exprs = [ expr.copy() for expr in exprs ]
    cexpr = mk_cexpr(*exprs, diagram_type_dict=diagram_type_dict)
    return cexpr

@q.timer
def contract_simplify_compile(*exprs, is_isospin_symmetric_limit=True, diagram_type_dict=None):
    """
    interface function
    Call `contract_simplify` and then `compile_expr`
    #
    This function can be used to construct the first argument of `cached_comipled_cexpr`.
    `cached_comipled_cexpr` will call `cexpr.optimize()`.
    #
    e.g. exprs = [ Qb("u", "x", s, c) * Qv("u", "x", s, c) + "u_bar*u", Qb("s", "x", s, c) * Qv("s", "x", s, c) + "s_bar*s", Qb("c", "x", s, c) * Qv("c", "x", s, c) + "c_bar*c", ]
    e.g. exprs = [ mk_pi_p("x2", True) * mk_pi_p("x1") + "(pi   * pi)", mk_j5pi_mu("x2", 3) * mk_pi_p("x1") + "(a_pi * pi)", mk_k_p("x2", True)  * mk_k_p("x1")  + "(k    * k )", mk_j5k_mu("x2", 3)  * mk_k_p("x1")  + "(a_k  * k )", ]
    """
    contracted_simplified_exprs = contract_simplify(
            *exprs,
            is_isospin_symmetric_limit=is_isospin_symmetric_limit,
            diagram_type_dict=diagram_type_dict)
    cexpr = compile_expr(*contracted_simplified_exprs, diagram_type_dict=diagram_type_dict)
    return cexpr

def show_variable_value(value):
    if isinstance(value, list):
        return "*".join(map(show_variable_value, value))
    elif isinstance(value, Var):
        return f"{value.name}"
    elif isinstance(value, G) and value.tag in [ 0, 1, 2, 3, 5, ]:
        tag = { 0: "x", 1: "y", 2: "z", 3: "t", 5: "5", }[value.tag]
        if value.s1 == "auto" and value.s2 == "auto":
            return f"gamma_{tag}"
        else:
            return f"gamma_{tag}({value.s1},{value.s2})"
    elif isinstance(value, G):
        if value.s1 == "auto" and value.s2 == "auto":
            return f"gamma({value.tag})"
        else:
            return f"gamma({value.tag},{value.s1},{value.s2})"
    elif isinstance(value, U):
        if value.c1 == "auto" and value.c2 == "auto":
            return f"U({value.tag},{value.p},{value.mu})"
        else:
            return f"U({value.tag},{value.p},{value.mu},{value.c1},{value.c2})"
    elif isinstance(value, S):
        if value.s1 == "auto" and value.s2 == "auto" and value.c1 == "auto" and value.c2 == "auto":
            return f"S_{value.f}({value.p1},{value.p2})"
        else:
            return f"S_{value.f}({value.p1},{value.p2},{value.s1},{value.s2},{value.c1},{value.c2})"
    elif isinstance(value, Tr):
        expr = "*".join(map(show_variable_value, value.ops))
        return f"tr({expr})"
    elif isinstance(value, Chain):
        expr = "*".join(map(show_variable_value, value.ops))
        return f"chain({expr})"
    elif isinstance(value, BS):
        elem_list = value.elem_list
        elem_list_expr = ",".join([ f"({v},{c!r},)" for v, c, in elem_list ])
        chain_list_expr = ",".join(map(show_variable_value, value.chain_list))
        return f"bs([{elem_list_expr}],{chain_list_expr})"
    elif isinstance(value, Term):
        if value.coef == 1:
            return "*".join(map(show_variable_value, value.c_ops + value.a_ops))
        else:
            return "*".join(map(show_variable_value, [ f"({value.coef})", ] + value.c_ops + value.a_ops))
    elif isinstance(value, tuple) and len(value) == 2:
        return f"{show_variable_value(value[0])}*{show_variable_value(value[1])}"
    else:
        return f"{value}"

def display_cexpr(cexpr:CExpr):
    """
    interface function
    return a string
    """
    lines = []
    lines.append(f"# Begin CExpr")
    if cexpr.diagram_types:
        lines.append(f"diagram_type_dict = dict()")
        for name, diagram_type in cexpr.diagram_types:
            lines.append(f"diagram_type_dict[{diagram_type}] = {name!r}")
    if cexpr.positions:
        position_vars = ", ".join(cexpr.positions)
        lines.append(f"# Positions:")
        lines.append(f"{position_vars} = {cexpr.positions}")
    if cexpr.variables_prop:
        lines.append(f"# Variables prop:")
    for name, value in cexpr.variables_prop:
        lines.append(f"{name:<30} = {show_variable_value(value)}")
    if cexpr.variables_color_matrix:
        lines.append(f"# Variables color matrix:")
    for name, value in cexpr.variables_color_matrix:
        lines.append(f"{name:<30} = {show_variable_value(value)}")
    if cexpr.variables_factor_intermediate:
        lines.append(f"# Variables factor intermediate:")
    for name, value in cexpr.variables_factor_intermediate:
        lines.append(f"{name:<30} = {show_variable_value(value)}")
    if cexpr.variables_factor:
        lines.append(f"# Variables factor:")
    for name, value in cexpr.variables_factor:
        lines.append(f"{name:<30} = {show_variable_value(value)}")
    if cexpr.variables_prod:
        lines.append(f"# Variables prod:")
    for name, value in cexpr.variables_prod:
        lines.append(f"{name:<30} = {show_variable_value(value)}")
    if cexpr.variables_chain:
        lines.append(f"# Variables chain:")
    for name, value in cexpr.variables_chain:
        lines.append(f"{name:<30} = {show_variable_value(value)}")
    if cexpr.variables_tr:
        lines.append(f"# Variables tr:")
    for name, value in cexpr.variables_tr:
        lines.append(f"{name:<30} = {show_variable_value(value)}")
    if cexpr.variables_baryon_prop:
        lines.append(f"# Variables baryon_prop:")
    for name_list, value_list in cexpr.variables_baryon_prop:
        for name, value in zip(name_list, value_list):
            lines.append(f"{name:<30} = {show_variable_value(value)}")
    if cexpr.diagram_types:
        lines.append(f"# Diagram type coef:")
    for name, diagram_type in cexpr.diagram_types:
        if name is not None:
            coef_name = f"coef_{name}"
            lines.append(f"{coef_name:<30} = 1")
    if cexpr.named_terms:
        lines.append(f"# Named terms:")
    for idx, (name, term,) in enumerate(cexpr.named_terms):
        name_type = "_".join([ "coef", ] + name.split("_")[1:-1])
        lines.append(f"{name:<30} = {name_type} * {show_variable_value(term)}")
    lines.append(f"terms = [ 0 for i in range({len(cexpr.named_terms)}) ]")
    for idx, (name, term,) in enumerate(cexpr.named_terms):
        lines.append(f"terms[{idx}] = {name}")
    if cexpr.named_exprs:
        lines.append(f"# Named exprs:")
    lines.append(f"exprs = [ 0 for i in range({len(cexpr.named_exprs)}) ]")
    for idx, (name, expr,) in enumerate(cexpr.named_exprs):
        lines.append(f"# {name}")
        for e in expr:
            lines.append(f"exprs[{idx}] += {show_variable_value(e)}")
    lines.append(f"# End CExpr")
    return "\n".join(lines)

@q.timer_verbose
def cexpr_code_gen_py(cexpr:CExpr, *, is_cython=True, is_distillation=False):
    """
    interface function
    return a string
    #
    if is_distillation:
        assert is_cython == False
    """
    gen = CExprCodeGenPy(cexpr,
                         is_cython=is_cython,
                         is_distillation=is_distillation)
    return gen.code_gen()

class CExprCodeGenPy:

    """
    self.cexpr
    self.is_cython
    self.is_distillation
    self.var_dict_for_factors
    self.lines
    self.indent
    self.total_sloppy_flops
    #
    flops per complex addition: 2
    flops per complex multiplication: 6
    flops per matrix multiplication: 6 M N L + 2 M L (N-1) ==> 13536 (sc * sc), 4320 (sc * s), 480 (s * s), 3168 (sc * c), 198 (c * c)
    flops per trace 2 (M-1) ==> 22 (sc)
    flops per trace2 6 M N + 2 (M N - 1) ==> 1150 (sc, sc)
    """

    def __init__(self, cexpr, *, is_cython=True, is_distillation=False):
        self.cexpr = cexpr
        self.is_cython = is_cython
        self.is_distillation = is_distillation
        self.var_dict_for_factors = {}
        self.lines = []
        self.indent = 0
        self.total_sloppy_flops = 0
        #
        if self.is_distillation:
            assert self.is_cython == False
        #
        # self.set_var_dict_for_factors();

    def code_gen(self):
        """
        main function
        """
        lines = self.lines
        append = self.append
        append_cy = self.append_cy
        append_py = self.append_py
        if self.is_distillation:
            append(f"from auto_contractor.runtime_distillation import *")
        else:
            append(f"from auto_contractor.runtime import *")
        append_cy(f"import cython")
        append_cy(f"cimport qlat_utils.everything as cc")
        append_cy(f"cimport qlat_utils.all as qu")
        append_cy(f"cimport libcpp.complex")
        append_cy(f"cimport numpy")
        self.sep()
        self.cexpr_function()
        self.sep()
        self.cexpr_function_get_prop()
        self.sep()
        self.cexpr_function_eval()
        self.sep()
        self.cexpr_function_eval_with_props()
        self.sep()
        self.cexpr_function_bs_eval()
        self.sep()
        self.total_flops()
        return "\n".join(lines)

    def set_var_dict_for_factors(self):
        """
        If this function is called (currently not),
        factor will be referred use array and index,
        instead of its variable name.
        """
        cexpr = self.cexpr
        self.var_dict_for_factors = {}
        var_dict_for_factors = self.var_dict_for_factors
        for idx, (name, value,) in enumerate(cexpr.variables_factor_intermediate):
            assert name.startswith("V_factor_")
            var_dict_for_factors[name] = f"factors_intermediate_view[{idx}]"
        for idx, (name, value,) in enumerate(cexpr.variables_factor):
            assert name.startswith("V_factor_")
            var_dict_for_factors[name] = f"factors_view[{idx}]"

    def gen_expr(self, x):
        """
        return code_str, type_str
        `type_str` follows convention of `get_var_name_type`
        """
        if isinstance(x, (int, float, complex)):
            return f"{x}", "V_a"
        elif isinstance(x, (ea.Expr, ea.Factor)):
            return f"({ea.compile_py(x, self.var_dict_for_factors)})", "V_a"
        elif isinstance(x, tuple) and len(x) == 2 and isinstance(x[1], list) and len(x[1]) > 0 and isinstance(x[1][0], BS):
            name_list, value_list, = x
            for name in name_list:
                assert isinstance(name, str)
            for value in value_list:
                assert isinstance(value, BS)
            bs_list = value_list
            factor_list = get_bs_factor_variable_list(bs_list)
            chain_list = [ ch.name for ch in bs_list[0].chain_list ]
            for bs in bs_list:
                assert chain_list == [ ch.name for ch in bs.chain_list ]
            chain_list_uniq = sorted(set(chain_list))
            arg_list_cy = []
            arg_list_py = []
            for name in name_list:
                arg_list_cy.append(f"&{name}")
            for c in chain_list_uniq:
                arg_list_cy.append(f"&{c}")
                arg_list_py.append(f"{c}")
            for f in factor_list:
                arg_list_cy.append(f"{f}")
                arg_list_py.append(f"{f}")
            arg_str_cy = ", ".join(arg_list_cy)
            arg_str_py = ", ".join(arg_list_py)
            c_cy = f"cexpr_function_bs_eval_{name_list[0]}({arg_str_cy})"
            c_py = f"cexpr_function_bs_eval_{name_list[0]}({arg_str_py})"
            if self.is_cython:
                return c_cy, "None"
            else:
                return c_py, "list[V_a]"
        assert isinstance(x, Op)
        if x.otype == "S":
            return f"get_prop('{x.f}', {x.p1}, {x.p2})", "V_S"
        elif x.otype == "U":
            return f"get_prop('U', '{x.tag}', {x.p}, {x.mu})", "V_U"
        elif x.otype == "G":
            assert x.s1 == "auto" and x.s2 == "auto"
            assert x.tag in [ 0, 1, 2, 3, 5, ] or isinstance(x.tag, str)
            if self.is_cython:
                if x.tag in [ 0, 1, 2, 3, 5, ]:
                    return f"qu.gamma_matrix_{x.tag}", "V_G"
                else:
                    return f"cc.get_gamma_matrix({x.tag})", "V_G"
            else:
                return f"get_gamma_matrix({x.tag})", "V_G"
        elif x.otype == "Chain":
            assert x.tag == "sc"
            if len(x.ops) == 0:
                assert False
            else:
                c, t = self.gen_expr_prod_list(x.ops)
                assert t == "V_S"
                return c, t
        elif x.otype == "Tr":
            assert x.tag == "sc"
            if len(x.ops) == 0:
                assert False
            elif len(x.ops) == 1:
                c, t = self.gen_expr(x.ops[0])
                assert t == "V_S"
                self.total_sloppy_flops += 22
                if self.is_cython:
                    return f"cc.matrix_trace({c})", "V_a"
                else:
                    return f"mat_tr_wm({c})", "V_a"
            else:
                c1, t1 = self.gen_expr_prod_list(x.ops[:-1])
                c2, t2 = self.gen_expr(x.ops[-1])
                if t1 == "V_S" and t2 == "V_S":
                    self.total_sloppy_flops += 1150
                    if self.is_cython:
                        return f"cc.matrix_trace({c1}, {c2})", "V_a"
                    else:
                        return f"mat_tr_wm_wm({c1}, {c2})", "V_a"
                elif t1 == "V_S" and t2 == "V_G":
                    if self.is_cython:
                        return f"cc.matrix_trace({c1}, {c2})", "V_a"
                    else:
                        return f"mat_tr_wm_sm({c1}, {c2})", "V_a"
                elif t1 == "V_G" and t2 == "V_S":
                    if self.is_cython:
                        return f"cc.matrix_trace({c1}, {c2})", "V_a"
                    else:
                        return f"mat_tr_sm_wm({c1}, {c2})", "V_a"
                elif t1 == "V_S" and t2 == "V_U":
                    if self.is_cython:
                        return f"cc.matrix_trace({c1}, {c2})", "V_a"
                    else:
                        return f"mat_tr_wm_cm({c1}, {c2})", "V_a"
                elif t1 == "V_U" and t2 == "V_S":
                    if self.is_cython:
                        return f"cc.matrix_trace({c1}, {c2})", "V_a"
                    else:
                        return f"mat_tr_cm_wm({c1}, {c2})", "V_a"
                else:
                    assert False
        elif x.otype == "Var":
            if x.name.startswith("V_S_") or x.name.startswith("V_U_"):
                if self.is_cython:
                    return f"p_{x.name}[0]", get_var_name_type(x.name)
                else:
                    return f"p_{x.name}", get_var_name_type(x.name)
            else:
                return f"{x.name}", get_var_name_type(x.name)
        raise Exception(f"gen_expr: x='{x}'")

    def gen_expr_prod(self, ct1, ct2):
        """
        return code_str, type_str
        `type_str` follows convention of `get_var_name_type`
        """
        c1, t1 = ct1
        c2, t2 = ct2
        if c1 == "(1+0j)":
            assert t1 == "V_a"
            return ct2
        elif c2 == "(1+0j)":
            assert t2 == "V_a"
            return ct1
        elif t1 == "V_S" and t2 == "V_S":
            self.total_sloppy_flops += 13536
            if self.is_cython:
                return f"{c1} * {c2}", "V_S"
            else:
                return f"mat_mul_wm_wm({c1}, {c2})", "V_S"
        #
        elif t1 == "V_S" and t2 == "V_a":
            if self.is_cython:
                return f"{c1} * {c2}", "V_S"
            else:
                return f"mat_mul_a_wm({c2}, {c1})", "V_S"
        elif t1 == "V_a" and t2 == "V_S":
            if self.is_cython:
                return f"{c1} * {c2}", "V_S"
            else:
                return f"mat_mul_a_wm({c1}, {c2})", "V_S"
        elif t1 == "V_a" and t2 == "V_a":
            return f"{c1} * {c2}", "V_a"
        #
        elif t1 == "V_S" and t2 == "V_G":
            self.total_sloppy_flops += 4320
            if self.is_cython:
                return f"{c1} * {c2}", "V_S"
            else:
                return f"mat_mul_wm_sm({c1}, {c2})", "V_S"
        elif t1 == "V_G" and t2 == "V_S":
            self.total_sloppy_flops += 4320
            if self.is_cython:
                return f"{c1} * {c2}", "V_S"
            else:
                return f"mat_mul_sm_wm({c1}, {c2})", "V_S"
        elif t1 == "V_G" and t2 == "V_G":
            self.total_sloppy_flops += 480
            if self.is_cython:
                return f"{c1} * {c2}", "V_G"
            else:
                return f"mat_mul_sm_sm({c1}, {c2})", "V_G"
        elif t1 == "V_G" and t2 == "V_a":
            if self.is_cython:
                return f"{c1} * {c2}", "V_G"
            else:
                return f"mat_mul_a_sm({c2}, {c1})", "V_G"
        elif t1 == "V_a" and t2 == "V_G":
            if self.is_cython:
                return f"{c1} * {c2}", "V_G"
            else:
                return f"mat_mul_a_sm({c1}, {c2})", "V_G"
        #
        elif t1 == "V_S" and t2 == "V_U":
            self.total_sloppy_flops += 3168
            if self.is_cython:
                return f"{c1} * {c2}", "V_S"
            else:
                return f"mat_mul_wm_cm({c1}, {c2})", "V_S"
        elif t1 == "V_U" and t2 == "V_S":
            self.total_sloppy_flops += 3168
            if self.is_cython:
                return f"{c1} * {c2}", "V_S"
            else:
                return f"mat_mul_cm_wm({c1}, {c2})", "V_S"
        elif t1 == "V_U" and t2 == "V_U":
            self.total_sloppy_flops += 198
            if self.is_cython:
                return f"{c1} * {c2}", "V_U"
            else:
                return f"mat_mul_cm_cm({c1}, {c2})", "V_U"
        elif t1 == "V_U" and t2 == "V_a":
            if self.is_cython:
                return f"{c1} * {c2}", "V_U"
            else:
                return f"mat_mul_a_cm({c2}, {c1})", "V_U"
        elif t1 == "V_a" and t2 == "V_U":
            if self.is_cython:
                return f"{c1} * {c2}", "V_U"
            else:
                return f"mat_mul_a_cm({c1}, {c2})", "V_U"
        #
        else:
            raise Exception(f"gen_expr_prod: ct1='{ct1}' ; ct1='{ct1}'")

    def gen_expr_prod_list(self, x_list):
        """
        return code_str, type_str
        `type_str` follows convention of `get_var_name_type`
        """
        if len(x_list) == 0:
            return f"1", "V_a"
        elif len(x_list) == 1:
            return self.gen_expr(x_list[0])
        else:
            assert len(x_list) > 1
            return self.gen_expr_prod(self.gen_expr_prod_list(x_list[:-1]), self.gen_expr(x_list[-1]))

    def append(self, line):
        lines = self.lines
        if line == "":
            lines.append("")
        else:
            lines.append(self.indent * ' ' + line)

    def append_py(self, line):
        if self.is_cython:
            return
        lines = self.lines
        if line == "":
            lines.append("# Python only")
        else:
            lines.append(self.indent * ' ' + line + " # Python only")

    def append_cy(self, line):
        if not self.is_cython:
            return
        lines = self.lines
        if line == "":
            lines.append("# Cython")
        else:
            lines.append(self.indent * ' ' + line + " # Cython")

    def sep(self):
        append = self.append
        append(f"")
        append(f"### ----")
        append(f"")

    def cexpr_function(self):
        append = self.append
        append_cy = self.append_cy
        append_py = self.append_py
        append(f"@timer")
        append(f"def cexpr_function(*, positions_dict, get_prop, is_ama_and_sloppy=False):")
        self.indent += 4
        append(f"# get_props")
        append(f"props, cms, factors = cexpr_function_get_prop(positions_dict, get_prop)")
        append(f"# eval")
        append(f"ama_val = cexpr_function_eval(positions_dict, props, cms, factors)")
        append(f"# extract sloppy val")
        append(f"val_sloppy = ama_extract(ama_val, is_sloppy=True)")
        append(f"# extract AMA val")
        append(f"val_ama = ama_extract(ama_val)")
        append(f"# return")
        append(f"if is_ama_and_sloppy:")
        append(f"    # return both AMA corrected results and sloppy results")
        append(f"    return val_ama, val_sloppy")
        append(f"else:")
        append(f"    # return AMA corrected results by default")
        append(f"    return val_ama")
        self.indent -= 4

    def cexpr_function_get_prop(self):
        append = self.append
        append_cy = self.append_cy
        append_py = self.append_py
        cexpr = self.cexpr
        append(f"@timer_flops")
        append_cy(f"@cython.boundscheck(False)")
        append_cy(f"@cython.wraparound(False)")
        append(f"def cexpr_function_get_prop(positions_dict, get_prop):")
        self.indent += 4
        append(f"# set positions")
        for position_var in cexpr.positions:
            if position_var in aff.auto_fac_funcs_list:
                append(f"{position_var} = aff.{position_var}")
            else:
                append(f"{position_var} = positions_dict['{position_var}']")
        append(f"# get prop")
        for name, value in cexpr.variables_prop:
            assert name.startswith("V_S_")
            x = value
            assert isinstance(x, Op)
            assert x.otype == "S"
            c, t = self.gen_expr(x)
            assert t == "V_S"
            # potentially be AMA prop (cannot use cdef in cython)
            append(f"{name} = {c}")
        append(f"# get color matrix")
        for name, value in cexpr.variables_color_matrix:
            assert name.startswith("V_U_")
            x = value
            assert isinstance(x, Op)
            assert x.otype == "U"
            c, t = self.gen_expr(x)
            assert t == "V_U"
            append_cy(f"cdef qu.ColorMatrix {name} = {c}")
            append_py(f"{name} = {c}")
        append(f"# set props for return")
        append(f"props = [")
        self.indent += 4
        for name, value in cexpr.variables_prop:
            append(f"{name},")
        append(f"]")
        self.indent -= 4
        append(f"# set color matrix for return")
        append(f"cms = [")
        self.indent += 4
        for name, value in cexpr.variables_color_matrix:
            append(f"{name},")
        append(f"]")
        self.indent -= 4
        append(f"# set intermediate factors")
        for idx, (name, value,) in enumerate(cexpr.variables_factor_intermediate):
            assert name.startswith("V_factor_")
            x = value
            assert isinstance(x, ea.Expr)
            c, t = self.gen_expr(x)
            assert t == "V_a"
            append(f"# {idx} {name}")
            append_cy(f"cdef cc.PyComplexD {name} = {c}")
            append_py(f"{name} = {c}")
        append(f"# declare factors")
        append_cy(f"cdef numpy.ndarray[numpy.complex128_t] factors")
        append(f"factors = np.zeros({len(cexpr.variables_factor)}, dtype=np.complex128)")
        append_cy(f"cdef cc.PyComplexD[:] factors_view = factors")
        append_py(f"factors_view = factors")
        append(f"# set factors")
        for idx, (name, value,) in enumerate(cexpr.variables_factor):
            assert name.startswith("V_factor_")
            x = value
            assert isinstance(x, ea.Expr)
            c, t = self.gen_expr(x)
            assert t == "V_a"
            append(f"# {name}")
            append_cy(f"cdef cc.PyComplexD {name} = {c}")
            append_py(f"{name} = {c}")
            append(f"factors_view[{idx}] = {name}")
        append(f"# set flops")
        append(f"total_flops = len(props) * 144 * 2 * 8 + len(cms) * 9 * 2 * 8 + len(factors) * 2 * 8")
        append(f"# return")
        append(f"return total_flops, (props, cms, factors,)")
        self.indent -= 4

    def cexpr_function_eval(self):
        append = self.append
        append_cy = self.append_cy
        append_py = self.append_py
        cexpr = self.cexpr
        append(f"@timer_flops")
        append(f"def cexpr_function_eval(positions_dict, props, cms, factors):")
        self.indent += 4
        append(f"# load AMA props with proper format")
        append(f"props = [ load_prop(p) for p in props ]")
        append(f"# join the AMA props")
        append(f"ama_props = ama_list(*props)")
        append(f"# apply eval to the factors and AMA props")
        append(f"ama_val = ama_apply1(lambda x_props: cexpr_function_eval_with_props(positions_dict, x_props, cms, factors), ama_props)")
        append(f"# set flops")
        append(f"total_flops = ama_counts(ama_val) * total_sloppy_flops")
        append(f"# return")
        append(f"return total_flops, ama_val")
        self.indent -= 4

    def cexpr_function_bs_eval(self):
        append = self.append
        append_cy = self.append_cy
        append_py = self.append_py
        cexpr = self.cexpr
        for name_list, value_list in cexpr.variables_baryon_prop:
            for name in name_list:
                assert isinstance(name, str)
            for value in value_list:
                assert isinstance(value, BS)
            bs_list = value_list
            factor_list = get_bs_factor_variable_list(bs_list)
            chain_list = [ ch.name for ch in bs_list[0].chain_list ]
            for bs in bs_list:
                assert chain_list == [ ch.name for ch in bs.chain_list ]
            chain_list_uniq = sorted(set(chain_list))
            arg_list_cy = []
            arg_list_py = []
            for name in name_list:
                arg_list_cy.append(f"cc.PyComplexD* {name}")
            for c in chain_list_uniq:
                arg_list_cy.append(f"cc.WilsonMatrix* {c}")
                arg_list_py.append(f"{c}")
            for f in factor_list:
                arg_list_cy.append(f"cc.PyComplexD {f}")
                arg_list_py.append(f"{f}")
            arg_str_cy = ", ".join(arg_list_cy)
            arg_str_py = ", ".join(arg_list_py)
            append_cy(f"@cython.boundscheck(False)")
            append_cy(f"@cython.wraparound(False)")
            append_cy(f"cdef void cexpr_function_bs_eval_{name_list[0]}({arg_str_cy}):")
            append_py(f"def cexpr_function_bs_eval_{name_list[0]}({arg_str_py}):")
            self.indent += 4
            append_cy(f"cdef cc.PyComplexD v")
            for name in name_list:
                append_cy(f"cdef cc.PyComplexD {name}_v = 0")
                append_py(f"{name} = 0")
            ch1, ch2, ch3, = chain_list
            def get_ss_tag(ss):
                return "_".join([ str(s) for s in ss ])
            ss_dict = dict()
            for bs in bs_list:
                elem_list = bs.get_spin_spin_tensor_elem_list_code()
                for ss, coef, in elem_list:
                    if ss in ss_dict:
                        continue
                    ss_name = f"ss_res_{get_ss_tag(ss)}"
                    ss_dict[ss] = ss_name
                    v_s1, b_s1, v_s2, b_s2, v_s3, b_s3, = ss
                    append_cy(f"cdef cc.PyComplexD {ss_name} = cc.pycc_d(cc.epsilon_contraction({v_s1}, {b_s1}, {v_s2}, {b_s2}, {v_s3}, {b_s3}, {ch1}[0], {ch2}[0], {ch3}[0]))")
                    append_py(f"{ss_name} = mat_epsilon_contraction_wm_wm_wm({v_s1}, {b_s1}, {v_s2}, {b_s2}, {v_s3}, {b_s3}, {ch1}, {ch2}, {ch3})")
                    self.total_sloppy_flops += 504 # 6 * 6 * (2 * 6 + 2)
            for name, bs in zip(name_list, bs_list):
                elem_list = bs.get_spin_spin_tensor_elem_list_code()
                for ss, coef, in elem_list:
                    ss_name = f"ss_res_{get_ss_tag(ss)}"
                    c, t, = self.gen_expr(coef)
                    assert t == "V_a"
                    append(f"v = {c}")
                    append_cy(f"{name}_v += v * {ss_name}")
                    append_py(f"{name} += v * {ss_name}")
            for name in name_list:
                append_cy(f"{name}[0] = {name}_v")
            name_list_str = ", ".join(name_list)
            append_py(f"return {name_list_str}")
            self.indent -= 4
            append(f"")

    def cexpr_function_eval_with_props(self):
        append = self.append
        append_cy = self.append_cy
        append_py = self.append_py
        cexpr = self.cexpr
        append(f"@timer_flops")
        append_cy(f"@cython.boundscheck(False)")
        append_cy(f"@cython.wraparound(False)")
        append_cy(f"def cexpr_function_eval_with_props(dict positions_dict, list props, list cms, cc.PyComplexD[:] factors_view):")
        append_py(f"def cexpr_function_eval_with_props(positions_dict, props, cms, factors_view):")
        self.indent += 4
        append(f"# set positions")
        for position_var in cexpr.positions:
            if position_var in aff.auto_fac_funcs_list:
                append(f"{position_var} = aff.{position_var}")
            else:
                append(f"{position_var} = positions_dict['{position_var}']")
        append(f"# set props")
        for idx, (name, value,) in enumerate(cexpr.variables_prop):
            append_cy(f"cdef cc.WilsonMatrix* p_{name} = &(<qu.WilsonMatrix>props[{idx}]).xx")
            append_py(f"p_{name} = props[{idx}]")
        append(f"# set cms")
        for idx, (name, value,) in enumerate(cexpr.variables_color_matrix):
            append_cy(f"cdef cc.ColorMatrix* p_{name} = &(<qu.ColorMatrix>cms[{idx}]).xx")
            append_py(f"p_{name} = cms[{idx}]")
        append(f"# set factors")
        for idx, (name, value,) in enumerate(cexpr.variables_factor):
            append_cy(f"cdef cc.PyComplexD {name} = factors_view[{idx}]")
            append_py(f"{name} = factors_view[{idx}]")
        append(f"# compute products")
        for name, value in cexpr.variables_prod:
            assert name.startswith("V_prod_")
            x = value
            assert isinstance(x, tuple)
            assert len(x) == 2
            c, t = self.gen_expr_prod_list(x)
            assert t == get_var_name_type(name)
            if t == "V_G":
                append_cy(f"cdef cc.SpinMatrix {name} = {c}")
                append_py(f"{name} = {c}")
            elif t == "V_S":
                append_cy(f"cdef cc.WilsonMatrix {name} = {c}")
                append_py(f"{name} = {c}")
            elif t == "V_U":
                append_cy(f"cdef cc.ColorMatrix {name} = {c}")
                append_py(f"{name} = {c}")
            else:
                assert False
        append(f"# compute chains")
        for name, value in cexpr.variables_chain:
            assert name.startswith("V_chain_")
            x = value
            assert isinstance(x, Op)
            assert x.otype == "Chain"
            c, t = self.gen_expr(x)
            assert t == "V_S"
            append_cy(f"cdef cc.WilsonMatrix {name} = {c}")
            append_py(f"{name} = {c}")
        append(f"# compute traces")
        for name, value in cexpr.variables_tr:
            assert name.startswith("V_tr_")
            x = value
            assert isinstance(x, Op)
            assert x.otype == "Tr"
            c, t = self.gen_expr(x)
            assert t == "V_a"
            append_cy(f"cdef cc.PyComplexD {name} = cc.pycc_d({c})")
            append_py(f"{name} = {c}")
        append(f"# compute baryon_props")
        for name_list, value_list in cexpr.variables_baryon_prop:
            for x in value_list:
                assert isinstance(x, Op)
                assert x.otype == "BS"
            c, t = self.gen_expr((name_list, value_list))
            if self.is_cython:
                assert t == "None"
            else:
                assert t == "list[V_a]"
            for name in name_list:
                append_cy(f"cdef cc.PyComplexD {name} = 0")
            append_cy(c)
            name_list_str = ", ".join(name_list)
            append_py(f"{name_list_str} = {c}")
        append(f"# set terms")
        term_type_dict = dict()
        for name, term in cexpr.named_terms:
            x = term
            if x.coef == 1:
                c_ops = x.c_ops
            else:
                c_ops = [ x.coef, ] + x.c_ops
            c, t = self.gen_expr_prod_list(c_ops)
            term_type_dict[name] = t
            if t == "V_a":
                append_cy(f"cdef cc.PyComplexD {name} = {c}")
                append_py(f"{name} = {c}")
            elif t == "V_G":
                append_cy(f"cdef cc.SpinMatrix {name} = {c}")
                append_py(f"{name} = {c}")
            elif t == "V_S":
                append_cy(f"cdef cc.WilsonMatrix {name} = {c}")
                append_py(f"{name} = {c}")
            elif t == "V_U":
                append_cy(f"cdef cc.ColorMatrix {name} = {c}")
                append_py(f"{name} = {c}")
            else:
                assert False
        is_terms_all_complex_number = all("V_a" == x for x in term_type_dict.values())
        append(f"# declare exprs")
        if is_terms_all_complex_number:
            # Return complex array if all the exprs return complex numbers
            append_cy(f"cdef numpy.ndarray[numpy.complex128_t] exprs")
            append(f"exprs = np.zeros({len(cexpr.named_exprs)}, dtype=np.complex128)")
            append_cy(f"cdef cc.PyComplexD[:] exprs_view = exprs")
            append_py(f"exprs_view = exprs")
        else:
            # Return object array if some exprs return Matrices
            append(f"exprs = np.empty({len(cexpr.named_exprs)}, dtype=object)")
            append(f"exprs_view = exprs")
        append(f"# set exprs")
        def show_coef_term(coef, tname, ttype):
            coef, t = self.gen_expr(coef)
            assert t == "V_a"
            if coef == "1":
                return f"{tname}"
            else:
                if self.is_cython:
                    if ttype == "V_a":
                        return f"({coef}) * {tname}"
                    else:
                        return f"cc.ccpy_d({coef}) * {tname}"
                else:
                    if ttype == "V_a":
                        return f"({coef}) * {tname}"
                    elif ttype == "V_S":
                        return f"mat_mul_a_wm({coef}, {tname})"
                    elif ttype == "V_G":
                        return f"mat_mul_a_sm({coef}, {tname})"
                    elif ttype == "V_U":
                        return f"mat_mul_a_cm({coef}, {tname})"
                    else:
                        raise Exception(f"coef={coef}; tname={tname} ttype={ttype}")
        append_cy(f"cdef cc.PyComplexD expr_V_a")
        if not is_terms_all_complex_number:
            append_cy(f"cdef cc.SpinMatrix expr_V_G")
            append_cy(f"cdef cc.WilsonMatrix expr_V_S")
            append_cy(f"cdef cc.ColorMatrix expr_V_U")
            append_cy(f"cdef qu.SpinMatrix expr_V_G_box")
            append_cy(f"cdef qu.WilsonMatrix expr_V_S_box")
            append_cy(f"cdef qu.ColorMatrix expr_V_U_box")
            append_py(f"expr_V_G = SpinMatrix()")
            append_py(f"expr_V_S = WilsonMatrix()")
            append_py(f"expr_V_U = ColorMatrix()")
        for idx, (name, expr,) in enumerate(cexpr.named_exprs):
            name = name.replace("\n", "  ")
            append(f"# {idx} name='{name}' ")
            if len(expr) == 0:
                append(f"exprs_view[{idx}] = 0")
            else:
                assert len(expr) >= 1
                term_type_set = set(term_type_dict[tname] for _, tname in expr)
                if len(term_type_set) != 1:
                    raise Exception(f"term_type_set={term_type_set}")
                expr_type, = term_type_set
                expr_var_name = f"expr_{expr_type}"
                if expr_type == "V_a":
                    append(f"{expr_var_name} = 0")
                else:
                    append_cy(f"cc.set_zero({expr_var_name})")
                    append_py(f"{expr_var_name}.set_zero()")
                for coef, tname in expr:
                    s = show_coef_term(coef, tname, expr_type)
                    append(f"{expr_var_name} += {s}")
                if expr_type == "V_a":
                    append_cy(f"exprs_view[{idx}] = {expr_var_name}")
                    append_py(f"exprs_view[{idx}] = {expr_var_name}")
                else:
                    if expr_type == "V_S":
                        append_cy(f"{expr_var_name}_box = WilsonMatrix()")
                    elif expr_type == "V_G":
                        append_cy(f"{expr_var_name}_box = SpinMatrix()")
                    elif expr_type == "V_U":
                        append_cy(f"{expr_var_name}_box = ColorMatrix()")
                    else:
                        raise Exception(f"{expr_type}")
                    append_cy(f"{expr_var_name}_box.xx = {expr_var_name}")
                    append_cy(f"exprs_view[{idx}] = {expr_var_name}_box")
                    append_py(f"exprs_view[{idx}] = {expr_var_name}.copy()")
        append(f"# set flops")
        append(f"total_flops = total_sloppy_flops")
        append(f"# return")
        append(f"return total_flops, exprs")
        self.indent -= 4

    def total_flops(self):
        append = self.append
        append(f"# Total flops per sloppy call is: {self.total_sloppy_flops}")
        append(f"total_sloppy_flops = {self.total_sloppy_flops}")

#### ----

def get_bs_factor_variable_list(bs_list:list[BS]) -> list[str]:
    """
    Get the factor variable names used in `bs_list`.
    """
    variables_set = set()
    for bs in bs_list:
        for v, c in bs.elem_list:
            variables_set |= ea.mk_fac(c).get_variable_set()
    return sorted(variables_set)

#### ----

def mk_test_expr_compile_01():
    expr = Qb("d", "x1", "s1", "c1") * G(5, "s1", "s2") * Qv("u", "x1", "s2", "c1") * Qb("u", "x2", "s3", "c2") * G(5, "s3", "s4") * Qv("d", "x2", "s4", "c2")
    return expr

if __name__ == "__main__":
    expr = mk_test_expr_compile_01()
    print(expr)
    print()
    expr = simplified(contract_expr(expr))
    print(expr)
    print()
    cexpr = mk_cexpr(expr).copy()
    print(cexpr)
    print()
    cexpr.optimize()
    print(cexpr)
    print()
    print(display_cexpr(cexpr))
    print()
    print(expr)
    print()
    cexpr = contract_simplify_compile(expr, is_isospin_symmetric_limit=True)
    print(display_cexpr(cexpr))
    print()
    cexpr.optimize()
    print(display_cexpr(cexpr))
    print(cexpr_code_gen_py(cexpr))
    print()
    print("mk_test_expr_wick")
    print()
    expr_list = mk_test_expr_wick_07()
    with q.TimerFork():
        cexpr = contract_simplify_compile(*expr_list, is_isospin_symmetric_limit=True)
        cexpr.optimize()
    print(display_cexpr(cexpr))
    print()
    is_cython = False
    base_positions_dict = dict()
    print(cexpr_code_gen_py(cexpr, is_cython=is_cython))
    print()
