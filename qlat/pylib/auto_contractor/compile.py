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

try:
    from .wick import *
except:
    from wick import *

from itertools import permutations

import qlat as q

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

### ----

def get_var_name_type(x):
    # types include: V_S (wilson matrix), V_G (spin matrix), V_Tr (AMA c-number), V_a (c-number)
    if x.startswith("V_S_"):
        return "V_S"
    elif x.startswith("V_prod_GG_"):
        return "V_G"
    elif x.startswith("V_prod_GS_"):
        return "V_S"
    elif x.startswith("V_prod_SG_"):
        return "V_S"
    elif x.startswith("V_prod_SS_"):
        return "V_S"
    elif x.startswith("V_tr_"):
        return "V_Tr"
    else:
        assert False

def get_op_type(x):
    # x should be an Op
    # return a string which indicate the type of the op
    if not isinstance(x, Op):
        return None
    if x.otype == "S":
        return "S"
    elif x.otype == "G":
        return "G"
    elif x.otype == "Tr":
        return "Tr"
    elif x.otype == "Var":
        return get_var_name_type(x.name)
    else:
        assert False
    return None

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

def get_positions(term):
    s = set()
    add_positions(s, term)
    return sorted(list(s))

def collect_position_in_cexpr(named_terms):
    s = set()
    for name, term in named_terms:
        add_positions(s, term)
    positions = sorted(list(s))
    return positions

def collect_factor_in_cexpr(named_exprs):
    # collect the factors in all ea_coef
    variables_factor = []
    var_counter = 0
    var_dataset = {} # var_dataset[factor_code] = factor_var
    var_nameset = set()
    for name, expr in named_exprs:
        for i, (ea_coef, term_name,) in enumerate(expr):
            expr[i] = (ea.simplified(ea_coef), term_name,)
    for name, expr in named_exprs:
        for ea_coef, term_name in expr:
            if not isinstance(ea_coef, ea.Expr):
                continue
            for t in ea_coef.terms:
                x = t.factors
                for i, f in enumerate(x):
                    if f.otype != "Var":
                        assert f.otype == "Expr"
                        if f.code in var_dataset:
                            x[i] = var_dataset[f.code]
                        else:
                            while True:
                                name = f"V_factor_{var_counter}"
                                var_counter += 1
                                if name not in var_nameset:
                                    break
                            var_nameset.add(name)
                            variables_factor.append((name, f))
                            var = ea.Factor(name)
                            x[i] = var
                            var_dataset[f.code] = var
    return variables_factor

def collect_prop_in_cexpr(named_terms):
    # collect the propagators
    # modify the named_terms in-place and return the prop variable definitions as variables
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
        elif isinstance(x, Term):
            add_prop_variables(x.c_ops)
        elif isinstance(x, Expr):
            for t in x.terms:
                add_prop_variables(t)
    for name, term in named_terms:
        add_prop_variables(term)
        term.sort()
    return variables_prop

def collect_tr_in_cexpr(named_terms):
    # collect common traces
    # modify named_terms in-place and return definitions as variables_tr
    # variables_tr = [ (name, value,), ... ]
    # possible (name, value,) includes
    # ("V_tr_0", op,) where op.otype == "Tr"
    var_nameset = set()
    variables_tr = []
    var_counter = 0
    var_dataset_tr = {} # var_dataset[op_repr] = op_var
    def add_tr_varibles(x):
        nonlocal var_counter
        if isinstance(x, Term):
            add_tr_varibles(x.c_ops)
        elif isinstance(x, Op) and x.otype == "Tr":
            add_tr_varibles(x.ops)
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

def find_common_subexpr_in_tr(variables_tr):
    # return None or [ op, op1, ]
    subexpr_count = {}
    def add(x, count_added):
        op_repr = repr(x)
        if op_repr in subexpr_count:
            c, op = subexpr_count[op_repr]
            assert x == op
            subexpr_count[op_repr] = (c + count_added, x)
        else:
            subexpr_count[op_repr] = (count_added, x)
    def find_op_pair(x):
        if len(x) < 2:
            return None
        for i, op in enumerate(x):
            factor = 1
            if len(x) > 2:
                factor = 2
            op_type = get_op_type(op)
            if op_type in [ "V_S", "V_G", "S", "G", ]:
                op1 = x[(i+1) % len(x)]
                op1_type = get_op_type(op1)
                if op1_type in [ "V_S", "V_G", "S", "G", ]:
                    prod = [op, op1]
                    if op_type in [ "V_G", "G", ] and op1_type in [ "V_G", "G", ]:
                        add(prod, 1.02 * factor)
                    elif op_type in [ "V_G", "G", ] or op1_type in [ "V_G", "G", ]:
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
            find_op_pair(x.ops)
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

def collect_common_subexpr_in_tr(variables_tr, op_common, var):
    op_repr = repr(op_common)
    def replace(x):
        if x is None:
            return None
        elif isinstance(x, list):
            # need to represent the product of the list of operators
            for op in x:
                replace(op)
            if len(x) < 2:
                return None
            for i, op in enumerate(x):
                if isinstance(op, Op) and op.otype in ["Var", "S", "G",]:
                    i1 = (i+1) % len(x)
                    op1 = x[i1]
                    if isinstance(op1, Op) and op1.otype in ["Var", "S", "G",]:
                        prod = [op, op1]
                        if repr(prod) == op_repr:
                            x[i1] = None
                            x[i] = var
        elif isinstance(x, Op) and x.otype == "Tr" and len(x.ops) >= 2:
            replace(x.ops)
        elif isinstance(x, Term):
            replace(x.c_ops)
        elif isinstance(x, Expr):
            for t in x.terms:
                replace(t)
    def remove_none(x):
        # return a None removed x
        # possibly modify in-place
        if isinstance(x, list):
            return [ remove_none(op) for op in x if op is not None ]
        elif isinstance(x, Op):
            if x.otype == "Tr":
                x.ops = remove_none(x.ops)
            return x
        elif isinstance(x, Term):
            x.c_ops = remove_none(x.c_ops)
            x.a_ops = remove_none(x.a_ops)
            return x
        elif isinstance(x, Expr):
            x.terms = [ remove_none(t) for t in x.terms]
            return x
        else:
            assert False
    for name, tr in variables_tr:
        replace(tr)
        remove_none(tr)

def collect_subexpr_in_cexpr(variables_tr):
    # collect common sub-expressions
    # modify variables_tr in-place and return definitions as variables_prod
    # variables_prod = [ (name, value,), ... ]
    # possible (name, value,) includes
    # ("V_prod_SG_0", [ op, op1, ],) where get_op_type(op) in [ "V_S", "S", ] and get_op_type(op1) in [ "V_G", "G", ]
    var_nameset = set()
    var_counter_dict = {}
    var_counter_dict["V_prod_GG_"] = 0
    var_counter_dict["V_prod_GS_"] = 0
    var_counter_dict["V_prod_SG_"] = 0
    var_counter_dict["V_prod_SS_"] = 0
    variables_prod = []
    while True:
        subexpr = find_common_subexpr_in_tr(variables_tr)
        if subexpr is None:
            break
        [ op, op1, ] = subexpr
        op_type = get_op_type(op)
        op1_type = get_op_type(op1)
        assert op_type in [ "V_S", "V_G", "S", "G", ]
        assert op1_type in [ "V_S", "V_G", "S", "G", ]
        if op_type in [ "V_G", "G", ] and op1_type in [ "V_G", "G", ]:
            name_prefix = "V_prod_GG_"
        elif op_type in [ "V_G", "G", ] and op1_type in [ "V_S", "S", ]:
            name_prefix = "V_prod_GS_"
        elif op_type in [ "V_S", "S", ] and op1_type in [ "V_G", "G", ]:
            name_prefix = "V_prod_SG_"
        elif op_type in [ "V_S", "S", ] and op1_type in [ "V_S", "S", ]:
            name_prefix = "V_prod_SS_"
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

    # self.diagram_types
    # self.positions
    # self.variables_prop
    # self.variables_factor
    # self.variables_prod
    # self.variables_tr
    # self.named_terms
    # self.named_exprs
    # self.function
    #
    # self.named_terms[i] = (term_name, Term(c_ops, [], 1),)
    # self.named_exprs[i] = (expr_name, [ (ea_coef, term_name,), ... ],)
    # self.positions == sorted(list(self.positions))

    def __init__(self):
        self.diagram_types = []
        self.positions = []
        self.variables_factor = []
        self.variables_prop = []
        self.variables_prod = []
        self.variables_tr = []
        self.named_terms = []
        self.named_exprs = []
        self.function = None

    @q.timer
    def optimize(self):
        # interface function
        self.collect_op()

    @q.timer
    def collect_op(self):
        # interface function
        # Performing common sub-expression elimination
        # Should be called after contract_simplify_compile(*exprs) or mk_cexpr(*exprs)
        # The cexpr cannot be evaluated before collect_op!!!
        # eval term factor
        for name, term in self.named_terms:
            assert term.coef == 1
            assert term.a_ops == []
        checks = [
                self.variables_factor == [],
                self.variables_prop == [],
                self.variables_tr == [],
                self.variables_prod == [],
                ]
        if not all(checks):
            # likely collect_op is already performed
            return
        # collect ea_coef factors into variables
        self.variables_factor = collect_factor_in_cexpr(self.named_exprs)
        # collect prop expr into variables
        self.variables_prop = collect_prop_in_cexpr(self.named_terms)
        # collect trace expr into variables
        self.variables_tr = collect_tr_in_cexpr(self.named_terms)
        # collect common prod into variables
        self.variables_prod = collect_subexpr_in_cexpr(self.variables_tr)

### ----

def inc(type_dict, key):
    if key in type_dict:
        type_dict[key] += 1
    else:
        type_dict[key] = 1

def drop_tag_last_subscript(tag):
    return tag.rsplit("_", 1)[0]

def mk_permuting_dicts(pos_list):
    for p_pos_list in permutations(pos_list):
        yield dict(zip(pos_list, p_pos_list))

def get_position_permutation_groups(term):
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
            inc(type_dict, (op.p1, op.p2,))
        elif op.otype == "Tr":
            loop_term_ops(type_dict, op.ops)

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

def filter_diagram_type(expr, diagram_type_dict = None, included_types = None):
    """
    drop diagrams with diagram_type_dict[diagram_type] == None
    if included_types is not None:
        only keep diagrams with diagram_type that diagram_type in included_types.
    """
    if diagram_type_dict is None:
        return expr
    new_terms = []
    for term in expr.terms:
        diagram_type = get_term_diagram_type_info(term)
        if diagram_type in diagram_type_dict:
            if diagram_type_dict[diagram_type] is None:
                continue
            if included_types is not None:
                diagram_type_name = diagram_type_dict.get(diagram_type)
                if diagram_type_name not in included_types:
                    continue
        new_terms.append(term)
    included_types_tag = ""
    if included_types is not None:
        included_types_tag = " (" + ','.join(included_types) + ")"
    return Expr(new_terms, expr.description + included_types_tag)

def mk_cexpr(*exprs, diagram_type_dict = None):
    """
    exprs already finished wick contraction,
    otherwise use contract_simplify_compile(*exprs, is_isospin_symmetric_limit, diagram_type_dict)
    !!!if diagram_type_dict[diagram_type] == None: this diagram_type should have already be dropped!!!
    interface function
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
                diagram_type_name = f"ADT{diagram_type_counter:0>2}" # ADT is short for "auto diagram type"
                diagram_type_counter += 1
                diagram_type_dict[diagram_type] = diagram_type_name
            diagram_type_name = diagram_type_dict[diagram_type]
            diagram_type_term_dict[repr_term] = diagram_type_name
            assert diagram_type_name is not None
            term_name_counter = 0
            while True:
                term_name = f"term_{diagram_type_name}_{term_name_counter:0>4}"
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
        typed_expr_list_dict = { name : [] for name, diagram_type in diagram_types }
        for j, term in enumerate(expr.terms):
            coef = term.coef
            term.coef = 1
            repr_term = repr(term)
            diagram_type_name = diagram_type_term_dict[repr_term]
            assert diagram_type_name is not None
            term_name = term_dict[repr_term]
            typed_expr_list_dict[diagram_type_name].append((coef, term_name,))
            expr_list.append((coef, term_name,))
        named_exprs.append((f"# {descriptions[i]}\nexprs[{i}]", expr_list,))
    # positions
    positions = collect_position_in_cexpr(named_terms)
    # cexpr
    cexpr = CExpr()
    cexpr.diagram_types = diagram_types
    cexpr.positions = positions
    cexpr.named_terms = named_terms
    cexpr.named_exprs = named_exprs
    return cexpr

@q.timer
def contract_simplify(*exprs, is_isospin_symmetric_limit = True, diagram_type_dict = None):
    """
    exprs = [ expr, (expr, *included_types,), ... ]\n
    In case diagram_type_dict is not None, perform the following filter
    If diagram_type_dict[diagram_type] is None: term is removed.
    For (expr, *included_types,), only terms (in expr) with diagram_type that is in included_types is kept
    included_types should be a list/tuple of string
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
        expr.simplify(is_isospin_symmetric_limit = is_isospin_symmetric_limit)
        expr = filter_diagram_type(expr,
                diagram_type_dict = diagram_type_dict,
                included_types = included_types)
        return expr
    return q.parallel_map(func, exprs)

@q.timer
def compile_expr(*exprs, diagram_type_dict = None):
    """
    interface function
    """
    exprs = copy.deepcopy(exprs)
    cexpr = mk_cexpr(*exprs, diagram_type_dict = diagram_type_dict)
    return cexpr

@q.timer
def contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True, diagram_type_dict = None):
    """
    Call ``contract_simplify`` and then ``compile_expr``\n
    This function can be used to construct the first argument of ``cached_comipled_cexpr``.
    e.g. exprs = [ Qb("u", "x", s, c) * Qv("u", "x", s, c) + "u_bar*u", Qb("s", "x", s, c) * Qv("s", "x", s, c) + "s_bar*s", Qb("c", "x", s, c) * Qv("c", "x", s, c) + "c_bar*c", ]
    e.g. exprs = [ mk_pi_p("x2", True) * mk_pi_p("x1") + "(pi   * pi)", mk_j5pi_mu("x2", 3) * mk_pi_p("x1") + "(a_pi * pi)", mk_k_p("x2", True)  * mk_k_p("x1")  + "(k    * k )", mk_j5k_mu("x2", 3)  * mk_k_p("x1")  + "(a_k  * k )", ]
    """
    contracted_simplified_exprs = contract_simplify(
            *exprs,
            is_isospin_symmetric_limit = is_isospin_symmetric_limit,
            diagram_type_dict = diagram_type_dict)
    cexpr = compile_expr(*contracted_simplified_exprs, diagram_type_dict = diagram_type_dict)
    return cexpr

def show_variable_value(value):
    if isinstance(value, list):
        return "*".join(map(show_variable_value, value))
    elif isinstance(value, Var):
        return f"{value.name}"
    elif isinstance(value, G) and value.tag in [0, 1, 2, 3, 5]:
        tag = { 0: "x", 1: "y", 2: "z", 3: "t", 5: "5", }[value.tag]
        return f"gamma_{tag}"
    elif isinstance(value, G):
        return f"gamma({value.tag})"
    elif isinstance(value, S):
        return f"S_{value.f}({value.p1},{value.p2})"
    elif isinstance(value, Tr):
        expr = "*".join(map(show_variable_value, value.ops))
        return f"tr({expr})"
    elif isinstance(value, Term):
        if value.coef == 1:
            return "*".join(map(show_variable_value, value.c_ops + value.a_ops))
        else:
            return "*".join(map(show_variable_value, [ f"({value.coef})", ] + value.c_ops + value.a_ops))
    elif isinstance(value, tuple) and len(value) == 2 and isinstance(value[1], str):
        return f"({value[0]})*{value[1]}"
    elif isinstance(value, ea.Factor):
        return f"({value.code})"
    else:
        return f"{value}"

def display_cexpr_raw(cexpr : CExpr):
    # return a string
    # interface function
    lines = []
    lines.append(f"Begin CExpr")
    if cexpr.diagram_types:
        lines.append(f"diagram_type_dict = dict()")
        for name, diagram_type in cexpr.diagram_types:
            lines.append(f"diagram_type_dict[{diagram_type}] = {name!r}")
    if cexpr.positions:
        lines.append(f"Positions: {cexpr.positions}")
    if cexpr.variables_prop:
        lines.append(f"Variables prop:")
    for name, value in cexpr.variables_prop:
        lines.append(f"{name:>20} : {value}")
    if cexpr.variables_factor:
        lines.append(f"Variables factor:")
    for name, value in cexpr.variables_factor:
        lines.append(f"{name:>20} : {value}")
    if cexpr.variables_prod:
        lines.append(f"Variables prod:")
    for name, value in cexpr.variables_prod:
        lines.append(f"{name:>20} : {value}")
    if cexpr.variables_tr:
        lines.append(f"Variables tr:")
    for name, value in cexpr.variables_tr:
        lines.append(f"{name:>20} : {value}")
    lines.append(f"Named terms:")
    for name, term in cexpr.named_terms:
        lines.append(f"{name:>20} : {term}")
    lines.append(f"Named exprs:")
    for name, expr in cexpr.named_exprs:
        lines.append(f"{name} :\n  {expr}")
    lines.append(f"End CExpr")
    return "\n".join(lines)

def display_cexpr(cexpr : CExpr):
    # return a string
    # interface function
    lines = []
    lines.append(f"Begin CExpr")
    if cexpr.diagram_types:
        lines.append(f"diagram_type_dict = dict()")
        for name, diagram_type in cexpr.diagram_types:
            lines.append(f"diagram_type_dict[{diagram_type}] = {name!r}")
    if cexpr.positions:
        position_vars = ", ".join(cexpr.positions)
        lines.append(f"Positions:\n{position_vars} = {cexpr.positions}")
    if cexpr.variables_prop:
        lines.append(f"Variables prop:")
    for name, value in cexpr.variables_prop:
        lines.append(f"{name:>20} = {show_variable_value(value)}")
    if cexpr.variables_factor:
        lines.append(f"Variables factor:")
    for name, value in cexpr.variables_factor:
        lines.append(f"{name:>20} = {show_variable_value(value)}")
    if cexpr.variables_prod:
        lines.append(f"Variables expr:")
    for name, value in cexpr.variables_prod:
        lines.append(f"{name:>20} = {show_variable_value(value)}")
    if cexpr.variables_tr:
        lines.append(f"Variables tr:")
    for name, value in cexpr.variables_tr:
        lines.append(f"{name:>20} = {show_variable_value(value)}")
    lines.append(f"terms = [")
    for name, term in cexpr.named_terms:
        lines.append(f"  {show_variable_value(term)}, # {name}")
    lines.append(f"]")
    for name, diagram_type in cexpr.diagram_types:
        if name is not None:
            lines.append(f"coef_{name} = 1")
    for idx, (name, term) in enumerate(cexpr.named_terms):
        name_type = "_".join([ "coef", ] + name.split("_")[1:-1])
        lines.append(f"{name} = {name_type} * terms[{idx}]")
    lines.append(f"exprs = [ None for i in range({len(cexpr.named_exprs)}) ]")
    for name, expr, in cexpr.named_exprs:
        s = "+".join(map(show_variable_value, expr))
        if s == "":
            s = 0
        lines.append(f"{name} = {s}")
    lines.append(f"End CExpr")
    return "\n".join(lines)

@q.timer_verbose
def cexpr_code_gen_py(cexpr : CExpr):
    # return a string
    # interface function
    return CExprCodeGenPy(cexpr).code_gen()

class CExprCodeGenPy:

    # self.cexpr
    # self.lines
    # self.indent
    # self.total_sloppy_flops

    # flops per complex multiplication: 6
    # flops per matrix multiplication: 6 M N L + 2 M L (N-1) ==> 13536 (sc * sc), 4320 (sc * s), 480 (s * s)
    # flops per trace 2 (M-1) ==> 22
    # flops per trace2 6 M N + 2 (M N - 1) ==> 1150

    def __init__(self, cexpr):
        self.cexpr = cexpr
        self.lines = []
        self.indent = 0
        self.total_sloppy_flops = 0

    def code_gen(self):
        lines = self.lines
        lines.append(f"from auto_contractor.runtime import *")
        lines.append(f"cimport qlat_utils.everything as cc")
        lines.append(f"cimport qlat_utils.cp as cp")
        lines.append(f"cimport numpy as np")
        self.sep()
        self.cexpr_function()
        self.sep()
        self.cexpr_function_get_prop()
        self.sep()
        self.cexpr_function_eval()
        self.sep()
        self.cexpr_function_get_factor()
        self.sep()
        self.cexpr_function_eval_with_props()
        self.sep()
        self.total_flops()
        return "\n".join(lines)

    def gen_expr(self, x):
        # return code_str, type_str
        if isinstance(x, (int, float, complex)):
            return f"{x}", "V_a"
        elif isinstance(x, (ea.Expr, ea.Factor)):
            return f"({ea.compile_py(x)})", "V_a"
        assert isinstance(x, Op)
        if x.otype == "S":
            return f"get_prop('{x.f}', {x.p1}, {x.p2})", "V_S"
        elif x.otype == "G":
            assert x.s1 == "auto" and x.s2 == "auto"
            assert x.tag in [0, 1, 2, 3, 5]
            return f"cp.gamma_matrix_{x.tag}", "V_G"
        elif x.otype == "Tr":
            if len(x.ops) == 0:
                assert False
            elif len(x.ops) == 1:
                c, t = self.gen_expr(x.ops[0])
                assert t == "V_S"
                self.total_sloppy_flops += 22
                return f"cc.mat_tr({c})", "V_Tr"
            else:
                c1, t1 = self.gen_expr_prod_list(x.ops[:-1])
                c2, t2 = self.gen_expr(x.ops[-1])
                if t1 == "V_S" and t2 == "V_S":
                    self.total_sloppy_flops += 1150
                    return f"cc.mat_tr({c1}, {c2})", "V_Tr"
                elif t1 == "V_S" and t2 == "V_G":
                    return f"cc.mat_tr({c1}, {c2})", "V_Tr"
                elif t1 == "V_G" and t2 == "V_S":
                    return f"cc.mat_tr({c1}, {c2})", "V_Tr"
                else:
                    assert False
        elif x.otype == "Var":
            if x.name.startswith("V_S_"):
                return f"p_{x.name}[0]", get_var_name_type(x.name)
            else:
                return f"{x.name}", get_var_name_type(x.name)

    def gen_expr_prod(self, ct1, ct2):
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
            return f"{c1} * {c2}", "V_S"
        elif t1 == "V_S" and t2 == "V_G":
            self.total_sloppy_flops += 4320
            return f"{c1} * {c2}", "V_S"
        elif t1 == "V_G" and t2 == "V_S":
            self.total_sloppy_flops += 4320
            return f"{c1} * {c2}", "V_S"
        elif t1 == "V_G" and t2 == "V_G":
            self.total_sloppy_flops += 480
            return f"{c1} * {c2}", "V_G"
        elif t1 == "V_G" and t2 == "V_a":
            return f"{c1} * {c2}", "V_G"
        elif t1 == "V_a" and t2 == "V_G":
            return f"{c1} * {c2}", "V_G"
        elif t1 == "V_S" and t2 == "V_a":
            return f"{c1} * {c2}", "V_S"
        elif t1 == "V_a" and t2 == "V_S":
            return f"{c1} * {c2}", "V_S"
        elif t1 == "V_a" and t2 == "V_a":
            return f"{c1} * {c2}", "V_a"
        elif t1 == "V_Tr" and t2 == "V_Tr":
            return f"{c1} * {c2}", "V_Tr"
        elif t1 == "V_a" and t2 == "V_Tr":
            return f"{c2} * {c1}", "V_Tr"
        elif t1 == "V_Tr" and t2 == "V_a":
            return f"{c1} * {c2}", "V_Tr"
        else:
            print(ct1, ct2)
            assert False

    def gen_expr_prod_list(self, x_list):
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

    def sep(self):
        append = self.append
        append(f"")
        append(f"### ----")
        append(f"")

    def cexpr_function(self):
        append = self.append
        append(f"@timer")
        append(f"def cexpr_function(*, positions_dict, get_prop, is_ama_and_sloppy = False):")
        self.indent += 4
        append(f"# get_props")
        append(f"props = cexpr_function_get_prop(positions_dict, get_prop)")
        append(f"# eval")
        append(f"ama_val = cexpr_function_eval(positions_dict, props)")
        append(f"# extract sloppy val")
        append(f"val_sloppy = ama_extract(ama_val, is_sloppy = True)")
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
        cexpr = self.cexpr
        append(f"@timer_flops")
        append(f"def cexpr_function_get_prop(positions_dict, get_prop):")
        self.indent += 4
        append(f"# set positions")
        for position_var in cexpr.positions:
            append(f"{position_var} = positions_dict['{position_var}']")
        append(f"# get prop")
        for name, value in cexpr.variables_prop:
            assert name.startswith("V_S_")
            x = value
            assert isinstance(x, Op)
            assert x.otype == "S"
            c, t = self.gen_expr(x)
            assert t == "V_S"
            append(f"{name} = {c}")
        append(f"# prop")
        append(f"# set props for return")
        append(f"props = [")
        self.indent += 4
        for name, value in cexpr.variables_prop:
            append(f"{name},")
        append(f"]")
        self.indent -= 4
        append(f"# set flops")
        append(f"total_flops = len(props) * 144 * 2 * 8")
        append(f"# return")
        append(f"return total_flops, props")
        self.indent -= 4

    def cexpr_function_get_factor(self):
        append = self.append
        cexpr = self.cexpr
        append(f"@timer")
        append(f"def cexpr_function_get_factor(positions_dict):")
        self.indent += 4
        append(f"# set positions")
        append(f"size = positions_dict.get('size')")
        for position_var in cexpr.positions:
            append(f"{position_var}_type, {position_var} = positions_dict['{position_var}']")
        append(f"# declare factors")
        append(f"cdef np.ndarray[np.complex128_t] factors")
        append(f"factors = numpy.zeros({len(cexpr.variables_factor)}, dtype = numpy.complex128)")
        append(f"cdef cc.Complex[:] factors_view = factors")
        append(f"# set factors")
        for idx, (name, value,) in enumerate(cexpr.variables_factor):
            assert name.startswith("V_factor_")
            x = value
            assert isinstance(x, ea.Factor)
            assert x.otype == "Expr"
            c, t = self.gen_expr(x)
            assert t == "V_a"
            append(f"# {name}")
            append(f"factors_view[{idx}] = {c}")
        append(f"# return factors")
        append(f"return factors")
        self.indent -= 4

    def cexpr_function_eval(self):
        append = self.append
        cexpr = self.cexpr
        append(f"@timer_flops")
        append(f"def cexpr_function_eval(positions_dict, props):")
        self.indent += 4
        append(f"# load AMA props with proper format")
        append(f"props = [ load_prop(p) for p in props ]")
        append(f"# join the AMA props")
        append(f"ama_props = ama_list(*props)")
        append(f"# get_factor")
        append(f"factors = cexpr_function_get_factor(positions_dict)")
        append(f"# apply eval to the factors and AMA props")
        append(f"ama_val = ama_apply2_r(cexpr_function_eval_with_props, factors, ama_props)")
        append(f"# set flops")
        append(f"total_flops = ama_counts(ama_val) * total_sloppy_flops")
        append(f"# return")
        append(f"return total_flops, ama_val")
        self.indent -= 4

    def cexpr_function_eval_with_props(self):
        append = self.append
        cexpr = self.cexpr
        append(f"@timer_flops")
        append(f"def cexpr_function_eval_with_props(cc.Complex[:] factors, list props):")
        self.indent += 4
        append(f"# set factors")
        for idx, (name, value,) in enumerate(cexpr.variables_factor):
            append(f"cdef cc.Complex {name} = factors[{idx}]")
        append(f"# set props")
        for idx, (name, value,) in enumerate(cexpr.variables_prop):
            append(f"cdef cc.WilsonMatrix* p_{name} = &(<cp.WilsonMatrix>props[{idx}]).xx")
        append(f"# compute products")
        for name, value in cexpr.variables_prod:
            assert name.startswith("V_prod_")
            x = value
            assert isinstance(x, list)
            c, t = self.gen_expr_prod_list(x)
            assert t == get_var_name_type(name)
            if t == "V_G":
                append(f"cdef cc.SpinMatrix {name} = {c}")
            elif t == "V_S":
                append(f"cdef cc.WilsonMatrix {name} = {c}")
            else:
                assert False
        append(f"# compute traces")
        for name, value in cexpr.variables_tr:
            assert name.startswith("V_tr_")
            x = value
            assert isinstance(x, Op)
            assert x.otype == "Tr"
            c, t = self.gen_expr(x)
            assert t == "V_Tr"
            append(f"cdef cc.Complex {name} = {c}")
        append(f"# set terms")
        for name, term in cexpr.named_terms:
            x = term
            if x.coef == 1:
                c_ops = x.c_ops
            else:
                c_ops = [ x.coef, ] + x.c_ops
            c, t = self.gen_expr_prod_list(c_ops)
            append(f"cdef cc.Complex {name} = {c}")
        append(f"# declare exprs")
        append(f"cdef np.ndarray[np.complex128_t] exprs")
        append(f"exprs = numpy.zeros({len(cexpr.named_exprs)}, dtype = numpy.complex128)")
        append(f"cdef cc.Complex[:] exprs_view = exprs")
        append(f"# set exprs")
        def show_coef_term(coef, tname):
            coef = ea.compile_py(coef)
            if coef == "1":
                return f"{tname}"
            else:
                return f"({coef}) * {tname}"
        append(f"cdef cc.Complex expr")
        for idx, (name, expr,) in enumerate(cexpr.named_exprs):
            name = name.replace("\n", "  ")
            append(f"# {idx} name='{name}' ")
            append(f"expr = 0")
            for coef, tname in expr:
                s = show_coef_term(coef, tname)
                append(f"expr += {s}")
            append(f"exprs_view[{idx}] = expr")
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

if __name__ == "__main__":
    expr = Qb("d", "x1", "s1", "c1") * G(5, "s1", "s2") * Qv("u", "x1", "s2", "c1") * Qb("u", "x2", "s3", "c2") * G(5, "s3", "s4") * Qv("d", "x2", "s4", "c2")
    print(expr)
    print()
    expr = simplified(contract_expr(expr))
    print(expr)
    print()
    cexpr = copy.deepcopy(mk_cexpr(expr))
    print(cexpr)
    print()
    cexpr.optimize()
    print(cexpr)
    print()
    print(display_cexpr(cexpr))
    print()
    print(expr)
    print()
    cexpr = contract_simplify_compile(expr, is_isospin_symmetric_limit = True)
    print(display_cexpr(cexpr))
    print()
    cexpr.optimize()
    print(display_cexpr(cexpr))
    print(cexpr_code_gen_py(cexpr))
