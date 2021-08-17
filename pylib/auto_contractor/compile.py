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

def collect_op_in_cexpr(variables, named_terms):
    var_counter = 0
    var_dataset = {} # var_dataset[op_repr] = op_var
    var_nameset = set()
    for name, value in variables:
        var_nameset.add(name)
    def add_prop_variables(x):
        nonlocal var_counter
        if isinstance(x, list):
            for i, op in enumerate(x):
                if op.otype in ["S",]:
                    op_repr = repr(op)
                    if op_repr in var_dataset:
                        x[i] = var_dataset[op_repr]
                    else:
                        while True:
                            var_counter += 1
                            name = f"V0_{var_counter}"
                            if name not in var_nameset:
                                break
                        variables.append((name, op,))
                        var = Var(name)
                        x[i] = var
                        var_dataset[op_repr] = var
                        var_nameset.add(name)
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

def find_common_subexpr_in_tr(variables_trs):
    subexpr_count = {}
    def add(x):
        op_repr = repr(x)
        if op_repr in subexpr_count:
            c, op = subexpr_count[op_repr]
            assert x == op
            subexpr_count[op_repr] = (c + 1, x)
        else:
            subexpr_count[op_repr] = (1, x)
    def find(x):
        if isinstance(x, list):
            # need to represent the product of the list of operators
            for op in x:
                find(op)
            if len(x) < 2:
                return None
            for i, op in enumerate(x):
                if isinstance(op, Op):
                    op1 = x[(i+1) % len(x)]
                    if op.otype in ["Var", "S",]:
                        if isinstance(op1, Op) and op1.otype in ["Var", "S", "G",]:
                            prod = [op, op1]
                            add(prod)
                    elif op.otype in ["G",]:
                        if isinstance(op1, Op) and op1.otype in ["Var", "S",]:
                            prod = [op, op1]
                            add(prod)
        elif isinstance(x, Op) and x.otype == "Tr" and len(x.ops) >= 2:
            find(x.ops)
        elif isinstance(x, Term):
            find(x.c_ops)
        elif isinstance(x, Expr):
            for t in x.terms:
                find(t)
    for name, tr in variables_trs:
        find(tr)
    max_num_repeat = 1
    best_match = None
    for num_repeat, op in subexpr_count.values():
        if num_repeat > max_num_repeat:
            max_num_repeat = num_repeat
            best_match = op
    return best_match

def collect_common_subexpr_in_tr(variables_trs, op_common, var):
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
    for name, tr in variables_trs:
        replace(tr)
        remove_none(tr)

def collect_subexpr_in_cexpr(variables, named_terms):
    var_nameset = set()
    for name, value in variables:
        var_nameset.add(name)
    variables_trs = []
    var_counter_tr = 0
    var_dataset_tr = {} # var_dataset[op_repr] = op_var
    def add_tr_varibles(x):
        nonlocal var_counter_tr
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
                            var_counter_tr += 1
                            name = f"V2_{var_counter_tr}"
                            if name not in var_nameset:
                                break
                        variables_trs.append((name, op,))
                        var = Var(name)
                        x[i] = var
                        var_dataset_tr[op_repr] = var
                        var_nameset.add(name)
    for name, term in named_terms:
        add_tr_varibles(term)
        term.sort()
    var_counter = 0
    while True:
        op = find_common_subexpr_in_tr(variables_trs)
        if op is None:
            break
        while True:
            var_counter += 1
            name = f"V1_{var_counter}"
            if name not in var_nameset:
                break
        variables.append((name, op,))
        var = Var(name)
        collect_common_subexpr_in_tr(variables_trs, op, var)
    variables += variables_trs

class CExpr:

    def __init__(self, diagram_types, variables, named_terms, named_typed_exprs, named_exprs, positions = None):
        self.diagram_types = diagram_types
        self.variables = variables
        self.named_terms = named_terms
        # typed_expr and expr are a collection of term names representing the sum of these terms
        self.named_typed_exprs = named_typed_exprs
        self.named_exprs = named_exprs
        if positions is not None:
            self.positions = positions
        else:
            s = set()
            for name, term in named_terms:
                add_positions(s, term)
            self.positions = sorted(list(s))

    def __repr__(self) -> str:
        return f"CExpr({self.diagram_types},{self.variables},{self.named_terms},{self.named_typed_exprs},{self.named_exprs},{self.positions})"

    def collect_op(self):
        # interface function
        collect_op_in_cexpr(self.variables, self.named_terms)
        collect_subexpr_in_cexpr(self.variables, self.named_terms)

def inc(type_dict, key):
    if key in type_dict:
        type_dict[key] += 1
    else:
        type_dict[key] = 1

def drop_tag_last_subscript(tag):
    return tag.rsplit("_", 1)[0]

def loop_term_ops(type_dict, ops):
    for op in ops:
        if op.otype == "S":
            # ADJUST DIAGRAM TYPE DEFINITION
            # inc(type_dict, (op.f, op.p1, op.p2,))
            # inc(type_dict, (op.p1, op.p2,))
            inc(type_dict, (drop_tag_last_subscript(op.p1), drop_tag_last_subscript(op.p2),))
            #
        elif op.otype == "Tr":
            loop_term_ops(type_dict, op.ops)

def get_term_diagram_type_info(term):
    type_dict = dict()
    loop_term_ops(type_dict, term.c_ops)
    return tuple(sorted(type_dict.items()))

def mk_cexpr(*exprs):
    # interface function
    diagram_type_dict = dict()
    diagram_type_counter = 0
    for i, expr in enumerate(exprs):
        for j, term in enumerate(expr.terms):
            diagram_type = get_term_diagram_type_info(term)
            if diagram_type not in diagram_type_dict:
                diagram_type_counter += 1
                diagram_type_name = f"ADT{diagram_type_counter}" # ADT is short for "auto diagram type"
                diagram_type_dict[diagram_type] = diagram_type_name
    diagram_types = []
    for diagram_type, diagram_type_name in diagram_type_dict.items():
        diagram_types.append((diagram_type_name, diagram_type,))
    named_terms = []
    named_typed_exprs = []
    named_exprs = []
    for i, expr in enumerate(exprs):
        expr_list = []
        typed_expr_list_dict = { name : [] for name, diagram_type in diagram_types }
        for j, term in enumerate(expr.terms):
            diagram_type = get_term_diagram_type_info(term)
            diagram_type_name = diagram_type_dict[diagram_type]
            name = f"T{i+1}_{j+1}_{diagram_type_name}"
            named_terms.append((name, term,))
            typed_expr_list_dict[diagram_type_name].append(name)
            expr_list.append(name)
        for diagram_type_name, typed_expr_list in typed_expr_list_dict.items():
            named_typed_exprs.append((f"E{i+1}_{diagram_type_name}", typed_expr_list,))
        named_exprs.append((f"E{i+1}", expr_list,))
    return CExpr(diagram_types, [], named_terms, named_typed_exprs, named_exprs)

def contract_simplify_round_compile(*exprs, is_isospin_symmetric_limit = True):
    # interface function
    exprs = list(exprs)
    for i in range(len(exprs)):
        expr = copy.deepcopy(exprs[i])
        expr = contract_expr(expr)
        expr.simplify(is_isospin_symmetric_limit = is_isospin_symmetric_limit)
        exprs[i] = expr.round()
    cexpr = mk_cexpr(*exprs)
    return cexpr

def display_cexpr(cexpr : CExpr):
    # interface function
    # return a string
    lines = []
    lines.append(f"Begin CExpr")
    lines.append(f"{'Positions':>10} : {cexpr.positions}")
    for name, diagram_type in cexpr.diagram_types:
        lines.append(f"{name:>10} : {diagram_type}")
    for name, value in cexpr.variables:
        lines.append(f"{name:>10} : {value}")
    for name, term in cexpr.named_terms:
        lines.append(f"{name:>20} : {term}")
    for name, typed_expr in cexpr.named_typed_exprs:
        lines.append(f"{name:>15} : {typed_expr}")
    for name, expr in cexpr.named_exprs:
        lines.append(f"{name:>10} : {expr}")
    lines.append(f"End CExpr")
    return "\n".join(lines)

if __name__ == "__main__":
    expr = Qb("d", "x1", "s1", "c1") * G(5, "s1", "s2") * Qv("u", "x1", "s2", "c1") * Qb("u", "x2", "s3", "c2") * G(5, "s3", "s4") * Qv("d", "x2", "s4", "c2")
    print(expr)
    expr = simplified(contract_expr(expr))
    print(expr.round())
    cexpr = mk_cexpr(expr.round())
    print(cexpr)
    cexpr.collect_op()
    print(cexpr)
    print(display_cexpr(cexpr))
    print(CExpr([('S_1', S('d','x2','x1')), ('S_2', S('u','x1','x2'))],[('T_1', Term([Tr([G(5), Var('S_1'), G(5), Var('S_2')],'sc')],[],(-1+0j)))],['x1', 'x2']))
    expr = Expr([Term([G(5,'a_s_203','a_s_204'), G(5,'a_s_207','a_s_208'), G(5,'a_s_219','a_s_220'), G(5,'a_s_223','a_s_224')],[Qb('u','x21','a_s_203','a_c_102'), Qv('u','x21','a_s_204','a_c_102'), Qb('u','x22','a_s_207','a_c_104'), Qv('u','x22','a_s_208','a_c_104'), Qb('u','x11','a_s_219','a_c_110'), Qv('u','x11','a_s_220','a_c_110'), Qb('u','x12','a_s_223','a_c_112'), Qv('u','x12','a_s_224','a_c_112')],(0.08333333333333333+0j)), Term([G(5,'a_s_203','a_s_204'), G(5,'a_s_207','a_s_208'), G(5,'a_s_219','a_s_220'), G(5,'a_s_225','a_s_226')],[Qb('u','x21','a_s_203','a_c_102'), Qv('u','x21','a_s_204','a_c_102'), Qb('u','x22','a_s_207','a_c_104'), Qv('u','x22','a_s_208','a_c_104'), Qb('u','x11','a_s_219','a_c_110'), Qv('u','x11','a_s_220','a_c_110'), Qb('d','x12','a_s_225','a_c_113'), Qv('d','x12','a_s_226','a_c_113')],(-0.08333333333333333+0j)), Term([G(5,'a_s_203','a_s_204'), G(5,'a_s_207','a_s_208'), G(5,'a_s_221','a_s_222'), G(5,'a_s_223','a_s_224')],[Qb('u','x21','a_s_203','a_c_102'), Qv('u','x21','a_s_204','a_c_102'), Qb('u','x22','a_s_207','a_c_104'), Qv('u','x22','a_s_208','a_c_104'), Qb('d','x11','a_s_221','a_c_111'), Qv('d','x11','a_s_222','a_c_111'), Qb('u','x12','a_s_223','a_c_112'), Qv('u','x12','a_s_224','a_c_112')],(-0.08333333333333333+0j)), Term([G(5,'a_s_203','a_s_204'), G(5,'a_s_207','a_s_208'), G(5,'a_s_221','a_s_222'), G(5,'a_s_225','a_s_226')],[Qb('u','x21','a_s_203','a_c_102'), Qv('u','x21','a_s_204','a_c_102'), Qb('u','x22','a_s_207','a_c_104'), Qv('u','x22','a_s_208','a_c_104'), Qb('d','x11','a_s_221','a_c_111'), Qv('d','x11','a_s_222','a_c_111'), Qb('d','x12','a_s_225','a_c_113'), Qv('d','x12','a_s_226','a_c_113')],(0.08333333333333333+0j)), Term([G(5,'a_s_203','a_s_204'), G(5,'a_s_207','a_s_208'), G(5,'a_s_227','a_s_228'), G(5,'a_s_229','a_s_230')],[Qb('u','x21','a_s_203','a_c_102'), Qv('u','x21','a_s_204','a_c_102'), Qb('u','x22','a_s_207','a_c_104'), Qv('u','x22','a_s_208','a_c_104'), Qb('d','x11','a_s_227','a_c_114'), Qv('u','x11','a_s_228','a_c_114'), Qb('u','x12','a_s_229','a_c_115'), Qv('d','x12','a_s_230','a_c_115')],(0.16666666666666669+0j)), Term([G(5,'a_s_203','a_s_204'), G(5,'a_s_207','a_s_208'), G(5,'a_s_231','a_s_232'), G(5,'a_s_233','a_s_234')],[Qb('u','x21','a_s_203','a_c_102'), Qv('u','x21','a_s_204','a_c_102'), Qb('u','x22','a_s_207','a_c_104'), Qv('u','x22','a_s_208','a_c_104'), Qb('u','x11','a_s_231','a_c_116'), Qv('d','x11','a_s_232','a_c_116'), Qb('d','x12','a_s_233','a_c_117'), Qv('u','x12','a_s_234','a_c_117')],(0.16666666666666669+0j)), Term([G(5,'a_s_203','a_s_204'), G(5,'a_s_209','a_s_210'), G(5,'a_s_219','a_s_220'), G(5,'a_s_223','a_s_224')],[Qb('u','x21','a_s_203','a_c_102'), Qv('u','x21','a_s_204','a_c_102'), Qb('d','x22','a_s_209','a_c_105'), Qv('d','x22','a_s_210','a_c_105'), Qb('u','x11','a_s_219','a_c_110'), Qv('u','x11','a_s_220','a_c_110'), Qb('u','x12','a_s_223','a_c_112'), Qv('u','x12','a_s_224','a_c_112')],(-0.08333333333333333+0j)), Term([G(5,'a_s_203','a_s_204'), G(5,'a_s_209','a_s_210'), G(5,'a_s_219','a_s_220'), G(5,'a_s_225','a_s_226')],[Qb('u','x21','a_s_203','a_c_102'), Qv('u','x21','a_s_204','a_c_102'), Qb('d','x22','a_s_209','a_c_105'), Qv('d','x22','a_s_210','a_c_105'), Qb('u','x11','a_s_219','a_c_110'), Qv('u','x11','a_s_220','a_c_110'), Qb('d','x12','a_s_225','a_c_113'), Qv('d','x12','a_s_226','a_c_113')],(0.08333333333333333-0j)), Term([G(5,'a_s_203','a_s_204'), G(5,'a_s_209','a_s_210'), G(5,'a_s_221','a_s_222'), G(5,'a_s_223','a_s_224')],[Qb('u','x21','a_s_203','a_c_102'), Qv('u','x21','a_s_204','a_c_102'), Qb('d','x22','a_s_209','a_c_105'), Qv('d','x22','a_s_210','a_c_105'), Qb('d','x11','a_s_221','a_c_111'), Qv('d','x11','a_s_222','a_c_111'), Qb('u','x12','a_s_223','a_c_112'), Qv('u','x12','a_s_224','a_c_112')],(0.08333333333333333-0j)), Term([G(5,'a_s_203','a_s_204'), G(5,'a_s_209','a_s_210'), G(5,'a_s_221','a_s_222'), G(5,'a_s_225','a_s_226')],[Qb('u','x21','a_s_203','a_c_102'), Qv('u','x21','a_s_204','a_c_102'), Qb('d','x22','a_s_209','a_c_105'), Qv('d','x22','a_s_210','a_c_105'), Qb('d','x11','a_s_221','a_c_111'), Qv('d','x11','a_s_222','a_c_111'), Qb('d','x12','a_s_225','a_c_113'), Qv('d','x12','a_s_226','a_c_113')],(-0.08333333333333333+0j)), Term([G(5,'a_s_203','a_s_204'), G(5,'a_s_209','a_s_210'), G(5,'a_s_227','a_s_228'), G(5,'a_s_229','a_s_230')],[Qb('u','x21','a_s_203','a_c_102'), Qv('u','x21','a_s_204','a_c_102'), Qb('d','x22','a_s_209','a_c_105'), Qv('d','x22','a_s_210','a_c_105'), Qb('d','x11','a_s_227','a_c_114'), Qv('u','x11','a_s_228','a_c_114'), Qb('u','x12','a_s_229','a_c_115'), Qv('d','x12','a_s_230','a_c_115')],(-0.16666666666666669+0j)), Term([G(5,'a_s_203','a_s_204'), G(5,'a_s_209','a_s_210'), G(5,'a_s_231','a_s_232'), G(5,'a_s_233','a_s_234')],[Qb('u','x21','a_s_203','a_c_102'), Qv('u','x21','a_s_204','a_c_102'), Qb('d','x22','a_s_209','a_c_105'), Qv('d','x22','a_s_210','a_c_105'), Qb('u','x11','a_s_231','a_c_116'), Qv('d','x11','a_s_232','a_c_116'), Qb('d','x12','a_s_233','a_c_117'), Qv('u','x12','a_s_234','a_c_117')],(-0.16666666666666669+0j)), Term([G(5,'a_s_205','a_s_206'), G(5,'a_s_207','a_s_208'), G(5,'a_s_219','a_s_220'), G(5,'a_s_223','a_s_224')],[Qb('d','x21','a_s_205','a_c_103'), Qv('d','x21','a_s_206','a_c_103'), Qb('u','x22','a_s_207','a_c_104'), Qv('u','x22','a_s_208','a_c_104'), Qb('u','x11','a_s_219','a_c_110'), Qv('u','x11','a_s_220','a_c_110'), Qb('u','x12','a_s_223','a_c_112'), Qv('u','x12','a_s_224','a_c_112')],(-0.08333333333333333+0j)), Term([G(5,'a_s_205','a_s_206'), G(5,'a_s_207','a_s_208'), G(5,'a_s_219','a_s_220'), G(5,'a_s_225','a_s_226')],[Qb('d','x21','a_s_205','a_c_103'), Qv('d','x21','a_s_206','a_c_103'), Qb('u','x22','a_s_207','a_c_104'), Qv('u','x22','a_s_208','a_c_104'), Qb('u','x11','a_s_219','a_c_110'), Qv('u','x11','a_s_220','a_c_110'), Qb('d','x12','a_s_225','a_c_113'), Qv('d','x12','a_s_226','a_c_113')],(0.08333333333333333-0j)), Term([G(5,'a_s_205','a_s_206'), G(5,'a_s_207','a_s_208'), G(5,'a_s_221','a_s_222'), G(5,'a_s_223','a_s_224')],[Qb('d','x21','a_s_205','a_c_103'), Qv('d','x21','a_s_206','a_c_103'), Qb('u','x22','a_s_207','a_c_104'), Qv('u','x22','a_s_208','a_c_104'), Qb('d','x11','a_s_221','a_c_111'), Qv('d','x11','a_s_222','a_c_111'), Qb('u','x12','a_s_223','a_c_112'), Qv('u','x12','a_s_224','a_c_112')],(0.08333333333333333-0j)), Term([G(5,'a_s_205','a_s_206'), G(5,'a_s_207','a_s_208'), G(5,'a_s_221','a_s_222'), G(5,'a_s_225','a_s_226')],[Qb('d','x21','a_s_205','a_c_103'), Qv('d','x21','a_s_206','a_c_103'), Qb('u','x22','a_s_207','a_c_104'), Qv('u','x22','a_s_208','a_c_104'), Qb('d','x11','a_s_221','a_c_111'), Qv('d','x11','a_s_222','a_c_111'), Qb('d','x12','a_s_225','a_c_113'), Qv('d','x12','a_s_226','a_c_113')],(-0.08333333333333333+0j)), Term([G(5,'a_s_205','a_s_206'), G(5,'a_s_207','a_s_208'), G(5,'a_s_227','a_s_228'), G(5,'a_s_229','a_s_230')],[Qb('d','x21','a_s_205','a_c_103'), Qv('d','x21','a_s_206','a_c_103'), Qb('u','x22','a_s_207','a_c_104'), Qv('u','x22','a_s_208','a_c_104'), Qb('d','x11','a_s_227','a_c_114'), Qv('u','x11','a_s_228','a_c_114'), Qb('u','x12','a_s_229','a_c_115'), Qv('d','x12','a_s_230','a_c_115')],(-0.16666666666666669+0j)), Term([G(5,'a_s_205','a_s_206'), G(5,'a_s_207','a_s_208'), G(5,'a_s_231','a_s_232'), G(5,'a_s_233','a_s_234')],[Qb('d','x21','a_s_205','a_c_103'), Qv('d','x21','a_s_206','a_c_103'), Qb('u','x22','a_s_207','a_c_104'), Qv('u','x22','a_s_208','a_c_104'), Qb('u','x11','a_s_231','a_c_116'), Qv('d','x11','a_s_232','a_c_116'), Qb('d','x12','a_s_233','a_c_117'), Qv('u','x12','a_s_234','a_c_117')],(-0.16666666666666669+0j)), Term([G(5,'a_s_205','a_s_206'), G(5,'a_s_209','a_s_210'), G(5,'a_s_219','a_s_220'), G(5,'a_s_223','a_s_224')],[Qb('d','x21','a_s_205','a_c_103'), Qv('d','x21','a_s_206','a_c_103'), Qb('d','x22','a_s_209','a_c_105'), Qv('d','x22','a_s_210','a_c_105'), Qb('u','x11','a_s_219','a_c_110'), Qv('u','x11','a_s_220','a_c_110'), Qb('u','x12','a_s_223','a_c_112'), Qv('u','x12','a_s_224','a_c_112')],(0.08333333333333333+0j)), Term([G(5,'a_s_205','a_s_206'), G(5,'a_s_209','a_s_210'), G(5,'a_s_219','a_s_220'), G(5,'a_s_225','a_s_226')],[Qb('d','x21','a_s_205','a_c_103'), Qv('d','x21','a_s_206','a_c_103'), Qb('d','x22','a_s_209','a_c_105'), Qv('d','x22','a_s_210','a_c_105'), Qb('u','x11','a_s_219','a_c_110'), Qv('u','x11','a_s_220','a_c_110'), Qb('d','x12','a_s_225','a_c_113'), Qv('d','x12','a_s_226','a_c_113')],(-0.08333333333333333+0j)), Term([G(5,'a_s_205','a_s_206'), G(5,'a_s_209','a_s_210'), G(5,'a_s_221','a_s_222'), G(5,'a_s_223','a_s_224')],[Qb('d','x21','a_s_205','a_c_103'), Qv('d','x21','a_s_206','a_c_103'), Qb('d','x22','a_s_209','a_c_105'), Qv('d','x22','a_s_210','a_c_105'), Qb('d','x11','a_s_221','a_c_111'), Qv('d','x11','a_s_222','a_c_111'), Qb('u','x12','a_s_223','a_c_112'), Qv('u','x12','a_s_224','a_c_112')],(-0.08333333333333333+0j)), Term([G(5,'a_s_205','a_s_206'), G(5,'a_s_209','a_s_210'), G(5,'a_s_221','a_s_222'), G(5,'a_s_225','a_s_226')],[Qb('d','x21','a_s_205','a_c_103'), Qv('d','x21','a_s_206','a_c_103'), Qb('d','x22','a_s_209','a_c_105'), Qv('d','x22','a_s_210','a_c_105'), Qb('d','x11','a_s_221','a_c_111'), Qv('d','x11','a_s_222','a_c_111'), Qb('d','x12','a_s_225','a_c_113'), Qv('d','x12','a_s_226','a_c_113')],(0.08333333333333333+0j)), Term([G(5,'a_s_205','a_s_206'), G(5,'a_s_209','a_s_210'), G(5,'a_s_227','a_s_228'), G(5,'a_s_229','a_s_230')],[Qb('d','x21','a_s_205','a_c_103'), Qv('d','x21','a_s_206','a_c_103'), Qb('d','x22','a_s_209','a_c_105'), Qv('d','x22','a_s_210','a_c_105'), Qb('d','x11','a_s_227','a_c_114'), Qv('u','x11','a_s_228','a_c_114'), Qb('u','x12','a_s_229','a_c_115'), Qv('d','x12','a_s_230','a_c_115')],(0.16666666666666669+0j)), Term([G(5,'a_s_205','a_s_206'), G(5,'a_s_209','a_s_210'), G(5,'a_s_231','a_s_232'), G(5,'a_s_233','a_s_234')],[Qb('d','x21','a_s_205','a_c_103'), Qv('d','x21','a_s_206','a_c_103'), Qb('d','x22','a_s_209','a_c_105'), Qv('d','x22','a_s_210','a_c_105'), Qb('u','x11','a_s_231','a_c_116'), Qv('d','x11','a_s_232','a_c_116'), Qb('d','x12','a_s_233','a_c_117'), Qv('u','x12','a_s_234','a_c_117')],(0.16666666666666669+0j)), Term([G(5,'a_s_211','a_s_212'), G(5,'a_s_213','a_s_214'), G(5,'a_s_219','a_s_220'), G(5,'a_s_223','a_s_224')],[Qb('u','x21','a_s_211','a_c_106'), Qv('d','x21','a_s_212','a_c_106'), Qb('d','x22','a_s_213','a_c_107'), Qv('u','x22','a_s_214','a_c_107'), Qb('u','x11','a_s_219','a_c_110'), Qv('u','x11','a_s_220','a_c_110'), Qb('u','x12','a_s_223','a_c_112'), Qv('u','x12','a_s_224','a_c_112')],(0.16666666666666669+0j)), Term([G(5,'a_s_211','a_s_212'), G(5,'a_s_213','a_s_214'), G(5,'a_s_219','a_s_220'), G(5,'a_s_225','a_s_226')],[Qb('u','x21','a_s_211','a_c_106'), Qv('d','x21','a_s_212','a_c_106'), Qb('d','x22','a_s_213','a_c_107'), Qv('u','x22','a_s_214','a_c_107'), Qb('u','x11','a_s_219','a_c_110'), Qv('u','x11','a_s_220','a_c_110'), Qb('d','x12','a_s_225','a_c_113'), Qv('d','x12','a_s_226','a_c_113')],(-0.16666666666666669+0j)), Term([G(5,'a_s_211','a_s_212'), G(5,'a_s_213','a_s_214'), G(5,'a_s_221','a_s_222'), G(5,'a_s_223','a_s_224')],[Qb('u','x21','a_s_211','a_c_106'), Qv('d','x21','a_s_212','a_c_106'), Qb('d','x22','a_s_213','a_c_107'), Qv('u','x22','a_s_214','a_c_107'), Qb('d','x11','a_s_221','a_c_111'), Qv('d','x11','a_s_222','a_c_111'), Qb('u','x12','a_s_223','a_c_112'), Qv('u','x12','a_s_224','a_c_112')],(-0.16666666666666669+0j)), Term([G(5,'a_s_211','a_s_212'), G(5,'a_s_213','a_s_214'), G(5,'a_s_221','a_s_222'), G(5,'a_s_225','a_s_226')],[Qb('u','x21','a_s_211','a_c_106'), Qv('d','x21','a_s_212','a_c_106'), Qb('d','x22','a_s_213','a_c_107'), Qv('u','x22','a_s_214','a_c_107'), Qb('d','x11','a_s_221','a_c_111'), Qv('d','x11','a_s_222','a_c_111'), Qb('d','x12','a_s_225','a_c_113'), Qv('d','x12','a_s_226','a_c_113')],(0.16666666666666669+0j)), Term([G(5,'a_s_211','a_s_212'), G(5,'a_s_213','a_s_214'), G(5,'a_s_227','a_s_228'), G(5,'a_s_229','a_s_230')],[Qb('u','x21','a_s_211','a_c_106'), Qv('d','x21','a_s_212','a_c_106'), Qb('d','x22','a_s_213','a_c_107'), Qv('u','x22','a_s_214','a_c_107'), Qb('d','x11','a_s_227','a_c_114'), Qv('u','x11','a_s_228','a_c_114'), Qb('u','x12','a_s_229','a_c_115'), Qv('d','x12','a_s_230','a_c_115')],(0.3333333333333334+0j)), Term([G(5,'a_s_211','a_s_212'), G(5,'a_s_213','a_s_214'), G(5,'a_s_231','a_s_232'), G(5,'a_s_233','a_s_234')],[Qb('u','x21','a_s_211','a_c_106'), Qv('d','x21','a_s_212','a_c_106'), Qb('d','x22','a_s_213','a_c_107'), Qv('u','x22','a_s_214','a_c_107'), Qb('u','x11','a_s_231','a_c_116'), Qv('d','x11','a_s_232','a_c_116'), Qb('d','x12','a_s_233','a_c_117'), Qv('u','x12','a_s_234','a_c_117')],(0.3333333333333334+0j)), Term([G(5,'a_s_215','a_s_216'), G(5,'a_s_217','a_s_218'), G(5,'a_s_219','a_s_220'), G(5,'a_s_223','a_s_224')],[Qb('d','x21','a_s_215','a_c_108'), Qv('u','x21','a_s_216','a_c_108'), Qb('u','x22','a_s_217','a_c_109'), Qv('d','x22','a_s_218','a_c_109'), Qb('u','x11','a_s_219','a_c_110'), Qv('u','x11','a_s_220','a_c_110'), Qb('u','x12','a_s_223','a_c_112'), Qv('u','x12','a_s_224','a_c_112')],(0.16666666666666669+0j)), Term([G(5,'a_s_215','a_s_216'), G(5,'a_s_217','a_s_218'), G(5,'a_s_219','a_s_220'), G(5,'a_s_225','a_s_226')],[Qb('d','x21','a_s_215','a_c_108'), Qv('u','x21','a_s_216','a_c_108'), Qb('u','x22','a_s_217','a_c_109'), Qv('d','x22','a_s_218','a_c_109'), Qb('u','x11','a_s_219','a_c_110'), Qv('u','x11','a_s_220','a_c_110'), Qb('d','x12','a_s_225','a_c_113'), Qv('d','x12','a_s_226','a_c_113')],(-0.16666666666666669+0j)), Term([G(5,'a_s_215','a_s_216'), G(5,'a_s_217','a_s_218'), G(5,'a_s_221','a_s_222'), G(5,'a_s_223','a_s_224')],[Qb('d','x21','a_s_215','a_c_108'), Qv('u','x21','a_s_216','a_c_108'), Qb('u','x22','a_s_217','a_c_109'), Qv('d','x22','a_s_218','a_c_109'), Qb('d','x11','a_s_221','a_c_111'), Qv('d','x11','a_s_222','a_c_111'), Qb('u','x12','a_s_223','a_c_112'), Qv('u','x12','a_s_224','a_c_112')],(-0.16666666666666669+0j)), Term([G(5,'a_s_215','a_s_216'), G(5,'a_s_217','a_s_218'), G(5,'a_s_221','a_s_222'), G(5,'a_s_225','a_s_226')],[Qb('d','x21','a_s_215','a_c_108'), Qv('u','x21','a_s_216','a_c_108'), Qb('u','x22','a_s_217','a_c_109'), Qv('d','x22','a_s_218','a_c_109'), Qb('d','x11','a_s_221','a_c_111'), Qv('d','x11','a_s_222','a_c_111'), Qb('d','x12','a_s_225','a_c_113'), Qv('d','x12','a_s_226','a_c_113')],(0.16666666666666669+0j)), Term([G(5,'a_s_215','a_s_216'), G(5,'a_s_217','a_s_218'), G(5,'a_s_227','a_s_228'), G(5,'a_s_229','a_s_230')],[Qb('d','x21','a_s_215','a_c_108'), Qv('u','x21','a_s_216','a_c_108'), Qb('u','x22','a_s_217','a_c_109'), Qv('d','x22','a_s_218','a_c_109'), Qb('d','x11','a_s_227','a_c_114'), Qv('u','x11','a_s_228','a_c_114'), Qb('u','x12','a_s_229','a_c_115'), Qv('d','x12','a_s_230','a_c_115')],(0.3333333333333334+0j)), Term([G(5,'a_s_215','a_s_216'), G(5,'a_s_217','a_s_218'), G(5,'a_s_231','a_s_232'), G(5,'a_s_233','a_s_234')],[Qb('d','x21','a_s_215','a_c_108'), Qv('u','x21','a_s_216','a_c_108'), Qb('u','x22','a_s_217','a_c_109'), Qv('d','x22','a_s_218','a_c_109'), Qb('u','x11','a_s_231','a_c_116'), Qv('d','x11','a_s_232','a_c_116'), Qb('d','x12','a_s_233','a_c_117'), Qv('u','x12','a_s_234','a_c_117')],(0.3333333333333334+0j))])
    cexpr = contract_simplify_round_compile(expr, is_isospin_symmetric_limit = True)
    print(display_cexpr(cexpr))
    cexpr.collect_op()
    print(display_cexpr(cexpr))
