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

from auto_contractor.compile import *
import gpt as g
import numpy as np
import qlat as q
import copy
import cmath
import math

def get_spin_matrix(op):
    assert op.otype == "G"
    assert op.s1 == "auto" and op.s2 == "auto"
    assert op.tag in [0, 1, 2, 3, 5]
    return g.gamma[op.tag]

def ascontiguoustensor(x):
    # isinstance(x, g.core.tensor)
    return g.tensor(np.ascontiguousarray(x.array), x.otype)

def as_mspincolor(x):
    if isinstance(x, g.core.tensor):
        return ascontiguoustensor(x)
    else:
        return g.tensor(np.ascontiguousarray(np.array(x)), g.ot_matrix_spin_color(4, 3))

def adj_msc(x):
    # isinstance(x, g.core.tensor)
    x = g.adj(x)
    return ascontiguoustensor(x)

def g5_herm(x):
    x_h = ascontiguoustensor(
            ascontiguoustensor(
                g.gamma[5]
                * adj_msc(
                    as_mspincolor(x)))
                * g.gamma[5])
    return x_h

class AmaVal:

    def __init__(self, val = None, corrections = None):
        # should have the following:
        # self.val: sloppy value
        # self.corrections: list [ (val, description_dict,), ... ]
        # description_dict:
        # { source_specification:
        #   (accuracy_level_relative_to_the_basic_accuracy,
        #    probability_of_having_this_accuracy,),
        # }
        self.val = val
        if corrections is None:
            self.corrections = []
        else:
            self.corrections = corrections

def ama_apply1(f, x):
    if not isinstance(x, AmaVal):
        return f(x)
    elif isinstance(x, AmaVal):
        val = f(x.val)
        corrections = [ (f(v), d,) for v, d in x.corrections ]
        return AmaVal(val, corrections)
    else:
        assert False

def merge_description_dict(d1, d2):
    sd1 = set(d1)
    sd2 = set(d2)
    common_keys = sd1 & sd2
    for key in common_keys:
        if d1[key] != d2[key]:
            return None
    new_keys = sd2 - sd1
    d = d1.copy()
    for key in new_keys:
        d[key] = d2[key]
    return d

def ama_apply2(f, x, y):
    if not isinstance(x, AmaVal) and not isinstance(y, AmaVal):
        return f(x, y)
    elif isinstance(x, AmaVal) and not isinstance(y, AmaVal):
        def f1(x1):
            return f(x1, y)
        return ama_apply1(f1, x)
    elif not isinstance(x, AmaVal) and isinstance(y, AmaVal):
        def f1(y1):
            return f(x, y1)
        return ama_apply1(f1, y)
    elif isinstance(x, AmaVal) and isinstance(y, AmaVal):
        val = f(x.val, y.val)
        corrections = []
        for v_x, d_x in x.corrections:
            for v_y, d_y in y.corrections:
                d = merge_description_dict(d_x, d_y)
                if d is not None:
                    corrections.append((f(v_x, v_y), d,))
        return AmaVal(val, corrections)
    else:
        assert False

def ama_extract(x):
    if not isinstance(x, AmaVal):
        return x
    elif isinstance(x, AmaVal):
        val = x.val
        corrections = x.corrections
        assert isinstance(corrections, list)
        assert corrections
        # keys = [ source_specification, ... ]
        keys = list(corrections[0][1].keys())
        def get_level_prob(key):
            s = set([ d[key] for v, d in corrections ])
            return sorted(list(s))
        dict_level_prob = { k: get_level_prob(k) for k in keys }
        for k, v in dict_level_prob.items():
            assert len(v) >= 2
        # dict_level[key] = list of accuracy levels for this key (source_specification)
        dict_level = { k: [ l for l, prob in v ] for k, v in dict_level_prob.items() }
        # dict_prob[key] = list of probability_of_having_this_accuracy
        dict_prob = { k: [ prob for l, prob in v ] for k, v in dict_level_prob.items() }
        # dict_val[(level, ...)] = v
        dict_val = { tuple([ d[k][0] for k in keys ]): v for v, d in corrections }
        def ama_corr(fixed_levels, remaining_keys):
            if not remaining_keys:
                return dict_val[tuple(fixed_levels)]
            else:
                key = remaining_keys[0]
                rest_keys = remaining_keys[1:]
                levels = dict_level[key]
                probs = dict_prob[key]
                vals = [ ama_corr(fixed_levels + [ l, ], rest_keys) for l in levels ]
                corr = vals[0]
                for i in range(1, len(vals)):
                    corr += (vals[i] - vals[i - 1]) / probs[i]
                return corr
        return ama_corr([], keys)
    else:
        assert False

def eval_op_term_expr(expr, variable_dict, positions_dict, get_prop):
    def l_eval(x):
        if isinstance(x, list):
            ans = l_eval(x[0])
            def f(x, y):
                return ascontiguoustensor(x * y)
            for op in x[1:]:
                ans = ama_apply2(f, ans, l_eval(op))
            return ans
        elif isinstance(x, Op):
            if x.otype == "S":
                flavor = x.f
                xg_snk = positions_dict[x.p1]
                xg_src = positions_dict[x.p2]
                return get_prop(flavor, xg_snk, xg_src)
            elif x.otype == "G":
                return get_spin_matrix(x)
            elif x.otype == "Tr" and len(x.ops) == 2:
                def f(x, y):
                    return g.trace(x * y)
                return ama_apply2(f, l_eval(x.ops[0]), l_eval(x.ops[1]))
            elif x.otype == "Tr":
                if not x.ops:
                    return 1
                start_idx = None
                for i in range(len(x.ops)):
                    if x.ops[i].otype != "G":
                        start_idx = i
                        break
                assert start_idx is not None
                ans = l_eval(x.ops[start_idx])
                def f(x, y):
                    return ascontiguoustensor(x * y)
                for op in x.ops[start_idx + 1:] + x.ops[:start_idx]:
                    ans = ama_apply2(f, ans, l_eval(op))
                return ama_apply1(g.trace, ans)
            elif x.otype == "Var":
                return variable_dict[x.name]
            else:
                q.displayln_info(f"eval_op_term_expr: ERROR: l_eval({x})")
                assert False
        elif isinstance(x, Term):
            assert not x.a_ops
            ans = x.coef
            def f(x, y):
                return x * y
            for op in x.c_ops:
                ans = ama_apply2(f, ans, l_eval(op))
            return ans
        elif isinstance(x, Expr):
            ans = 0
            def f(x, y):
                return x + y
            for term in x.terms:
                ans = ama_apply2(f, ans, l_eval(term))
            return ans
        else:
            q.displayln_info(f"eval_op_term_expr: ERROR: l_eval({x})")
            assert False
    return l_eval(expr)

@q.timer
def eval_cexpr(cexpr : CExpr, *, positions_dict, get_prop, is_only_total):
    # interface function
    # return 1 dimensional np.array
    for pos in cexpr.positions:
        assert pos in positions_dict
    variable_dict = {}
    for name, op in cexpr.variables:
        variable_dict[name] = eval_op_term_expr(op, variable_dict, positions_dict, get_prop)
    tvals = { name : ama_extract(eval_op_term_expr(term, variable_dict, positions_dict, get_prop)) for name, term in cexpr.named_terms }
    evals = { name : sum([ tvals[tname] for tname in expr ]) for name, expr in cexpr.named_exprs }
    if is_only_total in [ True, "total", ]:
        return np.array([ evals[name] for name, expr in cexpr.named_exprs])
    elif is_only_total in [ "typed_total", ]:
        tevals = { name : sum([ tvals[tname] for tname in expr ]) for name, expr in cexpr.named_typed_exprs }
        return np.array([ tevals[name] for name, expr in cexpr.named_typed_exprs ] + [ evals[name] for name, expr in cexpr.named_exprs])
    elif is_only_total in [ False, "term", ]:
        tevals = { name : sum([ tvals[tname] for tname in expr ]) for name, expr in cexpr.named_typed_exprs }
        return np.array([ tvals[name] for name, term in cexpr.named_terms ] + [ tevals[name] for name, expr in cexpr.named_typed_exprs ] + [ evals[name] for name, expr in cexpr.named_exprs])
    else:
        assert False

def sqr_component(x):
    return x.real * x.real + 1j * x.imag * x.imag

def sqrt_component(x):
    return math.sqrt(x.real) + 1j * math.sqrt(x.imag)

def sqr_component_array(arr):
    return np.array([ sqr_component(x) for x in arr ])

def sqrt_component_array(arr):
    return np.array([ sqrt_component(x) for x in arr ])

@q.timer
def get_cexpr_names(cexpr, is_only_total = "total"):
    if is_only_total in [ True, "total", ]:
        names = [ name for name, expr in cexpr.named_exprs ]
    elif is_only_total in [ "typed_total", ]:
        names = ([ name for name, expr in cexpr.named_typed_exprs ]
                + [ name for name, expr in cexpr.named_exprs ])
    elif is_only_total in [ False, "term", ]:
        names = ([ name for name, term in cexpr.named_terms ]
                + [ name for name, expr in cexpr.named_typed_exprs ]
                + [ name for name, expr in cexpr.named_exprs ])
    else:
        assert False
    return names

def get_mpi_chunk(total_list, *, rng_state = None):
    if rng_state is None:
        rng_state = q.RngState("get_mpi_chunk")
    total = len(total_list)
    id_worker = q.get_id_node()
    num_worker = q.get_num_node()
    size_max = (total - 1) // num_worker + 1;
    start = min(id_worker * size_max, total);
    stop = min(start + size_max, total);
    # size = stop - start;
    total_list = q.random_permute(total_list, rng_state)
    return total_list[start:stop]

@q.timer
def eval_cexpr_simulation(cexpr : CExpr, *, positions_dict_maker, trial_indices, get_prop, is_only_total = "total"):
    # interface function
    if len(trial_indices) == 0:
        return None
    names = get_cexpr_names(cexpr, is_only_total)
    num_value = len(names)
    results = None
    num_fac = None
    for idx in trial_indices:
        positions_dict, facs = positions_dict_maker(idx)
        if num_fac is None:
            num_fac = len(facs)
            results = [ [] for i in range(num_fac) ]
        else:
            assert num_fac == len(facs)
        res = eval_cexpr(cexpr, positions_dict = positions_dict, get_prop = get_prop, is_only_total = is_only_total)
        assert len(res) == num_value
        for i in range(num_fac):
            results[i].append(facs[i] * res)
    for i in range(num_fac):
        assert len(results[i]) == len(trial_indices)
    total_num_trials = q.glb_sum(len(trial_indices))
    ld_info = [
            [ "fac", num_fac ],
            [ "name", num_value, ],
            ]
    ld = q.mk_lat_data(ld_info)
    for i in range(num_fac):
        ld[(i,)] = sum(results[i])
    ld_avg = (1 / total_num_trials) * q.glb_sum(ld)
    ld = q.mk_lat_data(ld_info)
    for i in range(num_fac):
        ld[(i,)] = sum([ sqr_component_array(r - ld_avg[(i,)]) for r in results[i] ])
    ld_var = (1 / total_num_trials) * q.glb_sum(ld)
    summaries = []
    for i in range(num_fac):
        results_avg = ld_avg[(i,)]
        results_err = sqrt_component_array(ld_var[(i,)]) / math.sqrt(total_num_trials)
        summary = {}
        for i in range(len(names)):
            summary[names[i]] = [ results_avg[i], results_err[i], ]
        summaries.append(summary)
    return summaries

@q.timer
def positions_dict_maker_example_1(idx):
    total_site = [4, 4, 4, 16,]
    rs = q.RngState("positions_dict_maker").split(str(idx))
    t2 = 3
    x1 = rs.c_rand_gen(total_site)
    x2 = rs.c_rand_gen(total_site)
    x2[3] = (x1[3] + t2) % total_site[3]
    pd = {
            "x1" : x1,
            "x2" : x2,
            }
    facs = [1.0,]
    return pd, facs

@q.timer
def positions_dict_maker_example_2(idx):
    total_site = [4, 4, 4, 16,]
    rs = q.RngState("positions_dict_maker").split(str(idx))
    lmom1 = [0.0, 0.0, 1.0, 0.0,]
    lmom2 = [0.0, 0.0, -1.0, 0.0,]
    t2 = 2
    t3 = 5
    x1 = rs.c_rand_gen(total_site)
    x2 = rs.c_rand_gen(total_site)
    x3 = rs.c_rand_gen(total_site)
    x2[3] = (x1[3] + t2) % total_site[3]
    x3[3] = (x1[3] + t3) % total_site[3]
    pd = {
            "x1" : x1,
            "x2" : x2,
            "x3" : x3,
            }
    phase = 0.0
    for mu in range(4):
        mom1 = 2.0 * math.pi / total_site[mu] * lmom1[mu]
        mom2 = 2.0 * math.pi / total_site[mu] * lmom2[mu]
        phase += mom1 * x1[mu]
        phase += mom2 * x2[mu]
    facs = [1.0, cmath.rect(1.0, phase),]
    return pd, facs

if __name__ == "__main__":
    print(positions_dict_maker_example_1(0))
    print(positions_dict_maker_example_2(0))

