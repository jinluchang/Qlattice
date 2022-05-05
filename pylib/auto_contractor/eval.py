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
from auto_contractor.ama import *

from auto_contractor.eval_sc_np import *

import numpy as np
import qlat as q
import copy
import cmath
import math

def eval_op_term_expr(expr, variable_dict, positions_dict, get_prop):
    def l_eval(x):
        if isinstance(x, list):
            ans = l_eval(x[0])
            def f(x, y):
                return x * y
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
            elif x.otype == "Tr":
                if len(x.ops) == 2:
                    def f(x, y):
                        return msc_trace2(x, y)
                    return ama_apply2(f, l_eval(x.ops[0]), l_eval(x.ops[1]))
                elif len(x.ops) == 1:
                    return ama_apply1(msc_trace, l_eval(x.ops[0]))
                elif len(x.ops) == 0:
                    return 1
                else:
                    assert len(x.ops) > 2
                    def f(x, y):
                        return x * y
                    ans = l_eval(x.ops[0])
                    for op in x.ops[1:-1]:
                        ans = ama_apply2(f, ans, l_eval(op))
                    return ama_apply2(msc_trace2, ans, l_eval(x.ops[-1]))
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
def get_cexpr_names(cexpr, *, is_only_total = "total"):
    if is_only_total in [ True, "total", ]:
        names = [ name for name, expr in cexpr.named_exprs ]
    elif is_only_total in [ "typed_total", ]:
        names = (
                [ name for name, expr in cexpr.named_exprs ]
                + [ name for name, expr in cexpr.named_typed_exprs ]
                )
    elif is_only_total in [ False, "term", ]:
        names = (
                [ name for name, expr in cexpr.named_exprs ]
                + [ name for name, expr in cexpr.named_typed_exprs ]
                + [ name for name, term in cexpr.named_terms ]
                )
    else:
        assert False
    return names

@q.timer
def eval_cexpr(cexpr : CExpr, *, positions_dict, get_prop, is_only_total = "total"):
    # interface function
    # return 1 dimensional np.array
    # xg = positions_dict[position]
    # mat_mspincolor = get_prop(flavor, xg_snk, xg_src)
    # is_only_total = "total", "typed_total", "term"
    # e.g. ("point-snk", [ 1, 2, 3, 4, ]) = positions_dict["x_1"]
    # e.g. flavor = "l"
    # e.g. xg_snk = ("point-snk", [ 1, 2, 3, 4, ])
    for pos in cexpr.positions:
        assert pos in positions_dict
    variable_dict = {}
    tvals = {}
    @q.timer
    def eval_cexpr_set_vars():
        for name, op in cexpr.variables:
            variable_dict[name] = eval_op_term_expr(op, variable_dict, positions_dict, get_prop)
    @q.timer
    def eval_cexpr_set_terms():
        for name, term in cexpr.named_terms:
            tvals[name] = ama_extract(eval_op_term_expr(term, variable_dict, positions_dict, get_prop))
    @q.timer
    def eval_cexpr_return_exprs():
        if is_only_total in [ True, "total", ]:
            return np.array(
                    [ sum([ coef * tvals[tname] for coef, tname in expr ]) for name, expr in cexpr.named_exprs ]
                    )
        elif is_only_total in [ "typed_total", ]:
            return np.array(
                    [ sum([ coef * tvals[tname] for coef, tname in expr ]) for name, expr in cexpr.named_exprs ]
                    + [ sum([ coef * tvals[tname] for coef, tname in expr ]) for name, expr in cexpr.named_typed_exprs ]
                    )
        elif is_only_total in [ False, "term", ]:
            tevals = { name : sum([ coef * tvals[tname] for coef, tname in expr ]) for name, expr in cexpr.named_typed_exprs }
            return np.array(
                    [ sum([ coef * tvals[tname] for coef, tname in expr ]) for name, expr in cexpr.named_exprs ]
                    + [ sum([ coef * tvals[tname] for coef, tname in expr ]) for name, expr in cexpr.named_typed_exprs ]
                    + [ tvals[name] for name, term in cexpr.named_terms ]
                    )
        else:
            assert False
    eval_cexpr_set_vars()
    eval_cexpr_set_terms()
    return eval_cexpr_return_exprs()

def make_rand_spin_color_matrix(rng_state):
    rs = rng_state
    return SpinColorMatrix(np.array(
        [ rs.u_rand_gen() + 1j * rs.u_rand_gen() for i in range(144) ],
        dtype = complex))

@q.timer
def benchmark_eval_cexpr(cexpr : CExpr, *, is_only_total = "total", benchmark_size = 10, benchmark_num = 10, rng_state = None):
    if rng_state is None:
        rng_state = q.RngState("benchmark_eval_cexpr")
    positions_dict = {}
    prop_dict = {}
    for pos in cexpr.positions:
        positions_dict[pos] = pos
    for pos in cexpr.positions:
        prop_dict[(pos, pos,)] = make_rand_spin_color_matrix(rng_state.split(pos))
    @q.timer
    def get_prop(flavor, xg_snk, xg_src):
        return prop_dict[(xg_src, xg_src)]
    @q.timer_verbose
    def benchmark_eval_cexpr_run():
        for k in range(benchmark_size):
            eval_cexpr(cexpr, positions_dict = positions_dict, get_prop = get_prop, is_only_total = is_only_total)
    q.displayln_info(f"benchmark_eval_cexpr: benchmark_size={benchmark_size} is_only_total={is_only_total}")
    def run(*args):
        for i in range(benchmark_num):
            benchmark_eval_cexpr_run()
    q.parallel_map(1, run, [ None, ])

def sqr_component(x):
    return x.real * x.real + 1j * x.imag * x.imag

def sqrt_component(x):
    return math.sqrt(x.real) + 1j * math.sqrt(x.imag)

def sqr_component_array(arr):
    return np.array([ sqr_component(x) for x in arr ])

def sqrt_component_array(arr):
    return np.array([ sqrt_component(x) for x in arr ])

def get_mpi_chunk(total_list, *, rng_state = None):
    # e.g. rng_state = q.RngState("get_mpi_chunk")
    total = len(total_list)
    id_worker = q.get_id_node()
    num_worker = q.get_num_node()
    size_max = (total - 1) // num_worker + 1;
    start = min(id_worker * size_max, total);
    stop = min(start + size_max, total);
    # size = stop - start;
    if rng_state is not None:
        assert isinstance(rng_state, q.RngState)
        total_list = q.random_permute(total_list, rng_state)
    return total_list[start:stop]

@q.timer
def eval_cexpr_simulation(cexpr : CExpr, *, positions_dict_maker, trial_indices, get_prop, is_only_total = "total"):
    # interface function
    # 1. positions_dict_maker(idx) = positions_dict
    # 2. positions_dict_maker(idx) = (positions_dict, facs)
    # positions_dict[pos] = xg
    # get_prop(flavor, xg1, xg2) = prop
    # facs = [ fac1, fac2, ... ]
    # trial_indices = [ idx1, idx2, ... ]
    assert isinstance(trial_indices, list)
    total_num_trials = q.glb_sum(len(trial_indices))
    if total_num_trials == 0:
        return None
    if trial_indices:
        has_data = 1
    else:
        has_data = 0
    num_has_data = q.glb_sum(has_data)
    names = get_cexpr_names(cexpr, is_only_total = is_only_total)
    num_value = len(names)
    num_fac = None
    results = None
    for idx in trial_indices:
        pdi = positions_dict_maker(idx)
        if isinstance(pdi, dict):
            positions_dict = pdi
            facs = [ 1.0, ]
        else:
            positions_dict, facs = pdi
        if num_fac is None:
            num_fac = len(facs)
            results = [ [] for i in range(num_fac) ]
        else:
            assert num_fac == len(facs)
        res = eval_cexpr(cexpr, positions_dict = positions_dict, get_prop = get_prop, is_only_total = is_only_total)
        assert len(res) == num_value
        for i in range(num_fac):
            results[i].append(facs[i] * res)
    if trial_indices:
        assert num_fac is not None
    else:
        assert num_fac is None
        num_fac = 0
    g_num_fac = q.glb_sum(num_fac)
    if trial_indices:
        assert g_num_fac == num_fac * num_has_data
        for i in range(num_fac):
            assert len(results[i]) == len(trial_indices)
    else:
        num_fac = g_num_fac // num_has_data
    ld_info = [
            [ "fac", num_fac ],
            [ "name", num_value, ],
            ]
    ld = q.mk_lat_data(ld_info)
    if trial_indices:
        for i in range(num_fac):
            ld[(i,)] = sum(results[i])
    ld_avg = (1 / total_num_trials) * q.glb_sum(ld)
    ld = q.mk_lat_data(ld_info)
    if trial_indices:
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
