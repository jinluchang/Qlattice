from auto_contractor.compile import *
import gpt as g
import numpy as np
import qlat as q
import copy
import cmath
import math

@q.timer
def get_prop_psrc(prop_cache, flavor : str, xg_src):
    # prop_cache[flavor][src_p] = prop
    # call load_prop_psrc_all(flavor, path_s) first
    return prop_cache[flavor][f"xg=({xg_src[0]},{xg_src[1]},{xg_src[2]},{xg_src[3]})"]

@q.timer
def get_prop_psnk_psrc(cache, flavor : str, xg_snk, xg_src):
    wm = get_prop_psrc(cache, flavor, xg_src).get_elem(xg_snk)
    return g.tensor(np.ascontiguousarray(np.array(wm)), g.ot_matrix_spin_color(4, 3))

def get_spin_matrix(op):
    assert op.otype == "G"
    assert op.s1 == "auto" and op.s2 == "auto"
    assert op.tag in [0, 1, 2, 3, 5]
    return g.gamma[op.tag]

def ascontiguoustensor(tensor):
    return g.tensor(np.ascontiguousarray(tensor.array), tensor.otype)

def eval_op_term_expr(expr, variables_dict, positions_dict, prop_cache):
    def l_eval(x):
        if isinstance(x, Op):
            if x.otype == "S":
                flavor = x.f
                xg_snk = positions_dict[x.p1]
                xg_src = positions_dict[x.p2]
                return get_prop_psnk_psrc(prop_cache, flavor, xg_snk, xg_src)
            elif x.otype == "G":
                return get_spin_matrix(x)
            elif x.otype == "Tr":
                if not x.ops:
                    return 1
                ans = l_eval(x.ops[0])
                for op in x.ops[1:]:
                    ans = ascontiguoustensor(ans * l_eval(op))
                return g.trace(ans)
            elif x.otype == "Var":
                return variables_dict[x.name]
            else:
                q.displayln_info(f"eval_op_term: ERROR: l_eval({x})")
                assert False
        elif isinstance(x, Term):
            assert not x.a_ops
            ans = x.coef
            for op in x.c_ops:
                ans = ans * l_eval(op)
            return ans
        elif isinstance(x, Expr):
            ans = 0
            for term in x.terms:
                ans = ans + l_eval(term)
            return ans
        else:
            q.displayln_info(f"eval_op_term: ERROR: l_eval({x})")
            assert False
    return g.eval(l_eval(expr))

@q.timer
def eval_cexpr(cexpr : CExpr, *, positions_dict, prop_cache):
    # interface function
    # the last element is the sum
    for pos in cexpr.positions:
        assert pos in positions_dict
    variables_dict = {}
    for name, op in cexpr.variables:
        variables_dict[name] = eval_op_term_expr(op, variables_dict, positions_dict, prop_cache)
    tvals = [ eval_op_term_expr(term, variables_dict, positions_dict, prop_cache) for name, term in cexpr.named_terms ]
    return np.array(tvals + [sum(tvals),])

def sqr_component(x):
    return x.real * x.real + 1j * x.imag * x.imag

def sqrt_component(x):
    return math.sqrt(x.real) + 1j * math.sqrt(x.imag)

def sqr_component_array(arr):
    return np.array([ sqr_component(x) for x in arr ])

def sqrt_component_array(arr):
    return np.array([ sqrt_component(x) for x in arr ])

def contract_simpify_round_compile(expr : Expr, is_isospin_symmetric_limit = True):
    # interface function
    expr = copy.deepcopy(expr)
    expr = contract_expr(expr)
    expr.simplify(is_isospin_symmetric_limit = is_isospin_symmetric_limit)
    cexpr = mk_cexpr(expr.round())
    cexpr.collect_prop()
    return cexpr

@q.timer
def eval_cexpr_simulation(cexpr : CExpr, *, positions_dict_maker, rng_state, trial_indices, total_site, prop_cache):
    # interface function
    results = []
    for idx in trial_indices:
        rs = rng_state.split(str(idx))
        positions_dict, fac = positions_dict_maker(rs, total_site)
        results.append(fac * eval_cexpr(cexpr, positions_dict = positions_dict, prop_cache = prop_cache))
    results_avg = sum(results) / len(results)
    results_err = sqrt_component_array(sum([ sqr_component_array(r - results_avg) for r in results ])) / len(results)
    return results_avg, results_err

@q.timer
def positions_dict_maker_example_1(rs, total_site):
    t2 = 3
    x1 = rs.c_rand_gen(total_site)
    x2 = rs.c_rand_gen(total_site)
    x2[3] = (x1[3] + t2) % total_site[3]
    pd = {
            "x1" : x1,
            "x2" : x2,
            }
    fac = 1.0
    return pd, fac

@q.timer
def positions_dict_maker_example_2(rs, total_site):
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
    fac = cmath.rect(1.0, phase)
    return pd, fac

if __name__ == "__main__":
    rs = q.RngState("3")
    q.displayln_info(rs)
    total_site = [4, 4, 4, 16,]
    print(positions_dict_maker_example_1(rs, total_site))
    print(positions_dict_maker_example_2(rs, total_site))
    print(CExpr([('S_1', S('d','x2','x1')), ('S_2', S('u','x1','x2'))],[('T_1', Term([Tr([G(5), Var('S_1'), G(5), Var('S_2')],'sc')],[],(-1+0j)))],['x1', 'x2']))

