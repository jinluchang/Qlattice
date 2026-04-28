#!/usr/bin/env python3
"""
Pi0 decay lattice QCD data generation script.\n
This script computes meson correlation functions and pi0 decay observables
using lattice QCD propagators. It generates data for pi0 -> gamma gamma
decay analysis, including:\n
- Meson two-point correlators (wall-wall, wall-point, point-wall, point-point)
- Meson tensor currents (J_mu * M)
- Meson mass insertions (m_q * M)
- Meson two-current correlators (J * J * M)
- Pi0 -> J J parity-violating correlators (epsilon_{ijk} x_i J_j J_k * pi0)
- Tadpole current correlators (disconnected contributions)
- Pi0 current correlators (connected contributions)
- Pi0 gamma-gamma disconnected diagrams\n
The script supports multiple gauge ensembles (24D, 48I, 64I, 64I-pq) and
test configurations. It uses adaptive sampling for point-source/point-sink
measurements with probability-weighted estimators.\n
Usage:
    python gpt-qlat-data-gen-pi0-decay.py [--job_tag_list tag1,tag2] [--no-inversion] [--no-contract]
"""

import qlat_gpt as qg

import numpy as np
import sympy
import qlat as q

from qlat_scripts.v1 import (
    ama_extract,
    check_job,
    get_expr_names,
    get_job_seed,
    get_load_path,
    get_param,
    get_r_list,
    get_r_sq_interp_idx_coef_list,
    get_save_path,
    is_test,
    load_path_list,
    load_point_distribution,
    run_eig,
    run_eig_strange,
    run_f_rand_01,
    run_f_weight_from_wsrc_prop_full,
    run_fsel_from_fsel_prob,
    run_fsel_prob,
    run_fsel_split,
    run_field_rand_u1_dict,
    run_gf,
    run_gf_ape,
    run_get_prop,
    run_gt,
    run_params,
    run_prop_psrc,
    run_prop_rand_u1,
    run_prop_smear,
    run_prop_wsrc_full,
    run_prop_wsrc_sparse,
    run_psel_from_psel_prob,
    run_psel_prob,
    run_psel_smear,
    run_psel_smear_median,
    run_psel_split,
    run_r_list,
    run_wi,
    set_param,
)
from auto_contractor.operators import (
    contract_simplify_compile,
    mk_a0_p,
    mk_eta_l,
    mk_eta_s,
    mk_expr,
    mk_fac,
    mk_j5eta_l_mu,
    mk_j5eta_s_mu,
    mk_j5k_mu,
    mk_j5pi_mu,
    mk_j_mu,
    mk_jk_mu,
    mk_jl_mu,
    mk_js_mu,
    mk_k_0,
    mk_k_0_bar,
    mk_k_m,
    mk_k_p,
    mk_k_p_star_mu,
    mk_kappa_p,
    mk_m,
    mk_meson,
    mk_omega,
    mk_pi_0,
    mk_pi_m,
    mk_pi_p,
    mk_proton,
    mk_scalar,
    mk_sym,
    mk_vec5_mu,
    mk_vec_mu,
)
from auto_contractor.eval import (
    benchmark_eval_cexpr,
    cache_compiled_cexpr,
    eval_cexpr,
)

# ----

load_path_list[:] = [
    "results",
    "qcddata",
    "qcddata1",
    "qcddata2",
    "/lustre20/volatile/qcdqedta/qcddata",
    "/lustre20/volatile/decay0n2b/qcddata",
    "/lustre20/volatile/pqpdf/ljin/qcddata",
]

is_cython = not is_test()

# is_cython = True

pname = "pi0_decay"

# ----

@q.timer
def get_cexpr_meson_corr(is_both_prop=True):
    """
    Build compiled expressions for meson two-point correlation functions.\n
    Computes correlators of the form <O2(0) O1(-tsep)> with various meson
    operators including pi+, K+, eta_l, eta_s, kappa, omega, proton, and
    vector/axial currents (j_mu, jl_mu, js_mu, jk_mu, j5pi_mu, j5k_mu).\n
    Args:
        is_both_prop: If True, uses wall-source/wall-sink propagators for
            both operators (Type1 + Type2 diagrams). If False, uses
            wall-source/point-sink (Type1 only).\n
    Returns:
        Compiled expression object for use with eval_cexpr.
    """
    if is_both_prop:
        fn_base = "cache/auto_contract_cexpr/get_cexpr_meson_corr_src_src"
    else:
        fn_base = "cache/auto_contract_cexpr/get_cexpr_meson_corr_snk_src"
    #
    def calc_cexpr():
        diagram_type_dict = dict()
        diagram_type_dict[()] = "T1"
        diagram_type_dict[((("x_1", "x_2"), 1), (("x_2", "x_1"), 1))] = "Type1"
        diagram_type_dict[((("x_1", "x_1"), 1), (("x_2", "x_2"), 1))] = "Type2"
        diagram_type_dict[((("x_1", "x_2"), 3),)] = "TypeB1"
        exprs_1_list = [
            mk_fac(1) + "1",
        ]
        exprs_12_corr_list = []
        exprs_corr_list = []
        exprs_b_corr_list = []
        exprs_12_corr_list += [
            mk_j_mu("x_2", 3) * mk_j_mu("x_1", 3) + "j_t(0) * j_t(-tsep)",
            mk_jl_mu("x_2", 3) * mk_jl_mu("x_1", 3) + "jl_t(0) * jl_t(-tsep)",
            mk_js_mu("x_2", 3) * mk_js_mu("x_1", 3) + "js_t(0) * js_t(-tsep)",
            sum([mk_j_mu("x_2", mu) * mk_j_mu("x_1", mu) for mu in range(3)])
            + "j_i(0) * j_i(-tsep)",
            sum([mk_jl_mu("x_2", mu) * mk_jl_mu("x_1", mu) for mu in range(3)])
            + "jl_i(0) * jl_i(-tsep)",
            sum([mk_js_mu("x_2", mu) * mk_js_mu("x_1", mu) for mu in range(3)])
            + "js_i(0) * js_i(-tsep)",
        ]
        exprs_corr_list += [
            mk_a0_p("x_2", True) * mk_a0_p("x_1") + "a0+^dag(0) * a0+(-tsep)",
            mk_jk_mu("x_2", 3, True) * mk_jk_mu("x_1", 3) + "jk_t^dag(0) * jk_t(-tsep)",
            sum([mk_jk_mu("x_2", mu, True) * mk_jk_mu("x_1", mu) for mu in range(3)])
            + "jk_i^dag(0) * jk_i(-tsep)",
            sum(
                [mk_j5pi_mu("x_2", mu) * mk_j5pi_mu("x_1", mu, True) for mu in range(3)]
            )
            + "j5pi_i(0) * j5pi_i^dag(-tsep)",
            sum([mk_j5k_mu("x_2", mu) * mk_j5k_mu("x_1", mu, True) for mu in range(3)])
            + "j5k_i(0) * j5k_i^dag(-tsep)",
        ]
        pi_1_list = [
            mk_pi_p("x_1") + "pi+(-tsep)",
            mk_j5pi_mu("x_1", 3, True) + "j5pi_t^dag(-tsep)",
        ]
        pi_2_list = [
            mk_pi_p("x_2", True) + "pi+^dag(0)",
            mk_j5pi_mu("x_2", 3) + "j5pi_t(0)",
        ]
        exprs_corr_list += [pi_2 * pi_1 for pi_1 in pi_1_list for pi_2 in pi_2_list]
        kk_1_list = [
            mk_k_p("x_1") + "k+(-tsep)",
            mk_j5k_mu("x_1", 3, True) + "j5k_t^dag(-tsep)",
        ]
        kk_2_list = [
            mk_k_p("x_2", True) + "k+^dag(0)",
            mk_j5k_mu("x_2", 3) + "j5k_t(0)",
        ]
        exprs_corr_list += [kk_2 * kk_1 for kk_1 in kk_1_list for kk_2 in kk_2_list]
        eta_1_list = [
            mk_eta_l("x_1") + "eta_l(-tsep)",
            mk_eta_s("x_1") + "eta_s(-tsep)",
            mk_j5eta_l_mu("x_1", 3, True) + "j5eta_l_t^dag(-tsep)",
            mk_j5eta_s_mu("x_1", 3, True) + "j5eta_s_t^dag(-tsep)",
        ]
        eta_2_list = [
            mk_eta_l("x_2", True) + "eta_l^dag(0)",
            mk_eta_s("x_2", True) + "eta_s^dag(0)",
            mk_j5eta_l_mu("x_2", 3) + "j5eta_l_t(0)",
            mk_j5eta_s_mu("x_2", 3) + "j5eta_s_t(0)",
        ]
        exprs_12_corr_list += [
            eta_2 * eta_1 for eta_1 in eta_1_list for eta_2 in eta_2_list
        ]
        kappa_1_list = [
            mk_kappa_p("x_1") + "kappa+(-tsep)",
            mk_k_p_star_mu("x_1", 3) + "k_p_star_t+(-tsep)",
        ]
        kappa_2_list = [
            mk_kappa_p("x_2", True) + "kappa+^dag(0)",
            mk_k_p_star_mu("x_2", 3, True) + "k_p_star_t+^dag(0)",
        ]
        exprs_corr_list += [
            kappa_2 * kappa_1 for kappa_1 in kappa_1_list for kappa_2 in kappa_2_list
        ]
        for spin in [
            "u3",
            "u1",
            "d1",
            "d3",
        ]:
            omega_1_list = []
            omega_2_list = []
            for baryon_type in [
                "std3",
                "pos3",
            ]:
                omega_1_list += [
                    mk_omega("x_1", spin, baryon_type)
                    + f"omega[{baryon_type},{spin}](-tsep)",
                ]
                omega_2_list += [
                    mk_omega("x_2", spin, baryon_type, True)
                    + f"omega[{baryon_type},{spin}]^dag(0)",
                ]
            exprs_b_corr_list += [
                omega_2 * omega_1
                for omega_1 in omega_1_list
                for omega_2 in omega_2_list
            ]
        for spin in [
            "u",
            "d",
        ]:
            proton_1_list = []
            proton_2_list = []
            for baryon_type in [
                "std",
                "pos",
            ]:
                proton_1_list += [
                    mk_proton("x_1", spin, baryon_type)
                    + f"proton[{baryon_type},{spin}](-tsep)",
                ]
                proton_2_list += [
                    mk_proton("x_2", spin, baryon_type, True)
                    + f"proton[{baryon_type},{spin}]^dag(0)",
                ]
            exprs_b_corr_list += [
                proton_2 * proton_1
                for proton_1 in proton_1_list
                for proton_2 in proton_2_list
            ]
        if is_both_prop:
            exprs_12_corr_list = [
                (
                    expr,
                    "Type1",
                    "Type2",
                )
                for expr in exprs_12_corr_list
            ]
        else:
            exprs_12_corr_list = [
                (
                    expr,
                    "Type1",
                )
                for expr in exprs_12_corr_list
            ]
        exprs_corr_list = [
            (
                expr,
                None,
            )
            for expr in exprs_corr_list
        ]
        exprs_b_corr_list = [
            (
                expr,
                None,
            )
            for expr in exprs_b_corr_list
        ]
        exprs = exprs_1_list + exprs_12_corr_list + exprs_corr_list + exprs_b_corr_list
        cexpr = contract_simplify_compile(
            *exprs, is_isospin_symmetric_limit=True, diagram_type_dict=diagram_type_dict
        )
        return cexpr
    #
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

@q.timer(is_timer_fork=True)
def auto_contract_meson_corr(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob):
    """
    Compute meson two-point correlators with wall-source/wall-sink.\n
    Evaluates all meson correlation functions defined by get_cexpr_meson_corr()
    using wall-source and wall-sink propagators. Results are saved as lat data
    with dimensions (expr_name, t_sep).\n
    Args:
        job_tag: Gauge ensemble identifier (e.g., "24D", "48I").
        traj: Trajectory number.
        get_get_prop: Callable returning propagator accessor function.
        get_psel_prob: Callable returning point selection probability object.
        get_fsel_prob: Callable returning field selection probability object.
    """
    fname = q.get_fname()
    fn = f"{job_tag}/{pname}/traj-{traj}/meson_corr.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_corr()
    expr_names = get_expr_names(cexpr)
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    t_size = total_site[3]
    get_prop = get_get_prop()
    psel_prob = get_psel_prob()
    fsel_prob = get_fsel_prob()
    psel = psel_prob.psel
    fsel = fsel_prob.fsel
    if not fsel.is_containing(psel):
        q.displayln_info(
            -1,
            "WARNING: fsel is not containing psel. The probability weighting may be wrong.",
        )
    fsel_prob[:].ravel()
    psel_prob[:].ravel()
    fsel.to_psel_local()[:]
    q.Geometry(total_site)
    #
    def load_data():
        t_t_list = q.get_mpi_chunk(
            [
                (
                    t_src,
                    t_snk,
                )
                for t_snk in range(total_site[3])
                for t_src in range(total_site[3])
            ],
            rng_state=None,
        )
        for t_src, t_snk in t_t_list:
            yield t_src, t_snk
    #
    @q.timer
    def feval(args):
        t_src, t_snk = args
        t = (t_snk - t_src) % total_site[3]
        pd = {
            "x_2": (
                "wall",
                t_snk,
            ),
            "x_1": (
                "wall",
                t_src,
            ),
            "size": total_site,
        }
        val = eval_cexpr(cexpr, positions_dict=pd, get_prop=get_prop)
        return val, t
    #
    def sum_function(val_list):
        values = np.zeros(
            (
                total_site[3],
                len(expr_names),
            ),
            dtype=complex,
        )
        for val, t in val_list:
            values[t] += val
        return q.glb_sum(values.transpose(1, 0))
    #
    auto_contractor_chunk_size = get_param(
        job_tag, "measurement", "auto_contractor_chunk_size", default=128
    )
    q.timer_fork(0)
    res_sum = q.parallel_map_sum(
        feval,
        load_data(),
        sum_function=sum_function,
        chunksize=auto_contractor_chunk_size,
    )
    q.displayln_info(f"{fname}: timer_display for parallel_map_sum")
    q.timer_display()
    q.timer_merge()
    res_sum *= 1.0 / total_site[3]
    assert q.qnorm(res_sum[0] - 1.0) < 1e-10
    ld = q.mk_lat_data(
        [
            [
                "expr_name",
                len(expr_names),
                expr_names,
            ],
            [
                "t_sep",
                t_size,
                [str(q.rel_mod(t, t_size)) for t in range(t_size)],
            ],
        ]
    )
    ld.from_numpy(res_sum)
    ld.save(get_save_path(fn))
    q.json_results_append(f"{fname}: ld sig", q.get_data_sig(ld, q.RngState()))
    for i, en in enumerate(expr_names):
        q.json_results_append(
            f"{fname}: ld '{en}' sig", q.get_data_sig(ld[i], q.RngState())
        )

@q.timer(is_timer_fork=True)
def auto_contract_meson_corr_psnk(
    job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob
):
    """
    Compute meson two-point correlators with wall-source/point-sink.\n
    Similar to auto_contract_meson_corr but uses point-sink propagators
    for the sink operator. Applies probability weighting by dividing
    by fsel_prob for each sink point. Results normalized by total volume.\n
    Args:
        job_tag: Gauge ensemble identifier.
        traj: Trajectory number.
        get_get_prop: Callable returning propagator accessor function.
        get_psel_prob: Callable returning point selection probability object.
        get_fsel_prob: Callable returning field selection probability object.
    """
    fname = q.get_fname()
    fn = f"{job_tag}/{pname}/traj-{traj}/meson_corr_psnk.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_corr(False)
    expr_names = get_expr_names(cexpr)
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    t_size = total_site[3]
    get_prop = get_get_prop()
    psel_prob = get_psel_prob()
    fsel_prob = get_fsel_prob()
    psel = psel_prob.psel
    fsel = fsel_prob.fsel
    if not fsel.is_containing(psel):
        q.displayln_info(
            -1,
            "WARNING: fsel is not containing psel. The probability weighting may be wrong.",
        )
    fsel_n_elems = fsel.n_elems
    fsel_prob_arr = fsel_prob[:].ravel()
    psel_prob[:].ravel()
    xg_fsel_arr = fsel.to_psel_local()[:]
    geo = q.Geometry(total_site)
    total_volume = geo.total_volume
    #
    def load_data():
        for t_src in range(total_site[3]):
            for idx in range(fsel_n_elems):
                yield t_src, idx
    #
    @q.timer
    def feval(args):
        t_src, idx = args
        xg_snk = tuple(xg_fsel_arr[idx])
        prob_snk = fsel_prob_arr[idx]
        t = (xg_snk[3] - t_src) % total_site[3]
        pd = {
            "x_2": (
                "point-snk",
                xg_snk,
            ),
            "x_1": (
                "wall",
                t_src,
            ),
            "size": total_site,
        }
        val = eval_cexpr(cexpr, positions_dict=pd, get_prop=get_prop)
        return val / prob_snk, t
    #
    def sum_function(val_list):
        values = np.zeros(
            (
                total_site[3],
                len(expr_names),
            ),
            dtype=complex,
        )
        for val, t in val_list:
            values[t] += val
        return values.transpose(1, 0)
    #
    auto_contractor_chunk_size = get_param(
        job_tag, "measurement", "auto_contractor_chunk_size", default=128
    )
    q.timer_fork(0)
    res_sum = q.glb_sum(
        q.parallel_map_sum(
            feval,
            load_data(),
            sum_function=sum_function,
            chunksize=auto_contractor_chunk_size,
        )
    )
    q.displayln_info(f"{fname}: timer_display for parallel_map_sum")
    q.timer_display()
    q.timer_merge()
    res_sum *= 1.0 / total_volume
    q.displayln_info(0, res_sum[0])
    ld = q.mk_lat_data(
        [
            [
                "expr_name",
                len(expr_names),
                expr_names,
            ],
            [
                "t_sep",
                t_size,
                [str(q.rel_mod(t, t_size)) for t in range(t_size)],
            ],
        ]
    )
    ld.from_numpy(res_sum)
    ld.save(get_save_path(fn))
    q.json_results_append(f"{fname}: ld sig", q.get_data_sig(ld, q.RngState()))
    for i, en in enumerate(expr_names):
        q.json_results_append(
            f"{fname}: ld '{en}' sig", q.get_data_sig(ld[i], q.RngState())
        )

@q.timer(is_timer_fork=True)
def auto_contract_meson_corr_psrc(
    job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob
):
    """
    Compute meson two-point correlators with point-source/wall-sink.\n
    Uses point-source propagators for the source operator with probability
    weighting. Results normalized by total volume.\n
    Args:
        job_tag: Gauge ensemble identifier.
        traj: Trajectory number.
        get_get_prop: Callable returning propagator accessor function.
        get_psel_prob: Callable returning point selection probability object.
        get_fsel_prob: Callable returning field selection probability object.
    """
    fname = q.get_fname()
    fn = f"{job_tag}/{pname}/traj-{traj}/meson_corr_psrc.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_corr()
    expr_names = get_expr_names(cexpr)
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    t_size = total_site[3]
    get_prop = get_get_prop()
    psel_prob = get_psel_prob()
    fsel_prob = get_fsel_prob()
    psel = psel_prob.psel
    fsel = fsel_prob.fsel
    if not fsel.is_containing(psel):
        q.displayln_info(
            -1,
            "WARNING: fsel is not containing psel. The probability weighting may be wrong.",
        )
    fsel_prob[:].ravel()
    psel_prob_arr = psel_prob[:].ravel()
    fsel.to_psel_local()[:]
    xg_psel_arr = psel[:]
    geo = q.Geometry(total_site)
    total_volume = geo.total_volume
    #
    def load_data():
        x_t_list = q.get_mpi_chunk(
            [
                (
                    pidx,
                    t_snk,
                )
                for t_snk in range(total_site[3])
                for pidx in range(len(xg_psel_arr))
            ],
            rng_state=None,
        )
        for pidx, t_snk in x_t_list:
            yield pidx, t_snk
    #
    @q.timer
    def feval(args):
        pidx, t_snk = args
        xg_src = tuple(xg_psel_arr[pidx])
        prob_src = psel_prob_arr[pidx]
        t = (xg_src[3] - t_snk) % total_site[3]
        pd = {
            "x_2": (
                "point",
                xg_src,
            ),
            "x_1": (
                "wall",
                t_snk,
            ),
            "size": total_site,
        }
        val = eval_cexpr(cexpr, positions_dict=pd, get_prop=get_prop)
        return val / prob_src, t
    #
    def sum_function(val_list):
        values = np.zeros(
            (
                total_site[3],
                len(expr_names),
            ),
            dtype=complex,
        )
        for val, t in val_list:
            values[t] += val
        return values.transpose(1, 0)
    #
    auto_contractor_chunk_size = get_param(
        job_tag, "measurement", "auto_contractor_chunk_size", default=128
    )
    q.timer_fork(0)
    res_sum = q.glb_sum(
        q.parallel_map_sum(
            feval,
            load_data(),
            sum_function=sum_function,
            chunksize=auto_contractor_chunk_size,
        )
    )
    q.displayln_info(f"{fname}: timer_display for parallel_map_sum")
    q.timer_display()
    q.timer_merge()
    res_sum *= 1.0 / total_volume
    q.displayln_info(0, res_sum[0])
    ld = q.mk_lat_data(
        [
            [
                "expr_name",
                len(expr_names),
                expr_names,
            ],
            [
                "t_sep",
                t_size,
                [str(q.rel_mod(t, t_size)) for t in range(t_size)],
            ],
        ]
    )
    ld.from_numpy(res_sum)
    ld.save(get_save_path(fn))
    q.json_results_append(f"{fname}: ld sig", q.get_data_sig(ld, q.RngState()))
    for i, en in enumerate(expr_names):
        q.json_results_append(
            f"{fname}: ld '{en}' sig", q.get_data_sig(ld[i], q.RngState())
        )

@q.timer(is_timer_fork=True)
def auto_contract_meson_corr_psnk_psrc(
    job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob
):
    """
    Compute meson two-point correlators with point-source/point-sink.\n
    Uses both point-source and point-sink propagators with probability
    weighting. Results are binned by spatial distance r for analysis of
    position-dependent correlators. Normalized by total_volume^2/t_size.\n
    Args:
        job_tag: Gauge ensemble identifier.
        traj: Trajectory number.
        get_get_prop: Callable returning propagator accessor function.
        get_psel_prob: Callable returning point selection probability object.
        get_fsel_prob: Callable returning field selection probability object.
    """
    fname = q.get_fname()
    fn = f"{job_tag}/{pname}/traj-{traj}/meson_corr_psnk_psrc.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_corr(False)
    expr_names = get_expr_names(cexpr)
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    t_size = total_site[3]
    get_prop = get_get_prop()
    psel_prob = get_psel_prob()
    fsel_prob = get_fsel_prob()
    psel = psel_prob.psel
    fsel = fsel_prob.fsel
    if not fsel.is_containing(psel):
        q.displayln_info(
            -1,
            "WARNING: fsel is not containing psel. The probability weighting may be wrong.",
        )
    fsel_prob_arr = fsel_prob[:].ravel()
    psel_prob_arr = psel_prob[:].ravel()
    xg_fsel_arr = fsel.to_psel_local()[:]
    xg_psel_arr = psel[:]
    geo = q.Geometry(total_site)
    total_volume = geo.total_volume
    r_list = get_r_list(job_tag)
    r_sq_interp_idx_coef_list = get_r_sq_interp_idx_coef_list(job_tag)
    #
    def load_data():
        for pidx in range(len(xg_psel_arr)):
            yield pidx
    #
    @q.timer
    def feval(args):
        pidx = args
        xg_src = tuple(xg_psel_arr[pidx])
        prob_src = psel_prob_arr[pidx]
        res_list = []
        for idx in range(len(xg_fsel_arr)):
            xg_snk = tuple(xg_fsel_arr[idx])
            if xg_snk == xg_src:
                prob_snk = 1.0
            else:
                prob_snk = fsel_prob_arr[idx]
            prob = prob_src * prob_snk
            x_rel = [
                q.rel_mod(xg_snk[mu] - xg_src[mu], total_site[mu]) for mu in range(4)
            ]
            r_sq = q.get_r_sq(x_rel)
            r_idx_low, r_idx_high, coef_low, coef_high = r_sq_interp_idx_coef_list[r_sq]
            t = (xg_snk[3] - xg_src[3]) % total_site[3]
            pd = {
                "x_2": (
                    "point-snk",
                    xg_snk,
                ),
                "x_1": (
                    "point",
                    xg_src,
                ),
                "size": total_site,
            }
            val = eval_cexpr(cexpr, positions_dict=pd, get_prop=get_prop)
            res_list.append((val / prob, t, r_idx_low, r_idx_high, coef_low, coef_high))
        return res_list
    #
    def sum_function(val_list):
        values = np.zeros(
            (
                total_site[3],
                len(r_list),
                len(expr_names),
            ),
            dtype=complex,
        )
        for idx, res_list in enumerate(val_list):
            for val, t, r_idx_low, r_idx_high, coef_low, coef_high in res_list:
                values[t, r_idx_low] += coef_low * val
                values[t, r_idx_high] += coef_high * val
            q.displayln_info(f"{fname}: {idx + 1}/{len(xg_psel_arr)}")
        return values.transpose(2, 0, 1)
    #
    q.timer_fork(0)
    res_sum = q.glb_sum(
        q.parallel_map_sum(feval, load_data(), sum_function=sum_function, chunksize=1)
    )
    q.displayln_info(f"{fname}: timer_display for parallel_map_sum")
    q.timer_display()
    q.timer_merge()
    res_sum *= 1.0 / (total_volume**2 / total_site[3])
    q.displayln_info(res_sum[0].sum(1))
    ld = q.mk_lat_data(
        [
            [
                "expr_name",
                len(expr_names),
                expr_names,
            ],
            [
                "t_sep",
                t_size,
                [str(q.rel_mod(t, t_size)) for t in range(t_size)],
            ],
            [
                "r",
                len(r_list),
                [f"{r:.5f}" for r in r_list],
            ],
        ]
    )
    ld.from_numpy(res_sum)
    ld.save(get_save_path(fn))
    q.json_results_append(f"{fname}: ld sig", q.get_data_sig(ld, q.RngState()))
    for i, en in enumerate(expr_names):
        q.json_results_append(
            f"{fname}: ld '{en}' sig", q.get_data_sig(ld[i], q.RngState())
        )

# ----

@q.timer
def get_cexpr_meson_jt():
    """
    Build compiled expressions for meson tensor current correlators.\n
    Computes <J_mu(x) M(t1) M^dag(t2)> where M is a meson operator
    (pi+, K+) and J_mu is a vector current (ubar*gamma_mu*u, sbar*gamma_mu*s).
    Uses wall-source/wall-sink meson operators at time slices tsep apart.\n
    Returns:
        Compiled expression object for use with eval_cexpr.
    """
    fn_base = "cache/auto_contract_cexpr/get_cexpr_meson_jt"
    #
    def calc_cexpr():
        diagram_type_dict = dict()
        diagram_type_dict[()] = "T1"
        diagram_type_dict[
            ((("t_1", "t_2"), 1), (("t_2", "x"), 1), (("x", "t_1"), 1))
        ] = "Type1"
        diagram_type_dict[
            ((("t_1", "t_2"), 1), (("t_2", "t_1"), 1), (("x", "x"), 1))
        ] = None
        diagram_type_dict[
            ((("t_1p", "t_2p"), 1), (("t_2p", "x"), 1), (("x", "t_1p"), 1))
        ] = "Type2"
        diagram_type_dict[
            ((("t_1p", "t_2p"), 1), (("t_2p", "t_1p"), 1), (("x", "x"), 1))
        ] = None
        mm_list = [
            mk_pi_p("t_2", True) * mk_pi_p("t_1") + "pi+^dag(+tsep) * pi+(-tsep)",
            mk_k_p("t_2", True) * mk_k_p("t_1") + "K+^dag(+tsep) * K+(-tsep)",
            mk_pi_p("t_2p", True) * mk_pi_p("t_1p")
            + "pi+^dag(T/2+tsep) * pi+(T/2-tsep)",
            mk_k_p("t_2p", True) * mk_k_p("t_1p") + "K+^dag(T/2+tsep) * K+(T/2-tsep)",
        ]
        op_list = [
            mk_vec_mu("u", "u", "x", 3) + "ubar_gt_u(0)",
            mk_vec_mu("s", "s", "x", 3) + "sbar_gt_s(0)",
        ]
        exprs = [
            mk_expr(1) + "1",
            op_list[0] * mm_list[0],
            op_list[0] * mm_list[1],
            -op_list[1] * mm_list[1],
            op_list[0] * mm_list[2],
            op_list[0] * mm_list[3],
            -op_list[1] * mm_list[3],
        ]
        cexpr = contract_simplify_compile(
            *exprs, is_isospin_symmetric_limit=True, diagram_type_dict=diagram_type_dict
        )
        return cexpr
    #
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

@q.timer(is_timer_fork=True)
def auto_contract_meson_jt(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob):
    """
    Compute meson tensor current correlators with point-sink.\n
    Evaluates <J_mu(x) M(t1) M^dag(t2)> using point-sink propagators.
    Results are scalar values (not position-dependent) normalized by volume.\n
    Args:
        job_tag: Gauge ensemble identifier.
        traj: Trajectory number.
        get_get_prop: Callable returning propagator accessor function.
        get_psel_prob: Callable returning point selection probability object.
        get_fsel_prob: Callable returning field selection probability object.
    """
    fname = q.get_fname()
    fn = f"{job_tag}/{pname}/traj-{traj}/meson_jt.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_jt()
    expr_names = get_expr_names(cexpr)
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    total_site[3]
    get_prop = get_get_prop()
    psel_prob = get_psel_prob()
    fsel_prob = get_fsel_prob()
    psel = psel_prob.psel
    fsel = fsel_prob.fsel
    if not fsel.is_containing(psel):
        q.displayln_info(
            -1,
            "WARNING: fsel is not containing psel. The probability weighting may be wrong.",
        )
    fsel_prob_arr = fsel_prob[:].ravel()
    psel_prob[:].ravel()
    xg_fsel_arr = fsel.to_psel_local()[:]
    psel[:]
    tsep = get_param(job_tag, "meson_tensor_tsep")
    geo = q.Geometry(total_site)
    total_volume = geo.total_volume
    #
    def load_data():
        for idx in range(len(xg_fsel_arr)):
            yield idx
    #
    @q.timer
    def feval(args):
        idx = args
        xg_snk = tuple(xg_fsel_arr[idx])
        prob_snk = fsel_prob_arr[idx]
        t = xg_snk[3]
        t_1 = (t - tsep) % total_site[3]
        t_2 = (t + tsep) % total_site[3]
        t_1p = (t_1 + total_site[3] // 2) % total_site[3]
        t_2p = (t_2 + total_site[3] // 2) % total_site[3]
        pd = {
            "x": (
                "point-snk",
                xg_snk,
            ),
            "t_1": ("wall", t_1),
            "t_2": ("wall", t_2),
            "t_1p": ("wall", t_1p),
            "t_2p": ("wall", t_2p),
            "size": total_site,
        }
        val = eval_cexpr(cexpr, positions_dict=pd, get_prop=get_prop)
        return val / prob_snk
    #
    def sum_function(val_list):
        values = np.zeros(len(expr_names), dtype=complex)
        for val in val_list:
            values += val
        return values
    #
    auto_contractor_chunk_size = get_param(
        job_tag, "measurement", "auto_contractor_chunk_size", default=128
    )
    q.timer_fork(0)
    res_sum = q.glb_sum(
        q.parallel_map_sum(
            feval,
            load_data(),
            sum_function=sum_function,
            chunksize=auto_contractor_chunk_size,
        )
    )
    q.displayln_info(f"{fname}: timer_display for parallel_map_sum")
    q.timer_display()
    q.timer_merge()
    res_sum *= 1.0 / total_volume
    q.displayln_info(0, res_sum[0])
    ld_sum = q.mk_lat_data(
        [
            [
                "expr_name",
                len(expr_names),
                expr_names,
            ],
        ]
    )
    ld_sum.from_numpy(res_sum)
    ld_sum.save(get_save_path(fn))
    q.json_results_append(f"{fname}: ld_sum sig", q.get_data_sig(ld_sum, q.RngState()))
    for i, en in enumerate(expr_names):
        q.json_results_append(
            f"{fname}: ld_sum '{en}' sig", q.get_data_sig(ld_sum[i], q.RngState())
        )

# ----

@q.timer
def get_cexpr_meson_m():
    """
    Build compiled expressions for meson mass insertion correlators.\n
    Computes <m_q(x) M(t1) M^dag(t2)> where M is a meson (pi0, pi+, K0, K+)
    and m_q is a scalar density (ubar*u, dbar*d, sbar*s). These correlators
    are needed for computing matrix elements of the axial vector current.\n
    Returns:
        Compiled expression object for use with eval_cexpr.
    """
    fn_base = "cache/auto_contract_cexpr/get_cexpr_meson_m"
    #
    def calc_cexpr():
        diagram_type_dict = dict()
        diagram_type_dict[()] = "T1"
        diagram_type_dict[
            ((("t_1", "t_2"), 1), (("t_2", "t_1"), 1), (("x", "x"), 1))
        ] = None
        diagram_type_dict[
            ((("t_1", "t_2"), 1), (("t_2", "x"), 1), (("x", "t_1"), 1))
        ] = "Type1"
        mm_list = [
            mk_pi_0("t_2", True) * mk_pi_0("t_1") + "pi0^dag(+tsep) * pi0(-tsep)",
            mk_pi_p("t_2", True) * mk_pi_p("t_1") + "pi+^dag(+tsep) * pi+(-tsep)",
            mk_k_0("t_2", True) * mk_k_0("t_1") + "K0^dag(+tsep) * K0(-tsep)",
            mk_k_p("t_2", True) * mk_k_p("t_1") + "K+^dag(+tsep) * K+(-tsep)",
        ]
        m_list = [
            mk_m("u", "x") + "ubar_u(0)",
            mk_m("d", "x") + "dbar_d(0)",
            mk_m("s", "x") + "sbar_s(0)",
        ]
        exprs = [
            mk_expr(1) + "1",
        ]
        exprs += [m * mm for mm in mm_list for m in m_list]
        cexpr = contract_simplify_compile(
            *exprs, is_isospin_symmetric_limit=True, diagram_type_dict=diagram_type_dict
        )
        return cexpr
    #
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

@q.timer(is_timer_fork=True)
def auto_contract_meson_m(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob):
    """
    Compute meson mass insertion correlators with point-sink.\n
    Evaluates <m_q(x) M(t1) M^dag(t2)> using point-sink propagators.
    Results are scalar values normalized by volume.\n
    Args:
        job_tag: Gauge ensemble identifier.
        traj: Trajectory number.
        get_get_prop: Callable returning propagator accessor function.
        get_psel_prob: Callable returning point selection probability object.
        get_fsel_prob: Callable returning field selection probability object.
    """
    fname = q.get_fname()
    fn = f"{job_tag}/{pname}/traj-{traj}/meson_m.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_m()
    expr_names = get_expr_names(cexpr)
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    total_site[3]
    get_prop = get_get_prop()
    psel_prob = get_psel_prob()
    fsel_prob = get_fsel_prob()
    psel = psel_prob.psel
    fsel = fsel_prob.fsel
    if not fsel.is_containing(psel):
        q.displayln_info(
            -1,
            "WARNING: fsel is not containing psel. The probability weighting may be wrong.",
        )
    fsel_prob_arr = fsel_prob[:].ravel()
    psel_prob[:].ravel()
    xg_fsel_arr = fsel.to_psel_local()[:]
    psel[:]
    tsep = get_param(job_tag, "meson_tensor_tsep")
    geo = q.Geometry(total_site)
    total_volume = geo.total_volume
    #
    def load_data():
        for idx in range(len(xg_fsel_arr)):
            yield idx
    #
    @q.timer
    def feval(args):
        idx = args
        xg_snk = tuple(xg_fsel_arr[idx])
        prob_snk = fsel_prob_arr[idx]
        t = xg_snk[3]
        t_2 = (t + tsep) % total_site[3]
        t_1 = (t - tsep) % total_site[3]
        pd = {
            "x": (
                "point-snk",
                xg_snk,
            ),
            "t_1": (
                "wall",
                t_1,
            ),
            "t_2": (
                "wall",
                t_2,
            ),
            "size": total_site,
        }
        val = eval_cexpr(cexpr, positions_dict=pd, get_prop=get_prop)
        return val / prob_snk
    #
    def sum_function(val_list):
        values = np.zeros(len(expr_names), dtype=complex)
        for val in val_list:
            values += val
        return values
    #
    auto_contractor_chunk_size = get_param(
        job_tag, "measurement", "auto_contractor_chunk_size", default=128
    )
    q.timer_fork(0)
    res_sum = q.glb_sum(
        q.parallel_map_sum(
            feval,
            load_data(),
            sum_function=sum_function,
            chunksize=auto_contractor_chunk_size,
        )
    )
    q.displayln_info(f"{fname}: timer_display for parallel_map_sum")
    q.timer_display()
    q.timer_merge()
    total_volume = geo.total_volume
    res_sum *= 1.0 / total_volume
    q.displayln_info(0, res_sum[0])
    ld_sum = q.mk_lat_data(
        [
            [
                "expr_name",
                len(expr_names),
                expr_names,
            ],
        ]
    )
    ld_sum.from_numpy(res_sum)
    ld_sum.save(get_save_path(fn))
    q.json_results_append(f"{fname}: ld_sum sig", q.get_data_sig(ld_sum, q.RngState()))
    for i, en in enumerate(expr_names):
        q.json_results_append(
            f"{fname}: ld_sum '{en}' sig", q.get_data_sig(ld_sum[i], q.RngState())
        )

# ----

@q.timer
def get_cexpr_meson_jj():
    """
    Build compiled expressions for meson two-current correlators.\n
    Computes <J_mu(x) J_nu(0) M(t1) M^dag(t2)> with various current
    combinations (j_mu*j_mu, j_t*j_t, j_i*j_i, x_i*x_j*j_i*j_j, etc.)
    and mass insertions (m_u*m_u, m_u*m_d, etc.). Also includes OPE
    operator product expansion terms with u, d, s quark currents.\n
    Used for computing the pi0 -> gamma gamma transition form factor.\n
    Returns:
        Compiled expression object for use with eval_cexpr.
    """
    fn_base = "cache/auto_contract_cexpr/get_cexpr_meson_jj"
    #
    def calc_cexpr():
        diagram_type_dict = dict()
        diagram_type_dict[()] = "T1"
        diagram_type_dict[
            (
                (("t_1", "x_1"), 1),
                (("t_2", "x_2"), 1),
                (("x_1", "t_1"), 1),
                (("x_2", "t_2"), 1),
            )
        ] = "Type1"
        diagram_type_dict[
            (
                (("t_1", "x_1"), 1),
                (("t_2", "x_2"), 1),
                (("x_1", "t_2"), 1),
                (("x_2", "t_1"), 1),
            )
        ] = "Type2"
        diagram_type_dict[
            (
                (("t_1", "t_2"), 1),
                (("t_2", "x_1"), 1),
                (("x_1", "x_2"), 1),
                (("x_2", "t_1"), 1),
            )
        ] = "Type3"
        diagram_type_dict[
            (
                (("t_1", "t_2"), 1),
                (("t_2", "x_1"), 1),
                (("x_1", "t_1"), 1),
                (("x_2", "x_2"), 1),
            )
        ] = None
        diagram_type_dict[
            (
                (("t_1", "t_2"), 1),
                (("t_2", "t_1"), 1),
                (("x_1", "x_1"), 1),
                (("x_2", "x_2"), 1),
            )
        ] = None
        diagram_type_dict[
            (
                (("t_1", "t_2"), 1),
                (("t_2", "t_1"), 1),
                (("x_1", "x_2"), 1),
                (("x_2", "x_1"), 1),
            )
        ] = None
        jj_list = [
            sum([mk_j_mu("x_2", mu) * mk_j_mu("x_1", mu) for mu in range(4)])
            + "j_mu(x) * j_mu(0)",
            #
            mk_j_mu("x_2", 3) * mk_j_mu("x_1", 3) + "j_t(x) * j_t(0)",
            #
            sum([mk_j_mu("x_2", mu) * mk_j_mu("x_1", mu) for mu in range(3)])
            + "j_i(x) * j_i(0)",
            #
            sum(
                [
                    mk_fac(f"rel_mod_sym(x_2[1][{mu}] - x_1[1][{mu}], size[{mu}])")
                    * mk_fac(f"rel_mod_sym(x_2[1][{nu}] - x_1[1][{nu}], size[{nu}])")
                    * mk_j_mu("x_2", mu)
                    * mk_j_mu("x_1", nu)
                    for mu in range(3)
                    for nu in range(3)
                ]
            )
            + "x[i] * x[j] * j_i(x) * j_j(0)",
            #
            sum(
                [
                    mk_fac(f"rel_mod_sym(x_2[1][{mu}] - x_1[1][{mu}], size[{mu}])")
                    * mk_j_mu("x_2", mu)
                    * mk_j_mu("x_1", 3)
                    for mu in range(3)
                ]
            )
            + "x[i] * j_i(x) * j_t(0)",
            #
            sum(
                [
                    mk_fac(f"rel_mod_sym(x_1[1][{mu}] - x_2[1][{mu}], size[{mu}])")
                    * mk_j_mu("x_1", mu)
                    * mk_j_mu("x_2", 3)
                    for mu in range(3)
                ]
            )
            + "-x[i] * j_i(-x) * j_t(0)",
        ]
        assert len(jj_list) == 6
        m2_list = [
            mk_m("u", "x_2") + "ubar_u(x)",
            mk_m("d", "x_2") + "dbar_d(x)",
            mk_m("s", "x_2") + "sbar_s(x)",
        ]
        assert len(m2_list) == 3
        m1_list = [
            mk_m("u", "x_1") + "ubar_u(0)",
            mk_m("d", "x_1") + "dbar_d(0)",
            mk_m("s", "x_1") + "sbar_s(0)",
        ]
        assert len(m1_list) == 3
        m1m2_list = [m2 * m1 for m2 in m2_list for m1 in m1_list]
        assert len(m1m2_list) == 9
        op_list = jj_list + m1m2_list
        assert len(op_list) == 15
        mm_list = [
            mk_pi_0("t_2", True) * mk_pi_0("t_1") + "pi0^dag(x[t]+tsep) * pi0(-tsep)",
            #
            mk_sym(1)
            / 2
            * (
                mk_pi_p("t_2", True) * mk_pi_p("t_1")
                + mk_pi_m("t_2", True) * mk_pi_m("t_1")
            )
            + "pi+^dag(x[t]+tsep) * pi+(-tsep)",
            #
            mk_sym(1)
            / 2
            * (
                mk_k_0("t_2", True) * mk_k_0("t_1")
                + mk_k_0_bar("t_2", True) * mk_k_0_bar("t_1")
            )
            + "K0^dag(x[t]+tsep) * K0(-tsep)",
            #
            mk_sym(1)
            / 2
            * (
                mk_k_p("t_2", True) * mk_k_p("t_1")
                + mk_k_m("t_2", True) * mk_k_m("t_1")
            )
            + "K+^dag(x[t]+tsep) * K+(-tsep)",
        ]
        assert len(mm_list) == 4
        exprs_self_energy = [op * mm for mm in mm_list for op in op_list]
        assert len(exprs_self_energy) == 60
        #
        op_u_ope_list = [
            sum(
                [
                    mk_vec_mu("u", "u'", "x_2", mu) * mk_vec_mu("u'", "u", "x_1", mu)
                    for mu in range(4)
                ]
            )
            + "juu'_mu(x) * ju'u_mu(0)",
            #
            mk_vec_mu("u", "u'", "x_2", 3) * mk_vec_mu("u'", "u", "x_1", 3)
            + "juu'_t(x) * ju'u_t(0)",
            #
            sum(
                [
                    mk_vec_mu("u", "u'", "x_2", mu) * mk_vec_mu("u'", "u", "x_1", mu)
                    for mu in range(3)
                ]
            )
            + "juu'_i(x) * ju'u_i(0)",
        ]
        assert len(op_u_ope_list) == 3
        op_d_ope_list = [
            sum(
                [
                    mk_vec_mu("d", "d'", "x_2", mu) * mk_vec_mu("d'", "d", "x_1", mu)
                    for mu in range(4)
                ]
            )
            + "jdd'_mu(x) * jd'd_mu(0)",
            #
            mk_vec_mu("d", "d'", "x_2", 3) * mk_vec_mu("d'", "d", "x_1", 3)
            + "jdd'_t(x) * jd'd_t(0)",
            #
            sum(
                [
                    mk_vec_mu("d", "d'", "x_2", mu) * mk_vec_mu("d'", "d", "x_1", mu)
                    for mu in range(3)
                ]
            )
            + "jdd'_i(x) * jd'd_i(0)",
        ]
        assert len(op_d_ope_list) == 3
        op_s_ope_list = [
            sum(
                [
                    mk_vec_mu("s", "s'", "x_2", mu) * mk_vec_mu("s'", "s", "x_1", mu)
                    for mu in range(4)
                ]
            )
            + "jss'_mu(x) * js's_mu(0)",
            #
            mk_vec_mu("s", "s'", "x_2", 3) * mk_vec_mu("s'", "s", "x_1", 3)
            + "jss'_t(x) * js's_t(0)",
            #
            sum(
                [
                    mk_vec_mu("s", "s'", "x_2", mu) * mk_vec_mu("s'", "s", "x_1", mu)
                    for mu in range(3)
                ]
            )
            + "jss'_i(x) * js's_i(0)",
        ]
        assert len(op_s_ope_list) == 3
        mm_pi_ope_list = [
            mk_pi_p("t_2", True) * mk_pi_p("t_1") + "pi+^dag(x[t]+tsep) * pi+(-tsep)",
        ]
        assert len(mm_pi_ope_list) == 1
        mm_kk_ope_list = [
            mk_k_p("t_2", True) * mk_k_p("t_1") + "K+^dag(x[t]+tsep) * K+(-tsep)",
        ]
        assert len(mm_kk_ope_list) == 1
        exprs_ope = []
        exprs_ope += [op * mm for mm in mm_pi_ope_list for op in op_u_ope_list]
        exprs_ope += [op * mm for mm in mm_pi_ope_list for op in op_d_ope_list]
        exprs_ope += [op * mm for mm in mm_kk_ope_list for op in op_u_ope_list]
        exprs_ope += [op * mm for mm in mm_kk_ope_list for op in op_s_ope_list]
        assert len(exprs_ope) == 12
        #
        exprs = [
            mk_expr(1) + "1",
        ]
        exprs += exprs_self_energy
        exprs += exprs_ope
        assert len(exprs) == 1 + 60 + 12
        #
        cexpr = contract_simplify_compile(
            *exprs, is_isospin_symmetric_limit=True, diagram_type_dict=diagram_type_dict
        )
        return cexpr
    #
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

@q.timer(is_timer_fork=True)
def auto_contract_meson_jj(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob):
    """
    Compute meson two-current correlators with point-source/point-sink.\n
    Evaluates <J(x) J(0) M(t1) M^dag(t2)> using point-source and point-sink
    propagators. Results are binned by spatial separation r and time
    difference t. Uses adaptive probability weighting.\n
    Args:
        job_tag: Gauge ensemble identifier.
        traj: Trajectory number.
        get_get_prop: Callable returning propagator accessor function.
        get_psel_prob: Callable returning point selection probability object.
        get_fsel_prob: Callable returning field selection probability object.
    """
    fname = q.get_fname()
    fn = f"{job_tag}/{pname}/traj-{traj}/meson_jj.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_jj()
    expr_names = get_expr_names(cexpr)
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    t_size = total_site[3]
    get_prop = get_get_prop()
    psel_prob = get_psel_prob()
    fsel_prob = get_fsel_prob()
    psel = psel_prob.psel
    fsel = fsel_prob.fsel
    if not fsel.is_containing(psel):
        q.displayln_info(
            -1,
            "WARNING: fsel is not containing psel. The probability weighting may be wrong.",
        )
    fsel_prob_arr = fsel_prob[:].ravel()
    psel_prob_arr = psel_prob[:].ravel()
    xg_fsel_arr = fsel.to_psel_local()[:]
    xg_psel_arr = psel[:]
    tsep = get_param(job_tag, "meson_tensor_tsep")
    geo = q.Geometry(total_site)
    total_volume = geo.total_volume
    r_list = get_r_list(job_tag)
    r_sq_interp_idx_coef_list = get_r_sq_interp_idx_coef_list(job_tag)
    #
    def load_data():
        for pidx in range(len(xg_psel_arr)):
            yield pidx
    #
    @q.timer
    def feval(args):
        pidx = args
        xg_src = tuple(xg_psel_arr[pidx])
        prob_src = psel_prob_arr[pidx]
        res_list = []
        for idx in range(len(xg_fsel_arr)):
            xg_snk = tuple(xg_fsel_arr[idx])
            if xg_snk == xg_src:
                prob_snk = 1.0
            else:
                prob_snk = fsel_prob_arr[idx]
            prob = prob_src * prob_snk
            x_rel = [
                q.rel_mod(xg_snk[mu] - xg_src[mu], total_site[mu]) for mu in range(4)
            ]
            r_sq = q.get_r_sq(x_rel)
            r_idx_low, r_idx_high, coef_low, coef_high = r_sq_interp_idx_coef_list[r_sq]
            x_rel_t = x_rel[3]
            x_2_t = xg_src[3]
            x_1_t = x_2_t + x_rel_t
            t_2 = (max(x_1_t, x_2_t) + tsep) % total_site[3]
            t_1 = (min(x_1_t, x_2_t) - tsep) % total_site[3]
            pd = {
                "x_2": (
                    "point-snk",
                    xg_snk,
                ),
                "x_1": (
                    "point",
                    xg_src,
                ),
                "t_2": ("wall", t_2),
                "t_1": ("wall", t_1),
                "size": total_site,
            }
            t = x_rel_t % t_size
            val = eval_cexpr(cexpr, positions_dict=pd, get_prop=get_prop)
            res_list.append(
                (
                    val / prob,
                    t,
                    r_idx_low,
                    r_idx_high,
                    coef_low,
                    coef_high,
                )
            )
        return res_list
    #
    def sum_function(val_list):
        values = np.zeros(
            (
                t_size,
                len(r_list),
                len(expr_names),
            ),
            dtype=complex,
        )
        for idx, res_list in enumerate(val_list):
            for val, t, r_idx_low, r_idx_high, coef_low, coef_high in res_list:
                values[t, r_idx_low] += coef_low * val
                values[t, r_idx_high] += coef_high * val
            q.displayln_info(f"{fname}: {idx + 1}/{len(xg_psel_arr)}")
        return q.glb_sum(values.transpose(2, 0, 1))
    #
    q.timer_fork(0)
    res_sum = q.parallel_map_sum(
        feval, load_data(), sum_function=sum_function, chunksize=1
    )
    q.displayln_info(f"{fname}: timer_display for parallel_map_sum")
    q.timer_display()
    q.timer_merge()
    res_sum *= 1.0 / total_volume
    ld_sum = q.mk_lat_data(
        [
            [
                "expr_name",
                len(expr_names),
                expr_names,
            ],
            [
                "t",
                t_size,
                [str(q.rel_mod(t, t_size)) for t in range(t_size)],
            ],
            [
                "r",
                len(r_list),
                [f"{r:.5f}" for r in r_list],
            ],
        ]
    )
    ld_sum.from_numpy(res_sum)
    ld_sum.save(get_save_path(fn))
    q.json_results_append(f"{fname}: ld_sum sig", q.get_data_sig(ld_sum, q.RngState()))
    for i, en in enumerate(expr_names):
        q.json_results_append(
            f"{fname}: ld_sum '{en}' sig", q.get_data_sig(ld_sum[i], q.RngState())
        )

# ----

@q.timer
def get_cexpr_pi0_jjp():
    """
    Build compiled expressions for pi0 -> J J parity-violating correlators.\n
    Computes <epsilon_{ijk} x_i J_j(x) J_k(0) pi0(w)> where J is the
    electromagnetic current and pi0 is represented by either:
    - pp: pseudoscalar density (ubar*gamma5*u - dbar*gamma5*d)
    - j5: time-component of axial vector current\n
    Uses diagram types Type1 (connected) and Type3 (disconnected).\n
    Returns:
        Compiled expression object for use with eval_cexpr.
    """
    fn_base = "cache/auto_contract_cexpr/get_cexpr_pi0_jjp"
    #
    def calc_cexpr():
        diagram_type_dict = dict()
        diagram_type_dict[()] = "T1"
        diagram_type_dict[
            ((("w", "x_1"), 1), (("x_1", "x_2"), 1), (("x_2", "w"), 1))
        ] = "Type1"
        diagram_type_dict[
            ((("w", "x_1"), 1), (("x_1", "w"), 1), (("x_2", "x_2"), 1))
        ] = "Type2"
        diagram_type_dict[
            ((("w", "w"), 1), (("x_1", "x_2"), 1), (("x_2", "x_1"), 1))
        ] = "Type3"
        diagram_type_dict[
            ((("w", "w"), 1), (("x_1", "x_1"), 1), (("x_2", "x_2"), 1))
        ] = "Type4"
        jj_list = [
            sum(
                [
                    q.epsilon_tensor(mu, nu, rho)
                    * mk_fac(f"rel_mod_sym(x_2[1][{mu}] - x_1[1][{mu}], size[{mu}])")
                    * mk_j_mu("x_2", nu)
                    * mk_j_mu("x_1", rho)
                    for mu in range(3)
                    for nu in range(3)
                    for rho in range(3)
                ]
            )
            + "e(i,j,k) * x[i] * j_j(x) * j_k(0)",
        ]
        assert len(jj_list) == 1
        pp_u = mk_meson("u", "u", "w")
        pp_d = mk_meson("d", "d", "w")
        j5_u = (mk_fac(sympy.I) + "i") * mk_vec5_mu("u", "u", "w", 3)
        j5_d = (mk_fac(sympy.I) + "i") * mk_vec5_mu("d", "d", "w", 3)
        p_list = [
            pp_u,
            pp_d,
            j5_u,
            j5_d,
        ]
        assert len(p_list) == 4
        exprs_1_list = [
            mk_expr(1) + "1",
        ]
        exprs_corr_list = [jj * p for p in p_list for jj in jj_list]
        exprs_corr_list = [
            (
                expr,
                # None,
                "Type1",
                # "Type2",
                "Type3",
                # "Type4",
            )
            for expr in exprs_corr_list
        ]
        exprs = exprs_1_list + exprs_corr_list
        assert len(exprs) == 1 + 4
        cexpr = contract_simplify_compile(
            *exprs, is_isospin_symmetric_limit=True, diagram_type_dict=diagram_type_dict
        )
        return cexpr
    #
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

@q.timer(is_timer_fork=True)
def auto_contract_pi0_jjp(job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob):
    """
    Compute pi0 -> J J parity-violating correlators with adaptive sampling.\n
    Evaluates <epsilon_{ijk} x_i J_j(x) J_k(0) pi0(w)> using point-source
    propagators for all three points (x, 0, w). Uses adaptive importance
    sampling with threshold-based probability weighting to reduce variance.
    Results are binned by (t1, t2, r) where t1, t2 are time separations
    from the pi0 operator.\n
    Args:
        job_tag: Gauge ensemble identifier.
        traj: Trajectory number.
        get_get_prop: Callable returning propagator accessor function.
        get_psel_prob: Callable returning point selection probability object.
        get_fsel_prob: Callable returning field selection probability object.
    """
    fname = q.get_fname()
    fn = f"{job_tag}/{pname}/traj-{traj}/pi0_jjp.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_pi0_jjp()
    expr_names = get_expr_names(cexpr)
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    t_size = total_site[3]
    load_point_distribution(job_tag)
    np.array(total_site.to_list())
    get_prop = get_get_prop()
    psel_prob = get_psel_prob()
    fsel_prob = get_fsel_prob()
    psel = psel_prob.psel
    fsel = fsel_prob.fsel
    if not fsel.is_containing(psel):
        q.displayln_info(
            -1,
            "WARNING: fsel is not containing psel. The probability weighting may be wrong.",
        )
    fsel_prob_arr = fsel_prob[:].ravel()
    psel_prob_arr = psel_prob[:].ravel()
    fsel_prob_inv_avg = np.average(1.0 / fsel_prob_arr)
    psel_prob_inv_avg = np.average(1.0 / psel_prob_arr)
    xg_fsel_arr = fsel.to_psel_local()[:]
    xg_psel_arr = psel[:]
    n_points = len(xg_psel_arr)
    get_param(job_tag, "meson_tensor_tsep")
    geo = q.Geometry(total_site)
    total_volume = geo.total_volume
    r_list = get_r_list(job_tag)
    r_sq_interp_idx_coef_list = get_r_sq_interp_idx_coef_list(job_tag)
    n_elems = len(xg_fsel_arr)
    n_pairs = n_points * n_points
    total_site_arr = np.array(total_site.to_list())
    total_site_arr = np.broadcast_to(
        total_site_arr,
        (
            n_elems,
            4,
        ),
    )
    #
    threshold = get_param(job_tag, "pi0_jjp_threshold")
    u_rand_prob = q.SelectedFieldRealD(fsel, 1)
    u_rand_prob.set_rand(
        q.RngState(f"auto_contract_pi0_jjp,{get_job_seed(job_tag)},{traj}"), 1.0, 0.0
    )
    u_rand_prob_arr = u_rand_prob[:].ravel()
    fn_meson_corr = f"{job_tag}/{pname}/traj-{traj}/meson_corr.lat"
    if get_load_path(fn_meson_corr) is None:
        q.displayln_info(f"{fname}: '{fn_meson_corr}' does not exist. Skipping.")
        return
    ld_meson_corr = q.load_lat_data(get_load_path(fn_meson_corr))
    expr_idx = ld_meson_corr.dim_idx(0, "< pi+^dag(0) * pi+(-tsep) >  exprs")
    meson_corr_arr = ld_meson_corr[expr_idx]
    meson_corr_arr = meson_corr_arr / meson_corr_arr[0]
    #
    def get_prop_norm_sqrt(*args):
        is_sloppy = True
        return abs(ama_extract(get_prop(*args, is_norm_sqrt=True), is_sloppy=is_sloppy))
    #
    def load_psrc_psrc_prop_norm_sqrt(flavor, i):
        xg1_src = tuple(xg_psel_arr[i])
        x_1 = (
            "point",
            xg1_src,
        )
        v_list = []
        for j in range(n_points):
            xg2_src = tuple(xg_psel_arr[j])
            x_2 = (
                "point",
                xg2_src,
            )
            v = get_prop_norm_sqrt(flavor, x_1, x_2)
            v_list.append(v)
        return np.array(v_list)
    #
    psrc_psrc_prop_norm_sqrt = np.array(
        [
            q.parallel_map(
                lambda i: load_psrc_psrc_prop_norm_sqrt(flavor, i),
                range(n_points),
                verbose=1,
            )
            for flavor in [
                "l",
                "s",
            ]
        ],
        dtype=float,
    )
    #
    def load_psrc_psnk_prop_norm_sqrt(flavor, i):
        xg1_src = tuple(xg_psel_arr[i])
        x_1 = (
            "point",
            xg1_src,
        )
        v_list = []
        for j in range(n_elems):
            xg2_snk = tuple(xg_fsel_arr[j])
            x_2 = (
                "point-snk",
                xg2_snk,
            )
            v = get_prop_norm_sqrt(flavor, x_2, x_1)
            v_list.append(v)
        return np.array(v_list)
    #
    psrc_psnk_prop_norm_sqrt = np.array(
        [
            q.parallel_map(
                lambda i: load_psrc_psnk_prop_norm_sqrt(flavor, i),
                range(n_points),
                verbose=1,
            )
            for flavor in [
                "l",
                "s",
            ]
        ],
        dtype=float,
    )
    t_psel_arr = xg_psel_arr[:, 3]
    t_fsel_arr = xg_fsel_arr[:, 3]
    #
    def get_estimate(idx_w, idx_1, idx_2):
        flavor_l = 0
        prob_1 = psel_prob_arr[idx_1] * psel_prob_inv_avg
        prob_w = psel_prob_arr[idx_w] * psel_prob_inv_avg
        prob_2 = fsel_prob_arr[idx_2] * fsel_prob_inv_avg
        xg_1_t = t_psel_arr[idx_1]
        xg_w_t = t_psel_arr[idx_w]
        xg_2_t = t_fsel_arr[idx_2]
        corr1 = np.abs(meson_corr_arr[(xg_w_t - xg_1_t) % t_size])
        corr2 = np.abs(meson_corr_arr[(xg_w_t - xg_2_t) % t_size])
        p1p2 = psrc_psnk_prop_norm_sqrt[flavor_l, idx_1, idx_2]
        wp1 = psrc_psrc_prop_norm_sqrt[flavor_l, idx_1, idx_w]
        wp2 = psrc_psnk_prop_norm_sqrt[flavor_l, idx_w, idx_2]
        value = 0
        value += p1p2 * wp1 * wp2 / np.sqrt(corr1 * corr2)
        value /= prob_1 * prob_w * prob_2
        assert np.all(value > 0)
        return value
    #
    @q.timer
    def get_weight(idx_w, idx_1, idx_2):
        """
        return weight for point ``idx_2`` (1 / prob or zero)
        """
        est = get_estimate(idx_w, idx_1, idx_2)
        assert est.shape == idx_2.shape
        prob = est / threshold
        assert prob.shape == idx_2.shape
        rand = u_rand_prob_arr[idx_2]
        assert rand.shape == idx_2.shape
        weight = 1.0 / prob
        weight[prob >= 1] = 1.0
        weight[rand >= prob] = 0.0
        return weight
    #
    def load_data():
        idx_pair = 0
        for idx_1 in range(n_points):
            for idx_w in range(n_points):
                idx_pair += 1
                yield idx_1, idx_w
    #
    @q.timer
    def feval(args):
        idx_1, idx_w = args
        idx_2_arr = np.arange(n_elems)
        xg_1 = tuple(xg_psel_arr[idx_1])
        xg_w = tuple(xg_psel_arr[idx_w])
        #
        # Adaptive sampling method:
        prob_1 = psel_prob_arr[idx_1]
        prob_w = psel_prob_arr[idx_w]
        if idx_1 != idx_w:
            prob = total_volume * prob_1 * prob_w
        else:
            prob = total_volume * prob_1
        #
        weight_base = 1.0 / prob / total_volume
        xg_1_arr = np.broadcast_to(np.array(xg_1), total_site_arr.shape)
        xg_2_arr = xg_fsel_arr[idx_2_arr]
        xg_w_arr = np.broadcast_to(np.array(xg_w), total_site_arr.shape)
        t_size_arr = total_site_arr[:, 3]
        xg_1_t_arr = xg_1_arr[:, 3]
        xg_2_t_arr = xg_2_arr[:, 3]
        xg_w_t_arr = xg_w_arr[:, 3]
        x_rel_arr = q.rel_mod_arr(xg_2_arr - xg_1_arr, total_site_arr)
        r_sq_arr = q.sqr(x_rel_arr[:, :3]).sum(-1)
        xg_1_xg_t_arr = q.rel_mod_arr(xg_1_t_arr - xg_w_t_arr, t_size_arr)
        xg_2_xg_t_arr = q.rel_mod_arr(xg_2_t_arr - xg_w_t_arr, t_size_arr)
        weight_arr = weight_base * get_weight(idx_w, idx_1, idx_2_arr)
        weight_arr[np.abs(xg_2_xg_t_arr - xg_1_xg_t_arr) >= t_size_arr // 2] = 0.0
        results = []
        for idx_2 in idx_2_arr[weight_arr > 0]:
            xg_2 = tuple(xg_fsel_arr[idx_2])
            weight = weight_arr[idx_2]
            #
            # Adaptive sampling method:
            if xg_2 != xg_1 and xg_2 != xg_w:
                prob_2 = fsel_prob_arr[idx_2]
            else:
                prob_2 = 1
            #
            r_sq = r_sq_arr[idx_2]
            xg_1_xg_t = xg_1_xg_t_arr[idx_2]
            xg_2_xg_t = xg_2_xg_t_arr[idx_2]
            pd = {
                "w": (
                    "point",
                    xg_w,
                ),
                "x_1": (
                    "point",
                    xg_1,
                ),
                "x_2": (
                    "point-snk",
                    xg_2,
                ),
                "size": total_site.to_list(),
            }
            t_1 = xg_1_xg_t
            t_2 = xg_2_xg_t
            val = eval_cexpr(cexpr, positions_dict=pd, get_prop=get_prop)
            r_idx_low, r_idx_high, coef_low, coef_high = r_sq_interp_idx_coef_list[r_sq]
            results.append(
                (
                    weight * val / prob_2,
                    t_1,
                    t_2,
                    r_idx_low,
                    r_idx_high,
                    coef_low,
                    coef_high,
                )
            )
        return idx_1, idx_w, results
    #
    def sum_function(val_list):
        n_total = 0
        n_selected = 0
        idx_pair = 0
        values = np.zeros(
            (
                t_size,
                t_size,
                len(r_list),
                len(expr_names),
            ),
            dtype=complex,
        )
        for idx_1, idx_w, results in val_list:
            idx_pair += 1
            n_total += n_elems
            for val, t_1, t_2, r_idx_low, r_idx_high, coef_low, coef_high in results:
                n_selected += 1
                values[t_1, t_2, r_idx_low] += coef_low * val
                values[t_1, t_2, r_idx_high] += coef_high * val
            if idx_pair % (n_pairs // 1000 + 100) == 0:
                xg_1_src = tuple(xg_psel_arr[idx_1])
                xg_w_src = tuple(xg_psel_arr[idx_w])
                q.displayln_info(
                    1,
                    f"{fname}: {idx_pair}/{n_pairs} {xg_1_src} {xg_w_src} {len(results)}/{n_elems} n_total={n_total} n_selected={n_selected} ratio={n_selected / n_total}",
                )
        q.displayln_info(
            1,
            f"{fname}: Final: n_total={n_total} n_selected={n_selected} ratio={n_selected / n_total}",
        )
        n_total = q.glb_sum(n_total)
        n_selected = q.glb_sum(n_selected)
        q.displayln_info(
            1,
            f"{fname}: Final(glb_sum): n_total={n_total} n_selected={n_selected} ratio={n_selected / n_total}",
        )
        return values.transpose(3, 0, 1, 2)
    #
    q.timer_fork(0)
    res_sum = q.parallel_map_sum(
        feval, load_data(), sum_function=sum_function, chunksize=1
    )
    res_sum = q.glb_sum(res_sum)
    q.displayln_info(f"{fname}: timer_display")
    q.timer_display()
    q.timer_merge()
    res_sum *= 1.0 / (total_volume / t_size)
    ld_sum = q.mk_lat_data(
        [
            [
                "expr_name",
                len(expr_names),
                expr_names,
            ],
            [
                "t1",
                t_size,
                [str(q.rel_mod(t, t_size)) for t in range(t_size)],
            ],
            [
                "t2",
                t_size,
                [str(q.rel_mod(t, t_size)) for t in range(t_size)],
            ],
            [
                "r",
                len(r_list),
                [f"{r:.5f}" for r in r_list],
            ],
        ]
    )
    ld_sum.from_numpy(res_sum)
    ld_sum.save(get_save_path(fn))
    q.json_results_append(f"{fname}: ld_sum sig", q.get_data_sig(ld_sum, q.RngState()))
    for i, en in enumerate(expr_names):
        q.json_results_append(
            f"{fname}: ld_sum '{en}' sig", q.get_data_sig(ld_sum[i], q.RngState())
        )

### ------

@q.timer
def get_cexpr_tadpole_current():
    """
    Build compiled expressions for tadpole (disconnected) current correlators.\n
    Computes <J_mu(x_1)> with various quark flavor combinations (ls, sl, ll, ss).
    Results must be multiplied by (m_s - m_l) and summed over x_2 (space-time
    volume) to get the full tadpole contribution. The 1/3 charge factor is
    already included.\n
    Returns:
        Compiled expression object for use with eval_cexpr.
    """
    """
    Intend to calculating < J_mu(x_1) >
    !!! Results needs to be multiplied by (m_s - m_l) and sum over "x_2" over the entire space-time volume !!!
    After summing over "x_2", the first and second expr should be same (and so is the following exprs). (We can average these two.)
    1/3 charge factor is already multiplied.
    """
    fn_base = "cache/auto_contract_cexpr/get_cexpr_tadpole_current"
    #
    def calc_cexpr():
        diagram_type_dict = dict()
        diagram_type_dict[()] = "T1"
        diagram_type_dict[((("x_1", "x_2"), 1), (("x_2", "x_1"), 1))] = "TypeC"
        diagram_type_dict[((("x_1", "x_1"), 1), (("x_2", "x_2"), 1))] = None
        #
        jj_list = []
        jj_list += [
            1 / 3 * mk_vec_mu("l", "s", "x_1", mu) * mk_scalar("s", "l", "x_2")
            for mu in range(4)
        ]
        jj_list += [
            1 / 3 * mk_vec_mu("s", "l", "x_1", mu) * mk_scalar("l", "s", "x_2")
            for mu in range(4)
        ]
        jj_list += [
            1 / 3 * mk_vec_mu("l", "l", "x_1", mu) * mk_scalar("l", "l", "x_2")
            for mu in range(4)
        ]
        jj_list += [
            1 / 3 * mk_vec_mu("s", "s", "x_1", mu) * mk_scalar("s", "s", "x_2")
            for mu in range(4)
        ]
        #
        exprs = [
            mk_expr(1) + "1",
        ]
        exprs += jj_list
        #
        cexpr = contract_simplify_compile(
            *exprs, is_isospin_symmetric_limit=True, diagram_type_dict=diagram_type_dict
        )
        return cexpr
        #
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

@q.timer
def get_cexpr_pi0_current():
    """
    Build compiled expressions for pi0 current correlators.\n
    Computes <J_mu(x_1) pi0(t_1)> where J is the electromagnetic current
    and pi0 is represented by pseudoscalar density or axial vector current.\n
    Returns:
        Compiled expression object for use with eval_cexpr.
    """
    """
    Intend to calculating < J_mu(x_1) pi0(t_1) >
    """
    fn_base = "cache/auto_contract_cexpr/get_cexpr_pi0_current"
    #
    def calc_cexpr():
        diagram_type_dict = dict()
        diagram_type_dict[()] = "T1"
        diagram_type_dict[((("w", "x_1"), 1), (("x_1", "w"), 1))] = "TypeC"
        #
        jj_d_list = [mk_j_mu("x_1", mu) for mu in range(4)]
        assert len(jj_d_list) == 4
        #
        pp_u = mk_meson("u", "u", "w")
        pp_d = mk_meson("d", "d", "w")
        j5_u = (mk_fac(sympy.I) + "i") * mk_vec5_mu("u", "u", "w", 3)
        j5_d = (mk_fac(sympy.I) + "i") * mk_vec5_mu("d", "d", "w", 3)
        sqrt_2_inv = mk_fac(1 / sympy.sqrt(2)) + "1/sqrt(2)"
        pp_pi = sqrt_2_inv * (pp_u - pp_d)
        j5_pi = sqrt_2_inv * (j5_u - j5_d)
        pi0d_list = [
            pp_pi,
            j5_pi,
        ]
        assert len(pi0d_list) == 2
        #
        exprs_list_pi0_decay = [jj_d * pi0d for pi0d in pi0d_list for jj_d in jj_d_list]
        assert len(exprs_list_pi0_decay) == 2 * 4
        #
        exprs = [
            mk_expr(1) + "1",
        ]
        exprs += exprs_list_pi0_decay
        assert len(exprs) == 1 + 4 * 2
        #
        cexpr = contract_simplify_compile(
            *exprs, is_isospin_symmetric_limit=True, diagram_type_dict=diagram_type_dict
        )
        return cexpr
        #
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

@q.timer(is_timer_fork=True)
def auto_contract_tadpole_current(
    job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob
):
    """
    Compute tadpole (disconnected) current correlators.\n
    Evaluates <J_mu(x_1)> summed over x_2 for each field-selection point.
    Results are stored as a SelectedField with shape (n_elems, n_exprs).
    The first expression is the normalization factor.\n
    Args:
        job_tag: Gauge ensemble identifier.
        traj: Trajectory number.
        get_get_prop: Callable returning propagator accessor function.
        get_psel_prob: Callable returning point selection probability object.
        get_fsel_prob: Callable returning field selection probability object.
    """
    fname = q.get_fname()
    fn = f"{job_tag}/{pname}/traj-{traj}/tadpole-current.sfield"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_tadpole_current()
    expr_names = get_expr_names(cexpr)
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    total_site[3]
    get_prop = get_get_prop()
    psel_prob = get_psel_prob()
    fsel_prob = get_fsel_prob()
    psel = psel_prob.psel
    fsel = fsel_prob.fsel
    if not fsel.is_containing(psel):
        q.displayln_info(
            -1,
            "WARNING: fsel is not containing psel. The probability weighting may be wrong.",
        )
    fsel_prob_arr = fsel_prob[:].ravel()
    psel_prob_arr = psel_prob[:].ravel()
    xg_fsel_arr = fsel.to_psel_local()[:]
    xg_psel_arr = psel[:]
    q.Geometry(total_site)
    sf_tadpole_current = q.SelectedFieldComplexD(fsel, len(expr_names))
    q.set_zero(sf_tadpole_current)
    sf_tadpole_current_arr = sf_tadpole_current[:]
    #
    def load_data():
        for fidx in range(len(xg_fsel_arr)):
            yield fidx
    #
    @q.timer
    def feval(args):
        fidx = args
        xg_snk = tuple(xg_fsel_arr[fidx])
        prob_snk = fsel_prob_arr[fidx]
        val_sum = 0
        for pidx in range(len(xg_psel_arr)):
            xg_src = tuple(xg_psel_arr[pidx])
            prob_src = psel_prob_arr[pidx]
            if xg_snk == xg_src:
                prob_src = prob_src / prob_snk
            prob = prob_src * prob_snk
            pd = {
                "x_1": (
                    "point-snk",
                    xg_snk,
                ),
                "x_2": (
                    "point",
                    xg_src,
                ),
                "size": total_site,
            }
            val = eval_cexpr(cexpr, positions_dict=pd, get_prop=get_prop)
            val_sum = val_sum + val / prob
        res = (
            fidx,
            val_sum,
        )
        return res
    #
    def sum_function(val_list):
        for idx, res in enumerate(val_list):
            (
                fidx,
                val_sum,
            ) = res
            sf_tadpole_current_arr[fidx] = val_sum
            if (idx + 1) % (len(xg_fsel_arr) // 128 + 16) == 0:
                q.displayln_info(f"{fname}: {idx + 1}/{len(xg_fsel_arr)}")
        return None
    #
    auto_contractor_chunk_size = get_param(
        job_tag, "measurement", "auto_contractor_chunk_size", default=128
    )
    q.timer_fork(0)
    q.parallel_map_sum(
        feval,
        load_data(),
        sum_function=sum_function,
        chunksize=auto_contractor_chunk_size,
    )
    q.displayln_info(f"{fname}: timer_display for parallel_map_sum")
    q.timer_display()
    q.timer_merge()
    sf_tadpole_current.save_double(get_save_path(fn))
    q.json_results_append(
        f"{fname}: sf_tadpole_current sig",
        q.get_data_sig(sf_tadpole_current, q.RngState()),
    )
    for i, en in enumerate(expr_names):
        q.json_results_append(
            f"{fname}: sf_tadpole_current '{en}' sig",
            q.glb_sum(q.get_data_sig(sf_tadpole_current[:, i], q.RngState())),
        )

@q.timer(is_timer_fork=True)
def auto_contract_pi0_current(
    job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob
):
    """
    Compute pi0 current correlators <J_mu(x_1) pi0(t_1)>.\n
    Evaluates the correlator for each source time slice t_src, storing
    results as separate SelectedFields. Output is saved as a multi-field
    archive file indexed by t_src.\n
    Args:
        job_tag: Gauge ensemble identifier.
        traj: Trajectory number.
        get_get_prop: Callable returning propagator accessor function.
        get_psel_prob: Callable returning point selection probability object.
        get_fsel_prob: Callable returning field selection probability object.
    """
    fname = q.get_fname()
    fn = f"{job_tag}/{pname}/traj-{traj}/pi0-current"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_pi0_current()
    expr_names = get_expr_names(cexpr)
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    t_size = total_site[3]
    get_prop = get_get_prop()
    psel_prob = get_psel_prob()
    fsel_prob = get_fsel_prob()
    psel = psel_prob.psel
    fsel = fsel_prob.fsel
    if not fsel.is_containing(psel):
        q.displayln_info(
            -1,
            "WARNING: fsel is not containing psel. The probability weighting may be wrong.",
        )
    fsel_prob_arr = fsel_prob[:].ravel()
    psel_prob_arr = psel_prob[:].ravel()
    xg_fsel_arr = fsel.to_psel_local()[:]
    xg_psel_arr = psel[:]
    sf_pi0_current_list = []
    sf_pi0_current_arr_list = []
    for t_src in range(t_size):
        sf = q.SelectedFieldComplexD(fsel, len(expr_names))
        q.set_zero(sf)
        sf_pi0_current_list.append(sf)
        sf_pi0_current_arr_list.append(sf[:])
    #
    def load_data():
        for fidx in range(len(xg_fsel_arr)):
            yield fidx
    #
    @q.timer
    def feval(args):
        fidx = args
        xg_snk = tuple(xg_fsel_arr[fidx])
        prob_snk = fsel_prob_arr[fidx]
        val_arr = np.zeros(
            (
                t_size,
                len(expr_names),
            ),
            dtype=np.complex128,
        )
        for pidx in range(len(xg_psel_arr)):
            xg_src = tuple(xg_psel_arr[pidx])
            t_src = xg_src[3]
            prob_src = psel_prob_arr[pidx]
            if xg_snk == xg_src:
                prob = prob_src
            else:
                prob = prob_snk * prob_src
            pd = {
                "x_1": (
                    "point-snk",
                    xg_snk,
                ),
                "w": (
                    "point",
                    xg_src,
                ),
                "size": total_site,
            }
            val = eval_cexpr(cexpr, positions_dict=pd, get_prop=get_prop)
            val_arr[t_src] += val / prob
        res = (
            fidx,
            val_arr,
        )
        return res
    #
    def sum_function(val_list):
        for idx, res in enumerate(val_list):
            (
                fidx,
                val_arr,
            ) = res
            for t_src in range(t_size):
                sf_pi0_current_arr_list[t_src][fidx] = val_arr[t_src]
            if (idx + 1) % (len(xg_fsel_arr) // 128 + 16) == 0:
                q.displayln_info(f"{fname}: {idx + 1}/{len(xg_fsel_arr)}")
        return None
    #
    auto_contractor_chunk_size = get_param(
        job_tag, "measurement", "auto_contractor_chunk_size", default=128
    )
    q.timer_fork(0)
    q.parallel_map_sum(
        feval,
        load_data(),
        sum_function=sum_function,
        chunksize=auto_contractor_chunk_size,
    )
    q.displayln_info(f"{fname}: timer_display for parallel_map_sum")
    q.timer_display()
    q.timer_merge()
    sfw = q.open_fields(
        get_save_path(fn + ".acc"),
        "w",
        q.Coordinate(
            [
                2,
                2,
                2,
                4,
            ]
        ),
    )
    for t_src in range(t_size):
        sf = sf_pi0_current_list[t_src]
        tag = f"sf_pi0_current ; t_src={t_src}"
        sf.save_double(sfw, tag)
    sfw.close()
    q.qrename_info(get_save_path(fn + ".acc"), get_save_path(fn))
    sig_arr = np.zeros(
        (
            t_size,
            len(expr_names),
        ),
        dtype=np.complex128,
    )
    for t_src in range(t_size):
        sf = sf_pi0_current_list[t_src]
        q.json_results_append(
            f"{fname}: sf_pi0_current t_src={t_src} sig",
            q.get_data_sig(sf, q.RngState()),
        )
        for i, en in enumerate(expr_names):
            sig_arr[t_src, i] = q.glb_sum(
                q.get_data_sig(sf[:, i], q.RngState(f"t_src={t_src}"))
            )
    for i, en in enumerate(expr_names):
        q.json_results_append(
            f"{fname}: sf_pi0_current '{en}' sig", sig_arr[:, i].sum()
        )

@q.timer(is_timer_fork=True)
def auto_contract_pi0_gg_disc(
    job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob
):
    """
    Compute pi0 -> gamma gamma disconnected diagrams.\n
    Combines the tadpole current and pi0 current correlators to compute
    the disconnected contribution to pi0 -> gamma gamma decay. Uses field
    convolution to evaluate the spatial integral. Results are binned by
    (t1, t2, r) where t1 is the pi0 current separation, t2 is the time
    difference, and r is the spatial distance.\n
    Args:
        job_tag: Gauge ensemble identifier.
        traj: Trajectory number.
        get_get_prop: Callable returning propagator accessor function.
        get_psel_prob: Callable returning point selection probability object.
        get_fsel_prob: Callable returning field selection probability object.
    """
    fname = q.get_fname()
    fn = f"{job_tag}/{pname}/traj-{traj}/pi0-gg-disc.lat"
    if get_load_path(fn) is not None:
        return
    fn_tadpole_current = get_load_path(
        f"{job_tag}/{pname}/traj-{traj}/tadpole-current.sfield"
    )
    fn_pi0_current = get_load_path(f"{job_tag}/{pname}/traj-{traj}/pi0-current")
    if fn_tadpole_current is None or fn_pi0_current is None:
        return
    tadpole_current_expr_names = get_expr_names(get_cexpr_tadpole_current())
    pi0_current_expr_names = get_expr_names(get_cexpr_pi0_current())
    assert len(tadpole_current_expr_names) == 1 + 4 * 4
    assert len(pi0_current_expr_names) == 1 + 4 * 2
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    t_size = total_site[3]
    psel_prob = get_psel_prob()
    fsel_prob = get_fsel_prob()
    psel = psel_prob.psel
    fsel = fsel_prob.fsel
    if not fsel.is_containing(psel):
        q.displayln_info(
            -1,
            "WARNING: fsel is not containing psel. The probability weighting may be wrong.",
        )
    fsel_prob_arr = fsel_prob[:].ravel()
    psel_prob[:].ravel()
    m_l = get_param(job_tag, "m_l")
    m_h = get_param(job_tag, "m_h")
    geo = q.Geometry(total_site)
    total_volume = geo.total_volume
    r_list = get_r_list(job_tag)
    r_sq_interp_idx_coef_list = get_r_sq_interp_idx_coef_list(job_tag)
    sf_tadpole_current = q.SelectedFieldComplexD(fsel)
    sf_tadpole_current.load_double(fn_tadpole_current)
    assert sf_tadpole_current.multiplicity == len(tadpole_current_expr_names)
    sfr = q.open_fields(fn_pi0_current, "r")
    sf_pi0_current_list = []
    for t_src in range(t_size):
        sf = q.SelectedFieldComplexD(fsel)
        tag = f"sf_pi0_current ; t_src={t_src}"
        sf.load_double(sfr, tag)
        sf_pi0_current_list.append(sf)
        assert sf.multiplicity == len(pi0_current_expr_names)
    sfr.close()
    q.displayln_info(f"{len(sf_pi0_current_list)}")
    expr_names = [
        "< 1 [disc x] [con 0] >",
    ]
    for tadpole in [
        "ls",
        "sl",
        "ll",
        "ss",
    ]:
        for pi_op in [
            "pp",
            "j5",
        ]:
            expr_names += [
                f"< e(i,j,k) * (x_2-x_1)[i] * disc[ j_j(x_2) ] * j_k(x_1) * pi0(0) > (tadpole {tadpole}) (pi_op {pi_op})",
            ]
    # Some final modification for sf_tadpole_current
    sf_tadpole_current[:, 0] *= 1 / total_volume
    sf_tadpole_current[:, 1:] *= m_h - m_l
    # compute zero distance counts sum
    zero_dis_counts_sum_t_1_arr = np.zeros(t_size, dtype=np.float64)
    zero_dis_counts_sum_no_correction_t_1_arr = np.zeros(t_size, dtype=np.float64)
    sf_count = q.SelectedFieldRealD(fsel, 1)
    for t_src in range(t_size):
        sf_count[:] = (
            fsel_prob_arr[:, None]
            * sf_tadpole_current[:, 0, None].real
            * sf_pi0_current_list[t_src][:, 0, None].real
        )
        tslice_sum = sf_count.glb_sum_tslice()[:, 0]
        zero_dis_counts_sum_t_1_arr += np.roll(tslice_sum, -t_src)
        sf_count[:] = (
            sf_tadpole_current[:, 0, None].real
            * sf_pi0_current_list[t_src][:, 0, None].real
        )
        tslice_sum = sf_count.glb_sum_tslice()[:, 0]
        zero_dis_counts_sum_no_correction_t_1_arr += np.roll(tslice_sum, -t_src)
    # Convert tadpole to field
    f_tadpole_current = q.FieldComplexD(geo, sf_tadpole_current.multiplicity)
    q.set_zero(f_tadpole_current)
    f_tadpole_current @= sf_tadpole_current
    # Convert pi0 current list to field list
    f_pi0_current_list = []
    for t_src in range(t_size):
        sf = sf_pi0_current_list[t_src]
        f = q.FieldComplexD(geo, sf.multiplicity)
        q.set_zero(f)
        f @= sf
        f_pi0_current_list.append(f)
    # make new pi0 current field list based on pi0_current_sep (t_1)
    xg_arr = geo.xg_arr()
    t_arr = xg_arr[:, 3]
    f_pi0_current_sep_list = []
    for pi0_current_sep in range(t_size):
        f = q.FieldComplexD(geo, sf.multiplicity)
        q.set_zero(f)
        for t_src in range(t_size):
            f_src = f_pi0_current_list[t_src]
            sel = (t_arr - t_src) % t_size == pi0_current_sep
            f[sel] = f_src[sel]
        f_pi0_current_sep_list.append(f)
    # perform all needed convolution
    ff_idx = 0
    idx1 = [
        0,
    ]  # for con 0 (f_pi0_current_sep_list)
    idx2 = [
        0,
    ]  # for disc x (f_tadpole_current)
    info = [
        (ff_idx,),
    ]  # for (i, j, k,) tuple
    ff_idx += 1
    for i in range(3):
        for j in range(3):
            if i == j:
                continue
            for k in range(3):
                if k == i or k == j:
                    continue
                for tadpole_idx in range(4):
                    for pi_op_idx in range(2):
                        idx1.append(1 + pi_op_idx * 4 + k)
                        idx2.append(1 + tadpole_idx * 4 + j)
                        info.append(
                            (
                                ff_idx,
                                i,
                                j,
                                k,
                                tadpole_idx,
                                pi_op_idx,
                            )
                        )
                        ff_idx += 1
    idx1 = np.array(idx1, dtype=np.int32)
    idx2 = np.array(idx2, dtype=np.int32)
    info = np.array(info, dtype=object)
    assert ff_idx == len(idx1)
    assert ff_idx == len(idx2)
    assert ff_idx == len(info)
    ff_list = [
        q.field_convolution(
            f_pi0_current_sep_list[pi0_current_sep], f_tadpole_current, idx1, idx2
        )
        for pi0_current_sep in range(t_size)
    ]
    values = np.zeros(
        (
            t_size,
            t_size,
            len(r_list),
            len(expr_names),
        ),
        dtype=np.complex128,
    )
    for xg_idx, xg_rel in enumerate(geo.xg_arr()):
        x_rel = q.rel_mod_arr(xg_rel, total_site.to_numpy())
        r_sq = q.get_r_sq(x_rel)
        r_idx_low, r_idx_high, coef_low, coef_high = r_sq_interp_idx_coef_list[r_sq]
        x_rel_t = x_rel[3]
        t = x_rel_t % t_size
        for t_1 in range(t_size):
            # t_1 = pi0_current_sep
            t_2 = (t_1 + t) % t_size
            ff = ff_list[t_1]
            val = np.zeros(len(expr_names), dtype=np.complex128)
            if np.all(x_rel == 0):
                zero_dis_counts_sum = zero_dis_counts_sum_t_1_arr[t_1]
                zero_dis_counts_sum_no_correction = (
                    zero_dis_counts_sum_no_correction_t_1_arr[t_1]
                )
                assert abs(ff[xg_idx, 0] / zero_dis_counts_sum_no_correction - 1) < 1e-8
                val[0] = zero_dis_counts_sum
            else:
                val[0] = ff[xg_idx, 0]
            for ff_idx, i, j, k, tadpole_idx, pi_op_idx in info[1:]:
                x_i = q.rel_mod_sym(x_rel[i], total_site[i])  # x[i]
                eijk = q.epsilon_tensor(i, j, k)  # e(i,j,k)
                jjpi0_t1 = ff[xg_idx, ff_idx]  # disc[ j_j(x_2) ] * j_k(x_1) * pi0(0)
                val[1 + tadpole_idx * 2 + pi_op_idx] += eijk * x_i * jjpi0_t1
            values[t_1, t_2, r_idx_low] += coef_low * val
            values[t_1, t_2, r_idx_high] += coef_high * val
    res_sum = q.glb_sum(values.transpose(3, 0, 1, 2))
    res_sum *= 1.0 / total_volume
    ld_sum = q.mk_lat_data(
        [
            [
                "expr_name",
                len(expr_names),
                expr_names,
            ],
            [
                "t1",
                t_size,
                [str(q.rel_mod(t, t_size)) for t in range(t_size)],
            ],
            [
                "t2",
                t_size,
                [str(q.rel_mod(t, t_size)) for t in range(t_size)],
            ],
            [
                "r",
                len(r_list),
                [f"{r:.5f}" for r in r_list],
            ],
        ]
    )
    ld_sum.from_numpy(res_sum)
    ld_sum.save(get_save_path(fn))
    q.json_results_append(f"{fname}: ld_sum sig", q.get_data_sig(ld_sum, q.RngState()))
    for i, en in enumerate(expr_names):
        q.json_results_append(
            f"{fname}: ld_sum '{en}' sig", q.get_data_sig(ld_sum[i], q.RngState())
        )

### ------

@q.timer(is_timer_fork=True)
def run_job_inversion(job_tag, traj):
    """
    Run all quark propagator inversions for a given job_tag and trajectory.\n
    Performs the following steps:
    1. Load gauge field and compute gauge transformation
    2. Compute eigenvalues for deflated inversions (light and strange)
    3. Compute wall-source propagators (full time-slice)
    4. Compute field selection weights and probabilities
    5. Compute point-source, wall-source sparse, and smeared propagators
    6. Compute random U(1) propagators for all quark flavors\n
    Args:
        job_tag: Gauge ensemble identifier.
        traj: Trajectory number.
    """
    q.get_fname()
    #
    psel_split_num_piece = get_param(job_tag, "measurement", "psel_split_num_piece")
    fsel_psel_split_num_piece = get_param(
        job_tag, "measurement", "fsel_psel_split_num_piece"
    )
    #
    traj_gf = traj
    #
    if is_test():
        traj_gf = 1000
    #
    fns_produce = [
        f"{job_tag}/gauge-transform/traj-{traj_gf}.field",
        f"{job_tag}/points-selection/traj-{traj}.lati",
        f"{job_tag}/field-selection/traj-{traj}.field",
        #
        f"{job_tag}/field-selection-weight/traj-{traj}/weight.field",
        f"{job_tag}/field-selection-weight/traj-{traj}/f-rand-01.field",
        f"{job_tag}/field-selection-weight/traj-{traj}/fsel-prob.sfield",
        f"{job_tag}/field-selection-weight/traj-{traj}/psel-prob.lat",
        #
        f"{job_tag}/field-rand-u1/traj-{traj}/checkpoint.txt",
        #
        f"{job_tag}/points-selection-smear/traj-{traj}.lati",
        f"{job_tag}/psel_smear_median/traj-{traj}.lati",
        #
        f"{job_tag}/points-selection-split/traj-{traj}/num-piece-{psel_split_num_piece}/checkpoint.txt",
        f"{job_tag}/field-selection-split/traj-{traj}/num-piece-{fsel_psel_split_num_piece}/checkpoint.txt",
    ]
    for inv_type, quark_flavor in list(
        enumerate(get_param(job_tag, "quark_flavor_list"))
    )[:2]:
        fns_produce += [
            (
                f"{job_tag}/prop-psrc-{quark_flavor}/traj-{traj}.qar",
                f"{job_tag}/prop-psrc-{quark_flavor}/traj-{traj}/geon-info.txt",
            ),
            (
                f"{job_tag}/prop-wsrc-{quark_flavor}/traj-{traj}.qar",
                f"{job_tag}/prop-wsrc-{quark_flavor}/traj-{traj}/geon-info.txt",
            ),
            (
                f"{job_tag}/prop-smear-{quark_flavor}/traj-{traj}.qar",
                f"{job_tag}/prop-smear-{quark_flavor}/traj-{traj}/geon-info.txt",
            ),
            (
                f"{job_tag}/psel_smear_median-prop-smear-strange/traj-{traj}.qar",
                f"{job_tag}/psel_smear_median-prop-smear-strange/traj-{traj}/geon-info.txt.txt",
            ),
            f"{job_tag}/psel-prop-psrc-{quark_flavor}/traj-{traj}/checkpoint.txt",
            f"{job_tag}/psel-prop-wsrc-{quark_flavor}/traj-{traj}/checkpoint.txt",
            f"{job_tag}/psel-prop-smear-{quark_flavor}/traj-{traj}/checkpoint.txt",
        ]
    for inv_type, quark_flavor in list(
        enumerate(get_param(job_tag, "quark_flavor_list"))
    ):
        fns_produce += [
            (
                f"{job_tag}/prop-rand-u1-{quark_flavor}/traj-{traj}.qar",
                f"{job_tag}/prop-rand-u1-{quark_flavor}/traj-{traj}/geon-info.txt",
            ),
        ]
    fns_need = [
        (
            f"{job_tag}/configs/ckpoint_lat.{traj}",
            f"{job_tag}/configs/ckpoint_lat.IEEE64BIG.{traj}",
        ),
    ]
    if is_test():
        fns_need = []
    if not check_job(job_tag, traj, fns_produce, fns_need):
        return
    #
    get_gf = run_gf(job_tag, traj_gf)
    get_gt = run_gt(job_tag, traj_gf, get_gf)
    get_gf_ape = run_gf_ape(job_tag, get_gf)
    #
    get_wi = run_wi(job_tag, traj)
    #
    get_eig_light = run_eig(job_tag, traj_gf, get_gf)
    get_eig_strange = run_eig_strange(job_tag, traj_gf, get_gf)
    #
    def run_wsrc_full():
        get_eig = get_eig_light
        # run_get_inverter(job_tag, traj, inv_type=0, get_gf=get_gf, get_gt=get_gt, get_eig=get_eig)
        run_prop_wsrc_full(
            job_tag,
            traj,
            inv_type=0,
            get_gf=get_gf,
            get_eig=get_eig,
            get_gt=get_gt,
            get_wi=get_wi,
        )
        #
        get_eig = get_eig_strange
        # run_get_inverter(job_tag, traj, inv_type=1, get_gf=get_gf, get_gt=get_gt, get_eig=get_eig)
        run_prop_wsrc_full(
            job_tag,
            traj,
            inv_type=1,
            get_gf=get_gf,
            get_eig=get_eig,
            get_gt=get_gt,
            get_wi=get_wi,
        )
    #
    run_wsrc_full()
    #
    get_f_weight = run_f_weight_from_wsrc_prop_full(job_tag, traj)
    get_f_rand_01 = run_f_rand_01(job_tag, traj)
    # fsel should contain in psel (for old format, fsel from file will be combined with psel)
    get_fsel_prob = run_fsel_prob(
        job_tag, traj, get_f_rand_01=get_f_rand_01, get_f_weight=get_f_weight
    )
    get_psel_prob = run_psel_prob(
        job_tag, traj, get_f_rand_01=get_f_rand_01, get_f_weight=get_f_weight
    )
    get_psel_prob_median = run_psel_prob(
        job_tag,
        traj,
        get_f_rand_01=get_f_rand_01,
        get_f_weight=get_f_weight,
        tag="median",
    )
    get_fsel = run_fsel_from_fsel_prob(get_fsel_prob)
    get_psel = run_psel_from_psel_prob(get_psel_prob)
    run_psel_from_psel_prob(get_psel_prob_median)
    #
    run_psel_split(job_tag, traj, get_psel=get_psel, num_piece=psel_split_num_piece)
    run_fsel_split(
        job_tag, traj, get_fsel=get_fsel, num_piece=fsel_psel_split_num_piece
    )
    #
    run_field_rand_u1_dict(job_tag, traj)
    #
    get_psel_smear = run_psel_smear(job_tag, traj)
    get_psel_smear_median = run_psel_smear_median(job_tag, traj)
    #
    get_eig = get_eig_light
    run_prop_wsrc_sparse(
        job_tag,
        traj,
        inv_type=0,
        get_gf=get_gf,
        get_eig=get_eig,
        get_gt=get_gt,
        get_psel=get_psel,
        get_fsel=get_fsel,
        get_wi=get_wi,
    )
    get_eig = get_eig_strange
    run_prop_wsrc_sparse(
        job_tag,
        traj,
        inv_type=1,
        get_gf=get_gf,
        get_eig=get_eig,
        get_gt=get_gt,
        get_psel=get_psel,
        get_fsel=get_fsel,
        get_wi=get_wi,
    )
    #
    def run_with_eig():
        get_eig = get_eig_light
        # run_get_inverter(job_tag, traj, inv_type=0, get_gf=get_gf, get_eig=get_eig)
        run_prop_rand_u1(
            job_tag, traj, inv_type=0, get_gf=get_gf, get_fsel=get_fsel, get_eig=get_eig
        )
        run_prop_psrc(
            job_tag,
            traj,
            inv_type=0,
            get_gf=get_gf,
            get_eig=get_eig,
            get_gt=get_gt,
            get_psel=get_psel,
            get_fsel=get_fsel,
            get_f_rand_01=get_f_rand_01,
        )
        run_prop_smear(
            job_tag,
            traj,
            inv_type=0,
            get_gf=get_gf,
            get_gf_ape=get_gf_ape,
            get_eig=get_eig,
            get_gt=get_gt,
            get_psel=get_psel,
            get_fsel=get_fsel,
            get_psel_smear=get_psel_smear,
            get_psel_smear_median=get_psel_smear_median,
        )
        q.clean_cache(q.cache_inv)
    #
    def run_with_eig_strange():
        get_eig = get_eig_strange
        # run_get_inverter(job_tag, traj, inv_type=1, get_gf=get_gf, get_eig=get_eig)
        run_prop_rand_u1(
            job_tag, traj, inv_type=1, get_gf=get_gf, get_fsel=get_fsel, get_eig=get_eig
        )
        run_prop_psrc(
            job_tag,
            traj,
            inv_type=1,
            get_gf=get_gf,
            get_eig=get_eig,
            get_gt=get_gt,
            get_psel=get_psel,
            get_fsel=get_fsel,
            get_f_rand_01=get_f_rand_01,
        )
        run_prop_smear(
            job_tag,
            traj,
            inv_type=1,
            get_gf=get_gf,
            get_gf_ape=get_gf_ape,
            get_eig=get_eig,
            get_gt=get_gt,
            get_psel=get_psel,
            get_fsel=get_fsel,
            get_psel_smear=get_psel_smear,
            get_psel_smear_median=get_psel_smear_median,
        )
        q.clean_cache(q.cache_inv)
    #
    def run_charm():
        # run_get_inverter(job_tag, traj, inv_type=2, get_gf=get_gf)
        for inv_type, quark_flavor in list(
            enumerate(get_param(job_tag, "quark_flavor_list"))
        )[2:]:
            run_prop_rand_u1(
                job_tag, traj, inv_type=inv_type, get_gf=get_gf, get_fsel=get_fsel
            )
        q.clean_cache(q.cache_inv)
    #
    run_with_eig()
    run_with_eig_strange()
    run_charm()
    #
    q.clean_cache()
    if q.obtained_lock_history_list:
        q.timer_display()

@q.timer(is_timer_fork=True)
def run_job_contract(job_tag, traj):
    """
    Run all contraction measurements for a given job_tag and trajectory.\n
    Loads pre-computed propagators and performs all auto-contractor
    measurements including:
    - Meson correlation functions (wall-wall, wall-point, point-wall, point-point)
    - Pi0 -> J J parity-violating correlators
    - Meson two-current and tensor current correlators
    - Meson mass insertion correlators
    - Tadpole current and pi0 current correlators
    - Pi0 gamma-gamma disconnected diagrams\n
    Args:
        job_tag: Gauge ensemble identifier.
        traj: Trajectory number.
    """
    #
    traj_gf = traj
    if job_tag[:5] == "test-":
        # ADJUST ME
        traj_gf = 1000
        #
    fns_produce = [
        f"{job_tag}/{pname}/traj-{traj}/checkpoint.txt",
        #
    ]
    fns_need = [
        (
            f"{job_tag}/prop-psrc-light/traj-{traj}.qar",
            f"{job_tag}/prop-psrc-light/traj-{traj}/geon-info.txt",
        ),
        (
            f"{job_tag}/psel-prop-psrc-light/traj-{traj}.qar",
            f"{job_tag}/psel-prop-psrc-light/traj-{traj}/checkpoint.txt",
        ),
        (
            f"{job_tag}/prop-psrc-strange/traj-{traj}.qar",
            f"{job_tag}/prop-psrc-strange/traj-{traj}/geon-info.txt",
        ),
        (
            f"{job_tag}/psel-prop-psrc-strange/traj-{traj}.qar",
            f"{job_tag}/psel-prop-psrc-strange/traj-{traj}/checkpoint.txt",
        ),
        #
        (
            f"{job_tag}/prop-wsrc-light/traj-{traj}.qar",
            f"{job_tag}/prop-wsrc-light/traj-{traj}/geon-info.txt",
        ),
        (
            f"{job_tag}/psel-prop-wsrc-light/traj-{traj}.qar",
            f"{job_tag}/psel-prop-wsrc-light/traj-{traj}/checkpoint.txt",
        ),
        (
            f"{job_tag}/prop-wsrc-strange/traj-{traj}.qar",
            f"{job_tag}/prop-wsrc-strange/traj-{traj}/geon-info.txt",
        ),
        (
            f"{job_tag}/psel-prop-wsrc-strange/traj-{traj}.qar",
            f"{job_tag}/psel-prop-wsrc-strange/traj-{traj}/checkpoint.txt",
        ),
        f"{job_tag}/gauge-transform/traj-{traj_gf}.field",
        f"{job_tag}/points-selection/traj-{traj}.lati",
        f"{job_tag}/field-selection/traj-{traj}.field",
        # f"{job_tag}/wall-src-info-light/traj-{traj}.txt",
        # f"{job_tag}/wall-src-info-strange/traj-{traj}.txt",
        # (f"{job_tag}/configs/ckpoint_lat.{traj}", f"{job_tag}/configs/ckpoint_lat.IEEE64BIG.{traj}",),
    ]
    if not check_job(job_tag, traj, fns_produce, fns_need):
        return
    #
    get_gf = None
    get_gt = run_gt(job_tag, traj_gf, get_gf)
    #
    run_wi(job_tag, traj)
    #
    get_f_weight = run_f_weight_from_wsrc_prop_full(job_tag, traj)
    get_f_rand_01 = run_f_rand_01(job_tag, traj)
    get_fsel_prob = run_fsel_prob(
        job_tag, traj, get_f_rand_01=get_f_rand_01, get_f_weight=get_f_weight
    )
    get_psel_prob = run_psel_prob(
        job_tag, traj, get_f_rand_01=get_f_rand_01, get_f_weight=get_f_weight
    )
    get_fsel = run_fsel_from_fsel_prob(get_fsel_prob)
    get_psel = run_psel_from_psel_prob(get_psel_prob)
    #
    get_psel_smear = run_psel_smear(job_tag, traj)
    run_psel_smear_median(job_tag, traj)
    #
    get_get_prop = run_get_prop(
        job_tag,
        traj,
        get_gf=get_gf,
        get_gt=get_gt,
        get_psel=get_psel,
        get_fsel=get_fsel,
        get_psel_smear=get_psel_smear,
        prop_types=[
            "wsrc psel s",
            "wsrc psel l",
            "wsrc fsel s",
            "wsrc fsel l",
            "psrc psel s",
            "psrc psel l",
            "psrc fsel s",
            "psrc fsel l",
            # "rand_u1 fsel c",
            # "rand_u1 fsel s",
            # "rand_u1 fsel l",
        ],
    )
    #
    run_r_list(job_tag)
    #
    fn_checkpoint = f"{job_tag}/{pname}/traj-{traj}/checkpoint.txt"
    if get_load_path(fn_checkpoint) is None:
        if q.obtain_lock(f"locks/{job_tag}-{traj}-{pname}"):
            get_prop = get_get_prop()
            if get_prop is not None:
                q.timer_fork()
                # ADJUST ME
                auto_contract_meson_corr(
                    job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob
                )
                auto_contract_meson_corr_psnk(
                    job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob
                )
                auto_contract_meson_corr_psrc(
                    job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob
                )
                auto_contract_meson_corr_psnk_psrc(
                    job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob
                )
                auto_contract_pi0_jjp(
                    job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob
                )
                auto_contract_meson_jj(
                    job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob
                )
                auto_contract_meson_jt(
                    job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob
                )
                auto_contract_meson_m(
                    job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob
                )
                auto_contract_tadpole_current(
                    job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob
                )
                auto_contract_pi0_current(
                    job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob
                )
                auto_contract_pi0_gg_disc(
                    job_tag, traj, get_get_prop, get_psel_prob, get_fsel_prob
                )
                #
                q.qtouch_info(get_save_path(fn_checkpoint))
                q.displayln_info("timer_display for runjob")
                q.timer_display()
                q.timer_merge()
            q.release_lock()
    q.clean_cache()
    if q.obtained_lock_history_list:
        q.timer_display()

### ------

@q.timer(is_timer_fork=True)
def get_all_cexpr():
    """
    Pre-compile and benchmark all contraction expressions.\n
    Builds all compiled expressions used by the measurement functions
    and runs benchmark evaluations to ensure they are cached for later use.
    """
    benchmark_eval_cexpr(get_cexpr_meson_corr())
    benchmark_eval_cexpr(get_cexpr_meson_corr(False))
    benchmark_eval_cexpr(get_cexpr_meson_m())
    benchmark_eval_cexpr(get_cexpr_meson_jt())
    benchmark_eval_cexpr(get_cexpr_meson_jj())
    benchmark_eval_cexpr(get_cexpr_pi0_jjp())
    benchmark_eval_cexpr(get_cexpr_tadpole_current())
    benchmark_eval_cexpr(get_cexpr_pi0_current())

### ------

job_tag = "24D"
set_param(job_tag, "traj_list")(
    [
        2430,
        2550,
        2590,
        2610,
        2630,
        2940,
        2960,
    ]
)
set_param(job_tag, "meson_tensor_tsep")(8)
set_param(job_tag, "pi0_jjp_threshold")(0.02)
set_param(job_tag, "measurement", "auto_contractor_chunk_size")(128)

job_tag = "48I"
set_param(job_tag, "seed")("48I")
set_param(job_tag, "traj_list")(list(range(1000, 2000, 20)))
set_param(job_tag, "meson_tensor_tsep")(12)
set_param(job_tag, "pi0_jjp_threshold")(0.01)
set_param(job_tag, "measurement", "auto_contractor_chunk_size")(128)

job_tag = "64I"
set_param(job_tag, "seed")("64I")
set_param(job_tag, "traj_list")(list(range(1200, 3680, 20)))
set_param(job_tag, "meson_tensor_tsep")(18)
set_param(job_tag, "pi0_jjp_threshold")(0.0005)
set_param(job_tag, "measurement", "auto_contractor_chunk_size")(128)

job_tag = "64I-pq"
set_param(job_tag, "seed")("64I")
set_param(job_tag, "traj_list")(list(range(1200, 3680, 80)))
set_param(job_tag, "meson_tensor_tsep")(18)
set_param(job_tag, "pi0_jjp_threshold")(0.0005)
set_param(job_tag, "measurement", "auto_contractor_chunk_size")(128)

# ----

job_tag = "test-4nt8-checker"
#
set_param(job_tag, "seed")("test-4nt8")
set_param(job_tag, "traj_list")(list(range(1000, 1001)))
#
set_param(job_tag, "total_site")(
    [
        4,
        4,
        4,
        8,
    ]
)
set_param(job_tag, "load_config_params", "twist_boundary_at_boundary")(
    [
        0.0,
        0.0,
        0.0,
        -0.5,
    ]
)
#
set_param(job_tag, "mk_sample_gauge_field", "rand_n_step")(2)
set_param(job_tag, "mk_sample_gauge_field", "flow_n_step")(8)
set_param(job_tag, "mk_sample_gauge_field", "hmc_n_traj")(1)
#
set_param(job_tag, "quark_flavor_list")(
    [
        "light",
        "strange",
        "charm-1",
    ]
)
set_param(job_tag, "quark_mass_list")(
    [
        0.01,
        0.04,
        0.2,
    ]
)
set_param(job_tag, "fermion_params", 0, 0)(
    {
        "Ls": 8,
        "M5": 1.8,
        "b": 1.5,
        "c": 0.5,
        "boundary_phases": [
            1.0,
            1.0,
            1.0,
            1.0,
        ],
    }
)
for inv_type, mass in enumerate(get_param(job_tag, "quark_mass_list")):
    set_param(job_tag, "fermion_params", inv_type, 0)(
        get_param(job_tag, "fermion_params", 0, 0).copy()
    )
    set_param(job_tag, "fermion_params", inv_type, 0, "mass")(mass)
    for inv_acc in [
        0,
        1,
        2,
    ]:
        # set_param(job_tag, "fermion_params", inv_type, inv_acc)(get_param(job_tag, "fermion_params", inv_type, 0).copy())
        set_param(job_tag, f"cg_params-{inv_type}-{inv_acc}", "maxiter")(10)
        set_param(job_tag, f"cg_params-{inv_type}-{inv_acc}", "maxcycle")(1 + inv_acc)
#
set_param(job_tag, "lanc_params", 0, 0, "fermion_params")(
    get_param(job_tag, "fermion_params", 0, 0).copy()
)
set_param(job_tag, "lanc_params", 0, 0, "cheby_params")(
    {
        "low": 0.3,
        "high": 5.5,
        "order": 40,
    }
)
set_param(job_tag, "lanc_params", 0, 0, "irl_params")(
    {
        "Nstop": 50,
        "Nk": 80,
        "Nm": 100,
        "resid": 1e-8,
        "betastp": 0.0,
        "maxiter": 20,
        "Nminres": 0,
    }
)
set_param(job_tag, "lanc_params", 0, 0, "pit_params")(
    {
        "eps": 0.01,
        "maxiter": 500,
        "real": True,
    }
)
#
# set_param(job_tag, "clanc_params", 0, 0, "nbasis")(100)
# set_param(job_tag, "clanc_params", 0, 0, "block")([ 4, 4, 2, 2, ])
# set_param(job_tag, "clanc_params", 0, 0, "cheby_params")({ "low": 0.5, "high": 5.5, "order": 40, })
# set_param(job_tag, "clanc_params", 0, 0, "save_params")({ "nsingle": 100, "mpi": [ 1, 1, 1, 4, ], })
# set_param(job_tag, "clanc_params", 0, 0, "irl_params")({ "Nstop": 100, "Nk": 150, "Nm": 200, "resid": 1e-8, "betastp": 0.0, "maxiter": 20, "Nminres": 0, })
# set_param(job_tag, "clanc_params", 0, 0, "smoother_params")({'eps': 1e-08, 'maxiter': 10})
#
# set_param(job_tag, "clanc_params", 1, 0)(get_param(job_tag, "clanc_params", 0, 0).copy())
# set_param(job_tag, "lanc_params", 1, 0)(get_param(job_tag, "lanc_params", 0, 0).copy())
# set_param(job_tag, "lanc_params", 1, 0, "fermion_params")(get_param(job_tag, "fermion_params", 1, 0).copy())
#
set_param(job_tag, "field_selection_psel_rate")(1 / 32)
set_param(job_tag, "field_selection_psel_rate_median")(1 / 16)
set_param(job_tag, "field_selection_fsel_rate")(1 / 8)
set_param(job_tag, "field_selection_fsel_psrc_prop_norm_threshold")(1e-3)
#
set_param(job_tag, "prob_exact_wsrc")(1 / 4)
#
set_param(job_tag, "prob_acc_1_psrc")(1 / 4)
set_param(job_tag, "prob_acc_2_psrc")(1 / 16)
#
set_param(job_tag, "n_per_tslice_smear")(2)
set_param(job_tag, "n_per_tslice_smear_median")(8)
set_param(job_tag, "gf_ape_smear_coef")(0.5)
set_param(job_tag, "gf_ape_smear_step")(30)
set_param(job_tag, "prop_smear_coef")(0.9375)
set_param(job_tag, "prop_smear_step")(10)
set_param(job_tag, "prob_acc_1_smear")(1 / 4)
set_param(job_tag, "prob_acc_2_smear")(1 / 16)
#
set_param(job_tag, "measurement", "psel_split_num_piece")(2)
set_param(job_tag, "measurement", "fsel_psel_split_num_piece")(4)
set_param(job_tag, "prob_acc_1_rand_u1_sparse")(1 / 4)
set_param(job_tag, "prob_acc_2_rand_u1_sparse")(1 / 16)
#
set_param(job_tag, "n_rand_u1_fsel")(4)
set_param(job_tag, "prob_acc_1_rand_u1")(1 / 4)
set_param(job_tag, "prob_acc_2_rand_u1")(1 / 16)
#
set_param(job_tag, "m_l")(get_param(job_tag, "quark_mass_list")[0])
set_param(job_tag, "m_h")(get_param(job_tag, "quark_mass_list")[1])
#
set_param(job_tag, "meson_tensor_tsep")(1)
#
set_param(job_tag, "pi0_jjp_threshold")(0.001)
#
set_param(job_tag, "measurement", "auto_contract_meson_corr_wf", "sample_num")(32)
set_param(job_tag, "measurement", "auto_contract_meson_corr_wf", "sample_size")(2)
set_param(job_tag, "measurement", "auto_contract_meson_corr_wf", "t_sep_range")(6)
set_param(
    job_tag, "measurement", "auto_contract_meson_meson_i0_j0_corr_wf", "sample_num"
)(32)
set_param(
    job_tag, "measurement", "auto_contract_meson_meson_i0_j0_corr_wf", "sample_size"
)(2)
set_param(
    job_tag, "measurement", "auto_contract_meson_meson_i0_j0_corr_wf", "t_sep_range"
)(6)
set_param(job_tag, "measurement", "auto_contractor_chunk_size")(2)

# ----

##################### CMD options #####################

job_tag_list_default = [
    "test-4nt8-checker",
]
job_tag_list_str_default = ",".join(job_tag_list_default)
job_tag_list = q.get_arg("--job_tag_list", default=job_tag_list_str_default).split(",")

is_performing_inversion = q.get_arg("--no-inversion", default=None) is None

is_performing_contraction = q.get_arg("--no-contract", default=None) is None

#######################################################

def gracefully_finish():
    """
    Clean up and finalize the script execution.\n
    Displays timer information, checks log JSON for test mode, and
    terminates GPT/qlat MPI environment.
    """
    q.displayln_info("Begin to gracefully_finish.")
    q.timer_display()
    if is_test():
        q.json_results_append(
            f"q.obtained_lock_history_list={q.obtained_lock_history_list}"
        )
        q.check_log_json(__file__, check_eps=5e-5)
    qg.end_with_gpt()
    q.displayln_info("CHECK: finished successfully.")
    exit()

def try_gracefully_finish():
    """
    Call `gracefully_finish` if not test and if some work is done (q.obtained_lock_history_list != [])
    """
    if (not is_test()) and (len(q.obtained_lock_history_list) > 0):
        gracefully_finish()

if __name__ == "__main__":
    qg.begin_with_gpt()
    q.check_time_limit()
    get_all_cexpr()
    #
    job_tag_traj_list = []
    for job_tag in job_tag_list:
        run_params(job_tag)
        traj_list = get_param(job_tag, "traj_list")
        for traj in traj_list:
            job_tag_traj_list.append(
                (
                    job_tag,
                    traj,
                )
            )
    if not is_test():
        job_tag_traj_list = q.random_permute(
            job_tag_traj_list, q.RngState(f"{q.get_time()}")
        )
        job_tag_traj_list = q.get_comm().bcast(job_tag_traj_list)
    for job_tag, traj in job_tag_traj_list:
        if is_performing_inversion:
            q.check_time_limit()
            run_job_inversion(job_tag, traj)
            q.clean_cache()
            try_gracefully_finish()
    for job_tag, traj in job_tag_traj_list:
        if is_performing_contraction:
            q.check_time_limit()
            run_job_contract(job_tag, traj)
            q.clean_cache()
            try_gracefully_finish()
    #
    gracefully_finish()

# ----
