#!/usr/bin/env python3
"""
Auto-contractor lattice QCD data generation script.\n
This script computes meson correlation functions and related observables
using lattice QCD propagators. It generates data for auto-contractor
measurements including:\n
- Meson two-point correlators (wall-wall, wall-point, point-wall, point-point)
- Meson tensor currents (J_mu * M)
- Meson mass insertions (m_q * M)
- Meson two-current correlators (J * J * M)
- Meson weak-current correlators (J_w * J * M)
- Tadpole current correlators (disconnected contributions)
- Pi0 current correlators (connected contributions)
- Pi0 gamma-gamma disconnected diagrams\n
The script supports multiple gauge ensembles (24D, 48I, 64I, 64I-pq, 64I-pq2) and
test configurations. It uses adaptive sampling for point-source/point-sink
measurements with probability-weighted estimators.\n
Usage:
    python gpt-qlat-data-gen-auto.py [--job_tag_list tag1,tag2] [--no_inversion] [--no_contract]
"""

import argparse
import sys

import qlat_gpt as qg

import numpy as np
import qlat as q

from qlat_scripts.v1 import (
    benchmark_eval_cexpr,
    check_job,
    get_expr_names,
    get_load_path,
    get_param,
    get_r_sq_interp_idx_coef_list,
    get_save_path,
    is_test,
    load_path_list,
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
    mk_fac,
    mk_scalar5,
)
from auto_contractor.eval import (
    cache_compiled_cexpr,
    eval_cexpr,
)

# ----

load_path_list[:] = [
    "results",
    "qcddata",
    "qcddata-1",
    "qcddata-2",
    "qcddata-3",
    "qcddata-4",
    "qcddata-5",
    "qcddata-6",
]

is_cython = not is_test()

pname = "topo_sqr"

# ----

@q.timer
def get_cexpr_corr():
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
    fn_base = "cache/auto_contract_cexpr/get_cexpr_corr"
    #
    def calc_cexpr():
        diagram_type_dict = dict()
        diagram_type_dict[()] = "T1"
        diagram_type_dict[((("x_1", "x_2"), 1), (("x_2", "x_1"), 1))] = "Type1"
        diagram_type_dict[((("x_1", "x_1"), 1), (("x_2", "x_2"), 1))] = "Type2"
        exprs_1_list = [
            mk_fac(1) + "1",
        ]
        exprs_corr_list = []
        #
        def mk_op_light(p):
            return mk_scalar5("u", "u", p) + mk_scalar5("d", "d", p)
        #
        def mk_op_strange(p):
            return mk_scalar5("s", "s", p)
        #
        exprs_corr_list += [
            (
                mk_op_light("x_2") * mk_op_light("x_1"),
                "Type1",
                "Type2",
            ),
            (
                mk_op_strange("x_2") * mk_op_strange("x_1"),
                "Type1",
                "Type2",
            ),
            (
                mk_op_light("x_2") * mk_op_strange("x_1"),
                "Type1",
                "Type2",
            ),
            (
                mk_op_strange("x_2") * mk_op_light("x_1"),
                "Type1",
                "Type2",
            ),
        ]
        #
        exprs = exprs_1_list + exprs_corr_list
        cexpr = contract_simplify_compile(
            *exprs, is_isospin_symmetric_limit=True, diagram_type_dict=diagram_type_dict
        )
        return cexpr
    #
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

@q.timer(is_timer_fork=True)
def auto_contract_topo_corr(job_tag, traj, get_get_prop, get_psel_prob):
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
    """
    fname = q.get_fname()
    fn = f"{job_tag}/{pname}/traj-{traj}/topo_corr.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_corr()
    expr_names = get_expr_names(cexpr)
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    total_site[3]
    get_prop = get_get_prop()
    psel_prob = get_psel_prob()
    psel = psel_prob.psel
    psel_prob_arr = psel_prob[:].ravel()
    xg_psel_arr = psel[:]
    geo = q.Geometry(total_site)
    total_volume = geo.total_volume
    get_r_sq_interp_idx_coef_list(job_tag)
    #
    def load_data():
        for pidx in q.get_mpi_chunk(range(len(xg_psel_arr)), rng_state=None):
            yield pidx
    #
    @q.timer
    def feval(args):
        pidx = args
        xg_src = tuple(xg_psel_arr[pidx])
        prob_src = psel_prob_arr[pidx]
        res_list = []
        for idx in range(len(xg_psel_arr)):
            xg_snk = tuple(xg_psel_arr[idx])
            if xg_snk == xg_src:
                prob_snk = 1.0
            else:
                prob_snk = psel_prob_arr[idx]
            prob = prob_src * prob_snk
            x_rel = np.array(
                [q.rel_mod(xg_snk[mu] - xg_src[mu], total_site[mu]) for mu in range(4)],
                dtype=np.int32,
            )
            x_rel_abs = abs(x_rel)
            pd = {
                "x_2": (
                    "point",
                    xg_snk,
                ),
                "x_1": (
                    "point",
                    xg_src,
                ),
                "size": total_site,
            }
            val = eval_cexpr(cexpr, positions_dict=pd, get_prop=get_prop)
            res_list.append((val / prob, x_rel_abs))
        return res_list
    #
    def sum_function(val_list):
        values = np.zeros(
            (
                total_site[3] // 2 + 1,
                total_site[2] // 2 + 1,
                total_site[1] // 2 + 1,
                total_site[0] // 2 + 1,
                len(expr_names),
            ),
            dtype=complex,
        )
        for idx, res_list in enumerate(val_list):
            for val, x_rel_abs in res_list:
                x, y, z, t = x_rel_abs
                values[t, z, y, x] += val
            q.displayln_info(f"{fname}: {idx + 1}/{len(xg_psel_arr)}")
        return values.transpose(4, 0, 1, 2, 3)
    #
    q.timer_fork(0)
    res_sum = q.parallel_map_sum(
        feval, load_data(), sum_function=sum_function, chunksize=1
    )
    res_sum = q.glb_sum(res_sum)
    q.displayln_info(f"{fname}: timer_display for parallel_map_sum")
    q.timer_display()
    q.timer_merge()
    res_sum *= 1.0 / total_volume
    q.displayln_info(res_sum[0].sum(1))
    ld = q.mk_lat_data(
        [
            [
                "expr_name",
                len(expr_names),
                expr_names,
            ],
            [
                "abs(t)",
                total_site[3] // 2 + 1,
            ],
            [
                "abs(z)",
                total_site[2] // 2 + 1,
            ],
            [
                "abs(y)",
                total_site[1] // 2 + 1,
            ],
            [
                "abs(x)",
                total_site[0] // 2 + 1,
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
                f"{job_tag}/psel_smear_median-prop-smear-strange/traj-{traj}/geon-info.txt",
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
    - Meson two-current and tensor current correlators
    - Meson mass insertion correlators
    - Meson weak-current correlators
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
            f"{job_tag}/psel-prop-psrc-light/traj-{traj}.qar",
            f"{job_tag}/psel-prop-psrc-light/traj-{traj}/checkpoint.txt",
        ),
        (
            f"{job_tag}/psel-prop-psrc-strange/traj-{traj}.qar",
            f"{job_tag}/psel-prop-psrc-strange/traj-{traj}/checkpoint.txt",
        ),
        (
            f"{job_tag}/psel-prop-wsrc-light/traj-{traj}.qar",
            f"{job_tag}/psel-prop-wsrc-light/traj-{traj}/checkpoint.txt",
        ),
        (
            f"{job_tag}/psel-prop-wsrc-strange/traj-{traj}.qar",
            f"{job_tag}/psel-prop-wsrc-strange/traj-{traj}/checkpoint.txt",
        ),
        f"{job_tag}/gauge-transform/traj-{traj_gf}.field",
        f"{job_tag}/points-selection/traj-{traj}.lati",
        f"{job_tag}/field-selection/traj-{traj}.field",
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
    # get_psel_smear = run_psel_smear(job_tag, traj)
    #
    get_get_prop = run_get_prop(
        job_tag,
        traj,
        get_gf=get_gf,
        get_gt=get_gt,
        get_psel=get_psel,
        get_fsel=get_fsel,
        prop_types=[
            "wsrc psel s",
            "wsrc psel l",
            "psrc psel s",
            "psrc psel l",
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
                auto_contract_topo_corr(job_tag, traj, get_get_prop, get_psel_prob)
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
    benchmark_eval_cexpr(get_cexpr_corr())

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
set_param(job_tag, "meson_jwjj_threshold")(0.02)
set_param(job_tag, "measurement", "auto_contractor_chunk_size")(128)

job_tag = "48I"
set_param(job_tag, "seed")("48I")
set_param(job_tag, "n_per_tslice_smear_median")(512)
set_param(job_tag, "field_selection_psel_rate_median")(1 / 16)
set_param(job_tag, "meson_tensor_tsep")(13)  # Previously set to 12 in old runs
set_param(job_tag, "meson_jwjj_threshold")(0.01)
set_param(job_tag, "measurement", "auto_contractor_chunk_size")(128)
set_param(job_tag, "traj_list")(
    [
        2055,
        1025,
        1222,
        2085,
        1815,
        1715,
        1865,
        1102,
        1472,
        1845,
        2165,
        1352,
        1212,
        1312,
        2005,
        1242,
        1765,
        1885,
        1302,
        1015,
        1975,
        1432,
        2105,
        1095,
        1705,
        1595,
        1005,
        1505,
        1422,
        2135,
        1075,
        1402,
        1262,
        1492,
        1515,
        1985,
        1875,
        1322,
        1995,
        1585,
        1122,
        1555,
        1162,
        1452,
        1725,
        1342,
        1132,
        1855,
        1935,
        2065,
        1775,
        1412,
        2115,
        995,
        1382,
        1272,
        1392,
        2175,
        1055,
        1142,
        1482,
        1172,
        1202,
        1755,
        1152,
        1525,
        975,
        1605,
        1735,
        2045,
        1805,
        1065,
        1925,
        1192,
        1955,
        1085,
        1635,
        1292,
        2025,
        1565,
        1182,
        1945,
        1282,
        1252,
        985,
        1332,
        1625,
        1035,
        2095,
        1462,
        2015,
        1232,
        1905,
        1112,
        1045,
        2155,
        1965,
        1895,
        1615,
        1362,
        1745,
        1795,
        1785,
        1372,
        1545,
        2145,
        2125,
        1915,
        1575,
        1835,
        1825,
        1535,
        1442,
    ]
)

job_tag = "64I"
set_param(job_tag, "seed")("64I")
set_param(job_tag, "n_per_tslice_smear_median")(512)
set_param(job_tag, "field_selection_psel_rate_median")(1 / 16)
set_param(job_tag, "meson_tensor_tsep")(18)
set_param(job_tag, "meson_jwjj_threshold")(0.0005)
set_param(job_tag, "measurement", "auto_contractor_chunk_size")(128)
set_param(job_tag, "traj_list")(
    [
        2880,
        3600,
        1520,
        2000,
        3120,
        2800,
        3280,
        1440,
        2400,
        3360,
        3200,
        2560,
        2160,
        1600,
        1200,
        1280,
        1360,
        2080,
        3040,
        2320,
        2960,
        2240,
        2720,
        1840,
        2480,
        1920,
        3440,
        1680,
        2640,
        1760,
        3520,
    ]
    + [
        1740,
        2520,
        2060,
        2400,
        2580,
        3100,
        2880,
        3420,
        2140,
        1420,
        1860,
        2000,
        2640,
        2340,
        3320,
        1780,
        1200,
        3460,
        2660,
        1880,
        1480,
        3180,
        1620,
        2020,
        1800,
        2440,
        1920,
        3140,
        2280,
        1580,
        3200,
        2860,
        3160,
        3080,
        1540,
        3120,
        1720,
        3040,
        2200,
        2040,
        1380,
        1680,
        2320,
        1240,
        3360,
        1520,
        1560,
        1700,
        2120,
        2080,
        2720,
        1320,
        1960,
        3280,
        3340,
        2780,
        3500,
        2560,
        1980,
        2420,
        3220,
        3440,
        2460,
        1460,
        1220,
        3560,
        1360,
        2740,
        3300,
        2100,
        2900,
        3000,
        1940,
        3240,
        1760,
        2960,
        1500,
        1300,
        2500,
        2760,
        2480,
        3520,
        1640,
        2820,
        2380,
        3400,
        2220,
        1440,
        1340,
        1900,
        2620,
        3060,
        3020,
        2700,
        2980,
        1400,
        2360,
        2940,
        2800,
        2680,
        1280,
        2840,
        2300,
        1260,
        2260,
        2540,
        1820,
        1840,
        1600,
        3380,
        1660,
        3260,
        3480,
        2180,
        2600,
        3600,
        2160,
        2240,
        2920,
    ]
)

job_tag = "64I-pq"
set_param(job_tag, "seed")("64I")
set_param(job_tag, "n_per_tslice_smear_median")(512)
set_param(job_tag, "field_selection_psel_rate_median")(1 / 16)
set_param(job_tag, "meson_tensor_tsep")(18)
set_param(job_tag, "meson_jwjj_threshold")(0.0005)
set_param(job_tag, "measurement", "auto_contractor_chunk_size")(128)
set_param(job_tag, "traj_list")(
    [
        2880,
        3600,
        1520,
        2000,
        3120,
        2800,
        3280,
        1440,
        2400,
        3360,
        3200,
        2560,
        2160,
        1600,
        1200,
        1280,
        1360,
        2080,
        3040,
        2320,
        2960,
        2240,
        2720,
        1840,
        2480,
        1920,
        3440,
        1680,
        2640,
        1760,
        3520,
    ]
)

job_tag = "64I-pq2"
set_param(job_tag, "job_tag")("64I-pq2")
set_param(job_tag, "seed")("64I")
set_param(job_tag, "n_per_tslice_smear_median")(512)
set_param(job_tag, "field_selection_psel_rate_median")(1 / 16)
set_param(job_tag, "meson_tensor_tsep")(18)
set_param(job_tag, "meson_jwjj_threshold")(0.0005)
set_param(job_tag, "measurement", "auto_contractor_chunk_size")(128)
set_param(job_tag, "total_site")([64, 64, 64, 128])
set_param(job_tag, "load_config_params")(
    {"twist_boundary_at_boundary": [0.0, 0.0, 0.0, 0.0]}
)
set_param(job_tag, "fermion_params")(
    {
        0: {
            0: {
                "M5": 1.8,
                "boundary_phases": [1.0, 1.0, 1.0, 1.0],
                "b": 1.5,
                "c": 0.5,
                "mass": 0.0006203,
                "Ls": 12,
            },
            1: {
                "M5": 1.8,
                "boundary_phases": [1.0, 1.0, 1.0, 1.0],
                "b": 1.5,
                "c": 0.5,
                "mass": 0.0006203,
                "Ls": 12,
            },
            2: {
                "M5": 1.8,
                "boundary_phases": [1.0, 1.0, 1.0, 1.0],
                "b": 1.5,
                "c": 0.5,
                "mass": 0.0006203,
                "Ls": 12,
            },
        },
        1: {
            0: {
                "M5": 1.8,
                "boundary_phases": [1.0, 1.0, 1.0, 1.0],
                "b": 1.5,
                "c": 0.5,
                "mass": 0.02539,
                "Ls": 12,
            },
            1: {
                "M5": 1.8,
                "boundary_phases": [1.0, 1.0, 1.0, 1.0],
                "b": 1.5,
                "c": 0.5,
                "mass": 0.02539,
                "Ls": 12,
            },
            2: {
                "M5": 1.8,
                "boundary_phases": [1.0, 1.0, 1.0, 1.0],
                "b": 1.5,
                "c": 0.5,
                "mass": 0.02539,
                "Ls": 12,
            },
        },
        2: {
            0: {
                "M5": 1.8,
                "boundary_phases": [1.0, 1.0, 1.0, 1.0],
                "b": 1.5,
                "c": 0.5,
                "mass": 0.31,
                "Ls": 12,
            },
            1: {
                "M5": 1.8,
                "boundary_phases": [1.0, 1.0, 1.0, 1.0],
                "b": 1.5,
                "c": 0.5,
                "mass": 0.31,
                "Ls": 12,
            },
            2: {
                "M5": 1.8,
                "boundary_phases": [1.0, 1.0, 1.0, 1.0],
                "b": 1.5,
                "c": 0.5,
                "mass": 0.31,
                "Ls": 12,
            },
        },
    }
)
set_param(job_tag, "field_selection_fsel_rate")(1 / 32)
set_param(job_tag, "field_selection_psel_rate")(2048 / (64**3 * 128))
set_param(job_tag, "field_selection_fsel_psrc_prop_norm_threshold")(2e-5)
set_param(job_tag, "n_points_psel")(1024)
set_param(job_tag, "prob_exact_wsrc")(0.015625)
set_param(job_tag, "n_per_tslice_smear")(8)
set_param(job_tag, "prob_acc_1_smear")(0.03125)
set_param(job_tag, "prob_acc_2_smear")(0.0078125)
set_param(job_tag, "prob_acc_1_psrc")(0.03125)
set_param(job_tag, "prob_acc_2_psrc")(0.0078125)
set_param(job_tag, "n_rand_u1_fsel")(64)
set_param(job_tag, "prob_acc_1_rand_u1")(0.03125)
set_param(job_tag, "prob_acc_2_rand_u1")(0.0078125)
set_param(job_tag, "prop_smear_coef")(0.9375)
set_param(job_tag, "prop_smear_step")(54)
set_param(job_tag, "gf_ape_smear_coef")(0.5)
set_param(job_tag, "gf_ape_smear_step")(30)
set_param(job_tag, "meson_tsep_list")(list(range(6, 128, 4)))
set_param(job_tag, "a_inv_gev")(2.359)
set_param(job_tag, "zz_vv")(0.74293)
set_param(job_tag, "zz_aa")(0.74341)
set_param(job_tag, "m_res")(0.0003116)
set_param(job_tag, "m_l")(0.000678)
set_param(job_tag, "m_h")(0.02539)
set_param(job_tag, "zz_m_l")(
    2.997 / 2.198
)  # PHYSICAL REVIEW D 93, 074505 (2016) zz_m_l * m_l => m_l in MSbar scheme 3 GeV
set_param(job_tag, "zz_m_h")(
    81.64 / 60.62 * 0.9628
)  # PHYSICAL REVIEW D 93, 074505 (2016) zz_m_h * m_h => m_h in MSbar scheme 3 GeV
set_param(job_tag, "zz_ss_l")(1 / get_param(job_tag, "zz_m_l"))
set_param(job_tag, "zz_ss_h")(1 / get_param(job_tag, "zz_m_h"))
set_param(job_tag, "traj_list")(
    [
        2880,
        3600,
        1520,
        2000,
        3120,
        2800,
        3280,
        1440,
        2400,
        3360,
        3200,
        2560,
        2160,
        1600,
        1200,
        1280,
        1360,
        2080,
        3040,
        2320,
        2960,
        2240,
        2720,
        1840,
        2480,
        1920,
        3440,
        1680,
        2640,
        1760,
        3520,
    ]
)

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
set_param(job_tag, "meson_jwjj_threshold")(0.1)
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

def parse_args():
    """
    Parse command-line arguments using standard argparse.\n
    Returns:
        argparse.Namespace with the following attributes:
            job_tag_list: list of job tag strings to process
            no_inversion: bool, skip inversion if True
            no_contract: bool, skip contraction if True
            random_permute_job_tag_traj_list: bool, randomly permute job/traj list
    """
    parser = argparse.ArgumentParser(
        description="Auto-contractor lattice QCD data generation script."
    )
    parser.add_argument(
        "--job_tag_list",
        type=str,
        default="test-4nt8-checker",
        help="Comma-separated list of job tags to process (default: 'test-4nt8-checker')",
    )
    parser.add_argument(
        "--no_inversion",
        action="store_true",
        default=False,
        help="Skip the inversion step",
    )
    parser.add_argument(
        "--no_contract",
        action="store_true",
        default=False,
        help="Skip the contraction step",
    )
    parser.add_argument(
        "--random_permute_job_tag_traj_list",
        action="store_true",
        default=False,
        help="Randomly permute the job_tag/traj list before processing",
    )
    args, _ = parser.parse_known_args()
    args.job_tag_list = args.job_tag_list.split(",")
    return args

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
    sys.exit(0)

def try_gracefully_finish():
    """
    Call `gracefully_finish` if not test and if some work is done (q.obtained_lock_history_list != [])
    """
    if (not is_test()) and (len(q.obtained_lock_history_list) > 0):
        gracefully_finish()

if __name__ == "__main__":
    sys_args = parse_args()
    job_tag_list = sys_args.job_tag_list
    #
    if "64I-pq2" in job_tag_list:
        assert job_tag_list == [
            "64I-pq2",
        ]
        qlat_scripts.v1.load_data.is_source_specification_include_inv_type = True
    #
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
    if sys_args.random_permute_job_tag_traj_list:
        if not is_test():
            job_tag_traj_list = q.random_permute(
                job_tag_traj_list, q.RngState(f"{q.get_time()}")
            )
            job_tag_traj_list = q.get_comm().bcast(job_tag_traj_list)
    for job_tag, traj in job_tag_traj_list:
        if not sys_args.no_inversion:
            q.check_time_limit()
            run_job_inversion(job_tag, traj)
            q.clean_cache()
            try_gracefully_finish()
    for job_tag, traj in job_tag_traj_list:
        if not sys_args.no_contract:
            q.check_time_limit()
            run_job_contract(job_tag, traj)
            q.clean_cache()
            try_gracefully_finish()
    #
    gracefully_finish()

# ----
