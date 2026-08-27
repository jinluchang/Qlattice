#!/usr/bin/env python3

import argparse

import qlat_gpt as qg
import qlat as q

from qlat_scripts.v1 import (
    load_path_list,
    get_load_path,
    get_param,
    set_param,
    is_test,
    check_job,
    run_params,
    run_gf,
    run_gt,
    run_wi,
    run_eig,
    run_eig_strange,
    run_prop_wsrc_full,
    run_f_weight_from_wsrc_prop_full,
    run_f_rand_01,
    run_fsel_prob,
    run_psel_prob,
    run_fsel_from_fsel_prob,
    run_psel_from_psel_prob,
    run_prop_wsrc_sparse,
    run_prop_psrc,
    run_prop_rand_u1,
    run_gf_ape,
    run_psel_split,
    run_fsel_split,
    run_field_rand_u1_dict,
    run_psel_smear,
    run_psel_smear_median,
    run_prop_smear,
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
]

### ------

@q.timer_verbose
def run_job_inversion(job_tag, traj):
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

### ------

@q.timer_verbose
def run_job_load_data_demo(job_tag, traj):
    """
    Demonstrate how to load data generated by run_job_inversion.
    Shows the actual qlat loading calls without using run_* wrappers.
    Lists content of each prop type, loads one example, and verifies numpy view.
    """
    import os
    import numpy as np
    #
    fname = q.get_fname()
    #
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    #
    # Load gauge transform directly
    fn_gt = f"{job_tag}/gauge-transform/traj-{traj}.field"
    path_gt = get_load_path(fn_gt)
    assert path_gt is not None, f"{fname}: gt file not found: {fn_gt}"
    gt = q.GaugeTransform()
    gt.load(path_gt)
    gt_inv = gt.inv()
    q.json_results_append(f"{fname}: loaded gauge transform from {path_gt}")
    #
    # Load point selection directly
    fn_psel = f"{job_tag}/points-selection/traj-{traj}.lati"
    path_psel = get_load_path(fn_psel)
    assert path_psel is not None, f"{fname}: psel file not found: {fn_psel}"
    psel = q.PointsSelection()
    psel.load(path_psel)
    q.json_results_append(f"{fname}: loaded psel with {psel.n_points} points")
    #
    # Get xg_arr from psel and verify shape
    psel_xg_arr = psel.xg_arr
    q.json_results_append(
        f"{fname}: psel xg_arr shape = {psel_xg_arr.shape}, dtype = {psel_xg_arr.dtype}"
    )
    q.json_results_append(f"{fname}: psel n_points = {psel.n_points}")
    assert psel_xg_arr.shape[0] == psel.n_points
    assert psel_xg_arr.shape[1] == 4
    #
    # Load field selection directly and add psel to make union set
    fn_fsel = f"{job_tag}/field-selection/traj-{traj}.field"
    path_fsel = get_load_path(fn_fsel)
    assert path_fsel is not None, f"{fname}: fsel file not found: {fn_fsel}"
    fsel = q.FieldSelection()
    fsel.load(path_fsel)
    fsel.add_psel(psel)
    q.json_results_append(
        f"{fname}: loaded fsel with {fsel.n_elems} selected sites (including psel)"
    )
    #
    # Get xg_arr from fsel and verify shape
    fsel_xg_arr = np.array([xg.to_tuple() for xg in fsel], dtype=np.int32)
    q.json_results_append(
        f"{fname}: fsel xg_arr shape = {fsel_xg_arr.shape}, dtype = {fsel_xg_arr.dtype}"
    )
    q.json_results_append(f"{fname}: fsel n_elems = {fsel.n_elems}")
    assert fsel_xg_arr.shape[0] == fsel.n_elems
    assert fsel_xg_arr.shape[1] == 4
    #
    inv_type = 0
    flavor_tag = "light"
    prob_suffix = " ; fsel-prob-psrc-prop"
    #
    # List wall source psel propagator content
    path_sp = f"{job_tag}/psel-prop-wsrc-{flavor_tag}/traj-{traj}"
    wsrc_psel_entries = []
    if get_load_path(f"{path_sp}.qar", f"{path_sp}/checkpoint.txt") is not None:
        for tslice in range(total_site[3]):
            tag = f"tslice={tslice} ; type={inv_type} ; accuracy=1"
            fn_sp = os.path.join(path_sp, f"{tag}.lat")
            if get_load_path(fn_sp) is not None:
                wsrc_psel_entries.append(tag)
    q.json_results_append(f"{fname}: wsrc psel entries = {len(wsrc_psel_entries)}")
    #
    # List wall source fsel propagator content
    path_s = f"{job_tag}/prop-wsrc-{flavor_tag}/traj-{traj}"
    wsrc_fsel_all = []
    wsrc_fsel_entries = []
    if get_load_path(f"{path_s}.qar", f"{path_s}/geon-info.txt") is not None:
        sfr = q.open_fields(get_load_path(path_s + "/geon-info.txt"), "r")
        wsrc_fsel_all = [tag for tag in sfr.list() if f"type={inv_type}" in tag]
        wsrc_fsel_entries = [
            tag for tag in wsrc_fsel_all if not tag.endswith(prob_suffix)
        ]
        sfr.close()
    q.json_results_append(f"{fname}: wsrc fsel total entries = {len(wsrc_fsel_all)}")
    q.json_results_append(f"{fname}: wsrc fsel prop entries = {len(wsrc_fsel_entries)}")
    #
    # List point source psel propagator content
    path_sp_psrc = f"{job_tag}/psel-prop-psrc-{flavor_tag}/traj-{traj}"
    psrc_psel_entries = []
    if (
        get_load_path(f"{path_sp_psrc}.qar", f"{path_sp_psrc}/checkpoint.txt")
        is not None
    ):
        for xg in psel:
            xg_str = f"({xg[0]},{xg[1]},{xg[2]},{xg[3]})"
            tag = f"xg={xg_str} ; type={inv_type} ; accuracy=1"
            fn_sp = os.path.join(path_sp_psrc, f"{tag}.lat")
            if get_load_path(fn_sp) is not None:
                psrc_psel_entries.append(tag)
    q.json_results_append(f"{fname}: psrc psel entries = {len(psrc_psel_entries)}")
    #
    # List point source fsel propagator content
    # Note: psrc fsel dataset contains both propagator entries and probability entries
    # (with " ; fsel-prob-psrc-prop" suffix). Only propagator entries can be loaded as SelProp.
    path_s_psrc = f"{job_tag}/prop-psrc-{flavor_tag}/traj-{traj}"
    psrc_fsel_all = []
    psrc_fsel_entries = []
    psrc_fsel_prob_entries = []
    if get_load_path(f"{path_s_psrc}.qar", f"{path_s_psrc}/geon-info.txt") is not None:
        sfr = q.open_fields(get_load_path(path_s_psrc + "/geon-info.txt"), "r")
        psrc_fsel_all = [tag for tag in sfr.list() if f"type={inv_type}" in tag]
        psrc_fsel_entries = [
            tag for tag in psrc_fsel_all if not tag.endswith(prob_suffix)
        ]
        psrc_fsel_prob_entries = [
            tag for tag in psrc_fsel_all if tag.endswith(prob_suffix)
        ]
        sfr.close()
    q.json_results_append(f"{fname}: psrc fsel total entries = {len(psrc_fsel_all)}")
    q.json_results_append(f"{fname}: psrc fsel prop entries = {len(psrc_fsel_entries)}")
    q.json_results_append(
        f"{fname}: psrc fsel prob entries = {len(psrc_fsel_prob_entries)}"
    )
    #
    # Load one wsrc psel propagator and verify numpy view
    if wsrc_psel_entries:
        tag = wsrc_psel_entries[0]
        fn_sp = os.path.join(path_sp, f"{tag}.lat")
        sp_prop = q.PselProp(psel)
        sp_prop.load(get_load_path(fn_sp))
        sp_prop = gt_inv * sp_prop
        arr = np.asarray(sp_prop)
        q.json_results_append(
            f"{fname}: wsrc psel numpy shape = {arr.shape}, dtype = {arr.dtype}"
        )
        assert arr.shape[0] == psel_xg_arr.shape[0]
        elem = sp_prop.get_elem_wm(0)
        q.json_results_append(f"{fname}: wsrc psel sample norm", elem.qnorm(), 1e-8)
    #
    # Load one wsrc fsel propagator and verify numpy view
    if wsrc_fsel_entries:
        tag = wsrc_fsel_entries[0]
        sfr = q.open_fields(get_load_path(path_s + "/geon-info.txt"), "r")
        sc_prop = q.SelProp(fsel)
        sc_prop.load_double_from_float(sfr, tag)
        sc_prop = gt_inv * sc_prop
        sfr.close()
        arr = np.asarray(sc_prop)
        q.json_results_append(
            f"{fname}: wsrc fsel numpy shape = {arr.shape}, dtype = {arr.dtype}"
        )
        assert arr.shape[0] == fsel_xg_arr.shape[0]
        sp_check = q.PselProp(psel)
        sp_check @= sc_prop
        q.json_results_append(f"{fname}: wsrc fsel sample norm", sp_check.qnorm(), 1e-8)
    #
    # Load one psrc psel propagator and verify numpy view
    if psrc_psel_entries:
        tag = psrc_psel_entries[0]
        fn_sp = os.path.join(path_sp_psrc, f"{tag}.lat")
        sp_prop = q.PselProp(psel)
        sp_prop.load(get_load_path(fn_sp))
        arr = np.asarray(sp_prop)
        q.json_results_append(
            f"{fname}: psrc psel numpy shape = {arr.shape}, dtype = {arr.dtype}"
        )
        assert arr.shape[0] == psel_xg_arr.shape[0]
        elem = sp_prop.get_elem_wm(0)
        q.json_results_append(f"{fname}: psrc psel sample norm", elem.qnorm(), 1e-8)
    #
    # Load one psrc fsel propagator and verify numpy view
    # psrc fsel uses fsel_combine (original fsel + probabilistically selected sites)
    if psrc_fsel_entries:
        tag = psrc_fsel_entries[0]
        sfr = q.open_fields(get_load_path(path_s_psrc + "/geon-info.txt"), "r")
        sc_prop = q.SelProp(fsel)
        sc_prop.load_double_from_float(sfr, tag)
        sfr.close()
        arr = np.asarray(sc_prop)
        q.json_results_append(
            f"{fname}: psrc fsel numpy shape = {arr.shape}, dtype = {arr.dtype}"
        )
        assert arr.shape[0] == fsel_xg_arr.shape[0]
        sp_check = q.PselProp(psel)
        sp_check @= sc_prop
        q.json_results_append(f"{fname}: psrc fsel sample norm", sp_check.qnorm(), 1e-8)
    #
    q.clean_cache()
    q.json_results_append(f"{fname}: completed successfully")

### ------

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

##################### CMD options #####################

job_tag_list_default = [
    "test-4nt8-checker",
]
job_tag_list_str_default = ",".join(job_tag_list_default)

def parse_args():
    parser = argparse.ArgumentParser(
        description="Data generation demo script."
    )
    parser.add_argument(
        "--job_tag_list",
        type=str,
        default=job_tag_list_str_default,
        help="Comma-separated list of job tags",
    )
    parser.add_argument(
        "--random_order",
        action="store_true",
        default=False,
        help="Randomly permute the job_tag/traj list",
    )
    args, _ = parser.parse_known_args()
    args.job_tag_list = args.job_tag_list.split(",")
    return args

#######################################################

def gracefully_finish():
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
    sys_args = parse_args()
    job_tag_list = sys_args.job_tag_list
    is_random_order = sys_args.random_order

    qg.begin_with_gpt()
    q.check_time_limit()
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
    if is_random_order:
        job_tag_traj_list = q.random_permute(
            job_tag_traj_list, q.RngState(f"{q.get_time()}")
        )
        job_tag_traj_list = q.get_comm().bcast(job_tag_traj_list)
    for job_tag, traj in job_tag_traj_list:
        q.check_time_limit()
        run_job_inversion(job_tag, traj)
        run_job_load_data_demo(job_tag, traj)
        q.clean_cache()
        try_gracefully_finish()
    #
    gracefully_finish()

# ----
