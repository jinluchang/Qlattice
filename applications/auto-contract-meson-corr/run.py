#!/usr/bin/env python3

from auto_contractor.operators import *

import functools
import math
import os
import sys

from jobs import *
from load_data import *
from params import *

# ----

load_path_list[:] = [
        "results",
        "../qcddata",
        os.path.join(os.getenv("HOME"), "qcddata"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-gf-gt/results"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-sel/results"),
        # os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-fsel-self-loop/results"),
        # os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-selected-data/results"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-wsrc-prop/results"),
        # os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-psrc-prop/results"),
        # os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-smear-prop/results"),
        "/sdcc/u/jluchang/qcdqedta/luchang/data-gen/fill-wsnk-prop/results",
        ]

# ----

@q.timer
def get_cexpr_meson_corr():
    def calc_cexpr():
        exprs = [
                mk_pi_p("x_2", True)    * mk_pi_p("x_1")    + "pi     * pi   ",
                mk_k_p("x_2", True)     * mk_k_p("x_1")     + "k      * k    ",
                mk_a0_p("x_2", True)    * mk_a0_p("x_1")    + "a0     * a0   ",
                mk_kappa_p("x_2", True) * mk_kappa_p("x_1") + "kappa  * kappa",
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.collect_op()
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, f"cache/auto_contract_cexpr/meson_corr-cexpr")

@q.timer_verbose
def auto_contract_meson_corr(job_tag, traj, get_prop, get_psel, get_fsel):
    fn = f"{job_tag}/auto-contract/traj-{traj}/meson_corr.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_meson_corr()
    expr_names = get_cexpr_names(cexpr)
    total_site = rup.get_total_site(job_tag)
    fsel, fselc = get_fsel()
    xg_fsel_list = list(map(tuple, fsel.to_psel_local().to_list()))
    geo = q.Geometry(total_site, 1)
    total_volume = geo.total_volume()
    def fempty():
        values = np.zeros((len(expr_names), total_site[3],), dtype = complex)
        return q.Data([ 0, values, ])
    @q.timer
    def feval(t_snk):
        counts, values = fempty().get_val()
        counts += 1
        values = values.transpose()
        for t_src in range(total_site[3]):
            t = (t_snk - t_src) % total_site[3]
            pd = {
                    "x_2" : ("wall", t_snk,),
                    "x_1" : ("wall", t_src,),
                    }
            values[t] += eval_cexpr(cexpr, positions_dict = pd, get_prop = get_prop, is_only_total = "total")
        values = values.transpose() # res[expr_name, t_sep]
        return q.Data([ counts, values, ])
    t_snk_list = get_mpi_chunk(list(range(total_site[3])))
    res_count, res_sum = q.glb_sum(q.parallel_map_sum(q.get_q_mp_proc(), feval, t_snk_list, sum_start = fempty())).get_val()
    res_count *= 1.0 / total_site[3]
    res_sum *= 1.0 / total_site[3]
    assert q.qnorm(res_count - 1.0) < 1e-10
    ld = q.mk_lat_data([
        [ "expr_name", len(expr_names), expr_names, ],
        [ "t_sep", total_site[3], ],
        ])
    ld.from_numpy(res_sum)
    q.displayln_info(ld.show())
    ld.save(get_save_path(fn))

### ------

@q.timer_verbose
def run_job(job_tag, traj):
    fns_produce = [
            f"{job_tag}/auto-contract/traj-{traj}/checkpoint.txt",
            ]
    fns_need = [
            # (f"{job_tag}/configs/ckpoint_lat.{traj}", f"{job_tag}/configs/ckpoint_lat.IEEE64BIG.{traj}",),
            f"{job_tag}/point-selection/traj-{traj}.txt",
            f"{job_tag}/field-selection/traj-{traj}.field",
            f"{job_tag}/gauge-transform/traj-{traj}.field",
            f"{job_tag}/wall-src-info-light/traj-{traj}.txt",
            f"{job_tag}/wall-src-info-strange/traj-{traj}.txt",
            f"{job_tag}/psel-prop-wsrc-light/traj-{traj}/checkpoint.txt",
            f"{job_tag}/psel-prop-wsrc-strange/traj-{traj}/checkpoint.txt",
            # f"{job_tag}/psel-prop-psrc-light/traj-{traj}/checkpoint.txt",
            # f"{job_tag}/psel-prop-psrc-strange/traj-{traj}/checkpoint.txt",
            # f"{job_tag}/prop-wsrc-light/traj-{traj}/geon-info.txt",
            # f"{job_tag}/prop-wsrc-strange/traj-{traj}/geon-info.txt",
            # f"{job_tag}/prop-psrc-light/traj-{traj}/geon-info.txt",
            # f"{job_tag}/prop-psrc-strange/traj-{traj}/geon-info.txt",
            # f"{job_tag}/prop-rand-u1-light/traj-{traj}/geon-info.txt",
            # f"{job_tag}/prop-rand-u1-strange/traj-{traj}/geon-info.txt",
            # f"{job_tag}/prop-rand-u1-charm/traj-{traj}/geon-info.txt",
            ]
    if not check_job(job_tag, traj, fns_produce, fns_need):
        return
    #
    traj_gf = traj
    if job_tag[:5] == "test-":
        # ADJUST ME
        traj_gf = 1000
        #
    #
    get_gf = run_gf(job_tag, traj_gf)
    get_gt = run_gt(job_tag, traj_gf, get_gf)
    #
    get_psel = run_psel(job_tag, traj)
    get_fsel = run_fsel(job_tag, traj, get_psel)
    #
    get_wi = run_wi(job_tag, traj)
    get_psel_smear = run_psel_smear(job_tag, traj)
    #
    get_get_prop = run_get_prop(job_tag, traj,
            get_gt = get_gt,
            get_psel = get_psel,
            get_fsel = get_fsel,
            get_psel_smear = get_psel_smear,
            get_wi = get_wi,
            prop_types = [
                "wsrc psel s",
                "wsrc psel l",
                "wsrc fsel s",
                "wsrc fsel l",
                ],
            )
    #
    fn_checkpoint = f"{job_tag}/auto-contract/traj-{traj}/checkpoint.txt"
    if get_load_path(fn_checkpoint) is None:
        if q.obtain_lock(f"locks/{job_tag}-{traj}-auto-contract"):
            get_prop = get_get_prop()
            # ADJUST ME
            auto_contract_meson_corr(job_tag, traj, get_prop, get_psel, get_fsel)
            #
            q.qtouch_info(get_save_path(fn_checkpoint))
            q.release_lock()
    #
    q.clean_cache()
    q.timer_display()

def get_all_cexpr():
    benchmark_eval_cexpr(get_cexpr_meson_corr())

def test():
    q.qremove_all_info("locks")
    q.qremove_all_info("results")
    run_job("test-4nt8", 1000)
    # run_job("16IH2", 1000)

size_node_list = [
        [ 1, 1, 1, 1, ],
        [ 1, 1, 1, 2, ],
        [ 1, 1, 1, 4, ],
        [ 1, 1, 1, 8, ],
        ]

q.begin(sys.argv, size_node_list)

# ADJUST ME
q.qremove_all_info("cache")
get_all_cexpr()
test()

# ADJUST ME
job_tags = [
        # "test-4nt8", "test-4nt16",
        # "32IH1",
        # "32IH2",
        # "24IH1",
        # "24IH2",
        # "24IH3",
        # "24D",
        # "24DH",
        # "16IH2",
        # "32IfineH",
        # "32IcoarseH1",
        # "48I",
        ]

q.check_time_limit()

for job_tag in job_tags:
    q.displayln_info(pprint.pformat(rup.dict_params[job_tag]))
    for traj in rup.dict_params[job_tag]["trajs"]:
        run_job(job_tag, traj)

q.end()
