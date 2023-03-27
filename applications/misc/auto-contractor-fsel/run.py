#!/usr/bin/env python3

import qlat as q
#import gpt as g
#import qlat_gpt as qg
#import rbc_ukqcd as ru
import rbc_ukqcd_params as rup
import pprint

import os

from auto_contractor.eval import *
from auto_contractor.operators import *

from jobs import *
from load_data import *
from params import *

from cexpr import *

get_pi = None

load_path_list[:] = [
        "results",
        "../mk-gf-gt/results",
        "../mk-sel/results",
        "../mk-fsel-self-loop/results",
        "../mk-wsrc-prop/results",
        "../qcddata",
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/mk-gf-gt/results"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/mk-sel/results"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/mk-fsel-self-loop/results"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/mk-wsrc-prop/results"),
        os.path.join(os.getenv("HOME"), "qcddata"),
        "/work/2/gu19/share/ljin/data-gen/mk-fsel-self-loop/16IH2/results",
        "/work/2/gu19/share/ljin/data-gen/mk-wsrc-prop/16IH2/results",
        "/work/2/gu19/share/ljin/data-gen/mk-fsel-self-loop/32IfineH/results",
        "/work/2/gu19/share/ljin/data-gen/mk-wsrc-prop/32IfineH/results",
        ]

def mk_ama_val(val, source_specification, val_list, rel_acc_list, prob_list):
    # source_specification need to be unique for each propagator source to ensure proper AMA correction for final result
    # e.g. source_specification = ("point", (12, 2, 3, 4,),)
    assert len(val_list) == len(prob_list)
    corrections = []
    for val_i, rel_acc_i, prob_i in zip(val_list, rel_acc_list, prob_list):
        if val_i is not None:
            corrections.append((val_i, { source_specification: (rel_acc_i, prob_i), },))
    return AmaVal(val, corrections)


@q.timer_verbose
def auto_contractor_meson_corr_wsnk_wsrc(job_tag, traj, get_prop, get_fsel, get_pi, get_wi):
    fn = f"{job_tag}/auto-contractor-fsel/traj-{traj}/meson_corr/wsnk_wsrc.lat"
    if get_load_path(fn) is not None:
        return
    total_site = rup.get_total_site(job_tag)
    cexpr = get_cexpr_meson_corr()
    names_expr = get_cexpr_names(cexpr)
    names_fac = [ "rest", ]
    ld = q.mk_lat_data([
        [ "name_fac", len(names_fac), names_fac, ],
        [ "expr_name", len(names_expr), names_expr, ],
        [ "tsep", total_site[3], ],
        [ "val-err-n", 3, [ "val", "err", "n-trails", ] ],
        ])
    for tsep in range(total_site[3]):
        trial_indices = []
        for t1 in range(total_site[3]):
            for t2 in range(total_site[3]):
                if tsep == (t2 - t1) % total_site[3]:
                    pd = {
                            "x1" : ("wall", t1,),
                            "x2" : ("wall", t2,),
                            }
                    trial_indices.append(pd)
        if len(trial_indices) == 0:
            continue
        def positions_dict_maker(idx):
            pd = idx
            facs = [ 1.0, ]
            return pd, facs
        results_list = eval_cexpr_simulation(
                cexpr,
                positions_dict_maker = positions_dict_maker,
                trial_indices = get_mpi_chunk(trial_indices),
                get_prop = get_prop,
                is_only_total = "total"
                )
        for idx_name_fac, (name_fac, results,) in enumerate(zip(names_fac, results_list)):
            for i_k, (k, v,) in enumerate(results.items()):
                ld[(idx_name_fac, i_k, tsep,)] = v + [ complex(len(trial_indices)), ]
    q.displayln_info(ld.show())
    ld.save(get_save_path(fn))

@q.timer_verbose
def auto_contractor_various_corr_wsnk_wsrc(job_tag, traj, get_prop, get_fsel, get_pi, get_wi):
    fn = f"{job_tag}/auto-contractor-fsel/traj-{traj}/corr/wsnk_wsrc.lat"
    if get_load_path(fn) is not None:
        return
    total_site = rup.get_total_site(job_tag)
    cexpr = get_cexpr_various_corr()
    names_expr = get_cexpr_names(cexpr)
    names_fac = [ "rest", ]
    ld = q.mk_lat_data([
        [ "name_fac", len(names_fac), names_fac, ],
        [ "expr_name", len(names_expr), names_expr, ],
        [ "tsep", total_site[3], ],
        [ "val-err-n", 3, [ "val", "err", "n-trails", ] ],
        ])
    for tsep in range(total_site[3]):
        trial_indices = []
        for t1 in range(total_site[3]):
            for t2 in range(total_site[3]):
                if tsep == (t2 - t1) % total_site[3]:
                    pd = {
                            "t1" : ("wall", t1,),
                            "t2" : ("wall", t2,),
                            }
                    trial_indices.append(pd)
        if len(trial_indices) == 0:
            continue
        def positions_dict_maker(idx):
            pd = idx
            facs = [ 1.0, ]
            return pd, facs
        results_list = eval_cexpr_simulation(
                cexpr,
                positions_dict_maker = positions_dict_maker,
                trial_indices = get_mpi_chunk(trial_indices),
                get_prop = get_prop,
                is_only_total = "total"
                )
        results = results_list[0]
        if q.get_id_node() == 0:
            fn = get_save_path(f"{job_tag}/auto-contractor-fsel/traj-{traj}/corr/tsep{tsep}.bin")
            with open(fn, mode='wb') as f:
                for k, v in results.items():
                    if v[1].real == 0:
                        ratio_real = None
                    else:
                        ratio_real = v[0].real / v[1].real
                    if v[1].imag == 0:
                        ratio_imag = None
                    else:
                        ratio_imag = v[0].imag / v[1].imag
                    q.displayln_info(f"{k}:\n  {v}, ({ratio_real}, {ratio_imag})")
                    ###
                    f.write(v[0].real)
                    f.write(v[0].imag)
                    f.write(v[1].real)
                    f.write(v[1].imag)
    if q.get_id_node() == 0:
        def mk_key(info):
            def f(c):
                if c in "()<>/* ":
                    return "_"
                else:
                    return c
            info = "".join(map(f, info))
            while True:
                fn = info.replace("__", "_")
                if fn == info:
                    break
                info = fn
            if fn[-1] == "_":
                fn = fn[:-1]
            return fn
        metafn = get_save_path(f"{job_tag}/auto-contractor-fsel/traj-{traj}/corr/meta.txt")
        with open(metafn, mode='w') as metaf:
            for k, v in results.items():
                key = mk_key(f"{k}")
                metaf.write(f"{key}\n")

@q.timer_verbose
def auto_contractor_meson_corr_psnk_wsrc(job_tag, traj, get_prop, get_fsel, get_pi, get_wi):
    fn = f"{job_tag}/auto-contractor-fsel/traj-{traj}/meson_corr/psnk_wsrc.lat"
    if get_load_path(fn) is not None:
        return
    total_site = rup.get_total_site(job_tag)
    cexpr = get_cexpr_meson_corr()
    names_expr = get_cexpr_names(cexpr)
    names_fac = [ "rest", ]
    ld = q.mk_lat_data([
        [ "name_fac", len(names_fac), names_fac, ],
        [ "expr_name", len(names_expr), names_expr, ],
        [ "tsep", total_site[3], ],
        [ "val-err-n", 3, [ "val", "err", "n-trails", ] ],
        ])
    fsel, fselc = get_fsel()
    for tsep in range(total_site[3]):
        trial_indices = []
        for x2 in fsel.to_psel_local().to_list():
            for t1 in range(total_site[3]):
                t2 = x2[3]
                if tsep == (t2 - t1) % total_site[3]:
                    pd = {
                            "x1" : ("wall", t1,),
                            "x2" : ("point-snk", tuple(x2),),
                            }
                    trial_indices.append(pd)
        if len(trial_indices) == 0:
            continue
        def positions_dict_maker(idx):
            pd = idx
            facs = [ 1.0, ]
            return pd, facs
        results_list = eval_cexpr_simulation(
                cexpr,
                positions_dict_maker = positions_dict_maker,
                trial_indices = trial_indices,
                get_prop = get_prop,
                is_only_total = "total"
                )
        for idx_name_fac, (name_fac, results,) in enumerate(zip(names_fac, results_list)):
            for i_k, (k, v,) in enumerate(results.items()):
                ld[(idx_name_fac, i_k, tsep,)] = v + [ complex(len(trial_indices)), ]
    q.displayln_info(ld.show())
    ld.save(get_save_path(fn))

@q.timer_verbose
def auto_contractor_vev(job_tag, traj, get_prop, get_fsel, get_pi, get_wi):
    fn = f"{job_tag}/auto-contractor-fsel/traj-{traj}/vev.lat"
    if get_load_path(fn) is not None:
        return
    cexpr = get_cexpr_vev()
    names_expr = get_cexpr_names(cexpr)
    names_fac = [ "rest", ]
    ld = q.mk_lat_data([
        [ "name_fac", len(names_fac), names_fac, ],
        [ "expr_name", len(names_expr), names_expr, ],
        [ "val-err-n", 3, [ "val", "err", "n-trails", ] ],
        ])
    trial_indices = []
    fsel, fselc = get_fsel()
    for x in fsel.to_psel_local().to_list():
        pd = {
                "x" : ("point-snk", tuple(x),),
                }
        trial_indices.append(pd)
    def positions_dict_maker(idx):
        pd = idx
        facs = [ 1.0, ]
        return pd, facs
    assert trial_indices
    results_list = eval_cexpr_simulation(
            cexpr,
            positions_dict_maker = positions_dict_maker,
            trial_indices = trial_indices,
            get_prop = get_prop,
            is_only_total = "total"
            )
    for idx_name_fac, (name_fac, results,) in enumerate(zip(names_fac, results_list)):
        for i_k, (k, v,) in enumerate(results.items()):
            ld[(idx_name_fac, i_k,)] = v + [ complex(len(trial_indices)), ]
    q.displayln_info(ld.show())
    ld.save(get_save_path(fn))

@q.timer_verbose
def auto_contractor_3f4f_matching(job_tag, traj, get_prop, get_fsel, get_pi, get_wi):
    total_site = rup.get_total_site(job_tag)
    cexpr = get_cexpr_3f4f_matching()
    names_expr = get_cexpr_names(cexpr)
    src_snk_seps = [2,4,6,8]
    tsep_src2 = -2
    tsep_snk2 = 2
    tsep_src3 = -4
    tsep_snk3 = 4
    tsep_src4 = -6
    tsep_snk4 = 6
    q.mk_dirs_info(get_save_path(f"{job_tag}/auto-contractor-fsel/traj-{traj}/3f4f_b81"))
    fsel, fselc = get_fsel()
    for tsnk_tsrc in src_snk_seps:
        max_top_tsrc = tsnk_tsrc // 2
        min_top_tsrc = tsnk_tsrc // 2
        #
        for top_tsrc in range(min_top_tsrc,max_top_tsrc+1):
            tsrc1_top = - top_tsrc
            tsrc2_top = tsep_src2  + tsrc1_top
            tsrc3_top = tsep_src3  + tsrc1_top
            tsrc4_top = tsep_src4  + tsrc1_top
            tsnk1_top = tsnk_tsrc + tsrc1_top
            tsnk2_top = tsep_snk2  + tsnk1_top
            tsnk3_top = tsep_snk3  + tsnk1_top
            tsnk4_top = tsep_snk4  + tsnk1_top
            trial_indices = []
            for x in fsel.to_psel_local().to_list():
                t2_1 = ( tsrc1_top + x[3] + total_site[3] ) % total_site[3]
                t2_2 = ( tsrc2_top + x[3] + total_site[3] ) % total_site[3]
                t2_3 = ( tsrc3_top + x[3] + total_site[3] ) % total_site[3]
                t2_4 = ( tsrc4_top + x[3] + total_site[3] ) % total_site[3]
                t1_1 = ( tsnk1_top + x[3] + total_site[3] ) % total_site[3]
                t1_2 = ( tsnk2_top + x[3] + total_site[3] ) % total_site[3]
                t1_3 = ( tsnk3_top + x[3] + total_site[3] ) % total_site[3]
                t1_4 = ( tsnk4_top + x[3] + total_site[3] ) % total_site[3]
                pd = {
                    "t1_1" : ("wall", t1_1,),
                    "t1_2" : ("wall", t1_2,),
                    "t1_3" : ("wall", t1_3,),
                    "t1_4" : ("wall", t1_4,),
                    "x" : ("point-snk", tuple(x),),
                    "t2_1" : ("wall", t2_1,),
                    "t2_2" : ("wall", t2_2,),
                    "t2_3" : ("wall", t2_3,),
                    "t2_4" : ("wall", t2_4,),
                }
                trial_indices.append(pd)
            def positions_dict_maker(idx):
                pd = idx
                facs = [ 1.0, ]
                return pd, facs
            results_list = eval_cexpr_simulation(
                cexpr,
                positions_dict_maker = positions_dict_maker,
                trial_indices = trial_indices,
                get_prop = get_prop,
                is_only_total = "total"
            )
            results = results_list[0]
            if q.get_id_node() == 0:
                fn = get_save_path(f"{job_tag}/auto-contractor-fsel/traj-{traj}/3f4f_b81/tsnk_tsrc{tsnk_tsrc}_top_tsrc{top_tsrc}.bin")
                with open(fn, mode='wb') as f:
                    for k, v in results.items():
                        if v[1].real == 0:
                            ratio_real = None
                        else:
                            ratio_real = v[0].real / v[1].real
                        if v[1].imag == 0:
                            ratio_imag = None
                        else:
                            ratio_imag = v[0].imag / v[1].imag
                        q.displayln_info(f"{k}:\n  {v}, ({ratio_real}, {ratio_imag})")
                        ###
                        f.write(v[0].real)
                        f.write(v[0].imag)
                        f.write(v[1].real)
                        f.write(v[1].imag)
    if q.get_id_node() == 0:
        def mk_key(info):
            def f(c):
                if c in "()<>/* ":
                    return "_"
                else:
                    return c
            info = "".join(map(f, info))
            while True:
                fn = info.replace("__", "_")
                if fn == info:
                    break
                info = fn
            if fn[-1] == "_":
                fn = fn[:-1]
            return fn
        metafn = get_save_path(f"{job_tag}/auto-contractor-fsel/traj-{traj}/3f4f_b81/meta.txt")
        with open(metafn, mode='w') as metaf:
            for k, v in results.items():
                key = mk_key(f"{k}")
                metaf.write(f"{key}\n")

@q.timer_verbose
def auto_contractor_3f4f_matching_tslice(job_tag, traj, get_prop, get_fsel, get_pi, get_wi):
    total_site = rup.get_total_site(job_tag)
    cexpr = get_cexpr_3f4f_matching()
    names_expr = get_cexpr_names(cexpr)
    src_snk_seps = [2,4,6,8]
    tsep_src2 = -2
    tsep_snk2 = 2
    tsep_src3 = -4
    tsep_snk3 = 4
    tsep_src4 = -6
    tsep_snk4 = 6
    q.mk_dirs_info(get_save_path(f"{job_tag}/auto-contractor-fsel/traj-{traj}/3f4f_b81"))
    fsel, fselc = get_fsel()
    for tsnk_tsrc in src_snk_seps:
        max_top_tsrc = tsnk_tsrc // 2
        min_top_tsrc = tsnk_tsrc // 2
        #
        for top_tsrc in range(min_top_tsrc,max_top_tsrc+1):
            tsrc1_top = - top_tsrc
            tsrc2_top = tsep_src2  + tsrc1_top
            tsrc3_top = tsep_src3  + tsrc1_top
            tsrc4_top = tsep_src4  + tsrc1_top
            tsnk1_top = tsnk_tsrc + tsrc1_top
            tsnk2_top = tsep_snk2  + tsnk1_top
            tsnk3_top = tsep_snk3  + tsnk1_top
            tsnk4_top = tsep_snk4  + tsnk1_top
            for x3 in range(total_site[3]):
                trial_indices = []
                for x in fsel.to_psel_local().to_list():
                    if ( x[3] != x3 ):
                        continue
                    t2_1 = ( tsrc1_top + x[3] + total_site[3] ) % total_site[3]
                    t2_2 = ( tsrc2_top + x[3] + total_site[3] ) % total_site[3]
                    t2_3 = ( tsrc3_top + x[3] + total_site[3] ) % total_site[3]
                    t2_4 = ( tsrc4_top + x[3] + total_site[3] ) % total_site[3]
                    t1_1 = ( tsnk1_top + x[3] + total_site[3] ) % total_site[3]
                    t1_2 = ( tsnk2_top + x[3] + total_site[3] ) % total_site[3]
                    t1_3 = ( tsnk3_top + x[3] + total_site[3] ) % total_site[3]
                    t1_4 = ( tsnk4_top + x[3] + total_site[3] ) % total_site[3]
                    pd = {
                        "t1_1" : ("wall", t1_1,),
                        "t1_2" : ("wall", t1_2,),
                        "t1_3" : ("wall", t1_3,),
                        "t1_4" : ("wall", t1_4,),
                        "x" : ("point-snk", tuple(x),),
                        "t2_1" : ("wall", t2_1,),
                        "t2_2" : ("wall", t2_2,),
                        "t2_3" : ("wall", t2_3,),
                        "t2_4" : ("wall", t2_4,),
                    }
                    trial_indices.append(pd)
                def positions_dict_maker(idx):
                    pd = idx
                    facs = [ 1.0, ]
                    return pd, facs
                results_list = eval_cexpr_simulation(
                    cexpr,
                    positions_dict_maker = positions_dict_maker,
                    trial_indices = trial_indices,
                    get_prop = get_prop,
                    is_only_total = "total"
                )
                results = results_list[0]
                if q.get_id_node() == 0:
                    fn = get_save_path(f"{job_tag}/auto-contractor-fsel/traj-{traj}/3f4f_b81/tsnk_tsrc{tsnk_tsrc}_top_tsrc{top_tsrc}_t{x3}.bin")
                    with open(fn, mode='wb') as f:
                        for k, v in results.items():
                            if v[1].real == 0:
                                ratio_real = None
                            else:
                                ratio_real = v[0].real / v[1].real
                            if v[1].imag == 0:
                                ratio_imag = None
                            else:
                                ratio_imag = v[0].imag / v[1].imag
                            q.displayln_info(f"{k}:\n  {v}, ({ratio_real}, {ratio_imag})")
                            ###
                            f.write(v[0].real)
                            f.write(v[0].imag)
                            f.write(v[1].real)
                            f.write(v[1].imag)
    if q.get_id_node() == 0:
        def mk_key(info):
            def f(c):
                if c in "()<>/* ":
                    return "_"
                else:
                    return c
            info = "".join(map(f, info))
            while True:
                fn = info.replace("__", "_")
                if fn == info:
                    break
                info = fn
            if fn[-1] == "_":
                fn = fn[:-1]
            return fn
        metafn = get_save_path(f"{job_tag}/auto-contractor-fsel/traj-{traj}/3f4f_b81/meta.txt")
        with open(metafn, mode='w') as metaf:
            for k, v in results.items():
                key = mk_key(f"{k}")
                metaf.write(f"{key}\n")

@q.timer_verbose
def run_job(job_tag, traj):
    fns_produce = [
            f"{job_tag}/auto-contractor-fsel/traj-{traj}/checkpoint.txt",
            ]
    fns_need = [
            #(f"{job_tag}/configs/ckpoint_lat.{traj}", f"{job_tag}/configs/ckpoint_lat.IEEE64BIG.{traj}",),
            f"{job_tag}/point-selection/traj-{traj}.txt",
            f"{job_tag}/field-selection/traj-{traj}.field",
            f"{job_tag}/gauge-transform/traj-{traj}.field",
            f"{job_tag}/wall-src-info-light/traj-{traj}.txt",
            f"{job_tag}/wall-src-info-strange/traj-{traj}.txt",
            f"{job_tag}/prop-wsrc-strange/traj-{traj}/geon-info.txt",
            f"{job_tag}/psel-prop-wsrc-strange/traj-{traj}/checkpoint.txt",
            f"{job_tag}/prop-wsrc-light/traj-{traj}/geon-info.txt",
            f"{job_tag}/psel-prop-wsrc-light/traj-{traj}/checkpoint.txt",
            f"{job_tag}/prop-rand-u1-light/traj-{traj}/geon-info.txt",
            f"{job_tag}/prop-rand-u1-strange/traj-{traj}/geon-info.txt",
            f"{job_tag}/prop-rand-u1-charm/traj-{traj}/geon-info.txt",
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
            )
    #
    fn_checkpoint = f"{job_tag}/auto-contractor-fsel/traj-{traj}/checkpoint.txt"
    if get_load_path(fn_checkpoint) is None:
        if q.obtain_lock(f"locks/{job_tag}-{traj}-auto-contractor"):
            get_prop = get_get_prop()
            #get_prop = mk_get_prop(job_tag, traj,
            #        get_gt = get_gt,
            #        get_psel = get_psel,
            #        get_fsel = get_fsel,
            #        get_pi = get_pi,
            #        get_wi = get_wi,
            #        )
            # ADJUST ME
            auto_contractor_vev(job_tag, traj, get_prop, get_fsel, get_pi, get_wi)
            auto_contractor_meson_corr_wsnk_wsrc(job_tag, traj, get_prop, get_fsel, get_pi, get_wi)
            auto_contractor_meson_corr_psnk_wsrc(job_tag, traj, get_prop, get_fsel, get_pi, get_wi)
            auto_contractor_various_corr_wsnk_wsrc(job_tag, traj, get_prop, get_fsel, get_pi, get_wi)
            auto_contractor_3f4f_matching_tslice(job_tag, traj, get_prop, get_fsel, get_pi, get_wi)
            # auto_contractor_3f4f_matching(job_tag, traj, get_prop, get_fsel, get_pi, get_wi)
            #
            q.qtouch_info(get_save_path(fn_checkpoint))
            q.release_lock()
    #
    q.clean_cache()
    q.timer_display()

def get_all_cexpr():
    benchmark_eval_cexpr(get_cexpr_vev())
    benchmark_eval_cexpr(get_cexpr_meson_corr())
    benchmark_eval_cexpr(get_cexpr_meson_corr_with_env())
    benchmark_eval_cexpr(get_cexpr_various_corr())
    benchmark_eval_cexpr(get_cexpr_3f4f_matching())

def test():
    q.qremove_all_info("locks")
    q.qremove_all_info("results")
    run_job("test-4nt8", 1000)

def rel_mod(x, size):
    x = (x + 2 * size) % size
    assert x >= 0
    if 2 * x >= size:
        return x - size
    else:
        return x

size_node_list = [
        [1, 1, 1, 1],
        [1, 1, 2, 1],
        [1, 2, 2, 1],
        [2, 2, 2, 1],
        [2, 2, 4, 1],
        [2, 4, 4, 1],
        [4, 4, 4, 1],
        [4, 4, 8, 1],
        [4, 8, 8, 1],
        [8, 8, 8, 1],
        ]

q.begin_with_mpi()

# ADJUST ME
#q.qremove_all_info("cache")
get_all_cexpr()
#test()

# ADJUST ME
job_tags = [
        "test-4nt8", "test-4nt16",
        # "64I",
        # "48I",
        # "24D", "32D", "32Dfine","24DH",
        # "16IH2",
        # "32IfineH",
        ]

q.check_time_limit()

for job_tag in job_tags:
    q.displayln_info(pprint.pformat(rup.dict_params[job_tag]))
    for traj in rup.dict_params[job_tag]["trajs"]:
        run_job(job_tag, traj)

q.end_mpi()
