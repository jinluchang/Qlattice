#!/usr/bin/env python3

# Need --mpi X.X.X.X --mpi X.X.X runtime option

import qlat as q
import gpt as g
import qlat_gpt as qg
import rbc_ukqcd as ru
import rbc_ukqcd_params as rup
import pprint

import os

from auto_contractor.eval import *
from auto_contractor.operators import *

import jobs
from jobs import *

from cexpr import *

jobs.save_path_default = "results"

jobs.load_path_list = [
        "results",
        "../mk-gf-gt/results",
        "../mk-selected-data/results",
        "../mk-psel-self-loop/results",
        "/sdcc/u/jluchang/qcdqedta/hlbl-data-with-cache",
        "/sdcc/u/jluchang/qcdqedta/luchang/data-gen/fill-wsnk-prop/results"
        ]

@q.timer_verbose
def check_job(job_tag, traj):
    # return True if config is finished or unavailable
    fns_produce = [
            f"auto-contractor-psel/{job_tag}/traj={traj}/checkpoint.txt",
            ]
    is_job_done = True
    for fn in fns_produce:
        if get_load_path(fn) is None:
            q.displayln_info(f"check_job: {job_tag} {traj} to do as file '{fn}' does not exist.")
            is_job_done = False
            break
    if is_job_done:
        return True
    #
    fns_need = [
            # f"configs/{job_tag}/ckpoint_lat.{traj}",
            f"point-selection/{job_tag}/traj={traj}.txt",
            f"gauge-transform/{job_tag}/traj={traj}.field",
            f"psel-prop-rand-u1/{job_tag}/traj={traj}/checkpoint ; type=1.txt",
            f"psel-prop-rand-u1/{job_tag}/traj={traj}/checkpoint ; type=2.txt",
            f"psel-prop-wsrc-light/{job_tag}/traj={traj}/checkpoint.txt",
            f"psel-prop-wsrc-light/{job_tag}/traj={traj}/checkpoint ; wsnk.txt",
            f"psel-prop-wsrc-strange/{job_tag}/traj={traj}/checkpoint.txt",
            f"psel-prop-wsrc-strange/{job_tag}/traj={traj}/checkpoint ; wsnk.txt",
            f"psel-prop-psrc-light/{job_tag}/traj={traj}/checkpoint.txt",
            f"psel-prop-psrc-strange/{job_tag}/traj={traj}/checkpoint.txt",
            f"wall-src-info-light/{job_tag}/traj={traj}.txt",
            f"wall-src-info-strange/{job_tag}/traj={traj}.txt",
            f"point-src-info/{job_tag}/traj={traj}.txt",
            ]
    for fn in fns_need:
        if get_load_path(fn) is None:
            q.displayln_info(f"check_job: {job_tag} {traj} unavailable as {fn} does not exist.")
            return True
    #
    q.check_stop()
    q.check_time_limit()
    #
    return False

@q.timer
def load_prop_wsrc_all(job_tag, traj, flavor, wi, psel, gt):
    cache = q.mk_cache(f"prop_cache", "{job_tag}", "{traj}")
    total_site = ru.get_total_site(job_tag)
    psel_ts = q.get_psel_tslice(total_site)
    if flavor in ["l", "u", "d",]:
        flavor_inv_type = 0
        flavor_tag = "light"
    elif flavor in ["s",]:
        flavor_inv_type = 1
        flavor_tag = "strange"
    else:
        assert False
    path_sp = f"psel-prop-wsrc-{flavor_tag}/{job_tag}/traj={traj}"
    gt_inv = gt.inv()
    for idx, tslice, inv_type, inv_acc in wi:
        if inv_type != flavor_inv_type:
            continue
        q.displayln_info(f"load_prop_wsrc_all: idx={idx} tslice={tslice} inv_type={inv_type} path_sp={path_sp}")
        tag = f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc}"
        # load psel psnk prop
        fn_sp = os.path.join(path_sp, f"{tag}.lat")
        sp_prop = q.PselProp(psel)
        sp_prop.load(get_load_path(fn_sp))
        sp_prop = gt_inv * sp_prop
        # convert to GPT/Grid prop mspincolor order
        cache[f"{tag} ; psel"] = q.convert_mspincolor_from_wm(sp_prop)
        # load wsnk prop
        fn_spw = os.path.join(path_sp, f"{tag} ; wsnk.lat")
        spw_prop = q.PselProp(psel_ts)
        spw_prop.load(get_load_path(fn_spw))
        # convert to GPT/Grid prop mspincolor order
        cache[f"{tag} ; psel ; wsnk"] = q.convert_mspincolor_from_wm(spw_prop)

@q.timer
def load_prop_psrc_all(job_tag, traj, flavor, pi, psel):
    cache = q.mk_cache(f"prop_cache", "{job_tag}", "{traj}")
    total_site = ru.get_total_site(job_tag)
    psel_ts = q.get_psel_tslice(total_site)
    if flavor in ["l", "u", "d",]:
        flavor_inv_type = 0
        flavor_tag = "light"
    elif flavor in ["s",]:
        flavor_inv_type = 1
        flavor_tag = "strange"
    else:
        assert False
    path_sp = f"psel-prop-psrc-{flavor_tag}/{job_tag}/traj={traj}"
    for idx, xg, inv_type, inv_acc in pi:
        if inv_type != flavor_inv_type:
            continue
        q.displayln_info(f"load_prop_psrc_all: idx={idx} xg={xg} inv_type={inv_type} path_sp={path_sp}")
        tag = f"xg=({xg[0]},{xg[1]},{xg[2]},{xg[3]}) ; type={inv_type} ; accuracy={inv_acc}"
        # load psel psnk prop
        fn_sp = os.path.join(path_sp, f"{tag}.lat")
        sp_prop = q.PselProp(psel)
        sp_prop.load(get_load_path(fn_sp))
        # convert to GPT/Grid prop mspincolor order
        cache[f"{tag} ; psel"] = q.convert_mspincolor_from_wm(sp_prop)

@q.timer
def load_prop_rand_u1_all(job_tag, traj, flavor, psel):
    cache = q.mk_cache(f"prop_cache", "{job_tag}", "{traj}")
    total_site = ru.get_total_site(job_tag)
    psel_ts = q.get_psel_tslice(total_site)
    if flavor in [ "l", "u", "d", ]:
        flavor_inv_type = 0
        flavor_tag = "light"
    elif flavor in [ "s", ]:
        flavor_inv_type = 1
        flavor_tag = "strange"
    elif flavor in [ "c", ]:
        flavor_inv_type = 2
        flavor_tag = "charm"
    else:
        assert False
    inv_type = flavor_inv_type
    inv_acc = 2
    path_sp = f"psel-prop-rand-u1/{job_tag}/traj={traj}"
    n_rand_u1 = rup.dict_params[job_tag]["n_rand_u1"]
    sp_prop_avg = q.PselProp(psel)
    q.set_zero(sp_prop_avg)
    for idx_rand_u1 in range(n_rand_u1):
        tag = f"idx_rand_u1={idx_rand_u1} ; type={inv_type} ; accuracy={inv_acc}"
        fn_sp = os.path.join(path_sp, f"{tag}.lat")
        sp_prop = q.PselProp(psel)
        sp_prop.load(get_load_path(fn_sp))
        sp_prop_avg += sp_prop
    sp_prop_avg *= 1 / n_rand_u1
    cache[f"type={inv_type} ; accuracy={inv_acc} ; rand_u1 ; psel"] = q.convert_mspincolor_from_wm(sp_prop_avg)

@q.timer
def get_prop_rand_u1_psel(prop_cache, inv_type):
    inv_acc = 2
    tag = f"type={inv_type} ; accuracy={inv_acc} ; rand_u1 ; psel"
    return prop_cache.get(tag)

@q.timer
def get_prop_psrc_psel(prop_cache, inv_type, xg_src):
    inv_acc = 0
    xg = xg_src
    tag = f"xg=({xg[0]},{xg[1]},{xg[2]},{xg[3]}) ; type={inv_type} ; accuracy={inv_acc}"
    return prop_cache.get(f"{tag} ; psel")

@q.timer
def get_prop_wsrc_psel(prop_cache, inv_type, t_src):
    inv_acc = 1
    tslice = t_src
    tag = f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc}"
    return prop_cache.get(f"{tag} ; psel")

@q.timer
def get_prop_wsnk_wsrc(prop_cache, inv_type, t_snk, t_src):
    inv_acc = 1
    tslice = t_src
    tag = f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc}"
    return prop_cache[f"{tag} ; psel ; wsnk"].get_elem(t_snk)

@q.timer
def get_prop_snk_src(prop_cache, flavor, p_snk, p_src, *, psel_pos_dict):
    # psel_pos_dict[x] == idx
    # x == tuple(pos)
    # psel.to_list()[idx] == pos
    if flavor in [ "l", "u", "d", ]:
        flavor_inv_type = 0
    elif flavor in [ "s", ]:
        flavor_inv_type = 1
    elif flavor in [ "c", ]:
        flavor_inv_type = 2
    else:
        assert False
    assert isinstance(p_snk, tuple) and isinstance(p_src, tuple)
    assert 2 == len(p_snk)
    assert 2 == len(p_src)
    type_snk, pos_snk = p_snk
    type_src, pos_src = p_src
    if type_snk == "wall" and type_src == "wall":
        msc = get_prop_wsnk_wsrc(
                prop_cache, flavor_inv_type, pos_snk, pos_src)
    elif type_snk[:5] == "point" and type_src == "wall":
        pos_snk_tuple = tuple(pos_snk)
        assert pos_snk_tuple in psel_pos_dict
        msc = get_prop_wsrc_psel(
                prop_cache, flavor_inv_type, pos_src
                ).get_elem(psel_pos_dict[pos_snk_tuple])
    elif type_snk == "wall" and type_src[:5] == "point":
        pos_src_tuple = tuple(pos_src)
        assert pos_src_tuple in psel_pos_dict
        msc = g5_herm(
                get_prop_wsrc_psel(
                    prop_cache, flavor_inv_type, pos_snk
                    ).get_elem(psel_pos_dict[pos_src_tuple]))
    elif type_snk[:5] == "point" and type_src[:5] == "point":
        pos_snk_tuple = tuple(pos_snk)
        pos_src_tuple = tuple(pos_src)
        assert pos_snk_tuple in psel_pos_dict
        assert pos_src_tuple in psel_pos_dict
        if type_src == "point":
            sp_prop = get_prop_psrc_psel(prop_cache, flavor_inv_type, pos_src)
            msc = sp_prop.get_elem(psel_pos_dict[pos_snk_tuple])
        elif type_snk == "point":
            sp_prop = get_prop_psrc_psel(prop_cache, flavor_inv_type, pos_snk)
            msc = g5_herm(sp_prop.get_elem(psel_pos_dict[pos_src_tuple]))
        elif pos_snk_tuple == pos_src_tuple and flavor in [ "c", "s", ]:
            sp_prop = get_prop_rand_u1_psel(prop_cache, flavor_inv_type)
            msc = sp_prop.get_elem(psel_pos_dict[pos_snk_tuple])
        else:
            raise Exception("get_prop_snk_src unknown p_snk={p_snk} p_src={p_src}")
    else:
        raise Exception("get_prop_snk_src unknown p_snk={p_snk} p_src={p_src}")
    return as_mspincolor(msc)

@q.timer_verbose
def mk_get_prop(job_tag, traj, get_gt, get_psel, get_pi, get_wi):
    load_prop_psrc_all(job_tag, traj, "l", get_pi(), get_psel())
    load_prop_psrc_all(job_tag, traj, "s", get_pi(), get_psel())
    load_prop_wsrc_all(job_tag, traj, "l", get_wi(), get_psel(), get_gt())
    load_prop_wsrc_all(job_tag, traj, "s", get_wi(), get_psel(), get_gt())
    load_prop_rand_u1_all(job_tag, traj, "c", get_psel())
    load_prop_rand_u1_all(job_tag, traj, "s", get_psel())
    prop_cache = q.mk_cache(f"prop_cache", "{job_tag}", "{traj}")
    psel_pos_dict = dict([ (tuple(pos), i) for i, pos in enumerate(get_psel().to_list()) ])
    def get_prop(flavor, p_snk, p_src):
        return get_prop_snk_src(prop_cache, flavor, p_snk, p_src, psel_pos_dict = psel_pos_dict)
    return get_prop

@q.timer
def get_strange_psrc_psel(pi):
    coordinate_list = []
    for idx, xg, inv_type, inv_acc in pi:
        if inv_type == 1 and inv_acc == 0:
            coordinate_list.append(xg)
    return q.PointSelection(coordinate_list)

def rel_mod(x, size):
    x = (x + 2 * size) % size
    assert x >= 0
    if 2 * x >= size:
        return x - size
    else:
        return -x

@q.timer_verbose
def auto_contractor_meson_corr_wsnk_wsrc(job_tag, traj, get_prop, get_psel, get_pi, get_wi):
    total_site = ru.get_total_site(job_tag)
    cexpr = get_cexpr_meson_corr_wsnk_wsrc()
    names_expr = get_cexpr_names(cexpr)
    names_fac = [ "rest", ]
    ld = q.mk_lat_data([
        [ "tsep", total_site[3] // 2 + 1, ],
        [ "name_fac", len(names_fac), names_fac, ],
        [ "expr_name", len(names_expr), names_expr, ],
        [ "val-err-n", 3, [ "val", "err", "n-trails", ] ],
        ])
    for tsep in range(total_site[3] // 2 + 1):
        trial_indices = []
        for t1 in range(total_site[3]):
            t2 = (t1 + tsep) % total_site[3]
            pd = {
                    "t1" : ("wall", t1,),
                    "t2" : ("wall", t2,),
                    }
            trial_indices.append(pd)
            t2 = (t1 + total_site[3] - tsep) % total_site[3]
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
                trial_indices = trial_indices,
                get_prop = get_prop,
                is_only_total = "total"
                )
        for idx_name_fac, (name_fac, results,) in enumerate(zip(names_fac, results_list)):
            for i_k, (k, v,) in enumerate(results.items()):
                ld[(tsep, idx_name_fac, i_k,)] = v + [ complex(len(trial_indices)), ]
    q.displayln_info(ld.show())
    ld.save(f"results/meson_corr/{job_tag}/traj={traj}/wsnk_wsrc.lat")

@q.timer_verbose
def auto_contractor_meson_corr_psnk_wsrc(job_tag, traj, get_prop, get_psel, get_pi, get_wi):
    total_site = ru.get_total_site(job_tag)
    vol = total_site[0] * total_site[1] * total_site[2]
    cexpr = get_cexpr_meson_corr_psnk_wsrc(vol)
    names_expr = get_cexpr_names(cexpr)
    names_fac = [ "rest", ]
    ld = q.mk_lat_data([
        [ "tsep", total_site[3] // 2 + 1, ],
        [ "name_fac", len(names_fac), names_fac, ],
        [ "expr_name", len(names_expr), names_expr, ],
        [ "val-err-n", 3, [ "val", "err", "n-trails", ] ],
        ])
    for tsep in range(total_site[3] // 2 + 1):
        trial_indices = []
        for t1 in range(total_site[3]):
            for x2 in get_psel().to_list():
                if tsep == abs(rel_mod(x2[3] - t1, total_site[3])):
                    pd = {
                            "t1" : ("wall", t1,),
                            "x2" : ("point-snk", x2,),
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
                ld[(tsep, idx_name_fac, i_k,)] = v + [ complex(len(trial_indices)), ]
    q.displayln_info(ld.show())
    ld.save(f"results/meson_corr/{job_tag}/traj={traj}/psnk_wsrc.lat")

@q.timer_verbose
def auto_contractor_meson_corr_psnk_psrc(job_tag, traj, get_prop, get_psel, get_pi, get_wi):
    total_site = ru.get_total_site(job_tag)
    vol = total_site[0] * total_site[1] * total_site[2]
    cexpr = get_cexpr_meson_corr_psnk_psrc(vol)
    names_expr = get_cexpr_names(cexpr)
    names_fac = [ "rest", ]
    ld = q.mk_lat_data([
        [ "tsep", total_site[3] // 2 + 1, ],
        [ "name_fac", len(names_fac), names_fac, ],
        [ "expr_name", len(names_expr), names_expr, ],
        [ "val-err-n", 3, [ "val", "err", "n-trails", ] ],
        ])
    for tsep in range(total_site[3] // 2 + 1):
        trial_indices = []
        for x1 in get_strange_psrc_psel(get_pi()).to_list():
            for x2 in get_psel().to_list():
                if tsep == abs(rel_mod(x2[3] - x1[3], total_site[3])):
                    pd = {
                            "x1" : ("point", x1,),
                            "x2" : ("point-snk", x2,),
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
                ld[(tsep, idx_name_fac, i_k,)] = v + [ complex(len(trial_indices)), ]
    q.displayln_info(ld.show())
    ld.save(f"results/meson_corr/{job_tag}/traj={traj}/psnk_psrc.lat")

@q.timer_verbose
def auto_contractor_vev(job_tag, traj, get_prop, get_psel, get_pi, get_wi):
    cexpr = get_cexpr_vev()
    names_expr = get_cexpr_names(cexpr)
    names_fac = [ "rest", ]
    ld = q.mk_lat_data([
        [ "name_fac", len(names_fac), names_fac, ],
        [ "expr_name", len(names_expr), names_expr, ],
        [ "val-err-n", 3, [ "val", "err", "n-trails", ] ],
        ])
    trial_indices = []
    for x in get_psel().to_list():
        pd = {
                "x" : ("point", x,),
                "x_snk" : ("point-snk", x,),
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
    for idx_name_fac, (name_fac, results,) in enumerate(zip(names_fac, results_list)):
        for i_k, (k, v,) in enumerate(results.items()):
            ld[(idx_name_fac, i_k,)] = v + [ complex(len(trial_indices)), ]
    q.displayln_info(ld.show())
    ld.save(f"results/vev/{job_tag}/traj={traj}/vev.lat")

@q.timer_verbose
def run_job(job_tag, traj):
    if check_job(job_tag, traj):
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
    get_pi = run_pi(job_tag, traj, get_psel)
    get_wi = run_wi(job_tag, traj)
    #
    get_prop = mk_get_prop(job_tag, traj, get_gt, get_psel, get_pi, get_wi)
    #
    auto_contractor_vev(job_tag, traj, get_prop, get_psel, get_pi, get_wi)
    auto_contractor_meson_corr_wsnk_wsrc(job_tag, traj, get_prop, get_psel, get_pi, get_wi)
    auto_contractor_meson_corr_psnk_wsrc(job_tag, traj, get_prop, get_psel, get_pi, get_wi)
    auto_contractor_meson_corr_psnk_psrc(job_tag, traj, get_prop, get_psel, get_pi, get_wi)
    #
    q.clean_cache()
    q.timer_display()

rup.dict_params["test-4nt8"]["n_rand_u1"] = 4
rup.dict_params["test-4nt16"]["n_rand_u1"] = 4
rup.dict_params["48I"]["n_rand_u1"] = 4
rup.dict_params["64I"]["n_rand_u1"] = 4

rup.dict_params["test-4nt8"]["trajs"] = list(range(1000, 1400, 100))
rup.dict_params["test-4nt16"]["trajs"] = list(range(1000, 1400, 100))
rup.dict_params["48I"]["trajs"] = list(range(3000, 500, -5))
rup.dict_params["64I"]["trajs"] = list(range(3000, 500, -5))

# rup.dict_params["test-4nt8"]["fermion_params"][0][2]["Ls"] = 10
# rup.dict_params["test-4nt8"]["fermion_params"][1][2]["Ls"] = 10

# rup.dict_params["test-4nt16"]["fermion_params"][0][2]["Ls"] = 10
# rup.dict_params["test-4nt16"]["fermion_params"][1][2]["Ls"] = 10

qg.begin_with_gpt()

# ADJUST ME
job_tags = [
        "test-4nt8",
        "test-4nt16",
        # "test-8nt16",
        # "test-16nt32",
        # "test-32nt64",
        # "test-48nt96",
        # "test-64nt128",
        # "test-96nt192",
        # "test-128nt256",
        # "24D",
        # "48I",
        ]

q.check_time_limit()

for job_tag in job_tags:
    q.displayln_info(pprint.pformat(rup.dict_params[job_tag]))
    for traj in rup.dict_params[job_tag]["trajs"]:
        run_job(job_tag, traj)

qg.end_with_gpt()
