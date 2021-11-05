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
        "/sdcc/u/jluchang/qcdqedta/luchang/data-gen/mk-psel-self-loop/results",
        "/sdcc/u/jluchang/qcdqedta/luchang/data-gen/fill-wsnk-prop/results",
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
    cache = q.mk_cache(f"prop_cache", f"{job_tag}", f"{traj}")
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
    count = { 1: 0, 2: 0, }
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
        count[inv_acc] += 1
    assert count[1] == total_site[3]
    cache[f"type={flavor_inv_type} ; accuracy=1 ; wsrc ; prob"] = count[1] / total_site[3]
    cache[f"type={flavor_inv_type} ; accuracy=2 ; wsrc ; prob"] = count[2] / total_site[3]

@q.timer
def load_prop_psrc_all(job_tag, traj, flavor, pi, psel):
    cache = q.mk_cache(f"prop_cache", f"{job_tag}", f"{traj}")
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
def set_prop_psrc_prob(job_tag, traj, pi, psel):
    cache = q.mk_cache(f"prop_cache", f"{job_tag}", f"{traj}")
    n_point = len(psel.to_list())
    count_l = { 0: 0, 1: 0, 2: 0, }
    count_s = { 0: 0, 1: 0, 2: 0, }
    count = [ count_l, count_s, ]
    for idx, xg, inv_type, inv_acc in pi:
        count[inv_type][inv_acc] += 1
    assert count[0][0] == n_point
    assert count[0][1] == count[1][1]
    assert count[0][2] == count[1][2]
    for inv_type in [ 0, 1, ]:
        for inv_acc in [ 0, 1, 2, ]:
            cache[f"type={inv_type} ; accuracy={inv_acc} ; psrc ; prob"] = count[inv_type][inv_acc] / n_point
    sp_prop = q.PselProp(psel)
    q.set_zero(sp_prop)
    cache[f"psrc ; psel ; zero"] = sp_prop

@q.timer
def load_prop_rand_u1_all(job_tag, traj, flavor, psel):
    cache = q.mk_cache(f"prop_cache", f"{job_tag}", f"{traj}")
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

def mk_ama_val(val, source_specification, val_list, rel_acc_list, prob_list):
    assert len(val_list) == len(prob_list)
    corrections = []
    for val_i, rel_acc_i, prob_i in zip(val_list, rel_acc_list, prob_list):
        if val_i is not None:
            corrections.append((val_i, { source_specification: (rel_acc_i, prob_i), },))
    return AmaVal(val, corrections)

@q.timer
def get_prop_psrc_psel(prop_cache, inv_type, xg_src):
    inv_acc = 0
    xg = xg_src
    def mk_tag(inv_acc):
        return f"xg=({xg[0]},{xg[1]},{xg[2]},{xg[3]}) ; type={inv_type} ; accuracy={inv_acc} ; psel"
    tag = mk_tag(inv_acc = 0)
    tag1 = mk_tag(inv_acc = 1)
    tag2 = mk_tag(inv_acc = 2)
    prob = prop_cache[f"type={inv_type} ; accuracy=0 ; psrc ; prob"]
    # level light_accuracy strange_accuracy
    # 0     inv_acc=0      zero prop
    # 1     inv_acc=0      inv_acc=0
    # 2     inv_acc=1      inv_acc=1
    # 3     inv_acc=2      inv_acc=2
    if inv_type == 0:
        assert prob == 1
        val = prop_cache.get(tag)
        if tag1 not in prop_cache:
            return val
    else:
        assert prob < 1
        val = prop_cache.get(f"psrc ; psel ; zero")
        if tag not in prop_cache:
            return val
    source_specification = ("point", tuple(xg_src),)
    val_list = [
            val,
            prop_cache.get(tag),
            prop_cache.get(tag1),
            prop_cache.get(tag2),
            ]
    rel_acc_list = [ 0, 1, 2, 3, ]
    prob_list = [
            1.0,
            prop_cache[f"type=1 ; accuracy=0 ; psrc ; prob"],
            prop_cache[f"type=1 ; accuracy=1 ; psrc ; prob"],
            prop_cache[f"type=1 ; accuracy=2 ; psrc ; prob"],
            ]
    return mk_ama_val(val, source_specification, val_list, rel_acc_list, prob_list)

@q.timer
def get_prop_wsrc(prop_cache, inv_type, t_src, tag_snk_type):
    tslice = t_src
    def mk_tag(inv_acc):
        return f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc} ; {tag_snk_type}"
    tag = mk_tag(inv_acc = 1)
    tag1 = mk_tag(inv_acc = 2)
    prob = prop_cache[f"type={inv_type} ; accuracy=1 ; wsrc ; prob"]
    # level light_accuracy strange_accuracy
    # 0     inv_acc=1      inv_acc=1
    # 3     inv_acc=2      inv_acc=2
    assert prob == 1
    val = prop_cache.get(tag)
    if tag1 not in prop_cache:
        return val
    source_specification = ("wall", t_src,)
    val_list = [ val, prop_cache.get(tag1), ]
    rel_acc_list = [ 0, 3, ]
    prob_list = [ 1, prop_cache[f"type={inv_type} ; accuracy=2 ; wsrc ; prob"], ]
    return mk_ama_val(val, source_specification, val_list, rel_acc_list, prob_list)

@q.timer
def get_prop_wsrc_psel(prop_cache, inv_type, t_src):
    return get_prop_wsrc(prop_cache, inv_type, t_src, "psel")

@q.timer
def get_prop_wsnk_wsrc(prop_cache, inv_type, t_snk, t_src):
    sp_prop = get_prop_wsrc(prop_cache, inv_type, t_src, "psel ; wsnk")
    def f(x):
        return x.get_elem(t_snk)
    return ama_apply1(f, sp_prop)

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
        sp_prop = get_prop_wsrc_psel(prop_cache, flavor_inv_type, pos_src)
        def f(x):
            return x.get_elem(psel_pos_dict[pos_snk_tuple])
        msc = ama_apply1(f, sp_prop)
    elif type_snk == "wall" and type_src[:5] == "point":
        pos_src_tuple = tuple(pos_src)
        assert pos_src_tuple in psel_pos_dict
        sp_prop = get_prop_wsrc_psel(prop_cache, flavor_inv_type, pos_snk)
        def f(x):
            return g5_herm(x.get_elem(psel_pos_dict[pos_src_tuple]))
        msc = ama_apply1(f, sp_prop)
    elif type_snk[:5] == "point" and type_src[:5] == "point":
        # type can be "point" or "point-snk"
        pos_snk_tuple = tuple(pos_snk)
        pos_src_tuple = tuple(pos_src)
        assert pos_snk_tuple in psel_pos_dict
        assert pos_src_tuple in psel_pos_dict
        # ADJUST ME
        rand_u1_flavors = [ "c", "s", ]
        # rand_u1_flavors = [ "c", ]
        #
        if type_src == "point":
            # means we use point source at the source location
            sp_prop = get_prop_psrc_psel(prop_cache, flavor_inv_type, pos_src)
            def f(x):
                return x.get_elem(psel_pos_dict[pos_snk_tuple])
            msc = ama_apply1(f, sp_prop)
        elif type_snk == "point":
            # means we use point source at the sink location
            sp_prop = get_prop_psrc_psel(prop_cache, flavor_inv_type, pos_snk)
            def f(x):
                return g5_herm(x.get_elem(psel_pos_dict[pos_src_tuple]))
            msc = ama_apply1(f, sp_prop)
        elif pos_snk_tuple == pos_src_tuple and flavor in rand_u1_flavors:
            # use the rand_u1 source
            sp_prop = get_prop_rand_u1_psel(prop_cache, flavor_inv_type)
            def f(x):
                return x.get_elem(psel_pos_dict[pos_snk_tuple])
            msc = ama_apply1(f, sp_prop)
        else:
            # if nothing else work, try use point src propagator
            sp_prop = get_prop_psrc_psel(prop_cache, flavor_inv_type, pos_src)
            def f(x):
                return x.get_elem(psel_pos_dict[pos_snk_tuple])
            msc = ama_apply1(f, sp_prop)
    else:
        raise Exception("get_prop_snk_src unknown p_snk={p_snk} p_src={p_src}")
    return ama_apply1(as_mspincolor, msc)

@q.timer_verbose
def mk_get_prop(job_tag, traj, get_gt, get_psel, get_pi, get_wi):
    load_prop_psrc_all(job_tag, traj, "l", get_pi(), get_psel())
    load_prop_psrc_all(job_tag, traj, "s", get_pi(), get_psel())
    set_prop_psrc_prob(job_tag, traj, get_pi(), get_psel())
    load_prop_wsrc_all(job_tag, traj, "l", get_wi(), get_psel(), get_gt())
    load_prop_wsrc_all(job_tag, traj, "s", get_wi(), get_psel(), get_gt())
    load_prop_rand_u1_all(job_tag, traj, "c", get_psel())
    load_prop_rand_u1_all(job_tag, traj, "s", get_psel())
    prop_cache = q.mk_cache(f"prop_cache", f"{job_tag}", f"{traj}")
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
    cexpr = get_cexpr_meson_corr()
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
            for t2 in range(total_site[3]):
                if tsep == abs(rel_mod(t2 - t1, total_site[3])):
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
                ld[(tsep, idx_name_fac, i_k,)] = v + [ complex(len(trial_indices)), ]
    q.displayln_info(ld.show())
    ld.save(get_save_path(f"auto-contractor-psel/{job_tag}/traj={traj}/meson_corr/wsnk_wsrc.lat"))

@q.timer_verbose
def auto_contractor_meson_corr_psnk_wsrc(job_tag, traj, get_prop, get_psel, get_pi, get_wi):
    total_site = ru.get_total_site(job_tag)
    cexpr = get_cexpr_meson_corr()
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
                            "x1" : ("wall", t1,),
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
                trial_indices = get_mpi_chunk(trial_indices),
                get_prop = get_prop,
                is_only_total = "total"
                )
        for idx_name_fac, (name_fac, results,) in enumerate(zip(names_fac, results_list)):
            for i_k, (k, v,) in enumerate(results.items()):
                ld[(tsep, idx_name_fac, i_k,)] = v + [ complex(len(trial_indices)), ]
    q.displayln_info(ld.show())
    ld.save(get_save_path(f"auto-contractor-psel/{job_tag}/traj={traj}/meson_corr/psnk_wsrc.lat"))

@q.timer_verbose
def auto_contractor_meson_corr_psnk_psrc(job_tag, traj, get_prop, get_psel, get_pi, get_wi):
    total_site = ru.get_total_site(job_tag)
    cexpr = get_cexpr_meson_corr()
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
        for x1 in get_psel().to_list():
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
                trial_indices = get_mpi_chunk(trial_indices),
                get_prop = get_prop,
                is_only_total = "total"
                )
        for idx_name_fac, (name_fac, results,) in enumerate(zip(names_fac, results_list)):
            for i_k, (k, v,) in enumerate(results.items()):
                ld[(tsep, idx_name_fac, i_k,)] = v + [ complex(len(trial_indices)), ]
    q.displayln_info(ld.show())
    ld.save(get_save_path(f"auto-contractor-psel/{job_tag}/traj={traj}/meson_corr/psnk_psrc.lat"))

@q.timer_verbose
def auto_contractor_meson_corr_with_env_wsnk_wsrc(job_tag, traj, get_prop, get_psel, get_pi, get_wi):
    total_site = ru.get_total_site(job_tag)
    cexpr = get_cexpr_meson_corr_with_env()
    names_expr = get_cexpr_names(cexpr)
    names_fac = [ "rest", ]
    ld = q.mk_lat_data([
        [ "tsep", total_site[3] // 2 + 1, ],
        [ "name_fac", len(names_fac), names_fac, ],
        [ "expr_name", len(names_expr), names_expr, ],
        [ "val-err-n", 3, [ "val", "err", "n-trails", ] ],
        ])
    tsep_env = 6
    for tsep in range(total_site[3] // 2 + 1):
        trial_indices = []
        for t1 in range(total_site[3]):
            for t2 in range(total_site[3]):
                tsep_r = rel_mod(t2 - t1, total_site[3])
                if tsep == abs(tsep_r):
                    if tsep_r >= 0:
                        tsep_env_r = tsep_env
                    else:
                        tsep_env_r = -tsep_env
                    pd = {
                            "x1" : ("wall", t1,),
                            "x2" : ("wall", t2,),
                            "x1p" : ("wall", (t1 - tsep_env_r) % total_site[3],),
                            "x2p" : ("wall", (t2 + tsep_env_r) % total_site[3],),
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
                ld[(tsep, idx_name_fac, i_k,)] = v + [ complex(len(trial_indices)), ]
    q.displayln_info(ld.show())
    ld.save(get_save_path(f"auto-contractor-psel/{job_tag}/traj={traj}/meson_corr_with_env/wsnk_wsrc.lat"))

@q.timer_verbose
def auto_contractor_meson_corr_with_env_psnk_wsrc(job_tag, traj, get_prop, get_psel, get_pi, get_wi):
    total_site = ru.get_total_site(job_tag)
    cexpr = get_cexpr_meson_corr_with_env()
    names_expr = get_cexpr_names(cexpr)
    names_fac = [ "rest", ]
    ld = q.mk_lat_data([
        [ "tsep", total_site[3] // 2 + 1, ],
        [ "name_fac", len(names_fac), names_fac, ],
        [ "expr_name", len(names_expr), names_expr, ],
        [ "val-err-n", 3, [ "val", "err", "n-trails", ] ],
        ])
    tsep_env = 6
    for tsep in range(total_site[3] // 2 + 1):
        trial_indices = []
        for t1 in range(total_site[3]):
            for x2 in get_psel().to_list():
                t2 = x2[3]
                tsep_r = rel_mod(t2 - t1, total_site[3])
                if tsep == abs(tsep_r):
                    if tsep_r >= 0:
                        tsep_env_r = tsep_env
                    else:
                        tsep_env_r = -tsep_env
                    pd = {
                            "x1" : ("wall", t1,),
                            "x2" : ("point-snk", x2,),
                            "x1p" : ("wall", (t1 - tsep_env_r) % total_site[3],),
                            "x2p" : ("wall", (t2 + tsep_env_r) % total_site[3],),
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
                ld[(tsep, idx_name_fac, i_k,)] = v + [ complex(len(trial_indices)), ]
    q.displayln_info(ld.show())
    ld.save(get_save_path(f"auto-contractor-psel/{job_tag}/traj={traj}/meson_corr_with_env/psnk_wsrc.lat"))

@q.timer_verbose
def auto_contractor_meson_corr_with_env_psnk_psrc(job_tag, traj, get_prop, get_psel, get_pi, get_wi):
    total_site = ru.get_total_site(job_tag)
    cexpr = get_cexpr_meson_corr_with_env()
    names_expr = get_cexpr_names(cexpr)
    names_fac = [ "rest", ]
    ld = q.mk_lat_data([
        [ "tsep", total_site[3] // 2 + 1, ],
        [ "name_fac", len(names_fac), names_fac, ],
        [ "expr_name", len(names_expr), names_expr, ],
        [ "val-err-n", 3, [ "val", "err", "n-trails", ] ],
        ])
    tsep_env = 6
    for tsep in range(total_site[3] // 2 + 1):
        trial_indices = []
        for x1 in get_psel().to_list():
            t1 = x1[3]
            for x2 in get_psel().to_list():
                t2 = x2[3]
                tsep_r = rel_mod(t2 - t1, total_site[3])
                if tsep == abs(tsep_r):
                    if tsep_r >= 0:
                        tsep_env_r = tsep_env
                    else:
                        tsep_env_r = -tsep_env
                    pd = {
                            "x1" : ("point", x1,),
                            "x2" : ("point-snk", x2,),
                            "x1p" : ("wall", (t1 - tsep_env_r) % total_site[3],),
                            "x2p" : ("wall", (t2 + tsep_env_r) % total_site[3],),
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
                ld[(tsep, idx_name_fac, i_k,)] = v + [ complex(len(trial_indices)), ]
    q.displayln_info(ld.show())
    ld.save(get_save_path(f"auto-contractor-psel/{job_tag}/traj={traj}/meson_corr_with_env/psnk_psrc.lat"))

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
                "x" : ("point-snk", x,),
                }
        trial_indices.append(pd)
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
            ld[(idx_name_fac, i_k,)] = v + [ complex(len(trial_indices)), ]
    q.displayln_info(ld.show())
    ld.save(get_save_path(f"auto-contractor-psel/{job_tag}/traj={traj}/vev.lat"))

@q.timer_verbose
def auto_contractor_3f4f_matching(job_tag, traj, get_prop, get_psel, get_pi, get_wi):
    total_site = ru.get_total_site(job_tag)
    cexpr = get_cexpr_3f4f_matching()
    names_expr = get_cexpr_names(cexpr)
    src_snk_seps = [8,10,12,14,16]
    tsep_src = -4
    tsep_snk = 4
    q.mk_dirs_info(get_save_path(f"auto-contractor-psel/{job_tag}/traj={traj}/3f4f_b81"))
    for tsnk_tsrc in src_snk_seps:
        max_top_tsrc = tsnk_tsrc // 2
        min_top_tsrc = tsnk_tsrc // 2
        #
        for top_tsrc in range(min_top_tsrc,max_top_tsrc+1):
            tsrc1_top = - top_tsrc
            tsrc2_top = tsep_src  + tsrc1_top
            tsnk1_top = tsnk_tsrc + tsrc1_top
            tsnk2_top = tsep_snk  + tsnk1_top
            trial_indices = []
            for x in get_psel().to_list():
                t2_1 = ( tsrc1_top + x[3] + total_site[3] ) % total_site[3]
                t2_2 = ( tsrc2_top + x[3] + total_site[3] ) % total_site[3]
                t1_1 = ( tsnk1_top + x[3] + total_site[3] ) % total_site[3]
                t1_2 = ( tsnk2_top + x[3] + total_site[3] ) % total_site[3]
                pd = {
                    "t1_1" : ("wall", t1_1,),
                    "t1_2" : ("wall", t1_2,),
                    "x" : ("point-snk", x,),
                    "t2_1" : ("wall", t2_1,),
                    "t2_2" : ("wall", t2_2,),
                }
                trial_indices.append(pd)
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
                fn = get_save_path(f"auto-contractor-psel/{job_tag}/traj={traj}/3f4f_b81/tsnk_tsrc{tsnk_tsrc}_top_tsrc{top_tsrc}.bin")
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
        metafn = get_save_path(f"auto-contractor-psel/{job_tag}/traj={traj}/3f4f_b81/meta.txt")
        with open(metafn, mode='w') as metaf:
            for k, v in results.items():
                key = mk_key(f"{k}")
                metaf.write(f"{key}\n")

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
    fn_checkpoint = f"auto-contractor-psel/{job_tag}/traj={traj}/checkpoint.txt"
    if get_load_path(fn_checkpoint) is None:
        if q.obtain_lock(f"locks/{job_tag}-{traj}-auto-contractor"):
            get_prop = mk_get_prop(job_tag, traj, get_gt, get_psel, get_pi, get_wi)
            # ADJUST ME
            # auto_contractor_vev(job_tag, traj, get_prop, get_psel, get_pi, get_wi)
            # auto_contractor_meson_corr_wsnk_wsrc(job_tag, traj, get_prop, get_psel, get_pi, get_wi)
            # auto_contractor_meson_corr_psnk_wsrc(job_tag, traj, get_prop, get_psel, get_pi, get_wi)
            # auto_contractor_meson_corr_psnk_psrc(job_tag, traj, get_prop, get_psel, get_pi, get_wi)
            # auto_contractor_meson_corr_with_env_wsnk_wsrc(job_tag, traj, get_prop, get_psel, get_pi, get_wi)
            # auto_contractor_meson_corr_with_env_psnk_wsrc(job_tag, traj, get_prop, get_psel, get_pi, get_wi)
            # auto_contractor_meson_corr_with_env_psnk_psrc(job_tag, traj, get_prop, get_psel, get_pi, get_wi)
            auto_contractor_3f4f_matching(job_tag, traj, get_prop, get_psel, get_pi, get_wi)
            #
            # q.qtouch_info(get_save_path(fn_checkpoint))
            # q.release_lock()
    #
    q.clean_cache()
    q.timer_display()

rup.dict_params["test-4nt8"]["n_rand_u1"] = 4
rup.dict_params["test-4nt16"]["n_rand_u1"] = 4

rup.dict_params["test-4nt8"]["trajs"] = list(range(1000, 1400, 100))
rup.dict_params["test-4nt16"]["trajs"] = list(range(1000, 1400, 100))

rup.dict_params["48I"]["n_rand_u1"] = 2
rup.dict_params["64I"]["n_rand_u1"] = 2

for inv_acc in [ 0, 1, 2, ]:
    rup.dict_params["64I"]["fermion_params"][0][inv_acc]["mass"] = 0.0006203
    rup.dict_params["64I"]["fermion_params"][1][inv_acc]["mass"] = 0.02539

for inv_acc in [ 0, 1, 2, ]:
    rup.dict_params["48I"]["fermion_params"][0][inv_acc]["mass"] = 0.0006979
    rup.dict_params["48I"]["fermion_params"][1][inv_acc]["mass"] = 0.03580

rup.dict_params["48I"]["trajs"] = list(range(3000, 500, -5))
rup.dict_params["64I"]["trajs"] = list(range(3000, 500, -5))

rup.dict_params["24D"]["trajs"] = list(range(3000, 500, -5))

rup.dict_params["test-4nt8"]["fermion_params"][0][2]["Ls"] = 10
rup.dict_params["test-4nt8"]["fermion_params"][1][2]["Ls"] = 10
rup.dict_params["test-4nt8"]["fermion_params"][2][2]["Ls"] = 10

# rup.dict_params["test-4nt16"]["fermion_params"][0][2]["Ls"] = 10
# rup.dict_params["test-4nt16"]["fermion_params"][1][2]["Ls"] = 10
# rup.dict_params["test-4nt16"]["fermion_params"][2][2]["Ls"] = 10

qg.begin_with_gpt()

# ADJUST ME
job_tags = [
        "test-4nt8", "test-4nt16",
        # "64I",
        # "48I",
        # "24D",
        ]

q.check_time_limit()

for job_tag in job_tags:
    q.displayln_info(pprint.pformat(rup.dict_params[job_tag]))
    for traj in rup.dict_params[job_tag]["trajs"]:
        run_job(job_tag, traj)

qg.end_with_gpt()
