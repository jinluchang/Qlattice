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

from jobs import *

from cexpr import *

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

@q.timer
def load_prop_wsrc_all(job_tag, traj, flavor, *, wi, fsel, fselc, gt):
    # cache_fsel[f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc} ; wsrc ; fsel"]
    # cache_psel_ts[f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc} ; wsrc_wsnk ; psel_ts"]
    # cache_prob[f"type={inv_type} ; accuracy={inv_acc} ; wsrc ; prob"]
    cache_fsel = q.mk_cache(f"prop_cache", f"{job_tag}", f"{traj}", f"fsel")
    cache_psel_ts = q.mk_cache(f"prop_cache", f"{job_tag}", f"{traj}", f"psel_ts")
    cache_prob = q.mk_cache(f"prop_cache", f"{job_tag}", f"{traj}", f"prob")
    total_site = ru.get_total_site(job_tag)
    psel_ts = q.get_psel_tslice(total_site)
    if flavor in [ "l", "u", "d", ]:
        flavor_inv_type = 0
        flavor_tag = "light"
    elif flavor in [ "s", ]:
        flavor_inv_type = 1
        flavor_tag = "strange"
    else:
        assert False
    path_s = f"prop-wsrc-{flavor_tag}/{job_tag}/traj={traj}"
    sfr = q.open_fields(get_load_path(path_s), "r")
    path_sp = f"psel-prop-wsrc-{flavor_tag}/{job_tag}/traj={traj}"
    gt_inv = gt.inv()
    count = { 1: 0, 2: 0, }
    for idx, tslice, inv_type, inv_acc in wi:
        if inv_type != flavor_inv_type:
            continue
        q.displayln_info(f"load_prop_wsrc_all: idx={idx} tslice={tslice} inv_type={inv_type} path_sp={path_sp}")
        tag = f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc}"
        # load fsel psnk prop
        sc_prop = q.SelProp(fselc)
        sc_prop.load_double_from_float(sfr, tag)
        s_prop = q.SelProp(fsel)
        s_prop @= sc_prop
        s_prop = gt_inv * s_prop
        # convert to GPT/Grid prop mspincolor order
        cache_fsel[f"{tag} ; wsrc ; fsel"] = q.convert_mspincolor_from_wm(s_prop)
        # load wsnk prop
        fn_spw = os.path.join(path_sp, f"{tag} ; wsnk.lat")
        spw_prop = q.PselProp(psel_ts)
        spw_prop.load(get_load_path(fn_spw))
        # convert to GPT/Grid prop mspincolor order
        cache_psel_ts[f"{tag} ; wsrc_wsnk ; psel_ts"] = q.convert_mspincolor_from_wm(spw_prop)
        count[inv_acc] += 1
    sfr.close()
    assert count[1] == total_site[3]
    cache_prob[f"type={flavor_inv_type} ; accuracy=1 ; wsrc ; prob"] = 1
    cache_prob[f"type={flavor_inv_type} ; accuracy=2 ; wsrc ; prob"] = get_prob_exact_wsrc(job_tag)

@q.timer
def load_prop_rand_u1_all(job_tag, traj, flavor, *, fsel):
    # cache_fsel[f"type={inv_type} ; accuracy={inv_acc} ; rand_u1 ; fsel"]
    cache_fsel = q.mk_cache(f"prop_cache", f"{job_tag}", f"{traj}", f"fsel")
    total_site = ru.get_total_site(job_tag)
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
    s_prop_avg = q.SelProp(fsel)
    q.set_zero(s_prop_avg)
    path_s = f"prop-rand-u1-{flavor_tag}/{job_tag}/traj={traj}"
    sfr = q.open_fields(get_load_path(path_s), "r")
    tags = sfr.list()
    prob1 = rup.dict_params[job_tag]["prob_acc_1_rand_u1"]
    prob2 = rup.dict_params[job_tag]["prob_acc_2_rand_u1"]
    n_rand_u1_fsel = rup.dict_params[job_tag]["n_rand_u1_fsel"]
    def load(idx_rand_u1, inv_acc):
        tag = f"idx_rand_u1={idx_rand_u1} ; type={inv_type} ; accuracy={inv_acc}"
        if tag not in tags:
            return None
        s_prop = q.SelProp(fsel)
        total_bytes = s_prop.load_double_from_float(sfr, tag)
        assert total_bytes > 0
        return s_prop
    for idx_rand_u1 in range(n_rand_u1_fsel):
        sp0 = load(idx_rand_u1, inv_acc = 0)
        assert sp0 is not None
        sp1 = load(idx_rand_u1, inv_acc = 1)
        sp2 = load(idx_rand_u1, inv_acc = 2)
        if sp2 is not None:
            assert sp1 is not None
            sp2 -= sp1
            sp2 *= 1 / prob2
        if sp1 is not None:
            sp1 -= sp0
            sp1 *= 1 / prob1
        if sp1 is not None:
            sp0 += sp1
        if sp2 is not None:
            sp0 += sp2
        s_prop_avg += sp0
    s_prop_avg *= 1 / n_rand_u1_fsel
    inv_acc = 2
    cache_fsel[f"type={inv_type} ; accuracy={inv_acc} ; rand_u1 ; fsel"] = q.convert_mspincolor_from_wm(s_prop_avg)

def mk_ama_val(val, source_specification, val_list, rel_acc_list, prob_list):
    # source_specification need to be unique for each propagator source to ensure proper AMA correction for final result
    # e.g. source_specification = ("point", (12, 2, 3, 4,),)
    assert len(val_list) == len(prob_list)
    corrections = []
    for val_i, rel_acc_i, prob_i in zip(val_list, rel_acc_list, prob_list):
        if val_i is not None:
            corrections.append((val_i, { source_specification: (rel_acc_i, prob_i), },))
    return AmaVal(val, corrections)

@q.timer
def get_prop_rand_u1_fsel(prop_cache, inv_type):
    inv_acc = 2
    tag = f"type={inv_type} ; accuracy={inv_acc} ; rand_u1 ; fsel"
    return prop_cache["fsel"].get(tag)

@q.timer
def get_prop_wsrc(prop_cache, inv_type, t_src, tag_snk_type):
    cache_type_dict = {
            "wsrc_wsnk ; psel_ts": "psel_ts",
            "wsrc ; fsel": "fsel",
            }
    cache_type = cache_type_dict[tag_snk_type]
    tslice = t_src
    def mk_tag(inv_acc):
        return f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc} ; {tag_snk_type}"
    tag = mk_tag(inv_acc = 1)
    tag1 = mk_tag(inv_acc = 2)
    prob = prop_cache["prob"][f"type={inv_type} ; accuracy=1 ; wsrc ; prob"]
    # level light_accuracy strange_accuracy
    # 0     inv_acc=1      inv_acc=1
    # 3     inv_acc=2      inv_acc=2
    assert prob == 1
    val = prop_cache[cache_type].get(tag)
    if tag1 not in prop_cache[cache_type]:
        return val
    source_specification = ("wall", t_src,)
    val_list = [ val, prop_cache[cache_type].get(tag1), ]
    rel_acc_list = [ 0, 3, ]
    prob_list = [ 1, prop_cache["prob"][f"type={inv_type} ; accuracy=2 ; wsrc ; prob"], ]
    return mk_ama_val(val, source_specification, val_list, rel_acc_list, prob_list)

@q.timer
def get_prop_wsrc_fsel(prop_cache, inv_type, t_src):
    return get_prop_wsrc(prop_cache, inv_type, t_src, "wsrc ; fsel")

@q.timer
def get_prop_wsnk_wsrc(prop_cache, inv_type, t_snk, t_src):
    sp_prop = get_prop_wsrc(prop_cache, inv_type, t_src, "wsrc_wsnk ; psel_ts")
    def f(x):
        return x.get_elem(t_snk)
    return ama_apply1(f, sp_prop)

@q.timer
def get_prop_psnk_rand_u1_fsel(prop_cache, inv_type, xg_snk, fsel_pos_dict):
    assert isinstance(xg_snk, tuple) and len(xg_snk) == 4
    idx_snk = fsel_pos_dict[xg_snk]
    def f(x):
        return x.get_elem(idx_snk)
    return ama_apply1(f, get_prop_rand_u1_fsel(prop_cache, inv_type))

@q.timer
def get_prop_psnk_wsrc_fsel(prop_cache, inv_type, xg_snk, t_src, fsel_pos_dict):
    assert isinstance(xg_snk, tuple) and len(xg_snk) == 4
    idx_snk = fsel_pos_dict[xg_snk]
    def f(x):
        return x.get_elem(idx_snk)
    return ama_apply1(f, get_prop_wsrc_fsel(prop_cache, inv_type, t_src))

@q.timer
def get_prop_snk_src(prop_cache, flavor, p_snk, p_src, *, psel_pos_dict, fsel_pos_dict):
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
    inv_type = flavor_inv_type
    assert isinstance(p_snk, tuple) and isinstance(p_src, tuple)
    assert 2 == len(p_snk)
    assert 2 == len(p_src)
    type_snk, pos_snk = p_snk
    type_src, pos_src = p_src
    if type_snk[:5] == "point":
        pos_snk = tuple(pos_snk)
    if type_src[:5] == "point":
        pos_src = tuple(pos_src)
    if type_snk == "wall" and type_src == "wall":
        assert isinstance(pos_snk, int)
        assert isinstance(pos_src, int)
        msc = get_prop_wsnk_wsrc(
                prop_cache, inv_type, pos_snk, pos_src)
    elif type_snk[:5] == "point" and type_src == "wall":
        assert pos_snk in fsel_pos_dict
        assert isinstance(pos_src, int)
        msc = get_prop_psnk_wsrc_fsel(prop_cache, inv_type, pos_snk, pos_src, fsel_pos_dict)
    elif type_snk == "wall" and type_src[:5] == "point":
        assert isinstance(pos_snk, int)
        assert pos_src in fsel_pos_dict
        msc = ama_apply1(g5_herm,
                get_prop_psnk_wsrc_fsel(prop_cache, inv_type, pos_src, pos_snk, fsel_pos_dict))
    elif type_snk[:5] == "point" and type_src[:5] == "point":
        # type can be "point" or "point-snk"
        assert pos_snk in fsel_pos_dict
        assert pos_src in fsel_pos_dict
        # ADJUST ME
        rand_u1_flavors = [ "c", "s", "l", ]
        # rand_u1_flavors = [ "c", "s", ]
        # rand_u1_flavors = [ "c", ]
        #
        if type_src == "point":
            # means we use point source at the source location
            assert False
        elif type_snk == "point":
            # means we use point source at the sink location
            assert False
        elif pos_snk == pos_src and flavor in rand_u1_flavors:
            # use the rand_u1 source
            assert pos_snk in fsel_pos_dict
            msc = get_prop_psnk_rand_u1_fsel(prop_cache, inv_type, pos_snk, fsel_pos_dict)
        else:
            # if nothing else work, try use point src propagator
            assert False
    else:
        raise Exception("get_prop_snk_src unknown p_snk={p_snk} p_src={p_src}")
    return ama_apply1(as_mspincolor, msc)

@q.timer_verbose
def mk_get_prop(job_tag, traj, *, get_gt, get_psel, get_fsel, get_pi, get_wi):
    wi = get_wi()
    gt = get_gt()
    psel = get_psel()
    fsel, fselc = get_fsel()
    load_prop_wsrc_all(job_tag, traj, "l", wi = wi, fsel = fsel, fselc = fselc, gt = gt)
    load_prop_wsrc_all(job_tag, traj, "s", wi = wi, fsel = fsel, fselc = fselc, gt = gt)
    load_prop_rand_u1_all(job_tag, traj, "l", fsel = fsel)
    load_prop_rand_u1_all(job_tag, traj, "s", fsel = fsel)
    load_prop_rand_u1_all(job_tag, traj, "c", fsel = fsel)
    prop_cache = q.mk_cache(f"prop_cache", f"{job_tag}", f"{traj}")
    psel_pos_dict = dict([ (tuple(pos), i) for i, pos in enumerate(psel.to_list()) ])
    fsel_pos_dict = dict([ (tuple(pos), i) for i, pos in enumerate(fsel.to_psel_local().to_list()) ])
    def get_prop(flavor, p_snk, p_src):
        return get_prop_snk_src(prop_cache, flavor, p_snk, p_src, psel_pos_dict = psel_pos_dict, fsel_pos_dict = fsel_pos_dict)
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
        return x

@q.timer_verbose
def auto_contractor_meson_corr_wsnk_wsrc(job_tag, traj, get_prop, get_fsel, get_pi, get_wi):
    fn = f"auto-contractor-fsel/{job_tag}/traj={traj}/meson_corr/wsnk_wsrc.lat"
    if get_load_path(fn) is not None:
        return
    total_site = ru.get_total_site(job_tag)
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
def auto_contractor_meson_corr_psnk_wsrc(job_tag, traj, get_prop, get_fsel, get_pi, get_wi):
    fn = f"auto-contractor-fsel/{job_tag}/traj={traj}/meson_corr/psnk_wsrc.lat"
    if get_load_path(fn) is not None:
        return
    total_site = ru.get_total_site(job_tag)
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
                ld[(idx_name_fac, i_k, tsep,)] = v + [ complex(len(trial_indices)), ]
    q.displayln_info(ld.show())
    ld.save(get_save_path(fn))

@q.timer_verbose
def auto_contractor_vev(job_tag, traj, get_prop, get_fsel, get_pi, get_wi):
    fn = f"auto-contractor-fsel/{job_tag}/traj={traj}/vev.lat"
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
                "x" : ("point-snk", x,),
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
    total_site = ru.get_total_site(job_tag)
    cexpr = get_cexpr_3f4f_matching()
    names_expr = get_cexpr_names(cexpr)
    src_snk_seps = [2,4,6,8]
    tsep_src2 = -2
    tsep_snk2 = 2
    tsep_src3 = -4
    tsep_snk3 = 4
    tsep_src4 = -6
    tsep_snk4 = 6
    q.mk_dirs_info(get_save_path(f"auto-contractor-fsel/{job_tag}/traj={traj}/3f4f_b81"))
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
                    "x" : ("point-snk", x,),
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
                fn = get_save_path(f"auto-contractor-fsel/{job_tag}/traj={traj}/3f4f_b81/tsnk_tsrc{tsnk_tsrc}_top_tsrc{top_tsrc}.bin")
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
        metafn = get_save_path(f"auto-contractor-fsel/{job_tag}/traj={traj}/3f4f_b81/meta.txt")
        with open(metafn, mode='w') as metaf:
            for k, v in results.items():
                key = mk_key(f"{k}")
                metaf.write(f"{key}\n")

@q.timer_verbose
def run_job(job_tag, traj):
    fns_produce = [
            f"auto-contractor-fsel/{job_tag}/traj={traj}/checkpoint.txt",
            ]
    fns_need = [
            (f"configs/{job_tag}/ckpoint_lat.{traj}", f"configs/{job_tag}/ckpoint_lat.IEEE64BIG.{traj}",),
            f"point-selection/{job_tag}/traj={traj}.txt",
            f"field-selection/{job_tag}/traj={traj}.field",
            f"gauge-transform/{job_tag}/traj={traj}.field",
            f"prop-rand-u1-light/{job_tag}/traj={traj}",
            f"prop-rand-u1-strange/{job_tag}/traj={traj}",
            f"prop-rand-u1-charm/{job_tag}/traj={traj}",
            f"wall-src-info-light/{job_tag}/traj={traj}.txt",
            f"wall-src-info-strange/{job_tag}/traj={traj}.txt",
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
    get_pi = None 
    #
    fn_checkpoint = f"auto-contractor-fsel/{job_tag}/traj={traj}/checkpoint.txt"
    if get_load_path(fn_checkpoint) is None:
        if q.obtain_lock(f"locks/{job_tag}-{traj}-auto-contractor"):
            get_prop = mk_get_prop(job_tag, traj,
                    get_gt = get_gt,
                    get_psel = get_psel,
                    get_fsel = get_fsel,
                    get_pi = get_pi,
                    get_wi = get_wi,
                    )
            # ADJUST ME
            auto_contractor_vev(job_tag, traj, get_prop, get_fsel, get_pi, get_wi)
            auto_contractor_meson_corr_wsnk_wsrc(job_tag, traj, get_prop, get_fsel, get_pi, get_wi)
            auto_contractor_meson_corr_psnk_wsrc(job_tag, traj, get_prop, get_fsel, get_pi, get_wi)
            auto_contractor_3f4f_matching(job_tag, traj, get_prop, get_fsel, get_pi, get_wi)
            #
            q.qtouch_info(get_save_path(fn_checkpoint))
            q.release_lock()
    #
    q.clean_cache()
    q.timer_display()

rup.dict_params["test-4nt8"]["n_exact_wsrc"] = 2
rup.dict_params["48I"]["n_exact_wsrc"] = 2

tag = "prob_exact_wsrc"
rup.dict_params["test-4nt16"][tag] = 1/8
rup.dict_params["16IH2"][tag] = 1/16
rup.dict_params["32IfineH"][tag] = 1/32

tag = "n_rand_u1_fsel"
rup.dict_params["test-4nt8"][tag] = 16
rup.dict_params["test-4nt16"][tag] = 16
rup.dict_params["48I"][tag] = 16
rup.dict_params["64I"][tag] = 16
rup.dict_params["16IH2"][tag] = 16
rup.dict_params["32IfineH"][tag] = 64

tag = "prob_acc_1_rand_u1"
rup.dict_params["test-4nt8"][tag] = 1/4
rup.dict_params["test-4nt16"][tag] = 1/4
rup.dict_params["16IH2"][tag] = 1/16
rup.dict_params["32IfineH"][tag] = 1/32

tag = "prob_acc_2_rand_u1"
rup.dict_params["test-4nt8"][tag] = 1/16
rup.dict_params["test-4nt16"][tag] = 1/16
rup.dict_params["16IH2"][tag] = 1/64
rup.dict_params["32IfineH"][tag] = 1/128

tag = "trajs"
rup.dict_params["test-4nt8"][tag] = list(range(1000, 1400, 100))
rup.dict_params["test-4nt16"][tag] = list(range(1000, 1400, 100))
rup.dict_params["32Dfine"][tag] = list(range(500, 3000, 10))
rup.dict_params["16IH2"][tag] = list(range(500, 10000, 50))
rup.dict_params["32IfineH"][tag] = list(range(500, 10000, 50))
rup.dict_params["48I"][tag] = list(range(3000, 500, -5))
rup.dict_params["64I"][tag] = list(range(3000, 500, -5))

qg.begin_with_gpt()

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

qg.end_with_gpt()
