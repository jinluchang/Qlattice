import qlat as q
import rbc_ukqcd_params as rup
import rbc_ukqcd as ru

from auto_contractor.eval import *

import os

from jobs import *

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
def get_prop_wsrc(prop_cache, inv_type, t_src, tag_snk_type):
    cache_type_dict = {
            "wsrc_wsnk ; psel_ts": "psel_ts",
            "wsrc ; fsel": "fsel",
            "wsrc ; psel": "psel",
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
def get_prop_wsnk_wsrc(prop_cache, inv_type, t_snk, t_src):
    sp_prop = get_prop_wsrc(prop_cache, inv_type, t_src, "wsrc_wsnk ; psel_ts")
    def f(x):
        return x.get_elem(t_snk)
    return ama_apply1(f, sp_prop)

@q.timer
def get_prop_psnk_wsrc_fsel(prop_cache, inv_type, xg_snk, t_src, fsel_pos_dict):
    assert isinstance(xg_snk, tuple) and len(xg_snk) == 4
    idx_snk = fsel_pos_dict[xg_snk]
    def f(x):
        return x.get_elem(idx_snk)
    return ama_apply1(f, get_prop_wsrc(prop_cache, inv_type, t_src, "wsrc ; fsel"))

@q.timer
def get_prop_psnk_wsrc_psel(prop_cache, inv_type, xg_snk, t_src, psel_pos_dict):
    assert isinstance(xg_snk, tuple) and len(xg_snk) == 4
    idx_snk = psel_pos_dict[xg_snk]
    def f(x):
        return x.get_elem(idx_snk)
    return ama_apply1(f, get_prop_wsrc(prop_cache, inv_type, t_src, "wsrc ; psel"))

@q.timer
def get_prop_rand_u1_fsel(prop_cache, inv_type):
    inv_acc = 2
    tag = f"type={inv_type} ; accuracy={inv_acc} ; rand_u1 ; fsel"
    return prop_cache["fsel"].get(tag)

@q.timer
def get_prop_psnk_rand_u1_fsel(prop_cache, inv_type, xg_snk, fsel_pos_dict):
    assert isinstance(xg_snk, tuple) and len(xg_snk) == 4
    idx_snk = fsel_pos_dict[xg_snk]
    def f(x):
        return x.get_elem(idx_snk)
    return ama_apply1(f, get_prop_rand_u1_fsel(prop_cache, inv_type))

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
        msc = get_prop_wsnk_wsrc(prop_cache, inv_type, pos_snk, pos_src)
    elif type_snk[:5] == "point" and type_src == "wall":
        assert isinstance(pos_src, int)
        if type_snk == "point":
            assert pos_snk in psel_pos_dict
            msc = get_prop_psnk_wsrc_psel(prop_cache, inv_type, pos_snk, pos_src, psel_pos_dict)
        else:
            assert pos_snk in fsel_pos_dict
            msc = get_prop_psnk_wsrc_fsel(prop_cache, inv_type, pos_snk, pos_src, fsel_pos_dict)
    elif type_snk == "wall" and type_src[:5] == "point":
        assert isinstance(pos_snk, int)
        if type_snk == "point":
            assert pos_src in psel_pos_dict
            msc = get_prop_psnk_wsrc_psel(prop_cache, inv_type, pos_src, pos_snk, psel_pos_dict)
        else:
            assert pos_src in fsel_pos_dict
            msc = get_prop_psnk_wsrc_fsel(prop_cache, inv_type, pos_src, pos_snk, fsel_pos_dict)
        msc = ama_apply1(g5_herm, msc)
    elif type_snk[:5] == "point" and type_src[:5] == "point":
        # type can be "point" or "point-snk"
        # ADJUST ME
        rand_u1_flavors = [ "c", "s", "l", ]
        # rand_u1_flavors = [ "c", "s", ]
        # rand_u1_flavors = [ "c", ]
        #
        if type_src == "point" and type_snk == "point":
            # means we use point source at the source location and psel psnk at sink location
            assert pos_src in psel_pos_dict
            assert pos_snk in psel_pos_dict
            assert False
        elif type_src == "point":
            # means we use point source at the source location
            assert pos_src in psel_pos_dict
            assert pos_snk in fsel_pos_dict
            assert False
        elif type_snk == "point":
            # means we use point source at the sink location
            assert pos_snk in psel_pos_dict
            assert pos_src in fsel_pos_dict
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

@q.timer
def load_prop_wsrc_all(job_tag, traj, flavor, *, wi, psel, fsel, fselc, gt):
    # cache_fsel[f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc} ; wsrc ; fsel"]
    # cache_psel[f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc} ; wsrc ; psel"]
    # cache_psel_ts[f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc} ; wsrc_wsnk ; psel_ts"]
    # cache_prob[f"type={inv_type} ; accuracy={inv_acc} ; wsrc ; prob"]
    cache_fsel = q.mk_cache(f"prop_cache", f"{job_tag}", f"{traj}", f"fsel")
    cache_psel = q.mk_cache(f"prop_cache", f"{job_tag}", f"{traj}", f"psel")
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
    path_s = f"prop-wsrc-{flavor_tag}/{job_tag}/traj={traj}/geon-info.txt"
    sfr = q.open_fields(get_load_path(path_s), "r")
    path_sp = f"psel-prop-wsrc-{flavor_tag}/{job_tag}/traj={traj}"
    gt_inv = gt.inv()
    count = { 1: 0, 2: 0, }
    for idx, tslice, inv_type, inv_acc in wi:
        if inv_type != flavor_inv_type:
            continue
        q.displayln_info(0, f"load_prop_wsrc_all: idx={idx} tslice={tslice} inv_type={inv_type} path_sp={path_sp}")
        tag = f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc}"
        # load fsel psnk prop
        sc_prop = q.SelProp(fselc)
        sc_prop.load_double_from_float(sfr, tag)
        sc_prop = gt_inv * sc_prop
        s_prop = q.SelProp(fsel)
        s_prop @= sc_prop
        # convert to GPT/Grid prop mspincolor order
        cache_fsel[f"{tag} ; wsrc ; fsel"] = q.convert_mspincolor_from_wm(s_prop)
        # load psel psnk prop
        fn_sp = os.path.join(path_sp, f"{tag}.lat")
        sp_prop = q.PselProp(psel)
        sp_prop.load(get_load_path(fn_sp))
        sp_prop = gt_inv * sp_prop
        sp_prop_diff = q.PselProp(psel)
        sp_prop_diff @= sc_prop
        sp_prop_diff -= sp_prop
        assert sp_prop_diff.qnorm() <= 1e-14 * sp_prop.qnorm()
        # convert to GPT/Grid prop mspincolor order
        cache_psel[f"{tag} ; wsrc ; psel"] = q.convert_mspincolor_from_wm(sp_prop)
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
def load_prop_psrc_all(job_tag, traj, flavor, *, psel, fsel, fselc):
    # cache_fsel[f"xg=({xg[0]},{xg[1]},{xg[2]},{xg[3]}) ; type={inv_type} ; accuracy={inv_acc} ; psrc ; fsel"]
    # cache_psel[f"xg=({xg[0]},{xg[1]},{xg[2]},{xg[3]}) ; type={inv_type} ; accuracy={inv_acc} ; psrc ; psel"]
    # cache_psel_ts[f"xg=({xg[0]},{xg[1]},{xg[2]},{xg[3]}) ; type={inv_type} ; accuracy={inv_acc} ; psrc_wsnk ; psel_ts"]
    # cache_prob[f"type={inv_type} ; accuracy={inv_acc} ; psrc ; prob"]
    cache_fsel = q.mk_cache(f"prop_cache", f"{job_tag}", f"{traj}", f"fsel")
    cache_psel = q.mk_cache(f"prop_cache", f"{job_tag}", f"{traj}", f"psel")
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
    path_s = f"prop-psrc-{flavor_tag}/{job_tag}/traj={traj}/geon-info.txt"
    sfr = q.open_fields(get_load_path(path_s), "r")
    path_sp = f"psel-prop-psrc-{flavor_tag}/{job_tag}/traj={traj}"
    count = { 0: 0, 1: 0, 2: 0, }
    inv_type = flavor_inv_type
    idx = 0
    xg_list = psel.to_list()
    for xg, inv_acc in [ (xg, inv_acc) for xg in xg_list for inv_acc in (0, 1, 2,) ]:
        xg_str = f"({xg[0]},{xg[1]},{xg[2]},{xg[3]})"
        tag = f"xg={xg_str} ; type={inv_type} ; accuracy={inv_acc}"
        if not sfr.has(tag):
            continue
        q.displayln_info(0, f"load_prop_psrc_all: idx={idx} ; {tag} ; path_sp={path_sp}")
        idx += 1
        # load fsel psnk prop
        sc_prop = q.SelProp(fselc)
        sc_prop.load_double_from_float(sfr, tag)
        s_prop = q.SelProp(fsel)
        s_prop @= sc_prop
        # convert to GPT/Grid prop mspincolor order
        cache_fsel[f"{tag} ; psrc ; fsel"] = q.convert_mspincolor_from_wm(s_prop)
        # load psel psnk prop
        fn_sp = os.path.join(path_sp, f"{tag}.lat")
        sp_prop = q.PselProp(psel)
        sp_prop.load(get_load_path(fn_sp))
        sp_prop_diff = q.PselProp(psel)
        sp_prop_diff @= sc_prop
        sp_prop_diff -= sp_prop
        assert sp_prop_diff.qnorm() <= 1e-14 * sp_prop.qnorm()
        # convert to GPT/Grid prop mspincolor order
        cache_psel[f"{tag} ; psrc ; psel"] = q.convert_mspincolor_from_wm(sp_prop)
        # load wsnk prop
        fn_spw = os.path.join(path_sp, f"{tag} ; wsnk.lat")
        spw_prop = q.PselProp(psel_ts)
        spw_prop.load(get_load_path(fn_spw))
        # convert to GPT/Grid prop mspincolor order
        cache_psel_ts[f"{tag} ; psrc_wsnk ; psel_ts"] = q.convert_mspincolor_from_wm(spw_prop)
        count[inv_acc] += 1
    sfr.close()
    cache_prob[f"type={flavor_inv_type} ; accuracy=0 ; psrc ; prob"] = count[0] / len(xg_list)
    cache_prob[f"type={flavor_inv_type} ; accuracy=1 ; psrc ; prob"] = rup.dict_params[job_tag]["prob_acc_1_psrc"]
    cache_prob[f"type={flavor_inv_type} ; accuracy=2 ; psrc ; prob"] = rup.dict_params[job_tag]["prob_acc_2_psrc"]

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
    path_s = f"prop-rand-u1-{flavor_tag}/{job_tag}/traj={traj}/geon-info.txt"
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

@q.timer_verbose
def run_get_prop(job_tag, traj, *, get_gt, get_psel, get_fsel, get_psel_smear, get_wi):
    @q.timer_verbose
    def mk_get_prop():
        wi = get_wi()
        gt = get_gt()
        psel = get_psel()
        psel_smear = get_psel_smear()
        fsel, fselc = get_fsel()
        load_prop_wsrc_all(job_tag, traj, "l", wi = wi, psel = psel, fsel = fsel, fselc = fselc, gt = gt)
        load_prop_wsrc_all(job_tag, traj, "s", wi = wi, psel = psel, fsel = fsel, fselc = fselc, gt = gt)
        load_prop_psrc_all(job_tag, traj, "l", psel = psel, fsel = fsel, fselc = fselc)
        load_prop_psrc_all(job_tag, traj, "s", psel = psel, fsel = fsel, fselc = fselc)
        load_prop_rand_u1_all(job_tag, traj, "l", fsel = fsel)
        load_prop_rand_u1_all(job_tag, traj, "s", fsel = fsel)
        load_prop_rand_u1_all(job_tag, traj, "c", fsel = fsel)
        prop_cache = q.mk_cache(f"prop_cache", f"{job_tag}", f"{traj}")
        psel_pos_dict = dict([ (tuple(pos), i) for i, pos in enumerate(psel.to_list()) ])
        fsel_pos_dict = dict([ (tuple(pos), i) for i, pos in enumerate(fsel.to_psel_local().to_list()) ])
        def get_prop(flavor, p_snk, p_src):
            return get_prop_snk_src(prop_cache, flavor, p_snk, p_src, psel_pos_dict = psel_pos_dict, fsel_pos_dict = fsel_pos_dict)
        return get_prop
    return q.lazy_call(mk_get_prop)
